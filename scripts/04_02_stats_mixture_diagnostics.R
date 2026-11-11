# Rscript --no-save --no-restore scripts/04_02_stats_mixture_diagnostics.R
#
# Two descriptions of the deconvolution problem at each time and gene
# set, from tables and moments already on disk (no fit):
#
# 1. Agreement between the bulk and the reference at the observed
#    composition, standard-Wnt samples: Pearson and Spearman correlation
#    of log1p(bulk) with log1p(M p_observed), and the share of genes
#    whose bulk-to-reference ratio is more than 4-fold away from the
#    median ratio of the sample.
# 2. Error of every method at 96 h after the three paraxial mesoderm
#    types (caudal mesoderm, pre-somitic mesoderm, somite
#    differentiation front) are summed into one part: RMSE over the
#    three remaining parts.
#
# Inputs:
#   output/deconvolution/{proportions_long.csv, observed_proportions.csv}
#   output/grn/plugin/moments_<time>_<gene_set>.rds (means only)
#   data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz
#   output/feature_selection/tier_01/scored_slots.rds
# Outputs, output/deconvolution/:
#   bulk_reference_agreement.csv
#   lineage_collapsed_metrics.csv

root <- getwd()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "04_01_stats_deconvolution.R"))

out_dir <- file.path(root, "output", "deconvolution")
read_table <- function(name) {
  utils::read.csv(
    file.path(out_dir, name),
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
}
write_table <- function(x, name) {
  utils::write.csv(
    x,
    file.path(out_dir, name),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}
estimates <- read_table("proportions_long.csv")
observed <- read_table("observed_proportions.csv")
htseq <- read_bulk_htseq(
  file.path(root, "data", "raw", "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz")
)
design <- bulk_sample_design(colnames(htseq$counts))
slots <- readRDS(
  file.path(root, "output", "feature_selection", "tier_01", "scored_slots.rds")
)

# 1 · Bulk against the reference at the observed composition ----------------
agreement <- list()
for (tm in c("48h", "72h", "96h")) {
  for (gs in c("r2_05", "welch")) {
    moments <- readRDS(
      file.path(
        root,
        "output",
        "grn",
        "plugin",
        paste0("moments_", tm, "_", gs, ".rds")
      )
    )
    mu <- moments$mu
    rm(moments)
    obs <- observed[observed$time == tm, ]
    expected <- drop(mu %*% obs$observed[match(colnames(mu), obs$cell_type)])
    bulk <- bulk_symbol_counts(htseq, slots$shared, rownames(mu))
    ids <- design$sample_id[design$time == tm & design$arm == "standard"]
    for (id in ids) {
      y <- as.numeric(bulk[rownames(mu), id])
      ok <- y > 0 & expected > 0
      log_ratio <- log2(y[ok] / expected[ok])
      agreement[[length(agreement) + 1L]] <- data.frame(
        time = tm,
        gene_set = gs,
        genes = nrow(mu),
        sample_id = id,
        pearson_log = stats::cor(log1p(y), log1p(expected)),
        spearman = stats::cor(y, expected, method = "spearman"),
        share_over_4_fold = mean(
          abs(log_ratio - stats::median(log_ratio)) > 2
        ),
        genes_zero_in_bulk = sum(y == 0),
        stringsAsFactors = FALSE
      )
    }
  }
}
write_table(do.call(rbind, agreement), "bulk_reference_agreement.csv")

# 2 · Paraxial mesoderm summed into one part, 96 h --------------------------
lineage <- c(
  "Caudal mesoderm",
  "Pre-somitic mesoderm",
  "Somite differentiation front"
)
collapse <- function(cell_type) {
  ifelse(cell_type %in% lineage, "Paraxial mesoderm", cell_type)
}
late <- estimates[
  estimates$time == "96h" & estimates$arm == "standard",
  ,
  drop = FALSE
]
late$part <- collapse(late$cell_type)
obs_late <- observed[observed$time == "96h", ]
obs_part <- tapply(obs_late$observed, collapse(obs_late$cell_type), sum)
summed <- stats::aggregate(
  proportion ~ method + sample_id + part,
  data = late,
  FUN = sum
)
summed$observed <- as.numeric(obs_part[summed$part])
full <- late
full$observed <- obs_late$observed[match(full$cell_type, obs_late$cell_type)]
rmse <- function(d) sqrt(mean((d$proportion - d$observed)^2))
by_fit <- function(d) {
  do.call(
    rbind,
    lapply(split(d, list(d$method, d$sample_id), drop = TRUE), function(x) {
      data.frame(
        method = x$method[[1L]],
        sample_id = x$sample_id[[1L]],
        rmse = rmse(x),
        stringsAsFactors = FALSE
      )
    })
  )
}
collapsed <- merge(
  by_fit(full),
  by_fit(summed),
  by = c("method", "sample_id"),
  suffixes = c("_five_parts", "_three_parts")
)
within <- late[late$cell_type %in% lineage, ]
share <- do.call(
  rbind,
  lapply(split(within, list(within$method, within$sample_id), drop = TRUE), function(x) {
    data.frame(
      method = x$method[[1L]],
      sample_id = x$sample_id[[1L]],
      psm_share_of_lineage = x$proportion[x$cell_type == "Pre-somitic mesoderm"] /
        sum(x$proportion),
      stringsAsFactors = FALSE
    )
  })
)
collapsed <- merge(collapsed, share, by = c("method", "sample_id"))
descriptors <- unique(late[, c(
  "method",
  "family",
  "algorithm",
  "gene_set",
  "covariance",
  "platform_ratio",
  "optimiser"
)])
collapsed <- merge(collapsed, descriptors, by = "method")
collapsed$observed_psm_share <- obs_late$observed[
  obs_late$cell_type == "Pre-somitic mesoderm"
] / sum(obs_late$observed[obs_late$cell_type %in% lineage])
write_table(collapsed, "lineage_collapsed_metrics.csv")
message("Mixture diagnostics written.")
