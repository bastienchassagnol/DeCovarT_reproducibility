# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/03_05_stats_plugin_covariance.R \
#   > "logs/03_05_stats_plugin_covariance_$(date +%F).log" 2>&1 &
#
# Sample moments of the reference, per time, gene set and scored type.
# Nothing is regularised here: the deconvolution builds each regularised
# covariance from these moments with build_sigma()
# (R/03_01_stats_covariance_estimation.R), at the tuning values chosen by
# scripts/03_04_cv_covariance.R or by the bulk likelihood
# (scripts/04_01_stats_deconvolution.R).
#
# Input. Raw UMI counts of the scored cells (`scored_slots.rds`, tier
# (i)), no per-cell library division, restricted to one of two gene sets
# per time, both read from `selected_<time>.csv` (tier (ii)):
#   welch  every gene gated by the permuted Welch G at BH 0.05
#          (1,417 / 3,420 / 5,821 genes at 48 / 72 / 96 h);
#   r2_05  the gated genes with a deviance share R2_D >= 0.05
#   r2_01  the gated genes with R2_D >= 0.01, between the two
#          (17 / 140 / 635 genes).
#
# Outputs, output/grn/plugin/:
#   gene_sets.csv                  time_point, gene_set, gene
#   plugin_summary.csv             one row per time x gene set x type
#   moments_<time>_<gene_set>.rds  genes, types, mu (G x J), sigma_sample
#                                  (G x G x J, unbiased sample covariance),
#                                  n_cells, r_pool (G x G, pooled
#                                  within-type correlation) (not in Git)
#
# The closed-form shrinkage and the denoised covariance of the first pass
# are in scripts/archives/denoised_covariance/.

# ==========================================================================
# SECTION 0 · Dependencies, settings, paths ----
# ==========================================================================

root <- here::here()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))
check_blas_kernel()

time_levels <- c("48h", "72h", "96h")
r2_floors <- c(r2_05 = 0.05, r2_01 = 0.01) # deviance share of the smaller gene sets

slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
out_dir <- file.path(root, "output", "grn", "plugin")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# SECTION 1 · Inputs and gene sets ----
# ==========================================================================

if (!file.exists(slots_path)) {
  stop(
    "Missing ",
    slots_path,
    ". From the repository root run: dvc pull ",
    "output/feature_selection/tier_01/scored_slots.rds.dvc"
  )
}
slots <- readRDS(slots_path)

gene_sets <- lapply(time_levels, function(tm) {
  selected <- utils::read.csv(
    file.path(selected_dir, paste0("selected_", tm, ".csv")),
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  list(
    welch = selected$gene,
    r2_05 = selected$gene[selected$deviance_r2 >= r2_floors[["r2_05"]]],
    r2_01 = selected$gene[selected$deviance_r2 >= r2_floors[["r2_01"]]]
  )
})
names(gene_sets) <- time_levels

utils::write.csv(
  do.call(
    rbind,
    lapply(time_levels, function(tm) {
      do.call(
        rbind,
        lapply(names(gene_sets[[tm]]), function(gs) {
          data.frame(
            time_point = tm,
            gene_set = gs,
            gene = gene_sets[[tm]][[gs]],
            stringsAsFactors = FALSE
          )
        })
      )
    })
  ),
  file.path(out_dir, "gene_sets.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ==========================================================================
# SECTION 2 · Moments per time and gene set ----
# ==========================================================================

summary_rows <- list()
for (tm in time_levels) {
  # One dense read of the larger set; the smaller set is a subset of it.
  panel <- read_panel_slot(
    slots,
    selected_dir,
    tm,
    genes = gene_sets[[tm]]$welch
  )
  types <- sort(unique(panel$cell_type[!is.na(panel$cell_type)]))
  for (gs in names(gene_sets[[tm]])) {
    genes <- gene_sets[[tm]][[gs]]
    started <- proc.time()[["elapsed"]]
    moments <- sample_moments(
      panel$y[, genes, drop = FALSE],
      panel$cell_type,
      types
    )
    elapsed <- proc.time()[["elapsed"]] - started
    variances <- apply(moments$sigma_sample, 3L, diag)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      time_point = tm,
      gene_set = gs,
      cell_type = types,
      cells = as.integer(moments$n_cells[types]),
      genes = length(genes),
      constant_genes = colSums(variances <= 0),
      stringsAsFactors = FALSE
    )
    message(
      tm,
      " / ",
      gs,
      ": ",
      length(types),
      " types, ",
      length(genes),
      " genes, ",
      round(elapsed, 1),
      " s"
    )
    # Uncompressed: the arrays are dense doubles, and the deconvolution
    # script reads each file once per run.
    saveRDS(
      c(
        list(time_point = tm, gene_set = gs, genes = genes, types = types),
        moments
      ),
      file.path(out_dir, paste0("moments_", tm, "_", gs, ".rds")),
      compress = FALSE
    )
    rm(moments)
    invisible(gc())
  }
}

utils::write.csv(
  do.call(rbind, summary_rows),
  file.path(out_dir, "plugin_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
writeLines(
  utils::capture.output(utils::sessionInfo()),
  file.path(out_dir, "session_info.txt")
)
