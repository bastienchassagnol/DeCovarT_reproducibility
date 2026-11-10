# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_03_stats_cross_platform_consistency.R \
#   > "logs/02_03_stats_cross_platform_consistency_$(date +%F).log" 2>&1 &
#
# Tier (iii), two files: R/02_03_stats_influence.R (Cook's
# distance for cells and genes; the mean model inside cook_cells() is
# glmGamPoi::glm_gp()) and R/02_03_stats_tail_bounds.R (p-free range check,
# deMixSC plug-in, tail statements T1-T5, quantile band, drift across
# the three matched times, gene-wise platform ratio and transcript length).
#
# Three variants of the cross-platform model are run on the same panel;
# none removes a gene:
#   none   assumption A3, the same per-cell mean on both platforms;
#   length A3 relaxed by a ratio b_g = f(transcript length), cross-fitted;
#   gene   A3 relaxed by a free gene-wise ratio b_g, cross-fitted.
# Cross-fitted: the ratio used at one time is estimated on the bulk
# replicates of the two other times only.

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

root <- getwd()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "02_03_stats_influence.R"))
source(file.path(root, "R", "02_03_stats_tail_bounds.R"))

z_range <- 3
mad_k <- 3
level <- 0.01
c_reg <- 2
max_degree <- 10L
mean_degree <- 3
min_abs_cor <- 0.1
drift_statement <- "T4"
h_mult <- 3
r_max <- 3
ratio_n_iter <- 30L
ratio_tol <- 1e-3
ratio_pseudo <- 0.01
length_span <- 0.5
length_n_boot <- 1000L
length_seed <- 20261109L

bulk_path <- file.path(
  root,
  "data",
  "raw",
  "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz"
)
slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
rank_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_02",
  "cross_cell.rds"
)
out_dir <- file.path(root, "output", "feature_selection", "tier_03")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(bulk_path)) {
  stop(
    "Missing ",
    bulk_path,
    ". From the repository root run: dvc pull ",
    "data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz.dvc"
  )
}
if (!file.exists(slots_path)) {
  stop(
    "Missing ",
    slots_path,
    ". Run scripts/02_01_stats_quality_control.R first."
  )
}
if (!file.exists(rank_path)) {
  stop(
    "Missing ",
    rank_path,
    ". Run scripts/02_02_stats_cross_cell_consistency.R first."
  )
}

prepared <- readRDS(slots_path)
ranked <- readRDS(rank_path)
times <- prepared$time_levels
htseq <- read_bulk_htseq(bulk_path)

marker_genes <- prepared$symbol_map$mouse_symbol
marker_genes <- unique(marker_genes[!is.na(marker_genes)])
panel <- unique(c(ranked$panel, marker_genes))

dict <- utils::read.csv(
  file.path(
    root,
    "data",
    "dictionaries",
    "mouse_gastruloid_signaling_markers.csv"
  ),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
dict$gene <- prepared$symbol_map$mouse_symbol[
  match(dict$gene_marker, prepared$symbol_map$gene_marker)
]
modules <- lapply(split(dict$gene, dict$cell_type_excel), function(genes) {
  unique(stats::na.omit(genes))
})

coef_matrix <- function(composition, types) {
  samples <- unique(composition$sample)
  out <- matrix(
    NA_real_,
    nrow = length(types),
    ncol = length(samples),
    dimnames = list(types, samples)
  )
  for (i in seq_len(nrow(composition))) {
    out[composition$cell_type[i], composition$sample[i]] <- composition$coef[i]
  }
  out
}

# ==========================================================================
# SECTION 1 · Range, tails, Cook's distance ----
# ==========================================================================

range_flags <- list()
tail_res <- list()
outlier_cells <- list()
outlier_genes <- list()
inputs <- list()

for (tm in times) {
  slot <- prepared$slots[[tm]]
  genes <- intersect(panel, rownames(slot$counts))
  x <- slot$counts[genes, , drop = FALSE]
  ct <- droplevels(factor(slot$cell_type))
  sf <- slot$library_size / stats::median(slot$library_size)
  norm <- sweep(as.matrix(x), 2L, sf, "/")
  mu <- sapply(levels(ct), function(k) {
    rowMeans(norm[, ct == k, drop = FALSE])
  })
  sd <- sapply(levels(ct), function(k) {
    apply(norm[, ct == k, drop = FALSE], 1L, stats::sd)
  })
  samples <- standard_wnt_samples(colnames(htseq$counts), tm)
  bulk <- bulk_symbol_counts(htseq, prepared$shared, genes)
  bulk <- bulk[, samples, drop = FALSE]

  range_one <- bulk_consistency_filter(
    bulk = bulk,
    mu = mu,
    sd = sd,
    n_cells = table(ct)[levels(ct)],
    z_range = z_range,
    mad_k = mad_k
  )
  range_one$time_point <- tm
  range_flags[[tm]] <- range_one
  inputs[[tm]] <- list(
    bulk = bulk,
    x = x,
    cell_type = slot$cell_type,
    mu = mu,
    sd = sd,
    n_cells = table(ct)[levels(ct)]
  )

  tail_res[[tm]] <- bulk_tail_filter(
    bulk = bulk,
    counts = x,
    cell_type = slot$cell_type,
    modules = modules,
    level = level,
    c_reg = c_reg,
    max_degree = max_degree,
    mean_degree = mean_degree,
    min_abs_cor = min_abs_cor
  )

  outlier_cells[[tm]] <- cook_cells(
    x,
    slot$cell_type,
    size_factors = sf
  )
  plug <- tail_res[[tm]]$plugin
  outlier_genes[[tm]] <- cook_genes(
    bulk,
    tail_res[[tm]]$moments$mu,
    coef = coef_matrix(plug$composition, colnames(tail_res[[tm]]$moments$mu)),
    weights = plug$weights,
    h_mult = h_mult,
    r_max = r_max
  )
}

# Statement order per time (small, cited by the chapter); the full gene,
# module, band, range and Cook tables are in the RDS.
overlaps <- lapply(times, function(tm) {
  tail_flag_overlap(tail_res[[tm]]$gene, band = tail_res[[tm]]$quantile_band)
})
names(overlaps) <- times
for (tm in times) {
  utils::write.csv(
    overlaps[[tm]]$order,
    file.path(out_dir, paste0("tail_order_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

# Gene-wise bulk-to-reference ratios: are the discrepancies a stable
# property of the gene (high correlation across times)?
log_ratio <- platform_log_ratio(tail_res)
utils::write.csv(
  cbind(log_ratio$summary, as.data.frame(log_ratio$cor)[log_ratio$summary$time_point, ]),
  file.path(out_dir, "platform_log_ratio.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Drift: is a gene flagged at several matched times, with one sign?
drift <- platform_drift_summary(
  lapply(tail_res, function(r) r$gene),
  statement = drift_statement
)
utils::write.csv(
  drift$summary[drift$summary$n_flags > 0L, , drop = FALSE],
  file.path(out_dir, "platform_drift.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
utils::write.csv(
  data.frame(
    setting = c(
      "z_range",
      "mad_k",
      "level",
      "c_reg",
      "max_degree",
      "mean_degree",
      "min_abs_cor",
      "drift_statement",
      "h_mult",
      "r_max",
      "ratio_n_iter",
      "ratio_tol",
      "ratio_pseudo",
      "length_span",
      "length_n_boot",
      "length_seed"
    ),
    value = c(
      z_range,
      mad_k,
      level,
      c_reg,
      max_degree,
      mean_degree,
      min_abs_cor,
      drift_statement,
      h_mult,
      r_max,
      ratio_n_iter,
      ratio_tol,
      ratio_pseudo,
      length_span,
      length_n_boot,
      length_seed
    ),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "settings.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ==========================================================================
# SECTION 2 · Transcript length: does the platform ratio follow it? ----
# ==========================================================================

# Full-length bulk (Illumina Stranded mRNA) against 10x 3' v3 counts:
# the prediction is a slope of 1 of the log2 ratio on log2 length.
kept <- prepared$shared[prepared$shared$keep, , drop = FALSE]
lengths <- gene_transcript_length(panel, kept$row_id[match(panel, kept$gene)])
log2_length <- stats::setNames(log2(lengths$mean_tx_length), lengths$gene)
length_test <- platform_length_test(
  log_ratio$log_ratio,
  log2_length,
  n_boot = length_n_boot,
  seed = length_seed
)
length_test$mapped_genes <- sum(!is.na(lengths$route))
length_test$panel_genes <- length(panel)
utils::write.csv(
  length_test,
  file.path(out_dir, "platform_length.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ==========================================================================
# SECTION 3 · Cross-fitted platform ratio, three variants ----
# ==========================================================================

bulk_list <- lapply(inputs, function(z) z$bulk)
mu_list <- lapply(tail_res, function(r) r$moments$mu)
ratios <- list(
  length = cross_fit_platform_ratio(
    bulk_list,
    mu_list,
    method = "length",
    log2_length = log2_length,
    n_iter = ratio_n_iter,
    tol = ratio_tol,
    span = length_span,
    c_reg = c_reg,
    pseudo = ratio_pseudo
  ),
  gene = cross_fit_platform_ratio(
    bulk_list,
    mu_list,
    method = "gene",
    n_iter = ratio_n_iter,
    tol = ratio_tol,
    c_reg = c_reg,
    pseudo = ratio_pseudo
  )
)

run_variant <- function(ratio) {
  rng <- list()
  tails <- list()
  for (tm in times) {
    inp <- inputs[[tm]]
    est <- ratio[[tm]]
    b <- stats::setNames(exp(est$log_b), est$gene)
    v <- stats::setNames(est$var_log_b, est$gene)
    b_rows <- b[rownames(inp$bulk)]
    b_rows[!is.finite(b_rows)] <- 1
    rng[[tm]] <- bulk_consistency_filter(
      bulk = sweep(inp$bulk, 1L, b_rows, "/"),
      mu = inp$mu,
      sd = inp$sd,
      n_cells = inp$n_cells,
      z_range = z_range,
      mad_k = mad_k
    )
    rng[[tm]]$time_point <- tm
    tails[[tm]] <- bulk_tail_filter(
      bulk = inp$bulk,
      counts = inp$x,
      cell_type = inp$cell_type,
      modules = modules,
      level = level,
      c_reg = c_reg,
      max_degree = max_degree,
      mean_degree = mean_degree,
      min_abs_cor = min_abs_cor,
      platform_ratio = b,
      platform_var = v
    )
  }
  ov <- lapply(times, function(tm) {
    tail_flag_overlap(tails[[tm]]$gene, band = tails[[tm]]$quantile_band)
  })
  names(ov) <- times
  list(
    range = rng,
    tails = tails,
    overlaps = ov,
    drift = platform_drift_summary(
      lapply(tails, function(r) r$gene),
      statement = drift_statement
    ),
    log_ratio = platform_log_ratio(tails)
  )
}
variants <- lapply(ratios, run_variant)
variants <- c(
  list(
    none = list(
      range = range_flags,
      tails = tail_res,
      overlaps = overlaps,
      drift = drift,
      log_ratio = log_ratio
    )
  ),
  variants
)

a3_shift <- do.call(
  rbind,
  lapply(names(variants), function(v) {
    a3_shift_summary(variants[[v]]$tails, variants[[v]]$range, v, level = level)
  })
)
utils::write.csv(
  a3_shift,
  file.path(out_dir, "a3_shift_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
platform_ratio_table <- do.call(
  rbind,
  lapply(names(ratios), function(m) {
    do.call(
      rbind,
      lapply(times, function(tm) {
        r <- ratios[[m]][[tm]]
        data.frame(
          method = m,
          held_out_time = tm,
          gene = r$gene,
          log2_b = r$log_b / log(2),
          se_log2_b = sqrt(r$var_log_b) / log(2),
          n_obs = r$n_obs,
          iterations = attr(r, "iterations"),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
utils::write.csv(
  platform_ratio_table,
  file.path(out_dir, "platform_ratio.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

saveRDS(
  list(
    range = range_flags,
    tails = tail_res,
    drift = drift,
    log_ratio = log_ratio,
    overlaps = overlaps,
    outlier_cells = outlier_cells,
    outlier_genes = outlier_genes,
    panel = panel,
    lengths = lengths,
    length_test = length_test,
    platform_ratio = ratios,
    variants = variants,
    a3_shift = a3_shift
  ),
  file.path(out_dir, "cross_platform.rds")
)
message("Wrote ", out_dir)
