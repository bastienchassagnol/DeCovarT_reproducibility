# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_03_stats_multivariate_consistency.R \
#   > "logs/02_03_stats_multivariate_consistency_$(date +%F).log" 2>&1 &
#
# Tier (iii), multivariate bulk consistency test. For each matched time,
# the DeCovarT model selected in chapter 4 (Welch panel, raw shrinkage
# covariance at the level chosen by the bulk likelihood, free mean and
# variance scales c_i and phi_i) is fitted to each standard-Wnt bulk
# replicate, and every gene is tested on its conditional residual
# (multivariate_tail_filter()), without a gene-wise platform ratio. The
# flags are compared with the univariate n_eff variant of
# scripts/02_03_stats_cross_platform_consistency.R on the same genes.
#
# Depends on: output/grn/plugin/moments_<time>_welch.rds (03_05),
# output/deconvolution/lambda_selection.csv (04_01) and
# output/feature_selection/tier_03/cross_platform.rds (02_03, run first).

# ==========================================================================
# SECTION 0 · Dependencies, settings, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "02_03_stats_tail_bounds.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))
source(file.path(root, "R", "04_01_stats_deconvolution.R"))

level <- 0.01
shift <- 3
seed <- 20261118L
times <- c("48h", "72h", "96h")

bulk_path <- file.path(root, "data", "raw", "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz")
slots_path <- file.path(root, "output", "feature_selection", "tier_01", "scored_slots.rds")
moments_dir <- file.path(root, "output", "grn", "plugin")
selection_path <- file.path(root, "output", "deconvolution", "lambda_selection.csv")
out_dir <- file.path(root, "output", "feature_selection", "tier_03")
cross_path <- file.path(out_dir, "cross_platform.rds")
for (f in c(bulk_path, slots_path, selection_path, cross_path)) {
  if (!file.exists(f)) {
    stop("Missing ", f, ". See the dependencies in the header of this script.")
  }
}

# ==========================================================================
# SECTION 1 · Self-checks on simulated data ----
# ==========================================================================

withr::with_seed(seed, {
  # n_eff recovers n / kappa when the bulk is overdispersed by kappa.
  g <- 4000L
  v_bulk <- stats::rgamma(g, 2, 0.1)
  d <- stats::rnorm(g, sd = sqrt(v_bulk / (500 / 8)))
  n_eff <- estimate_effective_cells(d, v_bulk, rep(0, g), upper = 500)
  stopifnot(abs(n_eff / (500 / 8) - 1) < 0.1)

  # Under the model: about 2 % of genes beyond 2.33, module sums N(0, 1).
  g <- 300L
  a <- matrix(stats::rnorm(g * 5), g) / 3
  sig <- array(
    c(tcrossprod(a) + diag(g), tcrossprod(a[, 2:1]) + diag(g)),
    c(g, g, 2),
    dimnames = list(paste0("g", seq_len(g)), paste0("g", seq_len(g)), c("A", "B"))
  )
  mu <- matrix(stats::rgamma(2 * g, 2, 0.2), g, dimnames = list(dimnames(sig)[[1]], c("A", "B")))
  p <- c(A = 0.7, B = 0.3)
  s_mix <- 0.7 * sig[, , 1] + 0.3 * sig[, , 2]
  y <- 50 * drop(mu %*% p) + t(chol(4 * s_mix)) %*% matrix(stats::rnorm(g * 200), g)
  colnames(y) <- paste0("s", seq_len(200))
  sim <- multivariate_tail_filter(
    y,
    mu,
    sig,
    p_hat = matrix(p, 2, 200, dimnames = list(c("A", "B"), colnames(y))),
    c_hat = rep(50, 200),
    phi_hat = rep(4, 200),
    modules = list(m = paste0("g", 1:20)),
    robust = FALSE
  )
  stopifnot(
    abs(mean(sim$gene$flag_cond) - 2 * level) < 0.005,
    abs(stats::var(sim$module$s_cond) - 1) < 0.25,
    all(sim$gene$gain >= 1 - 1e-8)
  )
})
message("Self-checks passed.")

# ==========================================================================
# SECTION 2 · Fit the selected model and test, per time ----
# ==========================================================================

htseq <- read_bulk_htseq(bulk_path)
slots <- readRDS(slots_path)
selection <- utils::read.csv(selection_path, stringsAsFactors = FALSE)
cross <- readRDS(cross_path)
dict <- utils::read.csv(
  file.path(root, "data", "dictionaries", "mouse_gastruloid_signaling_markers.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
prepared <- readRDS(slots_path)
dict$gene <- prepared$symbol_map$mouse_symbol[
  match(dict$gene_marker, prepared$symbol_map$gene_marker)
]
modules <- lapply(split(dict$gene, dict$cell_type_excel), function(x) unique(stats::na.omit(x)))

res <- list()
for (tm in times) {
  moments <- readRDS(file.path(moments_dir, paste0("moments_", tm, "_welch.rds")))
  pick <- selection$time == tm & selection$gene_set == "welch" &
    selection$arm == "standard" & selection$family == "shrinkage" &
    selection$selected_profile %in% TRUE
  lambda <- selection$tuning[pick]
  stopifnot(length(lambda) == 1L)
  sigma <- build_sigma(moments, "shrinkage", lambda = lambda)
  genes <- moments$genes
  samples <- standard_wnt_samples(colnames(htseq$counts), tm)
  bulk <- bulk_symbol_counts(htseq, slots$shared, genes)[genes, samples, drop = FALSE]
  fits <- lapply(samples, function(s) {
    fit_decovart_once(bulk[, s], moments$mu, sigma, "lbfgsb")
  })
  names(fits) <- samples
  failed <- vapply(fits, function(f) is.null(f$p), logical(1))
  if (any(failed)) {
    stop("DeCovarT fit failed at ", tm, ": ", paste(samples[failed], collapse = ", "))
  }
  test <- multivariate_tail_filter(
    bulk,
    moments$mu,
    sigma,
    p_hat = sapply(fits, function(f) f$p),
    c_hat = vapply(fits, function(f) f$diagnostics$mean_scale, numeric(1)),
    phi_hat = vapply(fits, function(f) f$diagnostics$variance_scale, numeric(1)),
    modules = modules,
    level = level,
    shift = shift
  )
  test$scales$time_point <- tm
  test$scales$lambda <- lambda
  test$scales$loglik <- vapply(fits, function(f) f$diagnostics$loglik, numeric(1))
  uni <- cross$variants$gene_neff$tails[[tm]]
  test$scales$a_hat <- uni$plugin$scale$a_hat[match(test$scales$sample, uni$plugin$scale$sample)]
  test$scales$n_eff <- uni$plugin$scale$n_bulk[match(test$scales$sample, uni$plugin$scale$sample)]
  test$gene$time_point <- tm
  if (!is.null(test$module)) test$module$time_point <- tm
  test$composition <- data.frame(
    time_point = tm,
    sample = rep(samples, each = length(moments$types)),
    cell_type = moments$types,
    p_hat = as.numeric(sapply(fits, function(f) f$p)),
    stringsAsFactors = FALSE
  )

  # Genes flagged in at least one replicate, Welch genes only.
  flagged <- function(rows, col) unique(rows$gene[rows[[col]] %in% TRUE])
  ug <- uni$gene[uni$gene$gene %in% genes, , drop = FALSE]
  sets <- list(
    multivariate_conditional = flagged(test$gene, "flag_cond"),
    multivariate_diagonal = flagged(test$gene, "flag_diag"),
    univariate_T1_neff = flagged(ug[ug$statement == "T1", ], "flag"),
    univariate_T4_neff = flagged(ug[ug$statement == "T4", ], "flag")
  )
  pairs <- utils::combn(names(sets), 2L)
  test$overlap <- data.frame(
    time_point = tm,
    n_genes = length(genes),
    set_a = pairs[1, ],
    set_b = pairs[2, ],
    n_a = lengths(sets[pairs[1, ]]),
    n_b = lengths(sets[pairs[2, ]]),
    n_both = mapply(function(a, b) length(intersect(sets[[a]], sets[[b]])), pairs[1, ], pairs[2, ]),
    jaccard = mapply(
      function(a, b) length(intersect(sets[[a]], sets[[b]])) / max(1, length(union(sets[[a]], sets[[b]]))),
      pairs[1, ],
      pairs[2, ]
    ),
    stringsAsFactors = FALSE
  )
  res[[tm]] <- test
  message(tm, ": ", length(genes), " genes, lambda ", lambda, ", flags ",
          paste(names(sets), lengths(sets), sep = " = ", collapse = "; "))
  rm(moments, sigma)
  gc()
}

# ==========================================================================
# SECTION 3 · Tables ----
# ==========================================================================

stack <- function(el) do.call(rbind, lapply(res, function(r) r[[el]]))
scales <- stack("scales")
gene <- stack("gene")
module <- stack("module")
power <- do.call(rbind, lapply(times, function(tm) {
  g <- gene[gene$time_point == tm, ]
  m <- module[module$time_point == tm, ]
  data.frame(
    time_point = tm,
    unit = c("gene", "module"),
    n = c(nrow(g), nrow(m)),
    shift_sd = shift,
    median_gain = c(stats::median(g$gain), NA),
    q10_gain = c(stats::quantile(g$gain, 0.1), NA),
    q90_gain = c(stats::quantile(g$gain, 0.9), NA),
    mean_power_cond = c(mean(g$power_cond), mean(m$power_cond)),
    mean_power_diag = c(mean(g$power_diag), mean(m$power_diag)),
    share_cond_higher = c(mean(g$power_cond > g$power_diag + 1e-12), mean(m$power_cond > m$power_diag + 1e-12)),
    stringsAsFactors = FALSE
  )
}))
flags <- do.call(rbind, lapply(times, function(tm) {
  g <- gene[gene$time_point == tm, ]
  data.frame(
    time_point = tm,
    n_genes = length(unique(g$gene)),
    share_abs_z_cond_gt = mean(g$flag_cond),
    share_abs_z_diag_gt = mean(g$flag_diag),
    flagged_cond = length(unique(g$gene[g$flag_cond])),
    flagged_diag = length(unique(g$gene[g$flag_diag])),
    flagged_cond_modules = length(unique(module$module[module$time_point == tm & module$flag_cond])),
    flagged_diag_modules = length(unique(module$module[module$time_point == tm & module$flag_diag])),
    n_modules = length(unique(module$module[module$time_point == tm])),
    stringsAsFactors = FALSE
  )
}))
write <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, name), row.names = FALSE, fileEncoding = "UTF-8")
}
write(scales, "multivariate_scales.csv")
write(flags, "multivariate_flags.csv")
write(power, "multivariate_power.csv")
write(stack("overlap"), "multivariate_overlap.csv")
saveRDS(
  list(
    gene = gene,
    module = module,
    scales = scales,
    flags = flags,
    power = power,
    overlap = stack("overlap"),
    composition = stack("composition"),
    settings = list(level = level, shift = shift, seed = seed)
  ),
  file.path(out_dir, "multivariate_consistency.rds")
)
message("Wrote ", out_dir)
