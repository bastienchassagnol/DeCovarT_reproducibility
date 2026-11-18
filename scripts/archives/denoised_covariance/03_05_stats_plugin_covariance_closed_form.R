# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/03_05_stats_plugin_covariance.R \
#   > "logs/03_05_stats_plugin_covariance_$(date +%F).log" 2>&1 &
#
# Plug-in within-type moments for the quick deconvolution benchmark
# (scripts/04_01_stats_deconvolution.R). No cross-validation: the
# shrinkage intensity is the closed-form Schafer-Strimmer value of each
# type. The cross-validated comparison of the five arms stays in
# scripts/03_04_cv_covariance.R.
#
# Input. Raw UMI counts of the scored cells (`scored_slots.rds`, tier
# (i)), no per-cell library division, restricted to one of two gene sets
# per time, both read from `selected_<time>.csv` (tier (ii)):
#   welch  every gene gated by the permuted Welch G at BH 0.05
#          (1,417 / 3,420 / 5,821 genes at 48 / 72 / 96 h);
#   r2_05  the gated genes with a deviance share R2_D >= 0.05
#          (17 / 140 / 635 genes).
#
# For each time x gene set x scored type, two covariances with the same
# intensity lambda:
#   raw       target D on the sample covariance (corpcor::cov.shrink);
#   denoised  target D on the denoised covariance, kappa = 0: the Poisson
#             diagonal is removed and negative eigenvalues are clipped
#             (denoise_covariance()). This is Sigma_bio.
# lambda is the closed-form value of the raw counts. It is derived for
# the raw correlation; for the denoised arm it is a convenience value.
#
# Outputs, output/grn/plugin/:
#   gene_sets.csv                  time_point, gene_set, gene
#   plugin_summary.csv             one row per time x gene set x type
#   moments_<time>_<gene_set>.rds  genes, types, mu (G x J), sigma_raw and
#                                  sigma_bio (G x G x J), lambda, n_cells
#                                  (DVC, not in Git)

# ==========================================================================
# SECTION 0 · Dependencies, settings, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "03_02_stats_regularisation_path.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))
check_blas_kernel()

time_levels <- c("48h", "72h", "96h")
r2_floor <- 0.05 # deviance share of the small gene set
kappa <- 0 # all size variation treated as biological
# Target D needs lambda > 0 to make the clipped denoised matrix positive
# definite; the same floor as shrinkage_lambda_grid().
lambda_floor <- 1e-3

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
    r2_05 = selected$gene[selected$deviance_r2 >= r2_floor]
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
# SECTION 2 · Moments per time, gene set and type ----
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
    n_genes <- length(genes)
    mu <- matrix(
      NA_real_,
      n_genes,
      length(types),
      dimnames = list(genes, types)
    )
    sigma_raw <- array(
      NA_real_,
      c(n_genes, n_genes, length(types)),
      dimnames = list(genes, genes, types)
    )
    sigma_bio <- sigma_raw
    lambda <- stats::setNames(rep(NA_real_, length(types)), types)
    n_cells <- stats::setNames(rep(NA_integer_, length(types)), types)
    for (tp in types) {
      started <- proc.time()[["elapsed"]]
      X <- panel$y[which(panel$cell_type %in% tp), genes, drop = FALSE]
      lambda_star <- shrinkage_lambda_closed_form(X)
      lambda_used <- min(max(lambda_star, lambda_floor), 1)
      raw <- fit_cov_shrinkage(X, lambda_used)
      den <- fit_cov_shrinkage(
        X,
        lambda_used,
        preprocess = "denoised",
        kappa = kappa
      )
      if (isTRUE(raw$failed) || isTRUE(den$failed)) {
        stop("Shrinkage failed for ", tm, " / ", gs, " / ", tp, ".")
      }
      mu[, tp] <- raw$mu
      sigma_raw[,, tp] <- raw$sigma
      sigma_bio[,, tp] <- den$sigma
      lambda[[tp]] <- lambda_used
      n_cells[[tp]] <- nrow(X)
      elapsed <- proc.time()[["elapsed"]] - started
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        time_point = tm,
        gene_set = gs,
        cell_type = tp,
        cells = nrow(X),
        genes = n_genes,
        lambda_star = lambda_star,
        lambda_used = lambda_used,
        constant_genes = sum(apply(X, 2L, stats::var) <= 0),
        floored_genes = length(den$floored_genes),
        clipped_eigenvalues = den$n_clipped,
        elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      )
      message(
        tm,
        " / ",
        gs,
        " / ",
        tp,
        ": ",
        nrow(X),
        " cells, ",
        n_genes,
        " genes, lambda = ",
        signif(lambda_used, 3),
        ", ",
        length(den$floored_genes),
        " floored genes, ",
        round(elapsed, 1),
        " s"
      )
    }
    # Uncompressed: the arrays are dense doubles, and the deconvolution
    # script reads each file once per run.
    saveRDS(
      list(
        time_point = tm,
        gene_set = gs,
        genes = genes,
        types = types,
        mu = mu,
        sigma_raw = sigma_raw,
        sigma_bio = sigma_bio,
        lambda = lambda,
        n_cells = n_cells,
        kappa = kappa
      ),
      file.path(out_dir, paste0("moments_", tm, "_", gs, ".rds")),
      compress = FALSE
    )
    rm(sigma_raw, sigma_bio)
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
