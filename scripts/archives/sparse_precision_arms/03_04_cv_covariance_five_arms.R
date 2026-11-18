# mkdir -p logs output/grn/cv
# nohup Rscript --vanilla \
#   scripts/03_04_cv_covariance.R \
#   > "logs/03_04_cv_covariance_$(date +%F).log" 2>&1 &
#
# Cross-validated within-type covariance for DeCovarT, per scored type
# and time. Run from the repository root. Written, not yet run.
#
# Input. Raw UMI counts of the scored cells (`scored_slots.rds`, tier
# (i)) restricted to the gated panel of each time (`selected_<time>.csv`,
# tier (ii): permuted Welch G at BH 0.05, no cap): 1,417 genes at 48 h,
# 3,420 at 72 h, 5,821 at 96 h. `L_c` is the full RNA library of the
# cell (`slot$library_size`).
#
# Scale (`normalise`). "none" (default): raw per-cell counts, no per-cell
# library division. Dividing by the library removes the type-level
# transcriptome size and the biological size variation that the bulk
# accumulates; see 03-2, #sec-raw-vs-denoised. "library": the earlier
# `x_cg = y_cg / (L_c / median L)`, kept as a sensitivity variant.
# PLNnetwork always fits the raw counts with offset log L_c.
#
# Arms, one tuning parameter each, all scored by the same loss: the
# K-fold held-out Gaussian log-score of the cells on the scale of X.
# - shrinkage (the plug-in): Schäfer–Strimmer, correlation intensity
#   lambda, empirical variances (`corpcor::cov.shrink(lambda.var = 0)`);
#   the closed-form lambda* is one grid point;
# - glasso: graphical lasso on the correlation, unpenalised diagonal,
#   penalty lambda;
# - silggm_refit: SILGGM (de-sparsified graphical lasso on the
#   nonparanormal transform) FDR support at level alpha, then the
#   covariance refitted under those zeros;
# - plnnetwork: PLNnetwork with offset log L_c, penalty from the
#   package's own grid; latent moments mapped to the scale of X;
# - shrinkage_denoised: the shrinkage of the denoised covariance
#   (`denoise_covariance()`: V - diag(xbar), then the kappa term, then
#   eigenvalues clipped at 0). It is scored through `score_sigma`, the
#   count covariance rebuilt from the denoised estimate; the exported
#   `sigma` is the denoised estimate. Only with normalise = "none".
# The four arms after the first are sensitivity analyses.
#
# kappa. The score cannot choose kappa: every kappa rebuilds the same
# count moments up to regularisation. kappa is therefore a bracket, not
# a tuning parameter. Choice made here: the nested comparison uses
# kappa = 0 (all size variation treated as biological), which needs no
# estimate and so leaks nothing across folds. The other end,
# kappa_max = nu_hat / (1 + nu_hat) with nu_hat from
# `extrinsic_noise_nu()` on all cells of the unit, is a variant:
# `cv_select_tuning()` chooses its lambda and its fit is exported as
# `shrinkage_denoised_kappa_max`, outside the nested comparison. nu_hat
# and kappa_max are written per unit to `extrinsic_noise.csv`.
#
# PLN brackets. Beside the scored PLN moments (Poisson diagonal kept,
# depth variation kept), the PLN entry also exports the two denoised
# brackets of `pln_linear_moments()`: kappa = 0 (size term kept, no
# Poisson diagonal) and the all-technical end (size term dropped, no
# Poisson diagonal).
#
# Selection: `cv_select_tuning()` on all cells of the type (rule set by
# `selection_rule`). Comparison: `nested_cv_compare()` (inner loop
# chooses the tuning value, outer loop estimates the held-out risk).
# Fold assignments share one seed, so per-fold scores are paired across
# arms. Cost: SILGGM and PLNnetwork are refitted K (selection) and
# K_outer * (K_inner + 1) (comparison) times per type; at 96 h the panel
# has 5,821 genes, so these two arms dominate the run time. No cap on
# the panel is imposed to save time. The denoised arm costs one
# eigen-decomposition of a G x G matrix per training set.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "03_02_stats_regularisation_path.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))

time_levels <- c("48h", "72h", "96h")
normalise <- "none" # "none" (default) or "library" (sensitivity variant)
min_cells <- 20L
K <- 5L
K_outer <- 5L
K_inner <- 5L
seed <- 20261030L
selection_rule <- "min"
run_nested <- TRUE
n_lambda_glasso <- 12L
alpha_grid_silggm <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3)
n_penalties_pln <- 12L
min_ratio_pln <- 0.1
kappa_arm <- 0 # kappa of the shrinkage_denoised arm in the nested comparison
nu_min_mean <- 1 # genes entering extrinsic_noise_nu()

normalise <- match.arg(normalise, c("none", "library"))
slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
out_dir <- file.path(root, "output", "grn", "cv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# SECTION 1 · Inputs ----
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

#' Estimator specifications for one type.
#'
#' @param y Cells by genes, raw counts of the type.
#' @param lib Full library sizes of the same cells.
#' @param ref_lib Reference library of the time point.
#' @param pln_grid Penalties of the PLNnetwork family fitted on `y`.
#' @param pln_scale Scale of `X` for [pln_linear_moments()].
#' @param kappa `kappa` of the denoised shrinkage arm, or `NULL` to leave
#'   that arm out (library-divided `X`).
#' @return Named list for [cv_select_tuning()] and [nested_cv_compare()].
estimator_specs <- function(y, lib, ref_lib, pln_grid, pln_scale, kappa) {
  specs <- list(
    shrinkage = list(
      fit = fit_cov_shrinkage,
      grid = shrinkage_lambda_grid,
      regularisation = "increasing",
      path_fit = FALSE
    ),
    glasso = list(
      fit = fit_cov_glasso,
      grid = function(X, rows) {
        glasso_lambda_grid(X, n_lambda = n_lambda_glasso)
      },
      regularisation = "increasing",
      path_fit = FALSE
    ),
    silggm_refit = list(
      fit = fit_cov_silggm_refit,
      grid = alpha_grid_silggm,
      regularisation = "decreasing",
      path_fit = TRUE
    ),
    plnnetwork = list(
      fit = make_pln_fit_fun(y, lib, ref_lib, scale = pln_scale),
      grid = pln_grid,
      regularisation = "increasing",
      path_fit = TRUE
    )
  )
  if (!is.null(kappa)) {
    specs$shrinkage_denoised <- denoised_shrinkage_spec(kappa)
  }
  specs
}

#' Specification of the denoised shrinkage arm at one `kappa`.
#'
#' @param kappa Passed to [fit_cov_shrinkage()].
#' @return One estimator specification (path fit: the denoising is done
#'   once per training set).
denoised_shrinkage_spec <- function(kappa) {
  force(kappa)
  list(
    fit = function(X_train, grid, train_rows) {
      fit_cov_shrinkage(
        X_train,
        grid,
        train_rows,
        preprocess = "denoised",
        kappa = kappa,
        as_path = TRUE
      )
    },
    grid = shrinkage_lambda_grid,
    regularisation = "increasing",
    path_fit = TRUE
  )
}

# ==========================================================================
# SECTION 2 · Selection and comparison per type and time ----
# ==========================================================================

tuning_rows <- list()
selected_rows <- list()
nested_rows <- list()
nested_summary_rows <- list()
nu_rows <- list()

#' Select one arm's tuning value on all cells and refit it.
#'
#' @return A list with `row` (one row of `cv_selected_tuning.csv`),
#'   `table` (the CV table) and `entry` (the exported moments), or `NULL`.
select_and_refit <- function(spec, est_name, X, n, p, time_level, tp) {
  grid <- if (is.function(spec$grid)) spec$grid(X, seq_len(n)) else spec$grid
  message("  ", est_name, ": ", length(grid), " grid values")
  sel <- tryCatch(
    cv_select_tuning(
      X,
      fit_fun = spec$fit,
      grid = grid,
      K = K,
      seed = seed,
      regularisation = spec$regularisation,
      path_fit = spec$path_fit
    ),
    error = function(e) {
      message("  CV failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(sel)) {
    return(NULL)
  }
  chosen <- sel[[paste0("tuning_", selection_rule)]]
  final <- if (!is.finite(chosen)) {
    list(failed = TRUE)
  } else if (spec$path_fit) {
    spec$fit(X, chosen, seq_len(n))[[1L]]
  } else {
    spec$fit(X, chosen, seq_len(n))
  }
  moments <- if (isTRUE(final$failed)) {
    NULL
  } else {
    fit_moments(final, seq_len(n))
  }
  list(
    table = data.frame(
      time_point = time_level,
      cell_type = tp,
      estimator = est_name,
      sel$table,
      stringsAsFactors = FALSE
    ),
    row = data.frame(
      time_point = time_level,
      cell_type = tp,
      estimator = est_name,
      n_cells = n,
      n_genes = p,
      tuning_min = sel$tuning_min,
      tuning_1se = sel$tuning_1se,
      rule = selection_rule,
      tuning = chosen,
      n_edges = final$n_edges %||% NA_integer_,
      kappa = final$kappa %||% NA_real_,
      n_clipped = final$n_clipped %||% NA_integer_,
      n_floored = length(final$floored_genes),
      final_failed = is.null(moments),
      stringsAsFactors = FALSE
    ),
    entry = list(
      tuning = chosen,
      cv = sel,
      mu = moments$mu,
      # Exported covariance: the denoised estimate for a denoised arm.
      sigma = moments$sigma,
      # Count covariance rebuilt for scoring (denoised arms only).
      score_sigma = moments$score_sigma,
      kappa = final$kappa,
      floored_genes = final$floored_genes,
      final = final
    )
  )
}

results <- lapply(time_levels, function(time_level) {
  panel <- read_panel_slot(slots, selected_dir, time_level)
  labels <- panel$cell_type
  lib_size <- panel$library_size
  ref_lib <- stats::median(lib_size)
  types <- sort(unique(labels[!is.na(labels) & nzchar(labels)]))
  fits <- list()
  for (tp in types) {
    cells <- which(labels %in% tp)
    if (length(cells) < max(min_cells, K_outer * K_inner)) {
      message(
        "Skipping ",
        time_level,
        " / ",
        tp,
        " (n = ",
        length(cells),
        " cells)"
      )
      next
    }
    y <- panel$y[cells, , drop = FALSE]
    lib <- as.numeric(lib_size[cells])
    X <- if (identical(normalise, "none")) y else y / (lib / ref_lib)
    constant <- apply(y, 2L, function(v) all(v == v[[1L]]))
    excluded <- colnames(y)[constant]
    y <- y[, !constant, drop = FALSE]
    X <- X[, !constant, drop = FALSE]
    n <- nrow(X)
    p <- ncol(X)
    message(
      "CV covariance: ",
      time_level,
      " / ",
      tp,
      " (n = ",
      n,
      " cells, p = ",
      p,
      " genes, ",
      length(excluded),
      " constant genes left out; normalise = ",
      normalise,
      ")"
    )

    # Upper end of the kappa bracket, on raw counts of all cells of the
    # unit. Used only by the variant outside the nested comparison.
    nu_fit <- tryCatch(
      extrinsic_noise_nu(y, min_mean = nu_min_mean),
      error = function(e) {
        message("  nu estimate failed: ", conditionMessage(e))
        NULL
      }
    )
    nu_rows[[length(nu_rows) + 1L]] <<- data.frame(
      time_point = time_level,
      cell_type = tp,
      nu_hat = nu_fit$nu %||% NA_real_,
      kappa_max = nu_fit$kappa_max %||% NA_real_,
      n_genes_nu = nu_fit$n_genes %||% NA_integer_,
      n_pairs_nu = nu_fit$n_pairs %||% NA_integer_,
      stringsAsFactors = FALSE
    )

    pln_grid <- tryCatch(
      pln_penalty_grid(y, lib, n_penalties_pln, min_ratio_pln),
      error = function(e) {
        message("  PLNnetwork grid failed: ", conditionMessage(e))
        NULL
      }
    )
    pln_scale <- if (identical(normalise, "none")) "counts" else "normalised"
    specs <- estimator_specs(
      y,
      lib,
      ref_lib,
      pln_grid,
      pln_scale,
      kappa = if (identical(normalise, "none")) kappa_arm else NULL
    )
    if (is.null(pln_grid)) {
      specs$plnnetwork <- NULL
    }
    # The kappa_max variant is selected and exported, not compared.
    variants <- list()
    if (identical(normalise, "none") && !is.null(nu_fit)) {
      variants$shrinkage_denoised_kappa_max <-
        denoised_shrinkage_spec(nu_fit$kappa_max)
    }

    per_estimator <- list()
    for (est_name in c(names(specs), names(variants))) {
      spec <- if (est_name %in% names(specs)) {
        specs[[est_name]]
      } else {
        variants[[est_name]]
      }
      res <- select_and_refit(spec, est_name, X, n, p, time_level, tp)
      if (is.null(res)) {
        next
      }
      tuning_rows[[length(tuning_rows) + 1L]] <<- res$table
      selected_rows[[length(selected_rows) + 1L]] <<- res$row
      entry <- res$entry
      final <- entry$final
      entry$final <- NULL
      if (identical(est_name, "plnnetwork") && !isTRUE(final$failed)) {
        # The two denoised brackets of the PLN moments (no Poisson
        # diagonal), on the scale of X.
        entry$brackets <- list(
          kappa_0 = pln_linear_moments(
            final$mu_z,
            final$sigma_z,
            lib,
            ref_lib,
            scale = pln_scale,
            size_term = "keep",
            poisson_diag = FALSE
          ),
          all_technical = pln_linear_moments(
            final$mu_z,
            final$sigma_z,
            lib,
            ref_lib,
            scale = pln_scale,
            size_term = "drop",
            poisson_diag = FALSE
          )
        )
      }
      per_estimator[[est_name]] <- entry
    }

    if (run_nested && length(specs) > 0L) {
      nested <- nested_cv_compare(
        X,
        specs,
        K_outer = K_outer,
        K_inner = K_inner,
        seed = seed,
        rule = selection_rule
      )
      nested_rows[[length(nested_rows) + 1L]] <<- data.frame(
        time_point = time_level,
        cell_type = tp,
        nested$by_fold,
        stringsAsFactors = FALSE
      )
      nested_summary_rows[[length(nested_summary_rows) + 1L]] <<- data.frame(
        time_point = time_level,
        cell_type = tp,
        nested$summary,
        stringsAsFactors = FALSE
      )
    }

    fits[[tp]] <- list(
      genes = colnames(X),
      excluded_genes = excluded,
      n = n,
      ref_lib = ref_lib,
      normalise = normalise,
      scale = if (identical(normalise, "none")) {
        "raw UMI counts (no per-cell library division)"
      } else {
        "library-divided counts (full library / median)"
      },
      nu = nu_fit,
      estimators = per_estimator
    )
  }
  fits
})
names(results) <- time_levels

# ==========================================================================
# SECTION 3 · Outputs ----
# ==========================================================================

suffix <- if (identical(normalise, "none")) "" else paste0("_", normalise)
out_rds <- file.path(out_dir, paste0("cv_covariance_by_type_time", suffix, ".rds"))
saveRDS(results, out_rds)
message("Wrote ", out_rds)

write_table <- function(rows, filename) {
  if (length(rows) == 0L) {
    message("Nothing to write for ", filename)
    return(invisible(NULL))
  }
  path <- file.path(out_dir, sub("\\.csv$", paste0(suffix, ".csv"), filename))
  utils::write.csv(do.call(rbind, rows), path, row.names = FALSE)
  message("Wrote ", path)
}
write_table(tuning_rows, "cv_tuning_scores.csv")
write_table(selected_rows, "cv_selected_tuning.csv")
write_table(nested_rows, "nested_cv_by_fold.csv")
write_table(nested_summary_rows, "nested_cv_summary.csv")
write_table(nu_rows, "extrinsic_noise.csv")
