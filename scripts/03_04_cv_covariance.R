# mkdir -p logs output/grn/cv
# nohup Rscript --vanilla \
#   scripts/03_04_cv_covariance.R \
#   > "logs/03_04_cv_covariance_$(date +%F).log" 2>&1 &
#
# Cross-validated within-type covariance for DeCovarT, per labelled type
# and time. Run from the repository root. Written, not yet run.
#
# Four estimators, one tuning parameter each, all scored by the same
# loss: the K-fold held-out Gaussian log-score of the cells on the
# linear scale of the convolution (size-factor-normalised counts,
# `x_cg = y_cg / s_c`, `s_c = L_c / median L`, `L_c` the full RNA
# library from metadata `log_umi`).
#
# - shrinkage: Schäfer–Strimmer, correlation intensity lambda, empirical
#   variances (`corpcor::cov.shrink(lambda.var = 0)`); the closed-form
#   lambda* is one grid point;
# - glasso: graphical lasso on the correlation, unpenalised diagonal,
#   penalty lambda;
# - silggm_refit: SILGGM (de-sparsified graphical lasso on the
#   nonparanormal transform) FDR support at level alpha, then the
#   covariance refitted under those zeros;
# - plnnetwork: PLNnetwork with offset log L_c, penalty from the
#   package's own grid; latent moments mapped to the linear scale.
#
# Selection: `cv_select_tuning()` on all cells of the type (rule set by
# `selection_rule`). Comparison: `nested_cv_compare()` (inner loop
# chooses the tuning value, outer loop estimates the held-out risk).
# Fold assignments share one seed, so per-fold scores are paired across
# estimators. Cost: SILGGM and PLNnetwork are refitted K (selection) and
# K_outer * (K_inner + 1) (comparison) times per type.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "ebic_regularisation_path.R"))
source(file.path(root, "R", "cv_covariance_selection.R"))

time_levels <- c("48h", "72h", "96h")
cell_type_col <- "celltypeannotation"
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

input_rds <- file.path(
  root,
  "data",
  "intermediate",
  "equal_type_top500_by_time.rds"
)
out_dir <- file.path(root, "output", "grn", "cv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# SECTION 1 · Inputs ----
# ==========================================================================

if (!file.exists(input_rds)) {
  stop(
    "Missing ",
    input_rds,
    ". Run scripts/01_03_export_equal_type_seurat.R first."
  )
}
objects <- readRDS(input_rds)

#' Estimator specifications for one type.
#'
#' @param y Cells by genes, raw counts of the type.
#' @param lib Full library sizes of the same cells.
#' @param ref_lib Reference library of the time point.
#' @param pln_grid Penalties of the PLNnetwork family fitted on `y`.
#' @return Named list for [cv_select_tuning()] and [nested_cv_compare()].
estimator_specs <- function(y, lib, ref_lib, pln_grid) {
  list(
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
      fit = make_pln_fit_fun(y, lib, ref_lib),
      grid = pln_grid,
      regularisation = "increasing",
      path_fit = TRUE
    )
  )
}

# ==========================================================================
# SECTION 2 · Selection and comparison per type and time ----
# ==========================================================================

tuning_rows <- list()
selected_rows <- list()
nested_rows <- list()
nested_summary_rows <- list()

results <- lapply(time_levels, function(time_level) {
  seu <- objects[[time_level]]
  counts <- assay_matrix(seu, "RNA", "counts")
  labels <- seurat_metadata_chr(seu, cell_type_col)
  lib_size <- 10^as.numeric(seurat_metadata_chr(seu, "log_umi"))
  names(lib_size) <- colnames(seu)
  if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    stop("Metadata log_umi is missing or invalid at ", time_level, ".")
  }
  ref_lib <- stats::median(lib_size)
  types <- sort(unique(labels[!is.na(labels) & nzchar(labels)]))
  fits <- list()
  for (tp in types) {
    cells <- unique(colnames(seu)[labels %in% tp])
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
    y <- t(as.matrix(counts[, cells, drop = FALSE]))
    lib <- as.numeric(lib_size[cells])
    X <- y / (lib / ref_lib)
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
      " constant genes left out)"
    )

    pln_grid <- tryCatch(
      pln_penalty_grid(y, lib, n_penalties_pln, min_ratio_pln),
      error = function(e) {
        message("  PLNnetwork grid failed: ", conditionMessage(e))
        NULL
      }
    )
    specs <- estimator_specs(y, lib, ref_lib, pln_grid)
    if (is.null(pln_grid)) {
      specs$plnnetwork <- NULL
    }

    per_estimator <- list()
    for (est_name in names(specs)) {
      spec <- specs[[est_name]]
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
        next
      }
      tuning_rows[[length(tuning_rows) + 1L]] <<- data.frame(
        time_point = time_level,
        cell_type = tp,
        estimator = est_name,
        sel$table,
        stringsAsFactors = FALSE
      )
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
      selected_rows[[length(selected_rows) + 1L]] <<- data.frame(
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
        final_failed = is.null(moments),
        stringsAsFactors = FALSE
      )
      per_estimator[[est_name]] <- list(
        tuning = chosen,
        cv = sel,
        mu = moments$mu,
        sigma = moments$sigma
      )
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
      scale = "size-factor-normalised counts (full library / median)",
      estimators = per_estimator
    )
  }
  fits
})
names(results) <- time_levels

# ==========================================================================
# SECTION 3 · Outputs ----
# ==========================================================================

out_rds <- file.path(out_dir, "cv_covariance_by_type_time.rds")
saveRDS(results, out_rds)
message("Wrote ", out_rds)

write_table <- function(rows, filename) {
  if (length(rows) == 0L) {
    message("Nothing to write for ", filename)
    return(invisible(NULL))
  }
  utils::write.csv(
    do.call(rbind, rows),
    file.path(out_dir, filename),
    row.names = FALSE
  )
  message("Wrote ", file.path(out_dir, filename))
}
write_table(tuning_rows, "cv_tuning_scores.csv")
write_table(selected_rows, "cv_selected_tuning.csv")
write_table(nested_rows, "nested_cv_by_fold.csv")
write_table(nested_summary_rows, "nested_cv_summary.csv")
