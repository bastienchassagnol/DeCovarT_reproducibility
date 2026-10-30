#' Cross-validated selection and comparison of within-type covariances.
#'
#' Every estimator in this file returns a mean `mu` and a covariance
#' `sigma` on the scale of DeCovarT's convolution: size-factor-normalised
#' counts, `x_cg = y_cg / s_c` with `s_c = L_c / L_ref`, where `L_c` is
#' the full RNA library of cell `c` and `L_ref` a reference library (the
#' median over the cells of one time point). All four estimators are
#' therefore scored by the same loss: the held-out Gaussian log-score of
#' the test cells.
#'
#' Why the Gaussian log-score is a legitimate loss when cells are not
#' Gaussian: for any law of the cells with mean `m` and covariance
#' `Sigma_true`, the expected negative log-score of `N(mu, Sigma)` is
#' `0.5 * {log det Sigma + tr(Sigma^-1 Sigma_true)
#'   + (m - mu)' Sigma^-1 (m - mu)} + const`,
#' which is minimised at `mu = m`, `Sigma = Sigma_true` (Stein's loss).
#' DeCovarT's bulk model uses only these two moments, so the score ranks
#' covariance estimators by the quantity the deconvolution needs.
#'
#' Requires `R/utils_general.R` (`%||%`) and
#' `R/ebic_regularisation_path.R` (`make_glasso_solver()`,
#' `rescale_correlation_fit()`, `count_undirected_edges()`) to be
#' sourced first.

# ==========================================================================
# Held-out score and folds ----
# ==========================================================================

#' Held-out Gaussian log-score through a Cholesky factor.
#'
#' With `Sigma = R'R` (upper-triangular `R`), the log-density of a test
#' row `x` is `-0.5 * {p log(2 pi) + 2 sum(log diag R) + |R^-T (x - mu)|^2}`.
#' No inverse is formed.
#'
#' @param mu Mean vector, length `p`.
#' @param Sigma Covariance, `p` by `p`. It is symmetrised before the
#'   factorisation.
#' @param X_test Test cells in rows, genes in columns (`n_test` by `p`).
#' @param average If `TRUE` (default), return the mean log-score per
#'   test cell; otherwise one log-density per row.
#' @return A scalar (or a vector of length `n_test`). `-Inf` when
#'   `Sigma` is not numerically positive definite, so that such a fit
#'   can never be selected.
heldout_gaussian_logscore <- function(mu, Sigma, X_test, average = TRUE) {
  X_test <- as.matrix(X_test)
  Sigma <- as.matrix(Sigma)
  mu <- as.numeric(mu)
  p <- ncol(X_test)
  if (length(mu) != p || nrow(Sigma) != p || ncol(Sigma) != p) {
    stop("`mu`, `Sigma` and `X_test` must share the gene dimension.")
  }
  Sigma <- (Sigma + t(Sigma)) / 2
  R <- if (all(is.finite(Sigma))) {
    tryCatch(chol(Sigma), error = function(e) NULL)
  } else {
    NULL
  }
  if (is.null(R)) {
    return(if (average) -Inf else rep(-Inf, nrow(X_test)))
  }
  centred <- t(X_test) - mu
  z <- backsolve(R, centred, transpose = TRUE)
  log_det <- 2 * sum(log(diag(R)))
  per_cell <- -0.5 * (p * log(2 * pi) + log_det + colSums(z * z))
  if (average) mean(per_cell) else per_cell
}

#' Balanced random assignment of rows to K folds.
#'
#' @param n Number of rows (cells).
#' @param K Number of folds.
#' @param seed Optional integer seed. When `NULL`, the current RNG
#'   stream is used and left advanced.
#' @return Integer vector of length `n` with values in `1:K`.
make_cv_folds <- function(n, K = 5L, seed = NULL) {
  K <- as.integer(K)
  if (K < 2L || K > n) {
    stop("Need 2 <= K <= n; got K = ", K, " for n = ", n, ".")
  }
  assign_folds <- function() sample(rep_len(seq_len(K), n))
  if (is.null(seed)) {
    assign_folds()
  } else {
    withr::with_seed(seed, assign_folds())
  }
}

#' Mean and covariance that a fit predicts for a set of test cells.
#'
#' Most fits return fixed `mu` and `sigma`. A PLNnetwork fit returns a
#' `moments` function instead, because its linear-scale covariance
#' depends on the size factors of the cells it is evaluated on (see
#' [pln_linear_moments()]).
#'
#' @param fit List with `mu` and `sigma`, or with `moments`, a function
#'   of the test-row indices.
#' @param test_rows Indices of the test cells in the caller's full data.
#' @return A list with `mu` and `sigma`.
fit_moments <- function(fit, test_rows) {
  if (is.function(fit$moments)) {
    fit$moments(test_rows)
  } else {
    list(mu = fit$mu, sigma = fit$sigma)
  }
}

#' Held-out log-score of one fit, `NA` for a failed fit.
#'
#' @param fit A fit as described in [fit_moments()], possibly with
#'   `failed = TRUE`.
#' @param X_test Test cells by genes.
#' @param test_rows Indices of the test cells in the full data.
#' @return Scalar mean log-score per test cell.
score_fit <- function(fit, X_test, test_rows) {
  if (is.null(fit) || isTRUE(fit$failed)) {
    return(NA_real_)
  }
  m <- fit_moments(fit, test_rows)
  heldout_gaussian_logscore(m$mu, m$sigma, X_test)
}

# ==========================================================================
# K-fold selection of one tuning parameter ----
# ==========================================================================

#' K-fold cross-validated choice of one tuning parameter.
#'
#' For each fold, every grid value is fitted on the training cells and
#' scored on the held-out cells by [heldout_gaussian_logscore()]. The
#' `min` choice maximises the mean score; the `1se` choice is the most
#' regularised value whose mean score is within one standard error of
#' that maximum. A grid value with a failed fit in any fold is not
#' eligible.
#'
#' `fit_fun(X_train, tuning, train_rows)` returns one fit (a list with
#' `mu`/`sigma` or `moments`, or `failed = TRUE`). When `path_fit =
#' TRUE` it receives the whole grid and returns a list of fits in grid
#' order; SILGGM (one fit for all FDR levels) and PLNnetwork (one fit for
#' a whole penalty path) are run this way.
#'
#' @param X Cells by genes, on the linear scale used for scoring.
#' @param fit_fun Fit function, see above.
#' @param grid Numeric vector of candidate tuning values.
#' @param K Number of folds. Default `5`.
#' @param seed Optional seed for the fold assignment.
#' @param rows Indices of the rows of `X` in the caller's full data,
#'   passed on to `fit_fun` and to `moments`. Default `seq_len(nrow(X))`.
#' @param regularisation `"increasing"` (default) when a larger value
#'   regularises more (lasso penalty, shrinkage intensity);
#'   `"decreasing"` when a smaller value does (SILGGM FDR level).
#' @param path_fit Whether `fit_fun` takes the whole grid at once.
#' @param folds Optional precomputed fold labels (overrides `K`, `seed`).
#' @param verbose Print progress.
#' @return A list with `table` (`tuning`, `mean_score`, `se_score`,
#'   `n_folds_ok`), `fold_scores` (grid by fold), `tuning_min`,
#'   `tuning_1se`, `folds`, and `regularisation`.
cv_select_tuning <- function(
  X,
  fit_fun,
  grid,
  K = 5L,
  seed = NULL,
  rows = NULL,
  regularisation = c("increasing", "decreasing"),
  path_fit = FALSE,
  folds = NULL,
  verbose = FALSE
) {
  regularisation <- match.arg(regularisation)
  X <- as.matrix(X)
  n <- nrow(X)
  rows <- rows %||% seq_len(n)
  if (length(rows) != n) {
    stop("`rows` must have one entry per row of `X`.")
  }
  grid <- as.numeric(grid)
  if (length(grid) == 0L || any(!is.finite(grid))) {
    stop("`grid` must be a non-empty vector of finite values.")
  }
  folds <- folds %||% make_cv_folds(n, K, seed)
  n_folds <- max(folds)
  scores <- matrix(
    NA_real_,
    nrow = length(grid),
    ncol = n_folds,
    dimnames = list(NULL, paste0("fold_", seq_len(n_folds)))
  )
  for (k in seq_len(n_folds)) {
    test <- which(folds == k)
    train <- which(folds != k)
    if (verbose) {
      message("  fold ", k, "/", n_folds, " (", length(test), " test cells)")
    }
    X_train <- X[train, , drop = FALSE]
    fits <- if (path_fit) {
      fit_fun(X_train, grid, rows[train])
    } else {
      lapply(grid, function(g) fit_fun(X_train, g, rows[train]))
    }
    if (length(fits) != length(grid)) {
      stop("`fit_fun` returned ", length(fits), " fits for ", length(grid),
           " grid values.")
    }
    for (i in seq_along(grid)) {
      scores[i, k] <- score_fit(fits[[i]], X[test, , drop = FALSE], rows[test])
    }
  }

  n_ok <- rowSums(!is.na(scores))
  mean_score <- rowMeans(scores)
  se_score <- apply(scores, 1L, stats::sd) / sqrt(n_folds)
  mean_score[n_ok < n_folds] <- NA_real_
  se_score[n_ok < n_folds] <- NA_real_
  table <- data.frame(
    tuning = grid,
    mean_score = mean_score,
    se_score = se_score,
    n_folds_ok = n_ok
  )

  eligible <- which(is.finite(mean_score))
  if (length(eligible) == 0L) {
    warning("No grid value was fitted in every fold.", call. = FALSE)
    return(list(
      table = table,
      fold_scores = scores,
      tuning_min = NA_real_,
      tuning_1se = NA_real_,
      folds = folds,
      regularisation = regularisation
    ))
  }
  i_min <- eligible[[which.max(mean_score[eligible])]]
  se_min <- se_score[[i_min]]
  if (!is.finite(se_min)) {
    se_min <- 0
  }
  within <- eligible[mean_score[eligible] >= mean_score[[i_min]] - se_min]
  i_1se <- if (identical(regularisation, "increasing")) {
    within[[which.max(grid[within])]]
  } else {
    within[[which.min(grid[within])]]
  }
  list(
    table = table,
    fold_scores = scores,
    tuning_min = grid[[i_min]],
    tuning_1se = grid[[i_1se]],
    folds = folds,
    regularisation = regularisation
  )
}

# ==========================================================================
# Nested cross-validation across estimators ----
# ==========================================================================

#' Nested cross-validation of several covariance estimators.
#'
#' The outer loop estimates the held-out risk of each estimator; the
#' inner loop, run on the outer training cells only, chooses that
#' estimator's tuning value. Inner folds use the same seed for every
#' estimator within an outer fold, so the per-fold scores are paired.
#' Nested CV is needed only for this comparison; the tuning value that is
#' finally exported comes from [cv_select_tuning()] on all cells of the
#' type.
#'
#' Each element of `estimators` is a list with
#' - `fit`: a fit function as in [cv_select_tuning()];
#' - `grid`: a numeric grid, or a function `(X_train, train_rows)`
#'   returning one, so that a data-driven grid is built from the outer
#'   training cells only;
#' - `regularisation` (optional, default `"increasing"`);
#' - `path_fit` (optional, default `FALSE`).
#'
#' @param X Cells by genes, on the scoring scale.
#' @param estimators Named list of estimator specifications.
#' @param K_outer,K_inner Numbers of outer and inner folds. Default `5`.
#' @param seed Optional seed. The outer folds use `seed`; the inner
#'   folds of outer fold `k` use `seed + k`.
#' @param rule `"min"` (default) or `"1se"`, the inner-loop choice.
#' @param rows Indices of the rows of `X` in the caller's full data.
#' @param verbose Print progress.
#' @return A list with `by_fold` (one row per outer fold and estimator:
#'   `outer_fold`, `estimator`, `tuning`, `n_test`, `score`), `summary`
#'   (per estimator: `mean_score`, `se_score`, `n_folds_ok`) and
#'   `outer_folds`.
nested_cv_compare <- function(
  X,
  estimators,
  K_outer = 5L,
  K_inner = 5L,
  seed = NULL,
  rule = c("min", "1se"),
  rows = NULL,
  verbose = FALSE
) {
  rule <- match.arg(rule)
  if (is.null(names(estimators)) || any(!nzchar(names(estimators)))) {
    stop("`estimators` must be a named list.")
  }
  X <- as.matrix(X)
  n <- nrow(X)
  rows <- rows %||% seq_len(n)
  outer <- make_cv_folds(n, K_outer, seed)
  by_fold <- list()
  for (k in seq_len(max(outer))) {
    test <- which(outer == k)
    train <- which(outer != k)
    X_train <- X[train, , drop = FALSE]
    inner_seed <- if (is.null(seed)) NULL else seed + k
    inner_folds <- make_cv_folds(length(train), K_inner, inner_seed)
    for (est_name in names(estimators)) {
      est <- estimators[[est_name]]
      if (verbose) {
        message("Outer fold ", k, ": ", est_name)
      }
      grid <- if (is.function(est$grid)) {
        est$grid(X_train, rows[train])
      } else {
        est$grid
      }
      path_fit <- isTRUE(est$path_fit)
      sel <- tryCatch(
        cv_select_tuning(
          X_train,
          fit_fun = est$fit,
          grid = grid,
          rows = rows[train],
          regularisation = est$regularisation %||% "increasing",
          path_fit = path_fit,
          folds = inner_folds
        ),
        error = function(e) {
          message("  inner CV failed: ", conditionMessage(e))
          NULL
        }
      )
      tuning <- if (is.null(sel)) NA_real_ else sel[[paste0("tuning_", rule)]]
      score <- NA_real_
      if (is.finite(tuning)) {
        fit <- tryCatch(
          if (path_fit) {
            est$fit(X_train, tuning, rows[train])[[1L]]
          } else {
            est$fit(X_train, tuning, rows[train])
          },
          error = function(e) list(failed = TRUE)
        )
        score <- score_fit(fit, X[test, , drop = FALSE], rows[test])
      }
      by_fold[[length(by_fold) + 1L]] <- data.frame(
        outer_fold = k,
        estimator = est_name,
        tuning = tuning,
        n_test = length(test),
        score = score,
        stringsAsFactors = FALSE
      )
    }
  }
  by_fold <- do.call(rbind, by_fold)
  summary <- do.call(
    rbind,
    lapply(split(by_fold, by_fold$estimator), function(d) {
      ok <- is.finite(d$score)
      data.frame(
        estimator = d$estimator[[1L]],
        mean_score = if (all(ok)) mean(d$score) else NA_real_,
        se_score = if (all(ok)) stats::sd(d$score) / sqrt(nrow(d)) else NA_real_,
        n_folds_ok = sum(ok),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(summary) <- NULL
  list(by_fold = by_fold, summary = summary, outer_folds = outer)
}

# ==========================================================================
# Shared helpers for the correlation-scale estimators ----
# ==========================================================================

#' Means, standard deviations and correlation of the genes.
#'
#' A gene that is constant in the training cells has no correlation with
#' any other gene: its row and column of `R` are zero (diagonal one), and
#' its variance is floored so that the covariance stays positive
#' definite. The floor is a numerical guard, not a model; a gene constant
#' in training and variable in the test cells is still scored badly.
#'
#' @param X Cells by genes.
#' @param var_floor Variance floor. Default `NULL`: `1e-4` times the
#'   median positive variance (or `1` when no gene varies).
#' @return A list with `mu`, `sd`, `R`, `constant` (logical), `n`.
standardise_genes <- function(X, var_floor = NULL) {
  X <- as.matrix(X)
  p <- ncol(X)
  mu <- colMeans(X)
  v <- apply(X, 2L, stats::var)
  constant <- !is.finite(v) | v <= 0
  if (is.null(var_floor)) {
    var_floor <- if (any(!constant)) 1e-4 * stats::median(v[!constant]) else 1
  }
  v[constant] <- 0
  sd <- sqrt(pmax(v, var_floor))
  R <- diag(1, p)
  if (sum(!constant) >= 2L) {
    R[!constant, !constant] <- stats::cor(X[, !constant, drop = FALSE])
  }
  R <- (R + t(R)) / 2
  diag(R) <- 1
  dimnames(R) <- list(colnames(X), colnames(X))
  names(sd) <- colnames(X)
  list(mu = mu, sd = sd, R = R, constant = constant, n = nrow(X))
}

#' Covariance refitted under a fixed zero pattern of the precision.
#'
#' Maximises the Gaussian likelihood of the correlation `R` with the
#' precision entries outside `adjacency` constrained to zero, through
#' `glasso::glasso(zero = ...)`. A small `rho` keeps the problem well
#' posed when cells are fewer than genes; the diagonal is not penalised.
#' The result is rescaled to the data scale by `sd`.
#'
#' @param R Correlation matrix, `p` by `p`.
#' @param adjacency `p` by `p` matrix; non-zero off-diagonal entries are
#'   the edges kept free.
#' @param sd Standard deviations on the data scale.
#' @param rho Ridge-like penalty on the free off-diagonal entries.
#'   Default `1e-4`.
#' @param nobs Number of cells behind `R` (only fills `loglik`).
#' @param thr,maxit Passed to [glasso::glasso()].
#' @return A list with `sigma`, `omega` (data scale), `n_edges`, and
#'   `failed`.
refit_covariance_on_support <- function(
  R,
  adjacency,
  sd,
  rho = 1e-4,
  nobs = NULL,
  thr = 1e-4,
  maxit = 1e4
) {
  free <- as.matrix(adjacency) != 0
  free <- free | t(free)
  diag(free) <- TRUE
  zero <- which(!free, arr.ind = TRUE)
  if (nrow(zero) == 0L) {
    zero <- NULL
  }
  fit <- tryCatch(
    glasso::glasso(
      R,
      rho = rho,
      nobs = nobs,
      zero = zero,
      thr = thr,
      maxit = maxit,
      penalize.diagonal = FALSE
    ),
    error = function(e) e
  )
  if (inherits(fit, "error") || fit$errflag != 0L) {
    return(list(failed = TRUE))
  }
  omega <- (fit$wi + t(fit$wi)) / 2
  sigma <- (fit$w + t(fit$w)) / 2
  if (!all(is.finite(omega)) || !all(is.finite(sigma))) {
    return(list(failed = TRUE))
  }
  dimnames(omega) <- dimnames(R)
  dimnames(sigma) <- dimnames(R)
  scaled <- rescale_correlation_fit(omega, sigma, sd)
  list(
    sigma = scaled$sigma,
    omega = scaled$omega,
    n_edges = sum(free[upper.tri(free)]),
    failed = FALSE
  )
}

# ==========================================================================
# Estimator 1 · Schäfer–Strimmer shrinkage ----
# ==========================================================================

#' Linear shrinkage covariance with a fixed intensity.
#'
#' `corpcor::cov.shrink()` shrinks the correlations towards the
#' identity with intensity `lambda` and keeps the empirical variances
#' (`lambda.var = 0`), so the target is the diagonal of the empirical
#' covariance. Constant genes are kept on the diagonal with the floored
#' variance of [standardise_genes()].
#'
#' @param X Cells by genes (training cells).
#' @param lambda Correlation shrinkage intensity in `[0, 1]`.
#' @param train_rows Unused; kept for the common fit signature.
#' @param var_floor Passed to [standardise_genes()].
#' @return A list with `mu`, `sigma`, `failed`.
fit_cov_shrinkage <- function(X, lambda, train_rows = NULL, var_floor = NULL) {
  X <- as.matrix(X)
  st <- standardise_genes(X, var_floor)
  vary <- !st$constant
  sigma <- diag(st$sd^2, ncol(X))
  if (sum(vary) >= 2L) {
    s <- tryCatch(
      corpcor::cov.shrink(
        X[, vary, drop = FALSE],
        lambda = lambda,
        lambda.var = 0,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(s)) {
      return(list(failed = TRUE))
    }
    sigma[vary, vary] <- matrix(as.numeric(s), sum(vary), sum(vary))
  }
  dimnames(sigma) <- list(colnames(X), colnames(X))
  list(mu = st$mu, sigma = sigma, failed = FALSE)
}

#' Closed-form Schäfer–Strimmer shrinkage intensity.
#'
#' Used as one point of the CV grid (and as its natural start), not as
#' the final choice.
#'
#' @param X Cells by genes.
#' @return Scalar `lambda*` in `[0, 1]`.
shrinkage_lambda_closed_form <- function(X) {
  X <- as.matrix(X)
  vary <- apply(X, 2L, stats::var) > 0
  as.numeric(corpcor::estimate.lambda(X[, vary, drop = FALSE], verbose = FALSE))
}

#' Default shrinkage grid: fixed points in `(0, 1)` plus `lambda*`.
#'
#' @param X Cells by genes.
#' @param train_rows Unused; kept for the grid-function signature.
#' @return Sorted unique numeric grid.
shrinkage_lambda_grid <- function(X, train_rows = NULL) {
  fixed <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9)
  star <- tryCatch(shrinkage_lambda_closed_form(X), error = function(e) NULL)
  grid <- c(fixed, star)
  sort(unique(pmin(pmax(grid, 1e-3), 1)))
}

# ==========================================================================
# Estimator 2 · Graphical lasso, unpenalised diagonal ----
# ==========================================================================

#' Graphical lasso on the correlation, rescaled to covariance.
#'
#' Same solver as the EBIC path (`make_glasso_solver()`:
#' `penalize.diagonal = FALSE`), so the penalty is on the correlation
#' scale for every gene; the fit is rescaled by the gene standard
#' deviations with `rescale_correlation_fit()`.
#'
#' @param X Cells by genes (training cells).
#' @param lambda Lasso penalty on the correlation scale.
#' @param train_rows Unused; kept for the common fit signature.
#' @param var_floor Passed to [standardise_genes()].
#' @param thr,maxit Passed to [glasso::glasso()].
#' @return A list with `mu`, `sigma`, `omega`, `n_edges`, `failed`.
fit_cov_glasso <- function(
  X,
  lambda,
  train_rows = NULL,
  var_floor = NULL,
  thr = 1e-4,
  maxit = 1e4
) {
  st <- standardise_genes(X, var_floor)
  fit <- make_glasso_solver(st$R, st$n, thr = thr, maxit = maxit)(lambda)
  if (isTRUE(fit$failed)) {
    return(list(failed = TRUE))
  }
  scaled <- rescale_correlation_fit(fit$omega, fit$sigma, st$sd)
  list(
    mu = st$mu,
    sigma = scaled$sigma,
    omega = scaled$omega,
    n_edges = count_undirected_edges(fit$omega),
    failed = FALSE
  )
}

#' Log-spaced graphical-lasso grid on the correlation scale.
#'
#' Runs from `lambda_max = max |R_gh|` (empty graph) down to
#' `min_ratio * lambda_max`.
#'
#' @param X Cells by genes.
#' @param train_rows Unused; kept for the grid-function signature.
#' @param n_lambda Number of grid values. Default `12`.
#' @param min_ratio Smallest value as a fraction of `lambda_max`.
#'   Default `0.05`.
#' @return Decreasing numeric grid.
glasso_lambda_grid <- function(
  X,
  train_rows = NULL,
  n_lambda = 12L,
  min_ratio = 0.05
) {
  R <- standardise_genes(X)$R
  lambda_max <- max(abs(R[upper.tri(R)]))
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    stop("No non-zero correlation: the graphical-lasso grid is empty.")
  }
  exp(seq(log(lambda_max), log(min_ratio * lambda_max), length.out = n_lambda))
}

# ==========================================================================
# Estimator 3 · SILGGM support, then refit under the zeros ----
# ==========================================================================

#' SILGGM FDR-controlled support, refitted on the linear scale.
#'
#' The support comes from SILGGM's global FDR decision on the
#' nonparanormal (shrinkage) transform of the training cells. The
#' covariance is then refitted on the linear-scale correlation with the
#' precision constrained to zero off that support
#' ([refit_covariance_on_support()]). The tuning parameter is the FDR
#' level `alpha`; a smaller `alpha` gives a sparser graph. One SILGGM
#' run serves every `alpha`, so this wrapper always returns a list of
#' fits (use `path_fit = TRUE` in [cv_select_tuning()]).
#'
#' @param X Cells by genes (training cells), linear scale.
#' @param alpha Vector of FDR levels.
#' @param train_rows Unused; kept for the common fit signature.
#' @param method SILGGM method. Default `"D-S_GL"` (de-sparsified
#'   graphical lasso).
#' @param npn_func Passed to `huge::huge.npn(npn.func = ...)`. Default
#'   `"shrinkage"`.
#' @param rho_refit Passed to [refit_covariance_on_support()].
#' @param var_floor Passed to [standardise_genes()].
#' @return A list of fits, one per `alpha`, each with `mu`, `sigma`,
#'   `omega`, `n_edges`, `failed`.
fit_cov_silggm_refit <- function(
  X,
  alpha,
  train_rows = NULL,
  method = "D-S_GL",
  npn_func = "shrinkage",
  rho_refit = 1e-4,
  var_floor = NULL
) {
  X <- as.matrix(X)
  st <- standardise_genes(X, var_floor)
  vary <- !st$constant
  p <- ncol(X)
  failed_all <- rep(list(list(failed = TRUE)), length(alpha))
  if (sum(vary) < 2L) {
    return(failed_all)
  }
  x_npn <- huge::huge.npn(
    X[, vary, drop = FALSE],
    npn.func = npn_func,
    verbose = FALSE
  )
  sg <- tryCatch(
    SILGGM::SILGGM(
      x_npn,
      method = method,
      global = TRUE,
      alpha = alpha,
      cytoscape_format = FALSE,
      csv_save = FALSE
    ),
    error = function(e) e
  )
  if (inherits(sg, "error")) {
    message("  SILGGM failed (", conditionMessage(sg), ").")
    return(failed_all)
  }
  decisions <- sg$global_decision
  if (is.matrix(decisions)) {
    decisions <- list(decisions)
  }
  if (length(decisions) != length(alpha)) {
    stop("SILGGM returned ", length(decisions), " decisions for ",
         length(alpha), " FDR levels.")
  }
  lapply(decisions, function(decision) {
    adjacency <- matrix(0, p, p)
    adjacency[vary, vary] <- as.matrix(decision)
    refit <- refit_covariance_on_support(
      st$R,
      adjacency,
      st$sd,
      rho = rho_refit,
      nobs = st$n
    )
    if (isTRUE(refit$failed)) {
      return(refit)
    }
    c(list(mu = st$mu), refit)
  })
}

# ==========================================================================
# Estimator 4 · PLNnetwork, mapped to linear-scale moments ----
# ==========================================================================

#' Linear-scale moments of a Poisson log-normal model.
#'
#' Model: `y_cg | z_c ~ Poisson(L_c exp(z_cg))`, `z_c ~ N(mu_z, Sigma_z)`,
#' independent over genes given `z_c`. Write
#' `lambda = exp(mu_z + diag(Sigma_z) / 2)` and
#' `K = exp(Sigma_z) - 1` (entrywise). Then, given `L_c`,
#' `E y = L_c lambda` and
#' `Cov y = diag(L_c lambda) + L_c^2 (lambda lambda') * K`.
#'
#' - `scale = "normalised"`: `x_c = y_c L_ref / L_c`. With
#'   `nu = L_ref lambda`, `E x = nu` for every cell and, averaging over
#'   the test cells (law of total covariance, the conditional means do
#'   not vary), `Cov x = diag(nu) * mean(L_ref / L_c) + (nu nu') * K`.
#' - `scale = "counts"`: raw counts of a cell drawn from the test cells,
#'   `E y = mean(L) lambda` and
#'   `Cov y = diag(lambda) mean(L) + mean(L^2) (lambda lambda') * K
#'   + var(L) lambda lambda'`, with the population variance of `L`.
#'
#' @param mu_z Latent mean (length `p`), on the log scale per unit of
#'   library (the PLN intercept with offset `log L_c`).
#' @param sigma_z Latent covariance, `p` by `p`.
#' @param lib_size Library sizes `L_c` of the cells the moments are for.
#' @param ref_lib Reference library `L_ref`. Default the median of
#'   `lib_size`; pass the value used to normalise the data.
#' @param scale `"normalised"` (default) or `"counts"`.
#' @return A list with `mu` and `sigma`.
pln_linear_moments <- function(
  mu_z,
  sigma_z,
  lib_size,
  ref_lib = stats::median(lib_size),
  scale = c("normalised", "counts")
) {
  scale <- match.arg(scale)
  sigma_z <- as.matrix(sigma_z)
  mu_z <- as.numeric(mu_z)
  if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    stop("`lib_size` must be positive and finite.")
  }
  lambda <- exp(mu_z + diag(sigma_z) / 2)
  K <- exp(sigma_z) - 1
  out <- if (identical(scale, "normalised")) {
    nu <- ref_lib * lambda
    list(
      mu = nu,
      sigma = diag(nu * mean(ref_lib / lib_size), length(nu)) +
        outer(nu, nu) * K
    )
  } else {
    m1 <- mean(lib_size)
    m2 <- mean(lib_size^2)
    list(
      mu = m1 * lambda,
      sigma = diag(lambda * m1, length(lambda)) +
        m2 * outer(lambda, lambda) * K +
        (m2 - m1^2) * outer(lambda, lambda)
    )
  }
  genes <- rownames(sigma_z) %||% names(mu_z)
  names(out$mu) <- genes
  dimnames(out$sigma) <- list(genes, genes)
  out
}

#' PLNnetwork fits along a penalty path, one per requested penalty.
#'
#' The offset is the log of the full RNA library `L_c` (as in
#' `scripts/01_01`), not the row sum of the panel genes. `prepare_data()`
#' is called with `offset = "none"` so that it does not add its own TSS
#' offset.
#'
#' @param counts Cells by genes, integer counts of the panel genes.
#' @param lib_size Full library size of each cell (length `nrow(counts)`).
#' @param penalties Penalties to fit. A requested value that the fitted
#'   family does not contain (a truncated path) gives `failed = TRUE`.
#' @param control Optional `PLNmodels::PLNnetwork_param()` list. Default:
#'   unpenalised diagonal, `trace = 0`.
#' @return A list of fits, one per penalty, each with `mu_z`, `sigma_z`,
#'   `penalty`, `n_edges` (`NA` when the fit does not expose it), and
#'   `failed`.
fit_cov_plnnetwork <- function(counts, lib_size, penalties, control = NULL) {
  y <- round(as.matrix(counts))
  if (is.null(rownames(y))) {
    rownames(y) <- paste0("cell_", seq_len(nrow(y)))
  }
  if (length(lib_size) != nrow(y)) {
    stop("`lib_size` must have one entry per cell (row of `counts`).")
  }
  covariates <- data.frame(
    log_library = log(as.numeric(lib_size)),
    row.names = rownames(y)
  )
  dat <- PLNmodels::prepare_data(
    counts = y,
    covariates = covariates,
    offset = "none"
  )
  control <- control %||%
    PLNmodels::PLNnetwork_param(penalize_diagonal = FALSE, trace = 0L)
  family <- tryCatch(
    PLNmodels::PLNnetwork(
      Abundance ~ 1 + offset(log_library),
      data = dat,
      penalties = penalties,
      control = control
    ),
    error = function(e) e
  )
  if (inherits(family, "error")) {
    message("  PLNnetwork failed (", conditionMessage(family), ").")
    return(rep(list(list(failed = TRUE)), length(penalties)))
  }
  available <- family$penalties
  lapply(penalties, function(pen) {
    k <- which(abs(available - pen) <= 1e-8 * max(abs(pen), 1e-12))
    if (length(k) == 0L) {
      return(list(failed = TRUE, penalty = pen))
    }
    model <- PLNmodels::getModel(family, available[[k[[1L]]]])
    sigma_z <- as.matrix(stats::sigma(model))
    mu_z <- as.numeric(stats::coef(model))
    names(mu_z) <- colnames(y)
    dimnames(sigma_z) <- list(colnames(y), colnames(y))
    n_edges <- tryCatch(model$n_edges, error = function(e) NULL)
    list(
      mu_z = mu_z,
      sigma_z = sigma_z,
      penalty = pen,
      n_edges = if (length(n_edges) == 1L) as.integer(n_edges) else NA_integer_,
      failed = !(all(is.finite(mu_z)) && all(is.finite(sigma_z)))
    )
  })
}

#' Penalty grid of a PLNnetwork fit on all cells of one type.
#'
#' PLNmodels builds its own data-driven grid (and its scale changed
#' between versions: covariance scale up to 1.3.2, correlation scale in
#' later versions). Taking the grid from a fitted family keeps the CV
#' grid on whatever scale the installed version uses. Only the candidate
#' values come from all cells; the choice among them is cross-validated.
#'
#' @inheritParams fit_cov_plnnetwork
#' @param n_penalties,min_ratio Passed to
#'   `PLNmodels::PLNnetwork_param()`.
#' @return Numeric vector of penalties, decreasing.
pln_penalty_grid <- function(
  counts,
  lib_size,
  n_penalties = 12L,
  min_ratio = 0.1
) {
  y <- round(as.matrix(counts))
  if (is.null(rownames(y))) {
    rownames(y) <- paste0("cell_", seq_len(nrow(y)))
  }
  dat <- PLNmodels::prepare_data(
    counts = y,
    covariates = data.frame(
      log_library = log(as.numeric(lib_size)),
      row.names = rownames(y)
    ),
    offset = "none"
  )
  family <- PLNmodels::PLNnetwork(
    Abundance ~ 1 + offset(log_library),
    data = dat,
    control = PLNmodels::PLNnetwork_param(
      n_penalties = n_penalties,
      min_ratio = min_ratio,
      penalize_diagonal = FALSE,
      trace = 0L
    )
  )
  sort(family$penalties, decreasing = TRUE)
}

#' Fit function for [cv_select_tuning()] around [fit_cov_plnnetwork()].
#'
#' The CV helpers pass the normalised matrix `X`, which PLN cannot use;
#' this closure fits the raw counts of the same rows and returns, for
#' each penalty, a `moments` function that maps the latent fit to the
#' linear scale of `X` for the test cells' library sizes
#' ([pln_linear_moments()]). Use with `path_fit = TRUE`.
#'
#' @param counts Cells by genes, raw counts of every cell that the CV
#'   rows index.
#' @param lib_size Full library sizes of the same cells.
#' @param ref_lib Reference library used to build `X = counts * ref_lib
#'   / lib_size`.
#' @param control Passed to [fit_cov_plnnetwork()].
#' @return A function `(X_train, penalties, train_rows)` returning a list
#'   of fits.
make_pln_fit_fun <- function(counts, lib_size, ref_lib, control = NULL) {
  force(counts)
  force(lib_size)
  force(ref_lib)
  force(control)
  function(X_train, penalties, train_rows) {
    fits <- fit_cov_plnnetwork(
      counts[train_rows, , drop = FALSE],
      lib_size[train_rows],
      penalties,
      control = control
    )
    lapply(fits, function(fit) {
      if (isTRUE(fit$failed)) {
        return(fit)
      }
      mu_z <- fit$mu_z
      sigma_z <- fit$sigma_z
      fit$moments <- function(test_rows) {
        pln_linear_moments(mu_z, sigma_z, lib_size[test_rows], ref_lib)
      }
      fit
    })
  }
}
