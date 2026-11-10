#' Cross-validated selection and comparison of within-type covariances.
#'
#' Every estimator in this file returns a mean `mu` and a covariance
#' `sigma` on the linear scale of DeCovarT's convolution. By default
#' that scale is the raw per-cell UMI count of the panel genes, with no
#' per-cell library division (`normalise = "none"` in the drivers): the
#' division removes the type-level transcriptome size and the biological
#' size variation that the bulk accumulates. The earlier scale,
#' `x_cg = y_cg / s_c` with `s_c = L_c / L_ref` (`L_c` the full RNA
#' library of cell `c`, `L_ref` its median over the cells of one time
#' point), is kept as a sensitivity variant (`normalise = "library"`).
#' All estimators are scored by the same loss: the held-out Gaussian
#' log-score of the test cells on that scale.
#'
#' A denoised fit (`preprocess = "denoised"`) exports the biological
#' covariance `sigma`, without the Poisson diagonal and with part or all
#' of the size term removed ([denoise_covariance()]). The test cells
#' still carry those terms, so such a fit also returns `score_sigma`, the
#' count covariance rebuilt from `sigma` ([reconstruct_count_covariance()]),
#' and [score_fit()] scores `score_sigma` when it is present. Every
#' `kappa` rebuilds the same count moments up to regularisation, so the
#' score cannot choose `kappa`; it stays a bracket.
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
#' Requires `R/02_00_utils.R` (`%||%`) and
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
#' @param fit List with `mu` and `sigma` (and, for a denoised fit,
#'   `score_sigma`), or with `moments`, a function of the test-row
#'   indices.
#' @param test_rows Indices of the test cells in the caller's full data.
#' @return A list with `mu`, `sigma` (the exported covariance) and
#'   `score_sigma` (`NULL` unless the fit is denoised).
fit_moments <- function(fit, test_rows) {
  if (is.function(fit$moments)) {
    fit$moments(test_rows)
  } else {
    list(mu = fit$mu, sigma = fit$sigma, score_sigma = fit$score_sigma)
  }
}

#' Held-out log-score of one fit, `NA` for a failed fit.
#'
#' A denoised fit is scored through `score_sigma`, the count covariance
#' rebuilt from its exported `sigma`; every other fit through `sigma`.
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
  heldout_gaussian_logscore(m$mu, m$score_sigma %||% m$sigma, X_test)
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
# Denoised covariance: moment correction, then a PSD projection ----
# ==========================================================================

#' Denoised within-type covariance of raw UMI counts.
#'
#' Moment estimator of the covariance of the captured molecule content
#' (eq-denoised-within of vignette 05-1). With `V` the unbiased sample
#' covariance of the `n` cells and `xbar` their mean,
#' `V - diag(xbar) - kappa * (V - diag(xbar) + xbar xbar' - V / n)`.
#' The first subtraction removes the Poisson measurement diagonal. The
#' `kappa` term removes the technical share of the size term, with
#' `kappa = nu_tech / (1 + nu_tech)`. Without spike-ins `kappa` is not
#' identified and is bracketed: `kappa = 0` keeps the whole size term
#' (all size variation biological); `kappa = nu / (1 + nu)`, with `nu`
#' from [extrinsic_noise_nu()], removes it (all size variation
#' technical, which matches library normalisation for the covariance).
#'
#' The estimate can be indefinite. With `psd = "clip"` its eigenvalues
#' below `eig_floor` are raised to `eig_floor`; with `eig_floor = 0`
#' this is the nearest positive semi-definite matrix in Frobenius norm.
#' The eigen-decomposition costs of order `p^3` operations.
#'
#' @param X Cells by genes, raw counts of one type at one time.
#' @param kappa Share of the size term removed, in `[0, 1)`. Default `0`.
#' @param psd `"clip"` (default) or `"none"`.
#' @param eig_floor Smallest eigenvalue kept by `psd = "clip"`.
#'   Default `0`.
#' @return A list with `sigma` (denoised covariance, `p` by `p`), `mu`
#'   (`xbar`), `variance` (diagonal before the projection; a gene with
#'   `variance <= 0` is a Poisson gene at this depth), `n_clipped`
#'   (eigenvalues raised), `kappa`, `n`.
denoise_covariance <- function(
  X,
  kappa = 0,
  psd = c("clip", "none"),
  eig_floor = 0
) {
  psd <- match.arg(psd)
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < 2L) {
    stop("`X` needs at least two cells.")
  }
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) ||
      kappa < 0 || kappa >= 1) {
    stop("`kappa` must be one number in [0, 1); got ", kappa, ".")
  }
  xbar <- colMeans(X)
  V <- stats::cov(X)
  A <- V - diag(xbar, length(xbar))
  sigma <- A - kappa * (A + tcrossprod(xbar) - V / n)
  sigma <- (sigma + t(sigma)) / 2
  variance <- diag(sigma)
  n_clipped <- 0L
  if (identical(psd, "clip")) {
    eg <- eigen(sigma, symmetric = TRUE)
    low <- eg$values < eig_floor
    n_clipped <- sum(low)
    if (n_clipped > 0L) {
      values <- pmax(eg$values, eig_floor)
      sigma <- eg$vectors %*% (values * t(eg$vectors))
      sigma <- (sigma + t(sigma)) / 2
    }
  }
  genes <- colnames(X)
  dimnames(sigma) <- list(genes, genes)
  names(xbar) <- genes
  names(variance) <- genes
  list(
    sigma = sigma,
    mu = xbar,
    variance = variance,
    n_clipped = as.integer(n_clipped),
    kappa = kappa,
    n = n
  )
}

#' Count covariance rebuilt from a denoised covariance.
#'
#' Inverts the moment correction at reference depth, so that a denoised
#' fit can be scored on raw test cells:
#' `sigma_bio + diag(mu) + kappa / (1 - kappa) * (sigma_bio + mu mu')`,
#' with `nu_tech = kappa / (1 - kappa)`. With `kappa = 0` it is
#' `sigma_bio + diag(mu)`. The finite-sample term `V / n` of
#' [denoise_covariance()] is not added back: this is the population
#' identity, evaluated at the estimates.
#'
#' @param sigma_bio Denoised (and regularised) covariance.
#' @param mu Mean counts of the same genes.
#' @param kappa The `kappa` used to denoise.
#' @return The rebuilt covariance, `p` by `p`.
reconstruct_count_covariance <- function(sigma_bio, mu, kappa = 0) {
  mu <- as.numeric(mu)
  out <- sigma_bio + diag(mu, length(mu)) +
    kappa / (1 - kappa) * (sigma_bio + tcrossprod(mu))
  out <- (out + t(out)) / 2
  dimnames(out) <- dimnames(sigma_bio)
  out
}

#' Means, standard deviations and correlation of a denoised covariance.
#'
#' The counterpart of [standardise_genes()] for [denoise_covariance()].
#' A gene whose denoised variance (before the projection) is not
#' positive is a Poisson gene: its biological variance is set to
#' `var_floor`, its correlations to zero, and it is listed in
#' `floored_genes`. The floor is a numerical guard, not an estimate.
#'
#' @param den Output of [denoise_covariance()].
#' @param var_floor Variance floor. Default `NULL`: `1e-4` times the
#'   median positive denoised variance (or `1` when none is positive).
#' @return A list with `mu`, `sd`, `R`, `floored` (logical),
#'   `floored_genes`, `n`.
standardise_denoised <- function(den, var_floor = NULL) {
  v <- diag(den$sigma)
  floored <- !is.finite(den$variance) | den$variance <= 0 |
    !is.finite(v) | v <= 0
  if (is.null(var_floor)) {
    var_floor <- if (any(!floored)) 1e-4 * stats::median(v[!floored]) else 1
  }
  v[floored] <- 0
  sd <- sqrt(pmax(v, var_floor))
  p <- length(sd)
  R <- diag(1, p)
  keep <- !floored
  if (sum(keep) >= 2L) {
    s_keep <- sqrt(v[keep])
    R[keep, keep] <- den$sigma[keep, keep, drop = FALSE] / outer(s_keep, s_keep)
  }
  R <- (R + t(R)) / 2
  diag(R) <- 1
  genes <- colnames(den$sigma)
  dimnames(R) <- list(genes, genes)
  names(sd) <- genes
  list(
    mu = den$mu,
    sd = sd,
    R = R,
    floored = floored,
    floored_genes = genes[floored],
    n = den$n
  )
}

#' Extrinsic noise `nu` of one cell type, from the normalised covariances.
#'
#' After Fang and Pachter (2025): inside one type, the normalised
#' covariance `cov(x_g, x_h) / (mean_g * mean_h)` of two genes that are
#' intrinsically uncorrelated equals `nu = CV^2(s_c)`, the squared
#' coefficient of variation of the cell size factor. The estimate is the
#' mode of these values over gene pairs (maximum of a Gaussian kernel
#' density). The estimate assumes that most gene pairs are intrinsically
#' uncorrelated; co-regulated modules shift the distribution to the
#' right and bias `nu` upwards. Pairs of genes with small means have a
#' noisy ratio, so only genes with `mean >= min_mean` enter, at most
#' `max_genes` of them (largest means first); `min_mean` is a choice of
#' this repository, not of the paper.
#'
#' `nu` mixes the biological and technical size factors,
#' `1 + nu = (1 + nu_bio)(1 + nu_tech)`. Without spike-ins only the
#' upper end of the `kappa` bracket follows: `kappa_max = nu / (1 + nu)`.
#'
#' @param X Cells by genes, raw counts of one type at one time.
#' @param min_mean Smallest mean count of a gene that enters. Default `1`.
#' @param max_genes Largest number of genes. Default `2000`.
#' @param trim Share of pairs trimmed in each tail before the density,
#'   which keeps the bandwidth from being set by outliers. Default `0.01`.
#' @param adjust Bandwidth multiplier passed to [stats::density()].
#' @return A list with `nu`, `kappa_max`, `n_genes`, `n_pairs`,
#'   `bandwidth`.
extrinsic_noise_nu <- function(
  X,
  min_mean = 1,
  max_genes = 2000L,
  trim = 0.01,
  adjust = 1
) {
  X <- as.matrix(X)
  m <- colMeans(X)
  keep <- which(is.finite(m) & m >= min_mean)
  if (length(keep) > max_genes) {
    keep <- keep[order(m[keep], decreasing = TRUE)][seq_len(max_genes)]
  }
  if (length(keep) < 3L) {
    stop(
      "Fewer than 3 genes have a mean count >= ", min_mean,
      "; lower `min_mean` to estimate nu."
    )
  }
  V <- stats::cov(X[, keep, drop = FALSE])
  eta <- V / tcrossprod(m[keep])
  vals <- eta[upper.tri(eta)]
  vals <- vals[is.finite(vals)]
  bounds <- stats::quantile(vals, c(trim, 1 - trim), names = FALSE)
  vals_mid <- vals[vals >= bounds[[1L]] & vals <= bounds[[2L]]]
  dens <- stats::density(vals_mid, adjust = adjust, n = 2048L)
  nu <- dens$x[[which.max(dens$y)]]
  if (nu <= 0) {
    warning(
      "The mode of the normalised covariances is not positive (nu = ",
      signif(nu, 3), "); kappa_max is set to 0.",
      call. = FALSE
    )
  }
  nu_pos <- max(nu, 0)
  list(
    nu = nu,
    kappa_max = nu_pos / (1 + nu_pos),
    n_genes = length(keep),
    n_pairs = length(vals),
    bandwidth = dens$bw
  )
}

#' Attach the denoised export and its scoring covariance to a fit.
#'
#' @param fit A fit with `mu` and `sigma` on the denoised scale.
#' @param st Output of [standardise_denoised()].
#' @param den Output of [denoise_covariance()].
#' @return `fit` with `score_sigma`, `preprocess`, `kappa`, `n_clipped`
#'   and `floored_genes` added.
finish_denoised_fit <- function(fit, st, den) {
  fit$score_sigma <- reconstruct_count_covariance(fit$sigma, fit$mu, den$kappa)
  fit$preprocess <- "denoised"
  fit$kappa <- den$kappa
  fit$n_clipped <- den$n_clipped
  fit$floored_genes <- st$floored_genes
  fit
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
#' With `preprocess = "denoised"` the same target D is applied to the
#' denoised covariance of [denoise_covariance()] (positive semi-definite
#' after clipping): `D^{1/2} {(1 - lambda) R + lambda I} D^{1/2}`, where
#' `R` and `D` are the correlation and the variances of the denoised
#' matrix (as `corpcor` does on the sample covariance). Poisson genes
#' (denoised variance `<= 0`) get the floor of [standardise_denoised()].
#' The returned `sigma` is the denoised estimate (the export);
#' `score_sigma` is the count covariance used by [score_fit()].
#'
#' @param X Cells by genes (training cells), raw counts for
#'   `preprocess = "denoised"`.
#' @param lambda Correlation shrinkage intensity in `[0, 1]`; a vector
#'   when `as_path = TRUE`.
#' @param train_rows Unused; kept for the common fit signature.
#' @param var_floor Passed to [standardise_genes()] or
#'   [standardise_denoised()].
#' @param preprocess `"raw"` (default) or `"denoised"`.
#' @param kappa Passed to [denoise_covariance()]. Default `0`.
#' @param as_path If `TRUE`, `lambda` is a grid and a list of fits is
#'   returned, with the preprocessing done once (use with
#'   `path_fit = TRUE` in [cv_select_tuning()]). Default `FALSE`.
#' @return A list with `mu`, `sigma`, `failed` (and, when denoised,
#'   `score_sigma`, `kappa`, `n_clipped`, `floored_genes`); a list of
#'   such fits when `as_path = TRUE`.
fit_cov_shrinkage <- function(
  X,
  lambda,
  train_rows = NULL,
  var_floor = NULL,
  preprocess = c("raw", "denoised"),
  kappa = 0,
  as_path = FALSE
) {
  preprocess <- match.arg(preprocess)
  X <- as.matrix(X)
  if (identical(preprocess, "denoised")) {
    den <- tryCatch(denoise_covariance(X, kappa = kappa), error = function(e) e)
    if (inherits(den, "error")) {
      message("  denoising failed (", conditionMessage(den), ").")
      failed <- list(failed = TRUE)
      return(if (as_path) rep(list(failed), length(lambda)) else failed)
    }
    st <- standardise_denoised(den, var_floor)
    off <- st$R
    diag(off) <- 0
    one_fit <- function(lam) {
      R_shrunk <- (1 - lam) * off
      diag(R_shrunk) <- 1
      sigma <- R_shrunk * outer(st$sd, st$sd)
      dimnames(sigma) <- list(colnames(X), colnames(X))
      finish_denoised_fit(
        list(mu = st$mu, sigma = sigma, failed = FALSE),
        st,
        den
      )
    }
    return(if (as_path) lapply(lambda, one_fit) else one_fit(lambda))
  }
  if (as_path) {
    return(lapply(lambda, function(lam) {
      fit_cov_shrinkage(X, lam, train_rows, var_floor)
    }))
  }
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
#' With `preprocess = "denoised"` the input is the correlation of the
#' denoised covariance ([denoise_covariance()], clipped to positive
#' semi-definite) and the rescaling uses its standard deviations; Poisson
#' genes get the floor of [standardise_denoised()]. `sigma` and `omega`
#' are then the denoised estimate and its precision; `score_sigma` is the
#' count covariance used by [score_fit()]. The denoising is redone for
#' every `lambda`, which costs one eigen-decomposition per call.
#'
#' @param X Cells by genes (training cells).
#' @param lambda Lasso penalty on the correlation scale.
#' @param train_rows Unused; kept for the common fit signature.
#' @param var_floor Passed to [standardise_genes()] or
#'   [standardise_denoised()].
#' @param thr,maxit Passed to [glasso::glasso()].
#' @param preprocess `"raw"` (default) or `"denoised"`.
#' @param kappa Passed to [denoise_covariance()]. Default `0`.
#' @return A list with `mu`, `sigma`, `omega`, `n_edges`, `failed` (and,
#'   when denoised, `score_sigma`, `kappa`, `n_clipped`, `floored_genes`).
fit_cov_glasso <- function(
  X,
  lambda,
  train_rows = NULL,
  var_floor = NULL,
  thr = 1e-4,
  maxit = 1e4,
  preprocess = c("raw", "denoised"),
  kappa = 0
) {
  preprocess <- match.arg(preprocess)
  den <- NULL
  st <- if (identical(preprocess, "denoised")) {
    den <- tryCatch(denoise_covariance(X, kappa = kappa), error = function(e) e)
    if (inherits(den, "error")) {
      message("  denoising failed (", conditionMessage(den), ").")
      return(list(failed = TRUE))
    }
    standardise_denoised(den, var_floor)
  } else {
    standardise_genes(X, var_floor)
  }
  fit <- make_glasso_solver(st$R, st$n, thr = thr, maxit = maxit)(lambda)
  if (isTRUE(fit$failed)) {
    return(list(failed = TRUE))
  }
  scaled <- rescale_correlation_fit(fit$omega, fit$sigma, st$sd)
  out <- list(
    mu = st$mu,
    sigma = scaled$sigma,
    omega = scaled$omega,
    n_edges = count_undirected_edges(fit$omega),
    failed = FALSE
  )
  if (is.null(den)) out else finish_denoised_fit(out, st, den)
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
#' With `preprocess = "denoised"` only the refit changes: ranks cannot be
#' denoised, so the support still comes from the nonparanormal scores of
#' `X`, and [refit_covariance_on_support()] receives the correlation and
#' standard deviations of the denoised covariance
#' ([denoise_covariance()], [standardise_denoised()]). `sigma` and
#' `omega` are then the denoised estimate and its precision;
#' `score_sigma` is the count covariance used by [score_fit()].
#'
#' @param X Cells by genes (training cells), linear scale.
#' @param alpha Vector of FDR levels.
#' @param train_rows Unused; kept for the common fit signature.
#' @param method SILGGM method. Default `"D-S_GL"` (de-sparsified
#'   graphical lasso).
#' @param npn_func Passed to `huge::huge.npn(npn.func = ...)`. Default
#'   `"shrinkage"`.
#' @param rho_refit Passed to [refit_covariance_on_support()].
#' @param var_floor Passed to [standardise_genes()] or
#'   [standardise_denoised()].
#' @param preprocess `"raw"` (default) or `"denoised"`, for the refit.
#' @param kappa Passed to [denoise_covariance()]. Default `0`.
#' @return A list of fits, one per `alpha`, each with `mu`, `sigma`,
#'   `omega`, `n_edges`, `failed` (and, when denoised, `score_sigma`,
#'   `kappa`, `n_clipped`, `floored_genes`).
fit_cov_silggm_refit <- function(
  X,
  alpha,
  train_rows = NULL,
  method = "D-S_GL",
  npn_func = "shrinkage",
  rho_refit = 1e-4,
  var_floor = NULL,
  preprocess = c("raw", "denoised"),
  kappa = 0
) {
  preprocess <- match.arg(preprocess)
  X <- as.matrix(X)
  st <- standardise_genes(X, var_floor)
  vary <- !st$constant
  p <- ncol(X)
  failed_all <- rep(list(list(failed = TRUE)), length(alpha))
  if (sum(vary) < 2L) {
    return(failed_all)
  }
  # Scale of the refit: the linear-scale sample correlation, or that of
  # the denoised covariance. The support below uses the ranks of X either
  # way.
  den <- NULL
  st_refit <- st
  if (identical(preprocess, "denoised")) {
    den <- tryCatch(denoise_covariance(X, kappa = kappa), error = function(e) e)
    if (inherits(den, "error")) {
      message("  denoising failed (", conditionMessage(den), ").")
      return(failed_all)
    }
    st_refit <- standardise_denoised(den, var_floor)
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
      st_refit$R,
      adjacency,
      st_refit$sd,
      rho = rho_refit,
      nobs = st_refit$n
    )
    if (isTRUE(refit$failed)) {
      return(refit)
    }
    out <- c(list(mu = st_refit$mu), refit)
    if (is.null(den)) out else finish_denoised_fit(out, st_refit, den)
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
#' The offset `log L_c` removes every cell-level size variation from
#' `Sigma_z`, biological size included. Two arguments decide what the
#' export puts back (the `kappa` bracket of [denoise_covariance()]):
#' - `size_term = "keep"` (default) averages over the libraries of the
#'   cells, as above, so the depth variation is kept (`kappa = 0` with
#'   `scale = "counts"`); `size_term = "drop"` gives every cell the mean
#'   library, so `var(L) = 0` and `mean(L^2) = mean(L)^2` (all size
#'   variation treated as technical);
#' - `poisson_diag = TRUE` (default) keeps the Poisson measurement
#'   diagonal; `FALSE` omits it, for a denoised export. Scoring on raw
#'   test cells needs `TRUE`.
#'
#' @param mu_z Latent mean (length `p`), on the log scale per unit of
#'   library (the PLN intercept with offset `log L_c`).
#' @param sigma_z Latent covariance, `p` by `p`.
#' @param lib_size Library sizes `L_c` of the cells the moments are for.
#' @param ref_lib Reference library `L_ref`. Default the median of
#'   `lib_size`; pass the value used to normalise the data.
#' @param scale `"normalised"` (default) or `"counts"`.
#' @param size_term `"keep"` (default) or `"drop"`.
#' @param poisson_diag Keep the Poisson diagonal. Default `TRUE`.
#' @return A list with `mu` and `sigma`.
pln_linear_moments <- function(
  mu_z,
  sigma_z,
  lib_size,
  ref_lib = stats::median(lib_size),
  scale = c("normalised", "counts"),
  size_term = c("keep", "drop"),
  poisson_diag = TRUE
) {
  scale <- match.arg(scale)
  size_term <- match.arg(size_term)
  sigma_z <- as.matrix(sigma_z)
  mu_z <- as.numeric(mu_z)
  if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    stop("`lib_size` must be positive and finite.")
  }
  # `ref_lib` is evaluated before the libraries are equalised, so its
  # default stays the median of the libraries passed in.
  force(ref_lib)
  if (identical(size_term, "drop")) {
    lib_size <- rep(mean(lib_size), length(lib_size))
  }
  lambda <- exp(mu_z + diag(sigma_z) / 2)
  K <- exp(sigma_z) - 1
  diag_weight <- if (isTRUE(poisson_diag)) 1 else 0
  out <- if (identical(scale, "normalised")) {
    nu <- ref_lib * lambda
    list(
      mu = nu,
      sigma = diag_weight *
        diag(nu * mean(ref_lib / lib_size), length(nu)) +
        outer(nu, nu) * K
    )
  } else {
    m1 <- mean(lib_size)
    m2 <- mean(lib_size^2)
    list(
      mu = m1 * lambda,
      sigma = diag_weight * diag(lambda * m1, length(lambda)) +
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
#' `scripts/legacy/02_90_stats_prepare_panel.R`), not the row sum of
#' the panel genes. `prepare_data()`
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
#' The CV helpers pass the scoring matrix `X` (raw or library-divided
#' counts), which PLN does not fit directly; this closure fits the raw
#' counts of the same rows with the full-library offset and returns, for
#' each penalty, a `moments` function that maps the latent fit to the
#' linear scale of `X` for the test cells' library sizes
#' ([pln_linear_moments()]). Use with `path_fit = TRUE`.
#'
#' @param counts Cells by genes, raw counts of every cell that the CV
#'   rows index.
#' @param lib_size Full library sizes of the same cells.
#' @param ref_lib Reference library used to build `X = counts * ref_lib
#'   / lib_size` (ignored when `scale = "counts"`).
#' @param control Passed to [fit_cov_plnnetwork()].
#' @param scale Scale of `X`, passed to [pln_linear_moments()]:
#'   `"normalised"` (default, library-divided `X`) or `"counts"` (raw
#'   counts, `normalise = "none"` in the drivers).
#' @return A function `(X_train, penalties, train_rows)` returning a list
#'   of fits.
make_pln_fit_fun <- function(
  counts,
  lib_size,
  ref_lib,
  control = NULL,
  scale = c("normalised", "counts")
) {
  scale <- match.arg(scale)
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
        pln_linear_moments(
          mu_z,
          sigma_z,
          lib_size[test_rows],
          ref_lib,
          scale = scale
        )
      }
      fit
    })
  }
}

# ==========================================================================
# Inputs: scored counts restricted to the gated panel ----
# ==========================================================================

#' Panel counts of one bulk-matched time.
#'
#' Reads one slot of `scored_slots.rds` (tier (i) of the feature
#' selection: raw UMI counts of the scored cells, genes by cells, with
#' `cell_type`, `batch` and the full RNA `library_size`) and keeps the
#' genes of `selected_<time>.csv` (tier (ii): every gene that passes the
#' permuted Welch G at BH 0.05, no cap).
#'
#' @param slots Output of `readRDS("output/feature_selection/tier_01/scored_slots.rds")`.
#' @param selected_dir Directory of `selected_<time>.csv`
#'   (`output/feature_selection/tier_02`).
#' @param time_level One of `slots$time_levels`.
#' @return A list with `y` (cells by genes, dense, raw counts of the
#'   panel genes, not rounded), `cell_type`, `batch`, `library_size`
#'   (full library, one value per cell), `cells`, `genes`.
read_panel_slot <- function(slots, selected_dir, time_level) {
  slot <- slots$slots[[time_level]]
  if (is.null(slot)) {
    stop("No slot for time ", time_level, " in scored_slots.rds.")
  }
  if (methods::is(slot$counts, "Matrix")) {
    # Methods for `[` and `as.matrix()` on the sparse counts.
    requireNamespace("Matrix", quietly = TRUE)
  }
  selected_path <- file.path(
    selected_dir,
    paste0("selected_", time_level, ".csv")
  )
  if (!file.exists(selected_path)) {
    stop(
      "Missing ", selected_path,
      ". Run scripts/02_02_stats_cross_cell_consistency.R first."
    )
  }
  genes <- utils::read.csv(
    selected_path,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )$gene
  missing <- setdiff(genes, rownames(slot$counts))
  if (length(missing) > 0L) {
    stop(
      length(missing), " panel genes of ", time_level,
      " are absent from the slot counts (first: ", missing[[1L]], ")."
    )
  }
  cells <- slot$cells %||% colnames(slot$counts)
  y <- t(as.matrix(slot$counts[genes, , drop = FALSE]))
  dimnames(y) <- list(cells, genes)
  lib <- as.numeric(slot$library_size)
  if (length(lib) != nrow(y) || any(!is.finite(lib)) || any(lib <= 0)) {
    stop("Invalid library sizes in the ", time_level, " slot.")
  }
  names(lib) <- cells
  list(
    y = y,
    cell_type = as.character(slot$cell_type),
    batch = as.character(slot$batch),
    library_size = lib,
    cells = cells,
    genes = genes
  )
}
