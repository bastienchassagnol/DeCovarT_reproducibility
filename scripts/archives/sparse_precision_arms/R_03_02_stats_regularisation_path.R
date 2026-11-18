#' Theoretical endpoints of a graphical-lasso penalty path.
#'
#' `lambda_max` is the largest off-diagonal absolute covariance, at or
#' above which the estimated graph is empty. `lambda_min` follows the
#' high-dimensional rate `C * sqrt(log(p) / n)`, with a fallback
#' fraction of `lambda_max` when that bound is not smaller than the
#' upper endpoint.
#'
#' @param n Number of independent observations (cells).
#' @param p Number of variables (genes).
#' @param S Optional `p` by `p` sample covariance. If `NULL`,
#'   `lambda_max` is left as `NA`.
#' @param c_lambda Leading constant in the theory bound. Default `1`.
#' @param min_ratio Fallback `lambda_min / lambda_max` when the theory
#'   bound is not strictly below `lambda_max`. Default `0.01`.
#' @return A list with `lambda_min` and `lambda_max`.
lambda_endpoints <- function(
  n,
  p,
  S = NULL,
  c_lambda = 1,
  min_ratio = 0.01
) {
  lambda_max <- if (is.null(S)) {
    NA_real_
  } else {
    max(abs(S[upper.tri(S)]))
  }
  lambda_min <- c_lambda * sqrt(log(p) / n)
  if (is.finite(lambda_max) && lambda_min >= lambda_max) {
    lambda_min <- min_ratio * lambda_max
  }
  list(lambda_min = lambda_min, lambda_max = lambda_max)
}

#' Gaussian log-likelihood of a precision (up to 2 pi constants).
#'
#' @param omega Symmetric positive-definite precision.
#' @param S Sample covariance on the same scale.
#' @param n Sample size.
#' @return Scalar log-likelihood `n/2 * (log det Omega - tr(S Omega))`.
gaussian_precision_loglik <- function(omega, S, n) {
  omega <- as.matrix(omega)
  ld <- as.numeric(determinant(omega, logarithm = TRUE)$modulus)
  n / 2 * (ld - sum(diag(S %*% omega)))
}

#' Extended BIC of a Gaussian graphical model.
#'
#' @param loglik Gaussian log-likelihood from
#'   [gaussian_precision_loglik()].
#' @param n_edges Number of undirected off-diagonal edges.
#' @param n Sample size.
#' @param p Number of genes.
#' @param gamma EBIC extra-penalty in `[0, 1]`. Default `0.5`.
#' @return Scalar EBIC.
ebic_gaussian_ggm <- function(loglik, n_edges, n, p, gamma = 0.5) {
  -2 * loglik + n_edges * log(n) + 4 * gamma * n_edges * log(p)
}

#' Count undirected edges from a precision matrix.
count_undirected_edges <- function(omega, eps = 1e-8) {
  as.integer(sum(abs(omega[upper.tri(omega)]) > eps))
}

#' Frobenius relative gap between two precisions.
precision_path_distance <- function(omega_a, omega_b) {
  denom <- sqrt(sum(as.matrix(omega_a)^2))
  if (!is.finite(denom) || denom < .Machine$double.eps) {
    return(NA_real_)
  }
  sqrt(sum((as.matrix(omega_b) - as.matrix(omega_a))^2)) / denom
}

#' Interval priority for adaptive refinement.
#'
#' @param d2_ebic Second difference of EBIC, or first difference if
#'   the second is unavailable.
#' @param d_edges Change in edge count.
#' @param d_omega Relative Frobenius change of the precision.
#' @param n_possible `p * (p - 1) / 2`.
#' @param w1,w2,w3 Weights on EBIC curvature, support change, and
#'   precision distance. Defaults `1`.
interval_priority <- function(
  d2_ebic,
  d_edges,
  d_omega,
  n_possible,
  w1 = 1,
  w2 = 1,
  w3 = 1
) {
  w1 *
    abs(d2_ebic) +
    w2 * abs(d_edges) / max(n_possible, 1) +
    w3 * (d_omega %||% 0)
}

#' Safeguarded log-scale interval search for an EBIC-minimising penalty.
#'
#' The same search works for graphical lasso, `SILGGM` (its `lambda`),
#' and `PLNnetwork` (its sparsity penalty): only `fit_precision` changes.
#' This is **not** golden-section search. Neighbouring intervals of the
#' current EBIC minimum stay alive so discrete support changes cannot
#' trap the path in a single basin.
#'
#' `fit_precision` must accept a single positive `penalty` and return a
#' list with `omega` (precision) and optionally `sigma` and `loglik`. If
#' `loglik` is missing it is computed from `S`, so `omega` and `S` must
#' be on the same scale. A solve that returns `failed = TRUE` is kept out
#' of the selection and out of the returned path; its penalty is listed
#' in `failed_lambda`.
#'
#' Toy example (do not use as a real GRN):
#' ```
#' set.seed(1)
#' x <- matrix(rnorm(80 * 12), 80, 12)
#' R <- stats::cor(x)
#' ebic_adaptive_path(
#'   fit_precision = make_glasso_solver(R, n = 80),
#'   n = 80, p = 12, S = R, n_solves = 8
#' )
#' ```
#'
#' @param fit_precision Function of one penalty.
#' @param n Number of cells.
#' @param p Number of genes.
#' @param S Sample covariance (or correlation) used for the endpoints and
#'   the Gaussian likelihood, on the scale of `omega`.
#' @param n_solves Maximum number of precision solves. Default `40`.
#' @param gamma EBIC extra-penalty. Default `0.5`.
#' @param c_lambda,min_ratio Passed to [lambda_endpoints()].
#' @param lambda_min,lambda_max Optional overrides of the endpoints.
#' @param w1,w2,w3 Weights in [interval_priority()].
#' @param eps Numerical zero for edges.
#' @param verbose Print the penalty being solved.
#' @return A list with one entry per successful solve: `lambda`, `ebic`,
#'   `n_edges`, `density`, `priority`, `loglik`, `omega`, `sigma` (`NULL`
#'   entries when the solver returns none), plus `selected` (index of the
#'   EBIC minimum), `failed_lambda`, and the endpoints. Stops when every
#'   solve fails.
ebic_adaptive_path <- function(
  fit_precision,
  n,
  p,
  S,
  n_solves = 40L,
  gamma = 0.5,
  c_lambda = 1,
  min_ratio = 0.01,
  lambda_min = NULL,
  lambda_max = NULL,
  w1 = 1,
  w2 = 1,
  w3 = 1,
  eps = 1e-8,
  verbose = TRUE
) {
  n_solves <- as.integer(n_solves)
  ends <- lambda_endpoints(n, p, S, c_lambda, min_ratio)
  if (is.null(lambda_max)) {
    lambda_max <- ends$lambda_max
  }
  if (is.null(lambda_min)) {
    lambda_min <- ends$lambda_min
  }
  if (
    !(is.finite(lambda_min) && is.finite(lambda_max) && lambda_min < lambda_max)
  ) {
    stop("Need finite lambda_min < lambda_max.")
  }

  n_global <- min(6L, n_solves)
  n_possible <- p * (p - 1L) / 2L

  eval_one <- function(penalty) {
    if (verbose) {
      message("  penalty = ", signif(penalty, 4))
    }
    fit <- fit_precision(penalty)
    if (isTRUE(fit$failed)) {
      return(list(
        lambda = penalty,
        omega = NULL,
        sigma = NULL,
        loglik = NA_real_,
        n_edges = NA_integer_,
        density = NA_real_,
        ebic = NA_real_,
        priority = NA_real_,
        failed = TRUE
      ))
    }
    omega <- as.matrix(fit$omega)
    loglik <- fit$loglik %||% gaussian_precision_loglik(omega, S, n)
    n_edges <- count_undirected_edges(omega, eps)
    list(
      lambda = penalty,
      omega = omega,
      sigma = fit$sigma,
      loglik = loglik,
      n_edges = n_edges,
      density = n_edges / n_possible,
      ebic = ebic_gaussian_ggm(loglik, n_edges, n, p, gamma),
      priority = NA_real_,
      failed = FALSE
    )
  }

  anchors <- exp(seq(log(lambda_max), log(lambda_min), length.out = n_global))
  path <- lapply(anchors, eval_one)

  n_basin <- min(n_solves - length(path), max(0L, round(0.4 * n_solves)))
  n_adjacent <- min(
    n_solves - length(path) - n_basin,
    max(0L, round(0.2 * n_solves))
  )

  score_intervals <- function(path_sorted) {
    m <- length(path_sorted)
    scores <- rep(NA_real_, m - 1L)
    for (k in seq_len(m - 1L)) {
      d_ebic <- path_sorted[[k + 1L]]$ebic - path_sorted[[k]]$ebic
      d2 <- if (k > 1L) {
        d_ebic - (path_sorted[[k]]$ebic - path_sorted[[k - 1L]]$ebic)
      } else {
        d_ebic
      }
      if (path_sorted[[k]]$failed || path_sorted[[k + 1L]]$failed) {
        next
      }
      d_edges <- path_sorted[[k + 1L]]$n_edges - path_sorted[[k]]$n_edges
      d_omega <- precision_path_distance(
        path_sorted[[k]]$omega,
        path_sorted[[k + 1L]]$omega
      )
      scores[[k]] <- interval_priority(
        d2_ebic = d2,
        d_edges = d_edges,
        d_omega = d_omega,
        n_possible = n_possible,
        w1 = w1,
        w2 = w2,
        w3 = w3
      )
    }
    scores
  }

  phase_budget <- c(
    rep("basin", n_basin),
    rep("adjacent", n_adjacent),
    rep("verify", max(0L, n_solves - n_global - n_basin - n_adjacent))
  )

  for (phase in phase_budget) {
    if (length(path) >= n_solves) {
      break
    }
    ord <- order(vapply(path, `[[`, numeric(1), "lambda"), decreasing = TRUE)
    path <- path[ord]
    scores <- score_intervals(path)
    ebic <- vapply(path, `[[`, numeric(1), "ebic")
    if (all(is.na(ebic)) || all(is.na(scores))) {
      break
    }
    i_min <- which.min(ebic)
    eligible <- seq_along(scores)
    if (identical(phase, "basin")) {
      eligible <- unique(pmax(1L, pmin(length(scores), (i_min - 1L):i_min)))
    } else if (identical(phase, "adjacent")) {
      eligible <- unique(pmax(
        1L,
        pmin(length(scores), (i_min - 2L):(i_min + 1L))
      ))
    }
    eligible <- eligible[!is.na(scores[eligible])]
    if (length(eligible) == 0L) {
      eligible <- which(!is.na(scores))
    }
    pick <- eligible[[which.max(scores[eligible])]]
    lambda_new <- sqrt(path[[pick]]$lambda * path[[pick + 1L]]$lambda)
    if (
      any(
        abs(vapply(path, `[[`, numeric(1), "lambda") - lambda_new) <
          .Machine$double.eps^0.5
      )
    ) {
      next
    }
    added <- eval_one(lambda_new)
    added$priority <- scores[[pick]]
    path[[length(path) + 1L]] <- added
  }

  failed <- vapply(path, `[[`, logical(1), "failed")
  failed_lambda <- sort(
    vapply(path[failed], `[[`, numeric(1), "lambda"),
    decreasing = TRUE
  )
  path <- path[!failed]
  if (length(path) == 0L) {
    stop("All ", length(failed), " penalty solves failed.")
  }
  ord <- order(vapply(path, `[[`, numeric(1), "lambda"), decreasing = TRUE)
  path <- path[ord]
  ebic <- vapply(path, `[[`, numeric(1), "ebic")
  list(
    lambda = vapply(path, `[[`, numeric(1), "lambda"),
    ebic = ebic,
    n_edges = vapply(path, `[[`, integer(1), "n_edges"),
    density = vapply(path, `[[`, numeric(1), "density"),
    priority = vapply(path, `[[`, numeric(1), "priority"),
    loglik = vapply(path, `[[`, numeric(1), "loglik"),
    omega = lapply(path, `[[`, "omega"),
    sigma = lapply(path, `[[`, "sigma"),
    failed_lambda = failed_lambda,
    selected = which.min(ebic),
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    n = n,
    p = p,
    gamma = gamma
  )
}

#' Graphical-lasso solver for [ebic_adaptive_path()] via `glasso`.
#'
#' Fits the correlation matrix `R`, so the penalty is on the same scale
#' for every gene and the `lambda` endpoints of [lambda_endpoints()]
#' apply unchanged. `penalize.diagonal = FALSE` keeps
#' `sigma[g, g] = R[g, g] = 1`; the `glasso` default (`TRUE`) would add
#' `lambda` to every variance. The caller rescales `omega` and `sigma`
#' back to the data scale.
#'
#' A failed solve is reported as `failed = TRUE`, never replaced by a
#' diagonal precision: a diagonal stand-in would be scored by EBIC as
#' an empty graph and could be selected.
#'
#' @param R Correlation matrix of the genes, `p` by `p`.
#' @param n Number of cells used to compute `R`.
#' @param thr,maxit Passed to [glasso::glasso()].
#' @return A function of one `penalty` returning `omega`, `sigma`
#'   (the fitted covariance `W`, the inverse of `omega`), or
#'   `failed = TRUE`.
make_glasso_solver <- function(R, n, thr = 1e-4, maxit = 1e4) {
  function(penalty) {
    fit <- tryCatch(
      glasso::glasso(
        R,
        rho = penalty,
        nobs = n,
        thr = thr,
        maxit = maxit,
        penalize.diagonal = FALSE
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      message("  glasso failed (", conditionMessage(fit), ").")
      return(list(failed = TRUE))
    }
    omega <- (fit$wi + t(fit$wi)) / 2
    sigma <- (fit$w + t(fit$w)) / 2
    if (fit$errflag != 0L || !all(is.finite(omega)) || !all(is.finite(sigma))) {
      message("  glasso returned no finite estimate at this penalty.")
      return(list(failed = TRUE))
    }
    dimnames(omega) <- dimnames(R)
    dimnames(sigma) <- dimnames(R)
    list(omega = omega, sigma = sigma, failed = FALSE)
  }
}

#' Rescale a correlation-scale fit to the scale of the data.
#'
#' With `D = diag(sd^2)`, the covariance is `D^{1/2} W D^{1/2}` and the
#' precision `D^{-1/2} Omega D^{-1/2}`. Zeros of the precision, and
#' partial correlations, are unchanged.
#'
#' @param omega,sigma Precision and covariance of the standardised genes.
#' @param sd Standard deviations of the genes on the data scale.
#' @return A list with `omega` and `sigma` on the data scale.
rescale_correlation_fit <- function(omega, sigma, sd) {
  scale_outer <- outer(sd, sd)
  list(omega = omega / scale_outer, sigma = sigma * scale_outer)
}

#' Stop when a covariance is not a non-degenerate Gaussian scale.
#'
#' A non-positive diagonal variance makes the Hellinger denominator
#' undefined and sends `MixSim::overlap()` into a long failing solve.
#'
#' @param sigma Square covariance.
#' @param label Which covariance, used in the error message.
#' @return `TRUE`, invisibly, when every variance is positive.
require_positive_variances <- function(sigma, label) {
  variances <- diag(sigma)
  if (
    any(!is.finite(sigma)) || any(!is.finite(variances)) || any(variances <= 0)
  ) {
    stop(label, " covariance has a non-positive variance.")
  }
  invisible(TRUE)
}

#' Hellinger distance between two non-degenerate Gaussians.
gaussian_hellinger <- function(mu1, sigma1, mu2, sigma2) {
  require_positive_variances(sigma1, "First")
  require_positive_variances(sigma2, "Second")
  mu1 <- as.numeric(mu1)
  mu2 <- as.numeric(mu2)
  sigma_bar <- (sigma1 + sigma2) / 2
  ld1 <- as.numeric(determinant(sigma1, logarithm = TRUE)$modulus)
  ld2 <- as.numeric(determinant(sigma2, logarithm = TRUE)$modulus)
  ldb <- as.numeric(determinant(sigma_bar, logarithm = TRUE)$modulus)
  dmu <- mu1 - mu2
  quad <- as.numeric(crossprod(dmu, solve(sigma_bar, dmu)))
  bhattacharyya <- 0.125 * quad + 0.5 * ldb - 0.25 * (ld1 + ld2)
  h2 <- 1 - exp(-bhattacharyya)
  sqrt(pmax(0, pmin(1, h2)))
}

#' MixSim pairwise overlap of two Gaussians.
#'
#' @return Sum of the two misclassification probabilities.
gaussian_mixsim_overlap <- function(mu1, sigma1, n1, mu2, sigma2, n2) {
  require_positive_variances(sigma1, "First")
  require_positive_variances(sigma2, "Second")
  g <- length(mu1)
  s <- array(0, dim = c(g, g, 2L))
  s[,, 1L] <- sigma1
  s[,, 2L] <- sigma2
  ov <- MixSim::overlap(
    Pi = c(n1, n2) / (n1 + n2),
    Mu = rbind(mu1, mu2),
    S = s
  )
  as.numeric(ov$OmegaMap[1L, 2L] + ov$OmegaMap[2L, 1L])
}
