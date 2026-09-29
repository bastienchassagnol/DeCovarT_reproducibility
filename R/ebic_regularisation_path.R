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
  w1 * abs(d2_ebic) +
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
#' list with `omega` (precision) and optionally `loglik`. If `loglik` is
#' missing it is computed from `S`.
#'
#' Toy example (do not use as a real GRN):
#' ```
#' set.seed(1)
#' x <- matrix(rnorm(80 * 12), 80, 12)
#' S <- stats::cov(x)
#' fit <- function(penalty) {
#'   fit <- huge::huge(x, method = "glasso", lambda = penalty, verbose = FALSE)
#'   list(omega = as.matrix(fit$icov[[1]]))
#' }
#' ebic_adaptive_path(fit_precision = fit, n = 80, p = 12, S = S, n_solves = 8)
#' ```
#'
#' @param fit_precision Function of one penalty.
#' @param n Number of cells.
#' @param p Number of genes.
#' @param S Sample covariance used for the Gaussian likelihood.
#' @param n_solves Maximum number of precision solves. Default `40`.
#' @param gamma EBIC extra-penalty. Default `0.5`.
#' @param c_lambda,min_ratio Passed to [lambda_endpoints()].
#' @param lambda_min,lambda_max Optional overrides of the endpoints.
#' @param w1,w2,w3 Weights in [interval_priority()].
#' @param eps Numerical zero for edges.
#' @param verbose Print the penalty being solved.
#' @return A list with one entry per solve: `lambda`, `ebic`, `n_edges`,
#'   `density`, `priority`, `omega`, plus `selected` (index of the EBIC
#'   minimum) and the endpoints.
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
  if (!(is.finite(lambda_min) && is.finite(lambda_max) && lambda_min < lambda_max)) {
    stop("Need finite lambda_min < lambda_max.")
  }

  n_global <- min(6L, n_solves)
  n_possible <- p * (p - 1L) / 2L

  eval_one <- function(penalty) {
    if (verbose) {
      message("  penalty = ", signif(penalty, 4))
    }
    fit <- fit_precision(penalty)
    omega <- as.matrix(fit$omega)
    loglik <- fit$loglik %||% gaussian_precision_loglik(omega, S, n)
    n_edges <- count_undirected_edges(omega, eps)
    list(
      lambda = penalty,
      omega = omega,
      loglik = loglik,
      n_edges = n_edges,
      density = n_edges / n_possible,
      ebic = ebic_gaussian_ggm(loglik, n_edges, n, p, gamma),
      priority = NA_real_
    )
  }

  anchors <- exp(seq(log(lambda_max), log(lambda_min), length.out = n_global))
  path <- lapply(anchors, eval_one)

  n_basin <- min(n_solves - length(path), max(0L, round(0.4 * n_solves)))
  n_adjacent <- min(n_solves - length(path) - n_basin, max(0L, round(0.2 * n_solves)))

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
    i_min <- which.min(ebic)
    eligible <- seq_along(scores)
    if (identical(phase, "basin")) {
      eligible <- unique(pmax(1L, pmin(length(scores), (i_min - 1L):i_min)))
    } else if (identical(phase, "adjacent")) {
      eligible <- unique(pmax(1L, pmin(length(scores), (i_min - 2L):(i_min + 1L))))
    }
    pick <- eligible[[which.max(scores[eligible])]]
    lambda_new <- sqrt(path[[pick]]$lambda * path[[pick + 1L]]$lambda)
    if (any(abs(vapply(path, `[[`, numeric(1), "lambda") - lambda_new) <
      .Machine$double.eps^0.5)) {
      next
    }
    added <- eval_one(lambda_new)
    added$priority <- scores[[pick]]
    path[[length(path) + 1L]] <- added
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
    selected = which.min(ebic),
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    n = n,
    p = p,
    gamma = gamma
  )
}

#' Graphical-lasso solver for [ebic_adaptive_path()] via `huge`.
#'
#' The diagonal is left unpenalised (huge default).
#'
#' @param x Cells-by-genes numeric matrix.
#' @return A function of one `lambda`.
make_huge_solver <- function(x) {
  function(penalty) {
    fit <- huge::huge(
      x,
      method = "glasso",
      lambda = penalty,
      verbose = FALSE
    )
    list(omega = as.matrix(fit$icov[[1L]]))
  }
}

#' Hellinger distance between two non-degenerate Gaussians.
gaussian_hellinger <- function(mu1, sigma1, mu2, sigma2) {
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
  g <- length(mu1)
  s <- array(0, dim = c(g, g, 2L))
  s[, , 1L] <- sigma1
  s[, , 2L] <- sigma2
  ov <- MixSim::overlap(
    Pi = c(n1, n2) / (n1 + n2),
    Mu = rbind(mu1, mu2),
    S = s
  )
  as.numeric(ov$OmegaMap[1L, 2L] + ov$OmegaMap[2L, 1L])
}
