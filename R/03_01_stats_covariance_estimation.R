#' Regularised within-type covariances and their cross-validated choice.
#'
#' Every estimate keeps the sample variances of the raw per-cell UMI
#' counts and regularises the correlation only. With `R_j` the sample
#' correlation of type `j`, `Rbar` the pooled within-type correlation of
#' the time point, and `D_j` the sample variances:
#'
#' - diagonal: `I`;
#' - shrinkage: `(1 - lambda) R_j + lambda I`;
#' - two-parameter: `(1 - lambda) {(1 - alpha) R_j + alpha Rbar} + lambda I`;
#' - factor analysis: `B_j B_j' + Psi_j`, or `B B' + Psi_j` with loadings
#'   shared by the types of a time point.
#'
#' The covariance is `D_j^{1/2} (.) D_j^{1/2}`. Tuning values are chosen
#' by the held-out Gaussian log-score of cells, which is minimised in
#' expectation at the true mean and covariance whatever the law of the
#' cells (Stein's loss).
#'
#' No estimate is refitted per tuning value. For the two shrinkage forms
#' the eigenvectors do not depend on `lambda`, so one decomposition per
#' training set scores a whole grid. For factor analysis the score uses
#' the Woodbury identity.
#'
#' The graphical-lasso, SILGGM, PLNnetwork and denoised arms are in
#' `scripts/archives/sparse_precision_arms/` and
#' `scripts/archives/denoised_covariance/`.
#'
#' Requires `R/02_00_utils.R` (`%||%`).

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
#' @param genes Optional character vector of gene symbols. When supplied,
#'   these genes are kept, in this order, and `selected_<time>.csv` is not
#'   read. Default `NULL`: every gene of `selected_<time>.csv`.
#' @return A list with `y` (cells by genes, dense, raw counts of the
#'   panel genes, not rounded), `cell_type`, `batch`, `library_size`
#'   (full library, one value per cell), `cells`, `genes`.
read_panel_slot <- function(slots, selected_dir, time_level, genes = NULL) {
  slot <- slots$slots[[time_level]]
  if (is.null(slot)) {
    stop("No slot for time ", time_level, " in scored_slots.rds.")
  }
  if (methods::is(slot$counts, "Matrix")) {
    # Methods for `[` and `as.matrix()` on the sparse counts.
    requireNamespace("Matrix", quietly = TRUE)
  }
  if (is.null(genes)) {
    selected_path <- file.path(
      selected_dir,
      paste0("selected_", time_level, ".csv")
    )
    if (!file.exists(selected_path)) {
      stop(
        "Missing ",
        selected_path,
        ". Run scripts/02_02_stats_cross_cell_consistency.R first."
      )
    }
    genes <- utils::read.csv(
      selected_path,
      stringsAsFactors = FALSE,
      fileEncoding = "UTF-8"
    )$gene
  }
  genes <- as.character(genes)
  if (length(genes) == 0L || anyDuplicated(genes) > 0L) {
    stop("`genes` must be a non-empty vector of distinct symbols.")
  }
  missing <- setdiff(genes, rownames(slot$counts))
  if (length(missing) > 0L) {
    stop(
      length(missing),
      " panel genes of ",
      time_level,
      " are absent from the slot counts (first: ",
      missing[[1L]],
      ")."
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

# ==========================================================================
# Folds ----
# ==========================================================================

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

#' Folds of one type x time unit: leave-one-out when cells are few.
#'
#' With fewer than `loo_below` cells every cell is its own fold; a
#' training set then loses one cell instead of a fifth of them.
#'
#' @param n Number of cells.
#' @param seed Seed of the K-fold assignment.
#' @param K Number of folds above `loo_below`. Default `5`.
#' @param loo_below Cell count under which leave-one-out is used.
#' @return Integer fold labels of length `n`.
unit_folds <- function(n, seed = NULL, K = 5L, loo_below = 100L) {
  if (n < loo_below) seq_len(n) else make_cv_folds(n, K, seed)
}

# ==========================================================================
# Standardisation of a training set ----
# ==========================================================================

#' Standardise training cells, and test cells with the training moments.
#'
#' A gene that is constant in the training cells has no correlation with
#' any other gene. It is taken out of the correlation and scored as an
#' independent Gaussian with a floored variance. That part of the score
#' does not depend on any tuning value; it is returned apart (`const`)
#' so that it does not inflate the fold-to-fold spread used by the
#' one-standard-error rule.
#'
#' Every gene's variance is floored at its Poisson level plus one count
#' over the training cells, `mean + 1/n`: a count cannot vary less than a
#' Poisson count, and a variance below that level means too few cells
#' (or the fractional counts of ambient correction, e.g. a training mean
#' of 4e-7 at 96 h). The previous floor, `1e-4` times the median variance
#' for constant genes only, let one such gene cost 6e11 in a held-out
#' score.
#'
#' @param X_train Training cells by genes.
#' @param X_test Optional test cells by genes.
#' @return A list with `mu`, `sd`, `varying` (logical, per gene), `n`,
#'   `Z` (training cells by varying genes, scaled so that
#'   `crossprod(Z)` is the sample correlation, diagonal below 1 for a floored
#'   gene), and, with `X_test`: `Zt`
#'   (standardised test cells, varying genes), `base` (the part of the
#'   log-score that is the same for every correlation estimate: one
#'   number) and `const` (score of the constant genes, one value per test
#'   cell).
split_standardise <- function(X_train, X_test = NULL) {
  X_train <- as.matrix(X_train)
  n <- nrow(X_train)
  if (n < 3L) {
    stop("Need at least three training cells.")
  }
  mu <- colMeans(X_train)
  centred <- sweep(X_train, 2L, mu, check.margin = FALSE)
  v <- colSums(centred^2) / (n - 1)
  varying <- is.finite(v) & v > 0
  sd <- sqrt(pmax(ifelse(varying, v, 0), mu + 1 / n))
  names(sd) <- colnames(X_train)
  Z <- centred[, varying, drop = FALSE]
  Z <- sweep(Z, 2L, sd[varying] * sqrt(n - 1), "/", check.margin = FALSE)
  out <- list(mu = mu, sd = sd, varying = varying, n = n, Z = Z)
  if (!is.null(X_test)) {
    X_test <- as.matrix(X_test)
    zt <- sweep(X_test, 2L, mu, check.margin = FALSE)
    zt <- sweep(zt, 2L, sd, "/", check.margin = FALSE)
    out$Zt <- zt[, varying, drop = FALSE]
    out$base <- -0.5 *
      (sum(varying) * log(2 * pi) + 2 * sum(log(sd[varying])))
    out$const <- -0.5 *
      (sum(!varying) *
        log(2 * pi) +
        2 * sum(log(sd[!varying])) +
        rowSums(zt[, !varying, drop = FALSE]^2))
  }
  out
}

#' Non-zero spectrum and right singular vectors of a matrix.
#'
#' Works on the smaller of `Z Z'` and `Z' Z`.
#'
#' @param Z Matrix, `n` by `p`.
#' @param tol Relative threshold under which an eigenvalue is zero.
#' @return A list with `values` (eigenvalues of `Z' Z`, decreasing,
#'   zeros dropped) and `vectors` (`p` by `length(values)`).
right_spectrum <- function(Z, tol = 1e-10) {
  n <- nrow(Z)
  p <- ncol(Z)
  if (n < p) {
    eg <- eigen(tcrossprod(Z), symmetric = TRUE)
    keep <- eg$values > tol * max(eg$values, 1)
    values <- eg$values[keep]
    vectors <- crossprod(Z, eg$vectors[, keep, drop = FALSE])
    vectors <- sweep(vectors, 2L, sqrt(values), "/", check.margin = FALSE)
  } else {
    eg <- eigen(crossprod(Z), symmetric = TRUE)
    keep <- eg$values > tol * max(eg$values, 1)
    values <- eg$values[keep]
    vectors <- eg$vectors[, keep, drop = FALSE]
  }
  list(values = values, vectors = vectors)
}

# ==========================================================================
# Scores along a shrinkage path ----
# ==========================================================================

#' Log-score core of `(1 - lambda) M + lambda I` for a grid of `lambda`.
#'
#' With `M = V diag(values) V'` on its range and zero elsewhere, the
#' eigenvalues of the shrunk matrix are `(1 - lambda) values + lambda` on
#' the range and `lambda` on the null space. For a standardised test
#' cell `z`, the quadratic form is
#' `sum_i (v_i' z)^2 / {(1 - lambda) values_i + lambda}
#'   + (|z|^2 - sum_i (v_i' z)^2) / lambda`.
#'
#' @param Zt Standardised test cells by genes.
#' @param values,vectors Non-zero spectrum of `M`.
#' @param grid Values of `lambda` in `[0, 1]`.
#' @return Matrix, test cells by grid values, of
#'   `-0.5 * (log det + quadratic form)`; `-Inf` at `lambda = 0` when
#'   `M` is singular.
spectral_path_core <- function(Zt, values, vectors, grid) {
  p <- ncol(Zt)
  r <- length(values)
  w2 <- (Zt %*% vectors)^2
  rest <- pmax(rowSums(Zt^2) - rowSums(w2), 0)
  out <- matrix(NA_real_, nrow(Zt), length(grid))
  for (i in seq_along(grid)) {
    lam <- grid[[i]]
    e <- (1 - lam) * values + lam
    if (lam <= 0 && r < p) {
      out[, i] <- -Inf
      next
    }
    log_det <- sum(log(e)) + if (r < p) (p - r) * log(lam) else 0
    quad <- drop(w2 %*% (1 / e)) + if (r < p) rest / lam else 0
    out[, i] <- -0.5 * (log_det + quad)
  }
  out
}

#' Pooled within-type correlation.
#'
#' The correlation of `W = sum_j (n_j - 1) V_j / (N - J)`, the tied
#' covariance of linear discriminant analysis, over the cell types of one
#' time point and one gene panel. A gene with no pooled variance gets a
#' unit diagonal and no correlation.
#'
#' @param scatter Sum over types of the centred cross-products, genes by
#'   genes.
#' @return A correlation matrix.
pooled_correlation <- function(scatter) {
  d <- diag(scatter)
  ok <- is.finite(d) & d > 0
  s <- ifelse(ok, sqrt(d), 1)
  R <- scatter / tcrossprod(s)
  R[!ok, ] <- 0
  R[, !ok] <- 0
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R
}

#' Centred cross-product of a set of cells.
#'
#' @param X Cells by genes.
#' @return Genes by genes matrix `sum_c (x_c - xbar)(x_c - xbar)'`.
centred_scatter <- function(X) {
  crossprod(sweep(as.matrix(X), 2L, colMeans(X), check.margin = FALSE))
}

# ==========================================================================
# Gaussian factor analysis on the correlation scale ----
# ==========================================================================

#' Woodbury pieces of `B B' + Psi`.
#'
#' @param B Loadings, genes by factors.
#' @param psi Positive uniquenesses, one per gene.
#' @return A list with `Bp = Psi^-1 B`, `Minv = (I + B' Psi^-1 B)^-1` and
#'   `log_det`, the log-determinant of `B B' + Psi`.
factor_woodbury <- function(B, psi) {
  q <- ncol(B)
  Bp <- B / psi
  if (q == 0L) {
    return(list(Bp = Bp, Minv = matrix(0, 0L, 0L), log_det = sum(log(psi))))
  }
  M <- diag(q) + crossprod(B, Bp)
  R <- chol((M + t(M)) / 2)
  list(
    Bp = Bp,
    Minv = chol2inv(R),
    log_det = sum(log(psi)) + 2 * sum(log(diag(R)))
  )
}

#' Log-score core of `B B' + Psi` for standardised test cells.
#'
#' @param Zt Standardised test cells by genes.
#' @param B,psi Factor model on the same genes.
#' @return One value per test cell, `-0.5 * (log det + quadratic form)`.
factor_core <- function(Zt, B, psi) {
  w <- factor_woodbury(B, psi)
  quad <- drop(Zt^2 %*% (1 / psi))
  if (ncol(B) > 0L) {
    u <- Zt %*% w$Bp
    quad <- quad - rowSums((u %*% w$Minv) * u)
  }
  -0.5 * (w$log_det + quad)
}

#' Maximum-likelihood factor analysis of a correlation, by EM.
#'
#' Model `R = B B' + Psi` with `q` factors. The EM steps are those of
#' Rubin and Thayer, written on the cells so that the genes by genes
#' correlation is never formed: with `beta = (I + B' Psi^-1 B)^-1 B'
#' Psi^-1`, `B <- R beta' (I - beta B + beta R beta')^-1` and
#' `Psi <- diag(R - B beta R)`. The start is the probabilistic PCA
#' solution. Uniquenesses are floored: with fewer cells than genes the
#' likelihood is unbounded when one of them tends to zero.
#'
#' @param Z Cells by genes, scaled so that `crossprod(Z)` is the sample
#'   correlation.
#' @param q Number of factors; `0` returns the identity.
#' @param spectrum Optional output of [right_spectrum()] on `Z`.
#' @param psi_floor Smallest uniqueness. Default `0.01`.
#' @param maxit,tol Stop after `maxit` steps or when no uniqueness moves
#'   by more than `tol`.
#' @return A list with `B`, `psi`, `iterations`, `converged` and
#'   `n_floored`.
fit_factor_correlation <- function(
  Z,
  q,
  spectrum = NULL,
  psi_floor = 0.01,
  maxit = 200L,
  tol = 1e-5
) {
  p <- ncol(Z)
  q <- as.integer(q)
  if (q == 0L) {
    return(list(
      B = matrix(0, p, 0L),
      psi = rep(1, p),
      iterations = 0L,
      converged = TRUE,
      n_floored = 0L
    ))
  }
  spectrum <- spectrum %||% right_spectrum(Z)
  if (q > min(length(spectrum$values), p - 1L)) {
    stop("`q` exceeds the rank of the training cells.")
  }
  ell <- spectrum$values[seq_len(q)]
  noise <- max((p - sum(ell)) / (p - q), psi_floor)
  B <- sweep(
    spectrum$vectors[, seq_len(q), drop = FALSE],
    2L,
    sqrt(pmax(ell - noise, 1e-8)),
    "*",
    check.margin = FALSE
  )
  psi <- pmax(1 - rowSums(B^2), psi_floor)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    w <- factor_woodbury(B, psi)
    beta_t <- w$Bp %*% w$Minv
    r_beta <- crossprod(Z, Z %*% beta_t)
    e_zz <- w$Minv + crossprod(beta_t, r_beta)
    B <- r_beta %*% solve((e_zz + t(e_zz)) / 2)
    psi_new <- pmax(1 - rowSums(B * r_beta), psi_floor)
    moved <- max(abs(psi_new - psi))
    psi <- psi_new
    if (moved < tol) {
      converged <- TRUE
      break
    }
  }
  list(
    B = B,
    psi = psi,
    iterations = it,
    converged = converged,
    n_floored = sum(psi <= psi_floor)
  )
}

#' Factor analysis with loadings shared by the types of a time point.
#'
#' Model `R_j = B B' + Psi_j`: the same axes for every type, and one
#' uniqueness per gene and type. Conditional maximisation: the row of
#' `B` for gene `g` solves a `q` by `q` system weighted by
#' `n_j / psi_jg`; `Psi_j` follows in closed form. A gene constant in a
#' type has weight zero for that type.
#'
#' @param Z_list One matrix per type, cells by all genes, scaled as in
#'   [fit_factor_correlation()], with zero columns for the genes that are
#'   constant in that type.
#' @param varying_list One logical vector per type, `FALSE` for those
#'   genes.
#' @param q Number of factors.
#' @param psi_floor,maxit,tol As in [fit_factor_correlation()].
#' @return A list with `B` (genes by factors), `psi` (genes by types),
#'   `iterations` and `converged`.
fit_factor_shared <- function(
  Z_list,
  varying_list,
  q,
  psi_floor = 0.01,
  maxit = 100L,
  tol = 1e-4
) {
  n_types <- length(Z_list)
  p <- ncol(Z_list[[1L]])
  n_cells <- vapply(Z_list, nrow, integer(1))
  q <- as.integer(q)
  psi <- matrix(1, p, n_types, dimnames = list(NULL, names(Z_list)))
  if (q == 0L) {
    return(list(
      B = matrix(0, p, 0L),
      psi = psi,
      iterations = 0L,
      converged = TRUE
    ))
  }
  # Start: principal axes of the cell-weighted mean correlation.
  weight <- sqrt((n_cells - 1) / sum(n_cells - 1))
  stacked <- do.call(
    rbind,
    lapply(seq_len(n_types), function(j) Z_list[[j]] * weight[[j]])
  )
  spectrum <- right_spectrum(stacked)
  rm(stacked)
  if (q > min(length(spectrum$values), p - 1L)) {
    stop("`q` exceeds the rank of the training cells.")
  }
  ell <- spectrum$values[seq_len(q)]
  noise <- max((p - sum(ell)) / (p - q), psi_floor)
  B <- sweep(
    spectrum$vectors[, seq_len(q), drop = FALSE],
    2L,
    sqrt(pmax(ell - noise, 1e-8)),
    "*",
    check.margin = FALSE
  )
  for (j in seq_len(n_types)) {
    psi[, j] <- pmax(1 - rowSums(B^2), psi_floor)
  }
  converged <- FALSE
  for (it in seq_len(maxit)) {
    steps <- lapply(seq_len(n_types), function(j) {
      Bj <- B * varying_list[[j]]
      w <- factor_woodbury(Bj, psi[, j])
      beta_t <- w$Bp %*% w$Minv
      r_beta <- crossprod(Z_list[[j]], Z_list[[j]] %*% beta_t)
      list(r_beta = r_beta, e_zz = w$Minv + crossprod(beta_t, r_beta))
    })
    weights <- vapply(
      seq_len(n_types),
      function(j) n_cells[[j]] * varying_list[[j]] / psi[, j],
      numeric(p)
    )
    rhs <- Reduce(
      `+`,
      lapply(seq_len(n_types), function(j) weights[, j] * steps[[j]]$r_beta)
    )
    e_flat <- matrix(
      vapply(steps, function(s) as.numeric(s$e_zz), numeric(q * q)),
      nrow = q * q
    )
    lhs <- weights %*% t(e_flat) # genes by q^2
    for (g in seq_len(p)) {
      a <- matrix(lhs[g, ], q, q)
      if (sum(weights[g, ]) > 0) {
        B[g, ] <- solve((a + t(a)) / 2, rhs[g, ])
      } else {
        B[g, ] <- 0
      }
    }
    psi_new <- psi
    for (j in seq_len(n_types)) {
      fitted <- 1 -
        2 * rowSums(B * steps[[j]]$r_beta) +
        rowSums((B %*% steps[[j]]$e_zz) * B)
      psi_new[, j] <- ifelse(varying_list[[j]], pmax(fitted, psi_floor), 1)
    }
    moved <- max(abs(psi_new - psi))
    psi <- psi_new
    if (moved < tol) {
      converged <- TRUE
      break
    }
  }
  list(B = B, psi = psi, iterations = it, converged = converged)
}

# ==========================================================================
# Held-out scores of every arm for one training and test set ----
# ==========================================================================

#' Largest number of factors a correlation of `p` genes can identify.
#'
#' The model has `p + q p - q (q - 1) / 2` free parameters, which must
#' not exceed the `p (p + 1) / 2` entries of the matrix.
#'
#' @param p Genes.
#' @return The largest admissible `q`.
factor_rank_bound <- function(p) {
  floor(p + 0.5 * (1 - sqrt(1 + 8 * p)))
}

#' Usable factor ranks for a training set.
#'
#' @param ranks Candidate ranks.
#' @param n Training cells.
#' @param rank Rank of the training correlation.
#' @param p Genes.
#' @return The ranks no larger than `min(n - 2, rank)` and than
#'   [factor_rank_bound()].
usable_ranks <- function(ranks, n, rank, p = Inf) {
  ranks[ranks <= min(n - 2L, rank, factor_rank_bound(p))]
}

#' Mean held-out log-score of the per-type arms, at every tuning value.
#'
#' One spectral decomposition of the training correlation serves the
#' shrinkage path and the start of every factor fit. The two-parameter
#' arm needs one eigen-decomposition per `alpha`.
#'
#' @param X_train,X_test Cells by genes.
#' @param lambda_grid Shrinkage levels.
#' @param ranks Factor ranks; `NULL` skips factor analysis.
#' @param alpha_grid Pooling weights; `NULL` skips the two-parameter arm.
#' @param other_scatter Centred cross-product of the training cells of
#'   the other types of the time point (needed with `alpha_grid`).
#' @return A list with `shrinkage` (one mean score per `lambda`),
#'   `factor` (one per rank, `NA` when the rank is not usable),
#'   `two_param` (`alpha` by `lambda`), `const` (mean score of the genes
#'   constant in training; add it for the full log-score) and `n_test`.
score_unit_arms <- function(
  X_train,
  X_test,
  lambda_grid,
  ranks = NULL,
  alpha_grid = NULL,
  other_scatter = NULL
) {
  sp <- split_standardise(X_train, X_test)
  spectrum <- right_spectrum(sp$Z)
  out <- list(const = mean(sp$const), n_test = nrow(sp$Zt))
  out$shrinkage <- sp$base +
    colMeans(spectral_path_core(
      sp$Zt,
      spectrum$values,
      spectrum$vectors,
      lambda_grid
    ))
  if (!is.null(ranks)) {
    ok <- usable_ranks(ranks, sp$n, length(spectrum$values), ncol(sp$Z))
    out$factor <- vapply(
      ranks,
      function(q) {
        if (!q %in% ok) {
          return(NA_real_)
        }
        fit <- fit_factor_correlation(sp$Z, q, spectrum)
        sp$base + mean(factor_core(sp$Zt, fit$B, fit$psi))
      },
      numeric(1)
    )
  }
  if (!is.null(alpha_grid)) {
    if (is.null(other_scatter)) {
      stop("The two-parameter arm needs `other_scatter`.")
    }
    pooled <- pooled_correlation(other_scatter + centred_scatter(X_train))
    pooled <- pooled[sp$varying, sp$varying, drop = FALSE]
    own <- crossprod(sp$Z)
    out$two_param <- t(vapply(
      alpha_grid,
      function(alpha) {
        eg <- eigen((1 - alpha) * own + alpha * pooled, symmetric = TRUE)
        keep <- eg$values > 1e-10 * max(eg$values)
        sp$base +
          colMeans(spectral_path_core(
            sp$Zt,
            eg$values[keep],
            eg$vectors[, keep, drop = FALSE],
            lambda_grid
          ))
      },
      numeric(length(lambda_grid))
    ))
    dimnames(out$two_param) <- list(NULL, NULL)
  }
  out
}

#' Mean held-out log-score of the shared-loadings factor model.
#'
#' @param train_list,test_list One cells by genes matrix per type.
#' @param ranks Factor ranks.
#' @return A list with `scores` (types by ranks, mean scores without the
#'   part of the constant genes), `const` (that part, one value per
#'   type) and `n_test`.
score_shared_factor <- function(train_list, test_list, ranks) {
  splits <- Map(split_standardise, train_list, test_list)
  p <- ncol(train_list[[1L]])
  Z_list <- lapply(splits, function(sp) {
    Z <- matrix(0, sp$n, p)
    Z[, sp$varying] <- sp$Z
    Z
  })
  varying_list <- lapply(splits, `[[`, "varying")
  ok <- usable_ranks(
    ranks,
    sum(vapply(splits, `[[`, numeric(1), "n")) - length(splits),
    p,
    min(vapply(varying_list, sum, numeric(1)))
  )
  out <- vapply(
    ranks,
    function(q) {
      if (!q %in% ok) {
        return(rep(NA_real_, length(splits)))
      }
      fit <- fit_factor_shared(Z_list, varying_list, q)
      vapply(
        seq_along(splits),
        function(j) {
          sp <- splits[[j]]
          sp$base +
            mean(factor_core(
              sp$Zt,
              fit$B[sp$varying, , drop = FALSE],
              fit$psi[sp$varying, j]
            ))
        },
        numeric(1)
      )
    },
    numeric(length(splits))
  )
  list(
    scores = matrix(
      out,
      nrow = length(splits),
      dimnames = list(names(train_list), as.character(ranks))
    ),
    const = vapply(splits, function(sp) mean(sp$const), numeric(1)),
    n_test = vapply(splits, function(sp) nrow(sp$Zt), integer(1))
  )
}

#' Cross-validated scores of the per-type arms for one unit.
#'
#' @param X Cells by genes of one type at one time.
#' @param folds Fold label of each cell.
#' @param lambda_grid,ranks,alpha_grid,other_scatter See
#'   [score_unit_arms()]. `other_scatter` holds the other types only; the
#'   training cells of this type are added in each fold.
#' @return A list with one element per arm: a matrix, folds by tuning
#'   values, of mean held-out scores; and `const`, one value per fold.
cv_unit_arms <- function(
  X,
  folds,
  lambda_grid,
  ranks = NULL,
  alpha_grid = NULL,
  other_scatter = NULL
) {
  parts <- lapply(sort(unique(folds)), function(k) {
    score_unit_arms(
      X[folds != k, , drop = FALSE],
      X[folds == k, , drop = FALSE],
      lambda_grid,
      ranks,
      alpha_grid,
      other_scatter
    )
  })
  out <- list(
    shrinkage = do.call(rbind, lapply(parts, `[[`, "shrinkage")),
    const = vapply(parts, `[[`, numeric(1), "const")
  )
  if (!is.null(ranks)) {
    out$factor <- do.call(rbind, lapply(parts, `[[`, "factor"))
  }
  if (!is.null(alpha_grid)) {
    out$two_param <- do.call(
      rbind,
      lapply(parts, function(p) as.numeric(t(p$two_param)))
    )
  }
  out
}

#' Mean, standard error and the two usual choices along one path.
#'
#' `min` maximises the mean held-out score. `1se` is the most regularised
#' value whose mean is within one standard error of that maximum. A
#' value that could not be fitted in every fold is not eligible.
#'
#' @param fold_scores Folds by tuning values.
#' @param tuning The tuning values.
#' @param regularised `"larger"` when a larger value regularises more
#'   (shrinkage level), `"smaller"` when a smaller one does (factor
#'   rank).
#' @return A list with `table` (`tuning`, `mean_score`, `se_score`),
#'   `min` and `one_se`.
summarise_path <- function(
  fold_scores,
  tuning,
  regularised = c("larger", "smaller")
) {
  regularised <- match.arg(regularised)
  mean_score <- colMeans(fold_scores)
  se_score <- apply(fold_scores, 2L, stats::sd) / sqrt(nrow(fold_scores))
  ok <- which(is.finite(mean_score))
  if (length(ok) == 0L) {
    stop("No tuning value was scored in every fold.")
  }
  best <- ok[[which.max(mean_score[ok])]]
  within <- ok[mean_score[ok] >= mean_score[[best]] - se_score[[best]]]
  one_se <- if (regularised == "larger") {
    within[[which.max(tuning[within])]]
  } else {
    within[[which.min(tuning[within])]]
  }
  list(
    table = data.frame(
      tuning = tuning,
      mean_score = mean_score,
      se_score = se_score
    ),
    min = tuning[[best]],
    one_se = tuning[[one_se]]
  )
}

# ==========================================================================
# Reference score through a Cholesky factor (checks only) ----
# ==========================================================================

#' Held-out Gaussian log-score through a Cholesky factor.
#'
#' The direct computation, one factorisation per covariance. The path
#' functions above are checked against it.
#'
#' @param mu Mean vector, length `p`.
#' @param Sigma Covariance, `p` by `p`.
#' @param X_test Test cells in rows, genes in columns.
#' @param average If `TRUE` (default), the mean log-score per test cell.
#' @return A scalar (or one value per test cell); `-Inf` when `Sigma` is
#'   not numerically positive definite.
heldout_gaussian_logscore <- function(mu, Sigma, X_test, average = TRUE) {
  X_test <- as.matrix(X_test)
  Sigma <- (as.matrix(Sigma) + t(as.matrix(Sigma))) / 2
  p <- ncol(X_test)
  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(R)) {
    return(if (average) -Inf else rep(-Inf, nrow(X_test)))
  }
  z <- backsolve(R, t(X_test) - as.numeric(mu), transpose = TRUE)
  per_cell <- -0.5 *
    (p * log(2 * pi) + 2 * sum(log(diag(R))) + colSums(z * z))
  if (average) mean(per_cell) else per_cell
}

# ==========================================================================
# Plug-in covariances for the deconvolution ----
# ==========================================================================

#' Sample moments of every type of one time point and gene panel.
#'
#' @param y Cells by genes, raw counts.
#' @param cell_type Type label of each cell.
#' @param types Types to keep, in the order of the output.
#' @param library_size Optional library size of each cell. When given,
#'   the covariance (not the mean) is computed on cells divided by their
#'   size factor `L_c / mean(L)` within the type: the size term of the
#'   covariance is removed and the Poisson part kept (the "all technical"
#'   end of the kappa bracket of chapter 3, without denoising).
#' @return A list with `mu` (genes by types), `sigma_sample` (genes by
#'   genes by types, unbiased sample covariance, not regularised),
#'   `n_cells`, `r_pool` (pooled within-type correlation) and, with
#'   `library_size`, `mu_cov` (means of the size-normalised cells).
sample_moments <- function(y, cell_type, types, library_size = NULL) {
  genes <- colnames(y)
  mu <- matrix(NA_real_, length(genes), length(types), dimnames = list(genes, types))
  sigma <- array(
    NA_real_,
    c(length(genes), length(genes), length(types)),
    dimnames = list(genes, genes, types)
  )
  n_cells <- stats::setNames(integer(length(types)), types)
  scatter <- matrix(0, length(genes), length(genes))
  mu_cov <- if (!is.null(library_size)) mu
  for (tp in types) {
    rows <- which(cell_type %in% tp)
    X <- y[rows, , drop = FALSE]
    mu[, tp] <- colMeans(X)
    if (!is.null(library_size)) {
      X <- X / (library_size[rows] / mean(library_size[rows]))
      mu_cov[, tp] <- colMeans(X)
    }
    s <- centred_scatter(X)
    sigma[,, tp] <- s / (nrow(X) - 1)
    n_cells[[tp]] <- nrow(X)
    scatter <- scatter + s
  }
  r_pool <- pooled_correlation(scatter)
  dimnames(r_pool) <- list(genes, genes)
  list(mu = mu, sigma_sample = sigma, n_cells = n_cells, r_pool = r_pool, mu_cov = mu_cov)
}

#' Regularised covariances of every type, for one arm and tuning value.
#'
#' @param moments List with `sigma_sample` (genes by genes by types),
#'   `types`, and `r_pool` for the two-parameter arm.
#' @param arm `"shrinkage"`, `"two_param"` or `"factor"`.
#' @param lambda Shrinkage level: one value, or one per type (named).
#' @param alpha Pooling weight of the two-parameter arm: one value, or
#'   one per type (named).
#' @param factors For `arm = "factor"`: one list per type with `B` and
#'   `psi` on the correlation scale ([fit_factor_correlation()] or the
#'   columns of [fit_factor_shared()]).
#' @param var_floor_rel Shrinkage arms: a gene with no variance in a type
#'   gets this fraction of the median variance of the type. The factor arm
#'   uses the scale of its fit instead, the Poisson floor of
#'   [split_standardise()] (`mean + 1/n`), so that `B B' + Psi` maps back
#'   to the sample variance.
#' @return An array like `moments$sigma_sample`.
build_sigma <- function(
  moments,
  arm = c("shrinkage", "two_param", "factor"),
  lambda = NULL,
  alpha = 0,
  factors = NULL,
  var_floor_rel = 1e-4
) {
  arm <- match.arg(arm)
  types <- dimnames(moments$sigma_sample)[[3L]]
  p <- dim(moments$sigma_sample)[[1L]]
  if (arm != "factor") {
    if (length(lambda) == 1L) {
      lambda <- stats::setNames(rep(lambda, length(types)), types)
    }
    if (!all(types %in% names(lambda)) || any(lambda < 0 | lambda > 1)) {
      stop("`lambda` needs one value in [0, 1] per type.")
    }
    if (length(alpha) == 1L) {
      alpha <- stats::setNames(rep(alpha, length(types)), types)
    }
  }
  out <- moments$sigma_sample
  for (tp in types) {
    v <- diag(out[,, tp])
    sd <- if (arm == "factor") {
      sqrt(pmax(v, (moments$mu_cov %||% moments$mu)[, tp] + 1 / moments$n_cells[[tp]]))
    } else {
      sqrt(pmax(v, var_floor_rel * stats::median(v[v > 0])))
    }
    R <- out[,, tp] / tcrossprod(sd)
    if (arm == "shrinkage") {
      R <- (1 - lambda[[tp]]) * R
    } else if (arm == "two_param") {
      R <- (1 - lambda[[tp]]) *
        ((1 - alpha[[tp]]) * R + alpha[[tp]] * moments$r_pool)
    } else {
      f <- factors[[tp]]
      R <- tcrossprod(f$B) + diag(f$psi, p)
    }
    if (arm != "factor") {
      diag(R) <- 1
    }
    out[,, tp] <- R * tcrossprod(sd)
  }
  out
}

#' Free covariance parameters of the arms that have a count.
#'
#' Factor analysis with `q` factors on `p` genes has
#' `p + q p - q (q - 1) / 2` free covariance parameters (loadings up to
#' rotation, plus uniquenesses); `q = 0` is the diagonal model. With
#' shared loadings over `n_types` types the loadings are counted once.
#' The shrinkage arms have no such count.
#'
#' @param p Genes.
#' @param q Factors.
#' @param n_types Types sharing the loadings; `1` for a per-type model.
#' @return The number of free covariance parameters.
factor_parameter_count <- function(p, q, n_types = 1L) {
  q * p - q * (q - 1) / 2 + n_types * p
}

#' AIC and BIC from a log-likelihood and a parameter count.
#'
#' @param loglik Maximised log-likelihood of the cells.
#' @param k Free parameters, means included.
#' @param n Cells.
#' @return A named vector with `aic` and `bic`.
information_criteria <- function(loglik, k, n) {
  c(aic = -2 * loglik + 2 * k, bic = -2 * loglik + k * log(n))
}

#' Factor models of every type, on raw or size-normalised cells.
#'
#' The fits of the `fit` stage of `scripts/03_04_cv_covariance.R`, padded
#' to all genes (a gene constant in a type gets no loading and
#' uniqueness 1), with the largest cosine between a loading and the
#' standardised mean of the type (`align`; near 1 for a cell-size
#' factor).
#'
#' @param y Cells by genes, raw counts.
#' @param cell_type Type of each cell.
#' @param types Types to fit.
#' @param ranks Ranks; a rank above the identifiable bound of a type is
#'   skipped.
#' @param library_size Optional library sizes; when given, cells are
#'   divided by `L_c / mean(L)` within their type first.
#' @return One list per type, named by rank, of `B`, `psi`, `align`.
fit_factor_types <- function(y, cell_type, types, ranks, library_size = NULL) {
  genes <- colnames(y)
  lapply(stats::setNames(types, types), function(tp) {
    rows <- which(cell_type %in% tp)
    X <- y[rows, , drop = FALSE]
    if (!is.null(library_size)) {
      X <- X / (library_size[rows] / mean(library_size[rows]))
    }
    sp <- split_standardise(X)
    spectrum <- right_spectrum(sp$Z)
    ok <- usable_ranks(ranks, sp$n, length(spectrum$values), ncol(sp$Z))
    z_mean <- (sp$mu / sp$sd)[sp$varying]
    lapply(stats::setNames(ok, ok), function(q) {
      fit <- fit_factor_correlation(sp$Z, q, spectrum)
      B <- matrix(0, length(genes), q)
      B[sp$varying, ] <- fit$B
      psi <- rep(1, length(genes))
      psi[sp$varying] <- fit$psi
      list(
        B = B,
        psi = psi,
        align = max(abs(drop(crossprod(fit$B, z_mean))) /
          (sqrt(colSums(fit$B^2)) * sqrt(sum(z_mean^2))))
      )
    })
  })
}
