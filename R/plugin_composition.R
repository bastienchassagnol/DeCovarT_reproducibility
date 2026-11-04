#' Plug-in composition of a bulk sample for the tier (iii) tail bounds.
#'
#' The bulk consistency statistic `D_gi = y_gi / a_i - sum_j p_ji mu_gj`
#' needs a composition `p_i` and a scale `a_i`. They are plugged in from
#' an iterated weighted non-negative least-squares fit in the style of
#' deMixSC (Eq. 4): weights `1 / (yhat_g^2 + (y_g - yhat_g)^2 + c)`, with
#' `c = 2`. deMixSC also separates a gene set G2 calibrated on a
#' same-aliquot benchmark; no such benchmark exists for Suppinger, so G2
#' is empty and every gene enters the fit. `a_i` is the sum of the NNLS
#' coefficients and `p_i` the coefficients divided by `a_i`; the extra
#' overdispersion between platforms is fixed at 1. MuSiC is not used,
#' because its weights need at least two donors.
#'
#' `c` has the units of a squared expression value, so it is not
#' scale-free; bulk and signature must be on the same linear scale.
#' The package `nnls` is not in `renv.lock`, so the Lawson-Hanson active-
#' set algorithm is written inline ([nnls_lawson_hanson()]).
#'
#' Every tail statement built on this plug-in holds only with these
#' parameters treated as known. Written, not run.

#' Non-negative least squares, Lawson-Hanson active-set algorithm.
#'
#' Solves `min ||A x - b||_2` subject to `x >= 0` (Lawson and Hanson
#' 1974, chapter 23). Variables move from the active set (fixed at zero)
#' to the passive set while the gradient `A'(b - A x)` has a positive
#' entry; the inner loop steps back along the segment to keep the
#' passive variables positive.
#'
#' @param a Design matrix (rows: equations, columns: variables).
#' @param b Response, one value per row of `a`.
#' @param tol Tolerance on the gradient. Default
#'   `10 eps ||a||_1 max(dim(a))`, as in MATLAB's `lsqnonneg`.
#' @param max_iter Maximum number of outer iterations. Default
#'   `3 ncol(a)`.
#' @return The solution vector, with attributes `iterations` and
#'   `converged`.
nnls_lawson_hanson <- function(a, b, tol = NULL, max_iter = NULL) {
  a <- as.matrix(a)
  b <- as.numeric(b)
  if (nrow(a) != length(b)) {
    stop("`b` must have one value per row of `a`.")
  }
  if (any(!is.finite(a)) || any(!is.finite(b))) {
    stop("`a` and `b` must be finite.")
  }
  n <- ncol(a)
  if (is.null(tol)) {
    tol <- 10 * .Machine$double.eps * norm(a, "1") * max(dim(a))
  }
  if (is.null(max_iter)) {
    max_iter <- 3L * n
  }
  ls_passive <- function(passive) {
    z <- numeric(n)
    coef <- qr.coef(qr(a[, passive, drop = FALSE]), b)
    coef[is.na(coef)] <- 0
    z[passive] <- coef
    z
  }
  x <- numeric(n)
  passive <- logical(n)
  grad <- drop(crossprod(a, b - a %*% x))
  iter <- 0L
  converged <- TRUE
  while (any(!passive) && max(grad[!passive]) > tol) {
    iter <- iter + 1L
    if (iter > max_iter) {
      converged <- FALSE
      warning("nnls_lawson_hanson() reached `max_iter`.", call. = FALSE)
      break
    }
    j <- which(!passive)[which.max(grad[!passive])]
    passive[j] <- TRUE
    z <- ls_passive(passive)
    inner <- 0L
    while (any(z[passive] <= tol)) {
      inner <- inner + 1L
      if (inner > max_iter) {
        converged <- FALSE
        break
      }
      neg <- passive & z <= tol
      ratio <- x[neg] / (x[neg] - z[neg])
      ratio[!is.finite(ratio)] <- 0
      step <- min(ratio)
      x <- x + step * (z - x)
      passive[passive & x <= tol] <- FALSE
      x[!passive] <- 0
      if (!any(passive)) {
        z <- numeric(n)
        break
      }
      z <- ls_passive(passive)
    }
    x <- z
    grad <- drop(crossprod(a, b - a %*% x))
  }
  x <- pmax(x, 0)
  attr(x, "iterations") <- iter
  attr(x, "converged") <- converged
  x
}

#' deMixSC-style iterated weighted NNLS composition of bulk samples.
#'
#' For each bulk sample: start from the unweighted NNLS fit, then
#' iterate `w_g = 1 / (yhat_g^2 + (y_g - yhat_g)^2 + c_reg)`,
#' `beta = NNLS(sqrt(w) X, sqrt(w) y)` until the largest change of
#' `beta`, relative to `sum(beta)`, is below `tol`.
#'
#' @param bulk Bulk expression: named vector over genes, or genes by
#'   samples matrix (raw counts of the panel genes).
#' @param signature Genes by types reference means per cell (raw counts
#'   per cell, so that a bulk is `a_i sum_j p_ji mu_gj`), gene names as
#'   row names.
#' @param c_reg Constant `c` of the deMixSC weights. Default `2`.
#' @param max_iter Maximum number of reweighting steps. Default `100`.
#' @param tol Convergence tolerance. Default `1e-6`.
#' @return A list with
#'   - `composition`: one row per sample and type: `sample`, `cell_type`,
#'     `coef` (`a_i p_ji`), `p_hat`;
#'   - `scale`: one row per sample: `sample`, `a_hat`, `iterations`,
#'     `converged`;
#'   - `weights`: genes by samples matrix of final weights;
#'   - `fitted`: genes by samples matrix of `X beta`.
plugin_composition_demixsc <- function(
  bulk,
  signature,
  c_reg = 2,
  max_iter = 100L,
  tol = 1e-6
) {
  signature <- as.matrix(signature)
  genes <- rownames(signature)
  if (is.null(genes) || is.null(colnames(signature))) {
    stop("`signature` needs gene names (rows) and type names (columns).")
  }
  if (any(!is.finite(signature)) || any(signature < 0)) {
    stop("`signature` must be finite and non-negative.")
  }
  bulk <- if (is.null(dim(bulk))) {
    matrix(bulk, ncol = 1L, dimnames = list(names(bulk), "sample_1"))
  } else {
    as.matrix(bulk)
  }
  if (is.null(rownames(bulk)) || !all(genes %in% rownames(bulk))) {
    stop("`bulk` needs row names (or names) covering every signature gene.")
  }
  bulk <- bulk[genes, , drop = FALSE]
  if (is.null(colnames(bulk))) {
    colnames(bulk) <- paste0("sample_", seq_len(ncol(bulk)))
  }
  if (!is.numeric(c_reg) || length(c_reg) != 1L || c_reg <= 0) {
    stop("`c_reg` must be one positive number.")
  }

  types <- colnames(signature)
  weights <- matrix(NA_real_, nrow(bulk), ncol(bulk), dimnames = dimnames(bulk))
  fitted <- weights
  comp <- list()
  scale <- list()
  for (i in seq_len(ncol(bulk))) {
    y <- bulk[, i]
    beta <- nnls_lawson_hanson(signature, y)
    converged <- FALSE
    it <- 0L
    for (it in seq_len(max_iter)) {
      yhat <- drop(signature %*% beta)
      w <- 1 / (yhat^2 + (y - yhat)^2 + c_reg)
      beta_new <- nnls_lawson_hanson(signature * sqrt(w), y * sqrt(w))
      change <- max(abs(beta_new - beta)) /
        max(sum(abs(beta)), .Machine$double.eps)
      beta <- beta_new
      if (change < tol) {
        converged <- TRUE
        break
      }
    }
    if (!converged) {
      warning(
        "deMixSC reweighting did not converge for sample ",
        colnames(bulk)[i],
        " in ",
        max_iter,
        " steps.",
        call. = FALSE
      )
    }
    yhat <- drop(signature %*% beta)
    weights[, i] <- 1 / (yhat^2 + (y - yhat)^2 + c_reg)
    fitted[, i] <- yhat
    a_hat <- sum(beta)
    comp[[i]] <- data.frame(
      sample = colnames(bulk)[i],
      cell_type = types,
      coef = as.numeric(beta),
      p_hat = if (a_hat > 0) as.numeric(beta) / a_hat else NA_real_,
      stringsAsFactors = FALSE
    )
    scale[[i]] <- data.frame(
      sample = colnames(bulk)[i],
      a_hat = a_hat,
      iterations = it,
      converged = converged,
      stringsAsFactors = FALSE
    )
  }
  list(
    composition = do.call(rbind, comp),
    scale = do.call(rbind, scale),
    weights = weights,
    fitted = fitted
  )
}
