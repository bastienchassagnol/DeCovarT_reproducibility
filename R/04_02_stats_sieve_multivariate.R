#' Bulk covariance of the DeCovarT law and its likelihood at a composition.
#'
#' Notation (DeCovarT vignettes): `G` genes, `J` cell types, `M` the
#' `G x J` matrix of mean counts per cell (`mu`), `Sigma_j` the
#' within-type covariance, `p` the cellular proportions.
#'
#' The quadratic-programme chain, the between-type DeCovarT fit and the
#' Firth fit of the second pass are in
#' `scripts/archives/second_pass_variants/`.

# Covariance and its derivatives ---------------------------------------------

#' Bulk covariance for unnormalised weights.
#'
#' `sum_j w_j Sigma_j`, plus with `between = TRUE` the term
#' `sum_j w_j m_j m_j' - (M w)(M w)' / total`. With `sum(weights) = 1`
#' and `total = 1` this is the covariance of one random cell.
#'
#' @param weights Non-negative weights, one per cell type.
#' @param sigma Array `G x G x J`, or matrix `G x J` of variances
#'   (independence model).
#' @param m Mean signature, `G x J`.
#' @param between Add the between-type term.
#' @param total Number of cells the weights sum to.
#' @return A `G x G` matrix, or a vector of `G` variances.
bulk_covariance <- function(
  weights,
  sigma,
  m,
  between = TRUE,
  total = sum(weights)
) {
  u <- drop(m %*% weights)
  if (is.matrix(sigma)) {
    out <- drop(sigma %*% weights)
    if (between) {
      out <- out + drop(m^2 %*% weights) - u^2 / total
    }
    return(out)
  }
  n_genes <- nrow(m)
  out <- matrix(0, n_genes, n_genes)
  for (j in seq_along(weights)) {
    if (weights[[j]] != 0) {
      out <- out + weights[[j]] * sigma[,, j]
    }
  }
  if (between) {
    out <- out +
      tcrossprod(m * rep(sqrt(weights), each = n_genes)) -
      tcrossprod(u) / total
  }
  (out + t(out)) / 2
}

#' Factorise a bulk covariance.
#'
#' @param covariance Output of [bulk_covariance()].
#' @param inverse Also return the explicit inverse (needed for traces).
#' @return A list with `whiten(x)` (maps `x` to `R^-T x`, so squared
#'   norms are quadratic forms in the inverse), `solve(x)`, `log_det` and
#'   `inverse` (matrix, vector of reciprocals, or `NULL`).
factorise_covariance <- function(covariance, inverse = FALSE) {
  if (!is.matrix(covariance)) {
    if (any(!is.finite(covariance)) || any(covariance <= 0)) {
      stop("A bulk variance is not positive.")
    }
    root <- sqrt(covariance)
    return(list(
      whiten = function(x) x / root,
      solve = function(x) x / covariance,
      log_det = sum(log(covariance)),
      inverse = if (inverse) 1 / covariance else NULL
    ))
  }
  r <- chol(covariance)
  list(
    whiten = function(x) backsolve(r, x, transpose = TRUE),
    solve = function(x) backsolve(r, backsolve(r, x, transpose = TRUE)),
    log_det = 2 * sum(log(diag(r))),
    inverse = if (inverse) chol2inv(r) else NULL
  )
}

#' Quadratic form and trace of the derivative of the bulk covariance.
#'
#' For each cell type `j`, with `dU_j` the derivative of
#' [bulk_covariance()] in `w_j`, returns `a' dU_j a` and
#' `tr(U^-1 dU_j)`. With `tied = TRUE` the total `n = sum(w)` moves with
#' `w_j` (RNA-Sieve scale); with `tied = FALSE` the total is held at 1
#' and `w` is a composition treated as free (DeCovarT profile score).
#' Ledger P2.
#'
#' @param a Vector `U^-1 r`.
#' @param inverse Explicit inverse of the covariance, or vector of
#'   reciprocal variances.
#' @param weights,sigma,m,between As in [bulk_covariance()].
#' @param tied See above.

# DeCovarT likelihood at a composition -------------------------------------------

#' Profile and restricted log-likelihood of DeCovarT at a composition.
#'
#' Model `y ~ N(c M p, phi S(p))`, `S(p) = sum_j p_j Sigma_j`.
#'
#' Profile: `c = m' S^-1 y / m' S^-1 m` with `m = M p`,
#' `phi = r' S^-1 r / G`, value `-G/2 log(phi) - 1/2 log det S - G/2`
#' (ledger P1), the convention of `DeCovarT::deconvolute_single_cell()`.
#'
#' Restricted: with `S` held at `S(p)`, the `J` coefficients of the mean
#' `M beta` are integrated out (Patterson and Thompson). With `r_R` the
#' generalised least-squares residual of `y` on `M`,
#' `phi_R = r_R' S^-1 r_R / (G - J)` and the value is
#' `-(G - J)/2 log(phi_R) - 1/2 log det S - 1/2 log det(M' S^-1 M)
#'   - (G - J)/2` (ledger R4). It is exact for a covariance that does not
#' depend on the mean coefficients; here `S` depends on `p`, so it is the
#' restricted likelihood of the covariance at the fitted composition.
#'
#' @param p Composition (non-negative, sums to one).
#' @param y Bulk vector.
#' @param mu Mean signature.
#' @param sigma Covariances `G x G x J`.
#' @return A list with `profile`, `restricted`, `mean_scale`,
#'   `variance_scale`, `variance_scale_restricted`, `log_det` (of
#'   `S(p)`) and `log_det_information` (of `M' S^-1 M`).
decovart_likelihood_at <- function(p, y, mu, sigma) {
  n_genes <- length(y)
  n_types <- ncol(mu)
  fac <- factorise_covariance(
    bulk_covariance(as.numeric(p), sigma, mu, between = FALSE, total = 1)
  )
  w_y <- drop(fac$whiten(as.numeric(y)))
  w_mu <- fac$whiten(mu)
  w_m <- drop(w_mu %*% as.numeric(p))
  mean_scale <- sum(w_m * w_y) / sum(w_m^2)
  rss <- sum((w_y - mean_scale * w_m)^2)
  information <- crossprod(w_mu)
  beta <- solve(information, crossprod(w_mu, w_y))
  rss_restricted <- sum((w_y - drop(w_mu %*% beta))^2)
  log_det_information <- as.numeric(
    determinant(information, logarithm = TRUE)$modulus
  )
  df <- n_genes - n_types
  list(
    profile = -0.5 * n_genes * log(rss / n_genes) -
      0.5 * fac$log_det -
      0.5 * n_genes,
    restricted = -0.5 * df * log(rss_restricted / df) -
      0.5 * fac$log_det -
      0.5 * log_det_information -
      0.5 * df,
    mean_scale = mean_scale,
    variance_scale = rss / n_genes,
    variance_scale_restricted = rss_restricted / df,
    log_det = fac$log_det,
    log_det_information = log_det_information
  )
}

# Scale models across the replicates of a slot ----------------------------------

#' Whitened terms of one bulk at a composition.
#'
#' With `W = S(p)^-1` and `m = M p`: `N = m' W y`, `D = m' W m`,
#' `Q = y' W y`, the log-determinant of `S(p)`, and the total-count scale
#' `1'y / 1'm`.
#'
#' @param p Composition.
#' @param y Bulk vector.
#' @param mu Mean signature.
#' @param sigma Covariances `G x G x J`.
#' @return A list with `N`, `D`, `Q`, `log_det`, `total_scale`.
scale_terms <- function(p, y, mu, sigma) {
  fac <- factorise_covariance(
    bulk_covariance(as.numeric(p), sigma, mu, between = FALSE, total = 1)
  )
  m <- drop(mu %*% as.numeric(p))
  w_y <- drop(fac$whiten(as.numeric(y)))
  w_m <- drop(fac$whiten(m))
  list(
    N = sum(w_m * w_y),
    D = sum(w_m^2),
    Q = sum(w_y^2),
    log_det = fac$log_det,
    total_scale = sum(y) / sum(m)
  )
}

#' Mean and variance scales of the replicates under one scale model.
#'
#' Log-likelihood of sample `i`, without `-G/2 log(2 pi)`:
#' `-G/2 log(phi_i) - 1/2 log det S_i - (Q_i - 2 c_i N_i + c_i^2 D_i) / (2 phi_i)`.
#' - `free`: `c_i = N_i / D_i`, `phi_i` the residual over `G` (the
#'   DeCovarT profile); `c_i` can be negative.
#' - `positive`: `c_i = max(c_min, N_i / D_i)` with `c_min` one
#'   thousandth of the total-count scale, `phi_i` at that `c_i`.
#' - `linked`: `c_i = a n_i`, `phi_i = a^2 n_i = a c_i` with `a` shared by
#'   the replicates. At fixed `a`, `N_i` drops out and
#'   `c_i = (-G a + sqrt(G^2 a^2 + 4 D_i Q_i)) / (2 D_i) > 0` (ledger R5);
#'   `a` maximises the summed likelihood (Brent on `log a`).
#' - `shared_c`: one `c` for the replicates, `phi_i` free; coordinate
#'   ascent between `c = sum(N_i/phi_i) / sum(D_i/phi_i)` (at least
#'   `c_min`) and `phi_i`.
#'
#' @param terms List of [scale_terms()], one per replicate.
#' @param scale Scale model.
#' @param G Number of genes.
#' @return A list with `c`, `phi`, `a` (`NA` unless `linked`),
#'   `boundary` (logical per replicate) and `loglik` (per replicate).
scales_given_terms <- function(terms, scale, G) {
  get <- function(k) vapply(terms, `[[`, numeric(1), k)
  N <- get("N")
  D <- get("D")
  Q <- get("Q")
  c_min <- 1e-3 * get("total_scale")
  rss <- function(c) Q - 2 * c * N + c^2 * D
  ll <- function(c, phi) -0.5 * G * log(phi) - 0.5 * get("log_det") - rss(c) / (2 * phi)
  a <- NA_real_
  boundary <- rep(FALSE, length(N))
  if (scale == "free") {
    c <- N / D
  } else if (scale == "positive") {
    boundary <- N / D < c_min
    c <- pmax(N / D, c_min)
  } else if (scale == "linked") {
    c_of <- function(a) (-G * a + sqrt((G * a)^2 + 4 * D * Q)) / (2 * D)
    a0 <- stats::median(pmax(Q - N^2 / D, 1e-12) / G / get("total_scale"))
    opt <- stats::optimize(
      function(log_a) sum(ll(c_of(exp(log_a)), exp(log_a) * c_of(exp(log_a)))),
      log(a0) + c(-12, 12),
      maximum = TRUE,
      tol = 1e-10
    )
    a <- exp(opt$maximum)
    c <- c_of(a)
    return(list(c = c, phi = a * c, a = a, boundary = boundary, loglik = ll(c, a * c)))
  } else {
    floor_c <- stats::median(c_min)
    phi <- rss(max(stats::median(N / D), floor_c)) / G
    c_old <- Inf
    for (it in seq_len(200L)) {
      c_new <- max(sum(N / phi) / sum(D / phi), floor_c)
      phi <- rss(c_new) / G
      if (abs(c_new - c_old) <= 1e-12 * c_new) {
        break
      }
      c_old <- c_new
    }
    c <- rep(c_new, length(N))
    boundary <- rep(c_new <= floor_c, length(N))
  }
  phi <- rss(c) / G
  list(c = c, phi = phi, a = a, boundary = boundary, loglik = ll(c, phi))
}

#' DeCovarT on the replicates of one slot under a scale model.
#'
#' `free` is the package profile fit, one call per replicate. The other
#' models alternate a composition step (each replicate by
#' `DeCovarT::deconvolute_single_cell()` at fixed scales, L-BFGS-B) and a
#' scale step ([scales_given_terms()]); the summed likelihood increases
#' at every step. The start is the non-negative least-squares composition.
#'
#' @param Y Genes by replicates.
#' @param mu Mean signature, genes by types.
#' @param sigma Covariances `G x G x J`.
#' @param scale `"free"`, `"positive"`, `"linked"` or `"shared_c"`.
#' @param depth Library size of each replicate; `shared_c` first scales
#'   the replicates to their median depth.
#' @param maxit,tol Alternation limit and relative tolerance on the
#'   summed log-likelihood.
#' @return A list with `P` (types by replicates), `samples` (one row per
#'   replicate: scales, `a`, `n_cells_hat`, log-likelihood, boundary,
#'   whitened cosine `N / sqrt(D Q)`, `info_c = D`), `loglik`, `status`,
#'   `message`, `iterations`, `elapsed_sec`.
fit_decovart_scales <- function(
  Y,
  mu,
  sigma,
  scale = c("free", "positive", "linked", "shared_c"),
  depth = NULL,
  maxit = 25L,
  tol = 1e-9
) {
  scale <- match.arg(scale)
  Y <- as.matrix(Y)
  G <- nrow(Y)
  types <- colnames(mu)
  started <- proc.time()[["elapsed"]]
  if (scale == "shared_c") {
    if (length(depth) != ncol(Y)) {
      stop("`shared_c` needs one depth per replicate.")
    }
    Y <- sweep(Y, 2L, stats::median(depth) / depth, "*")
  }
  out <- function(P, sc, terms, status, message, iterations) {
    get <- function(k) vapply(terms, `[[`, numeric(1), k)
    a <- if (is.na(sc$a)) sc$phi / sc$c else rep(sc$a, ncol(Y))
    list(
      P = P,
      samples = data.frame(
        sample_id = colnames(Y),
        mean_scale = sc$c,
        variance_scale = sc$phi,
        platform_scale = a,
        n_cells_hat = sc$c / a,
        loglik = sc$loglik,
        boundary = sc$boundary,
        cos_w = get("N") / sqrt(get("D") * get("Q")),
        info_c = get("D"),
        stringsAsFactors = FALSE
      ),
      loglik = sum(sc$loglik),
      status = status,
      message = message,
      iterations = iterations,
      elapsed_sec = proc.time()[["elapsed"]] - started
    )
  }
  failed <- function(message) {
    list(P = NULL, samples = NULL, loglik = NA_real_, status = "failed",
         message = message, iterations = NA_integer_,
         elapsed_sec = proc.time()[["elapsed"]] - started)
  }
  if (scale == "free") {
    fits <- lapply(seq_len(ncol(Y)), function(i) {
      tryCatch(
        DeCovarT::deconvolute_single_cell(Y[, i], mu, sigma, optimiser = "lbfgsb"),
        error = function(e) e
      )
    })
    bad <- vapply(fits, inherits, logical(1), "error")
    if (any(bad)) {
      return(failed(conditionMessage(fits[[which(bad)[[1L]]]])))
    }
    P <- vapply(fits, function(f) as.numeric(f$p[types]), numeric(length(types)))
    dimnames(P) <- list(types, colnames(Y))
    terms <- lapply(seq_len(ncol(Y)), function(i) scale_terms(P[, i], Y[, i], mu, sigma))
    return(out(P, scales_given_terms(terms, "free", G), terms, "ok", "", NA_integer_))
  }
  P <- vapply(
    seq_len(ncol(Y)),
    function(i) {
      x <- nnls::nnls(mu, Y[, i])$x + 1e-3
      x / sum(x)
    },
    numeric(length(types))
  )
  dimnames(P) <- list(types, colnames(Y))
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    terms <- lapply(seq_len(ncol(Y)), function(i) scale_terms(P[, i], Y[, i], mu, sigma))
    sc <- scales_given_terms(terms, scale, G)
    ll <- sum(sc$loglik)
    if (abs(ll - ll_old) <= tol * abs(ll) || it == maxit) {
      break
    }
    ll_old <- ll
    for (i in seq_len(ncol(Y))) {
      f <- tryCatch(
        DeCovarT::deconvolute_single_cell(
          Y[, i], mu, sigma,
          mean_scale = sc$c[[i]],
          variance_scale = sc$phi[[i]],
          initial_p = P[, i],
          optimiser = "lbfgsb"
        ),
        error = function(e) e
      )
      if (inherits(f, "error")) {
        return(failed(conditionMessage(f)))
      }
      P[, i] <- as.numeric(f$p[types])
    }
  }
  status <- if (any(sc$boundary)) "boundary" else "ok"
  out(P, sc, terms, status, "", it)
}
