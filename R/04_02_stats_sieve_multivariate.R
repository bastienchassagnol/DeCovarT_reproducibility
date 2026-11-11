#' Estimators built on the Gaussian bulk law with a mixture covariance.
#'
#' Notation (DeCovarT vignettes): `G` genes, `J` cell types, `M` the
#' `G x J` matrix of mean counts per cell (`mu`), `Sigma_j` the
#' within-type covariance, `C_j` the number of reference cells, `p` the
#' cellular proportions, `w = n p` unnormalised weights with `n` the
#' number of cells behind the bulk.
#'
#' Three families share one covariance:
#'
#' - [fit_sieve_chain()]: a multivariate version of RNA-Sieve
#'   (Erdmann-Pham et al., bioRxiv 2020.10.01.322867), written here. One
#'   `n` scales the mean and the covariance. Proportions come from a
#'   quadratic programme on `w >= 0`, normalised afterwards.
#' - [fit_decovart_total()]: the DeCovarT profile likelihood (free mean
#'   scale `c` and variance scale `phi`) with the between-type term of
#'   the covariance added.
#' - [fit_decovart_firth()]: the DeCovarT profile likelihood with the
#'   Firth penalty, from interior starts.
#'
#' The covariance of one cell drawn at random from the tissue is, by the
#' law of total covariance,
#' `sum_j p_j Sigma_j + sum_j p_j m_j m_j' - (M p)(M p)'`:
#' within-type covariance, plus the covariance of the type the cell
#' belongs to. A bulk of `n` such cells has `n` times that covariance
#' (ledger C1, C2). DeCovarT's `S_W(p)` keeps the first term only, which
#' holds the number of cells of each type fixed. RNA-Sieve's eq. 7 is the
#' diagonal of the full expression; the matrix form is derived here and
#' is not in that paper.
#'
#' A `sigma` given as a `G x J` matrix of variances is the independence
#' model (genes independent within a type): every covariance is then a
#' vector and no factorisation is needed.
#'
#' Requires `quadprog`, `nnls` and `DeCovarT`.

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
#' @return A list with `quadratic` and `trace`, one value per type.
covariance_derivative_terms <- function(
  a,
  inverse,
  weights,
  sigma,
  m,
  between,
  tied
) {
  n_types <- ncol(m)
  total <- if (tied) sum(weights) else 1
  u <- drop(m %*% weights)
  quadratic <- numeric(n_types)
  trace <- numeric(n_types)
  if (is.matrix(sigma)) {
    for (j in seq_len(n_types)) {
      dv <- sigma[, j]
      if (between) {
        dv <- dv + m[, j]^2 - 2 * m[, j] * u / total
        if (tied) {
          dv <- dv + u^2 / total^2
        }
      }
      quadratic[[j]] <- sum(a^2 * dv)
      trace[[j]] <- sum(inverse * dv)
    }
    return(list(quadratic = quadratic, trace = trace))
  }
  ma <- drop(crossprod(m, a))
  ua <- sum(u * a)
  inv_m <- inverse %*% m
  inv_u <- drop(inverse %*% u)
  for (j in seq_len(n_types)) {
    slice <- sigma[,, j]
    quadratic[[j]] <- sum(a * drop(slice %*% a))
    trace[[j]] <- sum(inverse * slice)
    if (between) {
      quadratic[[j]] <- quadratic[[j]] + ma[[j]]^2 - 2 * ma[[j]] * ua / total
      trace[[j]] <- trace[[j]] +
        sum(m[, j] * inv_m[, j]) -
        2 * sum(m[, j] * inv_u) / total
      if (tied) {
        quadratic[[j]] <- quadratic[[j]] + ua^2 / total^2
        trace[[j]] <- trace[[j]] + sum(u * inv_u) / total^2
      }
    }
  }
  list(quadratic = quadratic, trace = trace)
}

# Quadratic programme ----------------------------------------------------------

#' Non-negative least squares by quadratic programming.
#'
#' Minimises `||b - a w||^2` over `w >= 0` with
#' [quadprog::solve.QP()] (Goldfarb-Idnani), on columns scaled to unit
#' norm. When the scaled cross-product is not positive definite to
#' working precision, the Lawson-Hanson solver [nnls::nnls()] is used:
#' same problem, same solution.
#'
#' @param a Design, rows by cell types.
#' @param b Response.
#' @return Non-negative weights.
solve_nonnegative_ls <- function(a, b) {
  scale <- sqrt(colSums(a^2))
  if (any(!is.finite(scale)) || any(scale == 0)) {
    stop("A column of the whitened signature is zero or not finite.")
  }
  a <- a * rep(1 / scale, each = nrow(a))
  n_types <- ncol(a)
  fit <- tryCatch(
    quadprog::solve.QP(
      Dmat = crossprod(a),
      dvec = drop(crossprod(a, b)),
      Amat = diag(n_types),
      bvec = rep(0, n_types)
    )$solution,
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- nnls::nnls(a, b)$x
  }
  pmax(fit, 0) / scale
}

#' Weights of one bulk by iterated generalised least squares.
#'
#' Freezes the bulk covariance at the current weights, solves the
#' quadratic programme `min (y - M w)' U^-1 (y - M w)`, `w >= 0`, and
#' repeats until the proportions move by less than `tol`. With `huber`,
#' each solve is itself reweighted: whitened residuals larger than
#' `huber * MAD` get the weight `huber * MAD / |e|`.
#'
#' The iteration is a fixed point, not a descent: with a nearly singular
#' covariance it can cycle. After `undamped` passes the new weights are
#' averaged with the previous ones. A solve that returns all-zero
#' weights (every column of the whitened signature negatively related
#' to the whitened bulk) stops the iteration at the previous weights.
#' Both cases are reported by `converged = FALSE`.
#'
#' @param y Bulk vector.
#' @param m Mean signature in use, `G x J`.
#' @param sigma Covariances or variances (see [bulk_covariance()]).
#' @param d Gene-wise variance scale, or `NULL` for 1.
#' @param between Add the between-type term.
#' @param start Starting weights; `NULL` uses non-negative least squares.
#' @param huber `NULL`, or the Huber constant (1.345).
#' @param tol,max_pass Stopping rule on the largest change of a
#'   proportion.
#' @param undamped Passes before the damping starts.
#' @return A list with `w`, `passes`, `converged`.
sieve_weights <- function(
  y,
  m,
  sigma,
  d = NULL,
  between = TRUE,
  start = NULL,
  huber = NULL,
  tol = 1e-4,
  max_pass = 30L,
  undamped = 5L
) {
  root_d <- if (is.null(d)) 1 else sqrt(d)
  w <- if (is.null(start)) solve_nonnegative_ls(m, y) else start
  if (sum(w) <= 0) {
    stop("The starting weights are all zero.")
  }
  floor_w <- 1e-8 * sum(w)
  converged <- FALSE
  passes <- 0L
  while (passes < max_pass && !converged) {
    passes <- passes + 1L
    # A type at zero keeps a trace weight in the covariance, so the
    # matrix stays positive definite when the others are rank deficient.
    fac <- factorise_covariance(
      bulk_covariance(pmax(w, floor_w), sigma, m, between)
    )
    a <- fac$whiten(m / root_d)
    b <- drop(fac$whiten(y / root_d))
    w_new <- solve_nonnegative_ls(a, b)
    if (!is.null(huber)) {
      for (inner in seq_len(10L)) {
        e <- b - drop(a %*% w_new)
        s <- stats::mad(e, center = 0)
        if (!is.finite(s) || s <= 0) {
          break
        }
        omega <- sqrt(pmin(1, huber * s / pmax(abs(e), 1e-300)))
        w_in <- solve_nonnegative_ls(a * omega, b * omega)
        if (sum(w_in) <= 0 || sum(w_new) <= 0) {
          break
        }
        moved <- max(abs(w_in / sum(w_in) - w_new / sum(w_new)))
        w_new <- w_in
        if (moved < tol) {
          break
        }
      }
    }
    if (sum(w_new) <= 0) {
      break
    }
    converged <- max(abs(w_new / sum(w_new) - w / sum(w))) < tol
    w <- if (passes > undamped && !converged) (w + w_new) / 2 else w_new
  }
  list(w = w, passes = passes, converged = converged)
}

#' Refined mean signature of one bulk.
#'
#' Minimises, with the bulk covariance `U` frozen,
#' `(y - M w)' U^-1 (y - M w) + sum_j C_j (mu_j - m_j)' Sigma_j^-1 (mu_j - m_j)`:
#' the bulk term plus the sampling error of the reference means. The
#' minimiser is `m_j = mu_j + (w_j / C_j) Sigma_j (U + Q)^-1 (y - mu w)`
#' with `Q = sum_j (w_j^2 / C_j) Sigma_j` (ledger S1); no `Sigma_j` is
#' inverted. Negative entries are set to zero, as the released RNA-Sieve
#' code does.
#'
#' @param y Bulk vector.
#' @param mu Reference mean signature.
#' @param m Mean signature at which `U` is evaluated.
#' @param w Current weights.
#' @param sigma,d,between As in [sieve_weights()].
#' @param n_cells Reference cells per type, `C_j`.
#' @return The refined `G x J` signature.
sieve_means <- function(y, mu, m, w, sigma, d, between, n_cells) {
  n_genes <- nrow(mu)
  n_types <- ncol(mu)
  w_pos <- pmax(w, 1e-8 * sum(w))
  u <- bulk_covariance(w_pos, sigma, m, between)
  share <- w^2 / n_cells
  resid <- y - drop(mu %*% w)
  if (is.matrix(sigma)) {
    if (!is.null(d)) {
      u <- u * d
    }
    z <- resid / (u + drop(sigma %*% share))
    out <- mu + sigma * rep(w / n_cells, each = n_genes) * z
    return(pmax(out, 0))
  }
  if (!is.null(d)) {
    u <- u * tcrossprod(sqrt(d))
  }
  for (j in seq_len(n_types)) {
    if (share[[j]] != 0) {
      u <- u + share[[j]] * sigma[,, j]
    }
  }
  z <- factorise_covariance(u)$solve(resid)
  out <- mu
  for (j in seq_len(n_types)) {
    if (w[[j]] != 0) {
      out[, j] <- mu[, j] + (w[[j]] / n_cells[[j]]) * drop(sigma[,, j] %*% z)
    }
  }
  pmax(out, 0)
}

#' Log-likelihood of one bulk under the tied scale, with its gradient.
#'
#' `y ~ N(M w, D^1/2 U(w) D^1/2)`, `U` from [bulk_covariance()]. The
#' value is `-1/2 r' U^-1 r - 1/2 log det U - 1/2 sum(log d)` with
#' `r = (y - M w) / sqrt(d)`; the constant in `2 pi` is dropped.
#'
#' @param w Weights (all positive).
#' @param y,m,sigma,d,between As in [sieve_weights()].
#' @return A list with `value` and `gradient`.
sieve_loglik <- function(w, y, m, sigma, d = NULL, between = TRUE) {
  root_d <- if (is.null(d)) rep(1, length(y)) else sqrt(d)
  fac <- factorise_covariance(
    bulk_covariance(w, sigma, m, between),
    inverse = TRUE
  )
  r <- (y - drop(m %*% w)) / root_d
  a <- drop(fac$solve(r))
  terms <- covariance_derivative_terms(
    a,
    fac$inverse,
    w,
    sigma,
    m,
    between,
    tied = TRUE
  )
  list(
    value = -0.5 * sum(r * a) - 0.5 * fac$log_det - sum(log(root_d)),
    gradient = drop(crossprod(m / root_d, a)) +
      0.5 * terms$quadratic -
      0.5 * terms$trace
  )
}

#' Maximise [sieve_loglik()] over non-negative weights.
#'
#' L-BFGS-B from `start`, on weights divided by `sum(start)` so the
#' parameters are of order one. The lower bound is a small positive
#' share of the total, which keeps the covariance positive definite.
#'
#' @param start Starting weights.
#' @param y,m,sigma,d,between As in [sieve_weights()].
#' @param itmax Largest number of iterations.
#' @return A list with `w`, `loglik`, `loglik_start`, `convergence`,
#'   `evaluations`.
sieve_full_likelihood <- function(
  start,
  y,
  m,
  sigma,
  d = NULL,
  between = TRUE,
  itmax = 50L
) {
  unit <- sum(start)
  lower <- 1e-8
  cache <- new.env(parent = emptyenv())
  evaluate <- function(x) {
    if (!identical(cache$x, x)) {
      cache$res <- tryCatch(
        sieve_loglik(x * unit, y, m, sigma, d, between),
        error = function(e) NULL
      )
      cache$x <- x
      cache$n <- (cache$n %||% 0L) + 1L
    }
    cache$res
  }
  x0 <- pmax(start / unit, lower)
  at_start <- evaluate(x0)
  if (is.null(at_start)) {
    stop("The likelihood cannot be evaluated at the starting weights.")
  }
  fit <- stats::optim(
    par = x0,
    fn = function(x) {
      res <- evaluate(x)
      if (is.null(res) || !is.finite(res$value)) -1e300 else res$value
    },
    gr = function(x) {
      res <- evaluate(x)
      if (is.null(res)) rep(0, length(x)) else res$gradient * unit
    },
    method = "L-BFGS-B",
    lower = lower,
    control = list(fnscale = -1, maxit = itmax, factr = 1e7)
  )
  list(
    w = fit$par * unit,
    loglik = fit$value,
    loglik_start = at_start$value,
    convergence = fit$convergence,
    evaluations = cache$n
  )
}

# Multivariate RNA-Sieve chain ---------------------------------------------------

#' Multivariate RNA-Sieve on the bulks of one time point.
#'
#' Model, for bulk `i`: `y_i ~ N(M_i w_i, D^1/2 U(w_i, M_i) D^1/2)` and
#' reference means `mu_j ~ N(m_j, Sigma_j / C_j)`. Stages, each stored:
#'
#' - `qp_within`: weights by [sieve_weights()], reference fixed,
#'   within-type covariance only;
#' - `qp`: the same with the between-type term;
#' - `qp_huber`: `qp` with Huber weights on the whitened residuals;
#' - `qp_means`: `qp`, alternated with [sieve_means()];
#' - `qp_means_var`: `qp_means`, alternated with the gene-wise variance
#'   scale `d_g = mean_i r_ig^2 / U_i,gg` over the bulks (the diagonal of
#'   the DECALS moment regression, Cai et al., JASA 2024; correlations
#'   are kept from the reference);
#' - `full_lik`: from that state, [sieve_full_likelihood()] in the
#'   weights, means and `d` fixed.
#'
#' Departures from the RNA-Sieve paper: the means are refined per bulk,
#' not shared by the bulks; no variance floor of `+1`; no cell trimming
#' and no detection filter (the panel is already gated); proportions are
#' weights divided by their sum instead of a simplex constraint.
#'
#' @param bulk Raw counts, genes by bulk samples, rows as `mu`.
#' @param mu Reference mean signature, `G x J`.
#' @param sigma Covariances `G x G x J`, or variances `G x J`.
#' @param n_cells Reference cells per type.
#' @param workers Forked workers over the bulks.
#' @param tol,max_pass Passed to [sieve_weights()].
#' @param max_outer Largest number of alternations.
#' @param d_floor Lower bound of `d_g`, as a share of its median.
#' @return A list with `estimates` (`stage`, `sample_id`, `cell_type`,
#'   `proportion`), `diagnostics` (one row per stage and sample) and `d`.
fit_sieve_chain <- function(
  bulk,
  mu,
  sigma,
  n_cells,
  workers = 1L,
  tol = 1e-4,
  max_pass = 30L,
  max_outer = 10L,
  d_floor = 0.01
) {
  types <- colnames(mu)
  samples <- colnames(bulk)
  n_cells <- as.numeric(n_cells[types])
  if (anyNA(n_cells) || any(n_cells <= 0)) {
    stop("`n_cells` must give a positive count for every type of `mu`.")
  }
  over_bulks <- function(fun) {
    out <- parallel::mclapply(
      seq_along(samples),
      function(i) tryCatch(fun(i), error = function(e) e),
      mc.cores = min(workers, length(samples)),
      mc.preschedule = FALSE
    )
    bad <- vapply(out, inherits, logical(1), what = "error")
    bad <- bad | vapply(out, inherits, logical(1), what = "try-error")
    if (any(bad)) {
      stop(
        "Sieve stage failed for ",
        paste(samples[bad], collapse = ", "),
        ": ",
        if (inherits(out[[which(bad)[[1L]]]], "error")) {
          conditionMessage(out[[which(bad)[[1L]]]])
        } else {
          as.character(out[[which(bad)[[1L]]]])
        }
      )
    }
    out
  }
  relative_shift <- function(m) sqrt(sum((m - mu)^2) / sum(mu^2))
  record <- function(stage, fits) {
    list(
      estimates = do.call(
        rbind,
        lapply(seq_along(samples), function(i) {
          data.frame(
            stage = stage,
            sample_id = samples[[i]],
            cell_type = types,
            proportion = fits[[i]]$w / sum(fits[[i]]$w),
            stringsAsFactors = FALSE
          )
        })
      ),
      diagnostics = do.call(
        rbind,
        lapply(seq_along(samples), function(i) {
          f <- fits[[i]]
          data.frame(
            stage = stage,
            sample_id = samples[[i]],
            n_hat = sum(f$w),
            passes = f$passes %||% NA_integer_,
            converged = f$converged %||% NA,
            mean_shift = if (is.null(f$m)) 0 else relative_shift(f$m),
            loglik = f$loglik %||% NA_real_,
            loglik_start = f$loglik_start %||% NA_real_,
            stringsAsFactors = FALSE
          )
        })
      )
    )
  }

  # Stages with the reference fixed, and the per-bulk mean refinement.
  first <- over_bulks(function(i) {
    y <- bulk[, i]
    within <- sieve_weights(
      y,
      mu,
      sigma,
      between = FALSE,
      tol = tol,
      max_pass = max_pass
    )
    qp <- sieve_weights(
      y,
      mu,
      sigma,
      start = within$w,
      tol = tol,
      max_pass = max_pass
    )
    hub <- sieve_weights(
      y,
      mu,
      sigma,
      start = qp$w,
      huber = 1.345,
      tol = tol,
      max_pass = max_pass
    )
    w <- qp$w
    m <- mu
    outer <- 0L
    done <- FALSE
    while (outer < max_outer && !done) {
      outer <- outer + 1L
      m <- sieve_means(y, mu, m, w, sigma, NULL, TRUE, n_cells)
      step <- sieve_weights(
        y,
        m,
        sigma,
        start = w,
        tol = tol,
        max_pass = 3L
      )
      done <- max(abs(step$w / sum(step$w) - w / sum(w))) < tol
      w <- step$w
    }
    list(
      within = within,
      qp = qp,
      huber = hub,
      means = list(w = w, m = m, passes = outer, converged = done)
    )
  })
  stages <- list(
    record("qp_within", lapply(first, `[[`, "within")),
    record("qp", lapply(first, `[[`, "qp")),
    record("qp_huber", lapply(first, `[[`, "huber")),
    record("qp_means", lapply(first, `[[`, "means"))
  )

  # Gene-wise variance scale, pooled over the bulks.
  state <- lapply(first, `[[`, "means")
  d <- rep(1, nrow(mu))
  outer <- 0L
  done <- FALSE
  while (outer < max_outer && !done) {
    outer <- outer + 1L
    ratios <- vapply(
      over_bulks(function(i) {
        s <- state[[i]]
        u <- bulk_covariance(
          pmax(s$w, 1e-8 * sum(s$w)),
          sigma,
          s$m,
          TRUE
        )
        (bulk[, i] - drop(s$m %*% s$w))^2 / (if (is.matrix(u)) diag(u) else u)
      }),
      identity,
      numeric(nrow(mu))
    )
    d <- rowMeans(ratios)
    d <- pmax(d, d_floor * stats::median(d))
    moved <- over_bulks(function(i) {
      s <- state[[i]]
      y <- bulk[, i]
      m <- sieve_means(y, mu, s$m, s$w, sigma, d, TRUE, n_cells)
      step <- sieve_weights(
        y,
        m,
        sigma,
        d = d,
        start = s$w,
        tol = tol,
        max_pass = 3L
      )
      list(
        w = step$w,
        m = m,
        change = max(abs(step$w / sum(step$w) - s$w / sum(s$w)))
      )
    })
    done <- max(vapply(moved, `[[`, numeric(1), "change")) < tol
    state <- lapply(moved, function(s) {
      list(w = s$w, m = s$m, passes = outer, converged = done)
    })
  }
  stages <- c(stages, list(record("qp_means_var", state)))

  full <- over_bulks(function(i) {
    s <- state[[i]]
    fit <- sieve_full_likelihood(
      pmax(s$w, 1e-6 * sum(s$w)),
      bulk[, i],
      s$m,
      sigma,
      d = d,
      between = TRUE
    )
    list(
      w = fit$w,
      m = s$m,
      passes = fit$evaluations,
      converged = fit$convergence == 0L,
      loglik = fit$loglik,
      loglik_start = fit$loglik_start
    )
  })
  stages <- c(stages, list(record("full_lik", full)))

  list(
    estimates = do.call(rbind, lapply(stages, `[[`, "estimates")),
    diagnostics = do.call(rbind, lapply(stages, `[[`, "diagnostics")),
    d = stats::setNames(d, rownames(mu))
  )
}

# DeCovarT profile likelihood with the between-type term ------------------------

#' Profile log-likelihood of DeCovarT at a composition, with its score.
#'
#' `y ~ N(c M p, phi S(p))`, `S` from [bulk_covariance()] at `total = 1`.
#' Both scales are in closed form: `c = m' S^-1 y / m' S^-1 m` with
#' `m = M p`, and `phi = r' S^-1 r / G` (ledger P1). The value is
#' `-G/2 log(phi) - 1/2 log det S - G/2`, the convention of
#' `DeCovarT::deconvolute_single_cell()`. The score treats `p` as free
#' (ledger P2); the caller projects it on the simplex.
#'
#' @param p Composition (positive, sums to one).
#' @param y Bulk vector.
#' @param mu Mean signature.
#' @param sigma Covariances `G x G x J`, or variances `G x J`.
#' @param between Add the between-type term.
#' @return A list with `value`, `score`, `mean_scale`, `variance_scale`
#'   and `log_det` (of `S(p)`).
decovart_profile_total <- function(p, y, mu, sigma, between = TRUE) {
  n_genes <- length(y)
  fac <- factorise_covariance(
    bulk_covariance(p, sigma, mu, between, total = 1),
    inverse = TRUE
  )
  mean_vec <- drop(mu %*% p)
  s_m <- drop(fac$solve(mean_vec))
  mean_scale <- sum(s_m * y) / sum(s_m * mean_vec)
  resid <- y - mean_scale * mean_vec
  a <- drop(fac$solve(resid))
  variance_scale <- sum(resid * a) / n_genes
  if (!is.finite(mean_scale) || mean_scale <= 0 || variance_scale <= 0) {
    stop("A profiled scale is not positive.")
  }
  terms <- covariance_derivative_terms(
    a,
    fac$inverse,
    p,
    sigma,
    mu,
    between,
    tied = FALSE
  )
  list(
    value = -0.5 * n_genes * log(variance_scale) -
      0.5 * fac$log_det -
      0.5 * n_genes,
    score = (mean_scale / variance_scale) *
      drop(crossprod(mu, a)) +
      0.5 * terms$quadratic / variance_scale -
      0.5 * terms$trace,
    mean_scale = mean_scale,
    variance_scale = variance_scale,
    log_det = fac$log_det
  )
}

#' DeCovarT profile fit with the between-type covariance term.
#'
#' Same optimiser as `optimiser = "lbfgsb"` of
#' `DeCovarT::deconvolute_single_cell()`: L-BFGS-B on non-negative
#' weights, renormalised inside the objective, from the non-negative
#' least-squares start. With `between = FALSE` the objective is the
#' package's own, which gives a regression test.
#'
#' @param y Bulk vector, genes as the rows of `mu`.
#' @param mu Mean signature.
#' @param sigma Covariances `G x G x J`, or variances `G x J`.
#' @param between Add the between-type term.
#' @param initial_p Optional start.
#' @param epsilon,itmax Tolerance and iteration cap, as in the package.
#' @return A list with `p`, `loglik`, `mean_scale`, `variance_scale`,
#'   `convergence`, `evaluations`, `message`.
fit_decovart_total <- function(
  y,
  mu,
  sigma,
  between = TRUE,
  initial_p = NULL,
  epsilon = 1e-8,
  itmax = 100L
) {
  types <- colnames(mu)
  floor_p <- 100 * .Machine$double.eps
  if (is.null(initial_p)) {
    initial_p <- nnls::nnls(mu, y)$x
    if (sum(initial_p) <= 0) {
      initial_p <- rep(1, ncol(mu))
    }
  }
  p0 <- pmax(initial_p / sum(initial_p), floor_p)
  p0 <- p0 / sum(p0)
  cache <- new.env(parent = emptyenv())
  evaluate <- function(w) {
    if (!identical(cache$w, w)) {
      p <- pmax(w / sum(w), floor_p)
      cache$res <- tryCatch(
        decovart_profile_total(p / sum(p), y, mu, sigma, between),
        error = function(e) NULL
      )
      cache$w <- w
      cache$n <- (cache$n %||% 0L) + 1L
    }
    cache$res
  }
  fit <- stats::optim(
    par = p0,
    fn = function(w) {
      res <- evaluate(w)
      if (is.null(res) || !is.finite(res$value)) -1e12 else res$value
    },
    gr = function(w) {
      res <- evaluate(w)
      if (is.null(res)) {
        return(rep(0, length(w)))
      }
      total <- sum(w)
      (res$score - sum(res$score * w / total)) / total
    },
    method = "L-BFGS-B",
    lower = 0,
    upper = Inf,
    control = list(
      fnscale = -1,
      maxit = itmax,
      factr = epsilon / .Machine$double.eps
    )
  )
  p <- fit$par / sum(fit$par)
  at <- evaluate(fit$par)
  list(
    p = stats::setNames(p, types),
    loglik = fit$value,
    mean_scale = if (is.null(at)) NA_real_ else at$mean_scale,
    variance_scale = if (is.null(at)) NA_real_ else at$variance_scale,
    convergence = fit$convergence,
    evaluations = cache$n,
    message = as.character(fit$message %||% "")
  )
}

# Firth-penalised DeCovarT from interior starts ----------------------------------

#' Firth-penalised DeCovarT profile fit.
#'
#' The objective of `DeCovarT::fit_firth_single_cell()`: the profile
#' log-likelihood plus half the log-determinant of the observed
#' information of the `J - 1` ILR coordinates. That function starts from
#' its own Newton fit and takes no start; a start on a face has a
#' singular information, an infinite penalty, and cannot move. Here the
#' same objective is rebuilt from the exported profile functions and
#' maximised by BFGS from each row of `starts`; the fit with the highest
#' penalised value is kept. The penalty gradient is a central
#' difference, as in the package.
#'
#' @param y Bulk vector.
#' @param mu Mean signature.
#' @param sigma Covariances `G x G x J`.
#' @param starts Matrix of interior compositions, one per row.
#' @param epsilon,itmax Tolerance and iteration cap.
#' @param step Step of the central difference on the ILR coordinates.
#' @return A list with `p`, `loglik`, `penalty`, `penalised_loglik`,
#'   `convergence`, `start` (row kept), `hessians` (evaluations).
fit_decovart_firth <- function(
  y,
  mu,
  sigma,
  starts,
  epsilon = 1e-8,
  itmax = 80L,
  step = 1e-4
) {
  types <- colnames(mu)
  count <- new.env(parent = emptyenv())
  count$n <- 0L
  penalty <- function(z) {
    count$n <- count$n + 1L
    h <- tryCatch(
      DeCovarT::hessian_profile_loglik_single_cell_ilr(z, y, mu, sigma),
      error = function(e) NULL
    )
    if (is.null(h) || any(!is.finite(h))) {
      return(-Inf)
    }
    ev <- eigen(-(h + t(h)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(ev)) || min(ev) <= 1e-8) {
      return(-Inf)
    }
    0.5 * sum(log(ev))
  }
  objective <- function(z) {
    ll <- tryCatch(
      DeCovarT::profile_loglik_single_cell_ilr(z, y, mu, sigma),
      error = function(e) -Inf
    )
    pen <- penalty(z)
    if (!is.finite(ll) || !is.finite(pen)) -1e10 else ll + pen
  }
  gradient <- function(z) {
    g_ll <- DeCovarT::gradient_profile_loglik_single_cell_ilr(z, y, mu, sigma)
    g_pen <- vapply(
      seq_along(z),
      function(k) {
        e <- replace(numeric(length(z)), k, step)
        (penalty(z + e) - penalty(z - e)) / (2 * step)
      },
      numeric(1)
    )
    if (any(!is.finite(g_pen))) {
      g_pen <- rep(0, length(z))
    }
    g_ll + g_pen
  }
  fits <- lapply(seq_len(nrow(starts)), function(k) {
    z0 <- DeCovarT::isometric_log_ratio(
      DeCovarT::repair_simplex(starts[k, ], open = TRUE)
    )
    tryCatch(
      stats::optim(
        par = z0,
        fn = objective,
        gr = gradient,
        method = "BFGS",
        control = list(fnscale = -1, reltol = epsilon, maxit = itmax)
      ),
      error = function(e) NULL
    )
  })
  values <- vapply(
    fits,
    function(f) if (is.null(f)) -Inf else f$value,
    numeric(1)
  )
  if (!any(is.finite(values)) || max(values) <= -1e10) {
    stop("The penalised objective is not finite from any start.")
  }
  best <- fits[[which.max(values)]]
  p <- DeCovarT::repair_simplex(
    DeCovarT::isometric_logistic(best$par),
    open = TRUE,
    nms = types
  )
  pen <- penalty(best$par)
  list(
    p = p,
    loglik = best$value - pen,
    penalty = pen,
    penalised_loglik = best$value,
    convergence = best$convergence,
    start = which.max(values),
    hessians = count$n
  )
}

# Shrinkage level ------------------------------------------------------------------

#' Move a shrinkage covariance to another shrinkage level.
#'
#' The Schafer-Strimmer estimate is
#' `D^1/2 {(1 - lambda) R + lambda I} D^1/2`: variances untouched,
#' off-diagonal terms multiplied by `1 - lambda`. Going from
#' `lambda_from` to `lambda_to` therefore multiplies every off-diagonal
#' term by `(1 - lambda_to) / (1 - lambda_from)` (ledger S2).
#'
#' @param sigma Covariances `G x G x J` estimated at `lambda_from`.
#' @param lambda_from Shrinkage of each type, below 1.
#' @param lambda_to Target shrinkage in `[0, 1]`, one value.
#' @return An array like `sigma`.
rescale_shrinkage <- function(sigma, lambda_from, lambda_to) {
  if (any(lambda_from >= 1)) {
    stop("`lambda_from` must be below 1.")
  }
  n_genes <- dim(sigma)[[1L]]
  on_diagonal <- (seq_len(n_genes) - 1L) * n_genes + seq_len(n_genes)
  for (j in seq_len(dim(sigma)[[3L]])) {
    idx <- (j - 1) * n_genes^2 + on_diagonal
    keep <- sigma[idx]
    slice <- (j - 1) * n_genes^2 + seq_len(n_genes^2)
    sigma[slice] <- sigma[slice] * (1 - lambda_to) / (1 - lambda_from[[j]])
    sigma[idx] <- keep
  }
  sigma
}
