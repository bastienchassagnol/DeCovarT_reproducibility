#' Tail statements for the bulk-versus-reference consistency filter.
#'
#' Tier (iii) of the gene selection, "Bulk consistency filter". For gene
#' `g`, bulk sample `i` and reference types `j = 1, ..., J`, the
#' statistic is
#'
#'   `D_gi = y_gi / a_i - sum_j p_ji mu_gj`,
#'
#' with `y_gi` the bulk count, `a_i` the plug-in scale (effective number
#' of cells), `p_ji` the plug-in composition ([plugin_composition_demixsc()])
#' and `mu_gj` the mean count per cell of type `j` over the `C_j`
#' reference cells. Under (A1) independent NB cells with mean `mu_gj`
#' and size `psi_gj`, (A2) `m_ji = p_ji n_i` bulk cells of type `j`,
#' (A3) the same law on both platforms and (A4) known parameters,
#'
#'   `Var D = sum_j p_j sigma_gj^2 / n_i + sum_j p_j^2 sigma_gj^2 / C_j`,
#'   `sigma_gj^2 = mu_gj + mu_gj^2 / psi_gj`,
#'
#' and the cumulant generating function of `D` is
#'
#'   `K_D(u) = sum_j m_ji log M_j(u / n_i) + sum_j C_j log M_j(-u p_j / C_j)`,
#'   `log M_j(l) = -psi_gj log(1 - mu_gj (e^l - 1) / psi_gj)`.
#'
#' Statements, at level `alpha / 2` per tail (`level`, default 0.01):
#'
#' | method    | gene   | module (sum of standardised `D_g`)          |
#' |-----------|--------|---------------------------------------------|
#' | `clt`     | T1     | T1, independent genes, `Var S = N`          |
#' | `clt_log` | T1-log | not defined                                 |
#' | `clt_dep` | --     | T2, dependency graph, `Var S = N + 2 sum Corr` |
#' | `poisson` | T3     | T3, Poisson CGFs, independent genes         |
#' | `nb`      | T4     | T4, NB CGFs, independent genes              |
#' | `nb_dep`  | --     | T5, NB CGFs after a proper colouring (Janson 2004) |
#'
#' T1, T1-log and T2 are approximations, T4 and T5 finite-sample bounds
#' under their assumptions; T3 is a bound only for Poisson counts. With
#' `a_i` and `p_i` plugged in, none of them accounts for estimated
#' parameters. The extra overdispersion between platforms is fixed at 1
#' (assumption A3). The module statements also assume (A5) a sparse
#' within-cell correlation on the edges of a dependency graph and (A6)
#' independence between gene sets with no edge between them.
#'
#' Port of `scripts/figures/fig_03_1_tail_bounds.R` of the DeCovarT
#' repository, generalised to type-specific `psi_gj` and to a Poisson
#' limit `psi = Inf`; deliberately duplicated, no cross-repository
#' dependency. Requires `R/02_00_utils.R`, `R/02_03_plugin_composition.R`.
#' `quantreg` (quantile band) is not in `renv.lock`. Written, not run.

#' Method-of-moments NB parameters per gene and type.
#'
#' Per type `j`, the mean `mu_gj` and unbiased variance `v_gj` of the raw
#' counts per cell, then `psi_gj = mu_gj^2 / (v_gj - mu_gj)` when
#' `v_gj > mu_gj`, `Inf` (Poisson) otherwise. Raw counts are used, not
#' size-normalised ones, because a bulk is a sum of raw cell profiles;
#' the library-size variation between cells is part of `v_gj`.
#'
#' @param counts Genes by cells raw counts, gene names as row names.
#' @param cell_type Cell-type label per column.
#' @return A list with `mu` and `psi` (genes by types matrices) and
#'   `n_cells` (named cells per type).
type_moments_nb <- function(counts, cell_type) {
  cell_type <- droplevels(factor(cell_type))
  if (length(cell_type) != ncol(counts)) {
    stop("`cell_type` must have one label per column of `counts`.")
  }
  if (is.null(rownames(counts))) {
    stop("`counts` needs gene names as row names.")
  }
  types <- levels(cell_type)
  n_cells <- as.numeric(table(cell_type)[types])
  names(n_cells) <- types
  if (any(n_cells < 2)) {
    stop("Every type needs at least two cells to estimate a variance.")
  }
  mu <- matrix(
    NA_real_,
    nrow(counts),
    length(types),
    dimnames = list(rownames(counts), types)
  )
  psi <- mu
  for (tp in types) {
    x_j <- counts[, cell_type == tp, drop = FALSE]
    n_j <- n_cells[[tp]]
    m1 <- Matrix::rowMeans(x_j)
    m2 <- Matrix::rowMeans(x_j * x_j)
    v <- pmax(m2 - m1^2, 0) * n_j / (n_j - 1)
    mu[, tp] <- m1
    psi[, tp] <- ifelse(v > m1, m1^2 / (v - m1), Inf)
  }
  list(mu = mu, psi = psi, n_cells = n_cells)
}

#' Greedy proper colouring of a graph (Welsh-Powell order).
#'
#' Vertices are coloured in decreasing degree; each takes the smallest
#' colour not used by a coloured neighbour. The number of colours `k`
#' satisfies `k <= max degree + 1` and `k <= (1 + sqrt(8 |E| + 1)) / 2`.
#'
#' @param adjacency Symmetric 0/1 (or logical) matrix, zero diagonal.
#' @return Integer colour per vertex, from 1 to `k`.
greedy_colouring <- function(adjacency) {
  adjacency <- as.matrix(adjacency) != 0
  if (nrow(adjacency) != ncol(adjacency)) {
    stop("`adjacency` must be a square matrix.")
  }
  if (!isSymmetric(unname(adjacency))) {
    stop("`adjacency` must be symmetric (undirected graph).")
  }
  diag(adjacency) <- FALSE
  n <- nrow(adjacency)
  deg <- rowSums(adjacency)
  colour <- integer(n)
  for (v in order(-deg, seq_len(n))) {
    used <- colour[adjacency[v, ]]
    colour[v] <- min(setdiff(seq_len(n), used))
  }
  colour
}

#' Sparse working dependency graph of a gene module.
#'
#' Candidate edges are gene pairs whose largest absolute within-type
#' correlation over the types, `max_j |r_ghj|`, is at least
#' `min_abs_cor`. They are added in decreasing order of that value while
#' both genes have degree below `max_degree` and the number of edges
#' stays at most `floor(mean_degree N / 2)`. The graph is a working
#' assumption (A6: genes with no edge are independent), not an estimate
#' of conditional independence; the residual correlations of the pairs
#' left out should be reported.
#'
#' @param cor_list List of `J` gene by gene correlation matrices (one per
#'   type), or one matrix.
#' @param max_degree Maximum degree. Default `10`.
#' @param mean_degree Maximum mean degree `2 |E| / N`. Default `3`.
#' @param min_abs_cor Smallest absolute correlation for an edge.
#'   Default `0.1`; a tuning choice, report it.
#' @return A symmetric 0/1 adjacency matrix with the gene names.
module_dependency_graph <- function(
  cor_list,
  max_degree = 10L,
  mean_degree = 3,
  min_abs_cor = 0.1
) {
  if (is.matrix(cor_list)) {
    cor_list <- list(cor_list)
  }
  abs_cor <- Reduce(
    pmax,
    lapply(cor_list, function(r) {
      r <- abs(as.matrix(r))
      r[!is.finite(r)] <- 0
      r
    })
  )
  n <- nrow(abs_cor)
  genes <- rownames(abs_cor) %||% as.character(seq_len(n))
  adjacency <- matrix(0, n, n, dimnames = list(genes, genes))
  if (n < 2L) {
    return(adjacency)
  }
  pairs <- which(upper.tri(abs_cor) & abs_cor >= min_abs_cor, arr.ind = TRUE)
  if (nrow(pairs) == 0L) {
    return(adjacency)
  }
  pairs <- pairs[order(-abs_cor[pairs]), , drop = FALSE]
  max_edges <- floor(mean_degree * n / 2)
  deg <- integer(n)
  n_edges <- 0L
  for (k in seq_len(nrow(pairs))) {
    if (n_edges >= max_edges) {
      break
    }
    g <- pairs[k, 1L]
    h <- pairs[k, 2L]
    if (deg[g] < max_degree && deg[h] < max_degree) {
      adjacency[g, h] <- 1
      adjacency[h, g] <- 1
      deg[c(g, h)] <- deg[c(g, h)] + 1L
      n_edges <- n_edges + 1L
    }
  }
  adjacency
}

# Per-type log-MGF of one cell and its derivative; psi = Inf is Poisson.
lmgf_cell <- function(l, mu, psi) {
  out <- mu * expm1(l)
  nb <- is.finite(psi)
  x <- mu[nb] * expm1(l[nb]) / psi[nb]
  out[nb] <- ifelse(x < 1, -psi[nb] * log1p(-x), Inf)
  out
}
dlmgf_cell <- function(l, mu, psi) {
  out <- mu * exp(l)
  nb <- is.finite(psi)
  out[nb] <- mu[nb] * exp(l[nb]) / (1 - mu[nb] * expm1(l[nb]) / psi[nb])
  out
}

#' Variance of the consistency statistic of one gene.
#'
#' @param mu Per-cell type means of the gene (length `J`).
#' @param p Composition (length `J`).
#' @param n_bulk Cells in the bulk sample (`n_i`, or the plug-in `a_i`).
#' @param n_ref Reference cells per type (`C_j`).
#' @param psi NB size per type (scalar or length `J`); `Inf` is Poisson.
#' @return `Var D`.
consistency_variance <- function(mu, p, n_bulk, n_ref, psi = Inf) {
  psi <- rep_len(psi, length(mu))
  s2 <- mu + ifelse(is.finite(psi), mu^2 / psi, 0)
  sum(p * s2) / n_bulk + sum(p^2 * s2 / n_ref)
}

#' Cumulant generating function of the consistency statistic of one gene.
#'
#' @inheritParams consistency_variance
#' @param family `"nb"` (exact NB CGF) or `"poisson"` (Poisson CGF, `psi`
#'   ignored).
#' @return A list with `K`, `dK` (functions of `u`) and the domain
#'   `lo < 0 < hi`, or `NULL` when no type has `p_j > 0` and `mu_j > 0`.
consistency_cgf <- function(
  mu,
  p,
  n_bulk,
  n_ref,
  psi = Inf,
  family = c("nb", "poisson")
) {
  family <- match.arg(family)
  psi <- rep_len(psi, length(mu))
  keep <- p > 0 & mu > 0
  if (!any(keep)) {
    return(NULL)
  }
  mu <- mu[keep]
  p <- p[keep]
  n_ref <- n_ref[keep]
  psi <- if (family == "poisson") rep(Inf, sum(keep)) else psi[keep]
  m <- p * n_bulk
  b <- ifelse(is.finite(psi), log1p(psi / mu), Inf) # -log pi_j
  list(
    K = function(u) {
      sum(m * lmgf_cell(rep(u / n_bulk, length(mu)), mu, psi)) +
        sum(n_ref * lmgf_cell(-u * p / n_ref, mu, psi))
    },
    dK = function(u) {
      sum(p * dlmgf_cell(rep(u / n_bulk, length(mu)), mu, psi)) -
        sum(p * dlmgf_cell(-u * p / n_ref, mu, psi))
    },
    lo = max(-n_ref * b / p),
    hi = n_bulk * min(b)
  )
}

# Legendre transform sup_{side u > 0} [side u t - K(u)], t >= 0, at the
# root of K'(u) = side t (K convex).
chernoff_rate <- function(cgf, t, side = 1) {
  if (t <= 0) {
    return(0)
  }
  edge <- if (side > 0) cgf$hi else cgf$lo
  if (is.finite(edge)) {
    edge <- edge * (1 - 1e-12)
  } else {
    edge <- side
    n_double <- 0L
    while (side * cgf$dK(edge) < t) {
      edge <- 2 * edge
      n_double <- n_double + 1L
      if (n_double > 200L) {
        stop("Could not bracket the saddle point; check the CGF inputs.")
      }
    }
  }
  g <- function(u) cgf$dK(u) - side * t
  if (!(side * g(edge) > 0)) {
    stop("The CGF derivative does not reach the deviation inside its domain.")
  }
  u_star <- stats::uniroot(
    g,
    sort(c(0, edge)),
    tol = 1e-13 * abs(edge),
    maxiter = 2000L
  )$root
  side * u_star * t - cgf$K(u_star)
}

# Smallest t with exp(-rate(t)) <= level.
chernoff_threshold <- function(cgf, side, t0, level) {
  target <- log(1 / level)
  hi <- t0
  n_double <- 0L
  while (chernoff_rate(cgf, hi, side) < target) {
    hi <- 2 * hi
    n_double <- n_double + 1L
    if (n_double > 200L) {
      stop("Could not bracket the Chernoff threshold.")
    }
  }
  stats::uniroot(
    function(t) chernoff_rate(cgf, t, side) - target,
    c(0, hi),
    tol = 1e-10 * hi
  )$root
}

#' Shared checks for one gene or one module tail statement.
#'
#' @return A list with the validated arguments, `z`, and either `s` and
#'   `mbar` (gene) or `sd_g` (module). `zero` is TRUE when a gene has
#'   no positive variance of `D`.
.tail_context <- function(mu, p, n_bulk, n_ref, psi, level) {
  if (
    !is.numeric(level) || length(level) != 1L || level <= 0 || level >= 0.5
  ) {
    stop(
      "`level` must be one number in (0, 0.5), ",
      "the tail probability per side."
    )
  }
  is_module <- is.matrix(mu) && nrow(mu) > 1L
  n_types <- if (is.matrix(mu)) ncol(mu) else length(mu)
  if (length(p) != n_types || length(n_ref) != n_types) {
    stop("`p` and `n_ref` must have one value per type (column of `mu`).")
  }
  if (any(p < 0) || abs(sum(p) - 1) > 1e-6) {
    stop("`p` must be non-negative and sum to one.")
  }
  if (!is.numeric(n_bulk) || length(n_bulk) != 1L || n_bulk <= 0) {
    stop("`n_bulk` must be one positive number (cells, or the plug-in a_i).")
  }
  if (any(n_ref[p > 0] <= 0)) {
    stop("`n_ref` must be positive for every type with p_j > 0.")
  }
  z <- stats::qnorm(1 - level)
  if (!is_module) {
    mu <- as.numeric(mu)
    psi <- rep_len(as.numeric(psi), n_types)
    s <- sqrt(consistency_variance(mu, p, n_bulk, n_ref, psi))
    return(list(
      is_module = FALSE,
      n_types = n_types,
      n_genes = 1L,
      mu = mu,
      psi = psi,
      p = p,
      n_bulk = n_bulk,
      n_ref = n_ref,
      level = level,
      z = z,
      s = s,
      mbar = sum(p * mu),
      zero = !is.finite(s) || s <= 0
    ))
  }
  n_genes <- nrow(mu)
  psi <- matrix(as.numeric(psi), n_genes, n_types)
  sd_g <- vapply(
    seq_len(n_genes),
    function(g) {
      sqrt(consistency_variance(mu[g, ], p, n_bulk, n_ref, psi[g, ]))
    },
    numeric(1)
  )
  if (any(!is.finite(sd_g)) || any(sd_g <= 0)) {
    stop("Every module gene needs a positive variance of D; drop the others.")
  }
  list(
    is_module = TRUE,
    n_types = n_types,
    n_genes = n_genes,
    mu = mu,
    psi = psi,
    p = p,
    n_bulk = n_bulk,
    n_ref = n_ref,
    level = level,
    z = z,
    sd_g = sd_g,
    zero = FALSE
  )
}

.tail_rows <- function(method, statement, lower, upper) {
  data.frame(
    method = method,
    statement = statement,
    side = c("lower", "upper"),
    threshold = c(lower, upper),
    stringsAsFactors = FALSE
  )
}

.tail_na <- function(method) {
  warning("Zero variance of D: thresholds are NA.", call. = FALSE)
  .tail_rows(method, NA_character_, NA_real_, NA_real_)
}

.tail_chernoff_rows <- function(method, statement, cgf, scale, level) {
  .tail_rows(
    method,
    statement,
    chernoff_threshold(cgf, -1, scale, level),
    chernoff_threshold(cgf, 1, scale, level)
  )
}

.module_gene_cgfs <- function(ctx, family) {
  n_genes <- ctx$n_genes
  lapply(seq_len(n_genes), function(g) {
    consistency_cgf(
      ctx$mu[g, ],
      ctx$p,
      ctx$n_bulk,
      ctx$n_ref,
      ctx$psi[g, ],
      family
    )
  })
}

.module_cgf_indep <- function(cl, sd_g) {
  n_genes <- length(sd_g)
  list(
    K = function(u) {
      sum(vapply(seq_len(n_genes), function(g) cl[[g]]$K(u / sd_g[g]), 0))
    },
    dK = function(u) {
      sum(vapply(
        seq_len(n_genes),
        function(g) cl[[g]]$dK(u / sd_g[g]) / sd_g[g],
        0
      ))
    },
    lo = max(vapply(seq_len(n_genes), function(g) cl[[g]]$lo * sd_g[g], 0)),
    hi = min(vapply(seq_len(n_genes), function(g) cl[[g]]$hi * sd_g[g], 0))
  )
}

# T5: L(u) = log(k^-1 sum_c exp(sum_{g in I_c} K_{Z_g}(k u))) >= log E e^{uS}.
.module_cgf_colour <- function(cl, sd_g, adjacency) {
  n_genes <- length(sd_g)
  colour <- greedy_colouring(adjacency)
  k <- max(colour)
  classes <- split(seq_len(n_genes), colour)
  a_c <- function(u) {
    vapply(
      classes,
      function(idx) {
        sum(vapply(idx, function(g) cl[[g]]$K(k * u / sd_g[g]), 0))
      },
      0
    )
  }
  da_c <- function(u) {
    vapply(
      classes,
      function(idx) {
        sum(vapply(
          idx,
          function(g) k * cl[[g]]$dK(k * u / sd_g[g]) / sd_g[g],
          0
        ))
      },
      0
    )
  }
  list(
    K = function(u) {
      a <- a_c(u)
      max(a) + log(mean(exp(a - max(a))))
    },
    dK = function(u) {
      a <- a_c(u)
      w <- exp(a - max(a))
      sum(w * da_c(u)) / sum(w)
    },
    lo = max(vapply(seq_len(n_genes), function(g) cl[[g]]$lo * sd_g[g], 0)) / k,
    hi = min(vapply(seq_len(n_genes), function(g) cl[[g]]$hi * sd_g[g], 0)) / k
  )
}

.require_module_graph <- function(mu, adjacency) {
  if (!(is.matrix(mu) && nrow(mu) > 1L)) {
    stop("This statement applies to a module (N by J, N >= 2).")
  }
  if (is.null(adjacency)) {
    stop("`adjacency` is required for clt_dep and nb_dep.")
  }
  adjacency <- as.matrix(adjacency) != 0
  if (!all(dim(adjacency) == nrow(mu))) {
    stop("`adjacency` must be N by N, N = nrow(mu).")
  }
  adjacency
}

#' T1, Gaussian approximation (Lindeberg-Feller).
#'
#' Gene: threshold `z * sd(D)` on both tails. Module of independent
#' genes: threshold `z * sqrt(N)` for `S = sum_g D_g / sd(D_g)`.
#' Arguments are those of [tail_threshold()].
tail_t1_clt <- function(
  mu,
  p,
  n_bulk,
  n_ref,
  psi = Inf,
  level = 0.01
) {
  ctx <- .tail_context(mu, p, n_bulk, n_ref, psi, level)
  if (ctx$zero) {
    return(.tail_na("clt"))
  }
  thr <- if (ctx$is_module) ctx$z * sqrt(ctx$n_genes) else ctx$z * ctx$s
  .tail_rows("clt", "T1", thr, thr)
}

#' T1-log, log-normal ratio, one gene only.
#'
#' Upper threshold `mbar * expm1(z * s / mbar)`, lower threshold
#' `-mbar * expm1(-z * s / mbar)`, with `mbar = sum_j p_j mu_j` and
#' `s = sd(D)`. Both stored values are the positive distances from
#' `mbar`; a lower flag is `D <= -threshold`.
tail_t1_log <- function(
  mu,
  p,
  n_bulk,
  n_ref,
  psi = Inf,
  level = 0.01
) {
  if (is.matrix(mu) && nrow(mu) > 1L) {
    stop("Method clt_log applies to one gene only.")
  }
  ctx <- .tail_context(mu, p, n_bulk, n_ref, psi, level)
  if (ctx$zero) {
    return(.tail_na("clt_log"))
  }
  upper <- ctx$mbar * expm1(ctx$z * ctx$s / ctx$mbar)
  lower <- -ctx$mbar * expm1(-ctx$z * ctx$s / ctx$mbar)
  .tail_rows("clt_log", "T1-log", lower, upper)
}

#' T2, dependency-graph normal approximation, module only.
#'
#' `Var S = sum(cov2cor(Cov D))` after within-cell correlations off the
#' edges of `adjacency` are set to zero. Threshold `z * sqrt(Var S)`.
tail_t2_clt_dep <- function(
  mu,
  p,
  n_bulk,
  n_ref,
  psi = Inf,
  level = 0.01,
  adjacency = NULL,
  cell_cor = NULL
) {
  if (!(is.matrix(mu) && nrow(mu) > 1L)) {
    stop("Methods clt_dep apply to modules only.")
  }
  adjacency <- .require_module_graph(mu, adjacency)
  if (is.null(cell_cor)) {
    stop("`cell_cor` is required for clt_dep.")
  }
  ctx <- .tail_context(mu, p, n_bulk, n_ref, psi, level)
  n_types <- ctx$n_types
  n_genes <- ctx$n_genes
  if (is.matrix(cell_cor)) {
    cell_cor <- rep(list(cell_cor), n_types)
  }
  if (length(cell_cor) != n_types) {
    stop("`cell_cor` must be one matrix or a list with one matrix per type.")
  }
  w_j <- ctx$p / ctx$n_bulk + ctx$p^2 / ctx$n_ref
  sigma <- sqrt(ctx$mu + ifelse(is.finite(ctx$psi), ctx$mu^2 / ctx$psi, 0))
  cov_d <- matrix(0, n_genes, n_genes)
  for (j in seq_len(n_types)) {
    r_j <- as.matrix(cell_cor[[j]])
    r_j[!adjacency] <- 0
    r_j[!is.finite(r_j)] <- 0
    diag(r_j) <- 1
    cov_d <- cov_d + w_j[j] * tcrossprod(sigma[, j]) * r_j
  }
  var_s_dep <- sum(stats::cov2cor(cov_d))
  thr <- ctx$z * sqrt(var_s_dep)
  .tail_rows("clt_dep", "T2", thr, thr)
}

#' T3, Poisson Chernoff bound.
#'
#' Gene: cumulant generating function of one gene with `psi` ignored
#' (`family = "poisson"`). Module: sum of standardised gene CGFs,
#' independent genes. The threshold scale is `sd(D)` for a gene and
#' `sqrt(N)` for a module.
tail_t3_poisson <- function(
  mu,
  p,
  n_bulk,
  n_ref,
  psi = Inf,
  level = 0.01
) {
  ctx <- .tail_context(mu, p, n_bulk, n_ref, psi, level)
  if (ctx$zero) {
    return(.tail_na("poisson"))
  }
  if (!ctx$is_module) {
    cgf <- consistency_cgf(
      ctx$mu, ctx$p, ctx$n_bulk, ctx$n_ref, ctx$psi, "poisson"
    )
    return(.tail_chernoff_rows("poisson", "T3", cgf, ctx$s, ctx$level))
  }
  cgf <- .module_cgf_indep(.module_gene_cgfs(ctx, "poisson"), ctx$sd_g)
  .tail_chernoff_rows("poisson", "T3", cgf, sqrt(ctx$n_genes), ctx$level)
}

#' T4, negative-binomial Chernoff bound.
#'
#' Same construction as [tail_t3_poisson()] with `family = "nb"`.
tail_t4_nb <- function(
  mu,
  p,
  n_bulk,
  n_ref,
  psi = Inf,
  level = 0.01
) {
  ctx <- .tail_context(mu, p, n_bulk, n_ref, psi, level)
  if (ctx$zero) {
    return(.tail_na("nb"))
  }
  if (!ctx$is_module) {
    cgf <- consistency_cgf(
      ctx$mu, ctx$p, ctx$n_bulk, ctx$n_ref, ctx$psi, "nb"
    )
    return(.tail_chernoff_rows("nb", "T4", cgf, ctx$s, ctx$level))
  }
  cgf <- .module_cgf_indep(.module_gene_cgfs(ctx, "nb"), ctx$sd_g)
  .tail_chernoff_rows("nb", "T4", cgf, sqrt(ctx$n_genes), ctx$level)
}

#' T5, negative-binomial Chernoff bound after a proper colouring.
#'
#' Module only. The bound is the colour-class mixture of Janson (2004)
#' implemented by [greedy_colouring()].
tail_t5_nb_dep <- function(
  mu,
  p,
  n_bulk,
  n_ref,
  psi = Inf,
  level = 0.01,
  adjacency = NULL
) {
  if (!(is.matrix(mu) && nrow(mu) > 1L)) {
    stop("Methods nb_dep apply to modules only.")
  }
  adjacency <- .require_module_graph(mu, adjacency)
  ctx <- .tail_context(mu, p, n_bulk, n_ref, psi, level)
  cgf <- .module_cgf_colour(
    .module_gene_cgfs(ctx, "nb"),
    ctx$sd_g,
    adjacency
  )
  .tail_chernoff_rows("nb_dep", "T5", cgf, sqrt(ctx$n_genes), ctx$level)
}

#' Thresholds of the tail statements for one gene or one module.
#'
#' Dispatches to [tail_t1_clt()], [tail_t1_log()], [tail_t2_clt_dep()],
#' [tail_t3_poisson()], [tail_t4_nb()] and [tail_t5_nb_dep()].
#' [quantile_band_flags()] is the model-free band and is not a method
#' here.
#'
#' For a gene (`mu` a vector of length `J`, or a one-row matrix) the
#' statistic is `D`; for a module (`mu` an `N` by `J` matrix, `N >= 2`)
#' it is `S = sum_g D_g / sd(D_g)`, with `sd(D_g)` the NB standard
#' deviation. A unit is flagged in the upper tail when the statistic is
#' at least `threshold`, in the lower tail when it is at most
#' `-threshold`.
#'
#' @param mu Per-cell type means: length-`J` vector (gene) or `N` by `J`
#'   matrix (module).
#' @param p Composition (length `J`, summing to one).
#' @param n_bulk Cells in the bulk sample, or the plug-in `a_i`.
#' @param n_ref Reference cells per type (length `J`).
#' @param psi NB size: scalar, length `J` (gene), or `N` by `J` matrix
#'   (module). `Inf` is Poisson.
#' @param method One or several of `"clt"`, `"clt_log"` (gene only),
#'   `"clt_dep"` (module only), `"poisson"`, `"nb"`, `"nb_dep"` (module
#'   only). Default: every statement defined for the unit.
#' @param level Tail probability per side (`alpha / 2`). Default `0.01`.
#' @param adjacency Module only: `N` by `N` 0/1 dependency graph,
#'   required by `clt_dep` and `nb_dep`.
#' @param cell_cor Module only, for `clt_dep`: within-cell correlation of
#'   the genes, one `N` by `N` matrix or a list of `J` (one per type).
#'   Entries off the edges of `adjacency` are set to zero (A6).
#' @return A data frame with `method`, `statement` (`T1`, `T1-log`, `T2`,
#'   `T3`, `T4`, `T5`), `side` (`lower`, `upper`) and `threshold`.
tail_threshold <- function(
  mu,
  p,
  n_bulk,
  n_ref,
  psi = Inf,
  method = c("clt", "clt_log", "clt_dep", "poisson", "nb", "nb_dep"),
  level = 0.01,
  adjacency = NULL,
  cell_cor = NULL
) {
  all_methods <- c("clt", "clt_log", "clt_dep", "poisson", "nb", "nb_dep")
  is_module <- is.matrix(mu) && nrow(mu) > 1L
  gene_methods <- c("clt", "clt_log", "poisson", "nb")
  module_methods <- c("clt", "clt_dep", "poisson", "nb", "nb_dep")
  method <- if (missing(method)) {
    if (is_module) module_methods else gene_methods
  } else {
    match.arg(method, all_methods, several.ok = TRUE)
  }
  if (!is_module) {
    bad <- intersect(method, c("clt_dep", "nb_dep"))
    if (length(bad) > 0L) {
      stop(
        "Methods ", paste(bad, collapse = ", "),
        " apply to modules only."
      )
    }
  } else if ("clt_log" %in% method) {
    stop("Method clt_log applies to one gene only.")
  }
  if (any(c("clt_dep", "nb_dep") %in% method)) {
    adjacency <- .require_module_graph(mu, adjacency)
  }
  if (!is_module) {
    ctx <- .tail_context(mu, p, n_bulk, n_ref, psi, level)
    if (ctx$zero) {
      warning("Zero variance of D: thresholds are NA.", call. = FALSE)
      rows <- lapply(method, function(m) {
        .tail_rows(m, NA_character_, NA_real_, NA_real_)
      })
      return(do.call(rbind, rows))
    }
  }
  call_one <- function(m) {
    switch(
      m,
      clt = tail_t1_clt(mu, p, n_bulk, n_ref, psi, level),
      clt_log = tail_t1_log(mu, p, n_bulk, n_ref, psi, level),
      clt_dep = tail_t2_clt_dep(
        mu, p, n_bulk, n_ref, psi, level, adjacency, cell_cor
      ),
      poisson = tail_t3_poisson(mu, p, n_bulk, n_ref, psi, level),
      nb = tail_t4_nb(mu, p, n_bulk, n_ref, psi, level),
      nb_dep = tail_t5_nb_dep(
        mu, p, n_bulk, n_ref, psi, level, adjacency
      )
    )
  }
  do.call(rbind, lapply(method, call_one))
}


#' Consistency statistic of every gene of one bulk sample.
#'
#' @param bulk Bulk counts of one sample (named vector over genes).
#' @param a_hat Plug-in scale `a_i`.
#' @param p_hat Plug-in composition (length `J`, named by type).
#' @param mu Genes by types per-cell means (gene names as row names).
#' @return Named vector `D_g = y_g / a_i - sum_j p_j mu_gj`.
consistency_statistic <- function(bulk, a_hat, p_hat, mu) {
  if (is.null(names(bulk)) || !all(rownames(mu) %in% names(bulk))) {
    stop("`bulk` must be named and cover every gene of `mu`.")
  }
  if (!is.null(names(p_hat)) && !is.null(colnames(mu))) {
    p_hat <- p_hat[colnames(mu)]
  }
  d <- bulk[rownames(mu)] / a_hat - drop(mu %*% p_hat)
  names(d) <- rownames(mu)
  d
}

#' Quantile-regression band of bulk against fitted bulk, across genes.
#'
#' A model-free reference for the tail statements: `quantreg::rq()` of
#' the observed bulk per cell `y_g / a_i` on the fitted value
#' `sum_j p_j mu_gj`, at `tau = level` and `1 - level`, pooled over the
#' genes of one sample. No count model is assumed; the band treats genes
#' as exchangeable given the fitted value.
#'
#' @param observed,fitted Numeric vectors over genes (same names).
#' @param tau Lower and upper quantile levels. Default `c(0.01, 0.99)`.
#' @return A data frame with `gene`, `observed`, `fitted`, `q_lower`,
#'   `q_upper` and `flag` (outside the band).
quantile_band_flags <- function(observed, fitted, tau = c(0.01, 0.99)) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop(
      "Package quantreg is required for the quantile band (not in ",
      "renv.lock). Install it with install.packages(\"quantreg\")."
    )
  }
  if (length(observed) != length(fitted)) {
    stop("`observed` and `fitted` must have the same length.")
  }
  if (length(tau) != 2L || tau[1L] >= tau[2L]) {
    stop("`tau` must be two increasing quantile levels.")
  }
  d <- data.frame(observed = observed, fitted = fitted)
  fit <- quantreg::rq(observed ~ fitted, tau = tau, data = d)
  q <- stats::predict(fit, newdata = d)
  q <- matrix(q, ncol = 2L)
  data.frame(
    gene = names(observed) %||% as.character(seq_along(observed)),
    observed = observed,
    fitted = fitted,
    q_lower = q[, 1L],
    q_upper = q[, 2L],
    flag = observed < q[, 1L] | observed > q[, 2L],
    stringsAsFactors = FALSE
  )
}

#' Bulk consistency filter by tail statements, for one matched pair.
#'
#' For one time slot: reference moments ([type_moments_nb()]), plug-in
#' composition and scale per bulk replicate
#' ([plugin_composition_demixsc()]), the statistic `D_gi`, the gene
#' statements T1, T1-log, T3, T4 and, per module, T1 to T5, plus the
#' quantile-regression band. Each statement flags at `level` per tail.
#'
#' @param bulk Genes by samples raw bulk counts (the replicates of the
#'   matched condition), gene names as row names.
#' @param counts Genes by cells raw reference counts of the same slot.
#' @param cell_type Cell-type label per column of `counts`.
#' @param modules Optional named list of gene vectors (for example the
#'   lineage groups of the signalling dictionary). Modules with fewer
#'   than two genes in the panel are skipped.
#' @param level Tail probability per side. Default `0.01`.
#' @param c_reg Constant of the deMixSC weights. Default `2`.
#' @param max_degree,mean_degree,min_abs_cor Passed to
#'   [module_dependency_graph()].
#' @param qr_tau Quantile levels of the band, or `NULL` to skip it.
#'   Default `c(level, 1 - level)`.
#' @return A list with `gene` (one row per sample, gene and statement:
#'   `sample`, `gene`, `statistic`, `statement`, `t_lower`, `t_upper`,
#'   `flag`), `module` (same columns with `module` instead of `gene`),
#'   `quantile_band`, `plugin` (output of
#'   [plugin_composition_demixsc()]), `moments` and `graphs` (adjacency
#'   per module).
bulk_tail_filter <- function(
  bulk,
  counts,
  cell_type,
  modules = NULL,
  level = 0.01,
  c_reg = 2,
  max_degree = 10L,
  mean_degree = 3,
  min_abs_cor = 0.1,
  qr_tau = c(level, 1 - level)
) {
  bulk <- as.matrix(bulk)
  if (is.null(rownames(bulk)) || is.null(rownames(counts))) {
    stop("`bulk` and `counts` need gene names as row names.")
  }
  if (is.null(colnames(bulk))) {
    colnames(bulk) <- paste0("sample_", seq_len(ncol(bulk)))
  }
  genes <- intersect(rownames(counts), rownames(bulk))
  if (length(genes) < 10L) {
    stop("Fewer than ten genes shared by `bulk` and `counts`.")
  }
  cell_type <- droplevels(factor(cell_type))
  counts <- counts[genes, , drop = FALSE]
  mom <- type_moments_nb(counts, cell_type)
  plug <- plugin_composition_demixsc(bulk[genes, , drop = FALSE], mom$mu, c_reg)
  types <- colnames(mom$mu)

  long_rows <- function(unit_col, unit, sample, statistic, thr) {
    wide_lo <- thr$threshold[thr$side == "lower"]
    wide_hi <- thr$threshold[thr$side == "upper"]
    stat_names <- thr$statement[thr$side == "lower"]
    out <- data.frame(
      sample = sample,
      unit = unit,
      statistic = statistic,
      statement = stat_names,
      t_lower = wide_lo,
      t_upper = wide_hi,
      stringsAsFactors = FALSE
    )
    out$flag <- out$statistic >= out$t_upper | out$statistic <= -out$t_lower
    names(out)[names(out) == "unit"] <- unit_col
    out
  }

  gene_rows <- list()
  module_rows <- list()
  band_rows <- list()
  graphs <- list()
  cor_by_module <- list()
  if (!is.null(modules)) {
    for (mod in names(modules)) {
      g_mod <- intersect(modules[[mod]], genes)
      if (length(g_mod) < 2L) {
        next
      }
      cor_by_module[[mod]] <- lapply(types, function(tp) {
        x <- as.matrix(counts[g_mod, cell_type == tp, drop = FALSE])
        r <- suppressWarnings(stats::cor(t(x)))
        r[!is.finite(r)] <- 0
        r
      })
      graphs[[mod]] <- module_dependency_graph(
        cor_by_module[[mod]],
        max_degree = max_degree,
        mean_degree = mean_degree,
        min_abs_cor = min_abs_cor
      )
    }
  }

  for (i in seq_len(ncol(bulk))) {
    sample <- colnames(bulk)[i]
    comp <- plug$composition[plug$composition$sample == sample, ]
    p_hat <- stats::setNames(comp$p_hat, comp$cell_type)[types]
    a_hat <- plug$scale$a_hat[plug$scale$sample == sample]
    y <- stats::setNames(bulk[genes, i], genes)
    d <- consistency_statistic(y, a_hat, p_hat, mom$mu)
    for (g in genes) {
      thr <- tail_threshold(
        mom$mu[g, ],
        p_hat,
        n_bulk = a_hat,
        n_ref = mom$n_cells,
        psi = mom$psi[g, ],
        level = level
      )
      gene_rows[[length(gene_rows) + 1L]] <- long_rows(
        "gene",
        g,
        sample,
        d[[g]],
        thr
      )
    }
    if (!is.null(qr_tau)) {
      band <- quantile_band_flags(
        y / a_hat,
        drop(mom$mu %*% p_hat),
        tau = qr_tau
      )
      band$sample <- sample
      band_rows[[i]] <- band
    }
    for (mod in names(graphs)) {
      g_mod <- rownames(graphs[[mod]])
      sd_g <- vapply(
        g_mod,
        function(g) {
          sqrt(consistency_variance(
            mom$mu[g, ],
            p_hat,
            a_hat,
            mom$n_cells,
            mom$psi[g, ]
          ))
        },
        numeric(1)
      )
      ok <- is.finite(sd_g) & sd_g > 0
      if (sum(ok) < 2L) {
        next
      }
      g_mod <- g_mod[ok]
      s_stat <- sum(d[g_mod] / sd_g[ok])
      thr <- tail_threshold(
        mom$mu[g_mod, , drop = FALSE],
        p_hat,
        n_bulk = a_hat,
        n_ref = mom$n_cells,
        psi = mom$psi[g_mod, , drop = FALSE],
        level = level,
        adjacency = graphs[[mod]][g_mod, g_mod, drop = FALSE],
        cell_cor = lapply(cor_by_module[[mod]], function(r) {
          r[g_mod, g_mod, drop = FALSE]
        })
      )
      module_rows[[length(module_rows) + 1L]] <- long_rows(
        "module",
        mod,
        sample,
        s_stat,
        thr
      )
    }
  }
  list(
    gene = do.call(rbind, gene_rows),
    module = if (length(module_rows) > 0L) {
      do.call(rbind, module_rows)
    } else {
      NULL
    },
    quantile_band = if (length(band_rows) > 0L) {
      do.call(rbind, band_rows)
    } else {
      NULL
    },
    plugin = plug,
    moments = mom,
    graphs = graphs
  )
}

#' Overlap and order of the flags of the tail statements.
#'
#' @param res_gene The `gene` table of [bulk_tail_filter()] (or the
#'   `module` table, with `unit = "module"`).
#' @param band Optional `quantile_band` table, added as statement `QR`.
#' @param unit Column naming the flagged unit. Default `"gene"`.
#' @return A list with `overlap` (output of [set_overlap()] on the units
#'   flagged in at least one sample, by statement) and `order` (one row
#'   per statement: `n_flagged`, `median_t_lower`, `median_t_upper`;
#'   sorted from the fewest flags, the most conservative, to the most).
tail_flag_overlap <- function(res_gene, band = NULL, unit = "gene") {
  flagged <- split(
    res_gene[[unit]][res_gene$flag %in% TRUE],
    res_gene$statement[res_gene$flag %in% TRUE]
  )
  statements <- unique(res_gene$statement)
  sets <- lapply(statements, function(s) unique(flagged[[s]] %||% character()))
  names(sets) <- statements
  if (!is.null(band)) {
    sets[["QR"]] <- unique(band$gene[band$flag %in% TRUE])
  }
  order_tab <- data.frame(
    statement = names(sets),
    n_flagged = lengths(sets),
    median_t_lower = vapply(
      names(sets),
      function(s) {
        v <- res_gene$t_lower[res_gene$statement == s]
        if (length(v) == 0L) NA_real_ else stats::median(v, na.rm = TRUE)
      },
      numeric(1)
    ),
    median_t_upper = vapply(
      names(sets),
      function(s) {
        v <- res_gene$t_upper[res_gene$statement == s]
        if (length(v) == 0L) NA_real_ else stats::median(v, na.rm = TRUE)
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  order_tab <- order_tab[order(order_tab$n_flagged, order_tab$statement), ]
  rownames(order_tab) <- NULL
  list(overlap = set_overlap(sets), order = order_tab)
}
