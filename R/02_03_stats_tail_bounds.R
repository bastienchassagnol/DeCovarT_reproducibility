#' Discrepancy and drift between the single-cell reference and the bulk.
#'
#' Tier (iii) of the gene selection, "Bulk consistency filter". One file
#' holds the whole cross-platform check; Cook's distances are in
#' `R/02_03_stats_influence.R`.
#'
#' 1. `p`-free range check ([bulk_consistency_filter()], RNA-Sieve range
#'    bound and MAD rules; a heuristic).
#' 2. Plug-in composition and scale per bulk replicate
#'    ([plugin_composition_demixsc()], deMixSC iterated weighted NNLS).
#' 3. Tail statements T1 to T5 per gene and per module, and the
#'    quantile-regression band ([bulk_tail_filter()]).
#' 4. Drift across the matched times ([platform_drift_summary()]).
#'
#' Why deMixSC and not plain NNLS for the plug-in. The variance of a
#' bulk count under the DeCovarT model is gene-specific and depends on
#' the composition and on the type means. Unweighted NNLS gives every
#' gene the same weight, so the most expressed genes decide `p`. The
#' deMixSC weights `1 / (yhat^2 + (y - yhat)^2 + c)` are inverse-variance
#' weights that follow the current fit, the diagonal plug-in analogue of
#' the generalised least squares of DeCovarT; their residual term
#' down-weights the discrepant genes, which are the genes this filter
#' must find, so that they do not bias `p`. NNLS is only the solver of
#' each reweighting step (`nnls::nnls()`, Lawson-Hanson).
#'
#' For gene `g`, bulk sample `i` and reference types `j = 1, ..., J`,
#' the statistic is
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
#' The tail statements are a port of
#' `scripts/figures/fig_03_1_tail_bounds.R` of the DeCovarT repository,
#' generalised to type-specific `psi_gj` and to a Poisson limit
#' `psi = Inf`; deliberately duplicated, no cross-repository dependency.
#' Requires `R/02_00_utils.R`. Uses `nnls` (plug-in) and `quantreg`
#' (quantile band).

# --------------------------------------------------------------------------
# Step 1 · p-free range check (RNA-Sieve heuristic) ----
# --------------------------------------------------------------------------

#' Flag genes whose bulk expression lies outside the reference range.
#'
#' **Heuristic.** This adapts the two gene filters of RNA-Sieve
#' (Erdmann-Pham et al. 2021, Genome Research, eq. 12 and 13) to a
#' DeCovarT panel. Neither rule is derived from DeCovarT's likelihood:
#' the range rule follows from the mixture constraint, the MAD rules are
#' empirical (the RNA-Sieve authors say their choice was "primarily
#' guided by empirical considerations"). It needs no composition, so it
#' does not share the plug-in's dependence on the genes it tests. A
#' flagged gene is a candidate for removal, to be reported, not a proven
#' platform artefact.
#'
#' Scale. Each type profile and each bulk sample is normalised to sum
#' one over the panel genes. If the bulk is a mixture of the types,
#' `b = sum_k alpha_k m_k`, the normalised bulk is a convex combination
#' of the normalised type profiles (weights proportional to
#' `alpha_k * sum_g m_gk`). So, gene by gene, the normalised bulk lies in
#' `[min_k m_gk, max_k m_gk]` up to sampling error, whatever the
#' proportions and whatever a global platform scale.
#'
#' Rule 1, range bound (RNA-Sieve eq. 12): the distance of the
#' normalised bulk `b_g` to the interval `[min_k m_gk, max_k m_gk]`,
#' divided by `se_g = sqrt(max_k s_gk^2 / C_k + v_g)`, where `s_gk^2 /
#' C_k` is the sampling variance of a type mean and `v_g` the variance
#' of `b_g` over bulk replicates divided by their number (`0` for one
#' sample). The gene is flagged when that ratio exceeds `z_range`. The
#' standard error is a working choice; RNA-Sieve bounds the tail
#' probability instead.
#'
#' Rule 2, MAD outliers (RNA-Sieve eq. 13), on four per-gene summaries:
#' `T_M = max_k m_gk`, `T_S = max_k s_gk`, `R_SM = max_k s_gk / m_gk`
#' (largest coefficient of variation) and `R_bM = b_g / max_k m_gk`.
#' A gene is flagged when any summary is more than `mad_k` scaled MADs
#' (`stats::mad()`, constant 1.4826) from that summary's median over the
#' panel. The summaries are compared on the log scale when
#' `log_scale = TRUE` (default), because they are positive and
#' right-skewed; RNA-Sieve does not state the scale.
#'
#' @param bulk Bulk expression of the panel genes: a vector, or a genes by
#'   samples matrix of replicates at one matched time (48, 72 or 96 h).
#'   Row names (or names) must be gene names.
#' @param mu Genes by types matrix of reference means on the linear
#'   scale (size-factor-normalised counts).
#' @param sd Genes by types matrix of reference within-type standard
#'   deviations, same scale as `mu`.
#' @param n_cells Number of cells per type (named by the columns of
#'   `mu`, or in column order).
#' @param z_range Threshold of the range rule. Default `3`.
#' @param mad_k Threshold of the MAD rule. Default `3`.
#' @param log_scale Compare the MAD summaries on the log scale. Default
#'   `TRUE`.
#' @return A data frame with one row per gene: `gene`, `bulk`,
#'   `ref_min`, `ref_max`, `range_distance`, `range_se`, `range_z`,
#'   `flag_range`, `T_M`, `T_S`, `R_SM`, `R_bM`, one `flag_<summary>`
#'   column each, `flag_mad`, and `keep` (no flag). The attribute
#'   `heuristic` repeats that the filter is heuristic.
bulk_consistency_filter <- function(
  bulk,
  mu,
  sd,
  n_cells,
  z_range = 3,
  mad_k = 3,
  log_scale = TRUE
) {
  mu <- as.matrix(mu)
  sd <- as.matrix(sd)
  if (!identical(dim(mu), dim(sd))) {
    stop("`mu` and `sd` must have the same dimensions.")
  }
  genes <- rownames(mu)
  if (is.null(genes)) {
    stop("`mu` needs gene names as row names.")
  }
  bulk <- if (is.null(dim(bulk))) {
    matrix(bulk, ncol = 1L, dimnames = list(names(bulk), NULL))
  } else {
    as.matrix(bulk)
  }
  if (is.null(rownames(bulk))) {
    stop("`bulk` needs gene names (names or row names).")
  }
  missing_genes <- setdiff(genes, rownames(bulk))
  if (length(missing_genes) > 0L) {
    stop(
      length(missing_genes),
      " reference genes are missing from `bulk`, for example ",
      paste(utils::head(missing_genes, 5L), collapse = ", "),
      "."
    )
  }
  bulk <- bulk[genes, , drop = FALSE]
  if (!is.null(names(n_cells)) && !is.null(colnames(mu))) {
    n_cells <- n_cells[colnames(mu)]
  }
  n_cells <- as.numeric(n_cells)
  if (length(n_cells) != ncol(mu) || any(!is.finite(n_cells))) {
    stop("`n_cells` must give one cell count per column of `mu`.")
  }

  # Sum-to-one normalisation over the panel, per type and per sample.
  type_total <- colSums(mu)
  m_norm <- sweep(mu, 2L, type_total, "/")
  s_norm <- sweep(sd, 2L, type_total, "/")
  b_norm <- sweep(bulk, 2L, colSums(bulk), "/")
  b <- rowMeans(b_norm)
  v_b <- if (ncol(b_norm) > 1L) {
    apply(b_norm, 1L, stats::var) / ncol(b_norm)
  } else {
    rep(0, length(b))
  }

  # Rule 1 · range bound.
  ref_min <- apply(m_norm, 1L, min)
  ref_max <- apply(m_norm, 1L, max)
  range_distance <- pmax(ref_min - b, b - ref_max, 0)
  se_mean <- sweep(s_norm^2, 2L, n_cells, "/")
  range_se <- sqrt(apply(se_mean, 1L, max) + v_b)
  range_z <- ifelse(range_se > 0, range_distance / range_se, Inf)
  range_z[range_distance == 0] <- 0
  flag_range <- range_z > z_range

  # Rule 2 · MAD outliers of four summaries.
  summaries <- data.frame(
    T_M = apply(m_norm, 1L, max),
    T_S = apply(s_norm, 1L, max),
    R_SM = apply(s_norm / m_norm, 1L, max),
    R_bM = b / apply(m_norm, 1L, max)
  )
  mad_flag <- function(v) {
    w <- if (log_scale) log(v) else v
    ok <- is.finite(w)
    if (sum(ok) < 3L) {
      return(rep(FALSE, length(v)))
    }
    centre <- stats::median(w[ok])
    spread <- stats::mad(w[ok])
    if (!is.finite(spread) || spread <= 0) {
      return(rep(FALSE, length(v)))
    }
    flag <- abs(w - centre) / spread > mad_k
    # A non-finite summary (zero bulk or zero mean) is an outlier.
    flag[!ok] <- TRUE
    flag
  }
  flags <- vapply(summaries, mad_flag, logical(nrow(summaries)))
  flags <- matrix(flags, nrow = nrow(summaries))
  colnames(flags) <- paste0("flag_", names(summaries))
  flag_mad <- rowSums(flags) > 0L

  out <- data.frame(
    gene = genes,
    bulk = b,
    ref_min = ref_min,
    ref_max = ref_max,
    range_distance = range_distance,
    range_se = range_se,
    range_z = range_z,
    flag_range = flag_range,
    summaries,
    flags,
    flag_mad = flag_mad,
    keep = !(flag_range | flag_mad),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "heuristic") <- paste(
    "RNA-Sieve range bound and MAD outlier flags, fitted per matched time;",
    "not derived from the DeCovarT likelihood."
  )
  out
}

# --------------------------------------------------------------------------
# Step 2 · Plug-in composition (deMixSC) ----
# --------------------------------------------------------------------------

#' deMixSC iterated weighted NNLS composition of bulk samples.
#'
#' For each bulk sample: start from the unweighted NNLS fit, then
#' iterate `w_g = 1 / (yhat_g^2 + (y_g - yhat_g)^2 + c_reg)`,
#' `beta = NNLS(sqrt(w) X, sqrt(w) y)` until the largest change of
#' `beta`, relative to `sum(beta)`, is below `tol` (deMixSC, Guo et al.,
#' Eq. 4). deMixSC also separates a gene set G2 calibrated on a
#' same-aliquot benchmark; no such benchmark exists for Suppinger, so G2
#' is empty and every gene enters the fit. `a_i` is the sum of the
#' coefficients and `p_i` the coefficients divided by `a_i`; the extra
#' overdispersion between platforms is fixed at 1. MuSiC is not used,
#' because its weights need at least two donors.
#'
#' `c_reg` has the units of a squared expression value, so it is not
#' scale-free; bulk and signature must be on the same linear scale.
#' Every tail statement built on this plug-in holds only with these
#' parameters treated as known.
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
  if (!requireNamespace("nnls", quietly = TRUE)) {
    stop("Package nnls is required for the deMixSC plug-in.")
  }
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
  solve_nnls <- function(a, b) nnls::nnls(a, b)$x

  types <- colnames(signature)
  weights <- matrix(NA_real_, nrow(bulk), ncol(bulk), dimnames = dimnames(bulk))
  fitted <- weights
  comp <- list()
  scale <- list()
  for (i in seq_len(ncol(bulk))) {
    y <- bulk[, i]
    beta <- solve_nnls(signature, y)
    converged <- FALSE
    it <- 0L
    for (it in seq_len(max_iter)) {
      yhat <- drop(signature %*% beta)
      w <- 1 / (yhat^2 + (y - yhat)^2 + c_reg)
      beta_new <- solve_nnls(signature * sqrt(w), y * sqrt(w))
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

# --------------------------------------------------------------------------
# Step 3 · Tail statements ----
# --------------------------------------------------------------------------

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
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 0.5) {
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
      ctx$mu,
      ctx$p,
      ctx$n_bulk,
      ctx$n_ref,
      ctx$psi,
      "poisson"
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
      ctx$mu,
      ctx$p,
      ctx$n_bulk,
      ctx$n_ref,
      ctx$psi,
      "nb"
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
        "Methods ",
        paste(bad, collapse = ", "),
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
        mu,
        p,
        n_bulk,
        n_ref,
        psi,
        level,
        adjacency,
        cell_cor
      ),
      poisson = tail_t3_poisson(mu, p, n_bulk, n_ref, psi, level),
      nb = tail_t4_nb(mu, p, n_bulk, n_ref, psi, level),
      nb_dep = tail_t5_nb_dep(
        mu,
        p,
        n_bulk,
        n_ref,
        psi,
        level,
        adjacency
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
    stop("Package quantreg is required for the quantile band.")
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
#' @param platform_ratio Optional named vector of plug-in gene-wise
#'   platform ratios `b_g` (bulk over reference, natural scale; see
#'   [cross_fit_platform_ratio()]). The bulk is divided by `b_g` before
#'   the plug-in fit, so the statistic becomes
#'   `D^b_gi = y_gi / (a_i b_g) - sum_j p_ji mu_gj`. Genes absent from
#'   the vector get `b_g = 1`. `NULL` (default) is the A3 model.
#' @param platform_var Optional named vector of `Var(log b_g)`, the
#'   sampling variance of the plug-in ratio. By the delta method it adds
#'   `(sum_j p_ji mu_gj)^2 Var(log b_g)` to `Var D^b`; that term enters
#'   `z` and the T1 threshold. T1-log, T3, T4 and the module statements
#'   treat `b_g` as known.
#' @return A list with `gene` (one row per sample, gene and statement:
#'   `sample`, `gene`, `statistic`, `statement`, `t_lower`, `t_upper`,
#'   `flag`, and `z = statistic / sd(D)`, the NB-standardised
#'   discrepancy used by [platform_drift_summary()]), `module` (same
#'   columns with `module` instead of `gene`, without `z`),
#'   `quantile_band`, `plugin` (output of
#'   [plugin_composition_demixsc()]), `moments`, `graphs` (adjacency per
#'   module), `platform` (`gene`, `b`, `var_log_b`; all ones and zeros
#'   without `platform_ratio`) and `settings` (`level`, `c_reg`,
#'   `max_degree`, `mean_degree`, `min_abs_cor`, `platform_ratio`).
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
  qr_tau = c(level, 1 - level),
  platform_ratio = NULL,
  platform_var = NULL
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
  b_g <- stats::setNames(rep(1, length(genes)), genes)
  v_log_b <- stats::setNames(rep(0, length(genes)), genes)
  if (!is.null(platform_ratio)) {
    if (is.null(names(platform_ratio)) || any(platform_ratio <= 0, na.rm = TRUE)) {
      stop("`platform_ratio` must be a named vector of positive ratios.")
    }
    known <- intersect(genes, names(platform_ratio))
    b_g[known] <- platform_ratio[known]
    b_g[!is.finite(b_g)] <- 1
  }
  if (!is.null(platform_var)) {
    if (is.null(names(platform_var))) {
      stop("`platform_var` must be a named vector.")
    }
    known <- intersect(genes, names(platform_var))
    v_log_b[known] <- platform_var[known]
    v_log_b[!is.finite(v_log_b) | v_log_b < 0] <- 0
  }
  bulk <- sweep(bulk[genes, , drop = FALSE], 1L, b_g, "/")
  plug <- plugin_composition_demixsc(bulk, mom$mu, c_reg)
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
      sd_g <- sqrt(consistency_variance(
        mom$mu[g, ],
        p_hat,
        a_hat,
        mom$n_cells,
        mom$psi[g, ]
      ))
      # Delta method for the plug-in ratio: D^b moves by -fit * d log b.
      sd_tot <- sqrt(sd_g^2 + sum(mom$mu[g, ] * p_hat)^2 * v_log_b[[g]])
      if (v_log_b[[g]] > 0 && is.finite(sd_g) && sd_g > 0) {
        is_t1 <- thr$statement %in% "T1"
        thr$threshold[is_t1] <- thr$threshold[is_t1] * sd_tot / sd_g
      }
      rows_g <- long_rows("gene", g, sample, d[[g]], thr)
      rows_g$z <- if (is.finite(sd_tot) && sd_tot > 0) {
        d[[g]] / sd_tot
      } else {
        NA_real_
      }
      gene_rows[[length(gene_rows) + 1L]] <- rows_g
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
    graphs = graphs,
    platform = data.frame(
      gene = genes,
      b = as.numeric(b_g),
      var_log_b = as.numeric(v_log_b),
      stringsAsFactors = FALSE
    ),
    settings = list(
      level = level,
      c_reg = c_reg,
      max_degree = max_degree,
      mean_degree = mean_degree,
      min_abs_cor = min_abs_cor,
      platform_ratio = !is.null(platform_ratio)
    )
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

# --------------------------------------------------------------------------
# Step 4 · Drift across the matched times ----
# --------------------------------------------------------------------------

#' Drift of the platform discrepancy across the matched times.
#'
#' The cross-platform filter assumes that a platform effect is a
#' property of the gene, constant over time. With three matched times
#' (48, 72 and 96 h) that assumption can be read from the flags: a gene
#' flagged at several times with the same sign of the standardised
#' discrepancy `z = D / sd(D)` drifts the same way on both platforms
#' whatever the composition; a gene flagged at one time only does not
#' support the transfer of its flag to the other times.
#'
#' @param gene_tables Named list of the `gene` tables of
#'   [bulk_tail_filter()], one per time (names are the times).
#' @param statement Gene statement whose flags are summarised: `"T4"`
#'   (default, the NB Chernoff bound), `"T1"`, `"T1-log"` or `"T3"`.
#' @return A list with
#'   - `by_time`: one row per gene and time: `gene`, `time_point`,
#'     `n_samples`, `n_flagged` (replicates flagged), `mean_z`;
#'   - `summary`: one row per gene: `gene`, `n_times` (times where the
#'     gene was tested), `n_times_flagged`, `n_flags` (time by replicate
#'     flags), `flagged_every_time`, `drift_class` (`"not flagged"`,
#'     `"one time"`, `"several times, same sign"`, `"several times, sign
#'     change"`) and one `mean_z_<time>` column per time.
platform_drift_summary <- function(gene_tables, statement = "T4") {
  times <- names(gene_tables)
  if (is.null(times) || any(!nzchar(times)) || anyDuplicated(times)) {
    stop("`gene_tables` must be a list named by distinct times.")
  }
  if (length(statement) != 1L) {
    stop("`statement` must be one statement label, for example \"T4\".")
  }
  by_time <- do.call(
    rbind,
    lapply(times, function(tm) {
      tab <- gene_tables[[tm]]
      tab <- tab[tab$statement %in% statement, , drop = FALSE]
      if (nrow(tab) == 0L) {
        stop("No row of statement ", statement, " at ", tm, ".")
      }
      by_gene <- split(tab, tab$gene)
      data.frame(
        gene = names(by_gene),
        time_point = tm,
        n_samples = vapply(by_gene, nrow, integer(1)),
        n_flagged = vapply(
          by_gene,
          function(r) sum(r$flag %in% TRUE),
          integer(1)
        ),
        mean_z = vapply(
          by_gene,
          function(r) {
            if (all(is.na(r$z))) NA_real_ else mean(r$z, na.rm = TRUE)
          },
          numeric(1)
        ),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(by_time) <- NULL
  summary <- do.call(
    rbind,
    lapply(split(by_time, by_time$gene), function(r) {
      flagged <- r$n_flagged > 0L
      signs <- sign(r$mean_z[flagged])
      signs <- unique(signs[!is.na(signs) & signs != 0])
      n_flag_t <- sum(flagged)
      drift_class <- if (n_flag_t == 0L) {
        "not flagged"
      } else if (n_flag_t == 1L) {
        "one time"
      } else if (length(signs) <= 1L) {
        "several times, same sign"
      } else {
        "several times, sign change"
      }
      data.frame(
        gene = r$gene[1L],
        n_times = nrow(r),
        n_times_flagged = n_flag_t,
        n_flags = sum(r$n_flagged),
        flagged_every_time = nrow(r) > 1L && n_flag_t == nrow(r),
        drift_class = drift_class,
        stringsAsFactors = FALSE
      )
    })
  )
  key <- paste(by_time$gene, by_time$time_point, sep = "\r")
  for (tm in times) {
    summary[[paste0("mean_z_", tm)]] <- by_time$mean_z[
      match(paste(summary$gene, tm, sep = "\r"), key)
    ]
  }
  summary <- summary[order(-summary$n_flags, summary$gene), , drop = FALSE]
  rownames(summary) <- NULL
  list(by_time = by_time, summary = summary)
}

#' Gene-wise bulk-to-reference ratio per matched time.
#'
#' For each gene and time, the mean over the bulk replicates of
#' `log2((y / a + pseudo) / (sum_j p_j mu_gj + pseudo))`, read from the
#' quantile-band table of [bulk_tail_filter()] (observed and fitted bulk
#' per cell). Under assumption A3 (same per-cell mean on both platforms)
#' the ratio is 0 up to sampling noise. A ratio that is far from 0 and
#' strongly correlated across times is a stable, gene-wise platform
#' effect, which the tail statements read as a discrepancy at every time.
#'
#' @param tails Named list (times) of [bulk_tail_filter()] outputs.
#' @param pseudo Pseudo-count on the per-cell scale. Default `0.01`.
#' @return A list with `log_ratio` (genes by times matrix), `cor`
#'   (Spearman correlation between times, on the genes present at every
#'   time) and `summary` (per time: genes, median, IQR and SD of the
#'   log2 ratio, share of genes with |log2 ratio| > 1).
platform_log_ratio <- function(tails, pseudo = 0.01) {
  times <- names(tails)
  per_time <- lapply(times, function(tm) {
    b <- tails[[tm]]$quantile_band
    if (is.null(b)) {
      stop("No quantile_band table at ", tm, "; rerun with qr_tau set.")
    }
    lr <- log2((b$observed + pseudo) / (b$fitted + pseudo))
    tapply(lr, b$gene, mean)
  })
  names(per_time) <- times
  genes <- Reduce(union, lapply(per_time, names))
  mat <- vapply(
    per_time,
    function(v) as.numeric(v[genes]),
    numeric(length(genes))
  )
  mat <- matrix(mat, nrow = length(genes), dimnames = list(genes, times))
  complete <- stats::complete.cases(mat)
  summary <- data.frame(
    time_point = times,
    genes = colSums(!is.na(mat)),
    median_log2_ratio = apply(mat, 2L, stats::median, na.rm = TRUE),
    iqr_log2_ratio = apply(mat, 2L, stats::IQR, na.rm = TRUE),
    sd_log2_ratio = apply(mat, 2L, stats::sd, na.rm = TRUE),
    share_abs_gt_1 = colMeans(abs(mat) > 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  rownames(summary) <- NULL
  list(
    log_ratio = mat,
    cor = stats::cor(mat[complete, , drop = FALSE], method = "spearman"),
    summary = summary
  )
}

# --------------------------------------------------------------------------
# Step 5 · Gene-wise platform ratio ----
# --------------------------------------------------------------------------
#
# Under A3 the bulk per cell and the reference mean of a gene agree up to
# sampling noise. A gene-wise platform ratio b_g (the b-tilde of the
# DeCovarT vignettes) relaxes A3 to
#
#   E[y_gi / a_i] = b_g sum_j p_ji mu_gj.
#
# Identifiability. b_g is diagonal (one scale per gene). A common factor
# of all b_g is absorbed by a_i, so the ratios are fixed by
# median_g log b_g = 0. Estimated on the sample under test, b_g absorbs
# any residual (MEAD, Theorem 1); it is therefore estimated on other
# samples whose compositions differ (cross-fitting over the matched
# times). A gene-by-type ratio b_gj is confounded with mu_gj and is not
# identified without paired samples.

#' Mean transcript length per gene, averaged over isoforms.
#'
#' Exonic length of every UCSC knownGene transcript of mm10
#' (`GenomicFeatures::transcriptLengths()` on
#' `TxDb.Mmusculus.UCSC.mm10.knownGene`), averaged over the isoforms of
#' each Entrez gene. Genes are matched by Ensembl id (version stripped)
#' through `org.Mm.eg.db`, then by symbol when the Ensembl id has no
#' Entrez gene with a transcript. Adapted from
#' `GastroDeconv2FateMap::annotate_genes_local()`, which uses the genomic
#' span (`genes()` width, introns included); a full-length library
#' samples the mature transcript, so the exonic length is used here.
#'
#' @param symbols Gene symbols (row names of the panel).
#' @param ensembl_ids Ensembl gene ids in the same order, with or without
#'   version suffix; `NA` allowed.
#' @return A data frame with `gene`, `ensembl_id`, `entrez_id`, `route`
#'   (`"ensembl"`, `"symbol"` or `NA`), `n_isoforms` and
#'   `mean_tx_length` (bases).
gene_transcript_length <- function(symbols, ensembl_ids = NULL) {
  pkgs <- c(
    "AnnotationDbi",
    "GenomicFeatures",
    "org.Mm.eg.db",
    "TxDb.Mmusculus.UCSC.mm10.knownGene"
  )
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0L) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), ".")
  }
  symbols <- as.character(symbols)
  if (is.null(ensembl_ids)) {
    ensembl_ids <- rep(NA_character_, length(symbols))
  }
  if (length(ensembl_ids) != length(symbols)) {
    stop("`ensembl_ids` and `symbols` must have the same length.")
  }
  ens <- sub("\\..*$", "", as.character(ensembl_ids))
  tx <- GenomicFeatures::transcriptLengths(
    TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  )
  tx <- tx[!is.na(tx$gene_id), , drop = FALSE]
  len <- tapply(tx$tx_len, tx$gene_id, mean)
  n_iso <- tapply(tx$tx_len, tx$gene_id, length)
  map_one <- function(keys, keytype) {
    keys <- unique(keys[!is.na(keys) & nzchar(keys)])
    if (length(keys) == 0L) {
      return(character())
    }
    out <- suppressMessages(AnnotationDbi::mapIds(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = keys,
      keytype = keytype,
      column = "ENTREZID",
      multiVals = "first"
    ))
    unlist(out)
  }
  by_ens <- map_one(ens, "ENSEMBL")
  by_sym <- map_one(symbols, "SYMBOL")
  entrez <- unname(by_ens[ens])
  route <- ifelse(!is.na(entrez) & entrez %in% names(len), "ensembl", NA)
  from_sym <- unname(by_sym[symbols])
  use_sym <- is.na(route) & !is.na(from_sym) & from_sym %in% names(len)
  entrez[use_sym] <- from_sym[use_sym]
  route[use_sym] <- "symbol"
  entrez[is.na(route)] <- NA_character_
  data.frame(
    gene = symbols,
    ensembl_id = ens,
    entrez_id = entrez,
    route = route,
    n_isoforms = as.integer(n_iso[entrez]),
    mean_tx_length = as.numeric(len[entrez]),
    stringsAsFactors = FALSE
  )
}

#' Log ratio of bulk to fitted reference, gene by replicate.
#'
#' Fits the deMixSC plug-in on `bulk / b` and returns, for each gene and
#' replicate, `log((y / a + pseudo) / (sum_j p_j mu_gj + pseudo))` with
#' the original counts `y`. At the fixed point of
#' [estimate_platform_ratio()] its median over replicates is `log b`.
#'
#' @param bulk Genes by replicates raw bulk counts (row names: genes).
#' @param mu Genes by types per-cell reference means (raw counts).
#' @param log_b Named vector of current `log b_g` (natural log); genes
#'   absent get 0.
#' @param c_reg,pseudo deMixSC constant and per-cell pseudo-count.
#' @return Genes by replicates matrix of natural-log ratios.
platform_ratio_residuals <- function(
  bulk,
  mu,
  log_b = NULL,
  c_reg = 2,
  pseudo = 0.01
) {
  genes <- rownames(mu)
  bulk <- as.matrix(bulk)[genes, , drop = FALSE]
  lb <- stats::setNames(rep(0, length(genes)), genes)
  if (!is.null(log_b)) {
    known <- intersect(genes, names(log_b))
    lb[known] <- log_b[known]
    lb[!is.finite(lb)] <- 0
  }
  plug <- plugin_composition_demixsc(sweep(bulk, 1L, exp(lb), "/"), mu, c_reg)
  a <- plug$scale$a_hat[match(colnames(bulk), plug$scale$sample)]
  fit <- sweep(plug$fitted, 2L, a, "/")
  obs <- sweep(bulk, 2L, a, "/")
  log((obs + pseudo) / (fit + pseudo))
}

#' Plug-in gene-wise platform ratio from a pool of matched samples.
#'
#' Alternates the deMixSC plug-in on `bulk / b` with an update of
#' `log b_g` (RCTD stage 1 in spirit: the ratio is learnt from many
#' samples, then held fixed), until the largest change is below `tol`.
#' After each update `median_g log b_g` is set to 0, because a common
#' factor is absorbed by the scale `a_i`.
#'
#' - `method = "gene"`: `log b_g` is the median over the pooled
#'   replicates of [platform_ratio_residuals()]; `var_log_b` is their
#'   variance divided by the number of replicates of the gene.
#' - `method = "length"`: `log b_g = f(log2 L_g)`, a robust loess
#'   (`family = "symmetric"`, `span`) of the per-gene median log ratio on
#'   the log2 mean transcript length; `var_log_b` is the squared
#'   standard error of the loess prediction. Genes without a length get
#'   `log b = 0` and `var_log_b = 0`.
#'
#' @param bulk_list,mu_list Named lists (times) of bulk counts and
#'   per-cell reference means, as in [platform_ratio_residuals()].
#' @param train_times Times whose replicates are pooled.
#' @param method `"gene"` or `"length"`.
#' @param log2_length Named vector of log2 mean transcript lengths
#'   (needed for `"length"`).
#' @param n_iter,tol Maximum alternations and tolerance on `log b`.
#' @param span Loess span for `"length"`. Default `0.5`.
#' @param c_reg,pseudo Passed to [platform_ratio_residuals()].
#' @return A data frame with `gene`, `log_b`, `var_log_b`, `n_obs`
#'   (pooled replicates), with attributes `train_times`, `method` and
#'   `iterations`.
estimate_platform_ratio <- function(
  bulk_list,
  mu_list,
  train_times,
  method = c("gene", "length"),
  log2_length = NULL,
  n_iter = 10L,
  tol = 1e-3,
  span = 0.5,
  c_reg = 2,
  pseudo = 0.01
) {
  method <- match.arg(method)
  if (!all(train_times %in% names(bulk_list))) {
    stop("Every `train_times` entry must name an element of `bulk_list`.")
  }
  if (method == "length" && (is.null(log2_length) || is.null(names(log2_length)))) {
    stop("`method = \"length\"` needs a named `log2_length`.")
  }
  genes <- Reduce(union, lapply(mu_list[train_times], rownames))
  log_b <- stats::setNames(rep(0, length(genes)), genes)
  var_log_b <- log_b
  n_obs <- log_b
  it <- 0L
  for (it in seq_len(n_iter)) {
    res <- lapply(train_times, function(tm) {
      platform_ratio_residuals(
        bulk_list[[tm]],
        mu_list[[tm]],
        log_b,
        c_reg = c_reg,
        pseudo = pseudo
      )
    })
    pooled <- matrix(
      NA_real_,
      nrow = length(genes),
      ncol = sum(vapply(res, ncol, integer(1))),
      dimnames = list(genes, NULL)
    )
    col <- 0L
    for (r in res) {
      idx <- col + seq_len(ncol(r))
      pooled[rownames(r), idx] <- r
      col <- col + ncol(r)
    }
    n_obs <- rowSums(!is.na(pooled))
    med <- apply(pooled, 1L, stats::median, na.rm = TRUE)
    if (method == "gene") {
      new <- med
      var_new <- apply(pooled, 1L, stats::var, na.rm = TRUE) / pmax(n_obs, 1)
    } else {
      x <- log2_length[genes]
      ok <- is.finite(x) & is.finite(med)
      fit <- stats::loess(
        y ~ x,
        data = data.frame(x = x[ok], y = med[ok]),
        span = span,
        degree = 1L,
        family = "symmetric"
      )
      new <- stats::setNames(rep(0, length(genes)), genes)
      var_new <- new
      inside <- is.finite(x) & x >= min(x[ok]) & x <= max(x[ok])
      pr <- stats::predict(fit, newdata = data.frame(x = x[inside]), se = TRUE)
      new[inside] <- pr$fit
      var_new[inside] <- pr$se.fit^2
    }
    new[!is.finite(new)] <- 0
    var_new[!is.finite(var_new)] <- 0
    new <- new - stats::median(new)
    change <- max(abs(new - log_b))
    log_b <- new
    var_log_b <- var_new
    if (change < tol) {
      break
    }
  }
  out <- data.frame(
    gene = genes,
    log_b = as.numeric(log_b),
    var_log_b = as.numeric(var_log_b),
    n_obs = as.integer(n_obs),
    stringsAsFactors = FALSE
  )
  attr(out, "train_times") <- train_times
  attr(out, "method") <- method
  attr(out, "iterations") <- it
  out
}

#' Cross-fitted platform ratio: each time from the other times only.
#'
#' For every time `t`, [estimate_platform_ratio()] on the replicates of
#' the other times. The ratio used at `t` therefore never sees a
#' replicate of `t`, so it cannot absorb the residual it is meant to
#' test, and its transfer across compositions (48 h epiblast, 72 h
#' primitive streak, 96 h presomitic and cardiac mesoderm) is what is
#' tested.
#'
#' @inheritParams estimate_platform_ratio
#' @return A named list (times) of [estimate_platform_ratio()] outputs.
cross_fit_platform_ratio <- function(
  bulk_list,
  mu_list,
  method = c("gene", "length"),
  log2_length = NULL,
  ...
) {
  method <- match.arg(method)
  times <- names(bulk_list)
  if (length(times) < 2L) {
    stop("Cross-fitting needs at least two matched times.")
  }
  out <- lapply(times, function(tm) {
    train <- setdiff(times, tm)
    stopifnot(!tm %in% train)
    estimate_platform_ratio(
      bulk_list,
      mu_list,
      train_times = train,
      method = method,
      log2_length = log2_length,
      ...
    )
  })
  names(out) <- times
  out
}

#' Does the gene-wise platform ratio follow transcript length?
#'
#' A full-length bulk library yields reads in proportion to molecules
#' times transcript length; a 3' UMI library counts molecules. If that is
#' the whole platform effect, the log2 ratio of bulk to reference grows
#' with slope 1 in log2 length. Per time: median regression
#' (`quantreg::rq`, `tau = 0.5`) of the log2 ratio on log2 mean
#' transcript length, with a gene-bootstrap percentile interval of the
#' slope (genes resampled as if independent), Spearman's correlation and
#' the share of variance explained by a loess on length.
#'
#' @param log2_ratio Genes by times matrix of log2 ratios
#'   ([platform_log_ratio()], element `log_ratio`).
#' @param log2_length Named vector of log2 mean transcript lengths.
#' @param n_boot Bootstrap replicates. Default `1000`.
#' @param seed Seed of the bootstrap.
#' @return A data frame with one row per time: `time_point`, `n_genes`,
#'   `slope`, `slope_lower`, `slope_upper`, `intercept`, `spearman`,
#'   `loess_r2`.
platform_length_test <- function(log2_ratio, log2_length, n_boot = 1000L, seed = 1L) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package quantreg is required.")
  }
  rows <- lapply(colnames(log2_ratio), function(tm) {
    y <- log2_ratio[, tm]
    x <- log2_length[rownames(log2_ratio)]
    ok <- is.finite(x) & is.finite(y)
    d <- data.frame(x = x[ok], y = y[ok])
    fit <- quantreg::rq(y ~ x, tau = 0.5, data = d)
    boot <- withr::with_seed(seed, vapply(
      seq_len(n_boot),
      function(b) {
        idx <- sample.int(nrow(d), replace = TRUE)
        stats::coef(quantreg::rq(y ~ x, tau = 0.5, data = d[idx, ]))[[2L]]
      },
      numeric(1)
    ))
    lo <- stats::loess(y ~ x, data = d, span = 0.5, degree = 1L)
    data.frame(
      time_point = tm,
      n_genes = nrow(d),
      slope = stats::coef(fit)[[2L]],
      slope_lower = stats::quantile(boot, 0.025, names = FALSE),
      slope_upper = stats::quantile(boot, 0.975, names = FALSE),
      intercept = stats::coef(fit)[[1L]],
      spearman = stats::cor(d$x, d$y, method = "spearman"),
      loess_r2 = 1 - stats::var(stats::residuals(lo)) / stats::var(d$y),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Summary of the A3 shift per time, for one platform-ratio variant.
#'
#' @param tails Named list (times) of [bulk_tail_filter()] outputs.
#' @param range Named list (times) of [bulk_consistency_filter()] outputs.
#' @param variant Label of the variant.
#' @param level Tail level of the statements. Default `0.01`.
#' @return A data frame, one row per time: `variant`, `time_point`,
#'   `n_genes`, `median_abs_z`, `share_abs_z_gt` (share of gene by
#'   replicate values with `|z| > qnorm(1 - level)`), the genes flagged
#'   in at least one replicate by T1, T4 and the quantile band, the range
#'   and MAD flags, and the IQR of the per-gene log2 ratio.
a3_shift_summary <- function(tails, range, variant, level = 0.01) {
  lr <- platform_log_ratio(tails)$log_ratio
  z_cut <- stats::qnorm(1 - level)
  rows <- lapply(names(tails), function(tm) {
    g <- tails[[tm]]$gene
    z <- g$z[g$statement == "T4"]
    flagged <- function(s) length(unique(g$gene[g$statement == s & g$flag %in% TRUE]))
    band <- tails[[tm]]$quantile_band
    data.frame(
      variant = variant,
      time_point = tm,
      n_genes = length(unique(g$gene)),
      median_abs_z = stats::median(abs(z), na.rm = TRUE),
      share_abs_z_gt = mean(abs(z) > z_cut, na.rm = TRUE),
      flagged_T1 = flagged("T1"),
      flagged_T4 = flagged("T4"),
      flagged_QR = length(unique(band$gene[band$flag %in% TRUE])),
      flagged_range = sum(range[[tm]]$flag_range),
      flagged_mad = sum(range[[tm]]$flag_mad),
      iqr_log2_ratio = stats::IQR(lr[, tm], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
