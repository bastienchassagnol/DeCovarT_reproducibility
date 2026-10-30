#' Cross-platform consistency filter between bulk and single-cell profiles.
#'
#' **Heuristic.** This adapts the two gene filters of RNA-Sieve
#' (Erdmann-Pham et al. 2021, Genome Research, eq. 12 and 13) to a
#' DeCovarT panel. Neither rule is derived from DeCovarT's likelihood:
#' the range rule follows from the mixture constraint, the MAD rules are
#' empirical (the RNA-Sieve authors say their choice was "primarily
#' guided by empirical considerations"). A flagged gene is a candidate
#' for removal, to be reported, not a proven platform artefact.
#'
#' In Suppinger, the single-cell reference matches the standard bulk
#' only at 48 h. The filter is therefore fitted on the 48 h bulk and
#' the 48 h reference, and the resulting gene list is applied at every
#' time. That transfer is an assumption: platform effects are taken to
#' be gene-specific and constant over time.
#'
#' Scale. Each type profile and each bulk sample is normalised to sum
#' one over the panel genes. If the bulk is a mixture of the types,
#' `b = sum_k alpha_k m_k`, the normalised bulk is a convex combination
#' of the normalised type profiles (weights proportional to
#' `alpha_k * sum_g m_gk`). So, gene by gene, the normalised bulk lies in
#' `[min_k m_gk, max_k m_gk]` up to sampling error, whatever the
#' proportions and whatever a global platform scale.
#'
#' Requires `R/utils_general.R` to be sourced first.

#' Flag genes whose bulk expression is inconsistent with the reference.
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
#'   samples matrix of replicates at the matched time (48 h). Row names
#'   (or names) must be gene names.
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
    "RNA-Sieve range bound and MAD outlier flags, fitted at 48 h only;",
    "not derived from the DeCovarT likelihood."
  )
  out
}
