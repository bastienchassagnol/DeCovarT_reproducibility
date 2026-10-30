#' Count-scale omnibus screen of cell-type effects.
#'
#' Replaces the between/within ratio of SCTransform Pearson residuals
#' (`R/rank_sctransform_naive.R`) by its count-scale counterpart:
#' a Gamma-Poisson (negative-binomial) GLM `y_cg ~ cell_type` with offset
#' `log s_c`, tested against `~ 1` by the quasi-likelihood F test of
#' `glmGamPoi::test_de()` (Lund et al. 2012). The test is a **gate**:
#' with one specimen, cells are not independent replicates, so the
#' p-values describe a screen inside that specimen, not population
#' differential expression. Genes that pass the gate are **ranked by an
#' effect size**, the deviance explained by the cell type with equal
#' type weights, so that a large type cannot dominate the rank.
#'
#' Requires `R/utils_general.R` to be sourced first.

#' Detection in every cell type.
#'
#' DeCovarT needs a non-singular within-type covariance for every type,
#' so a gene must be detected in a minimum fraction of the cells of
#' **every** type, not only in one. The gate also removes perfect
#' markers (zero in some type), which break identifiability once
#' gene-wise platform ratios are allowed.
#'
#' @param counts Genes by cells (dense or sparse `Matrix`).
#' @param cell_type Cell-type label per column of `counts`.
#' @param min_frac Minimum fraction of cells with a non-zero count, in
#'   each type. Default `0.1`; a tuning choice, report it.
#' @return A data frame with `gene`, `min_detection` (smallest fraction
#'   over types), `weakest_type`, and `pass`, plus the attribute
#'   `detection` (genes by types matrix of fractions).
presence_gate <- function(counts, cell_type, min_frac = 0.1) {
  cell_type <- droplevels(factor(cell_type))
  if (length(cell_type) != ncol(counts)) {
    stop("`cell_type` must have one label per column of `counts`.")
  }
  if (anyNA(cell_type)) {
    stop("`cell_type` has missing labels; drop those cells first.")
  }
  types <- levels(cell_type)
  detection <- vapply(
    types,
    function(tp) {
      idx <- which(cell_type == tp)
      Matrix::rowSums(counts[, idx, drop = FALSE] > 0) / length(idx)
    },
    numeric(nrow(counts))
  )
  detection <- matrix(
    detection,
    nrow = nrow(counts),
    dimnames = list(rownames(counts), types)
  )
  min_detection <- apply(detection, 1L, min)
  out <- data.frame(
    gene = rownames(counts) %||% as.character(seq_len(nrow(counts))),
    min_detection = min_detection,
    weakest_type = types[max.col(-detection, ties.method = "first")],
    pass = min_detection >= min_frac,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "detection") <- detection
  out
}

#' Negative-binomial unit deviances.
#'
#' Parameterisation of `glmGamPoi`: `Var y = mu + phi mu^2`. For
#' `phi = 0` the Poisson deviance is returned.
#'
#' @param y Counts (matrix or vector).
#' @param mu Means, same shape as `y`.
#' @param phi Overdispersion, recycled by row (one value per gene when
#'   `y` is genes by cells).
#' @return Unit deviances, same shape as `y`.
nb_unit_deviance <- function(y, mu, phi) {
  y <- as.matrix(y)
  mu <- as.matrix(mu)
  phi <- matrix(phi, nrow = nrow(y), ncol = ncol(y))
  term_y <- ifelse(y > 0, y * log(y / mu), 0)
  poisson <- 2 * (term_y - (y - mu))
  nb <- 2 *
    (term_y - (y + 1 / phi) * log((1 + phi * y) / (1 + phi * mu)))
  ifelse(phi > 0, nb, poisson)
}

#' Deviance explained by the cell type, with equal type weights.
#'
#' Let `r_gj` be the fitted rate of gene `g` in type `j` (fitted mean
#' divided by the size factor) and `r_g = J^-1 sum_j r_gj` the
#' equal-weight null rate. With `d` the NB unit deviance at the fitted
#' overdispersion and `D_gj(.)` its mean over the cells of type `j`,
#' `R2_eq(g) = 1 - sum_j D_gj(s_c r_gj) / sum_j D_gj(s_c r_g)`.
#' Each type contributes its mean deviance, so a type's weight does not
#' grow with its number of cells. The null rate is a plug-in, not the
#' weighted NB maximum-likelihood intercept, so `R2_eq` is an effect
#' size for ranking, not a test statistic.
#'
#' @param counts Genes by cells.
#' @param mu Fitted means of the full model, genes by cells.
#' @param size_factors Size factor `s_c` per cell.
#' @param cell_type Cell-type label per cell.
#' @param phi Overdispersion per gene.
#' @return Numeric vector, one value per gene.
deviance_explained_equal_type <- function(
  counts,
  mu,
  size_factors,
  cell_type,
  phi
) {
  cell_type <- droplevels(factor(cell_type))
  types <- levels(cell_type)
  rate <- vapply(
    types,
    function(tp) {
      idx <- which(cell_type == tp)
      rowMeans(sweep(mu[, idx, drop = FALSE], 2L, size_factors[idx], "/"))
    },
    numeric(nrow(mu))
  )
  rate <- matrix(rate, nrow = nrow(mu))
  null_rate <- rowMeans(rate)
  dev_full <- numeric(nrow(mu))
  dev_null <- numeric(nrow(mu))
  for (j in seq_along(types)) {
    idx <- which(cell_type == types[[j]])
    y_j <- as.matrix(counts[, idx, drop = FALSE])
    mu_full <- outer(rate[, j], size_factors[idx])
    mu_null <- outer(null_rate, size_factors[idx])
    dev_full <- dev_full + rowMeans(nb_unit_deviance(y_j, mu_full, phi))
    dev_null <- dev_null + rowMeans(nb_unit_deviance(y_j, mu_null, phi))
  }
  out <- rep(NA_real_, nrow(mu))
  ok <- dev_null > 0
  out[ok] <- 1 - dev_full[ok] / dev_null[ok]
  out
}

#' Omnibus cell-type screen on counts: gate by FDR, rank by effect size.
#'
#' Steps:
#' 1. [presence_gate()] in every type.
#' 2. `glmGamPoi::glm_gp(design = ~ cell_type, offset = log(s_c))`. The
#'    size factors are passed through `offset` with `size_factors` set
#'    to one, because `glm_gp()` adds `offset` to `log(size_factors)`.
#' 3. `glmGamPoi::test_de(reduced_design = ~ 1)`: quasi-likelihood F
#'    test of the cell-type effect, BH-adjusted.
#' 4. Genes with `adj_pval <= fdr` are ranked by
#'    [deviance_explained_equal_type()] (ties by the F statistic); the
#'    first `n_select` are selected.
#'
#' @param counts Genes by cells, raw UMI counts (dense matrix or an
#'   object `glm_gp()` accepts).
#' @param cell_type Cell-type label per column.
#' @param size_factors Size factor `s_c` per cell, for example the full
#'   library size divided by its median.
#' @param fdr BH level of the gate. Default `0.05`.
#' @param min_frac Passed to [presence_gate()].
#' @param n_select Number of genes to flag as selected. Default `500`.
#' @param overdispersion_shrinkage Passed to [glmGamPoi::glm_gp()].
#' @param verbose Passed to `glm_gp()` and `test_de()`.
#' @return A data frame with one row per gene of `counts`: `gene`,
#'   `min_detection`, `weakest_type`, `pass_presence`, `f_statistic`,
#'   `df1`, `df2`, `pval`, `adj_pval`, `pass_fdr`, `overdispersion`,
#'   `deviance_r2_equal_type`, `rank` (`NA` outside the gate),
#'   `selected`.
rank_genes_omnibus <- function(
  counts,
  cell_type,
  size_factors,
  fdr = 0.05,
  min_frac = 0.1,
  n_select = 500L,
  overdispersion_shrinkage = TRUE,
  verbose = FALSE
) {
  cell_type <- droplevels(factor(cell_type))
  if (nlevels(cell_type) < 2L) {
    stop("Need at least two cell types for an omnibus test.")
  }
  if (length(size_factors) != ncol(counts) || any(size_factors <= 0)) {
    stop("`size_factors` must be positive, one per column of `counts`.")
  }
  if (is.null(rownames(counts))) {
    stop("`counts` needs gene names as row names.")
  }
  gate <- presence_gate(counts, cell_type, min_frac)
  out <- data.frame(
    gene = gate$gene,
    min_detection = gate$min_detection,
    weakest_type = gate$weakest_type,
    pass_presence = gate$pass,
    f_statistic = NA_real_,
    df1 = NA_real_,
    df2 = NA_real_,
    pval = NA_real_,
    adj_pval = NA_real_,
    pass_fdr = FALSE,
    overdispersion = NA_real_,
    deviance_r2_equal_type = NA_real_,
    rank = NA_integer_,
    selected = FALSE,
    stringsAsFactors = FALSE
  )
  keep <- gate$pass
  if (sum(keep) == 0L) {
    warning("No gene passes the presence gate.", call. = FALSE)
    return(out)
  }
  y <- counts[keep, , drop = FALSE]
  col_data <- data.frame(cell_type = cell_type)
  fit <- glmGamPoi::glm_gp(
    y,
    design = ~cell_type,
    col_data = col_data,
    offset = log(size_factors),
    size_factors = rep(1, ncol(y)),
    overdispersion_shrinkage = overdispersion_shrinkage,
    verbose = verbose
  )
  de <- glmGamPoi::test_de(
    fit,
    reduced_design = ~1,
    pval_adjust_method = "BH",
    verbose = verbose
  )
  de <- de[match(rownames(y), de$name), , drop = FALSE]
  r2 <- deviance_explained_equal_type(
    y,
    fit$Mu,
    size_factors,
    cell_type,
    fit$overdispersions
  )
  idx <- which(keep)
  out$f_statistic[idx] <- de$f_statistic
  out$df1[idx] <- de$df1
  out$df2[idx] <- de$df2
  out$pval[idx] <- de$pval
  # The BH adjustment covers the genes that passed the presence gate.
  out$adj_pval[idx] <- de$adj_pval
  out$pass_fdr[idx] <- !is.na(de$adj_pval) & de$adj_pval <= fdr
  out$overdispersion[idx] <- fit$overdispersions
  out$deviance_r2_equal_type[idx] <- r2

  gated <- which(out$pass_fdr & is.finite(out$deviance_r2_equal_type))
  ord <- gated[order(
    -out$deviance_r2_equal_type[gated],
    -out$f_statistic[gated],
    out$gene[gated]
  )]
  out$rank[ord] <- seq_along(ord)
  out$selected <- !is.na(out$rank) & out$rank <= n_select
  out
}
