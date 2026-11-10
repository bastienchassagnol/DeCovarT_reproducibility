#' Leverage, standardised residuals and Cook's distance for cells and genes.
#'
#' Tier (iii) of the gene selection, "Outlier cells and genes". Two fits
#' are diagnosed:
#'
#' - **Cells**: the one-way negative-binomial GLM M1 of one gene,
#'   `log m_gc = beta_0g + beta_g,j(c) + log s_c`, fitted on the cells of
#'   one time slot ([cook_cells()]).
#' - **Genes**: the weighted least-squares regression of one bulk sample
#'   on the reference signature, across genes ([cook_genes()]).
#'
#' A cell or gene is **flagged** only when its leverage and its
#' standardised residual are both large. High leverage alone marks a
#' cell of a rare type or a highly expressed marker, which carries
#' information; a large residual alone marks an ordinary tail count.
#' Nothing is clipped or down-weighted here: flags are reported, and the
#' decision is taken in the chapter.
#'
#' This file and `R/02_03_stats_tail_bounds.R` (discrepancy and drift between
#' the reference and the bulk) are the two files of tier (iii).
#'
#' Port of the formulas of `scripts/figures/fig_03_1_influence.R` in the
#' DeCovarT repository (checked there against `hatvalues()`,
#' `rstandard()` and `cooks.distance()`); deliberately duplicated, no
#' cross-repository dependency. Requires `R/02_00_utils.R`.

#' Cell-level influence in the one-way NB GLM M1.
#'
#' For one gene, with fitted means `m_c` and overdispersion `phi`
#' (`V(m) = m + phi m^2`), the IRLS weight is `w_c = m_c^2 / V(m_c)`. In
#' a one-way design the hat matrix is block diagonal by type, so
#'
#' - leverage `h_c = w_c / sum_{c' in j(c)} w_c'` (mean leverage `J / C`);
#' - standardised Pearson residual
#'   `r_c = (x_c - m_c) / sqrt(dispersion V(m_c) (1 - h_c))`;
#' - Cook's distance `D_c = r_p^2 h_c / (J dispersion (1 - h_c)^2)`, with
#'   `r_p` the raw Pearson residual.
#'
#' `dispersion = 1` is the convention of `cooks.distance()` on a negbin
#' fit. A cell is flagged when `h_c > h_mult J / C` **and**
#' `|r_c| > r_max`.
#'
#' @param counts Genes by cells, raw counts, gene names as row names.
#' @param cell_type Cell-type label per column.
#' @param size_factors Size factor per cell (enters as `log` offset).
#' @param mu Optional genes by cells fitted means of M1 (for example
#'   `fit$Mu` of the fit in [compare_type_effect_tests()]). When `NULL`,
#'   M1 is fitted with `glmGamPoi::glm_gp()`.
#' @param overdispersion Optional overdispersion per gene, required when
#'   `mu` is given.
#' @param dispersion Scale in the residual and Cook formulas. Default 1.
#' @param h_mult Multiple of the mean leverage `J / C` above which the
#'   leverage is called large. Default `3`.
#' @param r_max Threshold on `|standardised Pearson residual|`.
#'   Default `3`.
#' @param block Genes per block, to bound memory. Default `200`.
#' @return A list with
#'   - `flagged`: one row per flagged (gene, cell): `gene`, `cell`,
#'     `cell_type`, `count`, `fitted`, `leverage`, `std_pearson`, `cook`;
#'   - `summary`: one row per gene: `gene`, `n_flagged`, `max_cook`,
#'     `cell_max_cook`, `max_leverage`, `n_high_leverage`,
#'     `n_large_residual`.
cook_cells <- function(
  counts,
  cell_type,
  size_factors,
  mu = NULL,
  overdispersion = NULL,
  dispersion = 1,
  h_mult = 3,
  r_max = 3,
  block = 200L
) {
  cell_type <- droplevels(factor(cell_type))
  if (length(cell_type) != ncol(counts)) {
    stop("`cell_type` must have one label per column of `counts`.")
  }
  if (length(size_factors) != ncol(counts) || any(size_factors <= 0)) {
    stop("`size_factors` must be positive, one per column of `counts`.")
  }
  if (is.null(rownames(counts))) {
    stop("`counts` needs gene names as row names.")
  }
  cells <- colnames(counts) %||% as.character(seq_len(ncol(counts)))
  if (is.null(mu)) {
    if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
      stop(
        "Package glmGamPoi is required to fit M1, ",
        "or pass `mu` and `overdispersion` from an existing fit."
      )
    }
    fit <- glmGamPoi::glm_gp(
      counts,
      design = ~cell_type,
      col_data = data.frame(cell_type = cell_type),
      offset = log(size_factors),
      size_factors = rep(1, ncol(counts)),
      overdispersion = TRUE,
      overdispersion_shrinkage = FALSE
    )
    mu <- fit$Mu
    overdispersion <- fit$overdispersions
  } else if (is.null(overdispersion)) {
    stop("`overdispersion` is required when `mu` is supplied.")
  }
  if (nrow(mu) != nrow(counts) || ncol(mu) != ncol(counts)) {
    stop("`mu` must have the dimensions of `counts` (genes by cells).")
  }
  if (length(overdispersion) != nrow(counts)) {
    stop("`overdispersion` must have one value per gene of `counts`.")
  }
  n_types <- nlevels(cell_type)
  n_cells <- ncol(counts)
  h_cut <- h_mult * n_types / n_cells
  type_id <- as.integer(cell_type)
  z <- Matrix::sparseMatrix(
    i = seq_len(n_cells),
    j = type_id,
    x = 1,
    dims = c(n_cells, n_types)
  )

  flagged <- list()
  summary_rows <- list()
  for (s in seq(1L, nrow(counts), by = block)) {
    rows <- s:min(s + block - 1L, nrow(counts))
    x <- as.matrix(counts[rows, , drop = FALSE])
    m <- as.matrix(mu[rows, , drop = FALSE])
    phi <- overdispersion[rows]
    v <- m + phi * m^2
    w <- m^2 / v
    w_type <- as.matrix(w %*% z)
    h <- w / w_type[, type_id, drop = FALSE]
    r_p <- (x - m) / sqrt(v)
    r_sp <- r_p / sqrt(dispersion * (1 - h))
    cook <- r_p^2 * h / (n_types * dispersion * (1 - h)^2)
    flag <- h > h_cut & abs(r_sp) > r_max
    flag[is.na(flag)] <- FALSE
    hit <- which(flag, arr.ind = TRUE)
    if (nrow(hit) > 0L) {
      flagged[[length(flagged) + 1L]] <- data.frame(
        gene = rownames(counts)[rows][hit[, 1L]],
        cell = cells[hit[, 2L]],
        cell_type = as.character(cell_type)[hit[, 2L]],
        count = x[hit],
        fitted = m[hit],
        leverage = h[hit],
        std_pearson = r_sp[hit],
        cook = cook[hit],
        stringsAsFactors = FALSE
      )
    }
    cook_max_col <- max.col(
      replace(cook, is.na(cook), -Inf),
      ties.method = "first"
    )
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      gene = rownames(counts)[rows],
      n_flagged = rowSums(flag),
      max_cook = apply(cook, 1L, max, na.rm = TRUE),
      cell_max_cook = cells[cook_max_col],
      max_leverage = apply(h, 1L, max, na.rm = TRUE),
      n_high_leverage = rowSums(h > h_cut, na.rm = TRUE),
      n_large_residual = rowSums(abs(r_sp) > r_max, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  empty <- data.frame(
    gene = character(),
    cell = character(),
    cell_type = character(),
    count = numeric(),
    fitted = numeric(),
    leverage = numeric(),
    std_pearson = numeric(),
    cook = numeric(),
    stringsAsFactors = FALSE
  )
  out <- list(
    flagged = if (length(flagged) > 0L) do.call(rbind, flagged) else empty,
    summary = do.call(rbind, summary_rows)
  )
  attr(out, "thresholds") <- c(leverage = h_cut, std_pearson = r_max)
  out
}

#' Gene-level influence in the bulk-on-reference weighted regression.
#'
#' One bulk sample `y` (genes) is regressed on the reference signature
#' `X` (genes by types) without intercept, with weights `w_g`. When the
#' composition comes from a non-negative fit, the diagnostics are those
#' of the unconstrained weighted fit restricted to the types with a
#' positive coefficient (the passive set, on which the NNLS solution is
#' the least-squares solution). With `P` the passive types, `p = |P|`,
#' residual `e_g = y_g - x_g' beta`,
#'
#' - leverage `h_g = w_g x_gP' (X_P' W X_P)^{-1} x_gP` (they sum to `p`);
#' - `sigma2 = sum_g w_g e_g^2 / (G - p)`;
#' - standardised residual `r_g = sqrt(w_g) e_g / sqrt(sigma2 (1 - h_g))`;
#' - Cook's distance `D_g = w_g e_g^2 h_g / (p sigma2 (1 - h_g)^2)`.
#'
#' A gene is flagged when `h_g > h_mult p / G` **and** `|r_g| > r_max`.
#'
#' @param bulk Bulk expression of one sample: a named numeric vector over
#'   genes, or a genes by samples matrix (one fit per column).
#' @param signature Genes by types reference means, same scale as `bulk`,
#'   gene names as row names.
#' @param coef Optional coefficients per type (a vector, or a types by
#'   samples matrix), for example `a_i p_i` from
#'   [plugin_composition_demixsc()]. When `NULL`, the unconstrained
#'   weighted least-squares fit is used.
#' @param weights Optional weights per gene (vector, or genes by samples
#'   matrix), for example the final deMixSC weights. Default: all one.
#' @param h_mult,r_max Flag thresholds. Defaults `3` and `3`.
#' @return A data frame with one row per gene and sample: `sample`,
#'   `gene`, `bulk`, `fitted`, `weight`, `leverage`, `std_residual`,
#'   `cook`, `flag`.
cook_genes <- function(
  bulk,
  signature,
  coef = NULL,
  weights = NULL,
  h_mult = 3,
  r_max = 3
) {
  signature <- as.matrix(signature)
  genes <- rownames(signature)
  if (is.null(genes)) {
    stop("`signature` needs gene names as row names.")
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
  n_genes <- nrow(signature)
  w_mat <- if (is.null(weights)) {
    matrix(1, n_genes, ncol(bulk))
  } else {
    matrix(as.numeric(weights), n_genes, ncol(bulk))
  }
  if (any(!is.finite(w_mat)) || any(w_mat <= 0)) {
    stop("`weights` must be positive and finite.")
  }
  coef_mat <- if (is.null(coef)) {
    NULL
  } else {
    matrix(as.numeric(coef), ncol(signature), ncol(bulk))
  }

  rows <- lapply(seq_len(ncol(bulk)), function(i) {
    y <- bulk[, i]
    w <- w_mat[, i]
    beta <- if (is.null(coef_mat)) {
      xw <- signature * sqrt(w)
      drop(qr.coef(qr(xw), y * sqrt(w)))
    } else {
      coef_mat[, i]
    }
    beta[is.na(beta)] <- 0
    passive <- beta > 0
    if (!any(passive)) {
      stop(
        "No type has a positive coefficient in sample ",
        colnames(bulk)[i],
        "."
      )
    }
    x_p <- signature[, passive, drop = FALSE]
    n_par <- sum(passive)
    fitted <- drop(x_p %*% beta[passive])
    e <- y - fitted
    xw <- x_p * sqrt(w)
    xtwx_inv <- solve(crossprod(xw))
    h <- rowSums((xw %*% xtwx_inv) * xw)
    sigma2 <- sum(w * e^2) / (n_genes - n_par)
    r <- sqrt(w) * e / sqrt(sigma2 * (1 - h))
    cook <- w * e^2 * h / (n_par * sigma2 * (1 - h)^2)
    data.frame(
      sample = colnames(bulk)[i],
      gene = genes,
      bulk = y,
      fitted = fitted,
      weight = w,
      leverage = h,
      std_residual = r,
      cook = cook,
      flag = h > h_mult * n_par / n_genes & abs(r) > r_max,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
