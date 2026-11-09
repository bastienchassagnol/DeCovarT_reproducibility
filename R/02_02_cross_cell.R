#' Count-scale screen of cell-type effects: four tests, one effect size.
#'
#' Tier (ii) of the gene selection (book chapter "Feature selection",
#' sections "Count-model fit per time" and "Four tests of the type
#' effect"). Per time slot and per gene, two negative-binomial GLMs are
#' fitted by `glmGamPoi::glm_gp()` with offset `log s_c`:
#' M0 = intercept only and M1 = intercept + cell type. No `SCTransform`
#' fit is involved: the overdispersion is estimated gene by gene and no
#' residual is clipped.
#'
#' Four tests of the type effect are computed side by side
#' ([compare_type_effect_tests()]): the likelihood-ratio test (analysis of
#' deviance), a one-way ANOVA F on the shifted logarithm of
#' size-normalised counts, the quasi-likelihood F test of
#' `glmGamPoi::test_de()` (Lund et al. 2012), and a Freedman-Lane label
#' permutation of Welch's G on the Pearson residuals of M0 (Winkler et
#' al. 2014). With one specimen, cells are not independent replicates,
#' so every p-value describes a screen inside that specimen, not
#' population differential expression. One test is used as a
#' Benjamini-Hochberg **gate**; genes that pass are **ranked by an effect
#' size**, the cell-weighted deviance explained `1 - D(M1) / D(M0)`.
#'
#' `glmGamPoi` is not in `renv.lock`; it is checked with
#' `requireNamespace()`. Requires `R/02_00_utils.R` and
#' `R/02_01_stats_quality_control.R` (`presence_gate()`) to be sourced
#' first. Written, not run.

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

#' Total NB deviance per gene, computed in blocks of genes.
#'
#' @param y Genes by cells counts (dense or sparse).
#' @param mu Fitted means, genes by cells.
#' @param phi Overdispersion per gene.
#' @param block Number of genes per block, to bound memory.
#' @return Numeric vector, the sum over cells of the unit deviances.
nb_deviance_rows <- function(y, mu, phi, block = 500L) {
  out <- numeric(nrow(y))
  starts <- seq(1L, nrow(y), by = block)
  for (s in starts) {
    rows <- s:min(s + block - 1L, nrow(y))
    out[rows] <- rowSums(nb_unit_deviance(
      as.matrix(y[rows, , drop = FALSE]),
      mu[rows, , drop = FALSE],
      phi[rows]
    ))
  }
  out
}

#' Deviance explained by the cell type, cell-weighted.
#'
#' `R2_D(g) = 1 - D_g(M1) / D_g(M0)`, with both deviances summed over
#' every cell of the slot (each cell has weight one) and evaluated at the
#' same overdispersion, so that `D_g(M0) - D_g(M1)` is the likelihood-
#' ratio statistic. A large type contributes more cells and therefore
#' more deviance; this is the weighting of the likelihood itself. The
#' type-equal variant of the legacy rank is not used any more.
#'
#' @param deviance_m0,deviance_m1 Deviances of M0 and M1 per gene.
#' @return Numeric vector, `NA` where `deviance_m0` is not positive.
deviance_explained <- function(deviance_m0, deviance_m1) {
  out <- rep(NA_real_, length(deviance_m0))
  ok <- is.finite(deviance_m0) & is.finite(deviance_m1) & deviance_m0 > 0
  out[ok] <- 1 - deviance_m1[ok] / deviance_m0[ok]
  out
}

#' Welch's G and the one-way F from per-type sums.
#'
#' Welch (1951) heteroscedastic one-way statistic, as written in
#' Winkler et al. (2014, Eq. 6-7) for a one-way design.
#'
#' @param s1,s2 Genes by types sums of `r` and `r^2` per type.
#' @param n_j Cells per type.
#' @param tot,tot2 Row totals of `r` and `r^2`.
#' @return Matrix with columns `G` and `F`; `G` is `NA` when a type has
#'   zero variance.
welch_g_stat <- function(s1, s2, n_j, tot, tot2) {
  s1 <- as.matrix(s1)
  s2 <- as.matrix(s2)
  n_types <- length(n_j)
  n <- sum(n_j)
  mean_j <- sweep(s1, 2L, n_j, "/")
  var_j <- sweep(s2 - sweep(s1^2, 2L, n_j, "/"), 2L, n_j - 1, "/")
  var_j[var_j <= 0] <- NA_real_
  w <- sweep(1 / var_j, 2L, n_j, "*")
  w_tot <- rowSums(w)
  r_w <- rowSums(w * mean_j) / w_tot
  num <- rowSums(w * (mean_j - r_w)^2) / (n_types - 1)
  a <- rowSums(sweep((1 - w / w_tot)^2, 2L, n_j - 1, "/"))
  den <- 1 + 2 * (n_types - 2) / (n_types^2 - 1) * a
  ss_b <- rowSums(sweep(s1^2, 2L, n_j, "/")) - tot^2 / n
  ss_w <- tot2 - rowSums(sweep(s1^2, 2L, n_j, "/"))
  cbind(
    G = num / den,
    F = (ss_b / (n_types - 1)) / (ss_w / (n - n_types))
  )
}

#' Four tests of the cell-type effect on the counts of one time slot.
#'
#' Per gene that passes [presence_gate()]:
#'
#' 1. `lr`: analysis of deviance. M1 is fitted with a gene-wise maximum-
#'    likelihood overdispersion (`overdispersion_shrinkage = FALSE`), M0
#'    is refitted at the same overdispersion, and
#'    `Lambda = D(M0) - D(M1)` is referred to `chi^2_{J-1}`.
#' 2. `anova_vst`: one-way ANOVA F on the shifted logarithm
#'    `log(x / s_c + y0)`, `y0 = 1 / (4 vst_alpha)` (Ahlmann-Eltze and
#'    Huber 2023), against `F_{J-1, C-J}`.
#' 3. `ql_f`: quasi-likelihood F test of `glmGamPoi::test_de()` on an M1
#'    fit with empirical-Bayes shrinkage of the overdispersion.
#' 4. `welch_perm`: Freedman-Lane permutation. The reduced model is M0
#'    (intercept and offset, own overdispersion); its Pearson residuals
#'    carry no type information, so permuting the residuals against fixed
#'    labels equals permuting the labels. Welch's G is the statistic
#'    (pivotal under unequal variances); `n_perm` permutations of the
#'    labels are drawn once and shared by every gene;
#'    `p = (1 + #{G* >= G}) / (n_perm + 1)`.
#'
#' Size factors enter through `offset = log(size_factors)` with
#' `size_factors = 1` in `glm_gp()`, because `glm_gp()` adds `offset` to
#' `log(size_factors)`. BH adjustment covers the genes that pass the
#' presence gate, test by test. `n_perm` is fixed before the run; its
#' Monte Carlo floor is `1 / (n_perm + 1)`.
#'
#' @param counts Genes by cells, raw UMI counts, gene names as row names.
#' @param cell_type Cell-type label per column. Labels must be fixed
#'   before any statistic is computed (no refinement on the same counts).
#' @param size_factors Size factor `s_c` per cell, for example the full
#'   library size divided by its median.
#' @param n_perm Number of shared label permutations. Default `2000`.
#' @param seed Integer seed for the permutations, or `NULL` to use the
#'   current RNG state. Applied with `withr::with_seed()`.
#' @param fdr BH level used for the `<test>_pass` columns. Default `0.05`.
#' @param min_frac Passed to [presence_gate()]. Default `0.05`.
#' @param tests Tests to compute; any of `"lr"`, `"anova_vst"`, `"ql_f"`,
#'   `"welch_perm"`. The M0 and M1 fits behind `deviance_r2` always run.
#' @param vst_alpha Overdispersion fixed in the shifted logarithm.
#'   Default `0.05`.
#' @param perm_block Permutations per matrix product. Default `250`.
#' @param verbose Passed to `glm_gp()` and `test_de()`.
#' @return A data frame with one row per gene of `counts`: `gene`,
#'   `min_detection`, `weakest_type`, `pass_presence`, `overdispersion`
#'   (M1, gene-wise), `deviance_m0`, `deviance_m1`, `deviance_r2`, then
#'   for each test `<test>_stat`, `<test>_pval`, `<test>_padj` and
#'   `<test>_pass`. Attributes: `n_perm`, `seed`, `vst_alpha`, `df1`
#'   (`J - 1`) and `df2` (`C - J`).
compare_type_effect_tests <- function(
  counts,
  cell_type,
  size_factors,
  n_perm = 2000L,
  seed = NULL,
  fdr = 0.05,
  min_frac = 0.05,
  tests = c("lr", "anova_vst", "ql_f", "welch_perm"),
  vst_alpha = 0.05,
  perm_block = 250L,
  verbose = FALSE
) {
  if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
    stop(
      "Package glmGamPoi is required (not in renv.lock). Install it with ",
      "BiocManager::install(\"glmGamPoi\") and record it with renv::snapshot()."
    )
  }
  tests <- match.arg(tests, several.ok = TRUE)
  cell_type <- droplevels(factor(cell_type))
  if (nlevels(cell_type) < 2L) {
    stop("Need at least two cell types for a test of the type effect.")
  }
  if (length(cell_type) != ncol(counts)) {
    stop("`cell_type` must have one label per column of `counts`.")
  }
  if (
    length(size_factors) != ncol(counts) ||
      any(!is.finite(size_factors)) ||
      any(size_factors <= 0)
  ) {
    stop("`size_factors` must be positive, one per column of `counts`.")
  }
  if (is.null(rownames(counts))) {
    stop("`counts` needs gene names as row names.")
  }
  n_perm <- as.integer(n_perm)
  if ("welch_perm" %in% tests && (is.na(n_perm) || n_perm < 1L)) {
    stop("`n_perm` must be a positive integer when `welch_perm` is requested.")
  }
  if (!is.numeric(vst_alpha) || length(vst_alpha) != 1L || vst_alpha <= 0) {
    stop("`vst_alpha` must be one positive number.")
  }

  gate <- presence_gate(counts, cell_type, min_frac)
  n_genes <- nrow(counts)
  out <- data.frame(
    gene = gate$gene,
    min_detection = gate$min_detection,
    weakest_type = gate$weakest_type,
    pass_presence = gate$pass,
    overdispersion = NA_real_,
    deviance_m0 = NA_real_,
    deviance_m1 = NA_real_,
    deviance_r2 = NA_real_,
    stringsAsFactors = FALSE
  )
  for (tt in tests) {
    out[[paste0(tt, "_stat")]] <- NA_real_
    out[[paste0(tt, "_pval")]] <- NA_real_
    out[[paste0(tt, "_padj")]] <- NA_real_
    out[[paste0(tt, "_pass")]] <- FALSE
  }
  n_types <- nlevels(cell_type)
  n_cells <- ncol(counts)
  df1 <- n_types - 1L
  df2 <- n_cells - n_types
  attr(out, "n_perm") <- if ("welch_perm" %in% tests) n_perm else NA_integer_
  attr(out, "seed") <- seed
  attr(out, "vst_alpha") <- vst_alpha
  attr(out, "df1") <- df1
  attr(out, "df2") <- df2

  keep <- which(gate$pass)
  if (length(keep) == 0L) {
    warning("No gene passes the presence gate.", call. = FALSE)
    return(out)
  }
  y <- counts[keep, , drop = FALSE]
  col_data <- data.frame(cell_type = cell_type)
  offset <- log(size_factors)
  unit_sf <- rep(1, n_cells)

  # M1 with gene-wise ML overdispersion, and M0 at the same overdispersion:
  # the pair behind the LR statistic and the deviance R^2.
  fit1 <- glmGamPoi::glm_gp(
    y,
    design = ~cell_type,
    col_data = col_data,
    offset = offset,
    size_factors = unit_sf,
    overdispersion = TRUE,
    overdispersion_shrinkage = FALSE,
    verbose = verbose
  )
  phi1 <- fit1$overdispersions
  fit0 <- glmGamPoi::glm_gp(
    y,
    design = ~1,
    col_data = col_data,
    offset = offset,
    size_factors = unit_sf,
    overdispersion = phi1,
    overdispersion_shrinkage = FALSE,
    verbose = verbose
  )
  dev1 <- nb_deviance_rows(y, fit1$Mu, phi1)
  dev0 <- nb_deviance_rows(y, fit0$Mu, phi1)
  out$overdispersion[keep] <- phi1
  out$deviance_m0[keep] <- dev0
  out$deviance_m1[keep] <- dev1
  out$deviance_r2[keep] <- deviance_explained(dev0, dev1)

  set_test <- function(out, test, stat, pval) {
    padj <- stats::p.adjust(pval, method = "BH")
    out[[paste0(test, "_stat")]][keep] <- stat
    out[[paste0(test, "_pval")]][keep] <- pval
    out[[paste0(test, "_padj")]][keep] <- padj
    out[[paste0(test, "_pass")]][keep] <- !is.na(padj) & padj <= fdr
    out
  }

  if ("lr" %in% tests) {
    lambda <- pmax(dev0 - dev1, 0)
    out <- set_test(
      out,
      "lr",
      lambda,
      stats::pchisq(lambda, df1, lower.tail = FALSE)
    )
  }

  z <- Matrix::sparseMatrix(
    i = seq_len(n_cells),
    j = as.integer(cell_type),
    x = 1,
    dims = c(n_cells, n_types)
  )
  n_j <- as.numeric(table(cell_type))

  if ("anova_vst" %in% tests) {
    y0 <- 1 / (4 * vst_alpha)
    f_vst <- numeric(length(keep))
    for (s in seq(1L, length(keep), by = 500L)) {
      rows <- s:min(s + 499L, length(keep))
      v <- log(
        sweep(as.matrix(y[rows, , drop = FALSE]), 2L, size_factors, "/") + y0
      )
      v <- v - rowMeans(v)
      sums <- as.matrix(v %*% z)
      ss_b <- rowSums(sweep(sums^2, 2L, n_j, "/"))
      ss_w <- rowSums(v^2) - ss_b
      f_vst[rows] <- (ss_b / df1) / (ss_w / df2)
    }
    out <- set_test(
      out,
      "anova_vst",
      f_vst,
      stats::pf(f_vst, df1, df2, lower.tail = FALSE)
    )
  }

  if ("ql_f" %in% tests) {
    fit_ql <- glmGamPoi::glm_gp(
      y,
      design = ~cell_type,
      col_data = col_data,
      offset = offset,
      size_factors = unit_sf,
      overdispersion = TRUE,
      overdispersion_shrinkage = TRUE,
      verbose = verbose
    )
    de <- glmGamPoi::test_de(
      fit_ql,
      reduced_design = ~1,
      pval_adjust_method = "BH",
      verbose = verbose
    )
    de <- de[match(rownames(y), de$name), , drop = FALSE]
    out <- set_test(out, "ql_f", de$f_statistic, de$pval)
  }

  if ("welch_perm" %in% tests) {
    fit0_own <- glmGamPoi::glm_gp(
      y,
      design = ~1,
      col_data = col_data,
      offset = offset,
      size_factors = unit_sf,
      overdispersion = TRUE,
      overdispersion_shrinkage = FALSE,
      verbose = verbose
    )
    m0 <- fit0_own$Mu
    phi0 <- fit0_own$overdispersions
    r0 <- (as.matrix(y) - m0) / sqrt(m0 + phi0 * m0^2)
    r0_sq <- r0^2
    tot <- rowSums(r0)
    tot2 <- rowSums(r0_sq)
    obs <- welch_g_stat(
      as.matrix(r0 %*% z),
      as.matrix(r0_sq %*% z),
      n_j,
      tot,
      tot2
    )[, "G"]
    labels <- as.integer(cell_type)
    draw <- function() replicate(n_perm, sample(labels))
    perm_lab <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
    perm_lab <- matrix(perm_lab, nrow = n_cells)
    exceed <- numeric(length(keep))
    for (b0 in seq(1L, n_perm, by = perm_block)) {
      bb <- b0:min(b0 + perm_block - 1L, n_perm)
      # Stacked indicators: column (k - 1) J + j flags type j in
      # permutation k of the block, so one product gives all type sums.
      zp <- Matrix::sparseMatrix(
        i = rep(seq_len(n_cells), length(bb)),
        j = as.vector(sweep(
          perm_lab[, bb, drop = FALSE],
          2L,
          (seq_along(bb) - 1L) * n_types,
          "+"
        )),
        x = 1,
        dims = c(n_cells, length(bb) * n_types)
      )
      s1 <- as.matrix(r0 %*% zp)
      s2 <- as.matrix(r0_sq %*% zp)
      for (k in seq_along(bb)) {
        cols <- (k - 1L) * n_types + seq_len(n_types)
        g_star <- welch_g_stat(
          s1[, cols, drop = FALSE],
          s2[, cols, drop = FALSE],
          n_j,
          tot,
          tot2
        )[, "G"]
        exceed <- exceed + (!is.na(g_star) & g_star >= obs)
      }
    }
    p_perm <- (1 + exceed) / (n_perm + 1)
    p_perm[is.na(obs)] <- NA_real_
    out <- set_test(out, "welch_perm", obs, p_perm)
  }
  out
}

#' Kendall's tau between the gene orders of the type-effect tests.
#'
#' Tau-b between the p-value rankings, on the log scale (no tie at
#' underflow for the analytic tests). The permutation p-value has a
#' floor of `1 / (n_perm + 1)`; its ties are kept.
#'
#' @param res Output of [compare_type_effect_tests()].
#' @param tests Tests to compare; columns `<test>_pval` must exist.
#' @param genes Optional logical or character selection of genes (for
#'   example the genes that pass the gate of one test).
#' @return A symmetric matrix of Kendall's tau, tests by tests, computed
#'   on the genes with a p-value for every test.
type_test_kendall <- function(
  res,
  tests = c("lr", "anova_vst", "ql_f", "welch_perm"),
  genes = NULL
) {
  cols <- paste0(tests, "_pval")
  missing_cols <- setdiff(cols, names(res))
  if (length(missing_cols) > 0L) {
    stop(
      "Columns missing from `res`: ",
      paste(missing_cols, collapse = ", "),
      ". Run compare_type_effect_tests() with these tests."
    )
  }
  sel <- if (is.null(genes)) {
    rep(TRUE, nrow(res))
  } else if (is.logical(genes)) {
    genes
  } else {
    res$gene %in% genes
  }
  p <- as.matrix(res[sel, cols, drop = FALSE])
  p <- p[stats::complete.cases(p), , drop = FALSE]
  if (nrow(p) < 3L) {
    stop("Fewer than three genes with a p-value for every test.")
  }
  tau <- stats::cor(log(p), method = "kendall")
  dimnames(tau) <- list(tests, tests)
  tau
}

#' Calibration of the type-effect tests under permuted labels.
#'
#' Run [compare_type_effect_tests()] once more with `cell_type` replaced
#' by one random permutation of the labels: the type effect is then null
#' by construction for every gene. A test is conservative when fewer than
#' the nominal fraction of its p-values fall below a level, and
#' anti-conservative when more do.
#'
#' @param res_null Output of [compare_type_effect_tests()] on permuted
#'   labels.
#' @param tests Tests to summarise.
#' @param levels Nominal levels. Default `c(0.01, 0.05)`.
#' @return A data frame with `test`, `level`, `n_genes`, `fraction_below`
#'   and `ratio_to_nominal`.
type_test_calibration <- function(
  res_null,
  tests = c("lr", "anova_vst", "ql_f", "welch_perm"),
  levels = c(0.01, 0.05)
) {
  rows <- lapply(tests, function(tt) {
    p <- res_null[[paste0(tt, "_pval")]]
    if (is.null(p)) {
      stop("Column ", tt, "_pval is missing from `res_null`.")
    }
    p <- p[!is.na(p)]
    frac <- vapply(levels, function(a) mean(p <= a), numeric(1))
    data.frame(
      test = tt,
      level = levels,
      n_genes = length(p),
      fraction_below = frac,
      ratio_to_nominal = frac / levels,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Type-effect screen: gate by one test, rank by deviance explained.
#'
#' Steps, per time slot:
#' 1. [compare_type_effect_tests()] with the gate test only (plus the M0
#'    and M1 fits it always runs).
#' 2. Genes with `<gate_test>_padj <= fdr` are ranked by `deviance_r2`,
#'    the cell-weighted `1 - D(M1) / D(M0)` (ties by the LR statistic
#'    `deviance_m0 - deviance_m1`, then by gene name); the first
#'    `n_select` are selected.
#'
#' @inheritParams compare_type_effect_tests
#' @param gate_test The test behind the BH gate: `"ql_f"` (default),
#'   `"lr"`, `"anova_vst"` or `"welch_perm"`. Choose it before the ranks
#'   are seen.
#' @param n_select Number of genes to flag as selected. Default `500`.
#' @param ... Further arguments to [compare_type_effect_tests()]
#'   (`n_perm`, `seed`, `vst_alpha`, `verbose`).
#' @return The data frame of [compare_type_effect_tests()] with, in
#'   addition, `gate_test`, `pval`, `adj_pval`, `pass_fdr` (copies of the
#'   gate test's columns), `rank` (`NA` outside the gate) and `selected`.
rank_genes_omnibus <- function(
  counts,
  cell_type,
  size_factors,
  fdr = 0.05,
  min_frac = 0.05,
  n_select = 500L,
  gate_test = c("ql_f", "lr", "anova_vst", "welch_perm"),
  ...
) {
  gate_test <- match.arg(gate_test)
  out <- compare_type_effect_tests(
    counts,
    cell_type,
    size_factors,
    fdr = fdr,
    min_frac = min_frac,
    tests = gate_test,
    ...
  )
  out$gate_test <- gate_test
  out$pval <- out[[paste0(gate_test, "_pval")]]
  out$adj_pval <- out[[paste0(gate_test, "_padj")]]
  out$pass_fdr <- out[[paste0(gate_test, "_pass")]]
  out$rank <- NA_integer_
  gated <- which(out$pass_fdr & is.finite(out$deviance_r2))
  ord <- gated[order(
    -out$deviance_r2[gated],
    -(out$deviance_m0[gated] - out$deviance_m1[gated]),
    out$gene[gated]
  )]
  out$rank[ord] <- seq_along(ord)
  out$selected <- !is.na(out$rank) & out$rank <= n_select
  out
}
