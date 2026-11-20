#' Which type drives the type effect: post hoc contrasts after the gate.
#'
#' Theory (hypotheses A and B, stage-wise gene-level OFDR, graphical
#' within-gene procedure, dropped strategies): DeCovarT vignette 3.1,
#' `#sec-fs-deg-posthoc`. This file holds the implementation.
#'
#' Per gene g and scored type j of one slot:
#'
#' - A_gj: type j is above the equal-weighted average of the other types
#'   on the log scale, Delta_gj = eta_gj - mean_{k != j} eta_gk > 0;
#' - B_gj: type j is above every other type, eta_gj - eta_gk > 0 for all
#'   k != j (intersection-union p-value: the largest of the pairwise
#'   one-sided p-values). B_gj implies A_gj.
#'
#' Two-sided variant (returned beside the one-sided values, same fits and
#' same permutation draw):
#'
#' - A2_gj: Delta_gj != 0 (two-sided p);
#' - B2_gj: type j is above every other type or below every other type,
#'   in the sign of Delta_gj (intersection-union of the two-sided pairwise
#'   p-values, set to 1 when a pairwise estimate has the other sign).
#'   Delta_gj is the mean of the pairwise contrasts eta_gj - eta_gk, so
#'   B2_gj implies A2_gj with the same sign and the graph carries over.
#'
#' Four families of one-sided p-values, one per omnibus test of tier (ii):
#'
#' - `lr`: Wald z on the NB M1 coefficients (cell-means design, gene-wise
#'   maximum-likelihood overdispersion), with the closed-form covariance
#'   of a one-way design, Var(beta_gj) = 1 / sum_{c in j} m_gc /
#'   (1 + phi_g m_gc); equal to `emmeans` on the same fit (checked by
#'   [emmeans_contrast_check()]);
#' - `ql_f`: `glmGamPoi::test_de()` with a contrast, on the fit with
#'   shrunk overdispersion; one-sided from the sign of the fold change;
#' - `anova_vst`: t on the shifted logarithm with the pooled residual
#'   variance (the contrasts `emmeans` gives on a one-way `lm`);
#' - `welch_perm`: Welch t on the Pearson residuals of M0 (own
#'   overdispersion), calibrated by the same within-stratum label
#'   permutations as the gate (same seed, same draw), so the gate's
#'   p-value is reproduced as a check.

#' Contrast matrices of one slot.
#'
#' @param types Type names (J >= 2).
#' @return List with `A` (J x J: row j is type j minus the mean of the
#'   others) and `pairs` (data frame of ordered pairs j, k with j != k)
#'   plus `P` (rows = ordered pairs, j minus k).
contrast_matrices <- function(types) {
  n_types <- length(types)
  a <- matrix(-1 / (n_types - 1), n_types, n_types, dimnames = list(types, types))
  diag(a) <- 1
  pairs <- expand.grid(j = seq_len(n_types), k = seq_len(n_types))
  pairs <- pairs[pairs$j != pairs$k, , drop = FALSE]
  p <- matrix(0, nrow(pairs), n_types, dimnames = list(NULL, types))
  p[cbind(seq_len(nrow(pairs)), pairs$j)] <- 1
  p[cbind(seq_len(nrow(pairs)), pairs$k)] <- -1
  list(A = a, pairs = pairs, P = p)
}

#' One-sided p-values for A and the intersection-union p for B.
#'
#' @param est_a,se_a Genes x J estimates and standard errors of A.
#' @param est_p,se_p Genes x (ordered pairs) estimates and standard errors.
#' @param pairs Data frame of ordered pairs from [contrast_matrices()].
#' @param df Degrees of freedom (`Inf` for z).
#' @return List of two genes x J matrices, `p_a` and `p_b`.
one_sided_ab <- function(est_a, se_a, est_p, se_p, pairs, df = Inf) {
  upper <- function(t) stats::pt(t, df = df, lower.tail = FALSE)
  p_a <- upper(est_a / se_a)
  p_pair <- upper(est_p / se_p)
  p_b <- vapply(
    seq_len(ncol(est_a)),
    function(j) apply(p_pair[, pairs$j == j, drop = FALSE], 1L, max),
    numeric(nrow(est_a))
  )
  dim(p_b) <- dim(p_a)
  list(p_a = p_a, p_b = p_b)
}

#' Sign-consistent intersection-union p for B2 from two-sided pairwise p.
#'
#' @param p_pair Genes x (ordered pairs) two-sided p-values.
#' @param est_a Genes x J estimates of A (their sign is the direction).
#' @param est_p Genes x (ordered pairs) estimates.
#' @param pairs Data frame of ordered pairs from [contrast_matrices()].
#' @return Genes x J matrix: the largest pairwise p of type j, or 1 when a
#'   pairwise estimate of j does not have the sign of its A estimate.
iut_two_sided <- function(p_pair, est_a, est_p, pairs) {
  p_b <- vapply(
    seq_len(ncol(est_a)),
    function(j) {
      cols <- pairs$j == j
      same <- rowSums(sign(est_p[, cols, drop = FALSE]) != sign(est_a[, j])) == 0
      pb <- apply(p_pair[, cols, drop = FALSE], 1L, max)
      ifelse(same %in% TRUE | is.na(pb), pb, 1)
    },
    numeric(nrow(est_a))
  )
  dim(p_b) <- dim(est_a)
  p_b
}

#' Two-sided p-values for A2 and the sign-consistent p for B2.
#'
#' @inheritParams one_sided_ab
#' @return List of genes x J matrices `p_a`, `p_b` and `sign_a`.
two_sided_ab <- function(est_a, se_a, est_p, se_p, pairs, df = Inf) {
  both <- function(t) 2 * stats::pt(-abs(t), df = df)
  list(
    p_a = both(est_a / se_a),
    p_b = iut_two_sided(both(est_p / se_p), est_a, est_p, pairs),
    sign_a = sign(est_a)
  )
}

#' Within-gene graphical procedure, A_j then B_j.
#'
#' Hypotheses A_1..A_J (weight 1/J each) and B_1..B_J (weight 0);
#' A_j passes all its weight to B_j (fixed sequence, B_j implies A_j);
#' B_j passes its weight to the other A_k equally. Tested with
#' `graphicalMCP::graph_test_shortcut()` at level `alpha`.
#'
#' @param n_types Number of types J.
#' @param types Names of the types.
#' @return A `graphicalMCP` initial graph.
graph_a_then_b <- function(n_types, types) {
  m <- 2L * n_types
  tr <- matrix(0, m, m)
  for (j in seq_len(n_types)) {
    tr[j, n_types + j] <- 1
    tr[n_types + j, setdiff(seq_len(n_types), j)] <- 1 / (n_types - 1)
  }
  graphicalMCP::graph_create(
    c(rep(1 / n_types, n_types), rep(0, n_types)),
    tr,
    c(paste0("A_", types), paste0("B_", types))
  )
}

#' Stage-wise decisions for one test of one slot.
#'
#' Stage 1: the gene passes its omnibus test at BH `fdr` (R of G genes).
#' Stage 2, per screened gene: level `alpha_I = fdr R / G`; three designs:
#' the graphical A-then-B procedure, Holm on the A family only, and (not
#' stage-wise) BH over all gene x type A p-values of the slot. The Holm
#' step uses the multipliers J - 1, J - 1, J - 2, ..., 1, as `stageR`'s
#' "holm" method: with delta = 0 the screen rejects the intersection of
#' the A nulls, so at most J - 1 of them can still be true (Shaffer's
#' logical restriction). The graph keeps the plain weights 1/J
#' (conservative).
#'
#' @param p_a,p_b Genes x J one-sided p-values.
#' @param omnibus_padj BH-adjusted omnibus p-value per gene (same rows).
#' @param fdr Target level. Default `0.05`.
#' @return Data frame, one row per gene x type: `screened`, `alpha_I`,
#'   `a_graph`, `b_graph`, `a_holm`, `a_flat_bh`.
stagewise_decisions <- function(p_a, p_b, omnibus_padj, fdr = 0.05) {
  n_genes <- nrow(p_a)
  n_types <- ncol(p_a)
  types <- colnames(p_a)
  screened <- !is.na(omnibus_padj) & omnibus_padj <= fdr
  alpha_i <- fdr * sum(screened) / sum(!is.na(omnibus_padj))
  a_graph <- b_graph <- a_holm <- matrix(FALSE, n_genes, n_types)
  g0 <- graph_a_then_b(n_types, types)
  for (i in which(screened)) {
    pv <- c(p_a[i, ], p_b[i, ])
    if (anyNA(pv)) {
      next
    }
    res <- graphicalMCP::graph_test_shortcut(g0, p = pv, alpha = alpha_i)
    rej <- res$outputs$rejected
    a_graph[i, ] <- rej[seq_len(n_types)]
    b_graph[i, ] <- rej[n_types + seq_len(n_types)]
    o <- order(p_a[i, ])
    mult <- c(n_types - 1, (n_types - 1):1)
    adj <- cummax(pmin(p_a[i, o] * mult, 1))
    a_holm[i, o] <- adj <= alpha_i
  }
  flat <- matrix(stats::p.adjust(as.vector(p_a), method = "BH"), n_genes) <= fdr
  data.frame(
    row = rep(seq_len(n_genes), times = n_types),
    type = rep(types, each = n_genes),
    screened = rep(screened, times = n_types),
    alpha_I = alpha_i,
    a_graph = as.vector(a_graph),
    b_graph = as.vector(b_graph),
    a_holm = as.vector(a_holm),
    a_flat_bh = as.vector(flat) %in% TRUE,
    stringsAsFactors = FALSE
  )
}

#' Post hoc contrasts of the type effect on one slot.
#'
#' @param counts Genes x cells raw counts of the genes to test (those that
#'   pass the presence gate), gene names as row names.
#' @param cell_type Type label per cell (fixed before any statistic).
#' @param size_factors Size factor per cell (as in tier (ii)).
#' @param strata Stratum per cell for the permutations (the batch).
#' @param n_perm,seed Permutations and seed; use the gate's values so the
#'   draw is the same.
#' @param vst_alpha Overdispersion of the shifted logarithm. Default `0.05`.
#' @param perm_block Permutations per matrix product. Default `250`.
#' @return List: `log2fc` (genes x J, A contrast of the NB M1 fit, log2),
#'   `log_mean` (genes x J, beta of M1), `p` (named list over tests of
#'   lists `p_a`, `p_b`), `p2` (the same for A2 and B2, two-sided, plus
#'   `sign_a`), `welch_gate_p` (gate p-value recomputed),
#'   `types`, `n_j`.
type_contrast_tests <- function(
  counts,
  cell_type,
  size_factors,
  strata,
  n_perm,
  seed,
  vst_alpha = 0.05,
  perm_block = 250L
) {
  cell_type <- droplevels(factor(cell_type))
  types <- levels(cell_type)
  n_types <- length(types)
  n_cells <- ncol(counts)
  n_j <- as.numeric(table(cell_type))
  cm <- contrast_matrices(types)
  col_data <- data.frame(cell_type = cell_type)
  offset <- log(size_factors)
  unit_sf <- rep(1, n_cells)
  z <- Matrix::sparseMatrix(
    i = seq_len(n_cells),
    j = as.integer(cell_type),
    x = 1,
    dims = c(n_cells, n_types)
  )
  p <- list()
  p2 <- list()

  # LR / Wald: NB M1, cell-means design, ML overdispersion (as the gate).
  fit1 <- glmGamPoi::glm_gp(
    counts,
    design = ~ 0 + cell_type,
    col_data = col_data,
    offset = offset,
    size_factors = unit_sf,
    overdispersion = TRUE,
    overdispersion_shrinkage = FALSE
  )
  beta <- fit1$Beta
  colnames(beta) <- types
  phi <- fit1$overdispersions
  w <- fit1$Mu / (1 + phi * fit1$Mu)
  info <- as.matrix(w %*% z)
  v_beta <- 1 / info
  est_a <- beta %*% t(cm$A)
  se_a <- sqrt(v_beta %*% t(cm$A^2))
  est_p <- beta %*% t(cm$P)
  se_p <- sqrt(v_beta %*% t(cm$P^2))
  p$lr <- one_sided_ab(est_a, se_a, est_p, se_p, cm$pairs)
  p2$lr <- two_sided_ab(est_a, se_a, est_p, se_p, cm$pairs)
  log2fc <- est_a / log(2)

  # QL F: contrasts on the fit with shrunk overdispersion.
  fit_ql <- glmGamPoi::glm_gp(
    counts,
    design = ~ 0 + cell_type,
    col_data = col_data,
    offset = offset,
    size_factors = unit_sf,
    overdispersion = TRUE,
    overdispersion_shrinkage = TRUE
  )
  ql_one <- function(cvec) {
    de <- glmGamPoi::test_de(fit_ql, contrast = cvec)
    de[match(rownames(counts), de$name), c("pval", "lfc"), drop = FALSE]
  }
  de_a <- lapply(seq_len(n_types), function(j) ql_one(cm$A[j, ]))
  de_p <- lapply(seq_len(nrow(cm$P)), function(r) ql_one(cm$P[r, ]))
  pick <- function(des, col) vapply(des, `[[`, numeric(nrow(counts)), col)
  qa2 <- pick(de_a, "pval")
  qp2 <- pick(de_p, "pval")
  la <- pick(de_a, "lfc")
  lp <- pick(de_p, "lfc")
  dim(qa2) <- dim(la) <- c(nrow(counts), n_types)
  dim(qp2) <- dim(lp) <- c(nrow(counts), nrow(cm$P))
  qa <- ifelse(la > 0, qa2 / 2, 1 - qa2 / 2)
  qp <- ifelse(lp > 0, qp2 / 2, 1 - qp2 / 2)
  qb <- vapply(
    seq_len(n_types),
    function(j) apply(qp[, cm$pairs$j == j, drop = FALSE], 1L, max),
    numeric(nrow(counts))
  )
  dim(qb) <- dim(qa)
  p$ql_f <- list(p_a = qa, p_b = qb)
  p2$ql_f <- list(
    p_a = qa2,
    p_b = iut_two_sided(qp2, la, lp, cm$pairs),
    sign_a = sign(la)
  )

  # ANOVA on the shifted logarithm, pooled residual variance.
  y0 <- 1 / (4 * vst_alpha)
  df2 <- n_cells - n_types
  means <- matrix(NA_real_, nrow(counts), n_types)
  mse <- numeric(nrow(counts))
  for (s in seq(1L, nrow(counts), by = 500L)) {
    rows <- s:min(s + 499L, nrow(counts))
    v <- log(sweep(as.matrix(counts[rows, , drop = FALSE]), 2L, size_factors, "/") + y0)
    sums <- as.matrix(v %*% z)
    means[rows, ] <- sweep(sums, 2L, n_j, "/")
    mse[rows] <- (rowSums(v^2) - rowSums(sweep(sums^2, 2L, n_j, "/"))) / df2
  }
  inv_n <- 1 / n_j
  est_a <- means %*% t(cm$A)
  se_a <- sqrt(outer(mse, as.vector(cm$A^2 %*% inv_n)))
  est_p <- means %*% t(cm$P)
  se_p <- sqrt(outer(mse, as.vector(cm$P^2 %*% inv_n)))
  p$anova_vst <- one_sided_ab(est_a, se_a, est_p, se_p, cm$pairs, df = df2)
  p2$anova_vst <- two_sided_ab(est_a, se_a, est_p, se_p, cm$pairs, df = df2)

  # Welch t on the Pearson residuals of M0, permutation-calibrated.
  fit0 <- glmGamPoi::glm_gp(
    counts,
    design = ~1,
    col_data = col_data,
    offset = offset,
    size_factors = unit_sf,
    overdispersion = TRUE,
    overdispersion_shrinkage = FALSE
  )
  m0 <- fit0$Mu
  r0 <- (as.matrix(counts) - m0) / sqrt(m0 + fit0$overdispersions * m0^2)
  r0_sq <- r0^2
  tot <- rowSums(r0)
  tot2 <- rowSums(r0_sq)
  welch_t <- function(s1, s2) {
    mean_j <- sweep(s1, 2L, n_j, "/")
    var_j <- sweep(s2 - sweep(s1^2, 2L, n_j, "/"), 2L, n_j - 1, "/")
    var_j[var_j <= 0] <- NA_real_
    vm <- sweep(var_j, 2L, n_j, "/")
    list(
      a = (mean_j %*% t(cm$A)) / sqrt(vm %*% t(cm$A^2)),
      p = (mean_j %*% t(cm$P)) / sqrt(vm %*% t(cm$P^2))
    )
  }
  s1 <- as.matrix(r0 %*% z)
  s2 <- as.matrix(r0_sq %*% z)
  obs <- welch_t(s1, s2)
  obs_g <- welch_g_stat(s1, s2, n_j, tot, tot2)[, "G"]
  # Same draw as compare_type_effect_tests(): same seed, same order.
  labels <- as.integer(cell_type)
  by_stratum <- split(seq_len(n_cells), strata)
  permute_once <- function() {
    out_lab <- labels
    for (idx in by_stratum) {
      out_lab[idx] <- labels[idx][sample.int(length(idx))]
    }
    out_lab
  }
  draw <- function() replicate(n_perm, permute_once())
  perm_lab <- matrix(withr::with_seed(seed, draw()), nrow = n_cells)
  ex_a <- matrix(0, nrow(counts), n_types)
  ex_p <- matrix(0, nrow(counts), nrow(cm$P))
  ex_a2 <- ex_a
  ex_p2 <- ex_p
  ex_g <- numeric(nrow(counts))
  for (b0 in seq(1L, n_perm, by = perm_block)) {
    bb <- b0:min(b0 + perm_block - 1L, n_perm)
    zp <- Matrix::sparseMatrix(
      i = rep(seq_len(n_cells), length(bb)),
      j = as.vector(sweep(perm_lab[, bb, drop = FALSE], 2L, (seq_along(bb) - 1L) * n_types, "+")),
      x = 1,
      dims = c(n_cells, length(bb) * n_types)
    )
    ps1 <- as.matrix(r0 %*% zp)
    ps2 <- as.matrix(r0_sq %*% zp)
    for (k in seq_along(bb)) {
      cols <- (k - 1L) * n_types + seq_len(n_types)
      st <- welch_t(ps1[, cols, drop = FALSE], ps2[, cols, drop = FALSE])
      ex_a <- ex_a + (!is.na(st$a) & st$a >= obs$a)
      ex_p <- ex_p + (!is.na(st$p) & st$p >= obs$p)
      ex_a2 <- ex_a2 + (!is.na(st$a) & abs(st$a) >= abs(obs$a))
      ex_p2 <- ex_p2 + (!is.na(st$p) & abs(st$p) >= abs(obs$p))
      g_star <- welch_g_stat(ps1[, cols, drop = FALSE], ps2[, cols, drop = FALSE], n_j, tot, tot2)[, "G"]
      ex_g <- ex_g + (!is.na(g_star) & g_star >= obs_g)
    }
  }
  pa <- (1 + ex_a) / (n_perm + 1)
  pa[is.na(obs$a)] <- NA_real_
  pp <- (1 + ex_p) / (n_perm + 1)
  pp[is.na(obs$p)] <- NA_real_
  pb <- vapply(
    seq_len(n_types),
    function(j) apply(pp[, cm$pairs$j == j, drop = FALSE], 1L, max),
    numeric(nrow(counts))
  )
  dim(pb) <- dim(pa)
  p$welch_perm <- list(p_a = pa, p_b = pb, t_a = obs$a)
  pa2 <- (1 + ex_a2) / (n_perm + 1)
  pa2[is.na(obs$a)] <- NA_real_
  pp2 <- (1 + ex_p2) / (n_perm + 1)
  pp2[is.na(obs$p)] <- NA_real_
  p2$welch_perm <- list(
    p_a = pa2,
    p_b = iut_two_sided(pp2, obs$a, obs$p, cm$pairs),
    sign_a = sign(obs$a),
    t_a = obs$a
  )
  gate_p <- (1 + ex_g) / (n_perm + 1)
  gate_p[is.na(obs_g)] <- NA_real_

  for (tt in names(p)) {
    dimnames(p[[tt]]$p_a) <- dimnames(p[[tt]]$p_b) <- list(rownames(counts), types)
    dimnames(p2[[tt]]$p_a) <- dimnames(p2[[tt]]$p_b) <-
      dimnames(p2[[tt]]$sign_a) <- list(rownames(counts), types)
  }
  dimnames(log2fc) <- list(rownames(counts), types)
  rownames(beta) <- rownames(counts)
  list(
    log2fc = log2fc,
    log_mean = beta,
    v_beta = v_beta,
    p = p,
    p2 = p2,
    welch_gate_p = stats::setNames(gate_p, rownames(counts)),
    types = types,
    n_j = stats::setNames(n_j, types)
  )
}

#' Cross-check of the Wald A and B p-values with emmeans.
#'
#' For a few genes, builds an `emmeans` reference grid from the M1
#' coefficients and their covariance with `emmeans::qdrg()`, and compares
#' `contrast(method = "del.eff", side = ">")` (A) and the one-sided
#' pairwise contrasts (B, intersection-union max) with the vectorised
#' values.
#'
#' @param res Output of [type_contrast_tests()].
#' @param genes Genes to check.
#' @param side `">"` (one-sided A and B) or `"two-sided"` (A2 and B2).
#' @return Data frame: `gene`, `side`, `max_abs_diff_a`, `max_abs_diff_b`.
emmeans_contrast_check <- function(res, genes, side = c(">", "two-sided")) {
  side <- match.arg(side)
  two <- identical(side, "two-sided")
  types <- res$types
  n_types <- length(types)
  grid <- data.frame(cell_type = factor(types, levels = types))
  do.call(rbind, lapply(genes, function(g) {
    rg <- emmeans::qdrg(
      ~ 0 + cell_type,
      data = grid,
      coef = unname(res$log_mean[g, ]),
      vcov = diag(unname(res$v_beta[g, ]), n_types),
      df = Inf
    )
    emm <- emmeans::emmeans(rg, "cell_type")
    a <- summary(
      emmeans::contrast(emm, method = "del.eff", adjust = "none"),
      side = if (two) "=" else ">",
      infer = c(FALSE, TRUE)
    )$p.value
    pw <- summary(
      emmeans::contrast(emm, method = "revpairwise", adjust = "none"),
      infer = c(FALSE, TRUE)
    )
    # Pairwise one-sided p for "j above k" from the two-sided z.
    zval <- pw$z.ratio
    lab <- strsplit(as.character(pw$contrast), " - ", fixed = TRUE)
    # emmeans wraps level names that contain "-" in parentheses.
    unwrap <- function(z) sub("^\\((.*)\\)$", "\\1", z)
    hi <- unwrap(vapply(lab, `[`, "", 1L))
    lo <- unwrap(vapply(lab, `[`, "", 2L))
    one <- c(stats::pnorm(zval, lower.tail = FALSE), stats::pnorm(zval))
    top <- c(hi, lo)
    if (two) {
      # B2: largest two-sided pairwise p, 1 unless every pairwise sign
      # matches the sign of the A estimate.
      est_a <- summary(emmeans::contrast(emm, method = "del.eff"))$estimate
      names(est_a) <- types
      sg <- c(sign(zval), -sign(zval))
      p2 <- rep(pw$p.value, 2L)
      b <- vapply(types, function(tp) {
        if (all(sg[top == tp] == sign(est_a[[tp]]))) max(p2[top == tp]) else 1
      }, numeric(1))
      ref <- res$p2$lr
    } else {
      b <- vapply(types, function(tp) max(one[top == tp]), numeric(1))
      ref <- res$p$lr
    }
    data.frame(
      gene = g,
      side = side,
      max_abs_diff_a = max(abs(a - ref$p_a[g, ])),
      max_abs_diff_b = max(abs(b - ref$p_b[g, ]))
    )
  }))
}
