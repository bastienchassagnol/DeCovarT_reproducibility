# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_02_vis_cross_cell.R \
#   > "logs/02_02_vis_cross_cell_$(date +%F).log" 2>&1 &
#
# Tier (ii) figures from output/feature_selection/tier_02/cross_cell.rds:
# null calibration, agreement of the four tests, gate and rank,
# overdispersion. PDFs in output/feature_selection/tier_02/figures/, PNG
# copies in docs/images/feature_selection/.

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_02_vis_cross_cell.R"))

tier_dir <- file.path(root, "output", "feature_selection", "tier_02")
fig_dir <- file.path(tier_dir, "figures")
res <- readRDS(file.path(tier_dir, "cross_cell.rds"))
prepared <- readRDS(
  file.path(root, "output", "feature_selection", "tier_01", "scored_slots.rds")
)
gene_means <- lapply(names(res$tests), function(tm) {
  slot <- prepared$slots[[tm]]
  Matrix::rowMeans(Matrix::t(Matrix::t(slot$counts) / res$sf[[tm]]))
})
names(gene_means) <- names(res$tests)

calib <- res$calibration
ratio_at <- function(tt, tm, lv = 0.01) {
  calib$ratio_to_nominal[
    calib$test == tt & calib$time_point == tm & calib$level == lv
  ]
}
analytic <- c("lr", "anova_vst", "ql_f")
r48 <- vapply(analytic, ratio_at, numeric(1), tm = "48h")
save_figure(
  plot_null_calibration(
    res$null_tests,
    title = sprintf(
      "At 48 h only the permuted Welch G holds its level under permuted labels (%.2f x nominal at 0.01)",
      ratio_at("welch_perm", "48h")
    ),
    subtitle = sprintf(
      "Analytic tests reject %.1f to %.1f times too often at 48 h (two cell lines); grey band: pointwise 95 %% interval under uniformity",
      min(r48),
      max(r48)
    )
  ),
  "fs_tier2_null_calibration",
  fig_dir,
  data = list(null_tests = res$null_tests, calibration = calib),
  width = 11,
  height = 8.4
)
tau_wl <- vapply(res$tau_all, function(m) m["lr", "welch_perm"], numeric(1))
tau_lq <- vapply(res$tau_all, function(m) m["lr", "ql_f"], numeric(1))
save_figure(
  plot_test_agreement(
    res$tau_all,
    title = sprintf(
      "LR and QL F order genes almost identically (tau %.2f); Welch G differs most (tau %.2f to %.2f with LR)",
      min(tau_lq),
      min(tau_wl),
      max(tau_wl)
    ),
    subtitle = "Kendall's tau between the p-value rankings of all tested genes\nLR: likelihood ratio; ANOVA: on the shifted log; QL F: quasi-likelihood F; Welch G: permuted Welch G"
  ),
  "fs_tier2_test_agreement",
  fig_dir,
  data = res$tau_all,
  width = 11,
  height = 4.6
)
n_perm <- attr(res$tests[[1L]], "n_perm")
if (is.null(n_perm)) {
  stop("cross_cell.rds has no n_perm attribute; rerun the stats script.")
}
n_pass <- vapply(res$ranked, function(r) sum(r$pass_fdr %in% TRUE), numeric(1))
n_floor <- vapply(
  res$ranked,
  function(r) sum(gate_class(r, n_perm) == "At the floor 1/(B+1)", na.rm = TRUE),
  numeric(1)
)
save_figure(
  plot_gate_rank(
    res$ranked,
    n_perm = n_perm,
    title = sprintf(
      "%s genes pass BH; %s beat all %s permutations (48 / 72 / 96 h)",
      paste(fmt_n(n_pass), collapse = " / "),
      paste(fmt_n(n_floor), collapse = " / "),
      fmt_n(n_perm)
    ),
    subtitle = "One point per tested gene; stripes are the discrete values k/(B+1) of the permutation p-value"
  ),
  "fs_tier2_gate_rank",
  fig_dir,
  data = lapply(res$ranked, function(r) {
    out <- r[r$pass_presence, c("gene", "pval", "adj_pval", "pass_fdr", "deviance_r2"), drop = FALSE]
    out$gate_class <- gate_class(r[r$pass_presence, , drop = FALSE], n_perm)
    out
  }),
  width = 11,
  height = 4.8
)
save_figure(
  plot_r2_vs_mean(
    res$ranked,
    gene_means,
    n_perm = n_perm,
    title = "Low-mean genes with a large R2_D fail the gate; floor genes mostly have a mean above 0.1",
    subtitle = "Between-type share of deviance against mean size-normalised count, per tested gene"
  ),
  "fs_tier2_r2_vs_mean",
  fig_dir,
  data = lapply(names(res$ranked), function(tm) {
    r <- res$ranked[[tm]]
    r <- r[r$pass_presence, , drop = FALSE]
    data.frame(gene = r$gene, mean = gene_means[[tm]][r$gene], r2 = r$deviance_r2, gate_class = gate_class(r, n_perm))
  }),
  width = 11,
  height = 4.6
)
save_figure(
  plot_pvalue_histograms(
    res$tests,
    title = "Every test has a peak at 0 over a flat remainder; the permuted Welch G is discrete",
    subtitle = "Unadjusted p-values of the tested genes; dashed: height of a uniform histogram of pi0 m genes"
  ),
  "fs_tier2_pvalue_histograms",
  fig_dir,
  data = lapply(res$tests, function(r) {
    r[r$pass_presence, c("gene", "lr_pval", "anova_vst_pval", "ql_f_pval", "welch_perm_pval"), drop = FALSE]
  }),
  width = 11,
  height = 8
)
tertiles <- lapply(names(res$tests), function(tm) {
  r <- res$tests[[tm]]
  m <- gene_means[[tm]][r$gene]
  q <- stats::quantile(m[r$pass_presence %in% TRUE], c(1, 2) / 3, na.rm = TRUE)
  factor(
    cut(m, c(-Inf, q, Inf), labels = c("low mean", "middle mean", "high mean")),
    levels = c("low mean", "middle mean", "high mean")
  )
})
names(tertiles) <- names(res$tests)
save_figure(
  plot_pvalue_histograms(
    res$tests,
    strata = tertiles,
    test = "welch_perm",
    title = "Low-expression genes carry less signal: the peak at 0 grows with the mean",
    subtitle = "Permuted Welch G p-values by tertile of the mean size-normalised count (tertiles per time)"
  ),
  "fs_tier2_pvalue_histograms_by_mean",
  fig_dir,
  data = list(tertiles = tertiles),
  width = 11,
  height = 6.4
)
save_figure(
  plot_overdispersion(
    res$tests,
    gene_means,
    title = "Overdispersion falls with the mean; low-mean genes sit at the Poisson boundary",
    subtitle = "Gene-wise maximum-likelihood phi of M1 against the mean size-normalised count; orange: loess trend of the overdispersed genes"
  ),
  "fs_tier2_overdispersion",
  fig_dir,
  data = list(
    overdispersion = lapply(res$tests, function(r) {
      r[r$pass_presence, c("gene", "overdispersion"), drop = FALSE]
    }),
    gene_means = gene_means
  ),
  width = 11,
  height = 4.2
)
message("Wrote ", fig_dir)
