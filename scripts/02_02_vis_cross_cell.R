# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_02_vis_cross_cell.R \
#   > "logs/02_02_vis_cross_cell_$(date +%F).log" 2>&1 &
#
# Tier (ii) figures from output/feature_selection/tier_02/cross_cell.rds:
# null calibration, agreement of the four tests, gate and rank,
# overdispersion. PDFs in output/feature_selection/tier_02/figures/, PNG
# copies in docs/images/feature_selection/.

root <- getwd()
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
  calib$ratio_to_nominal[calib$test == tt & calib$time_point == tm & calib$level == lv]
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
      "Analytic tests reject %.1f to %.1f times too often at 48 h (two cell lines); grey band: 95 %% interval under uniformity",
      min(r48),
      max(r48)
    )
  ),
  "fs_tier2_null_calibration",
  fig_dir,
  width = 11,
  height = 4.4
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
  width = 11,
  height = 4.6
)
n_pass <- vapply(res$ranked, function(r) sum(r$pass_fdr %in% TRUE), numeric(1))
small <- vapply(res$ranked, function(r) {
  mean(r$deviance_r2[r$pass_fdr %in% TRUE] < 0.05)
}, numeric(1))
save_figure(
  plot_gate_rank(
    res$ranked,
    n_perm = 2000L,
    title = sprintf(
      "%s genes pass at 48 / 72 / 96 h; %.0f to %.0f %% of them explain < 5 %% of the deviance",
      paste(fmt_n(n_pass), collapse = " / "),
      100 * min(small),
      100 * max(small)
    ),
    subtitle = "One point per tested gene; BH at 0.05 on the permuted Welch G, no cap and no effect-size floor"
  ),
  "fs_tier2_gate_rank",
  fig_dir,
  width = 11,
  height = 4.6
)
save_figure(
  plot_overdispersion(
    res$tests,
    gene_means,
    title = "Overdispersion falls with the mean; a group of genes sits at the Poisson limit",
    subtitle = "Gene-wise maximum-likelihood 1/psi of M1 against the mean size-normalised count; orange: loess trend"
  ),
  "fs_tier2_overdispersion",
  fig_dir,
  width = 11,
  height = 4.2
)
message("Wrote ", fig_dir)
