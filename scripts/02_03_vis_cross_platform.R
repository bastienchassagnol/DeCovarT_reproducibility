# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_03_vis_cross_platform.R \
#   > "logs/02_03_vis_cross_platform_$(date +%F).log" 2>&1 &
#
# Tier (iii) figures from output/feature_selection/tier_03/cross_platform.rds:
# Cook diagnostics, platform ratio against transcript length, A3 shift under
# three models, observed against fitted bulk, tail-statement flags, plug-in
# composition, stability of the ratio across times. Titles are computed
# from the outputs. PDFs in output/feature_selection/tier_03/figures/, PNG
# copies in docs/images/feature_selection/.

root <- getwd()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_03_vis_cross_platform.R"))

tier_dir <- file.path(root, "output", "feature_selection", "tier_03")
fig_dir <- file.path(tier_dir, "figures")
res <- readRDS(file.path(tier_dir, "cross_platform.rds"))
prepared <- readRDS(
  file.path(root, "output", "feature_selection", "tier_01", "scored_slots.rds")
)
cells <- lapply(prepared$slots, function(slot) table(slot$cell_type))
a3 <- res$a3_shift
pct_range <- function(x) sprintf("%.0f to %.0f %%", 100 * min(x), 100 * max(x))

# Cook ----------------------------------------------------------------------
n_gene_flags <- sum(vapply(res$outlier_genes, function(g) sum(g$flag), numeric(1)))
max_res <- max(vapply(res$outlier_genes, function(g) max(abs(g$std_residual), na.rm = TRUE), numeric(1)))
pairs <- vapply(res$outlier_cells, function(o) nrow(o$flagged), numeric(1))
save_figure(
  plot_cook(
    res$outlier_genes,
    res$outlier_cells,
    cells,
    titles = c(
      if (n_gene_flags == 0) {
        sprintf(
          "a. No gene is flagged in the weighted bulk fit: |standardised residual| stays below %.1f",
          max_res
        )
      } else {
        sprintf("a. %d gene flags in the weighted bulk fit", n_gene_flags)
      },
      sprintf(
        "b. %s to %s (gene, cell) pairs flagged per time, most per cell in the smallest types",
        fmt_n(min(pairs)),
        fmt_n(max(pairs))
      )
    )
  ),
  "fs_tier3_cook",
  fig_dir,
  width = 11,
  height = 10
)

# Transcript length -----------------------------------------------------------
lt <- res$length_test
save_figure(
  plot_platform_length(
    res$log_ratio,
    res$lengths,
    lt,
    title = sprintf(
      "The bulk / reference ratio grows with transcript length, slope %.2f to %.2f against a predicted 1",
      min(lt$slope),
      max(lt$slope)
    ),
    subtitle = sprintf(
      "Full-length bulk (Illumina Stranded mRNA) against 10x 3' v3; %s of %s panel genes have a length\nLength explains %.0f to %.0f %% of the variance of the log2 ratio (loess)",
      fmt_n(lt$mapped_genes[1L]),
      fmt_n(lt$panel_genes[1L]),
      100 * min(lt$loess_r2),
      100 * max(lt$loess_r2)
    )
  ),
  "fs_tier3_length",
  fig_dir,
  width = 11,
  height = 4.9
)

# A3 shift ------------------------------------------------------------------
share <- function(v) a3$share_abs_z_gt[a3$variant == v]
save_figure(
  plot_a3_shift(
    res$variants,
    level = 0.01,
    seed = 20261109L,
    title = sprintf(
      "A cross-fitted gene-wise ratio cuts the share of |z| > 2.33 from %s to %s; 2 %% is expected",
      pct_range(share("none")),
      pct_range(share("gene"))
    ),
    subtitle = "Distribution of z = D / sd(D) over panel genes and the three replicates, T4 standardisation; grey band: |z| <= 2.33"
  ),
  "fs_tier3_a3_shift",
  fig_dir,
  width = 11,
  height = 5
)

# Observed against fitted bulk ---------------------------------------------------
save_figure(
  plot_bulk_fit(
    res$variants,
    which = c("none", "gene"),
    title = "Divided by a cross-fitted gene-wise ratio, the bulk lies on the reference fit",
    subtitle = "Each bin counts gene x replicate values; dashed line: identity. Ratio for each time estimated on the two other times"
  ),
  "fs_tier3_bulk_fit",
  fig_dir,
  width = 10,
  height = 7.6
)

# Tail flags ------------------------------------------------------------------
t4 <- function(v) a3$flagged_T4[a3$variant == v]
save_figure(
  plot_tail_flags(
    res$variants,
    title = sprintf(
      "T4 flags %s genes per time under A3, %s with the gene-wise ratio",
      paste(fmt_n(range(t4("none"))), collapse = " to "),
      paste(fmt_n(range(t4("gene"))), collapse = " to ")
    ),
    subtitle = sprintf(
      "Genes flagged in at least one replicate, of %s panel genes; statements ranked by their count under A3",
      fmt_n(a3$n_genes[1L])
    )
  ),
  "fs_tier3_tail_flags",
  fig_dir,
  width = 11,
  height = 4.6
)

# Plug-in composition ---------------------------------------------------------
top_of <- function(v, tm) {
  cmp <- res$variants[[v]]$tails[[tm]]$plugin$composition
  m <- tapply(cmp$p_hat, cmp$cell_type, mean)
  list(type = names(which.max(m)), share = max(m))
}
last <- tail(names(res$variants$none$tails), 1L)
t_none <- top_of("none", last)
t_gene <- top_of("gene", last)
comp_title <- if (identical(t_none$type, t_gene$type)) {
  sprintf("At %s the largest plug-in type is the same with the gene-wise ratio", last)
} else {
  sprintf("At %s the largest plug-in type changes with the gene-wise ratio", last)
}
comp_subtitle <- sprintf(
  "%s: %s %.0f %% under A3, %s %.0f %% with the ratio\ndeMixSC weights (c = 2), one point per replicate; types ranked by mean share",
  last,
  t_none$type,
  100 * t_none$share,
  t_gene$type,
  100 * t_gene$share
)
save_figure(
  plot_plugin_composition(
    res$variants,
    which = c("none", "gene"),
    title = comp_title,
    subtitle = comp_subtitle
  ),
  "fs_tier3_plugin_composition",
  fig_dir,
  width = 11,
  height = 4.8
)

# Stability of the ratio across times -------------------------------------------
rho <- res$log_ratio$cor[upper.tri(res$log_ratio$cor)]
save_figure(
  plot_platform_ratio(
    res$log_ratio,
    res$drift,
    title = sprintf(
      "Each gene keeps its bulk / reference ratio across times (Spearman %.2f to %.2f)",
      min(rho),
      max(rho)
    ),
    subtitle = "Mean log2 ratio over the three replicates under A3; dashed line: identity; under A3 the points would gather at (0, 0)"
  ),
  "fs_tier3_drift",
  fig_dir,
  width = 11,
  height = 4.8
)
message("Wrote ", fig_dir)
