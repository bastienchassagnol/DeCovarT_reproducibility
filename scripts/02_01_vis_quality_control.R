# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_01_vis_quality_control.R \
#   > "logs/02_01_vis_quality_control_$(date +%F).log" 2>&1 &
#
# Tier (i) figures from output/feature_selection/tier_01/scored_slots.rds:
# low-quality cells, gene filters and fractional counts, and PCA/PLS-DA
# per time on the Pearson residuals of M0 (genes passing the presence
# gate), with the 48 h projections also coloured by cell line. PDFs in
# output/feature_selection/tier_01/figures/, PNG copies in
# docs/images/feature_selection/.

root <- getwd()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_vis_quality_control.R"))

tier_dir <- file.path(root, "output", "feature_selection", "tier_01")
fig_dir <- file.path(tier_dir, "figures")
slots_path <- file.path(tier_dir, "scored_slots.rds")
if (!file.exists(slots_path)) {
  stop(
    "Missing ", slots_path,
    ". Run scripts/02_01_stats_quality_control.R first, or dvc pull its .dvc."
  )
}
prepared <- readRDS(slots_path)
colours <- celltype_colours(
  file.path(root, "data", "dictionaries", "suppinger_celltype_colours.csv")
)
gene_summary <- utils::read.csv(
  file.path(tier_dir, "gene_filter_summary.csv"),
  stringsAsFactors = FALSE
)

# Low-quality cells (all bulk-matched cells, artefacts included).
save_figure(
  plot_low_quality_cells(prepared$cell_values, colours),
  "fs_tier1_low_quality_cells",
  fig_dir,
  width = 11,
  height = 8.5
)

# Gene filters and fractional counts.
save_figure(
  plot_gene_filters(gene_summary, prepared$fractional, prepared$gene_rules),
  "fs_tier1_gene_filters",
  fig_dir,
  width = 11,
  height = 8
)

# Projections: every time slot by cell type; 48 h also by cell line
# (batch 1 = B-S, batch 2 = SBR, the two lines pooled there).
line_label <- c(batch1 = "B-S (batch 1)", batch2 = "SBR (batch 2)")
line_colours <- c("B-S (batch 1)" = "#1d4e89", "SBR (batch 2)" = "#c45c26")
projections <- lapply(prepared$time_levels, function(tm) {
  slot <- prepared$slots[[tm]]
  gate <- prepared$presence[[tm]]
  genes <- gate$gene[gate$pass]
  message(tm, ": projections on ", length(genes), " genes")
  groups <- if (length(unique(slot$batch)) > 1L) {
    list("Cell line" = unname(line_label[slot$batch]))
  } else {
    NULL
  }
  drawn <- plot_type_pca_plsda(
    counts = slot$counts[genes, , drop = FALSE],
    cell_type = slot$cell_type,
    size_factors = slot$library_size / stats::median(slot$library_size),
    colours = colours,
    input = "pearson",
    groups = groups,
    group_colours = list("Cell line" = line_colours)
  )
  # Titles state the finding (#33); the PCA panel drops its legend, which
  # the PLS-DA panel of the same row repeats (#41).
  pca_share <- sum(drawn$pca_explained)
  drawn$pca <- drawn$pca +
    ggplot2::labs(
      title = sprintf("%s PCA", tm),
      subtitle = sprintf("PC1 + PC2: %.1f %% of the residual variance", pca_share)
    ) +
    ggplot2::theme(legend.position = "none")
  drawn$plsda <- drawn$plsda +
    ggplot2::labs(
      title = sprintf("%s PLS-DA", tm),
      subtitle = "Two components fitted to separate the types"
    ) +
    ggplot2::theme(legend.position = "right")
  if (!is.null(groups)) {
    line <- groups[["Cell line"]]
    pos <- tapply(drawn$scores$pc1 > 0, line, mean)
    drawn[["pca_Cell line"]] <- drawn[["pca_Cell line"]] +
      ggplot2::labs(
        title = sprintf("%s PCA: PC1 splits the cell lines", tm),
        subtitle = sprintf(
          "PC1 > 0: %.0f %% of B-S, %.0f %% of SBR cells",
          100 * pos[["B-S (batch 1)"]],
          100 * pos[["SBR (batch 2)"]]
        )
      ) +
      ggplot2::theme(legend.position = "none")
    drawn[["plsda_Cell line"]] <- drawn[["plsda_Cell line"]] +
      ggplot2::labs(
        title = sprintf("%s PLS-DA also splits the lines", tm),
        subtitle = "Fitted on the types; coloured by cell line"
      ) +
      ggplot2::theme(legend.position = "right")
  }
  drawn$time_point <- tm
  drawn
})
names(projections) <- prepared$time_levels
panels <- list()
for (tm in names(projections)) {
  d <- projections[[tm]]
  panels <- c(panels, list(d$pca, d$plsda))
  if (!is.null(d[["pca_Cell line"]])) {
    panels <- c(panels, list(d[["pca_Cell line"]], d[["plsda_Cell line"]]))
  }
}
save_figure(
  patchwork::wrap_plots(panels, ncol = 2, byrow = TRUE) &
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9.5)),
  "fs_tier1_pca_plsda",
  fig_dir,
  width = 11,
  height = 4 * length(panels) / 2
)

# Share of each cell line on either side of PC1, per type and time
# (summary cited by the chapter).
batch_split <- do.call(rbind, lapply(names(projections), function(tm) {
  sc <- projections[[tm]]$scores
  slot <- prepared$slots[[tm]]
  sc$batch <- slot$batch
  sc$pc1_positive <- sc$pc1 > 0
  tab <- as.data.frame(
    table(cell_type = sc$cell_type, batch = sc$batch, pc1_positive = sc$pc1_positive),
    stringsAsFactors = FALSE
  )
  tab$time_point <- tm
  tab[tab$Freq > 0, , drop = FALSE]
}))
utils::write.csv(
  batch_split,
  file.path(tier_dir, "projection_pc1_by_batch.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
message("Wrote ", fig_dir)
