# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_04_vis_deg_twosided.R \
#   > "logs/02_04_vis_deg_twosided_$(date +%F).log" 2>&1 &
#
# Figures of the two-sided post hoc contrasts (A2: the type differs from
# the geometric mean of the other types; B2: it is above, or below, every
# other type), appendix of chapter 2: volcano with the fold-change guides,
# confirmed genes per family, design and direction, UpSet of the DEG
# sets, heatmaps of the top up and down genes. PDFs and data in
# output/feature_selection/deg_twosided/figures/, PNGs in
# docs/images/feature_selection/ (prefix fs_deg_twosided_).

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_04_vis_deg.R"))

deg_dir <- file.path(root, "output", "feature_selection", "deg_twosided")
fig_dir <- file.path(deg_dir, "figures")
png_dir <- file.path(root, "docs", "images", "feature_selection")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
deg <- readRDS(file.path(deg_dir, "deg_twosided.rds"))
prepared <- readRDS(file.path(root, "output", "feature_selection", "tier_01", "scored_slots.rds"))
cc <- readRDS(file.path(root, "output", "feature_selection", "tier_02", "cross_cell.rds"))
type_cols <- celltype_colours(file.path(root, "data", "dictionaries", "suppinger_celltype_colours.csv"))
primary <- "welch_perm"
lfc_threshold <- deg$lfc_thresholds[[1L]]
dec <- deg$decisions
prim <- dec[dec$test == primary, , drop = FALSE]
times <- unique(prim$time_point)
seed <- 20261122L

# 1. Volcano.
n_up <- sum(prim$a_graph & prim$direction %in% "up")
n_down <- sum(prim$a_graph & prim$direction %in% "down")
save_figure(
  plot_deg_volcano(
    dec,
    colour_test = primary,
    lfc_threshold = lfc_threshold,
    title = sprintf("Two-sided: %s gene x type pairs above and %s below the average of the other types", fmt_n(n_up), fmt_n(n_down)),
    subtitle = sprintf("Colour: permuted Welch family, stage-wise graphical procedure; y: two-sided Wald p (NB fit); dashed: |log2FC| = %g (drawn, not tested)", lfc_threshold)
  ),
  "fs_deg_twosided_volcano",
  fig_dir,
  png_dir = png_dir,
  data = list(
    decisions = dec[dec$test %in% c("lr", primary), c("time_point", "test", "gene", "type", "log2fc", "direction", "p_a", "a_graph", "b_graph")],
    lfc_threshold = lfc_threshold
  ),
  width = 12,
  height = 8.5
)

# 2. Confirmed genes per family, design and direction.
save_figure(
  plot_deg_counts(
    deg$counts,
    title = "Two-sided contrasts add the genes below the average of the other types",
    subtitle = "Confirmed gene x type pairs at 0.05 per direction (flat BH: share of false pairs; stage-wise: gene-level OFDR)"
  ),
  "fs_deg_twosided_counts",
  fig_dir,
  png_dir = png_dir,
  data = deg$counts,
  width = 12,
  height = 7.5
)

# 3. UpSet of the DEG sets (primary family, A2 confirmed), one set per
# time, type and direction.
conf <- prim[prim$a_graph, , drop = FALSE]
sets <- split(conf$gene, paste0(conf$time_point, ": ", conf$type, ", ", conf$direction))
all_genes <- sort(unique(unlist(sets)))
up <- as.data.frame(lapply(sets, function(s) all_genes %in% s), check.names = FALSE)
up$gene <- all_genes
p_up <- ComplexUpset::upset(
  up,
  intersect = names(sets),
  min_size = 5,
  n_intersections = 25,
  width_ratio = 0.25,
  base_annotations = list("Genes" = ComplexUpset::intersection_size(fill = "#1d4e89", text = list(size = 2.5))),
  set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(fill = "#9a9286"))
) +
  patchwork::plot_annotation(
    title = "Two-sided: most confirmed genes belong to one type, time and direction",
    subtitle = "Genes confirmed against the average of the other types (A2; permuted Welch, graphical procedure); 25 largest intersections",
    theme = theme_slide()
  )
save_figure(p_up, "fs_deg_twosided_upset", fig_dir, png_dir = png_dir, data = up, width = 12, height = 8)

# 4. Heatmaps: per type, up to 6 genes up and 6 down (against every type
# first, then by |log2FC|).
tops <- list()
for (tm in times) {
  p <- prim[prim$time_point == tm & prim$a_graph, , drop = FALSE]
  p <- p[order(p$type, p$direction != "up", !p$b_graph, -abs(p$log2fc)), , drop = FALSE]
  top <- do.call(rbind, lapply(split(p, list(p$type, p$direction), drop = TRUE), utils::head, 6L))
  top <- top[order(top$type, top$direction != "up"), , drop = FALSE]
  top <- top[!duplicated(top$gene), , drop = FALSE]
  tops[[tm]] <- top[, c("gene", "type", "direction", "log2fc", "b_graph")]
  hm <- deg_heatmap(prepared$slots[[tm]], cc$sf[[tm]], top, type_cols, seed)
  nm <- paste0("fs_deg_twosided_heatmap_", tm)
  for (dev in c("pdf", "png")) {
    f <- if (dev == "pdf") file.path(fig_dir, paste0(nm, ".pdf")) else file.path(png_dir, paste0(nm, ".png"))
    if (dev == "pdf") grDevices::cairo_pdf(f, width = 13, height = 8.5) else ragg::agg_png(f, width = 13, height = 8.5, units = "in", res = 300, background = "white")
    ComplexHeatmap::draw(
      hm,
      column_title = sprintf("Two-sided: top confirmed genes per type at %s (up to 6 up, then 6 down), cells subsampled to 150 per type", tm),
      column_title_gp = grid::gpar(fontsize = 11, fontface = "bold", col = "#0e4d3a"),
      merge_legend = TRUE
    )
    grDevices::dev.off()
  }
}
saveRDS(list(genes = tops, seed = seed), file.path(fig_dir, "fs_deg_twosided_heatmap.rds"))
message("Wrote ", fig_dir)
