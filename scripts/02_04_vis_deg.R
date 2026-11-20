# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_04_vis_deg.R \
#   > "logs/02_04_vis_deg_$(date +%F).log" 2>&1 &
#
# Figures of the post hoc contrasts (appendix of chapter 2): volcano,
# log2FC distributions, confirmed genes per family and design, UpSet of
# the DEG sets, heatmap of the top confirmed genes, within-type
# bimodality, and the old (SCTransform) against the new panel. PDFs and
# data in output/feature_selection/deg/figures/, PNGs in
# docs/images/feature_selection/.

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_04_vis_deg.R"))

deg_dir <- file.path(root, "output", "feature_selection", "deg")
fig_dir <- file.path(deg_dir, "figures")
png_dir <- file.path(root, "docs", "images", "feature_selection")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
deg <- readRDS(file.path(deg_dir, "deg.rds"))
prepared <- readRDS(file.path(root, "output", "feature_selection", "tier_01", "scored_slots.rds"))
cc <- readRDS(file.path(root, "output", "feature_selection", "tier_02", "cross_cell.rds"))
type_cols <- celltype_colours(file.path(root, "data", "dictionaries", "suppinger_celltype_colours.csv"))
primary <- "welch_perm"
dec <- deg$decisions
prim <- dec[dec$test == primary, , drop = FALSE]
prim$status <- deg_status(prim)
seed <- 20261122L

# 1. Volcano.
n_ab <- sum(prim$b_graph)
n_a <- sum(prim$a_graph)
save_figure(
  plot_deg_volcano(
    dec,
    colour_test = primary,
    title = sprintf("%s gene x type pairs are confirmed above the average of the other types; %s above every type", fmt_n(n_a), fmt_n(n_ab)),
    subtitle = "Colour: permuted Welch family, stage-wise graphical procedure; y: one-sided Wald p (NB fit)"
  ),
  "fs_deg_volcano",
  fig_dir,
  png_dir = png_dir,
  data = dec[dec$test %in% c("lr", primary), c("time_point", "test", "gene", "type", "log2fc", "p_a", "a_graph", "b_graph")],
  width = 12,
  height = 8.5
)

# 2. log2FC distributions.
save_figure(
  plot_lfc_distributions(
    deg$results,
    type_cols,
    title = "One-versus-average fold changes are centred near 0 with long right tails",
    subtitle = "Density over all tested genes, with the middle 66 % and 95 % intervals"
  ),
  "fs_deg_lfc_distributions",
  fig_dir,
  png_dir = png_dir,
  data = lapply(deg$results, function(r) r$log2fc),
  width = 11,
  height = 4.8
)

# 3. Confirmed genes per family and design.
save_figure(
  plot_deg_counts(
    deg$counts,
    title = "The four families agree on which types carry most genes; B confirms far fewer than A",
    subtitle = "Confirmed gene x type pairs at 0.05 (flat BH: share of false pairs; stage-wise: gene-level OFDR)"
  ),
  "fs_deg_counts",
  fig_dir,
  png_dir = png_dir,
  data = deg$counts,
  width = 12,
  height = 5
)

# 4. UpSet of the DEG sets (primary family, A confirmed).
sets <- split(prim$gene[prim$a_graph], paste0(prim$time_point[prim$a_graph], ": ", prim$type[prim$a_graph]))
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
    title = "Most enriched genes are specific to one type and time",
    subtitle = "Genes confirmed above the average of the other types (permuted Welch, graphical procedure); 25 largest intersections",
    theme = theme_slide()
  )
save_figure(p_up, "fs_deg_upset", fig_dir, png_dir = png_dir, data = up, width = 12, height = 6.5)

# 5. Heatmap of the top confirmed genes per type, cells subsampled.
heat_one <- function(tm, k = 12L, max_cells = 150L) {
  slot <- prepared$slots[[tm]]
  p <- prim[prim$time_point == tm & prim$a_graph, , drop = FALSE]
  p <- p[order(p$type, !p$b_graph, -p$log2fc), , drop = FALSE]
  top <- do.call(rbind, lapply(split(p, p$type), utils::head, k))
  top <- top[!duplicated(top$gene), , drop = FALSE]
  ct <- as.character(slot$cell_type)
  keep <- unlist(withr::with_seed(seed, lapply(split(seq_along(ct), ct), function(i) {
    if (length(i) > max_cells) sample(i, max_cells) else i
  })))
  x <- as.matrix(slot$counts[top$gene, keep, drop = FALSE])
  v <- log1p(sweep(x, 2L, cc$sf[[tm]][keep], "/"))
  z <- t(scale(t(v)))
  z[!is.finite(z)] <- 0
  z <- pmin(pmax(z, -3), 3)
  types <- intersect(names(type_cols), unique(ct))
  col_split <- factor(ct[keep], levels = types)
  row_split <- factor(top$type, levels = types)
  ha <- ComplexHeatmap::HeatmapAnnotation(
    type = col_split,
    batch = factor(slot$batch[keep]),
    col = list(type = type_cols[types], batch = c(batch1 = "#d9d3c7", batch2 = "#3d4a42")),
    show_annotation_name = FALSE,
    annotation_legend_param = list(type = list(title = "Cell type"), batch = list(title = "Batch"))
  )
  ComplexHeatmap::Heatmap(
    z,
    name = "z",
    col = circlize::colorRamp2(c(-3, 0, 3), c("#1d4e89", "#f7f4ee", "#c45c26")),
    top_annotation = ha,
    column_split = col_split,
    row_split = row_split,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 6, fontface = "italic"),
    column_title_gp = grid::gpar(fontsize = 8),
    row_title_gp = grid::gpar(fontsize = 7),
    row_title_rot = 0,
    # Types are named by the colour annotation; slice titles of small
    # types would overlap.
    column_title = NULL,
    use_raster = TRUE,
    heatmap_legend_param = list(title = "z-score of\nlog(1 + x/s)")
  )
}
heat_files <- list()
for (tm in names(deg$results)) {
  hm <- heat_one(tm)
  nm <- paste0("fs_deg_heatmap_", tm)
  for (dev in c("pdf", "png")) {
    f <- if (dev == "pdf") file.path(fig_dir, paste0(nm, ".pdf")) else file.path(png_dir, paste0(nm, ".png"))
    if (dev == "pdf") grDevices::cairo_pdf(f, width = 13, height = 7.5) else ragg::agg_png(f, width = 13, height = 7.5, units = "in", res = 300, background = "white")
    ComplexHeatmap::draw(
      hm,
      column_title = sprintf("Top confirmed genes per type at %s (up to 12; above every type first), %s", tm, "cells subsampled to 150 per type"),
      column_title_gp = grid::gpar(fontsize = 11, fontface = "bold", col = "#0e4d3a"),
      merge_legend = TRUE
    )
    grDevices::dev.off()
  }
  heat_files[[tm]] <- nm
}
saveRDS(list(genes = lapply(names(deg$results), function(tm) prim[prim$time_point == tm & prim$a_graph, c("gene", "type", "log2fc", "b_graph")]), seed = seed),
        file.path(fig_dir, "fs_deg_heatmap.rds"))

# 6. Within-type bimodality of the genes confirmed above every type.
# Two cues of sub-states (CardamomOT basins: an NB mixture within a type):
# (a) excess zeros, the share of cells at zero minus the share the fitted
# NB (type mean, gene-wise phi of M1, cell size factors) predicts, an
# off-state; (b) the bimodality coefficient and dip test of the non-zero
# cells, two expressed modes.
bim <- list()
for (tm in names(deg$results)) {
  slot <- prepared$slots[[tm]]
  ct <- as.character(slot$cell_type)
  res <- deg$results[[tm]]
  phi <- stats::setNames(cc$tests[[tm]]$overdispersion, cc$tests[[tm]]$gene)
  p <- prim[prim$time_point == tm & prim$b_graph, , drop = FALSE]
  for (tp in unique(p$type)) {
    g <- p$gene[p$type == tp]
    cells <- which(ct == tp)
    sf <- cc$sf[[tm]][cells]
    x <- as.matrix(slot$counts[g, cells, drop = FALSE])
    mu <- exp(res$log_mean[g, tp])
    p0 <- vapply(seq_along(g), function(i) {
      m <- sf * mu[i]
      ph <- phi[[g[i]]]
      mean(if (ph > 0) (1 + ph * m)^(-1 / ph) else exp(-m))
    }, numeric(1))
    nz <- lapply(seq_along(g), function(i) log1p(x[i, x[i, ] > 0] / sf[x[i, ] > 0]))
    bim[[paste(tm, tp)]] <- data.frame(
      time_point = tm,
      type = tp,
      gene = g,
      mean_count = rowMeans(sweep(x, 2L, sf, "/")),
      share_zero = rowMeans(x == 0),
      expected_zero = p0,
      excess_zero = rowMeans(x == 0) - p0,
      bc_nonzero = vapply(nz, bimodality_coefficient, numeric(1)),
      dip_p_nonzero = vapply(nz, function(v) if (length(v) > 3L && stats::sd(v) > 0) diptest::dip.test(v)$p.value else NA_real_, numeric(1))
    )
  }
}
bim <- do.call(rbind, bim)
utils::write.csv(bim, file.path(deg_dir, "bimodality.csv"), row.names = FALSE, fileEncoding = "UTF-8")
# The six genes with the largest excess of zeros among those with a mean
# above one count per cell (where the NB predicts almost no zero).
cand <- bim[bim$mean_count > 1 & !is.na(bim$excess_zero), , drop = FALSE]
top6 <- utils::head(cand[order(-cand$excess_zero), , drop = FALSE], 6L)
expr <- do.call(rbind, lapply(seq_len(nrow(top6)), function(i) {
  slot <- prepared$slots[[top6$time_point[i]]]
  cells <- which(as.character(slot$cell_type) == top6$type[i])
  data.frame(
    time_point = top6$time_point[i],
    type = top6$type[i],
    gene = top6$gene[i],
    value = log1p(as.numeric(slot$counts[top6$gene[i], cells]) / cc$sf[[top6$time_point[i]]][cells])
  )
}))
n_off <- sum(bim$excess_zero > 0.1, na.rm = TRUE)
n_two <- sum(bim$bc_nonzero > 5 / 9 & bim$dip_p_nonzero < 0.05, na.rm = TRUE)
save_figure(
  plot_bimodality(
    top6,
    expr,
    title = sprintf("Within-type sub-states are rare: %s / %s genes above every type exceed the NB zero share by 10 points", fmt_n(n_off), fmt_n(nrow(bim))),
    subtitle = sprintf("Six largest excesses among genes with a mean above 1; two modes among non-zero cells (BC > 5/9, dip p < 0.05): %s genes", fmt_n(n_two))
  ),
  "fs_deg_bimodality",
  fig_dir,
  png_dir = png_dir,
  data = list(bimodality = bim, top6 = top6, expr = expr),
  width = 11,
  height = 6.4
)

# 7. Old (SCTransform top 500) against new (BH) panel.
old_path <- "output/naive_marker_selection/tables/ranks_global_equal_type.csv"
old <- utils::read.csv(text = system2("git", c("-C", root, "show", paste0("3b71ae0^:", old_path)), stdout = TRUE), stringsAsFactors = FALSE)
old <- old[old$selected %in% TRUE, , drop = FALSE]
cmp <- list()
ranks <- list()
for (tm in names(cc$ranked)) {
  r <- cc$ranked[[tm]]
  new <- r[r$selected %in% TRUE, , drop = FALSE]
  o <- old[old$time_point == tm, , drop = FALSE]
  both <- intersect(o$gene, new$gene)
  cmp[[tm]] <- data.frame(
    time_point = tm,
    class = c("old only", "both", "new only"),
    n = c(length(setdiff(o$gene, new$gene)), length(both), length(setdiff(new$gene, o$gene)))
  )
  ranks[[tm]] <- data.frame(
    time_point = tm,
    gene = both,
    old_rank = o$rank[match(both, o$gene)],
    new_rank = new$rank[match(both, new$gene)]
  )
}
cmp <- do.call(rbind, cmp)
ranks <- do.call(rbind, ranks)
utils::write.csv(cmp, file.path(deg_dir, "panel_comparison.csv"), row.names = FALSE, fileEncoding = "UTF-8")
shared <- cmp$n[cmp$class == "both"]
save_figure(
  plot_panel_comparison(
    cmp,
    ranks,
    title = sprintf("%s of the 500 SCTransform genes per time pass the new gate", paste(fmt_n(shared), collapse = " / ")),
    subtitle = "Old: top 500 by type-equal between/within ratio of SCTransform residuals (archived); new: BH 0.05 on the permuted Welch G"
  ),
  "fs_deg_panel_comparison",
  fig_dir,
  png_dir = png_dir,
  data = list(comparison = cmp, ranks = ranks),
  width = 11,
  height = 6
)
message("Wrote ", fig_dir)
