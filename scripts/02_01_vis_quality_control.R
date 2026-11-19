# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_01_vis_quality_control.R \
#   > "logs/02_01_vis_quality_control_$(date +%F).log" 2>&1 &
#
# Tier (i) figures from output/feature_selection/tier_01/scored_slots.rds:
# cells kept, gene filters, fractional counts, markers against the SCT
# genes, expression before/after the gate, and the projections (Pearson
# residuals of M0 in the chapter; Pearson and the authors' SCT residuals,
# with UMAP, in the appendix). PDFs and their regeneration RDS in
# output/feature_selection/tier_01/figures/, PNG copies in
# docs/images/feature_selection/.

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_vis_quality_control.R"))

tier_dir <- file.path(root, "output", "feature_selection", "tier_01")
fig_dir <- file.path(tier_dir, "figures")
png_dir <- file.path(root, "docs", "images", "feature_selection")
slots_path <- file.path(tier_dir, "scored_slots.rds")
if (!file.exists(slots_path)) {
  stop(
    "Missing ",
    slots_path,
    ". Run scripts/02_01_stats_quality_control.R first, or dvc pull its .dvc."
  )
}
prepared <- readRDS(slots_path)
colour_path <- file.path(root, "data", "dictionaries", "suppinger_celltype_colours.csv")
colours <- celltype_colours(colour_path)
lineages <- celltype_lineages(colour_path)
gene_summary <- utils::read.csv(
  file.path(tier_dir, "gene_filter_summary.csv"),
  stringsAsFactors = FALSE
)
save_fig <- function(plot, name, width, height, data) {
  save_figure(plot, name, fig_dir, png_dir, width = width, height = height, data = data)
}

# Cells kept by the authors (artefacts included).
cell_values <- prepared$cell_values[, c("time_point", "cell_type", "library_size", "detected_genes")]
save_fig(
  plot_low_quality_cells(cell_values, colours, lineages),
  "fs_tier1_low_quality_cells",
  width = 11,
  height = 8.5,
  data = list(cell_values = cell_values, colours = colours, lineages = lineages)
)

# Gene filters.
gene_rules <- prepared$gene_rules[, c("gene", "time_point", "enters_tier_ii", "in_sct_layer")]
save_fig(
  plot_gene_filters(gene_summary, gene_rules),
  "fs_tier1_gene_filters",
  width = 11,
  height = 6,
  data = list(gene_summary = gene_summary, gene_rules = gene_rules)
)

# Fractional counts by type and biotype.
frac_data <- fractional_figure_data(prepared$fractional_group, prepared$biotype)
save_fig(
  plot_fractional_genes(frac_data, colours, lineages),
  "fs_tier1_fractional",
  width = 11,
  height = 12,
  data = c(frac_data, list(colours = colours, lineages = lineages))
)

# Curated markers against the authors' SCT genes.
venn <- plot_marker_venn(
  prepared$symbol_map$mouse_symbol,
  prepared$sct_genes,
  universe = prepared$fractional$gene
)
save_fig(venn$plot, "fs_tier1_marker_venn", width = 7, height = 5, data = venn$data)

# Per-gene mean expression by type, before and after the gate.
means <- do.call(
  rbind,
  lapply(prepared$time_levels, function(tm) {
    gate <- prepared$presence[[tm]]
    gene_mean_by_type(prepared$slots[[tm]], gate$gene[gate$pass])
  })
)
save_fig(
  plot_density_before_after(means, colours),
  "fs_tier1_density_filter",
  width = 11,
  height = 7,
  data = list(means = means, colours = colours)
)

# Projections. Chapter: Pearson residuals, PCA and sPLS-DA per time, 48 h
# also by cell line (batch 1 = B-S, batch 2 = SBR). Appendix: PCA, UMAP
# and sPLS-DA on the Pearson and on the authors' SCT residuals.
line_label <- c(batch1 = "B-S (batch 1)", batch2 = "SBR (batch 2)")
line_colours <- c("B-S (batch 1)" = "#1d4e89", "SBR (batch 2)" = "#c45c26")
project <- function(tm, input) {
  slot <- prepared$slots[[tm]]
  if (input == "sct") {
    x <- prepared$sct_slots[[tm]]
  } else {
    gate <- prepared$presence[[tm]]
    x <- slot$counts[gate$gene[gate$pass], , drop = FALSE]
  }
  groups <- if (length(unique(slot$batch)) > 1L) {
    list("Cell line" = unname(line_label[slot$batch]))
  } else {
    NULL
  }
  message(tm, " / ", input, ": projections on ", nrow(x), " genes")
  drawn <- plot_type_pca_plsda(
    counts = x,
    cell_type = slot$cell_type,
    size_factors = slot$library_size / stats::median(slot$library_size),
    colours = colours,
    lineages = lineages,
    input = input,
    groups = groups,
    group_colours = list("Cell line" = line_colours)
  )
  drawn$scores$time_point <- tm
  drawn
}
title_panels <- function(d, tm) {
  titles <- c(pca = "PCA", umap = "UMAP", plsda = "sPLS-DA")
  for (m in names(titles)) {
    d[[m]] <- d[[m]] + ggplot2::labs(title = paste(tm, titles[[m]]))
  }
  d
}
inputs <- c("pearson", if (!is.null(prepared$sct_slots)) "sct")
if (is.null(prepared$sct_slots)) {
  message("No `sct_slots` in ", slots_path, ": SCT projections skipped.")
}
proj <- lapply(inputs, function(input) {
  out <- lapply(prepared$time_levels, function(tm) title_panels(project(tm, input), tm))
  names(out) <- prepared$time_levels
  out
})
names(proj) <- inputs
proj_data <- function(input) {
  lapply(proj[[input]], function(d) {
    d[c("scores", "pca_explained", "plsda_explained", "n_genes", "settings")]
  })
}

# Chapter figure: Pearson residuals, PCA and sPLS-DA.
panels <- list()
for (tm in prepared$time_levels) {
  d <- proj$pearson[[tm]]
  panels <- c(
    panels,
    list(
      d$pca +
        ggplot2::labs(
          subtitle = sprintf("PC1 + PC2: %.1f %% of the variance", sum(d$pca_explained))
        ) +
        ggplot2::theme(legend.position = "none"),
      d$plsda + ggplot2::theme(legend.position = "right")
    )
  )
  if (!is.null(d[["pca_Cell line"]])) {
    pos <- tapply(d$scores$pc1 > 0, d$scores[["Cell line"]], mean)
    panels <- c(
      panels,
      list(
        d[["pca_Cell line"]] +
          ggplot2::labs(
            title = sprintf("%s PCA: PC1 splits the cell lines", tm),
            subtitle = sprintf(
              "PC1 > 0: %.0f %% of B-S, %.0f %% of SBR cells",
              100 * pos[["B-S (batch 1)"]],
              100 * pos[["SBR (batch 2)"]]
            )
          ) +
          ggplot2::theme(legend.position = "none"),
        d[["plsda_Cell line"]] +
          ggplot2::labs(title = sprintf("%s sPLS-DA also splits the lines", tm)) +
          ggplot2::theme(legend.position = "right")
      )
    )
  }
}
save_fig(
  patchwork::wrap_plots(panels, ncol = 2, byrow = TRUE) &
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9.5)),
  "fs_tier1_pca_plsda",
  width = 11,
  height = 4 * length(panels) / 2,
  data = list(pearson = proj_data("pearson"), colours = colours, lineages = lineages)
)

# Appendix figures: PCA, UMAP, sPLS-DA per time, one figure per input.
for (input in inputs) {
  grid_panels <- unlist(
    lapply(prepared$time_levels, function(tm) {
      d <- proj[[input]][[tm]]
      list(
        d$pca + ggplot2::theme(legend.position = "none"),
        d$umap + ggplot2::theme(legend.position = "none"),
        d$plsda + ggplot2::theme(legend.position = "right")
      )
    }),
    recursive = FALSE
  )
  save_fig(
    patchwork::wrap_plots(grid_panels, ncol = 3, byrow = TRUE) &
      ggplot2::theme(plot.title = ggplot2::element_text(size = 9.5)),
    paste0("fs_tier1_projections_", input),
    width = 15,
    height = 4.5 * length(prepared$time_levels),
    data = list(scores = proj_data(input), colours = colours, lineages = lineages)
  )
}

# Appendix: one UMAP of the three slots on the authors' SCT residuals
# (same genes and model at every time), one panel per highlighted time.
if ("sct" %in% inputs) {
  pooled_x <- t(do.call(cbind, prepared$sct_slots))
  pooled_type <- unlist(lapply(prepared$slots, `[[`, "cell_type"), use.names = FALSE)
  pooled <- project_slot(pooled_x, pooled_type)
  pooled$scores$time_point <- rep(
    prepared$time_levels,
    vapply(prepared$sct_slots, ncol, integer(1))
  )
  save_fig(
    plot_pooled_umap(pooled$scores, colours, lineages),
    "fs_tier1_umap_pooled",
    width = 13,
    height = 6,
    data = list(scores = pooled$scores, settings = pooled$settings, colours = colours, lineages = lineages)
  )
}

# Share of each cell line on either side of PC1, per type and time
# (summary cited by the chapter).
batch_split <- do.call(
  rbind,
  lapply(prepared$time_levels, function(tm) {
    sc <- proj$pearson[[tm]]$scores
    tab <- as.data.frame(
      table(
        cell_type = sc$cell_type,
        batch = prepared$slots[[tm]]$batch,
        pc1_positive = sc$pc1 > 0
      ),
      stringsAsFactors = FALSE
    )
    tab$time_point <- tm
    tab[tab$Freq > 0, , drop = FALSE]
  })
)
utils::write.csv(
  batch_split,
  file.path(tier_dir, "projection_pc1_by_batch.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
message("Wrote ", fig_dir)
