# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/01_02_plot_naively_selected_genes.R \
#   > "logs/01_02_naive_plots_$(date +%F).log" 2>&1 &
#
# Reads tables written by
# scripts/01_01_prepare_and_filter_genes_sctransform.R
# and writes figures under output/naive_marker_selection/.

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  file.path(getwd(), "scripts")
}
root <- normalizePath(file.path(script_dir, ".."))
out_dir <- file.path(root, "output", "naive_marker_selection")
table_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(name) {
  path <- file.path(table_dir, name)
  if (!file.exists(path)) {
    stop("Missing table ", path, ". Run 01_01 first.")
  }
  utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}

read_dict <- function(name) {
  utils::read.csv(
    file.path(root, "data", "dictionaries", name),
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
}

save_pages <- function(grobs, filename, width, height) {
  pages <- gridExtra::marrangeGrob(
    grobs = grobs,
    nrow = 1L,
    ncol = 1L,
    top = NULL
  )
  ggplot2::ggsave(
    filename = file.path(out_dir, filename),
    plot = pages,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    limitsize = FALSE
  )
}

pretty_time <- function(x) sub("h$", " h", x)

colour_dict <- read_dict("suppinger_celltype_colours.csv")
type_colours <- stats::setNames(
  colour_dict$colour,
  colour_dict$celltypeannotation
)

# ==========================================================================
# SECTION 1 · Density of mean log1p CPM before and after the SCT list ----
# ==========================================================================

gene_means <- read_table("gene_mean_log_cpm.csv")
time_levels <- c("48h", "72h", "96h")
density_grobs <- list()

for (time_level in time_levels) {
  slice <- gene_means[gene_means$time_point == time_level, , drop = FALSE]
  present <- colour_dict$celltypeannotation[
    colour_dict$celltypeannotation %in% unique(slice$cell_type)
  ]
  slice$cell_type <- factor(slice$cell_type, levels = present)
  xmax <- max(slice$mean_log_cpm, na.rm = TRUE)
  panels <- list(
    Before = slice,
    After = slice[slice$in_sct_filter %in% TRUE, , drop = FALSE]
  )
  plots <- lapply(names(panels), function(panel_name) {
    ggplot2::ggplot(
      panels[[panel_name]],
      ggplot2::aes(
        x = mean_log_cpm,
        y = cell_type,
        fill = cell_type,
        colour = cell_type
      )
    ) +
      ggridges::geom_density_ridges(
        alpha = 0.65,
        linewidth = 0.35,
        scale = 1.15,
        rel_min_height = 0.01
      ) +
      ggplot2::scale_fill_manual(values = type_colours, drop = TRUE) +
      ggplot2::scale_colour_manual(values = type_colours, drop = TRUE) +
      ggplot2::scale_x_continuous(limits = c(0, xmax * 1.02)) +
      ggplot2::labs(
        x = "Mean log1p CPM",
        y = NULL,
        fill = "Cell type",
        colour = "Cell type"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
  })
  names(plots) <- names(panels)
  legend <- cowplot::get_legend(
    plots$Before +
      ggplot2::guides(colour = "none") +
      ggplot2::theme(legend.position = "bottom")
  )
  body <- cowplot::plot_grid(
    plots$Before + ggplot2::theme(legend.position = "none"),
    plots$After + ggplot2::theme(legend.position = "none"),
    ncol = 2L,
    labels = c("Before", "After"),
    label_size = 12
  )
  page <- cowplot::plot_grid(
    body,
    legend,
    ncol = 1L,
    rel_heights = c(1, 0.22)
  )
  page <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label(
        paste0(
          pretty_time(time_level),
          ": gene-wise mean log1p CPM, all RNA genes (left) versus ",
          "the SCT counts allow-list (right)"
        ),
        fontface = "bold",
        size = 12,
        x = 0.01,
        hjust = 0
      ),
    page,
    ncol = 1L,
    rel_heights = c(0.06, 1)
  )
  density_grobs[[time_level]] <- cowplot::as_grob(page)
}

save_pages(
  grobs = density_grobs,
  filename = "density_log_cpm_before_after.pdf",
  width = 14,
  height = 9
)
message("Wrote density_log_cpm_before_after.pdf")

# ==========================================================================
# SECTION 2 · Marker Venn diagrams ----
# ==========================================================================

venn_sets <- read_table("venn_gene_sets.csv")
venn_list <- split(venn_sets$gene, venn_sets$set)
venn_pages <- list(
  list(
    sets = list(
      `RNA counts` = unique(venn_list$rna_counts),
      `Marker genes` = unique(venn_list$markers)
    ),
    title = "Raw RNA symbols and signalling markers"
  ),
  list(
    sets = list(
      `SCT pre-filter` = unique(venn_list$sct_prefilter),
      `Marker genes` = unique(venn_list$markers)
    ),
    title = "SCT counts allow-list and signalling markers"
  )
)
venn_grobs <- lapply(venn_pages, function(page) {
    plot <- ggVennDiagram::ggVennDiagram(page$sets, label_alpha = 0) +
      ggplot2::labs(title = page$title) +
      ggplot2::scale_fill_gradient(low = "#f7f7f7", high = "#74a9cf") +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(12, 16, 12, 36)
      )
    cowplot::as_grob(plot)
})
save_pages(
  grobs = venn_grobs,
  filename = "venn_markers_vs_counts.pdf",
  width = 8,
  height = 7
)
message("Wrote venn_markers_vs_counts.pdf")

# ==========================================================================
# SECTION 3 · Upset of strategy-2 genes across cell types ----
# ==========================================================================

per_type <- read_table("ranks_per_cell_type.csv")
union_hits <- per_type[per_type$in_union %in% TRUE, c("gene", "cell_type")]
union_hits <- unique(union_hits)

if (nrow(union_hits) > 0L) {
  union_hits$present <- TRUE
  membership <- tidyr::pivot_wider(
    union_hits,
    names_from = cell_type,
    values_from = present,
    values_fill = FALSE
  )
  membership_cols <- setdiff(colnames(membership), "gene")
  upset_plot <- ComplexUpset::upset(
    membership,
    intersect = membership_cols,
    min_size = 1,
    n_intersections = 20,
    name = "Genes"
  ) +
    ggplot2::ggtitle(
      paste(
        "Strategy 2 union across bulk-matched times:",
        "scTransform rank or signalling marker"
      )
    )
  ggplot2::ggsave(
    filename = file.path(out_dir, "upset_strategy2_celltypes.pdf"),
    plot = upset_plot,
    width = 16,
    height = 9,
    units = "in",
    dpi = 320,
    limitsize = FALSE
  )
  message("Wrote upset_strategy2_celltypes.pdf")
} else {
  message("No strategy-2 genes; skipping the upset plot.")
}

# ==========================================================================
# SECTION 4 · Stability of strategy 2 across time, within a cell type ----
# ==========================================================================

# No annotated type is present, with enough cells, at all three
# bulk-matched slots. Diagrams are drawn for every type scored in
# two or more of those slots.
scored <- per_type[!is.na(per_type$rank) | per_type$in_union %in% TRUE, ]
type_time <- unique(scored[, c("cell_type", "time_point")])
n_times <- table(type_time$cell_type)
stable_types <- names(n_times)[n_times >= 2L]
stability_grobs <- list()

for (type_name in stable_types) {
  sets <- lapply(time_levels, function(time_level) {
    rows <- per_type$cell_type == type_name &
      per_type$time_point == time_level &
      per_type$in_union %in% TRUE
    unique(per_type$gene[rows])
  })
  names(sets) <- vapply(time_levels, pretty_time, character(1))
  # Keep a slot only when that type was actually scored there.
  scored_here <- type_time$time_point[type_time$cell_type == type_name]
  sets <- sets[sub(" h$", "h", names(sets)) %in% scored_here]
  sets <- sets[lengths(sets) > 0L]
  if (length(sets) < 2L) {
    next
  }
  plot <- ggVennDiagram::ggVennDiagram(sets, label_alpha = 0) +
    ggplot2::labs(
      title = type_name,
      subtitle = "Strategy 2 genes (rank union markers) by time slot"
    ) +
    ggplot2::scale_fill_gradient(low = "#f7f7f7", high = "#74a9cf") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(12, 16, 12, 36)
    )
  stability_grobs[[type_name]] <- cowplot::as_grob(plot)
}

if (length(stability_grobs) > 0L) {
  save_pages(
    grobs = stability_grobs,
    filename = "venn_strategy2_across_time.pdf",
    width = 8,
    height = 7
  )
  message("Wrote venn_strategy2_across_time.pdf")
} else {
  message("No cell type is scored in two time slots; skipping stability Venns.")
}

# ==========================================================================
# SECTION 5 · Global between versus within stacked bars ----
# ==========================================================================

global <- read_table("ranks_global.csv")
global <- global[is.finite(global$between_over_total) & !is.na(global$rank), ]
bar_grobs <- list()

for (time_level in time_levels) {
  df <- global[global$time_point == time_level, , drop = FALSE]
  if (nrow(df) == 0L) {
    next
  }
  n_cut <- max(df$rank[df$selected %in% TRUE], na.rm = TRUE)
  if (!is.finite(n_cut)) {
    n_cut <- min(500L, max(df$rank))
  }
  cross <- df$rank[df$between_over_within < 1]
  cross <- if (length(cross) == 0L) NA_integer_ else min(cross)
  show_until <- n_cut
  if (!is.na(cross)) {
    show_until <- max(show_until, min(cross + 40L, max(df$rank)))
  }
  show_until <- min(show_until, 1200L, max(df$rank))
  df <- df[df$rank <= show_until, , drop = FALSE]
  long <- rbind(
    data.frame(
      rank = df$rank,
      gene = df$gene,
      fraction = df$between_over_total,
      component = "between",
      stringsAsFactors = FALSE
    ),
    data.frame(
      rank = df$rank,
      gene = df$gene,
      fraction = 1 - df$between_over_total,
      component = "within",
      stringsAsFactors = FALSE
    )
  )
  long$component <- factor(long$component, levels = c("between", "within"))
  cut_row <- df[df$rank == n_cut, , drop = FALSE][1, , drop = FALSE]
  top10 <- df[df$rank <= 10L, , drop = FALSE]
  top10$label <- paste0(top10$gene, "\n", top10$top_cell_type)
  green_right <- if (is.na(cross)) {
    show_until + 0.5
  } else {
    min(cross, show_until + 1L) - 0.5
  }
  formula_y <- min(cut_row$between_over_total + 0.16, 1.28)
  plot <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = rank, y = fraction, fill = component)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = 0.5,
      xmax = green_right,
      ymin = 0,
      ymax = 1.65,
      fill = "#2ca25f",
      alpha = 0.08
    )
  if (!is.na(cross) && cross <= show_until) {
    plot <- plot +
      ggplot2::annotate(
        "rect",
        xmin = green_right,
        xmax = show_until + 0.5,
        ymin = 0,
        ymax = 1.65,
        fill = "#e34a33",
        alpha = 0.08
      )
  }
  plot <- plot +
    ggplot2::geom_col(width = 0.9, colour = NA) +
    ggplot2::geom_hline(
      yintercept = 0.5,
      linetype = "dashed",
      colour = "#e34a33",
      linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = cut_row$between_over_total,
      linetype = "solid",
      colour = "#222222",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = n_cut,
      linetype = "solid",
      colour = "#222222",
      linewidth = 0.35
    ) +
    ggplot2::annotate(
      "label",
      x = n_cut,
      y = formula_y,
      label = sprintf(
        "frac(SS[between], SS[within]) == %.3g",
        cut_row$between_over_within
      ),
      parse = TRUE,
      hjust = 1,
      size = 3.2,
      fill = "white",
      linewidth = 0.2
    ) +
    ggplot2::annotate(
      "label",
      x = n_cut,
      y = min(formula_y + 0.14, 1.48),
      label = paste0(cut_row$gene, " (", cut_row$top_cell_type, ")"),
      hjust = 1,
      size = 3,
      fill = "white",
      linewidth = 0.2
    ) +
    ggplot2::annotate(
      "label",
      x = max(1, show_until * 0.72),
      y = 0.5,
      label = "SS[within] > SS[between]",
      parse = TRUE,
      size = 3,
      fill = "#fee0d2",
      vjust = -0.4
    ) +
    ggrepel::geom_text_repel(
      data = top10,
      ggplot2::aes(x = rank, y = 1, label = label),
      inherit.aes = FALSE,
      nudge_y = 0.38,
      direction = "y",
      min.segment.length = 0,
      size = 2.6,
      max.overlaps = Inf,
      seed = 1,
      box.padding = 0.15
    ) +
    ggplot2::scale_fill_manual(
      values = c(between = "#2ca25f", within = "#e34a33"),
      labels = c(
        between = "Between cell types",
        within = "Within cell types"
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(c(1, 10, 50, 100, 250, n_cut, show_until)))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.65),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      x = "Gene rank (decreasing between / within)",
      y = "Share of residual sum of squares",
      fill = NULL,
      title = paste0(
        pretty_time(time_level),
        ": global screen, between-type share of Pearson-residual SS"
      ),
      subtitle = paste0(
        "Vertical line at rank ",
        n_cut,
        " (",
        cut_row$gene,
        "). Green background: between exceeds within."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  bar_grobs[[time_level]] <- cowplot::as_grob(plot)
}

if (length(bar_grobs) > 0L) {
  save_pages(
    grobs = bar_grobs,
    filename = "bar_between_within_global.pdf",
    width = 16,
    height = 8
  )
  message("Wrote bar_between_within_global.pdf")
}

# ==========================================================================
# SECTION 6 · Cell and gene counts ----
# ==========================================================================

cell_counts <- read_table("cell_counts.csv")
show_df <- cell_counts[, c(
  "time_point",
  "cell_type",
  "n_cells",
  "n_genes_union"
)]
names(show_df) <- c(
  "time_point",
  "cell_type",
  "n_cells_single_cell",
  "n_sctransform_plus_marker_genes"
)
tab <- tinytable::tt(
  show_df,
  caption = paste(
    "Strategy 2 at bulk-matched times.",
    "Gene counts are zero when a type was not scored",
    "(artefact, or fewer than the minimum number of cells)."
  )
)
html_path <- file.path(out_dir, "strategy2_celltype_counts.html")
tinytable::save_tt(tab, output = html_path, overwrite = TRUE)
utils::write.csv(
  show_df,
  file = file.path(out_dir, "strategy2_celltype_counts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
message("Wrote strategy2_celltype_counts.html")
