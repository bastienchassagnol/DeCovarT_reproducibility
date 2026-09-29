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
    panel <- panels[[panel_name]]
    null_share <- stats::aggregate(
      mean_log_cpm ~ cell_type,
      data = panel,
      FUN = function(z) mean(z <= 1e-8)
    )
    names(null_share)[2] <- "null_share"
    null_share$label <- sprintf("%d%% null", round(100 * null_share$null_share))
    ggplot2::ggplot(
      panel,
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
      ggplot2::geom_rug(
        sides = "b",
        alpha = 0.35,
        linewidth = 0.2,
        length = grid::unit(0.025, "npc")
      ) +
      ggplot2::geom_label(
        data = null_share,
        ggplot2::aes(x = xmax * 0.78, y = cell_type, label = label),
        inherit.aes = FALSE,
        size = 2.8,
        fill = "white",
        linewidth = 0.2,
        label.padding = ggplot2::unit(0.12, "lines")
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
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.y = ggplot2::element_text(face = "bold", size = 12)
      )
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
        x = 0.5,
        hjust = 0.5
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
sct_genes <- unique(venn_list$sct_prefilter)
marker_genes <- unique(venn_list$markers)
rna_genes <- unique(venn_list$rna_counts)
in_both <- intersect(marker_genes, sct_genes)
marker_only <- setdiff(marker_genes, sct_genes)
n_rna <- length(rna_genes)
n_sct <- length(sct_genes)
n_marker <- length(marker_genes)
n_both <- length(in_both)
# P(X >= n_both) under a uniform draw of the marker panel from the RNA universe.
hyper_p <- stats::phyper(
  n_both - 1L,
  n_sct,
  n_rna - n_sct,
  n_marker,
  lower.tail = FALSE
)
utils::write.csv(
  data.frame(
    n_rna = n_rna,
    n_sct = n_sct,
    n_marker = n_marker,
    n_in_both = n_both,
    hypergeometric_p = hyper_p,
    absent_markers = paste(sort(marker_only), collapse = ";"),
    stringsAsFactors = FALSE
  ),
  file = file.path(table_dir, "marker_sct_hypergeometric.csv"),
  row.names = FALSE
)
circle_df <- function(x0, y0, radius, id, n = 240L) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = x0 + radius * cos(theta),
    y = y0 + radius * sin(theta),
    id = id,
    stringsAsFactors = FALSE
  )
}
sct_r <- 1.85
marker_r <- 1.15
sct_x <- 0
marker_x <- 1.85
circles <- rbind(
  circle_df(sct_x, 0, sct_r, "sct"),
  circle_df(marker_x, 0, marker_r, "marker")
)
lens_x <- (
  (marker_x - sct_x)^2 + sct_r^2 - marker_r^2
) / (2 * (marker_x - sct_x))
both_sorted <- sort(in_both)
n_left <- ceiling(length(both_sorted) / 2)
marker_venn <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = circles,
    ggplot2::aes(x = x, y = y, group = id, fill = id),
    colour = "#333333",
    alpha = 0.28,
    linewidth = 0.4
  ) +
  ggplot2::annotate(
    "text",
    x = sct_x - 0.55,
    y = 0.15,
    label = format(n_sct - n_both, big.mark = ","),
    size = 4.5
  ) +
  ggplot2::annotate(
    "text",
    x = lens_x,
    y = 0.95,
    label = as.character(n_both),
    size = 4.5
  ) +
  ggplot2::annotate(
    "text",
    x = marker_x + 0.42,
    y = 0.72,
    label = as.character(length(marker_only)),
    size = 4.5
  ) +
  ggplot2::annotate(
    "text",
    x = lens_x - 0.22,
    y = -0.15,
    label = paste(both_sorted[seq_len(n_left)], collapse = "\n"),
    colour = "#2ca25f",
    size = 2.5,
    fontface = "bold",
    lineheight = 0.88
  ) +
  ggplot2::annotate(
    "text",
    x = lens_x + 0.18,
    y = -0.15,
    label = paste(both_sorted[-seq_len(n_left)], collapse = "\n"),
    colour = "#2ca25f",
    size = 2.5,
    fontface = "bold",
    lineheight = 0.88
  ) +
  ggplot2::annotate(
    "text",
    x = marker_x + 0.58,
    y = -0.05,
    label = paste(sort(marker_only), collapse = "\n"),
    colour = "#e34a33",
    size = 3,
    fontface = "bold",
    lineheight = 0.9
  ) +
  ggplot2::annotate(
    "text",
    x = sct_x - 0.15,
    y = sct_r + 0.18,
    label = "SCT pre-filter",
    fontface = "bold",
    size = 3.4
  ) +
  ggplot2::annotate(
    "text",
    x = marker_x,
    y = marker_r + 0.22,
    label = "Marker genes",
    fontface = "bold",
    size = 3.4
  ) +
  ggplot2::scale_fill_manual(
    values = c(sct = "#d9d9d9", marker = "#74a9cf"),
    guide = "none"
  ) +
  ggplot2::coord_fixed(
    xlim = c(sct_x - sct_r - 0.15, marker_x + marker_r + 0.35),
    ylim = c(-sct_r - 0.15, sct_r + 0.55),
    clip = "off"
  ) +
  ggplot2::labs(
    title = "SCT counts allow-list and signalling markers",
    subtitle = sprintf(
      paste0(
        "Hypergeometric P(X >= %d) = %.2g. ",
        "Uniform draw of %d markers from %s RNA genes; ",
        "%s of those genes are in the SCT allow-list."
      ),
      n_both,
      hyper_p,
      n_marker,
      format(n_rna, big.mark = ","),
      format(n_sct, big.mark = ",")
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.margin = ggplot2::margin(12, 28, 12, 16)
  )
save_pages(
  grobs = list(cowplot::as_grob(marker_venn)),
  filename = "venn_markers_vs_counts.pdf",
  width = 10,
  height = 8
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
# ==========================================================================
# SECTION 5 · Global between versus within stacked bars ----
# ==========================================================================

global_bar_pages <- function(ranks, weighting_label, share_title, y_label) {
  ranks <- ranks[
    is.finite(ranks$between_over_total) & !is.na(ranks$rank),
    ,
    drop = FALSE
  ]
  grobs <- list()
  for (time_level in time_levels) {
    df <- ranks[ranks$time_point == time_level, , drop = FALSE]
    if (nrow(df) == 0L) {
      next
    }
    n_cut <- max(df$rank[df$selected %in% TRUE], na.rm = TRUE)
    if (!is.finite(n_cut)) {
      n_cut <- min(500L, max(df$rank))
    }
    df <- df[df$rank <= n_cut, , drop = FALSE]
    df <- df[order(df$rank), , drop = FALSE]
    pieces <- rbind(
      data.frame(
        rank = df$rank,
        ymin = 0,
        ymax = pmin(pmax(df$between_over_total, 0), 1),
        component = "between",
        stringsAsFactors = FALSE
      ),
      data.frame(
        rank = df$rank,
        ymin = pmin(pmax(df$between_over_total, 0), 1),
        ymax = 1,
        component = "within",
        stringsAsFactors = FALSE
      )
    )
    pieces$component <- factor(
      pieces$component,
      levels = c("between", "within")
    )
    n_lab <- min(10L, nrow(df))
    top_n <- df[seq_len(n_lab), , drop = FALSE]
    x_lab <- seq(n_cut * 0.08, n_cut * 0.92, length.out = n_lab)
    leaders <- data.frame(
      x = top_n$rank,
      xend = x_lab,
      y = 1,
      yend = 1.22,
      label = sprintf(
        "%s\n%s\n%.3f",
        top_n$gene,
        top_n$top_cell_type,
        top_n$between_over_within
      ),
      stringsAsFactors = FALSE
    )
    line_key <- "Dashed line: between share equals within share"
    plot <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = pieces,
        ggplot2::aes(
          xmin = rank - 0.5,
          xmax = rank + 0.5,
          ymin = ymin,
          ymax = ymax,
          fill = component
        ),
        colour = NA
      ) +
      ggplot2::geom_hline(
        data = data.frame(y = 0.5, role = line_key),
        ggplot2::aes(yintercept = y, linetype = role),
        colour = "#e34a33",
        linewidth = 1.15
      ) +
      ggplot2::geom_segment(
        data = leaders,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        linewidth = 0.3,
        colour = "#333333"
      ) +
      ggplot2::geom_label(
        data = leaders,
        ggplot2::aes(x = xend, y = 1.48, label = label),
        size = 2.2,
        linewidth = 0.2,
        fill = "white",
        lineheight = 0.88,
        label.padding = ggplot2::unit(0.12, "lines")
      ) +
      ggplot2::scale_fill_manual(
        values = c(between = "#2ca25f", within = "#e34a33"),
        labels = c(
          between = "Between cell types",
          within = "Within cell types"
        ),
        name = NULL
      ) +
      ggplot2::scale_linetype_manual(
        values = stats::setNames("dashed", line_key),
        name = NULL
      ) +
      ggplot2::guides(
        linetype = ggplot2::guide_legend(
          override.aes = list(colour = "#e34a33", linewidth = 1.15)
        )
      ) +
      ggplot2::scale_x_continuous(
        breaks = sort(unique(c(1, 10, 50, 100, 250, n_cut)))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1.85),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        x = "Gene rank (decreasing between / within)",
        y = y_label,
        title = share_title,
        subtitle = paste0(
          pretty_time(time_level),
          " · ",
          weighting_label,
          "\nGreen: steadier between-type signal (marker-like). ",
          "Red: variation inside types (stress-like). ",
          "Callouts: gene, highest-mean cell type, between/within ratio."
        )
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        plot.margin = ggplot2::margin(8, 12, 8, 8)
      )
    grobs[[time_level]] <- cowplot::as_grob(plot)
  }
  grobs
}

global <- read_table("ranks_global.csv")
bar_grobs <- global_bar_pages(
  ranks = global,
  weighting_label = "Equal cell weight.",
  share_title = expression(
    "Share" ~ frac(SS[between], SS[between] + SS[within])
  ),
  y_label = "Share of residual sum of squares"
)
if (length(bar_grobs) > 0L) {
  save_pages(
    grobs = bar_grobs,
    filename = "bar_between_within_global.pdf",
    width = 16,
    height = 8
  )
  message("Wrote bar_between_within_global.pdf")
}

equal_path <- file.path(table_dir, "ranks_global_equal_type.csv")
if (file.exists(equal_path)) {
  equal <- utils::read.csv(
    equal_path,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  equal_grobs <- global_bar_pages(
    ranks = equal,
    weighting_label = "Equal type weight (1/J per scored type).",
    share_title = expression(
      "Share" ~ frac(MS[between], MS[between] + MS[within])
    ),
    y_label = "Share of residual mean squares"
  )
  if (length(equal_grobs) > 0L) {
    save_pages(
      grobs = equal_grobs,
      filename = "bar_between_within_equal_type.pdf",
      width = 16,
      height = 8
    )
    message("Wrote bar_between_within_equal_type.pdf")
  }

  venn_weight <- list()
  overlap_rows <- list()
  for (time_level in time_levels) {
    cell_genes <- global$gene[
      global$time_point == time_level & global$selected %in% TRUE
    ]
    type_genes <- equal$gene[
      equal$time_point == time_level & equal$selected %in% TRUE
    ]
    if (length(cell_genes) == 0L || length(type_genes) == 0L) {
      next
    }
    overlap_rows[[time_level]] <- data.frame(
      time_point = time_level,
      n_equal_cell = length(unique(cell_genes)),
      n_equal_type = length(unique(type_genes)),
      n_shared = length(intersect(cell_genes, type_genes)),
      stringsAsFactors = FALSE
    )
    venn_plot <- ggVennDiagram::ggVennDiagram(
      list(
        `Equal cell weight` = unique(cell_genes),
        `Equal type weight` = unique(type_genes)
      ),
      label_alpha = 0
    ) +
      ggplot2::labs(
        title = paste0(
          pretty_time(time_level),
          ": top ",
          length(unique(cell_genes)),
          " genes, equal cell weight versus equal type weight"
        ),
        subtitle = paste(
          length(intersect(cell_genes, type_genes)),
          "genes are shared."
        )
      ) +
      ggplot2::scale_fill_gradient(low = "#f7f7f7", high = "#74a9cf") +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        plot.margin = ggplot2::margin(12, 24, 12, 24)
      )
    venn_weight[[time_level]] <- cowplot::as_grob(venn_plot)
  }
  if (length(venn_weight) > 0L) {
    save_pages(
      grobs = venn_weight,
      filename = "venn_global_cell_vs_type_weight.pdf",
      width = 8,
      height = 7
    )
    utils::write.csv(
      do.call(rbind, overlap_rows),
      file = file.path(table_dir, "global_weighting_overlap.csv"),
      row.names = FALSE
    )
    message("Wrote venn_global_cell_vs_type_weight.pdf")
  }
} else {
  message(
    "ranks_global_equal_type.csv is absent; ",
    "skipping the equal-type bar and the weighting Venn."
  )
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
utils::write.csv(
  show_df,
  file = file.path(out_dir, "strategy2_celltype_counts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
