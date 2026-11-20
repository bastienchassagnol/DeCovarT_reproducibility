#' Figures of the post hoc contrasts (appendix of chapter 2).
#'
#' Plots only; `scripts/02_04_vis_deg.R` saves them. Inputs are elements
#' of `output/feature_selection/deg/deg.rds`. Requires `R/02_00_utils.R`.

deg_test_labels <- c(
  lr = "LR (Wald)",
  ql_f = "QL F",
  anova_vst = "ANOVA",
  welch_perm = "Permuted Welch"
)
deg_design_labels <- c(
  a_flat_bh = "Flat BH on A",
  a_holm = "Stage-wise, Holm on A",
  a_graph = "Stage-wise, graph: A",
  b_graph = "Stage-wise, graph: B"
)
deg_status_colours <- c(
  "A and B (above every type)" = "#c45c26",
  "A only (above the average)" = "#1d4e89",
  "Not confirmed" = "#c9c3b8"
)

#' Confirmation status of each gene x type for one family.
#'
#' @param d Rows of `deg$decisions` for one test.
#' @return Factor with the levels of `deg_status_colours`.
deg_status <- function(d) {
  factor(
    ifelse(d$b_graph, names(deg_status_colours)[1L],
      ifelse(d$a_graph, names(deg_status_colours)[2L], names(deg_status_colours)[3L])
    ),
    levels = names(deg_status_colours)
  )
}

#' Volcano plots, one panel per type and time.
#'
#' x: log2 fold change of the type against the equal-weighted geometric
#' mean of the other types (NB fit). y: -log10 one-sided p-value of the
#' Wald contrast (continuous; the permutation p-value is discrete with a
#' floor). Colour: decision of the `colour_test` family.
#'
#' @param dec `deg$decisions`.
#' @param colour_test Family whose decisions colour the points.
#' @param n_label Genes labelled per panel (largest log2FC among A and B).
#' @param y_cap Cap of the y axis (values above are drawn at the cap).
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_deg_volcano <- function(dec, colour_test = "welch_perm", n_label = 6L, y_cap = 300, title = NULL, subtitle = NULL) {
  y <- dec[dec$test == "lr", c("time_point", "gene", "type", "log2fc", "p_a")]
  col <- dec[dec$test == colour_test, ]
  col$status <- deg_status(col)
  d <- merge(y, col[, c("time_point", "gene", "type", "status")], by = c("time_point", "gene", "type"))
  d$mlog_p <- pmin(-log10(pmax(d$p_a, .Machine$double.xmin)), y_cap)
  d$panel <- paste0(d$time_point, ": ", d$type)
  d <- d[order(d$status, decreasing = TRUE), , drop = FALSE]
  lab <- do.call(rbind, lapply(split(d, d$panel), function(r) {
    r <- r[r$status == names(deg_status_colours)[1L], , drop = FALSE]
    utils::head(r[order(-r$log2fc), , drop = FALSE], n_label)
  }))
  ggplot2::ggplot(d, ggplot2::aes(x = log2fc, y = mlog_p, colour = status)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#1c2b24", linewidth = 0.3) +
    ggplot2::geom_point(size = 0.4, alpha = 0.5) +
    ggrepel::geom_text_repel(
      data = lab,
      ggplot2::aes(label = gene),
      size = 2.4,
      fontface = "italic",
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = deg_status_colours,
      name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))
    ) +
    ggplot2::facet_wrap(~panel, ncol = 4, scales = "free") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "log2 fold change, type vs geometric mean of the other types",
      y = sprintf("-log10 p, one-sided Wald (capped at %d)", y_cap)
    ) +
    theme_slide() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
}

#' Distribution of the one-versus-average log2 fold changes.
#'
#' Half-eye densities with the interval of the middle 66 % and 95 %, one
#' row per type, one panel per time (`ggdist`).
#'
#' @param results `deg$results`.
#' @param colours Named type colours.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_lfc_distributions <- function(results, colours, title = NULL, subtitle = NULL) {
  d <- do.call(rbind, lapply(names(results), function(tm) {
    l <- results[[tm]]$log2fc
    data.frame(
      time_point = tm,
      type = rep(colnames(l), each = nrow(l)),
      log2fc = as.vector(l)
    )
  }))
  d <- d[is.finite(d$log2fc), , drop = FALSE]
  d$type <- factor(d$type, levels = rev(intersect(names(colours), unique(d$type))))
  ggplot2::ggplot(d, ggplot2::aes(x = log2fc, y = type, fill = type)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#1c2b24", linewidth = 0.3) +
    ggdist::stat_halfeye(
      adjust = 0.6,
      .width = c(0.66, 0.95),
      point_size = 1.2,
      slab_alpha = 0.8,
      normalize = "groups"
    ) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-6, 6)) +
    ggplot2::facet_wrap(~time_point, ncol = 3, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "log2 fold change, type vs geometric mean of the other types (genes tested)",
      y = NULL
    ) +
    theme_slide()
}

#' Number of confirmed genes per type for each family and design.
#'
#' @param counts `deg$counts`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_deg_counts <- function(counts, title = NULL, subtitle = NULL) {
  long <- do.call(rbind, lapply(names(deg_design_labels), function(k) {
    data.frame(
      time_point = counts$time_point,
      type = counts$type,
      test = deg_test_labels[counts$test],
      design = deg_design_labels[[k]],
      n = counts[[k]]
    )
  }))
  long$design <- factor(long$design, levels = unname(deg_design_labels))
  long$test <- factor(long$test, levels = unname(deg_test_labels))
  long$row <- paste0(long$time_point, ": ", long$type)
  ord <- stats::aggregate(n ~ row, long[long$design == deg_design_labels[["a_graph"]], ], max)
  long$row <- factor(long$row, levels = ord$row[order(ord$n)])
  test_cols <- c("#1d4e89", "#5b8c5a", "#9a9286", "#c45c26")
  names(test_cols) <- unname(deg_test_labels)
  ggplot2::ggplot(long, ggplot2::aes(x = n, y = row, colour = test)) +
    ggplot2::geom_point(size = 1.8, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = test_cols, name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::label_comma()) +
    ggplot2::facet_wrap(~design, nrow = 1) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Confirmed genes", y = NULL) +
    theme_slide()
}

#' Bimodality coefficient of a sample.
#'
#' `BC = (g^2 + 1) / (k + 3 (n - 1)^2 / ((n - 2)(n - 3)))`, with `g` the
#' sample skewness and `k` the sample excess kurtosis; values above 5/9
#' (the value of a uniform law) suggest bimodality (Pfister et al. 2013).
#'
#' @param x Numeric vector (n > 3).
#' @return The coefficient, or `NA`.
bimodality_coefficient <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4L || stats::sd(x) == 0) {
    return(NA_real_)
  }
  z <- x - mean(x)
  m2 <- mean(z^2)
  g <- mean(z^3) / m2^1.5 * sqrt(n * (n - 1)) / (n - 2)
  k <- ((n + 1) * (mean(z^4) / m2^2 - 3) + 6) * (n - 1) / ((n - 2) * (n - 3))
  (g^2 + 1) / (k + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}

#' Within-type expression of the most bimodal confirmed genes.
#'
#' @param bim Data frame from the script: `time_point`, `type`, `gene`,
#'   `share_zero`, `expected_zero` (under the fitted NB), `bc_nonzero`.
#' @param expr Data frame: `time_point`, `type`, `gene`, `value`
#'   (log1p normalised counts of the cells of that type).
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_bimodality <- function(bim, expr, title = NULL, subtitle = NULL) {
  expr$panel <- paste0(expr$gene, " (", expr$type, ", ", expr$time_point, ")")
  lab <- merge(unique(expr[, c("time_point", "type", "gene", "panel")]), bim, by = c("time_point", "type", "gene"))
  lab$label <- sprintf(
    "zeros %.0f %% (NB: %.0f %%); non-zero BC %.2f",
    100 * lab$share_zero, 100 * lab$expected_zero, lab$bc_nonzero
  )
  ggplot2::ggplot(expr, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 40, fill = "#d9d3c7", colour = "white") +
    ggplot2::geom_text(
      data = lab,
      ggplot2::aes(x = Inf, y = Inf, label = label),
      hjust = 1.02,
      vjust = 1.4,
      size = 2.5,
      colour = "#1c2b24"
    ) +
    ggplot2::facet_wrap(~panel, scales = "free", ncol = 3) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "log(1 + count / size factor), cells of the type",
      y = "Cells"
    ) +
    theme_slide() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
}

#' Old (SCTransform top-500) against new (BH gate) panel, per time.
#'
#' @param cmp Data frame: `time_point`, `class` (`"old only"`, `"both"`,
#'   `"new only"`), `n`.
#' @param ranks Data frame of shared genes: `time_point`, `old_rank`,
#'   `new_rank`.
#' @param title,subtitle Plot title and subtitle.
#' @return A patchwork of two ggplots.
plot_panel_comparison <- function(cmp, ranks, title = NULL, subtitle = NULL) {
  cls <- c("old only" = "#9a9286", "both" = "#5b8c5a", "new only" = "#1d4e89")
  cmp$class <- factor(cmp$class, levels = names(cls))
  p1 <- ggplot2::ggplot(cmp, ggplot2::aes(x = n, y = time_point, fill = class)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(n)),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 2.8,
      colour = "white"
    ) +
    ggplot2::scale_fill_manual(values = cls, name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::label_comma()) +
    ggplot2::labs(x = "Genes", y = NULL) +
    theme_slide()
  rho <- do.call(rbind, lapply(split(ranks, ranks$time_point), function(r) {
    data.frame(
      time_point = r$time_point[1L],
      label = sprintf("Spearman %.2f, n = %d", stats::cor(r$old_rank, r$new_rank, method = "spearman"), nrow(r))
    )
  }))
  p2 <- ggplot2::ggplot(ranks, ggplot2::aes(x = old_rank, y = new_rank)) +
    ggplot2::geom_point(size = 0.6, alpha = 0.5, colour = "#5b8c5a") +
    ggplot2::geom_text(
      data = rho,
      ggplot2::aes(x = -Inf, y = Inf, label = label),
      hjust = -0.05,
      vjust = 1.3,
      size = 2.8,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(x = "Rank in the SCTransform panel (top 500)", y = "Rank by R2_D in the BH panel") +
    theme_slide()
  patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 1.6)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle, theme = theme_slide())
}
