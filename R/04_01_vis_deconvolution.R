#' Figures of the deconvolution benchmark.
#'
#' Requires `R/02_00_utils.R` (`theme_slide()`). The figure functions of
#' the second pass are in `scripts/archives/second_pass_variants/`.

benchmark_colours <- c(
  first = "#1d4e89",
  second = "#c45c26",
  ink = "#1c2b24",
  context = "#8a8f8c"
)
gene_set_colours <- c(
  "Welch-gated genes" = "#1d4e89",
  "deviance share \u2265 1%" = "#5a8f5a",
  "deviance share \u2265 5%" = "#c45c26"
)
gene_set_words <- c(
  welch = "Welch-gated genes",
  r2_01 = "deviance share \u2265 1%",
  r2_05 = "deviance share \u2265 5%"
)

#' One metric per method, time and gene set, with its replicates.
#'
#' Dot-and-range plot: one row per method, one panel per gene set (rows)
#' and time (columns). Point: mean of the three standard-Wnt replicates;
#' segment: their range. Colour: the group of the method (DeCovarT with
#' the level chosen by the bulk likelihood, DeCovarT with the level
#' chosen on held-out cells, other methods). Methods are ordered by their
#' mean over the Welch panels.
#'
#' @param metrics Rows of `benchmark_metrics.csv`, with `label` and
#'   `group` columns.
#' @param metric `"rmse"` or `"aitchison"`.
#' @param title,subtitle,x_label Text of the figure.
#' @return A ggplot.
plot_benchmark_metric <- function(metrics, metric, title, subtitle, x_label) {
  metrics$value <- metrics[[metric]]
  metrics <- metrics[is.finite(metrics$value), ]
  metrics$gene_set_label <- factor(gene_set_words[metrics$gene_set], gene_set_words)
  metrics$time <- sub("h$", " h", metrics$time)
  rank <- with(
    metrics[metrics$gene_set == "welch", ],
    sort(tapply(value, label, mean), decreasing = TRUE)
  )
  metrics$label <- factor(metrics$label, union(names(rank), unique(metrics$label)))
  summary <- do.call(rbind, lapply(
    split(metrics, list(metrics$label, metrics$time, metrics$gene_set_label), drop = TRUE),
    function(d) {
      data.frame(d[1L, c("label", "time", "gene_set_label", "group")],
                 mean = mean(d$value), low = min(d$value), high = max(d$value))
    }
  ))
  ggplot2::ggplot(summary, ggplot2::aes(y = label, colour = group)) +
    ggplot2::geom_segment(ggplot2::aes(x = low, xend = high, yend = label), linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = mean), size = 1.8) +
    ggplot2::facet_grid(gene_set_label ~ time, scales = "free_x") +
    ggplot2::scale_colour_manual(values = method_group_colours, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, NA)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    theme_slide() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(size = 10),
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0)
    )
}

method_group_colours <- c(
  "DeCovarT, level by bulk likelihood" = "#1d4e89",
  "DeCovarT, level by held-out cells" = "#7fa6c9",
  "Other methods" = "#c45c26"
)

#' Agreement between each bulk and the reference mixture at the observed
#' composition: one panel per measure, one row per gene set, one point
#' per standard-Wnt replicate, along time.
#'
#' @param agreement Rows of `bulk_reference_agreement.csv`.
#' @param title,subtitle Text of the figure.
#' @return A ggplot.
plot_bulk_agreement <- function(agreement, title, subtitle) {
  long <- rbind(
    data.frame(agreement[c("time", "gene_set")], measure = "Pearson correlation, log counts", value = agreement$pearson_log),
    data.frame(agreement[c("time", "gene_set")], measure = "Share of genes > 4-fold off the sample median ratio", value = agreement$share_over_4_fold)
  )
  long$gene_set <- factor(gene_set_words[long$gene_set], gene_set_words)
  long$time <- sub("h$", " h", long$time)
  means <- stats::aggregate(value ~ time + gene_set + measure, long, mean)
  ggplot2::ggplot(long, ggplot2::aes(time, value, colour = gene_set)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5), size = 1.2, alpha = 0.5) +
    ggplot2::geom_point(data = means, position = ggplot2::position_dodge(width = 0.5), size = 3) +
    ggplot2::facet_wrap(~measure, scales = "free_y") +
    ggplot2::scale_colour_manual(values = gene_set_colours, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_slide()
}

#' Error at 96 h on five parts against three parts (paraxial mesoderm
#' summed), one dumbbell per method.
#'
#' @param lineage Rows of `lineage_collapsed_metrics.csv`, with `label`.
#' @param title,subtitle,x_label Text of the figure.
#' @return A ggplot.
plot_lineage_collapse <- function(lineage, title, subtitle, x_label) {
  m <- stats::aggregate(cbind(rmse_five_parts, rmse_three_parts) ~ label + gene_set, lineage, mean)
  m$gene_set <- factor(gene_set_words[m$gene_set], gene_set_words)
  order <- with(m[m$gene_set == gene_set_words[["welch"]], ], label[order(-rmse_five_parts)])
  m$label <- factor(m$label, union(order, unique(m$label)))
  ggplot2::ggplot(m, ggplot2::aes(y = label)) +
    ggplot2::geom_segment(ggplot2::aes(x = rmse_three_parts, xend = rmse_five_parts, yend = label),
                          colour = "#d9d3c7", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(x = rmse_five_parts, colour = "five scored types"), size = 2) +
    ggplot2::geom_point(ggplot2::aes(x = rmse_three_parts, colour = "paraxial mesoderm summed"), size = 2) +
    ggplot2::facet_wrap(~gene_set, nrow = 1L) +
    ggplot2::scale_colour_manual(
      values = c("five scored types" = "#c45c26", "paraxial mesoderm summed" = "#1d4e89"),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(limits = c(0, NA)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    theme_slide() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
}

#' Selection criteria and error along one tuning value.
#'
#' One column per time, one row per quantity. Each criterion is drawn as
#' its difference to its own maximum, so the selected level is at zero.
#'
#' @param curves Data frame with `time`, `tuning`, `quantity` (factor,
#'   one level per row of panels), `value`, and `replicate` (`NA` for a
#'   summary line).
#' @param title,subtitle,x_label Text of the figure.
#' @return A ggplot.
plot_selection_curves <- function(curves, title, subtitle, x_label) {
  curves$time <- sub("h$", " h", curves$time)
  lines <- curves[is.na(curves$replicate), ]
  best <- do.call(
    rbind,
    lapply(split(lines, list(lines$time, lines$quantity), drop = TRUE), function(d) {
      lowest <- grepl("error", d$quantity[[1L]])
      d[if (lowest) which.min(d$value) else which.max(d$value), ]
    })
  )
  ggplot2::ggplot(curves, ggplot2::aes(tuning, value)) +
    ggplot2::geom_point(
      data = curves[!is.na(curves$replicate), ],
      colour = "#1d4e89",
      size = 0.9,
      alpha = 0.5
    ) +
    ggplot2::geom_line(data = lines, colour = "#1d4e89") +
    ggplot2::geom_point(data = best, colour = "#c45c26", size = 2.6) +
    ggplot2::facet_grid(quantity ~ time, scales = "free_y", switch = "y") +
    ggplot2::labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    theme_slide() +
    ggplot2::theme(strip.placement = "outside")
}

#' Correlation between the mean profiles of the cell types.
#'
#' One panel per time (columns) and coefficient (rows): Pearson on the
#' mean raw counts per cell, the scale the solvers see, and Kendall's
#' tau, which depends on the ranks only. Cell types are ordered within a
#' time by hierarchical clustering of `1 - Pearson` (average linkage)
#' and numbered in that order; columns carry the number only, so no
#' label is rotated. The coefficient is printed in each tile.
#'
#' @param correlations `profile_correlation.csv` as a data frame:
#'   `time`, `gene_set`, `type_a`, `type_b`, `pearson`, `kendall`.
#' @param gene_set `"welch"` or `"r2_05"`.
#' @param title,subtitle Text of the figure.
#' @return A ggplot.
plot_profile_correlation <- function(
  correlations,
  gene_set,
  title = NULL,
  subtitle = NULL
) {
  d <- correlations[correlations$gene_set == gene_set, , drop = FALSE]
  ordered <- lapply(split(d, d$time), function(x) {
    types <- unique(x$type_a)
    m <- matrix(
      1,
      length(types),
      length(types),
      dimnames = list(types, types)
    )
    m[cbind(x$type_a, x$type_b)] <- x$pearson
    types[stats::hclust(stats::as.dist(1 - m), method = "average")$order]
  })
  key <- function(type, time) paste(type, time, sep = "___")
  levels_key <- unlist(lapply(names(ordered), function(tm) {
    key(ordered[[tm]], tm)
  }))
  long <- rbind(
    data.frame(
      d[, c("time", "type_a", "type_b")],
      coefficient = "Pearson, mean counts",
      value = d$pearson
    ),
    data.frame(
      d[, c("time", "type_a", "type_b")],
      coefficient = "Kendall's tau",
      value = d$kendall
    )
  )
  long$coefficient <- factor(
    long$coefficient,
    levels = c("Pearson, mean counts", "Kendall's tau")
  )
  long$a <- factor(key(long$type_a, long$time), levels = levels_key)
  long$b <- factor(key(long$type_b, long$time), levels = rev(levels_key))
  long$time_label <- sub("h$", " h", long$time)
  lowest <- floor(min(long$value) * 10) / 10
  rank_of <- unlist(lapply(names(ordered), function(tm) {
    stats::setNames(seq_along(ordered[[tm]]), key(ordered[[tm]], tm))
  }))
  row_label <- function(x) paste0(sub("___.*$", "", x), "  ", rank_of[x])
  column_label <- function(x) as.character(rank_of[x])
  ggplot2::ggplot(long, ggplot2::aes(a, b, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = formatC(value, format = "f", digits = 2),
        colour = value > (lowest + 1) / 2
      ),
      size = 2.9
    ) +
    ggplot2::facet_wrap(
      coefficient ~ time_label,
      scales = "free",
      nrow = 2
    ) +
    ggplot2::scale_fill_gradient(
      low = "#f3f1ec",
      high = benchmark_colours[["first"]],
      limits = c(lowest, 1),
      name = NULL
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "white", `FALSE` = benchmark_colours[["ink"]]),
      guide = "none"
    ) +
    ggplot2::scale_x_discrete(labels = column_label, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(labels = row_label, expand = c(0, 0)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_slide(base_size = 10) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(2.2, "cm"),
      legend.key.height = ggplot2::unit(0.3, "cm")
    )
}

#' Error of each strategy of the light screen at the level its bulk
#' likelihood chooses.
#'
#' One row per strategy (covariance input, arm, scale model), ordered by
#' mean error; one point per replicate triplet, coloured by pseudo-bulk
#' version; one column per time.
#'
#' @param chosen Data frame with `strategy`, `time`, `version`, `rmse`
#'   (mean over the three replicates of a triplet).
#' @param title,subtitle,x_label Text of the figure.
#' @return A ggplot.
plot_scale_screen <- function(chosen, title, subtitle, x_label) {
  chosen$time <- sub("h$", " h", chosen$time)
  order <- sort(tapply(chosen$rmse, chosen$strategy, mean), decreasing = TRUE)
  chosen$strategy <- factor(chosen$strategy, names(order))
  ggplot2::ggplot(chosen, ggplot2::aes(rmse, strategy, colour = version)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::facet_wrap(~time, nrow = 1L) +
    ggplot2::scale_colour_manual(
      values = c("as is" = "#1d4e89", "with platform ratio" = "#c45c26"),
      name = NULL
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    theme_slide()
}
