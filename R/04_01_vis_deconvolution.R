#' Figures of the deconvolution benchmark.
#'
#' Every function takes the long table written by
#' `scripts/04_01_stats_deconvolution.R` (`proportions_long.csv`, or the
#' metrics derived from it) and returns a ggplot. Requires
#' `R/02_00_utils.R` (`theme_slide()`).

# Colour roles (data-visualisation rule): first series, second series,
# annotation.
benchmark_colours <- c(
  first = "#1d4e89",
  second = "#c45c26",
  ink = "#1c2b24",
  context = "#8a8f8c"
)

#' Display labels shared by the benchmark figures.
#'
#' Adds to a table of estimates or metrics:
#' - `method_label`: the algorithm, and for DeCovarT, its variants and
#'   the multivariate RNA-Sieve their covariance input (one facet or one
#'   row per value, in a fixed order);
#' - `ratio_label`: whether the bulk was divided by the gene-wise
#'   platform ratio;
#' - `gene_set_label`: the gene set, in words.
#'
#' @param x Data frame with `algorithm`, `covariance`, `platform_ratio`
#'   and `gene_set`.
#' @param wrap If `TRUE`, method labels are split over several lines
#'   (facet strips).
#' @return `x` with the three factor columns added.
add_benchmark_labels <- function(x, wrap = FALSE) {
  sep <- if (wrap) "\n" else " "
  covariance_words <- c(
    raw = "raw",
    denoised_bio = "denoised",
    denoised_plugin = paste0("denoised + Poisson,", sep, "plug-in a"),
    denoised_fixed = paste0("denoised + Poisson,", sep, "fixed-point a")
  )
  comparator_words <- c(
    "DWLS" = "DWLS",
    "SCDC" = "SCDC",
    "BayesPrism" = "BayesPrism",
    "RNA-Sieve" = "RNA-Sieve",
    "DeMixSC wNNLS" = "DeMixSC wNNLS",
    "MEAD (pseudo-individuals)" = paste0("MEAD,", sep, "pseudo-individuals"),
    "ReDeconv" = "ReDeconv",
    "RLM" = "Robust regression"
  )
  sieve_words <- c(
    diag_raw = "diagonal",
    full_raw = "full raw",
    full_bio = "full denoised"
  )
  shrinkage_levels <- c("0.25", "0.50", "0.75", "1.00")
  levels_method <- c(
    unname(comparator_words),
    paste0("RNA-Sieve (R),", sep, sieve_words),
    paste0("DeCovarT", sep, covariance_words),
    paste0("DeCovarT total", sep, "covariance, ", covariance_words[1:2]),
    paste0("DeCovarT Firth,", sep, covariance_words[1:2]),
    paste0("DeCovarT raw,", sep, "\u03bb = ", shrinkage_levels)
  )
  label <- unname(comparator_words[x$algorithm])
  is_decovart <- x$algorithm == "DeCovarT"
  label[is_decovart] <- paste0(
    "DeCovarT",
    sep,
    covariance_words[x$covariance[is_decovart]]
  )
  pick <- x$algorithm == "RNA-Sieve (R)"
  label[pick] <- paste0("RNA-Sieve (R),", sep, sieve_words[x$covariance[pick]])
  pick <- x$algorithm == "DeCovarT total covariance"
  label[pick] <- paste0(
    "DeCovarT total",
    sep,
    "covariance, ",
    covariance_words[x$covariance[pick]]
  )
  pick <- x$algorithm == "DeCovarT Firth"
  label[pick] <- paste0(
    "DeCovarT Firth,",
    sep,
    covariance_words[x$covariance[pick]]
  )
  pick <- x$algorithm == "DeCovarT shrinkage"
  label[pick] <- paste0(
    "DeCovarT raw,",
    sep,
    "\u03bb = ",
    sub("^lambda_", "", x$covariance[pick])
  )
  if (anyNA(label)) {
    stop(
      "No display label for: ",
      paste(unique(x$algorithm[is.na(label)]), collapse = ", ")
    )
  }
  x$method_label <- factor(label, levels = levels_method)
  x$ratio_label <- factor(
    ifelse(x$platform_ratio == "gene", "bulk / b_g", "bulk as is"),
    levels = c("bulk as is", "bulk / b_g")
  )
  gene_set_words <- c(
    welch = "Welch-gated\ngenes",
    r2_05 = "deviance\nshare \u2265 5%"
  )
  x$gene_set_label <- factor(
    gene_set_words[x$gene_set],
    levels = unname(gene_set_words)
  )
  x
}

#' Observed against estimated proportions, one panel per method and time.
#'
#' Standard-Wnt samples of one gene set. Each point is one cell type in
#' one bulk replicate; its observed value is the fraction of scored
#' cells of that type in the single-cell annotation at that time. For
#' DeCovarT the fits of one optimiser are drawn. The text above each
#' panel gives the Pearson correlation and the root mean squared error
#' over the points of the panel, separately for the two treatments of
#' the bulk.
#'
#' @param estimates `proportions_long.csv` as a data frame.
#' @param observed `observed_proportions.csv` as a data frame.
#' @param gene_set `"welch"` or `"r2_05"`.
#' @param colours Named vector of cell-type colours.
#' @param optimiser DeCovarT optimiser drawn: `"lbfgsb"` (default) or
#'   `"newton"`.
#' @param title,subtitle Text of the figure.
#' @return A ggplot.
plot_observed_vs_estimated <- function(
  estimates,
  observed,
  gene_set,
  colours,
  optimiser = c("lbfgsb", "newton"),
  title = NULL,
  subtitle = NULL
) {
  optimiser <- match.arg(optimiser)
  keep <- estimates$arm == "standard" &
    estimates$gene_set == gene_set &
    (estimates$algorithm != "DeCovarT" | estimates$optimiser %in% optimiser)
  d <- add_benchmark_labels(estimates[keep, , drop = FALSE], wrap = TRUE)
  d$observed <- observed$observed[match(
    paste(d$time, d$cell_type),
    paste(observed$time, observed$cell_type)
  )]
  d$time_label <- sub("h$", " h", d$time)
  d <- d[is.finite(d$proportion), , drop = FALSE]

  groups <- split(
    d,
    list(d$method_label, d$time_label, d$ratio_label),
    drop = TRUE
  )
  notes <- do.call(
    rbind,
    lapply(groups, function(g) {
      r <- if (stats::sd(g$proportion) > 0) {
        stats::cor(g$proportion, g$observed)
      } else {
        NA_real_
      }
      data.frame(
        method_label = g$method_label[[1L]],
        time_label = g$time_label[[1L]],
        ratio_label = g$ratio_label[[1L]],
        text = paste0(
          if (g$ratio_label[[1L]] == "bulk / b_g") "/ b_g: " else "as is: ",
          "r ",
          if (is.na(r)) "n/a" else formatC(r, format = "f", digits = 2),
          ", RMSE ",
          formatC(
            sqrt(mean((g$proportion - g$observed)^2)),
            format = "f",
            digits = 2
          )
        ),
        stringsAsFactors = FALSE
      )
    })
  )
  notes <- do.call(
    rbind,
    lapply(
      split(notes, list(notes$method_label, notes$time_label), drop = TRUE),
      function(n) {
        n <- n[order(n$ratio_label), , drop = FALSE]
        data.frame(
          method_label = n$method_label[[1L]],
          time_label = n$time_label[[1L]],
          text = paste(n$text, collapse = "\n"),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  ggplot2::ggplot(d, ggplot2::aes(observed, proportion)) +
    ggplot2::annotate(
      "segment",
      x = 0,
      y = 0,
      xend = 1,
      yend = 1,
      linetype = "dashed",
      colour = benchmark_colours[["context"]],
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = cell_type, shape = ratio_label),
      size = 1.7,
      stroke = 0.7,
      alpha = 0.85
    ) +
    ggplot2::geom_text(
      data = notes,
      ggplot2::aes(x = 0, y = 1.3, label = text),
      hjust = 0,
      vjust = 1,
      size = 2.2,
      lineheight = 0.95,
      colour = benchmark_colours[["ink"]]
    ) +
    ggplot2::facet_grid(time_label ~ method_label) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::scale_shape_manual(
      values = c("bulk as is" = 16, "bulk / b_g" = 2),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    ) +
    # The band above 1 holds the text, clear of the points.
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1.3), clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Observed proportion (scored cells of the single-cell annotation)",
      y = "Estimated proportion"
    ) +
    theme_slide(base_size = 9) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 7, lineheight = 0.95),
      strip.text.y = ggplot2::element_text(angle = 0),
      panel.spacing = ggplot2::unit(0.7, "lines"),
      legend.box = "vertical",
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, order = 1),
      shape = ggplot2::guide_legend(order = 2)
    )
}

# Axis title of each benchmark metric.
benchmark_metric_words <- c(
  rmse = "Root mean squared error of the proportions (per replicate)",
  pearson = "Pearson correlation with the observed proportions (per replicate)",
  aitchison = paste(
    "Aitchison distance to the observed composition",
    "(per replicate, zeros replaced at 0.001)"
  )
)

#' Axis of a benchmark metric: errors and distances start at zero, the
#' correlation spans -1 to 1.
#'
#' @param metric `"rmse"`, `"pearson"` or `"aitchison"`.
#' @param horizontal `TRUE` for an x-axis, `FALSE` for a y-axis.
#' @return A ggplot scale.
benchmark_metric_scale <- function(metric, horizontal = TRUE) {
  scale <- if (horizontal) {
    ggplot2::scale_x_continuous
  } else {
    ggplot2::scale_y_continuous
  }
  switch(
    metric,
    rmse = scale(
      limits = c(0, NA),
      breaks = scales::breaks_width(0.25),
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ),
    pearson = scale(limits = c(-1, 1), breaks = c(-1, 0, 1)),
    aitchison = scale(
      limits = c(0, NA),
      breaks = scales::breaks_width(2),
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    )
  )
}

#' One benchmark metric for every method, gene set and treatment.
#'
#' Rows of panels are the methods (comparators, then the four covariance
#' inputs of DeCovarT); columns are the bulk times. Inside a panel, one
#' line per gene set and treatment of the bulk. Small points are the
#' three standard-Wnt replicates, the large point their mean. For
#' DeCovarT the two optimisers are drawn side by side.
#'
#' @param metrics `benchmark_metrics.csv` as a data frame.
#' @param metric `"rmse"`, `"pearson"` or `"aitchison"`.
#' @param title,subtitle Text of the figure.
#' @return A ggplot.
plot_benchmark_metric <- function(
  metrics,
  metric = c("rmse", "pearson", "aitchison"),
  title = NULL,
  subtitle = NULL
) {
  metric <- match.arg(metric)
  d <- add_benchmark_labels(metrics)
  d$value <- d[[metric]]
  d <- d[is.finite(d$value), , drop = FALSE]
  d$time_label <- sub("h$", " h", d$time)
  d$line <- paste0(
    ifelse(d$gene_set == "welch", "Welch-gated", "share \u2265 5%"),
    ifelse(d$platform_ratio == "gene", ", bulk / b_g", "")
  )
  d$line <- factor(
    d$line,
    levels = rev(c(
      "Welch-gated",
      "Welch-gated, bulk / b_g",
      "share \u2265 5%",
      "share \u2265 5%, bulk / b_g"
    ))
  )
  d$solver <- factor(
    ifelse(
      is.na(d$optimiser),
      "comparator",
      ifelse(d$optimiser == "newton", "Newton", "L-BFGS-B or BFGS")
    ),
    levels = c("comparator", "Newton", "L-BFGS-B or BFGS")
  )
  means <- stats::aggregate(
    value ~ method_label + time_label + line + solver,
    data = d,
    FUN = mean
  )
  dodge <- ggplot2::position_dodge(width = 0.6)
  x_scale <- benchmark_metric_scale(metric)
  ggplot2::ggplot(d, ggplot2::aes(value, line, colour = solver)) +
    ggplot2::geom_point(position = dodge, size = 0.9, alpha = 0.45) +
    ggplot2::geom_point(data = means, position = dodge, size = 2.1) +
    ggplot2::facet_grid(
      method_label ~ time_label,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        comparator = benchmark_colours[["ink"]],
        Newton = benchmark_colours[["first"]],
        "L-BFGS-B or BFGS" = benchmark_colours[["second"]]
      ),
      name = NULL
    ) +
    x_scale +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = benchmark_metric_words[[metric]],
      y = NULL
    ) +
    theme_slide(base_size = 10) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
      panel.spacing.y = ggplot2::unit(0.5, "lines"),
      panel.spacing.x = ggplot2::unit(1.6, "lines"),
      panel.grid.major.y = ggplot2::element_blank()
    )
}

#' Estimated proportions of one time, standard against early Wnt.
#'
#' Columns of panels are the cell types of that time, rows the two gene
#' sets. One line per method (DeCovarT: one optimiser). Each point is one
#' bulk replicate. The dashed line is the observed fraction of scored
#' cells, which describes the standard-Wnt protocol only.
#'
#' @param estimates `proportions_long.csv` as a data frame.
#' @param observed `observed_proportions.csv` as a data frame.
#' @param time_level `"48h"`, `"72h"` or `"96h"`.
#' @param optimiser DeCovarT optimiser drawn: `"lbfgsb"` (default) or
#'   `"newton"`.
#' @param title,subtitle Text of the figure.
#' @return A ggplot.
plot_estimates_by_arm <- function(
  estimates,
  observed,
  time_level,
  optimiser = c("lbfgsb", "newton"),
  title = NULL,
  subtitle = NULL
) {
  optimiser <- match.arg(optimiser)
  keep <- estimates$time == time_level &
    (estimates$algorithm != "DeCovarT" | estimates$optimiser %in% optimiser)
  d <- add_benchmark_labels(estimates[keep, , drop = FALSE])
  d <- d[is.finite(d$proportion), , drop = FALSE]
  d$line <- paste0(
    as.character(d$method_label),
    ifelse(d$platform_ratio == "gene", ", bulk / b_g", "")
  )
  line_levels <- unique(d$line[order(d$method_label, d$ratio_label)])
  d$line <- factor(d$line, levels = rev(line_levels))
  d$arm_label <- factor(
    ifelse(d$arm == "standard", "standard Wnt", "early Wnt"),
    levels = c("standard Wnt", "early Wnt")
  )
  obs <- observed[observed$time == time_level, , drop = FALSE]
  type_levels <- obs$cell_type[order(-obs$observed)]
  d$cell_type <- factor(d$cell_type, levels = type_levels)
  obs$cell_type <- factor(obs$cell_type, levels = type_levels)
  dodge <- ggplot2::position_dodge(width = 0.6)
  ggplot2::ggplot(d, ggplot2::aes(proportion, line, colour = arm_label)) +
    ggplot2::geom_vline(
      data = obs,
      ggplot2::aes(xintercept = observed),
      linetype = "dashed",
      colour = benchmark_colours[["ink"]],
      linewidth = 0.4
    ) +
    ggplot2::geom_point(position = dodge, size = 1.4, alpha = 0.8) +
    ggplot2::facet_grid(
      gene_set_label ~ cell_type,
      labeller = ggplot2::labeller(
        cell_type = ggplot2::label_wrap_gen(width = 16)
      )
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "standard Wnt" = benchmark_colours[["first"]],
        "early Wnt" = benchmark_colours[["second"]]
      ),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      oob = scales::oob_squish
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Estimated proportion (dashed: observed, standard Wnt)",
      y = NULL
    ) +
    theme_slide(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(
        colour = "#f1eee8",
        linewidth = 0.3
      ),
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0)
    )
}

#' Stages of the multivariate RNA-Sieve, one metric.
#'
#' Rows of panels are the covariance inputs, columns the bulk times.
#' Inside a panel, one line per stage, in the order they are fitted;
#' the two gene sets are drawn side by side. Small points are the three
#' standard-Wnt replicates, the large point their mean.
#'
#' @param metrics `benchmark_metrics.csv` as a data frame.
#' @param metric `"rmse"`, `"pearson"` or `"aitchison"`.
#' @param title,subtitle Text of the figure.
#' @return A ggplot.
plot_sieve_stages <- function(
  metrics,
  metric = c("rmse", "pearson", "aitchison"),
  title = NULL,
  subtitle = NULL
) {
  metric <- match.arg(metric)
  d <- metrics[metrics$algorithm == "RNA-Sieve (R)", , drop = FALSE]
  d$value <- d[[metric]]
  d <- d[is.finite(d$value), , drop = FALSE]
  stage_words <- c(
    qp_within = "QP, within-type covariance",
    qp = "QP, total covariance",
    qp_huber = "QP, Huber weights",
    qp_means = "QP + refined means",
    qp_means_var = "QP + means + gene variances",
    full_lik = "then full likelihood"
  )
  input_words <- c(
    diag_raw = "Diagonal of the\nraw covariance",
    full_raw = "Full raw\ncovariance",
    full_bio = "Full denoised\ncovariance"
  )
  d$stage_label <- factor(
    stage_words[d$optimiser],
    levels = rev(unname(stage_words))
  )
  d$input_label <- factor(
    input_words[d$covariance],
    levels = unname(input_words)
  )
  d$time_label <- sub("h$", " h", d$time)
  d$panel <- factor(
    ifelse(d$gene_set == "welch", "Welch-gated genes", "deviance share \u2265 5%"),
    levels = c("Welch-gated genes", "deviance share \u2265 5%")
  )
  means <- stats::aggregate(
    value ~ stage_label + input_label + time_label + panel,
    data = d,
    FUN = mean
  )
  dodge <- ggplot2::position_dodge(width = 0.6)
  ggplot2::ggplot(d, ggplot2::aes(value, stage_label, colour = panel)) +
    ggplot2::geom_point(position = dodge, size = 0.9, alpha = 0.45) +
    ggplot2::geom_point(data = means, position = dodge, size = 2.1) +
    ggplot2::facet_grid(input_label ~ time_label) +
    ggplot2::scale_colour_manual(
      values = c(
        "Welch-gated genes" = benchmark_colours[["first"]],
        "deviance share \u2265 5%" = benchmark_colours[["second"]]
      ),
      name = NULL
    ) +
    benchmark_metric_scale(metric) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = benchmark_metric_words[[metric]],
      y = NULL
    ) +
    theme_slide(base_size = 10) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
      panel.spacing.x = ggplot2::unit(1.6, "lines"),
      panel.grid.major.y = ggplot2::element_blank()
    )
}

#' DeCovarT on the raw covariance along the shrinkage level.
#'
#' Columns of panels are the bulk times, rows the two metrics (root mean
#' squared error and Aitchison distance). The first level is the
#' closed-form shrinkage of each cell type; the others are one common
#' value. L-BFGS-B fits, bulk as is. Small points are the three
#' standard-Wnt replicates; the line joins their means.
#'
#' @param metrics `benchmark_metrics.csv` as a data frame.
#' @param title,subtitle Text of the figure.
#' @return A ggplot.
plot_shrinkage_sweep <- function(metrics, title = NULL, subtitle = NULL) {
  base <- metrics$algorithm == "DeCovarT" &
    metrics$covariance == "raw" &
    metrics$optimiser == "lbfgsb" &
    metrics$platform_ratio == "none"
  sweep <- metrics$algorithm == "DeCovarT shrinkage"
  d <- metrics[base | sweep, , drop = FALSE]
  d$level <- ifelse(
    d$algorithm == "DeCovarT",
    "closed\nform",
    sub("^lambda_", "", d$covariance)
  )
  d$level <- factor(
    d$level,
    levels = c("closed\nform", "0.25", "0.50", "0.75", "1.00")
  )
  d$time_label <- sub("h$", " h", d$time)
  d$panel <- factor(
    ifelse(d$gene_set == "welch", "Welch-gated genes", "deviance share \u2265 5%"),
    levels = c("Welch-gated genes", "deviance share \u2265 5%")
  )
  long <- rbind(
    data.frame(
      d[, c("level", "time_label", "panel")],
      metric = "Root mean squared error",
      value = d$rmse
    ),
    data.frame(
      d[, c("level", "time_label", "panel")],
      metric = "Aitchison distance",
      value = d$aitchison
    )
  )
  long$metric <- factor(
    long$metric,
    levels = c("Root mean squared error", "Aitchison distance")
  )
  long <- long[is.finite(long$value), , drop = FALSE]
  means <- stats::aggregate(
    value ~ level + time_label + panel + metric,
    data = long,
    FUN = mean
  )
  ggplot2::ggplot(long, ggplot2::aes(level, value, colour = panel)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.45) +
    ggplot2::geom_line(
      data = means,
      ggplot2::aes(group = panel),
      linewidth = 0.6
    ) +
    ggplot2::geom_point(data = means, size = 2.1) +
    ggplot2::facet_grid(metric ~ time_label, scales = "free_y", switch = "y") +
    ggplot2::scale_colour_manual(
      values = c(
        "Welch-gated genes" = benchmark_colours[["first"]],
        "deviance share \u2265 5%" = benchmark_colours[["second"]]
      ),
      name = NULL
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0.02, 0.08))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Shrinkage of the correlations towards zero (\u03bb)",
      y = NULL
    ) +
    theme_slide(base_size = 10) +
    ggplot2::theme(
      strip.placement = "outside",
      panel.spacing.x = ggplot2::unit(1.4, "lines"),
      panel.grid.major.x = ggplot2::element_blank()
    )
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
