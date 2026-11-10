#' Tier (iii) figures: influence, transcript length, A3 shift, bulk fit,
#' tail flags, composition, stability of the platform ratio.
#'
#' Plots only; nothing is saved here (see
#' `scripts/02_03_vis_cross_platform.R`). Inputs are elements of
#' `output/feature_selection/tier_03/cross_platform.rds`. Requires
#' `R/02_00_utils.R` (`theme_slide()`, `spread_labels()`, `fmt_n()`).
#' Titles are passed by the driver, which computes the finding from the
#' outputs (data-visualisation rule #33).

signed_log <- function(z) sign(z) * log10(1 + abs(z))

variant_labels <- c(
  none = "A3: no platform ratio",
  length = "Ratio from transcript length",
  gene = "Gene-wise ratio"
)
variant_colours <- c(
  "A3: no platform ratio" = "#9a9286",
  "Ratio from transcript length" = "#1d4e89",
  "Gene-wise ratio" = "#c45c26"
)

#' Cook diagnostics for genes (bulk fit) and cells (M1 fits).
#'
#' Panel a: leverage against standardised residual of each gene in the
#' weighted bulk-on-reference regression, one facet per time, all three
#' replicates; flagged genes (large leverage and |residual| > r_max), if
#' any, in orange. Panel b: (gene, cell) pairs flagged in the M1 fits
#' per cell of the type, types ranked, one colour (#18, #50).
#'
#' @param outlier_genes Named list (times) of [cook_genes()] tables.
#' @param outlier_cells Named list (times) of [cook_cells()] outputs.
#' @param cells Named list (times) of named cell counts per type.
#' @param r_max Residual threshold drawn as a line. Default `3`.
#' @param titles Two titles (panels a and b).
#' @return A patchwork object.
plot_cook <- function(outlier_genes, outlier_cells, cells, r_max = 3, titles = NULL) {
  g <- do.call(
    rbind,
    lapply(names(outlier_genes), function(tm) {
      d <- outlier_genes[[tm]]
      d$time_point <- tm
      d
    })
  )
  g <- g[order(g$flag), , drop = FALSE]
  a <- ggplot2::ggplot(
    g,
    ggplot2::aes(x = leverage, y = std_residual, colour = flag)
  ) +
    ggplot2::geom_point(size = 0.4, alpha = 0.3) +
    ggplot2::geom_hline(
      yintercept = c(-r_max, r_max),
      linetype = "dashed",
      colour = "#1c2b24",
      linewidth = 0.3
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = r_max,
      label = paste0("|residual| = ", r_max),
      hjust = 1.05,
      vjust = -0.4,
      size = 2.6,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "#1d4e89", "TRUE" = "#c45c26"),
      guide = "none"
    ) +
    ggplot2::scale_x_log10(
      labels = scales::label_log(),
      expand = ggplot2::expansion(mult = c(0.03, 0.08))
    ) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      x = "Leverage in the bulk-on-reference fit (log scale)",
      y = "Standardised residual",
      title = titles[1L]
    ) +
    theme_slide()
  cnt <- do.call(
    rbind,
    lapply(names(outlier_cells), function(tm) {
      f <- outlier_cells[[tm]]$flagged
      n <- cells[[tm]]
      k <- table(factor(f$cell_type, levels = names(n)))
      data.frame(
        time_point = tm,
        cell_type = names(n),
        pairs = as.integer(k),
        n_cells = as.integer(n),
        per_cell = as.integer(k) / as.numeric(n)
      )
    })
  )
  key <- function(tm, ct) paste(tm, ct, sep = "___")
  cnt <- cnt[order(cnt$time_point, cnt$per_cell), , drop = FALSE]
  cnt$row <- factor(
    key(cnt$time_point, cnt$cell_type),
    levels = key(cnt$time_point, cnt$cell_type)
  )
  row_labels <- stats::setNames(
    paste0(cnt$cell_type, " (", fmt_n(cnt$n_cells), " cells)"),
    as.character(cnt$row)
  )
  b <- ggplot2::ggplot(cnt, ggplot2::aes(x = per_cell, y = row)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = per_cell, yend = row),
      colour = "grey70"
    ) +
    ggplot2::geom_point(colour = "#1d4e89", size = 2) +
    ggplot2::geom_text(
      ggplot2::aes(label = fmt_n(pairs)),
      hjust = -0.1,
      size = 2.4,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_y_discrete(labels = row_labels) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::facet_wrap(~time_point, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      x = "Flagged (gene, cell) pairs per cell of the type (label: pairs)",
      y = NULL,
      title = titles[2L]
    ) +
    theme_slide()
  patchwork::wrap_plots(a, b, ncol = 1, heights = c(1, 1.6))
}

#' Platform ratio against transcript length.
#'
#' Per time, the per-gene log2 ratio of bulk to fitted reference (A3
#' fit, [platform_log_ratio()]) against the log2 mean transcript length,
#' as 2D bins (counts on a log scale) to handle overplotting; the median
#' regression line (orange) and, dashed, a line of slope 1 through the
#' medians of both axes, the prediction for a full-length bulk against
#' 3' UMI counts.
#'
#' @param log_ratio Output of [platform_log_ratio()].
#' @param lengths Output of [gene_transcript_length()].
#' @param length_test Output of [platform_length_test()].
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_platform_length <- function(log_ratio, lengths, length_test, title = NULL, subtitle = NULL) {
  m <- log_ratio$log_ratio
  l2 <- stats::setNames(log2(lengths$mean_tx_length), lengths$gene)
  d <- do.call(rbind, lapply(colnames(m), function(tm) {
    data.frame(time_point = tm, x = l2[rownames(m)], y = m[, tm])
  }))
  d <- d[is.finite(d$x) & is.finite(d$y), , drop = FALSE]
  centre <- do.call(rbind, lapply(split(d, d$time_point), function(r) {
    data.frame(
      time_point = r$time_point[1L],
      mx = stats::median(r$x),
      my = stats::median(r$y)
    )
  }))
  lt <- merge(length_test, centre, by = "time_point")
  lines <- do.call(rbind, lapply(seq_len(nrow(lt)), function(i) {
    xr <- range(d$x[d$time_point == lt$time_point[i]])
    rbind(
      data.frame(
        time_point = lt$time_point[i],
        kind = "fit",
        x = xr,
        y = lt$intercept[i] + lt$slope[i] * xr
      ),
      data.frame(
        time_point = lt$time_point[i],
        kind = "slope 1",
        x = xr,
        y = lt$my[i] + (xr - lt$mx[i])
      )
    )
  }))
  lab <- data.frame(
    time_point = lt$time_point,
    x = -Inf,
    y = Inf,
    label = sprintf(
      "slope %.2f (95 %% CI %.2f to %.2f)\nSpearman %.2f; loess R2 %.2f",
      lt$slope,
      lt$slope_lower,
      lt$slope_upper,
      lt$spearman,
      lt$loess_r2
    )
  )
  line_labels <- c(fit = "median regression", "slope 1" = "slope 1 (prediction)")
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bin_2d(bins = 60) +
    ggplot2::geom_hline(yintercept = 0, colour = "#1c2b24", linewidth = 0.3) +
    ggplot2::geom_line(
      data = lines,
      ggplot2::aes(x = x, y = y, linetype = kind, colour = kind),
      linewidth = 0.8
    ) +
    ggplot2::geom_text(
      data = lab,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = -0.05,
      vjust = 1.2,
      size = 2.6,
      colour = "#1c2b24",
      lineheight = 0.95
    ) +
    ggplot2::scale_fill_gradient(
      low = "#dfe7f1",
      high = "#1d4e89",
      transform = "log10",
      name = "Genes per bin"
    ) +
    ggplot2::scale_colour_manual(
      values = c(fit = "#c45c26", "slope 1" = "#1c2b24"),
      labels = line_labels,
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c(fit = "solid", "slope 1" = "dashed"),
      labels = line_labels,
      name = NULL
    ) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "log2 mean transcript length (bases, mm10 UCSC knownGene)",
      y = "log2 bulk / reference ratio (A3 fit)"
    ) +
    theme_slide()
}

#' Distribution of the standardised discrepancy under three models.
#'
#' Per time, one slab per variant of the cross-platform model: the
#' distribution of `z = D / sd(D)` over genes and replicates, on a
#' signed log10 scale, with a reference row drawn from N(0, 1) (the law
#' of `z` when the model holds) and the band `|z| <= qnorm(1 - level)`
#' shaded. The share of values beyond the band is written at the right
#' of each slab.
#'
#' @param variants Named list (`none`, `length`, `gene`) of variant
#'   results, each with `tails` (named list of [bulk_tail_filter()]
#'   outputs).
#' @param level Tail level per side. Default `0.01`.
#' @param seed Seed of the N(0, 1) reference draws.
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_a3_shift <- function(variants, level = 0.01, seed = 1L, title = NULL, subtitle = NULL) {
  z_cut <- stats::qnorm(1 - level)
  ref_label <- "Expected if the model holds: N(0, 1)"
  d <- do.call(rbind, lapply(names(variants), function(v) {
    do.call(rbind, lapply(names(variants[[v]]$tails), function(tm) {
      g <- variants[[v]]$tails[[tm]]$gene
      z <- g$z[g$statement == "T4"]
      data.frame(time_point = tm, variant = variant_labels[[v]], z = z[is.finite(z)])
    }))
  }))
  times <- unique(d$time_point)
  ref <- withr::with_seed(seed, data.frame(
    time_point = rep(times, each = 20000L),
    variant = ref_label,
    z = stats::rnorm(20000L * length(times))
  ))
  d <- rbind(d, ref)
  d$variant <- factor(d$variant, levels = rev(c(ref_label, unname(variant_labels))))
  d$sz <- signed_log(d$z)
  share <- stats::aggregate(
    z ~ time_point + variant,
    data = d,
    FUN = function(z) mean(abs(z) > z_cut)
  )
  x_max <- max(abs(d$sz))
  share$x <- x_max * 1.02
  share$label <- sprintf("%.0f %% beyond", 100 * share$z)
  cols <- c(variant_colours, stats::setNames("#3d4a42", ref_label))
  breaks_z <- c(-1000, -100, -10, 0, 10, 100, 1000)
  breaks_z <- breaks_z[abs(signed_log(breaks_z)) <= x_max]
  ggplot2::ggplot(d, ggplot2::aes(x = sz, y = variant, fill = variant)) +
    ggplot2::annotate(
      "rect",
      xmin = -signed_log(z_cut),
      xmax = signed_log(z_cut),
      ymin = -Inf,
      ymax = Inf,
      fill = "#e6e1d8",
      alpha = 0.7
    ) +
    ggplot2::geom_vline(xintercept = 0, colour = "#1c2b24", linewidth = 0.3) +
    ggdist::stat_slab(
      height = 0.9,
      normalize = "groups",
      colour = NA,
      alpha = 0.85
    ) +
    ggplot2::geom_text(
      data = share,
      ggplot2::aes(x = x, y = variant, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = -0.6,
      size = 2.5,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = signed_log(breaks_z),
      labels = scales::label_comma(accuracy = 1)(breaks_z),
      expand = ggplot2::expansion(mult = c(0.02, 0.32))
    ) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Standardised discrepancy z = D / sd(D) (signed log scale)",
      y = NULL
    ) +
    theme_slide()
}

#' Observed against fitted bulk, without and with the gene-wise ratio.
#'
#' Rows: variant of the cross-platform model; columns: time. Observed
#' bulk per cell `y / (a b)` against the fitted `sum_j p_j mu_j`, both
#' on `log10(1 + x)`, as 2D bins over genes and the three replicates,
#' with the identity line. Both axes share units, so the panels are
#' square (#58).
#'
#' @param variants Named list of variant results (elements `tails`).
#' @param which Variants to draw. Default `c("none", "gene")`.
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_bulk_fit <- function(variants, which = c("none", "gene"), title = NULL, subtitle = NULL) {
  d <- do.call(rbind, lapply(which, function(v) {
    do.call(rbind, lapply(names(variants[[v]]$tails), function(tm) {
      b <- variants[[v]]$tails[[tm]]$quantile_band
      data.frame(
        variant = variant_labels[[v]],
        time_point = tm,
        x = log10(1 + b$fitted),
        y = log10(1 + b$observed)
      )
    }))
  }))
  d$variant <- factor(d$variant, levels = unname(variant_labels[which]))
  stat <- do.call(rbind, lapply(split(d, list(d$variant, d$time_point), drop = TRUE), function(r) {
    data.frame(
      variant = r$variant[1L],
      time_point = r$time_point[1L],
      label = sprintf("median |log10 ratio| %.2f", stats::median(abs(r$y - r$x)))
    )
  }))
  top <- max(c(d$x, d$y))
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_bin_2d(bins = 60) +
    ggplot2::geom_abline(colour = "#1c2b24", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_text(
      data = stat,
      ggplot2::aes(x = 0, y = top, label = label),
      hjust = 0,
      vjust = 1,
      size = 2.6,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_fill_gradient(
      low = "#dfe7f1",
      high = "#1d4e89",
      transform = "log10",
      name = "Values per bin"
    ) +
    ggplot2::coord_equal(xlim = c(0, top), ylim = c(0, top)) +
    ggplot2::facet_wrap(
      ggplot2::vars(variant, time_point),
      ncol = length(unique(d$time_point)),
      labeller = ggplot2::label_wrap_gen(multi_line = FALSE)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Fitted bulk per cell, log10(1 + sum p mu)",
      y = "Observed bulk per cell, log10(1 + y / (a b))"
    ) +
    theme_slide()
}

#' Genes flagged by each statement, under the three models (dumbbell).
#'
#' Per time and statement, one point per variant on a common axis,
#' joined by a grey segment (#16); statements ranked by the count under
#' A3 (#50). One legend for the three variants: direct labels overlap
#' where the points of a row are close.
#'
#' @param variants Named list of variant results (elements `overlaps`,
#'   `range`).
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_tail_flags <- function(variants, title = NULL, subtitle = NULL) {
  statement_labels <- c(
    "T1" = "T1 (normal)",
    "T1-log" = "T1-log (normal, log scale)",
    "T3" = "T3 (Poisson Chernoff)",
    "T4" = "T4 (NB Chernoff)",
    "QR" = "Quantile band (1 % / 99 %)",
    "range" = "Range (RNA-Sieve, p-free)",
    "mad" = "MAD (RNA-Sieve, p-free)"
  )
  d <- do.call(rbind, lapply(names(variants), function(v) {
    do.call(rbind, lapply(names(variants[[v]]$overlaps), function(tm) {
      o <- variants[[v]]$overlaps[[tm]]$order
      r <- variants[[v]]$range[[tm]]
      rbind(
        data.frame(
          variant = variant_labels[[v]],
          time_point = tm,
          statement = o$statement,
          n = o$n_flagged
        ),
        data.frame(
          variant = variant_labels[[v]],
          time_point = tm,
          statement = c("range", "mad"),
          n = c(sum(r$flag_range), sum(r$flag_mad))
        )
      )
    }))
  }))
  d <- d[d$statement %in% names(statement_labels), , drop = FALSE]
  base <- stats::aggregate(
    n ~ statement,
    data = d[d$variant == variant_labels[["none"]], , drop = FALSE],
    FUN = mean
  )
  d$statement <- factor(
    statement_labels[d$statement],
    levels = statement_labels[base$statement[order(base$n)]]
  )
  d$variant <- factor(d$variant, levels = unname(variant_labels))
  seg <- stats::aggregate(n ~ time_point + statement, data = d, FUN = range)
  seg <- data.frame(
    time_point = seg$time_point,
    statement = seg$statement,
    lo = seg$n[, 1L],
    hi = seg$n[, 2L]
  )
  ggplot2::ggplot(d, ggplot2::aes(x = n, y = statement)) +
    ggplot2::geom_segment(
      data = seg,
      ggplot2::aes(x = lo, xend = hi, y = statement, yend = statement),
      inherit.aes = FALSE,
      colour = "grey75",
      linewidth = 1
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = variant), size = 2.2) +
    ggplot2::scale_colour_manual(values = variant_colours, name = NULL) +
    ggplot2::scale_x_continuous(
      labels = scales::label_comma(),
      breaks = scales::breaks_extended(4),
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0.02, 0.08))
    ) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Genes flagged in at least one replicate",
      y = NULL
    ) +
    theme_slide() +
    ggplot2::theme(legend.position = "top", legend.justification = "left")
}

#' Plug-in composition per bulk replicate, without and with the ratio.
#'
#' Dot plot (#57: no stacked bars): one point per replicate and type,
#' types ranked by their mean share, one facet per time, grey for the A3
#' fit and orange for the fit with the gene-wise ratio.
#'
#' @param variants Named list of variant results (elements `tails`).
#' @param which Variants to draw. Default `c("none", "gene")`.
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_plugin_composition <- function(variants, which = c("none", "gene"), title = NULL, subtitle = NULL) {
  d <- do.call(rbind, lapply(which, function(v) {
    do.call(rbind, lapply(names(variants[[v]]$tails), function(tm) {
      cmp <- variants[[v]]$tails[[tm]]$plugin$composition
      data.frame(
        variant = variant_labels[[v]],
        time_point = tm,
        cell_type = cmp$cell_type,
        p_hat = cmp$p_hat
      )
    }))
  }))
  key <- function(tm, ct) paste(tm, ct, sep = "___")
  ord <- stats::aggregate(p_hat ~ time_point + cell_type, data = d, FUN = mean)
  ord <- ord[order(ord$time_point, ord$p_hat), , drop = FALSE]
  d$row <- factor(key(d$time_point, d$cell_type), levels = key(ord$time_point, ord$cell_type))
  d$variant <- factor(d$variant, levels = unname(variant_labels[which]))
  ggplot2::ggplot(d, ggplot2::aes(x = p_hat, y = row, colour = variant)) +
    ggplot2::geom_point(
      size = 1.8,
      alpha = 0.8,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::scale_colour_manual(values = variant_colours, name = NULL) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(),
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      expand = ggplot2::expansion(mult = c(0.03, 0.08))
    ) +
    ggplot2::scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::facet_wrap(~time_point, nrow = 1, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Plug-in composition (deMixSC), one point per replicate",
      y = NULL
    ) +
    theme_slide()
}

#' Gene-wise bulk-to-reference ratio, one time against another.
#'
#' Each panel plots the mean log2 ratio of a gene (A3 fit,
#' [platform_log_ratio()]) at one matched time against another. Genes
#' flagged by T4 at two or three times with one sign are orange, the
#' others grey; both groups are counted in the first panel (#41). Under
#' A3 the points would gather at (0, 0); a stable gene-wise platform
#' effect puts them along the diagonal.
#'
#' @param log_ratio Output of [platform_log_ratio()].
#' @param drift Output of [platform_drift_summary()].
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_platform_ratio <- function(log_ratio, drift, title = NULL, subtitle = NULL) {
  m <- log_ratio$log_ratio
  times <- colnames(m)
  pairs <- utils::combn(times, 2L, simplify = FALSE)
  cls <- drift$summary$drift_class[match(rownames(m), drift$summary$gene)]
  stable <- cls %in% "several times, same sign"
  d <- do.call(
    rbind,
    lapply(pairs, function(pr) {
      rho <- log_ratio$cor[pr[1L], pr[2L]]
      data.frame(
        panel = sprintf("%s against %s (Spearman %.2f)", pr[2L], pr[1L], rho),
        x = m[, pr[1L]],
        y = m[, pr[2L]],
        stable = stable
      )
    })
  )
  d <- d[stats::complete.cases(d[, c("x", "y")]), , drop = FALSE]
  d <- d[order(d$stable), , drop = FALSE]
  lim <- range(c(d$x, d$y))
  first <- sort(unique(d$panel))[1L]
  lab <- data.frame(
    panel = first,
    x = lim[1L],
    y = lim[2L] - c(0, 0.09) * diff(lim),
    stable = c(TRUE, FALSE),
    label = c(
      paste0("T4 at 2-3 times, one sign: ", fmt_n(sum(stable)), " genes"),
      paste0("other genes: ", fmt_n(sum(!stable)))
    )
  )
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = stable)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#b8b1a5") +
    ggplot2::geom_vline(xintercept = 0, colour = "#b8b1a5") +
    ggplot2::geom_abline(colour = "#1c2b24", linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_point(size = 0.35, alpha = 0.3) +
    ggplot2::geom_text(
      data = lab,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 0,
      vjust = 1,
      size = 2.6,
      fontface = "bold"
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#c45c26", `FALSE` = "#9a9286"),
      guide = "none"
    ) +
    ggplot2::coord_equal(xlim = lim, ylim = lim) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "log2 bulk / reference ratio, earlier time",
      y = "log2 ratio, later time"
    ) +
    theme_slide()
}
