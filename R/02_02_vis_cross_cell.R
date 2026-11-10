#' Tier (ii) figures: calibration, agreement and the gate-and-rank screen.
#'
#' Plots only; nothing is saved here (see
#' `scripts/02_02_vis_cross_cell.R`). Inputs are elements of
#' `output/feature_selection/tier_02/cross_cell.rds`. Requires
#' `R/02_00_utils.R` (`theme_slide()`).

test_labels <- c(
  lr = "Likelihood ratio",
  anova_vst = "ANOVA, shifted log",
  ql_f = "Quasi-likelihood F",
  welch_perm = "Permuted Welch G"
)
test_short <- c(
  lr = "LR",
  anova_vst = "ANOVA",
  ql_f = "QL F",
  welch_perm = "Welch G"
)
time_colours <- c("48h" = "#1d4e89", "72h" = "#c45c26", "96h" = "#5b8c5a")

#' QQ plot of the null p-values of the four tests.
#'
#' Observed against expected `-log10 p` under one within-batch
#' permutation of the labels, one panel per test, one line per time,
#' labelled at its right end. The grey band is the pointwise 95 %
#' interval of uniform order statistics; a curve above it is
#' anti-conservative. Both axes are `-log10 p`, so the panels are square
#' with a common range.
#'
#' @param null_tests Named list (times) of data frames with `<test>_pval`.
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_null_calibration <- function(null_tests, title = NULL, subtitle = NULL) {
  rows <- list()
  bands <- list()
  for (tm in names(null_tests)) {
    for (tt in names(test_labels)) {
      p <- null_tests[[tm]][[paste0(tt, "_pval")]]
      p <- sort(p[!is.na(p)])
      if (length(p) == 0L) {
        next
      }
      n <- length(p)
      i <- seq_len(n)
      rows[[length(rows) + 1L]] <- data.frame(
        time_point = tm,
        test = test_labels[[tt]],
        expected = -log10(i / (n + 1)),
        observed = -log10(p)
      )
      if (tm == names(null_tests)[1L]) {
        bands[[length(bands) + 1L]] <- data.frame(
          test = test_labels[[tt]],
          expected = -log10(i / (n + 1)),
          lo = -log10(stats::qbeta(0.975, i, n - i + 1)),
          hi = -log10(stats::qbeta(0.025, i, n - i + 1))
        )
      }
    }
  }
  d <- do.call(rbind, rows)
  b <- do.call(rbind, bands)
  d$test <- factor(d$test, levels = unname(test_labels))
  b$test <- factor(b$test, levels = unname(test_labels))
  top <- max(c(d$expected, d$observed), na.rm = TRUE)
  ends <- do.call(rbind, lapply(split(d, list(d$test, d$time_point), drop = TRUE), function(r) {
    r[which.max(r$expected), , drop = FALSE]
  }))
  ends <- do.call(rbind, lapply(split(ends, ends$test), function(r) {
    r$label_y <- spread_labels(r$observed, min_gap = 0.07 * top)
    r
  }))
  ggplot2::ggplot(d, ggplot2::aes(x = expected, y = observed)) +
    ggplot2::geom_ribbon(
      data = b,
      ggplot2::aes(x = expected, ymin = lo, ymax = hi),
      inherit.aes = FALSE,
      fill = "#e6e1d8"
    ) +
    ggplot2::geom_abline(colour = "#1c2b24", linewidth = 0.3) +
    ggplot2::geom_line(ggplot2::aes(colour = time_point), linewidth = 0.6) +
    ggplot2::geom_text(
      data = ends,
      ggplot2::aes(y = label_y, label = time_point, colour = time_point),
      hjust = 0,
      nudge_x = 0.08 * top,
      size = 3
    ) +
    ggplot2::scale_colour_manual(values = time_colours, guide = "none") +
    ggplot2::coord_equal(xlim = c(0, top), ylim = c(0, top), clip = "off") +
    ggplot2::facet_wrap(~test, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Expected -log10 p (uniform)",
      y = "Observed -log10 p"
    ) +
    theme_slide() +
    ggplot2::theme(panel.spacing = ggplot2::unit(2.2, "lines"))
}

#' Kendall's tau between the gene orders of the four tests.
#'
#' Lower triangle only (the matrix is symmetric, the diagonal is 1);
#' every cell is labelled with its value, so the colour scale has no
#' legend.
#'
#' @param tau Named list (times) of test-by-test tau matrices.
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot (one tile grid per time).
plot_test_agreement <- function(tau, title = NULL, subtitle = NULL) {
  lev <- names(test_short)
  d <- do.call(
    rbind,
    lapply(names(tau), function(tm) {
      m <- tau[[tm]][lev, lev]
      keep <- row(m) > col(m)
      data.frame(
        time_point = tm,
        test_a = factor(test_short[colnames(m)[col(m)[keep]]], levels = unname(test_short)),
        test_b = factor(test_short[rownames(m)[row(m)[keep]]], levels = rev(unname(test_short))),
        tau = m[keep]
      )
    })
  )
  ggplot2::ggplot(d, ggplot2::aes(x = test_a, y = test_b, fill = tau)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", tau), colour = tau > 0.85),
      size = 3.2
    ) +
    ggplot2::scale_fill_gradient(
      low = "#f7f4ee",
      high = "#1d4e89",
      limits = c(0.5, 1),
      oob = scales::squish,
      guide = "none"
    ) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "#1c2b24", `TRUE` = "white"),
      guide = "none"
    ) +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::scale_y_discrete(drop = TRUE) +
    ggplot2::coord_equal() +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_slide() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank())
}

#' Gate and rank: deviance explained against the gate p-value.
#'
#' One point per gene that passes the presence gate: orange when it
#' passes Benjamini-Hochberg (and so enters the panel, `n_select = Inf`),
#' grey otherwise; a third class "passes BH, not selected" appears only
#' with a finite `n_select`. The dashed line is the largest raw p-value
#' that passes BH, the dotted line the permutation floor `1 / (B + 1)`;
#' both are labelled in the first panel, and so are the classes.
#'
#' @param ranked Named list (times) of [rank_genes_omnibus()] outputs.
#' @param n_perm Number of permutations (for the floor).
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_gate_rank <- function(ranked, n_perm = 2000L, title = NULL, subtitle = NULL) {
  d <- do.call(
    rbind,
    lapply(names(ranked), function(tm) {
      r <- ranked[[tm]]
      r <- r[r$pass_presence & !is.na(r$pval), , drop = FALSE]
      data.frame(
        time_point = tm,
        gene = r$gene,
        mlog_p = -log10(r$pval),
        r2 = r$deviance_r2,
        status = ifelse(
          r$selected,
          "Passes BH: in the panel",
          ifelse(r$pass_fdr, "Passes BH, not selected", "Fails BH")
        )
      )
    })
  )
  cuts <- do.call(
    rbind,
    lapply(names(ranked), function(tm) {
      r <- ranked[[tm]]
      p <- r$pval[r$pass_fdr %in% TRUE]
      data.frame(
        time_point = tm,
        cut = if (length(p) > 0L) -log10(max(p)) else NA_real_
      )
    })
  )
  class_colours <- c(
    "Fails BH" = "#9a9286",
    "Passes BH, not selected" = "#1d4e89",
    "Passes BH: in the panel" = "#c45c26"
  )
  d$status <- droplevels(factor(d$status, levels = names(class_colours)))
  # Draw the larger class first so the other stays visible (#55, #61).
  d <- d[order(-table(d$status)[d$status]), , drop = FALSE]
  floor_x <- -log10(1 / (n_perm + 1))
  first <- names(ranked)[1L]
  top <- max(d$r2, na.rm = TRUE)
  line_labels <- data.frame(
    time_point = first,
    x = c(cuts$cut[cuts$time_point == first], floor_x),
    y = top,
    label = c("BH 0.05", "floor 1/(B+1)"),
    hjust = c(1.05, 1.05)
  )
  class_labels <- do.call(rbind, lapply(levels(d$status), function(s) {
    r <- d[d$time_point == first & d$status == s, , drop = FALSE]
    data.frame(
      time_point = first,
      x = min(r$mlog_p) + 0.05,
      y = stats::quantile(r$r2, 0.995, names = FALSE),
      label = sub(": ", "\n", s, fixed = TRUE),
      status = s
    )
  }))
  ggplot2::ggplot(d, ggplot2::aes(x = mlog_p, y = r2, colour = status)) +
    ggplot2::geom_jitter(
      width = 0.06,
      height = 0,
      size = 0.35,
      alpha = 0.25
    ) +
    ggplot2::geom_vline(
      data = cuts,
      ggplot2::aes(xintercept = cut),
      linetype = "dashed",
      colour = "#1c2b24"
    ) +
    ggplot2::geom_vline(
      xintercept = floor_x,
      linetype = "dotted",
      colour = "#1c2b24"
    ) +
    ggplot2::geom_text(
      data = line_labels,
      ggplot2::aes(x = x, y = y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      size = 2.8,
      colour = "#1c2b24",
      vjust = 1
    ) +
    ggplot2::geom_label(
      data = class_labels,
      ggplot2::aes(x = x, y = y, label = label, colour = status),
      hjust = 0,
      lineheight = 0.9,
      size = 2.9,
      fill = "white",
      linewidth = 0,
      fontface = "bold"
    ) +
    ggplot2::scale_colour_manual(values = class_colours, guide = "none") +
    ggplot2::scale_y_sqrt() +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "-log10 p, permuted Welch G (jittered)",
      y = "Deviance explained R2_D (square-root scale)"
    ) +
    theme_slide()
}

#' Gene-wise overdispersion against mean expression.
#'
#' Small transparent points with a loess trend (orange); the genes at
#' the Poisson limit of the fit (overdispersion below `poisson_cut`) are
#' counted and labelled in each panel.
#'
#' @param tests Named list (times) of [compare_type_effect_tests()]
#'   outputs.
#' @param gene_means Named list (times) of named vectors: mean
#'   size-normalised count per gene.
#' @param poisson_cut Overdispersion below which a gene is reported as
#'   at the Poisson limit. Default `1e-4`.
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_overdispersion <- function(
  tests,
  gene_means,
  poisson_cut = 1e-4,
  title = NULL,
  subtitle = NULL
) {
  d <- do.call(
    rbind,
    lapply(names(tests), function(tm) {
      r <- tests[[tm]]
      r <- r[r$pass_presence & !is.na(r$overdispersion), , drop = FALSE]
      data.frame(
        time_point = tm,
        mean = gene_means[[tm]][r$gene],
        overdispersion = r$overdispersion
      )
    })
  )
  d <- d[d$overdispersion > 0 & d$mean > 0, , drop = FALSE]
  low <- do.call(rbind, lapply(split(d, d$time_point), function(r) {
    data.frame(
      time_point = r$time_point[1L],
      x = min(r$mean),
      y = poisson_cut,
      label = paste0("below: ", fmt_n(sum(r$overdispersion < poisson_cut)), " genes at the Poisson limit")
    )
  }))
  ggplot2::ggplot(d, ggplot2::aes(x = mean, y = overdispersion)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.2, colour = "#1d4e89") +
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      colour = "#c45c26",
      linewidth = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = poisson_cut,
      linetype = "dotted",
      colour = "#1c2b24"
    ) +
    ggplot2::geom_text(
      data = low,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 0,
      vjust = -0.5,
      size = 2.8,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_x_log10(labels = scales::label_log()) +
    ggplot2::scale_y_log10(labels = scales::label_log()) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Mean size-normalised count",
      y = "Overdispersion 1/psi (M1)"
    ) +
    theme_slide()
}
