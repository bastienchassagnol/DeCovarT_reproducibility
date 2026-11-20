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

#' QQ plot of the null p-values of the four tests (qqplotr).
#'
#' Under the null a p-value is uniform, so `-log10 p` is exponential with
#' rate `log(10)`. Each panel (test x time) plots the observed `-log10 p`
#' against that exponential with `qqplotr`: points, the identity line and
#' a pointwise 95 % band. Points above the band mean too many small
#' p-values (anti-conservative), points below it too few (conservative).
#' The band assumes independent genes, so it is too narrow when genes are
#' correlated.
#'
#' @param null_tests Named list (times) of data frames with `<test>_pval`.
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_null_calibration <- function(null_tests, title = NULL, subtitle = NULL) {
  d <- do.call(
    rbind,
    lapply(names(null_tests), function(tm) {
      do.call(rbind, lapply(names(test_labels), function(tt) {
        p <- null_tests[[tm]][[paste0(tt, "_pval")]]
        p <- p[!is.na(p) & p > 0]
        if (length(p) == 0L) {
          return(NULL)
        }
        data.frame(time_point = tm, test = test_labels[[tt]], mlog_p = -log10(p))
      }))
    })
  )
  d$test <- factor(d$test, levels = unname(test_labels))
  dp <- list(rate = log(10))
  # Square panels on the range of the expected values, x 1.6; the few
  # larger observed values are drawn at the top edge (one ANOVA gene).
  top <- 1.6 * -log10(1 / (max(table(d$test, d$time_point)) + 1))
  d$mlog_p <- pmin(d$mlog_p, top)
  ggplot2::ggplot(d, ggplot2::aes(sample = mlog_p)) +
    qqplotr::stat_qq_band(
      distribution = "exp",
      dparams = dp,
      bandType = "pointwise",
      fill = "#d9d3c7",
      alpha = 0.8
    ) +
    qqplotr::stat_qq_line(
      distribution = "exp",
      dparams = dp,
      identity = TRUE,
      colour = "#1c2b24",
      linewidth = 0.3
    ) +
    qqplotr::stat_qq_point(
      ggplot2::aes(colour = time_point),
      distribution = "exp",
      dparams = dp,
      size = 0.4
    ) +
    ggplot2::scale_colour_manual(values = time_colours, guide = "none") +
    ggplot2::coord_equal(xlim = c(0, top), ylim = c(0, top)) +
    ggplot2::facet_grid(time_point ~ test) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Expected -log10 p under uniformity",
      y = "Observed -log10 p, permuted labels"
    ) +
    theme_slide()
}

#' Histograms of the unadjusted p-values.
#'
#' One panel per test (rows) and time (columns), 20 bins of width 0.05.
#' The dashed line is the height a uniform (null) histogram would have,
#' `pi0 * m / 20`, with `pi0` from `qvalue::pi0est()`: bars above it near
#' 0 are the signal, a flat remainder is the null part; a peak near 1 or
#' gaps flag a discrete or conservative test.
#'
#' @param tests Named list (times) of [compare_type_effect_tests()] outputs.
#' @param strata Optional named list (times) of factors, one per gene of
#'   `tests[[tm]]` (for example an expression tertile). When given, only
#'   `test` is drawn and the rows are the strata.
#' @param test Test drawn when `strata` is given. Default `"welch_perm"`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_pvalue_histograms <- function(
  tests,
  strata = NULL,
  test = "welch_perm",
  title = NULL,
  subtitle = NULL
) {
  tt_set <- if (is.null(strata)) names(test_labels) else test
  d <- do.call(rbind, lapply(names(tests), function(tm) {
    r <- tests[[tm]]
    keep <- r$pass_presence %in% TRUE
    do.call(rbind, lapply(tt_set, function(tt) {
      p <- r[[paste0(tt, "_pval")]][keep]
      g <- if (is.null(strata)) rep(test_short[[tt]], sum(keep)) else as.character(strata[[tm]][keep])
      ok <- !is.na(p)
      data.frame(time_point = tm, row = g[ok], p = p[ok])
    }))
  }))
  row_levels <- if (is.null(strata)) unname(test_short) else levels(strata[[1L]])
  d$row <- factor(d$row, levels = row_levels)
  ref <- do.call(rbind, lapply(split(d, list(d$row, d$time_point), drop = TRUE), function(r) {
    pi0 <- tryCatch(qvalue::pi0est(r$p)$pi0, error = function(e) NA_real_)
    data.frame(
      row = r$row[1L],
      time_point = r$time_point[1L],
      height = pi0 * nrow(r) / 20,
      pi0 = pi0
    )
  }))
  ggplot2::ggplot(d, ggplot2::aes(x = p)) +
    ggplot2::geom_histogram(
      breaks = seq(0, 1, by = 0.05),
      fill = "#d9d3c7",
      colour = "white"
    ) +
    ggplot2::geom_hline(
      data = ref,
      ggplot2::aes(yintercept = height),
      linetype = "dashed",
      colour = "#1c2b24"
    ) +
    ggplot2::geom_text(
      data = ref,
      ggplot2::aes(x = 0.98, y = height, label = sprintf("pi0 = %.2f", pi0)),
      hjust = 1,
      vjust = -0.4,
      size = 2.6,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1)) +
    ggplot2::facet_grid(row ~ time_point, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Unadjusted p-value (bins of 0.05)",
      y = "Genes"
    ) +
    theme_slide()
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

#' Gate class of each tested gene: fails BH, passes BH, or at the floor.
#'
#' @param r One [rank_genes_omnibus()] output.
#' @param n_perm Number of permutations B.
#' @return Factor with levels `"Fails BH"`, `"Passes BH 0.05"`,
#'   `"At the floor 1/(B+1)"`.
gate_class <- function(r, n_perm) {
  floor_p <- 1 / (n_perm + 1)
  cls <- ifelse(
    r$pass_fdr %in% TRUE,
    ifelse(r$pval <= floor_p + 1e-12, "At the floor 1/(B+1)", "Passes BH 0.05"),
    "Fails BH"
  )
  factor(cls, levels = names(gate_colours))
}
gate_colours <- c(
  "Fails BH" = "#9a9286",
  "Passes BH 0.05" = "#c45c26",
  "At the floor 1/(B+1)" = "#1d4e89"
)

#' Gate and rank: deviance explained against the gate p-value.
#'
#' One point per gene that passes the presence gate, coloured by
#' [gate_class()]. The permutation p-value takes the values k / (B + 1),
#' so the genes form vertical stripes (jittered); the rightmost stripe is
#' the floor k = 1. The dashed line is the largest raw p-value that passes
#' BH at that time, the dotted line the floor. The formula of the p-value
#' is written in the first panel.
#'
#' @param ranked Named list (times) of [rank_genes_omnibus()] outputs.
#' @param n_perm Number of permutations B.
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_gate_rank <- function(ranked, n_perm, title = NULL, subtitle = NULL) {
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
        status = gate_class(r, n_perm)
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
  # Draw the larger class first so the other stays visible (#55, #61).
  d <- d[order(-table(d$status)[d$status]), , drop = FALSE]
  floor_x <- -log10(1 / (n_perm + 1))
  first <- names(ranked)[1L]
  top <- max(d$r2, na.rm = TRUE)
  line_labels <- data.frame(
    time_point = first,
    x = c(cuts$cut[cuts$time_point == first], floor_x),
    y = top * 0.9^2,
    label = c("BH 0.05", "floor"),
    hjust = c(1.05, 1.05)
  )
  counts <- stats::aggregate(gene ~ time_point + status, d, length)
  class_labels <- do.call(rbind, lapply(split(counts, counts$time_point), function(r) {
    r <- r[match(levels(d$status), r$status), , drop = FALSE]
    r <- r[!is.na(r$gene), , drop = FALSE]
    data.frame(
      time_point = r$time_point,
      x = 0,
      # Positions on the square-root axis: 80, 68 and 56 % of its height.
      y = top * c(0.8, 0.68, 0.56)[seq_len(nrow(r))]^2,
      label = paste0(r$status, ": ", fmt_n(r$gene)),
      status = r$status
    )
  }))
  formula_label <- data.frame(
    time_point = first,
    x = 0,
    y = top,
    label = latex2exp::TeX(
      "$p_g = (1 + |\\{b : G^*_{gb} \\geq G_g\\}|) / (B + 1)$",
      output = "character"
    )
  )
  ggplot2::ggplot(d, ggplot2::aes(x = mlog_p, y = r2, colour = status)) +
    ggplot2::geom_jitter(width = 0.04, height = 0, size = 0.35, alpha = 0.3) +
    ggplot2::geom_vline(
      data = cuts,
      ggplot2::aes(xintercept = cut),
      linetype = "dashed",
      colour = "#c45c26"
    ) +
    ggplot2::geom_vline(xintercept = floor_x, linetype = "dotted", colour = "#1d4e89") +
    ggplot2::geom_text(
      data = line_labels,
      ggplot2::aes(x = x, y = y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      size = 2.8,
      colour = "#1c2b24",
      vjust = 1
    ) +
    ggplot2::geom_text(
      data = formula_label,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 3,
      colour = "#1c2b24"
    ) +
    ggplot2::geom_text(
      data = class_labels,
      ggplot2::aes(x = x, y = y, label = label, colour = status),
      hjust = 0,
      size = 2.8,
      fontface = "bold"
    ) +
    ggplot2::scale_colour_manual(values = gate_colours, guide = "none") +
    ggplot2::scale_y_sqrt() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "-log10 p, permuted Welch G (jittered)",
      y = "Between-type share of deviance R2_D (sqrt scale)"
    ) +
    theme_slide() +
    ggplot2::theme(plot.margin = ggplot2::margin(18, 8, 6, 6))
}

#' Effect size against mean expression, coloured by gate class.
#'
#' @param ranked Named list (times) of [rank_genes_omnibus()] outputs.
#' @param gene_means Named list (times) of named vectors of mean
#'   size-normalised counts.
#' @param n_perm Number of permutations B.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_r2_vs_mean <- function(ranked, gene_means, n_perm, title = NULL, subtitle = NULL) {
  d <- do.call(rbind, lapply(names(ranked), function(tm) {
    r <- ranked[[tm]]
    r <- r[r$pass_presence & !is.na(r$pval), , drop = FALSE]
    data.frame(
      time_point = tm,
      mean = gene_means[[tm]][r$gene],
      r2 = r$deviance_r2,
      status = gate_class(r, n_perm)
    )
  }))
  d <- d[d$mean > 0, , drop = FALSE]
  d <- d[order(-table(d$status)[d$status]), , drop = FALSE]
  ggplot2::ggplot(d, ggplot2::aes(x = mean, y = r2, colour = status)) +
    ggplot2::geom_point(size = 0.35, alpha = 0.3) +
    ggplot2::scale_colour_manual(
      values = gate_colours,
      name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))
    ) +
    ggplot2::scale_x_log10(labels = scales::label_log()) +
    ggplot2::scale_y_sqrt() +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Mean size-normalised count",
      y = "R2_D (sqrt scale)"
    ) +
    theme_slide()
}

#' Gene-wise overdispersion against mean expression.
#'
#' `glmGamPoi` returns phi = 1/psi = 0 when the likelihood is maximal at
#' the Poisson boundary (no extra-Poisson variation detected) and very
#' large values when it diverges (counts in a handful of cells). Three
#' classes are coloured: Poisson-like genes (phi below `poisson_cut`,
#' zeros drawn on a labelled strip at the bottom), the upper boundary
#' (phi above `upper_cut`) and the others, with a loess trend on the
#' others. The NB variance is written in the first panel.
#'
#' @param tests Named list (times) of [compare_type_effect_tests()]
#'   outputs.
#' @param gene_means Named list (times) of named vectors: mean
#'   size-normalised count per gene.
#' @param poisson_cut,upper_cut Class limits on phi. Defaults `1e-4`,
#'   `1e10`.
#' @param title,subtitle Plot title (the finding) and subtitle.
#' @return A ggplot.
plot_overdispersion <- function(
  tests,
  gene_means,
  poisson_cut = 1e-4,
  upper_cut = 1e10,
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
  d <- d[d$mean > 0, , drop = FALSE]
  od_colours <- c(
    "Poisson-like" = "#5b8c5a",
    "Upper boundary" = "#9a9286",
    "Overdispersed" = "#1d4e89"
  )
  d$class <- factor(
    ifelse(
      d$overdispersion < poisson_cut,
      "Poisson-like",
      ifelse(d$overdispersion > upper_cut, "Upper boundary", "Overdispersed")
    ),
    levels = names(od_colours)
  )
  pos <- d$overdispersion[d$overdispersion > 0]
  zero_y <- 10^(floor(log10(min(pos))) - 1)
  d$y <- ifelse(d$overdispersion == 0, zero_y, d$overdispersion)
  counts <- do.call(rbind, lapply(split(d, d$time_point), function(r) {
    data.frame(
      time_point = r$time_point[1L],
      x = min(r$mean),
      y = c(poisson_cut, zero_y, upper_cut),
      label = c(
        paste0("Poisson-like (phi < 1e-4): ", fmt_n(sum(r$class == "Poisson-like"))),
        paste0("phi = 0: ", fmt_n(sum(r$overdispersion == 0))),
        paste0("upper boundary: ", fmt_n(sum(r$class == "Upper boundary")))
      ),
      class = c("Poisson-like", "Poisson-like", "Upper boundary")
    )
  }))
  formula_label <- data.frame(
    time_point = names(tests)[1L],
    x = max(d$mean),
    y = upper_cut,
    label = latex2exp::TeX(
      "$Var(x_{gc}) = m_{gc} + \\phi_g m_{gc}^2$, $\\phi_g = 1/\\psi_g$",
      output = "character"
    )
  )
  ggplot2::ggplot(d, ggplot2::aes(x = mean, y = y, colour = class)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.3) +
    ggplot2::geom_smooth(
      data = d[d$class == "Overdispersed", , drop = FALSE],
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      colour = "#c45c26",
      linewidth = 0.7
    ) +
    ggplot2::geom_hline(yintercept = poisson_cut, linetype = "dotted", colour = "#1c2b24") +
    ggplot2::geom_text(
      data = counts,
      ggplot2::aes(x = x, y = y, label = label, colour = class),
      hjust = 0,
      vjust = -0.5,
      size = 2.7
    ) +
    ggplot2::geom_text(
      data = formula_label,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      parse = TRUE,
      hjust = 1,
      vjust = -0.5,
      size = 3,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_colour_manual(values = od_colours, guide = "none") +
    ggplot2::scale_x_log10(labels = scales::label_log()) +
    ggplot2::scale_y_log10(labels = scales::label_log()) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Mean size-normalised count",
      y = "Overdispersion phi = 1/psi (M1, ML)"
    ) +
    theme_slide()
}
