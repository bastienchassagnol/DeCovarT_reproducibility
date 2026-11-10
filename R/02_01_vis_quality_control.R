#' Tier (i) figures: projections, low-quality cells, gene filters.
#'
#' Plots only; nothing is saved here (see
#' `scripts/02_01_vis_quality_control.R`). Requires `R/02_00_utils.R`
#' (`theme_slide()`).

#' Pearson residuals of the null count model M0.
#'
#' `log m_gc = beta_0g + log s_c`, negative binomial with a gene-wise
#' maximum-likelihood overdispersion, fitted with `glmGamPoi::glm_gp()`
#' as in [compare_type_effect_tests()]. The residual is
#' `(x - m) / sqrt(m + phi m^2)`; nothing is clipped. The type labels do
#' not enter the fit, so the type structure is left in the residuals.
#'
#' @param counts Genes by cells raw counts.
#' @param size_factors Size factor per cell.
#' @return Dense genes by cells matrix of Pearson residuals.
pearson_residuals_m0 <- function(counts, size_factors) {
  if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
    stop("Package glmGamPoi is required for the Pearson residuals.")
  }
  y <- as.matrix(counts)
  fit <- glmGamPoi::glm_gp(
    y,
    design = ~1,
    offset = log(size_factors),
    size_factors = rep(1, ncol(y)),
    overdispersion = TRUE,
    overdispersion_shrinkage = FALSE
  )
  m <- fit$Mu
  (y - m) / sqrt(m + fit$overdispersions * m^2)
}

#' PCA and PLS-DA of one time slot.
#'
#' Descriptive projections of the cells that enter tier (ii); not a
#' selection rule. Input (`input`):
#'
#' - `"pearson"` (default): Pearson residuals of M0
#'   ([pearson_residuals_m0()]), the variance-stabilised scale of the
#'   tier (ii) count model, without `SCTransform` and without clipping.
#'   Genes are centred, not scaled: scaling would undo the variance
#'   stabilisation and up-weight genes that carry mostly noise.
#' - `"shifted_log"`: `log1p(counts / size_factor)`, centred and scaled.
#'
#' PLS-DA is `mixOmics::plsda()` on the same matrix (same scaling
#' choice), two components, with the cell type as the response. The
#' projections are computed once; `groups` adds the same scatters
#' coloured by another labelling (for example the cell line).
#'
#' @param counts Genes by cells, symbols as row names.
#' @param cell_type Cell-type label per column.
#' @param size_factors Size factor per cell. Default: column sums over
#'   their median.
#' @param colours Optional named colours for the cell types.
#' @param ncomp Components requested from PLS-DA. Default `2`.
#' @param input `"pearson"` or `"shifted_log"`.
#' @param groups Optional named list of extra labellings (one value per
#'   cell), each drawn as `pca_<name>` and `plsda_<name>`.
#' @param group_colours Optional named list of named colour vectors, one
#'   per entry of `groups`.
#' @return A list with `pca`, `plsda` (ggplots coloured by type), one
#'   pair per entry of `groups`, `n_genes`, and `scores` (data frame of
#'   the PCA and PLS-DA coordinates per cell).
plot_type_pca_plsda <- function(
  counts,
  cell_type,
  size_factors = NULL,
  colours = NULL,
  ncomp = 2L,
  input = c("pearson", "shifted_log"),
  groups = NULL,
  group_colours = NULL
) {
  input <- match.arg(input)
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package mixOmics is required for PLS-DA.")
  }
  if (length(cell_type) != ncol(counts)) {
    stop("`cell_type` must have one label per column of `counts`.")
  }
  cell_type <- as.character(cell_type)
  if (anyNA(cell_type) || any(!nzchar(cell_type))) {
    stop("`cell_type` has missing labels.")
  }
  if (length(unique(cell_type)) < 2L) {
    stop("PLS-DA needs at least two cell types.")
  }
  if (ncomp < 2L) {
    stop("`ncomp` must be at least 2 so both scatters have two axes.")
  }
  if (is.null(size_factors)) {
    lib <- Matrix::colSums(counts)
    size_factors <- lib / stats::median(lib)
  }
  if (
    length(size_factors) != ncol(counts) ||
      any(!is.finite(size_factors) | size_factors <= 0)
  ) {
    stop("`size_factors` must be positive and finite, one per column.")
  }
  scale_genes <- input == "shifted_log"
  x <- if (input == "pearson") {
    t(pearson_residuals_m0(counts, size_factors))
  } else {
    t(log1p(sweep(as.matrix(counts), 2L, size_factors, "/")))
  }
  gene_sd <- apply(x, 2L, stats::sd)
  x <- x[, is.finite(gene_sd) & gene_sd > 0, drop = FALSE]
  if (nrow(x) < 3L || ncol(x) < 2L) {
    stop("PCA and PLS-DA need at least three cells and two varying genes.")
  }
  pc <- stats::prcomp(x, center = TRUE, scale. = scale_genes, rank. = 2L)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  fit <- mixOmics::plsda(
    X = x,
    Y = factor(cell_type),
    ncomp = ncomp,
    scale = scale_genes
  )
  pls <- as.matrix(fit$variates$X)
  explained <- fit$prop_expl_var$X
  scores <- data.frame(
    pc1 = pc$x[, 1L],
    pc2 = pc$x[, 2L],
    pls1 = pls[, 1L],
    pls2 = pls[, 2L],
    cell_type = cell_type,
    stringsAsFactors = FALSE
  )
  lab <- function(prefix, k, share) {
    paste0(prefix, k, " (", format(round(share, 1), nsmall = 1), "%)")
  }
  pca_x <- lab("PC", 1, pct[[1L]])
  pca_y <- lab("PC", 2, pct[[2L]])
  pls_x <- lab("Component ", 1, 100 * explained[[1L]])
  pls_y <- lab("Component ", 2, 100 * explained[[2L]])
  input_label <- if (input == "pearson") {
    "Pearson residuals of M0"
  } else {
    "shifted log"
  }
  draw <- function(x_col, y_col, labels, x_lab, y_lab, title, legend, pal) {
    df <- data.frame(
      axis1 = scores[[x_col]],
      axis2 = scores[[y_col]],
      label = labels
    )
    # Largest group first, so the rarer groups stay visible (#55, #61).
    df <- df[order(-table(df$label)[df$label]), , drop = FALSE]
    plot <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = axis1, y = axis2, colour = label)
    ) +
      ggplot2::geom_point(alpha = 0.35, size = 0.5) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        x = x_lab,
        y = y_lab,
        colour = legend,
        title = title,
        subtitle = input_label
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))
      ) +
      theme_slide()
    if (!is.null(pal)) {
      plot <- plot + ggplot2::scale_colour_manual(values = pal)
    }
    plot
  }
  out <- list(
    pca = draw("pc1", "pc2", cell_type, pca_x, pca_y, "PCA", "Cell type", colours),
    plsda = draw(
      "pls1", "pls2", cell_type, pls_x, pls_y, "PLS-DA", "Cell type", colours
    )
  )
  for (g in names(groups)) {
    pal <- group_colours[[g]]
    out[[paste0("pca_", g)]] <- draw(
      "pc1", "pc2", groups[[g]], pca_x, pca_y, "PCA", g, pal
    )
    out[[paste0("plsda_", g)]] <- draw(
      "pls1", "pls2", groups[[g]], pls_x, pls_y, "PLS-DA (fitted on types)", g, pal
    )
    scores[[g]] <- groups[[g]]
  }
  out$n_genes <- ncol(x)
  out$pca_explained <- pct[1:2]
  out$plsda_explained <- 100 * as.numeric(explained[1:2])
  out$scores <- scores
  out
}

#' Low-quality-cell diagnostics on the cells kept by the authors.
#'
#' Panel a: detected genes against library size (log10), one facet per
#' time, coloured by cell type (largest type drawn first), with the
#' Suppinger cut-off (1,800 detected genes) and the deMixSC cut-offs
#' (200 genes, 500 counts) as dashed lines labelled in empty space;
#' cells beyond 3 scaled MADs of the log library size of their type are
#' circled. Panel b: distribution of the log10 library size per type and
#' time, types ranked by median, with the median, the +-3 MAD bounds and
#' the number of cells.
#'
#' @param cell_values Per-cell table of the tier (i) RDS: `time_point`,
#'   `cell_type`, `library_size`, `detected_genes`, `mad_flag`.
#' @param colours Named colours per cell type.
#' @return A patchwork object.
plot_low_quality_cells <- function(cell_values, colours = NULL) {
  d <- cell_values
  d$log_lib <- log10(d$library_size)
  d <- d[order(-table(d$cell_type)[d$cell_type]), , drop = FALSE]
  flagged <- d[d$mad_flag, , drop = FALSE]
  n_below <- sum(d$detected_genes < 1800)
  cut_labels <- data.frame(
    x = c(Inf, Inf, 520),
    y = c(1800, 200, Inf),
    label = c(
      "1,800 genes (Suppinger)",
      "200 genes (deMixSC)",
      "500 counts (deMixSC)"
    ),
    hjust = c(1.02, 1.02, 0),
    vjust = c(1.5, -0.5, 1.5)
  )
  a <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = library_size, y = detected_genes, colour = cell_type)
  ) +
    ggplot2::geom_point(size = 0.25, alpha = 0.35) +
    ggplot2::geom_point(
      data = flagged,
      shape = 21,
      size = 1.2,
      colour = "#1c2b24",
      fill = NA,
      stroke = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = c(200, 1800),
      linetype = "dashed",
      colour = "#1c2b24",
      linewidth = 0.3
    ) +
    ggplot2::geom_vline(
      xintercept = 500,
      linetype = "dashed",
      colour = "#1c2b24",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = cut_labels,
      ggplot2::aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE,
      size = 2.6,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_x_log10(labels = scales::label_comma()) +
    ggplot2::scale_y_log10(labels = scales::label_comma()) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      x = "Library size (counts, log scale)",
      y = "Detected genes (log scale)",
      colour = "Cell type",
      title = paste0(
        "a. ",
        fmt_n(nrow(flagged)),
        " kept cells lie beyond 3 MAD of their type's library size (circled)"
      ),
      subtitle = paste0(
        fmt_n(n_below),
        " kept cells count fewer than 1,800 genes on the deposited matrix (genes removed after the authors' QC)"
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 3)
    ) +
    theme_slide()
  bounds <- stats::aggregate(
    log_lib ~ time_point + cell_type,
    data = d,
    FUN = function(v) c(med = stats::median(v), mad = stats::mad(v), n = length(v))
  )
  bounds <- data.frame(
    time_point = bounds$time_point,
    cell_type = bounds$cell_type,
    med = bounds$log_lib[, "med"],
    lo = bounds$log_lib[, "med"] - 3 * bounds$log_lib[, "mad"],
    hi = bounds$log_lib[, "med"] + 3 * bounds$log_lib[, "mad"],
    n = bounds$log_lib[, "n"]
  )
  # Rank the types by median within each time (#50): one factor level per
  # time and type, ordered by median; the time prefix is dropped on the
  # axis.
  key <- function(tm, ct) paste(tm, ct, sep = "___")
  bounds <- bounds[order(bounds$time_point, bounds$med), , drop = FALSE]
  bounds$row <- factor(key(bounds$time_point, bounds$cell_type), levels = key(bounds$time_point, bounds$cell_type))
  d$row <- factor(key(d$time_point, d$cell_type), levels = levels(bounds$row))
  row_labels <- stats::setNames(
    paste0(bounds$cell_type, " (n = ", fmt_n(bounds$n), ")"),
    as.character(bounds$row)
  )
  med_range <- range(10^bounds$med)
  b <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = log_lib, y = row, fill = cell_type)
  )
  b <- if (requireNamespace("ggdist", quietly = TRUE)) {
    b +
      ggdist::stat_halfeye(
        adjust = 0.8,
        width = 0.7,
        .width = 0,
        point_colour = NA,
        alpha = 0.8
      )
  } else {
    b + ggplot2::geom_violin(alpha = 0.8, colour = NA)
  }
  b <- b +
    ggplot2::geom_errorbar(
      data = bounds,
      ggplot2::aes(xmin = lo, xmax = hi, y = row),
      inherit.aes = FALSE,
      orientation = "y",
      width = 0.25,
      colour = "#1c2b24",
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      data = bounds,
      ggplot2::aes(x = med, y = row),
      inherit.aes = FALSE,
      size = 1,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_y_discrete(labels = row_labels) +
    ggplot2::facet_wrap(~time_point, nrow = 1, scales = "free_y") +
    ggplot2::labs(
      x = "log10 library size",
      y = NULL,
      title = paste0(
        "b. Median library sizes of the types span only ",
        fmt_n(med_range[1L]),
        " to ",
        fmt_n(med_range[2L]),
        " counts"
      ),
      subtitle = "Types ranked by median, n = cells; point: median; bar: median +- 3 MAD (log10 scale)"
    ) +
    theme_slide() +
    ggplot2::theme(legend.position = "none")
  if (!is.null(colours)) {
    a <- a + ggplot2::scale_colour_manual(values = colours)
    b <- b + ggplot2::scale_fill_manual(values = colours)
  }
  patchwork::wrap_plots(a, b, ncol = 1, heights = c(1.1, 1))
}

#' Gene-filter diagnostics of tier (i).
#'
#' Panel a: genes failing each recorded rule, per time. Panel b: share
#' of fractional non-zero entries per gene (genes detected in at least
#' `min_nonzero` cells), with the 0.9 cut. Panel c: genes entering tier
#' (ii) against the authors' 2,944 SCT genes, per time.
#'
#' @param gene_summary `gene_filter_summary` table (time, rule, genes).
#' @param fractional `fractional` table of the tier (i) RDS.
#' @param gene_rules `gene_rules` table of the tier (i) RDS.
#' @param cut Fractional-share cut. Default `0.9`.
#' @param min_nonzero Smallest number of non-zero cells for panel b.
#'   Default `50`.
#' @return A patchwork object.
plot_gene_filters <- function(
  gene_summary,
  fractional,
  gene_rules,
  cut = 0.9,
  min_nonzero = 50L
) {
  rules <- c(
    "zero in every cell of the slot" = "Zero in the slot",
    "detected in fewer than 5 cells of a batch (sctransform min_cells)" =
      "< 5 cells in a batch (SCTransform)",
    "detected in < 5 % of the cells of the slot (deMixSC)" =
      "< 5 % of cells (deMixSC)",
    "detected in < 20 % of the cells of every type (RNA-Sieve)" =
      "< 20 % in every type (RNA-Sieve)",
    "fractional share above 0.9" = "Fractional share > 0.9",
    "not shared one-to-one with the bulk" = "Not one-to-one with bulk",
    "fails the presence gate (>= 5 % in every scored type, variance > 0)" =
      "Fails presence gate (this book)"
  )
  gs <- gene_summary[gene_summary$rule %in% names(rules), , drop = FALSE]
  gs$label <- factor(rules[gs$rule], levels = rev(unname(rules)))
  gate_n <- gs$genes[gs$label == "Fails presence gate (this book)"]
  a <- ggplot2::ggplot(gs, ggplot2::aes(x = genes, y = label)) +
    ggplot2::geom_col(fill = "#1d4e89", width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(genes)),
      hjust = -0.1,
      size = 2.6,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 10000, 20000),
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.35))
    ) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      x = "Genes failing the rule (of 27,027)",
      y = NULL,
      title = paste0(
        "a. Most genes are removed for low detection: ",
        fmt_n(min(gate_n)),
        " to ",
        fmt_n(max(gate_n)),
        " fail the presence gate"
      ),
      subtitle = "Rules in the order of the pipeline, top to bottom; a gene can fail several rules"
    ) +
    theme_slide()
  fr <- fractional[
    !is.na(fractional$fractional_share) &
      fractional$n_nonzero >= min_nonzero,
    ,
    drop = FALSE
  ]
  b <- ggplot2::ggplot(fr, ggplot2::aes(x = fractional_share)) +
    ggplot2::geom_histogram(
      binwidth = 0.02,
      boundary = 0,
      fill = "#d9d3c7",
      colour = "white"
    ) +
    ggplot2::geom_vline(
      xintercept = cut,
      linetype = "dashed",
      colour = "#c45c26"
    ) +
    ggplot2::annotate(
      "text",
      x = cut,
      y = Inf,
      vjust = 1.5,
      hjust = 1.05,
      size = 2.6,
      colour = "#c45c26",
      label = paste0(
        scales::comma(sum(fr$fractional_share > cut)),
        " genes dropped"
      )
    ) +
    ggplot2::scale_y_continuous(
      trans = "log1p",
      breaks = c(0, 10, 100, 1000, 10000),
      labels = scales::label_comma()
    ) +
    ggplot2::labs(
      x = "Share of fractional non-zero entries",
      y = "Genes (log1p scale)",
      title = paste0(
        "b. ",
        scales::comma(sum(fr$fractional_share > cut)),
        " well-detected genes are mostly fractional counts"
      ),
      subtitle = paste0(
        scales::comma(nrow(fr)),
        " genes detected in >= ",
        min_nonzero,
        " cells; bins of 0.02. Panel a counts the rule on all genes"
      )
    ) +
    theme_slide()
  ov <- do.call(
    rbind,
    lapply(split(gene_rules, gene_rules$time_point), function(r) {
      data.frame(
        time_point = r$time_point[1L],
        set = c("Tier (ii) and SCT", "Tier (ii) only", "SCT only"),
        genes = c(
          sum(r$enters_tier_ii & r$in_sct_layer),
          sum(r$enters_tier_ii & !r$in_sct_layer),
          sum(!r$enters_tier_ii & r$in_sct_layer)
        ),
        stringsAsFactors = FALSE
      )
    })
  )
  ov$set <- factor(
    ov$set,
    levels = c("SCT only", "Tier (ii) only", "Tier (ii) and SCT")
  )
  shared <- ov$genes[ov$set == "Tier (ii) and SCT"]
  c_plot <- ggplot2::ggplot(ov, ggplot2::aes(x = genes, y = set)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = genes, yend = set),
      colour = "grey70"
    ) +
    ggplot2::geom_point(colour = "#1d4e89", size = 2.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(genes)),
      hjust = -0.25,
      size = 2.5,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_comma(),
      breaks = scales::breaks_extended(3),
      expand = ggplot2::expansion(mult = c(0, 0.45))
    ) +
    ggplot2::facet_wrap(~time_point, ncol = 1) +
    ggplot2::labs(
      x = "Genes",
      y = NULL,
      title = paste0(
        "c. ",
        fmt_n(min(shared)),
        " to ",
        fmt_n(max(shared)),
        " of the 2,944\nSCT genes enter tier (ii)"
      )
    ) +
    theme_slide()
  patchwork::wrap_plots(
    a,
    patchwork::wrap_plots(b, c_plot, nrow = 1, widths = c(1.4, 1)),
    ncol = 1,
    heights = c(1, 1)
  )
}
