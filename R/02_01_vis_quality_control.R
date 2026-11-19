#' Tier (i) figures: projections, cells, gene filters, fractional counts,
#' markers and expression before/after the gate.
#'
#' Plots only; nothing is saved here (see
#' `scripts/02_01_vis_quality_control.R`). Requires `R/02_00_utils.R`
#' (`theme_slide()`, `fmt_n()`).

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

#' Colour scale and grouped legend of the cell types.
#'
#' Types keep the order of `suppinger_celltype_colours.csv` and are
#' numbered in that order (the numbers drawn on the projections); the
#' legend groups them by lineage (`legendry::guide_legend_group()`).
#'
#' @param colours Named colours per type ([celltype_colours()]).
#' @param lineages Named lineage per type ([celltype_lineages()]).
#' @param present Types present in the plot.
#' @param aesthetics `"colour"` or `"fill"`.
#' @param numbered Prefix each label with the type number. Default `TRUE`.
#' @return A list of a scale and a guide, to add to a ggplot.
celltype_scale <- function(
  colours,
  lineages,
  present,
  aesthetics = "colour",
  numbered = TRUE
) {
  breaks <- names(colours)[names(colours) %in% present]
  number <- stats::setNames(seq_along(colours), names(colours))
  labels <- if (numbered) paste0(number[breaks], "  ", breaks) else breaks
  key <- legendry::key_group_lut(
    members = names(lineages),
    group = as.character(lineages)
  )
  guide <- legendry::guide_legend_group(
    key = key,
    override.aes = list(size = 2.5, alpha = 1)
  )
  list(
    ggplot2::scale_discrete_manual(
      aesthetics = aesthetics,
      values = colours,
      breaks = breaks,
      labels = labels
    ),
    do.call(ggplot2::guides, stats::setNames(list(guide), aesthetics))
  )
}

#' PCA, UMAP and sPLS-DA scores of one slot.
#'
#' PCA on the centred, unscaled matrix (`prcomp()`, first `n_pcs`
#' components); UMAP (`uwot::umap()`) on those components;
#' sPLS-DA (`mixOmics::splsda()`) with the cell type as response and
#' `keep_x` genes per component.
#'
#' @param x Cells by genes matrix, already on its analysis scale.
#' @param cell_type Cell-type label per row.
#' @param n_pcs Components passed to UMAP. Default `30`.
#' @param keep_x Genes kept per sPLS-DA component. Default `200`.
#' @param n_neighbors,min_dist UMAP settings. Defaults `30`, `0.3`.
#' @param seed Seed of the UMAP initialisation.
#' @return A list with `scores` (data frame: `pc1`, `pc2`, `umap1`,
#'   `umap2`, `pls1`, `pls2`, `cell_type`), `pca_explained` and
#'   `plsda_explained` (percent, two values each), `n_genes` and the
#'   settings.
project_slot <- function(
  x,
  cell_type,
  n_pcs = 30L,
  keep_x = 200L,
  n_neighbors = 30L,
  min_dist = 0.3,
  seed = 20261119L
) {
  for (pkg in c("mixOmics", "uwot")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package ", pkg, " is required for the projections.")
    }
  }
  cell_type <- as.character(cell_type)
  if (length(cell_type) != nrow(x) || anyNA(cell_type)) {
    stop("`cell_type` needs one non-missing label per row of `x`.")
  }
  if (length(unique(cell_type)) < 2L) {
    stop("sPLS-DA needs at least two cell types.")
  }
  gene_sd <- apply(x, 2L, stats::sd)
  x <- x[, is.finite(gene_sd) & gene_sd > 0, drop = FALSE]
  n_pcs <- min(n_pcs, nrow(x) - 1L, ncol(x))
  pc <- stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = n_pcs)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  um <- withr::with_seed(
    seed,
    uwot::umap(
      pc$x,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      n_threads = 1L
    )
  )
  keep_x <- min(keep_x, ncol(x))
  fit <- mixOmics::splsda(
    X = x,
    Y = factor(cell_type),
    ncomp = 2L,
    keepX = c(keep_x, keep_x),
    scale = FALSE
  )
  pls <- as.matrix(fit$variates$X)
  list(
    scores = data.frame(
      pc1 = pc$x[, 1L],
      pc2 = pc$x[, 2L],
      umap1 = um[, 1L],
      umap2 = um[, 2L],
      pls1 = pls[, 1L],
      pls2 = pls[, 2L],
      cell_type = cell_type,
      stringsAsFactors = FALSE
    ),
    pca_explained = pct[1:2],
    plsda_explained = 100 * as.numeric(fit$prop_expl_var$X[1:2]),
    n_genes = ncol(x),
    settings = list(
      n_pcs = n_pcs,
      keep_x = keep_x,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      seed = seed
    )
  )
}

#' One projection scatter.
#'
#' Largest group drawn first so rare groups stay visible; for the cell
#' types, each type's number is drawn at its median (as in Suppinger et
#' al. 2023, Fig. 1B), and PCA and sPLS-DA panels get dashed zero lines.
#'
#' @param scores `scores` of [project_slot()], with a `label` column (the
#'   colour grouping).
#' @param method `"pca"`, `"umap"` or `"plsda"`.
#' @param explained Percent of variance of the two axes (`NULL` for UMAP).
#' @param scale List from [celltype_scale()], or a manual colour scale.
#' @param number_labels Named type numbers to draw at the medians, or
#'   `NULL`.
#' @param title,subtitle Panel title and subtitle.
#' @return A ggplot.
plot_projection <- function(
  scores,
  method = c("pca", "umap", "plsda"),
  explained = NULL,
  scale = NULL,
  number_labels = NULL,
  title = NULL,
  subtitle = NULL
) {
  method <- match.arg(method)
  cols <- switch(
    method,
    pca = c("pc1", "pc2"),
    umap = c("umap1", "umap2"),
    plsda = c("pls1", "pls2")
  )
  df <- data.frame(
    axis1 = scores[[cols[1L]]],
    axis2 = scores[[cols[2L]]],
    label = scores$label
  )
  df <- df[order(-table(df$label)[df$label]), , drop = FALSE]
  axis_title <- switch(
    method,
    pca = sprintf("PC%d (%.1f %%)", 1:2, explained),
    umap = sprintf("UMAP %d (30 PCs, 30 neighbours)", 1:2),
    plsda = sprintf("sPLS-DA comp. %d (%.1f %%)", 1:2, explained)
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(axis1, axis2, colour = label))
  if (method != "umap") {
    p <- p +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "#1c2b24",
        linewidth = 0.3
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        colour = "#1c2b24",
        linewidth = 0.3
      )
  }
  p <- p +
    ggplot2::geom_point(alpha = 0.4, size = 0.4) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = axis_title[1L],
      y = axis_title[2L],
      colour = NULL,
      title = title,
      subtitle = subtitle
    ) +
    scale +
    theme_slide()
  if (!is.null(number_labels)) {
    med <- stats::aggregate(cbind(axis1, axis2) ~ label, data = df, FUN = stats::median)
    med$number <- number_labels[med$label]
    p <- p +
      ggplot2::geom_label(
        data = med,
        ggplot2::aes(axis1, axis2, label = number),
        inherit.aes = FALSE,
        size = 2.6,
        label.padding = ggplot2::unit(0.12, "lines"),
        label.size = 0.2,
        colour = "#1c2b24",
        fill = "white",
        alpha = 0.85
      )
  }
  p
}

#' Projection matrix of one slot on the analysis scale.
#'
#' @param counts Genes by cells raw counts, or, for `input = "sct"`, the
#'   authors' SCT residuals (genes by cells).
#' @param size_factors Size factor per cell (not used for `"sct"`).
#' @param input `"pearson"` (Pearson residuals of M0,
#'   [pearson_residuals_m0()]), `"shifted_log"` (`log1p(x / s)`) or
#'   `"sct"` (residuals passed as they are).
#' @return Cells by genes matrix.
projection_input <- function(
  counts,
  size_factors = NULL,
  input = c("pearson", "shifted_log", "sct")
) {
  input <- match.arg(input)
  if (input != "sct") {
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
  }
  switch(
    input,
    pearson = t(pearson_residuals_m0(counts, size_factors)),
    shifted_log = t(log1p(sweep(as.matrix(counts), 2L, size_factors, "/"))),
    sct = t(as.matrix(counts))
  )
}

#' PCA, UMAP and sPLS-DA of one time slot.
#'
#' Descriptive projections of the cells that enter tier (ii); not a
#' selection rule. The matrix is built by [projection_input()], centred
#' and not scaled (scaling would undo the variance stabilisation of the
#' residuals), then projected once by [project_slot()]. `groups` adds the
#' same scatters coloured by another labelling (for example the cell
#' line).
#'
#' @param counts,size_factors,input See [projection_input()].
#' @param cell_type Cell-type label per column.
#' @param colours,lineages Named colours and lineages per type
#'   ([celltype_colours()], [celltype_lineages()]).
#' @param groups Optional named list of extra labellings (one value per
#'   cell), each drawn as `pca_<name>`, `umap_<name>` and `plsda_<name>`.
#' @param group_colours Optional named list of named colour vectors.
#' @param ... Passed to [project_slot()].
#' @return A list with `pca`, `umap`, `plsda` (ggplots coloured by type),
#'   the extra panels of `groups`, and the elements of [project_slot()]
#'   (`scores` gains one column per entry of `groups`).
plot_type_pca_plsda <- function(
  counts,
  cell_type,
  size_factors = NULL,
  colours,
  lineages,
  input = c("pearson", "shifted_log", "sct"),
  groups = NULL,
  group_colours = NULL,
  ...
) {
  input <- match.arg(input)
  if (length(cell_type) != ncol(counts)) {
    stop("`cell_type` must have one label per column of `counts`.")
  }
  proj <- project_slot(
    projection_input(counts, size_factors, input),
    cell_type,
    ...
  )
  scores <- proj$scores
  number <- stats::setNames(seq_along(colours), names(colours))
  type_scale <- celltype_scale(colours, lineages, unique(scores$cell_type))
  input_label <- c(
    pearson = "Pearson residuals of M0",
    shifted_log = "shifted log",
    sct = "authors' SCT residuals"
  )[[input]]
  draw <- function(method, labels, scale, numbers, title) {
    s <- scores
    s$label <- labels
    explained <- switch(
      method,
      pca = proj$pca_explained,
      plsda = proj$plsda_explained,
      umap = NULL
    )
    plot_projection(s, method, explained, scale, numbers, title, input_label)
  }
  out <- list()
  titles <- c(pca = "PCA", umap = "UMAP", plsda = "sPLS-DA")
  for (m in names(titles)) {
    out[[m]] <- draw(m, scores$cell_type, type_scale, number, titles[[m]])
  }
  for (g in names(groups)) {
    pal <- ggplot2::scale_colour_manual(
      values = group_colours[[g]],
      guide = ggplot2::guide_legend(override.aes = list(size = 2.5, alpha = 1))
    )
    for (m in names(titles)) {
      out[[paste0(m, "_", g)]] <- draw(m, groups[[g]], pal, NULL, titles[[m]])
    }
    scores[[g]] <- groups[[g]]
  }
  proj$scores <- scores
  c(out, proj)
}

#' UMAP of several slots, one panel per highlighted slot.
#'
#' As Suppinger et al. (2023), Fig. S2A: one embedding of all the cells,
#' redrawn once per time with the cells of that time coloured by type
#' and the others in grey.
#'
#' @param scores Data frame with `umap1`, `umap2`, `cell_type`,
#'   `time_point`.
#' @param colours,lineages See [celltype_scale()].
#' @return A ggplot.
plot_pooled_umap <- function(scores, colours, lineages) {
  times <- unique(scores$time_point)
  bg <- do.call(
    rbind,
    lapply(times, function(tm) transform(scores[, c("umap1", "umap2")], panel = tm))
  )
  fg <- transform(scores, panel = time_point)
  fg <- fg[order(-table(fg$cell_type)[fg$cell_type]), , drop = FALSE]
  number <- stats::setNames(seq_along(colours), names(colours))
  med <- stats::aggregate(cbind(umap1, umap2) ~ cell_type + panel, data = fg, FUN = stats::median)
  med$number <- number[med$cell_type]
  ggplot2::ggplot(fg, ggplot2::aes(umap1, umap2)) +
    ggplot2::geom_point(data = bg, colour = "grey85", size = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = cell_type), size = 0.4, alpha = 0.5) +
    ggplot2::geom_label(
      data = med,
      ggplot2::aes(label = number),
      size = 2.6,
      label.padding = ggplot2::unit(0.12, "lines"),
      label.size = 0.2,
      colour = "#1c2b24",
      fill = "white",
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::coord_equal() +
    celltype_scale(colours, lineages, unique(scores$cell_type)) +
    ggplot2::labs(
      x = "UMAP 1 (30 PCs, 30 neighbours)",
      y = "UMAP 2",
      colour = NULL
    ) +
    theme_slide() +
    ggplot2::theme(legend.position = "right")
}

#' Cells kept by the authors at the bulk-matched times.
#'
#' Panel a: detected genes against library size (log10), one facet per
#' time, coloured by cell type (largest type drawn first), with the
#' Suppinger cut-off (1,800 detected genes). Panel b: distribution of the
#' log10 library size per type and time, types ranked by median, with
#' the median and the number of cells.
#'
#' @param cell_values Per-cell table of the tier (i) RDS: `time_point`,
#'   `cell_type`, `library_size`, `detected_genes`.
#' @param colours,lineages See [celltype_scale()].
#' @return A patchwork object.
plot_low_quality_cells <- function(cell_values, colours, lineages) {
  d <- cell_values
  d$log_lib <- log10(d$library_size)
  d <- d[order(-table(d$cell_type)[d$cell_type]), , drop = FALSE]
  n_below <- sum(d$detected_genes <= 1800)
  a <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = library_size, y = detected_genes, colour = cell_type)
  ) +
    ggplot2::geom_point(size = 0.25, alpha = 0.35) +
    ggplot2::geom_hline(
      yintercept = 1800,
      linetype = "dashed",
      colour = "#1c2b24",
      linewidth = 0.3
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = 1800,
      hjust = 1.02,
      vjust = 1.5,
      size = 2.6,
      colour = "#1c2b24",
      label = "1,800 genes (Suppinger et al.)"
    ) +
    ggplot2::scale_x_log10(labels = scales::label_comma()) +
    ggplot2::scale_y_log10(labels = scales::label_comma()) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      x = "Library size (counts, log scale)",
      y = "Detected genes (log scale)",
      title = paste0(
        "a. ",
        fmt_n(n_below),
        " kept cells count at most 1,800 genes on the deposited matrix"
      ),
      subtitle = "Gm, mitochondrial and ribosomal genes were removed after the authors' cut-off"
    ) +
    celltype_scale(colours, lineages, unique(d$cell_type), numbered = FALSE) +
    theme_slide() +
    ggplot2::theme(legend.position = "right")
  med <- stats::aggregate(
    log_lib ~ time_point + cell_type,
    data = d,
    FUN = function(v) c(med = stats::median(v), n = length(v))
  )
  med <- data.frame(
    time_point = med$time_point,
    cell_type = med$cell_type,
    med = med$log_lib[, "med"],
    n = med$log_lib[, "n"]
  )
  # One factor level per time and type, ordered by median (#50).
  key <- function(tm, ct) paste(tm, ct, sep = "___")
  med <- med[order(med$time_point, med$med), , drop = FALSE]
  med$row <- factor(key(med$time_point, med$cell_type), levels = key(med$time_point, med$cell_type))
  d$row <- factor(key(d$time_point, d$cell_type), levels = levels(med$row))
  row_labels <- stats::setNames(
    paste0(med$cell_type, " (n = ", fmt_n(med$n), ")"),
    as.character(med$row)
  )
  med_range <- range(10^med$med)
  b <- ggplot2::ggplot(d, ggplot2::aes(x = log_lib, y = row, fill = cell_type)) +
    ggdist::stat_halfeye(
      adjust = 0.8,
      width = 0.7,
      .width = 0,
      point_colour = NA,
      alpha = 0.85
    ) +
    ggplot2::geom_point(
      data = med,
      ggplot2::aes(x = med, y = row),
      inherit.aes = FALSE,
      size = 1,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_y_discrete(labels = row_labels) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::facet_wrap(~time_point, nrow = 1, scales = "free_y") +
    ggplot2::labs(
      x = "log10 library size",
      y = NULL,
      title = paste0(
        "b. Median library sizes of the types span ",
        fmt_n(med_range[1L]),
        " to ",
        fmt_n(med_range[2L]),
        " counts"
      ),
      subtitle = "Types ranked by median (point); n = cells"
    ) +
    theme_slide()
  patchwork::wrap_plots(a, b, ncol = 1, heights = c(1.1, 1))
}

#' Gene filters of tier (i).
#'
#' Panel a: genes removed by each applied rule, per time, bars ranked by
#' decreasing count. Panel b: genes entering tier (ii) against the
#' authors' SCT genes, per time.
#'
#' @param gene_summary `gene_filter_summary` table (time, rule, genes).
#' @param gene_rules `gene_rules` table of the tier (i) RDS.
#' @return A patchwork object.
plot_gene_filters <- function(gene_summary, gene_rules) {
  rules <- c(
    "zero in every cell of the slot" = "Zero in the slot",
    "not shared one-to-one with the bulk" = "Not one-to-one with the bulk",
    "fails the presence gate (>= 20 cells, >= 2 in every type, variance > 0)" =
      "Fails the presence gate"
  )
  gs <- gene_summary[gene_summary$rule %in% names(rules), , drop = FALSE]
  gs$label <- unname(rules[gs$rule])
  totals <- tapply(gs$genes, gs$label, sum)
  gs$label <- factor(gs$label, levels = names(sort(totals)))
  gate_n <- gs$genes[gs$label == "Fails the presence gate"]
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
      x = paste0("Genes (of ", fmt_n(gene_summary$genes[gene_summary$rule == "deposited genes"][1L]), ")"),
      y = NULL,
      title = paste0(
        "a. Low detection removes most genes: ",
        fmt_n(min(gate_n)),
        " to ",
        fmt_n(max(gate_n)),
        " fail the gate"
      ),
      subtitle = "A gene can fail several rules"
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
  ov$set <- factor(ov$set, levels = names(sort(tapply(ov$genes, ov$set, sum))))
  shared <- ov$genes[ov$set == "Tier (ii) and SCT"]
  b <- ggplot2::ggplot(ov, ggplot2::aes(x = genes, y = set)) +
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
      expand = ggplot2::expansion(mult = c(0, 0.3))
    ) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(
      x = "Genes",
      y = NULL,
      title = paste0(
        "b. ",
        fmt_n(min(shared)),
        " to ",
        fmt_n(max(shared)),
        " of the 2,944 SCT genes enter tier (ii)"
      )
    ) +
    theme_slide()
  patchwork::wrap_plots(a, b, ncol = 1, heights = c(1, 0.8))
}

#' Data behind the fractional-count figure.
#'
#' @param fractional_group `fractional_group` of the tier (i) RDS (one
#'   row per gene, type and time).
#' @param biotype `biotype` of the tier (i) RDS.
#' @param min_nonzero Genes with fewer non-zero cells in a type are left
#'   out of that type. Default `10`.
#' @param mostly Share above which a gene counts as mostly fractional.
#'   Default `0.5`.
#' @return A list with `shares` (gene, type, time, share, biotype) and
#'   `biotypes` (per type, time, biotype and class: genes and percent of
#'   the class).
fractional_figure_data <- function(
  fractional_group,
  biotype,
  min_nonzero = 10L,
  mostly = 0.5
) {
  grouping <- c(
    protein_coding = "Protein coding",
    pseudogene = "Pseudogene",
    lncRNA = "Non-coding RNA",
    other_non_coding_RNA = "Non-coding RNA",
    miRNA = "Small RNA",
    snRNA = "Small RNA",
    snoRNA = "Small RNA",
    rRNA = "Small RNA",
    tRNA = "Small RNA",
    other = "Other or unknown",
    not_annotated = "Not in org.Mm.eg.db"
  )
  s <- fractional_group[fractional_group$n_nonzero >= min_nonzero, , drop = FALSE]
  s$biotype <- unname(grouping[biotype$gene_biotype[match(s$gene, biotype$gene)]])
  s$biotype[is.na(s$biotype)] <- "Not in org.Mm.eg.db"
  s$class <- ifelse(
    s$n_fractional == 0,
    "Integer counts only",
    ifelse(s$fractional_share > mostly, "Mostly fractional", NA_character_)
  )
  cls <- s[!is.na(s$class), , drop = FALSE]
  tab <- stats::aggregate(
    gene ~ time_point + cell_type + class + biotype,
    data = cls,
    FUN = length
  )
  names(tab)[names(tab) == "gene"] <- "genes"
  tot <- stats::ave(tab$genes, tab$time_point, tab$cell_type, tab$class, FUN = sum)
  tab$percent <- 100 * tab$genes / tot
  list(
    shares = s[, c("gene", "cell_type", "time_point", "n_nonzero", "fractional_share", "biotype")],
    biotypes = tab,
    settings = list(min_nonzero = min_nonzero, mostly = mostly)
  )
}

#' Fractional counts by cell type and biotype.
#'
#' Panel a: ridges of the per-gene share of fractional non-zero entries,
#' one ridge per type, one facet per time. Panel b: per type and time,
#' the percent of each biotype among mostly fractional genes and among
#' genes with integer counts only (dumbbell).
#'
#' @param data Output of [fractional_figure_data()].
#' @param colours,lineages See [celltype_scale()].
#' @return A patchwork object.
plot_fractional_genes <- function(data, colours, lineages) {
  s <- data$shares
  types <- names(colours)[names(colours) %in% s$cell_type]
  s$cell_type <- factor(s$cell_type, levels = rev(types))
  med <- stats::aggregate(fractional_share ~ time_point + cell_type, data = s, FUN = stats::median)
  a <- ggplot2::ggplot(
    s,
    ggplot2::aes(x = fractional_share, y = cell_type, fill = cell_type)
  ) +
    ggridges::geom_density_ridges(
      scale = 1.1,
      rel_min_height = 0.005,
      colour = "white",
      linewidth = 0.3,
      bandwidth = 0.03
    ) +
    ggplot2::geom_point(
      data = med,
      ggplot2::aes(x = fractional_share, y = cell_type),
      inherit.aes = FALSE,
      size = 1,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(),
      breaks = c(0, 0.5, 1),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::facet_wrap(~time_point, nrow = 1, scales = "free_y") +
    ggplot2::labs(
      x = "Fractional entries among the gene's non-zero cells",
      y = NULL,
      title = "a. Most genes carry few fractional entries; a tail is mostly fractional",
      subtitle = sprintf(
        "Genes with >= %d non-zero cells in the type; point: median",
        data$settings$min_nonzero
      )
    ) +
    theme_slide()
  b_data <- data$biotypes
  wide <- stats::reshape(
    b_data[, c("time_point", "cell_type", "biotype", "class", "percent")],
    idvar = c("time_point", "cell_type", "biotype"),
    timevar = "class",
    direction = "wide"
  )
  names(wide) <- sub("^percent\\.", "", names(wide))
  wide[is.na(wide)] <- 0
  enrich <- tapply(
    wide[["Mostly fractional"]] - wide[["Integer counts only"]],
    wide$biotype,
    mean
  )
  b_data$biotype <- factor(b_data$biotype, levels = names(sort(enrich)))
  wide$biotype <- factor(wide$biotype, levels = levels(b_data$biotype))
  b_data$panel <- paste(b_data$time_point, b_data$cell_type, sep = " · ")
  wide$panel <- paste(wide$time_point, wide$cell_type, sep = " · ")
  order_panels <- unique(b_data$panel[order(
    b_data$time_point,
    match(b_data$cell_type, names(colours))
  )])
  b_data$panel <- factor(b_data$panel, levels = order_panels)
  wide$panel <- factor(wide$panel, levels = order_panels)
  class_colours <- c("Mostly fractional" = "#c45c26", "Integer counts only" = "#1d4e89")
  b <- ggplot2::ggplot(b_data, ggplot2::aes(x = percent, y = biotype)) +
    ggplot2::geom_segment(
      data = wide,
      ggplot2::aes(
        x = `Integer counts only`,
        xend = `Mostly fractional`,
        y = biotype,
        yend = biotype
      ),
      colour = "grey70"
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = class), size = 1.6) +
    ggplot2::scale_colour_manual(values = class_colours) +
    ggplot2::scale_x_sqrt(
      breaks = c(1, 5, 25, 50, 100),
      labels = function(x) paste0(x, " %")
    ) +
    ggplot2::facet_wrap(~panel, ncol = 4) +
    ggplot2::labs(
      x = "Percent of the genes of the class",
      y = NULL,
      colour = NULL,
      title = "b. Biotypes of mostly fractional and of integer-only genes",
      subtitle = sprintf(
        "Mostly fractional: > %.0f %% of non-zero entries; square-root scale; biotypes ranked by mean difference",
        100 * data$settings$mostly
      )
    ) +
    theme_slide()
  patchwork::wrap_plots(a, b, ncol = 1, heights = c(0.8, 1.6))
}

#' Curated markers against the authors' SCT genes.
#'
#' Two-set Venn (`ggVennDiagram`), with the markers outside the SCT
#' genes named and the hypergeometric tail of the overlap under a
#' uniform draw of the markers from all deposited symbols.
#'
#' @param markers Mouse symbols of the curated markers.
#' @param sct_genes Genes of the authors' SCT layer.
#' @param universe All deposited symbols.
#' @return A list with `plot` and `data` (sets, counts, p-value).
plot_marker_venn <- function(markers, sct_genes, universe) {
  markers <- unique(stats::na.omit(markers))
  in_both <- intersect(markers, sct_genes)
  absent <- sort(setdiff(markers, sct_genes))
  p_tail <- stats::phyper(
    length(in_both) - 1L,
    length(sct_genes),
    length(universe) - length(sct_genes),
    length(markers),
    lower.tail = FALSE
  )
  sets <- list(markers, sct_genes)
  names(sets) <- c(
    paste0("Markers (", length(markers), ")"),
    paste0("SCT genes (", fmt_n(length(sct_genes)), ")")
  )
  plot <- ggVennDiagram::ggVennDiagram(
    sets,
    label = "count",
    label_alpha = 0,
    set_size = 3.5,
    edge_size = 0.4
  ) +
    ggplot2::scale_fill_gradient(low = "#f4f1ea", high = "#9ab6d3", guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.4)) +
    ggplot2::labs(
      title = sprintf(
        "%d of %d markers are SCT genes (expected %.1f; P(X >= %d) = %s)",
        length(in_both),
        length(markers),
        length(markers) * length(sct_genes) / length(universe),
        length(in_both),
        format(signif(p_tail, 2))
      ),
      subtitle = paste0("Markers outside the SCT genes: ", paste(absent, collapse = ", "))
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", colour = "#0e4d3a", size = 11),
      plot.subtitle = ggplot2::element_text(colour = "#3d4a42", size = 10),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  list(
    plot = plot,
    data = list(
      markers = markers,
      sct_genes = sct_genes,
      n_universe = length(universe),
      in_both = in_both,
      absent = absent,
      p_tail = p_tail
    )
  )
}

#' Per-gene mean expression by type, before and after the gene gate.
#'
#' Mean over the cells of a type of `log1p(x_gc / s_c)`, `s_c` the
#' library size over its median.
#'
#' @param slot One slot of the tier (i) RDS (`counts`, `cell_type`,
#'   `library_size`, `time_point`).
#' @param kept Genes that pass the gate.
#' @return Data frame with `gene`, `cell_type`, `time_point`,
#'   `mean_log` and `gate` (`"Before the gate"` for every gene,
#'   `"After the gate"` for the kept genes, repeated).
gene_mean_by_type <- function(slot, kept) {
  sf <- slot$library_size / stats::median(slot$library_size)
  x <- methods::as(slot$counts %*% Matrix::Diagonal(x = 1 / sf), "CsparseMatrix")
  x@x <- log1p(x@x)
  types <- unique(slot$cell_type)
  out <- do.call(
    rbind,
    lapply(types, function(ty) {
      data.frame(
        gene = rownames(x),
        cell_type = ty,
        mean_log = Matrix::rowMeans(x[, slot$cell_type == ty, drop = FALSE]),
        stringsAsFactors = FALSE
      )
    })
  )
  out$time_point <- slot$time_point
  rbind(
    transform(out, gate = "Before the gate"),
    transform(out[out$gene %in% kept, , drop = FALSE], gate = "After the gate")
  )
}

#' Ridges of per-gene mean expression, before and after the gate.
#'
#' @param means Output of [gene_mean_by_type()], several slots bound.
#' @param colours Named colours per type.
#' @return A ggplot.
plot_density_before_after <- function(means, colours) {
  types <- names(colours)[names(colours) %in% means$cell_type]
  means$cell_type <- factor(means$cell_type, levels = rev(types))
  means$gate <- factor(means$gate, levels = c("Before the gate", "After the gate"))
  zero <- stats::aggregate(
    mean_log ~ time_point + gate + cell_type,
    data = means,
    FUN = function(v) mean(v <= 1e-8)
  )
  zero$label <- sprintf("%d %% zero", round(100 * zero$mean_log))
  n_genes <- tapply(means$gene, list(means$gate, means$time_point), function(g) length(unique(g)))
  ggplot2::ggplot(means, ggplot2::aes(x = mean_log, y = cell_type, fill = cell_type)) +
    ggridges::geom_density_ridges(
      scale = 1.1,
      rel_min_height = 0.005,
      colour = "white",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = zero,
      ggplot2::aes(x = Inf, y = cell_type, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = -0.6,
      size = 2.4,
      colour = "#1c2b24"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::facet_grid(gate ~ time_point, scales = "free_y") +
    ggplot2::labs(
      x = "Per-gene mean of log1p(count / size factor) in the type",
      y = NULL,
      title = "The gate removes the genes with a zero mean in some type",
      subtitle = paste0(
        "Genes before / after: ",
        paste(
          sprintf(
            "%s %s / %s",
            colnames(n_genes),
            fmt_n(n_genes["Before the gate", ]),
            fmt_n(n_genes["After the gate", ])
          ),
          collapse = "; "
        ),
        ". Label: share of genes with a zero mean in the type"
      )
    ) +
    theme_slide()
}
