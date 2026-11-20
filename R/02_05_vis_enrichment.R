#' Figures of the enrichment analysis (appendix of chapter 2).
#'
#' Plots only; `scripts/02_05_vis_enrichment.R` saves them. Inputs are
#' tables of `output/feature_selection/enrichment/`.

db_colours <- c("GO BP" = "#1d4e89", "Reactome" = "#c45c26", "Hallmark" = "#5b8c5a")
method_shapes <- c(ORA = 16, GSEA = 17)

#' Wrap long term names.
#'
#' @param x Character vector.
#' @param width Line width.
#' @return Character vector with line breaks.
wrap_terms <- function(x, width = 45) {
  vapply(x, function(s) paste(strwrap(s, width), collapse = "\n"), "")
}

#' Summary lollipop across databases and methods.
#'
#' One row per curated term (`themes`), ordered by significance within a
#' theme; x: -log10 adjusted p; colour: database; shape: ORA or GSEA;
#' the cell type and time are written in the row label.
#'
#' @param sig Significant terms (`enrichment_significant.csv`).
#' @param themes Data frame `theme`, `time_point`, `cell_type`, `method`,
#'   `database`, `id`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_enrichment_lollipop <- function(sig, themes, title = NULL, subtitle = NULL) {
  keys <- c("time_point", "cell_type", "method", "database", "id")
  d <- merge(themes[, c("theme", keys)], sig, by = keys)
  d$mlog <- -log10(d$padj)
  d$label <- paste0(wrap_terms(d$description, 50), "  (", d$cell_type, ", ", d$time_point, ")")
  d <- d[order(d$theme, d$mlog), , drop = FALSE]
  d$label <- factor(d$label, levels = unique(d$label))
  d$theme <- factor(d$theme, levels = unique(themes$theme))
  ggplot2::ggplot(d, ggplot2::aes(x = mlog, y = label, colour = database)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = mlog, yend = label), linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(shape = method), size = 2.4) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "#1c2b24") +
    ggplot2::scale_colour_manual(values = db_colours, name = NULL) +
    ggplot2::scale_shape_manual(values = method_shapes, name = NULL) +
    ggplot2::facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::labs(title = title, subtitle = subtitle, x = "-log10 BH-adjusted p", y = NULL) +
    theme_slide() +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7)
    )
}

#' Manhattan plot of the ORA results, one panel per type and time.
#'
#' Terms ordered by database then identifier on x; y: -log10 p; dashed
#' line at the largest p that passes BH 0.05 in the panel.
#'
#' @param enr All ORA terms (`enrichment.csv`, `method == "ORA"`).
#' @param n_label Terms labelled per panel.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_enrichment_manhattan <- function(enr, n_label = 2L, title = NULL, subtitle = NULL) {
  d <- enr[enr$method == "ORA", , drop = FALSE]
  d$panel <- paste0(d$time_point, ": ", d$cell_type)
  d$database <- factor(d$database, levels = names(db_colours))
  d <- d[order(d$database, d$id), , drop = FALSE]
  d$x <- stats::ave(seq_len(nrow(d)), d$panel, FUN = seq_along)
  d$mlog <- -log10(d$pval)
  cut <- do.call(rbind, lapply(split(d, d$panel), function(r) {
    p <- r$pval[r$padj <= 0.05]
    data.frame(panel = r$panel[1L], y = if (length(p)) -log10(max(p)) else NA_real_)
  }))
  lab <- do.call(rbind, lapply(split(d, d$panel), function(r) {
    utils::head(r[order(r$pval), , drop = FALSE], n_label)
  }))
  lab$label <- wrap_terms(lab$description, 28)
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = mlog, colour = database)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::geom_hline(data = cut, ggplot2::aes(yintercept = y), linetype = "dashed", colour = "#1c2b24") +
    ggrepel::geom_text_repel(
      data = lab, ggplot2::aes(label = label), size = 2.1, colour = "#1c2b24",
      max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, seed = 1
    ) +
    ggplot2::scale_colour_manual(values = db_colours, name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(size = 2))) +
    ggplot2::facet_wrap(~panel, ncol = 4, scales = "free") +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Terms, by database", y = "-log10 p (ORA)") +
    theme_slide() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 8))
}

#' Volcano of the GSEA results: NES against -log10 adjusted p.
#'
#' @param enr All GSEA terms.
#' @param n_label Terms labelled per panel (smallest adjusted p).
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_gsea_volcano <- function(enr, n_label = 3L, title = NULL, subtitle = NULL) {
  d <- enr[enr$method == "GSEA", , drop = FALSE]
  d$panel <- paste0(d$time_point, ": ", d$cell_type)
  d$mlog <- -log10(d$padj)
  d$sig <- d$padj <= 0.05
  lab <- do.call(rbind, lapply(split(d[d$sig, ], d$panel[d$sig]), function(r) {
    utils::head(r[order(r$padj, -abs(r$nes)), , drop = FALSE], n_label)
  }))
  lab$label <- wrap_terms(lab$description, 28)
  ggplot2::ggplot(d, ggplot2::aes(x = nes, y = mlog)) +
    ggplot2::geom_point(ggplot2::aes(colour = ifelse(sig, as.character(database), "n.s.")), size = 0.6, alpha = 0.6) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#1c2b24") +
    ggplot2::geom_vline(xintercept = 0, colour = "#1c2b24", linewidth = 0.3) +
    ggrepel::geom_text_repel(
      data = lab, ggplot2::aes(label = label), size = 2.1, colour = "#1c2b24",
      max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, seed = 1
    ) +
    ggplot2::scale_colour_manual(values = c(db_colours, "n.s." = "#d9d3c7"), name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::facet_wrap(~panel, ncol = 4) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Normalised enrichment score (one-versus-average Welch t)",
                  y = "-log10 BH-adjusted p") +
    theme_slide() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
}

#' Semantic map of the significant GO BP terms.
#'
#' Classical MDS of 1 - semantic similarity (`simplifyEnrichment`), points
#' sized by the number of type x time sets in which the term is
#' significant, coloured by cluster (the largest `n_colour` clusters; the
#' others grey), with the most frequent term of each coloured cluster
#' labelled.
#'
#' @param clusters `go_clusters.csv`.
#' @param sig Significant terms.
#' @param n_colour Clusters coloured.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_semantic_map <- function(clusters, sig, n_colour = 8L, title = NULL, subtitle = NULL) {
  go <- sig[sig$database == "GO BP", , drop = FALSE]
  u <- unique(go[, c("id", "description", "time_point", "cell_type")])
  freq <- stats::aggregate(list(n = rep(1L, nrow(u))), u[, c("id", "description")], sum)
  d <- merge(clusters, freq, by = "id")
  big <- names(sort(table(d$cluster), decreasing = TRUE))[seq_len(min(n_colour, length(unique(d$cluster))))]
  d$cl <- ifelse(as.character(d$cluster) %in% big, paste("cluster", d$cluster), "other")
  pal <- c(stats::setNames(colorspace::qualitative_hcl(length(big), palette = "Dark 3"), paste("cluster", big)), other = "#d9d3c7")
  lab <- do.call(rbind, lapply(split(d[d$cl != "other", ], d$cl[d$cl != "other"]), function(r) r[which.max(r$n), , drop = FALSE]))
  lab$label <- wrap_terms(lab$description, 26)
  ggplot2::ggplot(d, ggplot2::aes(x = mds1, y = mds2, colour = cl)) +
    ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.7) +
    ggrepel::geom_label_repel(
      data = lab, ggplot2::aes(label = label), size = 2.4, fill = "white", label.size = 0,
      max.overlaps = Inf, min.segment.length = 0, seed = 1, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_size_area(max_size = 5, name = "Type x time sets") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle, x = "MDS 1 (1 - semantic similarity)", y = "MDS 2") +
    theme_slide()
}
