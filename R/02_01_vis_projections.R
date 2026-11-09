#' PCA and PLS-DA of one time slot, coloured by cell type.
#'
#' Descriptive projections of the cells that enter tier (ii). They are
#' not a selection rule. Both use `log1p` of counts divided by the
#' size factor. PCA centres and scales genes. PLS-DA is
#' `mixOmics::plsda()` on the same log matrix, two components, with
#' the cell type as the response. `mixOmics` is not in `renv.lock`.
#'
#' @param counts Genes by cells, symbols as row names.
#' @param cell_type Cell-type label per column.
#' @param size_factors Size factor per cell. Default: the column sums.
#' @param colours Optional named vector of colours, names equal to the
#'   cell-type labels.
#' @param ncomp Components requested from PLS-DA. Default `2`.
#' @return A list with `pca` and `plsda` (ggplot objects) and `n_genes`
#'   (genes with positive variance that entered both fits).
plot_type_pca_plsda <- function(
  counts,
  cell_type,
  size_factors = NULL,
  colours = NULL,
  ncomp = 2L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for the PCA and PLS-DA scatters.")
  }
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop(
      "Package mixOmics is required for PLS-DA (not in renv.lock). ",
      "Install it with BiocManager::install(\"mixOmics\")."
    )
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
  if (is.null(size_factors)) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package Matrix is required to sum count columns.")
    }
    size_factors <- Matrix::colSums(counts)
  }
  if (length(size_factors) != ncol(counts)) {
    stop("`size_factors` must have one value per column of `counts`.")
  }
  if (any(!is.finite(size_factors) | size_factors <= 0)) {
    stop("`size_factors` must be positive and finite.")
  }
  log_counts <- log1p(sweep(as.matrix(counts), 2L, size_factors, "/"))
  x <- t(log_counts)
  gene_sd <- apply(x, 2L, stats::sd)
  keep <- is.finite(gene_sd) & gene_sd > 0
  x <- x[, keep, drop = FALSE]
  if (nrow(x) < 3L || ncol(x) < 2L) {
    stop("PCA and PLS-DA need at least three cells and two varying genes.")
  }
  if (ncomp < 2L) {
    stop("`ncomp` must be at least 2 so both scatters have two axes.")
  }
  pc <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  var <- pc$sdev^2
  pct <- 100 * var / sum(var)
  pca_df <- data.frame(
    axis1 = pc$x[, 1L],
    axis2 = pc$x[, 2L],
    type_label = cell_type,
    stringsAsFactors = FALSE
  )
  fit <- mixOmics::plsda(X = x, Y = factor(cell_type), ncomp = ncomp)
  scores <- as.matrix(fit$variates$X)
  if (ncol(scores) < 2L) {
    stop("plsda returned fewer than two components.")
  }
  plsda_df <- data.frame(
    axis1 = scores[, 1L],
    axis2 = scores[, 2L],
    type_label = cell_type,
    stringsAsFactors = FALSE
  )
  explained <- fit$prop_expl_var$X
  plsda_x <- "Component 1"
  plsda_y <- "Component 2"
  if (!is.null(explained) && length(explained) >= 2L) {
    plsda_x <- paste0(
      "Component 1 (",
      format(round(100 * explained[[1L]], 1), nsmall = 1),
      "%)"
    )
    plsda_y <- paste0(
      "Component 2 (",
      format(round(100 * explained[[2L]], 1), nsmall = 1),
      "%)"
    )
  }
  draw <- function(df, x_lab, y_lab, title) {
    plot <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = axis1, y = axis2, colour = type_label)
    ) +
      ggplot2::geom_point(alpha = 0.45, size = 0.7) +
      ggplot2::labs(
        x = x_lab,
        y = y_lab,
        colour = "Cell type",
        title = title
      ) +
      ggplot2::theme_minimal(base_size = 11)
    if (!is.null(colours)) {
      plot <- plot + ggplot2::scale_colour_manual(values = colours)
    }
    plot
  }
  list(
    pca = draw(
      pca_df,
      paste0(
        "PC1 (",
        format(round(pct[[1L]], 1), nsmall = 1),
        "%)"
      ),
      paste0(
        "PC2 (",
        format(round(pct[[2L]], 1), nsmall = 1),
        "%)"
      ),
      "PCA of size-factor log counts, by cell type"
    ),
    plsda = draw(
      plsda_df,
      plsda_x,
      plsda_y,
      "PLS-DA of size-factor log counts, by cell type"
    ),
    n_genes = ncol(x)
  )
}
