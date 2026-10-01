#' Largest absolute off-diagonal entry of a precision.
#'
#' @param omega Square precision matrix.
#' @return A non-negative scalar. Zero when `omega` has no off-diagonal
#'   entries.
offdiag_abs_max <- function(omega) {
  off <- omega[upper.tri(omega)]
  if (length(off) == 0L || !any(is.finite(off))) {
    return(0)
  }
  max(abs(off), na.rm = TRUE)
}

#' Whether any precision on the path has an off-diagonal entry.
#'
#' @param fit One cell-type fit, with `fit$path$omega`.
#' @param eps Absolute value below which an entry counts as zero.
#'   Default `1e-8`.
#' @return `TRUE` when at least one penalty has a non-zero edge.
fit_has_offdiag <- function(fit, eps = 1e-8) {
  any(vapply(
    fit$path$omega,
    function(omega) offdiag_abs_max(omega) > eps,
    logical(1)
  ))
}

#' Drop cell types whose precision path is purely diagonal.
#'
#' @param fits Named list of cell-type fits for one time point.
#' @param eps Passed to [fit_has_offdiag()].
#' @return `fits` with diagonal-only types removed.
drop_diagonal_fits <- function(fits, eps = 1e-8) {
  fits[vapply(fits, fit_has_offdiag, logical(1), eps = eps)]
}

#' Strongest direct edges along a regularisation path.
#'
#' An edge qualifies when its precision entry is non-zero for at least
#' one penalty. The score is the largest absolute entry of that edge
#' along the path. The returned edges are the `n_edges` largest scores.
#'
#' @param path List with `path$omega`, one precision per penalty.
#' @param n_edges Number of edges to keep. Default `10`.
#' @param eps Absolute value below which an entry counts as zero.
#'   Default `1e-8`.
#' @return A data frame with columns `i` and `j` (upper-triangle
#'   indices). Zero rows when every precision is diagonal.
pick_strongest_edges <- function(path, n_edges = 10L, eps = 1e-8) {
  omega0 <- path$omega[[1L]]
  ut <- which(upper.tri(omega0), arr.ind = TRUE)
  strength <- rep(0, nrow(ut))
  for (omega in path$omega) {
    strength <- pmax(strength, abs(omega[ut]))
  }
  keep <- which(strength > eps)
  if (length(keep) == 0L) {
    return(data.frame(i = integer(0), j = integer(0)))
  }
  keep <- keep[order(strength[keep], decreasing = TRUE)]
  keep <- keep[seq_len(min(n_edges, length(keep)))]
  data.frame(
    i = ut[keep, 1L],
    j = ut[keep, 2L],
    stringsAsFactors = FALSE
  )
}

path_edge_long <- function(path, edges, genes = NULL) {
  rows <- vector("list", nrow(edges) * length(path$lambda))
  k <- 1L
  for (e in seq_len(nrow(edges))) {
    ii <- edges$i[[e]]
    jj <- edges$j[[e]]
    lab <- if (is.null(genes)) {
      paste0(ii, "-", jj)
    } else {
      paste(genes[[ii]], genes[[jj]], sep = "–")
    }
    for (t in seq_along(path$lambda)) {
      rows[[k]] <- data.frame(
        lambda = path$lambda[[t]],
        partial = path$omega[[t]][ii, jj],
        edge = lab,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

save_paged_ggplot <- function(plots, filename, width = 11, height = 8.5) {
  grDevices::pdf(filename, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in plots) {
    print(p)
  }
  invisible(filename)
}

#' Regularisation paths of the strongest direct edges.
#'
#' Diagonal-only cell types are omitted. Each panel shows the ten edges
#' with the largest absolute precision anywhere on the path, among
#' edges that are non-zero for at least one penalty.
#'
#' @inheritParams plot_ebic_vs_lambda results filename colours
#' @param n_edges Number of strongest edges per cell type. Default `10`.
#' @return `filename`, invisibly.
plot_regularisation_paths <- function(
  results,
  filename,
  colours,
  n_edges = 10L
) {
  pages <- list()
  for (time_level in names(results)) {
    pieces <- list()
    fits <- drop_diagonal_fits(results[[time_level]])
    types <- names(fits)
    for (tp in types) {
      path <- fits[[tp]]$path
      genes <- fits[[tp]]$genes
      edges <- pick_strongest_edges(path, n_edges)
      if (nrow(edges) == 0L) {
        next
      }
      dat <- path_edge_long(path, edges, genes)
      dat$cell_type <- tp
      pieces[[tp]] <- dat
    }
    if (length(pieces) == 0L) {
      next
    }
    dat <- do.call(rbind, pieces)
    dat$cell_type <- factor(dat$cell_type, levels = types)
    pages[[time_level]] <- ggplot2::ggplot(
      dat,
      ggplot2::aes(lambda, partial, group = edge, colour = edge)
    ) +
      ggplot2::geom_line() +
      ggplot2::scale_x_log10() +
      ggplot2::geom_rug(sides = "b", alpha = 0.4) +
      ggplot2::facet_wrap(~cell_type, scales = "free_y") +
      ggplot2::labs(
        title = paste("Strongest direct edges,", time_level),
        x = expression(lambda),
        y = "Precision entry",
        colour = NULL
      ) +
      ggplot2::theme_minimal()
  }
  save_paged_ggplot(pages, filename)
}

#' EBIC against log penalty, one page per time.
#'
#' Diagonal-only cell types are omitted. Both axes use a logarithmic
#' scale. EBIC changes sign, so the vertical axis is a signed log
#' (`asinh`), which keeps the minimum visible. The score is
#' \eqn{\mathrm{EBIC} = -2\ell + |E|\log n + 4\gamma|E|\log p}
#' with \eqn{\ell = n/2 (\log\det\Omega - \mathrm{tr}(S\Omega))} and
#' \eqn{\gamma = 0.5}. Smaller EBIC is preferred.
#'
#' @param results Named list of time points. Each time point is a named
#'   list of cell-type fits. Every fit has a `path` element with
#'   `lambda`, `ebic`, `density`, `priority`, `n_edges`, `omega`, and
#'   `selected`, as returned by `ebic_adaptive_path()`.
#' @param filename PDF path. One page per time point.
#' @param colours Named character vector of colours. Names are the
#'   cell-type labels used in `results`.
#' @return `filename`, invisibly.
plot_ebic_vs_lambda <- function(results, filename, colours) {
  if (!requireNamespace("gridmicrotex", quietly = TRUE)) {
    stop("Package gridmicrotex is required for the EBIC formula.")
  }
  ebic_formula <- paste0(
    "$\\mathrm{EBIC} = -2\\ell + |E|\\log n + 4\\gamma|E|\\log p$, ",
    "$\\ell = \\frac{n}{2}(\\log\\det\\Omega - \\mathrm{tr}(S\\Omega))$, ",
    "$\\gamma = 0.5$.<br>**Minimise EBIC.**"
  )
  pages <- list()
  for (time_level in names(results)) {
    fits <- drop_diagonal_fits(results[[time_level]])
    types <- names(fits)
    if (length(types) == 0L) {
      next
    }
    pieces <- lapply(types, function(tp) {
      path <- fits[[tp]]$path
      data.frame(
        lambda = path$lambda,
        ebic = path$ebic,
        cell_type = tp,
        stringsAsFactors = FALSE
      )
    })
    dat <- do.call(rbind, pieces)
    dat <- dat[is.finite(dat$ebic) & dat$lambda > 0, , drop = FALSE]
    dat$is_min <- FALSE
    for (tp in types) {
      idx <- which(dat$cell_type == tp)
      if (length(idx) == 0L) {
        next
      }
      dat$is_min[idx[which.min(dat$ebic[idx])]] <- TRUE
    }
    dat$cell_type <- factor(dat$cell_type, levels = types)
    dat$point_label <- sprintf(
      "\u03bb = %s\nEBIC = %s",
      signif(dat$lambda, 3),
      signif(dat$ebic, 4)
    )
    pages[[time_level]] <- ggplot2::ggplot(
      dat,
      ggplot2::aes(lambda, ebic, colour = cell_type)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = dat[!dat$is_min, ], size = 1.2) +
      ggplot2::geom_point(data = dat[dat$is_min, ], size = 3.5) +
      ggplot2::geom_text(
        data = dat[!dat$is_min, ],
        ggplot2::aes(label = point_label),
        size = 1.8,
        vjust = -0.4,
        lineheight = 0.85,
        show.legend = FALSE
      ) +
      ggplot2::geom_label(
        data = dat[dat$is_min, ],
        ggplot2::aes(label = point_label),
        size = 2.6,
        vjust = 1.2,
        lineheight = 0.9,
        show.legend = FALSE
      ) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_continuous(trans = "pseudo_log") +
      ggplot2::geom_rug(sides = "b", alpha = 0.4) +
      ggplot2::scale_colour_manual(values = colours, guide = "none") +
      ggplot2::facet_wrap(~cell_type, scales = "free") +
      ggplot2::labs(
        title = paste("EBIC against log penalty,", time_level),
        subtitle = ebic_formula,
        x = expression(lambda),
        y = "EBIC (signed log)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.subtitle = gridmicrotex::element_markdown(
          fontsize = 9,
          width = grid::unit(1, "npc")
        )
      )
  }
  grDevices::cairo_pdf(filename, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in pages) {
    print(p)
  }
  invisible(filename)
}

#' Undirected edge density against EBIC, one page per time.
#'
#' Diagonal-only cell types are omitted. Density is the number of
#' edges divided by `G (G - 1) / 2`. The larger point is the selected
#' penalty.
#'
#' @inheritParams plot_ebic_vs_lambda
#' @return `filename`, invisibly.
plot_edge_count_vs_ebic <- function(results, filename, colours) {
  pages <- list()
  for (time_level in names(results)) {
    fits <- drop_diagonal_fits(results[[time_level]])
    types <- names(fits)
    if (length(types) == 0L) {
      next
    }
    pieces <- lapply(types, function(tp) {
      path <- fits[[tp]]$path
      data.frame(
        density = path$density,
        ebic = path$ebic,
        selected = seq_along(path$lambda) == path$selected,
        cell_type = tp,
        stringsAsFactors = FALSE
      )
    })
    dat <- do.call(rbind, pieces)
    dat$cell_type <- factor(dat$cell_type, levels = types)
    pages[[time_level]] <- ggplot2::ggplot(
      dat,
      ggplot2::aes(density, ebic, colour = cell_type)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = dat[dat$selected, ], size = 2) +
      ggplot2::scale_colour_manual(values = colours, guide = "none") +
      ggplot2::facet_wrap(~cell_type, scales = "free") +
      ggplot2::labs(
        title = paste("Edge density against EBIC,", time_level),
        x = "Undirected density (edges / [G(G-1)/2])",
        y = "EBIC"
      ) +
      ggplot2::theme_minimal()
  }
  save_paged_ggplot(pages, filename)
}

#' One covariance or precision heatmap, with null edges in grey.
#'
#' The fill limits are the minimum and maximum of the diagonal
#' (marginal variances, or partial variances on the precision).
#' Off-diagonal zeros are drawn in grey and do not use that scale.
#'
#' @param mat Square covariance or precision.
#' @param panel_title Panel title.
#' @param eps Absolute value below which an off-diagonal entry is null.
#' @return A ggplot.
matrix_heatmap <- function(mat, panel_title, eps = 1e-8) {
  p <- nrow(mat)
  dat <- data.frame(
    i = as.vector(row(mat)),
    j = as.vector(col(mat)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  dat$null_off <- dat$i != dat$j &
    (!is.finite(dat$value) | abs(dat$value) <= eps)
  diag_values <- diag(mat)
  diag_values <- diag_values[is.finite(diag_values)]
  limits <- range(diag_values)
  if (!all(is.finite(limits)) || diff(limits) == 0) {
    limits <- c(0, 1)
  }
  shown <- dat[!dat$null_off, , drop = FALSE]
  hidden <- dat[dat$null_off, , drop = FALSE]
  ggplot2::ggplot(shown, ggplot2::aes(j, i, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::geom_tile(
      data = hidden,
      mapping = ggplot2::aes(j, i),
      fill = "grey35",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_viridis_c(
      limits = limits,
      oob = scales::squish,
      na.value = "grey35"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = panel_title, x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

#' Covariance and precision heatmaps, one page per time.
#'
#' Every cell type is included, including those whose precision is
#' diagonal. Each cell type is one row: covariance beside precision.
#' A null off-diagonal entry is grey. The colour scale of each panel
#' runs from the smallest to the largest diagonal entry.
#'
#' @inheritParams plot_ebic_vs_lambda results filename
#' @return `filename`, invisibly.
plot_covariance_precision_heatmaps <- function(results, filename) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package patchwork is required for the matrix heatmaps.")
  }
  n_max <- max(1L, vapply(results, length, integer(1)))
  grDevices::cairo_pdf(
    filename,
    width = 8,
    height = 3.8 * n_max,
    onefile = TRUE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  drew <- FALSE
  for (time_level in names(results)) {
    fits <- results[[time_level]]
    rows <- list()
    for (tp in names(fits)) {
      item <- fits[[tp]]
      cov_plot <- matrix_heatmap(item$sigma, paste(tp, "covariance"))
      prec_plot <- matrix_heatmap(item$omega, paste(tp, "precision"))
      rows[[tp]] <- cov_plot + prec_plot
    }
    if (length(rows) == 0L) {
      next
    }
    page <- patchwork::wrap_plots(rows, ncol = 1) +
      patchwork::plot_annotation(
        title = paste("Covariance and precision,", time_level)
      )
    print(page)
    drew <- TRUE
  }
  if (!drew) {
    plot.new()
  }
  invisible(filename)
}

#' One row per selected network.
#'
#' The selected penalty is `path$selected`. Degree treats a precision
#' entry as an edge when its absolute value exceeds `1e-8`.
#'
#' @inheritParams plot_ebic_vs_lambda results
#' @param estimator Character label written to the `estimator` column,
#'   for example `"huge"`.
#' @return A data frame with one row per time point and cell type.
selected_metric_table <- function(results, estimator) {
  rows <- list()
  for (time_level in names(results)) {
    for (tp in names(results[[time_level]])) {
      item <- results[[time_level]][[tp]]
      path <- item$path
      sel <- path$selected
      omega <- path$omega[[sel]]
      deg <- rowSums(abs(omega) > 1e-8) - 1
      rows[[length(rows) + 1L]] <- data.frame(
        time_point = time_level,
        cell_type = tp,
        estimator = estimator,
        ebic = path$ebic[[sel]],
        n_edges = path$n_edges[[sel]],
        density = path$density[[sel]],
        mean_degree = mean(deg),
        hellinger_empirical = item$hellinger_empirical,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Heat map of the selected-network metrics.
#'
#' Uses `funkyheatmap` when that package is installed. Otherwise writes
#' a faceted tile plot of the same columns.
#'
#' @inheritParams plot_ebic_vs_lambda filename
#' @param metric_table Data frame returned by [selected_metric_table()].
#' @return `filename`, invisibly.
plot_selected_funkyheatmap <- function(metric_table, filename) {
  if (!requireNamespace("funkyheatmap", quietly = TRUE)) {
    warning("funkyheatmap is not installed; writing a tile table instead.")
    metrics <- c(
      "ebic",
      "n_edges",
      "density",
      "mean_degree",
      "hellinger_empirical"
    )
    long <- do.call(
      rbind,
      lapply(metrics, function(m) {
        data.frame(
          time_point = metric_table$time_point,
          cell_type = metric_table$cell_type,
          estimator = metric_table$estimator,
          metric = m,
          value = metric_table[[m]],
          stringsAsFactors = FALSE
        )
      })
    )
    p <- ggplot2::ggplot(
      long,
      ggplot2::aes(metric, cell_type, fill = value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_grid(time_point ~ estimator, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(
      filename,
      p,
      width = 10,
      height = 8,
      dpi = 320,
      limitsize = FALSE
    )
    return(invisible(filename))
  }
  pal <- funkyheatmap::funky_heatmap(
    data = metric_table,
    expand = list(xmax = 4)
  )
  ggplot2::ggsave(
    filename,
    pal,
    width = 12,
    height = 8,
    dpi = 320,
    limitsize = FALSE
  )
}

#' Pairwise Gaussian scores among cell types, one page per time.
#'
#' Each page is one time point. Rows and columns are cell types. The
#' fill is one column of `pair_table`.
#'
#' @inheritParams plot_ebic_vs_lambda filename
#' @param pair_table Data frame with `time_point`, `estimator`,
#'   `type_a`, `type_b`, and the column named by `score_col`.
#' @param score_col Name of the numeric column to map to the tile fill,
#'   either `"hellinger"` or `"mixsim_overlap"`.
#' @param title Stem of the page title. The time point is appended.
#' @return `filename`, invisibly.
plot_pairwise_gaussian_tiles <- function(
  pair_table,
  filename,
  score_col,
  title
) {
  pages <- list()
  times <- unique(pair_table$time_point)
  for (time_level in times) {
    dat <- pair_table[pair_table$time_point == time_level, , drop = FALSE]
    pages[[time_level]] <- ggplot2::ggplot(
      dat,
      ggplot2::aes(type_a, type_b, fill = .data[[score_col]])
    ) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~estimator) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      ) +
      ggplot2::labs(
        title = paste(title, time_level),
        x = NULL,
        y = NULL,
        fill = score_col
      )
  }
  save_paged_ggplot(pages, filename)
}
