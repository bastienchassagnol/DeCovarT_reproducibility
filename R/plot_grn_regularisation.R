#' Pick pedagogical edges for a regularisation-path plot.
#'
#' Five edges that vanish at the selected penalty, and five that remain.
pick_path_edges <- function(path, n_keep = 5L, n_drop = 5L, seed = 1L, eps = 1e-8) {
  sel <- path$selected
  omega_sel <- path$omega[[sel]]
  omega_dense <- path$omega[[which.min(path$lambda)]]
  ut <- which(upper.tri(omega_sel), arr.ind = TRUE)
  mag_sel <- abs(omega_sel[ut])
  mag_dense <- abs(omega_dense[ut])
  keep_id <- which(mag_sel > eps)
  drop_id <- which(mag_dense > eps & mag_sel <= eps)
  withr::with_seed(seed, {
    keep_id <- if (length(keep_id) > n_keep) sample(keep_id, n_keep) else keep_id
    drop_id <- if (length(drop_id) > n_drop) sample(drop_id, n_drop) else drop_id
  })
  rbind(
    data.frame(
      i = ut[keep_id, 1],
      j = ut[keep_id, 2],
      status = "retained",
      stringsAsFactors = FALSE
    ),
    data.frame(
      i = ut[drop_id, 1],
      j = ut[drop_id, 2],
      status = "driven to zero",
      stringsAsFactors = FALSE
    )
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
        status = edges$status[[e]],
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

#' Regularisation paths: five retained and five vanishing edges.
#'
#' One page per time; all cell types of that time are faceted.
plot_regularisation_paths <- function(
  results,
  filename,
  colours,
  n_keep = 5L,
  n_drop = 5L,
  seed = 1L
) {
  pages <- list()
  for (time_level in names(results)) {
    pieces <- list()
    types <- names(results[[time_level]])
    for (tp in types) {
      path <- results[[time_level]][[tp]]$path
      genes <- results[[time_level]][[tp]]$genes
      edges <- pick_path_edges(path, n_keep, n_drop, seed)
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
      ggplot2::aes(lambda, partial, group = edge, colour = status)
    ) +
      ggplot2::geom_line() +
      ggplot2::scale_x_log10() +
      ggplot2::geom_rug(sides = "b", alpha = 0.4) +
      ggplot2::facet_wrap(~cell_type, scales = "free_y") +
      ggplot2::labs(
        title = paste("Regularisation paths,", time_level),
        x = expression(lambda),
        y = "Precision entry",
        colour = NULL
      ) +
      ggplot2::theme_minimal()
  }
  save_paged_ggplot(pages, filename)
}

#' EBIC against log lambda, one page per time.
plot_ebic_vs_lambda <- function(results, filename, colours) {
  pages <- list()
  for (time_level in names(results)) {
    types <- names(results[[time_level]])
    pieces <- lapply(types, function(tp) {
      path <- results[[time_level]][[tp]]$path
      data.frame(
        lambda = path$lambda,
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
      ggplot2::aes(lambda, ebic, colour = cell_type)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = dat[dat$selected, ], size = 2) +
      ggplot2::scale_x_log10() +
      ggplot2::geom_rug(sides = "b", alpha = 0.4) +
      ggplot2::scale_colour_manual(values = colours, guide = "none") +
      ggplot2::facet_wrap(~cell_type, scales = "free_y") +
      ggplot2::labs(
        title = paste("EBIC against log lambda,", time_level),
        x = expression(lambda),
        y = "EBIC"
      ) +
      ggplot2::theme_minimal()
  }
  save_paged_ggplot(pages, filename)
}

#' Normalised support (undirected density) against EBIC.
plot_edge_count_vs_ebic <- function(results, filename, colours) {
  pages <- list()
  for (time_level in names(results)) {
    types <- names(results[[time_level]])
    pieces <- lapply(types, function(tp) {
      path <- results[[time_level]][[tp]]$path
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

#' Weighted interval priority along the path.
plot_priority_score <- function(results, filename, colours) {
  pages <- list()
  for (time_level in names(results)) {
    types <- names(results[[time_level]])
    pieces <- lapply(types, function(tp) {
      path <- results[[time_level]][[tp]]$path
      data.frame(
        lambda = path$lambda,
        priority = path$priority,
        cell_type = tp,
        stringsAsFactors = FALSE
      )
    })
    dat <- do.call(rbind, pieces)
    dat <- dat[is.finite(dat$priority), , drop = FALSE]
    dat$cell_type <- factor(dat$cell_type, levels = types)
    pages[[time_level]] <- ggplot2::ggplot(
      dat,
      ggplot2::aes(lambda, priority, colour = cell_type)
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_x_log10() +
      ggplot2::geom_rug(sides = "b", alpha = 0.4) +
      ggplot2::scale_colour_manual(values = colours, guide = "none") +
      ggplot2::facet_wrap(~cell_type, scales = "free_y") +
      ggplot2::labs(
        title = paste("Interval priority (EBIC + support),", time_level),
        x = expression(lambda),
        y = "Priority score"
      ) +
      ggplot2::theme_minimal()
  }
  save_paged_ggplot(pages, filename)
}

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

#' Folder-1 table: one selected network per type and time.
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

#' Pairwise Hellinger or MixSim overlap among types, one page per time.
plot_pairwise_gaussian_tiles <- function(pair_table, filename, score_col, title) {
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
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(
        title = paste(title, time_level),
        x = NULL,
        y = NULL,
        fill = score_col
      )
  }
  save_paged_ggplot(pages, filename)
}
