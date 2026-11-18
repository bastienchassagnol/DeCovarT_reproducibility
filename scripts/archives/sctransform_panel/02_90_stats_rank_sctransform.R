#' Cell-weighted between/within sums of squares of Pearson residuals.
#'
#' This is the one-way ANOVA split used in the naive-screen chapter:
#' each barcode has weight one, so a large type dominates `ss_between`
#' and `ss_within`.
#'
#' @param residual Numeric matrix, genes in rows and cells in columns.
#' @param groups Factor or character vector of group labels, one per
#'   column of `residual`.
#' @return A list with `ss_between`, `ss_within` (length `nrow(residual)`),
#'   and `means` (genes by groups).
#' @details
#'   For gene \(g\) and groups \(j=1,\ldots,J\) with sizes \(C_j\),
#'   \[
#'   \mathrm{SS}_{\mathrm{between},g}
#'   =
#'   \sum_{j} C_j (\bar r_{gj}-\bar r_{g})^{2},
#'   \qquad
#'   \mathrm{SS}_{\mathrm{within},g}
#'   =
#'   \sum_{j}\sum_{c_j}(r_{c_j g}-\bar r_{gj})^{2}.
#'   \]
#' @seealso [anova_ms_equal_type()]
anova_ss <- function(residual, groups) {
  groups <- droplevels(factor(groups))
  if (nlevels(groups) < 2L) {
    stop("Need at least two groups for a between/within split.")
  }
  grand <- rowMeans(residual)
  ss_between <- numeric(nrow(residual))
  ss_within <- numeric(nrow(residual))
  means <- matrix(
    NA_real_,
    nrow = nrow(residual),
    ncol = nlevels(groups),
    dimnames = list(rownames(residual), levels(groups))
  )
  for (group in levels(groups)) {
    idx <- which(groups == group)
    group_mean <- rowMeans(residual[, idx, drop = FALSE])
    means[, group] <- group_mean
    ss_between <- ss_between + length(idx) * (group_mean - grand)^2
    deviation <- residual[, idx, drop = FALSE] - group_mean
    ss_within <- ss_within + rowSums(deviation * deviation)
  }
  list(
    ss_between = ss_between,
    ss_within = ss_within,
    means = means
  )
}

#' Type-equal between/within mean squares of Pearson residuals.
#'
#' Each of the \(J\) types has weight \(1/J\), so a large type cannot
#' dominate the rank. The returned `ss_between` and `ss_within` are
#' mean squares, not ANOVA sums of squares.
#'
#' @param residual Numeric matrix, genes in rows and cells in columns.
#' @param groups Factor or character vector of group labels, one per
#'   column of `residual`.
#' @return A list with `ss_between` (variance of the \(J\) type means),
#'   `ss_within` (mean of the \(J\) within-type variances), and `means`.
#' @details
#'   Let \(\bar r_{gj}\) be the mean residual of gene \(g\) in type
#'   \(j\) and \(v_{gj}\) the mean squared deviation inside that type.
#'   Then
#'   \[
#'   \bar r_{g}^{\mathrm{eq}} = J^{-1}\sum_{j}\bar r_{gj},
#'   \quad
#'   \mathrm{MS}_{\mathrm{between},g}
#'   =
#'   J^{-1}\sum_{j}(\bar r_{gj}-\bar r_{g}^{\mathrm{eq}})^{2},
#'   \quad
#'   \mathrm{MS}_{\mathrm{within},g}
#'   =
#'   J^{-1}\sum_{j} v_{gj}.
#'   \]
#' @seealso [anova_ss()]
anova_ms_equal_type <- function(residual, groups) {
  groups <- droplevels(factor(groups))
  if (nlevels(groups) < 2L) {
    stop("Need at least two groups for a between/within split.")
  }
  group_levels <- levels(groups)
  means <- matrix(
    NA_real_,
    nrow = nrow(residual),
    ncol = length(group_levels),
    dimnames = list(rownames(residual), group_levels)
  )
  variances <- means
  for (group in group_levels) {
    idx <- which(groups == group)
    group_mean <- rowMeans(residual[, idx, drop = FALSE])
    means[, group] <- group_mean
    deviation <- residual[, idx, drop = FALSE] - group_mean
    variances[, group] <- rowSums(deviation * deviation) / length(idx)
  }
  grand <- rowMeans(means)
  list(
    ss_between = rowMeans((means - grand)^2),
    ss_within = rowMeans(variances),
    means = means
  )
}

#' Between/within ratio and share for a pair of residual scores.
#'
#' @param ss_between Numeric vector of between-type scores (sums of
#'   squares or mean squares).
#' @param ss_within Numeric vector of within-type scores, same length.
#' @return A list with `between_over_within` and `between_over_total`.
#'   Both are `NA` when `ss_between + ss_within` is not positive.
rank_ratio <- function(ss_between, ss_within) {
  total <- ss_between + ss_within
  keep <- total > 0
  ratio <- rep(NA_real_, length(ss_between))
  share <- rep(NA_real_, length(ss_between))
  ratio[keep] <- ss_between[keep] /
    pmax(
      ss_within[keep],
      .Machine$double.eps
    )
  share[keep] <- ss_between[keep] / total[keep]
  list(between_over_within = ratio, between_over_total = share)
}

#' Cell type with the largest mean residual for each gene.
#'
#' @param means Numeric matrix, genes in rows and types in columns.
#' @return Character vector of column names, one per gene. Ties take
#'   the first type in column order.
top_mean_group <- function(means) {
  idx <- max.col(means, ties.method = "first")
  colnames(means)[idx]
}

#' Empty rank table with the columns written by the naive screen.
#'
#' @return A zero-row data frame with `time_point`, `cell_type`,
#'   `gene`, the between/within scores, `top_cell_type`,
#'   `passes_between_gt_within`, `rank`, and `selected`.
empty_rank_table <- function() {
  data.frame(
    time_point = character(),
    cell_type = character(),
    gene = character(),
    ss_between = numeric(),
    ss_within = numeric(),
    between_over_within = numeric(),
    between_over_total = numeric(),
    top_cell_type = character(),
    passes_between_gt_within = logical(),
    rank = integer(),
    selected = logical(),
    stringsAsFactors = FALSE
  )
}

#' Record the SCTransform v2 arguments used for one time slot.
#'
#' @param seu_fit Seurat object after `Seurat::SCTransform()`.
#' @param time_point Label stored in the model card, typically `"48h"`.
#' @param n_cells_fit Number of cells passed to the fit.
#' @return A one-row data frame of flavour, method, offset variable,
#'   cell counts, residual gene count, and residual clip.
capture_sct_model <- function(seu_fit, time_point, n_cells_fit) {
  fallback <- data.frame(
    time_point = time_point,
    vst_flavor = NA_character_,
    method = NA_character_,
    latent_var = NA_character_,
    n_cells_param = NA_integer_,
    n_cells_fit = n_cells_fit,
    n_genes_residual = NA_integer_,
    clip = NA_character_,
    stringsAsFactors = FALSE
  )
  model <- tryCatch(
    seu_fit[["SCT"]]@SCTModel.list[[1]],
    error = function(e) NULL
  )
  if (is.null(model)) {
    return(fallback)
  }
  args <- tryCatch(model@arguments, error = function(e) list())
  residual <- tryCatch(
    assay_matrix(seu_fit, "SCT", "scale.data"),
    error = function(e) NULL
  )
  data.frame(
    time_point = time_point,
    vst_flavor = as.character(args$vst.flavor %||% NA),
    method = as.character(
      args$sct.method %||% args$method %||% NA
    ),
    latent_var = paste(args$latent_var %||% character(), collapse = ","),
    n_cells_param = as.integer(args$n_cells %||% NA_integer_),
    n_cells_fit = as.integer(n_cells_fit),
    n_genes_residual = if (is.null(residual)) {
      NA_integer_
    } else {
      nrow(residual)
    },
    clip = paste(tryCatch(model@clips, error = function(e) NA), collapse = ","),
    stringsAsFactors = FALSE
  )
}

#' Flag ranked genes that also sit in the signalling-marker panel.
#'
#' Adds `in_marker_panel`, `in_union`, and `source`. Marker genes that
#' were never scored (absent from the residual matrix) are appended
#' with missing between/within columns.
#'
#' @param ranks Rank table from the global or per-type screen.
#' @param marker_df Marker rows with `time_point`, `gene`, and, when
#'   `by_cell_type` is `TRUE`, `celltypeannotation`.
#' @param by_cell_type Logical. Match on type as well as time and gene
#'   when `TRUE` (per-type table); match on time and gene only when
#'   `FALSE` (global tables).
#' @return `ranks` with the three union columns. `source` is
#'   `"both"`, `"sctransform"`, `"marker"`, or `"screened"`.
attach_union <- function(ranks, marker_df, by_cell_type) {
  ranks$in_marker_panel <- FALSE
  if (by_cell_type) {
    id_rank <- paste(ranks$time_point, ranks$cell_type, ranks$gene, sep = "\r")
    id_mark <- paste(
      marker_df$time_point,
      marker_df$celltypeannotation,
      marker_df$gene,
      sep = "\r"
    )
    key_df <- unique(marker_df[, c(
      "time_point",
      "celltypeannotation",
      "gene"
    )])
  } else {
    id_rank <- paste(ranks$time_point, ranks$gene, sep = "\r")
    id_mark <- paste(marker_df$time_point, marker_df$gene, sep = "\r")
    key_df <- unique(marker_df[, c("time_point", "gene")])
  }
  if (nrow(key_df) > 0L) {
    ranks$in_marker_panel <- id_rank %in% unique(id_mark)
    missing_id <- unique(id_mark[!id_mark %in% id_rank])
    if (length(missing_id) > 0L) {
      if (by_cell_type) {
        bits <- strsplit(missing_id, "\r", fixed = TRUE)
        extra <- data.frame(
          time_point = vapply(bits, `[`, character(1), 1L),
          cell_type = vapply(bits, `[`, character(1), 2L),
          gene = vapply(bits, `[`, character(1), 3L),
          stringsAsFactors = FALSE
        )
      } else {
        bits <- strsplit(missing_id, "\r", fixed = TRUE)
        extra <- data.frame(
          time_point = vapply(bits, `[`, character(1), 1L),
          cell_type = NA_character_,
          gene = vapply(bits, `[`, character(1), 2L),
          stringsAsFactors = FALSE
        )
      }
      extra$ss_between <- NA_real_
      extra$ss_within <- NA_real_
      extra$between_over_within <- NA_real_
      extra$between_over_total <- NA_real_
      extra$top_cell_type <- extra$cell_type
      extra$passes_between_gt_within <- NA
      extra$rank <- NA_integer_
      extra$selected <- FALSE
      extra$in_marker_panel <- TRUE
      ranks <- rbind(ranks, extra[colnames(ranks)])
    }
  }
  ranks$in_union <- ranks$selected %in% TRUE | ranks$in_marker_panel %in% TRUE
  ranks$source <- ifelse(
    ranks$selected %in% TRUE & ranks$in_marker_panel %in% TRUE,
    "both",
    ifelse(
      ranks$selected %in% TRUE,
      "sctransform",
      ifelse(ranks$in_marker_panel %in% TRUE, "marker", "screened")
    )
  )
  ranks
}

#' Rank genes from one SCTransform v2 fit per split of the object.
#'
#' The negative-binomial mean is an intercept plus a library-size
#' offset. Cell type is not a covariate. After the fit, residuals are
#' scored twice globally (cell-weighted sums of squares, then
#' type-equal mean squares) and optionally once per type
#' (one-versus-rest).
#'
#' @param seu Seurat object. RNA `counts` are the UMI table.
#'   Metadata column `log_umi`, when present, is the v2 offset
#'   (log10 library size) and is not recomputed from a gene subset.
#' @param strategy `"global"`, `"per_cell_type"`, or `"both"`.
#'   `"both"` fits SCTransform once and scores both rules.
#' @param cell_type_col Metadata column of cell-type labels.
#' @param split_by Optional metadata column. The negative-binomial
#'   fit and the sum-of-squares split are redone in each level.
#' @param n_genes_global Cap for the global between/within ranking.
#' @param n_genes_celltype Cap, per cell type, on the one-versus-rest
#'   between/within ranking. The inequality between > within is stored
#'   and is not a hard gate: for Pearson residuals that cut is
#'   R-squared above one half, which almost no gene meets here.
#' @param exclude_cell_types Labels dropped before the fit (artefacts).
#' @param min_cells Types with fewer cells in a split are dropped
#'   from the sum-of-squares split only after the fit.
#' @param seed Passed to `SCTransform(seed.use)`.
#' @return A list with `global` (cell-weighted), `global_equal_type`
#'   (type-equal mean squares), `per_cell_type`, and `model`.
rank_genes_sctransform <- function(
  seu,
  strategy = c("global", "per_cell_type", "both"),
  cell_type_col = "celltypeannotation",
  split_by = "time_point",
  n_genes_global = 500L,
  n_genes_celltype = 50L,
  exclude_cell_types = character(),
  min_cells = 20L,
  seed = 1L
) {
  strategy <- match.arg(strategy)
  do_global <- strategy %in% c("global", "both")
  do_type <- strategy %in% c("per_cell_type", "both")
  n_genes_global <- as.integer(n_genes_global)
  n_genes_celltype <- as.integer(n_genes_celltype)
  min_cells <- as.integer(min_cells)
  meta <- seu@meta.data
  if (!cell_type_col %in% colnames(meta)) {
    stop("Column ", cell_type_col, " is missing from the Seurat metadata.")
  }
  if (!is.null(split_by) && !split_by %in% colnames(meta)) {
    stop("Column ", split_by, " is missing from the Seurat metadata.")
  }
  if ("log_umi" %in% colnames(meta)) {
    message(
      "Using metadata log_umi as the v2 library-size offset ",
      "(median log10 depth ",
      format(stats::median(meta$log_umi), digits = 4),
      ")."
    )
  } else {
    message(
      "metadata log_umi is absent; sctransform will take library ",
      "size from the genes still in the RNA assay."
    )
  }

  if (is.null(split_by)) {
    seu$sct_split <- "all"
    split_by <- "sct_split"
  }
  split_values <- as.character(seu@meta.data[[split_by]])
  levels_split <- unique(split_values)
  global_parts <- list()
  global_equal_parts <- list()
  type_parts <- list()
  model_parts <- list()

  for (level in levels_split) {
    message("SCTransform v2 at ", level, " ...")
    cells_level <- colnames(seu)[split_values == level]
    sub <- subset(x = seu, cells = cells_level)
    type_level <- as.character(sub@meta.data[[cell_type_col]])
    keep_cell <- !type_level %in% exclude_cell_types
    if (!any(keep_cell)) {
      warning(
        "No cells left at ",
        level,
        " after dropping excluded types.",
        call. = FALSE
      )
      next
    }
    sub <- subset(x = sub, cells = colnames(sub)[keep_cell])
    n_fit <- ncol(sub)
    # v2 replaces a NULL n_cells with 2,000. Pass the full split so
    # the offset model is estimated on every retained cell.
    sub <- Seurat::SCTransform(
      object = sub,
      assay = "RNA",
      new.assay.name = "SCT",
      vst.flavor = "v2",
      method = "glmGamPoi",
      vars.to.regress = NULL,
      residual.features = rownames(sub),
      return.only.var.genes = FALSE,
      do.correct.umi = FALSE,
      do.scale = FALSE,
      do.center = TRUE,
      ncells = n_fit,
      conserve.memory = FALSE,
      seed.use = as.integer(seed),
      verbose = TRUE
    )
    model_parts[[level]] <- capture_sct_model(sub, level, n_fit)
    residual <- as.matrix(assay_matrix(sub, "SCT", "scale.data"))
    type_fit <- as.character(sub@meta.data[[cell_type_col]])
    names(type_fit) <- colnames(sub)
    type_fit <- type_fit[colnames(residual)]
    n_by_type <- table(type_fit)
    rare <- names(n_by_type)[n_by_type < min_cells]
    if (length(rare) > 0L) {
      message(
        "  Dropping types with fewer than ",
        min_cells,
        " cells from the sum of squares at ",
        level,
        ": ",
        paste(rare, collapse = ", ")
      )
      keep_col <- !type_fit %in% rare
      residual <- residual[, keep_col, drop = FALSE]
      type_fit <- type_fit[keep_col]
    }
    if (length(unique(type_fit)) < 2L) {
      warning(
        "Fewer than two scored cell types at ",
        level,
        "; skipping the rank.",
        call. = FALSE
      )
      rm(sub, residual)
      gc()
      next
    }

    if (do_global) {
      # Cell-weighted global rank: each barcode has weight one.
      ss_cell <- anova_ss(residual, type_fit)
      scored_cell <- rank_ratio(ss_cell$ss_between, ss_cell$ss_within)
      order_cell <- order(
        -scored_cell$between_over_within,
        -ss_cell$ss_between,
        rownames(ss_cell$means)
      )
      rank_cell <- integer(nrow(ss_cell$means))
      rank_cell[order_cell] <- seq_along(order_cell)
      global_parts[[level]] <- data.frame(
        time_point = level,
        cell_type = NA_character_,
        gene = rownames(ss_cell$means),
        ss_between = ss_cell$ss_between,
        ss_within = ss_cell$ss_within,
        between_over_within = scored_cell$between_over_within,
        between_over_total = scored_cell$between_over_total,
        top_cell_type = top_mean_group(ss_cell$means),
        passes_between_gt_within = ss_cell$ss_between > ss_cell$ss_within,
        rank = rank_cell,
        selected = rank_cell <= n_genes_global &
          is.finite(scored_cell$between_over_within),
        stringsAsFactors = FALSE
      )

      # Type-equal global rank: each scored type has weight 1/J.
      ss_eq <- anova_ms_equal_type(residual, type_fit)
      scored_eq <- rank_ratio(ss_eq$ss_between, ss_eq$ss_within)
      order_eq <- order(
        -scored_eq$between_over_within,
        -ss_eq$ss_between,
        rownames(ss_eq$means)
      )
      rank_eq <- integer(nrow(ss_eq$means))
      rank_eq[order_eq] <- seq_along(order_eq)
      global_equal_parts[[level]] <- data.frame(
        time_point = level,
        cell_type = NA_character_,
        gene = rownames(ss_eq$means),
        ss_between = ss_eq$ss_between,
        ss_within = ss_eq$ss_within,
        between_over_within = scored_eq$between_over_within,
        between_over_total = scored_eq$between_over_total,
        top_cell_type = top_mean_group(ss_eq$means),
        passes_between_gt_within = ss_eq$ss_between > ss_eq$ss_within,
        rank = rank_eq,
        selected = rank_eq <= n_genes_global &
          is.finite(scored_eq$between_over_within),
        stringsAsFactors = FALSE
      )
    }

    if (do_type) {
      type_rows <- lapply(unique(type_fit), function(type_name) {
        contrast <- ifelse(type_fit == type_name, type_name, "rest")
        ss <- anova_ss(residual, contrast)
        scored <- rank_ratio(ss$ss_between, ss$ss_within)
        order_gene <- order(
          -scored$between_over_within,
          -ss$ss_between,
          rownames(residual)
        )
        rank <- integer(nrow(residual))
        rank[order_gene] <- seq_along(order_gene)
        data.frame(
          time_point = level,
          cell_type = type_name,
          gene = rownames(residual),
          ss_between = ss$ss_between,
          ss_within = ss$ss_within,
          between_over_within = scored$between_over_within,
          between_over_total = scored$between_over_total,
          top_cell_type = type_name,
          passes_between_gt_within = ss$ss_between > ss$ss_within,
          rank = rank,
          selected = rank <= n_genes_celltype &
            is.finite(scored$between_over_within),
          stringsAsFactors = FALSE
        )
      })
      type_parts[[level]] <- do.call(rbind, type_rows)
    }

    rm(sub, residual)
    gc()
  }

  list(
    global = if (do_global) {
      if (length(global_parts) == 0L) {
        empty_rank_table()
      } else {
        do.call(rbind, global_parts)
      }
    } else {
      NULL
    },
    global_equal_type = if (do_global) {
      if (length(global_equal_parts) == 0L) {
        empty_rank_table()
      } else {
        do.call(rbind, global_equal_parts)
      }
    } else {
      NULL
    },
    per_cell_type = if (do_type) {
      if (length(type_parts) == 0L) {
        empty_rank_table()
      } else {
        do.call(rbind, type_parts)
      }
    } else {
      NULL
    },
    model = if (length(model_parts) == 0L) {
      data.frame(
        time_point = character(),
        vst_flavor = character(),
        method = character(),
        latent_var = character(),
        n_cells_param = integer(),
        n_cells_fit = integer(),
        n_genes_residual = integer(),
        clip = character(),
        stringsAsFactors = FALSE
      )
    } else {
      do.call(rbind, model_parts)
    }
  )
}
