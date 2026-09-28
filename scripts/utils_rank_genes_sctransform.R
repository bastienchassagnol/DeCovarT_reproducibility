# Rank genes from SCTransform v2 Pearson residuals.
#
# The fit is an intercept plus a library-size offset. Cell type and
# cell line are not covariates. Pass `split_by` (here, time) to fit
# once within each level, on all retained types together. Do not
# subset by cell type before the fit: that would absorb the
# between-type mean the ranking is meant to keep.
#
# Caller responsibility: if the count matrix has already been
# restricted to a gene allow-list, set metadata column `log_umi` to
# log10 of the full-transcriptome library size. sctransform reuses
# that column as the v2 offset and otherwise recomputes depth from
# the genes still in the assay.

# ==========================================================================
# Helpers ----
# ==========================================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

assay_matrix <- function(seu, assay, layer) {
  got <- tryCatch(
    SeuratObject::GetAssayData(
      object = seu,
      assay = assay,
      layer = layer
    ),
    error = function(e) NULL
  )
  if (is.null(got)) {
    got <- SeuratObject::GetAssayData(
      object = seu,
      assay = assay,
      slot = layer
    )
  }
  got
}

# Between / within sums of squares on a genes x cells residual matrix.
# `groups` aligns with columns. This is the ANOVA split in eq-fs-ss.
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

# Equal type weight. Each of the J types has weight 1/J, so a large
# type cannot dominate the rank. `ss_between` and `ss_within` here are
# mean squares: the variance of the J type means, and the mean of the
# J within-type variances.
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

rank_from_ss <- function(ss, time_level, n_genes, cell_type = NA_character_) {
  scored <- rank_ratio(ss$ss_between, ss$ss_within)
  order_gene <- order(
    -scored$between_over_within,
    -ss$ss_between,
    rownames(ss$means)
  )
  rank <- integer(nrow(ss$means))
  rank[order_gene] <- seq_along(order_gene)
  data.frame(
    time_point = time_level,
    cell_type = cell_type,
    gene = rownames(ss$means),
    ss_between = ss$ss_between,
    ss_within = ss$ss_within,
    between_over_within = scored$between_over_within,
    between_over_total = scored$between_over_total,
    top_cell_type = if (length(cell_type) == 1L && is.na(cell_type)) {
      top_mean_group(ss$means)
    } else {
      cell_type
    },
    passes_between_gt_within = ss$ss_between > ss$ss_within,
    rank = rank,
    selected = rank <= n_genes & is.finite(scored$between_over_within),
    stringsAsFactors = FALSE
  )
}

rank_ratio <- function(ss_between, ss_within) {
  total <- ss_between + ss_within
  keep <- total > 0
  ratio <- rep(NA_real_, length(ss_between))
  share <- rep(NA_real_, length(ss_between))
  ratio[keep] <- ss_between[keep] / pmax(
    ss_within[keep],
    .Machine$double.eps
  )
  share[keep] <- ss_between[keep] / total[keep]
  list(between_over_within = ratio, between_over_total = share)
}

top_mean_group <- function(means) {
  idx <- max.col(means, ties.method = "first")
  colnames(means)[idx]
}

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

# Map dictionary symbols onto mouse symbols in `universe`.
# Offline: case fold, then org.Mm.eg.db SYMBOL / ALIAS. No biomaRt.
map_marker_symbols <- function(markers, universe) {
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop(
      "Package org.Mm.eg.db is required for offline mouse symbol ",
      "mapping. Install it with BiocManager::install(\"org.Mm.eg.db\")."
    )
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package AnnotationDbi is required for offline mouse symbol mapping.")
  }
  markers <- unique(as.character(markers))
  universe <- unique(as.character(universe))
  upper_hits <- split(universe, toupper(universe))
  symbols <- AnnotationDbi::keys(
    org.Mm.eg.db::org.Mm.eg.db,
    keytype = "SYMBOL"
  )
  alias_tbl <- suppressMessages(
    AnnotationDbi::select(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = symbols,
      columns = "ALIAS",
      keytype = "SYMBOL"
    )
  )
  alias_tbl <- alias_tbl[
    !is.na(alias_tbl$ALIAS) & nzchar(alias_tbl$ALIAS),
    ,
    drop = FALSE
  ]
  alias_upper <- split(alias_tbl$SYMBOL, toupper(alias_tbl$ALIAS))
  symbol_upper <- split(symbols, toupper(symbols))

  resolve_one <- function(marker) {
    key <- toupper(marker)
    in_matrix <- upper_hits[[key]] %||% character()
    if (length(unique(in_matrix)) == 1L) {
      return(list(mouse = unique(in_matrix), method = "casefold_universe"))
    }
    if (length(unique(in_matrix)) > 1L) {
      return(list(mouse = NA_character_, method = "ambiguous_universe"))
    }
    db_hit <- unique(symbol_upper[[key]] %||% character())
    in_universe <- intersect(db_hit, universe)
    if (length(in_universe) == 1L) {
      return(list(mouse = in_universe, method = "org.Mm.eg.db_symbol"))
    }
    if (length(db_hit) == 1L && length(in_universe) == 0L) {
      return(list(mouse = db_hit, method = "org.Mm.eg.db_symbol_absent"))
    }
    alias_hit <- unique(alias_upper[[key]] %||% character())
    alias_in <- intersect(alias_hit, universe)
    if (length(alias_in) == 1L) {
      return(list(mouse = alias_in, method = "org.Mm.eg.db_alias"))
    }
    if (length(alias_hit) == 1L && length(alias_in) == 0L) {
      return(list(
        mouse = alias_hit,
        method = "org.Mm.eg.db_alias_absent"
      ))
    }
    if (length(alias_in) > 1L || length(db_hit) > 1L) {
      return(list(mouse = NA_character_, method = "ambiguous_alias"))
    }
    list(mouse = NA_character_, method = "unmapped")
  }

  resolved <- lapply(markers, resolve_one)
  data.frame(
    gene_marker = markers,
    mouse_symbol = vapply(resolved, function(x) x$mouse, character(1)),
    map_method = vapply(resolved, function(x) x$method, character(1)),
    stringsAsFactors = FALSE
  )
}

# ==========================================================================
# Ranker ----
# ==========================================================================

#' Rank genes by between/within Pearson-residual sums of squares.
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
#' @return A list with `global`, `per_cell_type`, and `model`.
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
      global_parts[[level]] <- rank_from_ss(
        ss = anova_ss(residual, type_fit),
        time_level = level,
        n_genes = n_genes_global
      )
      global_equal_parts[[level]] <- rank_from_ss(
        ss = anova_ms_equal_type(residual, type_fit),
        time_level = level,
        n_genes = n_genes_global
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

  bind <- function(parts) {
    if (length(parts) == 0L) {
      return(empty_rank_table())
    }
    do.call(rbind, parts)
  }

  list(
    global = if (do_global) bind(global_parts) else NULL,
    global_equal_type = if (do_global) bind(global_equal_parts) else NULL,
    per_cell_type = if (do_type) bind(type_parts) else NULL,
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
