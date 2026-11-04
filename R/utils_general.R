#' Default to the second argument when the first is empty.
#'
#' @param x Primary value.
#' @param y Fallback used when `x` is `NULL`, has length zero, or is
#'   all `NA`.
#' @return `x` when it is usable, otherwise `y`.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

#' One character label per cell from a Seurat metadata column.
#'
#' `seu[[column]]` is a one-column data frame. Coercing that frame
#' with `as.character()` deparses every cell, so callers must take
#' the column vector before `unique()` or a type loop.
#'
#' @param seu A Seurat object (or anything whose `[[` returns the column).
#' @param column Metadata column name.
#' @return An unnamed character vector of length `ncol(seu)`.
seurat_metadata_chr <- function(seu, column) {
  values <- seu[[column]]
  if (is.data.frame(values)) {
    values <- values[[1L]]
  }
  as.character(values)
}

#' Pull a dense or sparse assay layer from a Seurat object.
#'
#' Tries the Seurat v5 `layer` argument first, then the v4 `slot`
#' argument, so the same call works on both object generations.
#'
#' @param seu A Seurat object.
#' @param assay Assay name, for example `"RNA"` or `"SCT"`.
#' @param layer Layer or slot name, for example `"counts"` or
#'   `"scale.data"`.
#' @return The requested matrix (typically genes by cells).
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

#' Pairwise Jaccard indices and membership patterns of named sets.
#'
#' Used to compare gene sets: the flags of the tail statements, or the
#' panels of the omnibus rank and of the classifiers.
#'
#' @param sets Named list of character vectors.
#' @return A list with `jaccard` (sets by sets matrix,
#'   `|A n B| / |A u B|`, `NA` when both are empty), `intersection`
#'   (sets by sets matrix of `|A n B|`), and `patterns`, an UpSet-style
#'   table with one logical column per set and `n`, the number of
#'   elements with exactly that membership pattern (largest first).
set_overlap <- function(sets) {
  if (!is.list(sets) || is.null(names(sets)) || any(names(sets) == "")) {
    stop("`sets` must be a named list of character vectors.")
  }
  sets <- lapply(sets, function(s) unique(as.character(s)))
  k <- length(sets)
  inter <- matrix(
    0L,
    k,
    k,
    dimnames = list(names(sets), names(sets))
  )
  jac <- matrix(NA_real_, k, k, dimnames = dimnames(inter))
  for (a in seq_len(k)) {
    for (b in seq_len(k)) {
      n_int <- length(intersect(sets[[a]], sets[[b]]))
      n_uni <- length(union(sets[[a]], sets[[b]]))
      inter[a, b] <- n_int
      jac[a, b] <- if (n_uni > 0L) n_int / n_uni else NA_real_
    }
  }
  universe <- unique(unlist(sets, use.names = FALSE))
  if (length(universe) == 0L) {
    patterns <- data.frame(n = integer())
  } else {
    member <- vapply(
      sets,
      function(s) universe %in% s,
      logical(length(universe))
    )
    member <- matrix(
      member,
      nrow = length(universe),
      dimnames = list(NULL, names(sets))
    )
    key <- apply(member, 1L, function(r) paste(as.integer(r), collapse = ""))
    counts <- table(key)
    first <- member[match(names(counts), key), , drop = FALSE]
    patterns <- data.frame(first, n = as.integer(counts), check.names = FALSE)
    patterns <- patterns[order(-patterns$n), , drop = FALSE]
    rownames(patterns) <- NULL
  }
  list(jaccard = jac, intersection = inter, patterns = patterns)
}
