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
