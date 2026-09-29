#' Default to the second argument when the first is empty.
#'
#' @param x Primary value.
#' @param y Fallback used when `x` is `NULL`, has length zero, or is
#'   all `NA`.
#' @return `x` when it is usable, otherwise `y`.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
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
