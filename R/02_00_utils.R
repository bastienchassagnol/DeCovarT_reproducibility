#' Default to the second argument when the first is empty.
#'
#' @param x Primary value.
#' @param y Fallback used when `x` is `NULL`, has length zero, or is
#'   all `NA`.
#' @return `x` when it is usable, otherwise `y`.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

#' Stop unless the OpenBLAS kernel gives correct factorisations.
#'
#' On this workstation (Ubuntu 22.04, OpenBLAS 0.3.20, AVX-512 CPU) the
#' default kernel returns wrong `chol()`, `solve()` and `eigen()` results
#' without an error; `OPENBLAS_CORETYPE=Haswell` is correct
#' (`DeCovarT/scripts/failed_reports/openblas_launchpad_comment.md`). The
#' variable is read when the library is loaded, so it must be set in the
#' shell that starts R. The check is the test of that report: the
#' Cholesky factor and the eigen-decomposition of the AR(1) matrix
#' `0.5^|i - j|`. Passing it does not exclude the rarer `eigen()`
#' failures that the report found under other kernels than Haswell.
#'
#' @param n Size of the test matrix. Default `64`.
#' @param tol Largest accepted reconstruction error. Default `1e-10`.
#' @return The reconstruction error, invisibly.
check_blas_kernel <- function(n = 64L, tol = 1e-10) {
  a <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
  r <- tryCatch(chol(a), error = function(e) NULL)
  eg <- tryCatch(eigen(a, symmetric = TRUE), error = function(e) NULL)
  err <- if (is.null(r) || is.null(eg)) {
    Inf
  } else {
    max(
      abs(crossprod(r) - a),
      abs(eg$vectors %*% (eg$values * t(eg$vectors)) - a)
    )
  }
  if (!is.finite(err) || err > tol) {
    stop(
      "chol() or eigen() is wrong under this BLAS kernel (error ",
      signif(err, 3),
      "; OPENBLAS_CORETYPE = '",
      Sys.getenv("OPENBLAS_CORETYPE"),
      "'). Start R with OPENBLAS_CORETYPE=Haswell.",
      call. = FALSE
    )
  }
  invisible(err)
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

#' Plot theme of the feature-selection figures.
#'
#' `theme_minimal()` with the grid, type and white background of the
#' book's other figures (`scripts/06_01_plot_early_wnt_compositions.R`).
#'
#' @param base_size Base font size. Default `11` (multi-panel figures).
#' @return A ggplot2 theme.
theme_slide <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e1d8"),
      plot.title = ggplot2::element_text(face = "bold", colour = "#0e4d3a"),
      plot.subtitle = ggplot2::element_text(colour = "#3d4a42"),
      axis.title = ggplot2::element_text(colour = "#1c2b24"),
      strip.text = ggplot2::element_text(face = "bold", colour = "#0e4d3a"),
      legend.position = "bottom",
      plot.title.position = "plot",
      panel.spacing = ggplot2::unit(1.2, "lines"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' Spread direct labels so that they do not overlap.
#'
#' Moves the positions `y` apart, keeping their order, so that two
#' neighbours are at least `min_gap` apart; used to place the labels of
#' several lines at their right end (data-visualisation rule #41).
#'
#' @param y Numeric positions (same scale as the axis).
#' @param min_gap Smallest distance between two labels.
#' @return The adjusted positions, in the input order.
spread_labels <- function(y, min_gap) {
  o <- order(y)
  z <- y[o]
  for (k in seq_along(z)[-1L]) {
    z[k] <- max(z[k], z[k - 1L] + min_gap)
  }
  # Re-centre so the block keeps the mean position of the input.
  z <- z - (mean(z) - mean(y[o]))
  out <- y
  out[o] <- z
  out
}

#' Format a count with a thousands separator.
#'
#' @param x Numbers.
#' @return Character vector, for example `"6,604"`.
fmt_n <- function(x) format(round(x), big.mark = ",", trim = TRUE)

#' Cell-type colours of the Suppinger annotation.
#'
#' @param path CSV with columns `celltypeannotation` and `colour`.
#' @return Named character vector of colours.
celltype_colours <- function(
  path = file.path("data", "dictionaries", "suppinger_celltype_colours.csv")
) {
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  stats::setNames(tab$colour, tab$celltypeannotation)
}

#' Save one figure as PDF (archive) and PNG (book).
#'
#' The PDF goes to `pdf_dir`, the PNG (`ragg`, white background) to
#' `png_dir`, both named `name`.
#'
#' @param plot A ggplot or patchwork object.
#' @param name File name without extension.
#' @param pdf_dir,png_dir Output folders, created when missing.
#' @param width,height Size in inches.
#' @param dpi PNG resolution. Default `300`.
#' @return The two paths, invisibly.
save_figure <- function(
  plot,
  name,
  pdf_dir,
  png_dir = file.path("docs", "images", "feature_selection"),
  width = 8,
  height = 6,
  dpi = 300
) {
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(pdf_dir, paste0(name, ".pdf"))
  png_path <- file.path(png_dir, paste0(name, ".png"))
  ggplot2::ggsave(
    pdf_path,
    plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    png_path,
    plot,
    width = width,
    height = height,
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white"
  )
  invisible(c(pdf = pdf_path, png = png_path))
}
