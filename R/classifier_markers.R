#' Classifier-based marker selection with stability selection.
#'
#' Perspective of the book chapter "Feature selection" (section
#' "Comparison with classifier-based marker selection"). Three
#' classifiers of the cell type are fitted on subsamples of the cells of
#' one time slot:
#'
#' - `glmnet`: penalised multinomial logistic regression with a grouped
#'   lasso penalty (`type.multinomial = "grouped"`), so that a gene is in
#'   or out for every type at once;
#' - `splsda`: sparse PLS discriminant analysis (`mixOmics::splsda()`);
#' - `ranger`: random forest (`ranger::ranger()`), genes ranked by
#'   impurity importance.
#'
#' Each subsample keeps a fraction `frac` of the cells of **each** type,
#' so that rare types stay represented. A gene's selection frequency is
#' the share of subsamples in which the classifier selects it among its
#' top `n_top`; the stable set is the genes with frequency at least
#' `threshold` (Meinshausen and Buehlmann 2010). The cells come from one
#' specimen, so the resampling unit is the cell, and the frequencies
#' measure how the selection generalises to new cells of this specimen,
#' not to new specimens. The purpose is to check that the omnibus panel
#' is biologically meaningful, robust and able to discriminate the types,
#' not to replace it.
#'
#' `glmnet`, `mixOmics` and `ranger` are not in `renv.lock`; each is
#' checked with `requireNamespace()`. Requires `R/utils_general.R`.
#' Written, not run.

#' Selection frequencies of three classifiers under subsampling.
#'
#' @param counts Genes by cells raw counts, gene names as row names.
#' @param cell_type Cell-type label per column.
#' @param size_factors Size factor per cell. Features are
#'   `log1p(counts / size_factors)`, cells in rows.
#' @param method Any of `"glmnet"`, `"splsda"`, `"ranger"`.
#' @param n_subsample Number of subsamples. Default `100`.
#' @param frac Fraction of the cells of each type per subsample.
#'   Default `0.5`.
#' @param n_top Genes selected per fit. Default `50`.
#' @param threshold Selection frequency of the stable set. Default `0.6`.
#' @param seed Integer seed, or `NULL` for the current RNG state.
#' @param num_threads Threads for `ranger`. Default `2`.
#' @return A data frame with one row per gene and method: `method`,
#'   `gene`, `frequency`, `stable`.
classifier_markers <- function(
  counts,
  cell_type,
  size_factors,
  method = c("glmnet", "splsda", "ranger"),
  n_subsample = 100L,
  frac = 0.5,
  n_top = 50L,
  threshold = 0.6,
  seed = NULL,
  num_threads = 2L
) {
  method <- match.arg(method, several.ok = TRUE)
  pkg <- c(glmnet = "glmnet", splsda = "mixOmics", ranger = "ranger")[method]
  missing_pkg <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkg) > 0L) {
    stop(
      "Packages not installed (not in renv.lock): ",
      paste(unique(missing_pkg), collapse = ", "),
      ". Install them, or drop the matching `method`."
    )
  }
  cell_type <- droplevels(factor(cell_type))
  if (length(cell_type) != ncol(counts)) {
    stop("`cell_type` must have one label per column of `counts`.")
  }
  if (nlevels(cell_type) < 2L) {
    stop("Need at least two cell types.")
  }
  if (length(size_factors) != ncol(counts) || any(size_factors <= 0)) {
    stop("`size_factors` must be positive, one per column of `counts`.")
  }
  if (is.null(rownames(counts))) {
    stop("`counts` needs gene names as row names.")
  }
  if (frac <= 0 || frac >= 1) {
    stop("`frac` must lie strictly between 0 and 1.")
  }
  genes <- rownames(counts)
  n_top <- min(as.integer(n_top), length(genes))
  x <- t(log1p(sweep(as.matrix(counts), 2L, size_factors, "/")))
  colnames(x) <- genes
  by_type <- split(seq_along(cell_type), cell_type)
  n_types <- nlevels(cell_type)

  select_one <- function(m, idx) {
    xs <- x[idx, , drop = FALSE]
    ys <- droplevels(cell_type[idx])
    keep <- apply(xs, 2L, stats::var) > 0
    xs <- xs[, keep, drop = FALSE]
    switch(
      m,
      glmnet = {
        fit <- glmnet::glmnet(
          xs,
          ys,
          family = "multinomial",
          type.multinomial = "grouped"
        )
        # Largest model on the path with at most n_top genes.
        n_sel <- fit$df
        step <- max(which(n_sel <= n_top))
        beta <- do.call(
          cbind,
          lapply(fit$beta, function(b) as.matrix(b[, step, drop = FALSE]))
        )
        colnames(xs)[rowSums(abs(beta)) > 0]
      },
      splsda = {
        ncomp <- n_types - 1L
        keep_x <- rep(max(1L, ceiling(n_top / ncomp)), ncomp)
        fit <- mixOmics::splsda(xs, ys, ncomp = ncomp, keepX = keep_x)
        load <- fit$loadings$X
        rownames(load)[rowSums(abs(load)) > 0]
      },
      ranger = {
        fit <- ranger::ranger(
          x = xs,
          y = ys,
          num.trees = 500L,
          importance = "impurity",
          num.threads = num_threads
        )
        imp <- sort(fit$variable.importance, decreasing = TRUE)
        names(imp)[seq_len(min(n_top, length(imp)))]
      }
    )
  }

  run <- function() {
    hits <- lapply(method, function(m) {
      stats::setNames(integer(length(genes)), genes)
    })
    names(hits) <- method
    for (b in seq_len(n_subsample)) {
      idx <- unlist(lapply(by_type, function(i) {
        i[sample.int(length(i), max(2L, floor(frac * length(i))))]
      }))
      for (m in method) {
        sel <- select_one(m, idx)
        hits[[m]][sel] <- hits[[m]][sel] + 1L
      }
    }
    hits
  }
  hits <- if (is.null(seed)) run() else withr::with_seed(seed, run())
  out <- do.call(
    rbind,
    lapply(method, function(m) {
      data.frame(
        method = m,
        gene = genes,
        frequency = hits[[m]] / n_subsample,
        stable = hits[[m]] / n_subsample >= threshold,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  attr(out, "settings") <- list(
    n_subsample = n_subsample,
    frac = frac,
    n_top = n_top,
    threshold = threshold,
    seed = seed
  )
  out
}
