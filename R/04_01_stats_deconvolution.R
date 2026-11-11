#' Deconvolution of the Suppinger bulk samples: DeCovarT flavours and
#' comparators, on one gene panel and one single-cell reference.
#'
#' Every function takes the bulk on raw counts (genes by samples) and a
#' reference of one time point: either the per-type moments `mu` (genes by
#' types, mean raw counts per cell) and `Sigma` (genes by genes by types),
#' or the raw counts of the reference cells. Every wrapper returns one row
#' per sample and cell type: `sample_id`, `cell_type`, `proportion`,
#' `lower`, `upper` (`NA` when the method gives no interval).
#'
#' Requires `R/02_00_utils.R` and, for the DeMixSC weights,
#' `R/02_03_stats_tail_bounds.R` (`plugin_composition_demixsc()`).

# ==========================================================================
# Inputs: design, observed proportions, platform ratio ----
# ==========================================================================

#' Design of the Suppinger bulk samples, from their column names.
#'
#' Names follow `<time>_<CTL|Early>_rep<k>`: `CTL` is the standard Wnt
#' pulse, `Early` the early pulse.
#'
#' @param sample_names Column names of the HTSeq count matrix.
#' @return A data frame with `sample_id`, `time`, `arm` (`"standard"` or
#'   `"early"`) and `replicate`.
bulk_sample_design <- function(sample_names) {
  sample_names <- as.character(sample_names)
  pattern <- "^(\\d+h)_(CTL|Early)_rep(\\d+)$"
  if (!all(grepl(pattern, sample_names))) {
    stop(
      "Unexpected bulk column names: ",
      paste(sample_names[!grepl(pattern, sample_names)], collapse = ", "),
      "."
    )
  }
  data.frame(
    sample_id = sample_names,
    time = sub(pattern, "\\1", sample_names),
    arm = ifelse(
      sub(pattern, "\\2", sample_names) == "CTL",
      "standard",
      "early"
    ),
    replicate = as.integer(sub(pattern, "\\3", sample_names)),
    stringsAsFactors = FALSE
  )
}

#' Observed proportions: fractions of the scored cells at each time.
#'
#' The proportions are those of the single-cell annotation, among the
#' scored types only (`scored == TRUE` in `cell_counts.csv`), with the
#' cells of the cell lines pooled. They are not flow-cytometry
#' measurements, and they come from the specimen that also gives the
#' reference.
#'
#' @param path Path to `output/feature_selection/tier_01/cell_counts.csv`.
#' @return A data frame with `time`, `cell_type`, `cells`, `observed`.
observed_scored_fractions <- function(path) {
  counts <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  counts <- counts[counts$scored, , drop = FALSE]
  total <- tapply(counts$cells, counts$time_point, sum)
  data.frame(
    time = counts$time_point,
    cell_type = counts$cell_type,
    cells = counts$cells,
    observed = counts$cells / as.numeric(total[counts$time_point]),
    stringsAsFactors = FALSE
  )
}

#' Cross-fitted gene-wise platform ratio `b_g` of one time.
#'
#' Reads `platform_ratio.csv` of tier (iii). The ratio used at one time
#' was estimated on the bulk replicates of the two other times only, and
#' its median log is zero.
#'
#' @param path Path to `output/feature_selection/tier_03/platform_ratio.csv`.
#' @param time_level Time whose bulk samples the ratio is applied to.
#' @param genes Genes to return, in this order.
#' @param method `"gene"` (default, a free ratio per gene) or `"length"`.
#' @return Named numeric vector `b_g` (linear scale).
read_platform_ratio <- function(
  path,
  time_level,
  genes,
  method = c("gene", "length")
) {
  method <- match.arg(method)
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  tab <- tab[tab$method == method & tab$held_out_time == time_level, ]
  missing <- setdiff(genes, tab$gene)
  if (length(missing) > 0L) {
    stop(
      length(missing),
      " genes have no platform ratio at ",
      time_level,
      " (first: ",
      missing[[1L]],
      ")."
    )
  }
  b <- 2^tab$log2_b[match(genes, tab$gene)]
  if (any(!is.finite(b)) || any(b <= 0)) {
    stop("The platform ratio must be positive and finite.")
  }
  stats::setNames(b, genes)
}

# ==========================================================================
# DeCovarT: single-cell model, four covariance inputs ----
# ==========================================================================

#' Denoised covariance plus the bulk Poisson term, per type.
#'
#' The bulk law with a denoised within-type covariance is
#' `Cov y = phi * {S_bio(p) + diag(mu p) / a}`, with `a = phi / c` the
#' bulk depth per cell (DeCovarT vignette 3.2, eq-bulk-plugin). The
#' diagonal term is linear in `p`, so it can be carried by the per-type
#' matrices: `Sigma_j = Sigma_bio_j + diag(mu_j) / a`.
#'
#' When the bulk was divided by a gene-wise platform ratio `b`, the mean
#' is unchanged and the Poisson variance of gene `g` is divided by `b_g`:
#' `Sigma_j = Sigma_bio_j + diag(mu_j / b) / a`. This second form is
#' derived here from the first; it is not in the vignette.
#'
#' @param sigma_bio Denoised covariances, genes by genes by types.
#' @param mu Mean raw counts per cell, genes by types.
#' @param a Positive bulk depth per cell.
#' @param b Optional gene-wise platform ratio the bulk was divided by.
#' @return An array like `sigma_bio`.
add_bulk_poisson_term <- function(sigma_bio, mu, a, b = NULL) {
  n_genes <- nrow(mu)
  if (!is.finite(a) || a <= 0) {
    stop("`a` must be one positive number.")
  }
  if (is.null(b)) {
    b <- rep(1, n_genes)
  }
  out <- sigma_bio
  on_diagonal <- (seq_len(n_genes) - 1L) * n_genes + seq_len(n_genes)
  for (j in seq_len(ncol(mu))) {
    idx <- (j - 1) * n_genes^2 + on_diagonal
    out[idx] <- out[idx] + mu[, j] / (b * a)
  }
  out
}

#' One DeCovarT single-cell fit, with its Wald interval and diagnostics.
#'
#' Profile fit of [DeCovarT::deconvolute_single_cell()] (both technical
#' scales profiled out, cellular estimand). The interval is
#' `p +/- z * se` from the delta-method covariance of the fit, clipped
#' to `[0, 1]`; it treats the reference moments as known. A fit that
#' stops with an error returns `NA` proportions and the message.
#'
#' @param y Bulk vector, genes named as the rows of `mu`.
#' @param mu Mean signature, genes by types.
#' @param sigma Covariances, genes by genes by types.
#' @param optimiser `"newton"` or `"lbfgsb"`.
#' @param initial_p Optional start on the simplex.
#' @param level Confidence level. Default `0.95`.
#' @return A list with `estimates` (one row per type), `diagnostics` (one
#'   row) and `p` (named vector, or `NULL` after an error).
fit_decovart_once <- function(
  y,
  mu,
  sigma,
  optimiser,
  initial_p = NULL,
  level = 0.95
) {
  types <- colnames(mu)
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    DeCovarT::deconvolute_single_cell(
      y,
      mu,
      sigma,
      optimiser = optimiser,
      initial_p = initial_p
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(fit, "error")) {
    return(list(
      estimates = data.frame(
        cell_type = types,
        proportion = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        stringsAsFactors = FALSE
      ),
      diagnostics = data.frame(
        mean_scale = NA_real_,
        variance_scale = NA_real_,
        platform_scale = NA_real_,
        n_cells_hat = NA_real_,
        loglik = NA_real_,
        convergence = NA_integer_,
        iterations = NA_integer_,
        elapsed_sec = elapsed,
        message = conditionMessage(fit),
        stringsAsFactors = FALSE
      ),
      p = NULL
    ))
  }
  p <- fit$p[types]
  se <- rep(NA_real_, length(types))
  if (is.matrix(fit$vcov) && all(is.finite(diag(fit$vcov)))) {
    se <- sqrt(pmax(diag(fit$vcov)[types], 0))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  list(
    estimates = data.frame(
      cell_type = types,
      proportion = as.numeric(p),
      lower = pmax(as.numeric(p) - z * se, 0),
      upper = pmin(as.numeric(p) + z * se, 1),
      stringsAsFactors = FALSE
    ),
    diagnostics = data.frame(
      mean_scale = fit$mean_scale,
      variance_scale = fit$variance_scale,
      platform_scale = fit$platform_scale,
      n_cells_hat = fit$n_cells,
      loglik = fit$loglik,
      convergence = as.integer(fit$convergence),
      iterations = as.integer(sum(fit$iterations)),
      elapsed_sec = elapsed,
      message = as.character(fit$message %||% ""),
      stringsAsFactors = FALSE
    ),
    p = p
  )
}

#' The four covariance inputs of DeCovarT for one bulk sample.
#'
#' Fits, in this order:
#' - `raw`: the regularised raw covariance, which assumes a bulk depth
#'   per cell `a = 1`;
#' - `denoised_bio`: the denoised covariance alone, with no bulk
#'   sampling term;
#' - `denoised_plugin`: the denoised covariance plus the bulk Poisson
#'   term ([add_bulk_poisson_term()]) at `a = phi / c` of the `raw` fit;
#' - `denoised_fixed`: the same form, with `a` replaced by `phi / c` of
#'   the previous fit until `|log(a_new / a)| < tol` or `max_passes`
#'   fits. Its first pass is the `denoised_plugin` fit. The limit is a
#'   fixed point, not the joint maximum-likelihood estimate.
#'
#' With `b`, the bulk is divided by the gene-wise platform ratio before
#' every fit; the reference moments are unchanged.
#'
#' @param y Bulk vector on raw counts, named by gene.
#' @param mu Mean signature, genes by types.
#' @param sigma_raw,sigma_bio Covariances, genes by genes by types.
#' @param optimiser `"newton"` or `"lbfgsb"`.
#' @param b Optional named gene-wise platform ratio.
#' @param max_passes Largest number of fits of the fixed point. Default 5.
#' @param tol Tolerance on `|log(a_new / a)|`. Default `0.01`.
#' @return A list with `estimates` (`covariance`, `cell_type`,
#'   `proportion`, `lower`, `upper`) and `diagnostics` (one row per
#'   covariance input, with `a_used`, `passes` and `fixed_point_reached`).
fit_decovart_chain <- function(
  y,
  mu,
  sigma_raw,
  sigma_bio,
  optimiser,
  b = NULL,
  max_passes = 5L,
  tol = 0.01
) {
  genes <- rownames(mu)
  if (
    is.null(genes) ||
      !identical(dimnames(sigma_raw)[[1L]], genes) ||
      !identical(dimnames(sigma_bio)[[1L]], genes) ||
      !all(genes %in% names(y))
  ) {
    stop("`y`, `mu` and the covariances must share their gene names.")
  }
  y <- as.numeric(y[genes])
  if (!is.null(b)) {
    b <- as.numeric(b[genes])
    y <- y / b
  }
  # One entry per covariance input: the fit, the depth used in the
  # Poisson term, the number of fits and whether the fixed point was
  # reached.
  entry <- function(fit, a_used, passes, reached) {
    list(fit = fit, a_used = a_used, passes = passes, reached = reached)
  }
  skipped <- list(
    estimates = data.frame(
      cell_type = colnames(mu),
      proportion = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      stringsAsFactors = FALSE
    ),
    diagnostics = data.frame(
      mean_scale = NA_real_,
      variance_scale = NA_real_,
      platform_scale = NA_real_,
      n_cells_hat = NA_real_,
      loglik = NA_real_,
      convergence = NA_integer_,
      iterations = NA_integer_,
      elapsed_sec = 0,
      message = "No bulk depth per cell from the raw fit.",
      stringsAsFactors = FALSE
    ),
    p = NULL
  )

  raw <- fit_decovart_once(y, mu, sigma_raw, optimiser)
  bio <- fit_decovart_once(y, mu, sigma_bio, optimiser)
  fits <- list(
    raw = entry(raw, NA_real_, 1L, NA),
    denoised_bio = entry(bio, NA_real_, 1L, NA)
  )

  a <- raw$diagnostics$platform_scale
  if (!is.finite(a) || a <= 0) {
    fits$denoised_plugin <- entry(skipped, NA_real_, 0L, NA)
    fits$denoised_fixed <- entry(skipped, NA_real_, 0L, NA)
  } else {
    current <- fit_decovart_once(
      y,
      mu,
      add_bulk_poisson_term(sigma_bio, mu, a, b),
      optimiser
    )
    fits$denoised_plugin <- entry(current, a, 1L, NA)
    passes <- 1L
    reached <- FALSE
    repeat {
      a_new <- current$diagnostics$platform_scale
      if (!is.finite(a_new) || a_new <= 0) {
        break
      }
      if (abs(log(a_new / a)) < tol) {
        reached <- TRUE
        break
      }
      if (passes >= max_passes) {
        break
      }
      a <- a_new
      current <- fit_decovart_once(
        y,
        mu,
        add_bulk_poisson_term(sigma_bio, mu, a, b),
        optimiser,
        initial_p = current$p
      )
      passes <- passes + 1L
    }
    fits$denoised_fixed <- entry(current, a, passes, reached)
  }
  estimates <- lapply(names(fits), function(covariance) {
    cbind(
      covariance = covariance,
      fits[[covariance]]$fit$estimates,
      stringsAsFactors = FALSE
    )
  })
  diagnostics <- lapply(names(fits), function(covariance) {
    cbind(
      covariance = covariance,
      a_used = fits[[covariance]]$a_used,
      passes = fits[[covariance]]$passes,
      fixed_point_reached = fits[[covariance]]$reached,
      fits[[covariance]]$fit$diagnostics,
      stringsAsFactors = FALSE
    )
  })
  out <- list(
    estimates = do.call(rbind, estimates),
    diagnostics = do.call(rbind, diagnostics)
  )
  rownames(out$estimates) <- NULL
  rownames(out$diagnostics) <- NULL
  out
}

# ==========================================================================
# Comparators ----
# ==========================================================================

#' Call a function; return its value or its error, and its warnings.
#'
#' Notices raised when a package is loaded (replaced imports) and the
#' command-line parsing of `snowfall` are not kept.
#'
#' @param fun A function of no argument.
#' @return A list with `value` (the result, or the error condition) and
#'   `warnings` (distinct warning messages, in order of appearance).
run_with_conditions <- function(fun) {
  seen <- new.env(parent = emptyenv())
  seen$warnings <- character()
  value <- withCallingHandlers(
    tryCatch(fun(), error = function(e) e),
    warning = function(w) {
      seen$warnings <- c(seen$warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  noise <- paste(
    "^replacing previous import",
    "^remplacement de l'importation",
    "^Unknown option on commandline",
    sep = "|"
  )
  kept <- unique(seen$warnings)
  list(value = value, warnings = kept[!grepl(noise, kept)])
}

#' Long table of proportions from a samples by types matrix.
#'
#' @param proportions Matrix, samples in rows, cell types in columns.
#' @param lower,upper Optional matrices of the same shape.
#' @return A data frame with `sample_id`, `cell_type`, `proportion`,
#'   `lower`, `upper`.
proportions_long <- function(proportions, lower = NULL, upper = NULL) {
  proportions <- as.matrix(proportions)
  grid <- expand.grid(
    sample_id = rownames(proportions),
    cell_type = colnames(proportions),
    stringsAsFactors = FALSE
  )
  pick <- function(m) {
    if (is.null(m)) {
      return(rep(NA_real_, nrow(grid)))
    }
    as.matrix(m)[cbind(grid$sample_id, grid$cell_type)]
  }
  grid$proportion <- pick(proportions)
  grid$lower <- pick(lower)
  grid$upper <- pick(upper)
  grid
}

#' DWLS: dampened weighted least squares on a supplied signature.
#'
#' One call of `DWLS::solveDampenedWLS()` per sample, which is what
#' `omnideconv::deconvolute_dwls(dwls_submethod = "DampenedWLS")` does.
#' The signature is the mean raw count per cell of each type on the
#' panel; the MAST signature step of DWLS is not used. The solver draws
#' gene subsets with its own fixed seeds; the RNG state of the caller is
#' restored.
#'
#' @param bulk Bulk raw counts, genes by samples.
#' @param mu Mean signature, genes by types.
#' @return Long table, see [proportions_long()].
deconvolute_dwls <- function(bulk, mu) {
  genes <- rownames(mu)
  bulk <- as.matrix(bulk)[genes, , drop = FALSE]
  est <- withr::with_preserve_seed(
    vapply(
      colnames(bulk),
      function(s) {
        sol <- DWLS::solveDampenedWLS(S = mu, B = bulk[, s])
        as.numeric(sol[colnames(mu)])
      },
      numeric(ncol(mu))
    )
  )
  rownames(est) <- colnames(mu)
  proportions_long(t(est))
}

#' SCDC with one single-cell reference.
#'
#' `SCDC::SCDC_prop_ONE()`: every reference cell is treated as coming
#' from one subject, so the cell line is ignored. The basis is built by
#' SCDC from the raw counts of the panel genes; the bulk is converted to
#' counts per million of the panel inside SCDC. Settings are the
#' defaults of upstream SCDC (`iter.max = 2000`, `nu = 1e-10`,
#' `epsilon = 0.001`, basis weighted by the cross-cell variance).
#'
#' @param bulk Bulk raw counts, genes by samples (at least two samples).
#' @param counts Reference raw counts, cells by genes.
#' @param cell_type Cell-type label of each reference cell.
#' @return Long table, see [proportions_long()].
deconvolute_scdc <- function(bulk, counts, cell_type) {
  genes <- colnames(counts)
  bulk <- as.matrix(bulk)[genes, , drop = FALSE]
  types <- sort(unique(cell_type))
  sc_pheno <- data.frame(
    cell_type = cell_type,
    subject = "reference",
    row.names = rownames(counts),
    stringsAsFactors = FALSE
  )
  sc_eset <- Biobase::ExpressionSet(
    assayData = t(counts),
    phenoData = Biobase::AnnotatedDataFrame(sc_pheno)
  )
  bulk_eset <- Biobase::ExpressionSet(
    assayData = bulk,
    phenoData = Biobase::AnnotatedDataFrame(
      data.frame(sample_id = colnames(bulk), row.names = colnames(bulk))
    )
  )
  fit <- SCDC::SCDC_prop_ONE(
    bulk.eset = bulk_eset,
    sc.eset = sc_eset,
    ct.varname = "cell_type",
    sample = "subject",
    ct.sub = types,
    iter.max = 2000,
    nu = 1e-10,
    epsilon = 0.001,
    weight.basis = TRUE
  )
  proportions_long(fit$prop.est.mvw[colnames(bulk), types, drop = FALSE])
}

#' BayesPrism on raw counts.
#'
#' `new.prism()` on the cell-level count matrix with cell states equal
#' to cell types and no malignant type, then `run.prism()` with the
#' second Gibbs pass on the updated reference, then the final cell-type
#' fractions.
#'
#' The outlier filter of `new.prism()` is switched off (`outlier.cut = 1`,
#' `outlier.fraction = 1`). Its default drops a gene that exceeds 1% of
#' the mixture counts in more than 10% of the samples, a rule written for
#' a whole transcriptome: on a panel of 17 genes it removed 13 of them.
#' Every method therefore sees the same panel.
#'
#' @param bulk Bulk raw counts, genes by samples.
#' @param counts Reference raw counts, cells by genes.
#' @param cell_type Cell-type label of each reference cell.
#' @param n_cores Cores passed to `run.prism()`. Default `2`.
#' @param seed Seed of the Gibbs sampler.
#' @return Long table, see [proportions_long()].
deconvolute_bayesprism <- function(
  bulk,
  counts,
  cell_type,
  n_cores = 2L,
  seed = 20261110L
) {
  genes <- colnames(counts)
  mixture <- t(as.matrix(bulk)[genes, , drop = FALSE])
  prism <- BayesPrism::new.prism(
    reference = counts,
    mixture = mixture,
    input.type = "count.matrix",
    cell.type.labels = cell_type,
    cell.state.labels = cell_type,
    key = NULL,
    outlier.cut = 1,
    outlier.fraction = 1
  )
  fit <- BayesPrism::run.prism(
    prism = prism,
    n.cores = n_cores,
    update.gibbs = TRUE,
    gibbs.control = list(seed = seed)
  )
  theta <- BayesPrism::get.fraction(
    bp = fit,
    which.theta = "final",
    state.or.type = "type"
  )
  types <- sort(unique(cell_type))
  proportions_long(theta[rownames(mixture), types, drop = FALSE])
}

#' RNA-Sieve through its Python implementation.
#'
#' Writes the bulk and one reference file per cell type, runs
#' `scripts/python/run_rnasieve.py` with `uv run --script`, and reads the
#' estimates. The intervals of `rnasieve` are not computed (the runner
#' has an `--intervals` switch), so `lower` and `upper` are `NA`. The
#' preprocessing
#' of `rnasieve` keeps its defaults for the cell trimming (2% per type)
#' and the sparsity filter (0.2). Its library-size normalisation is
#' switched off, because the column sum of a gene panel is not the
#' library of the cell.
#'
#' @param bulk Bulk raw counts, genes by samples.
#' @param counts Reference raw counts, cells by genes.
#' @param cell_type Cell-type label of each reference cell.
#' @param work_dir Folder for the exchange files (created, kept).
#' @param script Path to the Python runner.
#' @return Long table, see [proportions_long()].
deconvolute_rnasieve <- function(
  bulk,
  counts,
  cell_type,
  work_dir,
  script = file.path("scripts", "python", "run_rnasieve.py")
) {
  genes <- colnames(counts)
  types <- sort(unique(cell_type))
  in_dir <- file.path(work_dir, "input")
  out_dir <- file.path(work_dir, "output")
  unlink(c(in_dir, out_dir), recursive = TRUE)
  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_table <- function(m, path) {
    data.table::fwrite(
      data.table::data.table(gene = rownames(m), m),
      path
    )
  }
  bulk <- as.matrix(bulk)[genes, , drop = FALSE]
  write_table(bulk, file.path(in_dir, "bulk.csv"))
  utils::write.csv(
    data.frame(celltype = types),
    file.path(in_dir, "celltypes.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  for (k in seq_along(types)) {
    ref <- t(counts[cell_type == types[[k]], , drop = FALSE])
    colnames(ref) <- paste0("cell_", seq_len(ncol(ref)))
    write_table(ref, file.path(in_dir, paste0("reference_", k, ".csv")))
  }
  status <- system2(
    "uv",
    c("run", "--script", shQuote(script), shQuote(in_dir), shQuote(out_dir))
  )
  if (status != 0L) {
    stop("RNA-Sieve runner failed with status ", status, ".")
  }
  read_matrix <- function(name) {
    tab <- utils::read.csv(
      file.path(out_dir, name),
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fileEncoding = "UTF-8"
    )
    m <- as.matrix(tab[, -1L, drop = FALSE])
    rownames(m) <- tab[[1L]]
    t(m)[colnames(bulk), types, drop = FALSE]
  }
  proportions_long(
    read_matrix("estimates.csv"),
    read_matrix("lower.csv"),
    read_matrix("upper.csv")
  )
}

#' DeMixSC weights without the benchmark adjustment.
#'
#' The iterated weighted non-negative least squares of DeMixSC
#' ([plugin_composition_demixsc()]) on the mean signature. DeMixSC
#' calibrates a gene-wise platform adjustment on bulk and pseudo-bulk of
#' the same samples; no such pair exists for Suppinger, so that step is
#' absent. With `b`, the bulk is divided by the cross-fitted gene-wise
#' platform ratio instead.
#'
#' @param bulk Bulk raw counts, genes by samples.
#' @param mu Mean signature, genes by types.
#' @param b Optional named gene-wise platform ratio.
#' @param c_reg Constant of the DeMixSC weights. Default `2`.
#' @return Long table, see [proportions_long()].
deconvolute_demixsc_wnnls <- function(bulk, mu, b = NULL, c_reg = 2) {
  genes <- rownames(mu)
  bulk <- as.matrix(bulk)[genes, , drop = FALSE]
  if (!is.null(b)) {
    bulk <- bulk / as.numeric(b[genes])
  }
  fit <- plugin_composition_demixsc(bulk, mu, c_reg = c_reg)
  data.frame(
    sample_id = fit$composition$sample,
    cell_type = fit$composition$cell_type,
    proportion = fit$composition$p_hat,
    lower = NA_real_,
    upper = NA_real_,
    stringsAsFactors = FALSE
  )
}

#' MEAD with pseudo-individuals.
#'
#' `MEAD::MEAD()` on the cell-level raw counts of the panel, as in the
#' DeCovarT competitor runner (`get_R01(alpha = 0.1)` for the gene-gene
#' dependence set, falling back to `R01 = NULL`; gene filter off, so
#' every method sees the same panel; other settings at their defaults).
#'
#' MEAD weights genes by their variance between reference individuals.
#' The Suppinger reference has one cell line at 72 h and 96 h, so the
#' cells of each type are split at random into `n_pseudo` groups. Those
#' groups differ by sampling noise only: the weights measure the
#' precision of the mean, not donor variation.
#'
#' MEAD divides each cell by its total over the supplied genes, so its
#' estimate is the share of panel transcripts, `q`. It is converted to
#' cell proportions with the mean panel total per cell `s_j`:
#' `p_j = (q_j / s_j) / sum_k (q_k / s_k)`. MEAD itself drops the cells
#' with no count on the panel.
#'
#' @param bulk Bulk raw counts, genes by samples.
#' @param counts Reference raw counts, cells by genes.
#' @param cell_type Cell-type label of each reference cell.
#' @param mu Mean signature, genes by types (gives `s_j`).
#' @param n_pseudo Number of pseudo-individuals. Default `4`.
#' @param seed Seed of the random split.
#' @return Long table, see [proportions_long()], with the attribute
#'   `note` (route taken for `R01`, cells without counts).
deconvolute_mead <- function(
  bulk,
  counts,
  cell_type,
  mu,
  n_pseudo = 4L,
  seed = 20261110L
) {
  keep <- !is.na(cell_type)
  counts <- counts[keep, , drop = FALSE]
  cell_type <- cell_type[keep]
  genes <- colnames(counts)
  types <- colnames(mu)
  bulk <- as.matrix(bulk)[genes, , drop = FALSE]
  empty <- sum(rowSums(counts) == 0)
  individual <- withr::with_seed(seed, {
    out <- integer(length(cell_type))
    for (ty in types) {
      idx <- which(cell_type == ty)
      out[idx] <- sample(rep_len(seq_len(n_pseudo), length(idx)))
    }
    out
  })
  ref_counts <- t(counts)
  colnames(ref_counts) <- paste0("cell_", seq_len(ncol(ref_counts)))
  ref <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = ref_counts),
    colData = S4Vectors::DataFrame(
      cell_type = factor(cell_type, levels = types),
      individual = factor(paste0("pseudo_", individual))
    )
  )
  r01 <- tryCatch(MEAD::get_R01(bulk, alpha = 0.1), error = function(e) e)
  route <- "R01 from get_R01(alpha = 0.1)"
  if (inherits(r01, "error")) {
    route <- paste0("R01 = NULL (get_R01: ", conditionMessage(r01), ")")
    r01 <- NULL
  }
  fit <- suppressMessages(MEAD::MEAD(
    bulk,
    ref,
    cell_types = types,
    preprocessing_control = list(filter.gene = FALSE),
    R01 = r01
  ))
  transcript <- t(fit$p_hat[types, colnames(bulk), drop = FALSE])
  size <- colSums(mu)[types]
  cellular <- transcript / rep(size, each = nrow(transcript))
  cellular <- cellular / rowSums(cellular)
  out <- proportions_long(cellular)
  attr(out, "note") <- paste0(
    route,
    "; ",
    empty,
    " reference cells without a count on the panel"
  )
  out
}

#' ReDeconv through its Python implementation.
#'
#' Writes the bulk and the reference moments (mean and standard
#' deviation of the raw counts per cell, per type), runs
#' `scripts/python/run_redeconv.py` with `uv run --script`, and reads the
#' percentages. Supplying the moments bypasses ReDeconv's own
#' signature-gene selection and its transcriptome-size normalisation of
#' the reference: every panel gene is used, as for the other methods.
#' ReDeconv clips negative coefficients and rescales the positive ones;
#' the estimates are rescaled here to sum to one.
#'
#' @param bulk Bulk raw counts, genes by samples.
#' @param mu Mean signature, genes by types.
#' @param sd Standard deviations, genes by types.
#' @param work_dir Folder for the exchange files (created, kept).
#' @param script Path to the Python runner.
#' @return Long table, see [proportions_long()].
deconvolute_redeconv <- function(
  bulk,
  mu,
  sd,
  work_dir,
  script = file.path("scripts", "python", "run_redeconv.py")
) {
  genes <- rownames(mu)
  types <- colnames(mu)
  in_dir <- file.path(work_dir, "input")
  out_dir <- file.path(work_dir, "output")
  unlink(c(in_dir, out_dir), recursive = TRUE)
  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  bulk <- as.matrix(bulk)[genes, , drop = FALSE]
  # ReDeconv reads tab-separated files and takes the type name from the
  # column header; plain codes avoid spaces and accents there.
  codes <- paste0("type", seq_along(types))
  write_tsv <- function(m, path) {
    data.table::fwrite(
      data.table::data.table(gene = rownames(m), m),
      path,
      sep = "\t"
    )
  }
  write_tsv(bulk, file.path(in_dir, "bulk.tsv"))
  signature <- matrix(
    0,
    length(genes),
    2L * length(types),
    dimnames = list(
      genes,
      as.vector(rbind(paste0(codes, "_mean"), paste0(codes, "_std")))
    )
  )
  signature[, paste0(codes, "_mean")] <- mu[genes, types]
  signature[, paste0(codes, "_std")] <- sd[genes, types]
  write_tsv(signature, file.path(in_dir, "signature.tsv"))
  status <- system2(
    "uv",
    c("run", "--script", shQuote(script), shQuote(in_dir), shQuote(out_dir))
  )
  if (status != 0L) {
    stop("ReDeconv runner failed with status ", status, ".")
  }
  tab <- utils::read.csv(
    file.path(out_dir, "estimates.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  est <- as.matrix(tab[, -1L, drop = FALSE])
  rownames(est) <- tab[[1L]]
  est <- t(est)[colnames(bulk), codes, drop = FALSE]
  colnames(est) <- types
  if (any(rowSums(est) <= 0)) {
    stop("ReDeconv returned no positive coefficient for a sample.")
  }
  proportions_long(est / rowSums(est))
}

#' Robust linear regression of the bulk on the signature.
#'
#' `DeCovarT::deconvolute_ratios_rlm()`: Huber M-estimation
#' (`MASS::rlm()`), no intercept, no covariance, no sign constraint;
#' the coefficients are then divided by their sum. That function stops
#' when a coefficient is negative. In that case the same `MASS::rlm()`
#' fit is repeated here, negative coefficients are set to zero and the
#' rest rescaled; the attribute `note` counts those samples.
#'
#' @param bulk Bulk raw counts, genes by samples.
#' @param mu Mean signature, genes by types.
#' @return Long table, see [proportions_long()], with the attribute
#'   `note`.
deconvolute_rlm <- function(bulk, mu) {
  genes <- rownames(mu)
  bulk <- as.matrix(bulk)[genes, , drop = FALSE]
  fits <- lapply(colnames(bulk), function(s) {
    y <- bulk[, s]
    p <- tryCatch(
      DeCovarT::deconvolute_ratios_rlm(y, mu),
      error = function(e) NULL
    )
    was_clipped <- is.null(p)
    if (was_clipped) {
      p <- pmax(MASS::rlm(y ~ mu + 0, method = "M")$coefficients, 0)
      if (sum(p) <= 0) {
        stop("Robust regression gave no positive coefficient for ", s, ".")
      }
      p <- p / sum(p)
    }
    list(p = as.numeric(p), clipped = was_clipped)
  })
  est <- vapply(fits, `[[`, numeric(ncol(mu)), "p")
  dimnames(est) <- list(colnames(mu), colnames(bulk))
  clipped <- sum(vapply(fits, `[[`, logical(1), "clipped"))
  out <- proportions_long(t(est))
  attr(out, "note") <- if (clipped > 0L) {
    paste0(
      clipped,
      " of ",
      ncol(bulk),
      " samples had a negative coefficient, set to zero"
    )
  } else {
    ""
  }
  out
}

# ==========================================================================
# Benchmark metrics ----
# ==========================================================================

#' Aitchison distance between two compositions.
#'
#' Euclidean distance between the centred log-ratio vectors, which is
#' also the distance between ILR coordinates (ledger M1). Zeros are
#' handled by multiplicative replacement: every part below `floor` is
#' set to `floor` and the other parts are scaled so the total stays one.
#' The distance depends on the number of parts, so it compares
#' compositions of one time point only.
#'
#' @param p,q Compositions of the same parts, in the same order.
#' @param floor Smallest part after replacement.
#' @return One number, or `NA` when a part is not finite.
aitchison_distance <- function(p, q, floor = 0.001) {
  if (any(!is.finite(p)) || any(!is.finite(q))) {
    return(NA_real_)
  }
  replace_zeros <- function(x) {
    x <- x / sum(x)
    low <- x < floor
    if (any(low)) {
      x[!low] <- x[!low] * (1 - floor * sum(low)) / sum(x[!low])
      x[low] <- floor
    }
    x
  }
  clr <- function(x) log(x) - mean(log(x))
  sqrt(sum((clr(replace_zeros(p)) - clr(replace_zeros(q)))^2))
}

#' Agreement of estimated and observed proportions, per sample.
#'
#' For each method and bulk sample: Pearson correlation, root mean
#' squared error, mean absolute error and Aitchison distance
#' ([aitchison_distance()], zeros replaced at 0.001, and at 0.01 as a
#' sensitivity check) over the cell types of that
#' time. With three cell types the correlation rests on three points.
#' The correlation is `NA` when the estimates do not vary.
#'
#' @param estimates Long table with `method`, `sample_id`, `time`,
#'   `cell_type`, `proportion`.
#' @param observed Output of [observed_scored_fractions()].
#' @return One row per method and sample: `method`, `sample_id`, `time`,
#'   `n_types`, `pearson`, `rmse`, `mae`, `aitchison`, `aitchison_01`.
benchmark_metrics <- function(estimates, observed) {
  key <- paste(estimates$time, estimates$cell_type)
  estimates$observed <- observed$observed[
    match(key, paste(observed$time, observed$cell_type))
  ]
  if (anyNA(estimates$observed)) {
    stop("An estimated cell type has no observed proportion.")
  }
  groups <- split(
    estimates,
    list(estimates$method, estimates$sample_id),
    drop = TRUE
  )
  out <- lapply(groups, function(d) {
    ok <- is.finite(d$proportion)
    err <- d$proportion[ok] - d$observed[ok]
    varies <- sum(ok) >= 3L && stats::sd(d$proportion[ok]) > 0
    data.frame(
      method = d$method[[1L]],
      sample_id = d$sample_id[[1L]],
      time = d$time[[1L]],
      n_types = sum(ok),
      pearson = if (varies) {
        stats::cor(d$proportion[ok], d$observed[ok])
      } else {
        NA_real_
      },
      rmse = if (any(ok)) sqrt(mean(err^2)) else NA_real_,
      mae = if (any(ok)) mean(abs(err)) else NA_real_,
      aitchison = aitchison_distance(d$proportion, d$observed, 0.001),
      aitchison_01 = aitchison_distance(d$proportion, d$observed, 0.01),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}
