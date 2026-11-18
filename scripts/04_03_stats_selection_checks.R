# OPENBLAS_NUM_THREADS=1 Rscript --no-save --no-restore \
#   scripts/04_03_stats_selection_checks.R [times=72h,96h]
#
# Two checks of the rule that chooses the shrinkage level by the bulk
# likelihood (scripts/04_01_stats_deconvolution.R). Welch panel, DeCovarT
# with shrinkage of the raw covariance. Neither check selects anything.
#
# 1. Pseudo-bulks with a known composition. The reference moments come
#    from 80% of the cells of each type; pseudo-bulks are sums of 500
#    cells drawn from the other 20% (with replacement only when a type
#    has too few held-out cells), at three compositions. Two versions: as is, and multiplied gene by gene by
#    the bulk-to-reference ratio of tier (iii), which puts back a
#    platform shift. Reported: the level the rule picks and the error at
#    every level.
# 2. Leave-one-module-out predictive likelihood, real standard-Wnt bulks
#    at 72 h. Ten gene modules (Ward clustering on 1 - |pooled
#    correlation|); the model is fitted without one module and scores it
#    by its Gaussian law given the other genes.
#
# Output: output/deconvolution/lambda_checks.csv

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))
source(file.path(root, "R", "04_01_stats_deconvolution.R"))
source(file.path(root, "R", "04_02_stats_sieve_multivariate.R"))
check_blas_kernel()

cli <- commandArgs(trailingOnly = TRUE)
times <- strsplit(
  sub(
    "^times=",
    "",
    c(grep("^times=", cli, value = TRUE), "times=72h,96h")[[1L]]
  ),
  ","
)[[1L]]
lambda_levels <- c(0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 0.9, 1)
n_cells_bulk <- 500L
n_modules <- 10L
module_time <- "72h"
seed <- 20261117L
n_workers <- 16L

tier_01 <- file.path(root, "output", "feature_selection", "tier_01")
slots <- readRDS(file.path(tier_01, "scored_slots.rds"))
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
ratio_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_03",
  "platform_ratio.csv"
)
observed <- observed_scored_fractions(file.path(tier_01, "cell_counts.csv"))
out_dir <- file.path(root, "output", "deconvolution")

run <- function(tasks, fun) {
  out <- parallel::mclapply(
    seq_len(nrow(tasks)),
    function(k) fun(tasks[k, ]),
    mc.cores = n_workers,
    mc.preschedule = FALSE
  )
  bad <- !vapply(out, is.data.frame, logical(1))
  if (any(bad)) {
    stop(sum(bad), " task(s) failed: ", as.character(out[[which(bad)[[1L]]]]))
  }
  do.call(rbind, out)
}

# 1 · Pseudo-bulks --------------------------------------------------------------
pseudo <- do.call(
  rbind,
  lapply(times, function(tm) {
    panel <- read_panel_slot(slots, selected_dir, tm)
    types <- sort(unique(panel$cell_type[!is.na(panel$cell_type)]))
    held_out <- withr::with_seed(seed, {
      unlist(lapply(types, function(tp) {
        rows <- which(panel$cell_type %in% tp)
        sample(rows, ceiling(0.2 * length(rows)))
      }))
    })
    train <- setdiff(which(panel$cell_type %in% types), held_out)
    moments <- sample_moments(panel$y[train, ], panel$cell_type[train], types)
    at_time <- observed[observed$time == tm, ]
    obs <- at_time$observed[match(types, at_time$cell_type)]
    stopifnot(!anyNA(obs), abs(sum(obs) - 1) < 1e-8)
    compositions <- withr::with_seed(seed + 1L, {
      draw <- stats::rgamma(length(types), 1)
      list(
        observed = obs,
        uniform = rep(1 / length(types), length(types)),
        random = draw / sum(draw)
      )
    })
    ratio <- read_platform_ratio(ratio_path, tm, panel$genes)
    bulks <- withr::with_seed(
      seed + 2L,
      lapply(compositions, function(p) {
        counts <- as.integer(stats::rmultinom(1L, n_cells_bulk, p))
        rows <- unlist(lapply(seq_along(types), function(j) {
          pool <- held_out[panel$cell_type[held_out] %in% types[[j]]]
          # shortcut: a pool smaller than the draw is sampled with replacement,
          # so those cells are not independent; use more cells to upgrade.
          pool[sample.int(
            length(pool),
            counts[[j]],
            replace = counts[[j]] > length(pool)
          )]
        }))
        list(
          y = colSums(panel$y[rows, , drop = FALSE]),
          truth = counts / n_cells_bulk
        )
      })
    )
    tasks <- expand.grid(
      composition = names(bulks),
      version = c("as is", "with platform ratio"),
      lambda = lambda_levels,
      stringsAsFactors = FALSE
    )
    message("Pseudo-bulks ", tm, ": ", nrow(tasks), " fits.")
    run(tasks, function(task) {
      b <- bulks[[task$composition]]
      y <- if (task$version == "as is") b$y else b$y * ratio
      sigma <- build_sigma(moments, "shrinkage", lambda = task$lambda)
      fit <- fit_decovart_once(y, moments$mu, sigma, "lbfgsb")
      data.frame(
        check = "pseudo-bulk",
        time = tm,
        sample = task$composition,
        version = task$version,
        lambda = task$lambda,
        loglik = fit$diagnostics$loglik,
        rmse = sqrt(mean((fit$estimates$proportion - b$truth)^2)),
        stringsAsFactors = FALSE
      )
    })
  })
)

# 2 · Leave-one-module-out predictive likelihood ---------------------------------
moments <- readRDS(file.path(
  root,
  "output",
  "grn",
  "plugin",
  paste0("moments_", module_time, "_welch.rds")
))
genes <- moments$genes
module <- stats::cutree(
  stats::hclust(stats::as.dist(1 - abs(moments$r_pool)), "ward.D2"),
  n_modules
)
htseq <- read_bulk_htseq(file.path(
  root,
  "data",
  "raw",
  "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz"
))
bulk <- bulk_symbol_counts(htseq, slots$shared, genes)[genes, , drop = FALSE]
design <- bulk_sample_design(colnames(bulk))
samples <- design$sample_id[
  design$time == module_time & design$arm == "standard"
]
tasks <- expand.grid(
  sample = samples,
  lambda = lambda_levels,
  module = seq_len(n_modules),
  stringsAsFactors = FALSE
)
message(
  "Modules ",
  module_time,
  ": sizes ",
  paste(table(module), collapse = " "),
  "; ",
  nrow(tasks),
  " fits."
)
by_module <- run(tasks, function(task) {
  out <- module == task$module
  y <- bulk[, task$sample]
  sigma <- build_sigma(moments, "shrinkage", lambda = task$lambda)
  fit <- fit_decovart_once(
    y[!out],
    moments$mu[!out, ],
    sigma[!out, !out, ],
    "lbfgsb"
  )
  s <- fit$diagnostics$variance_scale *
    bulk_covariance(
      as.numeric(fit$p),
      sigma,
      moments$mu,
      between = FALSE,
      total = 1
    )
  resid <- y -
    fit$diagnostics$mean_scale * drop(moments$mu %*% as.numeric(fit$p))
  # Gaussian law of the left-out genes given the others.
  gain <- s[out, !out] %*%
    solve(s[!out, !out], cbind(resid[!out], s[!out, out]))
  data.frame(
    sample = task$sample,
    lambda = task$lambda,
    loglik = heldout_gaussian_logscore(
      gain[, 1L],
      s[out, out] - gain[, -1L, drop = FALSE],
      matrix(resid[out], 1L)
    ),
    stringsAsFactors = FALSE
  )
})
modules <- stats::aggregate(loglik ~ sample + lambda, by_module, sum)
modules <- data.frame(
  check = "leave-one-module-out",
  time = module_time,
  sample = modules$sample,
  version = "real bulk",
  lambda = modules$lambda,
  loglik = modules$loglik,
  rmse = NA_real_,
  stringsAsFactors = FALSE
)

checks <- rbind(pseudo, modules)
checks$selected <- as.logical(stats::ave(
  checks$loglik,
  checks$check,
  checks$time,
  checks$sample,
  checks$version,
  FUN = function(v) seq_along(v) == which.max(v)
))
utils::write.csv(
  checks,
  file.path(out_dir, "lambda_checks.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
message("Checks written.")
