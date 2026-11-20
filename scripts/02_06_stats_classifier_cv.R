# mkdir -p logs
# OPENBLAS_NUM_THREADS=1 nohup Rscript --no-save --no-restore \
#   scripts/02_06_stats_classifier_cv.R \
#   [times=48h,72h,96h] [methods=glmnet,msgl,plnlda,ranger,xgboost,splsda] \
#   [stages=nested,tables] [variants=primary,covariate,combat] \
#   [workers=20] [subsample=1] \
#   > "logs/02_06_stats_classifier_cv_$(date +%F).log" 2>&1 &
#
# Classifier-based marker selection against the two-stage panel
# (chapter 2, #sec-fs-classifiers). Each time point is analysed on its
# own: outer 5-fold cross-validation stratified by type, inner 3-fold
# tuning on the outer training cells only. The resampling unit is the
# cell (one specimen, no donor column). Methods run one after another;
# within a method the outer tasks (time x fold x variant) run in parallel.
#
# Variants: `primary` (every time); at 48 h only, `covariate` (cell line
# as an unpenalised covariate where the method allows) and `combat`
# (counts adjusted by sva::ComBat_seq(batch, group = NULL) on all cells
# of the slot: no labels used).
#
# Universe: the genes tested in tier (ii) (presence gate) for glmnet,
# ranger, xgboost and sPLS-DA; the global top 2,000 highly variable genes
# of the slot (Seurat v3 vst, select_hvg_vst()) for msgl and PLN-LDA,
# whose cost grows too fast with the number of genes.
#
# Inputs: output/feature_selection/tier_01/scored_slots.rds,
# output/feature_selection/tier_02/cross_cell.rds.
# Outputs: output/feature_selection/classifiers/ (fits/ git-ignored).
# msgl and sglOptim are archived on CRAN: installed from GitHub
# (nielsrhansen/sglOptim, nielsrhansen/msgl, recorded in renv.lock) with a
# user Makevars holding `CXX = g++ -std=gnu++17` and
# `CXXFLAGS += -fpermissive` (R_MAKEVARS_USER); the default C++20 rejects
# their templated constructors.
# `subsample` < 1 keeps that share of the cells of each type (smoke test,
# written to classifiers_subsample/).

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))
source(file.path(root, "R", "02_06_stats_classifier_cv.R"))

cli <- commandArgs(trailingOnly = TRUE)
arg <- function(key, default) {
  hit <- grep(paste0("^", key, "="), cli, value = TRUE)
  if (length(hit) == 0L) default else strsplit(sub(paste0("^", key, "="), "", hit[[1L]]), ",")[[1L]]
}
times <- arg("times", c("48h", "72h", "96h"))
methods <- arg("methods", classifier_methods)
stages <- arg("stages", c("nested", "tables"))
variants <- arg("variants", c("primary", "covariate", "combat"))
n_workers <- as.integer(arg("workers", "20"))
subsample <- as.numeric(arg("subsample", "1"))
stopifnot(all(methods %in% classifier_methods))

seed <- 20261121L
K_outer <- 5L
K_inner <- 3L
stable_threshold <- 0.6
sensitivity_time <- "48h"
n_hvg <- 2000L
hvg_methods <- c("msgl", "plnlda")
settings <- list(
  gamma = 0.5,
  n_top = 50L,
  glmnet_alpha = c(0.1, 0.25, 0.5, 0.75, 1),
  n_lambda = 40L,
  lambda_min_ratio = 0.01,
  msgl_alpha = c(0.25, 0.5, 0.75),
  n_modules = 50L,
  pln_covariance = "diagonal",
  rf_trees = 500L,
  rf_mtry_mult = c(0.5, 1, 2),
  rf_min_n = c(1L, 10L),
  xgb_trees = 300L,
  xgb_depth = c(3L, 6L),
  xgb_learn_rate = c(0.05, 0.2),
  xgb_mtry_mult = 2,
  splsda_keepx = c(10, 25, 50, 100, 200),
  seed = seed
)

tier1 <- file.path(root, "output", "feature_selection", "tier_01")
tier2 <- file.path(root, "output", "feature_selection", "tier_02")
out_dir <- file.path(
  root, "output", "feature_selection",
  if (subsample < 1) "classifiers_subsample" else "classifiers"
)
fits_dir <- file.path(out_dir, "fits")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(out_dir, ".gitignore"))) {
  writeLines(c("/fits/", "/predictions.rds"), file.path(out_dir, ".gitignore"))
}
progress <- progress_path("02_06_stats_classifier_cv")
write_table <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, name), row.names = FALSE, fileEncoding = "UTF-8")
}
checkpoint <- function(method, variant, tm, k) {
  file.path(fits_dir, sprintf("%s__%s__%s__fold%d.rds", method, variant, tm, k))
}

# Data of one slot: log1p(x/s) (sparse, cells x genes), raw counts of the
# HVGs (dense), labels, folds.
load_slot <- function(tm, variant) {
  prepared <- readRDS(file.path(tier1, "scored_slots.rds"))
  cc <- readRDS(file.path(tier2, "cross_cell.rds"))
  slot <- prepared$slots[[tm]]
  gate <- cc$tests[[tm]]
  genes <- gate$gene[gate$pass_presence %in% TRUE]
  cells <- which(!is.na(slot$cell_type))
  if (subsample < 1) {
    ct0 <- as.character(slot$cell_type[cells])
    cells <- sort(unlist(withr::with_seed(seed, lapply(split(cells, ct0), function(i) {
      i[sort(sample.int(length(i), max(15L, ceiling(subsample * length(i)))))]
    }))))
  }
  counts <- slot$counts[genes, cells, drop = FALSE]
  if (variant == "combat") {
    cache <- file.path(fits_dir, sprintf("combat_%s.rds", tm))
    if (!file.exists(cache)) {
      adj <- sva::ComBat_seq(as.matrix(counts), batch = as.character(slot$batch[cells]), group = NULL)
      saveRDS(adj, cache)
    }
    counts <- Matrix::Matrix(readRDS(cache), sparse = TRUE)
  }
  lib <- as.numeric(slot$library_size[cells])
  sf <- lib / stats::median(lib)
  hvg <- select_hvg_vst(counts, n_top = n_hvg)
  X <- Matrix::t(counts %*% Matrix::Diagonal(x = 1 / sf))
  X@x <- log1p(X@x)
  colnames(X) <- genes
  y <- droplevels(factor(slot$cell_type[cells]))
  hv <- hvg$gene[hvg$hvg]
  list(
    X = X,
    Y = t(as.matrix(counts[hv, , drop = FALSE])),
    hvg = hvg,
    lib = lib,
    y = y,
    batch = as.character(slot$batch[cells]),
    cells = (slot$cells %||% colnames(slot$counts))[cells],
    outer = stratified_fold_ids(y, K_outer, seed)
  )
}

run_task <- function(d, method, variant, k) {
  started <- proc.time()[["elapsed"]]
  train <- which(d$outer != k)
  test <- which(d$outer == k)
  inner <- stratified_fold_ids(droplevels(d$y[train]), K_inner, seed + 1000L + k)
  cols <- if (method %in% hvg_methods) d$hvg$gene[d$hvg$hvg] else colnames(d$X)
  covariate <- if (variant == "covariate") as.numeric(d$batch == sort(unique(d$batch))[[2L]]) else NULL
  r <- fit_classifier(
    method,
    d$X[, cols, drop = FALSE],
    d$Y,
    d$lib,
    d$y,
    train,
    test,
    inner,
    covariate = covariate,
    settings = settings
  )
  types <- levels(d$y)
  pred <- if (is.null(r$prob)) r$pred else factor(types[max.col(r$prob, ties.method = "first")], levels = types)
  r$pred <- data.frame(cell = d$cells[test], truth = d$y[test], pred = pred)
  if (!is.null(r$prob)) r$pred <- cbind(r$pred, as.data.frame(r$prob, optional = TRUE))
  r$prob <- NULL
  r$universe <- length(cols)
  r$minutes <- (proc.time()[["elapsed"]] - started) / 60
  r
}

# 1 · Nested cross-validation ----------------------------------------------------
# Methods run one after another; all outer tasks of a method (time x
# variant x fold) share one mclapply. The `covariate` variant reuses the
# primary data of its time.
if ("nested" %in% stages) {
  data_key <- function(tm, variant) paste(tm, ifelse(variant == "combat", "combat", "primary"))
  slot_data <- list()
  for (method in methods) {
    tasks <- expand.grid(k = seq_len(K_outer), variant = variants, time = times, stringsAsFactors = FALSE)
    tasks <- tasks[tasks$variant == "primary" | tasks$time == sensitivity_time, , drop = FALSE]
    tasks$file <- checkpoint(method, tasks$variant, tasks$time, tasks$k)
    tasks <- tasks[!file.exists(tasks$file), , drop = FALSE]
    if (nrow(tasks) == 0L) next
    for (key in unique(data_key(tasks$time, tasks$variant))) {
      if (is.null(slot_data[[key]])) {
        parts <- strsplit(key, " ")[[1L]]
        message("Loading ", key)
        slot_data[[key]] <- load_slot(parts[[1L]], parts[[2L]])
        if (parts[[2L]] == "primary") write_table(slot_data[[key]]$hvg, sprintf("hvg_%s.csv", parts[[1L]]))
      }
    }
    workers <- min(n_workers, nrow(tasks))
    log_progress(progress, "plan", method, nrow(tasks), workers)
    res <- parallel::mclapply(
      seq_len(nrow(tasks)),
      function(i) {
        t0 <- proc.time()[["elapsed"]]
        label <- paste(method, tasks$variant[[i]], tasks$time[[i]])
        d <- slot_data[[data_key(tasks$time[[i]], tasks$variant[[i]])]]
        out <- try(run_task(d, method, tasks$variant[[i]], tasks$k[[i]]), silent = TRUE)
        failed <- inherits(out, "try-error")
        if (!failed) saveRDS(out, tasks$file[[i]])
        log_progress(
          progress, "done", label, tasks$k[[i]],
          round(proc.time()[["elapsed"]] - t0, 1),
          if (failed) "failed" else "ok",
          if (failed) as.character(out) else ""
        )
        !failed
      },
      mc.cores = workers,
      mc.preschedule = FALSE
    )
    ok <- vapply(res, isTRUE, logical(1))
    if (!all(ok)) {
      message(method, ": ", sum(!ok), " task(s) failed or were killed; see ", progress, ". Rerun to retry them.")
    }
    gc()
  }
}

# 2 · Tables ---------------------------------------------------------------------
if ("tables" %in% stages) {
  files <- list.files(fits_dir, "__fold[0-9]+\\.rds$", full.names = TRUE)
  meta <- do.call(rbind, strsplit(sub("\\.rds$", "", basename(files)), "__"))
  meta <- data.frame(
    method = meta[, 1L],
    variant = meta[, 2L],
    time_point = meta[, 3L],
    fold = as.integer(sub("fold", "", meta[, 4L])),
    file = files,
    stringsAsFactors = FALSE
  )
  fits <- lapply(meta$file, readRDS)
  tag <- function(i, x) {
    if (is.null(x) || nrow(x) == 0L) return(NULL)
    cbind(meta[i, c("time_point", "method", "variant", "fold")], x, row.names = NULL)
  }
  metrics <- do.call(rbind, lapply(seq_along(fits), function(i) {
    p <- fits[[i]]$pred
    types <- levels(p$truth)
    prob <- if (all(types %in% names(p))) as.matrix(p[, types]) else NULL
    tag(i, cbind(
      classification_metrics(p$truth, p$pred, prob),
      n_selected = length(unique(fits[[i]]$selected$gene)),
      universe = fits[[i]]$universe,
      minutes = fits[[i]]$minutes
    ))
  }))
  write_table(metrics, "metrics.csv")
  recall <- do.call(rbind, lapply(seq_along(fits), function(i) {
    p <- fits[[i]]$pred
    tag(i, data.frame(
      type = levels(p$truth),
      n_test = as.integer(table(p$truth)),
      recall = vapply(levels(p$truth), function(k) mean(p$pred[p$truth == k] == k), numeric(1))
    ))
  }))
  write_table(recall, "per_class_recall.csv")
  confusion <- do.call(rbind, lapply(seq_along(fits), function(i) {
    tag(i, as.data.frame(table(truth = fits[[i]]$pred$truth, pred = fits[[i]]$pred$pred), responseName = "n"))
  }))
  confusion <- stats::aggregate(n ~ time_point + method + variant + truth + pred, confusion, sum)
  write_table(confusion, "confusion.csv")
  for (m in unique(meta$method)) {
    tu <- do.call(rbind, lapply(which(meta$method == m), function(i) tag(i, fits[[i]]$tuning)))
    if (!is.null(tu)) write_table(tu, sprintf("tuning_%s.csv", m))
  }
  selected <- do.call(rbind, lapply(seq_along(fits), function(i) tag(i, fits[[i]]$selected)))
  freq <- stats::aggregate(
    fold ~ time_point + method + variant + type + gene,
    unique(selected[, c("time_point", "method", "variant", "type", "gene", "fold")]),
    length
  )
  names(freq)[names(freq) == "fold"] <- "n_folds"
  freq$frequency <- freq$n_folds / K_outer
  # Direction: majority sign over the folds that selected the gene.
  dir_tab <- stats::aggregate(
    (direction == "up") ~ time_point + method + variant + type + gene,
    selected,
    mean
  )
  names(dir_tab)[ncol(dir_tab)] <- "share_up"
  freq <- merge(freq, dir_tab, by = c("time_point", "method", "variant", "type", "gene"))
  freq$stable <- freq$frequency >= stable_threshold
  write_table(freq, "selection_frequency.csv")
  write_table(freq[freq$stable, , drop = FALSE], "stable_sets.csv")
  importance <- do.call(rbind, lapply(seq_along(fits), function(i) {
    im <- fits[[i]]$importance
    im <- im[order(-im$score), , drop = FALSE]
    tag(i, utils::head(im[im$score > 0, , drop = FALSE], 30L))
  }))
  write_table(importance, "importance_top30.csv")
  sens <- metrics[metrics$time_point == sensitivity_time, , drop = FALSE]
  write_table(sens, "metrics_sensitivity_48h.csv")
  st <- freq[freq$stable & freq$time_point == sensitivity_time, , drop = FALSE]
  jac <- do.call(rbind, lapply(split(st, list(st$method, st$type), drop = TRUE), function(s) {
    base <- s$gene[s$variant == "primary"]
    do.call(rbind, lapply(setdiff(unique(s$variant), "primary"), function(v) {
      alt <- s$gene[s$variant == v]
      data.frame(
        method = s$method[[1L]],
        type = s$type[[1L]],
        variant = v,
        n_primary = length(base),
        n_variant = length(alt),
        jaccard = if (length(union(base, alt)) == 0L) NA_real_ else length(intersect(base, alt)) / length(union(base, alt))
      )
    }))
  }))
  if (!is.null(jac)) write_table(jac, "jaccard_sensitivity_48h.csv")
  saveRDS(
    list(
      predictions = lapply(seq_along(fits), function(i) tag(i, fits[[i]]$pred)),
      diag = lapply(seq_along(fits), function(i) c(as.list(meta[i, 1:4]), fits[[i]]$diag))
    ),
    file.path(out_dir, "predictions.rds")
  )
  utils::capture.output(
    list(settings = settings, K_outer = K_outer, K_inner = K_inner, seed = seed, n_hvg = n_hvg, hvg_methods = hvg_methods, stable_threshold = stable_threshold),
    file = file.path(out_dir, "settings.txt")
  )
  utils::capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
  message("Wrote ", out_dir)
}
