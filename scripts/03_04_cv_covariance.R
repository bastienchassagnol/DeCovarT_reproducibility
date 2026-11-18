# mkdir -p logs
# OPENBLAS_NUM_THREADS=1 nohup Rscript --no-save --no-restore \
#   scripts/03_04_cv_covariance.R \
#   > "logs/03_04_cv_covariance_$(date +%F).log" 2>&1 &
#
# Optional arguments, `key=value`, comma-separated values; each defaults
# to everything:
#   times=48h,72h,96h   gene_sets=r2_05,r2_01,welch
#   stages=select,nested,fit,tables
#
# Cell-wise choice of the covariance regularisation of each scored type,
# on the held-out Gaussian log-score of cells. The bulk samples are never
# read. Arms (R/03_01_stats_covariance_estimation.R), all on the
# correlation scale with the sample variances kept:
#   shrinkage   (1 - lambda) R_j + lambda I
#   two_param   (1 - lambda) {(1 - alpha) R_j + alpha Rbar} + lambda I,
#               Rbar the pooled within-type correlation of the time point
#   factor      B_j B_j' + Psi_j, rank q_j
#   factor_shared  B B' + Psi_j, loadings shared by the types of a time
#   diagonal    lambda = 1 (also rank 0)
#
# Stages:
#   select  cross-validated score of every tuning value, per type x time
#           x gene set. Leave-one-out under 100 cells, 5-fold otherwise
#           (5-fold for the two-parameter arm, 5 joint folds stratified by
#           type for the shared loadings).
#   nested  outer 5 folds stratified by type; each rule chooses its tuning
#           value on the outer training cells only and is scored on the
#           outer test cells. Scores are paired by fold and type.
#   fit     factor models on all cells (used by the deconvolution), with
#           AIC and BIC from their free-parameter counts.
#   tables  CSV tables from the checkpoints.
#
# Inputs:
#   output/feature_selection/tier_01/scored_slots.rds
#   output/grn/plugin/gene_sets.csv
# Outputs, output/grn/cv/:
#   cv_tuning_scores.csv      one row per unit x arm x tuning value
#   cv_selected_tuning.csv    one row per unit x rule
#   nested_cv_by_fold.csv     one row per unit x outer fold x rule
#   nested_cv_summary.csv     one row per unit x rule
#   information_criteria.csv  one row per unit x factor rank
#   factors_<time>_<gene_set>.rds  factor fits on all cells (not in Git)
#   fits/                     checkpoints (not in Git)

# ==========================================================================
# SECTION 0 · Dependencies, settings, paths ----
# ==========================================================================

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))
check_blas_kernel()

time_levels <- c("48h", "72h", "96h")
gene_set_levels <- c("r2_05", "r2_01", "welch")
stage_levels <- c("select", "nested", "fit", "tables")

# 13 levels on a log scale up to 0.1, then steps of 1/32, and 0.9.
lambda_grid <- sort(unique(c(
  10^seq(-3, -1, length.out = 13L),
  seq(0.125, 1, by = 1 / 32),
  0.9
)))
alpha_grid <- c(0.25, 0.5, 0.75)
ranks <- c(0L, 1L, 2L, 3L, 5L, 8L, 12L, 20L, 30L)
ranks_shared <- c(0L, 1L, 2L, 3L, 5L, 8L, 12L, 20L)
K <- 5L
K_outer <- 5L
loo_below <- 100L
seed <- 20261116L
n_workers <- 24L
# A task holds a few genes x genes matrices and the cells of its time.
memory_budget_gb <- 150

cli <- commandArgs(trailingOnly = TRUE)
cli_value <- function(key, default) {
  hit <- grep(paste0("^", key, "="), cli, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  strsplit(sub(paste0("^", key, "="), "", hit[[1L]]), ",", fixed = TRUE)[[1L]]
}
run_times <- match.arg(cli_value("times", time_levels), time_levels, TRUE)
run_gene_sets <- match.arg(
  cli_value("gene_sets", gene_set_levels),
  gene_set_levels,
  TRUE
)
run_stages <- match.arg(cli_value("stages", stage_levels), stage_levels, TRUE)

slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
gene_sets_path <- file.path(root, "output", "grn", "plugin", "gene_sets.csv")
out_dir <- file.path(root, "output", "grn", "cv")
fits_dir <- file.path(out_dir, "fits")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
progress <- progress_path("03_04_cv_covariance")
if (!file.exists(file.path(out_dir, ".gitignore"))) {
  writeLines(c("/fits/", "/factors_*.rds"), file.path(out_dir, ".gitignore"))
}

if (!file.exists(gene_sets_path)) {
  stop("Missing ", gene_sets_path, ". Run scripts/03_05 first.")
}
gene_sets <- utils::read.csv(
  gene_sets_path,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
slots <- if (any(run_stages != "tables")) readRDS(slots_path) else NULL

#' Cells of every scored type of one time, on one gene set.
read_units <- function(time_level, gene_set) {
  genes <- gene_sets$gene[
    gene_sets$time_point == time_level & gene_sets$gene_set == gene_set
  ]
  panel <- read_panel_slot(slots, selected_dir, time_level, genes = genes)
  types <- sort(unique(panel$cell_type[!is.na(panel$cell_type)]))
  X <- lapply(types, function(tp) {
    panel$y[which(panel$cell_type %in% tp), , drop = FALSE]
  })
  names(X) <- types
  X
}

workers_for <- function(n_tasks, n_genes, copies = 10) {
  max(
    1L,
    min(
      n_workers,
      n_tasks,
      floor(memory_budget_gb * 1e9 / (copies * 8 * n_genes^2 + 2e8))
    )
  )
}

checkpoint <- function(stage, time_level, gene_set) {
  file.path(
    fits_dir,
    paste0(paste(stage, time_level, gene_set, sep = "__"), ".rds")
  )
}

#' Stop on the first failed task of an mclapply() result.
check_tasks <- function(results, what) {
  failed <- vapply(results, inherits, logical(1), "try-error")
  empty <- vapply(results, is.null, logical(1))
  if (any(failed | empty)) {
    stop(
      what,
      ": ",
      sum(failed | empty),
      " task(s) failed. First message: ",
      if (any(failed)) {
        as.character(results[[which(failed)[[1L]]]])
      } else {
        "worker killed (memory?)"
      }
    )
  }
  results
}

#' Joint folds: every type is split into K folds.
stratified_folds <- function(X, K, seed) {
  lapply(seq_along(X), function(j) make_cv_folds(nrow(X[[j]]), K, seed + j))
}

#' Wrap a task function so that each task writes a progress line.
logged <- function(slot, n_tasks, workers, fun) {
  log_progress(progress, "plan", slot, n_tasks, workers)
  function(i) {
    started <- proc.time()[["elapsed"]]
    out <- try(fun(i), silent = TRUE)
    failed <- inherits(out, "try-error")
    log_progress(
      progress,
      "done",
      slot,
      i,
      round(proc.time()[["elapsed"]] - started, 1),
      if (failed) "failed" else "ok",
      if (failed) as.character(out) else ""
    )
    out
  }
}

which_max <- function(x) {
  ok <- which(is.finite(x))
  ok[[which.max(x[ok])]]
}

# ==========================================================================
# SECTION 1 · Selection on all cells of each unit ----
# ==========================================================================

if ("select" %in% run_stages) {
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      path <- checkpoint("select", tm, gs)
      if (file.exists(path)) {
        next
      }
      started <- proc.time()[["elapsed"]]
      X <- read_units(tm, gs)
      types <- names(X)
      n_genes <- ncol(X[[1L]])
      scatter <- lapply(X, centred_scatter)
      total <- Reduce(`+`, scatter)
      joint <- stratified_folds(X, K, seed)
      tasks <- rbind(
        data.frame(kind = "unit", j = seq_along(X), k = NA_integer_),
        data.frame(kind = "two_param", j = seq_along(X), k = NA_integer_),
        data.frame(kind = "shared", j = NA_integer_, k = seq_len(K))
      )
      message(
        "Select ",
        tm,
        " / ",
        gs,
        ": ",
        n_genes,
        " genes, ",
        nrow(tasks),
        " tasks."
      )
      results <- parallel::mclapply(
        seq_len(nrow(tasks)),
        logged(
          paste("select", tm, gs),
          nrow(tasks),
          workers_for(nrow(tasks), n_genes),
          function(i) {
            task <- tasks[i, ]
            if (task$kind == "unit") {
              folds <- unit_folds(nrow(X[[task$j]]), seed, K, loo_below)
              cv_unit_arms(X[[task$j]], folds, lambda_grid, ranks)
            } else if (task$kind == "two_param") {
              folds <- make_cv_folds(nrow(X[[task$j]]), K, seed + 100L)
              cv_unit_arms(
                X[[task$j]],
                folds,
                lambda_grid,
                ranks = NULL,
                alpha_grid = alpha_grid,
                other_scatter = total - scatter[[task$j]]
              )[c("two_param", "const")]
            } else {
              score_shared_factor(
                Map(function(x, f) x[f != task$k, , drop = FALSE], X, joint),
                Map(function(x, f) x[f == task$k, , drop = FALSE], X, joint),
                ranks_shared
              )
            }
          }
        ),
        mc.cores = workers_for(nrow(tasks), n_genes),
        mc.preschedule = FALSE
      )
      check_tasks(results, paste("select", tm, gs))
      names(results) <- paste(tasks$kind, tasks$j, tasks$k)
      saveRDS(
        list(
          time = tm,
          gene_set = gs,
          types = types,
          n_cells = vapply(X, nrow, integer(1)),
          n_genes = n_genes,
          unit = results[tasks$kind == "unit"],
          two_param = results[tasks$kind == "two_param"],
          shared = results[tasks$kind == "shared"],
          elapsed_sec = proc.time()[["elapsed"]] - started
        ),
        path
      )
      message(
        "  done in ",
        round(proc.time()[["elapsed"]] - started),
        " s."
      )
      rm(X, scatter, total, results)
      invisible(gc())
    }
  }
}

#' Curves and choices of one `select` checkpoint.
#'
#' @return A list with `curves` (one row per unit, arm and tuning value)
#'   and `selected` (one row per unit and rule).
read_selection <- function(path) {
  res <- readRDS(path)
  types <- res$types
  curve <- function(tp, arm, alpha, tuning, scores, const, scheme) {
    data.frame(
      time = res$time,
      gene_set = res$gene_set,
      cell_type = tp,
      arm = arm,
      alpha = alpha,
      tuning = tuning,
      mean_score = colMeans(scores) + mean(const),
      se_score = apply(scores, 2L, stats::sd) / sqrt(nrow(scores)),
      folds = nrow(scores),
      scheme = scheme,
      stringsAsFactors = FALSE
    )
  }
  pick <- function(tp, rule, arm, alpha, tuning, score) {
    data.frame(
      time = res$time,
      gene_set = res$gene_set,
      cell_type = tp,
      rule = rule,
      arm = arm,
      alpha = alpha,
      tuning = tuning,
      mean_score = score,
      stringsAsFactors = FALSE
    )
  }
  curves <- list()
  selected <- list()
  for (j in seq_along(types)) {
    tp <- types[[j]]
    unit <- res$unit[[j]]
    scheme <- if (res$n_cells[[j]] < loo_below) "leave-one-out" else "5-fold"
    full <- function(scores, at) mean(scores[, at]) + mean(unit$const)
    shr <- summarise_path(unit$shrinkage, lambda_grid, "larger")
    curves[[length(curves) + 1L]] <- curve(
      tp,
      "shrinkage",
      0,
      lambda_grid,
      unit$shrinkage,
      unit$const,
      scheme
    )
    at <- function(value) match(value, lambda_grid)
    selected[[length(selected) + 1L]] <- rbind(
      pick(
        tp,
        "shrinkage, per type",
        "shrinkage",
        0,
        shr$min,
        full(unit$shrinkage, at(shr$min))
      ),
      pick(
        tp,
        "shrinkage, per type, 1se",
        "shrinkage",
        0,
        shr$one_se,
        full(unit$shrinkage, at(shr$one_se))
      ),
      pick(tp, "diagonal", "shrinkage", 0, 1, full(unit$shrinkage, at(1)))
    )
    fac <- summarise_path(unit$factor, ranks, "smaller")
    curves[[length(curves) + 1L]] <- curve(
      tp,
      "factor",
      NA_real_,
      ranks,
      unit$factor,
      unit$const,
      scheme
    )
    selected[[length(selected) + 1L]] <- rbind(
      pick(
        tp,
        "factor, per type",
        "factor",
        NA_real_,
        fac$min,
        full(unit$factor, match(fac$min, ranks))
      ),
      pick(
        tp,
        "factor, per type, 1se",
        "factor",
        NA_real_,
        fac$one_se,
        full(unit$factor, match(fac$one_se, ranks))
      )
    )
    two <- res$two_param[[j]]
    pairs <- expand.grid(lambda = lambda_grid, alpha = alpha_grid)
    curves[[length(curves) + 1L]] <- curve(
      tp,
      "two_param",
      pairs$alpha,
      pairs$lambda,
      two$two_param,
      two$const,
      "5-fold"
    )
    best <- which_max(colMeans(two$two_param))
    selected[[length(selected) + 1L]] <- pick(
      tp,
      "two-parameter, per type",
      "two_param",
      pairs$alpha[[best]],
      pairs$lambda[[best]],
      mean(two$two_param[, best]) + mean(two$const)
    )
  }
  # Shared loadings: folds in the list, types by ranks in each.
  shared_scores <- simplify2array(lapply(res$shared, `[[`, "scores"))
  shared_const <- vapply(res$shared, `[[`, numeric(length(types)), "const")
  for (j in seq_along(types)) {
    scores <- t(shared_scores[j, , ])
    curves[[length(curves) + 1L]] <- curve(
      types[[j]],
      "factor_shared",
      NA_real_,
      ranks_shared,
      scores,
      matrix(shared_const, nrow = length(types))[j, ],
      "5 joint folds"
    )
  }
  curves <- do.call(rbind, curves)
  selected <- do.call(rbind, selected)
  # Common choices: the value with the highest sum over types of the
  # mean held-out score (one vote per type).
  common <- function(arm, rule) {
    d <- curves[curves$arm == arm, ]
    key <- paste(d$alpha, d$tuning)
    total <- tapply(d$mean_score, key, sum)
    best <- names(total)[[which_max(total)]]
    hit <- d[key == best, ]
    pick(hit$cell_type, rule, arm, hit$alpha, hit$tuning, hit$mean_score)
  }
  selected <- rbind(
    selected,
    common("shrinkage", "shrinkage, common"),
    common("two_param", "two-parameter, common"),
    common("factor_shared", "factor, shared loadings")
  )
  rownames(curves) <- NULL
  rownames(selected) <- NULL
  list(curves = curves, selected = selected)
}

# ==========================================================================
# SECTION 2 · Nested comparison of the rules ----
# ==========================================================================

if ("nested" %in% run_stages) {
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      path <- checkpoint("nested", tm, gs)
      if (file.exists(path)) {
        next
      }
      started <- proc.time()[["elapsed"]]
      X <- read_units(tm, gs)
      types <- names(X)
      n_genes <- ncol(X[[1L]])
      outer <- stratified_folds(X, K_outer, seed + 1000L)
      tasks <- rbind(
        expand.grid(kind = "unit", j = seq_along(X), k = seq_len(K_outer)),
        data.frame(kind = "shared", j = NA_integer_, k = seq_len(K_outer))
      )
      message(
        "Nested ",
        tm,
        " / ",
        gs,
        ": ",
        n_genes,
        " genes, ",
        nrow(tasks),
        " tasks."
      )
      results <- parallel::mclapply(
        seq_len(nrow(tasks)),
        logged(
          paste("nested", tm, gs),
          nrow(tasks),
          workers_for(nrow(tasks), n_genes, copies = 14),
          function(i) {
            task <- tasks[i, ]
            k <- task$k
            train <- Map(function(x, f) x[f != k, , drop = FALSE], X, outer)
            test <- Map(function(x, f) x[f == k, , drop = FALSE], X, outer)
            if (task$kind == "shared") {
              inner <- stratified_folds(train, K, seed + 2000L + k)
              inner_scores <- lapply(seq_len(K), function(m) {
                score_shared_factor(
                  Map(function(x, f) x[f != m, , drop = FALSE], train, inner),
                  Map(function(x, f) x[f == m, , drop = FALSE], train, inner),
                  ranks_shared
                )$scores
              })
              return(list(
                inner = Reduce(`+`, inner_scores) / K,
                outer = score_shared_factor(train, test, ranks_shared)
              ))
            }
            j <- task$j
            other <- Reduce(`+`, lapply(train[-j], centred_scatter))
            n_train <- nrow(train[[j]])
            inner <- cv_unit_arms(
              train[[j]],
              unit_folds(n_train, seed + 3000L + k, K, loo_below),
              lambda_grid,
              ranks
            )
            inner_two <- cv_unit_arms(
              train[[j]],
              make_cv_folds(n_train, K, seed + 4000L + k),
              lambda_grid,
              ranks = NULL,
              alpha_grid = alpha_grid,
              other_scatter = other
            )
            list(
              inner = list(
                shrinkage = colMeans(inner$shrinkage),
                factor = colMeans(inner$factor),
                two_param = colMeans(inner_two$two_param)
              ),
              outer = score_unit_arms(
                train[[j]],
                test[[j]],
                lambda_grid,
                ranks,
                alpha_grid,
                other
              )
            )
          }
        ),
        mc.cores = workers_for(nrow(tasks), n_genes, copies = 14),
        mc.preschedule = FALSE
      )
      check_tasks(results, paste("nested", tm, gs))
      saveRDS(
        list(
          time = tm,
          gene_set = gs,
          types = types,
          tasks = tasks,
          results = results,
          elapsed_sec = proc.time()[["elapsed"]] - started
        ),
        path
      )
      message(
        "  done in ",
        round(proc.time()[["elapsed"]] - started),
        " s."
      )
      rm(X, results)
      invisible(gc())
    }
  }
}

#' Outer-fold scores of every rule from one `nested` checkpoint.
read_nested <- function(path) {
  res <- readRDS(path)
  types <- res$types
  pairs <- expand.grid(lambda = lambda_grid, alpha = alpha_grid)
  # A local environment, so the nested adder updates this list with `<-`
  # and does not search enclosing frames.
  rows <- new.env(parent = emptyenv())
  rows$items <- list()
  add <- function(tp, k, rule, alpha, tuning, score, n_test) {
    rows$items[[length(rows$items) + 1L]] <- data.frame(
      time = res$time,
      gene_set = res$gene_set,
      cell_type = tp,
      outer_fold = k,
      rule = rule,
      alpha = alpha,
      tuning = tuning,
      score = score,
      n_test = n_test,
      stringsAsFactors = FALSE
    )
  }
  for (k in seq_len(K_outer)) {
    units <- lapply(seq_along(types), function(j) {
      res$results[[which(
        res$tasks$kind == "unit" & res$tasks$j == j & res$tasks$k == k
      )]]
    })
    shared <- res$results[[which(
      res$tasks$kind == "shared" & res$tasks$k == k
    )]]
    common_lambda <- which_max(Reduce(
      `+`,
      lapply(units, function(u) u$inner$shrinkage)
    ))
    common_pair <- which_max(Reduce(
      `+`,
      lapply(units, function(u) u$inner$two_param)
    ))
    shared_rank <- which_max(colSums(shared$inner))
    for (j in seq_along(types)) {
      u <- units[[j]]
      tp <- types[[j]]
      const <- u$outer$const
      n_test <- u$outer$n_test
      two_outer <- as.numeric(t(u$outer$two_param))
      own <- which_max(u$inner$shrinkage)
      add(
        tp,
        k,
        "diagonal",
        0,
        1,
        u$outer$shrinkage[[match(1, lambda_grid)]] + const,
        n_test
      )
      add(
        tp,
        k,
        "shrinkage, per type",
        0,
        lambda_grid[[own]],
        u$outer$shrinkage[[own]] + const,
        n_test
      )
      add(
        tp,
        k,
        "shrinkage, common",
        0,
        lambda_grid[[common_lambda]],
        u$outer$shrinkage[[common_lambda]] + const,
        n_test
      )
      own <- which_max(u$inner$two_param)
      add(
        tp,
        k,
        "two-parameter, per type",
        pairs$alpha[[own]],
        pairs$lambda[[own]],
        two_outer[[own]] + const,
        n_test
      )
      add(
        tp,
        k,
        "two-parameter, common",
        pairs$alpha[[common_pair]],
        pairs$lambda[[common_pair]],
        two_outer[[common_pair]] + const,
        n_test
      )
      own <- which_max(u$inner$factor)
      add(
        tp,
        k,
        "factor, per type",
        NA_real_,
        ranks[[own]],
        u$outer$factor[[own]] + const,
        n_test
      )
      add(
        tp,
        k,
        "factor, shared loadings",
        NA_real_,
        ranks_shared[[shared_rank]],
        shared$outer$scores[j, shared_rank] + shared$outer$const[[j]],
        n_test
      )
    }
  }
  do.call(rbind, rows$items)
}

# ==========================================================================
# SECTION 3 · Factor models on all cells, AIC and BIC ----
# ==========================================================================

if ("fit" %in% run_stages) {
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      path <- file.path(out_dir, paste0("factors_", tm, "_", gs, ".rds"))
      if (file.exists(path)) {
        next
      }
      started <- proc.time()[["elapsed"]]
      X <- read_units(tm, gs)
      types <- names(X)
      n_genes <- ncol(X[[1L]])
      splits <- lapply(X, function(x) split_standardise(x, x))
      per_type <- parallel::mclapply(
        seq_along(X),
        logged(
          paste("fit per type", tm, gs),
          length(X),
          workers_for(length(X), n_genes),
          function(j) {
            sp <- splits[[j]]
            spectrum <- right_spectrum(sp$Z)
            ok <- usable_ranks(ranks, sp$n, length(spectrum$values), ncol(sp$Z))
            lapply(stats::setNames(ok, ok), function(q) {
              fit <- fit_factor_correlation(sp$Z, q, spectrum)
              # In-sample log-likelihood of the cells of the type.
              fit$loglik <- sp$n *
                (sp$base +
                  mean(factor_core(sp$Zt, fit$B, fit$psi)) +
                  mean(sp$const))
              B <- matrix(0, n_genes, ncol(fit$B))
              B[sp$varying, ] <- fit$B
              psi <- rep(1, n_genes)
              psi[sp$varying] <- fit$psi
              fit$B <- B
              fit$psi <- psi
              fit
            })
          }
        ),
        mc.cores = workers_for(length(X), n_genes),
        mc.preschedule = FALSE
      )
      check_tasks(per_type, paste("fit", tm, gs))
      names(per_type) <- types
      Z_list <- lapply(splits, function(sp) {
        Z <- matrix(0, sp$n, n_genes)
        Z[, sp$varying] <- sp$Z
        Z
      })
      varying_list <- lapply(splits, `[[`, "varying")
      shared_ok <- usable_ranks(
        ranks_shared,
        sum(vapply(X, nrow, integer(1))) - length(X),
        n_genes,
        min(vapply(varying_list, sum, numeric(1)))
      )
      shared <- parallel::mclapply(
        shared_ok,
        logged(
          paste("fit shared", tm, gs),
          length(shared_ok),
          workers_for(length(shared_ok), n_genes),
          function(q) {
            fit <- fit_factor_shared(Z_list, varying_list, q)
            fit$loglik <- vapply(
              seq_along(splits),
              function(j) {
                sp <- splits[[j]]
                sp$n *
                  (sp$base +
                    mean(factor_core(
                      sp$Zt,
                      fit$B[sp$varying, , drop = FALSE],
                      fit$psi[sp$varying, j]
                    )) +
                    mean(sp$const))
              },
              numeric(1)
            )
            colnames(fit$psi) <- types
            fit
          }
        ),
        mc.cores = workers_for(length(shared_ok), n_genes),
        mc.preschedule = FALSE
      )
      check_tasks(shared, paste("fit shared", tm, gs))
      names(shared) <- shared_ok
      saveRDS(
        list(
          time = tm,
          gene_set = gs,
          types = types,
          n_cells = vapply(X, nrow, integer(1)),
          n_genes = n_genes,
          per_type = per_type,
          shared = shared
        ),
        path
      )
      message(
        "Factor fits ",
        tm,
        " / ",
        gs,
        ": ",
        round(proc.time()[["elapsed"]] - started),
        " s."
      )
      rm(X, splits, Z_list, per_type, shared)
      invisible(gc())
    }
  }
}

#' AIC and BIC of the factor fits of one slot. The mean adds one
#' parameter per gene and type.
read_criteria <- function(path) {
  res <- readRDS(path)
  p <- res$n_genes
  n_types <- length(res$types)
  row <- function(tp, arm, q, loglik, k, n, fit) {
    ic <- information_criteria(loglik, k, n)
    data.frame(
      time = res$time,
      gene_set = res$gene_set,
      cell_type = tp,
      arm = arm,
      rank = q,
      loglik = loglik,
      parameters = k,
      cells = n,
      aic = ic[["aic"]],
      bic = ic[["bic"]],
      iterations = fit$iterations,
      converged = fit$converged,
      stringsAsFactors = FALSE
    )
  }
  per_type <- do.call(
    rbind,
    lapply(res$types, function(tp) {
      do.call(
        rbind,
        lapply(names(res$per_type[[tp]]), function(q) {
          fit <- res$per_type[[tp]][[q]]
          q <- as.integer(q)
          row(
            tp,
            "factor",
            q,
            fit$loglik,
            p + factor_parameter_count(p, q),
            res$n_cells[[tp]],
            fit
          )
        })
      )
    })
  )
  shared <- do.call(
    rbind,
    lapply(names(res$shared), function(q) {
      fit <- res$shared[[q]]
      q <- as.integer(q)
      row(
        "all types",
        "factor_shared",
        q,
        sum(fit$loglik),
        n_types * p + factor_parameter_count(p, q, n_types),
        sum(res$n_cells),
        fit
      )
    })
  )
  rbind(per_type, shared)
}

# ==========================================================================
# SECTION 4 · Tables ----
# ==========================================================================

if ("tables" %in% run_stages) {
  write_table <- function(x, name) {
    utils::write.csv(
      x,
      file.path(out_dir, name),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }
  stage_files <- function(stage) {
    list.files(fits_dir, paste0("^", stage, "__.*\\.rds$"), full.names = TRUE)
  }
  selections <- lapply(stage_files("select"), read_selection)
  if (length(selections) > 0L) {
    write_table(
      do.call(rbind, lapply(selections, `[[`, "curves")),
      "cv_tuning_scores.csv"
    )
    write_table(
      do.call(rbind, lapply(selections, `[[`, "selected")),
      "cv_selected_tuning.csv"
    )
  }
  nested <- do.call(rbind, lapply(stage_files("nested"), read_nested))
  if (!is.null(nested)) {
    write_table(nested, "nested_cv_by_fold.csv")
    # Paired difference to the diagonal model, per unit and rule.
    reference <- nested[nested$rule == "diagonal", ]
    key <- function(d) paste(d$time, d$gene_set, d$cell_type, d$outer_fold)
    nested$gain <- nested$score -
      reference$score[match(key(nested), key(reference))]
    groups <- split(
      nested,
      list(nested$time, nested$gene_set, nested$cell_type, nested$rule),
      drop = TRUE
    )
    write_table(
      do.call(
        rbind,
        lapply(groups, function(d) {
          data.frame(
            time = d$time[[1L]],
            gene_set = d$gene_set[[1L]],
            cell_type = d$cell_type[[1L]],
            rule = d$rule[[1L]],
            mean_score = mean(d$score),
            se_score = stats::sd(d$score) / sqrt(nrow(d)),
            mean_gain_over_diagonal = mean(d$gain),
            se_gain = stats::sd(d$gain) / sqrt(nrow(d)),
            tuning_values = paste(unique(signif(d$tuning, 3)), collapse = " "),
            alpha_values = paste(unique(d$alpha), collapse = " "),
            outer_folds = nrow(d),
            stringsAsFactors = FALSE
          )
        })
      ),
      "nested_cv_summary.csv"
    )
  }
  factor_files <- list.files(out_dir, "^factors_.*\\.rds$", full.names = TRUE)
  if (length(factor_files) > 0L) {
    write_table(
      do.call(rbind, lapply(factor_files, read_criteria)),
      "information_criteria.csv"
    )
  }
  writeLines(
    utils::capture.output(utils::sessionInfo()),
    file.path(out_dir, "session_info.txt")
  )
  message("Tables written to ", out_dir, ".")
}
