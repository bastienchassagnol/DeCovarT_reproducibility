# mkdir -p logs
# OPENBLAS_NUM_THREADS=1 nohup Rscript --no-save --no-restore \
#   scripts/04_01_stats_deconvolution.R \
#   > "logs/04_01_stats_deconvolution_$(date +%F).log" 2>&1 &
#
# Optional arguments, `key=value`, comma-separated values; each defaults
# to everything:
#   gene_sets=r2_05,r2_01,welch   times=48h,72h,96h   samples=48h_CTL_rep1
#   stages=comparators,sweep,tables
#   families=shrinkage,two_param,factor,cv,scale
#
# Deconvolution benchmark of the Suppinger bulk samples (GSE229386: 3
# standard-Wnt and 3 early-Wnt replicates at 48, 72 and 96 h), each
# against the single-cell reference of its own time, on three gene sets
# (welch, r2_01, r2_05; scripts/03_05_stats_plugin_covariance.R).
#
# DeCovarT is one model throughout: single-cell law, raw-count
# covariance, mean scale and variance scale profiled out, L-BFGS-B, bulk
# as is. Only the regularisation of the within-type covariance moves
# (build_sigma(), R/03_01_stats_covariance_estimation.R):
#   shrinkage  one level for all types, 8 values
#   two_param  pooling weight x level, 3 x 3 values
#   factor     per-type factor analysis, 6 ranks
#   cv         the tuning values chosen on held-out cells by
#              scripts/03_04_cv_covariance.R (6 rules)
#   scale      the five strategies of lowest error in the light screen
#              (scripts/04_04_stats_scale_screen.R): covariance input
#              (raw or size-normalised cells) x arm x scale model
#              (positive, linked or shared c), fitted jointly on the three
#              replicates of a protocol (fit_decovart_scales()); those
#              equal to a family above (raw cells, free scales) are not
#              repeated
# Selection without the observed composition: per time, protocol, gene
# set and family, the level with the highest profile likelihood summed
# over the three replicates, and the same with the restricted
# likelihood (decovart_likelihood_at()).
#
# Comparators, on the same gene sets and reference cells: DWLS, SCDC,
# BayesPrism, RNA-Sieve (simplex-constrained update repaired), the DeMixSC
# weights, unaligned and aligned on pseudo-bulks of the other times, MEAD
# (pseudo-individuals), ReDeconv, DECALS, robust linear regression and
# non-negative least squares.
#
# The observed proportions are the fractions of scored cells in the
# single-cell annotation at each time, not flow cytometry. Metrics use
# the standard-Wnt samples only. The variants of the second pass
# (denoised covariance, platform ratio, Newton, between-type term, Firth,
# quadratic programme) are in scripts/archives/second_pass_variants/;
# their tables are in output/deconvolution/second_pass/.
#
# Inputs:
#   data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz
#   output/feature_selection/tier_01/{scored_slots.rds, cell_counts.csv}
#   output/grn/plugin/moments_<time>_<gene_set>.rds
#   output/grn/cv/{cv_selected_tuning.csv, factors_<time>_<gene_set>.rds}
#   output/deconvolution/scale_screen_summary.csv
# Outputs, output/deconvolution/:
#   proportions_long.csv      one row per method x sample x cell type
#   sweep_diagnostics.csv     one row per DeCovarT fit
#   lambda_selection.csv      one row per time x gene set x protocol x level
#   comparator_status.csv     one row per comparator x time x gene set
#   observed_proportions.csv  scored-cell fractions per time
#   benchmark_metrics.csv     per method x standard-Wnt sample
#   benchmark_metrics_pooled.csv  per method, all standard-Wnt samples
#   profile_correlation.csv   correlation between cell-type mean profiles
#   fits/<time>/<strategy>.rds   fit store (not in Git): comparators,
#                             shrinkage, factor, two_parameter; each a list
#                             of tidy tables (estimates, diagnostics, status)

# ==========================================================================
# SECTION 0 · Dependencies, settings, paths ----
# ==========================================================================

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "02_03_stats_tail_bounds.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))
source(file.path(root, "R", "04_01_stats_deconvolution.R"))
source(file.path(root, "R", "04_02_stats_sieve_multivariate.R"))
check_blas_kernel()

time_levels <- c("48h", "72h", "96h")
gene_set_levels <- c("r2_05", "r2_01", "welch")
stage_levels <- c("comparators", "sweep", "tables")
family_levels <- c("shrinkage", "two_param", "factor", "cv", "scale")
scale_screen_top <- 5L
lambda_levels <- c(0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 0.9, 1)
two_param_levels <- expand.grid(
  alpha = c(0.25, 0.5, 0.75),
  tuning = c(0.25, 0.5, 0.75)
)
factor_ranks <- c(1L, 2L, 3L, 5L, 8L, 12L)
cv_rules <- c(
  "shrinkage, per type",
  "shrinkage, common",
  "two-parameter, per type",
  "two-parameter, common",
  "factor, per type",
  "factor, shared loadings"
)
mead_pseudo_individuals <- 4L
n_workers <- 24L
# A fit holds about seven copies of one G x G x J array at its peak.
memory_budget_gb <- 160
bayesprism_cores <- 6L
seed <- 20261110L

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
run_families <- match.arg(
  cli_value("families", family_levels),
  family_levels,
  TRUE
)
run_samples <- cli_value("samples", NULL)

bulk_path <- file.path(
  root,
  "data",
  "raw",
  "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz"
)
tier_01 <- file.path(root, "output", "feature_selection", "tier_01")
slots_path <- file.path(tier_01, "scored_slots.rds")
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
moments_dir <- file.path(root, "output", "grn", "plugin")
cv_dir <- file.path(root, "output", "grn", "cv")
out_dir <- file.path(root, "output", "deconvolution")
fits_dir <- file.path(out_dir, "fits")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
progress <- progress_path("04_01_stats_deconvolution")

# ==========================================================================
# SECTION 1 · Inputs ----
# ==========================================================================

if (!file.exists(slots_path)) {
  stop(
    "Missing ",
    slots_path,
    ". From the repository root run: dvc pull ",
    "output/feature_selection/tier_01/scored_slots.rds.dvc"
  )
}
htseq <- read_bulk_htseq(bulk_path)
design <- bulk_sample_design(colnames(htseq$counts))
observed <- observed_scored_fractions(file.path(tier_01, "cell_counts.csv"))
slots <- readRDS(slots_path)

read_moments <- function(time_level, gene_set) {
  path <- file.path(
    moments_dir,
    paste0("moments_", time_level, "_", gene_set, ".rds")
  )
  if (!file.exists(path)) {
    stop(
      "Missing ",
      path,
      ". Run scripts/03_05_stats_plugin_covariance.R first."
    )
  }
  readRDS(path)
}
samples_of <- function(time_level) {
  ids <- design$sample_id[design$time == time_level]
  if (!is.null(run_samples)) {
    ids <- intersect(ids, run_samples)
  }
  ids
}

# ==========================================================================
# SECTION 2 · Comparators ----
# ==========================================================================

if ("comparators" %in% run_stages) {
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      moments <- read_moments(tm, gs)
      genes <- moments$genes
      mu <- moments$mu
      sd_ref <- sqrt(apply(moments$sigma_sample, 3, diag))
      dimnames(sd_ref) <- dimnames(mu)
      rm(moments)
      panel <- read_panel_slot(slots, selected_dir, tm, genes = genes)
      bulk <- bulk_symbol_counts(htseq, slots$shared, genes)
      bulk <- bulk[, samples_of(tm), drop = FALSE]
      runners <- list(
        dwls = function() deconvolute_dwls(bulk, mu),
        scdc = function() deconvolute_scdc(bulk, panel$y, panel$cell_type),
        bayesprism = function() {
          deconvolute_bayesprism(
            bulk,
            panel$y,
            panel$cell_type,
            n_cores = bayesprism_cores,
            seed = seed
          )
        },
        rnasieve = function() {
          deconvolute_rnasieve(bulk, panel$y, panel$cell_type)
        },
        demixsc_wnnls = function() deconvolute_demixsc_wnnls(bulk, mu),
        mead = function() {
          deconvolute_mead(
            bulk,
            panel$y,
            panel$cell_type,
            mu,
            n_pseudo = mead_pseudo_individuals,
            seed = seed
          )
        },
        redeconv = function() {
          deconvolute_redeconv(bulk, mu, sd_ref)
        },
        rlm = function() deconvolute_rlm(bulk, mu),
        nnls = function() deconvolute_nnls(bulk, mu),
        decals = function() {
          deconvolute_decals(bulk, mu, read_moments(tm, gs)$sigma_sample)
        },
        demixsc_aligned = function() {
          pairs <- demixsc_benchmark_pairs(
            bulk_symbol_counts(htseq, slots$shared, genes),
            slots,
            selected_dir,
            observed,
            held_out_time = tm,
            seed = seed
          )
          deconvolute_demixsc_aligned(bulk, mu, pairs$bulk, pairs$pseudo)
        }
      )
      path <- fit_store_path(fits_dir, tm, "comparators")
      for (name in names(runners)) {
        # A fit stored for fewer samples, or failed, is fitted again.
        stored <- if (file.exists(path)) readRDS(path)
        fitted <- stored$estimates$sample_id[
          stored$estimates$comparator == name & stored$estimates$gene_set == gs
        ]
        ok <- any(
          stored$status$comparator == name &
            stored$status$gene_set == gs &
            stored$status$status == "ok"
        )
        if (ok && all(colnames(bulk) %in% fitted)) {
          next
        }
        started <- proc.time()[["elapsed"]]
        run <- run_with_conditions(runners[[name]])
        est <- run$value
        elapsed <- proc.time()[["elapsed"]] - started
        failed <- inherits(est, "error")
        status <- data.frame(
          comparator = name,
          time = tm,
          gene_set = gs,
          genes = length(genes),
          samples = ncol(bulk),
          status = if (failed) "error" else "ok",
          message = if (failed) {
            conditionMessage(est)
          } else {
            note <- attr(est, "note")
            paste(
              c(if (!is.null(note) && nzchar(note)) note, run$warnings),
              collapse = " | "
            )
          },
          elapsed_sec = elapsed,
          stringsAsFactors = FALSE
        )
        message(
          "Comparator ",
          name,
          " ",
          tm,
          " / ",
          gs,
          ": ",
          status$status,
          ", ",
          round(elapsed, 1),
          " s",
          if (nzchar(status$message)) paste0(" (", status$message, ")")
        )
        # A failed fit stores its status only, so a rerun tries again.
        update_fit_store(
          path,
          list(
            estimates = if (!failed) {
              data.frame(comparator = name, time = tm, gene_set = gs, est)
            },
            status = status
          ),
          key = c("comparator", "gene_set")
        )
      }
    }
  }
}

# ==========================================================================
# SECTION 3 · DeCovarT along the regularisation of the covariance ----
# ==========================================================================

# Keys of the DeCovarT fits already stored for one time and gene set:
# "<id>\r<sample_id>". A fit whose likelihood is missing counts as stored,
# as a failed checkpoint did before.
stored_sweep_keys <- function(time_level, gene_set) {
  d <- read_fit_store(
    file.path(fits_dir, time_level),
    c("shrinkage", "factor", "two_parameter"),
    "diagnostics"
  )
  d <- d[d$gene_set == gene_set, , drop = FALSE]
  paste(d$id, d$sample_id, sep = "\r")
}

# Merges the fits returned by the workers into the store, one file per
# strategy. Each task catches its own errors (tryCatch), so a failed fit
# never loses the others; only a killed R process loses the fits of the
# current time x gene set loop.
store_sweep <- function(results, time_level) {
  results <- Filter(is.list, results)
  if (length(results) == 0L) {
    return(invisible(NULL))
  }
  estimates <- bind_fill(lapply(results, `[[`, "estimates"))
  diagnostics <- bind_fill(lapply(results, `[[`, "diagnostics"))
  for (st in unique(fit_strategy(diagnostics$id))) {
    update_fit_store(
      fit_store_path(fits_dir, time_level, st),
      list(
        estimates = estimates[fit_strategy(estimates$id) == st, , drop = FALSE],
        diagnostics = diagnostics[
          fit_strategy(diagnostics$id) == st,
          ,
          drop = FALSE
        ]
      ),
      key = c("gene_set", "id", "sample_id")
    )
  }
  invisible(NULL)
}

cv_table <- function() {
  path <- file.path(cv_dir, "cv_selected_tuning.csv")
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run scripts/03_04_cv_covariance.R first.")
  }
  utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}

#' Levels of one slot: `id`, `family`, `alpha`, `tuning` (level or rank;
#' `NA` for a rule chosen per type).
sweep_levels <- function() {
  slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
  rbind(
    data.frame(
      id = sprintf("shrinkage_%.3f", lambda_levels),
      family = "shrinkage",
      alpha = 0,
      tuning = lambda_levels
    ),
    data.frame(
      id = sprintf(
        "two_param_%.2f_%.2f",
        two_param_levels$alpha,
        two_param_levels$tuning
      ),
      family = "two_param",
      two_param_levels
    ),
    data.frame(
      id = sprintf("factor_%02d", factor_ranks),
      family = "factor",
      alpha = NA_real_,
      tuning = factor_ranks
    ),
    data.frame(
      id = paste0("cv_", slug(cv_rules)),
      family = "cv",
      alpha = NA_real_,
      tuning = NA_real_
    )
  )
}

#' Covariance array of one level.
build_level <- function(level, moments, factors, cv) {
  types <- moments$types
  # Largest fitted rank of the type that does not exceed `q`.
  factor_of <- function(tp, q) {
    fits <- factors$per_type[[tp]]
    ranks <- as.integer(names(fits))
    fits[[as.character(max(ranks[ranks <= q]))]]
  }
  if (level$family == "shrinkage") {
    return(build_sigma(moments, "shrinkage", lambda = level$tuning))
  }
  if (level$family == "two_param") {
    return(build_sigma(
      moments,
      "two_param",
      lambda = level$tuning,
      alpha = level$alpha
    ))
  }
  if (level$family == "factor") {
    return(build_sigma(
      moments,
      "factor",
      factors = lapply(
        stats::setNames(types, types),
        factor_of,
        q = level$tuning
      )
    ))
  }
  rule <- cv_rules[[match(
    level$id,
    sweep_levels()$id[sweep_levels()$family == "cv"]
  )]]
  chosen <- cv[
    cv$time == moments$time_point &
      cv$gene_set == moments$gene_set &
      cv$rule == rule,
  ]
  chosen <- chosen[match(types, chosen$cell_type), ]
  if (anyNA(chosen$tuning)) {
    stop("No cross-validated value for rule '", rule, "'.")
  }
  tuning <- stats::setNames(chosen$tuning, types)
  if (rule == "factor, per type") {
    return(build_sigma(
      moments,
      "factor",
      factors = lapply(stats::setNames(types, types), function(tp) {
        factor_of(tp, tuning[[tp]])
      })
    ))
  }
  if (rule == "factor, shared loadings") {
    fit <- factors$shared[[as.character(tuning[[1L]])]]
    return(build_sigma(
      moments,
      "factor",
      factors = lapply(stats::setNames(types, types), function(tp) {
        list(B = fit$B, psi = fit$psi[, tp])
      })
    ))
  }
  if (grepl("^two-parameter", rule)) {
    return(build_sigma(
      moments,
      "two_param",
      lambda = tuning,
      alpha = stats::setNames(chosen$alpha, types)
    ))
  }
  build_sigma(moments, "shrinkage", lambda = tuning)
}

if ("sweep" %in% run_stages) {
  levels_all <- sweep_levels()
  levels_all <- levels_all[levels_all$family %in% run_families, ]
  cv <- if ("cv" %in% run_families) cv_table()
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      tasks <- merge(levels_all, data.frame(sample_id = samples_of(tm)))
      tasks$file <- paste(tasks$id, tasks$sample_id, sep = "__")
      tasks <- tasks[
        !paste(tasks$id, tasks$sample_id, sep = "\r") %in%
          stored_sweep_keys(tm, gs),
        ,
        drop = FALSE
      ]
      if (nrow(tasks) == 0L) {
        next
      }
      moments <- read_moments(tm, gs)
      genes <- moments$genes
      factors <- if (any(c("factor", "cv") %in% tasks$family)) {
        readRDS(file.path(cv_dir, paste0("factors_", tm, "_", gs, ".rds")))
      }
      bulk <- bulk_symbol_counts(htseq, slots$shared, genes)[
        genes,
        ,
        drop = FALSE
      ]
      fit_bytes <- 7 * 8 * length(genes)^2 * length(moments$types)
      workers <- max(
        1L,
        min(n_workers, nrow(tasks), floor(memory_budget_gb * 1e9 / fit_bytes))
      )
      message(
        "Sweep ",
        tm,
        " / ",
        gs,
        ": ",
        nrow(tasks),
        " fits, ",
        length(genes),
        " genes, ",
        workers,
        " workers."
      )
      log_progress(
        progress,
        "plan",
        paste("sweep", tm, gs),
        nrow(tasks),
        workers
      )
      # Forked workers share the sample moments read above.
      done <- parallel::mclapply(
        seq_len(nrow(tasks)),
        function(k) {
          task <- tasks[k, ]
          res <- tryCatch(
            {
              sigma <- build_level(task, moments, factors, cv)
              y <- bulk[, task$sample_id]
              fit <- fit_decovart_once(y, moments$mu, sigma, "lbfgsb")
              lik <- if (is.null(fit$p)) {
                list(profile = NA_real_, restricted = NA_real_)
              } else {
                decovart_likelihood_at(fit$p, y, moments$mu, sigma)
              }
              meta <- data.frame(
                time = tm,
                gene_set = gs,
                genes = length(genes),
                task[c("id", "family", "alpha", "tuning", "sample_id")],
                stringsAsFactors = FALSE
              )
              list(
                estimates = cbind(meta, fit$estimates),
                diagnostics = cbind(
                  meta,
                  fit$diagnostics,
                  profile = lik$profile,
                  restricted = lik$restricted
                )
              )
            },
            error = function(e) conditionMessage(e)
          )
          if (is.character(res)) {
            log_progress(
              progress,
              "done",
              paste("sweep", tm, gs),
              task$file,
              NA,
              "failed",
              res
            )
            return(paste0(task$file, ": ", res))
          }
          log_progress(
            progress,
            "done",
            paste("sweep", tm, gs),
            task$file,
            round(res$diagnostics$elapsed_sec, 1),
            if (is.na(res$diagnostics$loglik)) "failed" else "ok",
            res$diagnostics$message
          )
          message(
            "  ",
            task$file,
            ": ",
            round(res$diagnostics$elapsed_sec),
            " s"
          )
          res
        },
        mc.cores = workers,
        mc.preschedule = FALSE
      )
      store_sweep(done, tm)
      for (p in Filter(Negate(is.list), done)) {
        message("  FAILED fit: ", if (is.character(p)) p else format(p))
      }
      rm(moments, factors)
      invisible(gc())
    }
  }
}

#' Strategies carried from the light screen: the `scale_screen_top`
#' distinct strategies of lowest error, less those equal to a per-sample
#' family (raw cells, free scales). A `positive` strategy that never hit
#' its bound in the screen is the `free` estimator (same c, same fit), so
#' it counts as its `free` counterpart, which is faster.
scale_strategies <- function() {
  path <- file.path(out_dir, "scale_screen_summary.csv")
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run scripts/04_04_stats_scale_screen.R first.")
  }
  s <- utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  s <- s[order(s$rank), , drop = FALSE]
  same_as_free <- grepl("\\| positive$", s$strategy) & s$boundary == 0
  s$strategy[same_as_free] <- sub("positive$", "free", s$strategy[same_as_free])
  s <- s[!duplicated(s$strategy), , drop = FALSE][
    seq_len(scale_screen_top),
    ,
    drop = FALSE
  ]
  parts <- do.call(rbind, strsplit(s$strategy, " | ", fixed = TRUE))
  out <- data.frame(
    input = parts[, 1L],
    arm = parts[, 2L],
    scale = parts[, 3L],
    stringsAsFactors = FALSE
  )
  out <- out[!(out$input == "raw" & out$scale == "free"), , drop = FALSE]
  out$family <- sprintf("%s %s, %s scales", out$input, out$arm, out$scale)
  out
}

if ("sweep" %in% run_stages && "scale" %in% run_families) {
  strategies <- scale_strategies()
  message("Scale strategies: ", paste(strategies$family, collapse = "; "))
  depth <- colSums(htseq$counts)
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      levels <- do.call(
        rbind,
        lapply(seq_len(nrow(strategies)), function(k) {
          st <- strategies[k, ]
          tuning <- if (st$arm == "factor") factor_ranks else lambda_levels
          data.frame(
            st,
            tuning = tuning,
            id = paste0(gsub("[^a-z0-9]+", "_", st$family), "_", tuning),
            row.names = NULL
          )
        })
      )
      tasks <- merge(
        levels,
        data.frame(protocol = unique(design$arm[design$time == tm]))
      )
      tasks$file <- paste(tasks$id, tasks$protocol, sep = "__")
      # A protocol is fitted jointly, so its first replicate stands for it.
      first_of <- design$sample_id[design$time == tm][match(
        tasks$protocol,
        design$arm[design$time == tm]
      )]
      tasks <- tasks[
        !paste(tasks$id, first_of, sep = "\r") %in% stored_sweep_keys(tm, gs),
        ,
        drop = FALSE
      ]
      if (nrow(tasks) == 0L) {
        next
      }
      moments <- read_moments(tm, gs)
      genes <- moments$genes
      types <- moments$types
      panel <- if (any(tasks$input != "raw")) {
        read_panel_slot(slots, selected_dir, tm, genes = genes)
      }
      moments_size <- if (!is.null(panel)) {
        sample_moments(
          panel$y,
          panel$cell_type,
          types,
          library_size = panel$library_size
        )
      }
      fits <- list(
        raw = if (any(tasks$arm == "factor" & tasks$input == "raw")) {
          readRDS(file.path(
            cv_dir,
            paste0("factors_", tm, "_", gs, ".rds")
          ))$per_type
        },
        `size-normalised` = if (
          any(tasks$arm == "factor" & tasks$input != "raw")
        ) {
          fit_factor_types(
            panel$y,
            panel$cell_type,
            types,
            factor_ranks,
            panel$library_size
          )
        }
      )
      rm(panel)
      bulk <- bulk_symbol_counts(htseq, slots$shared, genes)[
        genes,
        ,
        drop = FALSE
      ]
      fit_bytes <- 7 * 8 * length(genes)^2 * length(types)
      workers <- max(
        1L,
        min(n_workers, nrow(tasks), floor(memory_budget_gb * 1e9 / fit_bytes))
      )
      slot <- paste("scale", tm, gs)
      message(
        "Scale ",
        tm,
        " / ",
        gs,
        ": ",
        nrow(tasks),
        " triplet fits, ",
        workers,
        " workers."
      )
      log_progress(progress, "plan", slot, nrow(tasks), workers)
      done <- parallel::mclapply(
        seq_len(nrow(tasks)),
        function(k) {
          task <- tasks[k, ]
          ids <- design$sample_id[
            design$time == tm & design$arm == task$protocol
          ]
          # A failed fit is stored as such (NA estimates, its message), so
          # one error does not lose the other fits of the loop.
          fit <- tryCatch(
            {
              mom <- if (task$input == "raw") moments else moments_size
              sigma <- if (task$arm == "shrinkage") {
                build_sigma(mom, "shrinkage", lambda = task$tuning)
              } else {
                f <- fits[[task$input]]
                build_sigma(
                  mom,
                  "factor",
                  factors = lapply(stats::setNames(types, types), function(tp) {
                    r <- as.integer(names(f[[tp]]))
                    f[[tp]][[as.character(max(r[r <= task$tuning]))]]
                  })
                )
              }
              fit_decovart_scales(
                bulk[, ids, drop = FALSE],
                moments$mu,
                sigma,
                task$scale,
                depth = depth[ids]
              )
            },
            error = function(e) {
              list(
                status = "failed",
                message = conditionMessage(e),
                P = NULL,
                elapsed_sec = NA_real_
              )
            }
          )
          meta <- data.frame(
            time = tm,
            gene_set = gs,
            genes = length(genes),
            id = task$id,
            family = task$family,
            alpha = NA_real_,
            tuning = task$tuning,
            sample_id = ids,
            stringsAsFactors = FALSE
          )
          ok <- !is.null(fit$P)
          estimates <- do.call(
            rbind,
            lapply(seq_along(ids), function(i) {
              data.frame(
                meta[i, ],
                cell_type = types,
                proportion = if (ok) fit$P[types, i] else NA_real_,
                lower = NA_real_,
                upper = NA_real_,
                row.names = NULL
              )
            })
          )
          diagnostics <- if (ok) {
            data.frame(
              meta,
              fit$samples[, -1L],
              convergence = NA_integer_,
              iterations = fit$iterations,
              elapsed_sec = fit$elapsed_sec,
              message = fit$status,
              profile = fit$samples$loglik,
              restricted = NA_real_,
              input = task$input,
              scale = task$scale
            )
          } else {
            data.frame(
              meta,
              loglik = NA_real_,
              elapsed_sec = fit$elapsed_sec,
              message = fit$message,
              profile = NA_real_,
              restricted = NA_real_,
              input = task$input,
              scale = task$scale
            )
          }
          log_progress(
            progress,
            "done",
            slot,
            task$file,
            round(fit$elapsed_sec, 1),
            fit$status,
            fit$message
          )
          list(estimates = estimates, diagnostics = diagnostics)
        },
        mc.cores = workers,
        mc.preschedule = FALSE
      )
      store_sweep(done, tm)
      rm(moments, moments_size, fits)
      invisible(gc())
    }
  }
}

# ==========================================================================
# SECTION 4 · Selection, tables and metrics, from every checkpoint ----
# ==========================================================================

if ("tables" %in% run_stages) {
  comparator_label <- c(
    dwls = "DWLS",
    scdc = "SCDC",
    bayesprism = "BayesPrism",
    rnasieve = "RNA-Sieve",
    demixsc_wnnls = "DeMixSC wNNLS",
    mead = "MEAD (pseudo-individuals)",
    redeconv = "ReDeconv",
    rlm = "RLM",
    nnls = "NNLS",
    decals = "DECALS",
    demixsc_aligned = "DeMixSC aligned (other times)"
  )
  family_label <- c(
    shrinkage = "shrinkage",
    two_param = "two-parameter",
    factor = "factor analysis",
    cv = "cell CV"
  )
  write_table <- function(x, name) {
    utils::write.csv(
      x,
      file.path(out_dir, name),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }

  # Rows keep the order of the former one-file-per-fit checkpoints, sorted
  # by name, so the tables below do not change with the storage.
  checkpoint_order <- function(...) {
    order(paste0(paste(..., sep = "__"), ".rds"))
  }
  comparator_status <- read_fit_store(fits_dir, "comparators", "status")
  comparator_status <- comparator_status[
    comparator_status$status == "ok",
    ,
    drop = FALSE
  ]
  comparator_status <- comparator_status[
    checkpoint_order(
      "comparator",
      comparator_status$comparator,
      comparator_status$time,
      comparator_status$gene_set
    ),
    ,
    drop = FALSE
  ]
  est <- read_fit_store(fits_dir, "comparators", "estimates")
  st <- comparator_status[
    match(
      paste(est$comparator, est$time, est$gene_set),
      paste(
        comparator_status$comparator,
        comparator_status$time,
        comparator_status$gene_set
      )
    ),
  ]
  keep <- !is.na(st$comparator)
  est <- est[keep, , drop = FALSE]
  st <- st[keep, , drop = FALSE]
  comparator_estimates <- data.frame(
    time = st$time,
    gene_set = st$gene_set,
    sample_id = est$sample_id,
    genes = st$genes,
    covariance = NA_character_,
    cell_type = est$cell_type,
    proportion = est$proportion,
    lower = est$lower,
    upper = est$upper,
    family = "comparator",
    algorithm = unname(comparator_label[st$comparator]),
    method = paste0(
      unname(comparator_label[st$comparator]),
      " | ",
      st$gene_set
    ),
    stringsAsFactors = FALSE
  )[
    checkpoint_order("comparator", st$comparator, st$time, st$gene_set),
    ,
    drop = FALSE
  ]
  rownames(comparator_status) <- rownames(comparator_estimates) <- NULL

  # A scale-family fit covers the replicates of one protocol.
  sweep_order <- function(d) {
    unit <- ifelse(
      grepl("_scales_", d$id),
      design$arm[match(d$sample_id, design$sample_id)],
      d$sample_id
    )
    checkpoint_order("sweep", d$time, d$gene_set, d$id, unit)
  }
  sweep_strategies <- c("shrinkage", "factor", "two_parameter")
  sweep_estimates <- read_fit_store(fits_dir, sweep_strategies, "estimates")
  sweep_diagnostics <- read_fit_store(fits_dir, sweep_strategies, "diagnostics")
  if (!is.null(sweep_diagnostics)) {
    sweep_estimates <- sweep_estimates[
      sweep_order(sweep_estimates),
      ,
      drop = FALSE
    ]
    sweep_diagnostics <- sweep_diagnostics[
      sweep_order(sweep_diagnostics),
      ,
      drop = FALSE
    ]
    rownames(sweep_estimates) <- rownames(sweep_diagnostics) <- NULL
  }
  label_of <- function(f) ifelse(f %in% names(family_label), family_label[f], f)
  selection <- NULL
  if (!is.null(sweep_diagnostics)) {
    sweep_diagnostics$arm <- design$arm[
      match(sweep_diagnostics$sample_id, design$sample_id)
    ]
    write_table(sweep_diagnostics, "sweep_diagnostics.csv")
    # Per time, gene set, protocol and level: likelihoods summed over the
    # replicates. A level with a failed fit is not eligible.
    by_level <- c("time", "gene_set", "arm", "family", "id", "alpha", "tuning")
    key <- interaction(sweep_diagnostics[by_level[1:5]], drop = TRUE)
    first <- !duplicated(key)
    selection <- data.frame(
      sweep_diagnostics[first, by_level],
      replicates = as.integer(table(key)[key[first]]),
      sum_profile = as.numeric(tapply(sweep_diagnostics$profile, key, sum)[key[
        first
      ]]),
      sum_restricted = as.numeric(
        tapply(sweep_diagnostics$restricted, key, sum)[key[first]]
      ),
      stringsAsFactors = FALSE
    )
    group <- interaction(
      selection[c("time", "gene_set", "arm", "family")],
      drop = TRUE
    )
    is_best <- function(x) {
      as.logical(stats::ave(x, group, FUN = function(v) {
        seq_along(v) == if (all(is.na(v))) 0L else which.max(v)
      }))
    }
    selection$selected_profile <- selection$family != "cv" &
      is_best(selection$sum_profile)
    selection$selected_restricted <- selection$family != "cv" &
      is_best(selection$sum_restricted)
    # One shrinkage level for the three times of a gene set and protocol.
    shr <- selection$family == "shrinkage"
    selection$sum_profile_all_times <- NA_real_
    across <- interaction(
      selection[shr, c("gene_set", "arm", "id")],
      drop = TRUE
    )
    complete <- tapply(selection$time[shr], across, length) ==
      length(time_levels)
    total <- tapply(selection$sum_profile[shr], across, sum)
    total[!complete] <- NA_real_
    selection$sum_profile_all_times[shr] <- as.numeric(total[across])

    at <- match(
      paste(sweep_estimates$time, sweep_estimates$gene_set, sweep_estimates$id),
      paste(selection$time, selection$gene_set, selection$id)
    )
    label <- ifelse(
      sweep_estimates$family == "cv",
      cv_rules[match(
        sweep_estimates$id,
        sweep_levels()$id[sweep_levels()$family == "cv"]
      )],
      ifelse(
        sweep_estimates$family == "two_param",
        sprintf(
          "alpha %.2f, lambda %.2f",
          sweep_estimates$alpha,
          sweep_estimates$tuning
        ),
        ifelse(
          grepl("factor", sweep_estimates$family),
          sprintf("rank %d", as.integer(sweep_estimates$tuning)),
          sprintf("lambda %.3f", sweep_estimates$tuning)
        )
      )
    )
    sweep_rows <- data.frame(
      sweep_estimates[c(
        "time",
        "gene_set",
        "sample_id",
        "genes",
        "cell_type",
        "proportion",
        "lower",
        "upper"
      )],
      covariance = label,
      family = ifelse(
        sweep_estimates$family == "cv",
        "DeCovarT, cell CV",
        "DeCovarT sweep"
      ),
      algorithm = paste("DeCovarT", label_of(sweep_estimates$family)),
      stringsAsFactors = FALSE
    )
    sweep_rows$method <- paste0(
      sweep_rows$algorithm,
      ": ",
      label,
      " | ",
      sweep_rows$gene_set
    )
    # The level chosen by the bulk likelihood, per sample's own protocol.
    chosen <- function(flag, rule) {
      arm <- design$arm[match(sweep_estimates$sample_id, design$sample_id)]
      hit <- match(
        paste(
          sweep_estimates$time,
          sweep_estimates$gene_set,
          arm,
          sweep_estimates$id
        ),
        paste(selection$time, selection$gene_set, selection$arm, selection$id)
      )
      keep <- selection[[flag]][hit]
      out <- sweep_rows[keep, ]
      out$family <- "DeCovarT, bulk likelihood"
      out$method <- paste0(out$algorithm, ", ", rule, " | ", out$gene_set)
      out
    }
    sweep_rows <- rbind(
      sweep_rows,
      chosen("selected_profile", "bulk likelihood"),
      chosen("selected_restricted", "restricted bulk likelihood")
    )
  }

  columns <- c(
    "method",
    "sample_id",
    "time",
    "cell_type",
    "proportion",
    "lower",
    "upper",
    "family",
    "algorithm",
    "gene_set",
    "covariance",
    "genes"
  )
  estimates <- rbind(
    comparator_estimates[columns],
    if (!is.null(selection)) sweep_rows[columns]
  )
  at <- match(estimates$sample_id, design$sample_id)
  estimates$arm <- design$arm[at]
  estimates$replicate <- design$replicate[at]
  estimates <- estimates[
    order(
      estimates$family,
      estimates$method,
      estimates$time,
      estimates$sample_id,
      estimates$cell_type
    ),
  ]
  write_table(estimates, "proportions_long.csv")
  write_table(observed, "observed_proportions.csv")
  write_table(comparator_status, "comparator_status.csv")

  # Metrics: standard-Wnt samples, whose protocol the reference follows.
  standard <- estimates[estimates$arm == "standard", , drop = FALSE]
  metrics <- benchmark_metrics(standard, observed)
  # A selected method can sit at a different level at each time: the
  # level is a descriptor of method x time, not of the method.
  by_time <- unique(standard[, c(
    "method",
    "time",
    "family",
    "algorithm",
    "gene_set",
    "covariance"
  )])
  descriptors <- unique(by_time[c("method", "family", "algorithm", "gene_set")])
  stopifnot(
    !anyDuplicated(by_time[c("method", "time")]),
    !anyDuplicated(descriptors$method)
  )
  metrics <- merge(metrics, by_time, by = c("method", "time"), sort = FALSE)
  metrics$replicate <- design$replicate[
    match(metrics$sample_id, design$sample_id)
  ]
  write_table(
    metrics[order(metrics$method, metrics$time, metrics$sample_id), ],
    "benchmark_metrics.csv"
  )
  if (!is.null(selection)) {
    # Error of each level, beside the likelihoods that select it. The
    # error is not used by the selection.
    level_method <- paste0(
      "DeCovarT ",
      label_of(selection$family),
      ": ",
      sweep_rows$covariance[match(
        paste(selection$time, selection$gene_set, selection$id),
        paste(
          sweep_estimates$time,
          sweep_estimates$gene_set,
          sweep_estimates$id
        )
      )],
      " | ",
      selection$gene_set
    )
    mean_of <- function(column) {
      m <- tapply(metrics[[column]], paste(metrics$method, metrics$time), mean)
      ifelse(
        selection$arm == "standard",
        as.numeric(m[paste(level_method, selection$time)]),
        NA_real_
      )
    }
    selection$rmse <- mean_of("rmse")
    selection$aitchison <- mean_of("aitchison")
    write_table(
      selection[
        order(selection$gene_set, selection$time, selection$arm, selection$id),
      ],
      "lambda_selection.csv"
    )
  }

  standard$observed <- observed$observed[match(
    paste(standard$time, standard$cell_type),
    paste(observed$time, observed$cell_type)
  )]
  pooled <- do.call(
    rbind,
    lapply(split(standard, standard$method), function(d) {
      ok <- is.finite(d$proportion)
      err <- d$proportion[ok] - d$observed[ok]
      data.frame(
        method = d$method[[1L]],
        n_points = sum(ok),
        n_missing = sum(!ok),
        n_times = length(unique(d$time[ok])),
        pearson = if (sum(ok) >= 3L && stats::sd(d$proportion[ok]) > 0) {
          stats::cor(d$proportion[ok], d$observed[ok])
        } else {
          NA_real_
        },
        rmse = if (any(ok)) sqrt(mean(err^2)) else NA_real_,
        mae = if (any(ok)) mean(abs(err)) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
  )
  # Compositions of different times have different parts: the pooled
  # Aitchison distance is the mean of the per-sample distances.
  pooled$aitchison <- as.numeric(
    tapply(metrics$aitchison, metrics$method, mean)[pooled$method]
  )
  pooled <- merge(pooled, descriptors, by = "method", sort = FALSE)
  write_table(pooled[order(pooled$rmse), ], "benchmark_metrics_pooled.csv")

  # Correlation between the mean profiles of the cell types, per time
  # and gene set (Pearson and cosine on mean counts, and Kendall).
  profile_correlation <- do.call(
    rbind,
    lapply(time_levels, function(tm) {
      do.call(
        rbind,
        lapply(gene_set_levels, function(gs) {
          mu <- read_moments(tm, gs)$mu
          grid <- expand.grid(
            type_a = colnames(mu),
            type_b = colnames(mu),
            stringsAsFactors = FALSE
          )
          at <- cbind(grid$type_a, grid$type_b)
          data.frame(
            time = tm,
            gene_set = gs,
            genes = nrow(mu),
            grid,
            pearson = stats::cor(mu)[at],
            kendall = stats::cor(mu, method = "kendall")[at],
            cosine = (crossprod(mu) / tcrossprod(sqrt(colSums(mu^2))))[at],
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
  write_table(profile_correlation, "profile_correlation.csv")

  writeLines(
    utils::capture.output(utils::sessionInfo()),
    file.path(out_dir, "session_info.txt")
  )
  message(
    "Tables written: ",
    nrow(estimates),
    " rows, ",
    length(unique(estimates$method)),
    " methods."
  )
}
