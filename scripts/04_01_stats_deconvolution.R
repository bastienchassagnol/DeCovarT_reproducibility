# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/04_01_stats_deconvolution.R \
#   > "logs/04_01_stats_deconvolution_$(date +%F).log" 2>&1 &
#
# Optional arguments, `key=value`, comma-separated values; each defaults
# to everything:
#   gene_sets=r2_05,welch   times=48h,72h,96h   samples=48h_CTL_rep1
#   stages=decovart,comparators,sieve,lambda,total,firth,tables
#   covariances=diag_raw,full_raw,full_bio   (stage sieve only)
#   ratios=none,gene       platform ratios of the stages decovart, total
#                          and firth
#   firth_budget=3600      seconds; a Firth slot whose projected time per
#                          fit is above it is recorded as not run
# A partial run writes checkpoints; `stages=tables` rebuilds the tables
# from every checkpoint on disk.
#
# Quick deconvolution benchmark of the Suppinger bulk samples
# (GSE229386: 3 standard-Wnt and 3 early-Wnt replicates at 48, 72 and
# 96 h). Each bulk sample is deconvolved against the single-cell
# reference of its own time, on two gene sets (welch, r2_05; see
# scripts/03_05_stats_plugin_covariance.R). Hyperparameters are set
# here and in R/04_01_stats_deconvolution.R; there is no configuration
# file.
#
# DeCovarT, single-cell model, profile fit, per sample and gene set:
#   4 covariance inputs  raw, denoised_bio, denoised_plugin, denoised_fixed
#                        (fit_decovart_chain())
#   2 platform ratios    none, or the bulk divided by the cross-fitted
#                        gene-wise ratio b_g of tier (iii)
#   2 optimisers         trust-region Newton, L-BFGS-B
# Comparators, on the same gene sets and reference cells: DWLS, SCDC,
# BayesPrism, RNA-Sieve, the DeMixSC weights with and without b_g, MEAD
# (pseudo-individuals), ReDeconv and robust linear regression.
#
# Diagnostic estimators (R/04_02_stats_sieve_multivariate.R), L-BFGS-B
# or quadratic programming, bulk as is unless stated:
#   sieve   multivariate RNA-Sieve in six stored stages, on three
#           covariance inputs (diagonal of the raw covariance, full raw,
#           full denoised)
#   lambda  DeCovarT on the raw covariance moved to a shrinkage level of
#           0.25, 0.5, 0.75 and 1
#   total   DeCovarT with the between-type covariance term, raw and
#           denoised covariance, with and without b_g
#   firth   Firth-penalised DeCovarT from interior starts, same grid;
#           skipped with its timing when one fit is projected over 1 h
#
# The observed proportions are the fractions of scored cells in the
# single-cell annotation at each time (observed_scored_fractions()), not
# flow cytometry. Metrics use the standard-Wnt samples only.
#
# Inputs:
#   data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz
#   output/feature_selection/tier_01/{scored_slots.rds, cell_counts.csv}
#   output/feature_selection/tier_03/platform_ratio.csv
#   output/grn/plugin/moments_<time>_<gene_set>.rds
# Outputs, output/deconvolution/:
#   proportions_long.csv      one row per method x sample x cell type
#   fit_diagnostics.csv       one row per DeCovarT fit
#   comparator_status.csv     one row per comparator x time x gene set
#   observed_proportions.csv  scored-cell fractions per time
#   benchmark_metrics.csv     per method x standard-Wnt sample
#   benchmark_metrics_pooled.csv  per method, all standard-Wnt samples
#   profile_correlation.csv   correlation between cell-type mean profiles
#   session_info.txt
#   extra_diagnostics.csv     one row per fit of the diagnostic stages
#   fits/, rnasieve/, redeconv/  checkpoints and exchange files (not in
#                             Git)

# ==========================================================================
# SECTION 0 · Dependencies, settings, paths ----
# ==========================================================================

root <- getwd()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "02_03_stats_tail_bounds.R"))
source(file.path(root, "R", "03_02_stats_regularisation_path.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))
source(file.path(root, "R", "04_01_stats_deconvolution.R"))
source(file.path(root, "R", "04_02_stats_sieve_multivariate.R"))
check_blas_kernel()

time_levels <- c("48h", "72h", "96h")
gene_set_levels <- c("r2_05", "welch")
stage_levels <- c(
  "decovart",
  "comparators",
  "sieve",
  "lambda",
  "total",
  "firth",
  "tables"
)
sieve_inputs <- c("diag_raw", "full_raw", "full_bio")
lambda_levels <- c(0.25, 0.5, 0.75, 1)
variant_covariances <- c("raw", "denoised_bio")
firth_budget_default <- 3600 # projected seconds of one Firth fit
mead_pseudo_individuals <- 4L
optimisers <- c("newton", "lbfgsb")
platform_ratios <- c("none", "gene")
max_passes <- 5L # fits of the fixed point on the bulk depth per cell
fixed_point_tol <- 0.01 # on |log(a_new / a)|
n_workers <- 16L # forked DeCovarT chains
# A chain holds about six copies of one G x G x J array at its peak.
memory_budget_gb <- 120
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
run_samples <- cli_value("samples", NULL)
platform_ratios <- match.arg(
  cli_value("ratios", platform_ratios),
  platform_ratios,
  TRUE
)
firth_budget_sec <- as.numeric(cli_value("firth_budget", firth_budget_default))
run_sieve_inputs <- match.arg(
  cli_value("covariances", sieve_inputs),
  sieve_inputs,
  TRUE
)

bulk_path <- file.path(
  root,
  "data",
  "raw",
  "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz"
)
tier_01 <- file.path(root, "output", "feature_selection", "tier_01")
slots_path <- file.path(tier_01, "scored_slots.rds")
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
ratio_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_03",
  "platform_ratio.csv"
)
moments_dir <- file.path(root, "output", "grn", "plugin")
out_dir <- file.path(root, "output", "deconvolution")
fits_dir <- file.path(out_dir, "fits")
rnasieve_dir <- file.path(out_dir, "rnasieve")
redeconv_dir <- file.path(out_dir, "redeconv")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

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
# SECTION 2 · DeCovarT: 4 covariance inputs x 2 ratios x 2 optimisers ----
# ==========================================================================

chain_file <- function(time_level, gene_set, sample_id, ratio, optimiser) {
  file.path(
    fits_dir,
    paste0(
      paste(
        "decovart",
        time_level,
        gene_set,
        sample_id,
        ratio,
        optimiser,
        sep = "__"
      ),
      ".rds"
    )
  )
}

if ("decovart" %in% run_stages) {
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      tasks <- expand.grid(
        sample_id = samples_of(tm),
        ratio = platform_ratios,
        optimiser = optimisers,
        stringsAsFactors = FALSE
      )
      tasks$file <- mapply(
        chain_file,
        tm,
        gs,
        tasks$sample_id,
        tasks$ratio,
        tasks$optimiser
      )
      tasks <- tasks[!file.exists(tasks$file), , drop = FALSE]
      if (nrow(tasks) == 0L) {
        message("DeCovarT ", tm, " / ", gs, ": every chain is on disk.")
        next
      }
      moments <- read_moments(tm, gs)
      genes <- moments$genes
      bulk <- bulk_symbol_counts(htseq, slots$shared, genes)
      ratio_b <- read_platform_ratio(ratio_path, tm, genes, method = "gene")
      chain_bytes <- 6 * 8 * length(genes)^2 * length(moments$types)
      workers <- max(
        1L,
        min(n_workers, nrow(tasks), floor(memory_budget_gb * 1e9 / chain_bytes))
      )
      message(
        "DeCovarT ",
        tm,
        " / ",
        gs,
        ": ",
        nrow(tasks),
        " chains, ",
        length(genes),
        " genes, ",
        length(moments$types),
        " types, ",
        workers,
        " workers."
      )
      # Forked workers share the moment arrays read above.
      done <- parallel::mclapply(
        seq_len(nrow(tasks)),
        function(k) {
          task <- tasks[k, ]
          started <- proc.time()[["elapsed"]]
          res <- tryCatch(
            fit_decovart_chain(
              y = bulk[, task$sample_id],
              mu = moments$mu,
              sigma_raw = moments$sigma_raw,
              sigma_bio = moments$sigma_bio,
              optimiser = task$optimiser,
              b = if (task$ratio == "gene") ratio_b else NULL,
              max_passes = max_passes,
              tol = fixed_point_tol
            ),
            error = function(e) conditionMessage(e)
          )
          elapsed <- proc.time()[["elapsed"]] - started
          if (is.character(res)) {
            return(paste0(basename(task$file), ": ", res))
          }
          res$meta <- data.frame(
            time = tm,
            gene_set = gs,
            sample_id = task$sample_id,
            platform_ratio = task$ratio,
            optimiser = task$optimiser,
            genes = length(genes),
            stringsAsFactors = FALSE
          )
          saveRDS(res, task$file)
          message(
            "  ",
            basename(task$file),
            ": ",
            round(elapsed, 1),
            " s"
          )
          NULL
        },
        mc.cores = workers,
        mc.preschedule = FALSE
      )
      problems <- Filter(Negate(is.null), done)
      for (p in problems) {
        message("  FAILED chain: ", if (is.character(p)) p else format(p))
      }
      rm(moments)
      invisible(gc())
    }
  }
}

# ==========================================================================
# SECTION 3 · Comparators ----
# ==========================================================================

comparator_file <- function(name, time_level, gene_set) {
  file.path(
    fits_dir,
    paste0(paste("comparator", name, time_level, gene_set, sep = "__"), ".rds")
  )
}

if ("comparators" %in% run_stages) {
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      moments <- read_moments(tm, gs)
      genes <- moments$genes
      mu <- moments$mu
      sd_ref <- sqrt(apply(moments$sigma_raw, 3, diag))
      dimnames(sd_ref) <- dimnames(mu)
      rm(moments)
      panel <- read_panel_slot(slots, selected_dir, tm, genes = genes)
      bulk <- bulk_symbol_counts(htseq, slots$shared, genes)
      bulk <- bulk[, samples_of(tm), drop = FALSE]
      ratio_b <- read_platform_ratio(ratio_path, tm, genes, method = "gene")
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
          deconvolute_rnasieve(
            bulk,
            panel$y,
            panel$cell_type,
            work_dir = file.path(rnasieve_dir, paste0(tm, "_", gs))
          )
        },
        demixsc_wnnls = function() deconvolute_demixsc_wnnls(bulk, mu),
        demixsc_wnnls_b = function() {
          deconvolute_demixsc_wnnls(bulk, mu, b = ratio_b)
        },
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
          deconvolute_redeconv(
            bulk,
            mu,
            sd_ref,
            work_dir = file.path(redeconv_dir, paste0(tm, "_", gs))
          )
        },
        rlm = function() deconvolute_rlm(bulk, mu)
      )
      for (name in names(runners)) {
        path <- comparator_file(name, tm, gs)
        # A checkpoint from a run on fewer samples is fitted again.
        if (
          file.exists(path) &&
            all(colnames(bulk) %in% readRDS(path)$estimates$sample_id)
        ) {
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
        # An error is kept under another name, so a rerun tries the
        # method again.
        error_path <- sub("\\.rds$", "__error.rds", path)
        if (failed) {
          saveRDS(list(status = status), error_path)
        } else {
          saveRDS(list(estimates = est, status = status), path)
          unlink(error_path)
        }
      }
    }
  }
}

# ==========================================================================
# SECTION 4 · Diagnostic estimators ----
# ==========================================================================

extra_file <- function(stage, ...) {
  file.path(fits_dir, paste0(paste("extra", stage, ..., sep = "__"), ".rds"))
}

#' Rows of an `extra` checkpoint: the columns of the DeCovarT chains,
#' plus `family` and `algorithm`.
extra_rows <- function(
  est,
  time_level,
  gene_set,
  n_genes,
  family,
  algorithm,
  covariance,
  ratio,
  optimiser
) {
  data.frame(
    time = time_level,
    gene_set = gene_set,
    sample_id = est$sample_id,
    platform_ratio = ratio,
    optimiser = optimiser,
    genes = n_genes,
    covariance = covariance,
    cell_type = est$cell_type,
    proportion = est$proportion,
    lower = est$lower %||% NA_real_,
    upper = est$upper %||% NA_real_,
    family = family,
    algorithm = algorithm,
    stringsAsFactors = FALSE
  )
}

#' Workers for a stage that factorises one G x G matrix per task and
#' holds about `copies` such matrices.
stage_workers <- function(n_tasks, n_genes, copies = 12) {
  max(
    1L,
    min(
      n_workers,
      n_tasks,
      floor(memory_budget_gb * 1e9 / (copies * 8 * n_genes^2))
    )
  )
}

# 4a · Multivariate RNA-Sieve, six stages per covariance input --------------
if ("sieve" %in% run_stages) {
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      todo <- run_sieve_inputs[
        !file.exists(extra_file("sieve", tm, gs, run_sieve_inputs))
      ]
      if (length(todo) == 0L) {
        next
      }
      moments <- read_moments(tm, gs)
      genes <- moments$genes
      bulk <- bulk_symbol_counts(htseq, slots$shared, genes)
      bulk <- bulk[, samples_of(tm), drop = FALSE]
      for (input in todo) {
        sigma <- switch(
          input,
          diag_raw = {
            v <- apply(moments$sigma_raw, 3, diag)
            dimnames(v) <- dimnames(moments$mu)
            v
          },
          full_raw = moments$sigma_raw,
          full_bio = moments$sigma_bio
        )
        started <- proc.time()[["elapsed"]]
        chain <- tryCatch(
          fit_sieve_chain(
            bulk,
            moments$mu,
            sigma,
            moments$n_cells,
            workers = stage_workers(ncol(bulk), length(genes))
          ),
          error = function(e) e
        )
        elapsed <- proc.time()[["elapsed"]] - started
        if (inherits(chain, "error")) {
          message(
            "Sieve ",
            tm,
            " / ",
            gs,
            " / ",
            input,
            ": FAILED, ",
            conditionMessage(chain)
          )
          next
        }
        diag_rows <- chain$diagnostics
        diag_rows$optimiser <- diag_rows$stage
        diag_rows$stage <- "sieve"
        saveRDS(
          list(
            estimates = do.call(
              rbind,
              lapply(
                split(chain$estimates, chain$estimates$stage),
                function(est) {
                  extra_rows(
                    est,
                    tm,
                    gs,
                    length(genes),
                    family = "sieve",
                    algorithm = "RNA-Sieve (R)",
                    covariance = input,
                    ratio = "none",
                    optimiser = est$stage[[1L]]
                  )
                }
              )
            ),
            diagnostics = cbind(
              time = tm,
              gene_set = gs,
              covariance = input,
              platform_ratio = "none",
              diag_rows,
              d_median = stats::median(chain$d),
              d_q05 = stats::quantile(chain$d, 0.05, names = FALSE),
              d_q95 = stats::quantile(chain$d, 0.95, names = FALSE),
              elapsed_sec = elapsed,
              stringsAsFactors = FALSE
            )
          ),
          extra_file("sieve", tm, gs, input)
        )
        message(
          "Sieve ",
          tm,
          " / ",
          gs,
          " / ",
          input,
          ": ",
          round(elapsed, 1),
          " s"
        )
        rm(sigma)
        invisible(gc())
      }
      rm(moments)
      invisible(gc())
    }
  }
}

# 4b · DeCovarT on the raw covariance at other shrinkage levels -------------
if ("lambda" %in% run_stages) {
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      levels_todo <- lambda_levels[
        !file.exists(extra_file("lambda", tm, gs, lambda_levels))
      ]
      if (length(levels_todo) == 0L) {
        next
      }
      moments <- read_moments(tm, gs)
      genes <- moments$genes
      bulk <- bulk_symbol_counts(htseq, slots$shared, genes)
      bulk <- bulk[genes, samples_of(tm), drop = FALSE]
      for (lam in levels_todo) {
        sigma <- rescale_shrinkage(
          moments$sigma_raw,
          moments$lambda[moments$types],
          lam
        )
        started <- proc.time()[["elapsed"]]
        fits <- parallel::mclapply(
          colnames(bulk),
          function(s) {
            fit <- fit_decovart_once(bulk[, s], moments$mu, sigma, "lbfgsb")
            fit$estimates$sample_id <- s
            fit$diagnostics$sample_id <- s
            fit
          },
          mc.cores = stage_workers(ncol(bulk), length(genes)),
          mc.preschedule = FALSE
        )
        label <- paste0("lambda_", format(lam, nsmall = 2))
        saveRDS(
          list(
            estimates = extra_rows(
              do.call(rbind, lapply(fits, `[[`, "estimates")),
              tm,
              gs,
              length(genes),
              family = "DeCovarT variant",
              algorithm = "DeCovarT shrinkage",
              covariance = label,
              ratio = "none",
              optimiser = "lbfgsb"
            ),
            diagnostics = cbind(
              time = tm,
              gene_set = gs,
              stage = "lambda",
              covariance = label,
              platform_ratio = "none",
              optimiser = "lbfgsb",
              do.call(rbind, lapply(fits, `[[`, "diagnostics")),
              stringsAsFactors = FALSE
            )
          ),
          extra_file("lambda", tm, gs, lam)
        )
        message(
          "Shrinkage ",
          lam,
          " ",
          tm,
          " / ",
          gs,
          ": ",
          round(proc.time()[["elapsed"]] - started, 1),
          " s"
        )
        rm(sigma)
        invisible(gc())
      }
      rm(moments)
      invisible(gc())
    }
  }
}

# 4c, 4d · DeCovarT with the between-type term; Firth-penalised DeCovarT ----
#' L-BFGS-B estimate of the main DeCovarT chain, for the Firth start.
chain_estimate <- function(time_level, gene_set, sample_id, ratio, covariance) {
  path <- chain_file(time_level, gene_set, sample_id, ratio, "lbfgsb")
  if (!file.exists(path)) {
    return(NULL)
  }
  est <- readRDS(path)$estimates
  est <- est[est$covariance == covariance, , drop = FALSE]
  if (nrow(est) == 0L || anyNA(est$proportion)) {
    return(NULL)
  }
  stats::setNames(est$proportion, est$cell_type)
}

for (stage in intersect(c("total", "firth"), run_stages)) {
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      tasks <- expand.grid(
        sample_id = samples_of(tm),
        covariance = variant_covariances,
        ratio = platform_ratios,
        stringsAsFactors = FALSE
      )
      tasks$file <- extra_file(
        stage,
        tm,
        gs,
        tasks$sample_id,
        tasks$covariance,
        tasks$ratio
      )
      tasks <- tasks[!file.exists(tasks$file), , drop = FALSE]
      gate_path <- extra_file(stage, tm, gs, "not_run")
      # A slot recorded as not run is tried again only with a larger
      # budget than the default.
      if (
        nrow(tasks) == 0L ||
          (file.exists(gate_path) && firth_budget_sec <= firth_budget_default)
      ) {
        next
      }
      moments <- read_moments(tm, gs)
      genes <- moments$genes
      types <- moments$types
      mu <- moments$mu
      bulk <- bulk_symbol_counts(htseq, slots$shared, genes)
      ratio_b <- read_platform_ratio(ratio_path, tm, genes, method = "gene")
      sigma_of <- function(covariance) {
        if (covariance == "raw") moments$sigma_raw else moments$sigma_bio
      }
      bulk_of <- function(task) {
        y <- as.numeric(bulk[genes, task$sample_id])
        if (task$ratio == "gene") y / as.numeric(ratio_b[genes]) else y
      }
      if (stage == "firth") {
        # One Hessian of the profile likelihood, timed. An objective
        # value costs one, a gradient 2 (J - 1); about 20 iterations
        # from each of two starts.
        task <- tasks[1L, ]
        z0 <- DeCovarT::isometric_log_ratio(
          rep(1 / length(types), length(types))
        )
        hessian_sec <- system.time(
          DeCovarT::hessian_profile_loglik_single_cell_ilr(
            z0,
            bulk_of(task),
            mu,
            sigma_of(task$covariance)
          )
        )[["elapsed"]]
        projected <- hessian_sec * (2 * (length(types) - 1) + 2) * 20 * 2
        message(
          "Firth ",
          tm,
          " / ",
          gs,
          ": one Hessian ",
          round(hessian_sec, 2),
          " s, one fit projected at ",
          round(projected),
          " s."
        )
        if (projected > firth_budget_sec) {
          saveRDS(
            list(
              estimates = NULL,
              diagnostics = data.frame(
                time = tm,
                gene_set = gs,
                stage = "firth",
                message = "not run: projected time over the budget",
                hessian_sec = hessian_sec,
                projected_sec = projected,
                stringsAsFactors = FALSE
              )
            ),
            gate_path
          )
          rm(moments)
          invisible(gc())
          next
        }
        unlink(gate_path)
      }
      done <- parallel::mclapply(
        seq_len(nrow(tasks)),
        function(k) {
          task <- tasks[k, ]
          y <- bulk_of(task)
          sigma <- sigma_of(task$covariance)
          started <- proc.time()[["elapsed"]]
          fit <- tryCatch(
            if (stage == "total") {
              fit_decovart_total(y, mu, sigma, between = TRUE)
            } else {
              start <- chain_estimate(
                tm,
                gs,
                task$sample_id,
                task$ratio,
                task$covariance
              )
              starts <- rbind(
                if (!is.null(start)) {
                  floored <- pmax(start[types], 0.001)
                  floored / sum(floored)
                },
                rep(1 / length(types), length(types))
              )
              fit_decovart_firth(y, mu, sigma, starts)
            },
            error = function(e) e
          )
          elapsed <- proc.time()[["elapsed"]] - started
          if (inherits(fit, "error")) {
            return(paste0(basename(task$file), ": ", conditionMessage(fit)))
          }
          saveRDS(
            list(
              estimates = extra_rows(
                data.frame(
                  sample_id = task$sample_id,
                  cell_type = types,
                  proportion = as.numeric(fit$p[types]),
                  stringsAsFactors = FALSE
                ),
                tm,
                gs,
                length(genes),
                family = "DeCovarT variant",
                algorithm = if (stage == "total") {
                  "DeCovarT total covariance"
                } else {
                  "DeCovarT Firth"
                },
                covariance = task$covariance,
                ratio = task$ratio,
                optimiser = if (stage == "total") "lbfgsb" else "bfgs"
              ),
              diagnostics = data.frame(
                time = tm,
                gene_set = gs,
                stage = stage,
                covariance = task$covariance,
                platform_ratio = task$ratio,
                optimiser = if (stage == "total") "lbfgsb" else "bfgs",
                sample_id = task$sample_id,
                loglik = fit$loglik,
                mean_scale = fit$mean_scale %||% NA_real_,
                variance_scale = fit$variance_scale %||% NA_real_,
                penalty = fit$penalty %||% NA_real_,
                start_kept = fit$start %||% NA_integer_,
                convergence = as.integer(fit$convergence),
                evaluations = fit$evaluations %||% fit$hessians,
                elapsed_sec = elapsed,
                stringsAsFactors = FALSE
              )
            ),
            task$file
          )
          message("  ", basename(task$file), ": ", round(elapsed, 1), " s")
          NULL
        },
        mc.cores = stage_workers(nrow(tasks), length(genes)),
        mc.preschedule = FALSE
      )
      for (p in Filter(Negate(is.null), done)) {
        message("  FAILED ", stage, ": ", if (is.character(p)) p else format(p))
      }
      rm(moments)
      invisible(gc())
    }
  }
}

# ==========================================================================
# SECTION 5 · Tables and metrics, from every checkpoint on disk ----
# ==========================================================================

if ("tables" %in% run_stages) {
  covariance_label <- c(
    raw = "raw",
    denoised_bio = "denoised (bio only)",
    denoised_plugin = "denoised + Poisson (plug-in a)",
    denoised_fixed = "denoised + Poisson (fixed-point a)"
  )
  comparator_label <- c(
    dwls = "DWLS",
    scdc = "SCDC",
    bayesprism = "BayesPrism",
    rnasieve = "RNA-Sieve",
    demixsc_wnnls = "DeMixSC wNNLS",
    demixsc_wnnls_b = "DeMixSC wNNLS",
    mead = "MEAD (pseudo-individuals)",
    redeconv = "ReDeconv",
    rlm = "RLM"
  )

  chain_files <- list.files(fits_dir, "^decovart__.*\\.rds$", full.names = TRUE)
  chains <- lapply(chain_files, readRDS)
  decovart_estimates <- do.call(
    rbind,
    lapply(chains, function(res) {
      cbind(res$meta, res$estimates, stringsAsFactors = FALSE)
    })
  )
  diagnostics <- do.call(
    rbind,
    lapply(chains, function(res) {
      cbind(res$meta, res$diagnostics, stringsAsFactors = FALSE)
    })
  )
  if (!is.null(decovart_estimates)) {
    decovart_estimates$family <- "DeCovarT"
    decovart_estimates$algorithm <- "DeCovarT"
    decovart_estimates$method <- paste0(
      "DeCovarT ",
      covariance_label[decovart_estimates$covariance],
      " | ",
      decovart_estimates$optimiser,
      " | b_g ",
      decovart_estimates$platform_ratio,
      " | ",
      decovart_estimates$gene_set
    )
  }

  comparator_files <- list.files(
    fits_dir,
    "^comparator__.*\\.rds$",
    full.names = TRUE
  )
  comparators <- lapply(comparator_files, readRDS)
  comparator_status <- do.call(rbind, lapply(comparators, `[[`, "status"))
  comparator_estimates <- do.call(
    rbind,
    lapply(comparators, function(res) {
      if (is.null(res$estimates)) {
        return(NULL)
      }
      st <- res$status
      with_ratio <- st$comparator == "demixsc_wnnls_b"
      data.frame(
        time = st$time,
        gene_set = st$gene_set,
        sample_id = res$estimates$sample_id,
        platform_ratio = if (with_ratio) "gene" else "none",
        optimiser = NA_character_,
        genes = st$genes,
        covariance = NA_character_,
        cell_type = res$estimates$cell_type,
        proportion = res$estimates$proportion,
        lower = res$estimates$lower,
        upper = res$estimates$upper,
        family = "comparator",
        algorithm = comparator_label[[st$comparator]],
        method = paste0(
          comparator_label[[st$comparator]],
          if (with_ratio) " + b_g",
          " | ",
          st$gene_set
        ),
        stringsAsFactors = FALSE
      )
    })
  )

  extras <- lapply(
    list.files(fits_dir, "^extra__.*\\.rds$", full.names = TRUE),
    readRDS
  )
  extra_estimates <- do.call(rbind, lapply(extras, `[[`, "estimates"))
  if (!is.null(extra_estimates)) {
    extra_estimates$method <- paste0(
      extra_estimates$algorithm,
      " ",
      extra_estimates$covariance,
      " | ",
      extra_estimates$optimiser,
      " | b_g ",
      extra_estimates$platform_ratio,
      " | ",
      extra_estimates$gene_set
    )
  }
  extra_diagnostics <- data.table::rbindlist(
    lapply(extras, `[[`, "diagnostics"),
    fill = TRUE
  )

  estimates <- rbind(decovart_estimates, comparator_estimates, extra_estimates)
  if (is.null(estimates)) {
    stop("No checkpoint in ", fits_dir, ".")
  }
  at <- match(estimates$sample_id, design$sample_id)
  estimates$arm <- design$arm[at]
  estimates$replicate <- design$replicate[at]
  estimates$reference_slot <- estimates$time
  estimates <- estimates[
    order(
      estimates$family,
      estimates$method,
      estimates$time,
      estimates$sample_id,
      estimates$cell_type
    ),
    c(
      "method",
      "sample_id",
      "arm",
      "time",
      "replicate",
      "reference_slot",
      "cell_type",
      "proportion",
      "lower",
      "upper",
      "family",
      "algorithm",
      "gene_set",
      "covariance",
      "platform_ratio",
      "optimiser",
      "genes"
    )
  ]
  write_table <- function(x, name) {
    utils::write.csv(
      x,
      file.path(out_dir, name),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }
  write_table(estimates, "proportions_long.csv")
  write_table(observed, "observed_proportions.csv")
  if (!is.null(diagnostics)) {
    write_table(diagnostics, "fit_diagnostics.csv")
  }
  if (!is.null(comparator_status)) {
    write_table(comparator_status, "comparator_status.csv")
  }
  if (nrow(extra_diagnostics) > 0L) {
    write_table(as.data.frame(extra_diagnostics), "extra_diagnostics.csv")
  }

  # Metrics: standard-Wnt samples, whose protocol the reference follows.
  standard <- estimates[estimates$arm == "standard", , drop = FALSE]
  metrics <- benchmark_metrics(standard, observed)
  descriptors <- unique(standard[, c(
    "method",
    "family",
    "algorithm",
    "gene_set",
    "covariance",
    "platform_ratio",
    "optimiser"
  )])
  metrics <- merge(metrics, descriptors, by = "method", sort = FALSE)
  metrics$replicate <- design$replicate[
    match(metrics$sample_id, design$sample_id)
  ]
  write_table(
    metrics[order(metrics$method, metrics$time, metrics$sample_id), ],
    "benchmark_metrics.csv"
  )

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
