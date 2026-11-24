# mkdir -p logs
# OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1 nohup Rscript --no-save \
#   --no-restore scripts/04_04_stats_scale_screen.R \
#   > "logs/04_04_stats_scale_screen_$(date +%F).log" 2>&1 &
#
# Optional arguments: times=72h,96h  stages=fit,tables
#
# Light screen of the scale and covariance strategies of DeCovarT on
# pseudo-bulks of known composition, before the full pipeline. A
# strategy is a covariance input x arm x scale model:
#   input  raw cells, or cells divided by their size factor before the
#          covariance (size term removed, Poisson part kept;
#          sample_moments(library_size = ))
#   arm    shrinkage, lambda in {0.25, 0.5, 0.75, 0.9}; factor analysis
#          per type, rank in {1, 2, 3, 5, 8}
#   scale  free (c_i, phi_i per replicate, the DeCovarT profile),
#          positive (c_i >= c_min), linked (phi_i = a c_i, a shared by
#          the replicates, n_i per replicate), shared_c (one c after
#          scaling the replicates to a common depth)
#          (fit_decovart_scales(), R/04_02_stats_sieve_multivariate.R)
# Within a strategy the level is chosen by the bulk likelihood summed
# over the three replicates (rule A); the composition is used only to
# score the choice.
#
# Design, per time (72 h, 96 h): the 1,000 Welch genes of highest
# deviance R^2; reference moments from 80% of the cells of each type;
# triplets of pseudo-bulks from the other 20%, compositions
# Dirichlet(50 x observed) and 300, 500 and 800 cells (truth c_i = n_i,
# phi_i = n_i, a = 1 on raw cells); two seeds; as is and multiplied by
# the tier (iii) platform ratio. Draws are without replacement unless a
# held-out pool is smaller than the draw; their share is recorded.
#
# Outputs, output/deconvolution/:
#   scale_screen_fits.csv     one row per replicate x triplet x strategy x level
#   scale_screen_summary.csv  one row per strategy: error at the chosen
#                             level, failures, scale recovery, rank
#   scale_screen_factors_<time>.rds  per-type factor fits on the
#                             training cells (80%), raw and size-normalised:
#                             per_input[[input]][[type]][[rank]] with B, psi,
#                             align (not in Git)
#   fits/<time>/scale_screen.rds  fit store (not in Git): rows (one per
#                             replicate x fit) and compositions (estimate
#                             and truth per replicate x cell type)
# Progress: logs/04_04_stats_scale_screen_<date>.progress.tsv
#   (scripts/utils/progress_report.R).

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
for (f in c(
  "02_00_utils",
  "02_01_stats_quality_control",
  "03_01_stats_covariance_estimation",
  "04_01_stats_deconvolution",
  "04_02_stats_sieve_multivariate"
)) {
  source(file.path(root, "R", paste0(f, ".R")))
}
check_blas_kernel()

time_levels <- c("72h", "96h")
cli <- commandArgs(trailingOnly = TRUE)
cli_value <- function(key, default) {
  hit <- grep(paste0("^", key, "="), cli, value = TRUE)
  if (length(hit) == 0L) {
    default
  } else {
    strsplit(sub(paste0("^", key, "="), "", hit[[1L]]), ",")[[1L]]
  }
}
run_times <- match.arg(cli_value("times", time_levels), time_levels, TRUE)
run_stages <- match.arg(
  cli_value("stages", c("fit", "tables")),
  c("fit", "tables"),
  TRUE
)

n_genes_light <- 1000L
lambda_levels <- c(0.25, 0.5, 0.75, 0.9)
factor_ranks <- c(1L, 2L, 3L, 5L, 8L)
scale_models <- c("free", "positive", "linked", "shared_c")
inputs <- c("raw", "size-normalised")
bulk_cells <- c(300L, 500L, 800L)
concentration <- 50
seeds <- 1:2
seed <- 20261117L
n_workers <- 24L

tier_01 <- file.path(root, "output", "feature_selection", "tier_01")
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
ratio_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_03",
  "platform_ratio.csv"
)
out_dir <- file.path(root, "output", "deconvolution")
fits_dir <- file.path(out_dir, "fits")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
progress <- progress_path("04_04_stats_scale_screen")
observed <- observed_scored_fractions(file.path(tier_01, "cell_counts.csv"))

# 1 · Reference, covariances and pseudo-bulks of one time ----------------------
prepare_time <- function(tm, slots) {
  selected <- utils::read.csv(file.path(
    selected_dir,
    paste0("selected_", tm, ".csv")
  ))
  genes <- selected$gene[order(-selected$deviance_r2)][seq_len(min(
    n_genes_light,
    nrow(selected)
  ))]
  panel <- read_panel_slot(slots, selected_dir, tm, genes = genes)
  types <- sort(unique(panel$cell_type[!is.na(panel$cell_type)]))
  held_out <- withr::with_seed(
    seed,
    unlist(lapply(types, function(tp) {
      rows <- which(panel$cell_type %in% tp)
      sample(rows, ceiling(0.2 * length(rows)))
    }))
  )
  train <- setdiff(which(panel$cell_type %in% types), held_out)
  moments <- list(
    raw = sample_moments(panel$y[train, ], panel$cell_type[train], types),
    `size-normalised` = sample_moments(
      panel$y[train, ],
      panel$cell_type[train],
      types,
      library_size = panel$library_size[train]
    )
  )
  # Factor fits per type and input, padded to all genes (as scripts/03_04).
  factor_fits <- lapply(stats::setNames(inputs, inputs), function(input) {
    fit_factor_types(
      panel$y[train, ],
      panel$cell_type[train],
      types,
      factor_ranks,
      library_size = if (input == "size-normalised") panel$library_size[train]
    )
  })
  # Kept for figures (loadings, alignment with the mean); a few MB.
  saveRDS(
    list(
      time = tm,
      genes = genes,
      types = types,
      train = panel$cells[train],
      per_input = factor_fits
    ),
    file.path(out_dir, paste0("scale_screen_factors_", tm, ".rds"))
  )
  covariances <- c(
    lapply(lambda_levels, function(l) {
      stats::setNames(
        lapply(inputs, function(input) {
          list(
            arm = "shrinkage",
            tuning = l,
            align = NA_real_,
            sigma = build_sigma(moments[[input]], "shrinkage", lambda = l)
          )
        }),
        paste0(sprintf("shrinkage_%.2f", l), ".", inputs)
      )
    }),
    lapply(factor_ranks, function(q) {
      stats::setNames(
        lapply(inputs, function(input) {
          chosen <- lapply(factor_fits[[input]], function(fits) {
            r <- as.integer(names(fits))
            fits[[as.character(max(r[r <= q]))]]
          })
          list(
            arm = "factor",
            tuning = q,
            align = max(vapply(chosen, `[[`, numeric(1), "align")),
            sigma = build_sigma(moments[[input]], "factor", factors = chosen)
          )
        }),
        paste0(sprintf("factor_%02d", q), ".", inputs)
      )
    })
  )
  covariances <- unlist(covariances, recursive = FALSE)
  at_time <- observed[observed$time == tm, ]
  obs <- at_time$observed[match(types, at_time$cell_type)]
  stopifnot(!anyNA(obs), abs(sum(obs) - 1) < 1e-8)
  ratio <- read_platform_ratio(ratio_path, tm, genes)
  triplets <- lapply(seeds, function(s) {
    withr::with_seed(seed + 100L * s + match(tm, time_levels), {
      draws <- lapply(bulk_cells, function(n) {
        g <- stats::rgamma(length(types), concentration * obs)
        counts <- as.integer(stats::rmultinom(1L, n, g / sum(g)))
        replaced <- 0L
        rows_list <- vector("list", length(types))
        for (j in seq_along(types)) {
          pool <- held_out[panel$cell_type[held_out] %in% types[[j]]]
          with_replacement <- counts[[j]] > length(pool)
          if (with_replacement) {
            replaced <- replaced + counts[[j]]
          }
          rows_list[[j]] <- pool[sample.int(
            length(pool),
            counts[[j]],
            replace = with_replacement
          )]
        }
        rows <- unlist(rows_list)
        list(
          y = colSums(panel$y[rows, , drop = FALSE]),
          truth = counts / n,
          depth = sum(panel$library_size[rows]),
          replaced = replaced / n
        )
      })
      Y <- vapply(draws, `[[`, numeric(length(genes)), "y")
      dimnames(Y) <- list(genes, paste0("rep", seq_along(bulk_cells)))
      list(
        Y = Y,
        truth = vapply(draws, `[[`, numeric(length(types)), "truth"),
        depth = vapply(draws, `[[`, numeric(1), "depth"),
        replaced = vapply(draws, `[[`, numeric(1), "replaced")
      )
    })
  })
  list(
    types = types,
    moments = moments,
    covariances = covariances,
    triplets = triplets,
    ratio = ratio
  )
}

# Name of a fit, as its former checkpoint file; also the row order of the tables.
task_file <- function(task) {
  paste0(
    paste(
      task$time,
      paste0("s", task$seed),
      gsub(" ", "_", task$version),
      task$cov,
      task$scale,
      sep = "__"
    ),
    ".rds"
  )
}

# 2 · Fits ----------------------------------------------------------------------
if ("fit" %in% run_stages) {
  slots <- readRDS(file.path(tier_01, "scored_slots.rds"))
  for (tm in run_times) {
    prep <- prepare_time(tm, slots)
    cov_ids <- names(prep$covariances)
    tasks <- expand.grid(
      time = tm,
      seed = seeds,
      version = c("as is", "with platform ratio"),
      cov = cov_ids,
      scale = scale_models,
      stringsAsFactors = FALSE
    )
    tasks$file <- vapply(
      seq_len(nrow(tasks)),
      function(k) task_file(tasks[k, ]),
      ""
    )
    stored <- read_fit_store(file.path(fits_dir, tm), "scale_screen", "rows")
    tasks <- tasks[!tasks$file %in% task_file(stored), , drop = FALSE]
    message(
      "Screen ",
      tm,
      ": ",
      length(prep$types),
      " types, ",
      nrow(tasks),
      " fits."
    )
    log_progress(progress, "plan", tm, nrow(tasks), n_workers)
    done <- parallel::mclapply(
      seq_len(nrow(tasks)),
      function(k) {
        task <- tasks[k, ]
        cv <- prep$covariances[[task$cov]]
        input <- sub("^.*\\.", "", task$cov)
        mom <- prep$moments[[input]]
        trip <- prep$triplets[[task$seed]]
        Y <- if (task$version == "as is") trip$Y else trip$Y * prep$ratio
        fit <- tryCatch(
          fit_decovart_scales(
            Y,
            mom$mu,
            cv$sigma,
            task$scale,
            depth = trip$depth
          ),
          error = function(e) {
            list(
              status = "failed",
              message = conditionMessage(e),
              P = NULL,
              samples = NULL,
              loglik = NA_real_,
              iterations = NA_integer_,
              elapsed_sec = NA_real_
            )
          }
        )
        rows <- data.frame(
          time = tm,
          seed = task$seed,
          version = task$version,
          input = input,
          arm = cv$arm,
          tuning = cv$tuning,
          scale = task$scale,
          sample_id = colnames(Y),
          n_true = bulk_cells,
          replaced = trip$replaced,
          max_alignment = cv$align,
          status = fit$status,
          message = fit$message,
          iterations = fit$iterations,
          elapsed_sec = fit$elapsed_sec,
          loglik_total = fit$loglik,
          stringsAsFactors = FALSE
        )
        if (!is.null(fit$P)) {
          rows <- cbind(rows, fit$samples[, -1L])
          rows$rmse <- sqrt(colMeans((fit$P - trip$truth)^2))
          rows$aitchison <- vapply(
            seq_len(ncol(Y)),
            function(i) aitchison_distance(fit$P[, i], trip$truth[, i], 0.001),
            numeric(1)
          )
        }
        rows$cov <- task$cov
        log_progress(
          progress,
          "done",
          tm,
          task$file,
          round(fit$elapsed_sec %||% NA_real_, 1),
          fit$status,
          fit$message
        )
        key <- data.frame(
          task[c("seed", "version", "cov", "scale")],
          row.names = NULL
        )
        list(
          rows = rows,
          compositions = if (!is.null(fit$P)) {
            data.frame(
              key,
              sample_id = rep(colnames(Y), each = nrow(fit$P)),
              cell_type = rownames(fit$P),
              estimate = as.vector(fit$P),
              truth = as.vector(trip$truth),
              row.names = NULL
            )
          }
        )
      },
      mc.cores = n_workers,
      mc.preschedule = FALSE
    )
    # Each fit catches its own errors (tryCatch above); only a killed R
    # process loses the fits of the current time.
    ok <- Filter(is.list, done)
    if (length(ok) > 0L) {
      update_fit_store(
        fit_store_path(fits_dir, tm, "scale_screen"),
        list(
          rows = bind_fill(lapply(ok, `[[`, "rows")),
          compositions = bind_fill(lapply(ok, `[[`, "compositions"))
        ),
        key = c("seed", "version", "cov", "scale")
      )
    }
    errs <- Filter(Negate(is.list), done)
    if (length(errs) > 0L) {
      message("Worker errors: ", length(errs), "; first: ", format(errs[[1L]]))
    }
    rm(prep)
    invisible(gc())
  }
}

# 3 · Tables --------------------------------------------------------------------
if ("tables" %in% run_stages) {
  fits <- read_fit_store(fits_dir, "scale_screen", "rows")
  fits <- fits[order(task_file(fits)), setdiff(names(fits), "cov")]
  rownames(fits) <- NULL
  utils::write.csv(
    fits,
    file.path(out_dir, "scale_screen_fits.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  per_level <- unique(fits[c(
    "time",
    "seed",
    "version",
    "input",
    "arm",
    "tuning",
    "scale",
    "status",
    "loglik_total"
  )])
  per_level$triplet <- paste(per_level$time, per_level$seed, per_level$version)
  per_level$strategy <- paste(
    per_level$input,
    per_level$arm,
    per_level$scale,
    sep = " | "
  )
  eligible <- per_level$status != "failed" & is.finite(per_level$loglik_total)
  chosen <- do.call(
    rbind,
    lapply(
      split(
        per_level[eligible, ],
        interaction(
          per_level$triplet[eligible],
          per_level$strategy[eligible],
          drop = TRUE
        )
      ),
      function(d) d[which.max(d$loglik_total), ]
    )
  )
  key <- c("time", "seed", "version", "input", "arm", "tuning", "scale")
  at_choice <- merge(chosen[c(key, "triplet", "strategy")], fits, by = key)
  by_triplet <- stats::aggregate(
    cbind(rmse, aitchison) ~ strategy + triplet + tuning,
    at_choice,
    mean
  )
  as_is <- at_choice[at_choice$version == "as is", ]
  as_is$n_ratio <- as_is$n_cells_hat / as_is$n_true
  recovery <- stats::aggregate(n_ratio ~ strategy, as_is, stats::median)
  per_level$failed <- per_level$status == "failed"
  per_level$boundary <- per_level$status == "boundary"
  rates <- stats::aggregate(cbind(failed, boundary) ~ strategy, per_level, mean)
  strategies <- sort(unique(per_level$strategy))
  n_triplets <- length(unique(per_level$triplet))
  summary <- do.call(
    rbind,
    lapply(strategies, function(s) {
      d <- by_triplet[by_triplet$strategy == s, ]
      data.frame(
        strategy = s,
        triplets_chosen = nrow(d),
        rmse = mean(d$rmse),
        aitchison = mean(d$aitchison),
        levels_chosen = paste(sort(unique(d$tuning)), collapse = ", ")
      )
    })
  )
  summary <- merge(
    merge(summary, rates, by = "strategy", all.x = TRUE),
    recovery,
    by = "strategy",
    all.x = TRUE
  )
  # Standard error of the paired difference to the lowest-error strategy.
  complete <- summary$strategy[summary$triplets_chosen == n_triplets]
  best <- complete[which.min(summary$rmse[match(complete, summary$strategy)])]
  wide <- stats::reshape(
    by_triplet[c("strategy", "triplet", "rmse")],
    idvar = "triplet",
    timevar = "strategy",
    direction = "wide"
  )
  summary$se_vs_best <- vapply(
    summary$strategy,
    function(s) {
      d <- wide[[paste0("rmse.", s)]] - wide[[paste0("rmse.", best)]]
      d <- d[is.finite(d)]
      if (length(d) < 2L) NA_real_ else stats::sd(d) / sqrt(length(d))
    },
    numeric(1)
  )
  elapsed <- stats::aggregate(
    elapsed_sec ~ strategy,
    transform(
      unique(fits[c(key, "elapsed_sec")]),
      strategy = paste(input, arm, scale, sep = " | ")
    ),
    stats::median
  )
  summary <- merge(summary, elapsed, by = "strategy", all.x = TRUE)
  # Strategies chosen on every triplet first, then by error. Ties within
  # one SE (se_vs_best) are read off the table: fewer failures, then
  # less time.
  summary <- summary[order(-summary$triplets_chosen, summary$rmse), ]
  summary$rank <- seq_len(nrow(summary))
  utils::write.csv(
    summary,
    file.path(out_dir, "scale_screen_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  print(summary, digits = 3, row.names = FALSE)
}
