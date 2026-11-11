# mkdir -p logs
# Rscript --no-save --no-restore \
#   scripts/04_02_stats_likelihood_terms.R \
#   > "logs/04_02_stats_likelihood_terms_$(date +%F).log" 2>&1
#
# Optional arguments, `key=value`: gene_sets=r2_05,welch  times=48h,72h,96h
#
# Where does the DeCovarT profile log-likelihood prefer its estimate to
# the observed composition? For every standard-Wnt bulk sample, gene set
# and covariance input (raw, denoised), the profile log-likelihood of
# `DeCovarT::deconvolute_single_cell()` is evaluated at three
# compositions and split into its two terms:
#
#   loglik = -(G/2) log(phi_hat) - (1/2) log det S_W(p) - G/2
#            `quadratic_term`       `log_det_term`
#
#   observed   fraction of scored cells at that time
#   decovart   DeCovarT L-BFGS-B estimate, bulk as is, same covariance
#   sieve_qp   quadratic-programme estimate of the multivariate
#              RNA-Sieve, same covariance, stage `qp_within`
#
# Compositions are floored at 1e-6 and rescaled, so the covariance stays
# positive definite on a face. No fit is run here.
#
# Inputs:
#   output/deconvolution/{proportions_long.csv, observed_proportions.csv}
#   output/grn/plugin/moments_<time>_<gene_set>.rds
#   data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz
#   output/feature_selection/tier_01/scored_slots.rds
# Output:
#   output/deconvolution/likelihood_terms.csv

root <- getwd()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "04_01_stats_deconvolution.R"))
source(file.path(root, "R", "04_02_stats_sieve_multivariate.R"))
check_blas_kernel()

cli <- commandArgs(trailingOnly = TRUE)
cli_value <- function(key, default) {
  hit <- grep(paste0("^", key, "="), cli, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  strsplit(sub(paste0("^", key, "="), "", hit[[1L]]), ",", fixed = TRUE)[[1L]]
}
run_times <- cli_value("times", c("48h", "72h", "96h"))
run_gene_sets <- cli_value("gene_sets", c("r2_05", "welch"))
n_workers <- 9L
floor_p <- 1e-6

out_dir <- file.path(root, "output", "deconvolution")
estimates <- utils::read.csv(
  file.path(out_dir, "proportions_long.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
observed <- utils::read.csv(
  file.path(out_dir, "observed_proportions.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
htseq <- read_bulk_htseq(
  file.path(root, "data", "raw", "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz")
)
slots <- readRDS(
  file.path(root, "output", "feature_selection", "tier_01", "scored_slots.rds")
)
standard <- estimates[
  estimates$arm == "standard" & estimates$platform_ratio == "none",
  ,
  drop = FALSE
]

#' Composition of one method for one sample, in the order of `types`.
composition_of <- function(rows, types) {
  if (nrow(rows) != length(types)) {
    return(NULL)
  }
  p <- rows$proportion[match(types, rows$cell_type)]
  if (anyNA(p)) NULL else p
}

sieve_input <- c(raw = "full_raw", denoised_bio = "full_bio")
terms <- list()
for (tm in run_times) {
  for (gs in run_gene_sets) {
    moments <- readRDS(
      file.path(
        root,
        "output",
        "grn",
        "plugin",
        paste0("moments_", tm, "_", gs, ".rds")
      )
    )
    types <- moments$types
    genes <- moments$genes
    bulk <- bulk_symbol_counts(htseq, slots$shared, genes)
    in_slot <- standard[standard$time == tm & standard$gene_set == gs, ]
    obs <- observed[observed$time == tm, ]
    tasks <- expand.grid(
      sample_id = unique(in_slot$sample_id),
      covariance = names(sieve_input),
      composition = c("observed", "decovart", "sieve_qp"),
      stringsAsFactors = FALSE
    )
    rows <- parallel::mclapply(
      seq_len(nrow(tasks)),
      function(k) {
        task <- tasks[k, ]
        mine <- in_slot[in_slot$sample_id == task$sample_id, ]
        p <- switch(
          task$composition,
          observed = obs$observed[match(types, obs$cell_type)],
          decovart = composition_of(
            mine[
              mine$algorithm == "DeCovarT" &
                mine$covariance == task$covariance &
                mine$optimiser == "lbfgsb",
            ],
            types
          ),
          sieve_qp = composition_of(
            mine[
              mine$algorithm == "RNA-Sieve (R)" &
                mine$covariance == sieve_input[[task$covariance]] &
                mine$optimiser == "qp_within",
            ],
            types
          )
        )
        if (is.null(p)) {
          return(NULL)
        }
        p <- pmax(p, floor_p)
        p <- p / sum(p)
        sigma <- if (task$covariance == "raw") {
          moments$sigma_raw
        } else {
          moments$sigma_bio
        }
        at <- tryCatch(
          decovart_profile_total(
            p,
            as.numeric(bulk[genes, task$sample_id]),
            moments$mu,
            sigma,
            between = FALSE
          ),
          error = function(e) NULL
        )
        if (is.null(at)) {
          return(NULL)
        }
        data.frame(
          time = tm,
          gene_set = gs,
          genes = length(genes),
          sample_id = task$sample_id,
          covariance = task$covariance,
          composition = task$composition,
          loglik = at$value,
          quadratic_term = -0.5 * length(genes) * log(at$variance_scale),
          log_det_term = -0.5 * at$log_det,
          mean_scale = at$mean_scale,
          variance_scale = at$variance_scale,
          stringsAsFactors = FALSE
        )
      },
      mc.cores = min(n_workers, nrow(tasks)),
      mc.preschedule = FALSE
    )
    terms <- c(terms, rows)
    message("Likelihood terms ", tm, " / ", gs, ": ", nrow(tasks), " points.")
    rm(moments)
    invisible(gc())
  }
}
terms <- do.call(rbind, Filter(is.data.frame, terms))
utils::write.csv(
  terms,
  file.path(out_dir, "likelihood_terms.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
message("Likelihood terms written: ", nrow(terms), " rows.")
