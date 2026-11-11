# mkdir -p logs
# Rscript --no-save --no-restore \
#   scripts/04_01_vis_deconvolution.R \
#   > "logs/04_01_vis_deconvolution_$(date +%F).log" 2>&1
#
# Figures of the quick deconvolution benchmark, from the tables of
# scripts/04_01_stats_deconvolution.R. PDF in
# output/deconvolution/figures/, PNG in docs/images/deconvolution/.
#
#   observed_vs_estimated_<gene_set>   observed against estimated
#                                      proportions, standard-Wnt samples,
#                                      one panel per method and time:
#                                      comparators and RNA-Sieve (R)
#   observed_vs_estimated_decovart_<gene_set>   the same for the
#                                      DeCovarT flavours
#   benchmark_rmse, benchmark_pearson, benchmark_aitchison
#                                      one metric for every method, gene
#                                      set, platform ratio and optimiser
#   sieve_stages_rmse, sieve_stages_aitchison
#                                      the stages of the multivariate
#                                      RNA-Sieve
#   shrinkage_sweep                    DeCovarT along the shrinkage level
#   profile_correlation_<gene_set>     correlation between the mean
#                                      profiles of the cell types
#   deconvolution_estimates_<time>     estimates of the standard and
#                                      early Wnt samples, per cell type
#
# The multivariate RNA-Sieve enters the shared figures at its stage `qp`
# (quadratic programme, reference fixed); its other stages and the
# shrinkage sweep have their own figures.
#
# Titles are built from the tables, so they follow a new run.

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "04_01_vis_deconvolution.R"))

out_dir <- file.path(root, "output", "deconvolution")
pdf_dir <- file.path(out_dir, "figures")
png_dir <- file.path(root, "docs", "images", "deconvolution")

read_table <- function(name) {
  path <- file.path(out_dir, name)
  if (!file.exists(path)) {
    stop(
      "Missing ",
      path,
      ". Run scripts/04_01_stats_deconvolution.R first."
    )
  }
  utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}
estimates <- read_table("proportions_long.csv")
observed <- read_table("observed_proportions.csv")
metrics <- read_table("benchmark_metrics.csv")
pooled <- read_table("benchmark_metrics_pooled.csv")
correlations <- read_table("profile_correlation.csv")

#' Rows drawn in the figures shared by every method: one stage of the
#' multivariate RNA-Sieve (drawn as a comparator), no shrinkage sweep.
shared_rows <- function(x) {
  keep <- x$algorithm != "DeCovarT shrinkage" &
    (x$algorithm != "RNA-Sieve (R)" | x$optimiser %in% "qp")
  x <- x[keep, , drop = FALSE]
  x$optimiser[x$algorithm == "RNA-Sieve (R)"] <- NA_character_
  x
}
estimates_shared <- shared_rows(estimates)
metrics_shared <- shared_rows(metrics)
colours <- celltype_colours(
  file.path(root, "data", "dictionaries", "suppinger_celltype_colours.csv")
)
colours <- colours[unique(observed$cell_type)]

gene_set_words <- c(
  welch = "Welch-gated genes",
  r2_05 = "genes with a deviance share of at least 5%"
)
fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

#' Lowest pooled RMSE of a subset of methods, as "<rmse> (<method>)".
#' Only methods fitted at the three times enter.
best_of <- function(tab) {
  tab <- tab[is.finite(tab$rmse) & tab$n_times == 3L, , drop = FALSE]
  if (nrow(tab) == 0L) {
    return("not available")
  }
  best <- tab[which.min(tab$rmse), ]
  name <- sub(" \\| (welch|r2_05)$", "", best$method)
  name <- gsub(" | ", ", ", sub("^DeCovarT ", "", name), fixed = TRUE)
  name <- sub("lbfgsb", "L-BFGS-B", name, fixed = TRUE)
  name <- sub("newton", "Newton", name, fixed = TRUE)
  name <- sub("b_g gene", "bulk / b_g", name, fixed = TRUE)
  name <- sub("b_g none", "bulk as is", name, fixed = TRUE)
  paste0(fmt(best$rmse), " (", name, ")")
}

# ==========================================================================
# SECTION 1 · Observed against estimated proportions ----
# ==========================================================================

scatter_groups <- list(
  list(
    suffix = "",
    families = c("comparator", "sieve"),
    note = "RNA-Sieve (R): quadratic programme, reference fixed."
  ),
  list(
    suffix = "decovart_",
    families = c("DeCovarT", "DeCovarT variant"),
    note = "DeCovarT: L-BFGS-B fits; Firth: BFGS."
  )
)
for (gs in intersect(names(gene_set_words), unique(estimates$gene_set))) {
  in_set <- pooled[pooled$gene_set == gs, , drop = FALSE]
  for (group in scatter_groups) {
    shown <- estimates_shared[
      estimates_shared$family %in% group$families,
      ,
      drop = FALSE
    ]
    n_methods <- length(unique(
      add_benchmark_labels(shown[shown$gene_set == gs, ])$method_label
    ))
    p <- plot_observed_vs_estimated(
      shown,
      observed,
      gene_set = gs,
      colours = colours,
      title = paste0(
        "Lowest RMSE over the three times: ",
        best_of(in_set[in_set$family == "comparator", ]),
        " among comparators, ",
        best_of(in_set[in_set$family == "DeCovarT", ]),
        " for DeCovarT"
      ),
      subtitle = paste0(
        "Panel: ",
        gene_set_words[[gs]],
        ". Standard-Wnt bulk replicates, one point per cell type and ",
        "replicate. Above each panel: Pearson r and RMSE. ",
        group$note
      )
    )
    save_figure(
      p,
      paste0("observed_vs_estimated_", group$suffix, gs),
      pdf_dir = pdf_dir,
      png_dir = png_dir,
      width = max(15, 1.55 * n_methods + 1.2),
      height = 7
    )
  }
}

# ==========================================================================
# SECTION 2 · Metrics of the full factorial ----
# ==========================================================================

decovart <- metrics[metrics$family == "DeCovarT", , drop = FALSE]
# Largest gap between the two optimisers over matched fits.
solver_gap <- {
  key <- c("sample_id", "gene_set", "covariance", "platform_ratio")
  wide <- merge(
    decovart[decovart$optimiser == "newton", c(key, "rmse")],
    decovart[decovart$optimiser == "lbfgsb", c(key, "rmse")],
    by = key,
    suffixes = c("_newton", "_lbfgsb")
  )
  max(abs(wide$rmse_newton - wide$rmse_lbfgsb), na.rm = TRUE)
}
ratio_effect <- {
  mean_rmse <- tapply(
    decovart$rmse,
    decovart$platform_ratio,
    mean,
    na.rm = TRUE
  )
  paste0(
    "DeCovarT mean RMSE ",
    fmt(mean_rmse[["none"]]),
    " on the bulk as is, ",
    fmt(mean_rmse[["gene"]]),
    " on the bulk divided by b_g"
  )
}
metric_subtitle <- paste0(
  "Small points: the three standard-Wnt replicates; large point: their ",
  "mean. Largest RMSE gap between Newton and L-BFGS-B: ",
  fmt(solver_gap, 3),
  "."
)
metric_titles <- c(
  rmse = ratio_effect,
  pearson = paste(
    "Pearson correlation per replicate rests on 3 cell types",
    "at 48 and 72 h, 5 at 96 h"
  ),
  aitchison = paste(
    "Aitchison distance compares within a time only:",
    "3 parts at 48 and 72 h, 5 at 96 h"
  )
)
for (metric in names(metric_titles)) {
  save_figure(
    plot_benchmark_metric(
      metrics_shared,
      metric,
      title = metric_titles[[metric]],
      subtitle = metric_subtitle
    ),
    paste0("benchmark_", metric),
    pdf_dir = pdf_dir,
    png_dir = png_dir,
    width = 11.5,
    height = 15
  )
}

# Stages of the multivariate RNA-Sieve.
sieve <- metrics[metrics$algorithm == "RNA-Sieve (R)", , drop = FALSE]
if (nrow(sieve) > 0L) {
  stage_mean <- function(stage, input, column) {
    pick <- sieve$optimiser == stage &
      sieve$covariance == input &
      sieve$gene_set == "welch" &
      sieve$time == "96h"
    if (any(pick)) fmt(mean(sieve[[column]][pick])) else "not run"
  }
  sieve_subtitle <- paste(
    "Small points: the three standard-Wnt replicates; large point: their",
    "mean. Stages are fitted from top to bottom; each starts from the",
    "previous one."
  )
  save_figure(
    plot_sieve_stages(
      sieve,
      "rmse",
      title = paste0(
        "Welch-gated genes, 96 h, full raw covariance: RMSE ",
        stage_mean("qp", "full_raw", "rmse"),
        " by QP, ",
        stage_mean("qp_means", "full_raw", "rmse"),
        " with refined means"
      ),
      subtitle = sieve_subtitle
    ),
    "sieve_stages_rmse",
    pdf_dir = pdf_dir,
    png_dir = png_dir,
    width = 11.5,
    height = 6.5
  )
  save_figure(
    plot_sieve_stages(
      sieve,
      "aitchison",
      title = paste0(
        "Welch-gated genes, 96 h, full raw covariance: distance ",
        stage_mean("qp", "full_raw", "aitchison"),
        " by QP, ",
        stage_mean("qp_means", "full_raw", "aitchison"),
        " with refined means"
      ),
      subtitle = sieve_subtitle
    ),
    "sieve_stages_aitchison",
    pdf_dir = pdf_dir,
    png_dir = png_dir,
    width = 11.5,
    height = 6.5
  )
}

# DeCovarT along the shrinkage level of the raw covariance.
if (any(metrics$algorithm == "DeCovarT shrinkage")) {
  sweep_mean <- function(level) {
    pick <- metrics$gene_set == "welch" &
      metrics$time == "96h" &
      if (is.null(level)) {
        metrics$algorithm == "DeCovarT" &
          metrics$covariance == "raw" &
          metrics$optimiser == "lbfgsb" &
          metrics$platform_ratio == "none"
      } else {
        metrics$covariance %in% paste0("lambda_", level)
      }
    if (any(pick)) fmt(mean(metrics$rmse[pick])) else "not run"
  }
  save_figure(
    plot_shrinkage_sweep(
      metrics,
      title = paste0(
        "Welch-gated genes, 96 h: RMSE ",
        sweep_mean(NULL),
        " at the closed-form \u03bb, ",
        sweep_mean("0.50"),
        " at 0.5, ",
        sweep_mean("1.00"),
        " at 1"
      ),
      subtitle = paste0(
        "DeCovarT, raw covariance, L-BFGS-B, bulk as is. Small points: the ",
        "three standard-Wnt replicates; line: their mean.\n",
        "The closed-form level differs by cell type (0.006 to 0.74); ",
        "\u03bb = 1 keeps the variances only."
      )
    ),
    "shrinkage_sweep",
    pdf_dir = pdf_dir,
    png_dir = png_dir,
    width = 10,
    height = 6.5
  )
}

# Correlation between the mean profiles of the cell types.
for (gs in intersect(names(gene_set_words), unique(correlations$gene_set))) {
  in_set <- correlations[
    correlations$gene_set == gs & correlations$type_a != correlations$type_b,
    ,
    drop = FALSE
  ]
  top <- in_set[which.max(in_set$pearson), ]
  save_figure(
    plot_profile_correlation(
      correlations,
      gs,
      title = paste0(
        "Closest pair of mean profiles: ",
        top$type_a,
        " and ",
        top$type_b,
        " at ",
        sub("h$", " h", top$time),
        " (Pearson ",
        fmt(top$pearson),
        ", Kendall ",
        fmt(top$kendall),
        ")"
      ),
      subtitle = paste0(
        "Mean raw counts per cell of each scored cell type, ",
        gene_set_words[[gs]],
        " (",
        paste(
          unique(correlations$genes[correlations$gene_set == gs]),
          collapse = ", "
        ),
        " genes at 48, 72 and 96 h). Types ordered by clustering."
      )
    ),
    paste0("profile_correlation_", gs),
    pdf_dir = pdf_dir,
    png_dir = png_dir,
    width = 12,
    height = 8
  )
}

# ==========================================================================
# SECTION 3 · Estimates per time, standard against early Wnt ----
# ==========================================================================

for (tm in unique(estimates$time)) {
  obs_tm <- observed[observed$time == tm, , drop = FALSE]
  top <- obs_tm[which.max(obs_tm$observed), ]
  p <- plot_estimates_by_arm(
    estimates_shared,
    observed,
    time_level = tm,
    title = paste0(
      "Estimated composition at ",
      sub("h$", " h", tm),
      " (observed: ",
      fmt(top$observed),
      " ",
      top$cell_type,
      ")"
    ),
    subtitle = paste(
      "One point per bulk replicate (3 per arm). DeCovarT: L-BFGS-B fits",
      "(Firth: BFGS). RNA-Sieve (R): quadratic programme, reference fixed.",
      "The observed fraction comes from standard-Wnt gastruloids."
    )
  )
  save_figure(
    p,
    paste0("deconvolution_estimates_", tm),
    pdf_dir = pdf_dir,
    png_dir = png_dir,
    width = 4.6 + 2.4 * nrow(obs_tm),
    height = 13
  )
}
