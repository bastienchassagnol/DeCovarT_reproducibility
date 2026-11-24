# Rscript --no-save --no-restore scripts/04_01_vis_deconvolution.R
#
# Figures of the deconvolution benchmark, from the tables of
# scripts/04_01_stats_deconvolution.R and scripts/03_04_cv_covariance.R.
# PDF in output/deconvolution/figures/, PNG in docs/images/deconvolution/:
#   bulk_reference_agreement, lineage_collapse_96h   mixture diagnostics
#                                         (scripts/04_02)
#   benchmark_rmse, benchmark_aitchison   comparators, and DeCovarT at the
#                                         levels chosen without the
#                                         observed composition
#   selection_curves                      held-out cell score, bulk
#                                         likelihood and error along the
#                                         shrinkage level, Welch panel
#   scale_screen                          error of the strategies of the
#                                         light screen (scripts/04_04)
#   profile_correlation_welch             correlation of the mean profiles
# The benchmark, agreement and lineage figures are also written for the
# Welch panel alone, suffix _welch (main text of docs/04).
# The figures of the second pass (scripts/archives/second_pass_variants/)
# stay in docs/images/deconvolution/ under the prefix second_pass_.

root <- here::here()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "04_01_vis_deconvolution.R"))

out_dir <- file.path(root, "output", "deconvolution")
pdf_dir <- file.path(out_dir, "figures")
png_dir <- file.path(root, "docs", "images", "deconvolution")
read_table <- function(dir, name) {
  utils::read.csv(
    file.path(dir, name),
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
}
metrics <- read_table(out_dir, "benchmark_metrics.csv")
selection <- read_table(out_dir, "lambda_selection.csv")
correlations <- read_table(out_dir, "profile_correlation.csv")
cv <- read_table(file.path(root, "output", "grn", "cv"), "cv_tuning_scores.csv")
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

# 1 · Benchmark ---------------------------------------------------------------
shown <- metrics[metrics$family != "DeCovarT sweep", ]
shown <- shown[!grepl("restricted", shown$method), ]
#' Short row label: no gene set, no "DeCovarT" prefix, no repeated rule.
short_label <- function(method) {
  x <- sub(" \\| (welch|r2_01|r2_05)$", "", method)
  x <- sub(", bulk likelihood$", "", x)
  x <- sub("^DeCovarT cell CV: ", "DeCovarT, ", x)
  x <- sub("^RLM$", "Robust regression", x)
  x <- gsub("shared_c", "shared c", x, fixed = TRUE)
  x
}
shown$label <- short_label(shown$method)
shown$group <- ifelse(
  shown$family == "DeCovarT, bulk likelihood",
  "DeCovarT, level by bulk likelihood",
  ifelse(
    shown$family == "DeCovarT, cell CV",
    "DeCovarT, level by held-out cells",
    "Other methods"
  )
)
at_96 <- shown[shown$time == "96h" & shown$gene_set == "welch", ]
# Each figure twice: all three gene panels (appendix of the book) and the
# Welch panel alone (main text), suffix _welch.
panel_versions <- list(all = c("welch", "r2_01", "r2_05"), welch = "welch")
for (version in names(panel_versions)) for (metric in c("rmse", "aitchison")) {
  best <- sort(tapply(at_96[[metric]], at_96$label, mean))[1:2]
  save_figure(
    plot_benchmark_metric(
      shown[shown$gene_set %in% panel_versions[[version]], ],
      metric,
      title = paste0(
        "Lowest error at 96 h on the Welch panel: ",
        names(best)[[1L]],
        if (metric == "rmse") " (RMSE " else " (Aitchison ",
        fmt(best[[1L]]),
        ")"
      ),
      subtitle = paste0(
        "Point: mean of the three standard-Wnt replicates; segment: their range. ",
        "No DeCovarT row uses the observed composition to choose its level."
      ),
      x_label = if (metric == "rmse") {
        "Root mean squared error of the proportions"
      } else {
        "Aitchison distance to the observed composition (zeros at 0.001)"
      }
    ),
    paste0("benchmark_", metric, if (version == "welch") "_welch"),
    pdf_dir,
    png_dir,
    width = 15,
    height = if (version == "welch") 6 else 13
  )
}

# 1b · Mixture diagnostics (scripts/04_02_stats_mixture_diagnostics.R) ---------
agreement <- read_table(out_dir, "bulk_reference_agreement.csv")
off <- tapply(agreement$share_over_4_fold, agreement$time, mean)
for (version in names(panel_versions)) save_figure(
  plot_bulk_agreement(
    agreement[agreement$gene_set %in% panel_versions[[version]], ],
    title = paste0(
      "The bulk departs most from the reference mixture at ",
      sub("h$", " h", names(which.max(off))),
      ": ",
      round(100 * max(off)),
      "% of genes more than 4-fold off"
    ),
    subtitle = paste0(
      "Standard-Wnt bulk against the reference mean at the observed composition. ",
      "Small points: replicates; large: mean."
    )
  ),
  paste0("bulk_reference_agreement", if (version == "welch") "_welch"),
  pdf_dir,
  png_dir,
  width = 12,
  height = 6
)
lineage <- read_table(out_dir, "lineage_collapsed_metrics.csv")
lineage <- lineage[
  !grepl("restricted|DeCovarT sweep", paste(lineage$method, lineage$family)),
]
lineage$label <- short_label(lineage$method)
gain <- with(
  lineage[lineage$gene_set == "welch", ],
  mean(rmse_five_parts - rmse_three_parts)
)
for (version in names(panel_versions)) save_figure(
  plot_lineage_collapse(
    lineage[lineage$gene_set %in% panel_versions[[version]], ],
    title = paste0(
      "At 96 h most error is a split inside paraxial mesoderm: summing it removes ",
      fmt(gain, 2),
      " RMSE on average (Welch)"
    ),
    subtitle = paste0(
      "Mean over the standard-Wnt replicates; caudal mesoderm, pre-somitic mesoderm ",
      "and somite front summed into one part."
    ),
    x_label = "Root mean squared error at 96 h"
  ),
  paste0("lineage_collapse_96h", if (version == "welch") "_welch"),
  pdf_dir,
  png_dir,
  width = if (version == "welch") 12 else 15,
  height = 10
)

# 2 · Criteria along the shrinkage level, Welch panel --------------------------
shr <- selection[
  selection$family == "shrinkage" &
    selection$gene_set == "welch" &
    selection$arm == "standard",
]
relative <- function(x, by) x - stats::ave(x, by, FUN = max)
bulk_curve <- data.frame(
  time = shr$time,
  tuning = shr$tuning,
  quantity = "Bulk log-likelihood,\ndifference to its maximum",
  value = relative(shr$sum_profile, shr$time),
  replicate = NA_integer_
)
cells <- cv[cv$arm == "shrinkage" & cv$gene_set == "welch", ]
cell_sum <- stats::aggregate(mean_score ~ time + tuning, cells, sum)
cell_curve <- data.frame(
  time = cell_sum$time,
  tuning = cell_sum$tuning,
  quantity = "Held-out cell log-score, summed\nover types, difference to its maximum",
  value = relative(cell_sum$mean_score, cell_sum$time),
  replicate = NA_integer_
)
sweep <- metrics[
  metrics$family == "DeCovarT sweep" &
    metrics$algorithm == "DeCovarT shrinkage" &
    metrics$gene_set == "welch",
]
sweep$tuning <- as.numeric(sub("^lambda ", "", sweep$covariance))
error_points <- data.frame(
  time = sweep$time,
  tuning = sweep$tuning,
  quantity = "Root mean squared error\n(not used by the selection)",
  value = sweep$rmse,
  replicate = sweep$replicate
)
error_mean <- stats::aggregate(
  value ~ time + tuning + quantity,
  error_points,
  mean
)
error_mean$replicate <- NA_integer_
curves <- rbind(cell_curve, bulk_curve, error_points, error_mean)
# The cell score falls by thousands at small levels; the window keeps
# the region where the three curves can be compared.
curves <- curves[curves$tuning >= 0.1, ]
curves$quantity <- factor(curves$quantity, unique(curves$quantity))
chosen <- shr[shr$selected_profile, ]
save_figure(
  plot_selection_curves(
    curves,
    title = paste0(
      "Level chosen by the bulk likelihood: ",
      paste0(
        fmt(chosen$tuning, 3),
        " at ",
        sub("h$", " h", chosen$time),
        " (RMSE ",
        fmt(chosen$rmse),
        ")",
        collapse = ", "
      )
    ),
    subtitle = paste0(
      "DeCovarT, shrinkage of the raw covariance, Welch panel, standard-Wnt ",
      "replicates. Orange: best level of each curve. Levels of 0.1 and above."
    ),
    x_label = "Shrinkage level λ (0: sample correlations; 1: variances only)"
  ),
  "selection_curves",
  pdf_dir,
  png_dir,
  width = 11,
  height = 8
)

# 3 · Light screen of the scale and covariance strategies --------------------
screen_path <- file.path(out_dir, "scale_screen_fits.csv")
if (file.exists(screen_path)) {
  screen <- read_table(out_dir, "scale_screen_fits.csv")
  screen$strategy <- paste(
    screen$input,
    screen$arm,
    paste(screen$scale, "scales"),
    sep = ", "
  )
  screen$triplet <- paste(screen$time, screen$seed, screen$version)
  levels_tried <- unique(screen[c(
    "triplet",
    "strategy",
    "tuning",
    "loglik_total",
    "status"
  )])
  levels_tried <- levels_tried[
    levels_tried$status != "failed" & is.finite(levels_tried$loglik_total),
  ]
  best <- do.call(
    rbind,
    lapply(
      split(levels_tried, paste(levels_tried$triplet, levels_tried$strategy)),
      function(d) {
        d[which.max(d$loglik_total), c("triplet", "strategy", "tuning")]
      }
    )
  )
  at_choice <- merge(best, screen, by = c("triplet", "strategy", "tuning"))
  chosen <- stats::aggregate(
    rmse ~ strategy + time + version + seed,
    at_choice,
    mean
  )
  ranking <- sort(tapply(chosen$rmse, chosen$strategy, mean))
  save_figure(
    plot_scale_screen(
      chosen,
      title = paste0(
        "Lowest error on pseudo-bulks: ",
        names(ranking)[[1L]],
        " (mean RMSE ",
        fmt(ranking[[1L]]),
        ")"
      ),
      subtitle = paste0(
        "Each strategy at the level its summed bulk likelihood chooses; ",
        "one point per triplet of replicates (2 seeds). 1,000 genes."
      ),
      x_label = "Root mean squared error of the proportions (log scale)"
    ),
    "scale_screen",
    pdf_dir,
    png_dir,
    width = 13,
    height = 8
  )
}

# 4 · Mean profiles -----------------------------------------------------------
welch <- correlations[
  correlations$gene_set == "welch" & correlations$type_a < correlations$type_b,
]
top <- welch[which.max(welch$pearson), ]
save_figure(
  plot_profile_correlation(
    correlations,
    "welch",
    title = paste0(
      "Closest pair of mean profiles: ",
      top$type_a,
      " and ",
      top$type_b,
      " at ",
      sub("h$", " h", top$time),
      " (Pearson ",
      fmt(top$pearson, 2),
      ")"
    ),
    subtitle = "Scored cell types, Welch panel, mean raw counts per cell."
  ),
  "profile_correlation_welch",
  pdf_dir,
  png_dir,
  width = 11,
  height = 8
)
message("Figures written to ", png_dir, ".")
