# mkdir -p logs output/grn/huge
# nohup Rscript --vanilla \
#   scripts/03_01_grn_glasso.R \
#   > "logs/03_01_glasso_$(date +%F)_ebic-path.log" 2>&1 &
#
# Independent graphical lasso per scored type and time, EBIC path
# without diagonal penalisation. Run from the repository root. A
# diagnostic path: the tuning value used downstream is chosen by
# cross-validation in scripts/03_04_cv_covariance.R.
#
# The solver is `glasso::glasso(penalize.diagonal = FALSE)` on the
# correlation matrix, not `huge`, which always penalises the diagonal
# (the script was called 03_01_grn_huge.R, and the output folder keeps
# that name). Genes constant within a type are left out of that type's
# graph instead of emptying it.
#
# Input. Raw UMI counts of the scored cells (`scored_slots.rds`)
# restricted to the gated panel of each time (`selected_<time>.csv`):
# 1,417, 3,420 and 5,821 genes at 48, 72 and 96 h.
#
# Scale. The exported mean `mu` and covariance `sigma` are on DeCovarT's
# linear scale, set by `normalise`: "none" (default) for raw per-cell
# counts, "library" for the earlier `y_cg / s_c`, `s_c = L_c / L_ref`
# (`L_c` the full RNA library, `L_ref` its median over the cells of the
# time point), kept as a sensitivity variant. The support is detected on
# log1p(y / s_c) whatever `normalise` is (EBIC path on its correlation);
# the covariance is then refitted on the linear-scale correlation under
# that zero pattern (`refit_covariance_on_support()`). `path$omega` and
# `path$sigma` stay on the log1p correlation scale. Outputs written
# before these changes used the 500-gene legacy panel and are stale.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "ebic_regularisation_path.R"))
source(file.path(root, "R", "cv_covariance_selection.R"))
source(file.path(root, "R", "plot_grn_regularisation.R"))

time_levels <- c("48h", "72h", "96h")
normalise <- "none" # "none" (default) or "library" (sensitivity variant)
min_cells <- 20L
n_solves <- 40L
ebic_gamma <- 0.5
estimator <- "glasso"
rho_refit <- 1e-4

normalise <- match.arg(normalise, c("none", "library"))
slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
colour_path <- file.path(
  root,
  "data",
  "dictionaries",
  "suppinger_celltype_colours.csv"
)
out_dir <- file.path(root, "output", "grn", "huge")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
slots <- readRDS(slots_path)
colour_dict <- utils::read.csv(
  colour_path,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
type_colours <- stats::setNames(
  colour_dict$colour,
  colour_dict$celltypeannotation
)

# ==========================================================================
# SECTION 2 · EBIC path per type and time ----
# ==========================================================================

results <- lapply(time_levels, function(time_level) {
  panel <- read_panel_slot(slots, selected_dir, time_level)
  labels <- panel$cell_type
  # Full RNA library of each cell, computed before any gene subset.
  lib_size <- panel$library_size
  ref_lib <- stats::median(lib_size)
  types <- sort(unique(labels[!is.na(labels) & nzchar(labels)]))
  fits <- list()
  for (tp in types) {
    # Cells of one type at one time, treated as i.i.d. draws.
    cells <- which(labels == tp)
    if (length(cells) < min_cells) {
      message(
        "Skipping ",
        time_level,
        " / ",
        tp,
        " (n = ",
        length(cells),
        " cells)"
      )
      next
    }
    y <- panel$y[cells, , drop = FALSE]
    size_factor <- lib_size[cells] / ref_lib
    # Linear scale of the export (see `normalise`); log1p of the
    # library-divided counts only for the support.
    x_lin <- if (identical(normalise, "none")) y else y / size_factor
    x <- log1p(y / size_factor)
    n <- nrow(x)
    # A constant gene has no correlation with anything. Fit the graph on
    # the other genes of this type and record the excluded ones.
    constant <- vapply(
      seq_len(ncol(x)),
      function(j) all(x[, j] == x[1L, j]),
      logical(1)
    )
    excluded <- colnames(x)[constant]
    if (length(excluded) > 0L) {
      shown <- if (length(excluded) > 8L) {
        paste0(paste(excluded[seq_len(8L)], collapse = ", "), ", ...")
      } else {
        paste(excluded, collapse = ", ")
      }
      message(
        "  ",
        time_level,
        " / ",
        tp,
        ": ",
        length(excluded),
        " genes constant across these ",
        n,
        " cells are left out of this graph (",
        shown,
        ")."
      )
      x <- x[, !constant, drop = FALSE]
      x_lin <- x_lin[, !constant, drop = FALSE]
    }
    p <- ncol(x)
    if (p < 2L) {
      message("Skipping ", time_level, " / ", tp, " (fewer than 2 genes vary)")
      next
    }
    S <- stats::cov(x)
    sd_genes <- sqrt(diag(S))
    # Penalty, endpoints and EBIC all live on the correlation scale.
    R <- stats::cov2cor(S)
    ends <- lambda_endpoints(n, p, R)
    message(
      "Fitting glasso: ",
      time_level,
      " / ",
      tp,
      " (n = ",
      n,
      " cells, p = ",
      p,
      " genes); lambda_min = ",
      signif(ends$lambda_min, 4),
      "; lambda_max = max |R_gg'| = ",
      signif(ends$lambda_max, 4)
    )
    path <- tryCatch(
      ebic_adaptive_path(
        fit_precision = make_glasso_solver(R, n),
        n = n,
        p = p,
        S = R,
        n_solves = n_solves,
        gamma = ebic_gamma,
        lambda_min = ends$lambda_min,
        lambda_max = ends$lambda_max,
        verbose = TRUE
      ),
      error = function(e) e
    )
    if (inherits(path, "error")) {
      message(
        "Skipping ",
        time_level,
        " / ",
        tp,
        ": ",
        conditionMessage(path)
      )
      next
    }
    if (length(path$failed_lambda) > 0L) {
      message(
        "  ",
        length(path$failed_lambda),
        " penalties failed and were left out of the EBIC selection."
      )
    }
    sel <- path$selected
    # Support from the log1p fit, covariance refitted on the linear scale.
    lin <- standardise_genes(x_lin)
    refit <- refit_covariance_on_support(
      lin$R,
      adjacency = path$omega[[sel]],
      sd = lin$sd,
      rho = rho_refit,
      nobs = n
    )
    if (isTRUE(refit$failed)) {
      message(
        "Skipping ",
        time_level,
        " / ",
        tp,
        ": linear-scale refit on the EBIC support failed."
      )
      next
    }
    omega <- refit$omega
    sigma <- refit$sigma
    mu <- lin$mu
    S_lin <- stats::cov(x_lin)
    sigma_emp <- if (n > p) {
      S_lin
    } else {
      diag(diag(S_lin), p)
    }
    hell <- tryCatch(
      gaussian_hellinger(mu, sigma, mu, sigma_emp),
      error = function(e) NA_real_
    )
    # path$omega and path$sigma stay on the log1p correlation scale;
    # mu, sigma, omega and sd are on the linear scale set by `normalise`.
    fits[[tp]] <- list(
      path = path,
      genes = colnames(x),
      excluded_genes = excluded,
      scale = if (identical(normalise, "none")) {
        "raw UMI counts (no per-cell library division)"
      } else {
        "library-divided counts (full library / median)"
      },
      normalise = normalise,
      ref_lib = ref_lib,
      size_factor = size_factor,
      sd_log1p = sd_genes,
      sd = lin$sd,
      n = n,
      mu = mu,
      sigma = sigma,
      omega = omega,
      hellinger_empirical = hell
    )
  }
  fits
})
names(results) <- time_levels

out_rds <- file.path(out_dir, "huge_ebic_path_by_type_time.rds")
saveRDS(results, out_rds)
message("Wrote ", out_rds)

# ==========================================================================
# SECTION 3 · Path visualisations ----
# ==========================================================================

plot_regularisation_paths(
  results,
  file.path(out_dir, "regularisation_paths.pdf"),
  type_colours
)
plot_ebic_vs_lambda(
  results,
  file.path(out_dir, "ebic_vs_lambda.pdf"),
  type_colours
)
plot_edge_count_vs_ebic(
  results,
  file.path(out_dir, "edge_density_vs_ebic.pdf"),
  type_colours
)
plot_covariance_precision_heatmaps(
  results,
  file.path(out_dir, "covariance_precision_heatmaps.pdf")
)

# ==========================================================================
# SECTION 4 · Selected-network summaries ----
# ==========================================================================

metrics <- selected_metric_table(results, estimator)
utils::write.csv(
  metrics,
  file.path(out_dir, "selected_network_metrics.csv"),
  row.names = FALSE
)
results_linked <- lapply(results, drop_diagonal_fits)
plot_selected_funkyheatmap(
  selected_metric_table(results_linked, estimator),
  file.path(out_dir, "selected_network_funkyheatmap.pdf")
)

pair_rows <- list()
for (time_level in time_levels) {
  types <- names(results_linked[[time_level]])
  if (length(types) < 2L) {
    next
  }
  for (a in seq_along(types)) {
    for (b in seq_along(types)) {
      ia <- results[[time_level]][[types[[a]]]]
      ib <- results[[time_level]][[types[[b]]]]
      # Marginals of each Gaussian on the genes both types kept.
      shared <- intersect(ia$genes, ib$genes)
      hell <- tryCatch(
        gaussian_hellinger(
          ia$mu[shared],
          ia$sigma[shared, shared],
          ib$mu[shared],
          ib$sigma[shared, shared]
        ),
        error = function(e) NA_real_
      )
      ov <- tryCatch(
        gaussian_mixsim_overlap(
          ia$mu[shared],
          ia$sigma[shared, shared],
          ia$n,
          ib$mu[shared],
          ib$sigma[shared, shared],
          ib$n
        ),
        error = function(e) NA_real_
      )
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        time_point = time_level,
        estimator = estimator,
        type_a = types[[a]],
        type_b = types[[b]],
        n_shared_genes = length(shared),
        hellinger = hell,
        mixsim_overlap = ov,
        stringsAsFactors = FALSE
      )
    }
  }
}
pair_table <- do.call(rbind, pair_rows)
utils::write.csv(
  pair_table,
  file.path(out_dir, "pairwise_gaussian_scores.csv"),
  row.names = FALSE
)
plot_pairwise_gaussian_tiles(
  pair_table,
  file.path(out_dir, "pairwise_hellinger.pdf"),
  "hellinger",
  "Pairwise Hellinger"
)
plot_pairwise_gaussian_tiles(
  pair_table,
  file.path(out_dir, "pairwise_mixsim_overlap.pdf"),
  "mixsim_overlap",
  "Pairwise MixSim overlap"
)

message("Wrote figures under ", out_dir)
