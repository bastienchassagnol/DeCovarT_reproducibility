# mkdir -p logs output/grn/huge
# nohup Rscript --vanilla \
#   scripts/03_01_grn_huge.R \
#   > "logs/03_01_huge_$(date +%F)_ebic-path.log" 2>&1 &
#
# Independent graphical lasso per labelled type and time, EBIC path
# without diagonal penalisation. Run from the repository root.
#
# The file keeps its historical name; the solver is
# `glasso::glasso(penalize.diagonal = FALSE)` on the correlation matrix,
# not `huge`, which always penalises the diagonal. Genes constant within
# a type are left out of that type's graph instead of emptying it.
#
# Scale. The exported mean `mu` and covariance `sigma` are on DeCovarT's
# linear scale: size-factor-normalised counts `x_cg = y_cg / s_c`, with
# `s_c = L_c / L_ref`, `L_c` the full RNA library (metadata `log_umi`,
# log10 of the full column sum, set in 01_01) and `L_ref` its median over
# the cells of the time point. `log1p(x)` is used only to detect the
# support (EBIC path on its correlation); the covariance is then refitted
# on the linear-scale correlation under that zero pattern
# (`refit_covariance_on_support()`). `path$omega` and `path$sigma` stay on
# the log1p correlation scale. Outputs written before this change are on
# the log1p scale and are stale.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "ebic_regularisation_path.R"))
source(file.path(root, "R", "cv_covariance_selection.R"))
source(file.path(root, "R", "plot_grn_regularisation.R"))

time_levels <- c("48h", "72h", "96h")
cell_type_col <- "celltypeannotation"
min_cells <- 20L
n_solves <- 40L
ebic_gamma <- 0.5
n_genes_global <- 500L
estimator <- "glasso"
rho_refit <- 1e-4

input_rds <- file.path(
  root,
  "data",
  "intermediate",
  "equal_type_top500_by_time.rds"
)
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

objects <- readRDS(input_rds)
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
  seu <- objects[[time_level]]
  counts <- assay_matrix(seu, "RNA", "counts")
  labels <- seurat_metadata_chr(seu, cell_type_col)
  # Full RNA library, not the 500-gene column sum (`nCount_RNA` of this
  # object is recomputed on the panel by CreateSeuratObject).
  lib_size <- 10^as.numeric(seurat_metadata_chr(seu, "log_umi"))
  names(lib_size) <- colnames(seu)
  if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    stop("Metadata log_umi is missing or invalid at ", time_level, ".")
  }
  ref_lib <- stats::median(lib_size)
  types <- sort(unique(labels[!is.na(labels) & nzchar(labels)]))
  fits <- list()
  for (tp in types) {
    # Cells of one type at one time, treated as i.i.d. draws.
    cells <- unique(colnames(seu)[labels == tp])
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
    y <- t(as.matrix(counts[, cells, drop = FALSE]))
    size_factor <- lib_size[cells] / ref_lib
    # Linear scale of the convolution; log1p only for the support.
    x_lin <- y / size_factor
    x <- log1p(x_lin)
    n <- nrow(x)
    if (ncol(x) != n_genes_global) {
      stop("Expected ", n_genes_global, " genes, found ", ncol(x), ".")
    }
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
    # mu, sigma, omega and sd are on the linear (normalised-count) scale.
    fits[[tp]] <- list(
      path = path,
      genes = colnames(x),
      excluded_genes = excluded,
      scale = "size-factor-normalised counts (full library / median)",
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
