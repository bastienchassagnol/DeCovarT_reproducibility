# mkdir -p logs output/grn/huge
# nohup Rscript --vanilla \
#   scripts/03_01_grn_huge.R \
#   > "logs/03_01_huge_$(date +%F)_ebic-path.log" 2>&1 &
#
# Independent graphical lasso per labelled type and time, EBIC path
# without diagonal penalisation. Run from the repository root.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "utils_general.R"))
source(file.path(root, "R", "ebic_regularisation_path.R"))
source(file.path(root, "R", "plot_grn_regularisation.R"))

time_levels <- c("48h", "72h", "96h")
cell_type_col <- "celltypeannotation"
min_cells <- 20L
n_solves <- 40L
ebic_gamma <- 0.5
n_genes_global <- 500L
seed <- 1L
estimator <- "huge"

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
    x <- t(as.matrix(counts[, cells, drop = FALSE]))
    x <- log1p(x)
    n <- nrow(x)
    p <- ncol(x)
    if (p != n_genes_global) {
      stop("Expected ", n_genes_global, " genes, found ", p, ".")
    }
    S <- stats::cov(x)
    ends <- lambda_endpoints(n, p, S)
    message(
      "Fitting huge glasso: ",
      time_level,
      " / ",
      tp,
      " (n = ",
      n,
      " cells); lambda_min = ",
      signif(ends$lambda_min, 4),
      "; lambda_max = ",
      signif(ends$lambda_max, 4)
    )
    path <- ebic_adaptive_path(
      fit_precision = make_huge_solver(x),
      n = n,
      p = p,
      S = S,
      n_solves = n_solves,
      gamma = ebic_gamma,
      lambda_min = ends$lambda_min,
      lambda_max = ends$lambda_max,
      verbose = TRUE
    )
    sel <- path$selected
    omega <- path$omega[[sel]]
    mu <- colMeans(x)
    sigma <- tryCatch(
      solve(omega),
      error = function(e) {
        omega_r <- omega + diag(1e-4, p)
        solve(omega_r)
      }
    )
    sigma_emp <- if (n > p) {
      S
    } else {
      diag(diag(S), p)
    }
    hell <- tryCatch(
      gaussian_hellinger(mu, sigma, mu, sigma_emp),
      error = function(e) NA_real_
    )
    fits[[tp]] <- list(
      path = path,
      genes = colnames(x),
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
  type_colours,
  seed = seed
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
plot_priority_score(
  results,
  file.path(out_dir, "priority_score_vs_lambda.pdf"),
  type_colours
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
plot_selected_funkyheatmap(
  metrics,
  file.path(out_dir, "selected_network_funkyheatmap.pdf")
)

pair_rows <- list()
for (time_level in time_levels) {
  types <- names(results[[time_level]])
  if (length(types) < 2L) {
    next
  }
  for (a in seq_along(types)) {
    for (b in seq_along(types)) {
      ia <- results[[time_level]][[types[[a]]]]
      ib <- results[[time_level]][[types[[b]]]]
      hell <- tryCatch(
        gaussian_hellinger(ia$mu, ia$sigma, ib$mu, ib$sigma),
        error = function(e) NA_real_
      )
      ov <- tryCatch(
        gaussian_mixsim_overlap(
          ia$mu,
          ia$sigma,
          ia$n,
          ib$mu,
          ib$sigma,
          ib$n
        ),
        error = function(e) NA_real_
      )
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        time_point = time_level,
        estimator = estimator,
        type_a = types[[a]],
        type_b = types[[b]],
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
