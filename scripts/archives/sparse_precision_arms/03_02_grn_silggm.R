# mkdir -p logs output/grn/silggm
# nohup Rscript --vanilla \
#   scripts/03_02_grn_silggm.R \
#   > "logs/03_02_silggm_$(date +%F)_dsgl.log" 2>&1 &
#
# Independent SILGGM fit per scored type and time. Nonparanormal
# shrinkage, de-sparsified graphical lasso, global FDR. The R object is
# the only export (no Cytoscape / CSV dump). Run from the repository
# root. This script stops at the support; the linear-scale refit under
# that support, and its cross-validated FDR level, are in
# `fit_cov_silggm_refit()` (scripts/03_04_cv_covariance.R).
#
# Input. Raw UMI counts of the scored cells (`scored_slots.rds`)
# restricted to the gated panel of each time (`selected_<time>.csv`):
# 1,417, 3,420 and 5,821 genes at 48, 72 and 96 h. With 5,821 genes at
# 96 h, D-S_GL tests about 1.7e7 gene pairs per type; no cap is imposed.
#
# `normalise`: "none" (default) passes the raw counts to the
# nonparanormal map; "library" passes `y_cg / s_c`, `s_c = L_c / median
# L` (sensitivity variant). The ranks, hence the support, depend on that
# choice only through the cell-to-cell size variation.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))

time_levels <- c("48h", "72h", "96h")
normalise <- "none" # "none" (default) or "library" (sensitivity variant)
min_cells <- 20L
npn_fun <- "shrinkage"
silggm_method <- "D-S_GL"
global_flag <- TRUE
alpha_level <- 0.05
seed <- 1L

normalise <- match.arg(normalise, c("none", "library"))
slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
out_dir <- file.path(root, "output", "grn", "silggm")
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

# ==========================================================================
# SECTION 2 · Nonparanormal + SILGGM per type and time ----
# ==========================================================================

withr::with_seed(seed, {
  results <- lapply(time_levels, function(time_level) {
    panel <- read_panel_slot(slots, selected_dir, time_level)
    labels <- panel$cell_type
    ref_lib <- stats::median(panel$library_size)
    types <- sort(unique(labels[!is.na(labels) & nzchar(labels)]))
    fits <- list()
    for (tp in types) {
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
      x <- panel$y[cells, , drop = FALSE]
      if (identical(normalise, "library")) {
        x <- x / (panel$library_size[cells] / ref_lib)
      }
      x_npn <- huge::huge.npn(x, npn.func = npn_fun, verbose = FALSE)
      message("Fitting SILGGM: ", time_level, " / ", tp)
      fits[[tp]] <- SILGGM::SILGGM(
        x_npn,
        method = silggm_method,
        global = global_flag,
        alpha = alpha_level,
        cytoscape_format = FALSE,
        csv_save = FALSE
      )
      fits[[tp]]$genes <- colnames(x)
      fits[[tp]]$n <- nrow(x)
      fits[[tp]]$normalise <- normalise
    }
    fits
  })
  names(results) <- time_levels
})

out_rds <- file.path(out_dir, "silggm_dsgl_by_type_time.rds")
saveRDS(results, out_rds)
message("Wrote ", out_rds)
