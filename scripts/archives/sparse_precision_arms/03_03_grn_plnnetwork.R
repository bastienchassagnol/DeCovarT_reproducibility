# mkdir -p logs output/grn/plnnetwork
# nohup Rscript --vanilla \
#   scripts/03_03_grn_plnnetwork.R \
#   > "logs/03_03_plnnetwork_$(date +%F)_offset.log" 2>&1 &
#
# Independent PLNnetwork fit per scored type and time, with a
# library-size offset. Run from the repository root.
#
# Input. Raw UMI counts of the scored cells (`scored_slots.rds`)
# restricted to the gated panel of each time (`selected_<time>.csv`):
# 1,417, 3,420 and 5,821 genes at 48, 72 and 96 h. With 5,821 genes at
# 96 h each fit updates a 5,821 x 5,821 latent covariance; no cap is
# imposed.
#
# The offset is log of the full RNA library L_c (`slot$library_size`,
# the column sum over all deposited genes, computed before any gene
# subset), not the row sum of the panel genes: a panel sum moves with
# the panel's own expression and absorbs part of the cell-type signal.
# `prepare_data(offset = "none")` stops PLNmodels from adding its TSS
# (panel-sum) offset. There is no `normalise` switch: the offset is the
# model's own treatment of the size factor, and what is put back at
# export is decided by `pln_linear_moments(size_term, poisson_diag)`.
# The control list uses `PLNnetwork_param()`; the older `control_init` /
# `control_main` arguments were removed in PLNmodels 1.0.0. Outputs
# written before these changes used the 500-gene legacy panel and are
# stale.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))

time_levels <- c("48h", "72h", "96h")
min_cells <- 20L
n_penalties <- 40L
min_ratio <- 0.1
penalize_diagonal <- FALSE

slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
out_dir <- file.path(root, "output", "grn", "plnnetwork")
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
# SECTION 2 · PLNnetwork per type and time ----
# ==========================================================================

results <- lapply(time_levels, function(time_level) {
  panel <- read_panel_slot(slots, selected_dir, time_level)
  labels <- panel$cell_type
  lib_size <- panel$library_size
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
    y <- round(panel$y[cells, , drop = FALSE])
    storage.mode(y) <- "integer"
    lib <- as.numeric(lib_size[cells])
    dat <- PLNmodels::prepare_data(
      counts = y,
      covariates = data.frame(log_library = log(lib), row.names = rownames(y)),
      offset = "none"
    )
    message("Fitting PLNnetwork: ", time_level, " / ", tp)
    family <- PLNmodels::PLNnetwork(
      Abundance ~ 1 + offset(log_library),
      data = dat,
      control = PLNmodels::PLNnetwork_param(
        n_penalties = n_penalties,
        min_ratio = min_ratio,
        penalize_diagonal = penalize_diagonal
      )
    )
    # The family is an R6 object (possibly locked), so the metadata go
    # next to it in a plain list rather than into it.
    fits[[tp]] <- list(
      family = family,
      genes = colnames(y),
      n = nrow(y),
      lib_size = lib
    )
  }
  fits
})
names(results) <- time_levels

out_rds <- file.path(out_dir, "plnnetwork_by_type_time.rds")
saveRDS(results, out_rds)
message("Wrote ", out_rds)
