# mkdir -p logs output/grn/plnnetwork
# nohup Rscript --vanilla \
#   scripts/03_03_grn_plnnetwork.R \
#   > "logs/03_03_plnnetwork_$(date +%F)_offset.log" 2>&1 &
#
# Independent PLNnetwork fit per labelled type and time, with a
# library-size offset. Run from the repository root.
#
# The offset is log of the full RNA library L_c (metadata `log_umi`,
# log10 of the full column sum, as in
# scripts/legacy/02_90_stats_prepare_panel.R), not the row sum of
# the 500 panel genes: a panel sum moves with the panel's own expression
# and absorbs part of the cell-type signal. `prepare_data(offset =
# "none")` stops PLNmodels from adding its TSS (panel-sum) offset. The
# control list uses `PLNnetwork_param()`; the older `control_init` /
# `control_main` arguments were removed in PLNmodels 1.0.0. Outputs
# written before this change used the panel-sum offset and are stale.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))

time_levels <- c("48h", "72h", "96h")
cell_type_col <- "celltypeannotation"
min_cells <- 20L
n_penalties <- 40L
min_ratio <- 0.1
penalize_diagonal <- FALSE

input_rds <- file.path(
  root,
  "data",
  "intermediate",
  "equal_type_top500_by_time.rds"
)
out_dir <- file.path(root, "output", "grn", "plnnetwork")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# SECTION 1 · Inputs ----
# ==========================================================================

if (!file.exists(input_rds)) {
  stop(
    "Missing ",
    input_rds,
    ". Run scripts/01_03_export_equal_type_seurat.R first."
  )
}
objects <- readRDS(input_rds)

# ==========================================================================
# SECTION 2 · PLNnetwork per type and time ----
# ==========================================================================

results <- lapply(time_levels, function(time_level) {
  seu <- objects[[time_level]]
  counts <- assay_matrix(seu, "RNA", "counts")
  labels <- seurat_metadata_chr(seu, cell_type_col)
  lib_size <- 10^as.numeric(seurat_metadata_chr(seu, "log_umi"))
  names(lib_size) <- colnames(seu)
  if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    stop("Metadata log_umi is missing or invalid at ", time_level, ".")
  }
  types <- sort(unique(labels[!is.na(labels) & nzchar(labels)]))
  fits <- list()
  for (tp in types) {
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
    storage.mode(y) <- "integer"
    lib <- as.numeric(lib_size[cells])
    dat <- PLNmodels::prepare_data(
      counts = y,
      covariates = data.frame(log_library = log(lib), row.names = cells),
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
