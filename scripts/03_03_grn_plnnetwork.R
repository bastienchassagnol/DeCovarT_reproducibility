# mkdir -p logs output/grn/plnnetwork
# nohup Rscript --vanilla \
#   scripts/03_03_grn_plnnetwork.R \
#   > "logs/03_03_plnnetwork_$(date +%F)_offset.log" 2>&1 &
#
# Independent PLNnetwork fit per labelled type and time, with a
# library-size offset. Run from the repository root.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "utils_general.R"))

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
    lib <- pmax(as.numeric(Matrix::rowSums(y)), 1)
    dat <- PLNmodels::prepare_data(
      counts = y,
      covariates = data.frame(log_library = log(lib))
    )
    message("Fitting PLNnetwork: ", time_level, " / ", tp)
    fits[[tp]] <- PLNmodels::PLNnetwork(
      Abundance ~ 1 + offset(log_library),
      data = dat,
      control_init = list(nPenalties = n_penalties, min.ratio = min_ratio),
      control_main = list(penalize_diagonal = penalize_diagonal)
    )
    fits[[tp]]$genes <- colnames(y)
    fits[[tp]]$n <- nrow(y)
  }
  fits
})
names(results) <- time_levels

out_rds <- file.path(out_dir, "plnnetwork_by_type_time.rds")
saveRDS(results, out_rds)
message("Wrote ", out_rds)
