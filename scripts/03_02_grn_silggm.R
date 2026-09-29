# mkdir -p logs output/grn/silggm
# nohup Rscript --vanilla \
#   scripts/03_02_grn_silggm.R \
#   > "logs/03_02_silggm_$(date +%F)_dsgl.log" 2>&1 &
#
# Independent SILGGM fit per labelled type and time. Nonparanormal
# shrinkage, de-sparsified graphical lasso, global FDR. The R object is
# the only export (no Cytoscape / CSV dump). Run from the repository
# root.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "utils_general.R"))

time_levels <- c("48h", "72h", "96h")
cell_type_col <- "celltypeannotation"
min_cells <- 20L
npn_fun <- "shrinkage"
silggm_method <- "D-S_GL"
global_flag <- TRUE
alpha_level <- 0.05
seed <- 1L

input_rds <- file.path(
  root,
  "data",
  "intermediate",
  "equal_type_top500_by_time.rds"
)
out_dir <- file.path(root, "output", "grn", "silggm")
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
# SECTION 2 · Nonparanormal + SILGGM per type and time ----
# ==========================================================================

withr::with_seed(seed, {
  results <- lapply(time_levels, function(time_level) {
    seu <- objects[[time_level]]
    counts <- assay_matrix(seu, "RNA", "counts")
    types <- sort(unique(as.character(seu[[cell_type_col]])))
    fits <- list()
    for (tp in types) {
      cells <- colnames(seu)[as.character(seu[[cell_type_col]]) == tp]
      if (length(cells) < min_cells) {
        message("Skipping ", time_level, " ", tp, " (n = ", length(cells), ")")
        next
      }
      x <- t(as.matrix(counts[, cells, drop = FALSE]))
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
    }
    fits
  })
  names(results) <- time_levels
})

out_rds <- file.path(out_dir, "silggm_dsgl_by_type_time.rds")
saveRDS(results, out_rds)
message("Wrote ", out_rds)
