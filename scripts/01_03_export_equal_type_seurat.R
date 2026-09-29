# mkdir -p logs
# nohup Rscript --vanilla \
#   scripts/01_03_export_equal_type_seurat.R \
#   > "logs/01_03_export_top500_$(date +%F).log" 2>&1 &
#
# Subsets the SCTransform allow-list Seurat object to the 500 genes
# selected by the equal-type rank at each bulk-matched time.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "utils_general.R"))

time_levels <- c("48h", "72h", "96h")
n_genes_global <- 500L
cell_type_col <- "celltypeannotation"

slim_path <- file.path(
  root,
  "data",
  "intermediate",
  "suppinger_sct_allowlist_seurat.rds"
)
rank_path <- file.path(
  root,
  "output",
  "naive_marker_selection",
  "tables",
  "ranks_global_equal_type.csv"
)
crosswalk_path <- file.path(
  root,
  "data",
  "dictionaries",
  "suppinger_celltype_crosswalk.csv"
)
out_rds <- file.path(
  root,
  "data",
  "intermediate",
  "equal_type_top500_by_time.rds"
)

# ==========================================================================
# SECTION 1 · Inputs ----
# ==========================================================================

if (!file.exists(slim_path)) {
  stop("Missing ", slim_path, ". Run scripts/01_01 first.")
}
if (!file.exists(rank_path)) {
  stop("Missing ", rank_path, ". Run scripts/01_01 first.")
}

eq <- utils::read.csv(rank_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
crosswalk <- utils::read.csv(
  crosswalk_path,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
exclude_types <- unique(
  crosswalk$celltypeannotation[!crosswalk$include_in_ranking]
)

message("Reading ", slim_path)
seu_sct <- readRDS(slim_path)
rna_counts <- assay_matrix(seu_sct, "RNA", "counts")
keep_export <- !as.character(seu_sct$celltypeannotation) %in% exclude_types

# ==========================================================================
# SECTION 2 · One Seurat object per time ----
# ==========================================================================

objects <- vector("list", length(time_levels))
names(objects) <- time_levels

for (time_level in time_levels) {
  genes <- unique(eq$gene[eq$time_point == time_level & eq$selected %in% TRUE])
  if (length(genes) != n_genes_global) {
    stop(
      "Expected ",
      n_genes_global,
      " selected genes at ",
      time_level,
      ", found ",
      length(genes),
      "."
    )
  }
  missing_genes <- setdiff(genes, rownames(rna_counts))
  if (length(missing_genes) > 0L) {
    stop(
      "Selected genes missing from the allow-list object at ",
      time_level,
      ": ",
      paste(utils::head(missing_genes, 10L), collapse = ", "),
      "."
    )
  }
  cells <- colnames(seu_sct)[
    seu_sct$time_point == time_level & keep_export
  ]
  objects[[time_level]] <- Seurat::CreateSeuratObject(
    counts = rna_counts[genes, cells, drop = FALSE],
    meta.data = seu_sct@meta.data[cells, , drop = FALSE],
    min.cells = 0L,
    min.features = 0L,
    project = paste0("equal_type_top500_", time_level)
  )
  message(
    time_level,
    ": ",
    nrow(objects[[time_level]]),
    " genes x ",
    ncol(objects[[time_level]]),
    " cells"
  )
}

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(objects, out_rds)
message("Wrote ", out_rds)
