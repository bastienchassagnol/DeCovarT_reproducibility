# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_01_stats_quality_control.R \
#   > "logs/02_01_stats_quality_control_$(date +%F).log" 2>&1 &
#
# Written, not run. Tier (i) cell and gene filters. Does not call
# Seurat::SCTransform(). The mean model of later tiers is
# glmGamPoi::glm_gp().

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))

min_cells <- 20L
min_frac <- 0.05
time_levels <- c("48h", "72h", "96h")
cell_type_col <- "celltypeannotation"

seu_path <- file.path(root, "data", "raw", "GSE229513_gastruloidsobject.rds")
bulk_path <- file.path(
  root,
  "data",
  "raw",
  "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz"
)
dict_dir <- file.path(root, "data", "dictionaries")
out_dir <- file.path(root, "output", "feature_selection", "tier_01")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(seu_path)) {
  stop(
    "Missing ", seu_path,
    ". From the repository root run: dvc pull ",
    "data/raw/GSE229513_gastruloidsobject.rds.dvc."
  )
}
if (!file.exists(bulk_path)) {
  stop(
    "Missing ", bulk_path,
    ". Download GSE229386_AllHTSeqCountsWithGeneNames.txt.gz ",
    "into data/raw/ (see data/raw/README.md)."
  )
}

# ==========================================================================
# SECTION 1 · Cells, shared genes, mitochondrial symbols ----
# ==========================================================================

message("Reading ", seu_path)
seu <- readRDS(seu_path)
counts <- assay_matrix(seu, "RNA", "counts")
meta <- seu@meta.data
rm(seu)
gc()

hours <- as.integer(gsub("[^0-9]", "", as.character(meta$timepoints)))
time_point <- ifelse(
  hours %in% c(48L, 72L, 96L),
  paste0(hours, "h"),
  NA_character_
)
cell_type <- as.character(meta[[cell_type_col]])
library_size <- Matrix::colSums(counts)

crosswalk <- utils::read.csv(
  file.path(dict_dir, "suppinger_celltype_crosswalk.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
allowed <- unique(
  crosswalk$celltypeannotation[crosswalk$include_in_ranking]
)

htseq <- read_bulk_htseq(bulk_path)
prepared <- prepare_scored_counts(
  counts = counts,
  time_point = time_point,
  cell_type = cell_type,
  library_size = library_size,
  allowed_types = allowed,
  bulk_genes = htseq$genes$gene_name,
  bulk_ensembl = htseq$genes$ensembl_id,
  bulk_row_id = htseq$genes$row_id,
  time_levels = time_levels,
  min_cells = min_cells
)
rm(counts)
gc()

utils::write.csv(
  prepared$shared,
  file.path(out_dir, "shared_genes.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
utils::write.csv(
  prepared$mito,
  file.path(out_dir, "mito_fraction.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
utils::write.csv(
  data.frame(gene = prepared$dropped_mt, stringsAsFactors = FALSE),
  file.path(out_dir, "dropped_mt_genes.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

slot_counts <- do.call(rbind, lapply(time_levels, function(tm) {
  slot <- prepared$slots[[tm]]
  tab <- table(slot$cell_type)
  data.frame(
    time_point = tm,
    cell_type = names(tab),
    cells = as.integer(tab),
    scored = TRUE,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  slot_counts,
  file.path(out_dir, "cell_counts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ==========================================================================
# SECTION 2 · Presence gate and marker map ----
# ==========================================================================

presence <- lapply(time_levels, function(tm) {
  slot <- prepared$slots[[tm]]
  gate <- presence_gate(slot$counts, slot$cell_type, min_frac = min_frac)
  gate$time_point <- tm
  utils::write.csv(
    gate,
    file.path(out_dir, paste0("presence_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  gate
})
names(presence) <- time_levels

markers <- utils::read.csv(
  file.path(dict_dir, "mouse_gastruloid_signaling_markers.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
symbol_map <- map_marker_symbols(
  markers = markers$gene_marker,
  universe = rownames(prepared$slots[["48h"]]$counts)
)
utils::write.csv(
  symbol_map,
  file.path(out_dir, "marker_symbol_map.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

saveRDS(
  list(
    slots = prepared$slots,
    shared = prepared$shared,
    symbol_map = symbol_map,
    min_cells = min_cells,
    min_frac = min_frac,
    time_levels = time_levels
  ),
  file.path(out_dir, "scored_slots.rds")
)
message("Wrote ", out_dir)
