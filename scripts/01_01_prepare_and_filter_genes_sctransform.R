# mkdir -p logs
# nohup Rscript --vanilla \
#   scripts/01_01_prepare_and_filter_genes_sctransform.R \
#   > "logs/01_01_sctransform_$(date +%F)_naive.log" 2>&1 &
#
# Place GSE229513_gastruloidsobject.rds under data/raw/
# (dvc pull data/raw/GSE229513_gastruloidsobject.rds.dvc).
# GSE229386 is the bulk HTSeq table and is not read here.

# ==========================================================================
# SECTION 0 · Dependencies, hyperparameters, paths ----
# ==========================================================================

root <- getwd()
r_dir <- file.path(root, "R")
source(file.path(r_dir, "utils_general.R"))
source(file.path(r_dir, "rank_sctransform_naive.R"))
source(file.path(r_dir, "map_marker_symbols.R"))

n_genes_global <- 500L
n_genes_celltype <- 50L
min_cells <- 20L
sct_seed <- 1L
time_levels <- c("48h", "72h", "96h")
cell_type_col <- "celltypeannotation"

seu_path <- file.path(root, "data", "raw", "GSE229513_gastruloidsobject.rds")
table_dir <- file.path(root, "output", "naive_marker_selection", "tables")
slim_path <- file.path(
  root,
  "data",
  "intermediate",
  "suppinger_sct_allowlist_seurat.rds"
)
result_dir <- file.path(root, "data", "intermediate")
dict_dir <- file.path(root, "data", "dictionaries")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(slim_path), recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# SECTION 1 · Load counts, annotate, keep bulk-matched times ----
# ==========================================================================

if (!file.exists(seu_path)) {
  stop(
    "Missing ",
    seu_path,
    ". From the repository root run: dvc pull ",
    "data/raw/GSE229513_gastruloidsobject.rds.dvc."
  )
}
message("Reading ", seu_path)
seu <- readRDS(seu_path)

counts <- assay_matrix(seu, "RNA", "counts")
sct_genes <- rownames(assay_matrix(seu, "SCT", "counts"))
sct_genes <- intersect(sct_genes, rownames(counts))
meta <- seu@meta.data
rm(seu)
gc()

hours <- as.integer(gsub("[^0-9]", "", as.character(meta$timepoints)))
meta$time_point <- ifelse(
  hours %in% c(48L, 72L, 96L),
  paste0(hours, "h"),
  NA_character_
)
keep_cell <- !is.na(meta$time_point) & !is.na(meta$celltypeannotation)
lib_size <- Matrix::colSums(counts)
keep_cell <- keep_cell & lib_size > 0
counts <- counts[, keep_cell, drop = FALSE]
meta <- meta[keep_cell, , drop = FALSE]
lib_size <- lib_size[keep_cell]
message(
  "Bulk-matched cells: ",
  ncol(counts),
  " ; RNA genes: ",
  nrow(counts),
  " ; SCT allow-list: ",
  length(sct_genes)
)

batch_dict <- utils::read.csv(
  file.path(dict_dir, "suppinger_batch_shapes.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
colour_dict <- utils::read.csv(
  file.path(dict_dir, "suppinger_celltype_colours.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
crosswalk <- utils::read.csv(
  file.path(dict_dir, "suppinger_celltype_crosswalk.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
markers <- utils::read.csv(
  file.path(dict_dir, "mouse_gastruloid_signaling_markers.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

meta$batch_raw <- as.character(meta$batch)
recoded <- batch_dict$batch[match(meta$batch_raw, batch_dict$batch_raw)]
if (anyNA(recoded)) {
  unknown <- unique(meta$batch_raw[is.na(recoded)])
  stop(
    "Batch labels missing from suppinger_batch_shapes.csv: ",
    paste(unknown, collapse = ", "),
    "."
  )
}
meta$batch <- recoded

observed_types <- unique(as.character(meta$celltypeannotation))
missing_types <- setdiff(observed_types, colour_dict$celltypeannotation)
missing_cross <- setdiff(observed_types, crosswalk$celltypeannotation)
if (length(missing_types) > 0L || length(missing_cross) > 0L) {
  stop(
    "Cell types missing from the colour dictionary or the crosswalk: ",
    paste(union(missing_types, missing_cross), collapse = " | "),
    "."
  )
}

exclude_types <- unique(
  crosswalk$celltypeannotation[!crosswalk$include_in_ranking]
)

# ==========================================================================
# SECTION 2 · Gene-wise mean log1p CPM, marker symbols ----
# ==========================================================================

message("Gene-wise mean log1p CPM for the ridge plots ...")
lib_named <- lib_size
names(lib_named) <- colnames(counts)
mean_parts <- list()
for (time_level in time_levels) {
  cells_tp <- colnames(counts)[meta$time_point == time_level]
  types_tp <- unique(as.character(meta[cells_tp, "celltypeannotation"]))
  for (type_name in types_tp) {
    cells_ct <- cells_tp[
      meta[cells_tp, "celltypeannotation"] == type_name
    ]
    message(
      "  ",
      time_level,
      " ",
      type_name,
      " (",
      length(cells_ct),
      " cells)"
    )
    scaled <- counts[, cells_ct, drop = FALSE]
    scaled <- scaled %*%
      Matrix::Diagonal(
        x = 10000 / lib_named[cells_ct]
      )
    scaled@x <- log1p(scaled@x)
    mean_parts[[paste(time_level, type_name)]] <- data.frame(
      time_point = time_level,
      cell_type = type_name,
      gene = rownames(counts),
      mean_log_cpm = as.numeric(Matrix::rowSums(scaled) / ncol(scaled)),
      in_sct_filter = rownames(counts) %in% sct_genes,
      stringsAsFactors = FALSE
    )
    rm(scaled)
  }
}
gene_means <- do.call(rbind, mean_parts)
utils::write.csv(
  gene_means,
  file = file.path(table_dir, "gene_mean_log_cpm.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
rm(gene_means, mean_parts)
gc()

symbol_map <- map_marker_symbols(
  markers = markers$gene_marker,
  universe = rownames(counts)
)
symbol_map$in_rna <- symbol_map$mouse_symbol %in% rownames(counts)
symbol_map$in_sct <- symbol_map$mouse_symbol %in% sct_genes
utils::write.csv(
  symbol_map,
  file = file.path(table_dir, "marker_symbol_map.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
message("Marker symbol map:")
print(symbol_map[, c("gene_marker", "mouse_symbol", "map_method")])

marker_symbol <- ifelse(
  is.na(symbol_map$mouse_symbol),
  symbol_map$gene_marker,
  symbol_map$mouse_symbol
)
venn_sets <- rbind(
  data.frame(set = "rna_counts", gene = rownames(counts)),
  data.frame(set = "sct_prefilter", gene = sct_genes),
  data.frame(set = "markers", gene = unique(marker_symbol))
)
utils::write.csv(
  venn_sets,
  file = file.path(table_dir, "venn_gene_sets.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ==========================================================================
# SECTION 3 · SCTransform v2 ranks and marker union ----
# ==========================================================================

# log_umi is log10 of the full RNA library, not the 2,944-gene subset.
# sctransform v2 uses it as an offset with slope fixed at ln(10).
meta$log_umi <- log10(lib_named[rownames(meta)])
counts_sct <- counts[sct_genes, , drop = FALSE]
rm(counts)
gc()

seu_sct <- Seurat::CreateSeuratObject(
  counts = counts_sct,
  meta.data = meta,
  min.cells = 0L,
  min.features = 0L
)
seu_sct$log_umi <- meta[colnames(seu_sct), "log_umi"]
rm(counts_sct)
gc()

saveRDS(seu_sct, file = slim_path)

ranked <- rank_genes_sctransform(
  seu = seu_sct,
  strategy = "both",
  cell_type_col = cell_type_col,
  split_by = "time_point",
  n_genes_global = n_genes_global,
  n_genes_celltype = n_genes_celltype,
  exclude_cell_types = exclude_types,
  min_cells = min_cells,
  seed = sct_seed
)

marker_links <- merge(
  crosswalk[crosswalk$use_markers, c("celltypeannotation", "cell_type_excel")],
  markers[, c("cell_type_excel", "gene_marker")],
  by = "cell_type_excel",
  all = FALSE
)
marker_links <- merge(
  marker_links,
  symbol_map[, c("gene_marker", "mouse_symbol")],
  by = "gene_marker",
  all.x = TRUE
)
marker_links <- marker_links[!is.na(marker_links$mouse_symbol), , drop = FALSE]

scored_types <- unique(ranked$per_cell_type[, c("time_point", "cell_type")])
names(scored_types)[names(scored_types) == "cell_type"] <- "celltypeannotation"
marker_at_time <- merge(
  scored_types,
  marker_links,
  by = "celltypeannotation",
  all = FALSE
)
marker_at_time$gene <- marker_at_time$mouse_symbol

ranked$global <- attach_union(
  ranks = ranked$global,
  marker_df = marker_at_time,
  by_cell_type = FALSE
)
ranked$global_equal_type <- attach_union(
  ranks = ranked$global_equal_type,
  marker_df = marker_at_time,
  by_cell_type = FALSE
)
ranked$per_cell_type <- attach_union(
  ranks = ranked$per_cell_type,
  marker_df = marker_at_time,
  by_cell_type = TRUE
)

utils::write.csv(
  ranked$global,
  file = file.path(table_dir, "ranks_global.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
utils::write.csv(
  ranked$global_equal_type,
  file = file.path(table_dir, "ranks_global_equal_type.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
utils::write.csv(
  ranked$per_cell_type,
  file = file.path(table_dir, "ranks_per_cell_type.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
utils::write.csv(
  ranked$model,
  file = file.path(table_dir, "model_card.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cell_counts <- as.data.frame(
  table(
    time_point = meta$time_point,
    cell_type = meta$celltypeannotation
  ),
  stringsAsFactors = FALSE,
  responseName = "n_cells"
)
cell_counts <- cell_counts[cell_counts$n_cells > 0L, , drop = FALSE]
union_n <- aggregate(
  in_union ~ time_point + cell_type,
  data = ranked$per_cell_type,
  FUN = sum
)
names(union_n)[names(union_n) == "in_union"] <- "n_genes_union"
cell_counts <- merge(
  cell_counts,
  union_n,
  by = c("time_point", "cell_type"),
  all.x = TRUE
)
cell_counts$n_genes_union[is.na(cell_counts$n_genes_union)] <- 0L
cell_counts <- cell_counts[
  order(cell_counts$time_point, cell_counts$cell_type),
]
utils::write.csv(
  cell_counts,
  file = file.path(table_dir, "cell_counts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

writeLines(
  c(
    "SCTransform v2 design",
    "formula: log(mu_cg) = beta_0g + offset(log(library size_c))",
    "slope: fixed (ln 10 on the log10-depth scale; glmGamPoi offset)",
    "vars.to.regress: none",
    "cell type in the fit: no",
    "cell line / batch in the fit: no",
    "split: one fit per bulk-matched time (48h, 72h, 96h)",
    "gene allow-list: rownames of the author SCT counts layer",
    "library size: column sums of the full RNA count matrix",
    paste("n_genes_global:", n_genes_global),
    paste("n_genes_celltype:", n_genes_celltype),
    paste("min_cells:", min_cells),
    paste("excluded types:", paste(exclude_types, collapse = " | "))
  ),
  con = file.path(table_dir, "model_card.txt")
)
writeLines(
  capture.output(sessionInfo()),
  con = file.path(table_dir, "session_info.txt")
)

# ==========================================================================
# SECTION 4 · Equal-type top 500 raw-count Seurat objects ----
# ==========================================================================

rna_counts <- assay_matrix(seu_sct, "RNA", "counts")
objects <- vector("list", length(time_levels))
names(objects) <- time_levels
eq <- ranked$global_equal_type
keep_export <- !as.character(seu_sct$celltypeannotation) %in% exclude_types

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
out_rds <- file.path(result_dir, "equal_type_top500_by_time.rds")
saveRDS(objects, out_rds)
rm(seu_sct, rna_counts, objects)
gc()

message("Wrote tables under ", table_dir)
message("Wrote ", out_rds)
