# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/01_01_prepare_and_filter_genes_sctransform.R \
#   > "logs/01_01_sctransform_$(date +%F)_naive.log" 2>&1 &
#
# If this repository's renv library does not contain Seurat, the
# project .Rprofile hides the user library. Run with --vanilla so
# the user library is used:
#   Rscript --vanilla scripts/01_01_prepare_and_filter_genes_sctransform.R
#
# Optional key=value arguments:
#   seu_path=...  n_genes_global=500  n_genes_celltype=50  min_cells=20
#
# The Seurat RDS is the DVC target
# GastroDeconv2FateMap/data/raw/GSE229513_gastruloidsobject.rds
# (about 14 GB). Pull it with `dvc pull` in that repository. Do not
# copy it into a cloud-synced data/raw. GSE229386 is the bulk HTSeq
# table; this script uses it only to choose the shared times 48 h,
# 72 h and 96 h.

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  file.path(getwd(), "scripts")
}
root <- normalizePath(file.path(script_dir, ".."))
source(file.path(script_dir, "utils_rank_genes_sctransform.R"))

Sys.setenv(
  OMP_NUM_THREADS = "8",
  OPENBLAS_NUM_THREADS = "8",
  MKL_NUM_THREADS = "8"
)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

parse_kv_args <- function(defaults) {
  out <- defaults
  for (arg in commandArgs(trailingOnly = TRUE)) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    if (length(parts) != 2L) {
      stop("Expected key=value arguments, received ", arg, ".")
    }
    key <- parts[[1]]
    if (!key %in% names(out)) {
      stop("Unknown argument ", key, ".")
    }
    out[[key]] <- parts[[2]]
  }
  out
}

opts <- parse_kv_args(list(
  seu_path = "",
  n_genes_global = "500",
  n_genes_celltype = "50",
  min_cells = "20"
))
n_genes_global <- as.integer(opts$n_genes_global)
n_genes_celltype <- as.integer(opts$n_genes_celltype)
min_cells <- as.integer(opts$min_cells)

read_dict <- function(name) {
  path <- file.path(root, "data", "dictionaries", name)
  if (!file.exists(path)) {
    stop("Missing dictionary ", path, ".")
  }
  utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}

write_table <- function(x, name) {
  utils::write.csv(
    x,
    file = file.path(table_dir, name),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

table_dir <- file.path(root, "output", "naive_marker_selection", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# SECTION 1 · Load counts, annotate, keep bulk-matched times ----
# ==========================================================================

candidates <- c(
  opts$seu_path,
  file.path(root, "data", "raw", "GSE229513_gastruloidsobject.rds"),
  file.path(
    root,
    "..",
    "GastroDeconv2FateMap",
    "data",
    "raw",
    "GSE229513_gastruloidsobject.rds"
  ),
  paste0(
    "/mnt/DATA_11TB/projects/dtoo_project/",
    "GastroDeconv2FateMap/data/raw/",
    "GSE229513_gastruloidsobject.rds"
  )
)
candidates <- candidates[nzchar(candidates)]
seu_path <- candidates[file.exists(candidates)][1]
if (is.na(seu_path) || !nzchar(seu_path)) {
  stop(
    "GSE229513_gastruloidsobject.rds was not found. ",
    "From GastroDeconv2FateMap run: ",
    "dvc pull data/raw/GSE229513_gastruloidsobject.rds.dvc"
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

if (!all(c("timepoints", "batch", "celltypeannotation") %in% colnames(meta))) {
  stop(
    "Seurat metadata must contain timepoints, batch and ",
    "celltypeannotation."
  )
}

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

batch_dict <- read_dict("suppinger_batch_shapes.csv")
colour_dict <- read_dict("suppinger_celltype_colours.csv")
crosswalk <- read_dict("suppinger_celltype_crosswalk.csv")
markers <- read_dict("mouse_gastruloid_signaling_markers.csv")

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

batch_by_time <- as.data.frame(
  table(time_point = meta$time_point, batch = meta$batch),
  stringsAsFactors = FALSE,
  responseName = "n_cells"
)
write_table(batch_by_time, "batch_by_time.csv")

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
for (time_level in c("48h", "72h", "96h")) {
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
    scaled <- scaled %*% Matrix::Diagonal(
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
write_table(gene_means, "gene_mean_log_cpm.csv")
rm(gene_means, mean_parts)
gc()

symbol_map <- map_marker_symbols(
  markers = markers$gene_marker,
  universe = rownames(counts)
)
symbol_map$in_rna <- symbol_map$mouse_symbol %in% rownames(counts)
symbol_map$in_sct <- symbol_map$mouse_symbol %in% sct_genes
write_table(symbol_map, "marker_symbol_map.csv")
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
write_table(venn_sets, "venn_gene_sets.csv")

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

ranked <- rank_genes_sctransform(
  seu = seu_sct,
  strategy = "both",
  cell_type_col = "celltypeannotation",
  split_by = "time_point",
  n_genes_global = n_genes_global,
  n_genes_celltype = n_genes_celltype,
  exclude_cell_types = exclude_types,
  min_cells = min_cells,
  seed = 1L
)
rm(seu_sct)
gc()

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

attach_union <- function(ranks, marker_df, by_cell_type) {
  ranks$in_marker_panel <- FALSE
  if (by_cell_type) {
    id_rank <- paste(ranks$time_point, ranks$cell_type, ranks$gene, sep = "\r")
    id_mark <- paste(
      marker_df$time_point,
      marker_df$celltypeannotation,
      marker_df$gene,
      sep = "\r"
    )
    key_df <- unique(marker_df[, c(
      "time_point",
      "celltypeannotation",
      "gene"
    )])
  } else {
    id_rank <- paste(ranks$time_point, ranks$gene, sep = "\r")
    id_mark <- paste(marker_df$time_point, marker_df$gene, sep = "\r")
    key_df <- unique(marker_df[, c("time_point", "gene")])
  }
  if (nrow(key_df) > 0L) {
    ranks$in_marker_panel <- id_rank %in% unique(id_mark)
    missing_id <- unique(id_mark[!id_mark %in% id_rank])
    if (length(missing_id) > 0L) {
      if (by_cell_type) {
        bits <- strsplit(missing_id, "\r", fixed = TRUE)
        extra <- data.frame(
          time_point = vapply(bits, `[`, character(1), 1L),
          cell_type = vapply(bits, `[`, character(1), 2L),
          gene = vapply(bits, `[`, character(1), 3L),
          stringsAsFactors = FALSE
        )
      } else {
        bits <- strsplit(missing_id, "\r", fixed = TRUE)
        extra <- data.frame(
          time_point = vapply(bits, `[`, character(1), 1L),
          cell_type = NA_character_,
          gene = vapply(bits, `[`, character(1), 2L),
          stringsAsFactors = FALSE
        )
      }
      extra$ss_between <- NA_real_
      extra$ss_within <- NA_real_
      extra$between_over_within <- NA_real_
      extra$between_over_total <- NA_real_
      extra$top_cell_type <- extra$cell_type
      extra$passes_between_gt_within <- NA
      extra$rank <- NA_integer_
      extra$selected <- FALSE
      extra$in_marker_panel <- TRUE
      ranks <- rbind(ranks, extra[colnames(ranks)])
    }
  }
  ranks$in_union <- ranks$selected %in% TRUE | ranks$in_marker_panel %in% TRUE
  ranks$source <- ifelse(
    ranks$selected %in% TRUE & ranks$in_marker_panel %in% TRUE,
    "both",
    ifelse(
      ranks$selected %in% TRUE,
      "sctransform",
      ifelse(ranks$in_marker_panel %in% TRUE, "marker", "screened")
    )
  )
  ranks
}

ranked$global <- attach_union(
  ranks = ranked$global,
  marker_df = marker_at_time,
  by_cell_type = FALSE
)
ranked$per_cell_type <- attach_union(
  ranks = ranked$per_cell_type,
  marker_df = marker_at_time,
  by_cell_type = TRUE
)

write_table(ranked$global, "ranks_global.csv")
write_table(ranked$per_cell_type, "ranks_per_cell_type.csv")
write_table(ranked$model, "model_card.csv")

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
cell_counts <- cell_counts[order(cell_counts$time_point, cell_counts$cell_type), ]
write_table(cell_counts, "cell_counts.csv")

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

message("Wrote tables under ", table_dir)
