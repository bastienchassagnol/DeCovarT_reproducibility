# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_01_stats_quality_control.R \
#   > "logs/02_01_stats_quality_control_$(date +%F).log" 2>&1 &
#
# Tier (i) cell and gene filters. Writes three summary CSVs and one RDS
# (tracked with DVC) holding every table. Does not call
# Seurat::SCTransform(): the authors' SCT residuals are read only for the
# appendix projections, and their gene names for comparison. The mean
# model of later tiers is glmGamPoi::glm_gp(). Fractional (multi-mapping)
# counts are kept as they are.

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))

min_cells <- 20L # smallest scored type, and Luecken & Theis' gene rule
min_cells_type <- 2L # non-zero cells per type: non-singular covariance
max_fractional_share <- 1 # keep every gene, fractional counts included
time_levels <- c("48h", "72h", "96h")
cell_type_col <- "celltypeannotation"
batch_col <- "batch"

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

for (path in c(seu_path, bulk_path)) {
  if (!file.exists(path)) {
    stop(
      "Missing ",
      path,
      ". From the repository root run: dvc pull ",
      "data/raw/",
      basename(path),
      ".dvc"
    )
  }
}

write_out <- function(x, name) {
  utils::write.csv(
    x,
    file.path(out_dir, name),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

# ==========================================================================
# SECTION 1 · Read the deposited object ----
# ==========================================================================

message("Reading ", seu_path)
seu <- readRDS(seu_path)
counts <- assay_matrix(seu, "RNA", "counts")
meta <- seu@meta.data
# The authors' SCT layer holds their CCA integration features. Only the
# gene names are used, to compare with the filters below.
sct_genes <- rownames(assay_matrix(seu, "SCT", "counts"))
anchor_features <- seu@commands$FindIntegrationAnchors@params$anchor.features
# Authors' SCT Pearson residuals, for the appendix projections only.
sct_layer <- "scale.data"
sct_resid <- assay_matrix(seu, "SCT", sct_layer)
if (is.null(sct_resid) || nrow(sct_resid) == 0L) {
  sct_layer <- "data"
  sct_resid <- assay_matrix(seu, "SCT", sct_layer)
}
message("SCT layer for the projections: ", sct_layer, " (", nrow(sct_resid), " genes)")
rm(seu)
gc()

hours <- as.integer(gsub("[^0-9]", "", as.character(meta$timepoints)))
time_point <- ifelse(
  hours %in% c(48L, 72L, 96L),
  paste0(hours, "h"),
  NA_character_
)
cell_type <- as.character(meta[[cell_type_col]])
batch <- as.character(meta[[batch_col]])
# Full library size over the 27,027 deposited genes. It equals 10^log_umi
# of the authors' SCT models; meta$nCount_RNA describes another matrix.
library_size <- Matrix::colSums(counts)

# Cell types and curated markers (one row per type x marker group x gene).
type_markers <- utils::read.csv(
  file.path(dict_dir, "suppinger_celltype_markers.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA"),
  fileEncoding = "UTF-8"
)
allowed <- unique(
  type_markers$celltypeannotation[type_markers$include_in_ranking %in% TRUE]
)

# ==========================================================================
# SECTION 2 · Literature cell-QC rules on the authors' cells ----
# ==========================================================================

matched <- keep_bulk_matched_cells(time_point, cell_type, library_size)
cell_qc <- cell_qc_rules(
  counts[, matched, drop = FALSE],
  time_point[matched],
  cell_type[matched]
)
write_out(cell_qc, "author_qc_audit.csv")

# Per-cell QC values, kept in the RDS for the tier (i) figures.
cell_values <- data.frame(
  barcode = colnames(counts)[matched],
  time_point = time_point[matched],
  cell_type = cell_type[matched],
  batch = batch[matched],
  library_size = library_size[matched],
  detected_genes = Matrix::colSums(counts[, matched, drop = FALSE] > 0),
  stringsAsFactors = FALSE
)
cell_values$mad_flag <- stats::ave(
  log(cell_values$library_size),
  cell_values$time_point,
  cell_values$cell_type,
  FUN = function(v) {
    s <- stats::mad(v)
    if (is.finite(s) && s > 0) abs(v - stats::median(v)) > 3 * s else 0
  }
) ==
  1

# ==========================================================================
# SECTION 3 · Cells, shared genes, fractional counts ----
# ==========================================================================

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
  batch = batch,
  time_levels = time_levels,
  min_cells = min_cells,
  max_fractional_share = max_fractional_share
)

fractional <- prepared$fractional
fractional$dropped <- fractional$gene %in% prepared$dropped_fractional
fractional$shared_with_bulk <- fractional$gene %in%
  prepared$shared$gene[prepared$shared$keep]

slot_counts <- do.call(
  rbind,
  lapply(time_levels, function(tm) {
    slot <- prepared$slots[[tm]]
    tab <- table(slot$cell_type)
    data.frame(
      time_point = tm,
      cell_type = names(tab),
      cells = as.integer(tab),
      scored = TRUE,
      stringsAsFactors = FALSE
    )
  })
)
write_out(slot_counts, "cell_counts.csv")

# ==========================================================================
# SECTION 4 · Gene rules per slot, presence gate, marker map ----
# ==========================================================================

# Gene-level facts that do not depend on the slot.
gene_facts <- data.frame(
  gene = rownames(counts),
  n_nonzero = fractional$n_nonzero,
  fractional_share = fractional$fractional_share,
  shared_with_bulk = fractional$shared_with_bulk,
  in_sct_layer = rownames(counts) %in% sct_genes,
  in_anchor_features = rownames(counts) %in% anchor_features,
  stringsAsFactors = FALSE
)

gene_rules <- do.call(
  rbind,
  lapply(time_levels, function(tm) {
    slot <- prepared$slots[[tm]]
    rules <- gene_filter_rules(
      counts[, slot$cells, drop = FALSE],
      slot$cell_type,
      min_cells_total = min_cells,
      min_cells_type = min_cells_type
    )
    rules$time_point <- tm
    rules$enters_tier_ii <- rules$pass_presence &
      rules$gene %in% rownames(slot$counts)
    merge(gene_facts, rules, by = "gene", sort = FALSE)
  })
)

gene_rule_summary <- do.call(
  rbind,
  lapply(time_levels, function(tm) {
    r <- gene_rules[gene_rules$time_point == tm, , drop = FALSE]
    data.frame(
      time_point = tm,
      rule = c(
        "deposited genes",
        "zero in every cell of the slot",
        "not shared one-to-one with the bulk",
        "fails the presence gate (>= 20 cells, >= 2 in every type, variance > 0)",
        "in the authors' SCT layer (integration features)",
        "enters tier (ii)",
        "enters tier (ii) and in the SCT layer"
      ),
      genes = c(
        nrow(r),
        sum(r$all_zero),
        sum(!r$shared_with_bulk),
        sum(!r$pass_presence),
        sum(r$in_sct_layer),
        sum(r$enters_tier_ii),
        sum(r$enters_tier_ii & r$in_sct_layer)
      ),
      stringsAsFactors = FALSE
    )
  })
)
write_out(gene_rule_summary, "gene_filter_summary.csv")
rm(counts)
gc()

presence <- lapply(time_levels, function(tm) {
  slot <- prepared$slots[[tm]]
  gate <- presence_gate(
    slot$counts,
    slot$cell_type,
    min_cells_total = min_cells,
    min_cells_type = min_cells_type
  )
  gate$time_point <- tm
  gate
})
names(presence) <- time_levels

symbol_map <- map_marker_symbols(
  markers = unique(stats::na.omit(type_markers$gene_marker)),
  universe = rownames(prepared$slots[["48h"]]$counts)
)

# Fractional entries per gene, type and time (scored cells), and the
# biotype of every gene, for the fractional-count figure.
fractional_group <- do.call(
  rbind,
  lapply(time_levels, function(tm) {
    slot <- prepared$slots[[tm]]
    out <- fractional_by_group(slot$counts, slot$cell_type)
    names(out)[names(out) == "group"] <- "cell_type"
    out$time_point <- tm
    out
  })
)
biotype <- gene_biotypes(rownames(prepared$slots[[1L]]$counts))

# Authors' SCT residuals of the scored cells, one matrix per slot.
sct_slots <- lapply(prepared$slots, function(slot) {
  as.matrix(sct_resid[, slot$cells, drop = FALSE])
})
rm(sct_resid)

saveRDS(
  list(
    slots = prepared$slots,
    shared = prepared$shared,
    symbol_map = symbol_map,
    dropped_fractional = prepared$dropped_fractional,
    dropped_mt = prepared$dropped_mt,
    mito = prepared$mito,
    fractional = fractional,
    fractional_group = fractional_group,
    biotype = biotype,
    sct_slots = sct_slots,
    sct_layer = sct_layer,
    gene_rules = gene_rules,
    presence = presence,
    cell_qc = cell_qc,
    cell_values = cell_values,
    sct_genes = sct_genes,
    anchor_missing = setdiff(anchor_features, sct_genes),
    min_cells = min_cells,
    min_cells_type = min_cells_type,
    max_fractional_share = max_fractional_share,
    time_levels = time_levels
  ),
  file.path(out_dir, "scored_slots.rds")
)
message("Wrote ", out_dir)
