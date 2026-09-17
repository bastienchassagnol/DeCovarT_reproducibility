# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/01_01_single_cell_exploration_suppinger.R \
#   > "logs/01_01_single_cell_$(date +%F)_explore.log" 2>&1 &

# Template only — port from GastroDeconv2FateMap
# https://github.com/bastienchassagnol/GastroDeconv2FateMap/blob/main/scripts/01_01_single_cell_exploration_suppinger.R
# Paper: Suppinger et al., Cell Stem Cell 2023
# https://www.cell.com/cell-stem-cell/pdfExtended/S1934-5909(23)00170-4

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

source(here::here("R", "paths.R"))

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- path_output("single-cell")
# dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# SECTION 1 · Load Seurat object (GSE229513) ----
# ==========================================================================
# TODO: readRDS(path_raw("GSE229513_gastruloidsobject.rds"))
# TODO: inspect dim(), Idents(), meta.data (batch, timepoints,
#       celltypeannotation)
# TODO: recode batch labels (B-S / SBR); attach cell-type colour map

# ==========================================================================
# SECTION 2 · PCA + UMAP per time point ----
# ==========================================================================
# TODO: subset by ident; Seurat::DimPlot PCA and UMAP; patchwork; ggsave PDF

# ==========================================================================
# SECTION 3 · Marker genes (optional diagnostic) ----
# ==========================================================================
# TODO: FindAllMarkers on a chosen time point; tinytable of top genes
# Full DE for signatures belongs in scripts/04_01_*.R (pseudo-bulk DESeq2)

# ==========================================================================
# SECTION 4 · Barcode table (optional) ----
# ==========================================================================
# TODO: parse GSE229513_barcodes.tsv.gz if needed for FateMap-style IDs

stop("Template only: implement exploration; do not invent plots.")
