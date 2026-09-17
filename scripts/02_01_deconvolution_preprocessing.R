# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_01_deconvolution_preprocessing.R \
#   > "logs/02_01_preprocess_$(date +%F)_suppinger.log" 2>&1 &

# Template only — port from GastroDeconv2FateMap
# https://github.com/bastienchassagnol/GastroDeconv2FateMap/blob/main/scripts/02_01_deconvolution_preprocessing.R

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

# Run from the repository root. Paths are relative to that root.

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")

# ==========================================================================
# SECTION 1 · Load single-cell and bulk ----
# ==========================================================================
# TODO: Seurat RDS GSE229513; HTSeq table GSE229386

# ==========================================================================
# SECTION 2 · Harmonise gene names ----
# ==========================================================================
# TODO: overlap gene symbols; drop ambiguous Ensembl↔symbol maps
# TODO: optional Venn (ggVennDiagram); annotate biotype locally
# TODO: write a codebook row per retained gene → data/dictionaries/

# ==========================================================================
# SECTION 3 · SummarizedExperiment (bulk) ----
# ==========================================================================
# TODO: integer count matrix; colData (time, batch, treatment)
# TODO: GEO metadata list (GSE229386); saveRDS under data/intermediate/

# ==========================================================================
# SECTION 4 · SingleCellExperiment (RNA assay, shared genes) ----
# ==========================================================================
# TODO: subset features; factor timepoints; cell-type colours
# TODO: Seurat::as.SingleCellExperiment(); saveRDS

stop("Template only: implement preprocessing; do not invent objects.")
