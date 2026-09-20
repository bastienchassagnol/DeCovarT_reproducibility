# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/04_01_feature_selection_pseudobulk_deseq2.R \
#   > "logs/04_01_deseq2_$(date +%F)_markers.log" 2>&1 &

# Template only — pseudo-bulk DESeq2 markers per time point / cell type.
# Pattern: GastroDeconv2FateMap
# https://github.com/bastienchassagnol/GastroDeconv2FateMap/blob/main/scripts/06_02_deseq2_biological_exploration_transcriptomic_differential_analyses.R
# Use unnormalised counts; aggregate by sample × cell type; filterByExpr.

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

# Run from the repository root. Paths are relative to that root.

stopifnot(requireNamespace("sctransform", quietly = TRUE))

study <- "suppinger"
technique <- "pseudobulk_deseq2"
today <- format(Sys.Date(), "%Y-%m-%d")

# ==========================================================================
# SECTION 1 · Subset and aggregate ----
# ==========================================================================
# TODO: load processed Seurat/SCE; choose time point(s)
# TODO: aggregate counts by barcode/sample × cell type
# TODO: min cells per pseudo-bulk; design ~ celltype (or one-vs-rest)

# ==========================================================================
# SECTION 2 · DESeq2 / edgeR filter ----
# ==========================================================================
# TODO: DESeqDataSetFromMatrix; LRT or Wald; BH within type and globally
# TODO: write long and wide CSV under output/feature-selection/

# ==========================================================================
# SECTION 3 · Diagnostics (optional) ----
# ==========================================================================
# TODO: volcano / p-value histograms / upset of marker overlap

stop("Template only: implement pseudo-bulk DE; do not invent DEGs.")
