# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/03_01_grn_graphical_lasso.R \
#   > "logs/03_01_grn_$(date +%F)_glasso.log" 2>&1 &

# Infer a sparse precision matrix (GGM) per cell type via graphical lasso
# (`huge`), using the processed single-cell reference from 02_01_*.
# Align with DeCovarT covariance vignettes (huge / glasso path).

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

source(here::here("R", "paths.R"))

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
# seed <- 1L  # set explicitly when implementing

# ==========================================================================
# SECTION 1 · Inputs ----
# ==========================================================================
# TODO: load SingleCellExperiment / Seurat from data/intermediate/
# TODO: restrict to the gene panel from 04_* (or a HVG pre-filter)
# TODO: split cells by celltypeannotation (and optionally time point)

# ==========================================================================
# SECTION 2 · Graphical lasso per cell type ----
# ==========================================================================
# TODO: log-normalise / scale genes within type
# TODO: huge::huge(..., method = "glasso") + huge::huge.select()
# TODO: store Omega_j (precision) and Sigma_j = solve(Omega_j)
# TODO: record lambda, sparsity, n_cells, n_genes per type

# ==========================================================================
# SECTION 3 · Outputs ----
# ==========================================================================
# TODO: save RDS/CSV under data/processed/grn/
# TODO: optional igraph/ggplot of the sparsest selected graph

stop("Template only: implement glasso per cell type; do not invent GRNs.")
