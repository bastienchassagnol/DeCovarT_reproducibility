# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_01_vis_projections.R \
#   > "logs/02_01_vis_projections_$(date +%F).log" 2>&1 &
#
# Written, not run. One PCA scatter and one PLS-DA scatter per time,
# coloured by scored cell type. No barcode distribution.

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_01_vis_projections.R"))

seu_path <- file.path(root, "data", "raw", "GSE229513_gastruloidsobject.rds")
slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
fig_dir <- file.path(root, "output", "feature_selection", "tier_01", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(seu_path)) {
  stop(
    "Missing ", seu_path,
    ". From the repository root run: dvc pull ",
    "data/raw/GSE229513_gastruloidsobject.rds.dvc."
  )
}
if (!file.exists(slots_path)) {
  stop(
    "Missing ", slots_path,
    ". Run scripts/02_01_stats_quality_control.R first."
  )
}

colours <- utils::read.csv(
  file.path(root, "data", "dictionaries", "suppinger_celltype_colours.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
colour_map <- stats::setNames(colours$colour, colours$celltypeannotation)

prepared <- readRDS(slots_path)

# ==========================================================================
# SECTION 1 · PCA and PLS-DA ----
# ==========================================================================

for (tm in names(prepared$slots)) {
  slot <- prepared$slots[[tm]]
  size_factors <- slot$library_size / stats::median(slot$library_size)
  drawn <- plot_type_pca_plsda(
    counts = slot$counts,
    cell_type = slot$cell_type,
    size_factors = size_factors,
    colours = colour_map
  )
  ggplot2::ggsave(
    file.path(fig_dir, paste0("pca_cell_type_", tm, ".pdf")),
    drawn$pca,
    width = 8,
    height = 6,
    dpi = 320,
    limitsize = FALSE
  )
  ggplot2::ggsave(
    file.path(fig_dir, paste0("plsda_cell_type_", tm, ".pdf")),
    drawn$plsda,
    width = 8,
    height = 6,
    dpi = 320,
    limitsize = FALSE
  )
  message(tm, ": ", drawn$n_genes, " varying genes in the projections.")
}
