# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_03_cross_platform_consistency.R \
#   > "logs/02_03_cross_platform_consistency_$(date +%F).log" 2>&1 &
#
# Written, not run. Tier (iii): range check, deMixSC-style plug-in,
# tail statements T1-T5, quantile band, Cook's distance. The mean model
# inside cook_cells() is glmGamPoi::glm_gp().

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "02_03_influence.R"))
source(file.path(root, "R", "02_03_bulk_consistency.R"))
source(file.path(root, "R", "02_03_plugin_composition.R"))
source(file.path(root, "R", "02_03_tail_bounds.R"))

z_range <- 3
mad_k <- 3
level <- 0.01
c_reg <- 2
max_degree <- 10L
mean_degree <- 3
h_mult <- 3
r_max <- 3

seu_path <- file.path(root, "data", "raw", "GSE229513_gastruloidsobject.rds")
bulk_path <- file.path(
  root,
  "data",
  "raw",
  "GSE229386_AllHTSeqCountsWithGeneNames.txt.gz"
)
slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
rank_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_02",
  "cross_cell.rds"
)
out_dir <- file.path(root, "output", "feature_selection", "tier_03")
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
if (!file.exists(slots_path)) {
  stop(
    "Missing ", slots_path,
    ". Run scripts/02_01_stats_quality_control.R first."
  )
}
if (!file.exists(rank_path)) {
  stop(
    "Missing ", rank_path,
    ". Run scripts/02_02_cross_cell_consistency.R first."
  )
}

prepared <- readRDS(slots_path)
ranked <- readRDS(rank_path)
times <- prepared$time_levels
htseq <- read_bulk_htseq(bulk_path)

marker_genes <- prepared$symbol_map$mouse_symbol
marker_genes <- unique(marker_genes[!is.na(marker_genes)])
panel <- unique(c(ranked$panel, marker_genes))

dict <- utils::read.csv(
  file.path(
    root,
    "data",
    "dictionaries",
    "mouse_gastruloid_signaling_markers.csv"
  ),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
dict$gene <- prepared$symbol_map$mouse_symbol[
  match(dict$gene_marker, prepared$symbol_map$gene_marker)
]
modules <- lapply(split(dict$gene, dict$cell_type_excel), function(genes) {
  unique(stats::na.omit(genes))
})

coef_matrix <- function(composition, types) {
  samples <- unique(composition$sample)
  out <- matrix(
    NA_real_,
    nrow = length(types),
    ncol = length(samples),
    dimnames = list(types, samples)
  )
  for (i in seq_len(nrow(composition))) {
    out[composition$cell_type[i], composition$sample[i]] <- composition$coef[i]
  }
  out
}

# ==========================================================================
# SECTION 1 · Range, tails, Cook's distance ----
# ==========================================================================

range_flags <- list()
tail_res <- list()
outlier_cells <- list()
outlier_genes <- list()

for (tm in times) {
  slot <- prepared$slots[[tm]]
  genes <- intersect(panel, rownames(slot$counts))
  x <- slot$counts[genes, , drop = FALSE]
  ct <- droplevels(factor(slot$cell_type))
  sf <- slot$library_size / stats::median(slot$library_size)
  norm <- sweep(as.matrix(x), 2L, sf, "/")
  mu <- sapply(levels(ct), function(k) {
    rowMeans(norm[, ct == k, drop = FALSE])
  })
  sd <- sapply(levels(ct), function(k) {
    apply(norm[, ct == k, drop = FALSE], 1L, stats::sd)
  })
  samples <- standard_wnt_samples(colnames(htseq$counts), tm)
  bulk <- bulk_symbol_counts(htseq, prepared$shared, genes)
  bulk <- bulk[, samples, drop = FALSE]

  range_one <- bulk_consistency_filter(
    bulk = bulk,
    mu = mu,
    sd = sd,
    n_cells = table(ct)[levels(ct)],
    z_range = z_range,
    mad_k = mad_k
  )
  range_one$time_point <- tm
  range_flags[[tm]] <- range_one

  tail_res[[tm]] <- bulk_tail_filter(
    bulk = bulk,
    counts = x,
    cell_type = slot$cell_type,
    modules = modules,
    level = level,
    c_reg = c_reg,
    max_degree = max_degree,
    mean_degree = mean_degree
  )

  outlier_cells[[tm]] <- cook_cells(
    x,
    slot$cell_type,
    size_factors = sf
  )
  plug <- tail_res[[tm]]$plugin
  outlier_genes[[tm]] <- cook_genes(
    bulk,
    tail_res[[tm]]$moments$mu,
    coef = coef_matrix(plug$composition, colnames(tail_res[[tm]]$moments$mu)),
    weights = plug$weights,
    h_mult = h_mult,
    r_max = r_max
  )
}

for (tm in times) {
  utils::write.csv(
    range_flags[[tm]],
    file.path(out_dir, paste0("range_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  utils::write.csv(
    tail_res[[tm]]$gene,
    file.path(out_dir, paste0("tail_gene_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  if (!is.null(tail_res[[tm]]$module) && nrow(tail_res[[tm]]$module) > 0L) {
    utils::write.csv(
      tail_res[[tm]]$module,
      file.path(out_dir, paste0("tail_module_", tm, ".csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }
  overlap <- tail_flag_overlap(
    tail_res[[tm]]$gene,
    band = tail_res[[tm]]$quantile_band
  )
  utils::write.csv(
    as.data.frame(overlap$overlap$jaccard),
    file.path(out_dir, paste0("tail_jaccard_", tm, ".csv")),
    fileEncoding = "UTF-8"
  )
  utils::write.csv(
    overlap$overlap$patterns,
    file.path(out_dir, paste0("tail_patterns_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  utils::write.csv(
    overlap$order,
    file.path(out_dir, paste0("tail_order_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  utils::write.csv(
    outlier_genes[[tm]],
    file.path(out_dir, paste0("cook_genes_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  utils::write.csv(
    outlier_cells[[tm]]$flagged,
    file.path(out_dir, paste0("cook_cells_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  utils::write.csv(
    outlier_cells[[tm]]$summary,
    file.path(out_dir, paste0("cook_cell_summary_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

saveRDS(
  list(
    range = range_flags,
    tails = tail_res,
    outlier_cells = outlier_cells,
    outlier_genes = outlier_genes,
    panel = panel
  ),
  file.path(out_dir, "cross_platform.rds")
)
message("Wrote ", out_dir)
