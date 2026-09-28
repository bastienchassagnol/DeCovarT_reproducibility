# mkdir -p logs
# nohup Rscript --vanilla \
#   scripts/01_03_sensitivity_48h_cell_line.R \
#   > "logs/01_03_sensitivity_48h_$(date +%F).log" 2>&1 &
#
# Refits SCTransform v2 inside each cell line at 48 h, the only
# bulk-matched slot that contains both B-S and SBR. Reads the slim
# Seurat object written by
# scripts/01_01_prepare_and_filter_genes_sctransform.R.

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

slim_path <- file.path(
  root,
  "data",
  "intermediate",
  "suppinger_sct_allowlist_seurat.rds"
)
if (!file.exists(slim_path)) {
  stop(
    "Missing ",
    slim_path,
    ". Run scripts/01_01_prepare_and_filter_genes_sctransform.R first."
  )
}

out_dir <- file.path(root, "output", "naive_marker_selection", "cell_line_48h")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

crosswalk <- utils::read.csv(
  file.path(root, "data", "dictionaries", "suppinger_celltype_crosswalk.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
exclude_types <- unique(
  crosswalk$celltypeannotation[!crosswalk$include_in_ranking]
)

save_pages <- function(grobs, filename, width, height) {
  pages <- gridExtra::marrangeGrob(
    grobs = grobs,
    nrow = 1L,
    ncol = 1L,
    top = NULL
  )
  ggplot2::ggsave(
    filename = file.path(out_dir, filename),
    plot = pages,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    limitsize = FALSE
  )
}

# ==========================================================================
# SECTION 1 · Inputs ----
# ==========================================================================

message("Reading ", slim_path)
seu <- readRDS(slim_path)
seu <- subset(seu, cells = colnames(seu)[seu$time_point == "48h"])
lines <- intersect(c("B-S", "SBR"), unique(as.character(seu$batch)))
if (length(lines) < 2L) {
  stop("48 h does not contain both B-S and SBR.")
}

# ==========================================================================
# SECTION 2 · Method ----
# ==========================================================================

fit_line <- function(line) {
  message("SCTransform v2 at 48 h, line ", line, " ...")
  cells <- colnames(seu)[as.character(seu$batch) == line]
  sub <- subset(seu, cells = cells)
  ranked <- rank_genes_sctransform(
    seu = sub,
    strategy = "both",
    cell_type_col = "celltypeannotation",
    split_by = NULL,
    n_genes_global = 500L,
    n_genes_celltype = 50L,
    exclude_cell_types = exclude_types,
    min_cells = 20L,
    seed = 1L
  )
  ranked$global$time_point <- line
  ranked$per_cell_type$time_point <- line
  ranked
}

fits <- lapply(lines, fit_line)
names(fits) <- lines
rm(seu)
gc()

global <- do.call(rbind, lapply(fits, `[[`, "global"))
per_type <- do.call(rbind, lapply(fits, `[[`, "per_cell_type"))
utils::write.csv(
  global,
  file = file.path(out_dir, "ranks_global_48h_by_line.csv"),
  row.names = FALSE
)
utils::write.csv(
  per_type,
  file = file.path(out_dir, "ranks_per_cell_type_48h_by_line.csv"),
  row.names = FALSE
)

selected_genes <- function(df, line, cell_type = NULL) {
  keep <- df$time_point == line & df$selected %in% TRUE
  if (!is.null(cell_type)) {
    keep <- keep & df$cell_type == cell_type
  }
  unique(df$gene[keep])
}

# ==========================================================================
# SECTION 3 · Outputs ----
# ==========================================================================

bs_global <- selected_genes(global, "B-S")
sbr_global <- selected_genes(global, "SBR")
shared_global <- intersect(bs_global, sbr_global)
overlap <- data.frame(
  contrast = "global_equal_cell_top500",
  cell_type = NA_character_,
  n_bs = length(bs_global),
  n_sbr = length(sbr_global),
  n_shared = length(shared_global),
  jaccard = length(shared_global) /
    length(union(bs_global, sbr_global)),
  stringsAsFactors = FALSE
)

types_bs <- unique(per_type$cell_type[per_type$time_point == "B-S"])
types_sbr <- unique(per_type$cell_type[per_type$time_point == "SBR"])
shared_types <- intersect(types_bs, types_sbr)
one_sided <- lapply(
  setdiff(union(types_bs, types_sbr), shared_types),
  function(type_name) {
    a <- selected_genes(per_type, "B-S", type_name)
    b <- selected_genes(per_type, "SBR", type_name)
    data.frame(
      contrast = "per_cell_type_top50",
      cell_type = type_name,
      n_bs = length(a),
      n_sbr = length(b),
      n_shared = length(intersect(a, b)),
      jaccard = NA_real_,
      stringsAsFactors = FALSE
    )
  }
)
type_rows <- lapply(shared_types, function(type_name) {
  a <- selected_genes(per_type, "B-S", type_name)
  b <- selected_genes(per_type, "SBR", type_name)
  data.frame(
    contrast = "per_cell_type_top50",
    cell_type = type_name,
    n_bs = length(a),
    n_sbr = length(b),
    n_shared = length(intersect(a, b)),
    jaccard = length(intersect(a, b)) / length(union(a, b)),
    stringsAsFactors = FALSE
  )
})
overlap <- rbind(overlap, do.call(rbind, c(type_rows, one_sided)))
utils::write.csv(
  overlap,
  file = file.path(out_dir, "overlap_48h_cell_line.csv"),
  row.names = FALSE
)

venn_page <- function(set_a, set_b, title, subtitle) {
  plot <- ggVennDiagram::ggVennDiagram(
    list(`B-S` = set_a, SBR = set_b),
    label_alpha = 0
  ) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::scale_fill_gradient(low = "#f7f7f7", high = "#74a9cf") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(12, 24, 12, 24)
    )
  cowplot::as_grob(plot)
}

grobs <- list(
  global = venn_page(
    bs_global,
    sbr_global,
    title = "48 h global top 500, equal cell weight",
    subtitle = sprintf(
      "%d genes shared (Jaccard %.2f)",
      length(shared_global),
      overlap$jaccard[1]
    )
  )
)
for (type_name in shared_types) {
  row <- overlap[overlap$cell_type %in% type_name, , drop = FALSE]
  grobs[[type_name]] <- venn_page(
    selected_genes(per_type, "B-S", type_name),
    selected_genes(per_type, "SBR", type_name),
    title = paste0("48 h, ", type_name, ", top 50"),
    subtitle = sprintf(
      "%d genes shared (Jaccard %.2f)",
      row$n_shared,
      row$jaccard
    )
  )
}
save_pages(
  grobs = grobs,
  filename = "venn_48h_cell_line.pdf",
  width = 8,
  height = 7
)
message("Wrote ", file.path(out_dir, "venn_48h_cell_line.pdf"))
message("Wrote ", file.path(out_dir, "overlap_48h_cell_line.csv"))
