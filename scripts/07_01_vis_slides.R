# mkdir -p logs
# Rscript --no-save --no-restore \
#   scripts/07_01_vis_slides.R \
#   > "logs/07_01_vis_slides_$(date +%F).log" 2>&1
#
# PNG figures of the deck (slides/index.qmd), on a transparent background,
# in slides/images/. The book keeps its white figures in docs/images/.
#
# No plot is defined here. Each figure script of the book is run again with
# `save_figure()` replaced by a version that keeps the figures of the deck
# only, removes the white backgrounds, and writes the PNG. Titles and
# sizes therefore follow the book.
#
# A figure script whose input is missing is skipped, and the deck keeps
# its current PNG for those figures (two inputs are DVC-tracked .rds
# files: run `dvc pull` first to regenerate the feature-selection figures).

root <- getwd()
requireNamespace("patchwork", quietly = TRUE) # the `&` operator below

deck <- list(
  list(
    script = "02_01_vis_quality_control.R",
    folder = "feature_selection",
    figures = c("fs_tier1_low_quality_cells", "fs_tier1_gene_filters")
  ),
  list(
    script = "02_02_vis_cross_cell.R",
    folder = "feature_selection",
    figures = "fs_tier2_null_calibration"
  ),
  list(
    script = "04_01_vis_deconvolution.R",
    folder = "deconvolution",
    figures = c(
      "benchmark_rmse",
      "benchmark_aitchison",
      "profile_correlation_r2_05",
      "profile_correlation_welch"
    )
  )
)

clear <- ggplot2::element_rect(fill = NA, colour = NA)
transparent <- ggplot2::theme(
  plot.background = clear,
  panel.background = clear,
  legend.background = clear,
  legend.box.background = clear,
  legend.key = clear
)

for (job in deck) {
  png_dir <- file.path(root, "slides", "images", job$folder)
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
  env <- new.env(parent = globalenv())
  # Same signature as save_figure() in R/02_00_utils.R; found first because
  # the figure script is evaluated in `env`.
  env$save_figure <- function(
    plot,
    name,
    pdf_dir,
    png_dir = NULL,
    width = 8,
    height = 6,
    dpi = 300
  ) {
    if (!name %in% job$figures) {
      return(invisible(NULL))
    }
    path <- file.path(
      root, "slides", "images", job$folder, paste0(name, ".png")
    )
    ggplot2::ggsave(
      path,
      plot & transparent,
      width = width,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      bg = "transparent"
    )
    message("Wrote ", path)
    invisible(path)
  }
  tryCatch(
    {
      sys.source(file.path(root, "scripts", job$script), envir = env)
      TRUE
    },
    error = function(e) {
      message("Skipped ", job$script, ": ", conditionMessage(e))
      FALSE
    }
  )
}

# Per-time PLS-DA and the detection UpSet. The wrapped book page is not
# cropped. Stop when the scored cells are not on disk.
slots_path <- file.path(
  root, "output", "feature_selection", "tier_01", "scored_slots.rds"
)
if (!file.exists(slots_path)) {
  stop(
    "Missing ", slots_path, ". ",
    "The per-time PLS-DA and the detection UpSet are not drawn. ",
    "The DVC cache does not hold md5 9b1b1bf40efa4c5bedc15e04706e7d41."
  )
}

source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "02_01_vis_quality_control.R"))
if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
  stop("Package ComplexUpset is required for the detection UpSet.")
}

prepared <- readRDS(slots_path)
colours <- celltype_colours(
  file.path(root, "data", "dictionaries", "suppinger_celltype_colours.csv")
)
png_dir <- file.path(root, "slides", "images", "feature_selection")
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
line_label <- c(batch1 = "B-S (batch 1)", batch2 = "SBR (batch 2)")
line_colours <- c("B-S (batch 1)" = "#1d4e89", "SBR (batch 2)" = "#c45c26")

upset_panels <- lapply(prepared$time_levels, function(tm) {
  slot <- prepared$slots[[tm]]
  rules <- gene_filter_rules(slot$counts, slot$cell_type)
  flags <- data.frame(
    `Outside the presence gate` = !rules$pass_presence,
    `Under 20% in every type` = rules$fail_rnasieve,
    check.names = FALSE
  )
  ComplexUpset::upset(
    flags,
    c("Outside the presence gate", "Under 20% in every type"),
    name = "Detection rule",
    min_size = 0,
    width_ratio = 0.25
  ) +
    ggplot2::ggtitle(tm)
})
ggplot2::ggsave(
  file.path(png_dir, "fs_tier1_detection_upset.png"),
  patchwork::wrap_plots(upset_panels, ncol = 1),
  width = 10,
  height = 4 * length(upset_panels),
  dpi = 200,
  limitsize = FALSE,
  bg = "transparent"
)

for (tm in prepared$time_levels) {
  slot <- prepared$slots[[tm]]
  gate <- prepared$presence[[tm]]
  genes <- gate$gene[gate$pass]
  groups <- if (length(unique(slot$batch)) > 1L) {
    list("Cell line" = unname(line_label[slot$batch]))
  } else {
    NULL
  }
  drawn <- plot_type_pca_plsda(
    counts = slot$counts[genes, , drop = FALSE],
    cell_type = slot$cell_type,
    size_factors = slot$library_size / stats::median(slot$library_size),
    colours = colours,
    input = "pearson",
    groups = groups,
    group_colours = list("Cell line" = line_colours)
  )
  panels <- list(drawn$pca, drawn$plsda)
  if (!is.null(drawn[["pca_Cell line"]])) {
    panels <- c(
      panels,
      list(drawn[["pca_Cell line"]], drawn[["plsda_Cell line"]])
    )
  }
  stem <- gsub("[^0-9]", "", tm)
  ggplot2::ggsave(
    file.path(png_dir, paste0("fs_tier1_plsda_", stem, "h.png")),
    patchwork::wrap_plots(panels, ncol = 2) & transparent,
    width = 11,
    height = 4.2 * length(panels) / 2,
    dpi = 200,
    limitsize = FALSE,
    bg = "transparent"
  )
}
