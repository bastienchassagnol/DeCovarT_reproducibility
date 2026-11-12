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
