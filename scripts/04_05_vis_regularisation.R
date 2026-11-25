# Rscript --no-save --no-restore scripts/04_05_vis_regularisation.R
#
# Figures of the shrinkage regularisation of the per-type covariances,
# Welch gene panel, standard-Wnt protocol. The regularised covariance is
# the one of build_sigma() (R/03_01_stats_covariance_estimation.R): per
# type, the sample correlation R_j (genes with no variance get a unit
# diagonal) is replaced by (1 - lambda) R_j + lambda I.
#
# Inputs:
#   output/feature_selection/tier_01/scored_slots.rds   raw counts (DVC)
#   output/feature_selection/tier_02/selected_<time>.csv  Welch genes
#   output/deconvolution/lambda_selection.csv   bulk profile likelihood
#                                               and RMSE along lambda
#   output/grn/cv/cv_selected_tuning.csv        held-out cell choices
#   output/grn/plugin/plugin_summary.csv        cells and constant genes
#   output/grn/cv/factors_<time>_welch.rds      optional (section 7)
#
# Outputs, output/deconvolution/:
#   regularisation_spectra.csv   time, cell_type, cells, genes, lambda,
#                                kappa, eff_rank (lambda = 0 included;
#                                kappa is Inf there when R_j is singular)
#   figures/regularisation/      one SVG per page, the PDF made from it
#                                with Inkscape (pages merged with pdfunite),
#                                and an RDS of the data of each figure:
#     correlation_heatmaps_welch   3 pages, one per time
#     regularisation_spectra_welch condition number and effective rank
#     regularisation_map_welch     3 pages, one per time
#     lambda_path_welch            bulk likelihood and RMSE along lambda
#     size_term_<type>_96h_welch   only when the factor fits exist
# PNG for the book in docs/images/deconvolution/regularisation/.

# ==========================================================================
# SECTION 0 · Dependencies, settings, paths ----
# ==========================================================================

root <- here::here()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "03_01_stats_covariance_estimation.R"))

time_levels <- c("48h", "72h", "96h")
var_floor_rel <- 1e-4 # as in build_sigma()
n_heatmap_genes <- 50L
# Logit-spaced grid on [0.001, 0.999]: log-spaced near 0, where the
# sample correlation dominates, and near 1, where the levels are chosen.
lambda_grid <- stats::plogis(seq(
  stats::qlogis(1e-3),
  stats::qlogis(1 - 1e-3),
  length.out = 61L
))

slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
selected_dir <- file.path(root, "output", "feature_selection", "tier_02")
out_dir <- file.path(root, "output", "deconvolution")
fig_dir <- file.path(out_dir, "figures", "regularisation")
png_dir <- file.path(root, "docs", "images", "deconvolution", "regularisation")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
time_label <- function(x) sub("h$", " h", x)
lam <- "λ"
kap <- "κ"

# Colours: one per cell type, distinct within a time, the same for a type
# seen at two times (Okabe-Ito; the house lineage colours give three
# near-identical greens at 96 h).
type_colours <- c(
  "Naïve pluripotency" = "#CC79A7",
  "Exiting naïve pluripotency" = "#56B4E9",
  "Epiblast" = "#E69F00",
  "Ectopic pluripotency" = "#0072B2",
  "Primitive streak" = "#D55E00",
  "Caudal mesoderm" = "#009E73",
  "Gut" = "#E69F00",
  "Pre-somitic mesoderm" = "#D55E00",
  "Somite differentiation front" = "#56B4E9"
)
colour_series <- "#1d4e89"
colour_marked <- "#c45c26"
colour_ink <- "#1c2b24"

# Name Noto Sans in the SVG (same as ragg PNGs); svglite would else write Arial.
sans_family <- if ("Noto Sans" %in% systemfonts::system_fonts()$family) {
  "Noto Sans"
} else {
  systemfonts::font_info(family = "sans")$family[[1L]]
}
svg_device <- function(filename, ...) {
  svglite::svglite(filename, ..., system_fonts = list(sans = sans_family))
}

#' Save one page as SVG, then PDF through Inkscape (text as paths).
#'
#' @param plot A ggplot or patchwork object.
#' @param stem Path without extension.
#' @param width,height Size in inches.
#' @return The PDF path, invisibly.
save_page <- function(plot, stem, width, height) {
  svg <- paste0(stem, ".svg")
  pdf <- paste0(stem, ".pdf")
  ggplot2::ggsave(
    svg,
    plot,
    width = width,
    height = height,
    units = "in",
    device = svg_device,
    bg = "white"
  )
  status <- system2(
    "inkscape",
    c(
      shQuote(svg),
      "--export-text-to-path",
      "--export-type=pdf",
      "-o",
      shQuote(pdf)
    ),
    stdout = FALSE,
    stderr = FALSE
  )
  if (!identical(status, 0L) || !file.exists(pdf)) {
    stop("Inkscape failed on ", svg)
  }
  invisible(pdf)
}

#' Merge single-page PDFs with pdfunite and delete the pages.
#'
#' @param pages PDF paths, in page order.
#' @param out Path of the merged PDF.
#' @return `out`, invisibly.
merge_pages <- function(pages, out) {
  # pdfunite (poppler) ships with the system; qpdf is not installed here.
  status <- system2("pdfunite", c(shQuote(pages), shQuote(out)))
  if (!identical(status, 0L)) {
    stop("pdfunite failed on ", out)
  }
  unlink(pages)
  invisible(out)
}

#' Save the book PNG of a figure.
save_png <- function(plot, name, width, height) {
  ggplot2::ggsave(
    file.path(png_dir, paste0(name, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 200,
    device = ragg::agg_png,
    bg = "white"
  )
}

#' Non-zero eigenvalues of `Z' Z`, from the smaller Gram matrix.
#'
#' Values-only version of right_spectrum(), same zero threshold.
#'
#' @param Z Matrix, cells by genes.
#' @param tol Relative threshold under which an eigenvalue is zero.
#' @return Eigenvalues, decreasing, zeros dropped.
nonzero_values <- function(Z, tol = 1e-10) {
  gram <- if (nrow(Z) < ncol(Z)) tcrossprod(Z) else crossprod(Z)
  ev <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
  ev[ev > tol * max(ev, 1)]
}

#' Condition number and effective rank of `(1 - lambda) R + lambda I`.
#'
#' The shrunk matrix has the eigenvectors of `R` and the eigenvalues
#' `(1 - lambda) ell + lambda`. The condition number is the ratio of the
#' extreme eigenvalues; the effective rank is
#' `(sum ell)^2 / sum ell^2` of the shrunk eigenvalues (participation
#' ratio: `G` when all are equal, 1 when one dominates).
#'
#' @param ell All `G` eigenvalues of `R`, zeros included.
#' @param lambda Shrinkage levels in `[0, 1]`.
#' @return Data frame with `lambda`, `kappa`, `eff_rank`.
shrunk_conditioning <- function(ell, lambda) {
  out <- vapply(
    lambda,
    function(l) {
      s <- (1 - l) * ell + l
      c(max(s) / min(s), sum(s)^2 / sum(s^2))
    },
    numeric(2L)
  )
  data.frame(lambda = lambda, kappa = out[1L, ], eff_rank = out[2L, ])
}

# Closed form against base R, on a small full-rank and a small singular
# correlation matrix (the singular one checks the Gram path).
local({
  set.seed(20261124L)
  for (n in c(40L, 7L)) {
    X <- matrix(stats::rnorm(n * 12L), n, 12L)
    Xc <- sweep(X, 2L, colMeans(X))
    Z <- sweep(Xc, 2L, sqrt(colSums(Xc^2)), "/")
    R <- crossprod(Z)
    ell <- c(nonzero_values(Z), rep(0, 12L - min(n - 1L, 12L)))
    M <- (1 - 0.3) * R + 0.3 * diag(12L)
    direct <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    closed <- shrunk_conditioning(ell, 0.3)
    stopifnot(
      isTRUE(all.equal(nonzero_values(Z), right_spectrum(Z)$values)),
      isTRUE(all.equal(sort(direct), sort((1 - 0.3) * ell + 0.3))),
      isTRUE(all.equal(closed$kappa, kappa(M, exact = TRUE))),
      isTRUE(all.equal(closed$eff_rank, sum(direct)^2 / sum(direct^2)))
    )
  }
})

# ==========================================================================
# SECTION 1 · Tuning values ----
# ==========================================================================

selection <- read_table(file.path(out_dir, "lambda_selection.csv"))
path <- selection[
  selection$family == "shrinkage" &
    selection$arm == "standard" &
    selection$gene_set == "welch",
  c("time", "tuning", "sum_profile", "selected_profile", "rmse")
]
path <- path[order(path$time, path$tuning), ]
bulk_lambda <- with(
  path[path$selected_profile, ],
  stats::setNames(tuning, time)
)
stopifnot(setequal(names(bulk_lambda), time_levels))

cv <- read_table(file.path(
  root,
  "output",
  "grn",
  "cv",
  "cv_selected_tuning.csv"
))
cv <- cv[cv$gene_set == "welch", ]
cv_type <- cv[
  cv$rule == "shrinkage, per type",
  c("time", "cell_type", "tuning")
]
cv_common <- cv[
  cv$rule == "shrinkage, common",
  c("time", "cell_type", "tuning")
]
names(cv_type)[3L] <- "lambda_cv"
names(cv_common)[3L] <- "lambda_common"

plugin <- read_table(file.path(
  root,
  "output",
  "grn",
  "plugin",
  "plugin_summary.csv"
))
plugin <- plugin[plugin$gene_set == "welch", ]

grid_all <- sort(unique(c(
  0,
  1,
  lambda_grid,
  path$tuning,
  cv$tuning[cv$tuning > 0]
)))

# ==========================================================================
# SECTION 2 · Spectra and 50-gene correlations, one type at a time ----
# ==========================================================================

if (!file.exists(slots_path)) {
  stop(
    "Missing ",
    slots_path,
    ". From the repository root run: dvc pull ",
    "output/feature_selection/tier_01/scored_slots.rds.dvc"
  )
}
slots <- readRDS(slots_path)
requireNamespace("Matrix", quietly = TRUE) # rownames() of the sparse counts
# A local scored_slots.rds older than its DVC pointer lacks some genes of
# the current Welch panel. The run stops then, unless
# DECOVART_ALLOW_STALE_SLOTS=1 asks for provisional figures without them.
allow_stale <- identical(Sys.getenv("DECOVART_ALLOW_STALE_SLOTS"), "1")
missing_by_time <- stats::setNames(integer(length(time_levels)), time_levels)

spectra <- list()
spectrum_summary <- list()
heat <- list()
size_inputs <- list() # standardised means and correlations, section 7
for (tm in time_levels) {
  selected <- read_table(file.path(
    selected_dir,
    paste0("selected_", tm, ".csv")
  ))
  absent <- setdiff(selected$gene, rownames(slots$slots[[tm]]$counts))
  if (length(absent) > 0L) {
    msg <- paste0(
      length(absent),
      " Welch genes of ",
      tm,
      " are absent from scored_slots.rds (first: ",
      absent[[1L]],
      "): the local copy is older than its DVC pointer. Run ",
      "dvc pull output/feature_selection/tier_01/scored_slots.rds.dvc"
    )
    if (!allow_stale) {
      stop(
        msg,
        ", or set DECOVART_ALLOW_STALE_SLOTS=1 for provisional figures."
      )
    }
    warning(msg, ". Provisional run without these genes.")
    missing_by_time[[tm]] <- length(absent)
    selected <- selected[!selected$gene %in% absent, ]
  }
  panel <- read_panel_slot(slots, selected_dir, tm, genes = selected$gene)
  genes <- panel$genes
  G <- length(genes)
  top <- selected$gene[order(-selected$deviance_r2)][seq_len(n_heatmap_genes)]
  top_idx <- match(top, genes)
  types <- sort(unique(panel$cell_type[!is.na(panel$cell_type)]))
  scatter_top <- matrix(0, n_heatmap_genes, n_heatmap_genes)
  r_top <- list()
  size_inputs[[tm]] <- list()
  for (tp in types) {
    started <- proc.time()[["elapsed"]]
    rows <- which(panel$cell_type %in% tp)
    C <- length(rows)
    Xc <- panel$y[rows, , drop = FALSE]
    mu <- colMeans(Xc)
    Xc <- sweep(Xc, 2L, mu, check.margin = FALSE)
    v <- colSums(Xc^2) / (C - 1)
    floor_v <- var_floor_rel * stats::median(v[v > 0])
    sd <- sqrt(pmax(v, floor_v))
    n_const <- sum(v <= 0)
    n_floored <- sum(v > 0 & v < floor_v)
    ref <- plugin[plugin$time_point == tm & plugin$cell_type == tp, ]
    stopifnot(
      nrow(ref) == 1L,
      ref$cells == C,
      ref$genes == G + missing_by_time[[tm]],
      ref$constant_genes == n_const
    )
    scatter_top <- scatter_top + crossprod(Xc[, top_idx, drop = FALSE])
    # Z' Z is R_j of build_sigma() before its diagonal is set to 1.
    Z <- sweep(Xc, 2L, sd * sqrt(C - 1), "/", check.margin = FALSE)
    rm(Xc)
    r_top[[tp]] <- crossprod(Z[, top_idx, drop = FALSE])
    size_inputs[[tm]][[tp]] <- list(
      r_top = r_top[[tp]],
      m_top = (mu / sd)[top_idx],
      m_all = mu / sd
    )
    # A constant gene (zero column of Z) or a gene under the variance
    # floor has a diagonal below 1 in Z' Z, which build_sigma() resets to
    # 1. One extra row sqrt(1 - |z_g|^2) e_g' per such gene adds exactly
    # that deficit, so the Gram spectrum of the augmented Z is the
    # spectrum of R_j.
    deficit <- 1 - colSums(Z^2)
    fix <- which(deficit > 1e-8)
    if (length(fix) > 0L) {
      E <- matrix(0, length(fix), G)
      E[cbind(seq_along(fix), fix)] <- sqrt(deficit[fix])
      Z <- rbind(Z, E)
    }
    vals <- nonzero_values(Z)
    rm(Z)
    gc(verbose = FALSE)
    ell <- c(vals, rep(0, G - length(vals)))
    stopifnot(abs(sum(ell) - G) < 1e-6 * G)
    cond <- shrunk_conditioning(ell, grid_all)
    spectra[[length(spectra) + 1L]] <- data.frame(
      time = tm,
      cell_type = tp,
      cells = C,
      genes = G,
      genes_missing = missing_by_time[[tm]],
      cond
    )
    at <- function(l) cond[match(l, cond$lambda), ]
    l_cv <- cv_type$lambda_cv[cv_type$time == tm & cv_type$cell_type == tp]
    spectrum_summary[[length(spectrum_summary) + 1L]] <- data.frame(
      time = tm,
      cell_type = tp,
      cells = C,
      genes = G,
      rank = length(vals),
      floored_genes = n_floored,
      kappa_min_grid = at(min(lambda_grid))$kappa,
      lambda_bulk = bulk_lambda[[tm]],
      kappa_bulk = at(bulk_lambda[[tm]])$kappa,
      lambda_cv = l_cv,
      kappa_cv = at(l_cv)$kappa,
      eff_rank_0 = at(0)$eff_rank,
      eff_rank_bulk = at(bulk_lambda[[tm]])$eff_rank,
      eff_rank_cv = at(l_cv)$eff_rank
    )
    message(
      tm,
      " / ",
      tp,
      ": ",
      C,
      " cells, ",
      G,
      " genes, rank ",
      length(vals),
      ", ",
      round(proc.time()[["elapsed"]] - started, 1),
      " s"
    )
  }
  # Gene order of the page: pooled within-type correlation, 1 - |r|.
  r_pool <- pooled_correlation(scatter_top)
  dimnames(r_pool) <- list(top, top)
  gene_order <- top[
    stats::hclust(stats::as.dist(1 - abs(r_pool)), method = "average")$order
  ]
  heat[[tm]] <- list(r_top = r_top, genes = top, gene_order = gene_order)
  rm(panel)
  gc(verbose = FALSE)
}
rm(slots)
gc(verbose = FALSE)

spectra <- do.call(rbind, spectra)
spectrum_summary <- do.call(rbind, spectrum_summary)
utils::write.csv(
  spectra,
  file.path(out_dir, "regularisation_spectra.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
at_one <- spectra[spectra$lambda == 1, ]
stopifnot(
  all(abs(at_one$kappa - 1) < 1e-10),
  all(abs(at_one$eff_rank - at_one$genes) < 1e-6 * at_one$genes)
)
print(spectrum_summary, digits = 4)
# Appended to every subtitle of a provisional run.
provisional <- if (any(missing_by_time > 0L)) {
  paste0(
    "\nProvisional: ",
    paste0(missing_by_time, collapse = " / "),
    " Welch genes at 48 / 72 / 96 h are absent from the local ",
    "scored_slots.rds and left out."
  )
} else {
  ""
}

# Labels of the types: name and cells, rows ordered by cells.
types_tab <- unique(spectra[, c("time", "cell_type", "cells", "genes")])
types_tab <- merge(types_tab, cv_type, by = c("time", "cell_type"))
types_tab <- merge(types_tab, cv_common, by = c("time", "cell_type"))
types_tab$lambda_bulk <- bulk_lambda[types_tab$time]
types_tab <- types_tab[order(types_tab$time, -types_tab$cells), ]

# ==========================================================================
# SECTION 3 · (i) Correlation heatmaps, top 50 Welch genes ----
# ==========================================================================

column_levels <- function(tm) {
  c(
    "Sample correlation",
    paste0("Bulk likelihood, ", lam, " = ", fmt(bulk_lambda[[tm]], 2)),
    paste0("Held-out cells, ", lam, " per type")
  )
}
heat_long <- do.call(
  rbind,
  lapply(time_levels, function(tm) {
    h <- heat[[tm]]
    tt <- types_tab[types_tab$time == tm, ]
    cols <- column_levels(tm)
    do.call(
      rbind,
      lapply(seq_len(nrow(tt)), function(k) {
        tp <- tt$cell_type[[k]]
        r <- h$r_top[[tp]]
        diag(r) <- NA
        factors <- c(1, 1 - tt$lambda_bulk[[k]], 1 - tt$lambda_cv[[k]])
        do.call(
          rbind,
          lapply(1:3, function(m) {
            data.frame(
              time = tm,
              cell_type = tp,
              row_label = paste0(
                tp,
                "\nn = ",
                fmt_n(tt$cells[[k]]),
                "\nheld-out ",
                lam,
                " = ",
                fmt(tt$lambda_cv[[k]])
              ),
              column = cols[[m]],
              gene_x = rep(h$genes, times = length(h$genes)),
              gene_y = rep(h$genes, each = length(h$genes)),
              value = factors[[m]] * as.vector(r),
              stringsAsFactors = FALSE
            )
          })
        )
      })
    )
  })
)

#' Heatmap page of one time.
plot_heatmap_page <- function(d, tm) {
  h <- heat[[tm]]
  tt <- types_tab[types_tab$time == tm, ]
  d$column <- factor(d$column, levels = column_levels(tm))
  d$row_label <- factor(d$row_label, levels = unique(d$row_label))
  d$gene_x <- factor(d$gene_x, levels = h$gene_order)
  d$gene_y <- factor(d$gene_y, levels = rev(h$gene_order))
  biggest <- tt[1L, ]
  ggplot2::ggplot(d, ggplot2::aes(gene_x, gene_y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(row_label ~ column) +
    ggplot2::scale_fill_gradientn(
      colours = colorspace::diverging_hcl(11, palette = "Blue-Red 3"),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.5),
      na.value = "grey80",
      name = "Correlation",
      guide = ggplot2::guide_colourbar(
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    ggplot2::scale_x_discrete(labels = NULL, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
      title = paste0(
        "Shrinkage multiplies every off-diagonal correlation by 1 − ",
        lam,
        ": the pattern stays, its contrast fades"
      ),
      subtitle = paste0(
        "Welch panel, ",
        time_label(tm),
        ", standard Wnt: the ",
        n_heatmap_genes,
        " genes of highest deviance share, ordered by clustering of the ",
        "pooled within-type correlation (average linkage on 1 − |r|).\n",
        "Diagonal (grey) set to 1. Columns: sample correlation, then the ",
        "levels chosen by the bulk likelihood (common) and by held-out cells ",
        "(per type, in the row label).",
        provisional
      ),
      x = NULL,
      y = NULL
    ) +
    theme_slide(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.5, "lines"),
      axis.text.y = ggplot2::element_text(size = 4.5, colour = colour_ink),
      strip.text.y.right = ggplot2::element_text(angle = 0, hjust = 0),
      legend.key.width = ggplot2::unit(2.5, "cm"),
      legend.key.height = ggplot2::unit(0.3, "cm")
    )
}

heat_pages <- character()
heat_size <- function(tm) {
  n <- sum(types_tab$time == tm)
  c(width = 13, height = 2.4 + 3.3 * n)
}
for (tm in time_levels) {
  p <- plot_heatmap_page(heat_long[heat_long$time == tm, ], tm)
  size <- heat_size(tm)
  heat_pages[[tm]] <- save_page(
    p,
    file.path(fig_dir, paste0("correlation_heatmaps_welch_", tm)),
    size[["width"]],
    size[["height"]]
  )
  if (tm == "96h") {
    save_png(
      p,
      "correlation_heatmap_96h_welch",
      size[["width"]],
      size[["height"]]
    )
  }
}
merge_pages(heat_pages, file.path(fig_dir, "correlation_heatmaps_welch.pdf"))
saveRDS(
  list(
    data = heat_long,
    gene_order = lapply(heat, `[[`, "gene_order"),
    types = types_tab
  ),
  file.path(fig_dir, "correlation_heatmaps_welch.rds")
)

# ==========================================================================
# SECTION 4 · (ii) Conditioning along the shrinkage level ----
# ==========================================================================

# Plotted on a logit axis, so lambda = 0 and 1 (in the CSV) are left out;
# at lambda = 1 every line ends at kappa = 1 and r_eff = G.
spec_plot <- spectra[spectra$lambda > 0 & spectra$lambda < 1, ]
spec_plot$time_label <- factor(
  paste0(time_label(spec_plot$time), " · ", fmt_n(spec_plot$genes), " genes"),
  levels = unique(paste0(
    time_label(types_tab$time),
    " · ",
    fmt_n(types_tab$genes),
    " genes"
  ))
)
spec_plot$eff_share <- spec_plot$eff_rank / spec_plot$genes
spec_bulk <- unique(spec_plot[, c("time", "time_label")])
spec_bulk$lambda <- bulk_lambda[spec_bulk$time]
spec_cv <- merge(
  spec_plot,
  cv_type,
  by.x = c("time", "cell_type", "lambda"),
  by.y = c("time", "cell_type", "lambda_cv")
)
# Direct labels in the empty lower-left corner of the kappa panels,
# stacked in the order of the lines at lambda = 0.001 (the lines meet
# at the right end, so labels there would collide).
spec_left <- spec_plot[spec_plot$lambda == min(lambda_grid), ]
spec_left <- spec_left[order(spec_left$time, -spec_left$kappa), ]
spec_left$y_kappa <- 10^stats::ave(
  spec_left$kappa,
  spec_left$time,
  FUN = function(z) log10(80) - 0.38 * (seq_along(z) - 1)
)
x_label_pos <- 1.3e-3

kappa_max <- spectrum_summary[which.max(spectrum_summary$kappa_min_grid), ]
kappa_bulk_max <- max(spectrum_summary$kappa_bulk)
lambda_breaks <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
spectra_plot <- function() {
  base <- function(y, d_left = NULL) {
    p <- ggplot2::ggplot(
      spec_plot,
      ggplot2::aes(lambda, .data[[y]], colour = cell_type)
    ) +
      ggplot2::geom_vline(
        data = spec_bulk,
        ggplot2::aes(xintercept = lambda),
        linetype = "dashed",
        colour = colour_marked
      ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(data = spec_cv, size = 2.4) +
      ggplot2::facet_wrap(~time_label, nrow = 1L) +
      ggplot2::scale_colour_manual(values = type_colours, guide = "none") +
      ggplot2::scale_x_continuous(
        transform = "logit",
        breaks = lambda_breaks,
        labels = format(lambda_breaks, drop0trailing = TRUE)
      ) +
      theme_slide(base_size = 12) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
    if (!is.null(d_left)) {
      p <- p +
        ggplot2::geom_label(
          data = d_left,
          ggplot2::aes(x = x_label_pos, y = y_kappa, label = cell_type),
          hjust = 0,
          size = 3.4,
          fill = "white",
          linewidth = 0,
          label.padding = ggplot2::unit(0.1, "lines")
        )
    }
    p
  }
  p_kappa <- base("kappa", spec_left) +
    ggplot2::scale_y_log10(
      labels = scales::label_number(big.mark = ",", accuracy = 1),
      breaks = 10^(0:7)
    ) +
    ggplot2::labs(
      x = NULL,
      y = paste0("Condition number ", kap, " (log scale)")
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  p_rank <- base("eff_share") +
    ggplot2::scale_y_log10(
      breaks = c(0.001, 0.01, 0.1, 1),
      labels = c("0.1%", "1%", "10%", "100%")
    ) +
    ggplot2::labs(
      x = paste0(
        "Shrinkage level ",
        lam,
        " (logit scale; 0 = sample correlation, 1 = diagonal)"
      ),
      y = "Effective rank, share of the genes (log scale)"
    ) +
    ggplot2::theme(strip.text = ggplot2::element_blank())
  patchwork::wrap_plots(p_kappa, p_rank, ncol = 1L) +
    patchwork::plot_annotation(
      title = paste0(
        "Shrinkage brings ",
        kap,
        " from up to ",
        fmt(kappa_max$kappa_min_grid / 1e6, 1),
        " million (",
        lam,
        " = 0.001) to at most ",
        fmt_n(kappa_bulk_max),
        " at the level chosen by the bulk likelihood"
      ),
      subtitle = paste0(
        "Welch panel, standard Wnt. Lines: (1 − ",
        lam,
        ") R + ",
        lam,
        " I for each cell type; dashed: level chosen by the bulk ",
        "likelihood; point: level chosen per type by held-out cells.\n",
        "Effective rank: (Σ ℓ)² / Σ ℓ² of the ",
        "eigenvalues ℓ, as a share of the genes.",
        provisional
      ),
      theme = theme_slide(base_size = 12)
    )
}
p <- spectra_plot()
save_page(p, file.path(fig_dir, "regularisation_spectra_welch"), 15, 9.5)
save_png(p, "regularisation_spectra_welch", 15, 9.5)
saveRDS(
  list(
    data = spec_plot,
    bulk = spec_bulk,
    cv = spec_cv,
    labels = spec_left,
    summary = spectrum_summary
  ),
  file.path(fig_dir, "regularisation_spectra_welch.rds")
)

# ==========================================================================
# SECTION 5 · (iii) Regularisation map ----
# ==========================================================================

map_d <- types_tab
map_d$time_label <- time_label(map_d$time)
map_d$row_label <- paste0(map_d$cell_type, "  (n = ", fmt_n(map_d$cells), ")")
map_d$row_label <- factor(map_d$row_label, levels = rev(map_d$row_label))
map_top <- do.call(
  rbind,
  lapply(split(map_d, map_d$time), function(d) d[which.max(d$cells), ])
)
map_lines <- unique(map_d[, c("time", "time_label", "lambda_bulk")])
map_common <- unique(map_d[, c("time", "time_label", "lambda_common")])

plot_map <- function(d, top, lines, common, facet = FALSE) {
  p <- ggplot2::ggplot(d, ggplot2::aes(y = row_label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = lambda_cv, xend = lambda_bulk, yend = row_label),
      colour = "grey75",
      linewidth = 1.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = lambda_common),
      shape = 21,
      fill = NA,
      colour = colour_series,
      size = 5,
      stroke = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = lambda_cv),
      colour = colour_series,
      size = 3.2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = lambda_bulk),
      colour = colour_marked,
      size = 3.2
    ) +
    ggplot2::geom_text(
      data = lines,
      ggplot2::aes(
        x = lambda_bulk,
        y = Inf,
        label = paste0("bulk likelihood, ", lam, " = ", fmt(lambda_bulk, 2))
      ),
      inherit.aes = FALSE,
      colour = colour_marked,
      hjust = 0.5,
      vjust = 1.3,
      size = 3.6
    ) +
    ggplot2::geom_text(
      data = common,
      ggplot2::aes(
        x = lambda_common,
        y = Inf,
        label = paste0(
          "○ held-out cells, common ",
          lam,
          " = ",
          fmt(lambda_common)
        )
      ),
      inherit.aes = FALSE,
      colour = colour_series,
      hjust = 1,
      vjust = 1.3,
      size = 3.6
    ) +
    ggplot2::geom_text(
      data = top,
      ggplot2::aes(x = lambda_cv, label = "● held-out cells, per type"),
      colour = colour_series,
      hjust = 1,
      vjust = -1.1,
      size = 3.6
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = ggplot2::expansion(add = 0.02)
    ) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.6, 1.4))) +
    ggplot2::labs(
      x = paste0(
        "Shrinkage level ",
        lam,
        " (0 = sample correlation, 1 = diagonal)"
      ),
      y = NULL
    ) +
    theme_slide(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  if (facet) {
    p <- p +
      ggplot2::facet_grid(time_label ~ ., scales = "free_y", space = "free_y") +
      ggplot2::theme(strip.text.y.right = ggplot2::element_text(angle = 0))
  }
  p
}
map_title <- paste0(
  "Held-out cells choose heavier shrinkage (",
  lam,
  " ",
  fmt(min(map_d$lambda_cv), 2),
  "–",
  fmt(max(map_d$lambda_cv), 2),
  ") than the bulk likelihood (",
  fmt(min(bulk_lambda), 2),
  "–",
  fmt(max(bulk_lambda), 2),
  ")"
)
map_subtitle <- paste0(
  "Welch panel, standard Wnt. Rows: cell types, ordered by cells. ",
  "Filled blue: level chosen per type by held-out cells; open blue: one ",
  "level for all types;\norange: level chosen by the summed bulk profile ",
  "likelihood, common to all types.",
  provisional
)
map_pages <- character()
for (tm in time_levels) {
  n <- sum(map_d$time == tm)
  p <- plot_map(
    map_d[map_d$time == tm, ],
    map_top[map_top$time == tm, ],
    map_lines[map_lines$time == tm, ],
    map_common[map_common$time == tm, ]
  ) +
    ggplot2::labs(
      title = map_title,
      subtitle = paste0(time_label(tm), ". ", map_subtitle)
    )
  map_pages[[tm]] <- save_page(
    p,
    file.path(fig_dir, paste0("regularisation_map_welch_", tm)),
    12,
    2.6 + 0.6 * n
  )
}
merge_pages(map_pages, file.path(fig_dir, "regularisation_map_welch.pdf"))
p <- plot_map(
  map_d,
  map_top[map_top$time == "48h", ],
  map_lines,
  map_common,
  facet = TRUE
) +
  ggplot2::labs(title = map_title, subtitle = map_subtitle)
save_png(p, "regularisation_map_welch", 12, 8)
saveRDS(
  list(data = map_d, bulk = map_lines, common = map_common),
  file.path(fig_dir, "regularisation_map_welch.rds")
)

# ==========================================================================
# SECTION 6 · (v) Bulk likelihood and RMSE along lambda ----
# ==========================================================================

# A step between two grid values is "decoupled" when the bulk likelihood
# rises while the composition error rises too.
steps <- do.call(
  rbind,
  lapply(split(path, path$time), function(d) {
    k <- seq_len(nrow(d) - 1L)
    data.frame(
      time = d$time[k],
      xmin = d$tuning[k],
      xmax = d$tuning[k + 1L],
      d_loglik = diff(d$sum_profile),
      d_rmse = diff(d$rmse)
    )
  })
)
steps$decoupled <- steps$d_loglik > 0 & steps$d_rmse > 0
# Merge consecutive decoupled steps into intervals.
decoupled <- do.call(
  rbind,
  lapply(split(steps, steps$time), function(d) {
    d <- d[d$decoupled, ]
    if (nrow(d) == 0L) {
      return(NULL)
    }
    block <- cumsum(c(TRUE, d$xmin[-1L] != d$xmax[-nrow(d)]))
    do.call(
      rbind,
      lapply(split(d, block), function(b) {
        data.frame(
          time = b$time[[1L]],
          xmin = min(b$xmin),
          xmax = max(b$xmax),
          d_loglik = sum(b$d_loglik),
          d_rmse = sum(b$d_rmse)
        )
      })
    )
  })
)
print(steps, digits = 4)
print(decoupled, digits = 4)

quantity_levels <- c(
  "Summed profile log-likelihood, 3 replicates",
  "RMSE of the proportions"
)
path_long <- rbind(
  data.frame(
    path[, c("time", "tuning", "selected_profile")],
    quantity = quantity_levels[[1L]],
    value = path$sum_profile
  ),
  data.frame(
    path[, c("time", "tuning", "selected_profile")],
    quantity = quantity_levels[[2L]],
    value = path$rmse
  )
)
path_long$quantity <- factor(path_long$quantity, levels = quantity_levels)
path_long$time_label <- time_label(path_long$time)
if (!is.null(decoupled)) {
  decoupled$time_label <- time_label(decoupled$time)
}
widest <- if (is.null(decoupled)) {
  NULL
} else {
  decoupled[which.max(decoupled$xmax - decoupled$xmin), ]
}

#' One panel of the path figure: one quantity at one time, own y scale.
#' Separate plots (not free facets) so that patchwork aligns the panels.
plot_path_panel <- function(q, tm, y_label, top_row) {
  d <- path_long[path_long$quantity == q & path_long$time == tm, ]
  chosen <- d[d$selected_profile, ]
  band <- if (is.null(decoupled)) NULL else decoupled[decoupled$time == tm, ]
  p <- ggplot2::ggplot(d, ggplot2::aes(tuning, value))
  if (!is.null(band) && nrow(band) > 0L) {
    p <- p +
      ggplot2::geom_rect(
        data = band,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "#e6e1d8",
        alpha = 0.8
      )
    if (top_row) {
      p <- p +
        ggplot2::geom_text(
          data = band,
          ggplot2::aes(
            x = (xmin + xmax) / 2,
            y = -Inf,
            label = paste0(
              "decoupled\nRMSE +",
              formatC(d_rmse, format = "fg", digits = 2)
            )
          ),
          inherit.aes = FALSE,
          vjust = -0.8,
          colour = colour_ink,
          size = 3.6
        )
    }
  }
  p <- p +
    ggplot2::geom_vline(
      data = chosen,
      ggplot2::aes(xintercept = tuning),
      linetype = "dashed",
      colour = colour_marked
    ) +
    ggplot2::geom_line(colour = colour_series) +
    ggplot2::geom_point(colour = colour_series, size = 1.8) +
    ggplot2::geom_point(data = chosen, colour = colour_marked, size = 3) +
    ggplot2::scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 0.75, 1)) +
    ggplot2::labs(x = NULL, y = y_label) +
    theme_slide(base_size = 12)
  if (top_row) {
    p <- p +
      ggplot2::facet_wrap(~time_label) +
      ggplot2::scale_y_continuous(labels = scales::label_comma()) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  } else {
    p <- p + ggplot2::expand_limits(y = 0)
    if (diff(range(d$value)) < 1e-12) {
      p <- p +
        ggplot2::annotate(
          "text",
          x = mean(range(d$tuning)),
          y = d$value[[1L]] / 2,
          label = paste0("Same RMSE at every ", lam),
          colour = colour_ink,
          size = 3.6
        )
    }
  }
  p
}

plot_path <- function() {
  panels <- c(
    lapply(time_levels, function(tm) {
      plot_path_panel(
        quantity_levels[[1L]],
        tm,
        if (tm == time_levels[[1L]]) {
          "Summed profile log-likelihood\n(3 replicates)"
        } else {
          NULL
        },
        top_row = TRUE
      )
    }),
    lapply(time_levels, function(tm) {
      plot_path_panel(
        quantity_levels[[2L]],
        tm,
        if (tm == time_levels[[1L]]) "RMSE of the proportions" else NULL,
        top_row = FALSE
      ) +
        ggplot2::labs(
          x = if (tm == time_levels[[2L]]) {
            paste0("Shrinkage level ", lam)
          } else {
            NULL
          }
        )
    })
  )
  patchwork::wrap_plots(panels, ncol = 3L) +
    patchwork::plot_annotation(
      title = if (is.null(widest)) {
        "The bulk likelihood and the composition error move together"
      } else {
        paste0(
          "At ",
          time_label(widest$time),
          " the bulk likelihood rises from ",
          lam,
          " = ",
          widest$xmin,
          " to ",
          widest$xmax,
          " while the RMSE rises from ",
          fmt(path$rmse[path$time == widest$time & path$tuning == widest$xmin]),
          " to ",
          fmt(path$rmse[path$time == widest$time & path$tuning == widest$xmax])
        )
      },
      subtitle = paste0(
        "Welch panel, standard Wnt, shrinkage of the raw covariance with free ",
        "scales. Orange: level chosen by the bulk likelihood. ",
        "RMSE is not used by the selection.\nShaded: steps where the ",
        "likelihood rises and the RMSE rises too (decoupled). ",
        "Each panel has its own vertical scale.",
        provisional
      ),
      theme = theme_slide(base_size = 12)
    )
}
p <- plot_path()
save_page(p, file.path(fig_dir, "lambda_path_welch"), 13, 7.5)
save_png(p, "lambda_path_welch", 13, 7.5)
saveRDS(
  list(data = path_long, steps = steps, decoupled = decoupled),
  file.path(fig_dir, "lambda_path_welch.rds")
)

# ==========================================================================
# SECTION 7 · (iv) Size term and factor fit at 96 h (optional) ----
# ==========================================================================

# For the largest type at 96 h, on the same 50 genes: the sample
# correlation R, a rank-one size term nu u u' with u the standardised mean
# (mu / sd, unit length) and nu = u' R u (least-squares fit of R by a
# multiple of u u'), the remainder R - nu u u', and the factor fit
# B B' + Psi at the rank chosen on held-out cells ("factor, per type";
# largest fitted rank not above it), with the alignment A_k of each
# loading with the standardised mean over all genes.
factor_path <- file.path(root, "output", "grn", "cv", "factors_96h_welch.rds")
if (!file.exists(factor_path)) {
  message(
    "Section 7 skipped: ",
    factor_path,
    " does not exist (factor fits of scripts/03_04_cv_covariance.R)."
  )
} else {
  fits <- readRDS(factor_path)
  tt <- types_tab[types_tab$time == "96h", ]
  tp <- tt$cell_type[which.max(tt$cells)]
  inp <- size_inputs[["96h"]][[tp]]
  genes96 <- read_table(file.path(selected_dir, "selected_96h.csv"))$gene
  top <- heat[["96h"]]$genes
  u <- inp$m_top / sqrt(sum(inp$m_top^2))
  nu <- drop(crossprod(u, inp$r_top %*% u))
  size_term <- nu * tcrossprod(u)
  ranks <- as.integer(names(fits$per_type[[tp]]))
  cv_rank <- read_table(file.path(
    root,
    "output",
    "grn",
    "cv",
    "cv_selected_tuning.csv"
  ))
  cv_rank <- cv_rank$tuning[
    cv_rank$time == "96h" &
      cv_rank$gene_set == "welch" &
      cv_rank$rule == "factor, per type" &
      cv_rank$cell_type == tp
  ]
  stopifnot(length(cv_rank) == 1L)
  rank_q <- as.character(max(ranks[ranks <= cv_rank]))
  f <- fits$per_type[[tp]][[rank_q]]
  # A_k = |b_k' m| / (|b_k| |m|), m the standardised mean (genes with a
  # loading only, as fit_factor_types()).
  loaded <- rowSums(f$B != 0) > 0
  m_all <- inp$m_all[loaded]
  align <- abs(drop(crossprod(f$B[loaded, , drop = FALSE], m_all))) /
    (sqrt(colSums(f$B[loaded, , drop = FALSE]^2)) * sqrt(sum(m_all^2)))
  idx <- match(top, genes96)
  # A factor file of an older panel would not hold these genes.
  stopifnot(!anyNA(idx), nrow(f$B) == length(genes96))
  factor_fit <- tcrossprod(f$B[idx, , drop = FALSE]) + diag(f$psi[idx])
  panels <- stats::setNames(
    list(inp$r_top, size_term, inp$r_top - size_term, factor_fit),
    c(
      "Sample correlation",
      "Size term ν u u′",
      "Remainder R − ν u u′",
      paste0("Factor fit B B′ + Ψ, rank ", rank_q)
    )
  )
  size_long <- do.call(
    rbind,
    lapply(names(panels), function(nm) {
      m <- panels[[nm]]
      diag(m) <- NA
      data.frame(
        panel = nm,
        gene_x = rep(top, times = length(top)),
        gene_y = rep(top, each = length(top)),
        value = as.vector(m)
      )
    })
  )
  size_long$panel <- factor(size_long$panel, levels = names(panels))
  size_long$gene_x <- factor(
    size_long$gene_x,
    levels = heat[["96h"]]$gene_order
  )
  size_long$gene_y <- factor(
    size_long$gene_y,
    levels = rev(heat[["96h"]]$gene_order)
  )
  p <- ggplot2::ggplot(size_long, ggplot2::aes(gene_x, gene_y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~panel, nrow = 1L) +
    ggplot2::scale_fill_gradientn(
      colours = colorspace::diverging_hcl(11, palette = "Blue-Red 3"),
      limits = c(-1, 1),
      oob = scales::squish,
      na.value = "grey80",
      name = "Correlation"
    ) +
    ggplot2::scale_x_discrete(labels = NULL, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
      title = paste0(
        "A rank-one size term carries part of the correlation of ",
        tp,
        " at 96 h"
      ),
      subtitle = paste0(
        "Welch panel, ",
        n_heatmap_genes,
        " genes of highest deviance share; u: standardised mean, ",
        "ν = u′ R u = ",
        fmt(nu, 2),
        "; alignment of the loadings with the mean, A = ",
        paste(fmt(sort(align, decreasing = TRUE), 2), collapse = ", "),
        ". Diagonal in grey.",
        provisional
      ),
      x = NULL,
      y = NULL
    ) +
    theme_slide(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 4.5, colour = colour_ink),
      legend.key.width = ggplot2::unit(2.5, "cm"),
      legend.key.height = ggplot2::unit(0.3, "cm")
    )
  stem <- "size_term_96h_welch"
  save_page(p, file.path(fig_dir, stem), 16, 5.5)
  save_png(p, stem, 16, 5.5)
  saveRDS(
    list(
      data = size_long,
      cell_type = tp,
      nu = nu,
      rank = rank_q,
      alignment = align
    ),
    file.path(fig_dir, paste0(stem, ".rds"))
  )
}

writeLines(
  capture.output(utils::sessionInfo()),
  file.path(fig_dir, "session_info.txt")
)
