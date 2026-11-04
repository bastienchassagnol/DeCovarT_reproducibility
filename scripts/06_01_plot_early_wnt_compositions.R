# mkdir -p logs
# Rscript --no-save --no-restore \
#   scripts/06_01_plot_early_wnt_compositions.R \
#   > "logs/06_01_early_wnt_$(date +%F).log" 2>&1
#
# Written, not run. The input table does not exist yet: it is the output
# of the comparator runs planned in docs/04-deconvolution-pipeline.qmd
# (#sec-deconv-models). Once it exists, this script writes three vector
# panels for the TikZ composite docs/images/src/early_wnt/early_wnt.tex,
# which replaces docs/images/deconvolution_guided_early_treatment-1.png
# in docs/05-biological-insights.qmd (#fig-bio-early).
#
# Input: output/deconvolution/proportions_long.csv, UTF-8, one row per
# method x bulk sample x cell type.
#
#   method          chr  comparator (BayesPrism, SCDC, DWLS, DeCovarT, ...)
#   sample_id       chr  bulk column of GSE229386_AllHTSeqCountsWithGeneNames
#   arm             chr  "standard" (Wnt pulse at 48 h) or "early" (24 h)
#   time            chr  bulk time point: "48h", "72h" or "96h"
#   replicate       int  biological replicate within arm x time (1 to 3)
#   reference_slot  chr  single-cell slot whose reference was used
#   cell_type       chr  celltypeannotation, as in
#                        data/dictionaries/suppinger_celltype_crosswalk.csv
#   proportion      dbl  estimated proportion in [0, 1]; sums to one over
#                        cell_type within method x sample_id
#   lower, upper    dbl  optional; the method's own interval for that
#                        proportion (e.g. BayesPrism posterior), else NA
#
# Outputs:
#   docs/images/src/early_wnt/panel_b_levels.pdf     proportions per method
#   docs/images/src/early_wnt/panel_c_contrasts.pdf  early minus standard
#   docs/images/src/early_wnt/panel_d_agreement.pdf  sign agreement
#   output/deconvolution/early_wnt/contrasts.csv
#   output/deconvolution/early_wnt/agreement.csv
#   output/deconvolution/early_wnt/session_info.txt
#
# Statistical choices, fixed before any estimate is seen:
#   - Contrasts are on the centred log-ratio (CLR) scale within each
#     sample, because proportions are closed: one type going up forces
#     others down. A pseudo-count is added before the log; rerun with
#     pseudo_count = 1e-4 as a sensitivity check.
#   - The unit of the arm effect is the bulk biological replicate (three
#     per arm and time). Intervals resample replicates within arm x time,
#     independently in the two arms. With three replicates a resample has
#     only 10 distinct multisets per arm, so the percentile interval is
#     coarse; the nine pairwise early-minus-standard differences are drawn
#     beside it.
#   - Both arms at one time must be deconvolved with the same reference;
#     the script stops otherwise.

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

# Run from the repository root. Paths are relative to that root.

stopifnot(requireNamespace("ggplot2", quietly = TRUE))

root <- getwd()
in_file <- file.path(root, "output", "deconvolution", "proportions_long.csv")
out_dir <- file.path(root, "output", "deconvolution", "early_wnt")
fig_dir <- file.path(root, "docs", "images", "src", "early_wnt")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20261101L
n_boot <- 2000L
pseudo_count <- 1e-3
sum_tolerance <- 1e-6

arm_levels <- c("standard", "early")
arm_labels <- c(
  standard = "Standard Wnt (pulse at 48 h)",
  early = "Early Wnt (pulse at 24 h)"
)
arm_colours <- c(standard = "#1d4e89", early = "#c45c26")
time_levels <- c("48h", "72h", "96h")

# Developmental order of the scored types (docs/02, #tbl-fs-cell-counts).
# Any other label is appended in alphabetical order.
type_order <- c(
  "Naïve pluripotency",
  "Exiting naïve pluripotency",
  "Epiblast",
  "Ectopic pluripotency",
  "Primitive streak",
  "Caudal mesoderm",
  "Pre-somitic mesoderm",
  "Somite differentiation front",
  "Gut"
)

ink <- "#1c2b24"
neutral <- "#d9d3c7"

pretty_time <- function(x) sub("h$", " h", x)

# theme_slide() of the r-development rule (the repository has no theme of
# its own), with the smaller type that the rule asks for in multi-panel
# figures.
theme_slide <- function() {
  ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e1d8"),
      plot.title = ggplot2::element_text(
        face = "bold",
        colour = "#0e4d3a",
        size = 18
      ),
      plot.subtitle = ggplot2::element_text(colour = "#3d4a42", size = 13),
      axis.title = ggplot2::element_text(colour = "#1c2b24"),
      strip.text = ggplot2::element_text(face = "bold", colour = "#0e4d3a"),
      legend.position = "bottom",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

theme_panel <- function() {
  theme_slide() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10, hjust = 0),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11)
    )
}

# ==========================================================================
# SECTION 1 · Read and validate the proportions ----
# ==========================================================================

read_proportions <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing ",
      path,
      ". Run the comparator calls of docs/04-deconvolution-pipeline.qmd ",
      "and write one row per method x bulk sample x cell type.",
      call. = FALSE
    )
  }
  x <- utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  required <- c(
    "method",
    "sample_id",
    "arm",
    "time",
    "replicate",
    "reference_slot",
    "cell_type",
    "proportion"
  )
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0L) {
    stop(
      "Columns missing from ",
      path,
      ": ",
      toString(missing_cols),
      ".",
      call. = FALSE
    )
  }
  if (!all(x$arm %in% arm_levels)) {
    stop("`arm` must be one of ", toString(arm_levels), ".", call. = FALSE)
  }
  if (!all(x$time %in% time_levels)) {
    stop("`time` must be one of ", toString(time_levels), ".", call. = FALSE)
  }
  if (anyNA(x$proportion) || any(x$proportion < 0 | x$proportion > 1)) {
    stop(
      "`proportion` must lie in [0, 1] and have no missing value.",
      call. = FALSE
    )
  }

  sums <- stats::aggregate(
    proportion ~ method + sample_id,
    data = x,
    FUN = sum
  )
  off <- abs(sums$proportion - 1) > sum_tolerance
  if (any(off)) {
    stop(
      "Proportions do not sum to one for ",
      sum(off),
      " method x sample pairs, for example ",
      sums$method[off][1L],
      " on ",
      sums$sample_id[off][1L],
      ".",
      call. = FALSE
    )
  }

  refs <- unique(x[, c("method", "time", "arm", "reference_slot")])
  n_refs <- stats::aggregate(
    reference_slot ~ method + time,
    data = refs,
    FUN = function(v) length(unique(v))
  )
  if (any(n_refs$reference_slot > 1L)) {
    stop(
      "Both arms at one time must be deconvolved with the same single-cell ",
      "reference; otherwise the early-versus-standard contrast mixes ",
      "references.",
      call. = FALSE
    )
  }

  reps <- unique(x[, c("method", "arm", "time", "replicate")])
  n_reps <- stats::aggregate(
    replicate ~ method + arm + time,
    data = reps,
    FUN = length
  )
  if (any(n_reps$replicate != 3L)) {
    warning(
      "Expected three bulk replicates per method, arm and time; found ",
      toString(sort(unique(n_reps$replicate))),
      ".",
      call. = FALSE
    )
  }

  if (!"lower" %in% names(x)) {
    x$lower <- NA_real_
  }
  if (!"upper" %in% names(x)) {
    x$upper <- NA_real_
  }
  x
}

props <- read_proportions(in_file)

type_levels <- c(
  intersect(type_order, unique(props$cell_type)),
  sort(setdiff(unique(props$cell_type), type_order))
)
method_levels <- sort(unique(props$method))

props$arm <- factor(props$arm, levels = arm_levels)
props$time <- factor(props$time, levels = time_levels)
props$cell_type <- factor(props$cell_type, levels = type_levels)

# One facet per time x type, ordered by time, then by development.
facet_of <- function(time, cell_type) {
  key <- unique(data.frame(time = time, cell_type = cell_type))
  key <- key[order(key$time, key$cell_type), ]
  labels <- paste0(pretty_time(as.character(key$time)), " · ", key$cell_type)
  factor(
    paste0(pretty_time(as.character(time)), " · ", cell_type),
    levels = labels
  )
}

# ==========================================================================
# SECTION 2 · Centred log-ratio contrasts, early minus standard ----
# ==========================================================================

# CLR within each method x sample, over the types of that sample's
# reference.
sample_key <- paste(props$method, props$sample_id, sep = "\r")
log_p <- log(props$proportion + pseudo_count)
props$clr <- log_p - stats::ave(log_p, sample_key, FUN = mean)

resample <- function(v) v[sample.int(length(v), replace = TRUE)]

contrast_one <- function(clr_early, clr_standard) {
  boot <- replicate(
    n_boot,
    mean(resample(clr_early)) - mean(resample(clr_standard))
  )
  pairs <- outer(clr_early, clr_standard, "-")
  c(
    estimate = mean(clr_early) - mean(clr_standard),
    boot_lower = stats::quantile(boot, 0.025, names = FALSE),
    boot_upper = stats::quantile(boot, 0.975, names = FALSE),
    pair_min = min(pairs),
    pair_max = max(pairs),
    n_early = length(clr_early),
    n_standard = length(clr_standard)
  )
}

set.seed(seed)
groups <- split(
  props,
  list(props$method, props$time, props$cell_type),
  drop = TRUE
)
contrasts <- do.call(
  rbind,
  lapply(groups, function(g) {
    clr_early <- g$clr[g$arm == "early"]
    clr_standard <- g$clr[g$arm == "standard"]
    if (length(clr_early) == 0L || length(clr_standard) == 0L) {
      return(NULL)
    }
    data.frame(
      method = g$method[1L],
      time = g$time[1L],
      cell_type = g$cell_type[1L],
      t(contrast_one(clr_early, clr_standard)),
      row.names = NULL
    )
  })
)
contrasts$direction <- ifelse(
  contrasts$estimate > 0,
  "higher in early Wnt",
  "lower in early Wnt"
)
contrasts$facet <- facet_of(contrasts$time, contrasts$cell_type)
contrasts$method <- factor(contrasts$method, levels = rev(method_levels))

agreement <- do.call(
  rbind,
  lapply(
    split(contrasts, list(contrasts$time, contrasts$cell_type), drop = TRUE),
    function(g) {
      n_methods <- nrow(g)
      n_up <- sum(g$estimate > 0)
      consensus <- if (n_up == n_methods) {
        "higher in early Wnt, all methods"
      } else if (n_up == 0L) {
        "lower in early Wnt, all methods"
      } else {
        "methods disagree"
      }
      data.frame(
        time = g$time[1L],
        cell_type = g$cell_type[1L],
        n_methods = n_methods,
        n_up = n_up,
        n_pairs_same_sign = sum(g$pair_min > 0 | g$pair_max < 0),
        consensus = consensus,
        row.names = NULL
      )
    }
  )
)

utils::write.csv(
  contrasts[, setdiff(names(contrasts), "facet")],
  file.path(out_dir, "contrasts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
utils::write.csv(
  agreement,
  file.path(out_dir, "agreement.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ==========================================================================
# SECTION 3 · Panel b, estimated proportions per method ----
# ==========================================================================

# Methods on the y axis; the two arms are offset vertically by hand, so
# that the plot does not depend on vertical dodging.
arm_offset <- c(standard = -0.17, early = 0.17)
method_index <- stats::setNames(seq_along(method_levels), rev(method_levels))

props$facet <- facet_of(props$time, props$cell_type)
props$y <- method_index[props$method] + arm_offset[as.character(props$arm)]

level_means <- stats::aggregate(
  proportion ~ method + arm + time + cell_type,
  data = props,
  FUN = mean
)
level_means$facet <- facet_of(level_means$time, level_means$cell_type)
level_means$y <- method_index[level_means$method] +
  arm_offset[as.character(level_means$arm)]

panel_b <- ggplot2::ggplot(props, ggplot2::aes(y = y, colour = arm)) +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = lower, xmax = upper),
    data = props[!is.na(props$lower) & !is.na(props$upper), ],
    linewidth = 0.4,
    alpha = 0.6
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = proportion),
    size = 1.3,
    alpha = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = proportion),
    data = level_means,
    shape = 124,
    size = 4.5
  ) +
  ggplot2::scale_colour_manual(
    values = arm_colours,
    labels = arm_labels,
    name = NULL
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  ggplot2::scale_y_continuous(
    breaks = unname(method_index),
    labels = names(method_index),
    expand = ggplot2::expansion(add = 0.5)
  ) +
  ggplot2::facet_wrap(~facet, ncol = 4L) +
  ggplot2::labs(
    x = "Estimated proportion (dots: replicates; bar: mean)",
    y = NULL
  ) +
  theme_panel() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

# ==========================================================================
# SECTION 4 · Panel c, CLR contrast early minus standard ----
# ==========================================================================

panel_c <- ggplot2::ggplot(contrasts, ggplot2::aes(y = method)) +
  ggplot2::geom_vline(xintercept = 0, colour = ink, linewidth = 0.4) +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = pair_min, xmax = pair_max),
    colour = neutral,
    linewidth = 2.6
  ) +
  ggplot2::geom_linerange(
    ggplot2::aes(xmin = boot_lower, xmax = boot_upper),
    colour = ink,
    linewidth = 0.6
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = estimate, fill = direction),
    shape = 21,
    size = 2.4,
    colour = ink,
    stroke = 0.4
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "higher in early Wnt" = arm_colours[["early"]],
      "lower in early Wnt" = arm_colours[["standard"]]
    ),
    name = NULL
  ) +
  ggplot2::facet_wrap(~facet, ncol = 4L) +
  ggplot2::labs(
    x = paste0(
      "CLR contrast, early minus standard ",
      "(grey: nine replicate pairs; line: bootstrap 95 %)"
    ),
    y = NULL
  ) +
  theme_panel() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

# ==========================================================================
# SECTION 5 · Panel d, sign agreement across methods ----
# ==========================================================================

consensus_fill <- c(
  "higher in early Wnt, all methods" = arm_colours[["early"]],
  "lower in early Wnt, all methods" = arm_colours[["standard"]],
  "methods disagree" = neutral
)
agreement$label <- paste0(agreement$n_up, "/", agreement$n_methods)
agreement$text_colour <- ifelse(
  agreement$consensus == "methods disagree",
  ink,
  "white"
)
agreement$time_label <- factor(
  pretty_time(as.character(agreement$time)),
  levels = pretty_time(time_levels)
)

panel_d <- ggplot2::ggplot(
  agreement,
  ggplot2::aes(x = time_label, y = cell_type, fill = consensus)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 1) +
  ggplot2::geom_text(
    ggplot2::aes(label = label, colour = text_colour),
    size = 3.6
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_manual(values = consensus_fill, name = NULL) +
  ggplot2::scale_y_discrete(limits = rev(type_levels), drop = TRUE) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1L)) +
  ggplot2::labs(
    x = "Bulk time",
    y = NULL,
    subtitle = "Methods with a positive contrast / methods"
  ) +
  theme_panel() +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank())

# ==========================================================================
# SECTION 6 · Write panels and session information ----
# ==========================================================================

n_facets <- nlevels(droplevels(contrasts$facet))
facet_rows <- ceiling(n_facets / 4)
row_height <- 0.35 * length(method_levels) + 0.9

save_panel <- function(plot, filename, width, height) {
  ggplot2::ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

save_panel(panel_b, "panel_b_levels.pdf", 9.5, facet_rows * row_height + 1.2)
save_panel(
  panel_c,
  "panel_c_contrasts.pdf",
  9.5,
  facet_rows * row_height + 1.2
)
save_panel(panel_d, "panel_d_agreement.pdf", 4.2, 4.0)

writeLines(
  c(
    paste("seed:", seed),
    paste("n_boot:", n_boot),
    paste("pseudo_count:", pseudo_count),
    utils::capture.output(utils::sessionInfo())
  ),
  file.path(out_dir, "session_info.txt")
)
