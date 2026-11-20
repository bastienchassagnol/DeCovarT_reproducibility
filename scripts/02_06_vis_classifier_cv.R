# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_06_vis_classifier_cv.R [dir=classifiers] \
#   > "logs/02_06_vis_classifier_cv_$(date +%F).log" 2>&1 &
#
# Figures of the classifier benchmark (chapter 2, #sec-fs-classifiers and
# the diagnostics appendix): held-out metrics, confusion, recall,
# reliability, penalised-logit checks (CV curves, paths, events per
# variable, VIF, separation, linearity of the log-odds, influence), tree
# importance and tuning, sPLS-DA and PLN-LDA latent maps, HVGs, Jaccard,
# 48 h sensitivities and one UpSet per type x time against the two-stage
# two-sided panel. Reads output/feature_selection/<dir>/ (written by
# 02_06_stats_classifier_cv.R) and
# output/feature_selection/deg_twosided/panel_twosided_<time>.csv.
# PDFs and data in output/feature_selection/<dir>/figures/, PNGs in
# docs/images/feature_selection/ (prefix fs_clf_; with dir other than
# `classifiers`, PNGs stay in the figures folder).

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_06_stats_classifier_cv.R"))
source(file.path(root, "R", "02_06_vis_classifier_cv.R"))

cli <- commandArgs(trailingOnly = TRUE)
dir_name <- sub("^dir=", "", c(grep("^dir=", cli, value = TRUE), "dir=classifiers")[[1L]])
clf_dir <- file.path(root, "output", "feature_selection", dir_name)
fig_dir <- file.path(clf_dir, "figures")
png_dir <- if (dir_name == "classifiers") file.path(root, "docs", "images", "feature_selection") else fig_dir
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
read <- function(name) utils::read.csv(file.path(clf_dir, name), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
type_cols <- celltype_colours(file.path(root, "data", "dictionaries", "suppinger_celltype_colours.csv"))
lfc_column <- "abs_log2fc_ge_1.000"
save <- function(plot, name, data, width, height) {
  save_figure(plot, name, fig_dir, png_dir = png_dir, data = data, width = width, height = height)
}

metrics_all <- read("metrics.csv")
metrics <- metrics_all[metrics_all$variant == "primary", , drop = FALSE]
recall <- read("per_class_recall.csv")
recall <- recall[recall$variant == "primary", , drop = FALSE]
confusion <- read("confusion.csv")
confusion <- confusion[confusion$variant == "primary", , drop = FALSE]
freq <- read("selection_frequency.csv")
stable <- freq[freq$stable & freq$variant == "primary", , drop = FALSE]
imp <- read("importance_top30.csv")
imp <- imp[imp$variant == "primary", , drop = FALSE]
pr <- readRDS(file.path(clf_dir, "predictions.rds"))
pred <- Filter(function(p) !is.null(p) && p$variant[[1L]] == "primary", pr$predictions)
diag <- Filter(function(d) d$variant == "primary", pr$diag)
times <- sort(unique(metrics$time_point))
methods <- unique(metrics$method)
has <- function(m) m %in% methods

# 1 · Held-out metrics, confusion, recall -----------------------------------------
mf1 <- stats::aggregate(macro_f1 ~ method, metrics, mean)
best <- mf1$method[which.max(mf1$macro_f1)]
save(
  plot_clf_metrics(
    metrics,
    title = sprintf("%s has the highest mean macro-F1 (%.2f over the times shown)", clf_method_labels[[best]], max(mf1$macro_f1)),
    subtitle = "Outer 5-fold CV stratified by type; grey: folds and their range; blue: mean. Cells of one specimen"
  ),
  "fs_clf_metrics", metrics, 12, 7
)
save(
  plot_clf_confusion(confusion, title = "Held-out confusion: where each type is misread", subtitle = "Cells pooled over the outer folds; each row sums to 1"),
  "fs_clf_confusion", confusion, 14, 8
)
save(
  plot_clf_recall(recall, type_cols, title = "Recall of each type on held-out cells", subtitle = "Small points: outer folds; large points: mean"),
  "fs_clf_recall", recall, 13, 5.5
)

# 2 · Reliability and per-cell deviance -------------------------------------------
prob_pred <- Filter(function(p) all(levels(p$truth) %in% names(p)), pred)
save(
  plot_clf_reliability(
    prob_pred, type_cols,
    title = "Calibration of the probabilistic classifiers, one type against the rest",
    subtitle = "Held-out cells binned by predicted probability (10 bins); dashed: perfect calibration. sPLS-DA gives no probability"
  ),
  "fs_clf_reliability", lapply(prob_pred, function(p) p[, setdiff(names(p), "cell")]), 14, 8
)
dev <- do.call(rbind, lapply(prob_pred, function(p) {
  types <- levels(p$truth)
  pm <- as.matrix(p[, types])
  data.frame(
    time_point = p$time_point,
    method = p$method,
    fold = p$fold,
    truth = p$truth,
    dev = -log(pmax(pm[cbind(seq_len(nrow(p)), as.integer(p$truth))], 1e-15)),
    pmax = apply(pm, 1L, max)
  )
}))

# 3 · Penalised logistic models (glmnet, msgl) ------------------------------------
logit_methods <- intersect(c("glmnet", "msgl"), methods)
if (length(logit_methods) > 0L) {
  for (m in logit_methods) {
    tu <- read(sprintf("tuning_%s.csv", m))
    tu <- tu[tu$variant == "primary", , drop = FALSE]
    save(
      plot_clf_penalty_cv(
        tu,
        ylab = if (m == "glmnet") "Inner-CV multinomial deviance" else "Inner-CV negative log-likelihood",
        title = sprintf("%s: inner cross-validation over lambda and alpha", clf_method_labels[[m]]),
        subtitle = "One curve per outer fold and alpha; open circle: chosen (alpha, lambda)"
      ),
      sprintf("fs_clf_%s_cv", m), tu, 12, 4.5
    )
  }
  if (has("glmnet")) {
    path <- do.call(rbind, lapply(Filter(function(d) d$method == "glmnet" && d$fold == 1L, diag), function(d) {
      cbind(time_point = d$time_point, d$path)
    }))
    save(
      plot_clf_paths(path, title = "glmnet: coefficient paths of the 8 genes with the largest |beta| (outer fold 1)", subtitle = "Ungrouped multinomial elastic net at the chosen alpha; one panel per time and type"),
      "fs_clf_glmnet_paths", path, 13, 8
    )
  }
  # Checks per method, time and fold: EPV, separation; VIF on the stable
  # set of each time over all cells of the slot.
  slots <- readRDS(file.path(root, "output", "feature_selection", "tier_01", "scored_slots.rds"))
  cc <- readRDS(file.path(root, "output", "feature_selection", "tier_02", "cross_cell.rds"))
  tot <- stats::aggregate(n_test ~ time_point + method + variant + type, recall, sum)
  checks <- do.call(rbind, lapply(split(metrics[metrics$method %in% logit_methods, ], seq_len(sum(metrics$method %in% logit_methods))), function(r) {
    tr <- recall[recall$time_point == r$time_point & recall$method == r$method & recall$fold == r$fold, ]
    tt <- tot[tot$time_point == r$time_point & tot$method == r$method, ]
    n_train <- tt$n_test[match(tr$type, tt$type)] - tr$n_test
    dd <- dev[dev$time_point == r$time_point & dev$method == r$method & dev$fold == r$fold, ]
    data.frame(
      time_point = r$time_point, method = r$method, fold = r$fold,
      epv = min(n_train) / max(1, r$n_selected),
      share_separated = mean(dd$pmax >= 0.999)
    )
  }))
  vif <- do.call(rbind, lapply(logit_methods, function(m) {
    do.call(rbind, lapply(times, function(tm) {
      g <- unique(stable$gene[stable$method == m & stable$time_point == tm])
      if (length(g) < 2L) return(data.frame(time_point = tm, method = m, n_genes = length(g), max_vif = NA_real_, condition_number = NA_real_))
      sl <- slots$slots[[tm]]
      x <- log1p(t(as.matrix(sl$counts[g, , drop = FALSE])) / cc$sf[[tm]])
      r <- stats::cor(x)
      ev <- eigen(r, symmetric = TRUE, only.values = TRUE)$values
      inv <- tryCatch(solve(r), error = function(e) NULL)
      data.frame(
        time_point = tm, method = m, n_genes = length(g),
        max_vif = if (is.null(inv) || min(ev) <= 1e-10) Inf else max(diag(inv)),
        condition_number = sqrt(max(ev) / max(min(ev), .Machine$double.eps))
      )
    }))
  }))
  checks <- merge(checks, vif[, c("time_point", "method", "max_vif")], by = c("time_point", "method"))
  utils::write.csv(merge(checks, vif), file.path(clf_dir, "logit_checks.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  save(
    plot_clf_logit_checks(
      checks,
      title = "Penalised logit checks: few events per selected gene, collinear selected genes",
      subtitle = "Dashed: EPV 10; VIF 5 and 10. VIF of the stable genes over all cells of the slot (Inf: singular correlation)"
    ),
    "fs_clf_logit_checks", list(checks = checks, vif = vif), 12, 4.5
  )
  save(
    plot_clf_influence(
      dev[dev$method %in% logit_methods, ], type_cols,
      title = "Held-out deviance per cell: the poorly predicted cells",
      subtitle = "Each point is a cell of an outer test fold; dashed: p(true type) = 1/2"
    ),
    "fs_clf_logit_influence", dev[dev$method %in% logit_methods, ], 12, 7
  )
  # Linearity of the log-odds: two glmnet-stable genes per type and time.
  lin_method <- logit_methods[[1L]]
  lin <- do.call(rbind, lapply(times, function(tm) {
    sl <- slots$slots[[tm]]
    ct <- as.character(sl$cell_type)
    st <- stable[stable$method == lin_method & stable$time_point == tm & stable$share_up >= 0.5, ]
    mi <- stats::aggregate(score ~ gene, imp[imp$method == lin_method & imp$time_point == tm, ], mean)
    st$score <- mi$score[match(st$gene, mi$gene)]
    st$score[is.na(st$score)] <- 0
    do.call(rbind, lapply(split(st, st$type), function(s) {
      g <- utils::head(s$gene[order(-s$score, -s$frequency)], 2L)
      do.call(rbind, lapply(g, function(gg) {
        x <- log1p(as.numeric(sl$counts[gg, ]) / cc$sf[[tm]])
        bin <- ifelse(x == 0, 0L, 1L + dplyr::ntile(replace(x, x == 0, NA), 9L))
        y <- ct == s$type[[1L]]
        a <- stats::aggregate(cbind(x = x, y = y) ~ bin, FUN = mean)
        n <- as.numeric(table(bin)[as.character(a$bin)])
        k <- a$y * n
        data.frame(time_point = tm, type = s$type[[1L]], gene = gg, x = a$x, logit = log((k + 0.5) / (n - k + 0.5)), n = n)
      }))
    }))
  }))
  save(
    plot_clf_linearity(
      lin,
      title = sprintf("Linearity of the log-odds for two %s-stable genes per type", clf_method_labels[[lin_method]]),
      subtitle = "Empirical logit of the type in the zero bin and 9 quantile bins of the non-zero cells; orange: ns(x, 3) fit; dashed blue: straight line"
    ),
    "fs_clf_logit_linearity", lin, 13, 9
  )
}

# 4 · Trees ------------------------------------------------------------------------
tree_methods <- intersect(c("ranger", "xgboost"), methods)
if (length(tree_methods) > 0L) {
  vi <- imp[imp$method %in% tree_methods, , drop = FALSE]
  save(
    plot_clf_vip(vi, title = "Tree importance of the top genes, across outer folds", subtitle = "Random forest: permutation importance; XGBoost: gain; scaled by the largest of each fold; orange: mean"),
    "fs_clf_tree_vip", vi, 13, 10
  )
  tu <- do.call(rbind, lapply(tree_methods, function(m) {
    t <- read(sprintf("tuning_%s.csv", m))
    t <- t[t$variant == "primary", , drop = FALSE]
    t$setting <- if (m == "ranger") sprintf("mtry %d, min_n %d", t$mtry, t$min_n) else sprintf("depth %d, rate %.2f, mtry %d", t$tree_depth, t$learn_rate, t$mtry)
    t[, c("time_point", "method", "fold", "setting", "mean")]
  }))
  save(
    plot_clf_tree_tuning(tu, title = "Tree tuning: macro-F1 on the inner folds per grid point", subtitle = "Random forest: 500 trees; XGBoost: 300 rounds"),
    "fs_clf_tree_tuning", tu, 12, 5.5
  )
}

# 5 · sPLS-DA and PLN-LDA ------------------------------------------------------------
latent <- function(m, slot_name) {
  do.call(rbind, lapply(Filter(function(d) d$method == m && d$fold == 1L, diag), function(d) {
    v <- d[[slot_name]]
    num <- names(v)[vapply(v, is.numeric, logical(1))]
    data.frame(time_point = d$time_point, type = as.character(v$type), axis1 = v[[num[[1L]]]], axis2 = if (length(num) > 1L) v[[num[[2L]]]] else 0)
  }))
}
if (has("splsda")) {
  tu <- read("tuning_splsda.csv")
  tu <- tu[tu$variant == "primary", , drop = FALSE]
  save(
    plot_clf_splsda_tuning(tu, title = "sPLS-DA: inner balanced error rate against keepX", subtitle = "mixOmics::tune.splsda(), max distance; one curve per outer fold and component"),
    "fs_clf_splsda_tuning", tu, 12, 4.5
  )
  lv <- latent("splsda", "variates")
  save(
    plot_clf_latent(lv, type_cols, "Component 1", "Component 2", title = "sPLS-DA variates of the training cells (outer fold 1)", subtitle = "Up to 1,500 cells; 95 % normal ellipses"),
    "fs_clf_splsda_variates", lv, 12, 4.5
  )
}
if (has("plnlda")) {
  lv <- latent("plnlda", "scores")
  save(
    plot_clf_latent(lv, type_cols, "LDA axis 1", "LDA axis 2", title = "PLN-LDA scores of the training cells (outer fold 1)", subtitle = "Raw counts, offset log library size, diagonal covariance; up to 1,500 cells; 95 % ellipses"),
    "fs_clf_plnlda_scores", lv, 12, 4.5
  )
}

# 6 · HVGs --------------------------------------------------------------------------
hv <- do.call(rbind, lapply(times, function(tm) cbind(time_point = tm, read(sprintf("hvg_%s.csv", tm)))))
save(
  plot_clf_hvg(hv, title = "Global top 2,000 highly variable genes per slot (Seurat v3 vst)", subtitle = "Universe of msgl and PLN-LDA; black: loess of log variance on log mean (span 0.3)"),
  "fs_clf_hvg", hv, 12, 4.5
)

# 7 · Comparison with the two-stage panel -------------------------------------------
panel <- do.call(rbind, lapply(times, function(tm) {
  p <- utils::read.csv(file.path(root, "output", "feature_selection", "deg_twosided", sprintf("panel_twosided_%s.csv", tm)), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  p <- p[p[[lfc_column]] %in% TRUE, , drop = FALSE]
  data.frame(time_point = tm, method = "two_stage", type = p$type, gene = p$gene)
}))
sets <- rbind(panel, stable[, c("time_point", "method", "type", "gene")])
utils::write.csv(sets, file.path(clf_dir, "comparison_sets.csv"), row.names = FALSE, fileEncoding = "UTF-8")
save(
  plot_clf_jaccard(sets, title = "Overlap of the gene x type sets between methods", subtitle = "Two-stage: A2 confirmed and |log2FC| >= 1; classifiers: stable sets (>= 60 % of outer folds)"),
  "fs_clf_jaccard", sets, 13, 5.5
)
upset_files <- list()
for (tm in times) {
  for (tp in sort(unique(sets$type[sets$time_point == tm]))) {
    s <- sets[sets$time_point == tm & sets$type == tp, , drop = FALSE]
    groups <- split(s$gene, clf_method_labels[s$method])
    groups <- groups[lengths(groups) > 0L]
    if (length(groups) < 2L) next
    genes <- sort(unique(unlist(groups)))
    up <- as.data.frame(lapply(groups, function(g) genes %in% g), check.names = FALSE)
    up$gene <- genes
    p_up <- ComplexUpset::upset(
      up,
      intersect = names(groups),
      n_intersections = 20,
      width_ratio = 0.25,
      base_annotations = list("Genes" = ComplexUpset::intersection_size(fill = "#1d4e89", text = list(size = 2.5))),
      set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(fill = "#9a9286"))
    ) +
      patchwork::plot_annotation(
        title = sprintf("%s, %s: genes of the two-stage panel and of the classifiers", tm, tp),
        subtitle = "Two-stage: A2 confirmed, |log2FC| >= 1 (both directions); classifiers: stable sets; 20 largest intersections",
        theme = theme_slide()
      )
    nm <- sprintf("fs_clf_upset_%s_%s", tm, gsub("[^a-z0-9]+", "_", tolower(iconv(tp, to = "ASCII//TRANSLIT"))))
    save(p_up, nm, up, 11, 5.5)
    upset_files[[length(upset_files) + 1L]] <- data.frame(time_point = tm, type = tp, file = nm)
  }
}
utils::write.csv(do.call(rbind, upset_files), file.path(clf_dir, "upset_files.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 8 · 48 h sensitivities --------------------------------------------------------------
sens <- metrics_all[metrics_all$time_point == "48h", , drop = FALSE]
if (length(unique(sens$variant)) > 1L) {
  save(
    plot_clf_sensitivity(sens, title = "48 h: cell-line covariate and ComBat-seq against the uncorrected fit", subtitle = "Small points: outer folds; large points: mean"),
    "fs_clf_sensitivity_48h", sens, 11, 4.5
  )
}
message("Wrote ", fig_dir)
