#' Figures of the classifier benchmark (chapter 2, #sec-fs-classifiers and
#' #sec-fs-classifier-diagnostics).
#'
#' Plots only; `scripts/02_06_vis_classifier_cv.R` reads the tables of
#' `output/feature_selection/classifiers/` and saves them. Requires
#' `R/02_00_utils.R` (theme_slide()).

clf_method_labels <- c(
  two_stage = "Two-stage (two-sided)",
  glmnet = "glmnet, multinomial EN",
  msgl = "msgl, sparse-group lasso",
  plnlda = "PLN-LDA",
  ranger = "Random forest",
  ranger_balanced = "Balanced random forest",
  xgboost = "XGBoost",
  splsda = "sPLS-DA"
)
clf_method_colours <- c(
  two_stage = "#1c2b24",
  glmnet = "#1d4e89",
  msgl = "#5b8c5a",
  plnlda = "#8e6c8a",
  ranger = "#c45c26",
  ranger_balanced = "#e39a6c",
  xgboost = "#c9a227",
  splsda = "#4f9fa6"
)
clf_metric_labels <- c(
  macro_f1 = "Macro-F1",
  balanced_accuracy = "Balanced accuracy",
  accuracy = "Accuracy",
  log_loss = "Log loss (lower is better)"
)

#' Held-out metrics per method, one point per outer fold and the mean.
#'
#' @param metrics `metrics.csv` (one variant).
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_metrics <- function(metrics, title = NULL, subtitle = NULL) {
  long <- do.call(rbind, lapply(names(clf_metric_labels), function(k) {
    data.frame(
      time_point = metrics$time_point,
      method = metrics$method,
      fold = metrics$fold,
      metric = clf_metric_labels[[k]],
      value = metrics[[k]]
    )
  }))
  long <- long[is.finite(long$value), , drop = FALSE]
  long$metric <- factor(long$metric, levels = unname(clf_metric_labels))
  mf1 <- stats::aggregate(macro_f1 ~ method, metrics, mean)
  long$method <- factor(clf_method_labels[long$method], levels = clf_method_labels[mf1$method[order(mf1$macro_f1)]])
  summ <- stats::aggregate(value ~ time_point + method + metric, long, function(v) c(mean = mean(v), lo = min(v), hi = max(v)))
  summ <- cbind(summ[, 1:3], as.data.frame(summ$value))
  ggplot2::ggplot(long, ggplot2::aes(x = value, y = method)) +
    ggplot2::geom_linerange(data = summ, ggplot2::aes(xmin = lo, xmax = hi, x = NULL), colour = "#9a9286", linewidth = 0.6) +
    ggplot2::geom_point(colour = "#9a9286", size = 1, alpha = 0.7) +
    ggplot2::geom_point(data = summ, ggplot2::aes(x = mean), colour = "#1d4e89", size = 2.2) +
    # Free y drops a method with no fit at that time (msgl at 48 h).
    ggplot2::facet_grid(time_point ~ metric, scales = "free") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_slide() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 9))
}

#' Row-normalised confusion matrices.
#'
#' @param confusion `confusion.csv` (one variant).
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_confusion <- function(confusion, title = NULL, subtitle = NULL) {
  tot <- stats::ave(confusion$n, confusion$time_point, confusion$method, confusion$truth, FUN = sum)
  confusion$share <- ifelse(tot > 0, confusion$n / tot, NA_real_)
  confusion$method <- factor(clf_method_labels[confusion$method], levels = clf_method_labels)
  ggplot2::ggplot(confusion, ggplot2::aes(x = pred, y = truth, fill = share)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", share)), size = 2.2, colour = "#1c2b24") +
    ggplot2::scale_fill_gradient(low = "#f7f4ee", high = "#1d4e89", limits = c(0, 1), name = "Share of\ntrue type") +
    ggplot2::facet_grid(time_point ~ method, scales = "free") +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Predicted type", y = "True type") +
    theme_slide() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7),
      strip.text = ggplot2::element_text(size = 8)
    )
}

#' Per-type recall, one point per outer fold and the mean.
#'
#' @param recall `per_class_recall.csv` (one variant).
#' @param colours Named type colours.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_recall <- function(recall, colours, title = NULL, subtitle = NULL) {
  recall$method <- factor(clf_method_labels[recall$method], levels = rev(clf_method_labels))
  m <- stats::aggregate(recall ~ time_point + method + type, recall, mean)
  ggplot2::ggplot(recall, ggplot2::aes(x = recall, y = method, colour = type)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.4, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(data = m, size = 2.2, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::facet_wrap(~time_point, nrow = 1) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Recall on held-out cells", y = NULL) +
    theme_slide()
}

#' One-versus-rest reliability curves (calibration) of the probabilistic
#' methods.
#'
#' @param preds List of held-out prediction frames (one per time, method
#'   and fold) with `time_point`, `method`, `truth` and one probability
#'   column per type; frames without probabilities are skipped.
#' @param colours Named type colours.
#' @param n_bins Bins of predicted probability. Default `10`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_reliability <- function(preds, colours, n_bins = 10L, title = NULL, subtitle = NULL) {
  preds <- Filter(function(p) all(levels(p$truth) %in% names(p)), preds)
  key <- vapply(preds, function(p) paste(p$time_point[[1L]], p$method[[1L]]), "")
  d <- do.call(rbind, lapply(split(preds, key), function(group) {
    p <- do.call(rbind, group)
    do.call(rbind, lapply(levels(p$truth), function(k) {
      pr <- p[[k]]
      bin <- cut(pr, seq(0, 1, length.out = n_bins + 1L), include.lowest = TRUE)
      a <- stats::aggregate(cbind(pred = pr, obs = as.numeric(p$truth == k)) ~ bin, FUN = mean)
      a$n <- as.integer(table(bin)[as.character(a$bin)])
      cbind(time_point = p$time_point[[1L]], method = p$method[[1L]], type = k, a)
    }))
  }))
  d$method <- factor(clf_method_labels[d$method], levels = clf_method_labels)
  ggplot2::ggplot(d, ggplot2::aes(x = pred, y = obs, colour = type)) +
    ggplot2::geom_abline(colour = "#1c2b24", linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.7) +
    ggplot2::scale_size_area(max_size = 3, name = "Cells") +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::facet_grid(time_point ~ method) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Mean predicted probability (bin)", y = "Observed share of the type") +
    theme_slide() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8), axis.text = ggplot2::element_text(size = 7))
}

#' Inner cross-validation curve of the penalised logistic models.
#'
#' @param tuning Rows of `tuning_glmnet.csv` or `tuning_msgl.csv`.
#' @param ylab Axis label of the criterion.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_penalty_cv <- function(tuning, ylab, title = NULL, subtitle = NULL) {
  tuning$alpha <- factor(tuning$alpha)
  tuning$curve <- paste(tuning$fold, tuning$alpha)
  ch <- tuning[tuning$chosen %in% TRUE, , drop = FALSE]
  ggplot2::ggplot(tuning, ggplot2::aes(x = log(lambda), y = cvm, colour = alpha, group = curve)) +
    ggplot2::geom_line(linewidth = 0.35, alpha = 0.7) +
    ggplot2::geom_point(data = ch, colour = "#1c2b24", size = 1.6, shape = 21, fill = "white") +
    ggplot2::scale_colour_manual(values = grDevices::colorRampPalette(c("#c9c3b8", "#1d4e89"))(nlevels(tuning$alpha)), name = "alpha") +
    ggplot2::facet_wrap(~time_point, scales = "free") +
    ggplot2::labs(title = title, subtitle = subtitle, x = "log lambda", y = ylab) +
    theme_slide()
}

#' Coefficient paths (glmnet) of the 8 genes with the largest |beta|.
#'
#' @param path Data frame `time_point`, `type`, `gene`, `lambda`, `coef`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_paths <- function(path, title = NULL, subtitle = NULL) {
  path$panel <- paste0(path$time_point, ": ", path$type)
  ends <- path[path$lambda == stats::ave(path$lambda, path$panel, FUN = min), , drop = FALSE]
  ends <- ends[ends$coef != 0, , drop = FALSE]
  ggplot2::ggplot(path, ggplot2::aes(x = log(lambda), y = coef, group = gene)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#1c2b24", linewidth = 0.3) +
    ggplot2::geom_line(colour = "#1d4e89", linewidth = 0.4, alpha = 0.8) +
    ggrepel::geom_text_repel(data = ends, ggplot2::aes(label = gene), size = 2.2, fontface = "italic", direction = "y", hjust = 1, segment.size = 0.2, max.overlaps = Inf) +
    ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 4) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "log lambda (penalty decreases to the left)", y = "Coefficient") +
    ggplot2::scale_x_reverse() +
    theme_slide() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
}

#' Two-tier logit checks summarised per method, time and fold: events
#' per variable, largest VIF, and the share of near-separated cells.
#'
#' @param checks Data frame `time_point`, `method`, `fold`, `epv`,
#'   `max_vif`, `share_separated`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_logit_checks <- function(checks, title = NULL, subtitle = NULL) {
  lab <- c(epv = "Events per variable (rarest type / genes)", max_vif = "Largest VIF of the selected genes (log scale)", share_separated = "Held-out cells with p >= 0.999")
  ref <- data.frame(check = factor(lab[c("epv", "max_vif", "max_vif")], levels = lab), x = c(10, 5, 10))
  long <- do.call(rbind, lapply(names(lab), function(k) {
    data.frame(time_point = checks$time_point, method = clf_method_labels[checks$method], check = lab[[k]], value = checks[[k]])
  }))
  long$check <- factor(long$check, levels = lab)
  long <- long[is.finite(long$value), , drop = FALSE]
  ggplot2::ggplot(long, ggplot2::aes(x = value, y = time_point, colour = method)) +
    ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = x), colour = "#1c2b24", linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.8, alpha = 0.8, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_colour_manual(values = stats::setNames(clf_method_colours, clf_method_labels), name = NULL) +
    ggplot2::facet_wrap(~check, scales = "free_x") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_slide()
}

#' Linearity of the log-odds: binned empirical logits of one type against
#' log1p expression of its top genes, with a natural spline (3 df) fit.
#'
#' @param d Data frame `time_point`, `type`, `gene`, `x` (bin mean of
#'   log1p(x/s)), `logit` (empirical logit), `n`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_linearity <- function(d, title = NULL, subtitle = NULL) {
  d$panel <- paste0(d$time_point, ": ", d$type, ", ", d$gene)
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = logit)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ splines::ns(x, 3), colour = "#c45c26", fill = "#e8a77f", linewidth = 0.5) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#1d4e89", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(size = n), colour = "#1c2b24", alpha = 0.6) +
    ggplot2::scale_size_area(max_size = 2.5, name = "Cells") +
    ggplot2::facet_wrap(~panel, scales = "free", ncol = 4) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "log(1 + x/s), bin mean", y = "Empirical logit of the type") +
    theme_slide() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 7.5))
}

#' Held-out deviance contribution per cell, by true type.
#'
#' @param pred Held-out predictions with `dev` = -log p(true type).
#' @param colours Named type colours.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_influence <- function(pred, colours, title = NULL, subtitle = NULL) {
  pred$method <- factor(clf_method_labels[pred$method], levels = clf_method_labels)
  ggplot2::ggplot(pred, ggplot2::aes(x = pmax(dev, 1e-6), y = truth, colour = truth)) +
    ggplot2::geom_jitter(height = 0.2, size = 0.3, alpha = 0.3) +
    ggplot2::geom_vline(xintercept = log(2), colour = "#1c2b24", linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::facet_grid(time_point ~ method, scales = "free_y") +
    ggplot2::labs(title = title, subtitle = subtitle, x = "-log p(true type), held out (log scale)", y = NULL) +
    theme_slide() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
}

#' Variable importance of the top genes across outer folds (trees).
#'
#' @param imp `importance_top30.csv` rows of `ranger` and `xgboost`.
#' @param n Genes per method and time. Default `12`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_vip <- function(imp, n = 12L, title = NULL, subtitle = NULL) {
  imp$scaled <- imp$score / stats::ave(imp$score, imp$time_point, imp$method, imp$fold, FUN = max)
  m <- stats::aggregate(scaled ~ time_point + method + gene, imp, mean)
  m <- do.call(rbind, lapply(split(m, list(m$time_point, m$method), drop = TRUE), function(s) utils::head(s[order(-s$scaled), ], n)))
  d <- merge(imp, m[, c("time_point", "method", "gene")])
  d$panel <- paste0(d$time_point, ": ", clf_method_labels[d$method])
  m$panel <- paste0(m$time_point, ": ", clf_method_labels[m$method])
  key <- paste(m$panel, m$gene)
  m$key <- factor(key, levels = key[order(m$panel, m$scaled)])
  d$key <- factor(paste(d$panel, d$gene), levels = levels(m$key))
  ggplot2::ggplot(d, ggplot2::aes(x = scaled, y = key)) +
    ggplot2::geom_point(colour = "#9a9286", size = 0.9, alpha = 0.7) +
    ggplot2::geom_point(data = m, colour = "#c45c26", size = 2) +
    ggplot2::scale_y_discrete(labels = function(x) sub("^.* ", "", x)) +
    ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Importance / largest importance of the fold", y = NULL) +
    theme_slide() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7, face = "italic"), strip.text = ggplot2::element_text(size = 8))
}

#' Inner tuning of the tree models: macro-F1 per grid point.
#'
#' @param tuning Rows of `tuning_ranger.csv` and `tuning_xgboost.csv`
#'   (from `tune::collect_metrics()`), with a `setting` label column.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_tree_tuning <- function(tuning, title = NULL, subtitle = NULL) {
  tuning$method <- clf_method_labels[tuning$method]
  ggplot2::ggplot(tuning, ggplot2::aes(x = mean, y = setting)) +
    ggplot2::geom_point(colour = "#c45c26", size = 1.4, alpha = 0.7) +
    ggplot2::facet_grid(method ~ time_point, scales = "free_y", space = "free_y") +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Inner-CV macro-F1 (one point per outer fold)", y = NULL) +
    theme_slide() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
}

#' Latent coordinates of a discriminant method (sPLS-DA variates or
#' PLN-LDA scores) on training cells of one fold, with 95 % ellipses.
#'
#' @param d Data frame `time_point`, `type`, `axis1`, `axis2`.
#' @param colours Named type colours.
#' @param xlab,ylab Axis labels.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_latent <- function(d, colours, xlab, ylab, title = NULL, subtitle = NULL) {
  ggplot2::ggplot(d, ggplot2::aes(x = axis1, y = axis2, colour = type)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::stat_ellipse(level = 0.95, linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::facet_wrap(~time_point, scales = "free") +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_slide()
}

#' sPLS-DA inner tuning: balanced error rate against keepX per component.
#'
#' @param tuning `tuning_splsda.csv`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_splsda_tuning <- function(tuning, title = NULL, subtitle = NULL) {
  tuning$comp <- factor(paste("Component", tuning$comp))
  tuning$curve <- paste(tuning$fold, tuning$comp)
  ggplot2::ggplot(tuning, ggplot2::aes(x = keepX, y = ber, colour = comp, group = curve)) +
    ggplot2::geom_line(linewidth = 0.4, alpha = 0.7) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c("#1d4e89", "#5b8c5a", "#c45c26", "#9a9286"), name = NULL) +
    ggplot2::facet_wrap(~time_point) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "keepX (genes per component, log scale)", y = "Inner-CV balanced error rate") +
    theme_slide()
}

#' Mean-variance relation and the highly variable genes (Seurat v3 vst).
#'
#' @param hvg Rows of `hvg_<time>.csv` with `time_point`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_hvg <- function(hvg, title = NULL, subtitle = NULL) {
  hvg <- hvg[hvg$variance > 0, , drop = FALSE]
  hvg$class <- ifelse(hvg$hvg, "Top 2,000 (vst)", "Other genes")
  hvg <- hvg[order(hvg$hvg), , drop = FALSE]
  ggplot2::ggplot(hvg, ggplot2::aes(x = mean, y = variance, colour = class)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = variance_expected), data = hvg[order(hvg$mean), ], colour = "#1c2b24", linewidth = 0.5) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = c("Top 2,000 (vst)" = "#c45c26", "Other genes" = "#c9c3b8"), name = NULL) +
    ggplot2::facet_wrap(~time_point) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Mean count (log scale)", y = "Variance (log scale)") +
    theme_slide()
}

#' Jaccard similarity between the gene x type sets of the methods.
#'
#' @param sets Data frame `time_point`, `method`, `type`, `gene`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_jaccard <- function(sets, title = NULL, subtitle = NULL) {
  d <- do.call(rbind, lapply(split(sets, sets$time_point), function(s) {
    key <- split(paste(s$type, s$gene), s$method)
    ms <- intersect(names(clf_method_labels), names(key))
    g <- expand.grid(a = ms, b = ms, stringsAsFactors = FALSE)
    g$jaccard <- mapply(function(a, b) {
      u <- length(union(key[[a]], key[[b]]))
      if (u == 0L) NA_real_ else length(intersect(key[[a]], key[[b]])) / u
    }, g$a, g$b)
    cbind(time_point = s$time_point[[1L]], g)
  }))
  d$a <- factor(clf_method_labels[d$a], levels = clf_method_labels)
  d$b <- factor(clf_method_labels[d$b], levels = rev(clf_method_labels))
  ggplot2::ggplot(d, ggplot2::aes(x = a, y = b, fill = jaccard)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", jaccard)), size = 2.4) +
    ggplot2::scale_fill_gradient(low = "#f7f4ee", high = "#1d4e89", limits = c(0, 1), name = "Jaccard") +
    ggplot2::facet_wrap(~time_point) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_slide() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7), axis.text.y = ggplot2::element_text(size = 7))
}

#' Primary against the 48 h sensitivities (cell-line covariate, ComBat-seq).
#'
#' @param metrics `metrics_sensitivity_48h.csv`.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot.
plot_clf_sensitivity <- function(metrics, title = NULL, subtitle = NULL) {
  long <- do.call(rbind, lapply(c("macro_f1", "balanced_accuracy"), function(k) {
    data.frame(method = clf_method_labels[metrics$method], variant = metrics$variant, metric = clf_metric_labels[[k]], value = metrics[[k]])
  }))
  m <- stats::aggregate(value ~ method + variant + metric, long, mean)
  ggplot2::ggplot(long, ggplot2::aes(x = value, y = method, colour = variant)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.4, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(data = m, size = 2.2, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_colour_manual(values = c(primary = "#1d4e89", covariate = "#c45c26", combat = "#5b8c5a"), name = NULL) +
    ggplot2::facet_wrap(~metric, scales = "free_x") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_slide()
}
