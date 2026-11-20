#' Classifier-based marker selection: nested cross-validation per slot.
#'
#' Six classifiers of the cell type, compared with the two-stage panel
#' (omnibus gate, then two-sided post hoc contrasts; chapter 2,
#' `#sec-fs-classifiers`). Theory and the choices not taken: DeCovarT
#' vignette 3.1, `#sec-fs-classifier-theory`. Each time point is analysed
#' on its own. The resampling unit is the cell (one specimen, no donor
#' column): the scores measure generalisation to new cells of this
#' specimen.
#'
#' - `glmnet`: multinomial elastic net, `type.multinomial = "ungrouped"`;
#' - `msgl`: multinomial sparse-group lasso, groups = co-expression modules;
#' - `plnlda`: PLN-LDA on raw counts, offset log library size;
#' - `ranger`: random forest (`parsnip`), permutation importance;
#' - `xgboost`: boosted trees (`parsnip`), gain importance;
#' - `splsda`: sparse PLS discriminant analysis (`mixOmics`).
#'
#' Inputs other than PLN-LDA are log1p(x / s) (s = library size / median);
#' `glmnet`, `msgl` and `splsda` standardise the genes internally on the
#' training cells, and trees are invariant to monotone scaling.

classifier_methods <- c("glmnet", "msgl", "plnlda", "ranger", "xgboost", "splsda")

#' Highly variable genes, Seurat v3 `vst` (reimplemented, no Seurat).
#'
#' Per gene, mean and variance of the raw counts; a loess of log10
#' variance on log10 mean gives the expected variance; each count is
#' standardised with the expected standard deviation and clipped at
#' sqrt(C); genes are ranked by the variance of the standardised values.
#'
#' @param counts Genes x cells raw counts (dense or sparse).
#' @param n_top Genes kept. Default `2000`.
#' @param span Loess span. Default `0.3`.
#' @return Data frame: `gene`, `mean`, `variance`, `variance_expected`,
#'   `variance_standardised`, `rank`, `hvg` (rank <= n_top).
select_hvg_vst <- function(counts, n_top = 2000L, span = 0.3) {
  n_cells <- ncol(counts)
  mu <- Matrix::rowMeans(counts)
  v <- (Matrix::rowSums(counts^2) - n_cells * mu^2) / (n_cells - 1)
  ok <- v > 0
  fit <- stats::loess(log10(v[ok]) ~ log10(mu[ok]), span = span, degree = 2)
  v_exp <- rep(NA_real_, length(mu))
  v_exp[ok] <- 10^stats::fitted(fit)
  sd_exp <- sqrt(v_exp)
  clip <- sqrt(n_cells)
  # Standardised variance from sums of clipped values, gene by gene block.
  v_std <- rep(0, length(mu))
  for (s in seq(1L, nrow(counts), by = 500L)) {
    rows <- s:min(s + 499L, nrow(counts))
    rows <- rows[ok[rows]]
    if (length(rows) == 0L) next
    z <- (as.matrix(counts[rows, , drop = FALSE]) - mu[rows]) / sd_exp[rows]
    z <- pmin(z, clip)
    v_std[rows] <- rowSums((z - rowMeans(z))^2) / (n_cells - 1)
  }
  rk <- rank(-v_std, ties.method = "first")
  data.frame(
    gene = rownames(counts),
    mean = mu,
    variance = v,
    variance_expected = v_exp,
    variance_standardised = v_std,
    rank = rk,
    hvg = rk <= n_top,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Per-cell weights proportional to n_c^-gamma, mean one.
#'
#' @param y Factor of types.
#' @param gamma Exponent. Default `0.5`.
#' @return Numeric vector, one weight per cell.
class_weights <- function(y, gamma = 0.5) {
  n_c <- table(y)
  w <- as.numeric(n_c[as.character(y)])^(-gamma)
  w / mean(w)
}

#' Folds stratified by type.
#'
#' @param y Factor of types.
#' @param K Folds.
#' @param seed Seed; type j uses `seed + j`.
#' @return Integer fold per cell.
stratified_fold_ids <- function(y, K, seed) {
  f <- integer(length(y))
  for (j in seq_along(levels(y))) {
    idx <- which(y == levels(y)[[j]])
    f[idx] <- make_cv_folds(length(idx), K, seed + j)
  }
  f
}

#' Mean difference of log1p(x/s) of each type against the average of the
#' other types (genes x types), used to give each gene a type.
#'
#' @param X Cells x genes log1p(x/s) (sparse or dense).
#' @param y Factor of types.
#' @return Genes x types matrix.
type_contrast_lognorm <- function(X, y) {
  z <- Matrix::sparseMatrix(
    i = seq_along(y),
    j = as.integer(y),
    x = 1,
    dims = c(length(y), nlevels(y))
  )
  m <- as.matrix(Matrix::crossprod(X, z))
  m <- sweep(m, 2L, as.numeric(table(y)), "/")
  k <- ncol(m)
  a <- matrix(-1 / (k - 1), k, k)
  diag(a) <- 1
  out <- m %*% t(a)
  dimnames(out) <- list(colnames(X), levels(y))
  out
}

#' Give each gene the type with the largest |contrast|.
#'
#' @param genes Genes.
#' @param contrast Output of [type_contrast_lognorm()].
#' @return Data frame `gene`, `type`, `direction`.
assign_type <- function(genes, contrast) {
  cm <- contrast[genes, , drop = FALSE]
  j <- max.col(abs(cm), ties.method = "first")
  data.frame(
    gene = genes,
    type = colnames(cm)[j],
    direction = ifelse(cm[cbind(seq_along(genes), j)] > 0, "up", "down"),
    stringsAsFactors = FALSE
  )
}

#' Co-expression modules: Ward clustering on 1 - |cor|.
#'
#' Same construction as `scripts/04_03_stats_selection_checks.R`.
#'
#' @param X Cells x genes.
#' @param n_modules Modules.
#' @return Integer module per gene.
coexpression_modules <- function(X, n_modules) {
  r <- stats::cor(as.matrix(X))
  r[!is.finite(r)] <- 0
  stats::cutree(stats::hclust(stats::as.dist(1 - abs(r)), "ward.D2"), n_modules)
}

#' Held-out classification metrics.
#'
#' @param truth,pred Factors with the same levels.
#' @param prob Cells x types probabilities, or `NULL`.
#' @return One-row data frame: `macro_f1`, `balanced_accuracy`,
#'   `accuracy`, `log_loss`, `brier`.
classification_metrics <- function(truth, pred, prob = NULL) {
  lv <- levels(truth)
  rec <- vapply(lv, function(k) mean(pred[truth == k] == k), numeric(1))
  prec <- vapply(lv, function(k) {
    if (any(pred == k)) mean(truth[pred == k] == k) else 0
  }, numeric(1))
  f1 <- ifelse(rec + prec > 0, 2 * rec * prec / (rec + prec), 0)
  ll <- br <- NA_real_
  if (!is.null(prob)) {
    pt <- prob[cbind(seq_along(truth), match(as.character(truth), colnames(prob)))]
    ll <- -mean(log(pmax(pt, 1e-15)))
    onehot <- outer(as.character(truth), colnames(prob), "==") * 1
    br <- mean(rowSums((prob - onehot)^2))
  }
  data.frame(
    macro_f1 = mean(f1),
    balanced_accuracy = mean(rec),
    accuracy = mean(pred == truth),
    log_loss = ll,
    brier = br
  )
}

#' Inner-CV tuning and outer fit of one method on one outer split.
#'
#' @param method One of [classifier_methods].
#' @param X Cells x genes log1p(x/s) of all cells of the slot (sparse),
#'   restricted to the universe of the method.
#' @param Y Cells x genes raw counts (dense, same genes as `X`), for
#'   `plnlda` only.
#' @param lib Library size per cell, for `plnlda` only.
#' @param y Factor of types.
#' @param train,test Row indices of the outer split.
#' @param inner Inner fold per training cell (1..K_inner).
#' @param covariate Optional numeric column (cell line), unpenalised where
#'   the method allows; never selected.
#' @param settings List of grids and sizes (see the script).
#' @return List: `prob` (test cells x types, or `NULL`), `pred`,
#'   `selected` (gene, type, direction), `importance` (gene, score),
#'   `tuning` (data frame), `diag` (method-specific list).
fit_classifier <- function(method, X, Y, lib, y, train, test, inner, covariate = NULL, settings) {
  types <- levels(y)
  y_tr <- droplevels(y[train])
  stopifnot(identical(levels(y_tr), types))
  w <- class_weights(y_tr, settings$gamma)
  contrast <- type_contrast_lognorm(X[train, , drop = FALSE], y_tr)
  genes <- colnames(X)
  Xtr <- X[train, , drop = FALSE]
  Xte <- X[test, , drop = FALSE]
  if (!is.null(covariate)) {
    Xtr <- cbind(Xtr, cell_line = covariate[train])
    Xte <- cbind(Xte, cell_line = covariate[test])
  }
  top_assign <- function(score) {
    score <- score[names(score) %in% genes]
    g <- names(utils::head(sort(score, decreasing = TRUE), settings$n_top))
    assign_type(g, contrast)
  }
  switch(
    method,
    glmnet = {
      pf <- c(rep(1, length(genes)), if (!is.null(covariate)) 0)
      cvs <- lapply(settings$glmnet_alpha, function(a) {
        glmnet::cv.glmnet(
          Xtr, y_tr,
          family = "multinomial",
          type.multinomial = "ungrouped",
          alpha = a,
          weights = w,
          foldid = inner,
          penalty.factor = pf,
          nlambda = settings$n_lambda,
          type.measure = "deviance"
        )
      })
      best <- which.min(vapply(cvs, function(f) min(f$cvm), numeric(1)))
      cv <- cvs[[best]]
      fit <- cv$glmnet.fit
      il <- which(fit$lambda == cv$lambda.min)
      prob <- stats::predict(fit, Xte, s = cv$lambda.min, type = "response")[, , 1L]
      beta <- vapply(types, function(k) as.numeric(fit$beta[[k]][, il]), numeric(ncol(Xtr)))
      rownames(beta) <- colnames(Xtr)
      beta <- beta[genes, , drop = FALSE]
      nz <- which(beta != 0, arr.ind = TRUE)
      selected <- data.frame(
        gene = genes[nz[, 1L]],
        type = types[nz[, 2L]],
        direction = ifelse(beta[nz] > 0, "up", "down"),
        stringsAsFactors = FALSE
      )
      imp <- apply(abs(beta), 1L, max)
      # Coefficient paths of the 8 genes with the largest |beta| at lambda.min.
      top8 <- names(utils::head(sort(imp[imp > 0], decreasing = TRUE), 8L))
      path <- do.call(rbind, lapply(types, function(k) {
        b <- fit$beta[[k]][top8, , drop = FALSE]
        data.frame(
          type = k,
          gene = rep(top8, times = ncol(b)),
          lambda = rep(fit$lambda, each = length(top8)),
          coef = as.vector(as.matrix(b))
        )
      }))
      list(
        prob = prob,
        selected = selected,
        importance = data.frame(gene = genes, score = imp),
        tuning = do.call(rbind, Map(function(f, a) {
          data.frame(alpha = a, lambda = f$lambda, cvm = f$cvm, cvsd = f$cvsd, nzero = f$nzero, chosen = a == settings$glmnet_alpha[[best]] & f$lambda == cv$lambda.min)
        }, cvs, settings$glmnet_alpha)),
        diag = list(path = path, alpha = settings$glmnet_alpha[[best]], lambda = cv$lambda.min)
      )
    },
    msgl = {
      groups <- coexpression_modules(X[train, , drop = FALSE], settings$n_modules)
      k <- length(types)
      gw <- NULL
      pw <- NULL
      if (!is.null(covariate)) {
        groups <- c(groups, settings$n_modules + 1L)
        sizes <- as.numeric(table(factor(groups, levels = seq_len(settings$n_modules + 1L))))
        gw <- sqrt(k * sizes)
        gw[length(gw)] <- 0
        pw <- matrix(1, k, length(groups))
        pw[, length(groups)] <- 0
      }
      cv_idx <- lapply(sort(unique(inner)), function(f) which(inner == f))
      cvs <- lapply(settings$msgl_alpha, function(a) {
        lam <- msgl::lambda(
          Xtr, y_tr,
          sampleWeights = w,
          grouping = groups,
          groupWeights = gw,
          parameterWeights = pw,
          alpha = a,
          d = settings$n_lambda,
          lambda.min = settings$lambda_min_ratio,
          lambda.min.rel = TRUE
        )
        f <- msgl::cv(
          Xtr, y_tr,
          sampleWeights = w,
          grouping = groups,
          groupWeights = gw,
          parameterWeights = pw,
          alpha = a,
          lambda = lam,
          cv.indices = cv_idx
        )
        list(lambda = lam, err = as.numeric(msgl:::Err.msgl(f, type = "loglike")), fit = f)
      })
      best_a <- which.min(vapply(cvs, function(o) min(o$err), numeric(1)))
      best_l <- which.min(cvs[[best_a]]$err)
      lam <- cvs[[best_a]]$lambda
      fit <- msgl::fit(
        Xtr, y_tr,
        sampleWeights = w,
        grouping = groups,
        groupWeights = gw,
        parameterWeights = pw,
        alpha = settings$msgl_alpha[[best_a]],
        lambda = lam
      )
      pr <- stats::predict(fit, Xte)
      prob <- t(as.matrix(pr$response[[best_l]]))
      beta <- t(as.matrix(stats::coef(fit, best_l)))
      beta <- beta[intersect(genes, rownames(beta)), , drop = FALSE]
      nz <- which(beta != 0, arr.ind = TRUE)
      selected <- data.frame(
        gene = rownames(beta)[nz[, 1L]],
        type = colnames(beta)[nz[, 2L]],
        direction = ifelse(beta[nz] > 0, "up", "down"),
        stringsAsFactors = FALSE
      )
      imp <- stats::setNames(rep(0, length(genes)), genes)
      imp[rownames(beta)] <- apply(abs(beta), 1L, max)
      list(
        prob = prob,
        selected = selected,
        importance = data.frame(gene = genes, score = as.numeric(imp)),
        tuning = do.call(rbind, Map(function(o, a) {
          data.frame(alpha = a, lambda = o$lambda, cvm = o$err, nzero = NA_real_, chosen = a == settings$msgl_alpha[[best_a]] & seq_along(o$lambda) == if (a == settings$msgl_alpha[[best_a]]) best_l else 0L)
        }, cvs, settings$msgl_alpha)),
        diag = list(alpha = settings$msgl_alpha[[best_a]], lambda = lam[[best_l]], n_modules = settings$n_modules, module_sizes = as.numeric(table(groups)))
      )
    },
    plnlda = {
      mk <- function(rows) {
        d <- PLNmodels::prepare_data(
          counts = Y[rows, genes, drop = FALSE],
          covariates = data.frame(
            group = y[rows],
            cell_line = if (is.null(covariate)) 0 else covariate[rows],
            lib = lib[rows],
            row.names = rownames(Y)[rows]
          ),
          offset = "none"
        ) |> suppressWarnings()
        d$Offset <- lib[rows]
        d
      }
      d_tr <- mk(train)
      d_te <- mk(test)
      form <- if (is.null(covariate)) {
        Abundance ~ 0 + offset(log(Offset))
      } else {
        Abundance ~ 0 + cell_line + offset(log(Offset))
      }
      # do.call: PLNLDA() re-evaluates its own call, so pass values.
      fit <- do.call(PLNmodels::PLNLDA, list(
        formula = form,
        data = d_tr,
        weights = w,
        grouping = d_tr$group,
        control = PLNmodels::PLNLDA_param(covariance = settings$pln_covariance, trace = 0)
      ))
      prob <- stats::predict(fit, newdata = d_te, type = "posterior", scale = "prob")
      prob <- as.matrix(prob)[, types, drop = FALSE]
      cmap <- fit$corr_map
      imp <- stats::setNames(apply(abs(cmap), 1L, max), rownames(cmap))
      sc <- fit$scores
      keep <- withr::with_seed(settings$seed, sample(nrow(sc), min(nrow(sc), 1500L)))
      list(
        prob = prob,
        selected = top_assign(imp),
        importance = data.frame(gene = names(imp), score = as.numeric(imp)),
        tuning = data.frame(),
        diag = list(
          scores = data.frame(type = y_tr[keep], sc[keep, , drop = FALSE]),
          corr_map = cmap[names(utils::head(sort(imp, decreasing = TRUE), 30L)), , drop = FALSE],
          loglik = fit$loglik,
          percent_var = fit$percent_var
        )
      )
    },
    ranger = ,
    xgboost = {
      df_tr <- data.frame(as.matrix(Xtr), check.names = TRUE)
      df_te <- data.frame(as.matrix(Xte), check.names = TRUE)
      gmap <- stats::setNames(colnames(Xtr), names(df_tr))
      df_tr$.y <- y_tr
      p <- length(genes)
      rec <- recipes::recipe(.y ~ ., data = df_tr)
      if (method == "ranger") {
        cw <- stats::setNames(as.numeric(table(y_tr))^(-settings$gamma), types)
        cw <- cw / mean(cw[as.character(y_tr)])
        spec <- parsnip::rand_forest(trees = settings$rf_trees, mtry = parsnip::tune(), min_n = parsnip::tune()) |>
          parsnip::set_engine("ranger", class.weights = cw, num.threads = 1L, importance = "none") |>
          parsnip::set_mode("classification")
        grid <- expand.grid(mtry = unique(pmax(1L, round(sqrt(p) * settings$rf_mtry_mult))), min_n = settings$rf_min_n)
      } else {
        df_tr$.w <- hardhat::importance_weights(w)
        rec <- recipes::recipe(.y ~ ., data = df_tr)
        spec <- parsnip::boost_tree(
          trees = settings$xgb_trees,
          tree_depth = parsnip::tune(),
          learn_rate = parsnip::tune(),
          mtry = parsnip::tune()
        ) |>
          parsnip::set_engine("xgboost", nthread = 1L) |>
          parsnip::set_mode("classification")
        grid <- expand.grid(
          tree_depth = settings$xgb_depth,
          learn_rate = settings$xgb_learn_rate,
          mtry = unique(pmax(1L, round(sqrt(p) * settings$xgb_mtry_mult)))
        )
      }
      wf <- workflows::workflow() |> workflows::add_recipe(rec) |> workflows::add_model(spec)
      if (method == "xgboost") wf <- workflows::add_case_weights(wf, .w)
      splits <- lapply(sort(unique(inner)), function(f) {
        rsample::make_splits(list(analysis = which(inner != f), assessment = which(inner == f)), df_tr)
      })
      rs <- rsample::manual_rset(splits, paste0("inner", seq_along(splits)))
      # with_seed: tune draws its own seeds from .Random.seed.
      tuned <- withr::with_seed(settings$seed, tune::tune_grid(
        wf,
        resamples = rs,
        grid = grid,
        metrics = yardstick::metric_set(yardstick::f_meas),
        control = tune::control_grid(verbose = FALSE)
      ))
      best <- tune::select_best(tuned, metric = "f_meas")
      tun <- tune::collect_metrics(tuned)
      if (method == "ranger") {
        spec_f <- parsnip::rand_forest(trees = settings$rf_trees, mtry = best$mtry, min_n = best$min_n) |>
          parsnip::set_engine("ranger", class.weights = cw, num.threads = 1L, importance = "permutation") |>
          parsnip::set_mode("classification")
        wf_f <- workflows::update_model(wf, spec_f)
      } else {
        wf_f <- tune::finalize_workflow(wf, best)
      }
      fitted <- withr::with_seed(settings$seed + 1L, parsnip::fit(wf_f, df_tr))
      prob <- as.matrix(stats::predict(fitted, df_te, type = "prob"))
      colnames(prob) <- sub("^\\.pred_", "", colnames(prob))
      prob <- prob[, types, drop = FALSE]
      eng <- workflows::extract_fit_engine(fitted)
      if (method == "ranger") {
        imp <- eng$variable.importance
      } else {
        xi <- xgboost::xgb.importance(model = eng)
        imp <- stats::setNames(xi$Gain, xi$Feature)
      }
      names(imp) <- gmap[names(imp)]
      imp <- imp[names(imp) %in% genes]
      full <- stats::setNames(rep(0, length(genes)), genes)
      full[names(imp)] <- imp
      list(
        prob = prob,
        selected = top_assign(full),
        importance = data.frame(gene = genes, score = as.numeric(full)),
        tuning = as.data.frame(tun),
        diag = list(best = as.data.frame(best))
      )
    },
    splsda = {
      k <- length(types)
      ncomp <- k - 1L
      tune <- mixOmics::tune.splsda(
        as.matrix(Xtr), y_tr,
        ncomp = ncomp,
        test.keepX = settings$splsda_keepx,
        validation = "Mfold",
        folds = max(inner),
        nrepeat = 1L,
        dist = "max.dist",
        measure = "BER",
        seed = settings$seed
      )
      keepx <- tune$choice.keepX
      fit <- mixOmics::splsda(as.matrix(Xtr), y_tr, ncomp = ncomp, keepX = keepx)
      pr <- stats::predict(fit, as.matrix(Xte), dist = "max.dist")
      pred <- factor(pr$class$max.dist[, ncomp], levels = types)
      sel <- unique(unlist(lapply(seq_len(ncomp), function(cc) mixOmics::selectVar(fit, comp = cc)$name)))
      sel <- intersect(sel, genes)
      load <- fit$loadings$X[genes, , drop = FALSE]
      imp <- apply(abs(load), 1L, max)
      var_x <- fit$variates$X
      keep <- withr::with_seed(settings$seed, sample(nrow(var_x), min(nrow(var_x), 1500L)))
      list(
        prob = NULL,
        pred = pred,
        selected = assign_type(sel, contrast),
        importance = data.frame(gene = genes, score = imp),
        tuning = data.frame(
          comp = rep(seq_len(ncol(tune$error.rate)), each = nrow(tune$error.rate)),
          keepX = rep(as.numeric(rownames(tune$error.rate)), times = ncol(tune$error.rate)),
          ber = as.vector(tune$error.rate)
        ),
        diag = list(
          keepX = keepx,
          variates = data.frame(type = y_tr[keep], var_x[keep, , drop = FALSE]),
          loadings = load[sel, , drop = FALSE]
        )
      )
    },
    stop("Unknown method ", method, ".")
  )
}
