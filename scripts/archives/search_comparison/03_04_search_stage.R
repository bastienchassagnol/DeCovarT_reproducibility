# Archived 2026-11-17 from scripts/03_04_cv_covariance.R (stage `search`).
# Compared, for the shrinkage level, a 10-value refit grid, the exact path on
# 200 values, optimize() on the path score and 10 random values; for the
# two-parameter arm a 5 x 8 grid against 40 random pairs. Result: the path
# (a grid scored from one spectral decomposition) was best at no loss, so the
# pipeline keeps grid search only. Output kept: output/grn/cv/search_comparison.csv
# (read by docs/03-covariance-estimation.qmd). To rerun, paste this block back
# into scripts/03_04_cv_covariance.R before SECTION 5, restore "search" in
# stage_levels and `search_two_param_max_genes <- 1500L`, and the two lines of
# the tables stage that write search_comparison.csv from the search checkpoints.

# ==========================================================================
# SECTION 4 · Grid, path, optimiser and random search ----
# ==========================================================================

if ("search" %in% run_stages) {
  coarse_grid <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9)
  fine_grid <- seq(0.005, 1, length.out = 200L)
  alpha_search <- c(0, 0.25, 0.5, 0.75, 0.9)
  lambda_search <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9)
  for (tm in run_times) {
    for (gs in run_gene_sets) {
      path <- checkpoint("search", tm, gs)
      if (file.exists(path)) {
        next
      }
      X <- read_units(tm, gs)
      types <- names(X)
      n_genes <- ncol(X[[1L]])
      total <- if (n_genes <= search_two_param_max_genes) {
        Reduce(`+`, lapply(X, centred_scatter))
      }
      results <- parallel::mclapply(
        seq_along(X),
        function(j) {
          x <- X[[j]]
          folds <- make_cv_folds(nrow(x), K, seed + 5000L)
          timed <- function(expr) {
            t0 <- proc.time()[["elapsed"]]
            value <- force(expr)
            list(value = value, sec = proc.time()[["elapsed"]] - t0)
          }
          # Exact path: the reference every choice is scored against.
          prepared <- timed(lapply(seq_len(K), function(k) {
            sp <- split_standardise(
              x[folds != k, , drop = FALSE],
              x[folds == k, , drop = FALSE]
            )
            spectrum <- right_spectrum(sp$Z)
            w2 <- (sp$Zt %*% spectrum$vectors)^2
            list(
              base = sp$base,
              values = spectrum$values,
              w2 = w2,
              rest = pmax(rowSums(sp$Zt^2) - rowSums(w2), 0),
              null = ncol(sp$Zt) - length(spectrum$values)
            )
          }))
          score_at <- function(lam) {
            mean(vapply(
              prepared$value,
              function(f) {
                e <- (1 - lam) * f$values + lam
                f$base -
                  0.5 *
                    (sum(log(e)) +
                      f$null * log(lam) +
                      mean(drop(f$w2 %*% (1 / e)) + f$rest / lam))
              },
              numeric(1)
            ))
          }
          path_200 <- timed(vapply(fine_grid, score_at, numeric(1)))
          best <- max(path_200$value)
          refit <- timed({
            scores <- vapply(
              seq_len(K),
              function(k) {
                sp <- split_standardise(x[folds != k, , drop = FALSE])
                test <- x[folds == k, sp$varying, drop = FALSE]
                R <- crossprod(sp$Z)
                vapply(
                  coarse_grid,
                  function(lam) {
                    S <- (1 - lam) * R
                    diag(S) <- 1
                    heldout_gaussian_logscore(
                      sp$mu[sp$varying],
                      S * tcrossprod(sp$sd[sp$varying]),
                      test
                    )
                  },
                  numeric(1)
                )
              },
              numeric(length(coarse_grid))
            )
            coarse_grid[[which.max(rowMeans(scores))]]
          })
          optimised <- timed(
            stats::optimize(score_at, c(1e-3, 1), maximum = TRUE, tol = 1e-3)
          )
          random <- timed({
            draws <- withr::with_seed(seed + j, stats::runif(10L, 1e-3, 1))
            draws[[which.max(vapply(draws, score_at, numeric(1)))]]
          })
          row <- function(strategy, values, chosen, sec, extra = 0) {
            data.frame(
              time = tm,
              gene_set = gs,
              cell_type = types[[j]],
              cells = nrow(x),
              genes = n_genes,
              arm = "shrinkage",
              strategy = strategy,
              values_scored = values,
              alpha = 0,
              lambda = chosen,
              score_loss = best - score_at(chosen),
              elapsed_sec = sec + extra,
              stringsAsFactors = FALSE
            )
          }
          out <- rbind(
            row("grid, one factorisation per value", length(coarse_grid), refit$value, refit$sec),
            row("exact path", length(fine_grid), fine_grid[[which.max(path_200$value)]], path_200$sec, prepared$sec),
            row("optimiser on the path", NA_integer_, optimised$value$maximum, optimised$sec, prepared$sec),
            row("random values on the path", 10L, random$value, random$sec, prepared$sec)
          )
          if (!is.null(total)) {
            other <- total - centred_scatter(x)
            two_param_at <- function(alphas, lambdas_of) {
              # One eigen-decomposition per fold and distinct alpha.
              scores <- lapply(seq_len(K), function(k) {
                lapply(seq_along(alphas), function(a) {
                  score_unit_arms(
                    x[folds != k, , drop = FALSE],
                    x[folds == k, , drop = FALSE],
                    lambda_grid = lambdas_of(a),
                    alpha_grid = alphas[[a]],
                    other_scatter = other
                  )$two_param[1L, ]
                })
              })
              lapply(
                Reduce(function(a, b) Map(`+`, a, b), scores),
                function(s) s / K
              )
            }
            grid_run <- timed(two_param_at(alpha_search, function(a) lambda_search))
            grid_scores <- unlist(grid_run$value)
            grid_pairs <- expand.grid(lambda = lambda_search, alpha = alpha_search)
            draws <- withr::with_seed(
              seed + 10L * j,
              data.frame(
                alpha = stats::runif(40L, 0, 0.9),
                lambda = stats::runif(40L, 0.05, 0.9)
              )
            )
            random_run <- timed(two_param_at(
              draws$alpha,
              function(a) draws$lambda[[a]]
            ))
            random_scores <- unlist(random_run$value)
            top <- max(grid_scores, random_scores)
            two_row <- function(strategy, pairs, scores, sec) {
              i <- which.max(scores)
              data.frame(
                time = tm,
                gene_set = gs,
                cell_type = types[[j]],
                cells = nrow(x),
                genes = n_genes,
                arm = "two_param",
                strategy = strategy,
                values_scored = length(scores),
                alpha = pairs$alpha[[i]],
                lambda = pairs$lambda[[i]],
                score_loss = top - scores[[i]],
                elapsed_sec = sec,
                stringsAsFactors = FALSE
              )
            }
            out <- rbind(
              out,
              two_row("grid, path in lambda", grid_pairs, grid_scores, grid_run$sec),
              two_row("random pairs", draws, random_scores, random_run$sec)
            )
          }
          out
        },
        mc.cores = workers_for(length(X), n_genes),
        mc.preschedule = FALSE
      )
      check_tasks(results, paste("search", tm, gs))
      saveRDS(do.call(rbind, results), path)
      message("Search ", tm, " / ", gs, ": done.")
      rm(X, results)
      invisible(gc())
    }
  }
}

