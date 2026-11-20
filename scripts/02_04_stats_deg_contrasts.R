# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_04_stats_deg_contrasts.R \
#   > "logs/02_04_stats_deg_contrasts_$(date +%F).log" 2>&1 &
#
# Post hoc contrasts after the tier (ii) gate: which type drives the type
# effect, per slot (appendix of chapter 2). Reads the scored slots
# (tier (i)) and cross_cell.rds (tier (ii): size factors, presence gate,
# omnibus BH p-values, permutation settings). Writes per-slot contrast
# tables, DEG counts and checks to output/feature_selection/deg/
# (one-sided superiority: A and B) and output/feature_selection/deg_twosided/
# (two-sided: A2 and B2, same fits and permutation draw), plus the
# two-stage panel of the classifier benchmark (panel_twosided_<time>.csv).

root <- here::here()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "02_02_stats_cross_cell.R"))
source(file.path(root, "R", "02_04_stats_deg_contrasts.R"))

fdr <- 0.05
tests_used <- c("lr", "ql_f", "anova_vst", "welch_perm")
n_check_emmeans <- 25L
seed_check <- 20261120L
# Fold-change thresholds of the two-stage panel (|log2FC| of the type
# against the geometric mean of the others); drawn on the volcano, not tested.
lfc_thresholds <- c(1, log2(1.3))

tier1 <- file.path(root, "output", "feature_selection", "tier_01")
tier2 <- file.path(root, "output", "feature_selection", "tier_02")
out_dir <- file.path(root, "output", "feature_selection", "deg")
out_dir2 <- file.path(root, "output", "feature_selection", "deg_twosided")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)

prepared <- readRDS(file.path(tier1, "scored_slots.rds"))
cc <- readRDS(file.path(tier2, "cross_cell.rds"))
times <- names(cc$tests)

# DEG_REUSE=1 reuses the contrast p-values of an earlier deg.rds and only
# recomputes the decisions and checks (no refit, no permutation).
reuse <- nzchar(Sys.getenv("DEG_REUSE")) && file.exists(file.path(out_dir, "deg.rds"))
old_results <- if (reuse) readRDS(file.path(out_dir, "deg.rds"))$results else NULL

results <- list()
decisions <- list()
decisions2 <- list()
checks <- list()
checks2 <- list()

# Stage-wise decisions of every test for one list of p-values (p or p2).
decide <- function(res, pl, gate, genes) {
  rows <- list()
  for (tt in tests_used) {
    padj <- gate[[paste0(tt, "_padj")]][match(genes, gate$gene)]
    d <- stagewise_decisions(pl[[tt]]$p_a, pl[[tt]]$p_b, padj, fdr = fdr)
    d$gene <- genes[d$row]
    d$test <- tt
    d$p_a <- as.vector(pl[[tt]]$p_a)
    d$p_b <- as.vector(pl[[tt]]$p_b)
    d$log2fc <- as.vector(res$log2fc)
    if (!is.null(pl[[tt]]$sign_a)) {
      d$direction <- ifelse(as.vector(pl[[tt]]$sign_a) > 0, "up", "down")
    }
    d$row <- NULL
    rows[[tt]] <- d
  }
  rows
}
for (tm in times) {
  t0 <- Sys.time()
  slot <- prepared$slots[[tm]]
  gate <- cc$tests[[tm]]
  genes <- gate$gene[gate$pass_presence %in% TRUE]
  n_perm <- attr(gate, "n_perm")
  seed <- attr(gate, "seed")
  message(tm, ": ", length(genes), " genes, ", nlevels(droplevels(factor(slot$cell_type))), " types")
  res <- if (reuse) old_results[[tm]] else type_contrast_tests(
    slot$counts[genes, , drop = FALSE],
    slot$cell_type,
    size_factors = cc$sf[[tm]],
    strata = slot$batch,
    n_perm = n_perm,
    seed = seed,
    vst_alpha = attr(gate, "vst_alpha")
  )
  results[[tm]] <- res

  # Check 1: the same draw reproduces the gate's permutation p-value.
  gate_p <- gate$welch_perm_pval[match(genes, gate$gene)]
  same_draw <- isTRUE(all.equal(unname(res$welch_gate_p), gate_p))
  # Check 2: Wald p-values equal emmeans on the same coefficients.
  check_genes <- withr::with_seed(seed_check, sample(genes, n_check_emmeans))
  em <- emmeans_contrast_check(res, check_genes)

  rows <- decide(res, res$p, gate, genes)
  dec <- do.call(rbind, rows)
  dec$time_point <- tm
  decisions[[tm]] <- dec

  # Check 3: Holm on A equals stageR's stage-wise Holm (Welch test).
  w <- rows$welch_perm
  p_screen <- stats::setNames(gate$welch_perm_pval[match(genes, gate$gene)], genes)
  p_conf <- res$p$welch_perm$p_a
  sr <- stageR::stageR(pScreen = p_screen, pConfirmation = p_conf, pScreenAdjusted = FALSE)
  sr <- stageR::stageWiseAdjustment(sr, method = "holm", alpha = fdr, allowNA = TRUE)
  adj <- stageR::getAdjustedPValues(sr, onlySignificantGenes = FALSE, order = FALSE)
  stager_rej <- as.vector(adj[genes, res$types, drop = FALSE] <= fdr) %in% TRUE
  stager_agree <- mean(stager_rej == w$a_holm)

  checks[[tm]] <- data.frame(
    time_point = tm,
    genes_tested = length(genes),
    same_permutation_draw = same_draw,
    emmeans_max_abs_diff_a = max(em$max_abs_diff_a),
    emmeans_max_abs_diff_b = max(em$max_abs_diff_b),
    stager_holm_agreement = stager_agree,
    minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))
  )
  print(checks[[tm]])
  utils::write.csv(
    dec[, c(
      "time_point", "test", "gene", "type", "log2fc", "p_a", "p_b",
      "screened", "alpha_I", "a_graph", "b_graph", "a_holm", "a_flat_bh"
    )],
    file.path(out_dir, paste0("contrasts_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  # Two-sided: A2 and B2 on the same fits and draw.
  if (is.null(res$p2)) {
    stop("deg.rds has no two-sided p-values: rerun without DEG_REUSE.")
  }
  dec2 <- do.call(rbind, decide(res, res$p2, gate, genes))
  dec2$time_point <- tm
  decisions2[[tm]] <- dec2
  em2 <- emmeans_contrast_check(res, check_genes, side = "two-sided")
  checks2[[tm]] <- data.frame(
    time_point = tm,
    genes_tested = length(genes),
    same_permutation_draw = same_draw,
    emmeans_max_abs_diff_a2 = max(em2$max_abs_diff_a),
    emmeans_max_abs_diff_b2 = max(em2$max_abs_diff_b)
  )
  print(checks2[[tm]])
  utils::write.csv(
    dec2[, c(
      "time_point", "test", "gene", "type", "log2fc", "direction", "p_a",
      "p_b", "screened", "alpha_I", "a_graph", "b_graph", "a_holm",
      "a_flat_bh"
    )],
    file.path(out_dir2, paste0("contrasts_twosided_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  # Two-stage panel: primary family, A2 confirmed by the graph.
  pan <- dec2[dec2$test == "welch_perm" & dec2$a_graph %in% TRUE, ]
  pan <- pan[order(pan$type, -abs(pan$log2fc)), c(
    "time_point", "type", "gene", "log2fc", "direction", "p_a", "b_graph"
  )]
  for (thr in lfc_thresholds) {
    pan[[sprintf("abs_log2fc_ge_%.3f", thr)]] <- abs(pan$log2fc) >= thr
  }
  utils::write.csv(
    pan,
    file.path(out_dir2, paste0("panel_twosided_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

all_dec <- do.call(rbind, decisions)
counts_tbl <- stats::aggregate(
  cbind(a_flat_bh, a_holm, a_graph, b_graph) ~ time_point + test + type,
  data = all_dec,
  FUN = sum
)
utils::write.csv(counts_tbl, file.path(out_dir, "deg_counts.csv"), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(do.call(rbind, checks), file.path(out_dir, "checks.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(
  list(results = results, decisions = all_dec, counts = counts_tbl, checks = do.call(rbind, checks), fdr = fdr),
  file.path(out_dir, "deg.rds")
)
utils::capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
message("Wrote ", out_dir)

all_dec2 <- do.call(rbind, decisions2)
counts2 <- stats::aggregate(
  cbind(a_flat_bh, a_holm, a_graph, b_graph) ~ time_point + test + type + direction,
  data = all_dec2,
  FUN = sum
)
utils::write.csv(counts2, file.path(out_dir2, "deg_counts_twosided.csv"), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(do.call(rbind, checks2), file.path(out_dir2, "checks_twosided.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(
  list(
    decisions = all_dec2,
    counts = counts2,
    checks = do.call(rbind, checks2),
    fdr = fdr,
    lfc_thresholds = lfc_thresholds
  ),
  file.path(out_dir2, "deg_twosided.rds")
)
utils::capture.output(sessionInfo(), file = file.path(out_dir2, "session_info.txt"))
message("Wrote ", out_dir2)
