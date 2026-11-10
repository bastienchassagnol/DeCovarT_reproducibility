# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_02_stats_cross_cell_consistency.R \
#   > "logs/02_02_stats_cross_cell_consistency_$(date +%F).log" 2>&1 &
#
# Tier (ii): glmGamPoi::glm_gp() fits of M0 and M1,
# four tests of the type effect, BH gate on the Freedman-Lane permuted
# Welch G (labels permuted within batch, which matters at 48 h where
# the two cell lines are pooled), rank by deviance R^2. Does not call
# Seurat::SCTransform().

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

root <- getwd()
suppressPackageStartupMessages(library(Matrix))
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "02_02_stats_cross_cell.R"))

n_perm <- 2000L
fdr <- 0.05
min_frac <- 0.05
n_select <- Inf # every gene passing the BH gate; no cap
gate_test <- "welch_perm"
seed_fit <- 20261101L
seed_null_labels <- 20261102L
seed_null_fit <- 20261103L

slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
out_dir <- file.path(root, "output", "feature_selection", "tier_02")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(slots_path)) {
  stop(
    "Missing ",
    slots_path,
    ". Run scripts/02_01_stats_quality_control.R first."
  )
}

prepared <- readRDS(slots_path)
times <- prepared$time_levels
sf <- lapply(prepared$slots, function(slot) {
  slot$library_size / stats::median(slot$library_size)
})

# Labels drawn within a stratum keep that stratum's type counts.
permute_within <- function(labels, strata) {
  out <- labels
  for (idx in split(seq_along(labels), strata)) {
    out[idx] <- labels[idx][sample.int(length(idx))]
  }
  out
}

# ==========================================================================
# SECTION 1 · Fits, agreement, calibration, rank ----
# ==========================================================================

tests <- lapply(times, function(tm) {
  slot <- prepared$slots[[tm]]
  compare_type_effect_tests(
    slot$counts,
    slot$cell_type,
    size_factors = sf[[tm]],
    n_perm = n_perm,
    seed = seed_fit,
    strata = slot$batch,
    fdr = fdr,
    min_frac = min_frac
  )
})
names(tests) <- times

gate_pass <- paste0(gate_test, "_pass")
tau_all <- lapply(tests, type_test_kendall)
tau_gated <- lapply(tests, function(res) {
  type_test_kendall(res, genes = res[[gate_pass]])
})

# Null by construction: one within-batch permutation of the labels.
null_tests <- lapply(times, function(tm) {
  slot <- prepared$slots[[tm]]
  null_labels <- withr::with_seed(
    seed_null_labels,
    permute_within(slot$cell_type, slot$batch)
  )
  res_null <- compare_type_effect_tests(
    slot$counts,
    null_labels,
    size_factors = sf[[tm]],
    n_perm = n_perm,
    seed = seed_null_fit,
    strata = slot$batch,
    fdr = fdr,
    min_frac = min_frac
  )
  keep <- c("gene", grep("_pval$", names(res_null), value = TRUE))
  res_null[res_null$pass_presence, keep, drop = FALSE]
})
names(null_tests) <- times
calibration <- lapply(times, function(tm) {
  out <- type_test_calibration(null_tests[[tm]])
  out$time_point <- tm
  out
})

# The rank reuses the fits and permutations above.
ranked <- lapply(times, function(tm) {
  rank_genes_omnibus(
    fdr = fdr,
    n_select = n_select,
    gate_test = gate_test,
    results = tests[[tm]]
  )
})
names(ranked) <- times

write_out <- function(x, name, row_names = FALSE) {
  utils::write.csv(
    x,
    file.path(out_dir, name),
    row.names = row_names,
    fileEncoding = "UTF-8"
  )
}

# Small summaries cited by the chapter; every full table is in the RDS.
for (tm in times) {
  r <- ranked[[tm]]
  sel <- r[r$selected %in% TRUE, , drop = FALSE]
  write_out(
    sel[
      order(sel$rank),
      c("gene", "rank", "deviance_r2", "pval", "adj_pval", "overdispersion")
    ],
    paste0("selected_", tm, ".csv")
  )
  write_out(
    as.data.frame(tau_all[[tm]]),
    paste0("kendall_", tm, ".csv"),
    row_names = TRUE
  )
}
calibration_tbl <- do.call(rbind, calibration)
write_out(calibration_tbl, "calibration.csv")

# Floor of the permutation p-value and the BH condition for the genes
# at the floor: all of them pass iff n_floor >= m / (fdr (n_perm + 1)).
floor_tbl <- do.call(
  rbind,
  lapply(times, function(tm) {
    r <- tests[[tm]]
    p <- r[[paste0(gate_test, "_pval")]]
    p <- p[!is.na(p)]
    data.frame(
      time_point = tm,
      genes_tested = length(p),
      p_floor = 1 / (n_perm + 1),
      genes_at_floor = sum(p <= 1 / (n_perm + 1) + 1e-12),
      floor_genes_needed = length(p) / (fdr * (n_perm + 1)),
      genes_passing = sum(r[[gate_pass]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
write_out(floor_tbl, "permutation_floor.csv")

panel <- unique(unlist(lapply(ranked, function(res) res$gene[res$selected])))
saveRDS(
  list(
    tests = tests,
    ranked = ranked,
    tau_all = tau_all,
    tau_gated = tau_gated,
    calibration = calibration_tbl,
    null_tests = null_tests,
    permutation_floor = floor_tbl,
    panel = panel,
    sf = sf
  ),
  file.path(out_dir, "cross_cell.rds")
)
message("Selected genes, any time: ", length(panel))
