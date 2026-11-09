# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/02_02_cross_cell_consistency.R \
#   > "logs/02_02_cross_cell_consistency_$(date +%F).log" 2>&1 &
#
# Written, not run. Tier (ii): glmGamPoi::glm_gp() fits of M0 and M1,
# four tests of the type effect, BH gate on the quasi-likelihood F,
# rank by deviance R^2. Does not call Seurat::SCTransform().

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

root <- getwd()
source(file.path(root, "R", "02_00_utils.R"))
source(file.path(root, "R", "02_01_stats_quality_control.R"))
source(file.path(root, "R", "02_02_cross_cell.R"))

n_perm <- 2000L
fdr <- 0.05
min_frac <- 0.05
n_select <- 500L
gate_test <- "ql_f"
seed_fit <- 20261101L
seed_null_labels <- 20261102L
seed_null_fit <- 20261103L

seu_path <- file.path(root, "data", "raw", "GSE229513_gastruloidsobject.rds")
slots_path <- file.path(
  root,
  "output",
  "feature_selection",
  "tier_01",
  "scored_slots.rds"
)
out_dir <- file.path(root, "output", "feature_selection", "tier_02")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(seu_path)) {
  stop(
    "Missing ", seu_path,
    ". From the repository root run: dvc pull ",
    "data/raw/GSE229513_gastruloidsobject.rds.dvc."
  )
}
if (!file.exists(slots_path)) {
  stop(
    "Missing ", slots_path,
    ". Run scripts/02_01_stats_quality_control.R first."
  )
}

prepared <- readRDS(slots_path)
times <- prepared$time_levels
sf <- lapply(prepared$slots, function(slot) {
  slot$library_size / stats::median(slot$library_size)
})

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
    fdr = fdr,
    min_frac = min_frac
  )
})
names(tests) <- times

tau_all <- lapply(tests, type_test_kendall)
tau_gated <- lapply(tests, function(res) {
  type_test_kendall(res, genes = res$ql_f_pass)
})

calibration <- lapply(times, function(tm) {
  slot <- prepared$slots[[tm]]
  null_labels <- withr::with_seed(
    seed_null_labels,
    sample(slot$cell_type)
  )
  res_null <- compare_type_effect_tests(
    slot$counts,
    null_labels,
    size_factors = sf[[tm]],
    n_perm = n_perm,
    seed = seed_null_fit,
    fdr = fdr,
    min_frac = min_frac
  )
  out <- type_test_calibration(res_null)
  out$time_point <- tm
  out
})

ranked <- lapply(times, function(tm) {
  slot <- prepared$slots[[tm]]
  rank_genes_omnibus(
    slot$counts,
    slot$cell_type,
    size_factors = sf[[tm]],
    fdr = fdr,
    min_frac = min_frac,
    n_select = n_select,
    gate_test = gate_test,
    n_perm = n_perm,
    seed = seed_fit
  )
})
names(ranked) <- times

for (tm in times) {
  utils::write.csv(
    tests[[tm]],
    file.path(out_dir, paste0("tests_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  utils::write.csv(
    ranked[[tm]],
    file.path(out_dir, paste0("ranked_", tm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  utils::write.csv(
    as.data.frame(tau_all[[tm]]),
    file.path(out_dir, paste0("kendall_all_", tm, ".csv")),
    fileEncoding = "UTF-8"
  )
  utils::write.csv(
    as.data.frame(tau_gated[[tm]]),
    file.path(out_dir, paste0("kendall_gated_", tm, ".csv")),
    fileEncoding = "UTF-8"
  )
}
calibration_tbl <- do.call(rbind, calibration)
utils::write.csv(
  calibration_tbl,
  file.path(out_dir, "calibration.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

panel <- unique(unlist(lapply(ranked, function(res) res$gene[res$selected])))
saveRDS(
  list(tests = tests, ranked = ranked, panel = panel, sf = sf),
  file.path(out_dir, "cross_cell.rds")
)
message("Selected genes, any time: ", length(panel))
