# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/04_02_feature_selection_genetic_algorithm.R \
#   > "logs/04_02_ga_$(date +%F)_panel.log" 2>&1 &

# Advanced panel search after the DESeq2 shortlist.
# See DeCovarT vignette supp-S6-feature-selection.qmd:
# AutoGeneS-style all-vs-all objectives; NSGA-II via mco::nsga2();
# optional scalar refinement with GA / GenSA.
# Minimise kappa_2(mu_G) and mean overlap; keep several Pareto panels.

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

# Run from the repository root. Paths are relative to that root.

stopifnot(requireNamespace("MixSim", quietly = TRUE))

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
# seed <- 1L

# ==========================================================================
# SECTION 1 · Candidate universe ----
# ==========================================================================
# TODO: load mu (and optional Sigma) on the DESeq2 shortlist
# TODO: encode a continuous mask in [0, 1]^|G0|, threshold to panel size

# ==========================================================================
# SECTION 2 · Multi-objective search ----
# ==========================================================================
# TODO: fitness <- function(mask) c(kappa2, mean_overlap)
# TODO: mco::nsga2(...)  # template only; do not invent a Pareto front
# TODO: optional greedy / Fedorov exchange as auditable baseline

# ==========================================================================
# SECTION 3 · Outputs ----
# ==========================================================================
# TODO: write selected gene panels + score table under data/processed/

stop(
  "Template only: implement NSGA-II / GA on a real shortlist; ",
  "do not invent gene panels."
)
