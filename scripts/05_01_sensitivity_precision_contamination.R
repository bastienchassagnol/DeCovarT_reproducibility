# mkdir -p logs
# nohup Rscript --no-save --no-restore \
#   scripts/05_01_sensitivity_precision_contamination.R \
#   > "logs/05_01_sensitivity_$(date +%F)_omega.log" 2>&1 &

# Contaminate sparse precision matrices Omega_j from 03_01_* and re-run
# DeCovarT deconvolution. Strategies: DeCovarT vignette
# supp-S5-misspecification.qmd (description-only factors to port):
#   S5a observation noise (Gaussian / t / PLN / contaminated Gaussian)
#   S5b missing reference component (0 / 2 / 10%)
#   S5c copula or PLN inter-gene dependence
#   S5d spillover + missing type
#   S5e outlier genes
# Plus: perturb Omega_j itself (wrong covariance reference, S5
# "Wrong covariance reference" / Appendix S4 sample-size factor).

# ==========================================================================
# SECTION 0 · Dependencies and paths ----
# ==========================================================================

# Run from the repository root. Paths are relative to that root.

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
# seed <- 1L
# cores <- 2L

# ==========================================================================
# SECTION 1 · Load reference moments ----
# ==========================================================================
# TODO: mu, Sigma_j / Omega_j, bulk y, true or proxy p if available

# ==========================================================================
# SECTION 2 · Contamination operators (precision / covariance) ----
# ==========================================================================
# TODO: drop / add edges in Omega_j; shrink towards diagonal; mix with
#       another cell type's precision; use a smaller donor subset
# TODO: keep operators as pure functions: Omega -> Omega_tilde

# ==========================================================================
# SECTION 3 · Deconvolution sweep ----
# ==========================================================================
# TODO: DeCovarT::deconvolute_ratios() under each contamination
# TODO: baselines (NNLS / others) on the same genes and mixtures
# TODO: ADEMP metrics; write output/sensitivity/; optional forest plots

stop(
  "Template only: implement contamination operators on fitted Omega_j; ",
  "do not invent sensitivity tables."
)
