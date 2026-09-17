# mkdir -p logs
# nohup Rscript --no-save --no-restore scripts/03_01_grn_silggm.R \
#   npn=shrinkage method=D-S_GL global=TRUE alpha=0.05,0.10 \
#   input=data/processed/example.tsv \
#   > "logs/03_01_silggm_$(date +%F)_shrinkage.log" 2>&1 &
#
# Run from the repository root. Flags are key=value.
#   npn       shrinkage | truncation | skeptic
#   method    B_NW_SL | D-S_NW_SL | D-S_GL | GFC_SL | GFC_L
#   global    TRUE | FALSE
#   alpha     comma-separated FDR levels
#   seed      integer
#   input     cells x genes TSV
#
# https://github.com/hammodesharif-glitch/Mohamad_Al_charif_DTOO/blob/main/SILGGM/R/01_benchmark_silggm_dsgl.R

# ==========================================================================
# SECTION 0 · CLI and paths ----
# ==========================================================================

args <- commandArgs(trailingOnly = TRUE)
cli_flag <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    default
  } else {
    sub(paste0("^", name, "="), "", hit[[1L]], fixed = FALSE)
  }
}

npn_fun <- cli_flag("npn", "shrinkage")
silggm_method <- cli_flag("method", "D-S_GL")
global_flag <- identical(cli_flag("global", "TRUE"), "TRUE")
alpha_str <- cli_flag("alpha", "0.05,0.10")
alpha_levels <- as.numeric(strsplit(alpha_str, ",", fixed = TRUE)[[1L]])
seed <- as.integer(cli_flag("seed", "1"))
input_path <- cli_flag("input", "")

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
out_dir <- file.path("output", "grn", "silggm")

# ==========================================================================
# SECTION 1 · Inputs ----
# ==========================================================================

if (!nzchar(input_path) || !file.exists(input_path)) {
  stop(
    "No input matrix. Pass input=data/processed/cells_by_genes.tsv ",
    "from the repository root."
  )
}

set.seed(seed)
x_raw <- as.matrix(utils::read.delim(
  input_path,
  row.names = 1,
  check.names = FALSE
))
storage.mode(x_raw) <- "double"

# ==========================================================================
# SECTION 2 · Nonparanormal + SILGGM ----
# ==========================================================================

allowed_npn <- c("shrinkage", "truncation", "skeptic")
if (!npn_fun %in% allowed_npn) {
  stop("Unknown npn=", npn_fun, " (use shrinkage, truncation, or skeptic).")
}

if (identical(npn_fun, "skeptic")) {
  stop(
    "npn=skeptic returns a correlation matrix, not an n x p table. ",
    "Use shrinkage or truncation before SILGGM::SILGGM()."
  )
}

x_npn <- huge::huge.npn(x_raw, npn.func = npn_fun, verbose = TRUE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit <- SILGGM::SILGGM(
  x_npn,
  method = silggm_method,
  global = global_flag,
  alpha = alpha_levels,
  cytoscape_format = TRUE,
  csv_save = TRUE,
  directory = out_dir
)

# ==========================================================================
# SECTION 3 · Outputs ----
# ==========================================================================

rds_file <- file.path(
  out_dir,
  paste0(
    study,
    "_silggm_",
    gsub("[^A-Za-z0-9]+", "_", silggm_method),
    "_npn-",
    npn_fun,
    "_",
    today,
    ".rds"
  )
)
saveRDS(fit, rds_file)
message("Wrote ", rds_file)
