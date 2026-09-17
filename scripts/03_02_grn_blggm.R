# mkdir -p logs
# nohup Rscript --no-save --no-restore scripts/03_02_grn_blggm.R \
#   n_celltype=4 num_iterations=2000 num_threads=2 \
#   input=data/processed/example.tsv \
#   > "logs/03_02_blggm_$(date +%F)_K4.log" 2>&1 &
#
# Run from the repository root. Flags are key=value.
#   n_celltype       mixture components K
#   num_iterations   MCMC iterations
#   num_threads      OpenMP threads (CI: at most 2)
#   seed             integer
#   input            genes x cells TSV
#   npn              shrinkage | truncation | none
#
# https://github.com/WgitU/BLGGM

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

n_celltype <- as.integer(cli_flag("n_celltype", "4"))
num_iterations <- as.integer(cli_flag("num_iterations", "2000"))
num_threads <- as.integer(cli_flag("num_threads", "2"))
seed <- as.integer(cli_flag("seed", "1"))
npn_fun <- cli_flag("npn", "none")
input_path <- cli_flag("input", "")

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
out_dir <- file.path("output", "grn", "blggm")

# ==========================================================================
# SECTION 1 · Inputs ----
# ==========================================================================

if (!nzchar(input_path) || !file.exists(input_path)) {
  stop(
    "No input matrix. Pass input=data/processed/genes_by_cells.tsv ",
    "from the repository root."
  )
}

set.seed(seed)
x_mat <- as.matrix(utils::read.delim(
  input_path,
  row.names = 1,
  check.names = FALSE
))
storage.mode(x_mat) <- "double"

if (nrow(x_mat) > ncol(x_mat)) {
  warning(
    "More rows than columns; assuming cells x genes and transposing ",
    "to genes x cells."
  )
  x_mat <- t(x_mat)
}

# ==========================================================================
# SECTION 2 · Optional NPN + BLGGM ----
# ==========================================================================

if (!identical(npn_fun, "none")) {
  x_cells <- t(x_mat)
  x_cells <- huge::huge.npn(x_cells, npn.func = npn_fun, verbose = TRUE)
  x_mat <- t(x_cells)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit <- BLGGM::BLGGM(
  x_mat,
  n_celltype,
  num_iterations = num_iterations,
  num_threads = num_threads
)

# ==========================================================================
# SECTION 3 · Outputs ----
# ==========================================================================

rds_file <- file.path(
  out_dir,
  paste0(
    study,
    "_blggm_K",
    n_celltype,
    "_iter",
    num_iterations,
    "_",
    today,
    ".rds"
  )
)
saveRDS(fit, rds_file)
message("Wrote ", rds_file)
