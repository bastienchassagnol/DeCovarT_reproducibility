# Run from the repository root:
#   Rscript --no-save --no-restore scripts/ci_sysreqs.R install
#   Rscript --no-save --no-restore scripts/ci_sysreqs.R check
#
# install — apt packages required by renv.lock (pak::pkg_sysreqs), before restore
# check   — pak::sysreqs_check_installed() after restore

# ==========================================================================
# SECTION 0 · CLI ----
# ==========================================================================

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[[1L]] else "install"
if (!mode %in% c("install", "check")) {
  stop("Usage: scripts/ci_sysreqs.R install|check")
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

repos <- "https://cloud.r-project.org"
ensure_pak <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite", repos = repos)
  }
  if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak", repos = repos)
  }
}

# ==========================================================================
# SECTION 1 · Install system packages (pre-renv restore) ----
# ==========================================================================

install_sysreqs <- function() {
  ensure_pak()
  lock_path <- file.path("renv.lock")
  if (!file.exists(lock_path)) {
    stop("renv.lock not found; run from the repository root.")
  }

  lock <- jsonlite::fromJSON(lock_path, simplifyVector = FALSE)
  records <- lock$Packages
  pkg_names <- names(records)

  skip_pkgs <- c(
    "BLGGM",
    "MASS",
    "Matrix",
    "lattice",
    "codetools"
  )

  platform <- Sys.getenv("PAK_SYSREQS_PLATFORM", unset = "ubuntu")
  install_cmds <- character()

  for (pkg in pkg_names) {
    if (pkg %in% skip_pkgs) {
      next
    }
    rec <- records[[pkg]]
    remote_type <- rec$RemoteType %||% ""
    if (identical(remote_type, "github")) {
      next
    }

    sr <- tryCatch(
      pak::pkg_sysreqs(
        pkg,
        dependencies = FALSE,
        sysreqs_platform = platform
      ),
      error = function(err) {
        message("Skipping sysreqs for ", pkg, ": ", conditionMessage(err))
        NULL
      }
    )
    if (is.null(sr) || !length(sr$install_scripts)) {
      next
    }
    install_cmds <- c(install_cmds, sr$install_scripts)
  }

  install_cmds <- unique(install_cmds)
  if (!length(install_cmds)) {
    message("No system requirements reported for renv.lock.")
    return(invisible(TRUE))
  }

  message("System install commands:")
  for (cmd in install_cmds) {
    message("  ", cmd)
  }

  if (!nzchar(Sys.getenv("GITHUB_ACTIONS", unset = ""))) {
    message(
      "Not on GitHub Actions; run the apt-get commands above with sudo."
    )
    return(invisible(TRUE))
  }

  prefix_sudo <- function(cmd) {
    if (grepl("^sudo[[:space:]]", cmd)) {
      cmd
    } else {
      paste("sudo", cmd)
    }
  }

  status <- system("sudo apt-get update -y")
  if (!identical(status, 0L)) {
    stop("sudo apt-get update failed.")
  }

  for (cmd in install_cmds) {
    status <- system(prefix_sudo(cmd))
    if (!identical(status, 0L)) {
      stop("System requirement install failed: ", cmd)
    }
  }

  invisible(TRUE)
}

# ==========================================================================
# SECTION 2 · Verify installed R packages (post-renv restore) ----
# ==========================================================================

check_sysreqs <- function() {
  ensure_pak()
  result <- pak::sysreqs_check_installed()
  missing <- result[!result$installed, , drop = FALSE]

  if (nrow(missing) > 0L) {
    message("Missing system packages:")
    print(missing)
    stop(
      "Installed R packages are missing system requirements ",
      "(see pak::sysreqs_check_installed)."
    )
  }

  message("All system requirements for installed R packages are satisfied.")
  invisible(TRUE)
}

# ==========================================================================
# SECTION 3 · Dispatch ----
# ==========================================================================

if (identical(mode, "install")) {
  install_sysreqs()
} else {
  check_sysreqs()
}
