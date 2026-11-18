# Rscript --no-save --no-restore scripts/utils/progress_report.R \
#   logs/<script>_<date>.progress.tsv [more files]
#
# Summary of background runs from their progress files (log_progress(),
# R/02_00_utils.R): per slot, tasks done out of planned, failures,
# median seconds per task and an estimated finish time from the median
# and the number of workers. Reads only.

files <- commandArgs(trailingOnly = TRUE)
if (length(files) == 0L) {
  files <- sort(list.files("logs", "\\.progress\\.tsv$", full.names = TRUE))
}
for (f in files) {
  lines <- Filter(
    function(x) length(x) >= 5L,
    strsplit(readLines(f, warn = FALSE), "\t", fixed = TRUE)
  )
  field <- function(kind, k) {
    vapply(Filter(function(x) x[[2L]] == kind, lines), function(x) c(x, "")[[k]], "")
  }
  plan <- data.frame(
    slot = field("plan", 3L),
    tasks = as.integer(field("plan", 4L)),
    workers = as.integer(field("plan", 5L)),
    started = field("plan", 1L)
  )
  plan <- plan[!duplicated(plan$slot, fromLast = TRUE), , drop = FALSE]
  done <- data.frame(
    time = field("done", 1L),
    slot = field("done", 3L),
    sec = as.numeric(field("done", 5L)),
    status = field("done", 6L)
  )
  cat("\n==", f, "\n")
  if (nrow(plan) == 0L) {
    cat("no plan line yet\n")
    next
  }
  rows <- lapply(seq_len(nrow(plan)), function(i) {
    d <- done[done$slot == plan$slot[[i]], , drop = FALSE]
    left <- max(plan$tasks[[i]] - nrow(d), 0L)
    med <- if (nrow(d) > 0L) stats::median(d$sec) else NA_real_
    eta <- if (left == 0L) "done" else if (is.na(med)) "?" else
      format(Sys.time() + ceiling(left / max(plan$workers[[i]], 1L)) * med, "%d %b %H:%M")
    data.frame(
      slot = plan$slot[[i]],
      done = paste0(nrow(d), "/", plan$tasks[[i]]),
      failed = sum(d$status != "ok"),
      median_sec = round(med),
      last = if (nrow(d) > 0L) substr(max(d$time), 6L, 16L) else "",
      eta = eta
    )
  })
  print(do.call(rbind, rows), row.names = FALSE)
}
