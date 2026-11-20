#!/usr/bin/env bash
# Detached chain for the two-sided post hoc contrasts and the classifier
# benchmark of chapter 2. Start from the repository root:
#
#   setsid nohup bash scripts/run_02_twosided_classifiers.sh \
#     > logs/run_02_chain_$(date +%F).log 2>&1 < /dev/null &
#
# Each step writes its own log in logs/ and one line to
# logs/run_02_chain_status.tsv (start, end, exit code). A failed step does
# not stop the chain; rerun the chain to resume: finished classifier folds
# are checkpointed in output/feature_selection/classifiers/fits/ and are
# skipped, and the DEG step is skipped when deg_twosided.rds exists.
# Progress of the classifier folds: logs/02_06_stats_classifier_cv_<date>.progress.tsv.

set -u
cd "$(dirname "$0")/.."
mkdir -p logs
day=$(date +%F)
status=logs/run_02_chain_status.tsv
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

step() {
  local name=$1
  shift
  local log="logs/${name}_${day}.log"
  printf '%s\tstart\t%s\t%s\n' "$(date '+%F %T')" "$name" "$log" >> "$status"
  "$@" > "$log" 2>&1
  local code=$?
  printf '%s\tend\t%s\texit=%s\n' "$(date '+%F %T')" "$name" "$code" >> "$status"
  return 0
}

wait_for() {
  # Wait until no process matches the pattern (a run started earlier).
  while pgrep -f "$1" > /dev/null; do sleep 30; done
}

# 1. Two-sided contrasts (a run may already be going).
wait_for "file=scripts/02_04_stats_deg_contrasts.R"
if [ ! -f output/feature_selection/deg_twosided/deg_twosided.rds ]; then
  step 02_04_stats_deg_contrasts Rscript --no-save --no-restore scripts/02_04_stats_deg_contrasts.R
fi
# Figures are redrawn only when missing (delete them to redraw).
[ -f docs/images/feature_selection/fs_deg_twosided_volcano.png ] ||
  step 02_04_vis_deg_twosided Rscript --no-save --no-restore scripts/02_04_vis_deg_twosided.R
[ -f output/feature_selection/deg/figures/fs_deg_superiority_heatmap.rds ] ||
  step 02_04_vis_deg Rscript --no-save --no-restore scripts/02_04_vis_deg.R

# 2. Full benchmark: every time, method and 48 h sensitivity.
step 02_06_stats_classifier_cv Rscript --no-save --no-restore scripts/02_06_stats_classifier_cv.R stages=nested workers=20

# 3. Smoke test (20 % of the 72 h cells, started separately): its figures
#    check the plotting code before the full figures.
wait_for "file=scripts/02_06_stats_classifier_cv.R.*subsample"
step 02_06_vis_classifier_cv_smoke Rscript --no-save --no-restore scripts/02_06_vis_classifier_cv.R dir=classifiers_subsample

# 4. Retry folds killed or failed, then the tables and figures.
step 02_06_stats_classifier_cv_retry Rscript --no-save --no-restore scripts/02_06_stats_classifier_cv.R stages=nested workers=10
step 02_06_stats_classifier_cv_tables Rscript --no-save --no-restore scripts/02_06_stats_classifier_cv.R stages=tables
step 02_06_vis_classifier_cv Rscript --no-save --no-restore scripts/02_06_vis_classifier_cv.R
printf '%s\tchain done\n' "$(date '+%F %T')" >> "$status"
