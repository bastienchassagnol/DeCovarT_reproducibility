#!/usr/bin/env bash
# Detached chain for the classifier benchmark of chapter 2, after msgl was
# stopped (2026-11-22; see output/feature_selection/classifiers/msgl_runtime.csv).
# glmnet (25 folds) and msgl (72 h and 96 h primary) are already
# checkpointed. PLN-LDA, random forest, XGBoost and sPLS-DA run as four
# separate processes of 14 workers each. Start from the repository root:
#
#   setsid nohup bash scripts/run_02_classifiers_parallel.sh \
#     >> logs/run_02_chain_$(date +%F).log 2>&1 < /dev/null &
#
# Each step writes its own log in logs/ and one line to
# logs/run_02_chain_status.tsv (start, end, exit code). A failed step does
# not stop the chain; rerun the chain to resume: finished folds are
# checkpointed in output/feature_selection/classifiers/fits/ and skipped.
# Progress of the folds: logs/02_06_stats_classifier_cv_<date>.progress.tsv.

set -u
cd "$(dirname "$0")/.."
mkdir -p logs
day=$(date +%F)
status=logs/run_02_chain_status.tsv
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cv=scripts/02_06_stats_classifier_cv.R
vis=scripts/02_06_vis_classifier_cv.R
methods=(plnlda ranger xgboost splsda)

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

# 1. The four methods in parallel, plus the smoke figures (a check of the
#    plotting code on the 20 % subsample).
for m in "${methods[@]}"; do
  step "02_06_clf_${m}" Rscript --no-save --no-restore "$cv" stages=nested methods="$m" workers=14 &
done
step 02_06_vis_classifier_cv_smoke Rscript --no-save --no-restore "$vis" dir=classifiers_subsample &
wait

# 2. Retry folds that failed or were killed, one method after another.
step 02_06_clf_retry Rscript --no-save --no-restore "$cv" stages=nested \
  methods="$(IFS=,; echo "${methods[*]}")" workers=14

# 3. Tables over every checkpoint (glmnet, msgl and the four methods), then figures.
step 02_06_stats_classifier_cv_tables Rscript --no-save --no-restore "$cv" stages=tables
step 02_06_vis_classifier_cv Rscript --no-save --no-restore "$vis"
printf '%s\tchain done\n' "$(date '+%F %T')" >> "$status"
