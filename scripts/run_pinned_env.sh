#!/bin/bash
# Phase 20 (T81-T82): the model step in an environment as close as possible to the authors' (scikit-learn 1.1.2).
#   scripts/setup_pinned_env.sh   builds the env (micromamba env ml_pinned)
#   scripts/run_pinned_env.sh     T2.1 (shipped inputs), T2.1 with the authors' folds, T3.1 (regenerated features); resumes
# Output: results/audit_runs/<tag>/{per_phage,predictions}.csv, logs/pinned_<tag>.log, logs/pinned_runs.log
cd "$(dirname "$0")/.." || exit 1
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba-phage}
AF=results/audit/authors_folds.csv
run() {   # tag, then VAR=value pairs
  tag=$1; shift
  [ -f results/audit_runs/$tag/predictions.csv ] && return
  s=$(date +%s)
  env TAG=$tag "$@" micromamba run -n ml_pinned python scripts/predict.py > logs/pinned_$tag.log 2>&1
  echo "$tag $(( $(date +%s) - s )) s  exit $?" >> logs/pinned_runs.log
}
mkdir -p results/audit_runs
run T2_1_shipped FEATURE_SET=shipped
run T2_1_shipped_authors_folds FEATURE_SET=shipped FOLDS=$AF
run T3_1_regenerated FEATURE_SET=regenerated
echo done > logs/pinned_runs.done
