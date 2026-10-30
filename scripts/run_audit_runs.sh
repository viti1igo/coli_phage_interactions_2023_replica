#!/bin/bash
# Phase 18 attribution runs (T72a, T72, T73): predict.py with the audit switches; each run writes results/audit_runs/<tag>/{per_phage,predictions}.csv. Resumes: finished runs are skipped.
#   scripts/run_audit_runs.sh        (about 5 minutes per run on 8 cores; 27 runs)
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
AF=results/audit/authors_folds.csv
ALL="O-type,ST,LPS,Klebsiella,ABC,UMAP,host"
run() {   # tag, then VAR=value pairs
  tag=$1; shift
  [ -f results/audit_runs/$tag/predictions.csv ] && return
  s=$(date +%s)
  env TAG=$tag "$@" $PY scripts/predict.py > logs/audit_run_$tag.log 2>&1
  echo "$tag $(( $(date +%s) - s )) s  exit $?" >> logs/audit_runs.log
}
mkdir -p results/audit_runs
# T72a: the authors' folds, nothing else changed (saga unseeded as in the authors' script)
run fixed_default FOLDS=$AF
# T72: model seeds with the folds fixed, and seeds with a shuffled group partition
for s in 0 1 2 3 4; do run fixed_s$s FOLDS=$AF SEED=$s; done
for s in 0 1 2 3 4; do run shuffle_s$s FOLDS=shuffle SEED=$s; done
# T73: one regenerated feature group at a time into the shipped inputs, authors' folds, seed 0; "none" and "all" are the controls
run swap_none FOLDS=$AF SEED=0 SWAP=
for g in O-type ST LPS Klebsiella ABC UMAP host; do run swap_$g FOLDS=$AF SEED=0 SWAP=$g; done
run swap_all FOLDS=$AF SEED=0 SWAP=$ALL
# T73b: everything regenerated except one group (the single-group effects do not add up to the all-groups effect)
for g in O-type ST LPS Klebsiella ABC UMAP host; do
  rest=$(echo ",$ALL," | sed "s/,$g,/,/; s/^,//; s/,$//")
  run drop_$g FOLDS=$AF SEED=0 SWAP=$rest
done
echo done > logs/audit_runs.done
