#!/bin/bash
# Re-run every model result with fixed seeds (predict.py + ablate.py on the three input sets); about 1 hour on 8 cores.
#   scripts/run_seeded_models.sh
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
for set in shipped regenerated regenerated+matrix; do
  FEATURE_SET=$set $PY scripts/predict.py > logs/predict_$set.log 2>&1
done
for set in shipped regenerated; do
  FEATURE_SET=$set $PY scripts/ablate.py > logs/ablate_$set.log 2>&1
done
echo done > logs/seeded_models.done
