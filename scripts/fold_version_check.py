"""Fold assignment of sklearn's GroupKFold(10) on our CV groups, per phage, for the installed scikit-learn (T81).
    python scripts/fold_version_check.py <out.csv> [shipped|regenerated]     (run once per environment, then scripts/audit_pinned.py compares them)"""
import os, sys
os.environ["FEATURE_SET"] = sys.argv[2] if len(sys.argv) > 2 else "shipped"
sys.path.insert(0, os.path.dirname(__file__))
import pandas as pd, sklearn
import predict as P
from sklearn.model_selection import GroupKFold
rows = []
for ph in P.Y.columns:
    X, y, g = P.design(ph)
    for k, (_, te) in enumerate(GroupKFold(10).split(X, y, g)):
        rows += [(ph, b, k) for b in y.index[te]]
pd.DataFrame(rows, columns=["phage", "bacteria", "fold"]).to_csv(sys.argv[1], index=False)
print(sklearn.__version__, len(rows), "pairs")
