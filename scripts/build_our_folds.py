"""Our own fold table: which GroupKFold(10) fold held out each (phage, strain) pair, on our rebuilt CV groups (shipped feature set).
One row per pair, in the same layout as results/audit/authors_folds.csv (phage, bacteria, fold) plus the group and both scikit-learn versions.
Inputs are the outputs of scripts/audit_gaps.py folds (scikit-learn 1.9.1, with the group) and scripts/fold_version_check.py (1.9.1 and 1.1.2).
    .venv/bin/python scripts/build_our_folds.py"""
import os
import pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..")
A = f"{ROOT}/results/audit"
V = f"{ROOT}/results/audit_runs/fold_version"
k = ["phage", "bacteria"]
f191 = pd.read_csv(f"{A}/folds_shipped.csv").rename(columns={"fold": "fold_sklearn_1_9_1"})
f112 = pd.read_csv(f"{V}/sk112_shipped.csv").rename(columns={"fold": "fold_sklearn_1_1_2"})
chk = pd.read_csv(f"{V}/sk191_shipped.csv").rename(columns={"fold": "c"})
out = f191.merge(f112, on=k, validate="1:1").merge(chk, on=k, validate="1:1")
assert (out.fold_sklearn_1_9_1 == out.c).all() and len(out) == len(f191) == len(f112), "the two 1.9.1 runs disagree"
out = out.drop(columns="c")[k + ["group", "fold_sklearn_1_9_1", "fold_sklearn_1_1_2"]].sort_values(k)
out.to_csv(f"{A}/our_folds.csv", index=False)
print(len(out), "pairs,", out.phage.nunique(), "phages,", out.bacteria.nunique(), "strains,", out.group.nunique(), "groups;",
      f"same fold index in both versions: {(out.fold_sklearn_1_9_1 == out.fold_sklearn_1_1_2).mean():.3f}")
