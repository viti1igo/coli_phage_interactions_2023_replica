"""Rebuild the lost CV-group file: strains closer than 1e-4 subst/site (core-genome tree) share a group.
Then check it against the authors' own fold logs (no group of ours may be split across their folds)."""
import ast, glob, os, re, sys
import numpy as np, pandas as pd
from scipy.sparse.csgraph import connected_components
from treeutil import patristic

ROOT = os.path.join(os.path.dirname(__file__), "..")
TREE = f"{ROOT}/data/genomics/bacteria/panacota/tree/370+host_tree_root=albertii.nwk"
THRESH = 1e-4

CUSTOM = len(sys.argv) > 2   # usage: make_groups.py <our_tree.nwk> <out.csv>  (regenerated tree: gembase ids -> strain names, no fold-log check)
names, D = patristic(sys.argv[1] if CUSTOM else TREE)
if CUSTOM:
    lst = pd.read_csv(f"{ROOT}/replication/T14_panacota/annot/LSTINFO-list.lst", sep="\t").assign(strain=lambda d: d.orig_name.str.replace(r"\.fasta$", "", regex=True)).set_index("gembase_name").strain
    names = [lst.get(x, x) for x in names]
n = len(names)

_, lab = connected_components(D < THRESH, directed=False)  # single linkage
groups = pd.Series(lab, index=names, name="group").rename_axis("bacteria")
os.makedirs(f"{ROOT}/data/derived", exist_ok=True)
groups.to_csv(sys.argv[2] if CUSTOM else f"{ROOT}/data/derived/cv_groups_1e-4.csv", sep=";")
print(f"{n} strains -> {groups.nunique()} groups (max size {groups.value_counts().max()})")

if CUSTOM:
    sys.exit(0)
# --- check vs original logs: rows were sorted by bacteria, NaN pairs dropped, test_idx indexes those rows
m = pd.read_csv(f"{ROOT}/data/interactions/interaction_matrix.csv", sep=";", index_col=0).sort_index()
bad = tot = used = 0
for f in sorted(glob.glob(f"{ROOT}/dev/predictions/results/logs/logs_*_Group10Fold_CV.csv")):
    phage = re.match(r"logs_(.*)_Group10Fold_CV", os.path.basename(f)).group(1)
    rows = m.index[m[phage].notna()]
    lg = pd.read_csv(f, sep=";")
    if len(rows) != lg["train_size"][0] + lg["test_size"][0]:
        continue  # public matrix differs from the authors' private one for this phage: row mapping not reliable
    fold_of = {}
    for _, r in lg.iterrows():
        idx = [int(x) for x in re.findall(r"\d+", r["test_idx"].replace("\n", " "))]
        fold_of.update({rows[i]: r["fold"] for i in idx})
    g = groups.loc[list(fold_of)].to_frame().assign(fold=pd.Series(fold_of))
    bad += (g.groupby("group")["fold"].nunique() > 1).sum()
    tot += g["group"].nunique()
    used += 1
print(f"check ({used} phages with identical row count): {bad}/{tot} of our groups are split across the authors' folds (0 = consistent)")
