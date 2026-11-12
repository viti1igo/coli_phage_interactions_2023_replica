"""Phase 21 (T96): add strain names to biophage-reproduce predictions (its files hold none).
Repeats steps 1-5 of its run_training.py with its own functions, then the same deterministic GroupKFold split,
so test rows map to strains in the order it wrote them.
Run inside the biophage env from the run folder:  cd replication/T93_biophage/<env> && python ../../../scripts/biophage_names.py"""
import glob, os, sys
import pandas as pd
from sklearn.model_selection import GroupKFold

sys.path.insert(0, os.getcwd())
from training.cv_groups import assign_cv_groups, compute_cophenetic_distances, load_lstinfo_mapping, parse_newick_tree
from training.data_loader import load_bacterial_features, load_interaction_matrix, merge_datasets
from training.features import prepare_features

P = "pipeline_results"
m = load_interaction_matrix("interaction_matrix.csv", sep=";")
raw = pd.read_csv("interaction_matrix.csv", sep=";", index_col="bacteria")   # to flag the NaN cells it turned into 0
feat = load_bacterial_features(f"{P}/picard_collection.csv", sep=";")
umap = pd.read_csv(f"{P}/umap_coordinates.tsv", sep="\t", index_col="bacteria")
mp = load_lstinfo_mapping(f"{P}/panacota/annotate/LSTINFO-file_list.lst")
d = compute_cophenetic_distances(parse_newick_tree(f"{P}/panacota/tree/ESCO.nucl.grp.aln.fasttree_tree.nwk"))
d.index = [mp.get(n, n) for n in d.index]; d.columns = d.index
groups = assign_cv_groups(d, threshold=1e-4)
X, names = prepare_features(merge_datasets(feat, umap), rare_threshold=3)

out = []
for f in sorted(glob.glob("training_results/predictions/predictions_*_core_features_Group10Fold_CV.csv")):
    pr = pd.read_csv(f, sep=";")
    phage = pr.phage.iloc[0]
    y = m[phage]
    common = y.index.intersection(X.index)
    y, g = y.loc[common].astype(int), groups.loc[common, "group"]
    test = {k: common[te] for k, (_, te) in enumerate(GroupKFold(n_splits=10).split(X.loc[common], y, g))}
    pr["bacteria"] = None
    for (k, model), idx in pr.groupby(["fold", "model"], sort=False).groups.items():
        assert len(idx) == len(test[k]), (phage, k, model)
        pr.loc[idx, "bacteria"] = list(test[k])
    assert (pr.y_true.values == y.loc[pr.bacteria].values).all(), phage   # names line up with the labels it saved
    pr["label_was_nan"] = raw.loc[pr.bacteria, phage].isna().values
    out.append(pr)
out = pd.concat(out)
groups.to_csv("training_results/cv_groups_named.csv")
out.to_csv("training_results/predictions_named.csv", index=False)
print(len(out), "rows,", out.phage.nunique(), "phages,", X.shape[1], "encoded features:", ", ".join(c for c in names if not c.startswith("UMAP"))[:300])
