"""T14 + T15: compare our PanACoTA/FastTree core-genome tree with the shipped one, then re-embed ours with UMAP.
 (a) normalised Robinson-Foulds distance (unrooted splits) on the common strains, (b) correlation of pairwise patristic distances,
 (c) CV groups (distance < 1e-4) of both trees: adjusted Rand index, (d) kNN overlap of our UMAP vs the shipped coordinates.
python scripts/compare_tree.py [our_tree.nwk]   ->  results/replication/tree_agreement.csv"""
import copy, glob, os, sys
import numpy as np, pandas as pd
from Bio import Phylo
from scipy.sparse.csgraph import connected_components
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import adjusted_rand_score
import umap_embed as U
from treeutil import patristic

ROOT = os.path.join(os.path.dirname(__file__), "..")
SHIP = U.SHIPPED_TREE
cands = sys.argv[1:] or sorted(glob.glob(f"{ROOT}/replication/T14_panacota/tree/*.nwk"))
ours_path = next(c for c in cands if os.path.isfile(c))
names = pd.read_csv(f"{ROOT}/replication/T14_panacota/annot/LSTINFO-list.lst", sep="\t").assign(strain=lambda d: d.orig_name.str.replace(r"\.fasta$", "", regex=True)).set_index("gembase_name").strain

A = Phylo.read(SHIP, "newick"); B = Phylo.read(ours_path, "newick")
for t in B.get_terminals():
    t.name = names.get(t.name, t.name)
common = sorted({l.name for l in A.get_terminals()} & {l.name for l in B.get_terminals()})
for t in (A, B):
    for l in [l for l in t.get_terminals() if l.name not in set(common)]:
        t.prune(l)


def splits(t):
    full = frozenset(l.name for l in t.get_terminals()); ref = min(full)
    out = set()
    for c in t.get_nonterminals():
        s = frozenset(l.name for l in c.get_terminals())
        if 1 < len(s) < len(full) - 1:
            out.add(s if ref not in s else full - s)
    return out


sa, sb = splits(A), splits(B)
rf = len(sa ^ sb) / (len(sa) + len(sb))
na, Da = patristic(A); nb, Db = patristic(B)
ia = {n: i for i, n in enumerate(na)}; order = [ia[n] for n in nb]
Da = Da[np.ix_(order, order)]
iu = np.triu_indices(len(nb), 1)
x, y = Da[iu], Db[iu]
grp = lambda D: connected_components(D < 1e-4, directed=False)[1]
ari = adjusted_rand_score(grp(Da), grp(Db))
Eo = U.embed(Db, 15, 0.1)  # our tree, default UMAP settings (calibrated: the paper's are not stated)
ship_umap = pd.read_csv(U.SHIPPED_UMAP, sep="\t", index_col=0).loc[nb].values
pd.DataFrame(Eo, index=pd.Index(nb, name="bacteria"), columns=[f"UMAP{i}" for i in range(8)]).to_csv(f"{ROOT}/replication/T14_panacota/umap_8_dims_our_tree.tsv", sep="\t")
res = pd.DataFrame([dict(common_strains=len(common), our_tree=os.path.basename(ours_path), shipped_splits=len(sa), our_splits=len(sb), shared_splits=len(sa & sb),
    normalised_RF=round(rf, 3), patristic_pearson=round(pearsonr(x, y)[0], 4), patristic_spearman=round(spearmanr(x, y)[0], 4),
    patristic_median_ratio_ours_over_shipped=round(float(np.median(y[x > 0] / x[x > 0])), 3), cv_groups_ours=len(set(grp(Db))), cv_groups_shipped=len(set(grp(Da))), cv_groups_ARI=round(ari, 4),
    umap_knn_overlap_with_shipped=round(U.overlap(Eo, ship_umap), 3))])
os.makedirs(f"{ROOT}/results/replication", exist_ok=True)
res.to_csv(f"{ROOT}/results/replication/tree_agreement.csv", index=False)
print(res.T.to_string())
