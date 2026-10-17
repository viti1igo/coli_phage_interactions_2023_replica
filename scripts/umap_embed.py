"""Tree -> patristic distances -> 8-D UMAP (the paper's 'in-house script' is not shipped; UMAP settings are unspecified).
  python scripts/umap_embed.py calibrate                 # grid over UMAP settings on the shipped tree, scored against the shipped coordinates
  python scripts/umap_embed.py embed <tree.nwk> <out.tsv> [n_neighbors min_dist]"""
import itertools, os, sys
import numpy as np, pandas as pd
import umap
from sklearn.manifold import trustworthiness
from sklearn.neighbors import NearestNeighbors
from treeutil import patristic

ROOT = os.path.join(os.path.dirname(__file__), "..")
SHIPPED_TREE = f"{ROOT}/data/genomics/bacteria/panacota/tree/370+host_tree_root=albertii.nwk"
SHIPPED_UMAP = f"{ROOT}/data/genomics/bacteria/umap_phylogeny/coli_umap_8_dims.tsv"
K = 10


def embed(D, n_neighbors, min_dist, seed=0):
    return umap.UMAP(n_components=8, metric="precomputed", n_neighbors=n_neighbors, min_dist=min_dist, random_state=seed).fit_transform(D)


def knn(X, k=K):
    return NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X, return_distance=False)[:, 1:]


def overlap(A, B):  # mean fraction of shared k nearest neighbours
    return float(np.mean([len(set(a) & set(b)) / K for a, b in zip(knn(A), knn(B))]))


if __name__ == "__main__":
    if sys.argv[1] == "calibrate":
        names, D = patristic(SHIPPED_TREE)
        ship = pd.read_csv(SHIPPED_UMAP, sep="\t", index_col=0).loc[names].values
        rows = [dict(n_neighbors="shipped", min_dist="-", seed="-", knn_overlap_with_shipped=1.0, trustworthiness_vs_tree=trustworthiness(D, ship, n_neighbors=K, metric="precomputed"))]
        for nn, md, seed in itertools.product([5, 10, 15, 30, 50], [0.0, 0.1, 0.5], [0, 1]):
            E = embed(D, nn, md, seed)
            rows.append(dict(n_neighbors=nn, min_dist=md, seed=seed, knn_overlap_with_shipped=overlap(E, ship),
                             trustworthiness_vs_tree=trustworthiness(D, E, n_neighbors=K, metric="precomputed")))
            print(rows[-1], flush=True)
        os.makedirs(f"{ROOT}/results/replication", exist_ok=True)
        pd.DataFrame(rows).to_csv(f"{ROOT}/results/replication/umap_calibration.csv", index=False)
    else:
        tree, out = sys.argv[2], sys.argv[3]
        nn, md = (int(sys.argv[4]), float(sys.argv[5])) if len(sys.argv) > 5 else (15, 0.1)
        names, D = patristic(tree)
        E = embed(D, nn, md)
        pd.DataFrame(E, index=pd.Index(names, name="bacteria"), columns=[f"UMAP{i}" for i in range(8)]).to_csv(out, sep="\t")
