"""T25: assemble the regenerated bacterial feature table (all Phase 6 steps) and diff it column by column with the shipped picard_collection.csv.
Needs the outputs of compare_features.py (T16-T21), run_lps_markers.py (T19), run_omp.py (T23) and run_proteome_tools.sh (T22).
python scripts/assemble_features.py   ->   replication/T25_features/picard_collection_regenerated.csv, results/replication/phase6_summary.csv"""
import glob, os
import numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = os.path.join(os.path.dirname(__file__), "..")
R, RES = f"{ROOT}/replication", f"{ROOT}/results/replication"
ship = pd.read_csv(f"{ROOT}/data/genomics/bacteria/picard_collection.csv", sep=";", index_col=0)
names = pd.read_csv(f"{R}/T14_panacota/annot/LSTINFO-list.lst", sep="\t").assign(strain=lambda d: d.orig_name.str.replace(r"\.fasta$", "", regex=True)).set_index("gembase_name").strain

tab = pd.read_csv(f"{R}/T25_features/picard_collection_regenerated.csv", sep=";", index_col=0)           # O/H, ST, phylogroup, Klebsiella, ABC model
tab["LPS_type (cross-validated markers)"] = pd.read_csv(f"{R}/T19_lps/lps_predictions.csv", index_col=0)["predicted_cv"]
omp = pd.read_csv(f"{R}/T23_omp/omp_clusters_regenerated.tsv", sep="\t", index_col=0)
tab = tab.join(omp.add_prefix("OMP_cluster_"))
nd = {}
for f in glob.glob(f"{R}/T22_defense/runs/*/*_defense_finder_systems.tsv"):
    gem = os.path.basename(f).replace("_defense_finder_systems.tsv", "")
    nd[names[gem]] = len(pd.read_csv(f, sep="\t"))
tab["n_defense_systems"] = pd.Series(nd)
tab = tab.reindex(ship.index)
tab.rename_axis("bacteria").to_csv(f"{R}/T25_features/picard_collection_regenerated.csv", sep=";")

rows = []
add = lambda task, feature, metric, value, note="": rows.append(dict(task=task, feature=feature, metric=metric, value=value, note=note))
t = pd.read_csv(f"{RES}/tree_agreement.csv").iloc[0]
add("T14", "core-genome tree", "normalised Robinson-Foulds (0 = identical)", t.normalised_RF, f"{int(t.shared_splits)}/{int(t.shipped_splits)} splits shared, {int(t.common_strains)} strains")
add("T14", "core-genome tree", "patristic distance Pearson", t.patristic_pearson)
add("T14", "CV groups (<1e-4)", "adjusted Rand index", t.cv_groups_ARI, f"{int(t.cv_groups_ours)} vs {int(t.cv_groups_shipped)} groups")
add("T15", "UMAP 8-D", "10-NN overlap with shipped coordinates", t.umap_knn_overlap_with_shipped, "UMAP run-to-run overlap is 0.77-0.80")
for _, r in pd.read_csv(f"{RES}/feature_agreement.csv").iterrows():
    add("T16-T21", r.feature, "agreement % (strains)", r.agreement_pct, f"{int(r.n_agree)}/{int(r.n_strains)}")
lp = pd.read_csv(f"{RES}/lps_markers.csv").set_index("type")
add("T19", "LPS_type (pangenome marker families, 5-fold CV)", "accuracy %", lp.loc["ALL", "cv_recall_pct"], "supervised by shipped labels; K-12 waaL BLAST alone finds only K12/R1/R4")
for _, r in pd.read_csv(f"{RES}/omp_agreement.csv").iterrows():
    add("T23", f"OMP {r.gene} clusters (99%)", "adjusted Rand index", r.ARI, f"{int(r.ours_clusters)} vs {int(r.shipped_clusters)} clusters")
d = pd.read_csv(f"{RES}/defense_agreement.csv").iloc[0]
add("T22", "defense subtypes (presence)", "F1", d.presence_F1, f"precision {d.presence_precision}, recall {d.presence_recall}")
x, y = ship["n_defense_systems"], tab["n_defense_systems"]
ok = x.notna() & y.notna()
add("T22", "n_defense_systems per strain", "Pearson vs shipped column", round(np.corrcoef(x[ok], y[ok])[0, 1], 3), f"mean {y.mean():.2f} ours vs {x.mean():.2f} shipped; exact {100 * (x[ok] == y[ok]).mean():.1f}%")
h = pd.read_csv(f"{RES}/host_table_agreement.csv")
for _, r in h.iterrows():
    add("T24", f"host strains: {r.column}", "agreement % (strains)", r.agreement_pct, f"{int(r.n_agree)}/{int(r.n_strains)}")
pd.DataFrame(rows).to_csv(f"{RES}/phase6_summary.csv", index=False)
print(pd.DataFrame(rows)[["task", "feature", "metric", "value"]].to_string(index=False))
print("\nregenerated table:", tab.shape, "columns:", list(tab.columns))
