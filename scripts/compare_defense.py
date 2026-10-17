"""T22: compare our DefenseFinder 1.2.0 (models 1.2.3) systems with the shipped per-strain subtype table
(data/genomics/bacteria/defense_finder/370+host_defense_systems_subtypes.csv).   python scripts/compare_defense.py
Compared on the strains we have finished: (a) presence of each subtype, (b) systems per strain, (c) exact subtype-count vectors."""
import glob, os
import numpy as np, pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..")
names = pd.read_csv(f"{ROOT}/replication/T14_panacota/annot/LSTINFO-list.lst", sep="\t").assign(strain=lambda d: d.orig_name.str.replace(r"\.fasta$", "", regex=True)).set_index("gembase_name").strain
ship = pd.read_csv(f"{ROOT}/data/genomics/bacteria/defense_finder/370+host_defense_systems_subtypes.csv", sep=";", index_col=0)
rows = {}
for f in glob.glob(f"{ROOT}/replication/T22_defense/runs/*/*_defense_finder_systems.tsv"):
    gem = os.path.basename(f).replace("_defense_finder_systems.tsv", "")
    d = pd.read_csv(f, sep="\t")
    rows[names[gem]] = d.subtype.value_counts() if len(d) else pd.Series(dtype=int)
ours = pd.DataFrame(rows).T.fillna(0).astype(int)
common = [s for s in ours.index if s in ship.index]
cols = sorted(set(ours.columns) | set(ship.columns))
A = ship.reindex(index=common, columns=cols).fillna(0).astype(int)
B = ours.reindex(index=common, columns=cols).fillna(0).astype(int)
pa, pb = A > 0, B > 0
tp = (pa & pb).values.sum(); fp = (~pa & pb).values.sum(); fn = (pa & ~pb).values.sum()
res = pd.DataFrame([dict(n_strains=len(common), subtypes_shipped=int((ship.sum() > 0).sum()), subtypes_ours=int((ours.sum() > 0).sum()),
    presence_F1=round(2 * tp / (2 * tp + fp + fn), 3), presence_precision=round(tp / (tp + fp), 3), presence_recall=round(tp / (tp + fn), 3),
    systems_per_strain_pearson=round(np.corrcoef(A.sum(axis=1), B.sum(axis=1))[0, 1], 3) if len(common) > 2 else np.nan,
    mean_systems_shipped=round(A.sum(axis=1).mean(), 2), mean_systems_ours=round(B.sum(axis=1).mean(), 2),
    exact_vector_match_pct=round(100 * (A.values == B.values).all(axis=1).mean(), 1))])
os.makedirs(f"{ROOT}/results/replication", exist_ok=True)
res.to_csv(f"{ROOT}/results/replication/defense_agreement.csv", index=False)
diff = (B - A).sum()
diff[diff != 0].sort_values().to_frame("ours_minus_shipped").to_csv(f"{ROOT}/results/replication/defense_subtype_differences.csv")
print(res.T.to_string()); print(diff[diff != 0].sort_values().head(6).to_dict(), "...", diff[diff != 0].sort_values().tail(4).to_dict())
