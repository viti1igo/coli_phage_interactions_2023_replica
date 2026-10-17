"""T26: rebuild the MLC interaction matrix from data/interactions/raw/raw_interactions.csv with the paper's rule and compare with the shipped matrix.
Rule (Methods): MLC = the lowest phage concentration at which lysis is seen: 1 at 5e8 (log_dilution 0), 2 at 5e7 (-1), 3-4 at 5e6 (-2);
the 5e4 point (-4, replicate 1 only) is ignored; the final score averages the replicates. The raw table is binary (0/1/n), so plaque (3) vs
clearing (4) at 5e6 cannot be told apart: lysis at 5e6 is scored 3 here, and the comparison also reports "3 and 4 merged".
Only the 369 non-isolation strains are in the raw table (the 33 host strains were assayed separately).   python scripts/build_matrix.py"""
import os
import numpy as np, pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..")
raw = pd.read_csv(f"{ROOT}/data/interactions/raw/raw_interactions.csv", sep=";")
ship = pd.read_csv(f"{ROOT}/data/interactions/interaction_matrix.csv", sep=";", index_col=0)
raw["lyt"] = raw.score.astype(str).map({"0": 0.0, "1": 1.0})  # 'n' (not interpretable) -> NaN
raw = raw[raw.log_dilution.isin([0, -1, -2])]

# per replicate: score from the lowest dilution with lysis
w = raw.pivot_table(index=["bacteria", "phage", "replicate"], columns="log_dilution", values="lyt", aggfunc="max")
def rep_score(r):
    if r.isna().all(): return np.nan
    return 3 if r.get(-2) == 1 else 2 if r.get(-1) == 1 else 1 if r.get(0) == 1 else 0
per_rep = w.apply(rep_score, axis=1).unstack("replicate")           # bacteria-phage x replicate

# variant A: majority vote per dilution over the replicates, then the rule
vote = raw.groupby(["bacteria", "phage", "log_dilution"]).lyt.mean().unstack() >= 0.5
has = raw.groupby(["bacteria", "phage", "log_dilution"]).lyt.count().unstack() > 0
A = vote.apply(lambda r: 3 if r.get(-2) else 2 if r.get(-1) else 1 if r.get(0) else 0, axis=1).where(has.any(axis=1))
B = np.floor(per_rep.mean(axis=1) + 0.5)                              # mean of replicate scores, rounded half up
C = per_rep.median(axis=1).round()                                    # median
variants = {"max of replicate scores": per_rep.max(axis=1), "mean of replicate scores (rounded up)": np.ceil(per_rep.mean(axis=1) - 1e-9),
            "majority vote per dilution": A, "mean of replicate scores (rounded)": B, "median of replicate scores": C}

S = ship.stack(); S.index.names = ["bacteria", "phage"]
rows, best = [], None
for name, v in variants.items():
    d = pd.concat([v.rename("ours"), S.rename("ship")], axis=1).dropna()
    exact = (d.ours == d.ship.clip(upper=3)).mean() if False else None
    merged = d.ship.where(d.ship < 4, 3)                              # shipped 4 -> 3: plaque and clearing merged
    rows.append(dict(rule=name, cells=len(d), exact_0_1_2_and_3or4_merged_pct=round(100 * (d.ours == merged).mean(), 2),
                     lytic_vs_not_agreement_pct=round(100 * ((d.ours > 0) == (d.ship > 0)).mean(), 2),
                     exact_if_lysis_at_5e6_called_3_pct=round(100 * (d.ours == d.ship).mean(), 2), cells_missing_in_ours=int(S.dropna().index.difference(v.dropna().index).size)))
    if best is None or rows[-1]["exact_0_1_2_and_3or4_merged_pct"] > best[0]: best = (rows[-1]["exact_0_1_2_and_3or4_merged_pct"], name, v)
res = pd.DataFrame(rows)
os.makedirs(f"{ROOT}/results/replication", exist_ok=True); os.makedirs(f"{ROOT}/replication/T26_matrix", exist_ok=True)
res.to_csv(f"{ROOT}/results/replication/matrix_agreement.csv", index=False)
mat = best[2].unstack("phage").reindex(index=ship.index, columns=ship.columns)
mat.rename_axis("bacteria").to_csv(f"{ROOT}/replication/T26_matrix/mlc_matrix_regenerated.csv", sep=";")
d = pd.concat([best[2].rename("ours"), S.rename("ship")], axis=1).dropna()
pd.crosstab(d.ship, d.ours, rownames=["shipped MLC"], colnames=[f"ours ({best[1]})"]).to_csv(f"{ROOT}/results/replication/matrix_confusion.csv")
print(res.to_string(index=False)); print("\nbest rule:", best[1]); print(pd.crosstab(d.ship, d.ours, rownames=["shipped"], colnames=["ours"]).to_string())
print(f"\nregenerated matrix: {mat.notna().sum().sum()} cells filled of {mat.size}; strains without raw data: {int(mat.isna().all(axis=1).sum())}")
