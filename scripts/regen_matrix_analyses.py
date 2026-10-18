"""Redo the matrix-only analyses (Fig 2f ANOVA, Fig 3a PERMANOVA, ED2 nestedness and modularity, ED5, ED7 coverage) on the matrix regenerated from the raw plaque scores
(replication/T26_matrix). That matrix covers 369 strains (the 33 isolation strains were assayed separately and are not in the raw table), so each
analysis is shown three ways: the authoritative matrix with all strains, the same matrix restricted to the 369 strains, and the regenerated matrix.
python scripts/regen_matrix_analyses.py -> results/comparison/regenerated_matrix_analyses.csv"""
import numpy as np, pandas as pd
from scipy.stats import f_oneway
from matrix import load_matrix
from matrix_stats import coverage_curve, nested_modular, phage_distance, phage_models

ROOT = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..")
ph = pd.read_csv(f"{ROOT}/data/genomics/phages/guelin_collection.csv", sep=";", index_col=0)
reg = pd.read_csv(f"{ROOT}/replication/T26_matrix/mlc_matrix_regenerated.csv", sep=";", index_col=0); reg.index = reg.index.astype(str)
M = load_matrix()
phages = list(M.columns); reg = reg[phages]
strains = [s for s in reg.index if reg.loc[s].notna().any()]          # strains with raw plaque data
sets = {"all strains (authoritative matrix)": M, f"same {len(strains)} strains (authoritative matrix)": M.loc[strains], f"regenerated matrix ({len(strains)} strains)": reg.loc[strains]}
feat = pd.read_csv(f"{ROOT}/replication/T25_features/picard_collection_regenerated.csv", sep=";", index_col=0); feat.index = feat.index.astype(str)
ship = pd.read_csv(f"{ROOT}/data/genomics/bacteria/picard_collection.csv", sep=";", index_col=0); ship.index = ship.index.astype(str)
meta = ph.loc[phages, ["Genus", "Phage_host", "Morphotype"]].astype(str).rename(columns={"Phage_host": "Host"})
rows = []
for name, X in sets.items():
    hr = (X.fillna(0) > 0).sum() / X.notna().sum()
    F = f_oneway(*[hr[meta.Morphotype == m].values for m in ("Myoviridae", "Podoviridae", "Siphoviridae")])
    res = phage_models(phage_distance(X), meta)
    nm = nested_modular((X.fillna(0) > 0).values.astype(float))
    cov, b3 = coverage_curve(X.fillna(0) > 0)
    x = X["411_P1"].fillna(0) > 0                                   # PDP110_P1 in the paper's names
    k55_regen = feat.reindex(X.index).Klebs_capsule_type.fillna("").str.contains("K55"); k55_ship = ship.reindex(X.index).Klebs_capsule_type.fillna("").str.contains("K55")
    rows.append(dict(matrix=name, strains=len(X), anova_F=F.statistic, anova_P=F.pvalue,
                     adjR2_genus_after_host=res["genus x isolation strain"][1] - res["isolation strain"][1], adjR2_host=res["isolation strain"][1], adjR2_both=res["genus x isolation strain"][1],
                     NODF=nm["T0"], NODF_null=nm["nulls"].mean(), P_NODF=nm["p_nodf"], BRIM_Q=nm["Q0"], BRIM_modules=nm["cbest"], BRIM_null=nm["qn"].mean(), P_Q=nm["p_q"],
                     ED7_coverage_k1=cov[1], ED7_coverage_k3=cov[3], ED7_coverage_k4=cov[4], ED7_coverage_k25=cov[25], ED7_best_triplet=" + ".join(sorted(X.columns[list(b3)])),
                     PDP110_P1_infected=int(x.sum()), PDP110_P1_infected_and_K55_shipped_typing=int((x & k55_ship).sum()), K55_strains_shipped_typing=int(k55_ship.sum()),
                     PDP110_P1_infected_and_K55_regenerated_typing=int((x & k55_regen).sum()), K55_strains_regenerated_typing=int(k55_regen.sum())))
out = pd.DataFrame(rows); out.to_csv(f"{ROOT}/results/comparison/regenerated_matrix_analyses.csv", index=False)
print(out.round(3).T.to_string())
