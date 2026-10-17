"""T24: regenerate the table of the 34 phage-isolation (host) strains from our own tool outputs and compare it with the shipped
data/genomics/bacteria/isolation_strains/isolation_strain_features.csv. The authors' notebook only merges tool outputs (ECTyper, mlst,
ClermonTyping, waaL typing), so the regenerated table is the merge of our T16-T19 results.   python scripts/build_host_table.py"""
import os, re
import numpy as np, pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..")
ship = pd.read_csv(f"{ROOT}/data/genomics/bacteria/isolation_strains/isolation_strain_features.csv", sep=";", index_col=0)
feat = pd.read_csv(f"{ROOT}/replication/T25_features/picard_collection_regenerated.csv", sep=";", index_col=0)
lps = pd.read_csv(f"{ROOT}/replication/T19_lps/lps_predictions.csv", index_col=0)["predicted_cv"]
ours = feat[["Clermont_Phylo", "ST_Warwick", "O-type", "H-type"]].join(lps.rename("LPS_type")).reindex(ship.index)
ours.rename_axis("bacteria").to_csv(f"{ROOT}/replication/T25_features/host_features_regenerated.csv", sep=";")
cn = lambda x: np.nan if pd.isna(x) else re.sub(r"^e\.|[\s.]", "", str(x).lower())
rows = []
for c in ours.columns:
    a, b = (x.astype("string").str.replace(r"\.0$", "", regex=True) for x in (ship[c], ours[c]))  # ST numbers arrive as floats
    if c == "Clermont_Phylo":
        a, b = a.map(cn), b.map(cn)
    ok = ((a == b) | (a.isna() & b.isna())).fillna(False).astype(bool)
    rows.append(dict(column=c, n_strains=len(ship), n_agree=int(ok.sum()), agreement_pct=round(100 * ok.mean(), 1), mismatches="; ".join(f"{i}: shipped {a[i]} / ours {b[i]}" for i in ship.index[~ok])))
os.makedirs(f"{ROOT}/results/replication", exist_ok=True)
pd.DataFrame(rows).to_csv(f"{ROOT}/results/replication/host_table_agreement.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
