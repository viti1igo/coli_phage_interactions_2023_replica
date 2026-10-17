"""T19 (part 2): the K-12 WaaL query only recognises type K12 (R2 and R3 waaL are not homologous to it, and no per-type reference is shipped).
Data-driven substitute: find, in the PanACoTA pangenome (79k gene families), the family whose presence best marks each shipped LPS type,
and classify a strain by which markers it carries. Evaluated by 5-fold cross-validation (markers chosen on the training strains only).
python scripts/run_lps_markers.py   ->   results/replication/lps_markers.csv, lps_confusion.csv"""
import os, re
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold

ROOT = os.path.join(os.path.dirname(__file__), "..")
PAN = f"{ROOT}/replication/T14_panacota/pan/PanGenome-ESCO.All.prt-clust-0.8-mode1-th4.lst"
names = pd.read_csv(f"{ROOT}/replication/T14_panacota/annot/LSTINFO-list.lst", sep="\t").assign(strain=lambda d: d.orig_name.str.replace(r"\.fasta$", "", regex=True)).set_index("gembase_name").strain
ship = pd.read_csv(f"{ROOT}/data/genomics/bacteria/picard_collection.csv", sep=";", index_col=0)["LPS_type"].dropna()
strains = list(names.values); sidx = {s: i for i, s in enumerate(strains)}

fams = [l.split() for l in open(PAN)]
X = np.zeros((len(strains), len(fams)), dtype=np.uint8)
for j, f in enumerate(fams):
    for m in f[1:]:
        g = re.match(r"^(ESCO\.\d+\.\d+)\.", m)
        if g:
            X[sidx[names[g.group(1)]], j] = 1
y = ship.reindex(strains).dropna()
rows = np.array([sidx[s] for s in y.index]); Xy = X[rows]; types = sorted(y.unique())


def markers(Xt, yt):  # best-F1 family per type
    out = {}
    for t in types:
        pos = (yt == t).values
        tp = Xt[pos].sum(0); fp = Xt[~pos].sum(0); fn = pos.sum() - tp
        f1 = 2 * tp / np.maximum(2 * tp + fp + fn, 1)
        j = int(np.argmax(f1)); out[t] = (j, float(f1[j]))
    return out


def predict(Xs, mk):
    pred = []
    for r in Xs:
        c = [(f1, t) for t, (j, f1) in mk.items() if r[j] and t != "No_waaL"]
        pred.append(max(c)[1] if c else "No_waaL")
    return pred


full = markers(Xy, y)
pred = pd.Series("", index=y.index)
for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(Xy, y):
    mk = markers(Xy[tr], y.iloc[tr])
    pred.iloc[te] = predict(Xy[te], mk)
ok = pred == y
res = [dict(type=t, n=int((y == t).sum()), cv_recall_pct=round(100 * ok[y == t].mean(), 1), best_marker_family=fams[full[t][0]][0], marker_F1_all_strains=round(full[t][1], 3)) for t in types]
res.append(dict(type="ALL", n=len(y), cv_recall_pct=round(100 * ok.mean(), 1), best_marker_family="", marker_F1_all_strains=np.nan))
os.makedirs(f"{ROOT}/results/replication", exist_ok=True)
pd.DataFrame(res).to_csv(f"{ROOT}/results/replication/lps_markers.csv", index=False)
pd.crosstab(y, pred, rownames=["shipped"], colnames=["predicted"]).to_csv(f"{ROOT}/results/replication/lps_confusion.csv")
pd.DataFrame({"shipped": y, "predicted_cv": pred}).rename_axis("bacteria").to_csv(f"{ROOT}/replication/T19_lps/lps_predictions.csv")
print(pd.DataFrame(res).to_string(index=False)); print(pd.crosstab(y, pred, rownames=["shipped"], colnames=["predicted"]).to_string())
