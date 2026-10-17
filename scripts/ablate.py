"""Which feature groups are unnecessary?  (1) drop-one / drop-combination under the same group-CV, (2) the authors'
published importances aggregated by group.   python scripts/ablate.py   -> results/ablation/{ablation,importance}.csv"""
import glob, os, re
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import roc_auc_score
import predict as P

OUT = f"{P.ROOT}/results/ablation{P.SUFFIX}"
os.makedirs(OUT, exist_ok=True)
GROUPS = {  # disjoint; regexes on column names (host_sharing = same_*_as_host, same_O_and_ST, ST_host)
    "UMAP_phylogeny": "^UMAP", "O_type": "^O-type", "ST": "^ST_Warwick", "LPS_outer_core": "^LPS_type",
    "Klebsiella_capsule": "^Klebs", "ABC_capsule": "^ABC_serotype", "host_sharing": "^same_|^ST_host",
}
# singles + the redundancy suspects (UMAP/O/ST all encode phylogeny) + the two weak groups together
VARIANTS = {g: (g,) for g in GROUPS} | {
    "UMAP+O_type": ("UMAP_phylogeny", "O_type"), "UMAP+ST": ("UMAP_phylogeny", "ST"), "O_type+ST": ("O_type", "ST"),
    "UMAP+O_type+ST": ("UMAP_phylogeny", "O_type", "ST"), "Klebsiella+host_sharing": ("Klebsiella_capsule", "host_sharing"),
}


def ablation():
    base = pd.read_csv(f"{P.RES}/per_phage.csv", index_col=0)
    bp = pd.read_csv(f"{P.RES}/predictions.csv")
    rows = [dict(dropped="none", auroc=base.auroc.mean(), aupr=base.aupr.mean(), pooled_auroc=roc_auc_score(bp.y, bp.p))]
    for name, gs in VARIANTS.items():
        perf, pred = P.run(drop=tuple(GROUPS[g] for g in gs))
        perf.to_csv(f"{OUT}/per_phage_drop_{name}.csv")
        d, dp = (perf.auroc - base.auroc).dropna(), (perf.aupr - base.aupr).dropna()
        rows.append(dict(dropped=name, auroc=perf.auroc.mean(), aupr=perf.aupr.mean(), pooled_auroc=roc_auc_score(pred.y, pred.p),
                         d_auroc=d.mean(), d_aupr=dp.mean(), wilcoxon_p=wilcoxon(d).pvalue, phages_hurt_gt_0p02=int((d < -0.02).sum())))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/ablation.csv", index=False)


def importance():
    out = []
    for f in glob.glob(f"{P.ROOT}/dev/predictions/results/feature_importances/*_feature_importance.csv"):
        d = pd.read_csv(f, sep=";")
        d["g"] = d["variable"].map(lambda v: next((g for g, rx in GROUPS.items() if re.search(rx, v)), "other"))
        share = d.assign(w=d["value"].abs()).groupby("g")["w"].sum()  # RF importance >=0; LR |coef|
        out.append((share / share.sum()).rename(d["phage"].iloc[0]))
    t = pd.concat(out, axis=1).T.fillna(0)
    s = pd.DataFrame({"mean_share": t.mean(), "median_share": t.median(), "phages_where_top": t.idxmax(axis=1).value_counts()}).fillna(0)
    s.sort_values("mean_share", ascending=False).to_csv(f"{OUT}/importance.csv")
    print(s.sort_values("mean_share", ascending=False).round(3))


if __name__ == "__main__":
    importance()
    ablation()
