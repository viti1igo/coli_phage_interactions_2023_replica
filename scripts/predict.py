"""Reproduce dev/predictions/predict_all_phages.py on the public data (one model set per phage, Group10Fold CV).
Faithful to the original incl. its quirks (see docs/FEATURES.md); run(drop=...) is reused by ablate.py.
    python scripts/predict.py            -> results/repro/{per_phage,predictions}.csv + comparison to published"""
import glob, os, re, sys, warnings
import numpy as np, pandas as pd, sklearn
from joblib import Parallel, delayed
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

ROOT = os.path.join(os.path.dirname(__file__), "..")
G = f"{ROOT}/data/genomics"
# FEATURE_SET=shipped (default) | regenerated (our features, tree-based UMAP and CV groups) | regenerated+matrix (also our MLC matrix where we have raw data)
SET = os.environ.get("FEATURE_SET", "shipped"); SUFFIX = "" if SET == "shipped" else "_" + SET.replace("+", "_")
# Audit switches (scripts/run_audit_runs.sh), all off by default: TAG writes to results/audit_runs/<TAG>; SEED seeds RF, Dummy and LogReg(saga);
# FOLDS=<csv phage,bacteria,fold> reuses the authors' folds on the pairs in that file, FOLDS=shuffle draws a random partition of the groups with SEED;
# SWAP=<O-type,ST,LPS,Klebsiella,ABC,UMAP,host> puts those regenerated feature groups into the shipped inputs.
TAG, SEED, FOLDS = os.environ.get("TAG"), os.environ.get("SEED"), os.environ.get("FOLDS")
SWAP = [g for g in os.environ.get("SWAP", "").split(",") if g]
RES = f"{ROOT}/results/audit_runs/{TAG}" if TAG else f"{ROOT}/results/repro{SUFFIX}"
RS = int(SEED) if SEED else 0
warnings.filterwarnings("ignore")

M = pd.read_csv(f"{ROOT}/data/interactions/interaction_matrix.csv", sep=";", index_col=0)
BACT = pd.read_csv(f"{G}/bacteria/picard_collection.csv", sep=";", index_col=0)
UMAP = pd.read_csv(f"{G}/bacteria/umap_phylogeny/coli_umap_8_dims.tsv", sep="\t", index_col=0)
GROUPS = pd.read_csv(f"{ROOT}/data/derived/cv_groups_1e-4.csv", sep=";", index_col=0)["group"]
def regenerated_tables():
    reg = pd.read_csv(f"{ROOT}/replication/T25_features/picard_collection_regenerated.csv", sep=";", index_col=0, dtype={"ST_Warwick": str})
    b = reg.rename(columns={"LPS_type (cross-validated markers)": "LPS_type", "ABC_model (best system)": "ABC_serotype"})
    b["ST_Warwick"] = b["ST_Warwick"].str.replace(r"^~|\*$", "", regex=True)   # mlst marks approximate / novel alleles
    return b, pd.read_csv(f"{ROOT}/replication/T14_panacota/umap_8_dims_our_tree.tsv", sep="\t", index_col=0)


if SET != "shipped":
    BACT, UMAP = regenerated_tables()
    GROUPS = pd.read_csv(f"{ROOT}/data/derived/cv_groups_1e-4_regenerated.csv", sep=";", index_col=0)["group"]
if SET == "regenerated+matrix":
    mine = pd.read_csv(f"{ROOT}/replication/T26_matrix/mlc_matrix_regenerated.csv", sep=";", index_col=0)
    M = mine.combine_first(M)  # our matrix where the raw table gives a value, the shipped one elsewhere (the 33 host strains, uninterpretable cells)
BACT_HOST = BACT   # table the host-sharing features are built from
if SWAP:
    rb, ru = regenerated_tables()
    rb["ST_Warwick"] = rb["ST_Warwick"].str.replace(r"\.0$", "", regex=True)   # shipped ST is '1163', regenerated '1163.0'
    BACT = BACT.copy()
    COLS = {"O-type": ["O-type"], "ST": ["ST_Warwick"], "LPS": ["LPS_type"], "Klebsiella": ["Klebs_capsule_type"], "ABC": ["ABC_serotype"]}
    for g in SWAP:
        if g == "UMAP": UMAP = ru
        elif g == "host": BACT_HOST = rb
        else: BACT[COLS[g]] = rb[COLS[g]].reindex(BACT.index)
    if "host" not in SWAP: BACT_HOST = BACT
Y = (M > 0).astype(float).where(M.notna())  # lytic = MLC > 0 (paper: MLC is null for non-lytic)
PHAGE_HOST = pd.read_csv(f"{G}/phages/guelin_collection.csv", sep=";", index_col=0)["Phage_host"]
FEAT = BACT[["O-type", "LPS_type", "ST_Warwick", "Klebs_capsule_type", "ABC_serotype"]].join(UMAP, how="inner")
UMAP_COLS = [c for c in FEAT if c.startswith("UMAP")]
# original regex dropped H-type and phylogroup before the host merge, so host features = ST_host, O_host only
# the authors' penalty="l1" is deprecated in scikit-learn >= 1.8 (l1_ratio=1 is the same solution); in 1.1.2 l1_ratio=1 without penalty="l1" silently fits L2
L1 = {"penalty": "l1"} if tuple(map(int, sklearn.__version__.split(".")[:2])) < (1, 8) else {"l1_ratio": 1}
MODELS = {
    "RF_250_3": lambda cw: RandomForestClassifier(max_depth=3, n_estimators=250, class_weight=cw, random_state=RS),
    "RF_250_6": lambda cw: RandomForestClassifier(max_depth=6, n_estimators=250, class_weight=cw, random_state=RS),
    "LogReg_l2": lambda cw: LogisticRegression(class_weight=cw, max_iter=10000),
    "LogReg_l1": lambda cw: LogisticRegression(class_weight=cw, **L1, solver="saga", max_iter=10000, random_state=int(SEED) if SEED else None),
    "Dummy": lambda cw: DummyClassifier(strategy="stratified", random_state=RS),
}


if FOLDS and FOLDS != "shuffle":
    _f = pd.read_csv(FOLDS)
    AUTHOR_FOLDS = {ph: g.set_index("bacteria").fold for ph, g in _f.groupby("phage")}


def class_weight(y):
    f = y.mean()
    return {0: 1, 1: 0.8 if f >= .6 else 1 if f >= .4 else 1.5 if f >= .3 else 2 if f >= .2 else 3}


def design(phage, drop=()):
    """Feature matrix for one phage, built as in the original script. drop = regexes on pre-dummy column names."""
    d = pd.DataFrame({"y": Y[phage]}).join(FEAT).join(GROUPS).dropna(subset=["y"])
    for c in ("O-type", "ST_Warwick"):  # lump categories seen <3 times (counted over all strains)
        rare = FEAT.groupby(c).filter(lambda x: len(x) < 3)[c].unique()
        d.loc[d[c].isin(rare), c] = "Other"
    if not phage.startswith("LF110"):  # no host genome for LF110
        h = BACT_HOST.loc[PHAGE_HOST[phage]]
        d["same_O_as_host"] = d["O-type"] == h["O-type"]
        d["same_ST_as_host"] = d["ST_Warwick"] == h["ST_Warwick"]
        d["same_ABC_as_host"] = d["ABC_serotype"] == d["ABC_serotype"]  # original bug: x==x -> "ABC serotype is not missing"
        d["same_O_and_ST_as_host"] = d["same_O_as_host"] * d["same_ST_as_host"]
        d["ST_host"] = h["ST_Warwick"]
    g, y = d.pop("group"), d.pop("y")
    d = d.drop(columns=[c for c in d if any(re.search(p, c) for p in drop)])
    num = [c for c in d if c in UMAP_COLS]
    X = pd.concat([(d[num] - d[num].mean()) / d[num].std(), pd.get_dummies(d.drop(columns=num))], axis=1)
    return X.astype(float), y, g


def run_phage(phage, drop=()):
    X, y, g = design(phage, drop)
    if FOLDS and FOLDS != "shuffle":
        f = AUTHOR_FOLDS[phage].reindex(y.index)
        keep = f.notna().to_numpy()   # pairs the authors' run had
        X, y, g, f = X[keep], y[keep], g[keep], f[keep].to_numpy()
        splits = [(np.where(f != k)[0], np.where(f == k)[0]) for k in sorted(set(f))]
    else:
        splits = GroupKFold(10, shuffle=True, random_state=RS).split(X, y, g) if FOLDS == "shuffle" else GroupKFold(10).split(X, y, g)
    cw = class_weight(y)
    perf, preds = [], []
    for k, (tr, te) in enumerate(splits):
        for name, make in MODELS.items():
            if y.iloc[tr].nunique() < 2:  # very rare phage: every lytic strain sits in the test fold, nothing to learn -> constant prediction
                p = np.full(len(te), float(y.iloc[tr].iloc[0]))
            else:
                p = make(cw).fit(X.iloc[tr], y.iloc[tr]).predict_proba(X.iloc[te])[:, 1]
            if y.iloc[te].nunique() > 1:
                perf.append((name, k, average_precision_score(y.iloc[te], p), roc_auc_score(y.iloc[te], p)))
            preds.append(pd.DataFrame({"model": name, "bacteria": y.index[te], "y": y.iloc[te].values, "p": p}))
    perf = pd.DataFrame(perf, columns=["model", "fold", "aupr", "auroc"])
    best = perf.groupby("model")["aupr"].mean().idxmax()  # original selection rule: mean test AUPR
    preds = pd.concat(preds)
    return (dict(phage=phage, n_pos=int(y.sum()), n=len(y), best_model=best,
                 aupr=perf[perf.model == best].aupr.mean(), auroc=perf[perf.model == best].auroc.mean()),
            preds[preds.model == best].assign(phage=phage))


def run(drop=(), phages=None, jobs=-1):
    out = Parallel(n_jobs=jobs)(delayed(run_phage)(p, drop) for p in (phages or Y.columns))
    return pd.DataFrame([o[0] for o in out]).set_index("phage"), pd.concat([o[1] for o in out])


def published():
    """Same two summaries from the authors' shipped performances/ and predictions/ (best model by mean test AUPR)."""
    rows, preds = [], []
    for f in glob.glob(f"{ROOT}/dev/predictions/results/performances/performance_*.csv"):
        d = pd.read_csv(f, sep=";")
        d = d[d.dataset == "test"].dropna(subset=["model"])
        best = d.groupby("model").avg_precision.mean().idxmax()
        b = d[d.model == best]
        rows.append(dict(phage=d.phage.iloc[0], best_model=best, aupr=b.avg_precision.mean(), auroc=b.roc_auc.mean()))
        p = pd.read_csv(glob.glob(f"{ROOT}/dev/predictions/results/predictions/predictions_{rows[-1]['phage']}_*.csv")[0], sep=";")
        preds.append(p[(p.dataset == "test") & (p.model == best)].rename(columns={"y_pred_proba": "p"}))
    return pd.DataFrame(rows).set_index("phage"), pd.concat(preds)


if __name__ == "__main__":
    out = RES
    os.makedirs(out, exist_ok=True)
    perf, preds = run()
    perf.to_csv(f"{out}/per_phage.csv")
    preds.to_csv(f"{out}/predictions.csv", index=False)
    pub, pub_preds = published()
    print(f"{'':10}{'mean per-phage AUROC':>22}{'mean per-phage AUPR':>22}{'pooled AUROC':>14}")
    print(f"{'ours':10}{perf.auroc.mean():22.3f}{perf.aupr.mean():22.3f}{roc_auc_score(preds.y, preds.p):14.3f}")
    print(f"{'published':10}{pub.auroc.mean():22.3f}{pub.aupr.mean():22.3f}{roc_auc_score(pub_preds.y, pub_preds.p):14.3f}")
    j = perf.join(pub, rsuffix="_pub")
    print(f"per-phage AUROC corr ours vs published: {j.auroc.corr(j.auroc_pub):.3f}; best-model agreement: {(j.best_model.str[:2] == j.best_model_pub.str[:2]).mean():.2f} (RF vs LogReg family)")
