"""Phase 18 gap audit (T64-T67): gap register, tier alignment, overall metrics, paper against repository.
Tiers: T1 authors' shipped predictions, T2 our run on shipped inputs, T3 regenerated features, T4 regenerated features and matrix.
Outputs go to results/audit/.   python scripts/audit_gaps.py [register|paper_repo|align|overall|folds|phage|pairs|runs|final|external|preprint|all]; folds needs FEATURE_SET=shipped|regenerated|regenerated+matrix"""
import os, sys
import numpy as np, pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
sys.path.insert(0, os.path.dirname(__file__))
import predict as P

ROOT = os.path.abspath(P.ROOT)
R, OUT = f"{ROOT}/results", f"{ROOT}/results/audit"
os.makedirs(OUT, exist_ok=True)
TIER_DIRS = {"T2": "repro", "T3": "repro_regenerated", "T4": "repro_regenerated_matrix"}
B, SEED = 500, 0   # bootstrap replicates over strains


def tiers():
    """pair-level predictions of every tier: columns phage, bacteria, y, p"""
    _, ship = P.published()
    t = {"T1": ship[["phage", "bacteria", "y", "p"]]}
    for k, d in TIER_DIRS.items():
        t[k] = pd.read_csv(f"{R}/{d}/predictions.csv")[["phage", "bacteria", "y", "p"]]
    return t


# ------------------------------------------------------------------ T64 / T67: gap register
# (figure prefix, metric substring) -> (tier pair, cause code, basis, evidence, note); every unresolved row of paper_vs_replication.csv gets exactly one entry.
# Codes: P published materials disagree (paper, supplement, Source Data, repository); D definition or unit not stated;
# U input or code not published; M input data differs; F feature regeneration differs; S seed; C model choice rule; I our implementation choice;
# OK flagged only because the paper gives no number or a qualitative statement, and the numbers agree; X unexplained.
# basis: "reported" = evidence already in reports 05, 10, 11; "to test" = refined by T69, T71, T72 or T73.
RP = "reports/11_repository_vs_paper.md"
CODES = [
    ("Fig 4d", "pMCMC", "T0-T2", "X", "to test", f"{RP} #14", "beta equal (-0.032 with the panel-d counts); only pMCMC differs (0.0005 against 0.004); the panel-d counts are inconsistent with Fig 4c (P) but do not change beta"),
    ("ED1", "phages infecting", "T0-T2", "OK", "reported", "reports/05_paper_comparison.md", "45 and 48 equal the paper's counts; row marked qualitative only"),
    ("ED2", "nestedness (paper", "T0-T2", "I", "reported", "reports/05_paper_comparison.md", "NODF used instead of BINMATNEST; the paper's T and ours are not comparable"),
    ("ED2", "nestedness temperature", "T0-T2", "I", "reported", "reports/05_paper_comparison.md", "BINMATNEST not run; NODF used"),
    ("ED4A", "infected strains that are ABC K1", "T0-T2", "D", "reported", "reports/05_paper_comparison.md", "paper states a qualitative association with K1, no count"),
    ("ED5", "carry K55", "T2-T3", "F", "reported", "reports/03_feature_regeneration.md", "shipped types give 3 of 3; regenerated Kaptive types give 2 of 2 (K55 strains total 2)"),
    ("ED7", "baseline triplet", "T0-T2", "D", "reported", f"{RP} #20", "the paper does not name the triplet; ours is the most-covering set"),
    ("ED7", "coverage", "T0-T1", "D", "to test", f"{RP} #20", "no set of MLC scores gives the Source Data curve (best single phage 172 of 403 strains, any MLC > 0 gives 264); the published text (63 %) and the Source Data (67.0 %) also disagree (P); the 68 % for four phages is in the preprint only"),
    ("ED9A", "Mann-Whitney P", "T0-T1", "P", "reported", f"{RP} #19", "the notebook counts each strain three times (P = 3.5e-12); with each strain once P = 6.4e-5; the test unit is not stated (D)"),
    ("ED6A", "genome size vs lytic", "T0-T2", "D", "reported", "reports/05_paper_comparison.md", "the paper gives no value for this model; regenerated inputs give p = 0.55 against 0.12 (F)"),
    ("Fig 1j", "defense subtypes", "T0-T2", "D", "reported", "reports/05_paper_comparison.md", "no number stated"),
    ("Fig 3a", "genus + isolation strain", "T0-T2", "X", "to test", "reports/05_paper_comparison.md", "0.65 against 0.69; the repository has no PERMANOVA code"),
    ("Fig 3a", "genus alone", "T0-T2", "D", "reported", f"{RP} #18", "genus alone is not stated; 0.13 is the increment after isolation strain"),
    ("Fig 3a", "R2 (unadjusted)", "T0-T2", "D", "reported", f"{RP} #18", "unadjusted values are not stated"),
    ("Fig 5c", "Spearman", "T0-T2", "D", "reported", "reports/02_model_reproduction.md", "no value stated"),
    ("Fig 3b", "mean entropy in RBP", "T0-T2", "D", "reported", "reports/10_repository_analyses.md", "qualitative claim; ours supports it (0.62 against 0.24 bits); annotation not shipped"),
    ("Fig 3b", "mean entropy of the alignment", "T0-T2", "D", "reported", "reports/10_repository_analyses.md", "no value printed"),
    ("Fig 4d", "defense count in the Source Data", "T0-T1", "P", "reported", f"{RP} #14", "Source Data counts differ from Supplementary Table 1 for 353 of 400 strains"),
    ("Fig 4a", "largest pMCMC", "T0-T1", "P", "reported", f"{RP} #10", "Source Data keeps terms up to pMCMC 0.0089; Methods say 1e-3"),
    ("Fig 4a", "pMCMC < 0.001", "T0-T1", "P", "reported", f"{RP} #10", "threshold mismatch; at 0.01 the counts are 51 and 3 against 52 and 1"),
    ("Fig 4b", "significant in both", "T0-T2", "X", "to test", f"{RP} #11, #12", "19 of 54 pairs shared with the paper's 53; candidates are the unreproducible minimum-observations filter (D), the outer-membrane merge (P) and MCMC variation (S)"),
    ("Fig 6", "strains with tailored", "T0-T1", "P", "reported", f"{RP} #22", "Supplementary Table 7 gives 28 and 72; the Fig 6 analysis gives 24 and 76 (three or four generic cocktails)"),
    ("ED9b", "haplotypes covered", "T0-T1", "P", "reported", f"{RP} #21", "Source Data plots strain counts (31, 21, 19, 4); haplotype counts are 20, 18, 18, 4"),
    ("Fig 6 recommender", "identical cocktails", "T0-T2", "U", "reported", f"{RP} #27", "recommender code not published"),
    ("Fig 6 recommender", "phages in common", "T0-T2", "U", "reported", f"{RP} #27", "recommender code and test-strain genomes not published"),
    ("SI Text 4", "without morphotype", "T0-T2", "D", "reported", "reports/10_repository_analyses.md", "no value given"),
    ("SI Text 6", "K1 strains in the Picard", "T0-T1", "P", "reported", f"{RP} #4", "42 K1 strains in the repository, 45 in the text; regenerated types give 45 across the whole collection and 41 among the 369 Picard strains (F)"),
    ("SI Text 6", "K1 rule", "T0-T2", "X", "to test", f"{RP} #4", "the K1 count (42 against 45), the rule definition and the strain set are the candidates"),
    ("SI Text 7", "phages with 0 / 1 / 2 RBP", "T0-T1", "P", "reported", f"{RP} #3", "RBP_list.csv gives 7 phages no RBP; the text says every phage has one or two"),
    ("SI Text 7", "LMM beta of the RBP count", "T0-T1", "P", "reported", f"{RP} #3", "sign differs when the 7 phages count as 0; +16 against +28 with the 89 phages that have an RBP, a residual that stays unexplained (X)"),
    ("SI Text 8", "adsorption factors: AIC", "T0-T2", "D", "reported", "reports/10_repository_analyses.md", "AIC depends on the feature encoding and distance, which the supplement does not state"),
    ("SI Text 8", "defense subfamilies", "T0-T2", "U", "reported", "reports/10_repository_analyses.md", "the authors' gene clustering into defense subfamilies is not published"),
]
LEVEL = lambda fig: "overall" if fig.startswith("Fig 5") else "other"


def register():
    d = pd.read_csv(f"{R}/comparison/paper_vs_replication.csv")
    d = d[~d.status.str.startswith(("match", "consistent"))].reset_index(drop=True)
    rows = []
    for i, r in d.iterrows():
        hit = [c for c in CODES if r.figure.startswith(c[0]) and c[1] in r.metric]
        assert len(hit) >= 1, f"no code for {r.figure} | {r.metric}"
        c = max(hit, key=lambda c: len(c[0]) + len(c[1]))
        rows.append(dict(id=f"G{i + 1:02d}", figure=r.figure, metric=r.metric, level=LEVEL(r.figure), tier_pair=c[2],
                         paper=r.paper, shipped_data=r.published_data_result, replication_shipped=r.replication_shipped_inputs,
                         replication_regenerated=r.replication_regenerated_inputs, status=r.status,
                         cause_code=c[3], basis=c[4], evidence=c[5], note=c[6]))
    g = pd.DataFrame(rows)
    assert len(g) == len(d), (len(g), len(d))
    g.to_csv(f"{OUT}/gap_register.csv", index=False)
    print(f"gap_register.csv: {len(g)} rows\n{g.cause_code.value_counts().to_string()}\nbasis: {g.basis.value_counts().to_dict()}")
    return g


# ------------------------------------------------------------------ T67: paper against repository (comparison a), from the 28 numbered items of report 11
REPO_CODES = {1: "P", 2: "P", 3: "P", 4: "P", 5: "P", 6: "P", 7: "P", 8: "P", 9: "U", 10: "P", 11: "U", 12: "P", 13: "P", 14: "P", 15: "C", 16: "D", 17: "P",
              18: "D", 19: "P", 20: "D", 21: "P", 22: "P", 23: "P", 24: "D", 25: "P", 26: "P", 27: "U", 28: "P"}


def paper_repo():
    rows = []
    for line in open(f"{ROOT}/publish/phage-host-feature-necessity/reports/11_repository_vs_paper.md", encoding="utf-8"):
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if line.startswith("|") and cells[0].replace("★", "").strip().isdigit():
            n = int(cells[0].replace("★", "").strip())
            rows.append(dict(item=n, changes_a_reported_number="★" in cells[0], topic=cells[1], repository=cells[2], paper=cells[3],
                             effect=cells[4] if len(cells) > 4 else "", cause_code=REPO_CODES[n], evidence=f"reports/11_repository_vs_paper.md #{n}"))
    d = pd.DataFrame(rows)
    assert len(d) == 28 and d.item.tolist() == list(range(1, 29)), d.item.tolist()
    d.to_csv(f"{OUT}/paper_vs_repo.csv", index=False)
    print(f"paper_vs_repo.csv: {len(d)} items, {int(d.changes_a_reported_number.sum())} change a reported number\n{d.cause_code.value_counts().to_string()}")
    return d


# ------------------------------------------------------------------ T68-T71 helper: which fold held out each pair (rebuilt with the same GroupKFold call as predict.py)
def folds():
    """one row per (phage, strain): the fold that held it out under the feature set chosen by FEATURE_SET; run once per set"""
    from sklearn.model_selection import GroupKFold
    rows = []
    for ph in P.Y.columns:
        X, y, g = P.design(ph)
        for k, (_, te) in enumerate(GroupKFold(10).split(X, y, g)):
            rows += [(ph, b, k, grp) for b, grp in zip(y.index[te], g.iloc[te])]
    f = pd.DataFrame(rows, columns=["phage", "bacteria", "fold", "group"])
    f.to_csv(f"{OUT}/folds_{P.SET.replace('+', '_')}.csv", index=False)
    print(P.SET, len(f), "pairs,", f.group.nunique(), "groups")


def authors_folds():
    """the authors' fold index per pair (their prediction files), so predict.py can reuse their partition: FOLDS=results/audit/authors_folds.csv"""
    _, ship = P.published()
    ship[["phage", "bacteria", "fold"]].astype({"fold": int}).to_csv(f"{OUT}/authors_folds.csv", index=False)
    print("authors_folds.csv:", len(ship), "pairs")


# ------------------------------------------------------------------ T65: tier alignment
def align(t):
    key = {k: set(zip(v.phage, v.bacteria)) for k, v in t.items()}
    lab = {k: v.set_index(["phage", "bacteria"]).y for k, v in t.items()}
    rows = []
    for k, v in t.items():
        common = key[k] & key["T2"]
        flips = int((lab[k].loc[list(common)] != lab["T2"].loc[list(common)]).sum())
        rows.append(dict(tier=k, strains=v.bacteria.nunique(), phages=v.phage.nunique(), pairs=len(v), lytic_fraction=round(v.y.mean(), 4),
                         pairs_not_in_T2=len(key[k] - key["T2"]), T2_pairs_missing=len(key["T2"] - key[k]), label_flips_vs_T2=flips))
    a = pd.DataFrame(rows)
    a.to_csv(f"{OUT}/tier_alignment.csv", index=False)
    M = P.M
    diff = []
    for k in ("T1", "T3", "T4"):
        for only, other, side in ((key[k] - key["T2"], k, "only in " + k), (key["T2"] - key[k], k, "missing from " + k)):
            for ph, s in sorted(only):
                diff.append(dict(tier=k, side=side, phage=ph, bacteria=s, shipped_matrix_value=M.loc[s, ph] if s in M.index and ph in M.columns else np.nan))
    s = pd.DataFrame(diff)
    s["likely_reason"] = s.side.map({
        "missing from T1": "NaN in the authors' original N/P/F/U matrix (predict_all_phages.py:41 turns F into NaN); that file is not shipped, so unverified",
        "only in T1": "NaN in the public MLC matrix but labelled in the authors' original matrix; unverified",
        "only in T4": "our matrix from the raw plaque scores fills a cell the shipped matrix leaves empty"})
    s.to_csv(f"{OUT}/pair_set_differences.csv", index=False)
    print(a.to_string(index=False))
    print(s.groupby(["tier", "side"]).agg(pairs=("phage", "size"), phages=("phage", "nunique"), strains=("bacteria", "nunique"), shipped_value_nan=("shipped_matrix_value", lambda x: int(x.isna().sum()))).to_string())
    return a


# ------------------------------------------------------------------ T66: overall metrics with bootstrap intervals over strains
def wauc(y, p, w):
    return roc_auc_score(y, p, sample_weight=w) if (w[y == 1].sum() > 0 and w[y == 0].sum() > 0) else np.nan


def waupr(y, p, w):
    return average_precision_score(y, p, sample_weight=w) if w[y == 1].sum() > 0 else np.nan


def metrics(d, w):
    """d: pair table with a strain weight column; returns the five overall metrics for one set of strain weights"""
    pw = d.bacteria.map(w).to_numpy(float)
    keep = pw > 0
    y, p, pw, ph = d.y.to_numpy()[keep], d.p.to_numpy()[keep], pw[keep], d.phage.to_numpy()[keep]
    per = [(wauc(y[ph == q], p[ph == q], pw[ph == q]), waupr(y[ph == q], p[ph == q], pw[ph == q])) for q in np.unique(ph)]
    return dict(pooled_auroc=wauc(y, p, pw), pooled_aupr=waupr(y, p, pw), mean_phage_auroc=np.nanmean([a for a, _ in per]), mean_phage_aupr=np.nanmean([b for _, b in per]))


def per_strain_auroc(d):
    return d.groupby("bacteria").apply(lambda x: roc_auc_score(x.y, x.p) if x.y.nunique() > 1 else np.nan, include_groups=False)


def overall(t):
    key = set.intersection(*[set(zip(v.phage, v.bacteria)) for v in t.values()])
    c = {k: v[[(a, b) in key for a, b in zip(v.phage, v.bacteria)]].reset_index(drop=True) for k, v in t.items()}   # pairs present in every tier
    strains = np.array(sorted(set(c["T1"].bacteria)))
    rng = np.random.default_rng(SEED)
    W = [pd.Series(np.bincount(rng.integers(0, len(strains), len(strains)), minlength=len(strains)), index=strains) for _ in range(B)]
    one = pd.Series(1, index=strains)
    ps = {k: per_strain_auroc(v) for k, v in c.items()}
    point, boot = {}, {}
    for k, d in c.items():
        point[k] = metrics(d, one) | dict(mean_strain_auroc=ps[k].mean())
        boot[k] = pd.DataFrame([metrics(d, w) | dict(mean_strain_auroc=np.nansum(ps[k] * w.reindex(ps[k].index)) / w.reindex(ps[k].dropna().index).sum()) for w in W])
    rows = []
    for k in c:
        for m in point[k]:
            lo, hi = np.percentile(boot[k][m], [2.5, 97.5])
            row = dict(tier=k, metric=m, value=point[k][m], ci_low=lo, ci_high=hi)
            for ref in ("T1", "T2"):
                if k != ref:
                    dd = boot[k][m] - boot[ref][m]
                    row[f"gap_to_{ref}"] = point[k][m] - point[ref][m]
                    row[f"gap_to_{ref}_ci"] = f"{np.percentile(dd, 2.5):+.4f} to {np.percentile(dd, 97.5):+.4f}"
            rows.append(row)
    o = pd.DataFrame(rows)
    # T0: the paper's values, read from the comparison table so they carry the same source as everywhere else
    cmp = pd.read_csv(f"{R}/comparison/paper_vs_replication.csv")
    get = lambda fig, met: float(cmp[(cmp.figure == fig) & cmp.metric.str.startswith(met)].iloc[0].paper)
    t0 = {"pooled_auroc": get("Fig 5", "pooled AUROC"), "pooled_aupr": get("Fig 5", "average precision"), "mean_phage_auroc": get("Fig 5a", "mean AUROC per phage"),
          "mean_strain_auroc": get("Fig 5a", "mean AUROC per strain (published text")}
    o["paper_T0"] = o.metric.map(t0)
    o["gap_to_T0"] = o.value - o.paper_T0
    o["n_pairs"] = len(c["T1"])
    o["definition"] = o.metric.map({  # the per-phage means here pool the held-out predictions of the 10 folds; the authors' script averages the 10 per-fold values (results/*/per_phage.csv)
        "pooled_auroc": "AUROC of all held-out pairs", "pooled_aupr": "average precision of all held-out pairs",
        "mean_phage_auroc": "mean over 96 phages of the AUROC of the phage's pooled held-out predictions",
        "mean_phage_aupr": "mean over 96 phages of the average precision of the phage's pooled held-out predictions",
        "mean_strain_auroc": "mean over strains of the AUROC across their phages"})
    o.to_csv(f"{OUT}/overall.csv", index=False)
    pv = o.pivot(index="metric", columns="tier", values="value").join(pd.Series(t0, name="T0 paper"))
    print(f"pairs present in all four tiers: {len(c['T1'])}; strains {len(strains)}; bootstrap {B} replicates over strains\n{pv.round(3).to_string()}")
    print(o[o.tier != "T1"][["tier", "metric", "gap_to_T1", "gap_to_T1_ci", "gap_to_T2", "gap_to_T2_ci"]].round(4).to_string(index=False))
    return o


# ------------------------------------------------------------------ T68-T71 shared pieces
FEATS = ["O-type", "LPS_type", "ST_Warwick", "Klebs_capsule_type", "ABC_serotype"]


def fold_metrics():
    """per phage and fold: AUROC and AUPR of the tier's best model on the held-out fold (as predict.py computes them)"""
    out = {}
    # T1: the authors' per-fold test metrics of the best model by mean test AUPR (same rule as predict.published)
    rows = []
    for f in sorted(glob_(f"{ROOT}/dev/predictions/results/performances/performance_*.csv")):
        d = pd.read_csv(f, sep=";")
        d = d[d.dataset == "test"].dropna(subset=["model"])
        best = d.groupby("model").avg_precision.mean().idxmax()
        rows.append(d[d.model == best][["phage", "fold", "roc_auc", "avg_precision"]].rename(columns={"roc_auc": "auroc", "avg_precision": "aupr"}))
    out["T1"] = pd.concat(rows)
    t = tiers()
    for k, sfx in (("T2", "shipped"), ("T3", "regenerated"), ("T4", "regenerated_matrix")):
        f = pd.read_csv(f"{OUT}/folds_{sfx}.csv")
        d = t[k].merge(f, on=["phage", "bacteria"])
        out[k] = pd.DataFrame([dict(phage=ph, fold=fo, auroc=roc_auc_score(g.y, g.p), aupr=average_precision_score(g.y, g.p))
                               for (ph, fo), g in d.groupby(["phage", "fold"]) if g.y.nunique() > 1])
    return out


def glob_(pattern):
    import glob
    return glob.glob(pattern)


def tier_table():
    """per phage and tier: authors' definition of the metrics (mean over folds that hold both classes), their standard error, best model, pooled AUROC"""
    fm = fold_metrics()
    t = tiers()
    ship_pp, _ = P.published()
    best = {"T1": ship_pp.best_model, **{k: pd.read_csv(f"{R}/{d}/per_phage.csv", index_col=0).best_model for k, d in TIER_DIRS.items()}}
    npos = {k: pd.read_csv(f"{R}/{d}/per_phage.csv", index_col=0).n_pos for k, d in TIER_DIRS.items()}
    rows = {}
    for k in ("T1", "T2", "T3", "T4"):
        g = fm[k].groupby("phage")
        m = pd.DataFrame({"auroc": g.auroc.mean(), "se": g.auroc.std() / np.sqrt(g.auroc.count()), "aupr": g.aupr.mean(), "folds": g.auroc.count(),
                          "best_model": best[k], "pooled_auroc": t[k].groupby("phage").apply(lambda x: roc_auc_score(x.y, x.p) if x.y.nunique() > 1 else np.nan, include_groups=False)})
        if k != "T1":
            m["n_lytic"] = npos[k]
        rows[k] = m
    return rows


# ------------------------------------------------------------------ T68: per-phage table, 96 rows by tier
def per_phage(thr=0.05):
    rows = tier_table()
    d = pd.DataFrame(index=rows["T2"].index)
    d["n_lytic_T2"], d["n_lytic_T4"] = rows["T2"].n_lytic, rows["T4"].n_lytic
    for k, m in rows.items():
        for c in ("auroc", "se", "aupr", "best_model", "pooled_auroc", "folds"):
            d[f"{c}_{k}"] = m[c]
    for a, b in (("T2", "T1"), ("T3", "T2"), ("T4", "T2")):
        d[f"d_{a}_{b}"] = d[f"auroc_{a}"] - d[f"auroc_{b}"]
        d[f"z_{a}_{b}"] = d[f"d_{a}_{b}"] / np.sqrt(d[f"se_{a}"] ** 2 + d[f"se_{b}"] ** 2)
        d[f"flag_{a}_{b}"] = d[f"d_{a}_{b}"].abs() > thr
    key = {k: set(zip(v.phage, v.bacteria)) for k, v in tiers().items()}
    pairs = pd.read_csv(f"{OUT}/pair_set_differences.csv")
    cnt = lambda tier, side: pairs[(pairs.tier == tier) & (pairs.side == side)].groupby("phage").size().reindex(d.index).fillna(0).astype(int)
    d["pairs_T2_not_in_T1"], d["pairs_T1_not_in_T2"], d["pairs_T4_not_in_T2"] = cnt("T1", "missing from T1"), cnt("T1", "only in T1"), cnt("T4", "only in T4")
    t = tiers()
    t2, t4 = t["T2"].set_index(["phage", "bacteria"]).y, t["T4"].set_index(["phage", "bacteria"]).y
    j = pd.concat([t2.rename("y2"), t4.rename("y4")], axis=1).dropna()
    d["label_flips_T4_vs_T2"] = (j.y2 != j.y4).groupby(level=0).sum().reindex(d.index).fillna(0).astype(int)
    d.rename_axis("phage").to_csv(f"{OUT}/per_phage.csv")
    for a, b in (("T2", "T1"), ("T3", "T2"), ("T4", "T2")):
        c = f"d_{a}_{b}"
        print(f"{a} vs {b}: mean |d| {d[c].abs().mean():.3f}, median {d[c].abs().median():.3f}, max {d[c].abs().max():.3f}; flagged (>{thr}) {int(d['flag_' + a + '_' + b].sum())}; "
              f"|z| >= 2: {int((d['z_' + a + '_' + b].abs() >= 2).sum())}; corr of AUROC {d['auroc_' + a].corr(d['auroc_' + b]):.3f}")
    return d


# ------------------------------------------------------------------ T69: explain the large movers
def host_feature_changes():
    """per phage: whether the regenerated features of its isolation host differ from the shipped table (host-sharing features are built from them)"""
    G = f"{ROOT}/data/genomics"
    ship = pd.read_csv(f"{G}/bacteria/picard_collection.csv", sep=";", index_col=0)[FEATS]
    reg = pd.read_csv(f"{ROOT}/replication/T25_features/picard_collection_regenerated.csv", sep=";", index_col=0, dtype={"ST_Warwick": str})
    reg = reg.rename(columns={"LPS_type (cross-validated markers)": "LPS_type", "ABC_model (best system)": "ABC_serotype"})
    reg["ST_Warwick"] = reg["ST_Warwick"].str.replace(r"^~|\*$", "", regex=True)
    r = reg[FEATS].reindex(ship.index)
    s = ship.astype("string").replace({"<NA>": pd.NA})
    s["ST_Warwick"] = ship.ST_Warwick.map(lambda v: pd.NA if pd.isna(v) else str(int(float(v))) if str(v).replace(".", "").isdigit() else str(v))
    rr = r.astype("string")
    rr["ST_Warwick"] = r.ST_Warwick.map(lambda v: pd.NA if pd.isna(v) else str(int(float(v))) if str(v).replace(".", "").isdigit() else str(v))
    changed = ~((s == rr).fillna(False).astype(bool) | (s.isna() & rr.isna()))
    host = pd.read_csv(f"{G}/phages/guelin_collection.csv", sep=";", index_col=0)["Phage_host"]
    return changed, host


def t1_model_spread():
    """per phage: range of the mean test AUROC across the four non-dummy models of the authors' run (how far model choice alone can move an AUROC)"""
    rows = []
    for f in glob_(f"{ROOT}/dev/predictions/results/performances/performance_*.csv"):
        d = pd.read_csv(f, sep=";")
        d = d[(d.dataset == "test") & ~d.model.fillna("").str.startswith("Dummy")].dropna(subset=["model"])
        m = d.groupby("model").roc_auc.mean()
        rows.append((d.phage.iloc[0], m.max() - m.min()))
    return pd.Series(dict(rows))


def causes(thr=0.05):
    import re
    d = pd.read_csv(f"{OUT}/per_phage.csv", index_col=0)
    changed, host = host_feature_changes()
    spread = t1_model_spread()
    norm = lambda m: re.sub(r"_weight=[0-9.]+", "", m)
    t = tiers()
    t2, t4 = t["T2"].set_index(["phage", "bacteria"]), t["T4"].set_index(["phage", "bacteria"])

    def label_effect(ph):
        """pooled AUROC of T4's predictions of one phage scored against T4's own labels and against T2's labels, on the pairs both runs have"""
        j = t4.loc[ph].join(t2.loc[ph], how="inner", lsuffix="4", rsuffix="2")
        return roc_auc_score(j.y4, j.p4), roc_auc_score(j.y2, j.p4)
    rows = []
    for a, b in (("T2", "T1"), ("T3", "T2"), ("T4", "T2")):
        for ph, r in d[d[f"flag_{a}_{b}"]].iterrows():
            h = host.get(ph)
            hc = [c for c in FEATS if h in changed.index and changed.loc[h, c]] if a in ("T3", "T4") and not ph.startswith("LF110") else []
            diff_pairs = int(r.pairs_T2_not_in_T1 + r.pairs_T1_not_in_T2) if a == "T2" else int(r.label_flips_T4_vs_T2 + r.pairs_T4_not_in_T2) if a == "T4" else 0
            models_differ = norm(r[f"best_model_{a}"]) != norm(r[f"best_model_{b}"])
            m_note = ""
            if a == "T4" and diff_pairs:
                own, on2 = label_effect(ph)
                m_ruled_out = abs(own - on2) < 0.25 * abs(dd)
                m_note = f"; label changes {'ruled out' if m_ruled_out else 'not ruled out'} (T4 predictions score {own:.3f} on T4 labels and {on2:.3f} on T2 labels)"
            else:
                m_ruled_out = False
            cand = (["M"] if diff_pairs and not m_ruled_out else []) + (["C"] if models_differ else []) + (["F"] if hc else []) + ["U"]   # U: the CV partition differs for every phage (folds_*.csv, adjusted Rand index about 0.02); its share is untested
            z, dd, sp, fo = r[f"z_{a}_{b}"], r[f"d_{a}_{b}"], spread.get(ph, np.nan), min(r[f"folds_{a}"], r[f"folds_{b}"])
            if fo < 3:
                primary, basis = "S", f"AUROC rests on {int(fo)} fold(s) with both classes (n lytic {int(r.n_lytic_T2)}), too few for a stable value"
            elif abs(z) < 2:
                primary, basis = "S", f"within fold-sampling noise (|z| = {abs(z):.1f}); partition and model randomness differ between the runs"
            elif models_differ and sp >= abs(dd):
                primary, basis = "C", f"best model differs; the authors' own four models differ by up to {sp:.2f} AUROC, enough to produce a shift of {abs(dd):.2f} (sufficient, not shown)"
            else:
                primary, basis = "X", f"no listed factor is sufficient on the evidence in hand (the authors' four models differ by only {sp:.2f} AUROC){m_note}"
            rows.append(dict(comparison=f"{a} vs {b}", phage=ph, n_lytic=int(r.n_lytic_T2), d_auroc=dd, z=z, folds_with_both_classes=int(fo),
                             best_model_a=r[f"best_model_{a}"], best_model_b=r[f"best_model_{b}"], pairs_differing=diff_pairs,
                             host_features_changed=";".join(hc), candidates="+".join(cand), cause_code=primary, basis=basis))
    c = pd.DataFrame(rows)
    c.to_csv(f"{OUT}/per_phage_causes.csv", index=False)
    print(c.groupby(["comparison", "cause_code"]).size().unstack(fill_value=0).to_string())
    return c


# ------------------------------------------------------------------ T70: pair level, authors' predictions against T2
def pairs_t1_t2(thr=0.5):
    from scipy.stats import spearmanr
    _, ship = P.published()
    a = ship.rename(columns={"p": "p1", "fold": "fold1"})   # the authors' prediction files carry their fold index
    b = pd.read_csv(f"{R}/repro/predictions.csv").merge(pd.read_csv(f"{OUT}/folds_shipped.csv"), on=["phage", "bacteria"]).rename(columns={"p": "p2", "fold": "fold2", "model": "model2"})
    j = a[["phage", "bacteria", "y", "p1", "fold1", "model"]].rename(columns={"model": "model1"}).merge(b[["phage", "bacteria", "p2", "fold2", "model2", "group"]], on=["phage", "bacteria"])
    j["call1"], j["call2"] = j.p1 >= thr, j.p2 >= thr
    j["flip"] = j.call1 != j.call2
    norm = lambda m: m.str.replace(r"_weight=[0-9.]+", "", regex=True)
    j["same_family"] = j.model1.str.split("_").str[0] == j.model2.str.split("_").str[0]
    j["same_model"] = norm(j.model1) == norm(j.model2)
    rho = j.groupby("phage").apply(lambda x: spearmanr(x.p1, x.p2)[0], include_groups=False)
    # does each of our CV groups sit inside one authors' fold? (the earlier 0.3 % statement)
    gg = j.groupby(["phage", "group"]).fold1.nunique()
    from sklearn.metrics import adjusted_rand_score
    ari = j.groupby("phage").apply(lambda x: adjusted_rand_score(x.fold1, x.fold2), include_groups=False)   # fold indices are arbitrary labels, so compare the partitions
    summ = pd.DataFrame([dict(pairs=len(j), spearman_all=spearmanr(j.p1, j.p2)[0], spearman_phage_median=rho.median(), spearman_phage_min=rho.min(),
                              call_agreement=1 - j.flip.mean(), flipped_pairs=int(j.flip.sum()), flips_to_lytic=int((~j.call1 & j.call2).sum()), flips_to_nonlytic=int((j.call1 & ~j.call2).sum()),
                              same_best_model_family_share=j.groupby("phage").same_family.first().mean(),
                              same_best_model_share=j.groupby("phage").same_model.first().mean(),
                              flip_rate_same_model=j[j.same_model].flip.mean(), flip_rate_same_family_other_model=j[j.same_family & ~j.same_model].flip.mean(),
                              flip_rate_different_family=j[~j.same_family].flip.mean(),
                              groups_split_across_authors_folds=float((gg > 1).mean()), fold_partition_adjusted_rand_mean=float(ari.mean()), fold_partition_adjusted_rand_max=float(ari.max()))])
    summ.T.rename(columns={0: "value"}).to_csv(f"{OUT}/pairs_t1_t2_summary.csv")
    j[j.flip].drop(columns=["call1", "call2"]).to_csv(f"{OUT}/pairs_t1_t2.csv", index=False)
    byphage = j.groupby("phage").agg(pairs=("flip", "size"), flips=("flip", "sum"), mean_abs_dp=("p1", lambda s: np.nan)).drop(columns="mean_abs_dp")
    byphage["spearman"] = rho
    byphage.to_csv(f"{OUT}/pairs_t1_t2_by_phage.csv")
    print(summ.T.to_string(header=False))
    return j


# ------------------------------------------------------------------ T71: pair level, T2 against T3 and T4, with the changes that could explain each flip
def pairs_regeneration(thr=0.5):
    import re
    from scipy.stats import spearmanr
    t = tiers()
    changed, host = host_feature_changes()
    strain_changed = changed.any(axis=1)
    pp = pd.read_csv(f"{OUT}/per_phage.csv", index_col=0)
    norm = lambda m: re.sub(r"_weight=[0-9.]+", "", m)
    res, summ = [], []
    for k in ("T3", "T4"):
        j = t["T2"].merge(t[k], on=["phage", "bacteria"], suffixes=("2", "x"))
        j["flip"] = (j.p2 >= thr) != (j.px >= thr)
        j["label_changed"] = j.y2 != j.yx
        j["strain_features_changed"] = j.bacteria.map(strain_changed).fillna(False).astype(bool)
        j["host_features_changed"] = j.phage.map(lambda ph: bool(changed.loc[host.get(ph)].any()) if host.get(ph) in changed.index and not ph.startswith("LF110") else False)
        changed_model = {ph: norm(pp.loc[ph, "best_model_T2"]) != norm(pp.loc[ph, f"best_model_{k}"]) for ph in pp.index}
        j["best_model_changed"] = j.phage.map(changed_model).astype(bool)
        j["none_of_these"] = ~(j.label_changed | j.strain_features_changed | j.host_features_changed | j.best_model_changed)
        for name in ("label_changed", "strain_features_changed", "host_features_changed", "best_model_changed"):
            for v in (False, True):
                s = j[j[name] == v]
                if len(s): summ.append(dict(tier=k, factor=name, value=v, pairs=len(s), flips=int(s.flip.sum()), flip_rate=s.flip.mean()))
        for (sc, mc), s in j.groupby(["strain_features_changed", "best_model_changed"]):
            summ.append(dict(tier=k, factor="strain_features_changed x best_model_changed", value=f"{sc} x {mc}", pairs=len(s), flips=int(s.flip.sum()), flip_rate=s.flip.mean()))
        res.append(dict(tier=k, pairs=len(j), flips=int(j.flip.sum()), flip_rate=j.flip.mean(), spearman=spearmanr(j.p2, j.px)[0],
                        flips_with_label_change=int((j.flip & j.label_changed).sum()), flips_with_strain_feature_change=int((j.flip & j.strain_features_changed).sum()),
                        flips_with_host_feature_change=int((j.flip & j.host_features_changed).sum()), flips_with_best_model_change=int((j.flip & j.best_model_changed).sum()),
                        flips_with_none=int((j.flip & j.none_of_these).sum())))
        j[j.flip].to_csv(f"{OUT}/pairs_regeneration_{k}.csv", index=False)
    pd.DataFrame(res).to_csv(f"{OUT}/pairs_regeneration_summary.csv", index=False)
    pd.DataFrame(summ).to_csv(f"{OUT}/pairs_regeneration_by_factor.csv", index=False)
    print(pd.DataFrame(res).T.to_string(header=False)); print(pd.DataFrame(summ).round(4).to_string(index=False))


# ------------------------------------------------------------------ T72a, T72, T73: attribution runs (scripts/run_audit_runs.sh)
def load_run(tag):
    d = f"{R}/audit_runs/{tag}"
    return pd.read_csv(f"{d}/per_phage.csv", index_col=0), pd.read_csv(f"{d}/predictions.csv")


def pooled(pr):
    return roc_auc_score(pr.y, pr.p)


def fixed_folds():
    """T72a: our models on shipped inputs with the authors' folds, against the authors' own per-phage numbers and pair probabilities"""
    import re
    a, apr = load_run("fixed_default")
    ship_pp, ship_pr = P.published()
    norm = lambda m: re.sub(r"_weight=[0-9.]+", "", m)
    d = pd.DataFrame({"auroc_authors": ship_pp.auroc, "auroc_fixed_folds": a.auroc, "aupr_authors": ship_pp.aupr, "aupr_fixed_folds": a.aupr,
                      "best_authors": ship_pp.best_model.map(norm), "best_fixed_folds": a.best_model.map(norm), "n_lytic": a.n_pos})
    d["d_auroc"], d["d_aupr"] = d.auroc_fixed_folds - d.auroc_authors, d.aupr_fixed_folds - d.aupr_authors
    d["same_best_model"] = d.best_authors == d.best_fixed_folds
    t2 = pd.read_csv(f"{OUT}/per_phage.csv", index_col=0)
    d["d_auroc_our_folds"] = t2.auroc_T2 - ship_pp.auroc   # the earlier gap, with our own folds
    d.to_csv(f"{OUT}/fixed_folds_vs_authors.csv")
    j = ship_pr[["phage", "bacteria", "p", "y", "model"]].merge(apr[["phage", "bacteria", "p", "y"]], on=["phage", "bacteria"], suffixes=("1", "2"))
    j["same_best"] = j.phage.map(d.same_best_model)
    j["dp"] = (j.p1 - j.p2).abs()
    flip = (j.p1 >= 0.5) != (j.p2 >= 0.5)
    summ = dict(phages=len(d), identical_auroc_1e6=int((d.d_auroc.abs() < 1e-6).sum()), identical_aupr_1e6=int((d.d_aupr.abs() < 1e-6).sum()),
                mean_abs_d_auroc=d.d_auroc.abs().mean(), max_abs_d_auroc=d.d_auroc.abs().max(), mean_abs_d_auroc_our_folds=d.d_auroc_our_folds.abs().mean(),
                same_best_model=int(d.same_best_model.sum()), phages_d_above_0_05=int((d.d_auroc.abs() > 0.05).sum()),
                mean_auroc_authors=d.auroc_authors.mean(), mean_auroc_fixed_folds=d.auroc_fixed_folds.mean(),
                pooled_auroc_authors=pooled(ship_pr), pooled_auroc_fixed_folds=pooled(apr),
                pairs=len(j), pairs_prob_identical_1e6=int((j.dp < 1e-6).sum()), pairs_same_best_model=int(j.same_best.sum()),
                prob_identical_share_among_same_model=float((j[j.same_best].dp < 1e-6).mean()), flips_at_0_5=int(flip.sum()),
                flip_rate_same_model=float(flip[j.same_best].mean()), flip_rate_other_model=float(flip[~j.same_best].mean()))
    pd.Series(summ).to_csv(f"{OUT}/fixed_folds_summary.csv", header=["value"])
    print(pd.Series(summ).to_string())
    return d


def seed_spread():
    """T72: spread of the per-phage AUROC across seeds with the folds fixed, and across seeds with a shuffled group partition"""
    fx = {s: load_run(f"fixed_s{s}") for s in range(5)}
    sh = {s: load_run(f"shuffle_s{s}") for s in range(5)}
    A = lambda runs: pd.DataFrame({s: r[0].auroc for s, r in runs.items()})
    Pq = lambda runs: pd.DataFrame({s: r[0].aupr for s, r in runs.items()})
    af, ash = A(fx), A(sh)
    d = pd.DataFrame({"n_lytic": fx[0][0].n_pos, "mean_fixed": af.mean(axis=1), "sd_fixed": af.std(axis=1), "range_fixed": af.max(axis=1) - af.min(axis=1),
                      "mean_shuffled": ash.mean(axis=1), "sd_shuffled": ash.std(axis=1), "range_shuffled": ash.max(axis=1) - ash.min(axis=1)})
    t = pd.read_csv(f"{OUT}/per_phage.csv", index_col=0)
    d["d_T2_T1"] = t.d_T2_T1
    d["z_T2_T1_vs_partition"] = d.d_T2_T1 / (np.sqrt(2) * d.sd_shuffled)
    d.to_csv(f"{OUT}/seed_spread.csv")
    per_run = lambda runs: pd.DataFrame({s: dict(mean_phage_auroc=r[0].auroc.mean(), mean_phage_aupr=r[0].aupr.mean(), pooled_auroc=pooled(r[1])) for s, r in runs.items()}).T
    rf, rs = per_run(fx), per_run(sh)
    summ = pd.DataFrame({"fixed folds, seeds 0-4": rf.agg(["mean", "std", "min", "max"]).stack(), "shuffled partition, seeds 0-4": rs.agg(["mean", "std", "min", "max"]).stack()})
    summ.to_csv(f"{OUT}/seed_spread_summary.csv")
    print(summ.round(4).to_string())
    print(f"per-phage sd of AUROC: fixed folds mean {d.sd_fixed.mean():.4f} (max {d.sd_fixed.max():.4f}); shuffled partition mean {d.sd_shuffled.mean():.4f} (max {d.sd_shuffled.max():.4f})")
    print(f"phages with |d T2 vs T1| > 2 sqrt(2) sd_shuffled: {int((d.z_T2_T1_vs_partition.abs() > 2).sum())} of {len(d)}; "
          f"phages with |d| > 0.05: {int((d.d_T2_T1.abs() > .05).sum())}, of which beyond partition noise: {int(((d.d_T2_T1.abs() > .05) & (d.z_T2_T1_vs_partition.abs() > 2)).sum())}")
    return d


GROUPS_SWAPPED = ["O-type", "ST", "LPS", "Klebsiella", "ABC", "UMAP", "host"]


def swap_groups():
    """T73: effect of each regenerated feature group on its own, with the authors' folds and seed 0 held fixed"""
    base_pp, base_pr = load_run("swap_none")
    sp = pd.read_csv(f"{OUT}/seed_spread.csv", index_col=0)
    noise = 2 * np.sqrt(2) * sp.sd_fixed.clip(lower=1e-3)   # per phage: change that a seed alone can produce (fixed folds); floor 0.003
    rows, per = [], {}
    keyb = base_pr.set_index(["phage", "bacteria"])
    for g in GROUPS_SWAPPED + ["all"]:
        pp, pr = load_run(f"swap_{g}")
        d = (pp.auroc - base_pp.auroc)
        per[g] = d
        k = pr.set_index(["phage", "bacteria"])
        common = k.index.intersection(keyb.index)
        flips = ((k.loc[common, "p"] >= 0.5) != (keyb.loc[common, "p"] >= 0.5)).mean()
        big = base_pp.n_pos >= 10   # phages with fewer than 10 lytic strains rest on a handful of folds
        rows.append(dict(group=g, mean_d_auroc=d.mean(), mean_d_auroc_n10=d[big].mean(), mean_abs_d_auroc=d.abs().mean(), max_abs_d_auroc=d.abs().max(), phages_abs_d_above_0_05=int((d.abs() > .05).sum()),
                         phages_beyond_seed_noise=int((d.abs() > noise).sum()), mean_phage_auroc=pp.auroc.mean(), pooled_auroc=pooled(pr), d_pooled_auroc=pooled(pr) - pooled(base_pr), call_flip_rate=flips))
    s = pd.DataFrame(rows)
    s.to_csv(f"{OUT}/swap_groups.csv", index=False)
    per = pd.DataFrame(per)
    per["main_group"] = per[GROUPS_SWAPPED].abs().idxmax(axis=1)
    per["main_d"] = [per.loc[i, g] for i, g in zip(per.index, per.main_group)]
    per["sum_of_groups"] = per[GROUPS_SWAPPED].sum(axis=1)
    per["noise"] = noise
    per.to_csv(f"{OUT}/swap_per_phage.csv")
    # T73b: marginal effect of a group when the others are already regenerated = (all regenerated) - (all but that group regenerated)
    drops = []
    all_pp, all_pr = load_run("swap_all")
    big = base_pp.n_pos >= 10
    for g in GROUPS_SWAPPED:
        pp, pr = load_run(f"drop_{g}")
        m = all_pp.auroc - pp.auroc
        per[f"marginal_{g}"] = m
        drops.append(dict(group=g, mean_marginal_d_auroc=m.mean(), mean_marginal_d_auroc_n10=m[big].mean(), mean_abs_marginal=m.abs().mean(), phages_abs_marginal_above_0_05=int((m.abs() > .05).sum()),
                          mean_phage_auroc_without_group=pp.auroc.mean(), pooled_auroc_without_group=pooled(pr)))
    dr = pd.DataFrame(drops)
    dr.to_csv(f"{OUT}/swap_marginal.csv", index=False)
    per.to_csv(f"{OUT}/swap_per_phage.csv")
    t = pd.read_csv(f"{OUT}/per_phage.csv", index_col=0)
    ch, host = host_feature_changes()
    n_host = int(sum(bool(ch.loc[h, ["O-type", "ST_Warwick"]].any()) for h in host if h in ch.index))
    pd.Series(dict(swap_none_vs_fixed_s0_max_abs_d_auroc=(base_pp.auroc - load_run("fixed_s0")[0].auroc).abs().max(), swap_all_vs_T3_per_phage_corr=all_pp.auroc.corr(t.auroc_T3),
                   phages_whose_host_o_or_st_changed=n_host)).to_csv(f"{OUT}/swap_controls.csv", header=["value"])
    print(s.round(4).to_string(index=False)); print(dr.round(4).to_string(index=False))
    print(f"sum of the seven single-group effects on mean AUROC {s[s.group != 'all'].mean_d_auroc.sum():+.4f}; all groups together {s[s.group == 'all'].mean_d_auroc.iloc[0]:+.4f}")
    return s, per


def causes_runs():
    """T69 refined with the attribution runs: T2 vs T1 rows get the change that remains with the authors' folds; T3 and T4 rows get what the feature swaps show"""
    c = causes()
    fx = pd.read_csv(f"{OUT}/fixed_folds_vs_authors.csv", index_col=0)
    sw = pd.read_csv(f"{OUT}/swap_per_phage.csv", index_col=0)
    sp = pd.read_csv(f"{OUT}/seed_spread.csv", index_col=0)
    c["z_vs_partition_noise"] = [r.d_auroc / (np.sqrt(2) * sp.loc[r.phage, "sd_shuffled"]) for r in c.itertuples()]   # change against the spread of 5 shuffled CV partitions (T72)
    c["d_auroc_with_authors_folds"] = np.nan
    c["evidence"] = ""
    c["swap_all_d_auroc"], c["swap_main_group"], c["swap_main_d_auroc"] = np.nan, "", np.nan
    for i, r in c.iterrows():
        if r.comparison == "T2 vs T1":
            c.loc[i, "d_auroc_with_authors_folds"] = fx.loc[r.phage, "d_auroc"]
            if abs(fx.loc[r.phage, "d_auroc"]) < 0.05:
                c.loc[i, ["cause_code", "basis", "evidence"]] = ["U", f"CV partition: with the authors' folds the change shrinks from {r.d_auroc:+.3f} to {fx.loc[r.phage, 'd_auroc']:+.3f}", "tested: authors' folds (T72a)"]
            else:
                c.loc[i, "basis"] = f"change stays at {fx.loc[r.phage, 'd_auroc']:+.3f} with the authors' folds; best model {'identical' if fx.loc[r.phage, 'same_best_model'] else 'differs'} (T72a)"
        else:
            row = sw.loc[r.phage]
            c.loc[i, "swap_all_d_auroc"], c.loc[i, "swap_main_group"], c.loc[i, "swap_main_d_auroc"] = row["all"], row.main_group, row.main_d
            d_, a_ = r.d_auroc, row["all"]
            zf, zp = r.z, c.loc[i, "z_vs_partition_noise"]
            seed_noise = 2 * np.sqrt(2) * max(sp.loc[r.phage, "sd_fixed"], 1e-3)
            half = np.sign(a_) == np.sign(d_) and abs(a_) >= 0.5 * abs(d_)
            if r.cause_code in ("S", "X") and half and abs(a_) >= 0.05:
                c.loc[i, ["cause_code", "basis", "evidence"]] = ["F", f"regenerated features reproduce the change with the folds fixed ({a_:+.3f}); largest single group {row.main_group} ({row.main_d:+.3f})", "tested: feature swap (T73)"]
            elif r.basis.startswith("AUROC rests on"):
                c.loc[i, "evidence"] = "weak: fewer than 3 folds hold both classes"
            elif half and abs(a_) > seed_noise and (abs(zf) >= 2 or abs(zp) >= 2):
                c.loc[i, ["cause_code", "basis", "evidence"]] = ["F", f"regenerated features reproduce at least half of the change with the folds fixed ({a_:+.3f}, above seed noise {seed_noise:.3f}); the change is also beyond the random-partition spread (|z| = {abs(zp):.1f}); largest single group {row.main_group} ({row.main_d:+.3f})", "tested: feature swap, partial (T73)"]
            elif abs(zf) < 2 and abs(zp) < 2:
                c.loc[i, ["cause_code", "basis", "evidence"]] = ["S", f"inside both noise estimates: fold standard error (|z| = {abs(zf):.1f}) and the spread of random CV partitions (|z| = {abs(zp):.1f}); the regenerated features alone change it by {a_:+.3f} with the folds fixed", "noise only: no cause shown"]
            else:
                c.loc[i, ["cause_code", "basis", "evidence"]] = ["X", f"beyond the noise estimates (fold |z| = {abs(zf):.1f}, partition |z| = {abs(zp):.1f}); the regenerated features change it by {a_:+.3f} with the folds fixed, less than half of the change", "none: no tested factor is sufficient"]
    c.to_csv(f"{OUT}/per_phage_causes.csv", index=False)
    print(c.groupby(["comparison", "cause_code"]).size().unstack(fill_value=0).to_string())
    return c


# ------------------------------------------------------------------ T74: tests that settle open register rows, then the finished register
def ed7_definitions():
    """ED7: coverage of the most-covering k-set under candidate definitions of 'productive lysis', against the Source Data curve (counts of 403 strains)"""
    import matrix_stats as MS
    sd = pd.read_excel(f"{ROOT}/docs/paper/supplement/MOESM12.xlsx", sheet_name=0).iloc[1:, :2].astype(float)
    target = {int(k): round(v / 100 * 403) for k, v in zip(sd.iloc[:, 0], sd.iloc[:, 1])}
    rows = [dict(family="Source Data", definition="paper's curve", strains=403, k1=target[1], k2=target[2], k3=target[3], k4=target[4], plateau=target[max(target)])]
    M3 = pd.read_csv(f"{ROOT}/data/derived/interaction_matrix_403.csv", sep=";", index_col=0)
    r = pd.read_csv(f"{ROOT}/data/interactions/raw/raw_interactions.csv", sep=";")
    r = r[r.score != "n"].assign(pos=lambda d: (d.score == "1").astype(int))
    cnt = r.groupby(["bacteria", "phage", "log_dilution", "replicate"]).pos.max().reset_index().groupby(["bacteria", "phage", "log_dilution"]).pos.sum().unstack("log_dilution").fillna(0)
    defs = [("MLC score", f"MLC >= {t}", (M3 >= t).fillna(False)) for t in (1, 2, 3, 4)]
    for dil in (0, -1, -2, -4):
        for rep in (1, 2, 3):
            cs = [c for c in cnt.columns if c <= dil] if dil < 0 else [0]
            if cs: defs.append(("raw plaque scores", f"positive in >= {rep} of 3 replicates at dilution {dil} or more dilute", (cnt[cs].max(axis=1) >= rep).unstack("phage").fillna(False)))
    for fam, name, L in defs:
        cov, _ = MS.coverage_curve(L, kmax=4)
        rows.append(dict(family=fam, definition=name, strains=len(L), k1=round(cov[1] / 100 * len(L)), k2=round(cov[2] / 100 * len(L)), k3=round(cov[3] / 100 * len(L)),
                         k4=round(cov[4] / 100 * len(L)), plateau=int(L.any(axis=1).sum())))
    d = pd.DataFrame(rows)
    t = d.iloc[0]
    d["matches_source_data"] = pd.Series([(r_.strains == 403) and all(abs(r_[c] - t[c]) <= 2 for c in ("k1", "k2", "k3", "k4", "plateau")) for _, r_ in d.iterrows()], dtype=object)
    d.loc[0, "matches_source_data"] = "-"
    d.to_csv(f"{OUT}/ed7_definitions.csv", index=False)
    print(d.to_string(index=False))
    return d


def permanova_checks():
    """Fig 3a: the two candidate causes named in the register (LF110, genus labels), against the paper's adjusted R2 of 0.53 (isolation strain) and 0.65 (genus + isolation strain)"""
    import matrix_stats as MS
    G = pd.read_csv(f"{ROOT}/data/genomics/phages/guelin_collection.csv", sep=";", index_col=0)
    mats = {"402 strains (public matrix)": pd.read_csv(f"{ROOT}/data/interactions/interaction_matrix.csv", sep=";", index_col=0),
            "403 strains (with LF110)": pd.read_csv(f"{ROOT}/data/derived/interaction_matrix_403.csv", sep=";", index_col=0)}
    rows = []
    for name, Mx in mats.items():
        for gcol in ("Genus", "Old_Genus"):
            meta = pd.DataFrame({"Genus": G.loc[Mx.columns, gcol].astype(str), "Host": G.loc[Mx.columns, "Phage_host"].astype(str)})
            res = MS.phage_models(MS.phage_distance(Mx), meta, perms=199)
            rows.append(dict(matrix=name, genus_labels=gcol, adjR2_genus=res["genus"][1], adjR2_isolation_strain=res["isolation strain"][1], adjR2_both=res["genus x isolation strain"][1]))
    d = pd.DataFrame(rows)
    d.loc[len(d)] = dict(matrix="paper", genus_labels="-", adjR2_genus=np.nan, adjR2_isolation_strain=0.53, adjR2_both=0.65)
    d.to_csv(f"{OUT}/permanova_checks.csv", index=False)
    print(d.round(3).to_string(index=False))
    return d


# (figure prefix, metric substring) -> (what was tested, what would settle it); rows not listed get the generic text of their cause code
SETTLE = {"U": "the unpublished file or code", "D": "the authors stating the definition or unit", "P": "the authors correcting the paper, the supplement or the repository",
          "I": "nothing; our deliberate choice", "OK": "nothing; the numbers agree", "F": "nothing; a measured effect of regeneration", "S": "nothing; within resampling noise", "C": "the authors' model-selection seed and library versions", "X": "see the row"}
REFINE = [
    ("Fig 4d", "pMCMC", "not testable here: one chain takes about 4.5 hours and the authors' seed and chain settings are not stated; beta agrees, so only the posterior tail differs",
     "the authors' chain seed, burn-in and the exact defense counts behind panel d"),
    ("ED7", "coverage", "tested 4 MLC thresholds on the 403-strain matrix and 12 raw-plaque definitions (replicates by dilution); none reproduces the Source Data counts 172, 232, 270, 297 and plateau 360 (`ed7_definitions.csv`); the raw scores cover only 369 strains. An independent preprint reports a promiscuity baseline of about 61 % for one phage, near our 65 % and far from 42.7 % (`external_check.csv`)",
     "the authors' definition of productive lysis and the matrix they used for ED7"),
    ("Fig 3a", "genus + isolation strain", "tested adding LF110 and the old genus labels: adjusted R2 stays 0.686 and 0.682, so neither explains 0.65 (`permanova_checks.csv`); isolation strain alone 0.554 against 0.53",
     "the authors' PERMANOVA call (adonis2 formula, which host and genus columns, number of phages)"),
    ("Fig 4b", "significant in both", "not testable here: needs the minimum-observation filter and the outer-membrane clustering the repository does not contain, and another 4.5 hours of sampling",
     "the authors' filtering code and clustering of defense subfamilies"),
    ("SI Text 6", "K1 rule", "tested with the regenerated K1 calls (45 strains in the whole collection): precision 89 %, recall 65 % (reports/10), so the K1 count does not explain the gap; the strain set and the definition of infection behind 73 % and 72 % are not stated",
     "the list of the 45 K1 strains, the strain set and the definition of infection used for precision and recall"),
]


def numbers(x):
    import re
    m = re.findall(r"-?\d+\.?\d*(?:e-?\d+)?", str(x).replace(",", ""))
    return float(m[0]) if m else np.nan


def finish_register():
    g = pd.read_csv(f"{OUT}/gap_register.csv")
    g["tested"], g["would_settle"], g["size"] = "", "", ""
    for i, r in g.iterrows():
        ref = [x for x in REFINE if r.figure.startswith(x[0]) and x[1] in r.metric]
        g.loc[i, "tested"] = ref[0][2] if ref else ("evidence in the cited report" if r.basis == "reported" else "")
        g.loc[i, "would_settle"] = ref[0][3] if ref else SETTLE[r.cause_code]
        a, b = numbers(r.paper), numbers(r.replication_shipped)
        if not np.isnan(a) and not np.isnan(b):
            g.loc[i, "size"] = f"ours {b:g} against paper {a:g} (difference {b - a:+.3g})"
    g["basis"] = np.where(g.cause_code.isin(["X"]) | g.tested.str.startswith("not testable"), g.basis.where(g.basis != "to test", "residual"), g.basis.replace({"to test": "tested"}))
    assert (g.cause_code.notna() & (g.cause_code != "")).all() and (g.evidence.notna()).all()
    g.to_csv(f"{OUT}/gap_register.csv", index=False)
    print(g.cause_code.value_counts().to_dict(), g.basis.value_counts().to_dict())
    return g


def prediction_register():
    """gaps in the prediction results themselves (overall, per phage, per pair), each with one cause code and its evidence file"""
    o = pd.read_csv(f"{OUT}/overall.csv")
    c = pd.read_csv(f"{OUT}/per_phage_causes.csv")
    ps = pd.read_csv(f"{OUT}/pairs_t1_t2_summary.csv", index_col=0)["value"]
    fx = pd.read_csv(f"{OUT}/fixed_folds_summary.csv", index_col=0)["value"]
    sw = pd.read_csv(f"{OUT}/swap_groups.csv")
    ra = pd.read_csv(f"{OUT}/pairs_regeneration_summary.csv")
    pr = pd.read_csv(f"{OUT}/pair_set_differences.csv")
    ov = lambda tier, m: o[(o.tier == tier) & (o.metric == m)].iloc[0]
    rows = []
    add = lambda level, comp, metric, a, b, code, ev, size, settle: rows.append(dict(level=level, comparison=comp, metric=metric, value_a=a, value_b=b, cause_code=code, evidence=ev, size=size, would_settle=settle))
    for m in ("pooled_auroc", "mean_phage_auroc"):
        r = ov("T2", m)
        add("overall", "T2 vs T1", m, ov("T1", m).value, r.value, "OK", "results/audit/overall.csv", f"gap {r.gap_to_T1:+.4f}, 95 % interval {r.gap_to_T1_ci}; inside the interval", "nothing")
    for tier in ("T3", "T4"):
        for m in ("mean_phage_auroc", "mean_phage_aupr"):
            r = ov(tier, m)
            extra = ""
            if m == "mean_phage_auroc":
                extra = f"; with the folds fixed the regenerated features move the per-fold mean AUROC by {sw.loc[sw.group == 'all'].mean_d_auroc.iloc[0]:+.4f} (the same predictions pooled differently)"
            if tier == "T4":
                extra += "; T4 also uses the regenerated matrix (72 labels change, 135 pairs added), a share that was not separated"
            add("overall", f"{tier} vs T2", m, ov("T2", m).value, r.value, "F", "results/audit/overall.csv; swap_groups.csv", f"gap {r.gap_to_T2:+.4f}, 95 % interval {r.gap_to_T2_ci}{extra}",
                "nothing; a measured effect of regeneration (O-type and Klebsiella capsule type carry it)")
    r = ov("T2", "mean_strain_auroc")
    add("overall", "T0 vs T2", "mean_strain_auroc", r.paper_T0, r.value, "P", "reports/11_repository_vs_paper.md #17", f"paper text 0.85, its own Source Data 0.872, authors' predictions 0.881, ours {r.value:.3f}", "the authors correcting the text")
    n = c.groupby(["comparison", "cause_code"]).size()
    for comp in ("T2 vs T1", "T3 vs T2", "T4 vs T2"):
        for code in ("U", "F", "S", "X"):
            if (comp, code) in n.index:
                add("phage", comp, "phages moving by more than 0.05 AUROC", int(c[c.comparison == comp].shape[0]), int(n[(comp, code)]), code, "results/audit/per_phage_causes.csv",
                    f"{int(n[(comp, code)])} of {int(c[c.comparison == comp].shape[0])} flagged phages", SETTLE[code] if code != "U" else "the authors' CV groups file")
    add("pair", "T2 vs T1", "pairs in one set only", int(pr[pr.tier == "T1"].shape[0]), 0, "U", "results/audit/tier_alignment.csv; pair_set_differences.csv",
        "58 pairs T2 has and T1 lacks (16 strains, NILS81 38), 1 pair T1 has and T2 lacks; labels agree on all 38,377 shared pairs; cause unverified (the authors' script turns F into NaN in a matrix file that is not shipped)", "the authors' original N/P/F/U matrix")
    add("pair", "T2 vs T1", "calls flipped at 0.5", int(ps["flipped_pairs"]), 0, "U", "results/audit/pairs_t1_t2_summary.csv; fixed_folds_summary.csv",
        f"{int(ps['flipped_pairs'])} of {int(ps['pairs'])} pairs ({1 - ps['call_agreement']:.1%}) with our folds, {int(fx['flips_at_0_5'])} ({fx['flips_at_0_5'] / fx['pairs']:.1%}) with the authors' folds", "the authors' CV groups file")
    add("pair", "T2 vs T1 with the authors' folds", "calls flipped at 0.5", int(fx["flips_at_0_5"]), 0, "C", "results/audit/fixed_folds_summary.csv",
        f"{fx['flip_rate_same_model']:.1%} of pairs flip where the best model is the same, {fx['flip_rate_other_model']:.1%} where it differs (20 phages pick another model in a close AUPR race)", "the authors' model-selection seed and library versions")
    for tier in ("T3", "T4"):
        r = ra[ra.tier == tier].iloc[0]
        add("pair", f"{tier} vs T2", "calls flipped at 0.5", int(r.flips), 0, "F", "results/audit/pairs_regeneration_summary.csv; swap_groups.csv",
            f"{r.flip_rate:.1%} of pairs; with the folds fixed the regenerated features alone flip {sw.loc[sw.group == 'all'].call_flip_rate.iloc[0]:.1%}, led by UMAP ({sw.loc[sw.group == 'UMAP'].call_flip_rate.iloc[0]:.1%}) and O-type ({sw.loc[sw.group == 'O-type'].call_flip_rate.iloc[0]:.1%})", "nothing; a measured effect of regeneration")
    d = pd.DataFrame(rows)
    d.insert(0, "id", [f"P{i + 1:02d}" for i in range(len(d))])
    d.to_csv(f"{OUT}/gap_register_prediction.csv", index=False)
    print(d[["id", "level", "comparison", "metric", "cause_code", "size"]].to_string(index=False))
    return d


# ------------------------------------------------------------------ T76: check against an independent preprint (Noonan et al., bioRxiv 10.1101/2025.11.15.688630, version posted 2026-04-30, not peer reviewed)
def external_check():
    """Numbers the preprint reports on the same E. coli matrix (Table 1, Table 2, Results) against what this audit finds on the shipped matrix.
    Preprint values are typed from its text; our side is computed here. Writes results/audit/external_check.csv"""
    M = pd.read_csv(f"{ROOT}/data/interactions/interaction_matrix.csv", sep=";", index_col=0)
    L = (M > 0).fillna(False).values
    n, m = L.shape
    lytic = int(L.sum())
    cnt = pd.Series(L.sum(axis=0), index=M.columns)
    pairs = [(a, b) for i, a in enumerate(cnt.index) for b in cnt.index[i + 1:] if cnt[a] + cnt[b] == lytic - 7813]
    rng = np.random.default_rng(0)
    base = {k: [] for k in (1, 5)}
    for _ in range(50):
        for f in np.array_split(rng.permutation(n), 20):
            tr = np.setdiff1d(np.arange(n), f)
            o = np.argsort(-L[tr].mean(axis=0))
            for k in base:
                base[k].append(L[np.ix_(f, o[:k])].any(axis=1).mean())
    pooled = rd_overall("pooled_auroc")
    rows = [
        ("Strains in the E. coli matrix", "402 (Table 1)", f"{n} (public matrix)", "agrees", "both lack strain LF110; the paper has 403"),
        ("Phages", "94 (Table 1)", f"{m}", "differs", f"the preprint drops 2 phages and does not name them; {len(pairs)} pairs of phages would give its positive count"),
        ("Interactions", "37,788 = 402 x 94 (Table 1)", f"{n * 94:,} cells for 94 phages", "agrees", "arithmetic only"),
        ("Positive interactions", "7,813 (20.7 %) (Table 1)", f"{lytic:,} of {n * m:,} cells ({lytic / (n * m):.1%}); 7,813 after removing two phages with 163 lytic pairs", "agrees", "the README's '19 % positive' is the odd one out (report 11, item 8)"),
        ("Definition of a positive", "any non-zero interaction (Methods)", "MLC > 0", "agrees", "the same convention as the authors' script"),
        ("AUROC for unseen strains", "0.869 (Table 2); 0.87 in the text against 0.86 for the original study", f"pooled AUROC {pooled['T1']:.3f} (authors), {pooled['T2']:.3f} (this run), 0.86 (paper)", "agrees within 0.013",
         "different method (CatBoost on protein families, 20-fold hierarchical-cluster CV with 10 % of strains held out); the text does not say how the AUROC is pooled, and Table 2 does not say whether 'E. coli' is the 402-strain or the 177-strain set"),
        ("Best fixed single phage, share of strains infected", "about 61 % implied (66.9 % for the model's pick, 1.1-fold over promiscuity; rounding gives 58 to 64 %)", f"{np.mean(base[1]):.1%} (promiscuity chosen on training strains, 20-fold CV, 50 repeats); ED7 Source Data says 42.7 %",
         "supports this audit, not ED7", "with 'any non-zero' as positive, one promiscuous phage covers 61 to 65 % of strains in both analyses, far above ED7's 42.7 %"),
        ("Promiscuity baseline, 5 phages", "about 81 % implied (97.5 % for the model, 1.2-fold; rounding gives 78 to 85 %)", f"{np.mean(base[5]):.1%} (top 5 by promiscuity, no clustering)", "consistent",
         "the preprint's baseline picks one phage per cluster, so it is expected to sit slightly lower"),
        ("Model-guided 5-phage cocktail", "97.5 % of E. coli strains", "no counterpart", "not comparable", "per-strain selection by a model; ED7 is a fixed set"),
        ("Experimental validation", "AUROC 0.84 on 1,240 interactions (52 BASEL phages x 25 ECOR strains)", "no counterpart", "not comparable", "new phages, a different laboratory"),
        ("Arithmetic inside the preprint", "totals 115,037 interactions, 949 strains, 27,081 and 4,304 development, 1,240 validation", "recomputed from its Table 1 and Methods", "reconciles, with one slip",
         "Klebsiella-3 has 25 phages in Methods and 52 in Table 1; only 52 gives the stated 518 phages")]
    d = pd.DataFrame(rows, columns=["item", "preprint", "this audit", "verdict", "note"])
    d.to_csv(f"{OUT}/external_check.csv", index=False)
    print(d.to_string(index=False))
    return d


def rd_overall(metric):
    o = pd.read_csv(f"{OUT}/overall.csv")
    return {t: o[(o.tier == t) & (o.metric == metric)].value.iloc[0] for t in ("T1", "T2")}


# ------------------------------------------------------------------ T77: the original preprint (v1, 2023-11-22) against the published paper (T0) and the authors' repository (T1)
def preprint_check():
    """Each row: what the preprint says, what the published paper says (T0) and what the authors' repository holds (T1).
    Preprint and paper values are typed from docs/paper/preprint.txt and published.txt (line numbers in brackets); T1 values are computed here."""
    from sklearn.metrics import f1_score
    M = pd.read_csv(f"{ROOT}/data/interactions/interaction_matrix.csv", sep=";", index_col=0)
    L = (M > 0)
    n, m = M.shape
    tested = int(M.notna().values.sum())
    pp, pr = P.published()
    y = pr.y.astype(int)
    call = pr.p >= 0.5
    tp, fp, fn = int((call & (y == 1)).sum()), int((call & (y == 0)).sum()), int((~call & (y == 1)).sum())
    strain_auc = pr.groupby("bacteria").apply(lambda x: roc_auc_score(x.y, x.p) if x.y.nunique() > 1 else np.nan, include_groups=False).mean()
    cmp = pd.read_csv(f"{ROOT}/results/comparison/paper_vs_replication.csv")
    pvr = lambda fig, met: str(cmp[(cmp.figure == fig) & cmp.metric.str.startswith(met)].iloc[0].replication_shipped_inputs)
    B = pd.read_csv(f"{ROOT}/data/genomics/bacteria/picard_collection.csv", sep=";", index_col=0)
    lps = B.LPS_type.value_counts()
    ps, pn = L.sum(axis=0), L.sum(axis=1)
    rows = [
        ("Matrix size", "403 x 96 = 38,688 [l.242]", "403 x 96 [l.570]", f"{n} x {m} public matrix: {n * m:,} cells, {tested:,} tested; the shipped predictions cover {len(pr):,} pairs", "preprint and paper agree; T1 lacks LF110"),
        ("Share of interactions that are lytic", "20.3 % [l.245]", "20.3 % [l.344]", f"{L.values.sum() / tested:.2%} of tested cells", "T1 is 0.5 points higher; the paper's own 20.6 % and 20.8 % imply 20.7 %"),
        ("Strains infected per phage", "83 +/- 70 (20.6 %) [l.246]", "83 +/- 70 (20.6 %) [l.99]", f"{ps.mean():.1f} +/- {ps.std():.1f} ({ps.mean() / n:.1%})", "same in all three"),
        ("Phages lysing a strain", "20 +/- 15 (20.8 %) [l.247]", "20 +/- 15 (20.8 %) [l.100]", f"{pn.mean():.1f} +/- {pn.std():.1f} ({pn.mean() / m:.1%})", "same in all three"),
        ("Strains lysed by no phage", "12 (3.0 %) [l.248]", "12 (3.0 %) [l.102]", f"{int((pn == 0).sum())}", "same in all three"),
        ("One-way ANOVA of host range by morphotype", "F = 114.7, P < 1e-3 [l.256]", "F = 114.7 [l.339]; P = 1.44e-26 in the Fig 2f legend [l.460]", "F = 113.9 on the public matrix, 114.66 with LF110", "same once LF110 is added"),
        ("PERMANOVA of phage host range (adjusted R2)", "genus 0.13, isolation strain 0.53, both 0.65, defense 0.02 [l.308-311, l.432]", "same [l.365-369, l.424]", "no code in the repository", "same in both texts; T1 has no code"),
        ("LPS types", "K-12 43, R1 206, R2 22, R3 61, R4 35 [l.206]", "not in the extracted text", f"K-12 {lps.get('K12', 0)}, R1 {lps.get('R1', 0)}, R2 {lps.get('R2', 0)}, R3 {lps.get('R3', 0)}, R4 {lps.get('R4', 0)}", "T1 matches the preprint"),
        ("Strains with a Klebsiella capsule", "22 of 403 (5.4 %) [l.149]", "not stated", f"{int(B.Klebs_capsule_type.notna().sum())} of {len(B)}", "match (22/403 is 5.5 %)"),
        ("Mean AUROC of Fig 5a", "77 %, per phage [l.521]", "85 %, per strain [l.623]", f"per-phage mean {pp.auroc.mean():.3f}, per-strain mean {strain_auc:.3f}", "T1 matches the preprint's per-phage 77 %; the paper changed the unit"),
        ("Overall AUROC", "86 % in the Results [l.522], 85 % in the Discussion [l.672]", "86 % in the Results and the Discussion [l.625, l.824]", f"pooled {roc_auc_score(pr.y, pr.p):.3f}", "the preprint contradicts itself; paper and T1 agree on 86 %"),
        ("Average precision and F1", "not reported", "60 % and 60 % [l.625]", f"{average_precision_score(pr.y, pr.p):.3f} and {f1_score(y, call):.3f}", "paper and T1 agree"),
        ("False negatives and false positives", "3,379 and 2,922 of 38,688 [l.529]", "3,379 and 2,922 of 38,688 [l.633]", f"{fn:,} and {fp:,} of {len(pr):,} pairs at p >= 0.5", f"within 2 %; the totals of errors are 6,301 (paper) and {fn + fp:,} (T1)"),
        ("Most covering triplet and 4-set", "63 % and 68 % [l.575-576]", "63 % [l.802]; 4-set not stated; ED7 Source Data 67.0 % and 73.7 %", f"no code; the matrix gives {pvr('ED7', 'coverage at k = 3')} % and {pvr('ED7', 'coverage at k = 4')} % (any lysis)", "the published text drops the 68 %; the Source Data disagrees with both texts"),
        ("K1 rule: the three phages, then the other Vectreviruses and Kayfunaviruses", "precision 73 %, recall 72 %; 7 % and 14 % [l.343-345]", "same, in Supplementary Text 6", f"{pvr('SI Text 6', 'K1 rule, the three phages: precision')} % and {pvr('SI Text 6', 'K1 rule, the three phages: recall')} %; {pvr('SI Text 6', 'K1 rule, the other 28 Vectreviruses and Kayfunaviruses: precision')} % and {pvr('SI Text 6', 'K1 rule, the other 28 Vectreviruses and Kayfunaviruses: recall')} %", "same in both texts; T1 differs"),
        ("Tailored and generic recommendations", "15 cocktails for 24 strains, 3 cocktails for 76 [l.607-609]", "same [l.783-786]; Supplementary Table 7 gives 28 and 72", pvr("Fig 6", "strains with tailored"), "both texts match T1; Supplementary Table 7 differs"),
        ("Cocktail success: tailored, generic, baseline", "91.67, 78.95 and 81.00 % [l.628-630]", "91.67 %; 'baseline and generic 81 % and 78 %' [l.806-811]; 78.95 % in Fig 6d", f"{pvr('Fig 6d', 'success rate, tailored')}, {pvr('Fig 6d', 'success rate, generic')} and {pvr('Fig 6d', 'success rate, baseline')} %", "same; the published text truncates 78.95"),
        ("Mann-Whitney tests", "U = 399.5, P = 0.02; U = 1171.5, P = 0.02 [l.631-632]", "same text [l.812-815]; Source Data P = 0.019 and 0.034", f"U = {pvr('Fig 6e', 'Mann-Whitney U, tailored vs baseline')} and {pvr('Fig 6e', 'Mann-Whitney U, tailored vs generic')}; P = {pvr('Fig 6e', 'P, tailored vs baseline')} and {pvr('Fig 6e', 'P, tailored vs generic')}", "the U values agree; the text's 0.02 for tailored against generic should read 0.03"),
        ("Phylogroups covered by the isolation strains", "75 % [l.176]", "66 % [l.73]", "not checkable", "changed between versions"),
        ("Agreement of annotated scores across replicates", "88.64 %, 9.71 %, 98.35 % [l.815-817]", "88.6 %, 9.7 %, 98.3 %", "the raw table is binary: 90.7 %, 9.3 %, 100 %", "not checkable: the finer scores are not shipped"),
        ("Per-phage models", "LogReg L2 and L1, RF depth 3 and 6 with 250 trees, class weights by prevalence, scikit-learn 1.1.2 [l.984-992]", "same [l.1526-1562]", "same four models plus a Dummy baseline (never best in T1); our rerun uses scikit-learn 1.9.1", "texts agree; the code has 5 models, the texts 4"),
        ("Recommender random forests", "depth 4, 80 trees, class weight 1.5 [l.1014, l.1021]", "same [l.1608, l.1629]", "README the same", "same in all three"),
        ("Rule that picks the best model", "mean AUROC [l.1002]", "mean AUROC (report 11, item 15)", "the script uses mean AUPR", "texts say AUROC, code uses AUPR")]
    d = pd.DataFrame(rows, columns=["item", "preprint v1", "published paper (T0)", "repository (T1)", "verdict"])
    d.to_csv(f"{OUT}/preprint_check.csv", index=False)
    print(d[["item", "repository (T1)", "verdict"]].to_string(index=False))
    return d


def fold_mechanism():
    """Why the fold assignments differ. (1) With the groups the authors' folds split taken apart, the authors' folds are exactly what greedy GroupKFold packing of our groups
    produces for some tie order, for every phage. (2) With our group labels the tie order does not reproduce them: neither the installed scikit-learn nor a copy of the 1.1.2 algorithm."""
    from sklearn.metrics import adjusted_rand_score as ari
    import collections
    f = pd.read_csv(f"{OUT}/authors_folds.csv")
    G = pd.read_csv(f"{ROOT}/data/derived/cv_groups_1e-4.csv", sep=";", index_col=0)["group"]
    f["group"] = f.bacteria.map(G)
    nf = f.groupby(["group", "phage"]).fold.nunique().groupby("group").max()
    split = set(nf[nf > 1].index)
    f["g2"] = np.where(f.group.isin(split), "S_" + f.bacteria, f.group.astype(str))

    def packing(col):
        ok_phage, blocks, ok_blocks = 0, 0, 0
        for _, d in f.groupby("phage"):
            gs = d.groupby(col).agg(size=("bacteria", "size"), fold=("fold", lambda s: s.mode().iloc[0]))
            loads, all_ok = np.zeros(10), True
            for S in sorted(gs["size"].unique(), reverse=True):
                blk = gs[gs["size"] == S]
                seq = []
                for _ in range(len(blk)):
                    a = int(np.argmin(loads)); seq.append(a); loads[a] += S
                blocks += 1
                if collections.Counter(seq) == collections.Counter(blk.fold.astype(int)): ok_blocks += 1
                else:
                    all_ok = False
                    loads = loads - np.array([seq.count(i) * S for i in range(10)]) + np.array([int((blk.fold == i).sum()) * S for i in range(10)])
            ok_phage += all_ok
        return ok_phage, blocks, ok_blocks

    def gkf_1_1_2(groups, n_splits=10):   # GroupKFold.split of scikit-learn 1.1.2: groups by size (np.argsort, reversed), each into the lightest fold
        _, gi = np.unique(groups, return_inverse=True)
        n_per = np.bincount(gi); idx = np.argsort(n_per)[::-1]
        load, g2f = np.zeros(n_splits), np.zeros(len(n_per))
        for k, w in enumerate(n_per[idx]):
            lf = np.argmin(load); load[lf] += w; g2f[idx[k]] = lf
        return g2f[gi].astype(int)
    ours = pd.read_csv(f"{OUT}/folds_shipped.csv")
    a_new, a_old = [], []
    for ph, d in f.groupby("phage"):
        d = d.sort_values("bacteria")
        a_old.append(ari(d.fold, gkf_1_1_2(d.group.values)))
        o = ours[ours.phage == ph].set_index("bacteria").fold.reindex(d.bacteria).values
        ok = ~np.isnan(o)
        a_new.append(ari(d.fold.values[ok], o[ok].astype(int)))
    p1, b1, k1 = packing("group"); p2, b2, k2 = packing("g2")
    out = pd.Series(dict(our_groups=int(G.nunique()), our_groups_split_by_authors_folds=len(split), strains_in_them=int(G.isin(split).sum()),
                         phages_packing_consistent_original_groups=p1, phages_packing_consistent_after_splitting=p2, phages=int(f.phage.nunique()),
                         size_blocks_consistent_after_splitting=k2, size_blocks=b2,
                         ari_installed_sklearn=float(np.mean(a_new)), ari_installed_sklearn_max=float(np.max(a_new)), ari_sklearn_1_1_2_copy=float(np.mean(a_old)), ari_sklearn_1_1_2_copy_max=float(np.max(a_old))))
    out.to_csv(f"{OUT}/fold_mechanism.csv", header=["value"])
    print(out.to_string())
    return out


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    if step in ("register", "all"):
        register()
        paper_repo()
    if step == "paper_repo":
        paper_repo()
    if step == "folds":
        folds()
    if step == "authors_folds":
        authors_folds()
    if step == "external":
        external_check()
    if step == "preprint":
        preprint_check(); fold_mechanism()
    if step == "final":
        ed7_definitions(); permanova_checks(); finish_register(); prediction_register()
    if step == "runs":
        fixed_folds(); seed_spread(); swap_groups(); causes_runs()
    if step in ("phage", "all"):
        per_phage(); causes()
    if step in ("pairs", "all"):
        pairs_t1_t2(); pairs_regeneration()
    if step in ("align", "overall", "all"):
        t = tiers()
        if step in ("align", "all"):
            align(t)
        if step in ("overall", "all"):
            overall(t)
