"""T81-T82 summary: the model step in the pinned environment (scikit-learn 1.1.2, Python 3.10) against the same step in the main environment (scikit-learn 1.9.1, Python 3.14).
    python scripts/audit_pinned.py   ->  results/audit/pinned_env_{overall,per_phage,folds}.csv   (needs scripts/run_pinned_env.sh and scripts/fold_version_check.py runs)"""
import os, re, sys
import numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score as ari, roc_auc_score, average_precision_score
os.environ["FEATURE_SET"] = "shipped"
sys.path.insert(0, os.path.dirname(__file__))
import predict as P

R, OUT = P.ROOT, f"{P.ROOT}/results/audit"
norm = lambda m: re.sub(r"_weight=[0-9.]+", "", m)[:2]   # model family RF / Lo / Du, as in the audit
ship_pp, ship_pr = P.published()
RUNS = {  # label -> (directory, what it is)
    "T1": (None, "authors' shipped predictions (scikit-learn 1.1.2, their folds)"),
    "T2": ("repro", "shipped inputs, scikit-learn 1.9.1, our folds"),
    "T2.1": ("audit_runs/T2_1_shipped", "shipped inputs, scikit-learn 1.1.2, our groups"),
    "T2 authors' folds": ("audit_runs/fixed_default", "shipped inputs, scikit-learn 1.9.1, the authors' folds"),
    "T2.1 authors' folds": ("audit_runs/T2_1_shipped_authors_folds", "shipped inputs, scikit-learn 1.1.2, the authors' folds"),
    "T3": ("repro_regenerated", "regenerated features, scikit-learn 1.9.1"),
    "T3.1": ("audit_runs/T3_1_regenerated", "regenerated features, scikit-learn 1.1.2"),
}
load = lambda k: (ship_pp, ship_pr) if RUNS[k][0] is None else (pd.read_csv(f"{R}/results/{RUNS[k][0]}/per_phage.csv", index_col=0), pd.read_csv(f"{R}/results/{RUNS[k][0]}/predictions.csv"))
D = {k: load(k) for k in RUNS}

# overall
rows = []
for k, (pp, pr) in D.items():
    j = pp.join(ship_pp, rsuffix="_T1")
    r = dict(run=k, what=RUNS[k][1], mean_phage_auroc=pp.auroc.mean(), mean_phage_aupr=pp.aupr.mean(), pooled_auroc=roc_auc_score(pr.y, pr.p), pooled_aupr=average_precision_score(pr.y, pr.p))
    if k != "T1":
        r.update(mean_abs_d_auroc_vs_T1=(j.auroc - j.auroc_T1).abs().mean(), corr_auroc_vs_T1=j.auroc.corr(j.auroc_T1),
                 same_model_family_as_T1=(j.best_model.map(norm) == j.best_model_T1.map(norm)).mean())
    rows.append(r)
ov = pd.DataFrame(rows).set_index("run")
# pairwise: what one change (the scikit-learn version, the folds, the features) does
pairs = [("T2.1", "T2", "scikit-learn 1.1.2 against 1.9.1, our groups, shipped inputs"), ("T2.1 authors' folds", "T2 authors' folds", "scikit-learn 1.1.2 against 1.9.1, the authors' folds"),
         ("T2 authors' folds", "T2", "the authors' folds against ours, scikit-learn 1.9.1"), ("T2.1 authors' folds", "T2.1", "the authors' folds against ours, scikit-learn 1.1.2"),
         ("T3.1", "T3", "scikit-learn 1.1.2 against 1.9.1, regenerated features"), ("T3.1", "T2.1", "regenerated against shipped features, scikit-learn 1.1.2"), ("T3", "T2", "regenerated against shipped features, scikit-learn 1.9.1")]
pw = []
for a, b, what in pairs:
    d = D[a][0].auroc - D[b][0].auroc
    pa, pb = D[a][1].set_index(["phage", "bacteria"]), D[b][1].set_index(["phage", "bacteria"])
    jj = pa.p.rename("pa").to_frame().join(pb.p.rename("pb"), how="inner")
    pw.append(dict(change=what, a=a, b=b, mean_d_auroc=d.mean(), mean_abs_d_auroc=d.abs().mean(), max_abs_d_auroc=d.abs().max(), phages_above_0_05=int((d.abs() > 0.05).sum()),
                   same_best_model_family=(D[a][0].best_model.map(norm) == D[b][0].best_model.map(norm)).mean(), pairs=len(jj), flips_at_0_5=int(((jj.pa >= .5) != (jj.pb >= .5)).sum()),
                   spearman_p=jj.pa.corr(jj.pb, method="spearman")))
pw = pd.DataFrame(pw)
pd.concat([ov.reset_index().assign(section="run"), pw.assign(section="change")]).to_csv(f"{OUT}/pinned_env_overall.csv", index=False)
print(ov.round(4).to_string()); print(pw.round(4).to_string())

# per phage
pph = pd.DataFrame({k: D[k][0].auroc for k in RUNS}); pph.insert(0, "n_lytic", ship_pp.reindex(pph.index).apply(lambda s: s) if False else D["T2"][0].n_pos)
for k in RUNS:
    if k != "T1": pph[f"d_{k}_T1"] = pph[k] - pph["T1"]
pph.rename_axis("phage").to_csv(f"{OUT}/pinned_env_per_phage.csv")

# folds: does the scikit-learn version change the partition?
af = pd.read_csv(f"{OUT}/authors_folds.csv")
fr = []
for s in ("shipped", "regenerated"):
    a, b = (pd.read_csv(f"{R}/results/audit_runs/fold_version/{n}_{s}.csv") for n in ("sk191", "sk112"))
    j = a.merge(b, on=["phage", "bacteria"], suffixes=("_191", "_112"))
    row = dict(feature_set=s, pairs=len(j), same_fold_index_share=(j.fold_191 == j.fold_112).mean(), phages_identical=int(j.groupby("phage").apply(lambda d: (d.fold_191 == d.fold_112).all()).sum()),
               mean_ari_191_vs_112=j.groupby("phage").apply(lambda d: ari(d.fold_191, d.fold_112)).mean())
    if s == "shipped":   # the authors' folds exist only for the shipped pairs
        for tag, c in (("191", a), ("112", b)):
            m = c.merge(af, on=["phage", "bacteria"], suffixes=("", "_authors"))
            row[f"mean_ari_{tag}_vs_authors"] = m.groupby("phage").apply(lambda d: ari(d.fold, d.fold_authors)).mean()
    fr.append(row)
fr = pd.DataFrame(fr); fr.to_csv(f"{OUT}/pinned_env_folds.csv", index=False); print(fr.round(4).to_string())

# pair-level agreement with the authors' own probabilities, for the phages whose best model family is the same: only the deterministic LogReg L2 can match exactly
# (the authors' RandomForest, Dummy and saga LogReg have no random_state, see dev/predictions/predict_all_phages.py)
out, per = [], []
for tag, k in (("T2.1 authors' folds, scikit-learn 1.1.2", "T2.1 authors' folds"), ("T2 authors' folds, scikit-learn 1.9.1", "T2 authors' folds")):
    pp, pr = D[k]
    same = ship_pp.best_model.map(norm) == pp.best_model.reindex(ship_pp.index).map(norm)
    j = ship_pr[["phage", "bacteria", "p", "model"]].merge(pr[["phage", "bacteria", "p"]], on=["phage", "bacteria"], suffixes=("1", "2"))
    j["fam"], j["same"], j["dp"] = j.model.map(lambda m: re.sub(r"_weight=[0-9.]+", "", m)), j.phage.map(same), (j.p1 - j.p2).abs()
    for fam, g in j[j.same].groupby("fam"):
        out.append(dict(run=tag, model=fam, phages=g.phage.nunique(), pairs=len(g), share_identical_1e6=(g.dp < 1e-6).mean(), share_within_1e3=(g.dp < 1e-3).mean(), median_abs_dp=g.dp.median(), max_abs_dp=g.dp.max()))
        if fam == "LogReg_l2" and k == "T2.1 authors' folds":
            per.append(g.groupby("phage").agg(pairs=("dp", "size"), share_identical_1e6=("dp", lambda x: (x < 1e-6).mean()), max_abs_dp=("dp", "max")).join(pp.n_pos.rename("n_lytic")))
pd.DataFrame(out).to_csv(f"{OUT}/pinned_env_pair_agreement.csv", index=False)
pd.concat(per).sort_values("share_identical_1e6").rename_axis("phage").to_csv(f"{OUT}/pinned_env_logreg_l2_phages.csv")
print(pd.DataFrame(out).round(4).to_string())

# component by component: what the authors state or the code shows, what the main runs used, what the pinned runs used, and why a mismatch remains
ENV = [
    ("scikit-learn", "1.1.2 (Methods, statistics paragraph)", "1.9.1", "1.1.2", "same in pinned", "", "paper l.1562"),
    ("Python", "not stated (the paper, the prediction script and the repository give no Python version for the model step; the only Python version in the repository is 3.12.2 in the kernel metadata of the plotting notebook visualize_predictions.ipynb, which does not say what trained the models)", "3.14.7", "3.10.22", "differs",
     "scikit-learn 1.1.2 has no Python 3.12 build (conda-forge osx-arm64 has py310 only; the release dates from 2022). The version that fixes the results is scikit-learn's, so Python follows it. The authors' Python for the model step is unknown, so 3.10 is a choice, not a match or a mismatch. The notebook's 3.12.2 cannot be one environment with scikit-learn 1.1.2, so it does not describe the model step if the paper's version is right (inference)", "ml_pinned solver; notebook metadata"),
    ("scipy", "1.11.2 (Methods, scipy.stats)", "1.18.1", "1.11.2", "same in pinned (stated for scipy.stats, assumed for the model step)", "", "paper l.1650"),
    ("numpy", "not stated", "2.5.4", "1.23.5", "not stated, chosen", "release contemporary with scikit-learn 1.1.2; the conda build of 1.1.2 needs numpy < 2", ""),
    ("pandas", "not stated", "3.0.6", "1.5.3", "not stated, chosen", "same reasoning; pandas 3 postdates scikit-learn 1.1.2", ""),
    ("joblib / threadpoolctl", "not stated", "1.6.0 / -", "1.2.0 / 3.6.0", "not stated, chosen", "joblib contemporary with scikit-learn 1.1.2; threadpoolctl is pulled in by the solver", ""),
    ("LogisticRegression L1 argument", 'penalty="l1", solver="saga", max_iter=10000', 'l1_ratio=1 (penalty is deprecated from scikit-learn 1.8)', 'penalty="l1"', "same in pinned",
     'In 1.1.2, l1_ratio=1 without penalty="l1" silently fits L2 (checked); predict.py now chooses the argument by version', "dev/predictions/predict_all_phages.py l.252; scripts/predict.py"),
    ("model randomness", "np.random.seed(0) once at the top; RandomForest, Dummy and saga without random_state, phages in sequence, so each draw depends on every earlier draw", "random_state=0 for RF and Dummy, saga unseeded; phages in parallel", "same as the main runs",
     "not replicable", "the stream state at any model depends on the draws of all earlier models, including DummyClassifier(stratified), whose draws depend on the authors' private labels and on row order. Measured consequence: probabilities of RandomForest-best phages differ even with their folds (median |dp| 0.008 to 0.018), mean AUROC noise 0.0011 (T72)", "dev/predictions/predict_all_phages.py l.23, 145-216, 247-253; pinned_env_pair_agreement.csv"),
    ("row and column order of the training matrix", "rows sorted by bacteria then phage; columns in the order of their feature file", "order of the interaction matrix", "same as the main runs", "differs, not matched",
     "bootstrap and feature draws of the forests depend on positions; could be matched for the shipped set, not attempted because the random stream (above) cannot be matched anyway", "predict_all_phages.py l.59-100"),
    ("interaction labels", "private N/P/F/U matrix (ADSORPTION_MATRIX_370+host_vs_96_less_stringent.csv), U counted as lytic, F dropped", "lytic = MLC > 0 on the public 0-4 matrix", "same as the main runs", "not available", "the authors' matrix is not published; 59 of 38,435 pairs differ in presence (58 only in ours, 1 only in theirs); label differences are not measurable", "config/config_diff.csv"),
    ("CV groups and folds", "370+host_cross_validation_groups_1e-4.csv, GroupKFold(10) of scikit-learn 1.1.2", "rebuilt groups, GroupKFold(10) of 1.9.1", "rebuilt groups, GroupKFold(10) of 1.1.2; the authors' own folds as a second run", "groups not available",
     "the groups file is not published. The same groups give different folds under 1.1.2 and 1.9.1 (adjusted Rand index 0.047), and neither is close to the authors' (0.030 and 0.022). Their own folds are available from their prediction files and are used in the 'authors folds' runs", "pinned_env_folds.csv"),
    ("operating system, CPU, BLAS", "Windows (paths D:\\\\ in the script); CPU and BLAS not stated", "macOS arm64 (Apple M2), conda OpenBLAS", "macOS arm64 (Apple M2), conda OpenBLAS", "not replicable", "no Windows machine; floating-point differences of about 1e-4 in the probabilities of deterministic LogReg L2 on a few phages are consistent with this (not isolated)", "pinned_env_logreg_l2_phages.csv"),
    ("R, MCMCglmm, vegan, statsmodels", "MCMCglmm 2.35, vegan 2.6-4, statsmodels 0.14.1, R version not stated", "R 4.3.0, MCMCglmm 2.36; vegan and statsmodels not used", "not part of this rerun", "unchanged", "the pinned environment covers only the model step; see config/config_diff.csv for the rest", "config/config_diff.csv"),
]
pd.DataFrame(ENV, columns=["component", "authors", "main_runs", "pinned_runs", "status", "reason", "source"]).to_csv(f"{OUT}/env_match.csv", index=False)
print(f"wrote env_match.csv ({len(ENV)} rows)")
