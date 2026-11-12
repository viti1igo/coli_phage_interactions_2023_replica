"""Phase 21: figures and test-fold metrics for the four runs on the authors' shipped inputs:
T1 (authors' predictions), T2.1 (my run, scikit-learn 1.1.2), B2 (biophage code, paper-like env), B1 (biophage as written).
Metrics follow Fig 5 of the paper (cross-validated test folds, best model per phage by mean test AUPR, threshold 0.5).
Needs scripts/compare_biophage.py run first.   .venv/bin/python scripts/fig_biophage.py"""
import glob
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import average_precision_score, f1_score, precision_recall_curve, precision_score, recall_score, roc_auc_score, roc_curve
from figlib import AXIS, BLUE, INK, MUTED, ORANGE, ROOT, SURF, save, style

AQUA = "#1baf7a"  # palette slot 3
COL = {"T1": INK, "T2.1": BLUE, "B2": AQUA, "B1": ORANGE}
LAB = {"T1": "T1 authors' predictions", "T2.1": "T2.1 my run, scikit-learn 1.1.2", "B2": "B2 biophage, paper-like env", "B1": "B1 biophage as written"}
A = f"{ROOT}/results/audit"
M = pd.read_csv(f"{ROOT}/data/interactions/interaction_matrix.csv", sep=";", index_col=0)


def authors():
    out = []
    for f in glob.glob(f"{ROOT}/dev/predictions/results/performances/performance_*.csv"):
        d = pd.read_csv(f, sep=";"); d = d[d.dataset == "test"].dropna(subset=["model"])
        best = d.groupby("model").avg_precision.mean().idxmax(); ph = d.phage.iloc[0]
        p = pd.read_csv(glob.glob(f"{ROOT}/dev/predictions/results/predictions/predictions_{ph}_*.csv")[0], sep=";")
        out.append(p[(p.dataset == "test") & (p.model == best)].rename(columns={"y_pred_proba": "p"})[["phage", "bacteria", "y", "p"]])
    return pd.concat(out)


def biophage(env):
    d = f"{ROOT}/replication/T93_biophage/{env}/training_results"
    s = pd.read_csv(f"{d}/per_phage_perf.csv")
    p = pd.read_csv(f"{d}/predictions_named.csv").merge(s[["phage", "model"]], on=["phage", "model"])
    return p.rename(columns={"y_true": "y", "y_pred_proba": "p"})


preds = {"T1": authors(), "T2.1": pd.read_csv(f"{ROOT}/results/audit_runs/T2_1_shipped/predictions.csv"),
         "B2": biophage("bp_pinned"), "B1": biophage("bp_asis")}
nan_pairs = {(b, ph) for ph in M for b in M.index[M[ph].isna()]}
pp = pd.read_csv(f"{A}/biophage_compare_per_phage.csv", index_col=0)   # per-phage AUROC = mean of fold AUROCs, best model
SD = f"{ROOT}/docs/paper/supplement/MOESM9.xlsx"   # Source Data of Fig 5 (plotted, not re-shipped in these CSVs)
SD5A = pd.read_excel(SD, "Figure_5a", header=1).set_index("bacteria").auc   # AUROC per strain, 315 strains
SD5B = pd.read_excel(SD, "Figure_5b", header=1)
# the ROC columns are labelled the wrong way round (TPR < FPR along the whole curve gives AUROC 0.15); read swapped
SD_FPR, SD_TPR = SD5B.TPR.values, SD5B.FPR.values
SD_AUROC = float(np.trapezoid(SD_TPR, SD_FPR))
PAPER = dict(per_strain=0.85, per_phage=0.77, pooled_auroc=0.86, ap=0.60, f1=0.60, fn=3379, fp=2922, pairs=38688)

rows, strain_auc = [], {}
for r, d in preds.items():
    all_pairs = len(d)
    d = d[[(b, ph) not in nan_pairs for b, ph in zip(d.bacteria, d.phage)]]   # paper scores testable pairs only
    preds[r] = d
    yhat = (d.p >= 0.5).astype(int)
    sa = d.groupby("bacteria").apply(lambda x: roc_auc_score(x.y, x.p) if x.y.nunique() > 1 else np.nan, include_groups=False).dropna()
    strain_auc[r] = sa
    rows.append(dict(run=r, pairs_scored=len(d), pairs_trained=all_pairs, mean_auroc_per_strain=sa.mean(), mean_auroc_per_phage=pp[r].mean(),
                     median_auroc_per_phage=pp[r].median(), pooled_auroc=roc_auc_score(d.y, d.p), average_precision=average_precision_score(d.y, d.p),
                     precision=precision_score(d.y, yhat), recall=recall_score(d.y, yhat), f1=f1_score(d.y, yhat),
                     tp=int(((yhat == 1) & (d.y == 1)).sum()), fp=int(((yhat == 1) & (d.y == 0)).sum()), fn=int(((yhat == 0) & (d.y == 1)).sum()),
                     tn=int(((yhat == 0) & (d.y == 0)).sum()), mean_auroc_per_strain_paper_315=sa.reindex(SD5A.index).mean()))
met = pd.DataFrame(rows)
met = pd.concat([pd.DataFrame([dict(run="T0 paper", pairs_scored=PAPER["pairs"], mean_auroc_per_strain=PAPER["per_strain"], mean_auroc_per_phage=PAPER["per_phage"],
                                    pooled_auroc=PAPER["pooled_auroc"], average_precision=PAPER["ap"], f1=PAPER["f1"], fn=PAPER["fn"], fp=PAPER["fp"]),
                              dict(run="T0 paper Source Data", pairs_scored=np.nan, mean_auroc_per_strain=SD5A.mean(), mean_auroc_per_strain_paper_315=SD5A.mean(),
                                   pooled_auroc=SD_AUROC)]), met])
met.to_csv(f"{A}/biophage_testfold_metrics.csv", index=False)
print(met.round(3).T)

# 1. per-phage AUROC against the authors'
fig, axs = plt.subplots(1, 3, figsize=(13, 4.3), sharex=True, sharey=True)
for ax, r in zip(axs, ["T2.1", "B2", "B1"]):
    ax.plot([0.3, 1], [0.3, 1], color=AXIS, lw=1, zorder=1)
    ax.scatter(pp.T1, pp[r], s=18, color=COL[r], edgecolor=SURF, lw=0.6, zorder=2)
    d = pp[r] - pp.T1
    ax.text(0.33, 0.97, f"mean |difference| {d.abs().mean():.3f}\nPearson r {pp[r].corr(pp.T1):.2f}\nphages off by > 0.05: {(d.abs() > 0.05).sum()} of 96",
            va="top", fontsize=8, color=INK)
    ax.set_title(LAB[r], loc="left"); ax.set_xlabel("Authors' AUROC per phage (T1)"); style(ax, "both")
axs[0].set_ylabel("AUROC per phage, this run")
fig.tight_layout(); save(fig, "biophage_1_per_phage_vs_T1", pp.reset_index())

# 2. Fig 5a style: distributions per strain and per phage
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for r in COL:
    for ax, v in [(axs[0], strain_auc[r]), (axs[1], pp[r])]:
        x = np.sort(v.values); ax.step(x, np.arange(1, len(x) + 1) / len(x), where="post", color=COL[r], lw=2, label=f"{LAB[r]} (mean {v.mean():.3f})")
x = np.sort(SD5A.values); axs[0].step(x, np.arange(1, len(x) + 1) / len(x), where="post", color=MUTED, lw=2, ls=(0, (3, 2)),
                                      label=f"T0 paper, Source Data Fig 5a, {len(x)} strains (mean {SD5A.mean():.3f})")
for ax, t, ref, lab in [(axs[0], "a  AUROC per strain (Fig 5a)", PAPER["per_strain"], "paper text 0.85"), (axs[1], "b  AUROC per phage", PAPER["per_phage"], "preprint 0.77")]:
    ax.axvline(ref, color=MUTED, lw=1, ls=(0, (4, 3))); ax.text(ref - 0.005, 0.02, lab, ha="right", fontsize=8, color=MUTED)
    ax.set_title(t, loc="left"); ax.set_xlabel("AUROC"); ax.set_ylabel("Cumulative share"); ax.set_xlim(0.3, 1.0); style(ax, "both")
    ax.legend(frameon=False, fontsize=7, labelcolor=INK, loc="upper left")
fig.tight_layout()
save(fig, "biophage_2_auroc_distributions", pd.concat({r: strain_auc[r] for r in COL}, names=["run", "bacteria"]).rename("auroc_per_strain").reset_index())

# 3. Fig 5b style: pooled ROC and precision-recall
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
curves = []
for r, d in preds.items():
    f, t, _ = roc_curve(d.y, d.p); pr, rc, _ = precision_recall_curve(d.y, d.p)
    axs[0].plot(f, t, color=COL[r], lw=2, label=f"{LAB[r]}: {roc_auc_score(d.y, d.p):.3f}")
    axs[1].plot(rc, pr, color=COL[r], lw=2, label=f"{LAB[r]}: {average_precision_score(d.y, d.p):.3f}")
    curves.append(pd.DataFrame({"run": r, "fpr": f, "tpr": t}).iloc[::20])
axs[0].plot(SD_FPR, SD_TPR, color=MUTED, lw=2, ls=(0, (3, 2)), label=f"T0 paper, Source Data Fig 5b: {SD_AUROC:.3f}")
axs[0].plot([0, 1], [0, 1], color=AXIS, lw=1)
axs[0].set_title("a  Pooled ROC (paper text 0.86)", loc="left"); axs[0].set_xlabel("False positive rate"); axs[0].set_ylabel("True positive rate")
axs[1].set_title("b  Pooled precision-recall (paper text 0.60)", loc="left"); axs[1].set_xlabel("Recall"); axs[1].set_ylabel("Precision")
for ax, loc in zip(axs, ["lower right", "upper right"]): style(ax, "both"); ax.legend(frameon=False, fontsize=7, labelcolor=INK, loc=loc, title="AUROC / average precision", title_fontsize=7)
fig.tight_layout(); save(fig, "biophage_3_pooled_curves", pd.concat(curves))

# 4. difference to the authors by number of lytic strains
fig, ax = plt.subplots(figsize=(7.5, 4.3))
ax.axhline(0, color=AXIS, lw=1)
for r, dx in [("T2.1", -0.15), ("B2", 0), ("B1", 0.15)]:
    ax.scatter(np.log10(pp.n_lytic) + dx * 0.05, pp[r] - pp.T1, s=16, color=COL[r], edgecolor=SURF, lw=0.5, label=LAB[r])
ax.set_xticks(np.log10([2, 5, 10, 20, 50, 100, 200])); ax.set_xticklabels([2, 5, 10, 20, 50, 100, 200])
ax.set_xlabel("Lytic strains of the phage (log scale)"); ax.set_ylabel("AUROC per phage minus authors' (T1)")
ax.set_title("Difference to the authors by phage size", loc="left"); style(ax, "y"); ax.legend(frameon=False, fontsize=7, labelcolor=INK)
fig.tight_layout(); save(fig, "biophage_4_difference_by_size", pp.reset_index())

# 0. headline metrics, paper against the four runs
K = [("mean_auroc_per_strain", "AUROC per strain"), ("mean_auroc_per_phage", "AUROC per phage"), ("pooled_auroc", "Pooled AUROC"),
     ("average_precision", "Average precision"), ("f1", "F1 at 0.5")]
mt = met.set_index("run")
fig, axs = plt.subplots(1, 2, figsize=(13, 5.0), gridspec_kw={"width_ratios": [1.4, 1]})
ax = axs[0]
for j, (k, lab) in enumerate(K):
    ax.axhline(j, color="#e6e5e0", lw=0.6, zorder=0)
    ax.scatter(mt.loc["T0 paper", k], j, marker="D", s=70, facecolor="none", edgecolor=INK, lw=1.4, zorder=3, label="T0 paper (text)" if j == 0 else None)
    if k in ("mean_auroc_per_strain", "pooled_auroc"):
        ax.scatter(mt.loc["T0 paper Source Data", k], j, marker="s", s=46, facecolor="none", edgecolor=MUTED, lw=1.4, zorder=3, label="T0 paper (Source Data)" if j == 0 else None)
    for r in COL:
        ax.scatter(mt.loc[r, k], j, s=46, color=COL[r], edgecolor=SURF, lw=0.8, zorder=4, label=LAB[r] if j == 0 else None)
ax.set_yticks(range(len(K))); ax.set_yticklabels([l for _, l in K]); ax.invert_yaxis()
ax.set_xlabel("Value"); ax.set_title("a  Fig 5 metrics: paper against the four runs", loc="left"); style(ax, "x")
ax.legend(frameon=False, fontsize=7, labelcolor=INK, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
ax = axs[1]
runs = ["T0 paper"] + list(COL); w = 0.38; xs = np.arange(len(runs))
ax.bar(xs - w / 2, [mt.loc[r, "fn"] for r in runs], w, color=[MUTED] + [COL[r] for r in COL], edgecolor=SURF, label="false negatives")
ax.bar(xs + w / 2, [mt.loc[r, "fp"] for r in runs], w, color=[MUTED] + [COL[r] for r in COL], edgecolor=SURF, alpha=0.45, hatch="//", label="false positives")
for x, r in zip(xs, runs):
    ax.text(x - w / 2, mt.loc[r, "fn"] + 40, f"{mt.loc[r, 'fn']:,.0f}", ha="center", fontsize=7, color=INK, rotation=90)
    ax.text(x + w / 2, mt.loc[r, "fp"] + 40, f"{mt.loc[r, 'fp']:,.0f}", ha="center", fontsize=7, color=INK, rotation=90)
ax.set_xticks(xs); ax.set_xticklabels(["T0 paper", "T1", "T2.1", "B2", "B1"]); ax.set_ylim(0, 4600)
ax.set_ylabel("Pairs"); ax.set_title("b  Errors at 0.5 (Fig 5d)", loc="left"); style(ax, "y")
ax.legend(frameon=False, fontsize=7, labelcolor=INK, loc="upper left")
fig.tight_layout(); save(fig, "biophage_0_headline_vs_paper", met)
