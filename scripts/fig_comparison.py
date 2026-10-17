"""Phase 9 comparison figures; each is drawn from the CSV saved beside it (results/comparison/cmp_*.csv).
python scripts/fig_comparison.py  (after scripts/build_comparison.py)"""
import re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import figlib
from figlib import *

figlib.OUT = f"{ROOT}/results/comparison"
C = figlib.OUT

# 1. headline numbers as ratio replication / paper (p-values and non-numeric rows left out)
d = pd.read_csv(f"{C}/paper_vs_replication.csv")
num = lambda v: pd.to_numeric(v, errors="coerce")
first = lambda v: num(re.match(r"\s*(-?[\d.]+)", str(v)).group(1)) if re.match(r"\s*-?[\d.]", str(v)) else np.nan
d["paper_n"], d["rep_n"], d["regen_n"] = d.paper.map(first), d.replication_shipped_inputs.map(first), d.replication_regenerated_inputs.map(first)
d = d[d.paper_n.notna() & d.rep_n.notna() & (d.paper_n != 0) & ~d.metric.str.contains(r"\bP\b|pMCMC|P,", regex=True)].copy()
d["ratio_replication"] = d.rep_n / d.paper_n; d["ratio_regenerated"] = d.regen_n / d.paper_n
d["label"] = d.figure + "  " + d.metric.str.slice(0, 62)
d = d.iloc[::-1].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(8.5, 0.28 * len(d) + 1.4))
ax.axvline(1, color=MUTED, lw=1); ax.axvspan(0.95, 1.05, color="#e6e5e0", alpha=0.6, lw=0)
ax.scatter(d.ratio_replication, d.index, s=28, color=BLUE, label="Replication, shipped inputs", zorder=3, edgecolor=SURF, lw=0.5)
m = d.ratio_regenerated.notna()
ax.scatter(d.ratio_regenerated[m], d.index[m], s=28, color=ORANGE, label="Replication, regenerated inputs", zorder=3, edgecolor=SURF, lw=0.5)
ax.set_yticks(d.index); ax.set_yticklabels(d.label, fontsize=7); ax.set_xlabel("Our value / paper value (grey band = within 5%)"); ax.set_title("Headline numbers against the paper", loc="left")
ax.legend(frameon=False, labelcolor=INK, fontsize=8, loc="lower right"); style(ax, "x")
save(fig, "cmp_headline", d[["figure", "metric", "paper_n", "rep_n", "regen_n", "ratio_replication", "ratio_regenerated"]])

# 2. per-phage AUROC against the authors' shipped results
p = pd.read_csv(f"{C}/per_phage_auroc.csv")
fig, ax = plt.subplots(figsize=(5.4, 5.2))
ax.plot([0.4, 1], [0.4, 1], color=MUTED, lw=1, ls=(0, (4, 3)))
for col, c, lab in [("auroc_replication", BLUE, "Replication, shipped inputs"), ("auroc_regenerated", ORANGE, "Replication, regenerated inputs")]:
    r = p.auroc_shipped.corr(p[col])
    ax.scatter(p.auroc_shipped, p[col], s=20, color=c, edgecolor=SURF, lw=0.5, label=f"{lab} (r = {r:.2f})")
ax.set_xlabel("AUROC, authors' shipped results"); ax.set_ylabel("AUROC, our run"); ax.set_xlim(0.45, 1.02); ax.set_ylim(0.45, 1.02); ax.set_title("Per-phage AUROC, 96 phages", loc="left")
ax.legend(frameon=False, labelcolor=INK, fontsize=8, loc="lower right"); style(ax, "both")
save(fig, "cmp_per_phage_auroc", p)

# 3. regenerated features against the shipped ones
f = pd.read_csv(f"{C}/feature_agreement.csv")
pct = f[f.metric.str.contains("agreement %|accuracy %")].copy(); pct["value"] = pct.agreement_value
oth = f[~f.metric.str.contains("agreement %|accuracy %") & f.agreement_value.between(0, 1) & ~f.metric.str.contains("Robinson")].copy()
fig, axs = plt.subplots(1, 2, figsize=(12, 7.5), gridspec_kw={"width_ratios": [1.2, 1]})
pct["l"] = pct.feature.str.slice(0, 45); axs[0].barh(pct.l[::-1], pct.value[::-1], color=BLUE, height=0.65); axs[0].set_xlim(0, 101); axs[0].set_xlabel("Agreement with shipped values (%)"); axs[0].set_title("Typing features (strain-level agreement)", loc="left")
short = {"adjusted Rand index": "ARI", "10-NN overlap with shipped coordinates": "10-NN overlap", "patristic distance Pearson": "r", "F1": "F1", "Pearson vs shipped column": "r"}
oth["l"] = oth.feature.str.slice(0, 38) + "  [" + oth.metric.map(lambda m: short.get(m, m[:14])) + "]"; axs[1].barh(oth.l[::-1], oth.agreement_value[::-1], color=BLUE, height=0.65); axs[1].set_xlim(0, 1.02); axs[1].set_xlabel("Adjusted Rand index / F1 / r / overlap (0-1)"); axs[1].set_title("Clusterings, defense systems, tree", loc="left")
for a in axs: style(a, "x"); a.tick_params(axis="y", labelsize=7)
fig.tight_layout()
save(fig, "cmp_feature_agreement", pd.concat([pct.drop(columns="l"), oth.drop(columns="l")]))

# 4. ablation: shipped inputs vs regenerated inputs
a = pd.read_csv(f"{C}/ablation_comparison.csv").iloc[::-1].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(7.5, 5))
ax.axvline(0, color=MUTED, lw=1)
for i, r in a.iterrows(): ax.plot([r.d_auroc_shipped_inputs, r.d_auroc_regenerated_inputs], [i, i], color="#c9c8c2", lw=2, zorder=1)
ax.scatter(a.d_auroc_shipped_inputs, a.index, s=34, color=BLUE, label="Shipped inputs", zorder=3, edgecolor=SURF, lw=0.5)
ax.scatter(a.d_auroc_regenerated_inputs, a.index, s=34, color=ORANGE, label="Regenerated inputs", zorder=3, edgecolor=SURF, lw=0.5)
ax.set_yticks(a.index); ax.set_yticklabels(a.feature_group_removed, fontsize=8); ax.set_xlabel("Change in per-phage AUROC when the feature group is removed"); ax.set_title("Feature ablation: shipped vs regenerated inputs", loc="left")
ax.legend(frameon=False, labelcolor=INK, fontsize=8, loc="lower left"); style(ax, "x")
save(fig, "cmp_ablation", a)

# 5. pooled ROC and precision-recall curves: authors' shipped predictions vs this run on shipped and on regenerated inputs
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve
import predict as P

AQUA = "#1baf7a"  # palette slot 3 (slots 1-3 validate all-pairs)
runs = [("Authors' shipped predictions", P.published()[1], ORANGE),
        ("This run, shipped inputs", pd.read_csv(f"{ROOT}/results/repro/predictions.csv"), BLUE),
        ("This run, regenerated inputs", pd.read_csv(f"{ROOT}/results/repro_regenerated/predictions.csv"), AQUA)]
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.8))
curves = []
for name, d, c in runs:
    fpr, tpr, _ = roc_curve(d.y, d.p); pr, rc, _ = precision_recall_curve(d.y, d.p)
    axs[0].plot(fpr, tpr, color=c, lw=2, label=f"{name} (AUROC {roc_auc_score(d.y, d.p):.3f})")
    axs[1].plot(rc, pr, color=c, lw=2, label=f"{name} (AP {average_precision_score(d.y, d.p):.3f})")
    curves += [pd.DataFrame({"run": name, "curve": "ROC", "x": fpr, "y": tpr}), pd.DataFrame({"run": name, "curve": "PR", "x": rc, "y": pr})]
axs[0].plot([0, 1], [0, 1], color=MUTED, lw=1, ls=(0, (4, 3)), label="Chance")
prev = runs[1][1].y.mean(); axs[1].axhline(prev, color=MUTED, lw=1, ls=(0, (4, 3)), label=f"Chance (prevalence {prev:.2f})")
for a, t, xl, yl, loc in [(axs[0], "Pooled ROC, 96 phages", "False positive rate", "True positive rate", "lower right"), (axs[1], "Pooled precision-recall", "Recall", "Precision", "upper right")]:
    a.set_title(t, loc="left"); a.set_xlabel(xl); a.set_ylabel(yl); a.set_xlim(0, 1); a.set_ylim(0, 1.02); style(a, "both"); a.legend(frameon=False, labelcolor=INK, fontsize=8, loc=loc)
fig.tight_layout()
save(fig, "cmp_roc_pr", pd.concat(curves))
