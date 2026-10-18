"""Fig 3a (+PERMANOVA), Fig 5 a/c/d, Fig 6 c-e from shipped/our results.
python scripts/fig_main_3_6.py -> results/figures/fig3_a.*, fig5_a_c_d.*, fig6_c_e.*, stats_fig3_6.csv"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.stats import mannwhitneyu
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import f1_score, average_precision_score, roc_auc_score
from figlib import *
from matrix import load_matrix
from matrix_stats import phage_distance, phage_models
import supp
from supp import sheet

D = f"{ROOT}/data"
ph = pd.read_csv(f"{D}/genomics/phages/guelin_collection.csv", sep=";", index_col=0)
M = load_matrix()  # 403 strains when the supplement is present (the public matrix lacks LF110)

# ---------------- Fig 3a: Jaccard distance on binary lytic vectors (NaN = no lysis), UMAP, PERMANOVA
Dm = phage_distance(M)
import umap
emb = umap.UMAP(n_neighbors=15, min_dist=0.1, metric="precomputed", random_state=0).fit_transform(Dm)  # paper: parameters not stated


meta = ph.loc[M.columns, ["Genus", "Phage_host"]].copy()
meta["Genus"] = meta.Genus.astype(str); meta["Host"] = meta.Phage_host.astype(str)
res = phage_models(Dm, meta)
for k, (r2, adj, p) in res.items():
    print(f"PERMANOVA {k}: R2 {r2:.3f} adj {adj:.3f} p {p:.4f}")
inc = res["genus x isolation strain"][1] - res["isolation strain"][1]  # genus added after isolation strain (sequential terms, as adonis2 reports them)
stat("Fig 3a", "PERMANOVA adj. R2, genus added after isolation strain", 0.13, round(inc, 2), "Results p.2851 (inferred: the 0.13 equals the sequential increment of genus after isolation strain)")
stat("Fig 3a", "PERMANOVA adj. R2, isolation strain", 0.53, round(res["isolation strain"][1], 2), "Results p.2851")
stat("Fig 3a", "PERMANOVA adj. R2, genus + isolation strain", 0.65, round(res["genus x isolation strain"][1], 2), "Results p.2851")
stat("Fig 3a", "PERMANOVA adj. R2, genus alone", "not stated", round(res["genus"][1], 2), "-")
stat("Fig 3a", "PERMANOVA R2 (unadjusted), genus", "not stated", round(res["genus"][0], 3), "-")
stat("Fig 3a", "PERMANOVA R2 (unadjusted), isolation strain", "not stated", round(res["isolation strain"][0], 3), "-")

top = meta.Genus.value_counts().index[:8]
cmap = plt.get_cmap("tab10")
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
for g in meta.Genus.unique():
    m = (meta.Genus == g).values
    axs[0].scatter(emb[m, 0], emb[m, 1], s=22, color=cmap(list(top).index(g)) if g in top else "#b8b7b0", edgecolor=SURF, lw=0.5, label=g if g in top else None)
axs[0].legend(frameon=False, fontsize=7, labelcolor=INK, title="Genus (top 8; others grey)", title_fontsize=7)
tophost = meta.Host.value_counts().index[:8]
for g in meta.Host.unique():
    m = (meta.Host == g).values
    axs[1].scatter(emb[m, 0], emb[m, 1], s=22, color=cmap(list(tophost).index(g)) if g in tophost else "#b8b7b0", edgecolor=SURF, lw=0.5, label=g if g in tophost else None)
axs[1].legend(frameon=False, fontsize=7, labelcolor=INK, title="Isolation strain (top 8; others grey)", title_fontsize=7)
for a, t in zip(axs, ["a  Phage host-range UMAP, by genus", "a  by isolation strain"]):
    a.set_title(t, loc="left"); a.set_xlabel("UMAP 1"); a.set_ylabel("UMAP 2"); style(a, None)
fig.tight_layout()
save(fig, "fig3_a", pd.DataFrame({"phage": M.columns, "umap1": emb[:, 0], "umap2": emb[:, 1], "genus": meta.Genus.values, "isolation_strain": meta.Host.values}))

# ---------------- Fig 5 a, c, d
pp = pd.read_csv(f"{ROOT}/results/repro/per_phage.csv", index_col=0)
pr = pd.read_csv(f"{ROOT}/results/repro/predictions.csv")
pr["yhat"] = (pr.p >= 0.5).astype(int)
strain_auc = pr.groupby("bacteria").apply(lambda d: roc_auc_score(d.y, d.p) if d.y.nunique() > 1 else np.nan, include_groups=False).dropna()
n_lytic = pr.groupby("bacteria").y.sum().reindex(strain_auc.index)
paper5 = sheet(9, "Figure_5c", text_first=False) if supp.available() else None   # Source Data: per-strain AUROC against lytic interactions
fig, axs = plt.subplots(1, 4, figsize=(18, 4.2), gridspec_kw={"width_ratios": [1, 1, 1, 1.5]})
axs[0].hist(strain_auc, bins=np.arange(0.3, 1.0001, 0.05), color=BLUE, edgecolor=SURF)
axs[0].axvline(strain_auc.mean(), color=INK, lw=1, ls=(0, (4, 3))); axs[0].text(strain_auc.mean() - 0.01, axs[0].get_ylim()[1] * 0.93, f"mean {strain_auc.mean():.3f}", fontsize=8, color=INK, ha="right")
axs[0].set_xlabel("AUROC per strain"); axs[0].set_ylabel("Strains"); axs[0].set_title("a  AUROC per strain", loc="left")
axs[1].hist(pp.auroc, bins=np.arange(0.3, 1.0001, 0.05), color=BLUE, edgecolor=SURF)
axs[1].axvline(pp.auroc.mean(), color=INK, lw=1, ls=(0, (4, 3))); axs[1].text(pp.auroc.mean() + 0.01, axs[1].get_ylim()[1] * 0.93, f"mean {pp.auroc.mean():.3f}", fontsize=8, color=INK)
axs[1].set_xlabel("AUROC per phage"); axs[1].set_ylabel("Phages"); axs[1].set_title("(preprint's panel) AUROC per phage", loc="left", fontsize=9)
if paper5 is not None:
    axs[2].scatter(paper5.n_infections, paper5.auc, s=22, facecolor="none", edgecolor=ORANGE, lw=0.9, label="Paper (Source Data)")
axs[2].scatter(n_lytic, strain_auc, s=16, color=BLUE, edgecolor=SURF, lw=0.4, alpha=0.8, label="This run")
axs[2].set_xlabel("Lytic interactions of the strain"); axs[2].set_ylabel("AUROC of the strain"); axs[2].set_title("c  AUROC vs lytic interactions", loc="left"); axs[2].legend(frameon=False, labelcolor=INK, fontsize=8, loc="lower right")
for a in axs[:3]: style(a)
# d: prediction matrix, strains in tree order x phages by host range; categories TP, FP, FN, TN/NA
leaf = open(f"{D}/genomics/bacteria/panacota/tree/370+host_leaf_order_root=albertii.txt").read().strip().split(",")
rows = [b for b in leaf if b in M.index]
cols = M.columns[np.argsort(-(M > 0).sum().values)]
cat = pd.DataFrame(0, index=rows, columns=cols)
for r in pr.itertuples():
    if r.bacteria in cat.index:
        cat.loc[r.bacteria, r.phage] = 1 if (r.y and r.yhat) else 2 if (not r.y and r.yhat) else 3 if (r.y and not r.yhat) else 0
axs[3].imshow(cat.values.T, aspect="auto", interpolation="nearest", cmap=ListedColormap(["#f6f5f1", BLUE, "#e9b8a3", ORANGE]), vmin=0, vmax=3)
axs[3].set_xticks([]); axs[3].set_yticks([]); axs[3].set_title("d  Predictions (blue TP, light orange FP, orange FN)", loc="left", fontsize=9); axs[3].set_xlabel("Strains (tree order)"); axs[3].set_ylabel("Phages (by host range)")
style(axs[3], None)
fig.tight_layout()
save(fig, "fig5_a_c_d", pd.DataFrame({"n_lytic": n_lytic, "auroc_strain": strain_auc}).rename_axis("bacteria").reset_index())
pp.reset_index().to_csv(f"{OUT}/fig5_per_phage_auroc.csv", index=False)
cat.reset_index().rename(columns={"index": "bacteria"}).to_csv(f"{OUT}/fig5_d_prediction_matrix.csv", index=False)
tp, fp, fn = [(cat.values == k).sum() for k in (1, 2, 3)]
prec, rec = tp / (tp + fp), tp / (tp + fn)
stat("Fig 5a", "mean AUROC per strain (published text; Fig 5a plots one AUROC per strain)", 0.85, round(strain_auc.mean(), 3), "Results p.2854 (published.txt l.623); Fig 5a legend")
stat("Fig 5a", "mean AUROC per phage (preprint)", 0.77, round(pp.auroc.mean(), 3), "Results (preprint.txt l.521)")
stat("Fig 5", "pooled AUROC", 0.86, round(roc_auc_score(pr.y, pr.p), 3), "Fig 5b")
stat("Fig 5", "average precision", 0.60, round(average_precision_score(pr.y, pr.p), 3), "Fig 5b")
stat("Fig 5d", "F1 at p >= 0.5 (all pairs pooled)", 0.60, round(f1_score(pr.y, pr.yhat), 3), "Results")
stat("Fig 5d", "TP / FP / FN (pairs)", "not stated", f"{tp} / {fp} / {fn}", "-")
stat("Fig 5c", "Spearman(AUROC, n lytic) per strain", "not stated", round(strain_auc.corr(n_lytic, method="spearman"), 3), "-")

# ---------------- Fig 6 c-e
C = f"{ROOT}/dev/cocktails/data"
cml = pd.read_csv(f"{C}/cml_cocktails.csv", sep=";")
comp = pd.read_csv(f"{C}/cocktails_composition.csv", sep=";")
GENERIC = ["55989_P2-DIJ07_P2-LF82_P8", "536_P7-LF82_P8-LF82_P9", "55989_P2-LF82_P9-DIJ07_P2", "536_P7-DIJ07_P2-LF82_P8"]
cml = cml.merge(comp[["bacteria", "cocktails_str"]], on="bacteria")
cml["kind"] = np.where(cml.cocktails_str.isin(GENERIC), "Generic", "Tailored")
nc = cml.groupby(["cocktails_str", "kind"]).size().rename("strains").reset_index().sort_values("strains", ascending=False)
fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={"width_ratios": [1.5, 1, 1.2]})
axs[0].bar(range(len(nc)), nc.strains, color=[ORANGE if k == "Generic" else BLUE for k in nc.kind], width=0.75)
axs[0].set_xticks([]); axs[0].set_xlabel(f"{len(nc)} distinct cocktails"); axs[0].set_ylabel("Strains for which it was recommended"); axs[0].set_title("c  Recommendations per cocktail", loc="left")
from matplotlib.patches import Patch
axs[0].legend(handles=[Patch(color=ORANGE, label="Generic"), Patch(color=BLUE, label="Tailored")], frameon=False, labelcolor=INK)
sets = {"Tailored": cml[cml.kind == "Tailored"]["p1+p2+p3"], "Generic": cml[cml.kind == "Generic"]["p1+p2+p3"], "Baseline": cml["b1+b2+b3"]}
cols3 = {"Tailored": BLUE, "Generic": ORANGE, "Baseline": "#8a8983"}
succ = {k: (v > 0).mean() * 100 for k, v in sets.items()}
axs[1].bar(list(succ), list(succ.values()), color=[cols3[k] for k in succ], width=0.6)
for i, (k, v) in enumerate(succ.items()): axs[1].text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=9, color=INK)
axs[1].set_ylim(0, 105); axs[1].set_ylabel("Strains lysed (%)"); axs[1].set_title("d  In vitro success rate", loc="left")
bins = [-0.5, 0.5, 2.5, 4.5, 6.5]; names = ["0", "1-2", "3-4", "5-6"]
w = 0.26
for i, (k, v) in enumerate(sets.items()):
    h = np.histogram(v, bins=bins)[0] / len(v) * 100
    axs[2].bar(np.arange(4) + (i - 1) * w, h, width=w, color=cols3[k], label=k)
axs[2].set_xticks(range(4)); axs[2].set_xticklabels(names); axs[2].set_xlabel("MLC of the 3-phage cocktail"); axs[2].set_ylabel("Strains (%)"); axs[2].legend(frameon=False, labelcolor=INK); axs[2].set_title("e  Binned MLC", loc="left")
for a in axs: style(a)
fig.tight_layout()
save(fig, "fig6_c_e", cml)
mw = lambda a, b: mannwhitneyu(a, b, alternative="two-sided")
t = cml[cml.kind == "Tailored"]
r1, r2 = mw(sets["Tailored"], t["b1+b2+b3"]), mw(sets["Tailored"], sets["Generic"])
stat("Fig 6d", "success rate, tailored (%)", 91.67, round(succ["Tailored"], 2), "Results p.2855")
stat("Fig 6d", "success rate, generic (%)", 78.95, round(succ["Generic"], 2), "Results p.2855")
stat("Fig 6d", "success rate, baseline (%)", 81, round(succ["Baseline"], 1), "Results p.2855")
stat("Fig 6e", "Mann-Whitney U, tailored vs baseline (tailored strains)", 399.5, r1.statistic, "Fig 6 legend")
stat("Fig 6e", "P, tailored vs baseline", 0.019, round(r1.pvalue, 3), "Fig 6 legend")
stat("Fig 6e", "Mann-Whitney U, tailored vs generic", 1171.5, r2.statistic, "Fig 6 legend")
stat("Fig 6e", "P, tailored vs generic", 0.034, round(r2.pvalue, 3), "Fig 6 legend")
flush_stats("fig3_6")
