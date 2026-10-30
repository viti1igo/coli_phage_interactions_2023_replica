"""Fig 1 b-j (collection descriptors) and Fig 2 c-f (interaction matrix, host-range histograms, morphotype ANOVA) from shipped tables.
python scripts/fig_main_1_2.py  ->  results/figures/fig1_b_j.*, fig2_c.*, fig2_d_e.*, fig2_f.*, stats_fig1_2.csv"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.stats import f_oneway
from figlib import *
from matrix import load_matrix

D = f"{ROOT}/data"
pic = pd.read_csv(f"{D}/genomics/bacteria/picard_collection.csv", sep=";", index_col=0)
ph = pd.read_csv(f"{D}/genomics/phages/guelin_collection.csv", sep=";", index_col=0)
M = load_matrix()  # 403 strains when the supplement is present (the public matrix lacks LF110)
omp = pd.read_csv(f"{D}/genomics/bacteria/outer_membrane_proteins/blast_results_cured_clusters=99_wide.tsv", sep="\t", index_col=0)
dsub = pd.read_csv(f"{D}/genomics/bacteria/defense_finder/370+host_defense_systems_subtypes.csv", sep=";", index_col=0)

# ---------------- Fig 1 b-j
def barh(ax, s, title, top=None, color=BLUE):
    s = s.value_counts().sort_values(ascending=False)
    if top:
        s = s.iloc[:top]
    ax.barh(s.index.astype(str)[::-1], s.values[::-1], color=color, height=0.7)
    ax.set_title(title, loc="left"); ax.set_xlabel("Strains"); style(ax, "x")
    return s

rows = []
fig, axs = plt.subplots(3, 3, figsize=(13, 11))
ax = axs.ravel()
for a, (lab, col, title, top) in zip(ax[:6], [
        ("b", "ST_Warwick", "b  Sequence type (top 20)", 20), ("c", "O-type", "c  O-antigen (top 20)", 20),
        ("d", "Pathotype", "d  Pathotype", None), ("e", "Origin", "e  Source", None),
        ("f", "LPS_type", "f  LPS outer-core type", None), ("g", "ABC_serotype", "g  ABC capsule serotype", None)]):
    s = barh(a, pic[col].fillna("none"), title, top)
    rows += [dict(panel=lab, category=str(k), n=int(v)) for k, v in s.items()]
# h: number of variants (clusters at 99% identity) of each outer-membrane protein
nv = omp.nunique().sort_values()
ax[6].barh(nv.index, nv.values, color=BLUE, height=0.7); ax[6].set_title("h  Variants per outer-membrane protein", loc="left"); ax[6].set_xlabel("Clusters (99% identity)"); style(ax[6], "x")
rows += [dict(panel="h", category=k, n=int(v)) for k, v in nv.items()]
# i: defense subtypes per genome
ds = dsub.drop(columns=[c for c in dsub.columns if dsub[c].dtype == object], errors="ignore")
per = (ds > 0).sum(axis=1)
vc = per.value_counts().sort_index()
ax[7].bar(vc.index, vc.values, color=BLUE, width=0.8); ax[7].set_title("i  Defense subtypes per genome", loc="left"); ax[7].set_xlabel("Subtypes"); ax[7].set_ylabel("Strains"); style(ax[7])
rows += [dict(panel="i", category=int(k), n=int(v)) for k, v in vc.items()]
# j: frequency of defense subtypes (top 25)
fr = (ds > 0).sum().sort_values(ascending=False).iloc[:25]
ax[8].barh(fr.index[::-1], fr.values[::-1], color=BLUE, height=0.7); ax[8].set_title("j  Most frequent defense subtypes (top 25)", loc="left"); ax[8].set_xlabel("Strains"); style(ax[8], "x")
rows += [dict(panel="j", category=k, n=int(v)) for k, v in fr.items()]
fig.tight_layout()
save(fig, "fig1_b_j", pd.DataFrame(rows))
stat("Fig 1", "strains in the Picard collection", 403, len(pic), "Fig. 1 / Results, collection description")
stat("Fig 1j", "defense subtypes detected (columns)", "not stated numerically", ds.shape[1], "Fig. 1j")

# ---------------- Fig 2 c: matrix ordered by phylogeny; phages by genus then host range
leaf = open(f"{D}/genomics/bacteria/panacota/tree/370+host_leaf_order_root=albertii.txt").read().strip().split(",")
order = [b for b in leaf if b in M.index] + [b for b in M.index if b not in leaf]
hr = (M > 0).sum() / M.notna().sum()
pord = ph.assign(hr=hr).sort_values(["Morphotype", "Genus", "hr"], ascending=[True, True, False]).index
X = M.loc[order, pord]
fig, ax = plt.subplots(figsize=(11, 7))
Z = X.values.copy(); Z = np.where(np.isnan(Z), -1, Z) + 1
ax.imshow(Z.T, aspect="auto", cmap=ListedColormap(["#ecebe6", "#f6f5f1", "#bcd3f2", "#8bb4ea", "#2a78d6", "#154f9b"]), vmin=0, vmax=5, interpolation="nearest")
ax.set_xlabel("Strains (core-genome tree order)"); ax.set_ylabel("Phages (morphotype, genus)"); ax.set_title("Fig 2c  Interaction matrix (MLC 0-4; grey = not tested)", loc="left")
ax.set_xticks([]); ax.set_yticks([]); style(ax, None)
save(fig, "fig2_c", X.reset_index().rename(columns={"index": "bacteria"}))

# ---------------- Fig 2 d, e
pct_strain = (M > 0).sum(axis=1) / M.notna().sum(axis=1) * 100
pct_phage = hr * 100
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
axs[0].hist(pct_strain, bins=20, color=BLUE, edgecolor=SURF); axs[0].set_xlabel("% of phages infecting the strain"); axs[0].set_ylabel("Strains"); axs[0].set_title("d  Strains", loc="left")
axs[1].hist(pct_phage, bins=20, color=BLUE, edgecolor=SURF); axs[1].set_xlabel("% of strains infected by the phage"); axs[1].set_ylabel("Phages"); axs[1].set_title("e  Phages", loc="left")
for a in axs: style(a)
fig.tight_layout()
save(fig, "fig2_d_e", pd.concat([pct_strain.rename("pct_phages_infecting_strain").reset_index().rename(columns={"index": "bacteria"}),
                                  pct_phage.rename("pct_strains_infected").reset_index().rename(columns={"index": "phage"})], axis=1))
stat("Fig 2", "mean % of strains infected per phage", 20.6, round(pct_phage.mean(), 1), "Results (published.txt l.99-100: 83 +/- 70 strains, 20.6 % of the collection)")

# ---------------- Fig 2 f: host range by morphotype, one-way ANOVA
g = {m: pct_phage[ph.loc[pct_phage.index, "Morphotype"] == m] for m in ["Myoviridae", "Podoviridae", "Siphoviridae"]}
F, P = f_oneway(*g.values())
fig, ax = plt.subplots(figsize=(4.8, 4))
bp = ax.boxplot(list(g.values()), tick_labels=list(g), widths=0.5, patch_artist=True, medianprops=dict(color=INK), flierprops=dict(markersize=3))
for b in bp["boxes"]: b.set(facecolor="#bcd3f2", edgecolor=BLUE)
ax.set_ylabel("% of strains infected"); ax.set_title(f"f  Host range by morphotype  (ANOVA F = {F:.1f}, P = {P:.1e})", loc="left", fontsize=9); style(ax)
save(fig, "fig2_f", pd.DataFrame({"morphotype": [m for m, v in g.items() for _ in v], "pct_strains_infected": np.concatenate([v.values for v in g.values()])}))
stat("Fig 2f", "ANOVA F (host range ~ morphotype)", 114.7, round(F, 1), "Results p.2850; Fig 2f legend")
stat("Fig 2f", "ANOVA P", "1.44e-26", f"{P:.2e}", "Fig 2f legend")
for m, paper in [("Myoviridae", 41.3), ("Podoviridae", 11.1), ("Siphoviridae", 9.1)]:
    stat("Fig 2f", f"mean % lytic interactions, {m}", paper, round(g[m].mean(), 1), "Results p.2850")
flush_stats("fig1_2")
