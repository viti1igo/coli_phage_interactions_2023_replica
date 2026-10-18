"""Fig 4 a-d and Extended Data 1, 2, 4A, 5, 7, 9 (+ genome sizes for ED6, see fig_ed6.R) from shipped data and our results.
python scripts/fig_extended.py -> results/figures/fig4_*.*, ed1.*, ed2.*, ed4a.*, ed5.*, ed7.*, ed9.*, stats_extended.csv, data/derived/genome_size.csv"""
import itertools, os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu
from figlib import *
from matrix import load_matrix
from matrix_stats import nested_modular

D = f"{ROOT}/data"
pic = pd.read_csv(f"{D}/genomics/bacteria/picard_collection.csv", sep=";", index_col=0)
ph = pd.read_csv(f"{D}/genomics/phages/guelin_collection.csv", sep=";", index_col=0)
M = load_matrix()  # 403 strains when the supplement is present (the public matrix lacks LF110)
L = (M.fillna(0) > 0)
leaf = open(f"{D}/genomics/bacteria/panacota/tree/370+host_leaf_order_root=albertii.txt").read().strip().split(",")
tree_rows = [b for b in leaf if b in M.index]


def strips(name, rows, labels, vmax, title, height, data):
    """Heat-strip figure: one row per label over the strains in tree order."""
    fig, ax = plt.subplots(figsize=(13, height))
    ax.imshow(np.vstack(rows), aspect="auto", cmap="Blues", vmin=0, vmax=vmax, interpolation="nearest")
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels); ax.set_xticks([]); ax.set_xlabel("Strains (tree order)"); ax.set_title(title, loc="left"); style(ax, None)
    save(fig, name, data)

# ---------------- Fig 4 a, b (our MCMC re-run vs shipped AF-only) and c, d (scatter + fitted beta)
S = f"{ROOT}/results/fig4_summary"
hits = pd.read_csv(f"{S}/af_hits_by_family.csv", index_col=0)
fam = hits.sort_values("rerun_AF+DS")
fig, axs = plt.subplots(1, 4, figsize=(17, 4.2), gridspec_kw={"width_ratios": [0.7, 1.3, 1, 1]})
a = pd.read_csv(f"{S}/af_vs_ds_hits.csv", index_col=0).iloc[:, 0].reindex(["Adsorption factor", "Defense system"]).fillna(0)
axs[0].bar(a.index, a.values, color=[BLUE if i.startswith("Adsorption") else ORANGE for i in a.index], width=0.6)
for i, v in enumerate(a.values): axs[0].text(i, v + 0.5, int(v), ha="center", fontsize=9, color=INK)
axs[0].set_ylabel("Significant terms (pMCMC < 1e-3)"); axs[0].set_title("a  Term class", loc="left")
y = np.arange(len(fam))
axs[1].barh(y + 0.19, fam["rerun_AF+DS"], height=0.36, color=BLUE, label="Re-run (AF + defense)")
axs[1].barh(y - 0.19, fam["shipped_AF_only"], height=0.36, color=ORANGE, label="Shipped (AF only)")
axs[1].set_yticks(y); axs[1].set_yticklabels(fam.index); axs[1].set_xlabel("Phages with a significant term"); axs[1].set_title("b  Hits per adsorption factor", loc="left"); axs[1].legend(frameon=False, labelcolor=INK, fontsize=8)
cd = pd.read_csv(f"{S}/fig4cd.csv")
d = pd.DataFrame({"bacteria": M.index, "n_lytic": L.sum(axis=1).values,
                  "avg_mlc": [M.loc[b][M.loc[b] > 0].mean() if (M.loc[b] > 0).any() else np.nan for b in M.index],
                  "n_defense": pic.loc[M.index, "n_defense_systems"].values})
for ax, col, lab, key in [(axs[2], "n_lytic", "Lytic interactions of the strain", "C"), (axs[3], "avg_mlc", "Average MLC upon infection", "D")]:
    dd = d.dropna(subset=[col, "n_defense"])
    ax.scatter(dd.n_defense, dd[col], s=14, color=BLUE, alpha=0.55, edgecolor=SURF, lw=0.4)
    b = cd[(cd.model.str.startswith(key)) & (cd.defence_count == "ds_shipped")].iloc[0]
    ax.set_xlabel("Number of defense systems"); ax.set_ylabel(lab)
    ax.set_title(f"{'c' if key == 'C' else 'd'}  beta = {b.beta:.3f}, pMCMC = {b.pMCMC:.4f}  (paper {b.paper_beta}, {b.paper_pMCMC})", loc="left", fontsize=8.5)
for ax in axs: style(ax)
fig.tight_layout()
save(fig, "fig4_a_d", pd.concat([hits.reset_index(), d], axis=1))
for _, r in cd[cd.defence_count == "ds_shipped"].iterrows():
    k = "c" if r.model.startswith("C") else "d"
    stat(f"Fig 4{k}", "beta (defense systems)", r.paper_beta, round(r.beta, 3), "Results; Fig 4 legend")
    stat(f"Fig 4{k}", "pMCMC", r.paper_pMCMC, round(r.pMCMC, 4), "Results; Fig 4 legend")
for _, r in cd[cd.defence_count == "ds_paper4d"].iterrows():  # the defense count of the paper's Source Data for panel d
    stat("Fig 4d (paper's panel-d defense counts)", "beta (defense systems)", r.paper_beta, round(r.beta, 3), "Source Data Fig 4d")
    stat("Fig 4d (paper's panel-d defense counts)", "pMCMC", r.paper_pMCMC, round(r.pMCMC, 4), "Source Data Fig 4d")
for _, r in cd[cd.defence_count == "ds_ours"].iterrows():
    k = "c" if r.model.startswith("C") else "d"
    stat(f"Fig 4{k} (regenerated defense counts)", "beta (defense systems)", r.paper_beta, round(r.beta, 3), "Results; Fig 4 legend")

# ---------------- ED1: BE and BL21 (MG1655 was tested at a single concentration and is not in the shipped matrix)
lab = M.loc[["BE", "BL21"]]
fig, ax = plt.subplots(figsize=(13, 3))
order = ph.sort_values(["Morphotype", "Genus"]).index
im = ax.imshow(lab[order].values, aspect="auto", cmap="Blues", vmin=0, vmax=4, interpolation="nearest")
ax.set_yticks([0, 1]); ax.set_yticklabels(lab.index); ax.set_xticks([]); ax.set_xlabel("96 phages (morphotype, genus)"); ax.set_title("ED1  Laboratory-adapted strains (MLC 0-4; MG1655 not in the shipped matrix)", loc="left")
fig.colorbar(im, ax=ax, pad=0.01, label="MLC"); style(ax, None)
save(fig, "ed1", lab[order].T.rename_axis("phage").reset_index())
for s, paper in [("BE", "45 (46.9%)"), ("BL21", "48 (50.0%)")]:
    stat("ED1", f"phages infecting {s}", paper, f"{int((M.loc[s] > 0).sum())} ({(M.loc[s] > 0).mean() * 100:.1f}%)", "ED Fig. 1 legend")

# ---------------- ED2: nestedness (NODF as stand-in for BINMATNEST temperature) and BRIM modularity, 1000 null matrices
A = L.values.astype(float)
nm = nested_modular(A)
T0, nulls, p_nodf, cbest, Q0, qn, p_q = (nm[k] for k in ("T0", "nulls", "p_nodf", "cbest", "Q0", "qn", "p_q"))
print(f"NODF {T0:.2f} (null {nulls.mean():.2f} +- {nulls.std():.2f}, p {p_nodf:.4f}); Q {Q0:.3f} at c={cbest} (null {qn.mean():.3f} +- {qn.std():.3f}, p {p_q:.4f})")
ro, co = np.argsort(-A.sum(1)), np.argsort(-A.sum(0))
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1.3, 1, 1]})
axs[0].imshow(A[np.ix_(ro, co)].T, aspect="auto", cmap="Blues", interpolation="nearest"); axs[0].set_title("A  Matrix sorted by degree (nestedness)", loc="left"); axs[0].set_xticks([]); axs[0].set_yticks([]); axs[0].set_xlabel("Strains"); axs[0].set_ylabel("Phages")
axs[1].hist(nulls, bins=30, color="#bcd3f2"); axs[1].axvline(T0, color=ORANGE, lw=2); axs[1].set_xlabel("NODF"); axs[1].set_ylabel("Null matrices"); axs[1].set_title(f"NODF {T0:.1f}, null {nulls.mean():.1f}, P = {p_nodf:.3f}", loc="left", fontsize=9)
axs[2].hist(qn, bins=30, color="#bcd3f2"); axs[2].axvline(Q0, color=ORANGE, lw=2); axs[2].set_xlabel("Modularity Q (BRIM)"); axs[2].set_title(f"B  Q {Q0:.3f}, null {qn.mean():.3f}, P = {p_q:.3f}", loc="left", fontsize=9)
for a in axs[1:]: style(a)
style(axs[0], None); fig.tight_layout()
save(fig, "ed2", pd.DataFrame({"null_nodf": nulls, "null_Q": qn, "NODF": T0, "Q": Q0, "c_best": cbest}))
stat("ED2", "nestedness (paper: BINMATNEST temperature T; ours: NODF, not comparable)", "nested pattern, P not printed", f"NODF {T0:.1f}, P = {p_nodf:.3f}", "Results p.2850; ED Fig. 2")
stat("ED2", "modularity Q (BRIM) significance", "weakly modular, P not printed", f"Q = {Q0:.3f} ({cbest} modules), P = {p_q:.3f}", "ED Fig. 2 legend")

# The paper uses new phage names (supplement); the public matrix keeps the old ones. Mapping inferred from genus + isolation strain
# (PDP110_P1 = 411_P1 is confirmed by its 3 infected strains, as in ED5).
NAME = {"AL505_P3": "AL505_Sd2", "T145_P2": "412_P2", "MT1B1_P3": "MT1B1_3A1", "PDP110_P1": "411_P1"}

# ---------------- ED4 A: host range of three phages with similar RBP; K1 ABC capsule marked
trio = [NAME[t] for t in ("AL505_P3", "T145_P2", "MT1B1_P3")]
k1 = (pic.loc[tree_rows, "ABC_serotype"] == "K1").values
strips("ed4a", [k1.astype(int) * 4, *M.loc[tree_rows, trio].fillna(0).values.T], ["ABC K1"] + trio, 4, "ED4A  Host range of three phages sharing a K1-endosialidase RBP", 2.6,
       pd.DataFrame({"bacteria": tree_rows, "K1": k1, **{t: M.loc[tree_rows, t].values for t in trio}}))
for t in trio:
    inf = M.loc[tree_rows, t].fillna(0).values > 0
    stat("ED4A", f"{t}: infected strains that are ABC K1", "all K1-associated (qualitative)", f"{(inf & k1).sum()} of {inf.sum()} infected; {(inf & k1).sum()} of {k1.sum()} K1 strains infected", "ED Fig. 4A legend")

# ---------------- ED5: PDP110_P1 and Klebsiella K55 capsule
x = M[NAME["PDP110_P1"]].fillna(0) > 0
k55 = pic.loc[M.index, "Klebs_capsule_type"].fillna("").str.contains("K55")
strips("ed5", [k55.loc[tree_rows].astype(int) * 2, x.loc[tree_rows].astype(int)], ["Klebsiella K55 locus", "PDP110_P1 (411_P1) lysis"], 2, "ED5A  PDP110_P1 infects exactly the K55 strains", 2.2,
       pd.DataFrame({"bacteria": M.index, "K55": k55.values, "PDP110_P1_lytic": x.values}))
stat("ED5", "strains infected by PDP110_P1", 3, int(x.sum()), "ED Fig. 5 legend")
stat("ED5", "infected strains that carry K55", "3 of 3", f"{int((x & k55).sum())} of {int(x.sum())}; K55 strains total {int(k55.sum())}", "ED Fig. 5 legend")

# ---------------- ED7: most-covering k-uplet of phages (exact for k <= 3, greedy above)
cols = L.values; n = cols.shape[0]; names = list(L.columns)
def cover(idx): return cols[:, idx].any(axis=1).sum() / n * 100
best3 = max(itertools.combinations(range(cols.shape[1]), 3), key=cover)
cov = {1: max(cover([i]) for i in range(cols.shape[1]))}
cov[2] = max(cover(c) for c in itertools.combinations(range(cols.shape[1]), 2)); cov[3] = cover(best3)
chosen = list(best3)
for k in range(4, 26):
    nxt = max((i for i in range(cols.shape[1]) if i not in chosen), key=lambda i: cover(chosen + [i])); chosen.append(nxt); cov[k] = cover(chosen)
fig, ax = plt.subplots(figsize=(6, 3.8))
ax.plot(list(cov), list(cov.values()), color=BLUE, lw=2, marker="o", ms=4); ax.axhline(cov[3], color=MUTED, lw=1, ls=(0, (4, 3)))
ax.text(14, cov[3] - 4, f"best triplet {cov[3]:.1f}%", color=INK, fontsize=8); ax.set_xlabel("k (phages in the cocktail)"); ax.set_ylabel("Strains lysed by at least one (%)"); ax.set_title("ED7  Coverage of the most-covering k-uplet", loc="left"); ax.set_ylim(0, 100); style(ax)
save(fig, "ed7", pd.DataFrame({"k": list(cov), "coverage_pct": list(cov.values())}))
stat("ED7", "coverage of the most-covering triplet (%)", 63, round(cov[3], 1), "ED Fig. 7 legend; Methods p.2855")
stat("ED7", "baseline triplet", "not named", " + ".join(names[i] for i in best3), "-")

# ---------------- ED9: step at which each phage is recommended, haplotypes per cocktail
C = f"{ROOT}/dev/cocktails/data"
st = pd.read_csv(f"{C}/cocktails_step.csv", sep=";"); comp = pd.read_csv(f"{C}/cocktails_composition.csv", sep=";")
GENERIC = ["55989_P2-DIJ07_P2-LF82_P8", "536_P7-LF82_P8-LF82_P9", "55989_P2-LF82_P9-DIJ07_P2", "536_P7-DIJ07_P2-LF82_P8"]
kind = comp.set_index("bacteria").cocktails_str.isin(GENERIC).map({True: "Generic", False: "Tailored"})
st["kind"] = st.bacteria.map(kind); st["step"] = st.chosen_at_step.str[0].astype(int)
hap = st.drop_duplicates("bacteria").set_index("bacteria")
hap = (hap.ST_Warwick.astype(str) + "|" + hap["O-type"].astype(str) + "|" + hap["H-type"].astype(str))
comp["kind"] = comp.bacteria.map(kind); comp["haplotype"] = comp.bacteria.map(hap)
nh = comp.groupby(["cocktails_str", "kind"]).haplotype.nunique().rename("n_haplotypes").reset_index()
# the authors' notebook tests the per-strain mean of the step number (1 same host, 2 neighbours, 3 trait matching, 4 rescue), tailored vs generic strains
smean = st.assign(n=st.chosen_at_step.map({"0_same_as_host": 1, "1_neighbours": 2, "3_misc": 3, "4_rescue": 4})).groupby("bacteria").agg(n=("n", "mean"), kind=("kind", "first"))
u = mannwhitneyu(smean[smean.kind == "Tailored"].n, smean[smean.kind == "Generic"].n, alternative="two-sided")
# the notebook then merges the per-strain means back onto the three-rows-per-strain table, so each strain enters the test three times
u3 = mannwhitneyu(np.repeat(smean[smean.kind == "Tailored"].n.values, 3), np.repeat(smean[smean.kind == "Generic"].n.values, 3), alternative="two-sided")
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
ct = pd.crosstab(st.kind, st.step, normalize="index") * 100
w = 0.35
for i, (k, c) in enumerate([("Tailored", BLUE), ("Generic", ORANGE)]):
    axs[0].bar(np.arange(len(ct.columns)) + (i - 0.5) * w, ct.loc[k].values, width=w, color=c, label=k)
axs[0].set_xticks(range(len(ct.columns))); axs[0].set_xticklabels([f"step {s}" for s in ct.columns]); axs[0].set_ylabel("Recommended phages (%)"); axs[0].set_title(f"A  Step of recommendation (MWU P = {u.pvalue:.1e})", loc="left", fontsize=9); axs[0].legend(frameon=False, labelcolor=INK)
nh = nh.sort_values("n_haplotypes", ascending=False)
axs[1].bar(range(len(nh)), nh.n_haplotypes, color=[ORANGE if k == "Generic" else BLUE for k in nh.kind], width=0.75); axs[1].set_xticks([]); axs[1].set_xlabel("18 cocktails"); axs[1].set_ylabel("Haplotypes (ST/O/H) covered"); axs[1].yaxis.set_major_locator(plt.MaxNLocator(integer=True)); axs[1].set_title("B  Haplotypes per cocktail", loc="left")
axs[1].legend(handles=[Patch(color=ORANGE, label="Generic"), Patch(color=BLUE, label="Tailored")], frameon=False, labelcolor=INK)
for a in axs: style(a)
fig.tight_layout()
save(fig, "ed9", nh)
stat("ED9A", "Mann-Whitney P (per-strain mean step, each strain counted three times as in the notebook)", "3.5e-12", f"{u3.pvalue:.1e}", "ED Fig. 9 legend")
stat("ED9A", "Mann-Whitney P (per-strain mean step, each strain once)", "not given", f"{u.pvalue:.1e}", "-")
stat("ED9", "cocktails", 18, len(nh), "Fig. 6 legend")

# ---------------- genome sizes for ED6
out = f"{ROOT}/data/derived/genome_size.csv"
if not os.path.exists(out):
    A_ = f"{D}/genomics/bacteria/assemblies"; rows = []
    for f in sorted(os.listdir(A_)):
        if f.startswith("._") or not f.endswith(".fasta"): continue
        rows.append((f[:-6], sum(len(l.strip()) for l in open(f"{A_}/{f}") if not l.startswith(">"))))
    pd.DataFrame(rows, columns=["bacteria", "genome_size"]).to_csv(out, index=False)

# ---------------- ED6 A-D (fits from fig_ed6.R: genome size from our assemblies, AMR gene counts from the paper's Source Data)
e6 = f"{ROOT}/results/fig4_summary/ed6.csv"
if os.path.exists(e6):
    fit = pd.read_csv(e6).set_index("panel"); dd = pd.read_csv(f"{ROOT}/results/fig4_summary/ed6_data.csv")
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.5))
    for ax, pn, col, x, lab, xl, unit in [(axs[0, 0], "A", "n_lytic", "mb", "Lytic interactions of the strain", "Genome size (Mb)", "Mb"), (axs[0, 1], "B", "n_lytic", "amr", "Lytic interactions of the strain", "AMR genes (paper's counts)", "gene"),
                                          (axs[1, 0], "C", "avg_mlc", "mb", "Average MLC upon infection", "Genome size (Mb)", "Mb"), (axs[1, 1], "D", "avg_mlc", "amr", "Average MLC upon infection", "AMR genes (paper's counts)", "gene")]:
        z = dd.dropna(subset=[col, x]); ax.scatter(z[x], z[col], s=14, color=BLUE, alpha=0.55, edgecolor=SURF, lw=0.4)
        f = fit.loc[pn]; ax.set_xlabel(xl); ax.set_ylabel(lab); style(ax)
        ax.set_title(f"{pn}  beta = {f.beta:.3g} per {unit}, pMCMC = {f.pMCMC:.3f}", loc="left", fontsize=9)
    fig.tight_layout(); save(fig, "ed6_a_d", dd)
    for pn, what in [("A", "genome size vs lytic interactions"), ("C", "genome size vs average MLC"), ("B", "AMR genes vs lytic interactions"), ("D", "AMR genes vs average MLC")]:
        f = fit.loc[pn]
        stat(f"ED6{pn}", f"{what} (pMCMC)", "significant" if what.startswith("genome size vs average") else ("not significant" if what.startswith("AMR") else "not stated"), f"beta {f.beta:.3g}, pMCMC {f.pMCMC:.4f}", "ED Fig. 6 title")
flush_stats("extended")
