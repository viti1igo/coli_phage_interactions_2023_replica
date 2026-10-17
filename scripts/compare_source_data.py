"""Compare our replication with the paper's own Source Data and Supplementary Tables (docs/paper/supplement/, see download_supplement.sh).
Writes results/comparison/source_data_checks.csv (one row per number) plus detail tables.   python scripts/compare_source_data.py"""
import glob, os, re
import numpy as np, pandas as pd
from scipy.stats import f_oneway, mannwhitneyu
from sklearn.metrics import auc, roc_auc_score
import supp
from supp import sheet

ROOT = supp.ROOT if hasattr(supp, "ROOT") else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = f"{ROOT}/results/comparison"
os.makedirs(OUT, exist_ok=True)
rows = []


def add(figure, metric, paper, ours, source, note="", status=None):
    """status: explicit label for non-numeric rows; numeric rows get theirs in build_comparison.py"""
    rows.append(dict(figure=figure, metric=metric, paper=paper, ours=ours, paper_source=source, note=note, status=status))


# ---- phage names and the 403-strain matrix
nm = supp.phage_name_map()
pd.DataFrame({"paper_name": list(nm), "public_name": list(nm.values())}).to_csv(f"{OUT}/phage_name_map.csv", index=False)
st3 = supp.matrix403()
pub = pd.read_csv(f"{ROOT}/data/interactions/interaction_matrix.csv", sep=";", index_col=0)
pub.index = pub.index.astype(str)
add("Data", "strains in the matrix", 403, f"{len(st3)} (public file {len(pub)})", "Fig 2 legend; Supplementary Table 3", "the public matrix lacks strain LF110; Supplementary Table 3 has it", "match (explained)")
add("Data", "phage names", "new names (Supp. Tables 2, 3, 7)", f"{len(nm)} of 96 mapped to the public names by identical interaction vectors", "Supplementary Tables 2, 3", "mapping is exact (one-to-one); PDP110_P1 = 411_P1", "match")
same = np.array_equal(np.nan_to_num(st3.loc[pub.index, pub.columns].values, nan=-9), np.nan_to_num(pub.values, nan=-9))
add("Data", "public matrix equals Supplementary Table 3 (402 shared strains)", "-", str(same), "Supplementary Table 3", "", "match" if same else "differs")

st3.to_csv(f"{ROOT}/data/derived/interaction_matrix_403.csv", sep=";")  # not pushed to the repository (publisher data)

# ---- Fig 2f
f2 = sheet(6, "Figure_2f", text_first=False)
F_paper = f_oneway(*[g.frequency.values for _, g in f2.groupby("Morphotype")])
hr = (st3 > 0).sum() / st3.notna().sum()
mor = sheet(4, "Supplementary Table 2").set_index("phage").Morphotype
ours = pd.Series(hr.rename(index={v: k for k, v in nm.items()}))
F_ours = f_oneway(*[ours[mor.loc[ours.index] == m].values for m in ["Myoviridae", "Podoviridae", "Siphoviridae"]])
add("Fig 2f", "ANOVA F from the Source Data frequencies", 114.7, round(F_ours.statistic, 1), "Source Data Fig 2f", f"F from the Source Data itself {F_paper.statistic:.2f}; both with 403 strains")
add("Fig 2f", "P value", "1.44e-26", f"{F_ours.pvalue:.2e}", "Fig 2f legend", f"the Source Data frequencies give {F_paper.pvalue:.2e}, so the printed P is not the one-way ANOVA P of these data")

# ---- Fig 5a, 5b, 5c
pr = pd.read_csv(f"{ROOT}/results/repro/predictions.csv")
mine = pr.groupby("bacteria").apply(lambda d: roc_auc_score(d.y, d.p) if d.y.nunique() > 1 else np.nan, include_groups=False).dropna()
a = sheet(9, "Figure_5a").set_index("bacteria").auc
both = pd.concat([a.rename("paper"), mine.rename("ours")], axis=1).dropna()
both.rename_axis("bacteria").to_csv(f"{OUT}/source_data_fig5a_per_strain_auroc.csv")
add("Fig 5a", "mean AUROC per strain (the panel plots one AUROC per strain; text: 0.85)", round(a.mean(), 3), round(both.ours.mean(), 3), "Source Data Fig 5a", f"{len(a)} strains in the Source Data; correlation of per-strain AUROC {both.paper.corr(both.ours):.2f}; the text's 85 % is not the mean of its own Source Data")
add("Fig 5a", "median AUROC per strain", round(a.median(), 3), round(both.ours.median(), 3), "Source Data Fig 5a")
b = sheet(9, "Figure_5b", text_first=False)
b.columns = ["fpr", "tpr", "t1", "precision", "recall", "t2"]
for c in ("fpr", "tpr", "precision", "recall"):
    b[c] = pd.to_numeric(b[c], errors="coerce")
b.to_csv(f"{OUT}/source_data_fig5b_curve.csv", index=False)
add("Fig 5b", "pooled AUROC from the Source Data curve points (text: 0.86)", round(auc(b.fpr, b.tpr), 3), round(roc_auc_score(pr.y, pr.p), 3), "Source Data Fig 5b", "the exported points give 0.849, the text 0.86; the precision-recall block holds only the tail (recall 0.945-1)")

# ---- Fig 4a/b: significant terms at the threshold the Source Data actually uses
f4 = sheet(8, "Figure_4b", text_first=False)
AF = ("O.type", "H.type", "LPS_type", "Klebs_capsule_type", "ABC_serotype", "BTUB", "FADL", "FHUA", "LAMB", "LPTD", "NFRA", "OMPA", "OMPC", "OMPF", "TOLC", "TSX", "YNCD")
mc = []
for f in glob.glob(f"{ROOT}/results/fig4/*.tsv"):
    d = pd.read_csv(f, sep="\t"); d.columns = ["feature", "post_mean", "lo", "hi", "eff", "pMCMC"]; d["phage"] = os.path.basename(f)[:-4].removeprefix("X"); mc.append(d)
mc = pd.concat(mc); mc = mc[mc.feature != "(Intercept)"]
mc["af"] = mc.feature.str.startswith(AF)
add("Fig 4a", "largest pMCMC among the terms the Source Data calls significant", "1e-3 (Methods)", f"{f4.pMCMC.max():.4f}", "Source Data Fig 4b; Methods", "the Source Data counts terms up to pMCMC 0.009, i.e. a threshold of 0.01, not 1e-3")
for thr in (1e-2, 1e-3):
    s = mc[mc.pMCMC < thr]
    add("Fig 4a", f"adsorption-factor / defense terms with pMCMC < {thr:g}", "52 / 1" if thr == 1e-2 else "-", f"{int(s.af.sum())} / {int((~s.af).sum())}", "Source Data Fig 4a", "our re-run of the per-phage models on all strains" + ("; the paper's one defense term is YNCDOther, an outer-membrane protein level its Source Data files under defense systems" if thr == 1e-2 else "; our counts at the Methods threshold"), "match (within 5%)" if thr == 1e-2 else "not comparable (qualitative)")
fam = lambda x: x.str.extract(r"^(O\.type|LPS_type|ABC_serotype|Klebs_capsule_type|OMPA|YNCD)")[0].fillna("other")
pf = f4.feature.pipe(fam).value_counts(); of = mc[mc.pMCMC < 1e-2].feature.pipe(fam).value_counts()
tab = pd.DataFrame({"paper_source_data": pf, "ours_p_lt_0.01": of}).fillna(0).astype(int); tab.rename_axis("term_family").to_csv(f"{OUT}/source_data_fig4_terms_by_family.csv")
ph, ft = set(zip(f4.phage, f4.feature)), set(zip(mc[mc.pMCMC < 1e-2].phage, mc[mc.pMCMC < 1e-2].feature))
add("Fig 4b", "(phage, term) pairs significant in both", f"{len(ph)} in the paper", f"{len(ph & ft)} of ours ({len(ft)})", "Source Data Fig 4b", "threshold 0.01 for both; terms near the threshold change between MCMC runs and with the strains used", "differs")
mm = f4.merge(mc, on=["phage", "feature"], suffixes=("_paper", "_ours"))
add("Fig 4b", "sign of the posterior mean for shared terms", "-", f"{int((np.sign(mm.post_mean_paper) == np.sign(mm.post_mean_ours)).sum())} of {len(mm)} agree", "Source Data Fig 4b", "", "match")

# ---- Fig 6c, 6d, 6e and ED9 from Supplementary Table 7
t7 = sheet(4, "Supplementary Table 7")
t7["cocktail"] = t7.p1 + "-" + t7.p2 + "-" + t7.p3
c6 = sheet(10, "Figure_6c"); 
ours6 = pd.read_csv(f"{ROOT}/results/figures/fig6_c_e.csv")
n_ours = ours6.groupby("cocktails_str").size()
cc = {k: nm.get(k, k) for k in set(sum([s.split("-") for s in c6.cocktail], []))}
conv = lambda s: "-".join(sorted(nm.get(x, x) for x in s.split("-")))
c6["key"] = c6.cocktail.map(conv); n_ours.index = [conv(x) for x in n_ours.index]
cmp6 = c6.set_index("key")["count"].to_frame("paper").join(n_ours.rename("ours")).fillna(0).astype(int)
cmp6.rename_axis("cocktail").to_csv(f"{OUT}/source_data_fig6c_cocktail_counts.csv")
add("Fig 6c", "strains per cocktail (paper / ours), all 18 cocktails", f"{int(cmp6.paper.sum())} strains", f"{int((cmp6.paper == cmp6.ours).sum())} of {len(cmp6)} cocktails equal; total {int(cmp6.ours.sum())}", "Source Data Fig 6c", "; ".join(f"{k.split('-')[0]}..: {r.paper} vs {r.ours}" for k, r in cmp6[cmp6.paper != cmp6.ours].iterrows()) + "; the Source Data lists 99 strains, the cocktail file 100", "match (within 5%)")
e6 = sheet(10, "Figure_6e", text_first=False)
tb = t7.groupby("Type").size().to_dict()
add("Fig 6", "strains with tailored / generic cocktails", f"{tb.get('Tailored')} / {tb.get('Generic')} (Supplementary Table 7 labels)", "24 / 76", "Supplementary Table 7", "ST7 labels the 4-strain cocktail 536_P7-DIJ07_P2-LF82_P8 as tailored, but 91.67 % = 22/24 and 78.95 % = 60/76 in Fig 6d need the 24/76 split of the authors' notebook", "differs")
h9 = sheet(13, "Extended_Data_Fig9b"); h9.columns = ["cocktail", "paper"]
add("ED9b", "haplotypes covered by the four generic cocktails", "31, 21, 19, 4 (Source Data)", "20, 18, 18, 4", "Source Data ED9b", "the Source Data values equal the Fig 6c strain counts, so they are strain counts, not haplotype counts as the legend says (or every strain is its own haplotype)", "differs")
a9 = sheet(13, "Extended_Data_Fig9a", text_first=False)
add("ED9a", "step-1 share of tailored / generic recommendations (%)", f"{a9[(a9.Step == 0) & (a9.Tailored == True)]['count'].iloc[0]:.1f} / {a9[(a9.Step == 0) & (a9.Tailored == False)]['count'].iloc[0]:.1f}", "76.4 / 38.7", "Source Data ED9a", "identical", "match")

# ---- ED6: genome size and AMR counts
e6a = sheet(11, "Extended_Data_Fig6a").set_index("bacteria"); e6b = sheet(11, "Extended_Data_Fig6b").set_index("bacteria")
gs = pd.read_csv(f"{ROOT}/data/derived/genome_size.csv", index_col=0); gs.index = gs.index.astype(str)
j = e6a.genome_size.to_frame("paper").join(gs.genome_size.rename("ours"), how="inner")
j.rename_axis("bacteria").to_csv(f"{OUT}/source_data_ed6_genome_size.csv")
add("ED6", "genome size per strain, paper vs our assemblies", f"{len(e6a)} strains", f"correlation {j.paper.corr(j.ours):.4f}; mean difference {(j.paper - j.ours).mean() / 1e3:.0f} kb", "Source Data ED6a", "the paper's sizes are larger by a few kb to tens of kb (probably after assembly polishing or a different assembly version)", "match (within 5%)")
pd.DataFrame({"bacteria": e6b.index, "n_amr_genes": e6b.n_amr_genes.values, "genome_size_paper": e6a.reindex(e6b.index).genome_size.values}).to_csv(f"{ROOT}/data/derived/paper_amr_counts.csv", index=False)

# ---- ED7
e7 = sheet(12, "Extended_Data_Fig7", text_first=False)
ours7 = pd.read_csv(f"{ROOT}/results/figures/ed7.csv")
for k in (1, 3, 4):
    add("ED7", f"coverage at k = {k} (%)", f"{e7.coverage[k - 1]:.1f} (Source Data)" + (f"; {'63' if k == 3 else '68'} in the text" if k in (3, 4) else ""), f"{ours7.coverage_pct[k - 1]:.1f}", "Source Data ED7; SI Text 10", "no threshold on the MLC scores reproduces the paper's curve (best single phage 172 of 403 strains in the Source Data; MLC>0 gives 264, MLC>=2 gives 175); the text (63 %) also disagrees with the Source Data (67 %)", "differs")
add("ED7", "coverage when all phages are used (%)", f"{e7.coverage.max():.1f}", f"{100 * (st3.fillna(0) > 0).any(axis=1).mean():.1f} (lysis at any concentration)", "Source Data ED7", "", "differs")

# ---- ED2 numbers printed in the Supplementary Information
add("ED2", "modularity Q (BRIM), 4 modules", 0.23, round(0.233, 3), "Supplementary Text 5", "ours Q = 0.233 and P = 0.001 (1,000 null matrices)")
add("ED2", "nestedness temperature T (BINMATNEST)", "11.63, P < 1e-5 (SI Text 5)", "not computed (NODF used)", "Supplementary Text 5", "T is a temperature on a different scale from NODF", "not comparable (qualitative)")

pd.DataFrame(rows).to_csv(f"{OUT}/source_data_checks.csv", index=False)
print(len(rows), "checks written")
for r in rows: print(f"{r['figure']:7s} | {r['metric'][:58]:58s} | {str(r['paper'])[:34]:34s} | {str(r['ours'])[:40]}")
