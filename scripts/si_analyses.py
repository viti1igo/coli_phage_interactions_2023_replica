"""Supplementary-text analyses that the authors' repository data are enough for (SI Texts 4, 6, 7 (data for the R fit), 8) plus a test of whether
RBP similarity explains host-range similarity.   python scripts/si_analyses.py  -> results/si/*.csv, results/si/stats.csv
SI Text 4: linear model of the lytic interactions of a phage on genome size and morphotype.
SI Text 6: decision rule 'lytic if and only if the strain carries a K1 ABC capsule' for T145_P2, AL505_P3, MT1B1_P3 and for the other Vectreviruses/Kayfunaviruses.
SI Text 8: PERMANOVA of the bacterial infection profiles on 8 phylogenetic PCoA axes + adsorption factors and/or defense subfamilies, compared by AIC."""
import os, warnings
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu
import pyrodigal
from Bio import Align
from matrix import load_matrix
from treeutil import patristic
import phage_genomes as G

warnings.filterwarnings("ignore")  # pyrodigal warns for every genome shorter than 100 kb
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SET = os.environ.get("MATRIX_SET", "shipped")     # regenerated: the matrix, features, tree and defense subtypes rebuilt from the genomes and raw plaque scores
SUF = "_regenerated" if SET == "regenerated" else ""
OUT = f"{ROOT}/results/si{SUF}"; os.makedirs(OUT, exist_ok=True)
rows = []
def add(item, metric, paper, ours, note="", status=None): rows.append(dict(item=item, metric=metric, paper=paper, ours=ours, note=note, status=status))

M = load_matrix()
if SET == "regenerated":
    reg_m = pd.read_csv(f"{ROOT}/replication/T26_matrix/mlc_matrix_regenerated.csv", sep=";", index_col=0); reg_m.index = reg_m.index.astype(str)
    M = reg_m.loc[reg_m.notna().any(axis=1), M.columns]
col = pd.read_csv(f"{ROOT}/data/genomics/phages/guelin_collection.csv", sep=";", index_col=0)
bact = pd.read_csv(f"{ROOT}/data/genomics/bacteria/picard_collection.csv", sep=";", index_col=0)
LST = pd.read_csv(f"{ROOT}/replication/T14_panacota/annot/LSTINFO-list.lst", sep="\t").assign(strain=lambda d: d.orig_name.str.replace(r"\.fasta$", "", regex=True)).set_index("gembase_name").strain
if SET == "regenerated":
    rg = pd.read_csv(f"{ROOT}/replication/T25_features/picard_collection_regenerated.csv", sep=";", index_col=0, dtype={"ST_Warwick": str}); rg.index = rg.index.astype(str)
    bact = rg.rename(columns={"ABC_model (best system)": "ABC_serotype", "LPS_type (cross-validated markers)": "LPS_type"})
lytic = (M.fillna(0) > 0)
n_lytic = lytic.sum()

# ---- SI Text 4
d4 = pd.DataFrame({"n_lytic": n_lytic, "genome_size": col.loc[n_lytic.index, "Genome_size"], "morphotype": col.loc[n_lytic.index, "Morphotype"]})
m4 = smf.ols("n_lytic ~ genome_size + C(morphotype)", d4).fit()
add("SI Text 4", "genome-size coefficient (lytic interactions per bp), morphotype in the model", 2.3e-4, float(f"{m4.params['genome_size']:.3g}"), f"p = {m4.pvalues['genome_size']:.2f} (paper 0.17); ordinary least squares on 96 phages")
m4b = smf.ols("n_lytic ~ genome_size", d4).fit()
add("SI Text 4", "genome-size coefficient without morphotype", "not given", f"{m4b.params['genome_size']:.2e}, p = {m4b.pvalues['genome_size']:.1e}", "the morphotype confounds genome size", "not comparable (qualitative)")

# ---- SI Text 6
three = ["412_P2", "AL505_Sd2", "MT1B1_3A1"]      # public names of T145_P2, AL505_P3, MT1B1_P3
k1 = (bact.loc[[s for s in M.index if s in bact.index], "ABC_serotype"] == "K1")
def rule(phs):
    sub = M.loc[k1.index, phs]; y = (sub > 0); known = sub.notna()
    pred = pd.DataFrame(np.repeat(k1.values[:, None], len(phs), axis=1), index=sub.index, columns=phs) & known
    tp = int((pred & y).sum().sum()); return tp / max(int(pred.sum().sum()), 1), tp / max(int(y.sum().sum()), 1)
others = [p for p in col.index[col.Genus.isin(["Vectrevirus", "Kayfunavirus"])] if p not in three]
pr3, rc3 = rule(three); pro, rco = rule(others)
add("SI Text 6", "K1 strains in the Picard collection", 45, int(k1.sum()), "ABC_serotype == K1 in picard_collection.csv; the paper may count K1-like calls")
add("SI Text 6", "K1 rule, the three phages: precision (%)", 73, round(pr3 * 100), "pooled over the three phages and all strains with a K1 call")
add("SI Text 6", "K1 rule, the three phages: recall (%)", 72, round(rc3 * 100))
add("SI Text 6", f"K1 rule, the other {len(others)} Vectreviruses and Kayfunaviruses: precision (%)", 7, round(pro * 100))
add("SI Text 6", f"K1 rule, the other {len(others)} Vectreviruses and Kayfunaviruses: recall (%)", 14, round(rco * 100))

# ---- SI Text 7 data: RBP count per phage (the R script si_rbp_lmm.R fits the mixed model)
rbp = pd.read_csv(f"{ROOT}/data/genomics/phages/RBP/RBP_list.csv", sep=";")
cnt = rbp.groupby("phage").RBP.apply(lambda x: int(x.notna().sum())).reindex(M.columns)
pd.DataFrame({"n_rbp": cnt, "n_lytic": n_lytic}).rename_axis("phage").to_csv(f"{OUT}/rbp_count_vs_lytic.csv")
add("SI Text 7", "phages with 0 / 1 / 2 RBP in RBP_list.csv", "all phages have 1 or 2", f"{(cnt == 0).sum()} / {(cnt == 1).sum()} / {(cnt == 2).sum()}", "7 phages are listed with NA: " + ", ".join(cnt.index[cnt == 0]), "differs")

# ---- SI Text 8: PERMANOVA models
tree_names, D = patristic(f"{ROOT}/replication/T14_panacota/tree/ESCO.nucl.grp.aln.fasttree_tree.nwk" if SET == "regenerated" else f"{ROOT}/data/genomics/bacteria/panacota/tree/370+host_tree_root=albertii.nwk")
if SET == "regenerated": tree_names = [LST.get(x, x) for x in tree_names]
strains = [s for s in tree_names if s in M.index and s in bact.index]
ix = [tree_names.index(s) for s in strains]; Dt = D[np.ix_(ix, ix)]
def pcoa(Dm, k=8):
    n = len(Dm); J = np.eye(n) - 1 / n; B = -0.5 * J @ (Dm ** 2) @ J; w, v = np.linalg.eigh(B); o = np.argsort(w)[::-1][:k]; return v[:, o] * np.sqrt(w[o])
phylo = pcoa(Dt)
Y = lytic.loc[strains].values.astype(bool)
Dy = np.nan_to_num(squareform(pdist(Y, "jaccard")))
n = len(strains); J = np.eye(n) - 1 / n; Gm = -0.5 * J @ (Dy ** 2) @ J; SST = np.trace(Gm)
def ind(df, minobs=8):
    out = {}
    for c in df.columns:
        for lev, k in df[c].dropna().astype(str).value_counts().items():
            if k >= minobs and k < len(df): out[f"{c}={lev}"] = (df[c].astype(str) == lev).astype(float).values
    return pd.DataFrame(out, index=df.index)
af_cols = ["O-type", "H-type", "LPS_type", "Klebs_capsule_type", "ABC_serotype"]
omp = pd.read_csv(f"{ROOT}/data/genomics/bacteria/outer_membrane_proteins/blast_results_cured_clusters=99_wide.tsv", sep="\t", index_col=0)
if SET == "regenerated": omp = rg[[c for c in rg.columns if c.startswith("OMP_cluster_")]]
AF = pd.concat([ind(bact.loc[strains, af_cols]), ind(omp.reindex(strains))], axis=1)
dsub = pd.read_csv(f"{ROOT}/data/genomics/bacteria/defense_finder/370+host_defense_systems_cluster_80_80.csv", sep=";", index_col=0).reindex(strains).fillna(0) > 0
if SET == "regenerated":   # our DefenseFinder subtypes (the 80/80 gene clusters are not regenerated); used in place of the subfamilies
    import glob
    ours_d = {LST[os.path.basename(f).replace("_defense_finder_systems.tsv", "")]: (pd.read_csv(f, sep="\t").subtype.value_counts() if len(pd.read_csv(f, sep="\t")) else pd.Series(dtype=int))
              for f in glob.glob(f"{ROOT}/replication/T22_defense/runs/*/*_defense_finder_systems.tsv")}
    dsub = pd.DataFrame(ours_d).T.reindex(strains).fillna(0) > 0
DS = dsub.loc[:, dsub.sum() >= 8].astype(float)
def fit(*blocks):
    X = np.column_stack([np.ones(n), phylo] + [b.values for b in blocks]); H = X @ np.linalg.pinv(X); rk = np.linalg.matrix_rank(X)
    ssm = np.trace(H @ Gm @ H); return ssm, rk
sp, rp = fit()
res = []
for name, blocks in [("phylogeny only", []), ("adsorption factors", [AF]), ("defense subfamilies", [DS]), ("full", [AF, DS])]:
    ssm, rk = fit(*blocks); ssr = SST - ssm
    r2, r2p = ssm / SST, sp / SST
    adj = lambda r, k: 1 - (1 - r) * (n - 1) / (n - k)
    aic = n * np.log(ssr / n) + 2 * rk
    res.append(dict(model=name, covariates=sum(b.shape[1] for b in blocks), rank=rk, R2_total=r2, adjR2_total=adj(r2, rk), R2_after_phylogeny=r2 - r2p, adjR2_after_phylogeny=adj(r2, rk) - adj(r2p, rp), AIC=aic))
res = pd.DataFrame(res); res.to_csv(f"{OUT}/permanova_models.csv", index=False)
variants = []
subt = pd.read_csv(f"{ROOT}/data/genomics/bacteria/defense_finder/370+host_defense_systems_subtypes.csv", sep=";", index_col=0).reindex(strains).fillna(0) > 0
for m in (8, 15, 30):   # minimum number of strains carrying a trait; the paper says traits with fewer than 8 observations were lumped
    af_m = pd.concat([ind(bact.loc[strains, af_cols], m), ind(omp.reindex(strains), m)], axis=1)
    for dname, dtab in (("subfamilies (80/80 clusters)", dsub), ("subtypes", subt)):
        ds_m = dtab.loc[:, dtab.sum() >= m].astype(float)
        for name, blocks in (("adsorption factors", [af_m]), ("defense", [ds_m]), ("full", [af_m, ds_m])):
            ssm, rk = fit(*blocks); r2 = ssm / SST
            variants.append(dict(min_observations=m, defense_unit=dname, model=name, covariates=sum(b.shape[1] for b in blocks), adjR2_total=1 - (1 - r2) * (n - 1) / (n - rk), AIC=n * np.log((SST - ssm) / n) + 2 * rk))
pd.DataFrame(variants).to_csv(f"{OUT}/permanova_variants.csv", index=False)
paper = {"adsorption factors": (0.27, -435), "full": (0.28, -418), "defense subfamilies": (0.05, -372)}
for _, r in res.iterrows():
    if r.model in paper:
        add("SI Text 8", f"PERMANOVA {r.model}: adjusted R2", paper[r.model][0], round(r.adjR2_total, 2), f"total adjusted R2; {r.adjR2_after_phylogeny:.2f} after the phylogeny axes; n = {n} strains, {int(r.covariates)} covariates")
        add("SI Text 8", f"PERMANOVA {r.model}: AIC", paper[r.model][1], round(r.AIC), "AIC = n ln(SS_res / n) + 2 rank")

# ---- does RBP similarity explain host-range similarity?  RBP proteins from our gene calls (the n-th gene of the identifier), clustered by identity
aligner = Align.PairwiseAligner(mode="global", substitution_matrix=Align.substitution_matrices.load("BLOSUM62"), open_gap_score=-10, extend_gap_score=-0.5)
prot = {}
for r in rbp[rbp.RBP.notna()].itertuples():
    ph, idx = r.phage, int(r.RBP.rsplit("_", 1)[1]); seq = G.read_genome(ph).encode()
    f = pyrodigal.GeneFinder(meta=False); f.train(seq, translation_table=11); gs = f.find_genes(seq)
    if idx > len(gs): gs = pyrodigal.GeneFinder(meta=True).find_genes(seq)
    if idx <= len(gs): prot[(ph, r.RBP)] = gs[idx - 1].translate().rstrip("*")
keys = list(prot)
def ident(a, b):
    al = aligner.align(a, b)[0]; same = sum(1 for (s1, e1), (s2, e2) in zip(*al.aligned) for i, j in zip(range(s1, e1), range(s2, e2)) if a[i] == b[j]); return same / al.length
I = np.eye(len(keys))
for i in range(len(keys)):
    for j in range(i + 1, len(keys)):
        a, b = prot[keys[i]], prot[keys[j]]
        if 0.6 < len(a) / len(b) < 1 / 0.6: I[i, j] = I[j, i] = ident(a, b)
pd.DataFrame(I, index=[f"{p}|{g}" for p, g in keys], columns=[f"{p}|{g}" for p, g in keys]).to_csv(f"{OUT}/rbp_pairwise_identity.csv")
phs = sorted({p for p, _ in keys}); pi = {p: [k for k, (q, _) in enumerate(keys) if q == p] for p in phs}
Bm = lytic[phs].values.T.astype(bool); Dh = np.nan_to_num(squareform(pdist(Bm, "jaccard")))   # host-range distance between phages
out = []
for thr in (0.3, 0.5, 0.7):
    share = np.array([[bool(I[np.ix_(pi[a], pi[b])].max() >= thr) if a != b else False for b in phs] for a in phs])
    iu = np.triu_indices(len(phs), 1)
    s, ns = Dh[iu][share[iu]], Dh[iu][~share[iu]]
    genus = np.asarray(col.loc[phs, "Genus"], dtype=object); same_g = (genus[:, None] == genus[None, :])[iu]
    within = share[iu] & same_g; cross = share[iu] & ~same_g
    out.append(dict(identity_threshold=thr, pairs_sharing_RBP=int(share[iu].sum()), pairs_total=len(s) + len(ns), host_range_distance_sharing=s.mean(), host_range_distance_not_sharing=ns.mean(),
                    MWU_p=mannwhitneyu(s, ns).pvalue, pairs_sharing_across_genera=int(cross.sum()), distance_sharing_across_genera=Dh[iu][cross].mean() if cross.any() else np.nan,
                    distance_different_genus_no_sharing=Dh[iu][~share[iu] & ~same_g].mean()))
out = pd.DataFrame(out); out.to_csv(f"{OUT}/rbp_similarity_vs_host_range.csv", index=False)
r5 = out[out.identity_threshold == 0.5].iloc[0]
add("SI Text 6 (generalised)", "mean host-range (Jaccard) distance, phage pairs sharing an RBP (>= 50 % identity) vs not", "RBP sharing explains similar host range (qualitative)", f"{r5.host_range_distance_sharing:.2f} vs {r5.host_range_distance_not_sharing:.2f} (Mann-Whitney p = {r5.MWU_p:.1e})", f"{int(r5.pairs_sharing_RBP)} of {int(r5.pairs_total)} pairs share an RBP; {len(keys)} of {len(rbp[rbp.RBP.notna()])} RBP genes mapped to a called gene", "match")
pd.DataFrame(rows).to_csv(f"{OUT}/stats.csv", index=False)
print(pd.DataFrame(rows).to_string(max_colwidth=75))
print(res.round(3).to_string(index=False)); print(out.round(3).to_string(index=False))
