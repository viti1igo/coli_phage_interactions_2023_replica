"""Phase 9: comparison tables (paper vs authors' shipped results vs our replication), all in results/comparison/ and config/config_diff.csv.
Every paper number carries the source in docs/paper/published.txt (line or section).   python scripts/build_comparison.py"""
import glob, os
import numpy as np, pandas as pd
import predict as P

ROOT = P.ROOT
OUT = f"{ROOT}/results/comparison"
os.makedirs(OUT, exist_ok=True)
R = f"{ROOT}/results"

# ---------------------------------------------------------------- 1. headline numbers: paper | shipped data | replication (shipped inputs) | replication (regenerated inputs)
stats = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{R}/figures/stats_*.csv"))], ignore_index=True)
stats = stats[~stats.figure.str.contains("regenerated defense counts")]  # shown in the regenerated column of the Fig 4c/d rows below
ship_pp, ship_pred = P.published()
pp = {k: pd.read_csv(f"{R}/{d}/per_phage.csv", index_col=0) for k, d in [("shipped_inputs", "repro"), ("regen", "repro_regenerated"), ("regen_matrix", "repro_regenerated_matrix")]}
from sklearn.metrics import roc_auc_score, average_precision_score
pool = lambda d: roc_auc_score(d.y, d.p)
pr = {k: pd.read_csv(f"{R}/{d}/predictions.csv") for k, d in [("shipped_inputs", "repro"), ("regen", "repro_regenerated"), ("regen_matrix", "repro_regenerated_matrix")]}
def strain_auc(d):
    return d.groupby("bacteria").apply(lambda x: roc_auc_score(x.y, x.p) if x.y.nunique() > 1 else np.nan, include_groups=False).mean()


extra = {  # (figure, metric) -> (shipped data result, regenerated inputs result)
    ("Fig 5a", "mean AUROC per strain (published text; Fig 5a plots one AUROC per strain)"): (round(strain_auc(ship_pred), 3), f"{strain_auc(pr['regen']):.3f} (+ our matrix {strain_auc(pr['regen_matrix']):.3f})"),
    ("Fig 5a", "mean AUROC per phage (preprint)"): (round(ship_pp.auroc.mean(), 3), f"{pp['regen'].auroc.mean():.3f} (+ our matrix {pp['regen_matrix'].auroc.mean():.3f})"),
    ("Fig 5", "pooled AUROC"): (round(pool(ship_pred), 3), f"{pool(pr['regen']):.3f} (+ our matrix {pool(pr['regen_matrix']):.3f})"),
    ("Fig 5", "average precision"): (round(average_precision_score(ship_pred.y, ship_pred.p), 3), f"{average_precision_score(pr['regen'].y, pr['regen'].p):.3f}"),
}
fig4 = pd.read_csv(f"{R}/fig4_summary/fig4cd.csv")
for _, r in fig4[fig4.defence_count == "ds_ours"].iterrows():
    extra[(f"Fig 4{'c' if r.model.startswith('C') else 'd'}", "beta (defense systems)")] = ("-", round(r.beta, 3))
stats["published_data_result"] = [extra.get((f, m), ("-", "-"))[0] for f, m in zip(stats.figure, stats.metric)]
stats["replication_regenerated_inputs"] = [extra.get((f, m), ("-", "-"))[1] for f, m in zip(stats.figure, stats.metric)]


def status(paper, ours):
    try:
        p, o = float(paper), float(ours)
    except (TypeError, ValueError):
        return "not comparable (qualitative)"
    if p == 0: return "match" if o == 0 else "differs"
    if 0 < abs(p) < 0.5 and "e-" in f"{paper}".lower():  # p-values: same order of magnitude
        return "match" if abs(np.log10(max(o, 1e-300)) - np.log10(p)) < 1 else "differs"
    return "match (within 5%)" if abs(o - p) / abs(p) <= 0.05 else "differs"


def fmt(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return v
    return int(x) if x == int(x) and abs(x) < 1e6 else float(f"{x:.6g}")


stats["ours"] = stats.ours.map(fmt)
stats["status"] = [status(p, o) for p, o in zip(stats.paper, stats.ours)]
stats["note"] = ""
NOTES = {  # (figure prefix, metric) -> (status override or None, note)
    ("Fig 5a", "mean AUROC per strain (published text; Fig 5a plots one AUROC per strain)"): (None, "the paper's own Source Data for Fig 5a averages 0.872, above the 0.85 in the text; the authors' shipped predictions give 0.881"),
    ("Fig 4c", "beta (defense systems)"): ("consistent (inside paper CI -0.74 to 0.34)", "paper prints 2 decimals"),
    ("Fig 4d", "beta (defense systems)"): ("consistent (inside paper CI -0.05 to -0.01)", "paper prints 2 decimals"),
    ("Fig 4d", "pMCMC"): (None, "both significant; ours smaller than the paper's 0.004"),
    ("Fig 3a", "PERMANOVA adj. R2, isolation strain"): (None, "our own PERMANOVA; the adjusted-R2 formula may differ from the paper's"),
        ("ED9A", "Mann-Whitney P (per-strain mean step, each strain counted three times as in the notebook)"): (None, "the authors' notebook merges the per-strain means onto a table with three rows per strain, so every strain enters the test three times (U = 3856.5); with each strain once the P is 6.4e-5"),
    ("ED6C", "genome size vs average MLC (pMCMC)"): ("match", "significant at 0.01, the threshold of the Fig 4 Source Data; genome size from our assemblies (with the paper's sizes: pMCMC 0.001)"),
    ("ED6D", "AMR genes vs average MLC (pMCMC)"): ("match", "not significant at 0.01 (0.0105, borderline); AMR counts are the paper's Source Data"),
    ("ED6B", "AMR genes vs lytic interactions (pMCMC)"): ("match", "not significant at 0.01 (0.047); AMR counts are the paper's Source Data"),
    ("ED6A", "genome size vs lytic interactions (pMCMC)"): (None, "the paper's title only says genome size correlates with virulence (average MLC)"),
}
for i, r in stats.iterrows():
    for (fig, m), (st, nt) in NOTES.items():
        if r.figure.startswith(fig) and r.metric == m:
            stats.loc[i, "note"] = nt
            if st: stats.loc[i, "status"] = st
rc = pd.read_csv(f"{ROOT}/results/recommender/comparison_with_authors.csv")
sd_extra = pd.DataFrame([
    dict(figure="Fig 6 recommender", metric="identical cocktails out of 100", paper="100 (the authors' recommendations)", ours=int(rc.exact.sum()), paper_source="Fig 6; Methods p.2859; dev/cocktails/data/cocktails_step.csv", note="re-implementation from the Methods; the test strains' ABC and Klebsiella types and UMAP are unavailable (see reports/05)", status="differs"),
    dict(figure="Fig 6 recommender", metric="phages in common per cocktail (of 3)", paper=3, ours=round(rc.overlap.mean(), 2), paper_source="same", note="75 of 100 cocktails share exactly one phage with the authors'", status="differs"),
    dict(figure="Fig 6 recommender", metric="phages recommended by the authors' step 2 that our step 2 also recommends", paper=11, ours=10, paper_source="same", note="the haplotype-neighbour rule is reproduced; the two random-forest steps are not", status="consistent (10 of 11)"),
])
si = pd.read_csv(f"{ROOT}/results/si/stats.csv").rename(columns={"item": "figure"}); si["paper_source"] = si.figure
lm = pd.read_csv(f"{ROOT}/results/si/rbp_lmm.csv")
si = pd.concat([si, pd.DataFrame([dict(figure="SI Text 7", metric=f"LMM beta of the RBP count on lytic interactions, {r.set}", paper=28.11, ours=round(r.beta, 1), paper_source="SI Text 7", note=f"95% CI {r.ci_low:.0f} to {r.ci_high:.0f}, pMCMC {r.pMCMC:.2f} (paper CI -33 to 84, pMCMC 0.33): not significant in both", status="differs") for r in lm.itertuples()])], ignore_index=True)
sd_extra = pd.concat([sd_extra, si[["figure", "metric", "paper", "ours", "paper_source", "note", "status"]]], ignore_index=True)
DROP = [("Fig 2f", "ANOVA F (host range ~ morphotype)"), ("ED7", "coverage of the most-covering triplet (%)"), ("ED2", "modularity Q (BRIM) significance")]
stats = stats[[(f, m) not in DROP for f, m in zip(stats.figure, stats.metric)]].rename(columns={"ours": "replication_shipped_inputs"})
sd = pd.read_csv(f"{OUT}/source_data_checks.csv")  # from scripts/compare_source_data.py (needs the downloaded supplement)
sd["ours"] = sd.ours.map(fmt)
sd = pd.concat([sd, sd_extra], ignore_index=True)
sd["ours"] = sd.ours.map(fmt)
sd["status"] = [st if isinstance(st, str) else status(p, o) for st, p, o in zip(sd.status, sd.paper, sd.ours)]
sd["published_data_result"] = "-"
sd["replication_regenerated_inputs"] = "-"
stats = pd.concat([stats, sd.rename(columns={"ours": "replication_shipped_inputs"})], ignore_index=True)
stats = stats[["figure", "metric", "paper", "published_data_result", "replication_shipped_inputs", "replication_regenerated_inputs", "status", "note", "paper_source"]]
stats.to_csv(f"{OUT}/paper_vs_replication.csv", index=False)

# ---------------------------------------------------------------- 2. per phage
d = pd.DataFrame({"n_lytic": pp["shipped_inputs"].n_pos, "auroc_shipped": ship_pp.auroc, "aupr_shipped": ship_pp.aupr, "best_model_shipped": ship_pp.best_model,
                  "auroc_replication": pp["shipped_inputs"].auroc, "aupr_replication": pp["shipped_inputs"].aupr, "best_model_replication": pp["shipped_inputs"].best_model,
                  "auroc_regenerated": pp["regen"].auroc, "aupr_regenerated": pp["regen"].aupr, "auroc_regenerated_matrix": pp["regen_matrix"].auroc})
d["d_auroc_replication_vs_shipped"] = d.auroc_replication - d.auroc_shipped
d["d_auroc_regenerated_vs_shipped"] = d.auroc_regenerated - d.auroc_shipped
d.rename_axis("phage").to_csv(f"{OUT}/per_phage_auroc.csv")

# ---------------------------------------------------------------- 3. features regenerated vs shipped (Phases 5-7)
f = pd.read_csv(f"{R}/replication/phase6_summary.csv")
f.rename(columns={"value": "agreement_value"}).to_csv(f"{OUT}/feature_agreement.csv", index=False)
mx = pd.read_csv(f"{R}/replication/matrix_agreement.csv").assign(table="interaction matrix (our scoring rules vs shipped)")
tx = pd.read_csv(f"{R}/replication/phage_taxonomy_agreement.csv").assign(table="phage taxonomy (BLAST similarity vs VIRIDIC-derived shipped groups)")
pd.concat([mx, tx]).to_csv(f"{OUT}/matrix_and_taxonomy_agreement.csv", index=False)

# ---------------------------------------------------------------- 4. ablation: shipped inputs vs regenerated inputs
def verdict(dlt, p):
    return "keep" if (dlt <= -0.01 and p < 0.05) else ("weak" if p < 0.05 else "removable")


a1, a2 = [pd.read_csv(f"{R}/{x}/ablation.csv").query("dropped != 'none'").set_index("dropped") for x in ("ablation", "ablation_regenerated")]
ab = pd.DataFrame({"d_auroc_shipped_inputs": a1.d_auroc, "p_shipped_inputs": a1.wilcoxon_p, "d_auroc_regenerated_inputs": a2.d_auroc, "p_regenerated_inputs": a2.wilcoxon_p})
ab["verdict_shipped_inputs"] = [verdict(x, p) for x, p in zip(ab.d_auroc_shipped_inputs, ab.p_shipped_inputs)]
ab["verdict_regenerated_inputs"] = [verdict(x, p) for x, p in zip(ab.d_auroc_regenerated_inputs, ab.p_regenerated_inputs)]
ab["same_verdict"] = ab.verdict_shipped_inputs == ab.verdict_regenerated_inputs
ab["paper_statement"] = "defense systems not required (Fig 4c-d); no feature ablation in the paper"
ab.rename_axis("feature_group_removed").to_csv(f"{OUT}/ablation_comparison.csv")

# ---------------------------------------------------------------- 5. Fig 4 hits by adsorption factor
pd.read_csv(f"{R}/fig4_summary/af_hits_by_family.csv").rename(columns={"Unnamed: 0": "term"}).assign(
    paper="counts only in Fig 4b (image), not printed in the text").to_csv(f"{OUT}/fig4_hits_by_factor.csv", index=False)

# ---------------------------------------------------------------- 6. file inventory
inv = [
    ("paper", "docs/paper/preprint (PDF + txt)", "provided", "bioRxiv version of the article"),
    ("paper", "docs/paper/published (PDF + txt)", "provided", "Nature Microbiology main text, Methods, ED legends, reporting summary"),
    ("paper", "Supplementary Tables 1-6", "NOT provided", "needed for exact per-strain features and per-phage values"),
    ("paper", "Supplementary Text 4-5 and Supplementary Figures 1-9", "NOT provided", "needed for ED figure numbers (ED7 63% definition, ED2 P values)"),
    ("paper", "Source Data (Fig 1-6, ED 1-9)", "NOT provided", "needed to compare every plotted point"),
    ("authors' repo", "data/interactions/interaction_matrix.csv (402 strains x 96 phages)", "downloaded", "paper states 403 strains / 38,688 interactions"),
    ("authors' repo", "dev/predictions (results, logs, importances), dev/cocktails/data, dev/inference_bacteria", "downloaded", "shipped model outputs and cocktail results"),
    ("authors' repo", "CV groups file 370+host_cross_validation_groups_1e-4.csv", "NOT public", "rebuilt from the tree"),
    ("authors' repo", "recommender code (cocktail algorithm)", "NOT in repo", "Fig 6 statistics re-run, algorithm not re-run"),
    ("figshare", "bacterial assemblies (403)", "downloaded", "data/genomics/bacteria/assemblies"),
    ("not obtained", "1,026 GenBank phage genomes, AMR annotation, test-strain genomes", "NOT obtained", "Fig 1k, 3b, 6b, ED3, ED6 B/D"),
]
pd.DataFrame(inv, columns=["origin", "file", "status", "note"]).to_csv(f"{OUT}/file_inventory.csv", index=False)

# ---------------------------------------------------------------- 7. config differences (paper_config.yaml vs what was run)
C = [  # component, parameter, paper, used, reason, expected_effect, paper_source
    ("prediction", "scikit-learn", "1.1.2", "1.9.1", "1.1.2 not installable on Python 3.14", "none expected; L1 penalty arg replaced by l1_ratio=1, sparse solution checked identical", "Methods (statistics)"),
    ("prediction", "random seeds", "unspecified", "RF random_state=0; saga solver random_state=0", "unseeded saga gave tiny run-to-run changes", "reruns identical", "-"),
    ("prediction", "interaction label", "N/P/F/U matrix used by the original script", "lytic = MLC > 0 on the public 0-4 matrix", "private matrix not public", "20.7% positives (README says 19%)", "upstream README"),
    ("prediction", "CV groups", "authors' file (1e-4 substitutions/site)", "rebuilt: single linkage on the shipped tree at 1e-4", "file no longer public", "0.3% of our groups split across the authors' folds", "Methods"),
    ("prediction", "recommender (Fig 6)", "four-step algorithm", "not re-run; statistics re-run on shipped cocktail results", "code not in the repo", "Fig 6a-b only as schematic", "Methods; ED8"),
    ("features", "ST typing tool", "SRST2 (version unspecified)", "mlst 2.33.1, scheme ecoli_achtman_4", "SRST2 needs reads; only assemblies are public", "ST agreement 95.3%", "Methods l.1346-1348"),
    ("features", "phylogroup tool", "ClermonTyping (version unspecified)", "ClermonTyping 24.02 --minimal", "version not stated", "agreement 99.0%", "Methods l.1348"),
    ("features", "ECTyper", "1.0.0", "1.0.0", "-", "raw output 100% identical", "Methods l.1371"),
    ("features", "Kaptive", "2.0.6, default parameters", "2.0.6 run as script from tag v2.0.6 (python 3.8, biopython<1.78, osx-64 under Rosetta); ROAR131 header sanitised", "setup.py incompatible; one assembly malformed", "97.3% agreement", "Methods l.1387-1389"),
    ("features", "DefenseFinder", "1.2.0, models 1.2.3", "1.2.0, models 1.2.3, CasFinder 3.1.0, MacSyFinder 2.1.4", "CasFinder and MacSyFinder versions not stated", "F1 0.976; RM_Type_IV not detected", "Methods l.1372-1373"),
    ("features", "ABC capsule", "CapsuleFinder v1 default ABC model + custom models, inter-gene max space 20", "authors' custom ABC models only, MacSyFinder 2.1.6", "default CapsuleFinder model not run", "ABC presence 98.3%, K-types 89.5%", "Methods l.1373-1386"),
    ("features", "LPS outer-core type", "blastp of waaL + phylogenetic clustering, manual", "pangenome marker families, 5-fold CV against shipped labels", "no per-type waaL reference shipped; K-12 waaL finds only K12/R1/R4", "97.0% accuracy (supervised, not independent)", "upstream README 1.4"),
    ("features", "outer-membrane protein variants", "clusters >99% identity, method unspecified", "BLASTP (Prokka UniProt entries) + mmseqs at 99% identity; YncD from UniProt P76115", "method not stated", "ARI 0.06-0.96 per protein", "Fig 1h legend"),
    ("features", "PanACoTA", "1.3.1 (annotate, pangenome, corepers, align, tree --fasttree)", "1.3.1, Prodigal 2.6.3, FastTree 2.2.0, QC --l90 1000 --nbcont 9999; ESCO.All.prt excluded", "default QC dropped 10 genomes present in the paper's set", "RF distance 0.023; patristic r 0.99999988", "Methods l.1349-1352"),
    ("features", "tree to 8-D UMAP", "umap_dims 8; other parameters unspecified; distance method unspecified", "patristic distance, n_neighbors 15, min_dist 0.1, random_state 0, umap-learn 0.5.12", "parameters not stated; calibrated against shipped coordinates", "10-NN overlap 0.748 (UMAP run-to-run 0.77-0.80)", "upstream README"),
    ("features", "phage taxonomy", "VIRIDIC 1r3.6 (genus 70%, species 95%)", "BLASTn-based similarity approximation, single linkage", "VIRIDIC not installed", "genus ARI 0.929, species ARI 0.241", "Methods l.1353-1358"),
    ("features", "phage network / RBP", "vConTACT2 0.11.3 + Cytoscape 3.9.1; PhageRBPdetect 2.1.3 + HMMER 3.3.2 + ESMFold", "not run", "1,026 GenBank genomes, annotation and structure tools not available", "Fig 1k, 3b, ED3 not replicated; ED4A uses shipped RBP list only", "Methods l.1361-1368, 1390, 1429-1433"),
    ("modelling", "per-phage MCMCglmm", "2.35; 500,000 it, burn-in 50,000, thin 250; prior V=1 nu=0.002", "MCMCglmm 2.36, R 4.3.0; same settings; N(0,10^2) fixed-effect prior for 5 phages that failed as singular", "5 phages gave singular equations", "91 of 96 phages with default prior, 96 of 96 with the extra prior", "Methods l.1460-1475"),
    ("modelling", "Fig 4c-d LMM", "1,000,000 it, burn-in 50,000, thin 250", "same, inverseA from the shipped ultrametric tree", "-", "beta -0.231 vs -0.22; -0.033 vs -0.04", "Methods l.1515-1523"),
    ("modelling", "PERMANOVA", "vegan 2.6-4 adonis2, adjusted R2", "own Python implementation on the Jaccard matrix, 999 permutations, adjusted R2 = 1-(1-R2)(n-1)/(n-p-1)", "R vegan not used for this step", "isolation strain and combined model close; genus 0.27 vs 0.13", "Methods l.1415-1421"),
    ("modelling", "ANOVA / Tukey", "scipy 1.11.2, statsmodels 0.14.1", "scipy 1.18.1 f_oneway; Tukey not run", "-", "F 113.9 vs 114.7", "Methods l.1403-1405"),
    ("modelling", "nestedness / modularity (ED2)", "BINMATNEST temperature T; BRIM Q; 1,000 null matrices", "NODF instead of temperature; BRIM implemented here; 1,000 null matrices (Bernoulli, same connectance)", "no BINMATNEST implementation available", "nestedness number not comparable", "Methods l.1391-1398"),
    ("data", "strains in the interaction matrix", "403 (38,688 interactions)", "402 rows in the public matrix (38,592 cells)", "one strain absent from the public file", "small shifts (e.g. ANOVA F)", "Fig 2 legend"),
    ("data", "interaction matrix used by the figure scripts", "403 strains (Supplementary Table 3)", "the 403-strain matrix when the supplement is downloaded, otherwise the public 402-strain matrix", "the public file lacks strain LF110", "Fig 2f ANOVA F 114.7 with 403 strains, 113.9 with 402; model runs use the public matrix", "Fig 2 legend"),
    ("modelling", "significance threshold of the per-phage models", "pMCMC < 1e-3 (Methods)", "reported at 1e-3 and at 0.01", "the paper's Source Data for Fig 4a-b counts terms up to pMCMC 0.009, i.e. a threshold of 0.01", "28 adsorption-factor terms at 1e-3, 51 at 0.01 (paper: 52)", "Methods; Source Data Fig 4b"),
    ("modelling", "recommender (Fig 6)", "four-step algorithm, code not published", "re-implemented from the Methods (scripts/recommender.py); all features per split; step 2 needs a complete ST/O/H haplotype", "test-strain genomes are not shipped: ABC and Klebsiella types unknown, UMAP recomputed from the shipped distances", "2 of 100 cocktails identical, 0.97 of 3 phages in common", "Methods p.2859; ED Fig. 8"),
    ("data", "AMR gene counts (ED6 B, D)", "annotation tool not stated", "the paper's counts from its Source Data", "no AMR annotation is shipped", "ED6 B and D drawn from the paper's counts", "ED Fig. 6; Source Data ED6"),
    ("figures", "phage annotation (Fig 3b, ED3)", "Prokka 1.14.5 --kingdom Viruses; RBP genes by PhageRBPdetect + curation", "pyrodigal (single mode, metagenomic mode when an RBP identifier exceeds the gene count); RBP genes located by the number in their identifier in RBP_list.csv", "the authors' annotation is not shipped", "gene boundaries and numbering may differ for a few genes (127 of 130 RBP identifiers mapped)", "Methods l.1339-1344, 1429-1433"),
    ("figures", "genome maps and alignment (Fig 3b, ED3)", "Clinker for the maps; a multiple whole-genome alignment (tool not stated)", "own matplotlib maps (reciprocal best BLOSUM62 identity >= 30 %); MAFFT 7.526 --auto; entropy over nucleotides, 200-nt mean", "Clinker and the alignment tool are not available or not named", "layout differs; entropy peaks at the RBP region as in the paper", "Fig 3b legend"),
    ("notebooks", "the authors' four notebooks", "-", "run on a scratch copy, as shipped and with path/compatibility patches (scripts/notebook_patches.py)", "Windows and private paths, a moved file, a missing package, newer pandas and matplotlib", "no analysis changed; one cell of the authors' still fails", "-"),
    ("modelling", "PERMANOVA covariates (SI Text 8)", "adsorption factors and defense subfamilies as presence/absence, traits with fewer than 8 observations lumped, 8 PCoA phylogeny axes first", "indicators of the levels with at least 8 strains (O, H, LPS, Klebsiella, ABC, 12 outer-membrane variants; defense subfamilies at 80/80), Jaccard distance of the lytic profiles, AIC = n ln(SS_res/n) + 2 rank", "covariate lists and the distance are not given", "AF and full adjusted R2 agree (0.26, 0.27 vs 0.27, 0.28); defense 0.17 vs 0.05", "Methods l.1480-1500"),
    ("data", "genome set for tree", "403 strains + references", "402 strains in common with the shipped tree", "ESCO.All.prt is a proteome file, not a genome", "-", "Methods l.1349"),
]
pd.DataFrame(C, columns=["component", "parameter", "paper", "used", "reason", "expected_effect", "paper_source"]).to_csv(f"{ROOT}/config/config_diff.csv", index=False)
print(stats.status.value_counts().to_string()); print(len(C), "config differences;", ab.same_verdict.sum(), "of", len(ab), "ablation verdicts unchanged")
