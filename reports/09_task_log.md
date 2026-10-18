# 9. Task log

The task list of the replication with the outcome of each step, as kept during the work. File names that no longer exist (`run_lps.py`, `fig6_cocktails.py`) were removed in the code audit ([08](08_code_audit.md)).

## Tasks
## Phase 1 – Data
- [x] T1 Clone repo into working dir; confirm tables join (402 strains in matrix/features/tree/UMAP)
- [x] T2 Read paper: what was already tested (defense systems out; AF+phylogeny sufficient)
- [x] T3 `scripts/download_data.sh` (repo + extract .rar + optional figshare genomes w/ md5)
## Phase 2 – Reproduce
- [x] T4 `scripts/make_groups.py` → `data/derived/cv_groups_1e-4.csv`; validate vs original fold logs
- [x] T5 `scripts/predict.py` → `results/repro/`; compare to published per-phage perf
### Checkpoint: pooled AUROC ≈ 0.86, per-phage mean ≈ 0.77 (or deviation explained)
## Phase 3 – Feature necessity
- [x] T6 `scripts/ablate.py`: drop-one-group + published-importance aggregation → `results/ablation/`
- [x] T7 Verdict table: keep / redundant per feature, with "already in paper?" column
## Phase 4 – Docs & audit
- [x] T8 Markdown: `README_REPRO.md`, `docs/FEATURES.md`, `docs/RESULTS.md`
- [x] T9 /ponytail:ponytail-audit on scripts; apply fixes
- [x] T10 Remote Control on

## Added after request to run remaining parts
- [x] Fig 6 statistics (`scripts/fig6_cocktails.py`)
- [~] Fig 4 MCMC: 91/96 phages; 5 failed (singular), Fig 4C/D not run
- [ ] Optional: re-fit the 5 failed phages (stronger fixed-effect prior or per-phage min-frequency filter) as a labelled deviation

---
# Extension (2026-10-16): replicate from preprocessing (IN PROGRESS; approved 2026-10-17)
_This file is updated after every step (user request). Latest update: see the Status log at the bottom._

## Phase 5 - Feasibility and config
- [x] T11 Environment spike: all 9 tools run (`results/env_check.csv`). Timing: PanACoTA annotate 18 s per 5 genomes with Prodigal vs ~10 min per 5 with Prokka (Rosetta); ECTyper 7 s/genome; mlst 4 s; ClermonTyping 6 s; Kaptive ~2 min; DefenseFinder ~104 s.
- [x] T12 Download inputs: 403 genome assemblies (figshare, md5-checked); phage genomes already in repo; plaque images (9.9 GB) skipped unless approved. Done when: 403 files, md5 match.
- [x] T13 Pin the paper's tool versions and parameters in `config/paper_config.yaml`, marking unspecified ones. Done when: every item cites a paper line.
### Checkpoint 5: tools run, inputs present, config pinned. Review before long runs.

## Phase 6 - Bacterial features (each compared with the shipped output)
- [x] T14 [DONE: our independent core-genome tree (PanACoTA 1.3.1 + FastTree, Prodigal annotation, relaxed QC; 1,295 core families, 403 genomes) vs the shipped tree on 402 shared strains: 390 of 399 splits shared (normalised Robinson-Foulds 0.023), pairwise patristic distances Pearson 0.99999988 and Spearman 0.99999, median absolute difference 4.4e-6, total tree length 0.8683 vs 0.8683; CV groups at <1e-4: 300 vs 301 groups, adjusted Rand index 0.984. Deviations: Prodigal instead of Prokka, QC thresholds relaxed (paper does not state them), core = families present in all 403 genomes (default)] Core-genome tree. Files: `replication/T14_panacota/tree/`, `results/replication/tree_agreement.csv`.
- [x] T15 [DONE: UMAP (8-D, n_neighbors 15, min_dist 0.1, settings not stated in the paper) of our tree vs shipped coordinates: 10-nearest-neighbour overlap 0.748, the same as the run-to-run overlap of UMAP itself (0.77-0.80); embedding in `replication/T14_panacota/umap_8_dims_our_tree.tsv`] Tree to distance matrix to 8-D UMAP.
- [x] T16 [DONE: 384/403 = 95.3% agree; mlst scheme ecoli_achtman_4; mismatches mostly novel STs now numbered in the current DB] Sequence type (assembly-based substitute for SRST2). Done when: % agreement with `ST_Warwick`.
- [x] T17 [DONE: 399/403 = 99.0% agree] Phylogroup (ClermonTyping 24.02). Done when: % agreement with `Clermont_Phylo`.
- [x] T18 [DONE: raw ECTyper output equals the authors' raw output for 402/402 strains (O and H 100%); vs curated `picard_collection.csv` columns: O 74.7%, H 91.8% because the authors curated multi-calls and "-" by hand] O and H serotype (ECTyper 1.0.0).
- [x] T19 [DONE, two parts] (1) K-12 WaaL (UniProt P27243; protein BLAST and tblastn on genomes give the same result, as in the authors' README "K12 waaL as a query") recognises only type K12 (97.6% identity) and weakly R1/R4 (~30%); R2/R3 waaL are not homologous to it (no hit at e-value 1), so presence agreement is only 78.7% and per-type references are not shipped. (2) Substitute: marker families from the PanACoTA pangenome (`scripts/run_lps_markers.py`), 5-fold cross-validated: 97.0% (390/402); R3, R4, No_waaL 100%, R1 98.1%, K12 95.3%, R2 72.7%. Caveat: supervised by the shipped labels, so it shows the pangenome carries the information, not an independent rediscovery. LPS outer core (blastp *waaL* plus clustering; the paper's manual step is not reproducible).
- [x] T20 [DONE: 403/403 genomes. ABC capsule present: 396/403 = 98.3% agree. Named K-types in both: 85/95 = 89.5% (K1 42/42, K15 6/6, K7 6/6, K96 4/4, K5 19/19 incl. K5_like, K4 10/10; K2 only 5/10). The shipped HP, CatB, Unknown and numeric labels map one-to-one onto our class_1_x models (CatB=class_1_5 11/11, Unknown=class_1_7 7/7, '1'=class_1_15 4/4, '3'=class_1_4 5/5, HP=class_1_19 11/13). Mismatches: 3 shipped-none strains we call K1; 2 HP strains with no system. CapsuleFinder v1 presence step not run (only the custom typing models)] ABC capsule. Files: `results/replication/capsule_crosstab.csv`.
- [x] T21 [DONE: 403/403 genomes (Kaptive 2.0.6). K-type vs `Klebs_capsule_type`: 392/403 = 97.3% agree; 20 of the 22 shipped calls reproduced (missing: NILS22 K10, PDP110 K55); we call 9 more confident types (K51 x3, K123, K16 x2, K31, K57 x2) that the authors' manually cured table does not list. Data note: the figshare ROAR131 assembly has sequence glued to 2 header lines, which makes makeblastdb (inside Kaptive) abort and Kaptive skip the genomes after it in its list; we used a sanitized copy for Kaptive only (+120 bases recovered, 182 contigs unchanged) and the script now retries missing genomes] Klebsiella-like capsule (Kaptive 2.0.6).
- [x] T22 [DONE: 403/403 genomes, DefenseFinder 1.2.0 + models 1.2.3 + CasFinder 3.1.0; 137 subtypes on both sides; presence F1 0.976 (precision 0.981, recall 0.970); systems per strain Pearson 0.972 (mean 9.03 ours vs 9.17 shipped); exact subtype-count vector 58.6%. Systematic differences: RM_Type_IV not called by our run (-123 systems), slightly more AbiE (+15), Thoeris_II (+14), RM_Type_II (+7), AbiJ (+4); cause not established (possible: model or CasFinder version, hmmsearch version, or the authors' own filtering). Input proteins differ too (Prodigal vs Prokka)] Defense systems, for Fig 4 only. Files: `results/replication/defense_agreement.csv`, `defense_subtype_differences.csv`.
- [x] T23 [DONE, 12/12 genes; clusters at 99% identity of the gene nucleotides; adjusted Rand index vs shipped partition: ompC 0.96, yncD 0.92, lptD 0.90, ompF 0.88, lamB 0.86, tsx 0.81, ompA 0.74, btuB 0.67, fhuA 0.59, nfrA 0.40, fadL 0.34, tolC 0.06; identity sweep: mean ARI is best at 0.99 (0.655) but some genes fit other thresholds (tolC 0.81 at 0.97, fadL 0.95 at 0.95), so the authors' exact clustering is not recovered; yncD uses UniProt P76115 (K-12)] Outer-membrane-protein variants (method unspecified), for Fig 4 only. Files: `results/replication/omp_agreement.csv`.
- [x] T24 [DONE: regenerated host table from our T16-T19 results vs `isolation_strain_features.csv`: 33 strains; Clermont 33/33, ST 33/33, LPS 33/33 (cross-validated predictions), O-type 32/33 and H-type 32/33 where the one difference is only how "no call" is written (`-` vs empty). The authors' notebook only merges tool outputs; its paths are Windows-only] Host (isolation-strain) feature table: `scripts/build_host_table.py` -> `results/replication/host_table_agreement.csv`.
- [x] T25 [DONE: `scripts/assemble_features.py` -> `replication/T25_features/picard_collection_regenerated.csv` (403 strains x 19 columns: O/H, ST, phylogroup, Klebsiella K-type, ABC model, LPS (cross-validated), 12 OMP cluster ids, n_defense_systems) and the 35-metric diff `results/replication/phase6_summary.csv`] Assemble the regenerated feature table.
### Checkpoint 6: REACHED 2026-10-17 - feature table regenerated; agreement per feature known (`results/replication/phase6_summary.csv`). Review before modelling.

## Phase 7 - Phage and interaction data
- [x] T26 [DONE: `scripts/build_matrix.py`. Rule = lowest phage concentration with lysis (1 at 5e8, 2 at 5e7, 3-4 at 5e6), replicates combined by the median. 35,268 comparable cells (369 strains x 96 phages; the 33 host strains are not in the raw table): lytic vs not agrees in 99.8%; exact on {0,1,2,"3 or 4"} in 96.6% (94.9% with a mean, 89.3% with a max). Limits: the raw table is binary, so plaque (3) vs clearing (4) cannot be told apart (shipped 3 and 4 both show as lysis at all three concentrations); the largest disagreement is 867 cells shipped 2 where we give 1] Plaque scores to MLC matrix. Files: `replication/T26_matrix/mlc_matrix_regenerated.csv`, `results/replication/matrix_agreement.csv`, `matrix_confusion.csv`.
- [x] T27 [DONE as an approximation: VIRIDIC is not installed, so BLASTn all-vs-all with a VIRIDIC-like similarity (`scripts/run_phage_taxonomy.py`); 97 phage genomes; genus (>=70%): 26 clusters vs 20 shipped genera, adjusted Rand index 0.929; species (>=95%): 82 vs 46, ARI 0.241 (our measure under-counts identity for partially similar genomes; the real VIRIDIC is needed for species)] Phage genus and species. Files: `results/replication/phage_taxonomy_agreement.csv`.
- [-] T28 [NOT DONE: needs PhageRBPdetect (trained language-model classifier), Pfam-A 34 for HMMER (large download) and ESMFold structure prediction, plus the authors' manual curation; ESMFold is not practical on a 16 GB Mac without a GPU. The shipped `RBP_list.csv` stays the reference] RBP detection, optional.

## Phase 8 - Re-run models on regenerated inputs
- [x] T29 [DONE: `FEATURE_SET=... scripts/predict.py` (results in `results/repro*`). Same models and cross-validation on (a) shipped inputs: mean AUROC per phage 0.787, AUPR 0.591, pooled AUROC 0.856; (b) regenerated features + tree-based UMAP and CV groups (labels shipped): 0.765, 0.569, pooled 0.855; (c) also our MLC matrix where raw data exist: 0.769, 0.565, pooled 0.853. Published: 0.781, 0.586, 0.857. So the pooled AUROC survives regeneration (within 0.004), the per-phage mean loses about 2 points. Two code fixes: a phage whose lytic strains all fall in one test fold now gets a constant prediction instead of an error]  CV groups and Fig 5 predictions with regenerated features.
- [x] T30 [DONE: 96/96 phages. The 5 phages that failed (409_P3, 411_P1, 411_P2, 412_P1, 423_P10: "mixed model equations singular") were refitted by `scripts/fig4_mcmc.R` with an N(0, 10^2) prior on the fixed effects (a labelled deviation; listed in `results/fig4/STRONG_PRIOR_PHAGES.txt`). Fig 4A: 28 significant adsorption-factor terms, 0 defence-system terms (pMCMC < 1e-3); by family LPS 14, ABC capsule 12, O-type 1, OMPA 1. Against the shipped AF-only results (7,776 shared phage-term pairs): both significant 16, shipped only 10, re-run only 12, Spearman of pMCMC 0.625. Fig 4C/D reproduced (`results/fig4_summary/fig4cd.csv`): C beta -0.231 [-0.761, 0.284], pMCMC 0.39 vs paper -0.22 [-0.74, 0.34], 0.41; D beta -0.033 [-0.049, -0.014], pMCMC 0.0005 vs paper -0.04 [-0.05, -0.01], 0.004; with our own defence counts C -0.238 (p 0.40), D -0.037 (p 0.0003)] Fig 4 MCMC and Fig 4C/D.
- [x] T31 [DONE: `FEATURE_SET=regenerated scripts/ablate.py` -> `results/ablation_regenerated/`; baseline AUROC 0.765 / AUPR 0.569. Verdicts vs the shipped-input ablation: ABC capsule still necessary (dAUROC -0.019, dAUPR -0.056, p < 1e-4) and LPS outer core (-0.016, p < 1e-4); O-type + ST together still redundant (+0.001, p = 0.92); Klebsiella capsule now individually significant (-0.004, p = 0.003); host sharing -0.003 (p = 0.06). Difference: removing all three phylogeny channels (UMAP + O-type + ST) costs only -0.003 (p = 0.43) here, versus -0.012 (p = 0.035) with the shipped inputs, so the earlier advice to keep at least one phylogeny channel is not supported on regenerated inputs. Note: `results/ablation/` (shipped inputs) comes from the first run, before the `saga` solver was seeded; the GitHub repo has the seeded run, with the same conclusions] Re-run the feature ablation on regenerated features.
### Checkpoint 8: REACHED 2026-10-17 - replicated numbers vs paper known (T29, T30, T31 lines above).

## Phase 9 - Comparison deliverables
- [x] T32 Comparison CSVs in `results/comparison/`: `paper_vs_replication.csv` (46 rows, each with paper source), `per_phage_auroc.csv`, `feature_agreement.csv`, `matrix_and_taxonomy_agreement.csv`, `ablation_comparison.csv`, `fig4_hits_by_factor.csv`, `file_inventory.csv` (scripts/build_comparison.py).
- [x] T33 Four comparison figures with CSVs: `cmp_headline`, `cmp_per_phage_auroc`, `cmp_feature_agreement`, `cmp_ablation` (scripts/fig_comparison.py). Dataviz: 2-series palette slots 1-2, legend, no dual axis, rendered and checked by eye; the palette validator could not run (node missing libsimdjson), so the documented validated slots were used.
- [x] T34 `config/config_diff.csv`: 24 differences with reason, expected effect and paper line. Not yet re-read against `config/paper_config.yaml` line by line for omissions.
- [x] T35 Comparison report with the file inventory (provided vs still needed): `reports/05_paper_comparison.md` in the GitHub repo (48 numbers compared, 7 differ, list of missing files). Finalised on the user's request 2026-10-18 and published; `docs/COMPARISON_REPORT.md` now points to it.
- [x] T36 `/ponytail-audit` run over `scripts/`: 11 ranked findings in `docs/PONYTAIL_AUDIT.md` (about -135 lines, -2 envs). Nothing applied; awaiting approval. The separate `/ponytail` write-mode pass was not run (the audit covers the same ground).
### Checkpoint 9: ready for review.

## Phase 10 - Figure replication (added 2026-10-16; see docs/FIGURE_MAP.md)
- [x] T37 Fig 1 panels b-j from shipped tables, then 1a from the regenerated tree. Done: `results/figures/fig1_b_j.*`, `fig1a_tree.*` (scripts/fig_main_1_2.py, fig1a_tree.py). 1k not feasible (see T42).
- [x] T38 Fig 2 c-f. Done: `fig2_c`, `fig2_d_e`, `fig2_f`. ANOVA F = 113.9 vs paper 114.7 (public matrix has 402 strains, paper 403); morphotype means 41.3/11.2/9.1 vs 41.3/11.1/9.1.
- [x] T39 Fig 3a UMAP + PERMANOVA. Done: `fig3_a`. Adjusted R2: isolation strain 0.55 (paper 0.53), genus + strain 0.69 (0.65), genus 0.27 (0.13, NOT matched; cause unknown).
- [x] T40 Fig 4 a-d, Fig 5 a/c/d, Fig 6 c-e. Done: `fig4_a_d`, `fig5_a_c_d`, `fig6_c_e` (scripts/fig_extended.py, fig_main_3_6.py).
- [x] T41 Extended Data 1, 2, 4A, 5, 6A/C, 7, 9. Done: `ed1`, `ed2`, `ed4a`, `ed5`, `ed6_a_c`, `ed7`, `ed9`. Not matched: ED7 triplet coverage 87.6% vs 63%; ED9A P 4e-9 vs 3.5e-12. Not run: ED6 B/D (no AMR annotation), ED1 MG1655 (not shipped), ED2 uses NODF not BINMATNEST.
- [x] T42 Schematics drawn: `schematic_fig2ab`, `schematic_fig6a`, `schematic_ed8` (box text from the legends only). Not feasible: 1k (1,026 GenBank genomes + vConTACT2), 3b and ED3 (phage genome annotation and alignments), 6b (test-strain genomes not shipped).
### Checkpoint 10: all replicable panels drawn, side by side with the paper's.

## T11 status (2026-10-16, hand-verified; `results/env_check.csv`)
- Runs natively: ECTyper 1.0.0, MacSyFinder 2.1.6, BLAST 2.17.0, ClermonTyping 24.02.
- Broken: DefenseFinder 1.2.0 (missing macsypy), PanACoTA 1.3.1 (needs Rosetta or a pip install), SRST2 (needs Rosetta), Kaptive 2.0.6 (needs biopython <1.78).
- `--override-channels` was missing, so micromamba also used the commercial Anaconda `defaults` channel; fix before rebuilding.
- `.mamba/` uses ~155 GB on this exFAT drive (small-file cluster waste, no hardlinks); decision pending on moving it to the internal disk.

## Update (2026-10-17)
- Rosetta installed; envs rebuilt on the internal disk (`~/micromamba-phage`, ~13 GB; internal disk has ~27 GB free, watch it). All 8 envs pass their check (`results/env_check.csv`). Fixes: PanACoTA needs `python-utils<3.0`; DefenseFinder 1.2.0 needs `macsyfinder<=2.1.4`; Kaptive 2.0.6 runs as the authors' script from tag v2.0.6 (`tools/Kaptive`).
- T15 UMAP calibrated (`results/replication/umap_calibration.csv`): settings not identifiable; defaults n_neighbors=15, min_dist=0.1 used.

## Phase 11 - Publish to GitHub (2026-10-18, requested by the user)
- [x] T43 `/ponytail-audit` findings applied (9 of 11; #8 and #11 deliberately not, reasons in `reports/08_code_audit.md`) and `code-simplification` pass; outputs of the Fig 1-2 scripts checked byte-identical after the change.
- [x] T44 All model runs repeated with fixed seeds (`scripts/run_seeded_models.sh`) and the R fits repeated after the refactor (`lmm_common.R`); figures and comparison tables regenerated. Fig 4c-d changed only within MCMC noise (beta -0.233 / -0.033). Found and fixed on the way: the comparison table used the preprint's 77% as the paper's average AUROC (the published article says 85%); a 4-digit rounding bug turned U = 1171.5 into 1172.
- [x] T45 Nine English reports written in `reports/` from the current results; Vietnamese notes (`README_REPRO.md`, `docs/FEATURES.md`, `docs/RESULTS.md`) superseded.
- [x] T46 Commit and push to the existing repository (author Axel). Pushed 2026-10-18 as commit bebecfd on `main`; local copy `publish/phage-host-feature-necessity/`.

## Phase 12 - Close the gaps with the publisher's supplement (2026-10-18, approved by the user)
- [x] T47 Downloaded the 13 publisher files (Supplementary Information, Reporting Summary, Peer Review File, Supplementary Tables 1-7, Source Data for Fig 1-6, ED6, ED7, ED9; 8.7 MB) to `docs/paper/supplement/` (`scripts/download_supplement.sh`; not pushed). Checked the authors' 14 commits and the Zenodo record: no CV-group file and no recommender code anywhere.
- [x] T48 `scripts/compare_source_data.py`: exact phage-name map; the 403rd strain is LF110; ANOVA F 114.66 with 403 strains; Fig 5a is per strain (paper Source Data mean 0.872, ours 0.881); genus PERMANOVA 0.13 = sequential increment; Fig 4 Source Data uses pMCMC 0.01 (paper 52/1 terms, ours 51/3); ED2 Q 0.23 matches; ED7 and ED9b/ST7 labels are inconsistent inside the paper. Comparison table now 73 rows (12 differ).
- [x] T49 ED6 B and D drawn with the paper's AMR counts (genome size vs mean MLC significant, AMR not at 0.01). No independent AMR annotation run.
- [x] T50 Recommender re-implemented from the Methods (`scripts/recommender.py`): 2 of 100 cocktails identical, 0.97 of 3 phages in common; the haplotype step is reproduced (10 of 11), the two random-forest steps are not. Test-strain ABC/Klebsiella types and genomes are unavailable.
- [x] T51 Reports 02, 05, 06, 07, 08 and the README updated; pushed to the existing repository (author Axel). Publisher data is not in the repository.

## Phase 13 - Everything that can be run from the authors' repository (2026-10-18, requested by the user)
- [x] T52 Execute the four notebooks of the repo (bokeh interaction matrix, visualize_predictions, cocktail analysis, defense-arsenal distance table) as shipped (27 of 61 code cells fail: Windows/private paths, moved file, missing package, newer pandas) and with patches (1 authors' cell still fails). Outputs in `notebooks/`; report 10.
- [x] T53 Supplementary analyses computable from the repo data: `scripts/si_analyses.py`, `si_rbp_lmm.R`: SI Text 4 matches (2.32e-4, P 0.17); SI 6 K1 rule 87/55 % vs 73/72 %; SI 7 RBP LMM not significant (beta -13 vs paper +28); SI 8 PERMANOVA ordering reproduced (AF -458, full -437, defense -427); RBP sharing vs host range (genus-confounded).
- [x] T54 Fig 3b and ED3: pyrodigal + MAFFT (internal micromamba env `aln`, removed afterwards); entropy 0.62 bits in RBP genes vs 0.24 elsewhere.
- [x] T55 28 differences in `reports/11_repository_vs_paper.md`; also found: ED9a P 3.5e-12 comes from counting each strain three times (6.4e-5 once); Fig 4d defense counts differ from Fig 4c; BCH857/BCH857-72 name mismatch explains 99 vs 100 strains. Comparison table now 97 rows.

## Phase 14 - Matrix-only analyses on the regenerated matrix (2026-10-18, requested by the user)
- [x] T56 `scripts/regen_matrix_analyses.py`: ANOVA, PERMANOVA, NODF/BRIM and ED5 on the matrix regenerated from raw plaque scores (369 strains, no isolation strains), shown against the authoritative matrix with all strains and with the same 369 strains. Regenerating the matrix changes almost nothing; the gap to the paper comes from the missing 33 strains. Helpers moved to `scripts/matrix_stats.py` (outputs of the two figure scripts checked identical).
- [x] T57 Report 05 now lists why each empty cell of the headline table is empty (not published by the authors / cannot be run: Fig 6 wet-lab results and recommender / not run).
- [ ] Not run on regenerated inputs (could be): ED6, ED7, ED9, the supplementary analyses (SI texts 4-8), the RBP-count LMM.

## Phase 15 - Remaining analyses on regenerated inputs (2026-10-18, requested by the user)
- [x] T58 Fig 4c-d and ED6 mixed models on the regenerated matrix (`MATRIX_SET=regenerated`): 4c beta -0.193 (p 0.51), 4d -0.034; ED6 genome size vs avg MLC p 0.014 (0.003 on shipped); ED7 coverage 66.1/87.0/90.5 %.
- [x] T59 Supplementary analyses (SI 4, 6, 7, 8, RBP sharing) on regenerated inputs (`MATRIX_SET=regenerated`): K1 count 45 = the paper's (shipped table 42); SI 8 AF/full/defense adj R2 0.25/0.28/0.14.
- [x] T60 ED9: not regenerable (test-strain data); state it. ED6 AMR counts stay the paper's.
- [x] T61 Update reports 05/10, push.

## Status log (newest last; updated after every step)
- 2026-10-17 T11 timing; PanACoTA annotate Prodigal vs Prokka compared on 5 genomes (Prodigal 2-7% more CDS, 90-97% of its proteins identical to Prokka's).
- 2026-10-17 T18 done (ECTyper 403 genomes, 12 min). T16/T17 done.
- 2026-10-17 T14 annotate relaxed QC (default dropped 10 genomes present in the published tree); pangenome/core done.
- 2026-10-17 T21 Kaptive switched to 3 parallel chunks (single run was ~15 h).
- 2026-10-17 T23 OMP done for 11 genes; threshold sweep done.
- 2026-10-17 Housekeeping: micromamba cache cleaned (internal root 15 -> 12 GB); `replication/timing` and `publish/.../upstream` deleted; old `.mamba` (155 GB, external) being deleted in background. Internal disk free 14.6 GiB (swap file 11.7 GB is the main use); external free ~550 GB.
- 2026-10-17 References downloaded to the external drive only: `replication/T19_lps/refs/` (waaL P27243, yncD P76115; see `config/references.yaml`).
- 2026-10-17 T19 first result written to `results/replication/lps_agreement.csv`; next: improve waaL detection for R2/R3.
- 2026-10-17 T19 done (see its line above); files: `results/replication/lps_agreement.csv`, `lps_markers.csv`, `lps_confusion.csv`.
- 2026-10-17 T23 done with yncD (ARI 0.916). Old `.mamba` (155 GB on the external drive) deletion finished.
- 2026-10-17 T20 done (see its line). T22 DefenseFinder started (3 genomes at a time).
- 2026-10-17 T19: tblastn of K-12 WaaL on all genomes confirms the protein result (K12 43/43, R1 206/206 at ~30% identity, R4 35/35, R2 0/22, R3 2/61); the repo README says per-type representative sequences would be needed (not shipped). Scratch files (2 GB merged genomes, BLAST index) deleted.
- 2026-10-17 T24 done (see its line).
- 2026-10-17 `scripts/compare_tree.py` (T14 + T15 comparison) written and smoke-tested against the shipped ultrametric tree (splits RF 0.04, so the script works); the smoke-test outputs were deleted, the real run waits for our tree.
- 2026-10-17 T14 core alignment finished; FastTree started (SPR round 1 of 2).
- 2026-10-17 T14 and T15 done (see their lines): tree matches the published one almost exactly.
- 2026-10-17 T22 done (see its line).
- 2026-10-17 T21 done (see its line).
- 2026-10-17 T25 done; Checkpoint 6 reached. micromamba cache cleaned again after the jobs finished.
- 2026-10-17 after an unexpected computer reset: integrity re-checked (403/403 genome md5, all result files readable, all 9 tool envs run, Rosetta works, GitHub repo intact); removed ~2,000 BLAST index files that ECTyper had left beside the genomes. Note: `phase6_summary.csv` has 33 metric rows (not 35).
- 2026-10-17 internal disk cleaned at the user's request: removed `~/micromamba-phage` (12 GB: 9 tool environments). Rebuild any time with `scripts/setup_envs.sh` (specs in `config/envs/*.yml`, ~30 min); the `micromamba` program itself (Homebrew) is kept.
- 2026-10-17 external drive cleaned at the user's request with `scripts/cleanup_intermediates.sh`: removed PanACoTA work files, per-genome HMMER outputs of the capsule step, ClermonTyping per-genome files, the OMP BLAST index, Kaptive temp dirs and ~150,000 macOS `._*` twin files. Drive used 287 -> 209.8 GiB (free 721.7 GiB); `replication/` 83 GB -> 6.6 GB; project 14 GB. Kept: all result tables, genomes, proteomes and genes, the pangenome family list, final core alignment, tree, per-genome key outputs. Verified: all comparison scripts re-run from the kept files and give identical result tables (only change: the regenerated table now also has the `ABC_model` column, after fixing a bug in `compare_features.py` where that column made the comparison loop fail).
- 2026-10-17 Phase 8 COMPLETE (T29, T30, T31); Checkpoint 8 reached. Next: Phase 9 (comparison CSVs, figures, config differences, report).
- 2026-10-17 Phase 7: T26 done, T27 done (approximation), T28 not done (reasons on its line). BLAST env rebuilt for T27 and removed again (internal disk).
- Phase 6 COMPLETE (T14-T25). Phase 7 COMPLETE except T28. Next: Phase 8.
- 2026-10-17 Phase 10 COMPLETE (T37-T42): 19 figures in `results/figures/` (PNG, SVG, CSV each) + `stats_*.csv` (paper vs ours). Unmatched numbers listed in T38-T41. Next: Phase 9 (comparison CSVs, figures, config diff, report).
- 2026-10-17 Phase 9 done except review (T32-T34, T36 done; T35 awaiting user review). Disk after cleanup: external 722 GiB free, internal 22 GiB free.
- 2026-10-18 Phase 11 COMPLETE (T43-T46): repo reorganised, nine English reports in `reports/`, pushed to the existing GitHub repository.
- 2026-10-18 T35 finalised. Time Machine local snapshots deleted on the user's instruction (24 snapshots; internal disk free 28 GiB -> 65 GiB).
- 2026-10-18 GitHub repo history rewritten at the user's request: one commit `ab85ee0` on `main` (old commits removed from the branch); README rewritten around the from-scratch replication. Local `git gc` fails on exFAT `._*` files in `.git` (cosmetic).
- 2026-10-18 Phase 12 COMPLETE (T47-T51).
- 2026-10-18 Phase 13 COMPLETE (T52-T55); reports 10 and 11 written.
- 2026-10-18 Phase 14 COMPLETE (T56-T57).
- 2026-10-18 Phase 15 COMPLETE (T58-T61). 38 rows of the comparison table now have a regenerated value; Fig 4a-b and ED9 are the only items not run on regenerated inputs.

## Phase 16 - Annotated figures (2026-10-18, requested by the user)
- [x] T62 `reports/06_figures.md` rewritten: 31 figures, each with an annotation table (layout, axes, marks and colours, data file, how to read, numbers, paper comparison, caveats); one-paragraph annotations under every figure in reports 05 and 10; the headline figure now clips outliers at the axis edge.
