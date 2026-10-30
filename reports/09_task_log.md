# 9. Task log

One line per task, with what came out of it. `tasks/todo.md` in the working copy holds the longer running log. File names that no longer exist (`run_lps.py`, `fig6_cocktails.py`) were removed in the code audit ([08](08_code_audit.md)).

## Phase 1: data

- T1: I cloned the authors' repository and confirmed that the tables join (402 strains across matrix, features, tree and UMAP).
- T2: I read the paper to see what it had already tested. Defense systems are out, and adsorption factors plus phylogeny are enough.
- T3: `scripts/download_data.sh` fetches the repository, unpacks the archive and optionally downloads the genomes with md5 checks.

## Phase 2: reproduce the models

- T4: `scripts/make_groups.py` rebuilds the cross-validation groups. Checked against the authors' fold logs, 0.3 % of groups are split across folds.
- T5: `scripts/predict.py` reproduces the per-phage models. Pooled AUROC is 0.856 against the paper's 0.86, and the per-phage mean is 0.787 against the preprint's 0.77.

## Phase 3: feature necessity

- T6: `scripts/ablate.py` drops each feature group in turn and aggregates the authors' importances.
- T7: The verdict table says which groups to keep and which are redundant, with a column for what the paper already tested.

## Phase 4: documents and audit

- T8: I wrote the first notes and result tables. The English reports replaced them later.
- T9: I ran the ponytail audit on the scripts and applied its fixes.
- T10: I switched Remote Control on.
- Extra: I re-ran the Fig 6 statistics and the Fig 4 mixed models. 91 of 96 phages fitted at first, and the other five needed a prior on the fixed effects (T30).

## Phase 5: feasibility and configuration

- T11: A spike installed all nine tools and timed them. Prodigal annotates five genomes in 18 s where Prokka needs about 10 minutes under Rosetta.
- T12: I downloaded the 403 assemblies and checked their md5 sums. I skipped the 9.9 GB of plaque images.
- T13: `config/paper_config.yaml` pins the paper's tool versions and parameters, each with its line in the text. Unstated settings are marked.

## Phase 6: bacterial features

- T14: The core-genome tree shares 390 of 399 splits with the shipped tree (normalised Robinson-Foulds 0.023), and the patristic distances correlate at 0.99999988.
- T15: The 8-D UMAP of our tree has a 10-nearest-neighbour overlap of 0.748 with the shipped coordinates. UMAP's own run-to-run overlap is 0.77 to 0.80.
- T16: mlst replaces SRST2 and agrees on 95.3 % of strains. Most mismatches are novel STs that the current database numbers differently.
- T17: ClermonTyping agrees on 99.0 % of strains.
- T18: ECTyper output is identical to the authors' raw output for 402 of 402 strains. Against their hand-curated columns it agrees on 74.7 % (O) and 91.8 % (H).
- T19: K-12 WaaL finds only types K12, R1 and R4. Pangenome marker families, trained on the shipped labels, predict the LPS type for 97.0 % of strains in cross-validation. R2 is the weak class at 72.7 %.
- T20: The ABC capsule presence agrees on 98.3 % of strains and the K-type on 89.5 %. K2 is the weak type, with 5 of 10.
- T21: Kaptive agrees on 97.3 % of strains and reproduces 20 of the 22 shipped calls. One malformed assembly (ROAR131) needed a cleaned copy.
- T22: DefenseFinder gives presence F1 0.976. It does not call RM_Type_IV, which accounts for 123 missing systems.
- T23: Outer-membrane protein variants agree with the shipped partition with adjusted Rand index between 0.06 (tolC) and 0.96 (ompC). The authors' exact clustering is not recovered.
- T24: The regenerated host table matches the shipped one for 33 of 33 strains on phylogroup, ST and LPS. O and H differ in one strain, only in how "no call" is written.
- T25: `scripts/assemble_features.py` builds the regenerated feature table and a 33-metric comparison with the shipped one.

## Phase 7: phage and interaction data

- T26: The interaction matrix rebuilt from raw plaque scores agrees on lytic or not for 99.8 % of 35,268 comparable cells, and on the exact score for 96.6 %.
- T27: Phage taxonomy uses a BLASTn approximation of VIRIDIC. Genus matches well (adjusted Rand index 0.929), and species does not (0.241).
- T28: Not done. RBP detection needs PhageRBPdetect, Pfam-A 34 and ESMFold, which do not fit a 16 GB Mac without a GPU. The shipped `RBP_list.csv` stays the reference.

## Phase 8: models on regenerated inputs

- T29: With regenerated features, tree and groups the pooled AUROC is 0.855, and the per-phage mean falls from 0.787 to 0.765. With our matrix as well it is 0.853 and 0.769.
- T30: The five phages that failed as singular were refitted with a prior on the fixed effects. At pMCMC < 1e-3 there are 28 adsorption-factor terms and 0 defense terms. Fig 4c-d gives β = -0.231 (paper -0.22) and -0.033 (paper -0.04).
- T31: The ablation on regenerated inputs keeps the ABC capsule and LPS as necessary and the O-type and ST pair as redundant. Removing all phylogeny channels costs only 0.003 here against 0.012 on shipped inputs.

## Phase 9: comparison

- T32: `results/comparison/` holds the comparison CSVs, each paper number with its source.
- T33: Four comparison figures with CSVs. The palette validator would not run, so I used the documented palette slots and checked by eye.
- T34: `config/config_diff.csv` lists the differences between the paper's settings and what I ran. It has 32 rows by the end of the work.
- T35: The comparison report (05) covers the files provided and the files still needed. I finalised it on request and published it.
- T36: The ponytail audit of `scripts/` gave 11 ranked findings.

## Phase 10: figures

- T37: Fig 1b to j from the shipped tables and Fig 1a from the regenerated tree. Fig 1k needs 1,026 GenBank genomes and vConTACT2, so I did not draw it.
- T38: Fig 2c to f. The ANOVA F is 113.9 with 402 strains and 114.7 with all 403.
- T39: Fig 3a. PERMANOVA adjusted R² is 0.55 for isolation strain (paper 0.53) and 0.69 for genus plus strain (paper 0.65). The genus value of 0.27 against 0.13 was explained in T48.
- T40: Fig 4a to d, Fig 5a, c, d and Fig 6c to e.
- T41: ED1, ED2, ED4A, ED5, ED6, ED7 and ED9. ED7 triplet coverage (87.6 % against 63 %) and the ED9a P value do not match.
- T42: Schematics for Fig 2a-b, Fig 6a and ED8, with box text taken from the legends. Fig 6b stays undrawn because the test-strain genomes are not shipped.

## Phase 11: publish to GitHub

- T43: I applied 9 of the 11 audit findings and ran the code-simplification pass. The Fig 1 to 2 outputs are byte-identical afterwards.
- T44: I repeated all model runs with fixed seeds and the R fits after the refactor. On the way I fixed two comparison bugs. One used the preprint's 77 % as the paper's average AUROC. The other was a four-digit rounding that turned U = 1171.5 into 1172.
- T45: I wrote the English reports from the current results and retired the Vietnamese notes.
- T46: I pushed to the existing repository as Axel. At the user's request I later replaced its history with one commit.

## Phase 12: the publisher's supplement

- T47: I downloaded the 13 publisher files (8.7 MB). Neither the authors' 14 commits nor the Zenodo record holds the cross-validation groups or the recommender code.
- T48: `scripts/compare_source_data.py` settled several differences. It gave the exact phage-name map, identified the 403rd strain as LF110, gave ANOVA F = 114.66, showed that Fig 5a is per strain, explained the genus PERMANOVA of 0.13, and showed that Fig 4 counts terms at pMCMC 0.01.
- T49: ED6 B and D use the paper's AMR counts. Genome size against average MLC is significant, and AMR is not at 0.01.
- T50: The recommender, re-implemented from the README, reproduces 2 of 100 cocktails. The haplotype step matches (10 of 11), and the two forest steps do not.
- T51: I updated reports 02, 05, 06, 07 and 08 and the README, and pushed. The publisher's data stays out of the repository.

## Phase 13: everything runnable from the authors' repository

- T52: I ran the four notebooks as shipped (27 of 61 code cells fail) and with patches (one cell of the authors' still fails).
- T53: Supplementary analyses. SI 4 matches (2.32e-4, P = 0.17). The SI 6 K1 rule gives 87 % precision and 55 % recall against the paper's 73 % and 72 %. The SI 7 RBP model is not significant. SI 8 reproduces the PERMANOVA ordering. RBP sharing tracks host range but cannot be separated from genus.
- T54: Fig 3b and ED3 use pyrodigal and MAFFT. Entropy is 0.62 bits in the RBP genes against 0.24 elsewhere.
- T55: I collected 28 differences between the repository and the paper in report 11. Among them, the ED9a P value of 3.5e-12 counts each strain three times (6.4e-5 with each strain once), the Fig 4d defense counts differ from Fig 4c, and a BCH857 name mismatch explains 99 against 100 strains.

## Phase 14: matrix-only analyses on the regenerated matrix

- T56: ANOVA, PERMANOVA, NODF, BRIM and ED5 on the regenerated matrix (369 strains). Regenerating changes almost nothing, and the gap to the paper comes from the 33 missing isolation strains.
- T57: Report 05 explains why each empty cell of the headline table is empty.

## Phase 15: remaining analyses on regenerated inputs

- T58: Fig 4c-d and ED6 on the regenerated matrix give β = -0.193 (p = 0.51) and -0.034. Genome size against average MLC has p = 0.014. ED7 coverage is 66.1, 87.0 and 90.5 %.
- T59: SI 4, 6, 7, 8 and RBP sharing on regenerated inputs. The K1 count is 45, the paper's number, where the shipped table has 42.
- T60: ED9 cannot be regenerated, because it needs the test-strain data. ED6 keeps the paper's AMR counts.
- T61: I updated reports 05 and 10 and pushed.

## Phase 16: annotated figures

- T62: Report 06 now annotates 31 figures, with layout, axes, marks, data file, reading guide, numbers, paper comparison and caveats for each. Reports 05 and 10 have a short annotation under every figure, and the headline figure clips outliers at the axis edge.

## Phase 17: plain-prose rewrite

- T63: I rewrote the README and reports 01 to 11 in plain prose. Numbers, tables, code blocks, paths and links are unchanged.

## Phase 18: gap audit

- T64: I built a register of the 47 rows of the paper comparison that did not match, each with a cause code.
- T65: All tiers have 402 strains and 96 phages. The authors' predictions lack 58 pairs that the shipped matrix has and hold 1 pair it lacks.
- T66: Overall metrics with bootstrap intervals for the authors' run and ours. The shipped-input run lies inside the interval on every metric.
- T67: The 28 differences between the repository and the paper each have a code: P 20, D 4, U 3, C 1.
- T68, T69: Per-phage table for 96 phages and a cause for every phage that moves by more than 0.05.
- T70, T71: Pair-level comparison. Calls agree for 94.5 % of pairs between the authors' run and ours.
- T72a: With the authors' fold assignment the mean change in per-phage AUROC falls from 0.026 to 0.0065, and 34 of 96 phages match to six decimals.
- T72: Seed noise with fixed folds is 0.0011 on the mean per-phage AUROC and 0.0031 with a random partition.
- T73, T73b: Regenerated features cost 0.018 mean per-phage AUROC with the folds fixed. The O-type and Klebsiella capsule type carry it, mostly through phages with few lytic strains.
- T74: ED7 (4 thresholds, 12 raw-plaque definitions) and Fig 3a (LF110, genus labels) tested; 8 rows stay residual.
- T75: Report 12 written from the audit tables.
- T76: I checked the audit against an independent bioRxiv preprint (Noonan et al.) on the same matrix. Dataset counts and overall AUROC (0.869 against 0.86) agree, and its promiscuity baseline supports my ED7 reading.
- T77: I compared the original preprint with the published paper and the repository (23 statements). Two register rows wrongly said the paper gives no number and are fixed. The 68 % four-phage coverage is from the preprint only. The authors' folds fit greedy group packing of our groups once one group is taken apart, and scikit-learn 1.1.2 against 1.9.1 does not explain the fold difference.

## Phase 20: pinned environment rerun

- T81: I built the conda environment `ml_pinned` (scikit-learn 1.1.2, scipy 1.11.2, Python 3.10.22, numpy 1.23.5, pandas 1.5.3, joblib 1.2.0). The authors state no Python version for the model step, so 3.10 follows scikit-learn 1.1.2. `predict.py` now passes `penalty="l1"` below scikit-learn 1.8, because 1.1.2 fits L2 when only `l1_ratio=1` is given.
- T82: Reruns T2.1 (shipped inputs), T2.1 with the authors' folds, and T3.1 (regenerated features). With the authors' folds the mean change in per-phage AUROC is 0.0042 and mean AUROC 0.7813 against 0.7812. LogReg L2 probabilities are identical for 72.2 % of pairs. With our own groups the version changes the folds and does not close the gap (0.0296 against 0.0261). Every unmatched component of the authors' environment is listed with its reason in `results/audit/env_match.csv`.
