# 11. Where the authors' repository differs from the paper

I found everything below by running the repository's code and notebooks, or by comparing its files with the paper, its Supplementary Information, its Supplementary Tables and its Source Data. "Repository" is `mdmparis/coli_phage_interactions_2023` (14 commits up to 2026-01-22; the paper's Zenodo record is described as a copy of it). A ★ marks an item that changes a number the paper reports.

## Data files

| # | Topic | Repository | Paper | Effect |
|---|---|---|---|---|
| 1 ★ | Strains in the interaction matrix | 402 rows; strain LF110 is missing | 403 strains, 38,688 interactions (Supplementary Table 3 has LF110) | ANOVA F is 113.9 on the public file and 114.7 with LF110, as printed |
| 2 | Phage names | old names (`AL505_Sd2`, `412_P2`, `411_P1`, ...) | new names (`AL505_P3`, `T145_P2`, `PDP110_P1`, ...) | an exact one-to-one map exists (`results/comparison/phage_name_map.csv`); Fig 4b and some Source Data sheets already use the old names |
| 3 ★ | RBP genes per phage | `RBP_list.csv` gives 7 phages no RBP (536_P1, 536_P6, 536_P7, 536_P9, DIJ07_P1, DIJ07_P2, LF31_P1), 48 phages one, 41 two | Supplementary Text 7: every phage encodes one or two | the RBP-count mixed model gives β = -13 (all 96) or +16 (89 phages) against the paper's +28; not significant in all versions |
| 4 | K1 capsule strains | 42 strains with `ABC_serotype = K1` (38 collection + 4 isolation strains) | 45 (Supplementary Text 6) | the K1 decision rule gives precision 87 %, recall 55 % here against 73 % and 72 % |
| 5 | Cocktail strain names | `BCH857` in the cocktail, composition and feature files, `BCH857-72` in `cocktails_step.csv` and Supplementary Table 6 | not stated | the authors' notebook merges on names and silently loses this strain, so some counts have 99 strains (Source Data Fig 6c) and the cocktail 536_P7-LF82_P8-LF82_P9 counts 21 instead of 22 |
| 6 | Strain name used twice | `AN03` is a Picard strain and a test strain | not stated | any join on strain names mixes them |
| 7 | Phage genome files | `T4LD.fna` starts with an empty record | not stated | readers that take the first record fail |
| 8 | Fraction of positive pairs | README: "81 % negative, 19 % positive" | text (both versions): 20.3 % | the shipped matrix has 20.7 % lytic pairs (MLC above 0); the paper's own 20.6 % per phage and 20.8 % per strain imply 20.7 %, so the 20.3 % does not match either |
| 9 | Cross-validation groups | README: the file `370+host_cross_validation_groups_1e-4.csv` "is not available anymore" | Methods describe the grouping | the groups have to be rebuilt from the tree (0.3 % of the rebuilt groups differ from the authors' folds) |

## Statistics and definitions

| # | Topic | Repository | Paper | Effect |
|---|---|---|---|---|
| 10 ★ | Significance threshold of the per-phage models | the Source Data of Fig 4a-b keeps terms up to pMCMC 0.0089 | Methods: 1e-3 | counts at 0.01: 52 adsorption-factor and 1 defense term in the paper, 51 and 3 here; at 1e-3 this re-run gives 28 and 0 |
| 11 | Minimum observations per trait | the R markdown has no such filter (it lumps rare levels at 15) | Methods: traits with fewer than 8 observations were dropped | not reproducible from the code; five phages failed as singular without a prior |
| 12 ★ | Merge of the outer-membrane table | names are not harmonised, so `merge` silently drops 136 of 403 strains | not stated | I corrected it here, and the corrected analysis gives counts close to the paper's at pMCMC 0.01 (51 and 3 against 52 and 1), so the paper's figure probably used a corrected table |
| 13 | Fig 4b classification | `YNCDOther` (an outer-membrane variant) is filed under defense systems in the Source Data | one defense term reported | the paper's only "defense" hit is not a defense system |
| 14 ★ | Defense count in Fig 4d | the Source Data counts differ from those of Fig 4c and Supplementary Table 1 for 353 of 400 strains (typically 2 lower, from 2 higher to 8 lower) | the legend does not say which systems are left out | with the panel-d counts the estimate is the same (β -0.032) and still stronger than the paper's pMCMC (0.0005 against 0.004) |
| 15 ★ | Best model per phage | README and script: mean AUPR | Methods: average AUROC | this replication follows the script |
| 16 | Number of input features | script gives about 114 columns for a sample phage | 119 | within the effect of merging rare values |
| 17 ★ | Average AUROC of Fig 5a | the panel plots one AUROC per strain (Source Data: mean 0.872) | text: 85 % average; the preprint: 77 % per phage | the text's 85 % is below the mean of its own data |
| 18 ★ | PERMANOVA of the phage host range | no code in the repository | text: genus R² = 0.13 | 0.13 is the increase of the adjusted R² when genus is added after isolation strain; genus alone gives 0.27 |
| 19 ★ | ED9a test | the notebook merges the per-strain mean step number back onto the three-rows-per-strain table, so each strain enters the Mann-Whitney test three times | P = 3.5e-12 | with each strain once, P = 6.4e-5 |
| 20 ★ | ED7 coverage | no code in the repository | published text: 63 % (triplet); preprint text: 68 % (four phages); Source Data: 67.0 % and 73.7 % | the curve of the Source Data cannot be obtained from the matrix with any set of MLC scores |
| 21 | ED9b | Source Data "number of haplotypes" equals the strains per cocktail (31, 21, 19, 4) | legend: number of haplotypes | the plotted values are strain counts; the haplotype counts are 20, 18, 18 and 4 |
| 22 ★ | Generic cocktails | Fig 6 analysis: four generic cocktails (24 tailored, 76 generic strains); ED9 cell of the same notebook: three (Supplementary Table 7: 28 and 72) | Fig 6d: 91.67 % = 22/24, 78.95 % = 60/76 | the two figures use different groups |
| 23 | Printed P of Fig 2f | no code in the repository | 1.44e-26 | the one-way ANOVA of its own Source Data gives 7.9e-26 |
| 24 | Overfitting | the notebook reports a mean train-test AUROC gap of 16.3 points and 27 of 96 phages with a significant gap | not reported; the paper says features were reduced to limit overfitting | information the paper leaves out |

## Algorithms described twice

| # | Topic | README of the repository | Methods of the paper |
|---|---|---|---|
| 25 | Recommender, step 1 | compares O-antigen, ST, ABC, Klebsiella and LPS with the isolation strain | compares ST, O, H, LPS and ABC |
| 26 | Recommender, step 2 | phages infecting at least half of the strains of the haplotype, ranked by coverage | phages that infect the majority |
| 27 | Recommender, code | not in the repository | not in the supplement |
| 28 | Original prediction script | `same_ABC_as_host` compares the ABC serotype with itself, so it only flags strains that have an ABC capsule; its regular expression drops the H-type, so there is no `H_host` | host-sharing features described as comparisons with the isolation strain |

## Usability of the repository itself

- The four notebooks fail as shipped: 15 of 18 cells in `visualize_predictions`, 10 of 29 in the cocktail notebook and 2 of 11 in the defense notebook. The causes are Windows and private paths, a moved file, a missing package, newer pandas and matplotlib, and two errors of the authors' own (a figure saved under a `.tsv` name, and a cell that uses the wrong table). After the path and compatibility patches all but one cell run ([10](10_repository_analyses.md)).
- The prediction script `predict_all_phages.py` needs private files and Windows paths. I rewrote it as `scripts/predict.py` with the same models, class weights and cross-validation.
- The README names SRST2 for the sequence type, but only assemblies are public, so SRST2 cannot run. The README and the paper agree on the tool.
- The repository has no licence file.

## What does match

The interaction matrix (402 shared strains), the feature table of Supplementary Table 1, the shipped predictions (AUROC 0.857, F1 0.596, average precision 0.601), the cocktail results (91.67 %, 78.95 %, 81 %; U = 399.5 and 1171.5), the step shares of ED9a, the term counts at pMCMC 0.01 (52 adsorption-factor and 1 defense term in the paper, 51 and 3 here), the supplementary genome-size model (2.3e-4, P = 0.17), the modularity Q = 0.23 and the ANOVA F = 114.7 once LF110 is included.
