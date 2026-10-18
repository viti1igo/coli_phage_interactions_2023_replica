# 4. Which genomic features does the model need?

## Features of the final model

The paper's final model (`scripts/predict.py`, sections 1.4 and 1.5 of the authors' README) uses these groups:

| Group (name in the ablation) | Columns | Tool | Note |
|---|---|---|---|
| `UMAP_phylogeny` | UMAP0 to UMAP7 | core-genome tree, distance matrix, 8-D UMAP | z-scored |
| `O_type` | `O-type` | ECTyper | 94 values; values seen in fewer than 3 strains are merged into "Other" |
| `ST` | `ST_Warwick` | SRST2 (MLST) | 163 values, merged the same way |
| `LPS_outer_core` | `LPS_type` | blastp of *waaL* and manual grouping | R1, R2, R3, R4, K12, No_waaL |
| `Klebsiella_capsule` | `Klebs_capsule_type` | Kaptive | 381 of 403 strains empty |
| `ABC_capsule` | `ABC_serotype` | MacSyFinder with the authors' own models | 252 of 403 empty |
| `host_sharing` | `same_O_as_host`, `same_ST_as_host`, `same_O_and_ST_as_host`, `same_ABC_as_host`, `ST_host` | comparison with the phage's isolation strain | not available for the four LF110 phages (no host genome) |

The original code has three quirks, which this work keeps so the numbers reproduce:
- `same_ABC_as_host` is computed as `ABC_serotype == ABC_serotype`. It never looks at the host. It is true exactly when the strain has an ABC serotype (a missing value does not equal itself), so it works as a "has an ABC capsule" flag.
- The column filter drops `H-type`, so there is no `H_host`. `O_host` is dropped too, while `ST_host` stays (many levels, so many dummy columns).
- The UMAP columns are z-scored on all of a phage's data before cross-validation, which is a mild leak.

The model ignores these columns of `picard_collection.csv`: `Clermont_Phylo` (UMAP and ST already carry that information), `H-type` (the original filter excludes it, though the Fig 4 models use it), `Pathotype`, `Host`, `Origin`, `Mouse_killed_10` and `Collection` (metadata), the raw binary capsule columns (replaced by `ABC_serotype` and Kaptive), `n_infections` (a summary of the label, so using it would leak the target), `n_defense_systems` with the defense table, and the outer-membrane-protein clusters (Fig 4 only).

## What the paper already tested

- The authors excluded defense systems. Fig 4 and the section "Defense systems are not required to predict phage-bacteria interactions in our dataset" report β = -0.22 (95 % CrI -0.74 to 0.34, pMCMC 0.41) on the number of infecting phages, and only a very small effect on MLC (β = -0.04). This work did not repeat the test, but [02](02_model_reproduction.md) reproduces the estimates.
- Adsorption factors plus phylogeny are enough (per-phage and pooled AUROC around 86 %, Fig 5).
- Their screening before the model used per-phage mixed models (Fig 4a-b), not an ablation of the final model.
- No part of the preprint, the published article or its supplement describes dropping each final feature group under group cross-validation. That is the analysis below.
- The paper reports 119 input features. A sample phage has 114 columns here. The count varies by phage, because rare values are merged and host sharing differs.
- The script and the README pick the best model per phage by mean AUPR, and the paper's text says AUROC. This work follows the script.

## Ablation

The whole pipeline stays as it is (5 models, group 10-fold cross-validation, best model per phage by mean AUPR). Each run removes one feature group or one combination, and every phage is compared with its full model (96 pairs, Wilcoxon signed-rank test). The p-values are not corrected for the 12 tests, where the Bonferroni threshold would be about 0.004. `scripts/ablate.py` runs it, and the results are in `results/ablation` (shipped inputs) and `results/ablation_regenerated` (inputs rebuilt from the genomes). Both runs use fixed seeds.

| Removed | ΔAUROC per phage (shipped inputs) | ΔAUPR | p | Phages worse by more than 0.02 AUROC | ΔAUROC (regenerated inputs) | p |
|---|---|---|---|---|---|---|
| none | baseline 0.787 | baseline 0.591 | - | - | baseline 0.765 | - |
| ABC capsule | -0.021 | -0.053 | 1e-4 | 34 | -0.019 | < 1e-4 |
| LPS outer core | -0.016 | -0.017 | 2e-4 | 32 | -0.016 | < 1e-4 |
| Klebsiella capsule | -0.005 | -0.013 | 0.15 | 9 | -0.004 | 0.003 |
| Host sharing | -0.004 | -0.004 | 0.046 | 11 | -0.003 | 0.059 |
| UMAP phylogeny | +0.001 | -0.014 | 0.93 | 21 | +0.006 | 0.48 |
| O-type | +0.001 | -0.002 | 0.68 | 15 | +0.001 | 0.96 |
| ST | +0.002 | -0.001 | 0.56 | 10 | -0.001 | 0.95 |
| O-type and ST | +0.001 | +0.002 | 0.92 | 18 | +0.001 | 0.92 |
| UMAP and O-type | +0.000 | -0.019 | 0.43 | 31 | +0.006 | 0.86 |
| UMAP and ST | -0.003 | -0.026 | 0.27 | 29 | -0.000 | 0.47 |
| UMAP, O-type and ST (all phylogeny) | -0.012 | -0.044 | 0.035 | 44 | -0.003 | 0.43 |
| Klebsiella capsule and host sharing | -0.009 | -0.014 | 0.003 | 16 | -0.006 | 0.007 |

`results/comparison/ablation_comparison.csv` holds both runs with a verdict for each, and `cmp_ablation` ([06](06_figures.md)) charts them.

## Reading the results

1. Keep the ABC capsule and the LPS outer-core type. Removing either lowers AUROC on both input sets, and more than 30 of the 96 phages get worse. That agrees with the paper's conclusion that adsorption factors are the main axis.
2. O-type and ST are redundant as a pair. Removing both changes nothing (ΔAUROC +0.001, p = 0.92), because the UMAP phylogeny and host sharing carry the same information. For a sample phage the input columns fall from 114 to 45.
3. On the shipped inputs the phylogeny cannot go entirely. Removing UMAP, O-type and ST together costs 0.012 AUROC and 0.044 AUPR (p = 0.035). On the regenerated inputs the same removal costs only 0.003 (p = 0.43), so the advice depends on which inputs you use.
4. Klebsiella capsule and host sharing are weak. Alone, each adds very little. Removed together they cost 0.009 AUROC (p = 0.003) on shipped inputs and 0.006 (p = 0.007) on regenerated ones, so keeping them is cheap. Their separate verdicts swap between the two input sets (the Klebsiella capsule goes from not significant to significant, host sharing the other way), which is within noise.
5. Importance does not equal necessity. In the authors' shipped feature importances (`results/ablation/importance.csv`) the UMAP phylogeny holds 24 % of the importance, O-type 24 % (the top feature for 41 phages) and ST 17 %. O-type and ST can still be removed. LPS holds 8.5 % and cannot. Forest and logistic-regression importance gets shared among correlated features.

Three of the twelve verdicts differ between the input sets (Klebsiella capsule, host sharing and the all-phylogeny removal). The other nine agree.

## Limits

- Each run uses one seed and one cross-validation partition, so there is no interval for the split. Run five seeds before removing any feature for real.
- Every variant re-selects the best model per phage on the same test folds. That is slightly optimistic, but it applies to all variants equally.
- The groups are rebuilt (0.3 % differ from the authors' folds) and the public matrix differs from the private one, so absolute numbers may shift a little.
- "O-type and ST are redundant" holds for this collection and for strains inside its phylogenetic range. Strains further from the Picard collection were not tested.
