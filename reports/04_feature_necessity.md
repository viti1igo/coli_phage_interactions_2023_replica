# 4. Which genomic features does the model need?

## Features of the final model

The final model of the paper (`scripts/predict.py`, authors' README §1.4–1.5) uses these groups:

| Group (name in the ablation) | Columns | Tool | Note |
|---|---|---|---|
| `UMAP_phylogeny` | UMAP0–7 | core-genome tree, distance matrix, 8-D UMAP | z-scored |
| `O_type` | `O-type` | ECTyper | 94 values, those with fewer than 3 strains merged into "Other" |
| `ST` | `ST_Warwick` | SRST2 (MLST) | 163 values, same merging |
| `LPS_outer_core` | `LPS_type` | blastp of *waaL* and manual grouping | R1, R2, R3, R4, K12, No_waaL |
| `Klebsiella_capsule` | `Klebs_capsule_type` | Kaptive | 381 of 403 strains empty |
| `ABC_capsule` | `ABC_serotype` | MacSyFinder with the authors' own models | 252 of 403 empty |
| `host_sharing` | `same_O_as_host`, `same_ST_as_host`, `same_O_and_ST_as_host`, `same_ABC_as_host`, `ST_host` | comparison with the phage's isolation strain | not available for the four LF110 phages (no host genome) |

Quirks of the original code, kept so that the numbers reproduce:
- `same_ABC_as_host` is computed as `ABC_serotype == ABC_serotype`. It does not compare with the host; it is true exactly when the strain has an ABC serotype (missing is not equal to missing), so in effect it is a "has an ABC capsule" flag.
- The column filter drops `H-type`, so there is no `H_host`. `O_host` is dropped too; `ST_host` is kept (many levels, many dummy columns).
- The UMAP columns are z-scored on all data of a phage before cross-validation (a mild leak).

Columns in `picard_collection.csv` that the model does not use: `Clermont_Phylo` (the information is in UMAP and ST), `H-type` (excluded by the original filter; used in the Fig 4 models), `Pathotype`, `Host`, `Origin`, `Mouse_killed_10`, `Collection` (metadata), the raw binary capsule columns (replaced by `ABC_serotype` and Kaptive), `n_infections` (a summary of the label, which would leak it), `n_defense_systems` and the defense table, and the outer-membrane-protein clusters (Fig 4 only).

## What the paper already tested

- **Defense systems** were excluded by the authors. Fig 4 and the section "Defense systems are not required to predict phage–bacteria interactions in our dataset" report β = −0.22 (95 % CrI −0.74 to 0.34, pMCMC 0.41) on the number of infecting phages and only a very small effect on MLC (β = −0.04). This was therefore not repeated; [02](02_model_reproduction.md) reproduces the estimates.
- Adsorption factors plus phylogeny are sufficient (per-phage AUROC and pooled AUROC 86 %, Fig 5).
- The screening before the model was per-phage mixed models (Fig 4a–b), not an ablation of the final model.
- Not found in the preprint or in the main text of the published article: dropping each of the final feature groups under group cross-validation. That is the analysis below. The Supplementary Text was not available and may contain related analyses.
- The paper reports 119 input features; 114 columns were counted for a sample phage (the count depends on the phage through the merging of rare values and the host-sharing group).
- Model selection per phage is by mean AUPR in the script and the README, by AUROC in the paper text; the script is followed.

## Ablation

Design: the whole pipeline is kept (5 models, group 10-fold cross-validation, best model per phage by mean AUPR). One feature group, or a combination, is removed each time, and each phage is compared with the full model (96 pairs, Wilcoxon signed-rank test). The p-values are not corrected for 12 tests (Bonferroni threshold about 0.004). Script: `scripts/ablate.py`; results in `results/ablation` (shipped inputs) and `results/ablation_regenerated` (inputs rebuilt from the genomes). Both runs use fixed seeds.

| Removed | ΔAUROC per phage (shipped inputs) | ΔAUPR | p | Phages worse by more than 0.02 AUROC | ΔAUROC (regenerated inputs) | p |
|---|---|---|---|---|---|---|
| none | baseline 0.787 | baseline 0.591 | – | – | baseline 0.765 | – |
| **ABC capsule** | **−0.021** | −0.053 | 1e-4 | 34 | **−0.019** | < 1e-4 |
| **LPS outer core** | **−0.016** | −0.017 | 2e-4 | 32 | **−0.016** | < 1e-4 |
| Klebsiella capsule | −0.005 | −0.013 | 0.15 | 9 | −0.004 | 0.003 |
| Host sharing | −0.004 | −0.004 | 0.046 | 11 | −0.003 | 0.059 |
| UMAP phylogeny | +0.001 | −0.014 | 0.93 | 21 | +0.006 | 0.48 |
| O-type | +0.001 | −0.002 | 0.68 | 15 | +0.001 | 0.96 |
| ST | +0.002 | −0.001 | 0.56 | 10 | −0.001 | 0.95 |
| **O-type and ST** | **+0.001** | **+0.002** | 0.92 | 18 | +0.001 | 0.92 |
| UMAP and O-type | +0.000 | −0.019 | 0.43 | 31 | +0.006 | 0.86 |
| UMAP and ST | −0.003 | −0.026 | 0.27 | 29 | −0.000 | 0.47 |
| UMAP, O-type and ST (all phylogeny) | −0.012 | −0.044 | 0.035 | 44 | −0.003 | 0.43 |
| Klebsiella capsule and host sharing | −0.009 | −0.014 | 0.003 | 16 | −0.006 | 0.007 |

Table of both runs with verdicts: `results/comparison/ablation_comparison.csv`; chart: `cmp_ablation` ([06](06_figures.md)).

## Reading the results

1. **Keep the ABC capsule and the LPS outer-core type.** Removing either costs AUROC consistently (more than 30 of 96 phages worse), on both input sets. This agrees with the paper's conclusion that adsorption factors are the main axis.
2. **O-type and ST are redundant as a pair.** Removing both changes nothing (ΔAUROC +0.001, p = 0.92) because the UMAP phylogeny and host sharing carry the same information. For a sample phage the number of input columns falls from 114 to 45.
3. **Phylogeny cannot be dropped entirely on the shipped inputs.** Removing UMAP, O-type and ST together costs 0.012 AUROC and 0.044 AUPR (p = 0.035). On the regenerated inputs the same removal costs only 0.003 (p = 0.43), so that advice is not robust: it depends on which inputs are used.
4. **Klebsiella capsule and host sharing are weak.** Each contributes very little alone; removed together they cost 0.009 AUROC (p = 0.003) on shipped inputs and 0.006 (p = 0.007) on regenerated ones. Keeping them is cheap. Their individual verdicts swap between the two input sets (Klebsiella capsule: not significant, then significant; host sharing: the reverse), which is within noise.
5. **Importance is not necessity.** By the authors' shipped feature importances (`results/ablation/importance.csv`) the UMAP phylogeny holds 24 % of the importance, O-type 24 % (top feature for 41 phages) and ST 17 %, yet O-type and ST can be removed. LPS holds 8.5 % and cannot. Importance of forests and logistic regression is shared among correlated features.

Three of twelve verdicts differ between the input sets (Klebsiella capsule, host sharing, and the all-phylogeny removal); the other nine agree.

## Limits

- One seed and one cross-validation partition per run; no interval for the split. Five seeds should be run before removing any feature in practice.
- The best model per phage is re-selected in every variant on the same test folds, which is slightly optimistic; it applies equally to all variants.
- The groups are rebuilt (0.3 % differ from the authors' folds), and the public matrix differs from the private one, so absolute numbers may shift slightly.
- "O-type and ST are redundant" holds for this collection and for strains within its phylogenetic range; strains further from the Picard collection were not tested.
