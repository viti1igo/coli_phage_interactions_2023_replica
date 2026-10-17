# 2. Model reproduction (Fig 5, Fig 6, Fig 4)

## Method

For each of the 96 phages five models are trained: logistic regression with L1 and L2 penalties, random forests of depth 3 and 6 (250 trees) and a constant dummy model. Training uses group 10-fold cross-validation, where strains closer than 1e-4 substitutions per site on the core-genome tree always share a fold. Class weights follow the prevalence of lytic interactions. The best model per phage is the one with the highest mean test AUPR, as in the authors' script. A pair is lytic when the MLC score is above 0.

Script: `scripts/predict.py`, with `FEATURE_SET=shipped` (default), `regenerated` or `regenerated+matrix`. Results go to `results/repro`, `results/repro_regenerated` and `results/repro_regenerated_matrix`. All random seeds are fixed (random forests and the `saga` solver), so reruns are identical.

| Input set | Features | Cross-validation groups | Labels |
|---|---|---|---|
| `shipped` | authors' `picard_collection.csv` and UMAP | rebuilt from the shipped tree | shipped matrix |
| `regenerated` | tables built from the genomes ([03](03_feature_regeneration.md)), UMAP of our tree | rebuilt from our tree | shipped matrix |
| `regenerated+matrix` | as above | as above | our matrix from the raw plaque scores, where raw data exist |

## Fig 5: per-phage and pooled performance

| | Mean AUROC per phage | Mean AUPR per phage | Pooled AUROC | Pooled AP | F1 (threshold 0.5) | Mean AUROC per strain |
|---|---|---|---|---|---|---|
| Published article | 0.85 | – | 0.86 | 0.60 | 0.60 | – |
| Preprint | 0.77 | – | 0.86 | – | – | – |
| Authors' shipped predictions, recomputed here | 0.781 | 0.586 | 0.857 | 0.601 | 0.596 | 0.881 |
| **This run, shipped inputs** | **0.787** | 0.591 | **0.856** | 0.598 | 0.599 | 0.878 |
| This run, regenerated inputs | 0.765 | 0.569 | 0.855 | 0.596 | 0.592 | – |
| This run, regenerated inputs and our matrix | 0.769 | 0.565 | 0.853 | 0.588 | 0.589 | – |

- The pooled numbers (AUROC, average precision, F1) agree with the paper within 0.004 on the shipped inputs and stay within 0.004 on regenerated inputs.
- The published article states 85 % as the average AUROC of Fig 5a. No aggregation tried reproduces it: 0.787 per phage and 0.878 per strain here, 0.781 and 0.881 from the authors' own shipped predictions. The preprint stated 77 % for the same panel, which matches the per-phage mean. The published text and legend do not say which aggregation gives 85 %.
- The per-phage mean loses about two points when the inputs are regenerated (0.787 to 0.765), and the correlation of per-phage AUROC with the authors' results falls from 0.92 to 0.79 (0.81 with our matrix). The model chosen per phage belongs to the same family (forest or logistic regression) as the authors' for 76 % of phages on shipped inputs.
- Fixed in the code: a phage whose lytic strains all fall in one test fold now gets a constant prediction instead of an error, and `penalty="l1"` was replaced by `l1_ratio=1` (the sparse solution was checked to be identical).

Panels: `fig5_a_c_d` and `results/repro/fig5_curves.png` ([06](06_figures.md)).

## Fig 6: cocktail experiment

The statistics of the authors' notebook were re-run on their shipped experimental results (24 strains received a tailored cocktail, 76 a generic one). The recommender itself was not re-run because its code is not in the repository.

| | Here | Paper |
|---|---|---|
| Lysis, tailored | 91.67 % | 91.67 % |
| Lysis, generic | 78.95 % | 78.95 % |
| Lysis, baseline (100 strains) | 81 % | 81 % |
| MLC tailored vs baseline (Mann–Whitney) | U = 399.5, P = 0.019 | U = 399.5, P = 0.02 (text) |
| MLC tailored vs generic | U = 1171.5, P = 0.034 | P = 0.02 in the text, P = 0.034 in the Fig 6e legend |

Everything matches. For the last test the paper contradicts itself (text 0.02, legend 0.034); the result agrees with the legend.

## Fig 4: per-phage mixed models and defense count

**Per-phage models** (`scripts/fig4_mcmc.R`, `scripts/fig4_summarize.py`). `MCMCglmm` with a threshold family; adsorption factors and defense systems as covariates; the phylogeny as a random effect through the inverse relatedness matrix; 500,000 iterations, burn-in 50,000, thin 250; prior V = 1, ν = 0.002; significance at pMCMC < 1e-3. 91 of 96 phages fitted. Five (409_P3, 411_P1, 411_P2, 412_P1, 423_P10) failed with "mixed model equations singular". They were refitted with an N(0, 10²) prior on the fixed effects, a labelled deviation (`results/fig4/STRONG_PRIOR_PHAGES.txt`); all five have very few lytic strains (3 to 43).

| Result over 96 phages | Count (pMCMC < 1e-3) |
|---|---|
| Adsorption-factor terms | 28 (LPS 14, ABC capsule 12, O-type 1, OMPA 1) |
| Defense-system terms | 0 |

The direction matches the paper: adsorption factors matter, defense systems do not. The paper's counts for Fig 4a–b exist only in the figure, so the counts were not compared. Against the authors' shipped results (adsorption factors only, 96 phages), over 7,776 shared (phage, term) pairs both are significant for 16, only the shipped result for 10 and only this re-run for 12 (Spearman correlation of pMCMC 0.625). The agreement is moderate; the two models differ (with and without defense systems, different covariate sets).

Two things found in the authors' code or text:
- Their Rmd keeps the outer-membrane protein table under the original strain names while the other tables went through `make.names()`, so `merge` silently drops 136 of 403 strains. This was fixed in the R script.
- The Methods say traits with fewer than 8 observations were dropped. The Rmd has no such filter (it lumps rare levels at 15 and filters defense systems at 15). No filter was added here.

**Defense count vs infections (Fig 4c–d)** (`scripts/fig4_cd.R`): Gaussian LMM, phylogeny as random effect, 1,000,000 iterations, burn-in 50,000, thin 250.

| | Here (shipped counts) | Here (our defense counts) | Paper |
|---|---|---|---|
| Lytic interactions per strain (c) | β −0.233, 95 % CI −0.796 to 0.274, pMCMC 0.40 | β −0.235, pMCMC 0.39 | β −0.22, CI −0.74 to 0.34, pMCMC 0.41 |
| Average MLC upon infection (d) | β −0.033, CI −0.050 to −0.015, pMCMC 0.0003 | β −0.037, pMCMC 0.0003 | β −0.04, CI −0.05 to −0.01, pMCMC 0.004 |

Both estimates fall inside the paper's intervals. Panel d is significant in both, more strongly here than in the paper.

## Deviations from the original pipeline

1. The cross-validation group file of the authors is no longer public. Groups were rebuilt from the tree. Against the authors' fold logs (50 phages whose row counts agree), 50 of 14,994 rebuilt groups (0.3 %) are split across their folds.
2. The public matrix holds MLC scores 0–4; the original script used a private N/P/F/U matrix. Lytic = score above 0 gives 20.7 % positives (the README says 19 %).
3. scikit-learn 1.9.1 instead of 1.1.2.
4. The original keeps Windows paths and private files, so it was rewritten as `scripts/predict.py`. Training-set metrics, pickled models and self-computed importances were dropped.
5. Model selection by mean AUPR follows the script and the README, although the paper text says AUROC.
6. Fig 6: statistics only; Fig 4: Rmd bugs fixed and a prior added for five phages (above).
