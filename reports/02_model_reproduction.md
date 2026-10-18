# 2. Model reproduction (Fig 5, Fig 6, Fig 4)

## Method

Each of the 96 phages gets five models: logistic regression with an L1 and with an L2 penalty, random forests of depth 3 and depth 6 (250 trees each), and a constant dummy model. Training uses group 10-fold cross-validation. Strains closer than 1e-4 substitutions per site on the core-genome tree always share a fold. Class weights follow the prevalence of lytic interactions, and the best model per phage is the one with the highest mean test AUPR, as in the authors' script. A pair is lytic when its MLC score is above 0.

`scripts/predict.py` runs this with `FEATURE_SET=shipped` (the default), `regenerated` or `regenerated+matrix`. Results go to `results/repro`, `results/repro_regenerated` and `results/repro_regenerated_matrix`. Every random seed is fixed, including the one for the `saga` solver, so a rerun gives identical numbers.

| Input set | Features | Cross-validation groups | Labels |
|---|---|---|---|
| `shipped` | the authors' `picard_collection.csv` and UMAP | rebuilt from the shipped tree | shipped matrix |
| `regenerated` | tables built from the genomes ([03](03_feature_regeneration.md)) and the UMAP of our tree | rebuilt from our tree | shipped matrix |
| `regenerated+matrix` | as above | as above | our matrix from the raw plaque scores, where raw data exist |

## Fig 5: per-phage and pooled performance

| | Mean AUROC per phage | Mean AUPR per phage | Pooled AUROC | Pooled AP | F1 (threshold 0.5) | Mean AUROC per strain |
|---|---|---|---|---|---|---|
| Published article | - | - | 0.86 | 0.60 | 0.60 | 0.85 in the text; 0.872 in its Source Data |
| Preprint | 0.77 | - | 0.86 | - | - | - |
| Authors' shipped predictions, recomputed here | 0.781 | 0.586 | 0.857 | 0.601 | 0.596 | 0.881 |
| This run, shipped inputs | 0.787 | 0.591 | 0.856 | 0.598 | 0.599 | 0.878 |
| This run, regenerated inputs | 0.765 | 0.569 | 0.855 | 0.596 | 0.592 | - |
| This run, regenerated inputs and our matrix | 0.769 | 0.565 | 0.853 | 0.588 | 0.589 | - |

- The pooled numbers (AUROC, average precision, F1) agree with the paper within 0.004 on the shipped inputs, and they stay within 0.004 on the regenerated inputs.
- Fig 5a of the published article plots one AUROC per strain, and the text gives 85 % as the average. The paper's own Source Data for that panel averages 0.872 (median 0.892). This run gives 0.881 (median 0.895) on the same 315 strains, and 0.878 over the 390 strains with both classes. The authors' shipped predictions give 0.881. Per-strain AUROC correlates with the Source Data at 0.84, and the 85 % in the text sits slightly below the mean of its own data. The preprint's 77 % is the per-phage mean, and that matches (0.787 here).
- The per-phage mean loses about two points when the inputs are regenerated (0.787 to 0.765). The correlation of per-phage AUROC with the authors' results falls from 0.92 to 0.79 (0.81 with our matrix). On shipped inputs the model chosen for a phage belongs to the same family as the authors' (forest or logistic regression) for 76 % of phages.
- Two code fixes: a phage whose lytic strains all fall in one test fold now gets a constant prediction instead of an error, and `penalty="l1"` became `l1_ratio=1` (the sparse solution was checked to be identical).

The panels are `fig5_a_c_d` and `results/repro/fig5_curves.png` ([06](06_figures.md)).

## Fig 6: cocktail experiment

The statistics of the authors' notebook were re-run on their shipped experimental results. 24 strains received a tailored cocktail and 76 a generic one. The recommender itself was not re-run here, because its code is not in the repository (a separate re-implementation is described in [05](05_paper_comparison.md)).

| | Here | Paper |
|---|---|---|
| Lysis, tailored | 91.67 % | 91.67 % |
| Lysis, generic | 78.95 % | 78.95 % |
| Lysis, baseline (100 strains) | 81 % | 81 % |
| MLC tailored vs baseline (Mann-Whitney) | U = 399.5, P = 0.019 | U = 399.5, P = 0.02 (text) |
| MLC tailored vs generic | U = 1171.5, P = 0.034 | P = 0.02 in the text, P = 0.034 in the Fig 6e legend |

Everything matches. For the last test the paper contradicts itself (0.02 in the text, 0.034 in the legend), and the result here agrees with the legend.

## Fig 4: per-phage mixed models and the defense count

The per-phage models (`scripts/fig4_mcmc.R`, `scripts/fig4_summarize.py`) use `MCMCglmm` with a threshold family. Adsorption factors and defense systems are the covariates, and the phylogeny enters as a random effect through the inverse relatedness matrix. Each fit runs 500,000 iterations with a burn-in of 50,000, thinning of 250 and a prior of V = 1, ν = 0.002. Significance is pMCMC < 1e-3 (the paper's Source Data uses 0.01, see below).

91 of 96 phages fitted. Five (409_P3, 411_P1, 411_P2, 412_P1, 423_P10) failed with "mixed model equations singular". They were refitted with an N(0, 10²) prior on the fixed effects, which is a labelled deviation (`results/fig4/STRONG_PRIOR_PHAGES.txt`). All five have very few lytic strains (3 to 43).

| Result over 96 phages | Count (pMCMC < 1e-3) |
|---|---|
| Adsorption-factor terms | 28 (LPS 14, ABC capsule 12, O-type 1, OMPA 1) |
| Defense-system terms | 0 |

The direction matches the paper: adsorption factors matter and defense systems do not. The paper's Source Data for Fig 4a-b counts 52 adsorption-factor terms and one defense term, with pMCMC up to 0.0089. Its threshold is therefore 0.01, not the 1e-3 stated in the Methods. At 0.01 this re-run gives 51 and 3. The paper's single defense term is the outer-membrane variant YNCDOther, which its Source Data files under defense systems. By family at 0.01 (paper, then here): ABC capsule 25 and 21, LPS 17 and 22, O-type 7 and 6, Klebsiella capsule 2 and 0, OMPA 1 and 1.

The signs of the shared terms agree (31 of 31), but only 19 of the paper's 53 (phage, term) pairs are significant in both, because terms near the threshold change from run to run. Against the authors' shipped results (adsorption factors only, 96 phages), 7,776 (phage, term) pairs are shared. Both are significant for 16, only the shipped result for 10 and only this re-run for 12. The Spearman correlation of pMCMC is 0.625.

Two problems turned up in the authors' code and text:
- Their Rmd keeps the outer-membrane protein table under the original strain names, while the other tables went through `make.names()`. The `merge` therefore drops 136 of 403 strains without a warning. The R script here fixes it.
- The Methods say traits with fewer than 8 observations were dropped. The Rmd has no such filter: it lumps rare levels at 15 and filters defense systems at 15. No filter was added here.

For the defense count against infections (Fig 4c-d, `scripts/fig4_cd.R`), a Gaussian mixed model with the phylogeny as random effect runs 1,000,000 iterations with a burn-in of 50,000 and thinning of 250.

| | Here, shipped matrix and counts | Here, our defense counts | Here, regenerated matrix (our counts) | Paper |
|---|---|---|---|---|
| Lytic interactions per strain (c) | β -0.223, 95 % CI -0.729 to 0.303, pMCMC 0.40 | β -0.231, pMCMC 0.41 | β -0.193, pMCMC 0.51 (369 strains) | β -0.22, CI -0.74 to 0.34, pMCMC 0.41 |
| Average MLC upon infection (d) | β -0.033, CI -0.049 to -0.015, pMCMC 0.0005 | β -0.037, pMCMC 0.0003 | β -0.034, pMCMC 0.0003 (357 strains) | β -0.04, CI -0.05 to -0.01, pMCMC 0.004 |

Every estimate falls inside the paper's intervals. Panel d is significant in all runs, and more strongly here than in the paper. With the smaller defense counts from the paper's own Source Data for panel d, the shipped matrix gives β -0.032 (pMCMC 0.0005) and the regenerated matrix gives β -0.028 (pMCMC 0.0021), the closest to the paper's 0.004.

## Deviations from the original pipeline

1. The authors' cross-validation group file is no longer public, so the groups were rebuilt from the tree. Against the authors' fold logs (50 phages whose row counts agree), 50 of 14,994 rebuilt groups (0.3 %) are split across their folds.
2. The public matrix holds MLC scores from 0 to 4, while the original script used a private N/P/F/U matrix. Calling a pair lytic when its score is above 0 gives 20.7 % positives (the README says 19 %).
3. scikit-learn is 1.9.1 instead of 1.1.2.
4. The original keeps Windows paths and private files, so it was rewritten as `scripts/predict.py`. Training-set metrics, pickled models and self-computed importances were dropped.
5. The best model per phage is chosen by mean AUPR, as in the script and the README, although the paper's text says AUROC.
6. For Fig 6 only the statistics were repeated. For Fig 4 the Rmd bugs were fixed and a prior was added for five phages.
