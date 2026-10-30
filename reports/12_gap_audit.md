# 12. Gap audit

This report gives the reason for each gap between the paper, the authors' repository and this replication. Five sources are compared, as tiers:

| Tier | Source | Where |
|---|---|---|
| T0 | The paper: text, supplement and Source Data | `docs/paper/` (local only, never pushed) |
| T1 | The authors' repository: shipped predictions and per-phage performance files | `dev/predictions/results/` |
| T2 | This replication on the shipped inputs | `results/repro/` |
| T3 | This replication on regenerated features (own tree, UMAP, groups) | `results/repro_regenerated/` |
| T4 | T3 plus the interaction matrix regenerated from the raw plaque scores | `results/repro_regenerated_matrix/` |

Three comparisons run through the report: the paper against the repository (T0 against T1), our run on the shipped inputs against both (T2), and the shipped inputs against the regenerated ones (T2 against T3 and T4). Each is made overall, per phage (96 phages) and per strain-phage pair (38,377 pairs present in every tier). Every gap gets one cause code and a pointer to the file that holds the evidence:

| Code | Meaning |
|---|---|
| P | Published materials disagree: the paper, supplement, Source Data or repository contradict each other or the Methods |
| D | A definition or unit is not stated, so the number cannot be reproduced |
| U | An input or code is not published (CV groups file, recommender code, test-strain genomes, defense-subfamily clustering) |
| M | The input data differs (a strain or label) |
| F | Feature regeneration differs from the shipped table, with a measured effect |
| S | Resampling noise: a different fold partition, model randomness or a small number of lytic strains |
| C | The model-choice rule or the model picked differs |
| I | A deliberate choice on this side, for example NODF in place of BINMATNEST |
| OK | Flagged only because the paper gives no number or a qualitative statement; the numbers agree |
| X | Unexplained; the row says what was tested and what is missing |

A gap is coded X when the evidence in hand does not decide it. The code is never left empty. The register has 45 rows from the paper comparison ([`results/audit/gap_register.csv`](../results/audit/gap_register.csv)), 19 rows for the prediction results ([`gap_register_prediction.csv`](../results/audit/gap_register_prediction.csv)) and 28 rows for the paper against the repository ([`paper_vs_repo.csv`](../results/audit/paper_vs_repo.csv)). All of it comes from `scripts/audit_gaps.py`.

## Short answer

- Overall performance agrees. On the shipped inputs our pooled AUROC is 0.856 against 0.857 for the authors' predictions and 0.86 in the paper, and every overall metric of T2 lies inside the bootstrap interval of T1.
- The per-phage differences between our run and the authors' come from the cross-validation partition. The authors' groups file is not published. Our rebuilt groups sit inside one of their folds in 99.7 % of cases, yet the fold assignments barely overlap (adjusted Rand index 0.02). With their folds the mean absolute change in per-phage AUROC falls from 0.026 to 0.0065, and 34 of 96 phages match to six decimals.
- Pinning the environment to the paper (scikit-learn 1.1.2, scipy 1.11.2, Python 3.10 because 1.1.2 has no 3.12 build) and using the authors' folds gives a mean absolute change of 0.0042 in per-phage AUROC, and 72 % of the LogReg L2 probabilities equal theirs to six decimals. Random forests cannot match (the authors left them unseeded). With our own groups the version does not close the gap. [Details](#pinned-environment-t21-and-t31).
- Regenerating the features costs 0.018 mean per-phage AUROC with the folds fixed and leaves pooled AUROC almost unchanged (0.855 and 0.853 against 0.856). Most of the loss sits in 7 phages with fewer than 10 lytic strains. The O-type and Klebsiella capsule type carry it once the other groups are regenerated.
- The 45 rows of the paper comparison that did not match split as D 18, P 10, X 8, U 4, OK 2, I 2, F 1. The paper and its materials disagree with themselves in 10 rows, and a definition or unit is missing in 18.
- Of the 28 differences between the repository and the paper, 20 are inconsistencies between published materials (code P), 4 lack a definition, 3 rest on unpublished code or data and 1 is the model-selection rule.
- An independent preprint on the same matrix ([check](#check-against-an-independent-preprint)) agrees on the dataset counts and on overall AUROC (0.869 against 0.86), and its baseline supports our reading of ED7.
- 8 rows of the paper comparison stay unexplained, and 7 of the 58 per-phage rows lie beyond the noise estimates without a tested cause. [Open rows](#open-rows-and-what-would-close-them) lists them and what each needs.

## Overall performance

All five metrics below pool the held-out predictions of the 10 folds and use the 38,377 pairs that every tier has. The authors' script averages the 10 per-fold values instead, which gives the 0.787 and 0.765 of the earlier reports. The bootstrap resamples strains 500 times and keeps each phage's best model fixed.

| Metric | Paper (T0) | T1 | T2 | T3 | T4 |
|---|---|---|---|---|---|
| Pooled AUROC | 0.860 | 0.857 | 0.856 | 0.855 | 0.853 |
| Pooled average precision | 0.60 | 0.601 | 0.599 | 0.597 | 0.589 |
| Mean per-phage AUROC | 0.770 | 0.775 | 0.778 | 0.757 | 0.752 |
| Mean per-phage AUPR | - | 0.487 | 0.486 | 0.474 | 0.465 |
| Mean per-strain AUROC | 0.850 | 0.881 | 0.879 | 0.877 | 0.877 |

Gaps against the authors' predictions, with the 95 % bootstrap interval:

| Metric | T2 against T1 | T3 against T1 | T4 against T1 |
|---|---|---|---|
| Pooled AUROC | -0.0010 (-0.0031 to +0.0011) | -0.0018 (-0.0061 to +0.0026) | -0.0034 (-0.0079 to +0.0012) |
| Pooled average precision | -0.0020 (-0.0090 to +0.0057) | -0.0041 (-0.0174 to +0.0077) | -0.0124 (-0.0254 to -0.0004) |
| Mean per-phage AUROC | +0.0033 (-0.0036 to +0.0098) | -0.0180 (-0.0320 to -0.0018) | -0.0229 (-0.0372 to -0.0078) |
| Mean per-phage AUPR | -0.0007 (-0.0086 to +0.0069) | -0.0134 (-0.0257 to +0.0011) | -0.0226 (-0.0355 to -0.0068) |
| Mean per-strain AUROC | -0.0021 (-0.0048 to +0.0009) | -0.0039 (-0.0096 to +0.0012) | -0.0039 (-0.0092 to +0.0015) |

T2 sits inside the interval of T1 on every metric, so the shipped-input replication has no overall gap to explain. T3 and T4 lose 0.018 and 0.023 mean per-phage AUROC, and the intervals exclude zero. Pooled AUROC stays inside its interval. The paper's mean per-strain AUROC of 0.85 is below its own Source Data (0.872) and below the authors' predictions (0.881), so that gap is in the text (P).

## Per phage

[`per_phage.csv`](../results/audit/per_phage.csv) holds the 96 phages by tier: AUROC, its standard error over folds, AUPR, best model and pooled AUROC. The AUROC computed from the fold-level predictions equals the earlier `per_phage.csv` of each run exactly.

- T2 against T1: mean absolute change 0.026, 13 phages above 0.05.
- T3 against T2: mean absolute change 0.040, 20 phages above 0.05.
- T4 against T2: mean absolute change 0.044, 25 phages above 0.05.
- Phages with fewer lytic strains move more (Spearman -0.37 between the number of lytic strains and the size of the change).

[`per_phage_causes.csv`](../results/audit/per_phage_causes.csv) codes every flagged phage (58 phage-comparison rows). Each row carries a `basis` sentence and an `evidence` column that says how strong the support is:

| Comparison | Code | Evidence | Rows |
|---|---|---|---|
| T2 vs T1 | U | tested: authors' folds (T72a) | 13 |
| T3 vs T2 | F | tested: feature swap (T73) | 9 |
| T3 vs T2 | F | tested: feature swap, partial (T73) | 2 |
| T3 vs T2 | S | noise only: no cause shown | 7 |
| T3 vs T2 | X | none: no tested factor is sufficient | 2 |
| T4 vs T2 | F | tested: feature swap (T73) | 8 |
| T4 vs T2 | F | tested: feature swap, partial (T73) | 2 |
| T4 vs T2 | S | noise only: no cause shown | 10 |
| T4 vs T2 | X | none: no tested factor is sufficient | 5 |

Read the levels this way. The first two are tests: removing the factor removes the change (the authors' folds), or putting the factor in reproduces at least half of it (the feature swaps). Rows coded S are weaker. It says only that the change lies inside both noise estimates (the standard error over folds and the spread of five random partitions), so it cannot be told from noise, and no cause was shown. The code S does not mean a measured seed effect. The rows coded X lie beyond both noise estimates and no tested factor is sufficient: 55989_P1 (T3 against T2); LI10_P5 (T3 against T2, T4 against T2); LI10_P4 (T4 against T2); LM02_P1 (T4 against T2); MT1B1_3A1 (T4 against T2); NAN33_P4 (T4 against T2). The spread of random partitions comes from only five runs per phage, so it is a rough estimate, and a few of these X rows may still be noise. For BCH953_P2 and MT1B1_3A1 I re-scored the regenerated-matrix predictions against the shipped labels and got the same AUROC, which rules out the relabelled pairs.

The phages coded F for T3 against T2:

| Phage | Lytic strains | Change, T3 against T2 | Change with folds fixed, all groups | Largest single group |
|---|---|---|---|---|
| 411_P1 | 3 | -0.500 | -0.654 | Klebsiella (-0.052) |
| 411_P2 | 4 | -0.222 | -0.250 | Klebsiella (-0.188) |
| 412_P4 | 18 | -0.070 | -0.088 | O-type (-0.129) |
| 412_P5 | 39 | -0.050 | -0.054 | O-type (-0.053) |
| 416_P5 | 18 | -0.059 | -0.031 | ABC (-0.038) |
| 55989_P2 | 214 | -0.053 | -0.037 | UMAP (-0.020) |
| AN24_P2 | 5 | -0.056 | -0.101 | UMAP (-0.075) |
| BCH953_P2 | 138 | -0.063 | -0.055 | O-type (-0.047) |
| BCH953_P3 | 9 | -0.105 | -0.136 | UMAP (-0.073) |
| LF31_P3 | 42 | -0.053 | -0.050 | ABC (-0.039) |
| LF7074_P3 | 34 | -0.087 | -0.100 | Klebsiella (+0.020) |

## Per pair

Pair-level comparisons use the 38,377 pairs present in T1 and T2. The authors' predictions lack 58 pairs that T2 has (16 strains, 38 of them for NILS81) and have 1 pair that T2 lacks. Labels agree on every shared pair. The authors' script turns the label F into NaN in an N/P/F/U matrix that the repository does not ship, so the cause of the 59 pairs is a reading I could not verify (code U). T4 adds 135 pairs and changes 72 labels.

- T2 against T1: Spearman 0.927 over all pairs (median 0.894 per phage, lowest 0.536). Calls at 0.5 agree for 94.5 % of pairs (2,113 flips, 1,044 to lytic and 1,069 to non-lytic).
- The flip rate is 3.7 % when the best model is identical, 4.6 % for another model of the same family and 10.3 % for another family. The first figure is the floor set by partition and randomness.
- With the authors' folds the flip count falls to 657 (1.7 %). Where the best model is the same only 0.5 % flip, and where it differs 6.3 %. 20 phages pick another best model, because the choice is a close race in AUPR.
- T3 and T4 against T2: 8.2 % and 8.1 % of calls flip (Spearman 0.88 and 0.88). The flip rate is the same whether or not the strain's own categorical features changed, so those five columns do not explain it. With the folds fixed, the regenerated features alone flip 8.0 % of calls, led by UMAP (5.8 %) and the O-type (4.4 %). The regenerated tree gives different UMAP coordinates for every strain, which fits UMAP moving the most calls.
- Of the 72 relabelled pairs in T4, 12 flip.

## The fold partition

The authors fix `np.random.seed(0)` in the prediction script and say nothing about seeds elsewhere. The cross-validation groups file `370+host_cross_validation_groups_1e-4.csv` is "not available anymore" according to their README. I rebuilt the groups from the tree. Our groups agree with their folds (a group falls in a single authors' fold in 99.7 % of cases), yet the two sets of folds overlap little (adjusted Rand index 0.022 on average, 0.030 at most).

Two tests say where the difference comes from ([`fold_mechanism.csv`](../results/audit/fold_mechanism.csv)). First, the authors' folds are what greedy `GroupKFold` packing of groups like ours produces. That rule sorts groups by size and puts each into the lightest fold, so only the order among equal-size groups is free, and most of our 301 groups are single strains. One of our groups (2 strains) is split by their folds. With those two strains taken apart, the authors' folds fit the packing rule for every phage (96 of 96, all 672 size blocks). Without that change, none of the 96 phages fits. Second, the preprint, the paper and the README all name scikit-learn 1.1.2, and this run uses 1.9.1. A copy of the 1.1.2 algorithm run on our group labels overlaps the authors' folds no better (adjusted Rand index 0.028, against 0.022 for the installed version). The library version is therefore not the cause. What is missing is the order among equal groups, which comes from the labels in the unpublished file. This shows the authors' folds are reachable from groups like ours, and it does not prove that the order alone explains every difference.

Test T72a: the authors' prediction files carry the fold index of every pair. I added switches to `predict.py` (`FOLDS`, `SEED`, `SWAP`, `TAG`; all off by default, and a run with none set reproduces the earlier T2) and ran the shipped inputs on their folds, restricted to the pairs both runs have.

- Mean absolute change in per-phage AUROC against the authors' values: 0.026 with our folds, 0.0065 with theirs.
- 34 of 96 phages equal the authors' AUROC to six decimals, 76 have the same best model and 2 differ by more than 0.05.
- The two that differ by more than 0.05 are 409_P5, LI10_P4. Both picked another best model than the authors did.
- Mean per-phage AUROC 0.7789 against 0.7812; pooled AUROC 0.8559 against 0.8569.

The remaining differences come from the best-model race, the unseeded forests and the scikit-learn version (see [the pinned environment](#pinned-environment-t21-and-t31)). Run-to-run noise measured on the same inputs (test T72, five seeds each):

| Source of variation | Mean per-phage AUROC | Sd over runs | Range over runs | Mean per-phage sd | Largest per-phage sd |
|---|---|---|---|---|---|
| Model seed, folds fixed (5 seeds) | 0.7801 | 0.0011 | 0.7789 to 0.7810 | 0.0050 | 0.048 |
| Seed and a random group partition (5 seeds) | 0.7815 | 0.0031 | 0.7772 to 0.7846 | 0.0274 | 0.143 |

A random partition moves a phage's AUROC about five times as much as a model seed. Of the 13 phages that move by more than 0.05 between T2 and T1, 2 lie beyond the spread of random partitions, and none stays above 0.05 with the authors' folds.

## Pinned environment (T2.1 and T3.1)

The runs above used scikit-learn 1.9.1 on Python 3.14. The paper states scikit-learn 1.1.2. I built a second conda environment (`ml_pinned`, [`config/envs/ml_pinned.yml`](../config/envs/ml_pinned.yml), built by `scripts/setup_pinned_env.sh`) with scikit-learn 1.1.2 and every other version the paper gives, and ran the same `predict.py` three times: on the shipped inputs with our groups (T2.1), on the shipped inputs with the authors' folds, and on the regenerated features with our groups (T3.1). The authors state no Python version for the model step (3.12.2 appears only in the kernel metadata of their plotting notebook, which does not show what trained the models). scikit-learn 1.1.2 has no Python 3.12 build, so the environment runs Python 3.10.22. The table lists every component and what could not be matched ([`env_match.csv`](../results/audit/env_match.csv)).

| component | authors | pinned_runs | status | reason |
|---|---|---|---|---|
| scikit-learn | 1.1.2 (Methods, statistics paragraph) | 1.1.2 | same in pinned |  |
| Python | not stated (the paper, the prediction script and the repository give no Python version for the model step; the only Python version in the repository is 3.12.2 in the kernel metadata of the plotting notebook visualize_predictions.ipynb, which does not say what trained the models) | 3.10.22 | differs | scikit-learn 1.1.2 has no Python 3.12 build (conda-forge osx-arm64 has py310 only; the release dates from 2022). The version that fixes the results is scikit-learn's, so Python follows it. The authors' Python for the model step is unknown, so 3.10 is a choice, not a match or a mismatch. The notebook's 3.12.2 cannot be one environment with scikit-learn 1.1.2, so it does not describe the model step if the paper's version is right (inference) |
| scipy | 1.11.2 (Methods, scipy.stats) | 1.11.2 | same in pinned (stated for scipy.stats, assumed for the model step) |  |
| numpy | not stated | 1.23.5 | not stated, chosen | release contemporary with scikit-learn 1.1.2; the conda build of 1.1.2 needs numpy < 2 |
| pandas | not stated | 1.5.3 | not stated, chosen | same reasoning; pandas 3 postdates scikit-learn 1.1.2 |
| joblib / threadpoolctl | not stated | 1.2.0 / 3.6.0 | not stated, chosen | joblib contemporary with scikit-learn 1.1.2; threadpoolctl is pulled in by the solver |
| LogisticRegression L1 argument | penalty="l1", solver="saga", max_iter=10000 | penalty="l1" | same in pinned | In 1.1.2, l1_ratio=1 without penalty="l1" silently fits L2 (checked); predict.py now chooses the argument by version |
| model randomness | np.random.seed(0) once at the top; RandomForest, Dummy and saga without random_state, phages in sequence, so each draw depends on every earlier draw | same as the main runs | not replicable | the stream state at any model depends on the draws of all earlier models, including DummyClassifier(stratified), whose draws depend on the authors' private labels and on row order. Measured consequence: probabilities of RandomForest-best phages differ even with their folds (median /dp/ 0.008 to 0.018), mean AUROC noise 0.0011 (T72) |
| row and column order of the training matrix | rows sorted by bacteria then phage; columns in the order of their feature file | same as the main runs | differs, not matched | bootstrap and feature draws of the forests depend on positions; could be matched for the shipped set, not attempted because the random stream (above) cannot be matched anyway |
| interaction labels | private N/P/F/U matrix (ADSORPTION_MATRIX_370+host_vs_96_less_stringent.csv), U counted as lytic, F dropped | same as the main runs | not available | the authors' matrix is not published; 59 of 38,435 pairs differ in presence (58 only in ours, 1 only in theirs); label differences are not measurable |
| CV groups and folds | 370+host_cross_validation_groups_1e-4.csv, GroupKFold(10) of scikit-learn 1.1.2 | rebuilt groups, GroupKFold(10) of 1.1.2; the authors' own folds as a second run | groups not available | the groups file is not published. The same groups give different folds under 1.1.2 and 1.9.1 (adjusted Rand index 0.047), and neither is close to the authors' (0.030 and 0.022). Their own folds are available from their prediction files and are used in the 'authors folds' runs |
| operating system, CPU, BLAS | Windows (paths D:\\ in the script); CPU and BLAS not stated | macOS arm64 (Apple M2), conda OpenBLAS | not replicable | no Windows machine; floating-point differences of about 1e-4 in the probabilities of deterministic LogReg L2 on a few phages are consistent with this (not isolated) |
| R, MCMCglmm, vegan, statsmodels | MCMCglmm 2.35, vegan 2.6-4, statsmodels 0.14.1, R version not stated | not part of this rerun | unchanged | the pinned environment covers only the model step; see config/config_diff.csv for the rest |

Results ([`pinned_env_overall.csv`](../results/audit/pinned_env_overall.csv), [`pinned_env_per_phage.csv`](../results/audit/pinned_env_per_phage.csv)):

| run | mean per-phage AUROC | pooled AUROC | mean abs change vs T1 | correlation with T1 | same model family as T1 |
|---|---|---|---|---|---|
| T1 | 0.7812 | 0.8569 |  |  |  |
| T2 | 0.7867 | 0.8556 | 0.0261 | 0.9222 | 0.7604 |
| T2.1 | 0.7835 | 0.8545 | 0.0296 | 0.905 | 0.7292 |
| T2 authors' folds | 0.7789 | 0.8559 | 0.0065 | 0.9903 | 0.875 |
| T2.1 authors' folds | 0.7813 | 0.8568 | 0.0042 | 0.9965 | 0.9271 |
| T3 | 0.7655 | 0.8548 | 0.0349 | 0.79 | 0.6875 |
| T3.1 | 0.7682 | 0.8548 | 0.0385 | 0.6903 | 0.6979 |

- With the authors' folds, scikit-learn 1.1.2 brings the mean absolute change in per-phage AUROC from 0.0065 (1.9.1) to 0.0042, the correlation with the authors' values to 0.9965, and the share of phages with the same best model family to 92.7 %. Mean per-phage AUROC is 0.7813 against 0.7812.
- The pinned environment reproduces the deterministic model exactly. For the 24 phages whose best model is LogReg L2 in both runs, 72.2 % of the pair probabilities equal the authors' to six decimals, against 1.1 % with scikit-learn 1.9.1 ([`pinned_env_pair_agreement.csv`](../results/audit/pinned_env_pair_agreement.csv)). 16 of those phages agree on every pair. The others have 398 to 400 pairs instead of 401 or 402 (the pair set differs, see T70) and differ by at most 0.0002, except LF73_P4 (0.063), which I did not trace.
- The random forests cannot match. The authors seed NumPy once (`np.random.seed(0)`) and leave the forests, the Dummy and the saga solver without `random_state`, so each model draws from a stream that every earlier model advanced. Replaying it needs their private labels and their row order. With their folds and scikit-learn 1.1.2 the median absolute probability difference is 0.008 for depth-3 forests and 0.018 for depth-6 forests, and no pair is identical. This is the part of the gap that stays after the version and the folds are matched. I did not try to replay the stream.
- With our own groups, the version does not close the gap: the mean absolute change against the authors is 0.0296 (T2.1) against 0.0261 (T2). The two scikit-learn versions assign different folds to the same groups (adjusted Rand index 0.047, the same fold index for 12 % of pairs, no phage identical), and neither assignment is close to the authors' (0.030 for 1.1.2, 0.022 for 1.9.1). Between T2 and T2.1 the AUROC of a phage moves by 0.027 on average, as much as the gap to the authors. The version and the missing groups file cannot be separated, and the fold partition stays the main cause of the per-phage gaps.
- Regenerating the features (T3.1 against T2.1) lowers the mean per-phage AUROC by 0.0153, against 0.0212 with 1.9.1; the pooled AUROC is 0.8548 in both (T3: 0.8548). The environment does not change the conclusion about the features.

## Which regenerated features matter

Test T73 puts one regenerated feature group at a time into the shipped inputs, with the authors' folds and seed 0 held fixed. The groups are the O-type, sequence type (ST), LPS type, Klebsiella capsule type, ABC capsule type, the 8-D UMAP and the host features (the ST of the isolation host and the same-as-host comparisons). A group that is swapped carries its host-derived columns with it, except for "host", which swaps only the host side. Controls: swapping nothing reproduces the fixed-fold seed-0 run (largest AUROC difference 0), and swapping everything tracks the earlier regenerated run (per-phage correlation 0.886).

| Group regenerated | Mean AUROC change | Same, phages with 10 or more lytic strains | Phages moving by more than 0.05 | Calls flipped at 0.5 |
|---|---|---|---|---|
| O-type | -0.0004 | -0.0006 | 4 | 4.4 % |
| ST | +0.0015 | +0.0001 | 2 | 2.3 % |
| LPS | -0.0032 | -0.0035 | 1 | 2.1 % |
| Klebsiella | -0.0043 | -0.0021 | 5 | 1.8 % |
| ABC | +0.0029 | +0.0022 | 6 | 3.3 % |
| UMAP | +0.0014 | +0.0017 | 9 | 5.8 % |
| host | +0.0008 | +0.0009 | 1 | 0.0 % |
| all seven together | -0.0182 | -0.0069 | 18 | 8.0 % |

The single effects are small and have both signs, so they cancel: they sum to -0.0013, while the seven together lose 0.0182. Without the 7 phages that have fewer than 10 lytic strains the combined loss is 0.0069. Phages 411_P1 (-0.65) and 411_P2 (-0.25) carry most of it, with 3 and 4 lytic strains. The host group changes only the 9 phages whose isolation host got a different O-type or ST. Test T73b regenerates everything except one group, which gives the cost of each group once the others are in place:

| Group | Cost once the others are regenerated | Same, phages with 10 or more lytic strains | Phages moving by more than 0.05 |
|---|---|---|---|
| O-type | -0.0115 | -0.0016 | 6 |
| ST | -0.0024 | -0.0010 | 2 |
| LPS | -0.0024 | -0.0013 | 1 |
| Klebsiella | -0.0124 | -0.0017 | 5 |
| ABC | +0.0001 | +0.0011 | 2 |
| UMAP | +0.0001 | +0.0006 | 12 |
| host | +0.0000 | +0.0000 | 0 |

The O-type and the Klebsiella capsule type account for the loss once the other groups are regenerated, mostly through phages with few lytic strains. The agreement of the regenerated O-type with the shipped, hand-curated column is 74.7 % ([report 3](03_feature_regeneration.md)), the lowest of the typing features, which fits.

## Paper against repository

Report [11](11_repository_vs_paper.md) lists 28 differences between the repository and the paper. [`paper_vs_repo.csv`](../results/audit/paper_vs_repo.csv) gives each one a code: P 20, D 4, U 3, C 1. The 11 that change a number the paper reports:

| Item in report 11 | Topic | Code |
|---|---|---|
| 1 | Strains in the interaction matrix | P |
| 3 | RBP genes per phage | P |
| 10 | Significance threshold of the per-phage models | P |
| 12 | Merge of the outer-membrane table | P |
| 14 | Defense count in Fig 4d | P |
| 15 | Best model per phage | C |
| 17 | Average AUROC of Fig 5a | P |
| 18 | PERMANOVA of the phage host range | D |
| 19 | ED9a test | P |
| 20 | ED7 coverage | D |
| 22 | Generic cocktails | P |

## Register of the paper comparison

Each row of [`paper_vs_replication.csv`](../results/comparison/paper_vs_replication.csv) that did not match or could not be compared has a code and a reason (D 18, P 10, X 8, U 4, OK 2, I 2, F 1). The basis column of the CSV says whether the reason comes from an earlier report, from a test in this audit, or is a residual. Four ED7 rows were tested here, and eight rows remain residual.

| Id | Row | Code | Reason |
|---|---|---|---|
| G01 | Fig 4d: pMCMC | X | beta equal (-0.032 with the panel-d counts); only pMCMC differs (0.0005 against 0.004); the panel-d counts are inconsistent with Fig 4c (P) but do not change beta. Not testable here: one chain takes about 4.5 hours and the authors' seed and chain settings are not stated; beta agrees, so only the posterior tail differs |
| G02 | Fig 4d (paper's panel-d defense counts): pMCMC | X | beta equal (-0.032 with the panel-d counts); only pMCMC differs (0.0005 against 0.004); the panel-d counts are inconsistent with Fig 4c (P) but do not change beta. Not testable here: one chain takes about 4.5 hours and the authors' seed and chain settings are not stated; beta agrees, so only the posterior tail differs |
| G03 | ED1: phages infecting BE | OK | 45 and 48 equal the paper's counts; row marked qualitative only |
| G04 | ED1: phages infecting BL21 | OK | 45 and 48 equal the paper's counts; row marked qualitative only |
| G05 | ED2: nestedness (paper: BINMATNEST temperature T; ours: NODF, not comparable) | I | NODF used instead of BINMATNEST; the paper's T and ours are not comparable |
| G06 | ED4A: AL505_Sd2: infected strains that are ABC K1 | D | paper states a qualitative association with K1, no count |
| G07 | ED4A: 412_P2: infected strains that are ABC K1 | D | paper states a qualitative association with K1, no count |
| G08 | ED4A: MT1B1_3A1: infected strains that are ABC K1 | D | paper states a qualitative association with K1, no count |
| G09 | ED5: infected strains that carry K55 | F | shipped types give 3 of 3; regenerated Kaptive types give 2 of 2 (K55 strains total 2) |
| G10 | ED7: baseline triplet | D | the paper does not name the triplet; ours is the most-covering set |
| G11 | ED9A: Mann-Whitney P (per-strain mean step, each strain once) | P | the notebook counts each strain three times (P = 3.5e-12); with each strain once P = 6.4e-5; the test unit is not stated (D) |
| G12 | ED6A: genome size vs lytic interactions (pMCMC) | D | the paper gives no value for this model; regenerated inputs give p = 0.55 against 0.12 (F) |
| G13 | Fig 1j: defense subtypes detected (columns) | D | no number stated |
| G14 | Fig 3a: PERMANOVA adj. R2, genus + isolation strain | X | 0.65 against 0.69; the repository has no PERMANOVA code. Tested adding LF110 and the old genus labels: adjusted R2 stays 0.686 and 0.682, so neither explains 0.65 (`permanova_checks.csv`); isolation strain alone 0.554 against 0.53 |
| G15 | Fig 3a: PERMANOVA adj. R2, genus alone | D | genus alone is not stated; 0.13 is the increment after isolation strain |
| G16 | Fig 3a: PERMANOVA R2 (unadjusted), genus | D | unadjusted values are not stated |
| G17 | Fig 3a: PERMANOVA R2 (unadjusted), isolation strain | D | unadjusted values are not stated |
| G18 | Fig 5c: Spearman(AUROC, n lytic) per strain | D | no value stated |
| G19 | Fig 3b: mean entropy in RBP genes vs other genes (bits) | D | qualitative claim; ours supports it (0.62 against 0.24 bits); annotation not shipped |
| G20 | Fig 3b: mean entropy of the alignment (bits) | D | no value printed |
| G21 | Fig 4d: defense count in the Source Data equals Supplementary Table 1 | P | Source Data counts differ from Supplementary Table 1 for 353 of 400 strains |
| G22 | Fig 4a: largest pMCMC among the terms the Source Data calls significant | P | Source Data keeps terms up to pMCMC 0.0089; Methods say 1e-3 |
| G23 | Fig 4a: adsorption-factor / defense terms with pMCMC < 0.001 | P | threshold mismatch; at 0.01 the counts are 51 and 3 against 52 and 1 |
| G24 | Fig 4b: (phage, term) pairs significant in both | X | 19 of 54 pairs shared with the paper's 53; candidates are the unreproducible minimum-observations filter (D), the outer-membrane merge (P) and MCMC variation (S). Not testable here: needs the minimum-observation filter and the outer-membrane clustering the repository does not contain, and another 4.5 hours of sampling |
| G25 | Fig 6: strains with tailored / generic cocktails | P | Supplementary Table 7 gives 28 and 72; the Fig 6 analysis gives 24 and 76 (three or four generic cocktails) |
| G26 | ED9b: haplotypes covered by the four generic cocktails | P | Source Data plots strain counts (31, 21, 19, 4); haplotype counts are 20, 18, 18, 4 |
| G27 | ED7: coverage at k = 1 (%) | D | no set of MLC scores gives the Source Data curve (best single phage 172 of 403 strains, any MLC > 0 gives 264); the published text (63 %) and the Source Data (67.0 %) also disagree (P); the 68 % for four phages is in the preprint only. Tested 4 MLC thresholds on the 403-strain matrix and 12 raw-plaque definitions (replicates by dilution); none reproduces the Source Data counts 172, 232, 270, 297 and plateau 360 (`ed7_definitions.csv`); the raw scores cover only 369 strains. An independent preprint reports a promiscuity baseline of about 61 % for one phage, near our 65 % and far from 42.7 % (`external_check.csv`) |
| G28 | ED7: coverage at k = 3 (%) | D | no set of MLC scores gives the Source Data curve (best single phage 172 of 403 strains, any MLC > 0 gives 264); the published text (63 %) and the Source Data (67.0 %) also disagree (P); the 68 % for four phages is in the preprint only. Tested 4 MLC thresholds on the 403-strain matrix and 12 raw-plaque definitions (replicates by dilution); none reproduces the Source Data counts 172, 232, 270, 297 and plateau 360 (`ed7_definitions.csv`); the raw scores cover only 369 strains. An independent preprint reports a promiscuity baseline of about 61 % for one phage, near our 65 % and far from 42.7 % (`external_check.csv`) |
| G29 | ED7: coverage at k = 4 (%) | D | no set of MLC scores gives the Source Data curve (best single phage 172 of 403 strains, any MLC > 0 gives 264); the published text (63 %) and the Source Data (67.0 %) also disagree (P); the 68 % for four phages is in the preprint only. Tested 4 MLC thresholds on the 403-strain matrix and 12 raw-plaque definitions (replicates by dilution); none reproduces the Source Data counts 172, 232, 270, 297 and plateau 360 (`ed7_definitions.csv`); the raw scores cover only 369 strains. An independent preprint reports a promiscuity baseline of about 61 % for one phage, near our 65 % and far from 42.7 % (`external_check.csv`) |
| G30 | ED7: coverage when all phages are used (%) | D | no set of MLC scores gives the Source Data curve (best single phage 172 of 403 strains, any MLC > 0 gives 264); the published text (63 %) and the Source Data (67.0 %) also disagree (P); the 68 % for four phages is in the preprint only. Tested 4 MLC thresholds on the 403-strain matrix and 12 raw-plaque definitions (replicates by dilution); none reproduces the Source Data counts 172, 232, 270, 297 and plateau 360 (`ed7_definitions.csv`); the raw scores cover only 369 strains. An independent preprint reports a promiscuity baseline of about 61 % for one phage, near our 65 % and far from 42.7 % (`external_check.csv`) |
| G31 | ED2: nestedness temperature T (BINMATNEST) | I | BINMATNEST not run; NODF used |
| G32 | Fig 6 recommender: identical cocktails out of 100 | U | recommender code not published |
| G33 | Fig 6 recommender: phages in common per cocktail (of 3) | U | recommender code and test-strain genomes not published |
| G34 | SI Text 4: genome-size coefficient without morphotype | D | no value given |
| G35 | SI Text 6: K1 strains in the Picard collection | P | 42 K1 strains in the repository, 45 in the text; regenerated types give 45 across the whole collection and 41 among the 369 Picard strains (F) |
| G36 | SI Text 6: K1 rule, the three phages: precision (%) | X | the K1 count (42 against 45), the rule definition and the strain set are the candidates. Tested with the regenerated K1 calls (45 strains in the whole collection): precision 89 %, recall 65 % (reports/10), so the K1 count does not explain the gap; the strain set and the definition of infection behind 73 % and 72 % are not stated |
| G37 | SI Text 6: K1 rule, the three phages: recall (%) | X | the K1 count (42 against 45), the rule definition and the strain set are the candidates. Tested with the regenerated K1 calls (45 strains in the whole collection): precision 89 %, recall 65 % (reports/10), so the K1 count does not explain the gap; the strain set and the definition of infection behind 73 % and 72 % are not stated |
| G38 | SI Text 6: K1 rule, the other 28 Vectreviruses and Kayfunaviruses: precision (%) | X | the K1 count (42 against 45), the rule definition and the strain set are the candidates. Tested with the regenerated K1 calls (45 strains in the whole collection): precision 89 %, recall 65 % (reports/10), so the K1 count does not explain the gap; the strain set and the definition of infection behind 73 % and 72 % are not stated |
| G39 | SI Text 6: K1 rule, the other 28 Vectreviruses and Kayfunaviruses: recall (%) | X | the K1 count (42 against 45), the rule definition and the strain set are the candidates. Tested with the regenerated K1 calls (45 strains in the whole collection): precision 89 %, recall 65 % (reports/10), so the K1 count does not explain the gap; the strain set and the definition of infection behind 73 % and 72 % are not stated |
| G40 | SI Text 7: phages with 0 / 1 / 2 RBP in RBP_list.csv | P | RBP_list.csv gives 7 phages no RBP; the text says every phage has one or two |
| G41 | SI Text 8: PERMANOVA adsorption factors: AIC | D | AIC depends on the feature encoding and distance, which the supplement does not state |
| G42 | SI Text 8: PERMANOVA defense subfamilies: adjusted R2 | U | the authors' gene clustering into defense subfamilies is not published |
| G43 | SI Text 8: PERMANOVA defense subfamilies: AIC | U | the authors' gene clustering into defense subfamilies is not published |
| G44 | SI Text 7: LMM beta of the RBP count on lytic interactions, all 96 phages (the 7 listed without RBP count as 0) | P | sign differs when the 7 phages count as 0; +16 against +28 with the 89 phages that have an RBP, a residual that stays unexplained (X) |
| G45 | SI Text 7: LMM beta of the RBP count on lytic interactions, 89 phages with at least one RBP | P | sign differs when the 7 phages count as 0; +16 against +28 with the 89 phages that have an RBP, a residual that stays unexplained (X) |

## Prediction results

| Id | Level | Comparison | Metric | Code | Size and evidence |
|---|---|---|---|---|---|
| P01 | overall | T2 vs T1 | pooled AUROC | OK | gap -0.0010, 95 % interval -0.0031 to +0.0011; inside the interval |
| P02 | overall | T2 vs T1 | mean per-phage AUROC | OK | gap +0.0033, 95 % interval -0.0036 to +0.0098; inside the interval |
| P03 | overall | T3 vs T2 | mean per-phage AUROC | F | gap -0.0213, 95 % interval -0.0322 to -0.0090; with the folds fixed the regenerated features move the per-fold mean AUROC by -0.0182 (the same predictions pooled differently) |
| P04 | overall | T3 vs T2 | mean per-phage AUPR | F | gap -0.0127, 95 % interval -0.0257 to +0.0006 |
| P05 | overall | T4 vs T2 | mean per-phage AUROC | F | gap -0.0262, 95 % interval -0.0371 to -0.0141; with the folds fixed the regenerated features move the per-fold mean AUROC by -0.0182 (the same predictions pooled differently); T4 also uses the regenerated matrix (72 labels change, 135 pairs added), a share that was not separated |
| P06 | overall | T4 vs T2 | mean per-phage AUPR | F | gap -0.0219, 95 % interval -0.0354 to -0.0070; T4 also uses the regenerated matrix (72 labels change, 135 pairs added), a share that was not separated |
| P07 | overall | T0 vs T2 | mean per-strain AUROC | P | paper text 0.85, its own Source Data 0.872, authors' predictions 0.881, ours 0.879 |
| P08 | phage | T2 vs T1 | phages moving by more than 0.05 AUROC | U | 13 of 13 flagged phages |
| P09 | phage | T3 vs T2 | phages moving by more than 0.05 AUROC | F | 11 of 20 flagged phages |
| P10 | phage | T3 vs T2 | phages moving by more than 0.05 AUROC | S | 7 of 20 flagged phages |
| P11 | phage | T3 vs T2 | phages moving by more than 0.05 AUROC | X | 2 of 20 flagged phages |
| P12 | phage | T4 vs T2 | phages moving by more than 0.05 AUROC | F | 10 of 25 flagged phages |
| P13 | phage | T4 vs T2 | phages moving by more than 0.05 AUROC | S | 10 of 25 flagged phages |
| P14 | phage | T4 vs T2 | phages moving by more than 0.05 AUROC | X | 5 of 25 flagged phages |
| P15 | pair | T2 vs T1 | pairs in one set only | U | 58 pairs T2 has and T1 lacks (16 strains, NILS81 38), 1 pair T1 has and T2 lacks; labels agree on all 38,377 shared pairs; cause unverified (the authors' script turns F into NaN in a matrix file that is not shipped) |
| P16 | pair | T2 vs T1 | calls flipped at 0.5 | U | 2113 of 38377 pairs (5.5 %) with our folds, 657 (1.7 %) with the authors' folds |
| P17 | pair | T2 vs T1 with the authors' folds | calls flipped at 0.5 | C | 0.5 % of pairs flip where the best model is the same, 6.3 % where it differs (20 phages pick another model in a close AUPR race) |
| P18 | pair | T3 vs T2 | calls flipped at 0.5 | F | 8.2 % of pairs; with the folds fixed the regenerated features alone flip 8.0 %, led by UMAP (5.8 %) and O-type (4.4 %) |
| P19 | pair | T4 vs T2 | calls flipped at 0.5 | F | 8.1 % of pairs; with the folds fixed the regenerated features alone flip 8.0 %, led by UMAP (5.8 %) and O-type (4.4 %) |

## Open rows and what would close them

ED7: coverage under four thresholds on the MLC score of the 403-strain matrix and twelve definitions from the raw plaque scores (replicates by dilution) was compared with the Source Data counts (172, 232, 270, 297 strains for one to four phages and a plateau of 360). None matches ([`ed7_definitions.csv`](../results/audit/ed7_definitions.csv)). The closest, MLC of at least 2, gives 175, 265, 310, 330 and a plateau of 376. The raw scores cover only 369 strains, so they cannot reproduce a 403-strain curve.

Fig 3a: adding strain LF110 or using the old genus labels leaves the adjusted R² at 0.686 and 0.682, against 0.65 in the paper. The isolation strain alone gives 0.554 against 0.53. Both are 0.02 to 0.04 above the paper, so some setting of the PERMANOVA differs from the paper's in a way the paper does not describe.

The 8 unexplained rows are G01, G02, G14, G24, G36, G37, G38, G39. G01 and G02 (Fig 4d pMCMC) agree in beta and differ only in the posterior tail. G15 is the PERMANOVA above. G26 (Fig 4b) needs the authors' filtering and clustering code. G38 to G41 (the K1 rule) are not explained by the K1 count: the regenerated calls give 45 K1 strains, as in the text, and the rule still gives 89 % and 65 % against 73 % and 72 %. They need the strain set and the definition of infection.

| What the authors could supply | Code | Register rows |
|---|---|---|
| the authors stating the definition or unit | D | G06, G07, G08, G10, G12, G13, G15, G16, G17, G18, G19, G20, G34, G41 |
| the authors' PERMANOVA call (adonis2 formula, which host and genus columns, number of phages) | X | G14 |
| the authors' chain seed, burn-in and the exact defense counts behind panel d | X | G01, G02 |
| the authors' definition of productive lysis and the matrix they used for ED7 | D | G27, G28, G29, G30 |
| the authors' filtering code and clustering of defense subfamilies | X | G24 |
| the list of the 45 K1 strains, the strain set and the definition of infection used for precision and recall | X | G36, G37, G38, G39 |
| the unpublished file or code | U | G32, G33, G42, G43 |

## The original preprint against the paper and the repository

The preprint behind the link you sent is version 1, posted 22 November 2023, which is the file already in `docs/paper/`. [`preprint_check.csv`](../results/audit/preprint_check.csv) sets 23 statements side by side: the preprint, the published paper (T0) and the authors' repository (T1). Numbers in brackets are line numbers in the extracted text. T1 values are computed from the shipped matrix, feature table and predictions.

| Item | Preprint v1 | Published paper (T0) | Repository (T1) | Verdict |
|---|---|---|---|---|
| Matrix size | 403 x 96 = 38,688 [l.242] | 403 x 96 [l.570] | 402 x 96 public matrix: 38,592 cells, 38,435 tested; the shipped predictions cover 38,378 pairs | preprint and paper agree; T1 lacks LF110 |
| Share of interactions that are lytic | 20.3 % [l.245] | 20.3 % [l.344] | 20.75 % of tested cells | T1 is 0.5 points higher; the paper's own 20.6 % and 20.8 % imply 20.7 % |
| Strains infected per phage | 83 +/- 70 (20.6 %) [l.246] | 83 +/- 70 (20.6 %) [l.99] | 83.1 +/- 69.7 (20.7 %) | same in all three |
| Phages lysing a strain | 20 +/- 15 (20.8 %) [l.247] | 20 +/- 15 (20.8 %) [l.100] | 19.8 +/- 15.2 (20.7 %) | same in all three |
| Strains lysed by no phage | 12 (3.0 %) [l.248] | 12 (3.0 %) [l.102] | 12 | same in all three |
| One-way ANOVA of host range by morphotype | F = 114.7, P < 1e-3 [l.256] | F = 114.7 [l.339]; P = 1.44e-26 in the Fig 2f legend [l.460] | F = 113.9 on the public matrix, 114.66 with LF110 | same once LF110 is added |
| PERMANOVA of phage host range (adjusted R2) | genus 0.13, isolation strain 0.53, both 0.65, defense 0.02 [l.308-311, l.432] | same [l.365-369, l.424] | no code in the repository | same in both texts; T1 has no code |
| LPS types | K-12 43, R1 206, R2 22, R3 61, R4 35 [l.206] | not in the extracted text | K-12 43, R1 206, R2 22, R3 61, R4 35 | T1 matches the preprint |
| Strains with a Klebsiella capsule | 22 of 403 (5.4 %) [l.149] | not stated | 22 of 403 | match (22/403 is 5.5 %) |
| Mean AUROC of Fig 5a | 77 %, per phage [l.521] | 85 %, per strain [l.623] | per-phage mean 0.781, per-strain mean 0.881 | T1 matches the preprint's per-phage 77 %; the paper changed the unit |
| Overall AUROC | 86 % in the Results [l.522], 85 % in the Discussion [l.672] | 86 % in the Results and the Discussion [l.625, l.824] | pooled 0.857 | the preprint contradicts itself; paper and T1 agree on 86 % |
| Average precision and F1 | not reported | 60 % and 60 % [l.625] | 0.601 and 0.596 | paper and T1 agree |
| False negatives and false positives | 3,379 and 2,922 of 38,688 [l.529] | 3,379 and 2,922 of 38,688 [l.633] | 3,318 and 2,980 of 38,378 pairs at p >= 0.5 | within 2 %; the totals of errors are 6,301 (paper) and 6,298 (T1) |
| Most covering triplet and 4-set | 63 % and 68 % [l.575-576] | 63 % [l.802]; 4-set not stated; ED7 Source Data 67.0 % and 73.7 % | no code; the matrix gives 87.6 % and 90.8 % (any lysis) | the published text drops the 68 %; the Source Data disagrees with both texts |
| K1 rule: the three phages, then the other Vectreviruses and Kayfunaviruses | precision 73 %, recall 72 %; 7 % and 14 % [l.343-345] | same, in Supplementary Text 6 | 87 % and 55 %; 18 % and 18 % | same in both texts; T1 differs |
| Tailored and generic recommendations | 15 cocktails for 24 strains, 3 cocktails for 76 [l.607-609] | same [l.783-786]; Supplementary Table 7 gives 28 and 72 | 24 / 76 | both texts match T1; Supplementary Table 7 differs |
| Cocktail success: tailored, generic, baseline | 91.67, 78.95 and 81.00 % [l.628-630] | 91.67 %; 'baseline and generic 81 % and 78 %' [l.806-811]; 78.95 % in Fig 6d | 91.67, 78.95 and 81 % | same; the published text truncates 78.95 |
| Mann-Whitney tests | U = 399.5, P = 0.02; U = 1171.5, P = 0.02 [l.631-632] | same text [l.812-815]; Source Data P = 0.019 and 0.034 | U = 399.5 and 1171.5; P = 0.019 and 0.034 | the U values agree; the text's 0.02 for tailored against generic should read 0.03 |
| Phylogroups covered by the isolation strains | 75 % [l.176] | 66 % [l.73] | not checkable | changed between versions |
| Agreement of annotated scores across replicates | 88.64 %, 9.71 %, 98.35 % [l.815-817] | 88.6 %, 9.7 %, 98.3 % | the raw table is binary: 90.7 %, 9.3 %, 100 % | not checkable: the finer scores are not shipped |
| Per-phage models | LogReg L2 and L1, RF depth 3 and 6 with 250 trees, class weights by prevalence, scikit-learn 1.1.2 [l.984-992] | same [l.1526-1562] | same four models plus a Dummy baseline (never best in T1); our rerun uses scikit-learn 1.9.1 | texts agree; the code has 5 models, the texts 4 |
| Recommender random forests | depth 4, 80 trees, class weight 1.5 [l.1014, l.1021] | same [l.1608, l.1629] | README the same | same in all three |
| Rule that picks the best model | mean AUROC [l.1002] | mean AUROC (report 11, item 15) | the script uses mean AUPR | texts say AUROC, code uses AUPR |

Where the three agree. The descriptive numbers (83 +/- 70 strains per phage, 20 +/- 15 phages per strain, 12 strains lysed by no phage, the LPS counts, 22 Klebsiella capsules) are the same in the preprint, the paper and the repository, once strain LF110 is added. Both texts describe the per-phage models, the class weights and the recommender forests in the same way. The cocktail results (91.67, 78.95 and 81 %, and the U statistics) agree with the repository.

Where the preprint differs from the paper. The 77 % AUROC is a per-phage mean in the preprint and a per-strain 85 % in the paper. The preprint gives 86 % in its Results and 85 % in its Discussion, and the paper gives 86 % in both. The coverage of 68 % for the most covering four phages appears only in the preprint. The share of E. coli phylogroups covered by the isolation strains falls from 75 % to 66 %. The paper truncates 78.95 % to 78 % in its text.

Where the repository differs from both texts. The matrix has 402 strains. The K1 rule (87 and 55 % against 73 and 72 %), the triplet coverage (87.6 % against 63 % and the Source Data's 67.0 %) and the AUPR rule for choosing the best model are the differences found before. Two are new. The paper's 20.3 % lytic share does not match its own 20.6 % and 20.8 %, and the matrix gives 20.75 %. The text's P = 0.02 for tailored against generic cocktails is 0.034 in the Source Data and in the recomputation.

What this changed in the audit. Two rows of the register said the paper gives no number, and it does: 83 +/- 70 strains per phage (20.6 %), and 3,379 false negatives and 2,922 false positives over 38,688 pairs. Both now compare as numbers and match within 5 % (20.7 % here; 3,307 and 2,940 in our run, 3,318 and 2,980 in the authors' predictions), so the paper comparison has 45 unresolved rows instead of 47. The 68 % for four phages in ED7 comes from the preprint and not from the published text, and the affected rows and reports now say so.

## Check against an independent preprint

After the audit was finished I checked it against a preprint from another group that models the same E. coli matrix: Noonan et al., [Phylogeny-agnostic strain-level prediction of phage-host interactions from genomes](https://doi.org/10.1101/2025.11.15.688630), bioRxiv, version posted 30 April 2026, not peer reviewed. I typed its numbers from the PDF (Table 1, Table 2 and the Results) and recomputed our side from the shipped matrix ([`external_check.csv`](../results/audit/external_check.csv)).

| Item | Preprint | This audit | Verdict |
|---|---|---|---|
| Strains in the E. coli matrix | 402 (Table 1) | 402 (public matrix) | agrees |
| Phages | 94 (Table 1) | 96 | differs |
| Interactions | 37,788 = 402 x 94 (Table 1) | 37,788 cells for 94 phages | agrees |
| Positive interactions | 7,813 (20.7 %) (Table 1) | 7,976 of 38,592 cells (20.7 %); 7,813 after removing two phages with 163 lytic pairs | agrees |
| Definition of a positive | any non-zero interaction (Methods) | MLC > 0 | agrees |
| AUROC for unseen strains | 0.869 (Table 2); 0.87 in the text against 0.86 for the original study | pooled AUROC 0.857 (authors), 0.856 (this run), 0.86 (paper) | agrees within 0.013 |
| Best fixed single phage, share of strains infected | about 61 % implied (66.9 % for the model's pick, 1.1-fold over promiscuity; rounding gives 58 to 64 %) | 65.4 % (promiscuity chosen on training strains, 20-fold CV, 50 repeats); ED7 Source Data says 42.7 % | supports this audit, not ED7 |
| Promiscuity baseline, 5 phages | about 81 % implied (97.5 % for the model, 1.2-fold; rounding gives 78 to 85 %) | 84.6 % (top 5 by promiscuity, no clustering) | consistent |
| Model-guided 5-phage cocktail | 97.5 % of E. coli strains | no counterpart | not comparable |
| Experimental validation | AUROC 0.84 on 1,240 interactions (52 BASEL phages x 25 ECOR strains) | no counterpart | not comparable |
| Arithmetic inside the preprint | totals 115,037 interactions, 949 strains, 27,081 and 4,304 development, 1,240 validation | recomputed from its Table 1 and Methods | reconciles, with one slip |

What it supports. The dataset counts agree (402 strains, 20.7 % positives), so both analyses read the public matrix the same way. That is a check of loading, not an independent measurement, since the file is the same. The AUROC for unseen strains (0.869) lies within 0.013 of the paper, the authors' predictions and this run, although the preprint uses another model family and another split design. That fits the finding that overall AUROC does not depend on the fold partition. For ED7 the check is more useful. The preprint's promiscuity baseline implies about 61 % of strains for one phage, and the same baseline on our matrix gives 65 %. Both are far above the 42.7 % of the ED7 Source Data. The preprint calls any non-zero interaction a positive, which is the definition used here, so the Source Data curve needs a stricter definition that nobody has stated.

What it cannot test. The preprint reports no per-phage numbers for E. coli in its main text, so it says nothing about the fold-partition effect on single phages, the feature swaps, the recommender or the K1 rule. Its cross-validation (20 folds on hierarchical clusters of strains, 10 % held out) differs from the authors' group 10-fold, and the text does not say how its AUROC is pooled. I therefore read it as a consistency check from independent code on the same public file, and not as a replication of the paper. One slip inside the preprint: Klebsiella-3 has 25 phages in the Methods and 52 in Table 1, and only 52 gives its stated total of 518.

## Limits

- The code S is a statistical statement about noise, not a measured seed effect. The z test treats the two runs as independent samples of folds, which overstates the noise for a paired comparison, so it flags fewer phages than a paired test would.
- The bootstrap intervals keep each phage's best model fixed and do not reselect it per replicate.
- Each swap run uses one seed. The seed noise measured with the folds fixed is 0.0011 on the mean per-phage AUROC, so effects of 0.003 or more on the mean are above it. Effects on single phages use that phage's own spread.
- The swaps cover the five typing columns, the UMAP and the host features. They do not include the regenerated interaction matrix, so the T4 loss beyond T3 is not attributed to features or to labels beyond the rescoring check for two phages.
- The decomposition by group is not additive, and a group's cost depends on which other groups are regenerated. The marginal table is one ordering, with everything else regenerated.
- The 59 pairs that differ between the authors' predictions and ours have a cause I could not verify.
- Eight rows of the register are residual, as listed above.
- The pinned environment matches scikit-learn and scipy to the paper, but the authors do not state Python, numpy, pandas, joblib, the operating system or the BLAS for the model step, so those were chosen, not matched. The forests' random stream was not replayed. The effect of the unmatched parts is not isolated.
- `results/audit_runs/` (88 MB of run outputs), `results/audit/authors_folds.csv` and the pair-level CSVs that hold the authors' predictions are not in this repository.

## Reproduce

```bash
.venv/bin/python scripts/audit_gaps.py register      # gap register, paper against repository
.venv/bin/python scripts/audit_gaps.py align         # tier alignment
.venv/bin/python scripts/audit_gaps.py overall       # overall metrics with bootstrap intervals (about 8 minutes)
for s in shipped regenerated regenerated+matrix; do FEATURE_SET=$s .venv/bin/python scripts/audit_gaps.py folds; done
.venv/bin/python scripts/audit_gaps.py phage         # per-phage table and causes
.venv/bin/python scripts/audit_gaps.py pairs         # pair-level comparisons
.venv/bin/python scripts/audit_gaps.py authors_folds
scripts/run_audit_runs.sh                            # 27 attribution runs, about 100 seconds each
.venv/bin/python scripts/audit_gaps.py runs          # T72a, T72, T73 summaries and refined causes
.venv/bin/python scripts/audit_gaps.py final         # ED7 and PERMANOVA tests, final registers
scripts/setup_pinned_env.sh                          # conda environment with scikit-learn 1.1.2 (T81)
scripts/run_pinned_env.sh                            # T2.1, T2.1 with the authors' folds, T3.1 (about 2 minutes each)
for s in shipped regenerated; do .venv/bin/python scripts/fold_version_check.py results/audit_runs/fold_version/sk191_$s.csv $s; done   # and with the pinned python, as sk112_$s.csv
.venv/bin/python scripts/audit_pinned.py             # T82 tables
.venv/bin/python scripts/audit_report.py             # this report
```
