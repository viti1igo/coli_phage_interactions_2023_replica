# Report 13: fold tables, the authors' and ours

Both tables give, for each (phage, strain) pair, the fold of the 10-fold group cross-validation that held the pair out, so the two partitions can be compared pair by pair. Our table, [`our_folds.csv`](../results/audit/our_folds.csv), is in this repository. The authors' table, `authors_folds.csv`, is derived from the authors' shipped prediction files and is not redistributed here; the section below says how to rebuild it.

## `authors_folds.csv` (the authors' folds)

- Columns: `phage`, `bacteria`, `fold` (0 to 9). 38,378 rows, 96 phages, 402 strains. Not in this repository.
- Source: the authors' own prediction files, `dev/predictions/results/predictions/predictions_<phage>_core_features_Group10Fold_CV.csv` (96 files, from https://github.com/mdmparis/coli_phage_interactions_2023). Each of those has one row per model, dataset (train or test), fold and strain. The pair's fold is the `fold` of its `test` rows, which the authors' script wrote with the strain name (`bacteria`).
- Code: `P.published()` in `scripts/predict.py` (lines 133 to 144) reads the files and keeps the rows of the best model by mean test AUPR; `authors_folds()` in `scripts/audit_gaps.py` (lines 123 to 127) writes the three columns. After you fetch the authors' repository (`scripts/download_data.sh`), rebuild it with `.venv/bin/python scripts/audit_gaps.py authors_folds`.
- The fold is read straight from the strain name. It does not depend on any row-index mapping. (The older check in `scripts/make_groups.py`, lines 27 to 44, uses the authors' `logs_*_Group10Fold_CV.csv` files instead, whose `test_idx` holds row numbers; only that check needs the public matrix rows to line up, and it uses the 50 phages where the row count agrees.)
- Use: `predict.py` takes it through `FOLDS=results/audit/authors_folds.csv` to rerun the models on the authors' partition (`scripts/run_audit_runs.sh`, `scripts/run_pinned_env.sh`).
- The authors' group file that produced these folds, `370+host_cross_validation_groups_1e-4.csv`, is not published, so the group of each strain is unknown.

## `our_folds.csv` (our folds)

- Columns: `phage`, `bacteria`, `group`, `fold_sklearn_1_9_1`, `fold_sklearn_1_1_2`. 38,435 rows, 96 phages, 402 strains, 301 groups..
- Source: our rebuilt groups `data/derived/cv_groups_1e-4.csv` (made by `scripts/make_groups.py`, single linkage at 1e-4 substitutions per site on the shipped tree `370+host_tree_root=albertii.nwk`) and the shipped feature set, with `GroupKFold(10)` called once per phage exactly as `scripts/predict.py` does (`P.design(ph)`).
- The two fold columns come from the same call in two environments: scikit-learn 1.9.1 (`.venv`, from `scripts/audit_gaps.py folds`) and scikit-learn 1.1.2 (`ml_pinned`, from `scripts/fold_version_check.py`). The 1.9.1 column equals the independent 1.9.1 run of `fold_version_check.py` for every pair (checked in the build).
- Code: `scripts/build_our_folds.py` merges them; run `.venv/bin/python scripts/build_our_folds.py` after `scripts/audit_gaps.py folds` and `scripts/fold_version_check.py` (both write the inputs).
- This table holds only our own results.

## What the comparison shows

| Comparison | Value |
|---|---|
| Pairs in both tables | 38,377 (1 pair only in the authors', 58 only in ours: the pairs the authors' predictions lack) |
| Same fold index, ours (1.9.1) against the authors' | 14.3 % |
| Same fold index, ours (1.1.2) against the authors' | 16.9 % |
| Same fold index, our 1.9.1 against our 1.1.2 | 12.3 % |
| Mean adjusted Rand index per phage, ours (1.9.1) against the authors' | 0.022 |
| Mean adjusted Rand index per phage, ours (1.1.2) against the authors' | 0.030 |
| Mean adjusted Rand index per phage, our 1.9.1 against our 1.1.2 | 0.047 |

- `GroupKFold` is not random. For equal-size groups the order comes from the order of the group labels, so the same partition with renumbered labels gives a different fold assignment (adjusted Rand index about 0.02 between renumberings). Our labels differ from the authors' and theirs are lost, so the folds differ even where the groups agree. The steps are the `np.unique` and `np.argsort` lines of `GroupKFold._iter_test_indices` (lines 521 and 534 in scikit-learn 1.1.2, 616 and 640 in 1.9.1); [report 12](12_gap_audit.md) discusses the effect.
- One of our groups is split by the authors' folds: group 112, `IAI15` and `IAI17`. Their tree distance is 0, yet the authors' folds put them in different folds (5 and 6) for all 96 phages. With that pair taken apart, their folds equal greedy `GroupKFold` packing of our groups for all 96 phages ([`fold_mechanism.csv`](../results/audit/fold_mechanism.csv)). The cause is unknown.
