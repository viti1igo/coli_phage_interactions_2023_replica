# 1. Overview and plan

## Goal

Reproduce the strain-level phage–host prediction of Gaborieau et al., *Nature Microbiology* 9, 2847–2861 (2024) (repository `mdmparis/coli_phage_interactions_2023`), starting from the raw data wherever the data allow it, and report every difference from the paper. Two further questions were added along the way: which of the model's genomic features is unnecessary, and which figures of the paper can be redrawn.

## What was done

| Phase | Content | Report |
|---|---|---|
| 1–2 | Download the authors' data, rebuild the lost cross-validation groups, re-run the per-phage models (Fig 5), the cocktail statistics (Fig 6) and the per-phage mixed models (Fig 4) | [02](02_model_reproduction.md) |
| 3 | Feature-necessity analysis: remove each feature group under the authors' own cross-validation | [04](04_feature_necessity.md) |
| 5–7 | Regenerate every model input from the genomes and raw plaque scores with the paper's tools, and compare with the shipped tables | [03](03_feature_regeneration.md) |
| 8 | Re-run the models and the ablation on the regenerated inputs | [02](02_model_reproduction.md), [04](04_feature_necessity.md) |
| 9 | Paper-vs-shipped-vs-replication tables and figures, configuration differences | [05](05_paper_comparison.md), [07](07_configuration.md) |
| 10 | Redraw the paper's figures and extended-data figures | [06](06_figures.md) |

The task list with the outcome of each step is in [09_task_log.md](09_task_log.md). The code audit and simplification are in [08_code_audit.md](08_code_audit.md).

## Facts established before any code was written

- The authors' prediction script uses private Windows paths and a cross-validation group file that is no longer public. Groups are rebuilt from the core-genome tree (single linkage at 1e-4 substitutions per site).
- The public interaction matrix holds scores 0–4 (minimum lytic concentration, MLC; missing = not testable). The original script used a private N/P/F/U matrix. A pair counts as lytic when the score is above 0 (20.7 % of tested pairs; the authors' README says 19 %).
- The final model uses UMAP coordinates of the core-genome tree, O-type, sequence type, LPS outer-core type, Klebsiella capsule type, ABC capsule serotype and host-sharing features. Defense systems are not in the final model.
- Already done by the authors (and therefore skipped here): defense systems are not required (Fig 4 and the section "Defense systems are not required to predict phage–bacteria interactions"), and adsorption factors plus phylogeny are sufficient (Fig 5). The pre-model screening used per-phage mixed models, not an ablation of the final model.
- Not found in the paper: dropping each of the six final feature groups under group cross-validation. That is the new analysis in [04](04_feature_necessity.md). The Supplementary Text and Source Data were not available, so the paper may contain related analyses there.
- The original code has a quirk that is reproduced faithfully: `same_ABC_as_host` compares the ABC serotype with itself, so it is true exactly when the strain has an ABC serotype. Its regular expression also drops the H-type, so there is no `H_host`.

## Decisions

- Small scripts, pandas and scikit-learn only, no framework; the ablation reuses the reproduction function with a feature-group blacklist.
- One configuration record per source: `config/paper_config.yaml` (from the paper, with line numbers) and `config/config_diff.csv` (all differences).
- Each preprocessing step writes its own result and is compared with the shipped table before any modelling, so a failed step is visible early.
- Tools are pinned to the paper's versions through micromamba. Where a tool could not run (SRST2 needs reads, VIRIDIC is not installed), an approximation is used and listed as a deviation.

## Risks that were identified and how they were handled

| Risk | Handling |
|---|---|
| Rebuilt CV groups differ from the lost original | Checked against the authors' fold logs: 0.3 % of groups are split across their folds |
| scikit-learn 1.9 instead of 1.1.2 | Compare distributions, not single numbers; seeds fixed |
| Public matrix differs from the private one | Documented; all positives defined as score > 0 |
| Tools do not install on an arm64 Mac | A spike installed and timed all nine tools first; osx-64 builds run under Rosetta |
| UMAP and tree-to-distance not specified | Calibrated against the shipped coordinates; neighbourhoods compared, not coordinates |
| Manual steps (LPS clustering, RBP curation, plaque annotation) | Replicate the automatic part, report agreement, mark the manual part |
