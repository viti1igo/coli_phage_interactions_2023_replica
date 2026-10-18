# 1. Overview and plan

## Goal

The work reproduces the strain-level phage-host prediction of Gaborieau et al., *Nature Microbiology* 9, 2847-2861 (2024) (repository `mdmparis/coli_phage_interactions_2023`). It starts from the raw data wherever the data allow it and records every difference from the paper. Two questions were added on the way: which of the model's genomic features it can do without, and which of the paper's figures can be redrawn.

## What was done

| Phase | Content | Report |
|---|---|---|
| 1 and 2 | Download the authors' data, rebuild the lost cross-validation groups, and re-run the per-phage models (Fig 5), the cocktail statistics (Fig 6) and the per-phage mixed models (Fig 4) | [02](02_model_reproduction.md) |
| 3 | Remove each feature group in turn under the authors' own cross-validation | [04](04_feature_necessity.md) |
| 5 to 7 | Regenerate every model input from the genomes and raw plaque scores with the paper's tools, and compare each with the shipped table | [03](03_feature_regeneration.md) |
| 8 | Re-run the models and the ablation on the regenerated inputs | [02](02_model_reproduction.md), [04](04_feature_necessity.md) |
| 9 | Tables and figures comparing paper, shipped results and replication; configuration differences | [05](05_paper_comparison.md), [07](07_configuration.md) |
| 10 | Redraw the paper's figures and extended-data figures | [06](06_figures.md) |
| 11 to 15 | Use the publisher's supplement and Source Data, run everything that can be run from the authors' repository (notebooks, supplementary analyses, genome maps), and redo the remaining analyses on the regenerated inputs | [05](05_paper_comparison.md), [10](10_repository_analyses.md), [11](11_repository_vs_paper.md) |
| 16 | Annotate every figure | [06](06_figures.md) |

[09_task_log.md](09_task_log.md) has one line per task. [08_code_audit.md](08_code_audit.md) covers the code audit and the simplification.

## What was known before the first line of code

- The authors' prediction script uses private Windows paths and a cross-validation group file that is no longer public. The groups are rebuilt from the core-genome tree with single linkage at 1e-4 substitutions per site.
- The public interaction matrix holds scores from 0 to 4 (minimum lytic concentration, MLC; a missing value means the pair could not be tested). The original script used a private N/P/F/U matrix. A pair counts as lytic when its score is above 0. That gives 20.7 % of tested pairs, while the authors' README says 19 %.
- The final model uses the UMAP coordinates of the core-genome tree, O-type, sequence type, LPS outer-core type, Klebsiella capsule type, ABC capsule serotype and host-sharing features. Defense systems are not in it.
- The authors had already tested two questions, so this work skipped them: defense systems are not required (Fig 4 and the section "Defense systems are not required to predict phage-bacteria interactions"), and adsorption factors plus phylogeny are enough (Fig 5). Their screening before the model used per-phage mixed models, not an ablation of the final model.
- Neither the paper, its preprint nor its supplement describes dropping each of the six final feature groups under group cross-validation. That is the new analysis in [04](04_feature_necessity.md).
- The original code has a quirk that this work keeps so the numbers match. `same_ABC_as_host` compares the ABC serotype with itself, so it is true exactly when the strain has an ABC serotype. Its regular expression also drops the H-type, so there is no `H_host`.

## Decisions

- The code is a set of small scripts built on pandas and scikit-learn, with no framework. The ablation reuses the reproduction function with a feature-group blacklist.
- Each source of settings has one record: `config/paper_config.yaml` (the paper, with line numbers) and `config/config_diff.csv` (every difference).
- Each preprocessing step writes its own result and is checked against the shipped table before any modelling, so a failed step shows up early.
- Tools are pinned to the paper's versions through micromamba. Where a tool could not run (SRST2 needs reads, VIRIDIC is not installed), an approximation stands in and is listed as a deviation.

## Risks and how they were handled

| Risk | Handling |
|---|---|
| The rebuilt CV groups differ from the lost original | Checked against the authors' fold logs: 0.3 % of groups are split across their folds |
| scikit-learn 1.9 instead of 1.1.2 | Compare distributions rather than single numbers; seeds are fixed |
| The public matrix differs from the private one | Documented; every positive is defined as a score above 0 |
| Tools do not install on an arm64 Mac | A first spike installed and timed all nine tools; osx-64 builds run under Rosetta |
| UMAP settings and the tree-to-distance step are not specified | Calibrated against the shipped coordinates; neighbourhoods are compared, not coordinates |
| Manual steps (LPS clustering, RBP curation, plaque annotation) | Replicate the automatic part, report the agreement, mark the manual part |
