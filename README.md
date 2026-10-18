# Replicating a phage–host prediction pipeline from the raw genomes

A replication, started from the raw data, of the pipeline of Gaborieau et al., *Prediction of strain-level phage–host interactions across the *Escherichia* genus using only genomic information*, Nature Microbiology 9, 2847–2861 (2024), [doi:10.1038/s41564-024-01832-5](https://doi.org/10.1038/s41564-024-01832-5). The authors' repository is https://github.com/mdmparis/coli_phage_interactions_2023.

## The task

The authors ship the final feature tables and models but, for most preprocessing steps, only the outputs and not the code. This repository rebuilds every model input from the 403 bacterial assemblies, the phage genomes and the raw plaque scores, with the tools and versions stated in the paper, and records for each step how closely the result agrees with the authors' table, which settings the paper does not state, and every difference between the paper's configuration and what was run.

## What was rebuilt, and how closely it agrees

| Step | Tool (version) | Agreement with the shipped table |
|---|---|---|
| Core-genome tree | PanACoTA 1.3.1, Prodigal, FastTree | 390 of 399 splits shared (normalised Robinson–Foulds 0.023); patristic-distance Pearson 0.99999988 |
| Cross-validation groups | single linkage at 1e-4 substitutions per site | adjusted Rand index 0.984 |
| 8-D UMAP | umap-learn 0.5.12 | 10-nearest-neighbour overlap 0.748 (UMAP's own run-to-run overlap: 0.77–0.80) |
| Sequence type | mlst 2.33.1 (SRST2 needs reads, only assemblies are public) | 95.3 % |
| Phylogroup | ClermonTyping 24.02 | 99.0 % |
| O and H type | ECTyper 1.0.0 | raw output identical for 402 of 402 strains; 74.7 % and 91.8 % against the hand-curated columns |
| LPS outer-core type | pangenome marker families, cross-validated | 97.0 % (supervised by the shipped labels) |
| ABC capsule | MacSyFinder with the authors' custom models | presence 98.3 %, K-type 89.5 % |
| Klebsiella capsule | Kaptive 2.0.6 | 97.3 % |
| Defense systems | DefenseFinder 1.2.0, models 1.2.3, CasFinder 3.1.0 | presence F1 0.976; RM_Type_IV not called |
| Outer-membrane protein variants | BLASTP + mmseqs at 99 % identity | adjusted Rand index 0.06 to 0.96 by protein |
| Interaction matrix | lowest concentration with lysis, median of replicates | lytic vs not 99.8 %; exact score 96.6 % |
| Phage genus / species | BLASTn similarity in place of VIRIDIC | genus ARI 0.929; species ARI 0.241 |

Details, caveats and the manual steps that cannot be reproduced exactly (LPS clustering, receptor-binding-protein curation, plaque annotation) are in [reports/03_feature_regeneration.md](reports/03_feature_regeneration.md).

## Comparison with the paper

- The publisher's Supplementary Information, Supplementary Tables 1–7 and Source Data workbooks (13 files, 8.7 MB, fetched by `scripts/download_supplement.sh`, not redistributed here) are compared with this replication by `scripts/compare_source_data.py`. They settled several open differences: the 403rd strain (LF110) missing from the public matrix, the exact mapping of the paper's new phage names, the aggregation behind one headline AUROC, the significance threshold actually used in the Source Data, and the AMR counts needed for two extended-data panels.
- 97 numbers are compared one by one in [`results/comparison/paper_vs_replication.csv`](results/comparison/paper_vs_replication.csv), each with its place in the article: 33 match within 5 %, 12 match, 1 matches once explained, 4 are consistent with the paper's confidence interval or item list, 23 differ, and 24 are qualitative. The ones that differ, and what is known about each, are in [reports/05_paper_comparison.md](reports/05_paper_comparison.md).
- The authors' repository itself was run: its four notebooks (as shipped, which fail for 15 + 10 + 2 cells, and with path and compatibility patches), the supplementary analyses whose data it contains, and the genome maps. Every difference found between the repository and the paper is listed in [reports/11_repository_vs_paper.md](reports/11_repository_vs_paper.md); the runs are in [reports/10_repository_analyses.md](reports/10_repository_analyses.md).
- 28 differences between the paper's settings and the settings used are listed with reason and expected effect in [`config/config_diff.csv`](config/config_diff.csv); the paper's own settings, with line numbers and the ones it leaves unspecified, are in [`config/paper_config.yaml`](config/paper_config.yaml).
- 18 figures and extended-data figures of the paper are redrawn from the data (and three schematics redrawn from the legends), each as PNG, SVG and CSV: [reports/06_figures.md](reports/06_figures.md). Figures that need data not available here (the 1,026-genome phage network, the test-strain tree) are listed with the reason.
- The paper's cocktail recommender has no published code. A re-implementation from the Methods is in `scripts/recommender.py`; it does not reproduce the authors' cocktails (reports/05).
- What is still missing and would have to come from the authors (cross-validation group file, recommender code, test-strain genomes, the phage-network genomes): [reports/05_paper_comparison.md](reports/05_paper_comparison.md).

## Reports

| File | Content |
|---|---|
| [01_overview_and_plan.md](reports/01_overview_and_plan.md) | Goal, phases, facts established, decisions, risks |
| [02_model_reproduction.md](reports/02_model_reproduction.md) | The models run on the shipped inputs and on the regenerated inputs |
| [03_feature_regeneration.md](reports/03_feature_regeneration.md) | Each input rebuilt from the genomes and compared with the shipped table |
| [04_feature_necessity.md](reports/04_feature_necessity.md) | Which feature groups the model needs, on both input sets |
| [05_paper_comparison.md](reports/05_paper_comparison.md) | Paper vs shipped results vs replication, mismatches, files still needed |
| [06_figures.md](reports/06_figures.md) | Coverage of the paper's figures and every figure drawn |
| [07_configuration.md](reports/07_configuration.md) | Paper settings, what was run, software versions |
| [08_code_audit.md](reports/08_code_audit.md) | Code audit and the simplifications applied |
| [09_task_log.md](reports/09_task_log.md) | Task list with the outcome of each step |
| [10_repository_analyses.md](reports/10_repository_analyses.md) | The authors' notebooks, the supplementary analyses and the genome maps, run from their repository |
| [11_repository_vs_paper.md](reports/11_repository_vs_paper.md) | Every difference found between the authors' repository and the paper |

## Layout

```
config/       paper_config.yaml, config_diff.csv, references.yaml, envs/*.yml (tool environments)
scripts/      all code (Python, R, shell)
replication/  regenerated tables: features, interaction matrix, tree, UMAP, typing outputs
notebooks/    the authors' notebooks executed (as shipped and patched), their errors and outputs
results/      replication/ (agreement tables), comparison/, figures/, model outputs
data/derived/ cross-validation groups, genome sizes
reports/      the eleven reports above
```

The authors' data, the genome assemblies and the paper's PDFs are not redistributed. The authors' repository has no licence file; check its terms before reuse.

## Reproduce from scratch

Python 3.14 with `requirements.txt`; R 4.3 with MCMCglmm; micromamba for the tool environments. About 10 GB of disk for the intermediate files (the PanACoTA run alone is about 4 GB), 2 GB for the assemblies and several GB for the tool environments. Run from the repository root.

```bash
# 0. Environment
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
scripts/setup_envs.sh                      # tool environments (micromamba; Rosetta on Apple silicon)

# 1. Data: the authors' repository, the 403 assemblies (md5-checked) and the publisher's supplement
scripts/download_data.sh --genomes
scripts/download_supplement.sh

# 2. Model inputs rebuilt from the genomes (hours; each step is resumable)
scripts/run_panacota.sh                    # core-genome tree
scripts/run_typing.sh                      # ST, phylogroup, Klebsiella capsule
scripts/run_proteome_tools.sh all          # defense systems, ABC capsule
.venv/bin/python scripts/run_omp.py        # outer-membrane protein variants
.venv/bin/python scripts/run_lps_markers.py
.venv/bin/python scripts/compare_tree.py   # tree agreement and our UMAP
.venv/bin/python scripts/compare_defense.py
.venv/bin/python scripts/build_host_table.py
.venv/bin/python scripts/assemble_features.py
.venv/bin/python scripts/compare_features.py
.venv/bin/python scripts/build_matrix.py
.venv/bin/python scripts/run_phage_taxonomy.py

# 3. Models on the shipped and the regenerated inputs (fixed seeds, about 1 hour)
.venv/bin/python scripts/make_groups.py
.venv/bin/python scripts/make_groups.py replication/T14_panacota/tree/ESCO.nucl.grp.aln.fasttree_tree.nwk data/derived/cv_groups_1e-4_regenerated.csv
scripts/run_seeded_models.sh

# 4. Mixed models (R)
Rscript scripts/fig4_mcmc.R results/fig4 500000 250 50000 6   # about 4.5 hours on 6 cores
.venv/bin/python scripts/fig4_summarize.py
Rscript scripts/fig4_cd.R
Rscript scripts/fig_ed6.R                 # needs data/derived/paper_amr_counts.csv, written by compare_source_data.py

# 5. The authors' notebooks, supplementary analyses and genome maps
.venv/bin/python scripts/run_notebooks.py shipped
.venv/bin/python scripts/run_notebooks.py patched
.venv/bin/python scripts/si_analyses.py
.venv/bin/python scripts/regen_matrix_analyses.py
Rscript scripts/si_rbp_lmm.R
MAFFT=/path/to/mafft .venv/bin/python scripts/fig_phage_maps.py

# 6. Figures and comparison (the supplement comparison and the recommender come first)
.venv/bin/python scripts/compare_source_data.py
.venv/bin/python scripts/recommender.py
for s in fig_main_1_2 fig1a_tree fig_main_3_6 fig_extended fig_schematics fig5_curves; do .venv/bin/python scripts/$s.py; done
.venv/bin/python scripts/build_comparison.py
.venv/bin/python scripts/fig_comparison.py
```

## Caveats

- One seed per model run; UMAP layouts differ from run to run.
- The public interaction matrix has 402 strains where the paper has 403 (strain LF110); the figure scripts use Supplementary Table 3 when it has been downloaded.
- The recommender of the paper's cocktail experiment was not re-run (its code is not in the authors' repository).
- SRST2, VIRIDIC and the receptor-binding-protein tools were replaced or skipped; each case is in `config/config_diff.csv`.
- The definition behind one extended-data curve (coverage by the most-covering phage sets) could not be derived from the matrix, and the cocktail recommender could not be reproduced; see reports/05.
