# Replicating a phage-host prediction pipeline from the raw genomes

This repository replicates, starting from the raw data, the pipeline of Gaborieau et al., *Prediction of strain-level phage-host interactions across the *Escherichia* genus using only genomic information*, Nature Microbiology 9, 2847-2861 (2024), [doi:10.1038/s41564-024-01832-5](https://doi.org/10.1038/s41564-024-01832-5). The authors' own repository is https://github.com/mdmparis/coli_phage_interactions_2023.

## What this repository does

The authors share their final feature tables and models. For most preprocessing steps they share the outputs but not the code. Here, every model input is rebuilt from the 403 bacterial assemblies, the phage genomes and the raw plaque scores, using the tools and versions named in the paper. For each step the repository records how closely the result matches the authors' table, which settings the paper leaves out, and where the settings used differ from the paper's.

## What was rebuilt and how closely it matches

| Step | Tool (version) | Agreement with the shipped table |
|---|---|---|
| Core-genome tree | PanACoTA 1.3.1, Prodigal, FastTree | 390 of 399 splits shared (normalised Robinson-Foulds 0.023); patristic-distance Pearson 0.99999988 |
| Cross-validation groups | single linkage at 1e-4 substitutions per site | adjusted Rand index 0.984 |
| 8-D UMAP | umap-learn 0.5.12 | 10-nearest-neighbour overlap 0.748 (UMAP's own run-to-run overlap is 0.77 to 0.80) |
| Sequence type | mlst 2.33.1 (SRST2 needs reads, and only assemblies are public) | 95.3 % |
| Phylogroup | ClermonTyping 24.02 | 99.0 % |
| O and H type | ECTyper 1.0.0 | raw output identical for 402 of 402 strains; 74.7 % and 91.8 % against the hand-curated columns |
| LPS outer-core type | pangenome marker families, cross-validated | 97.0 % (trained on the shipped labels) |
| ABC capsule | MacSyFinder with the authors' custom models | presence 98.3 %, K-type 89.5 % |
| Klebsiella capsule | Kaptive 2.0.6 | 97.3 % |
| Defense systems | DefenseFinder 1.2.0, models 1.2.3, CasFinder 3.1.0 | presence F1 0.976; RM_Type_IV not called |
| Outer-membrane protein variants | BLASTP and mmseqs at 99 % identity | adjusted Rand index 0.06 to 0.96, depending on the protein |
| Interaction matrix | lowest concentration with lysis, median of replicates | lytic or not 99.8 %; exact score 96.6 % |
| Phage genus and species | BLASTn similarity in place of VIRIDIC | genus ARI 0.929; species ARI 0.241 |

[reports/03_feature_regeneration.md](reports/03_feature_regeneration.md) gives the details, the caveats, and the manual steps that cannot be reproduced exactly: LPS clustering, receptor-binding-protein curation and plaque annotation.

## Comparison with the paper

- The publisher's Supplementary Information, Supplementary Tables 1 to 7 and Source Data workbooks (13 files, 8.7 MB) are fetched by `scripts/download_supplement.sh` and compared with this replication by `scripts/compare_source_data.py`. They are not redistributed here. They settled several open differences: strain LF110 is missing from the public matrix, the paper's new phage names map exactly to the public ones, one headline AUROC is an average per strain, the Source Data counts significant terms at a threshold of 0.01, and the AMR counts needed for two extended-data panels are in the Source Data.
- [`results/comparison/paper_vs_replication.csv`](results/comparison/paper_vs_replication.csv) compares 97 numbers one by one, each with its place in the article. 33 match within 5 %, 12 match, 1 matches once explained, 4 agree with the paper's confidence interval or item list, 23 differ, and 24 are qualitative. [reports/05_paper_comparison.md](reports/05_paper_comparison.md) discusses the ones that differ.
- I ran the authors' repository itself: its four notebooks (as shipped, where 27 of 61 code cells fail, and with path and compatibility patches), the supplementary analyses whose data it holds, and the genome maps. [reports/11_repository_vs_paper.md](reports/11_repository_vs_paper.md) lists every difference found between the repository and the paper. [reports/10_repository_analyses.md](reports/10_repository_analyses.md) describes the runs.
- [`config/config_diff.csv`](config/config_diff.csv) lists 32 differences between the paper's settings and the settings used, each with a reason and an expected effect. [`config/paper_config.yaml`](config/paper_config.yaml) holds the paper's own settings with line numbers and marks the ones it leaves unspecified.
- 18 figures and extended-data figures are redrawn from the data, plus three schematics redrawn from the legends. Each comes as PNG, SVG and CSV, and [reports/06_figures.md](reports/06_figures.md) annotates every one. Figures that need data nobody has shared (the 1,026-genome phage network, the test-strain tree) are listed with the reason.
- The paper's cocktail recommender has no published code. `scripts/recommender.py` re-implements it from the README and the Methods, and it does not reproduce the authors' cocktails (see reports/05).
- Four things would have to come from the authors: the cross-validation group file, the recommender code, the test-strain genomes and the genomes of the phage network. [reports/05_paper_comparison.md](reports/05_paper_comparison.md) says what each would settle.

## Reports

| File | Content |
|---|---|
| [01_overview_and_plan.md](reports/01_overview_and_plan.md) | Goal, phases, what was known beforehand, decisions, risks |
| [02_model_reproduction.md](reports/02_model_reproduction.md) | The models run on the shipped inputs and on the regenerated inputs |
| [03_feature_regeneration.md](reports/03_feature_regeneration.md) | Each input rebuilt from the genomes and compared with the shipped table |
| [04_feature_necessity.md](reports/04_feature_necessity.md) | Which feature groups the model needs, on both input sets |
| [05_paper_comparison.md](reports/05_paper_comparison.md) | Paper, shipped results and replication side by side; mismatches; missing files |
| [06_figures.md](reports/06_figures.md) | Which figures were drawn, with an annotation for each |
| [07_configuration.md](reports/07_configuration.md) | The paper's settings, what was run, software versions |
| [08_code_audit.md](reports/08_code_audit.md) | The code audit and the simplifications applied |
| [09_task_log.md](reports/09_task_log.md) | One line per task with its outcome |
| [10_repository_analyses.md](reports/10_repository_analyses.md) | The authors' notebooks, the supplementary analyses and the genome maps |
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

The authors' data, the genome assemblies and the paper's PDFs are not redistributed. The authors' repository has no licence file, so check its terms before reusing their material.

## Reproduce from scratch

You need Python 3.14 with `requirements.txt`, R 4.3 with MCMCglmm, and micromamba for the tool environments. Plan for about 10 GB of disk for the intermediate files (the PanACoTA run alone takes 4 GB), 2 GB for the assemblies and several GB for the tool environments. Run everything from the repository root.

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
MATRIX_SET=regenerated Rscript scripts/fig4_cd.R      # the same models on the regenerated matrix
MATRIX_SET=regenerated Rscript scripts/fig_ed6.R

# 5. The authors' notebooks, supplementary analyses and genome maps
.venv/bin/python scripts/run_notebooks.py shipped
.venv/bin/python scripts/run_notebooks.py patched
.venv/bin/python scripts/si_analyses.py
MATRIX_SET=regenerated .venv/bin/python scripts/si_analyses.py
.venv/bin/python scripts/regen_matrix_analyses.py
Rscript scripts/si_rbp_lmm.R
MATRIX_SET=regenerated Rscript scripts/si_rbp_lmm.R
MAFFT=/path/to/mafft .venv/bin/python scripts/fig_phage_maps.py

# 6. Figures and comparison (the supplement comparison and the recommender come first)
.venv/bin/python scripts/compare_source_data.py
.venv/bin/python scripts/recommender.py
for s in fig_main_1_2 fig1a_tree fig_main_3_6 fig_extended fig_schematics fig5_curves; do .venv/bin/python scripts/$s.py; done
.venv/bin/python scripts/build_comparison.py
.venv/bin/python scripts/fig_comparison.py
```

## Caveats

- Each model run uses one seed, and UMAP layouts differ from run to run.
- The public interaction matrix has 402 strains where the paper has 403 (strain LF110). The figure scripts use Supplementary Table 3 when you have downloaded it.
- The recommender of the paper's cocktail experiment could not be re-run because its code is not in the authors' repository.
- SRST2, VIRIDIC and the receptor-binding-protein tools were replaced or skipped. `config/config_diff.csv` records each case.
- The definition behind one extended-data curve (coverage by the most-covering phage sets) could not be derived from the matrix. See reports/05.
