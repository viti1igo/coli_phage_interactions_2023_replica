# 7. Configuration: what the paper specifies and what was run

`config/paper_config.yaml` records every tool version and parameter the paper states, each with its line in the published text. Settings the paper leaves out are marked `unspecified`. `config/config_diff.csv` lists every difference between that file and what I ran (37 rows, reproduced below). `config/references.yaml` holds the two reference sequences I downloaded separately (UniProt P27243 *waaL*, P76115 *yncD*). `config/envs/*.yml` are the exported specifications of the eight tool environments.

## Differences

| Component | Parameter | Paper | Used | Reason | Expected effect |
|---|---|---|---|---|---|
| prediction | scikit-learn | 1.1.2 | 1.9.1 in the main runs (T2, T3, T4); 1.1.2 on Python 3.10 in the pinned reruns (T2.1, T3.1, Phase 20) | 1.1.2 is not installable on Python 3.14; it installs on Python 3.10 only (conda-forge osx-arm64 has no 3.12 build) | pinned rerun: mean per-phage AUROC 0.7813 against 0.7812 (the authors' folds); the version alone changes the fold partition (adjusted Rand 0.047); penalty="l1" restored in the pinned run (l1_ratio=1 would fit L2 in 1.1.2) |
| prediction | random seeds | unspecified | RF random_state=0; saga solver random_state=0 | unseeded saga gave tiny run-to-run changes | reruns identical |
| prediction | interaction label | N/P/F/U matrix used by the original script | lytic = MLC > 0 on the public 0-4 matrix | private matrix not public | 20.7% positives (README says 19%) |
| prediction | CV groups | authors' file (1e-4 substitutions/site) | rebuilt: single linkage on the shipped tree at 1e-4 | file no longer public | 0.3% of our groups split across the authors' folds |
| prediction | recommender (Fig 6) | four-step algorithm | not re-run; statistics re-run on shipped cocktail results | code not in the repo | Fig 6a-b only as schematic |
| features | ST typing tool | SRST2 (version unspecified) | mlst 2.33.1, scheme ecoli_achtman_4 | SRST2 needs reads; only assemblies are public | ST agreement 95.3% |
| features | phylogroup tool | ClermonTyping (version unspecified) | ClermonTyping 24.02 --minimal | version not stated | agreement 99.0% |
| features | ECTyper | 1.0.0 | 1.0.0 | - | raw output 100% identical |
| features | Kaptive | 2.0.6, default parameters | 2.0.6 run as script from tag v2.0.6 (python 3.8, biopython<1.78, osx-64 under Rosetta); ROAR131 header sanitised | setup.py incompatible; one assembly malformed | 97.3% agreement |
| features | DefenseFinder | 1.2.0, models 1.2.3 | 1.2.0, models 1.2.3, CasFinder 3.1.0, MacSyFinder 2.1.4 | CasFinder and MacSyFinder versions not stated | F1 0.976; RM_Type_IV not detected |
| features | ABC capsule | CapsuleFinder v1 default ABC model + custom models, inter-gene max space 20 | authors' custom ABC models only, MacSyFinder 2.1.6 | default CapsuleFinder model not run | ABC presence 98.3%, K-types 89.5% |
| features | LPS outer-core type | blastp of waaL + phylogenetic clustering, manual | pangenome marker families, 5-fold CV against shipped labels | no per-type waaL reference shipped; K-12 waaL finds only K12/R1/R4 | 97.0% accuracy (supervised, not independent) |
| features | outer-membrane protein variants | clusters >99% identity, method unspecified | BLASTP (Prokka UniProt entries) + mmseqs at 99% identity; YncD from UniProt P76115 | method not stated | ARI 0.06-0.96 per protein |
| features | PanACoTA | 1.3.1 (annotate, pangenome, corepers, align, tree --fasttree) | 1.3.1, Prodigal 2.6.3, FastTree 2.2.0, QC --l90 1000 --nbcont 9999; ESCO.All.prt excluded | default QC dropped 10 genomes present in the paper's set | RF distance 0.023; patristic r 0.99999988 |
| features | tree to 8-D UMAP | umap_dims 8; other parameters unspecified; distance method unspecified | patristic distance, n_neighbors 15, min_dist 0.1, random_state 0, umap-learn 0.5.12 | parameters not stated; calibrated against shipped coordinates | 10-NN overlap 0.748 (UMAP run-to-run 0.77-0.80) |
| features | phage taxonomy | VIRIDIC 1r3.6 (genus 70%, species 95%) | BLASTn-based similarity approximation, single linkage | VIRIDIC not installed | genus ARI 0.929, species ARI 0.241 |
| features | phage network / RBP | vConTACT2 0.11.3 + Cytoscape 3.9.1; PhageRBPdetect 2.1.3 + HMMER 3.3.2 + ESMFold | not run | 1,026 GenBank genomes, annotation and structure tools not available | Fig 1k, 3b, ED3 not replicated; ED4A uses shipped RBP list only |
| modelling | per-phage MCMCglmm | 2.35; 500,000 it, burn-in 50,000, thin 250; prior V=1 nu=0.002 | MCMCglmm 2.36, R 4.3.0; same settings; N(0,10^2) fixed-effect prior for 5 phages that failed as singular | 5 phages gave singular equations | 91 of 96 phages with default prior, 96 of 96 with the extra prior |
| modelling | Fig 4c-d LMM | 1,000,000 it, burn-in 50,000, thin 250 | same, inverseA from the shipped ultrametric tree | - | beta -0.231 vs -0.22; -0.033 vs -0.04 |
| modelling | PERMANOVA | vegan 2.6-4 adonis2, adjusted R2 | own Python implementation on the Jaccard matrix, 999 permutations, adjusted R2 = 1-(1-R2)(n-1)/(n-p-1) | R vegan not used for this step | isolation strain and combined model close; genus 0.27 vs 0.13 |
| modelling | ANOVA / Tukey | scipy 1.11.2, statsmodels 0.14.1 | scipy 1.18.1 f_oneway; Tukey not run | - | F 113.9 vs 114.7 |
| modelling | nestedness / modularity (ED2) | BINMATNEST temperature T; BRIM Q; 1,000 null matrices | NODF instead of temperature; BRIM implemented here; 1,000 null matrices (Bernoulli, same connectance) | no BINMATNEST implementation available | nestedness number not comparable |
| data | strains in the interaction matrix | 403 (38,688 interactions) | 402 rows in the public matrix (38,592 cells) | one strain absent from the public file | small shifts (e.g. ANOVA F) |
| data | interaction matrix used by the figure scripts | 403 strains (Supplementary Table 3) | the 403-strain matrix when the supplement is downloaded, otherwise the public 402-strain matrix | the public file lacks strain LF110 | Fig 2f ANOVA F 114.7 with 403 strains, 113.9 with 402; model runs use the public matrix |
| modelling | significance threshold of the per-phage models | pMCMC < 1e-3 (Methods) | reported at 1e-3 and at 0.01 | the paper's Source Data for Fig 4a-b counts terms up to pMCMC 0.009, i.e. a threshold of 0.01 | 28 adsorption-factor terms at 1e-3, 51 at 0.01 (paper: 52) |
| modelling | recommender (Fig 6) | four-step algorithm, code not published | re-implemented from the Methods (scripts/recommender.py); all features per split; step 2 needs a complete ST/O/H haplotype | test-strain genomes are not shipped: ABC and Klebsiella types unknown, UMAP recomputed from the shipped distances | 2 of 100 cocktails identical, 0.97 of 3 phages in common |
| data | AMR gene counts (ED6 B, D) | annotation tool not stated | the paper's counts from its Source Data | no AMR annotation is shipped | ED6 B and D drawn from the paper's counts |
| figures | phage annotation (Fig 3b, ED3) | Prokka 1.14.5 --kingdom Viruses; RBP genes by PhageRBPdetect + curation | pyrodigal (single mode, metagenomic mode when an RBP identifier exceeds the gene count); RBP genes located by the number in their identifier in RBP_list.csv | the authors' annotation is not shipped | gene boundaries and numbering may differ for a few genes (127 of 130 RBP identifiers mapped) |
| figures | genome maps and alignment (Fig 3b, ED3) | Clinker for the maps; a multiple whole-genome alignment (tool not stated) | own matplotlib maps (reciprocal best BLOSUM62 identity >= 30 %); MAFFT 7.526 --auto; entropy over nucleotides, 200-nt mean | Clinker and the alignment tool are not available or not named | layout differs; entropy peaks at the RBP region as in the paper |
| notebooks | the authors' four notebooks | - | run on a scratch copy, as shipped and with path/compatibility patches (scripts/notebook_patches.py) | Windows and private paths, a moved file, a missing package, newer pandas and matplotlib | no analysis changed; one cell of the authors' still fails |
| modelling | PERMANOVA covariates (SI Text 8) | adsorption factors and defense subfamilies as presence/absence, traits with fewer than 8 observations lumped, 8 PCoA phylogeny axes first | indicators of the levels with at least 8 strains (O, H, LPS, Klebsiella, ABC, 12 outer-membrane variants; defense subfamilies at 80/80), Jaccard distance of the lytic profiles, AIC = n ln(SS_res/n) + 2 rank | covariate lists and the distance are not given | AF and full adjusted R2 agree (0.26, 0.27 vs 0.27, 0.28); defense 0.17 vs 0.05 |
| data | genome set for tree | 403 strains + references | 402 strains in common with the shipped tree | ESCO.All.prt is a proteome file, not a genome | - |
| prediction | Python | not stated (only the plotting notebook's kernel metadata says 3.12.2; the prediction script records none) | 3.14.7 (main runs); 3.10.22 (pinned reruns) | authors' Python for the model step is unknown; scikit-learn 1.1.2 has no Python 3.12 build |  so Python follows scikit-learn |
| prediction | scipy | 1.11.2 (Methods, scipy.stats) | 1.18.1 (main runs); 1.11.2 (pinned reruns) | - | none expected |
| prediction | numpy / pandas / joblib | not stated | 2.5.4 / 3.0.6 / 1.6.0 (main); 1.23.5 / 1.5.3 / 1.2.0 (pinned) | not stated; pinned versions are contemporary with scikit-learn 1.1.2 and compatible with it | unknown, small |
| prediction | model randomness | np.random.seed(0) once; RandomForest, Dummy and saga have no random_state; phages run in sequence | random_state=0 per RF and Dummy; saga unseeded; phages in parallel | the global random stream cannot be replayed: it depends on private labels and on row order | RF probabilities differ (median abs 0.008-0.018 with their folds); mean AUROC noise 0.0011 |
| prediction | operating system | Windows (script paths D:\...) | macOS arm64 (Apple M2) | no Windows machine | floating-point differences about 1e-4 on a few LogReg L2 phages |

## Settings the paper does not state

The paper leaves these out: UMAP hyperparameters, how the tree becomes a distance matrix, the SRST2 and ClermonTyping versions, PanACoTA options other than `--fasttree`, the method for outer-membrane-protein clustering, the LPS clustering parameters and the random seeds. I calibrated or fixed each of them, and the table above marks them.

## Software

| Purpose | Environment | Versions actually run |
|---|---|---|
| Tree | `pan` | PanACoTA 1.3.1, Prodigal 2.6.3, FastTree 2.2.0, BLAST 2.16.0, Python 3.6 (needs python-utils < 3) |
| ST | `mlst` | mlst 2.33.1, BLAST 2.17.0 |
| Phylogroup | `clermont` | ClermonTyping 24.02 |
| Serotype | `ectyper` | ECTyper 1.0.0 |
| Klebsiella capsule | `kaptive` | Kaptive 2.0.6 run as a script from tag v2.0.6, Python 3.8, biopython < 1.78, osx-64 under Rosetta |
| Defense systems | `defense` | DefenseFinder 1.2.0, MacSyFinder 2.1.4, models 1.2.3, CasFinder 3.1.0 |
| ABC capsule | `macsy` | MacSyFinder 2.1.6 with the authors' custom models |
| Protein and nucleotide search | `blast` | BLAST 2.17.0, mmseqs |
| Alignment (Fig 3b) | `aln` | MAFFT 7.526 |
| Models and figures | `.venv` | Python 3.14, scikit-learn 1.9.1, scipy 1.18.1, pandas 3.0.6, numpy 2.5.4, umap-learn 0.5.12, matplotlib 3.11.2, biopython 1.88, pyrodigal, statsmodels 0.15, bokeh 3.10, nbclient |
| Mixed models | R 4.3.0 | MCMCglmm 2.36, ape 5.8.1 |
| Pinned model rerun (T2.1, T3.1) | `ml_pinned` ([`environment.yml`](../environment.yml)) | Python 3.10.22, scikit-learn 1.1.2, scipy 1.11.2, numpy 1.23.5, pandas 1.5.3, joblib 1.2.0, OpenBLAS 0.3.25 |

`scripts/setup_envs.sh <env>` creates the environments with micromamba, and `scripts/setup_pinned_env.sh` creates `ml_pinned`. [`results/audit/env_match.csv`](../results/audit/env_match.csv) lists 13 components of the authors' environment with what was matched and the reason for each mismatch. The paper states scikit-learn 1.1.2 and scipy but no Python version for the model step (3.12.2 appears only in the kernel metadata of the authors' plotting notebook), so Python 3.10 follows scikit-learn 1.1.2, which has no build for 3.12. On an Apple-silicon Mac the osx-64 environments need Rosetta. The `conda` installation on my machine was broken, which is why I used micromamba.

## Libraries per tier against the authors'

The tiers are defined in [report 12](12_gap_audit.md): T0 is the paper, T1 the authors' shipped predictions, T2 our model run on the shipped inputs, T3 on the regenerated features, T4 on the regenerated features and matrix. T2.1 and T3.1 repeat T2 and T3 in the pinned environment. "Not stated" means neither the paper, the authors' README nor their files give the version. Sources: the paper (Methods, `docs/paper/published.txt` line 1562 for scikit-learn and lines 1649 to 1651 for scipy), the authors' README (section 1.6), their scripts and notebooks, and the environments in `config/envs/`. [`results/audit/env_match.csv`](../results/audit/env_match.csv) holds the same comparison with the reason for each mismatch.

### Model step (`scripts/predict.py`)

| Tier | Python | scikit-learn | numpy | pandas | scipy | joblib | Linear algebra | System |
|---|---|---|---|---|---|---|---|---|
| T0 paper | not stated | 1.1.2 | not stated | not stated | 1.11.2 (scipy.stats only) | not stated | not stated | not stated |
| T1 authors' predictions | not stated (3.12.2 is only the kernel metadata of their plotting notebook) | 1.1.2 (paper and README) | not stated | not stated | not stated | not stated | not stated | Windows (their scripts read `D:\These\...`) |
| T2 shipped inputs | 3.14.7 | 1.9.1 | 2.5.4 | 3.0.6 | 1.18.1 | 1.6.0 | Apple Accelerate | macOS arm64, Apple M2, 8 cores |
| T3 regenerated features | 3.14.7 | 1.9.1 | 2.5.4 | 3.0.6 | 1.18.1 | 1.6.0 | Apple Accelerate | same machine |
| T4 regenerated features and matrix | 3.14.7 | 1.9.1 | 2.5.4 | 3.0.6 | 1.18.1 | 1.6.0 | Apple Accelerate | same machine |
| T2.1 shipped inputs, pinned | 3.10.22 | 1.1.2 | 1.23.5 | 1.5.3 | 1.11.2 | 1.2.0 | OpenBLAS 0.3.25 (conda-forge) | same machine |
| T3.1 regenerated features, pinned | 3.10.22 | 1.1.2 | 1.23.5 | 1.5.3 | 1.11.2 | 1.2.0 | OpenBLAS 0.3.25 (conda-forge) | same machine |

- T2, T3 and T4 use the same libraries and differ in their inputs. T2.1 and T3.1 use the same libraries and differ in their inputs in the same way. T2.1 also has a variant with the authors' fold assignment.
- scikit-learn and scipy match the paper only in T2.1 and T3.1. numpy, pandas and joblib are not stated anywhere, so the pinned versions (1.23.5, 1.5.3, 1.2.0) are contemporary with scikit-learn 1.1.2, not matched. threadpoolctl is 3.7.0 in `.venv` and 3.6.0 in `ml_pinned`.
- Python 3.10 is the only Python with a scikit-learn 1.1.2 build on conda-forge for macOS arm64. The authors' Python for the model step is unknown, so it is a choice, not a mismatch.
- The authors' script holds the same estimators in every tier. Only `penalty="l1"` is handled by version: scikit-learn 1.1.2 needs it and 1.8 and later drop it for `l1_ratio=1`.
- T3 and T4 build their features with the tools below and then run the model step in `.venv`. T3.1 reuses the T3 feature files and only reruns the model step in `ml_pinned`, so the tools below did not change between T3 and T3.1. T4 adds `scripts/build_matrix.py`, which runs in `.venv` with numpy 2.5.4 and pandas 3.0.6.

### Feature and matrix building (T3, T3.1, T4)

| Step | Paper | Used | Environment | Same |
|---|---|---|---|---|
| Core-genome tree | PanACoTA 1.3.1, `--fasttree` | PanACoTA 1.3.1, Prodigal 2.6.3, FastTree 2.2.0, MAFFT 7.526, Mash 2.3, BLAST 2.16.0, Python 3.6 | `pan` | tool yes; the rest not stated |
| Tree to 8-D UMAP | umap, version and parameters not stated | umap-learn 0.5.12 with numba 0.68.0, pynndescent 0.6.0, scikit-learn 1.9.1 | `.venv` | not stated |
| Sequence type | SRST2, version not stated | mlst 2.33.1, scheme `ecoli_achtman_4`, BLAST 2.17.0 | `mlst` | no (SRST2 needs reads) |
| Phylogroup | ClermonTyping, version not stated | ClermonTyping 24.02, Mash 2.3, BLAST 2.17.0 | `clermont` | tool yes; version not stated |
| O and H type | ECTyper 1.0.0 | ECTyper 1.0.0 | `ectyper` | yes |
| Klebsiella capsule | Kaptive 2.0.6 | Kaptive 2.0.6 (tag v2.0.6), Python 3.8, BLAST 2.17.0 | `kaptive` | yes |
| Defense systems | DefenseFinder 1.2.0, models 1.2.3 | DefenseFinder 1.2.0, models 1.2.3, MacSyFinder 2.1.4, HMMER 3.4 | `defense` | yes |
| ABC capsule | CapsuleFinder v1 on MacSyFinder, mmseqs2 13.45111 for the custom models | MacSyFinder 2.1.6, MacSylib 1.0.4, HMMER 3.4, mmseqs2 13.45111 | `macsy`, `pan` | mmseqs2 yes; MacSyFinder version not stated |
| LPS outer-core type | blastp of waaL and manual clustering | BLAST 2.17.0 and cross-validated pangenome markers | `blast`, `.venv` | no (manual step) |
| Outer-membrane proteins | clusters above 99 % identity, tool not stated | BLAST 2.17.0 and mmseqs2 | `blast` | not stated |
| Interaction matrix (T4) | median of replicates, rule partly stated | `scripts/build_matrix.py` | `.venv` | not applicable |

### Statistics and figures (not part of the tiers)

| Component | Paper | Used |
|---|---|---|
| scipy.stats | 1.11.2 | 1.18.1 (`.venv`) |
| statsmodels | 0.14.1 | 0.15.0 (`.venv`) |
| MCMCglmm | 2.35 | 2.36 (R 4.3.0, local library `.Rlib`) |
| ape | not stated | 5.8.1 |
| vegan | 2.6-4 for PERMANOVA | not installed; PERMANOVA is implemented in Python (`scripts/matrix_stats.py`) |
| R | not stated | 4.3.0 |
| matplotlib, seaborn, bokeh | not stated (the authors' notebooks import them) | 3.11.2, 0.13.2, 3.10.1 |

Nothing in this table was rerun in the pinned environment, so the figures and mixed-model tables still use the newer packages.

### What the numbers say about the environment

- Only scikit-learn and scipy.stats are stated, and only T2.1 and T3.1 match them. The other libraries, the operating system, the linear-algebra library and the random stream are unknown or cannot be matched, so the effect of those parts is not isolated.
- The authors' README says the cross-validation group file `370+host_cross_validation_groups_1e-4.csv` is no longer available. The name most likely refers to an earlier working set of 370 collection strains plus the phage-host strains; the authors do not explain it. The files in the repository that carry this name hold 402 or 404 strains, and the paper counts 369 collection strains plus 34 host strains (403). The model runs here use the 402 strains of the public matrix.

## Platform notes

I ran everything on macOS arm64 with 16 GB of memory and kept all data on an exFAT external drive. exFAT creates a `._*` twin for every file. The twins broke DefenseFinder and MacSyFinder, which scan their model folders, and a Python package installed on the drive (matplotlib). The scripts therefore run DefenseFinder in a temporary directory on the internal disk and give MacSyFinder a clean models folder, and `scripts/cleanup_intermediates.sh` removes the twins.
