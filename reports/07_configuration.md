# 7. Configuration: what the paper specifies and what was run

`config/paper_config.yaml` records every tool version and parameter stated in the paper, each with its line in the published text; settings the paper does not state are marked `unspecified`. `config/config_diff.csv` lists every difference between that file and what was run (28 rows, reproduced below). `config/references.yaml` holds the two reference sequences downloaded separately (UniProt P27243 *waaL*, P76115 *yncD*). `config/envs/*.yml` are the exported environment specifications of the eight tool environments.

## Differences

| Component | Parameter | Paper | Used | Reason | Expected effect |
|---|---|---|---|---|---|
| prediction | scikit-learn | 1.1.2 | 1.9.1 | 1.1.2 not installable on Python 3.14 | none expected; L1 penalty arg replaced by l1_ratio=1, sparse solution checked identical |
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
| data | genome set for tree | 403 strains + references | 402 strains in common with the shipped tree | ESCO.All.prt is a proteome file, not a genome | - |

## Settings the paper does not state

UMAP hyperparameters; how the tree is turned into a distance matrix; SRST2 version; ClermonTyping version; PanACoTA options other than `--fasttree`; the method for outer-membrane-protein clustering; the LPS clustering parameters; random seeds. These were calibrated or fixed here and are marked as such above.

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
| Models and figures | `.venv` | Python 3.14, scikit-learn 1.9.1, scipy 1.18.1, pandas 3.0.6, numpy 2.5.4, umap-learn 0.5.12, matplotlib 3.11.2, biopython 1.88 |
| Mixed models | R 4.3.0 | MCMCglmm 2.36, ape 5.8.1 |

Environments are created with micromamba (`scripts/setup_envs.sh <env>`). On an Apple-silicon Mac the osx-64 environments need Rosetta. The original `conda` installation on the test machine was broken, hence micromamba.

## Platform notes

The work ran on macOS arm64 with 16 GB of memory, with all data on an exFAT external drive. exFAT creates a `._*` twin for every file. These twins broke DefenseFinder and MacSyFinder (they scan their model folders) and a Python package installed on the drive (matplotlib). The scripts therefore run DefenseFinder in a temporary directory on the internal disk, give MacSyFinder a clean models folder, and `scripts/cleanup_intermediates.sh` removes the twins.
