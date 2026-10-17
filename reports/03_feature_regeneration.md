# 3. Regenerating the model inputs from the genomes

The authors ship the final feature tables but, for most preprocessing steps, only the outputs and not the code. This step rebuilds each input from the raw data with the tools and versions stated in the paper, and compares the result with the shipped table. Tool versions and every deviation are listed in [07_configuration.md](07_configuration.md).

Inputs: 403 bacterial assemblies (figshare 25941691, md5-checked), the 96 phage genomes shipped in the authors' repository, and the raw plaque-assay table `raw_interactions.csv`. The 9.9 GB of plaque images were not downloaded.

## Summary

| Step | Tool (version) | Agreement with the shipped table |
|---|---|---|
| Core-genome tree | PanACoTA 1.3.1, Prodigal, FastTree | 390 of 399 splits shared (normalised Robinson–Foulds 0.023); patristic-distance Pearson 0.99999988 |
| Cross-validation groups (< 1e-4) | single linkage on the tree | 300 vs 301 groups, adjusted Rand index 0.984 |
| 8-D UMAP | umap-learn 0.5.12 (n_neighbors 15, min_dist 0.1) | 10-nearest-neighbour overlap 0.748; UMAP's own run-to-run overlap is 0.77–0.80 |
| Sequence type | mlst 2.33.1 (scheme `ecoli_achtman_4`) instead of SRST2 | 95.3 % (384 of 403) |
| Phylogroup | ClermonTyping 24.02 | 99.0 % (399 of 403) |
| O and H type | ECTyper 1.0.0 | raw tool output identical for 402 of 402 strains; against the hand-curated columns 74.7 % (O) and 91.8 % (H) |
| LPS outer-core type | pangenome marker families, 5-fold cross-validated | 97.0 % (390 of 402), supervised by the shipped labels |
| ABC capsule | MacSyFinder with the authors' custom models | presence 98.3 %; K-type 89.5 % (85 of 95) |
| Klebsiella capsule | Kaptive 2.0.6 | 97.3 % (392 of 403); 20 of 22 shipped calls reproduced |
| Defense systems | DefenseFinder 1.2.0, models 1.2.3, CasFinder 3.1.0 | presence F1 0.976; RM_Type_IV not called by our run |
| Outer-membrane protein variants | BLASTP + mmseqs at 99 % identity | adjusted Rand index 0.06 (tolC) to 0.96 (ompC) |
| Host (isolation-strain) table | merge of the above | 33 of 33 strains agree on phylogroup, ST, LPS; O and H differ in 1 strain only in how "no call" is written |
| Interaction matrix | lowest concentration with lysis; median of replicates | lytic vs not 99.8 %; exact score 96.6 % |
| Phage genus / species | BLASTn similarity (VIRIDIC-like) | genus ARI 0.929; species ARI 0.241 |

All rows are in `results/replication/phase6_summary.csv`, `matrix_agreement.csv` and `phage_taxonomy_agreement.csv`. A bar chart is in [06_figures.md](06_figures.md) (`cmp_feature_agreement`).

## Details and caveats

**Tree.** Prodigal replaces Prokka (Prokka under Rosetta was about 30 times slower). The paper does not state the quality-control thresholds. PanACoTA's defaults dropped 10 genomes that are in the paper's set, so the run uses `--l90 1000 --nbcont 9999`. A file `ESCO.All.prt` in the genome folder is a proteome, not a genome, and is excluded. Core genes are those present in all genomes (1,295 families). Total tree length is 0.8683 in both trees.

**UMAP.** The paper states 8 dimensions but not the other parameters, nor how the tree becomes a distance matrix. Patristic distance is used. The parameters were calibrated against the shipped coordinates (`results/replication/umap_calibration.csv`). Layouts differ run to run, so neighbourhoods rather than coordinates are compared.

**Sequence type.** SRST2 needs sequencing reads, and only assemblies are public. Most mismatches are novel STs that the current database numbers differently.

**Serotype.** The authors curated multi-calls and "-" calls by hand, which explains the gap between raw and curated agreement.

**LPS.** The paper uses blastp of *waaL* followed by manual phylogenetic clustering. K-12 WaaL (UniProt P27243) recognises only type K12 (97.6 % identity) and weakly R1 and R4; R2 and R3 waaL are not homologous to it, and no per-type references are shipped. The substitute finds, in the 79 k-family pangenome, the families that best mark each type. Because it is trained on the shipped labels it shows that the pangenome carries the information; it is not an independent rediscovery. R2 is the weak class (72.7 %).

**ABC capsule.** Named K-types agree for K1 (42 of 42), K15, K7, K96, K5 (19 of 19) and K4; K2 only 5 of 10. The shipped labels HP, CatB, Unknown and the numeric ones map one to one onto our class models. The default CapsuleFinder presence model was not run, only the custom typing models.

**Kaptive.** The figshare assembly ROAR131 has sequence glued to two header lines, which makes `makeblastdb` (inside Kaptive) abort and Kaptive skip every genome after it. A sanitised copy is used for Kaptive only (120 bases recovered, 182 contigs unchanged), and the script retries missing genomes. We call nine more confident types than the authors' manually curated table.

**Defense systems.** 137 subtypes on both sides, mean 9.04 systems per strain here vs 9.17 shipped, Pearson 0.972 per strain; exact count per strain 63.9 %, exact subtype-count vector 58.6 %. Systematic differences: RM_Type_IV (−123 systems), AbiE (+15), Thoeris_II (+14). The cause is not established (model or CasFinder version, hmmsearch version, the authors' own filtering, or Prodigal vs Prokka proteins).

**Outer-membrane proteins.** The paper says only "variants with > 99 % identity". The authors' exact clustering is not recovered: mean agreement is best at 99 % identity, but some genes fit other thresholds (tolC 0.81 at 97 %, fadL 0.95 at 95 %). YncD is not in Prokka's database; its K-12 reference (UniProt P76115) is listed in `config/references.yaml`.

**Interaction matrix.** The raw table is binary per replicate, so a plaque (3) cannot be told from full clearing (4). Comparable cells: 35,268 (369 strains × 96 phages; the 33 host strains are not in the raw table). The largest disagreement is 867 cells scored 2 in the shipped matrix and 1 here.

**Phage taxonomy.** VIRIDIC is not installed, so its similarity is approximated. The genus threshold (70 %) is reproduced well; the species threshold (95 %) is not, because the approximation under-counts identity between partly similar genomes.

**Not done.** Receptor-binding-protein detection (PhageRBPdetect, Pfam-A 34 and ESMFold are not practical on a 16 GB Mac without a GPU, and the authors curated manually). The shipped `RBP_list.csv` remains the reference. The phage network (vConTACT2) was not attempted.

## Effect on the models

The regenerated tables, tree-based UMAP and CV groups feed the models in [02_model_reproduction.md](02_model_reproduction.md) (pooled AUROC 0.855 vs 0.856 on the shipped inputs) and the ablation in [04_feature_necessity.md](04_feature_necessity.md).
