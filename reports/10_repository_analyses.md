# 10. Everything that can be run from the authors' repository

This report covers what I could run from the authors' repository, either as it is or completed from its own data: the four notebooks, the supplementary analyses whose inputs the repository holds, and the genome maps. [11_repository_vs_paper.md](11_repository_vs_paper.md) collects every difference between the repository and the paper.

## 10.1 The four notebooks

`scripts/run_notebooks.py` copies the authors' `data/` and `dev/` folders to a scratch directory, which keeps their relative paths valid and leaves the clone untouched. It runs each notebook with `nbclient` and records every cell that fails. The `shipped` mode runs the notebooks unchanged. The `patched` mode applies the edits in `scripts/notebook_patches.py`. Those edits swap a Windows or private path for its public equivalent or fix an incompatibility with current pandas and matplotlib, and none of them changes an analysis. `notebooks/` holds the executed notebooks, the per-cell error tables, the files the notebooks write and the figures from their cells.

| Notebook | Code cells | Failing cells as shipped | Failing cells after patches | Produces |
|---|---|---|---|---|
| `bokeh_draw_interaction_matrix` | 3 | 0 | 0 | interactive heat map of the 402 × 96 matrix (the hover images point to private plaque pictures) |
| `visualize_predictions` | 18 | 15 | 1 | Fig 5 panels from the authors' shipped predictions |
| `analyze_cocktail_interactions_..._baseline` | 29 | 10 | 0 | Fig 6 and ED9 panels and tests |
| `build_bacteria_arsenale_pairwise_distance_table` | 11 | 2 | 0 | defense-arsenal Jaccard distance table |

Why the cells fail as shipped (`notebooks/shipped_errors.csv`):
- Windows paths: backslashes or paths on the author's disk (`D:\These\...`, `C:/Users/...`). `visualize_predictions` fails from such paths onward, and so do the figure-saving cells of the cocktail notebook.
- A private input: the notebooks reference the matrix `ADSORPTION_MATRIX_370+host.csv` (the N/P/F/U matrix of the original script) and the plaque pictures, and neither ships. The notebooks run on the public MLC matrix because those lines are commented out.
- A moved file: `../panacota/370+host_leaf_order_root=albertii.txt` now lives in `panacota/tree/`.
- A missing package: `bokeh_sampledata` (needed since Bokeh 3).
- Newer library behaviour: a figure saved under a `.tsv` name, positional indexing on a labelled Series, `fillna("")` on a numeric column.

The one cell that still fails after the patches is the authors' own: `visualize_predictions` cell 12 plots `n_infections` from `averaged_per_phage`, a table that does not have that column (`to_plot` does).

### What the notebooks print

`visualize_predictions` (the authors' shipped predictions, best model per phage, test folds). True positives 4,648, false positives 2,980, true negatives 27,432, false negatives 3,318; precision 0.609, recall 0.583, F1 0.596, average precision 0.601, AUROC 0.857. These are the paper's 86 %, 60 % and 60 %, and close to this replication (4,669 / 2,940 / 3,307; AUROC 0.856). The notebook also compares the train and test AUROC of every phage model. The mean gap is 16.3 AUROC points, and 27 of the 96 phages show a significant train-test difference after correction. The paper does not report this analysis. It says only that the authors reduced the number of features to limit overfitting.

![Feature importance of the authors' models (their notebook)](../notebooks/patched_images/visualize_predictions/cell06_1.png)

Annotation: Feature importance of the authors' per-phage models, one box per feature, sorted. The highest are specific O types, sequence types and LPS types.

![ROC curve of the shipped predictions (their notebook)](../notebooks/patched_images/visualize_predictions/cell16_8.png)

Annotation: Pooled ROC of the authors' shipped predictions (AUROC 86 %). The dashed line is random.

Cocktail notebook. Success 81.0 % (baseline), 78.95 % (generic), 91.67 % (tailored) with 100, 76 and 24 strains; tailored against baseline on tailored strains U = 399.5, P = 0.019 (paired t-test P = 0.030); tailored against baseline on all strains U = 1543.5, P = 0.028; the step-of-recommendation test of ED9 gives U = 3856.5, P = 3.5e-12 (see [11](11_repository_vs_paper.md): each strain enters it three times). The first panel counts the recommended cocktails.

![Recommended cocktails (their notebook)](../notebooks/patched_images/analyze_cocktail_interactions_/cell04_1.png)

Annotation: Strains per recommended cocktail in the authors' notebook (31, 21, 19, 4, then 3 to 1). It counts 99 strains because their merge loses BCH857.

Defense-arsenal notebook. It builds the Jaccard distance between the defense arsenals of the 404 strains in the defense tables and draws it in tree order. It wrote no file in this run.

## 10.2 Supplementary analyses computable from the repository

`scripts/si_analyses.py` (Python) and `scripts/si_rbp_lmm.R` (R) redo the analyses of the Supplementary Texts whose data the repository holds. The paper numbers come from the Supplementary Information, and the results are in `results/si/`.

| Supplementary Text | Analysis | Paper | This run (shipped inputs) | This run (regenerated inputs) |
|---|---|---|---|---|
| 4 | Lytic interactions of a phage on genome size, morphotype in the model (96 phages) | genome-size coefficient 2.3e-4, P = 0.17 | 2.32e-4, P = 0.17 (without the morphotype the coefficient is 1.2e-3, P = 3e-19: morphotype confounds genome size) | 2.2e-4, P = 0.16 |
| 6 | Decision rule "lytic if and only if the strain carries a K1 ABC capsule", the three phages T145_P2, AL505_P3, MT1B1_P3 | 45 K1 strains; precision 73 %, recall 72 % | 42 K1 strains; precision 87 %, recall 55 % | 45 K1 strains in the whole collection (41 in the 369 strains); precision 89 %, recall 65 % |
| 6 | The same rule on the other 28 Vectreviruses and Kayfunaviruses | precision 7 %, recall 14 % | precision 18 %, recall 18 % | precision 17 %, recall 20 % |
| 7 | LMM of the lytic interactions on the number of RBPs, phage tree as random effect | β = 28.11, 95 % CI -33.09 to 83.58, pMCMC 0.33 | β = -13.0 (CI -56.7 to 27.7, pMCMC 0.52) for all 96 phages; β = 16.0 (CI -39.4 to 72.6, pMCMC 0.57) for the 89 phages with an RBP. Not significant in either, as in the paper | β = -15.4 (CI -55.0 to 19.8, pMCMC 0.41) for all 96; β = 10.4 (pMCMC 0.67) for the 89 |
| 8 | PERMANOVA on the bacterial infection profiles, 8 phylogenetic PCoA axes first, then adsorption factors (AF) and/or defense subfamilies; adjusted R² (AIC) | AF 0.27 (-435), full 0.28 (-418), defense 0.05 (-372) | AF 0.26 (-458), full 0.27 (-437), defense 0.17 (-427) | AF 0.25 (-411), full 0.28 (-409), defense 0.14 (-395); 369 strains, our defense subtypes instead of the 80/80 subfamilies |

Notes:
- Supplementary Text 8 comes out in the same direction as the paper. Adsorption factors fit best (lowest AIC), the full model adds almost nothing, and defense subfamilies explain much less than adsorption factors. The adjusted R² of the adsorption-factor and full models agrees with the paper. The defense-only R² is 0.17 against 0.05, and the AIC values sit on a lower scale. Covariate definitions are not fully specified in the paper (traits with fewer than 8 observations are lumped); `results/si/permanova_variants.csv` gives 18 variants (minimum 8, 15 or 30 strains per trait, defense subfamilies or subtypes). With a minimum of 30 strains the defense R² falls to 0.09 to 0.10, still above 0.05. The distance is the Jaccard distance between the binary lytic profiles of the 402 strains of the tree; the AIC is n ln(SS_res / n) + 2·rank.
- Supplementary Text 6 (K1): the repository table has 42 K1 strains where the text says 45, with the four isolation strains included in both. Pooling the three phages over all strains with a K1 call gives a higher precision and a lower recall than the paper. The paper's recall of 72 % would need about 65 strains per phage to count as positive, so its counting differs from the 45 strains it names.
- Supplementary Text 7: the repository's `RBP_list.csv` lists 7 phages with no RBP (536_P1, 536_P6, 536_P7, 536_P9, DIJ07_P1, DIJ07_P2, LF31_P1), while the text says every phage encodes one or two. I fitted both versions above.

### Does RBP similarity explain host-range similarity?

The repository has the RBP gene identifiers but not their sequences. I located each RBP gene as the n-th gene of its genome, the way Prokka numbers them. 127 of 130 identifiers mapped: 107 with single-genome gene calling, and the rest with metagenomic calling where the identifier exceeded the gene count. I then clustered the proteins by pairwise identity. These are the results over the 3,916 pairs of the 89 phages with an RBP:

| RBP identity threshold | Pairs sharing an RBP | Host-range distance, sharing / not sharing | Pairs sharing across genera |
|---|---|---|---|
| 30 % | 267 | 0.61 / 0.85 | 27 |
| 50 % | 217 | 0.56 / 0.85 | 13 |
| 70 % | 182 | 0.53 / 0.85 | 1 |

Phages that share an RBP have much more similar host ranges (Mann-Whitney P < 1e-60 at every threshold). Most of the sharing happens within a genus, so I cannot separate it from genus here. Across genera only 13 pairs share an RBP at 50 % identity (distance 0.78 against 0.85 for the other cross-genus pairs). That agrees with the paper's claim, but the test does not stand apart from genus (`results/si/rbp_similarity_vs_host_range.csv`). The regenerated matrix gives the same numbers (0.56 against 0.85).

Regenerated inputs. `MATRIX_SET=regenerated python scripts/si_analyses.py` (and `Rscript scripts/si_rbp_lmm.R`) redo the analyses with the matrix, features, tree and defense subtypes regenerated from the genomes (results in `results/si_regenerated/`); the RBP list and the phage genomes stay the authors'. Our ABC calls give 45 K1 strains in the whole collection, the paper's count, where the shipped table has 42.

## 10.3 Genome maps (Fig 3b and ED3)

The authors' annotation is not shipped, so I annotated the 96 phage genomes again with `pyrodigal` (`scripts/fig_phage_maps.py`), and aligned them with MAFFT 7.526 in a small micromamba environment. The RBP genes are the ones named in `RBP_list.csv`, located by their number in the gene order. Fig 3b uses the six Vectrevirus phages of the paper's Source Data (LI10_P4, LI10_P5, 409_P5, 409_P3, LF31_P3, 412_P2, isolated on four hosts). ED3 uses up to four phages from different hosts for each of the five genera.

![Fig 3b](../results/figures/fig3b_vectrevirus.png)

Annotation: Top left: genome maps of six Vectrevirus phages. Light blue arrows are RBP genes and blue lines join protein matches of at least 30 % identity. Bottom left: entropy of the genome alignment, in bits. Right: host range of the six phages over the strains in tree order.

Protein links join the reciprocal best matches with at least 30 % identity between neighbouring genomes. The paper's legend says similarity, and I used identity. Position-wise entropy is the Shannon entropy of the nucleotide distribution in the six-genome alignment, averaged over 200 nucleotides. The RBP genes stand out. Their mean entropy is 0.62 bits against 0.24 bits in the other genes (alignment mean 0.25), and the peak lies in the RBP region. That supports the legend's statement that RBP genes are hotspots of variability. Data: `results/figures/fig3b_*.csv`.

![ED3](../results/figures/ed3_genus_maps.png)

Annotation: One panel per genus with up to four phages from different hosts. Arrows are genes, light blue arrows are RBP genes, and lines join protein matches of at least 30 % identity.

In each of the five genera the genes are conserved and linked across the phages, except the RBP genes (light blue), which have few or no links. Data: `results/figures/ed3_*.csv`. Two caveats apply. My gene calls differ from the authors' Prokka calls, so the position of the highlighted RBP genes rests on the numbering convention. I also picked up to four phages per genus by name from different isolation hosts, which may not match the authors' choice.
