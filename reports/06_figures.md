# 6. Figures

Every figure has a PNG, an SVG and a CSV with the numbers behind it, in `results/figures/` (panels of the paper) and `results/comparison/` (comparison with the paper). Scripts: `fig_main_1_2.py`, `fig1a_tree.py`, `fig_main_3_6.py`, `fig_extended.py`, `fig_schematics.py`, `fig_comparison.py` (R fits: `fig4_cd.R`, `fig_ed6.R`). The class of each panel was decided from the legends of the published article:

- **A** replicable from the shipped tables, **B** after the preprocessing steps are re-run, **C** needs inputs not shipped, **D** schematic (redrawn as a diagram).

## Coverage

| Panel | Content | Class | Status | File |
|---|---|---|---|---|
| Fig 1a | Core-genome tree with typing rings | B | drawn from our regenerated tree | `fig1a_tree` |
| Fig 1b–j | Collection descriptors (ST, O, pathotype, source, LPS, ABC, outer-membrane proteins, defense) | A | drawn | `fig1_b_j` |
| Fig 1k | Network of 1,122 phage genomes | C | not drawn: needs 1,026 GenBank genomes and vConTACT2 | – |
| Fig 2a–b | Plaque-assay and MLC schematics | D | redrawn | `schematic_fig2ab` |
| Fig 2c | Interaction matrix in tree order | A | drawn | `fig2_c` |
| Fig 2d–e | Host-range histograms | A | drawn | `fig2_d_e` |
| Fig 2f | Host range by morphotype, ANOVA | A | drawn, F = 113.9 (paper 114.7) | `fig2_f` |
| Fig 3a | Host-range UMAP, PERMANOVA | A | drawn | `fig3_a` |
| Fig 3b | Genome map of six Vectrevirus phages | C | not drawn: needs annotation and alignments | – |
| Fig 4a–b | Significant terms and hits per adsorption factor | A | drawn from our MCMC re-run | `fig4_a_d` |
| Fig 4c–d | Defense count vs infections and vs mean MLC | A | drawn with our fits | `fig4_a_d` |
| Fig 5a, c, d | Per-phage AUROC, AUROC vs lytic count, prediction matrix | A | drawn | `fig5_a_c_d` |
| Fig 5b | Pooled ROC and precision–recall | A | drawn (Phase 4) | `results/repro/fig5_curves.png` |
| Fig 6a | Experiment schematic | D | redrawn | `schematic_fig6a` |
| Fig 6b | Tree of the 100 test strains | B | not drawn: test-strain genomes are not shipped | – |
| Fig 6c–e | Cocktails, success rate, binned MLC | A | drawn | `fig6_c_e` |
| ED1 | Lab-adapted strains | A | BE and BL21 drawn; MG1655 is not in the public matrix | `ed1` |
| ED2 | Nestedness and modularity | A | drawn with NODF and our BRIM (see below) | `ed2` |
| ED3 | RBP genes as variability hotspots | C | not drawn | – |
| ED4 | RBPs and host range | B | panel A drawn from the shipped matrix | `ed4a` |
| ED5 | PDP110_P1 and Klebsiella K55 | A | drawn | `ed5` |
| ED6 | Genome size and AMR genes vs MLC | A | panels A and C drawn; B and D need an AMR annotation | `ed6_a_c` |
| ED7 | Coverage of the most-covering k-uplet | A | drawn; value differs from the paper | `ed7` |
| ED8 | Four-step recommender | D | redrawn (box titles from the legend only) | `schematic_ed8` |
| ED9 | Cocktails split into generic and tailored | A | drawn | `ed9` |

The paper uses new names for four phages; the public matrix keeps the old ones. The mapping used (AL505_P3 = AL505_Sd2, T145_P2 = 412_P2, MT1B1_P3 = MT1B1_3A1, PDP110_P1 = 411_P1) was inferred from genus and isolation strain. Only the last is confirmed: 411_P1 infects exactly three strains, all with a K55 capsule, as the paper states.

ED2 uses NODF as the nestedness measure because the BINMATNEST temperature has no implementation here, and a BRIM modularity written for this work (1,000 null matrices with the same size and connectance). The values are therefore not comparable with the paper's.

## Main figures

### Fig 1: the collection

![Fig 1b-j](../results/figures/fig1_b_j.png)

Bar charts of the 403 strains by sequence type, O-antigen, pathotype, source, LPS outer-core type and ABC capsule serotype; the number of variants of each of the 12 outer-membrane proteins; defense subtypes per genome; and the 25 most frequent defense subtypes (137 subtypes in all). The strain count (403) matches the paper.

![Fig 1a](../results/figures/fig1a_tree.png)

Our own core-genome tree (PanACoTA and FastTree, rooted at the midpoint) with rings for phylogroup, sequence type, O-type and H-type taken from the regenerated feature table. The ring colours show the 12 most frequent values; everything else is grey.

### Fig 2: the interaction matrix

![Fig 2c](../results/figures/fig2_c.png)
![Fig 2d-e](../results/figures/fig2_d_e.png)
![Fig 2f](../results/figures/fig2_f.png)

The matrix of MLC scores 0–4 with strains in core-genome tree order and phages grouped by morphotype and genus; the share of phages that infect each strain and of strains each phage infects (mean 20.7 % of strains per phage); and host range by morphotype. Myoviruses infect 41.3 % of strains on average, podoviruses 11.2 % and siphoviruses 9.1 % (paper: 41.3, 11.1, 9.1). One-way ANOVA F = 113.9 (paper 114.7; the public matrix has 402 strains where the paper has 403).

### Fig 3: host-range similarity of the phages

![Fig 3a](../results/figures/fig3_a.png)

UMAP of the Jaccard distance between the binary lytic vectors of the 96 phages, coloured by genus (left) and by isolation strain (right). PERMANOVA adjusted R²: isolation strain 0.55 (paper 0.53), genus and isolation strain together 0.69 (0.65), genus alone 0.27 (0.13, not matched; cause not found).

### Fig 4: what drives infection

![Fig 4](../results/figures/fig4_a_d.png)

(a) Significant terms in the per-phage mixed models: 28 adsorption-factor terms and no defense-system term. (b) Terms per adsorption factor, this re-run against the authors' shipped adsorption-factor-only model (LPS 14 vs 9, ABC capsule 12 vs 13). (c, d) Number of defense systems against lytic interactions per strain (β −0.233, pMCMC 0.40; paper −0.22, 0.41) and against average MLC (β −0.033, pMCMC 0.0003; paper −0.04, 0.004).

### Fig 5: prediction performance

![Fig 5](../results/figures/fig5_a_c_d.png)
![Fig 5b](../results/repro/fig5_curves.png)

(a) AUROC per phage (mean 0.787); (c) AUROC against the number of lytic interactions of the phage (Spearman −0.35: phages with few lytic strains are predicted better); (d) predicted against observed lysis over the 38,436 predicted pairs at threshold 0.5 (4,669 true positives, 2,940 false positives, 3,307 false negatives). (b) Pooled ROC and precision–recall curves of this run against the authors' shipped predictions. Comparison with the paper's numbers: [05](05_paper_comparison.md).

### Fig 6: the cocktail experiment

![Fig 6a](../results/figures/schematic_fig6a.png)

(a) Schematic redrawn from the legend (box text from the legend only). Panel (b), the tree of the 100 test strains among the 403, is not drawn because the test-strain genomes are not shipped.

![Fig 6c-e](../results/figures/fig6_c_e.png)

(c) The 18 distinct cocktails and the number of test strains each was recommended for: four generic cocktails cover 76 strains (31, 22, 19 and 4), fourteen tailored cocktails cover 24. (d) In vitro success rate: tailored 91.67 %, generic 78.95 %, baseline (the most-covering triplet) 81 %, identical to the paper. (e) Distribution of the binned MLC of the three-phage cocktail: tailored against baseline Mann–Whitney U = 399.5, P = 0.019; tailored against generic U = 1171.5, P = 0.034, both identical to the paper's legend. The recommender itself was not re-run; these panels re-analyse the shipped experimental results.

## Extended Data

### ED1: laboratory-adapted strains

![ED1](../results/figures/ed1.png)

MLC of the 96 phages on BE and BL21: 45 phages infect BE (46.9 %) and 48 infect BL21 (50.0 %), as in the paper's legend. MG1655, tested at a single concentration, is not in the public matrix.

### ED2: nestedness and modularity

![ED2](../results/figures/ed2.png)

Left: the matrix sorted by degree. Middle: NODF 52.1 against 21.5 ± 0.25 for 1,000 random matrices of the same size and connectance (P = 0.001). Right: BRIM modularity Q = 0.231 with 5 modules against 0.124 ± 0.003 (P = 0.001). The paper uses the BINMATNEST temperature for nestedness, so the first number is not comparable; the paper does not print its P values.

### ED4: receptor-binding proteins and host range

![ED4A](../results/figures/ed4a.png)

Panel A: the host range of three phages with a similar receptor-binding protein, over the strains in tree order, with the K1 ABC capsule marked in the top row. AL505_Sd2 infects 77 strains (36 K1), 412_P2 infects 53 (36 K1) and MT1B1_3A1 infects 68 (37 K1); 42 strains carry a K1 capsule. The phage names are the public names of AL505_P3, T145_P2 and MT1B1_P3.

### ED5: PDP110_P1 and the Klebsiella K55 capsule

![ED5](../results/figures/ed5.png)

PDP110_P1 (public name 411_P1) infects three strains, and these are exactly the three strains with a K55 capsule locus, as in the paper.

### ED6: genome size and mean MLC

![ED6](../results/figures/ed6_a_c.png)

Genome size (from our assemblies) against lytic interactions per strain (A, β −5.9 per Mb, pMCMC 0.12) and against average MLC (C, β −0.35 per Mb, pMCMC 0.004). Larger genomes go with a lower average MLC, consistent with the paper's title. Panels B and D (AMR genes) are not drawn.

### ED7: coverage by the most-covering phages

![ED7](../results/figures/ed7.png)

Share of strains lysed by at least one phage of the most-covering k-set (exact for k ≤ 3, greedy above): 65.4 % for one phage, 80.3 % for two, 87.6 % for the triplet (LF82_P8, NIC06_P2, 536_P7), 96.8 % for ten, 97.0 % for 25. The paper gives 63 % for the triplet; its definition of "productive lysis" is not in the legend, so the difference is unresolved.

### ED8: the four-step recommender

![ED8](../results/figures/schematic_ed8.png)

Redrawn from the legend; only the four box titles come from the figure, since the legend gives no further detail.

### ED9: generic and tailored cocktails

![ED9](../results/figures/ed9.png)

(A) Step of the pipeline at which each recommended phage was chosen. For tailored cocktails 76 % come from step 1 (similar to the isolation strain); for generic cocktails 39 % do, and 59 % come from the last two steps. Mann–Whitney P = 4e-9 (paper 3.5e-12, same direction, weaker). (B) Haplotypes (sequence type, O-type and H-type) covered by each cocktail: the four generic cocktails cover 4 to 20, the fourteen tailored ones 1 or 2.

## Schematics

![Fig 2a-b](../results/figures/schematic_fig2ab.png)

Fig 2a–b: the plaque assay at three phage concentrations and the MLC score (lowest concentration with lysis), redrawn from the legend.

## Comparison with the paper and with the shipped results

![Headline numbers](../results/comparison/cmp_headline.png)
![Pooled ROC and precision-recall](../results/comparison/cmp_roc_pr.png)
![Per-phage AUROC](../results/comparison/cmp_per_phage_auroc.png)
![Feature agreement](../results/comparison/cmp_feature_agreement.png)
![Ablation](../results/comparison/cmp_ablation.png)

## Colour and checks

Two series use the blue `#2a78d6` and orange `#eb6834` of the reference palette; every two-series chart has a legend, no chart has two y-axes, and grid and axes are muted. The palette validator of the visualisation skill could not run on the test machine (its Node binary lacked a library), so the palette's documented, previously validated slots were used and each figure was inspected by eye.
