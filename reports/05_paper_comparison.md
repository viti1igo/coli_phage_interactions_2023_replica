# 5. Comparison with the paper

Three sources are compared: the numbers printed in the paper (Gaborieau et al., *Nature Microbiology* 2024, and its preprint) together with the publisher's Supplementary Information, Supplementary Tables 1–7 and Source Data files; the authors' shipped results in their repository; and this replication, run once on the shipped inputs and once on inputs regenerated from the genomes. Every paper number carries its place in the article in the `paper_source` column of [`results/comparison/paper_vs_replication.csv`](../results/comparison/paper_vs_replication.csv).

The publisher's files (13 files, 8.7 MB) are downloaded by `scripts/download_supplement.sh` and are not redistributed here. `scripts/compare_source_data.py` compares them with this replication.

## What the supplement settled

| Open question | Answer from the supplement |
|---|---|
| Why 402 strains here and 403 in the paper? | The public matrix lacks strain LF110; Supplementary Table 3 has it. The other 402 strains are identical in both files |
| New phage names in the paper (PDP21_P1, CLB_P2, …) | Matched one to one with the public names by identical interaction vectors (`results/comparison/phage_name_map.csv`). The earlier guess of the mapping was right, and PDP110_P1 = 411_P1 is confirmed |
| ANOVA F of Fig 2f (113.9 here, 114.7 in the paper) | With all 403 strains F = 114.66, as in the paper and in its Source Data |
| The 85 % average AUROC of Fig 5a | Fig 5a plots one AUROC per **strain** (Source Data). The paper's own Source Data averages 0.872 (median 0.892) and this run gives 0.881 (0.895) on the same 315 strains; the text's 85 % is below the mean of its own data. The preprint's 77 % is the per-phage mean (0.787 here) |
| Genus PERMANOVA R² (0.27 here, 0.13 in the paper) | The 0.13 matches the increase of the adjusted R² when genus is added after isolation strain (0.686 − 0.554 = 0.13). That is our reading; the paper does not say so |
| BRIM modularity (ED2) | Supplementary Text 5 gives Q = 0.23 (P < 1e-3) and a BINMATNEST temperature T = 11.63; this run gives Q = 0.233 with four modules. T has no counterpart here (NODF was computed) |
| Fig 4a counts (28 adsorption-factor terms here, 52 in the paper) | The paper's Source Data counts terms up to pMCMC 0.0089, i.e. a threshold of 0.01, not the 1e-3 of the Methods. At 0.01 this re-run gives 51 adsorption-factor terms and 3 defense terms (paper 52 and 1; its single defense term is the outer-membrane variant YNCDOther, which its Source Data files under defense systems). Signs of the shared terms agree (31 of 31); only 19 of the 53 (phage, term) pairs are significant in both, because terms close to the threshold change from run to run |
| ED6 B and D (AMR genes) | The AMR counts per strain are in the Source Data (not in the public repository). With them: genome size against average MLC pMCMC 0.003 (significant), AMR against average MLC 0.0105 and against lytic interactions 0.047 (not significant at 0.01), genome size against lytic interactions 0.12 |
| ED7, the definition of "productive lysis" | Not settled; see below |
| Recommender code | Not in the supplement; see below |

## Headline numbers

| Item | Paper | Authors' shipped | This run (shipped inputs) | This run (regenerated inputs) |
|---|---|---|---|---|
| Mean AUROC per strain (Fig 5a) | 0.85 text; 0.872 Source Data | 0.881 | 0.881 (0.878 over all 390 strains with both classes) | 0.877 |
| Mean AUROC per phage (preprint) | 0.77 | 0.781 | 0.787 | 0.765 (0.769 with our matrix) |
| Pooled AUROC (Fig 5b) | 0.86 text; 0.849 from the curve points | 0.857 | 0.856 | 0.855 (0.853) |
| Average precision, F1 | 0.60, 0.60 | 0.601, 0.596 | 0.598, 0.599 | 0.596, 0.592 |
| Cocktail success tailored / generic / baseline (Fig 6d) | 91.67 / 78.95 / 81 % | same data | 91.67 / 78.95 / 81 % | not applicable |
| Mann–Whitney (Fig 6e) | U 399.5, P 0.019; U 1171.5, P 0.034 | same data | identical | not applicable |
| Defense count vs lytic interactions, β (Fig 4c) | −0.22 (CI −0.74 to 0.34) | – | −0.233 | −0.235 |
| Defense count vs average MLC, β (Fig 4d) | −0.04 (CI −0.05 to −0.01) | – | −0.033 | −0.037 |
| ANOVA, host range by morphotype (Fig 2f) | F = 114.7 | – | F = 114.7 | – |
| PERMANOVA adj. R², genus added after isolation strain / isolation strain / both | 0.13 / 0.53 / 0.65 | – | 0.13 / 0.55 / 0.69 | – |
| BRIM modularity Q (ED2) | 0.23 | – | 0.233 | – |
| PDP110_P1 infected strains, all K55 (ED5) | 3, all K55 | – | 3, all K55 | – |

The table has 97 rows in all: 33 match within 5 %, 12 match, 1 matches once explained (the 403rd strain), 4 are consistent with the paper's confidence interval or item list, 23 differ, and 24 cannot be compared numerically (the paper gives no number, or the panel is qualitative). Rows 74 to 97 come from the supplementary analyses and the genome maps of [10](10_repository_analyses.md). The chart shows each comparable number as our value divided by the paper's.

![Headline numbers against the paper](../results/comparison/cmp_headline.png)

Data: `results/comparison/cmp_headline.csv`.

## Pooled ROC and precision–recall curves

The three curves of the authors' shipped predictions, this run on the shipped inputs and this run on the regenerated inputs lie on top of each other (AUROC 0.857, 0.856, 0.855; average precision 0.601, 0.598, 0.596). The dashed black line is the curve exported in the paper's Source Data, which gives an AUROC of 0.849 from its points (the text says 0.86; the exported precision–recall block holds only the tail, recall 0.945 to 1, so it is not drawn).

![Pooled ROC and precision-recall curves](../results/comparison/cmp_roc_pr.png)

Data: `results/comparison/cmp_roc_pr.csv` (without the paper's points, which are the publisher's data). Per-phage AUROC against the authors' shipped values (correlation 0.92 on shipped inputs, 0.79 on regenerated ones):

![Per-phage AUROC](../results/comparison/cmp_per_phage_auroc.png)

Redrawn panels next to their paper numbers (per-strain AUROC with the paper's Source Data points in orange, cocktails, mixed models and defense count, host range by morphotype, host-range UMAP):

![Per-strain AUROC, AUROC vs lytic interactions, prediction matrix](../results/figures/fig5_a_c_d.png)
![Cocktail experiment](../results/figures/fig6_c_e.png)
![Mixed-model results and defense count](../results/figures/fig4_a_d.png)
![Host range by morphotype](../results/figures/fig2_f.png)
![Host-range UMAP](../results/figures/fig3_a.png)

## What still does not match

| Item | Difference | What is known |
|---|---|---|
| ED7, coverage of the most-covering k-sets | Source Data: 42.7 % for one phage, 67.0 % for the triplet, 73.7 % for four, 89.3 % at the plateau; the text says 63 % and 68 % for three and four; this run gives 65.5, 87.6, 90.8 and 97.0 % with lysis at any concentration | The paper contradicts itself (text 63 % against Source Data 67 %). No set of MLC scores reproduces the Source Data (best single phage 172 of 403 strains there; any MLC > 0: 264; MLC 2–4: 175; plateau 360 strains against 391 with any lysis). The matrix used for ED7 may differ from the final one |
| ED9b, "haplotypes per cocktail" | Source Data 31, 21, 19, 4; this run 20, 18, 18, 4 | The Source Data values equal the strain counts of Fig 6c, so they count strains, not haplotypes as the legend says |
| ED9a Mann–Whitney P | 4e-9 here, 3.5e-12 in the paper | The step shares are identical (tailored 76.4 %, generic 38.7 % at step 1); the test unit is not stated |
| Cocktail labels (Fig 6) | Supplementary Table 7 labels 28 strains as tailored and 72 as generic; the 24/76 split used here gives the paper's 91.67 % (22/24) and 78.95 % (60/76) | Table 7 marks the four-strain cocktail 536_P7-DIJ07_P2-LF82_P8 as tailored; the figure treats it as generic. The Source Data of Fig 6c lists 99 strains where the cocktail file has 100 |
| Fig 4d pMCMC | 0.0003 here, 0.004 in the paper | Both significant; the estimates agree |
| Combined PERMANOVA (0.69 vs 0.65) | 6 % | Own implementation; the adjusted-R² formula may differ |
| Fig 2f P value | 7.9e-26 here, 1.44e-26 printed | The Source Data frequencies give 7.9e-26 too, so the printed P is not the one-way ANOVA P of these data |
| Genome sizes (ED6) | Our assemblies are on average 69 kb smaller than the paper's (correlation 0.99) | Probably a different assembly version |
| ED2 nestedness | not comparable | NODF here, BINMATNEST temperature (11.63) in the paper |

## Cocktail recommender

The recommender's code is not published and is not in the supplement. The authors' repository describes it in its README (sections 2.3 to 2.6), which differs slightly from the Methods: step 1 compares O-type, ST, ABC, Klebsiella and LPS (the Methods say ST, O, H, LPS, ABC), and step 2 keeps phages that infect at least half of the group (the Methods say more than half). It was re-implemented from the README (`scripts/recommender.py`, `SPEC=readme`; `SPEC=paper` follows the Methods and gives identical cocktails here): a random forest (depth 4, 80 trees, positive class weight 1.5) on the phage identity and five "same as the isolation strain" indicators; the strains of the Picard collection with the same ST, O and H and the phages that infect most of them; the same forest on all bacterial features; the most generalist phages; then the first three allowed phages with probability above 0.5, never two of the same genus and isolation host. The test strains' genomes are not shipped, so their ABC and Klebsiella types are unknown (Supplementary Table 6 gives ST, O, H, LPS and phylogroup only) and their UMAP coordinates were recomputed from the shipped distance matrix. One test strain (AN03) has the name of a Picard strain and has no UMAP coordinates.

The re-implementation does not reproduce the authors' cocktails: 2 of 100 are identical and on average 0.97 of the 3 phages are shared (75 cocktails share exactly one phage). The haplotype-neighbour step is reproduced (10 of the 11 phages the authors chose at step 2 are also recommended by our step 2), the two forest steps are not (the authors' forests give probabilities of 0.5 to 0.66; ours stay below 0.5 unless every feature is considered at every split, which was therefore adopted, and then only 8 % of phages come from step 3 against 26 %). The authors' step-4 phages are LF82_P8 and 55989_P2; ranking phages by lysis at plaque or clearing (MLC 3–4) puts them second and fourth, ranking by lysis at any concentration puts them seventh and sixth. Four variants (one-hot or integer phage identity, all or √n features per split) are in `results/recommender/variants.csv`; none gets close.

![Recommender against the authors' cocktails](../results/comparison/cmp_recommender.png)

Data: `results/comparison/cmp_recommender.csv`, `results/recommender/`.

## Regenerated inputs and necessity of features

How well each input was regenerated is in [03](03_feature_regeneration.md); the ablation on both input sets is in [04](04_feature_necessity.md). In short, the pooled AUROC survives regeneration, and the verdicts on the ABC capsule, LPS and the O-type/ST pair are the same on both input sets; three of twelve verdicts differ.

![Agreement of regenerated features with the shipped ones](../results/comparison/cmp_feature_agreement.png)
![Feature ablation on shipped and regenerated inputs](../results/comparison/cmp_ablation.png)

## Not drawn, and why

Fig 1k (1,026 GenBank genomes and vConTACT2), Fig 6b (test-strain genomes are not shipped), ED1 MG1655 (not in the public matrix). Supplementary Tables 1–7 contain no accession list for the phage network and no genome identifiers for the test strains beyond Supplementary Table 6.

## Repository analyses and differences

The notebooks of the authors' repository, the supplementary analyses and the genome maps (Fig 3b, ED3) are in [10](10_repository_analyses.md); every difference between the repository and the paper that turned up is in [11](11_repository_vs_paper.md). Most of the 23 rows that still differ are differences inside the paper or between the repository and the paper (ED7, ED9b, the cocktail labels, the K1 counts, the RBP list), not failures of this replication.

## Configuration differences

28 differences between the paper's settings and what was run, each with a reason and an expected effect, are in [07](07_configuration.md) and `config/config_diff.csv`.

## Files compared and files still needed

| Origin | File | Status |
|---|---|---|
| Paper | preprint and published article (PDF and text) | used; not redistributed here |
| Publisher | Supplementary Information, Reporting Summary, Peer Review File | downloaded (`scripts/download_supplement.sh`) |
| Publisher | Supplementary Tables 1–7 (one workbook) | downloaded |
| Publisher | Source Data for Fig 1–6, ED6, ED7 and ED9 (nine workbooks) | downloaded; none for ED1, ED2, ED4 or ED5 |
| Authors' repository | interaction matrix (402 × 96), predictions, cocktail data, MCMC results | downloaded; also checked: 14 commits, none holds the cross-validation groups or the recommender |
| Zenodo 10.5281/zenodo.13831957 | archive of the same repository at publication (146.8 MB) | not downloaded; described as a copy of the GitHub repository |
| figshare | 403 bacterial assemblies | downloaded |
| Not available | cross-validation group file; recommender code; test-strain genomes; 1,026 GenBank phage genomes | would have to come from the authors |

List: `results/comparison/file_inventory.csv`.

## Limits

Single seed per model run; UMAP layouts differ run to run; the recommender is a best-effort re-implementation; the paper's ED7 curve could not be derived from the matrix.
