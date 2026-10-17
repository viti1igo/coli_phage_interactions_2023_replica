# 5. Comparison with the paper

Three sources are compared: the numbers printed in the paper (Gaborieau et al., *Nature Microbiology* 2024, and its preprint), the authors' shipped results in their repository, and this replication, run once on the shipped inputs and once on inputs regenerated from the genomes. Every paper number carries its place in the article in the `paper_source` column of [`results/comparison/paper_vs_replication.csv`](../results/comparison/paper_vs_replication.csv).

## Headline numbers

| Item | Paper | Authors' shipped | This run (shipped inputs) | This run (regenerated inputs) |
|---|---|---|---|---|
| Mean per-phage AUROC (Fig 5a) | 0.85 published, 0.77 preprint | 0.781 | 0.787 | 0.765 (0.769 with our matrix) |
| Mean per-strain AUROC | – | 0.881 | 0.878 | – |
| Pooled AUROC (Fig 5b) | 0.86 | 0.857 | 0.856 | 0.855 (0.853) |
| Average precision, F1 | 0.60, 0.60 | 0.601, 0.596 | 0.598, 0.599 | 0.596, 0.592 |
| Cocktail success tailored / generic / baseline (Fig 6d) | 91.67 / 78.95 / 81 % | same data | 91.67 / 78.95 / 81 % | not applicable |
| Mann–Whitney (Fig 6e) | U 399.5, P 0.019; U 1171.5, P 0.034 | same data | identical | not applicable |
| Defense count vs lytic interactions, β (Fig 4c) | −0.22 (CI −0.74 to 0.34) | – | −0.233 | −0.235 |
| Defense count vs average MLC, β (Fig 4d) | −0.04 (CI −0.05 to −0.01) | – | −0.033 | −0.037 |
| ANOVA, host range by morphotype (Fig 2f) | F = 114.7 | – | F = 113.9 | – |
| PERMANOVA adj. R², isolation strain (Fig 3a) | 0.53 | – | 0.55 | – |
| PERMANOVA adj. R², genus and isolation strain | 0.65 | – | 0.69 | – |
| PERMANOVA adj. R², genus | 0.13 | – | **0.27** | – |
| Coverage of the best triplet (ED7) | 63 % | – | **87.6 %** | – |
| Strains infected by PDP110_P1, and with K55 (ED5) | 3, all K55 | – | 3, all K55 | – |

The table has 48 rows in total: 20 match within 5 %, 2 lie inside the paper's confidence interval, 1 p-value has the same order of magnitude, 7 differ, and 18 cannot be compared numerically (the paper gives no number, or the panel is qualitative). Chart: `cmp_headline` ([06](06_figures.md)).

## What does not match

| Item | Difference | What is known |
|---|---|---|
| Fig 5a, published 85 % | ours 0.787 per phage, 0.878 per strain | No aggregation tried reproduces 85 %; the authors' own shipped predictions give 0.781 and 0.881, so the gap is not caused by this run. The preprint's 77 % matches. The text and legend do not say how 85 % was computed |
| PERMANOVA, genus (0.27 vs 0.13) | twofold | Cause not found. The isolation-strain R² is close (0.55 vs 0.53). The implementation is my own (Python, Jaccard distance, 999 permutations) rather than vegan's `adonis2` |
| PERMANOVA, genus and isolation strain (0.69 vs 0.65) | 6 % | Same implementation caveat |
| ED7 triplet coverage | 87.6 % vs 63 % | "Productive lysis" is not defined in the legend. MLC > 0, ≥ 2 and ≥ 3 give 87.6, 76.9 and 54.7 %; none gives 63 %. The Supplementary Text is probably needed |
| ED9A Mann–Whitney P | 4e-9 vs 3.5e-12 | Same direction and much weaker; the unit of analysis (phage or cocktail) is not stated |
| Fig 4d pMCMC | 0.0003 vs 0.004 | Both significant; the estimates agree (−0.033 vs −0.04) |
| Fig 2f ANOVA F | 113.9 vs 114.7 (0.7 %, counted as a match) | The public matrix has 402 strains where the paper has 403 |
| ED2 | not comparable | NODF is used for nestedness (the BINMATNEST temperature has no implementation here) and a BRIM written for this work; the paper does not print its P values |
| ED6 B and D | not run | No AMR-gene annotation is shipped. Panels A and C were run: genome size vs average MLC is significant (pMCMC 0.004), vs lytic interactions it is not (0.12), consistent with the paper's title that genome size correlates with virulence |

## Not drawn, and why

Fig 1k (1,026 GenBank genomes and vConTACT2), Fig 3b and ED3 (phage annotation and alignments), Fig 6b (test-strain genomes not shipped), ED1 MG1655 (not in the public matrix). The paper uses new names for four phages; the mapping to the public names was inferred from genus and isolation strain and is confirmed only for PDP110_P1 = 411_P1, which infects exactly the three K55 strains ([06](06_figures.md)).

## Regenerated inputs and necessity of features

How well each input was regenerated is in [03](03_feature_regeneration.md) (`cmp_feature_agreement`); the ablation on both input sets is in [04](04_feature_necessity.md). In short, the pooled AUROC survives regeneration, and the verdicts on the ABC capsule, LPS and the O-type/ST pair are the same on both input sets; three of twelve verdicts differ.

## Configuration differences

24 differences between the paper's settings and what was run, each with a reason and an expected effect, are in [07](07_configuration.md) and `config/config_diff.csv`. The ones that matter for interpretation: scikit-learn 1.9.1 for 1.1.2; cross-validation groups rebuilt; lytic defined as MLC above 0 on the public 0–4 matrix; `mlst` for SRST2; custom ABC models without the default CapsuleFinder model; own methods for LPS and outer-membrane variants; a BLAST similarity for VIRIDIC; an own PERMANOVA for vegan; NODF for the BINMATNEST temperature; an extra N(0, 10²) prior for the five phages whose MCMC failed as singular.

## Files compared and files still needed

| Origin | File | Status |
|---|---|---|
| Paper | preprint and published article (PDF and text) | used; not redistributed here |
| Paper | Supplementary Tables 1–6 | **not available** |
| Paper | Supplementary Text 4–5 and Supplementary Figures 1–9 | **not available** |
| Paper | Source Data for Fig 1–6 and ED 1–9 | **not available** |
| Authors' repository | interaction matrix (402 × 96), predictions, cocktail data, MCMC results | downloaded |
| Authors' repository | cross-validation group file; recommender code | not public / not in the repository |
| figshare | 403 bacterial assemblies | downloaded |
| Not obtained | 1,026 GenBank phage genomes; AMR annotation; test-strain genomes | – |

List: `results/comparison/file_inventory.csv`. With the Source Data each panel could be compared point by point; with the Supplementary Text the ED7 and ED2 numbers could probably be resolved.

## Limits

Single seed per model run; UMAP layouts differ run to run; the public matrix lacks one strain; the recommender was not re-run; the genus PERMANOVA and the ED7 value are unexplained.
