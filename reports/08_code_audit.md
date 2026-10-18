# 8. Code audit and simplification

The scripts were audited for over-engineering (the `ponytail-audit` pass) and simplified without changing behaviour (the `code-simplification` pass). Correctness, security and performance were out of scope. The audit listed 11 findings, ranked by the size of the cut; 9 were applied and 2 were deliberately not applied.

## Applied

| # | Finding | Action |
|---|---|---|
| 1 | `run_lps.py` (60 lines) was superseded by `run_lps_markers.py` | Deleted. Its output, `results/replication/lps_agreement.csv`, stays as the record that K-12 *waaL* alone finds only K12, R1 and R4 |
| 2 | The Fig 6 statistics were computed twice (`fig6_cocktails.py` and `fig_main_3_6.py`) | Deleted `fig6_cocktails.py`; `fig_main_3_6.py` also feeds the comparison table |
| 3 | `fig4_cd.R` and `fig_ed6.R` repeated about 16 lines of set-up | Shared part moved to `lmm_common.R` |
| 4 | ED4A and ED5 drew the same heat-strip block | One `strips()` helper in `fig_extended.py` |
| 5 | `fig5_curves.py` redeclared the colours and axis styling that `figlib.py` provides | Imports `figlib` |
| 6 | Dead `Xgh`, `meta["GH"]` and an unused inner function in `fig_main_3_6.py` | Removed |
| 7 | The `srst2` environment was never run (SRST2 needs reads; `mlst` replaced it) | Removed from `setup_envs.sh`; `config/envs/srst2.yml` deleted |
| 9 | An unused colour map in Fig 2c | Removed |
| 10 | `status()` returned a string that a later `.replace` renamed | Returns the final string |

## Not applied

| # | Finding | Reason |
|---|---|---|
| 8 | The `macsy` environment (MacSyFinder 2.1.6) duplicates MacSyFinder 2.1.4 in the `defense` environment | The ABC-capsule results were produced with the `macsy` environment. Switching the script without re-running would break the link between code and result |
| 11 | `build_comparison.py` imports the whole `predict` module only for `published()` | Gain too small for the risk of touching the model code |

## Checks after the changes

There is no test suite, so behaviour was checked by output. The Fig 1–2 and tree scripts were re-run and their CSVs are byte-identical to the earlier outputs. The model runs were repeated with fixed seeds (see below) and the figure and comparison tables regenerated from them.

The PERMANOVA, NODF and BRIM functions of `fig_main_3_6.py` and `fig_extended.py` were moved into `matrix_stats.py` after the audit; the outputs of both scripts were checked to be identical before and after.

## Code added after the audit

`supp.py`, `matrix.py`, `matrix_stats.py`, `regen_matrix_analyses.py`, `compare_source_data.py`, `recommender.py`, `download_supplement.sh`, `run_notebooks.py`, `notebook_patches.py`, `phage_genomes.py`, `fig_phage_maps.py`, `si_analyses.py`, `si_rbp_lmm.R` and the extensions of `fig_ed6.R`, `fig_comparison.py` and `build_comparison.py` were written after the audit and were not audited.

## Not published on purpose

The paper's PDFs and text (copyright), the authors' data, the genome assemblies, and the multi-gigabyte intermediate files, all of which the scripts regenerate.
