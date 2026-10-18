# NADIA 0.99.2

Changes in response to the Bioconductor review (#186).

## Significant user-visible changes

* NADIA no longer sets the random seed; call `set.seed()` to reproduce results.
  The `seed` arguments are removed, and `pattern_profiler_analysis()` replaces
  `seeds` with `n_starts`.
* 61 building blocks of the metric and benchmark functions are no longer
  exported (from 111 to 51 exported functions); the orchestrators still return
  everything they produced.
* Protein score columns are `PG.Score` and `PG.Score.RunWise_*` in every format,
  instead of `PG.Cscore`; the report's metric is kept in the `nadia_labels`
  attribute. Older tables and `.nadia` files still load.
* `imp_method = "missForest"` now uses the proteins as observations, as
  `MsCoreUtils` does; `observations = "samples"` keeps NAguideR's orientation.
* `imp_method = "nbavg"` is removed, and `bpca` no longer accepts `nPcs`.

## New features

* `process_proteomics()` accepts a `SummarizedExperiment`, through the new
  reader `preprocess_se()`.
* `export_format` also applies to the exported matrices.

## Dependencies

* `limma`, `rrcovNA`, `pROC` and `MsCoreUtils` move to `Imports`;
  `highcharter`, `reactable` and `paletteer` move to `Suggests`; `readr` is
  dropped. Every suggested package is checked before use, with the call that
  installs it.
* Imputation methods `zero`, `min`, `with`, `MinDet`, `MinProb`, `QRILC`, `bpca`
  and `missForest`, and `vsn` normalization, are computed by `MsCoreUtils`,
  with identical results.

## Bug fixes

* Batch correction with protected covariates no longer discards measured
  values: proteins it cannot correct with the covariates are corrected for the
  batch only.
* Every plotting function resolves palettes the same way; an unknown palette is
  an error.
* Heatmaps accept partial row annotations, keep the condition colors in a
  targeted heatmap and no longer change ComplexHeatmap's options.
* Interactive tables: two tables on one page no longer interfere, the Excel
  export keeps text that starts with digits, and conditions with underscores
  are accepted.
* Exported tables quote fields that contain tabs, newlines or quotes.
* `.nadia` files restore the processing result exactly, except `limpa`'s
  precision weights.
* `im_prepare_se()` accepts `norm_method = "log2"`.
* A Parquet export without `arrow` fails before the analysis starts.

# NADIA 0.99.0

Initial Bioconductor submission.

* Import of protein-level reports from Spectronaut (`preprocess_spectronaut()`),
  DIA-NN (`preprocess_diann()`) and Proteome Discoverer, both TMT
  (`preprocess_tmt()`) and label-free (`preprocess_lfq()`).
* `process_proteomics()` runs normalization (13 methods), optional batch
  correction, imputation (20 methods, including the two-stage MAR/MNAR `combo`
  strategy) and differential abundance with limma or limpa.
* Each protein's missingness is reported next to its fold change: across the
  experiment, in the comparison and in each condition (`MissGlobal`,
  `MissComp`, `MissCND1`, `MissCND2`).
* Evaluation of normalization and imputation methods
  (`normalization_metrics()`, `imputation_metrics()`) and spike-in
  benchmarking with OpDEA metrics (`benchmarking_proteomics()`,
  `benchmarking_multiple()`).
* Soft clustering of abundance profiles with Mfuzz
  (`pattern_profiler_analysis()`).
* Interactive and static figures (volcano, boxplot, PCA, heatmap, cluster
  profiles) and interactive result tables.
* `write_nadia()` and `read_nadia()` store a complete analysis in a single
  DuckDB file.
* Example dataset `nadia_dia` and trimmed example reports in `inst/extdata/`.
