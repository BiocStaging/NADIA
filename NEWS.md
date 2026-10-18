# NADIA 0.99.2

Changes made in response to the Bioconductor review (#186).

## Significant user-visible changes

* NADIA no longer sets the random seed; call `set.seed()` to reproduce results.
  The `seed` arguments have been removed, and `pattern_profiler_analysis()` replaces
  `seeds` with `n_starts`.
* 61 helper functions used by the metric and benchmarking functions are no longer
  exported (reducing the number of exported functions from 111 to 50; 51 with
  the new `preprocess_se()`). The
  high-level functions continue to return all previously available outputs.
* Protein score columns are `PG.Score` and `PG.Score.RunWise_*` in every format,
  instead of `PG.Cscore`; the original report metric is preserved in the
  `nadia_labels` attribute. Older tables and `.nadia` files can still be loaded.
* `imp_method = "missForest"` now treats proteins as observations, as
  `MsCoreUtils` does; `observations = "samples"` retains the orientation used by NAguideR.
* `imp_method = "nbavg"` has been removed, and `bpca` no longer accepts `nPcs`.

## New features

* `process_proteomics()` accepts a `SummarizedExperiment` through the new
  import function `preprocess_se()`.
* `export_format` also applies to the exported matrices.

## Dependencies

* `limma`, `rrcovNA`, `pROC` and `MsCoreUtils` move to `Imports`;
  `highcharter`, `reactable` and `paletteer` move to `Suggests`; `readr` is
  removed. The availability of each suggested package is checked before use,
  and an installation command is provided if it is missing.
* Imputation methods `zero`, `min`, `with`, `MinDet`, `MinProb`, `QRILC`, `bpca`
  and `missForest`, along with `vsn` normalization, now use `MsCoreUtils`
  and produce the same results as before.

## Bug fixes

* Batch correction with protected covariates no longer discards measured
  values: proteins that cannot be corrected while protecting the covariates
  receive batch-only correction.
* All plotting functions resolve palettes consistently; an unknown palette
  name raises an error.
* Heatmaps accept partial row annotations, preserve condition colors in
  comparison-specific heatmaps and no longer modify ComplexHeatmap options.
* Interactive tables on the same page no longer interfere with each other.
  Excel export preserves text that starts with digits, and condition names
  containing underscores are accepted.
* Exported tables quote fields that contain tabs, newlines or quotes.
* `.nadia` files restore the processing result exactly, except for the
  precision weights produced by `limpa`.
* `im_prepare_se()` accepts `norm_method = "log2"`.
* A Parquet export without `arrow` fails before the analysis starts.

# NADIA 0.99.0

Initial Bioconductor submission.

* Import of protein-level reports from Spectronaut (`preprocess_spectronaut()`),
  DIA-NN (`preprocess_diann()`) and Proteome Discoverer, including both TMT
  (`preprocess_tmt()`) and label-free (`preprocess_lfq()`).
* `process_proteomics()` runs normalization (13 methods), optional batch
  correction, imputation (20 methods, including the two-stage MAR/MNAR `combo`
  strategy) and differential abundance analysis with limma or limpa.
* Each protein's missingness is reported alongside its fold change: across the
  experiment, within the comparison and within each condition (`MissGlobal`,
  `MissComp`, `MissCND1`, `MissCND2`).
* Evaluation of normalization and imputation methods
  (`normalization_metrics()`, `imputation_metrics()`) and spike-in
  benchmarking with OpDEA metrics (`benchmarking_proteomics()`,
  `benchmarking_multiple()`).
* Soft clustering of abundance profiles with Mfuzz
  (`pattern_profiler_analysis()`).
* Interactive and static figures (volcano plots, boxplots, PCA plots, heatmaps
  and cluster profiles) and interactive result tables.
* `write_nadia()` and `read_nadia()` save and load a complete analysis in a
  single DuckDB file.
* Example dataset `nadia_dia` and trimmed example reports in `inst/extdata/`.
