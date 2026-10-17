# NADIA 0.99.2

Changes made in response to the Bioconductor review.

* NADIA no longer sets the random seed. The functions that use random numbers
  draw them from the session's generator, so call `set.seed()` beforehand to
  reproduce their results. The `seed` argument is removed from
  `normalization_metrics()`, `nm_compute_metrics()`, `nm_plot_metrics()`,
  `imputation_metrics()`, `im_compute_metrics()` and
  `pattern_profiler_analysis()`. The `PI`, `MLE` and `mice` imputation methods
  no longer use a fixed default seed; `MLE` and `mice` still accept one through
  `method_args`.
* `pattern_profiler_analysis()` replaces `seeds` with `n_starts`, the number of
  random starts used to evaluate each candidate number of clusters (default 3,
  as before).
* `limma`, `rrcovNA` and `pROC` move from `Suggests` to `Imports`: the default
  workflow and the benchmark metrics now run with a standard installation.
* `highcharter` and `reactable` move to `Suggests`, so the interactive figures
  and tables, and with them the Highcharts licence, are an explicit choice.
  The functions that need them say so when they are missing, and
  `benchmarking_proteomics()` skips its interactive volcano plots.
* `paletteer` moves to `Suggests`: it is needed only for palettes given as
  `"package::palette"`, and its absence is now always an error that says what
  to install. The p-value annotation of `proteomics_heatmap()` used to ignore
  the requested palette silently, and the Pattern Profiler charts with a
  warning; both now stop instead of drawing a palette nobody asked for.
* `readr` is no longer a dependency: files are read and written with base R.
  Exported tables now quote a field that contains a tab, a newline or a quote,
  which `write_tsv()` wrote unquoted, breaking the table.
* Every package in `Suggests` is now checked before it is used, with an error
  that names it and gives the `BiocManager::install()` call. This adds the
  checks that were missing (`ComplexHeatmap`, `tidyHeatmap` and `circlize` in
  the heatmap functions, `arrow` in the Pattern Profiler, `ggsci` for its
  palettes) and gives the existing ones a single format.
* A Parquet export that needs `arrow` now fails before the analysis starts.
  `process_proteomics()` used to run to the end and then skip the Parquet
  files with a warning.
* Batch correction with protected covariates no longer discards measured
  values. BERT needs two values of a protein in every combination of batch and
  covariate level, and set the protein's values in a batch to missing when they
  were lacking: those values were then imputed, and proteins left without
  values in a condition dropped out of the analysis. With two replicates per
  condition in each batch this affected over a quarter of the proteins. Those
  proteins are now corrected for the batch only, and are listed in
  `metadata(se)$batch_correction`.
* Every plotting function resolves palettes in the same way: `"pkg::name"`
  (paletteer), `"brewer:Name"` or the bare RColorBrewer name, or a vector of
  colours, named after the groups or not. A palette with fewer colours than
  groups is interpolated instead of recycled, and an unknown palette is an
  error instead of a silent fallback to the default. `volcano_highchart_list()`
  now accepts RColorBrewer palettes, and the Pattern Profiler charts accept a
  vector of colours, which used to fail.
* `proteomics_heatmap()` no longer changes ComplexHeatmap's `message` option
  for the rest of the session.
* The four `*_widget()` tables share one script, shipped in `inst/js/`,
  instead of embedding their own JavaScript. Two tables of the same kind on
  one page no longer act on each other: a condition chip or the filter button
  of one table used to change the other. The label above the condition chips
  gives the number of columns a condition hides, which was fixed at 16.

# NADIA 0.99.0

Initial Bioconductor submission.

* Import of protein-level reports from Spectronaut (`preprocess_spectronaut()`),
  DIA-NN (`preprocess_diann()`) and Proteome Discoverer, both TMT
  (`preprocess_tmt()`) and label-free (`preprocess_lfq()`).
* `process_proteomics()` runs normalisation (13 methods), optional batch
  correction, imputation (20 methods, including the two-stage MAR/MNAR `combo`
  strategy) and differential abundance with limma or limpa.
* Each protein's missingness is reported next to its fold change: across the
  experiment, in the comparison and in each condition (`MissGlobal`,
  `MissComp`, `MissCND1`, `MissCND2`).
* Evaluation of normalisation and imputation methods
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
