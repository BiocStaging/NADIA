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
* `proteomics_heatmap()` accepts a row annotation that does not mention every
  protein shown: those proteins are annotated "NA", drawn white and placed last,
  as empty values already were. The heatmap used to fail.
* In a targeted heatmap the conditions keep their colours. tidyHeatmap gives
  the colours to the condition levels by position, and the level of the
  condition left out of the comparison was kept, so in a comparison such as
  `D-A` the D samples were drawn without colour.
* The building blocks of the metric and benchmark orchestrators are no
  longer exported: 61 functions, which takes the exported interface from 111
  to 50. These are the `nm_plot_*()`, `nm_rank_*()`, `nm_compute_metrics()`,
  `im_plot_*()`, `im_compute_metrics()`, `bm_*()` and `import_*_results()`
  functions, the `benchmark_*_gg()`, `benchmark_volcano_hc*()` and
  `compute_*_metrics()` functions of the single benchmark, `pvca_compute()`,
  `pvca_plot()` and `pca_highchart()`. Everything they produce is still
  returned by `normalization_metrics()`, `imputation_metrics()`,
  `benchmarking_proteomics()`, `benchmarking_multiple()`, `pvca_analysis()` and
  `pca_highchart_list()` (the first two and `benchmarking_multiple()` choose
  which plots to draw with `plots`), and the functions themselves remain
  reachable as `NADIA:::name()`.
  `import_norm_matrices()` and `import_imp_matrices()` stay exported: they are
  how matrices from other pipelines enter the metrics.
* `imp_method = "missForest"` now fits the random forests with the proteins as
  observations and the samples as variables, the default of
  `MsCoreUtils::impute_RF()`. It used to transpose the matrix, as NAguideR
  does, so that each forest learnt from at most n - 1 samples; with the few
  replicates of a differential-abundance design that recovered hidden values
  clearly worse. On the example data the error on hidden values falls from
  0.203 to 0.110 (NRMSE, `imputation_metrics()`), and the spike-in MCC of the
  whole pipeline rises from 0.58 to 0.67. NAguideR's orientation remains
  available as `method_args = list(missForest = list(observations = "samples"))`.
* The methods that give identical results in `MsCoreUtils` are now computed
  by it: imputation `zero`, `min`, `with`, `MinDet`, `MinProb`, `QRILC`,
  `bpca` and `missForest`, and normalisation `vsn`. `MsCoreUtils` moves into
  `Imports`. Every result is unchanged (checked on the four example datasets).
  `bpca` no longer accepts `nPcs`, which `MsCoreUtils` fixes at the number of
  samples minus one, NADIA's default. Where NADIA differs from `MsCoreUtils`
  (`knn`, `MLE`, `quantile.robust`) it keeps its own implementation; see the
  Bioconductor review (#186) for the reasons.
* `imp_method = "nbavg"` is removed. Neighbour averaging comes from MSnbase,
  where it is meant for fractions ordered along a separation gradient and is
  invalid otherwise; NADIA's variant averaged the nearest columns, which
  mixes replicates of different conditions at the edges of each block.
* `process_proteomics()` accepts a `SummarizedExperiment`, so data imported
  with other tools (a protein-level `QFeatures` assay, a `limpa` import, or a
  matrix of intensities) enter NADIA's pipeline. The condition is read from
  `colData()$Condition`; `log2_input = TRUE` declares intensities already on
  the log2 scale and `assay_name` chooses the assay. The conversion is
  available on its own as `preprocess_se()`, a fifth reader returning the same
  `proteomics_data` object as the other four, which `write_nadia()` stores like
  any other. The example data rebuilt as a `SummarizedExperiment` give results
  identical to those from the Spectronaut report.
* The protein score columns are named `PG.Score` and `PG.Score.RunWise_*` in
  every format, instead of Spectronaut's `PG.Cscore`, under which the TMT and
  LFQ readers stored Proteome Discoverer's `Sum PEP Score` and Mascot score.
  Each table records which metric of the report it holds (attribute
  `nadia_labels`), and the interactive tables use it as the column label
  instead of "Spectronaut Cscore". Tables exported and `.nadia` files written
  before keep loading. The input-formats vignette lists, for every format,
  which metric of the report fills each identification column.
* The Excel export of `protein_list_widget()` and `quant_list_widget()` no
  longer turns text that starts with digits into a number: a description such
  as "40S ribosomal protein S26" was written as 40 (114 of the 2,000
  descriptions of the example data). Identifiers, descriptions and genes are
  written as text, and other columns become numbers only when the whole value
  is one.
* A `.nadia` file now restores the processing result exactly. The metadata
  of the rebuilt `SummarizedExperiment` used to keep only the condition and
  label columns, losing the MAR/MNAR masks of the combo imputation and the
  batch-correction diagnostics (the proteins corrected for the batch only or
  left unadjusted); the comparisons came back as text instead of a factor, and
  the parameters left at `NULL` were dropped. `nadia_result()` now returns an
  object `identical()` to the one `process_proteomics()` returned, except, with
  `limpa`, for the precision weights it keeps for its own analysis. Files
  written before still read, without those components.
* `im_prepare_se()` accepts `norm_method = "log2"`, given directly or as the
  winner of `normalization_metrics()$pc1_rank`; it used to stop with "Unknown
  normalization method". The result then has the single `"log2"` assay.
* `export_format` applies to the normalised, batch-corrected and imputed
  matrices written by `process_proteomics()`, which were always TSV; with
  `"both"` or `"parquet"` they now come in the format asked for, like the plot
  inputs. TSV exports are unchanged.

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
