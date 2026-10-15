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
* `proteomics_heatmap()` no longer changes ComplexHeatmap's `message` option
  for the rest of the session.

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
