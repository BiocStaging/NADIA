# NADIA 0.99.2

Changes made in response to the Bioconductor review.

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
