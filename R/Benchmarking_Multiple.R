# =============================================================================
# Benchmarking_Multiple.R
# Multiple-method OpDEA benchmarking: import, rank, and visualize
#
# Compares opdea_metrics from multiple normalization/imputation combinations.
# Ranking follows OpDEA methodology (Peng et al., Nature Comms 2024):
#   rank_final = mean(rank_nMCC, rank_G_mean, rank_pAUC_001, rank_pAUC_005, rank_pAUC_010)
#
# License: GPL-3 | Copyright (c) 2025 Sergio Ciordia
# =============================================================================

# =============================================================================
# SECTION 1: INTERNAL HELPERS
# =============================================================================

#' Required columns in the combined opdea data.frame
#' @keywords internal
.BM_REQUIRED_COLS <- c("Assay", "Comparison", "nMCC", "G_mean",
                       "pAUC_001", "pAUC_005", "pAUC_010")

#' Required columns in the combined benchmark metrics data.frame
#' @keywords internal
.BM_METRICS_REQUIRED_COLS <- c("Assay", "Comparison", "Sensitivity", "Specificity",
                               "Precision", "NPV", "F1", "Accuracy", "MCC")

#' Required columns in the combined confusion data.frame
#' @keywords internal
.BM_CONFUSION_COLS <- c("Assay", "Comparison", "TP", "FP", "TN", "FN")

#' Minimum required columns in the combined classified data.frame
#' @keywords internal
.BM_CLASSIFIED_COLS <- c("Assay", "Comparison", "truth")

#' Default metrics for ranking
#' @keywords internal
.BM_METRICS <- c("nMCC", "G_mean", "pAUC_001", "pAUC_005", "pAUC_010")

#' Extended metrics for combined ranking (all higher = better)
#' @keywords internal
.BM_EXTENDED_METRICS <- c("Sensitivity", "Specificity", "Precision", "NPV",
                          "F1", "pAUC_005", "Accuracy", "MCC",
                          "Performance", "nMCC", "G_mean")


#' Read a TSV file
#' @param path File path
#' @return data.frame
#' @keywords internal
.bm_read_tsv <- function(path) {
  .read_delim_file(path)
}

#' Row-bind data.frames with heterogeneous columns
#'
#' Combines a list of data.frames by the union of their columns, filling
#' missing columns with NA. Unlike \code{do.call(rbind, ...)}, this tolerates
#' methods whose exported tables carry different column sets (e.g. an external
#' tool like Proteome Discoverer whose classified table keeps its own de_res
#' columns).
#' @keywords internal
.bm_rbind_fill <- function(df_list) {
  df_list <- Filter(function(d) !is.null(d) && nrow(d) > 0, df_list)
  if (length(df_list) == 0) return(NULL)
  all_cols <- unique(unlist(lapply(df_list, colnames)))
  df_list <- lapply(df_list, function(d) {
    for (m in setdiff(all_cols, colnames(d))) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, df_list)
}


#' Validate combined opdea data.frame
#' @param df data.frame to validate
#' @keywords internal
.bm_validate_opdea <- function(df) {
  if (!is.data.frame(df))
    stop("opdea_combined must be a data.frame.")

  missing <- setdiff(.BM_REQUIRED_COLS, colnames(df))
  if (length(missing) > 0)
    stop("Missing required columns in opdea_combined: ",
         paste(missing, collapse = ", "),
         "\nRequired: ", paste(.BM_REQUIRED_COLS, collapse = ", "))

  if (nrow(df) == 0)
    stop("opdea_combined has 0 rows.")

  invisible(TRUE)
}


#' Validate combined confusion data.frame
#' @param df data.frame to validate
#' @keywords internal
.bm_validate_confusion <- function(df) {
  if (!is.data.frame(df))
    stop("confusion_combined must be a data.frame.")

  missing <- setdiff(.BM_CONFUSION_COLS, colnames(df))
  if (length(missing) > 0)
    stop("Missing required columns in confusion_combined: ",
         paste(missing, collapse = ", "),
         "\nRequired: ", paste(.BM_CONFUSION_COLS, collapse = ", "))

  if (nrow(df) == 0)
    stop("confusion_combined has 0 rows.")

  invisible(TRUE)
}


#' Validate combined classified data.frame
#' @param df data.frame to validate
#' @param p_col Name of the p-value column
#' @keywords internal
.bm_validate_classified <- function(df, p_col = "adj.P.Val") {
  if (!is.data.frame(df))
    stop("classified_combined must be a data.frame.")

  missing <- setdiff(.BM_CLASSIFIED_COLS, colnames(df))
  if (length(missing) > 0)
    stop("Missing required columns in classified_combined: ",
         paste(missing, collapse = ", "),
         "\nRequired: ", paste(.BM_CLASSIFIED_COLS, collapse = ", "))

  if (!p_col %in% colnames(df))
    stop("P-value column '", p_col, "' not found in classified_combined.")

  if (nrow(df) == 0)
    stop("classified_combined has 0 rows.")

  invisible(TRUE)
}


#' Validate combined benchmark metrics data.frame
#' @param df data.frame to validate
#' @keywords internal
.bm_validate_bench_metrics <- function(df) {
  if (!is.data.frame(df))
    stop("bench_metrics_combined must be a data.frame.")

  missing <- setdiff(.BM_METRICS_REQUIRED_COLS, colnames(df))
  if (length(missing) > 0)
    stop("Missing required columns in bench_metrics_combined: ",
         paste(missing, collapse = ", "),
         "\nRequired: ", paste(.BM_METRICS_REQUIRED_COLS, collapse = ", "))

  if (nrow(df) == 0)
    stop("bench_metrics_combined has 0 rows.")

  invisible(TRUE)
}


#' Convert Performance (categorical) to numeric score
#'
#' Excellent=5, Very Good=4, Good=3, Acceptable=2, Poor=1, NA=0
#' @param x Character vector of Performance labels
#' @return Numeric vector
#' @keywords internal
.bm_performance_to_numeric <- function(x) {
  map <- c("Excellent" = 5, "Very Good" = 4, "Good" = 3,
           "Acceptable" = 2, "Poor" = 1)
  vals <- map[x]
  vals[is.na(vals)] <- 0
  unname(vals)
}


#' Aggregate metrics by Assay
#'
#' @param df Combined opdea data.frame (with Assay column)
#' @param metrics Character vector of metric column names
#' @param agg_fun Aggregation function (mean or median)
#' @return data.frame with one row per Assay
#' @keywords internal
.bm_aggregate_metrics <- function(df, metrics, agg_fun = mean) {
  assays <- unique(df$Assay)
  agg_list <- lapply(assays, function(a) {
    sub <- df[df$Assay == a, , drop = FALSE]
    vals <- vapply(metrics, function(m) {
      agg_fun(sub[[m]], na.rm = TRUE)
    }, numeric(1))
    row <- data.frame(Assay = a, stringsAsFactors = FALSE)
    for (i in seq_along(metrics)) row[[metrics[i]]] <- round(vals[i], 4)
    row
  })
  do.call(rbind, agg_list)
}


#' Rank methods by descending metric value
#'
#' @param agg_df Aggregated data.frame (one row per Assay)
#' @param metrics Character vector of metric column names
#' @return data.frame with rank columns and rank_final, sorted by rank_final
#' @keywords internal
.bm_rank_methods <- function(agg_df, metrics) {
  rank_df <- data.frame(Assay = agg_df$Assay, stringsAsFactors = FALSE)

  for (m in metrics) {
    vals <- agg_df[[m]]
    vals[is.na(vals)] <- -Inf
    rank_df[[paste0("rank_", m)]] <- rank(-vals, ties.method = "average")
  }

  rank_cols <- paste0("rank_", metrics)
  rank_df$rank_final <- round(rowMeans(rank_df[, rank_cols, drop = FALSE]), 2)
  rank_df <- rank_df[order(rank_df$rank_final), , drop = FALSE]
  rownames(rank_df) <- NULL
  rank_df
}


#' Export data.frame to TSV
#' @param data data.frame
#' @param filepath Output path
#' @keywords internal
.bm_export_data <- function(data, filepath) {
  .write_tsv(data, filepath)
  invisible(filepath)
}


#' Export ggplot2 to PNG
#' @param gg ggplot2 object
#' @param filepath Output path
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution
#' @keywords internal
.bm_export_gg_plot <- function(gg, filepath, width = 12, height = 8, dpi = 150) {
  if (is.null(gg)) return(invisible(NULL))

  tryCatch({
    ggplot2::ggsave(filepath, plot = gg, width = width, height = height,
                    dpi = dpi, bg = "white")
  }, error = function(e) {
    warning("Could not export the plot to ", filepath, ": ", e$message)
  })

  invisible(filepath)
}


# =============================================================================
# SECTION 2: IMPORT FUNCTION
# =============================================================================

#' Import one kind of benchmark table from every pipeline folder
#'
#' Shared by the four `import_*()` functions below, which only differ in the
#' file they look for and in how they validate the result.
#'
#' @param results_dir,pattern,method_names,recursive See `import_opdea_results()`.
#' @param validate Function that validates the combined table, or `NULL`.
#' @param label Name of the calling function, for the message.
#' @return Combined data.frame with an `Assay` column naming the pipeline.
#' @keywords internal
#' @noRd
.bm_import <- function(results_dir, pattern, method_names, recursive,
                       validate, label) {
  if (!dir.exists(results_dir))
    stop("Directory not found: ", results_dir)

  files <- list.files(results_dir, pattern = pattern,
                      full.names = TRUE, recursive = recursive)
  if (length(files) == 0)
    stop("No files matching pattern '", pattern, "' found in: ", results_dir)

  # --- Determine Assay names ---
  if (is.null(method_names)) {
    method_names <- basename(dirname(files))
  }
  if (length(method_names) != length(files))
    stop("Length of method_names (", length(method_names),
         ") must match number of files (", length(files), ").")

  # --- Read and combine ---
  df_list <- vector("list", length(files))
  for (i in seq_along(files)) {
    df <- .bm_read_tsv(files[i])
    df$Assay <- method_names[i]
    df_list[[i]] <- df
  }

  combined <- .bm_rbind_fill(df_list)
  rownames(combined) <- NULL

  if (!is.null(validate)) validate(combined)

  n_methods <- length(unique(combined$Assay))
  n_comps   <- length(unique(combined$Comparison))
  message(label, ": loaded ", n_methods, " method(s), ",
          n_comps, " comparison(s), ", nrow(combined), " rows.")

  combined
}


#' Import OpDEA metrics from multiple benchmark result folders
#'
#' Reads \code{benchmark_opdea_metrics.tsv} files from subdirectories,
#' adding an \code{Assay} column derived from the folder name.
#'
#' @param results_dir Parent directory containing benchmark subfolders
#' @param pattern Regex pattern for the opdea metrics filename
#'   (default: \code{"benchmark_opdea_metrics\\.tsv$"})
#' @param method_names Optional character vector of Assay names.
#'   If NULL, names are derived from the parent folder of each file.
#' @param recursive Logical. Search subdirectories recursively? (default: TRUE)
#'
#' @return data.frame in long format with columns:
#'   Assay, Comparison, nMCC, G_mean, pAUC_001, pAUC_005, pAUC_010
#' @keywords internal
#' @noRd
import_opdea_results <- function(results_dir,
                                 pattern      = "benchmark_opdea_metrics\\.tsv$",
                                 method_names = NULL,
                                 recursive    = TRUE) {
  .bm_import(results_dir, pattern, method_names, recursive,
             validate = .bm_validate_opdea, label = "import_opdea_results")
}



#' Import confusion matrices from multiple benchmark result folders
#'
#' Reads \code{benchmark_confusion_overall.tsv} files from subdirectories,
#' adding an \code{Assay} column derived from the folder name.
#'
#' @param results_dir Parent directory containing benchmark subfolders
#' @param pattern Regex pattern for the confusion filename
#' @param method_names Optional character vector of Assay names
#' @param recursive Search subdirectories recursively? (default: TRUE)
#'
#' @return data.frame with columns:
#'   Assay, Comparison, N, TP, FP, TN, FN, plus percentage columns
#' @keywords internal
#' @noRd
import_confusion_results <- function(results_dir,
                                     pattern      = "benchmark_confusion_overall\\.tsv$",
                                     method_names = NULL,
                                     recursive    = TRUE) {
  .bm_import(results_dir, pattern, method_names, recursive,
             validate = .bm_validate_confusion, label = "import_confusion_results")
}


#' Import classified results from multiple benchmark result folders
#'
#' Reads \code{benchmark_classified.tsv} files from subdirectories,
#' adding an \code{Assay} column derived from the folder name.
#'
#' @param results_dir Parent directory containing benchmark subfolders
#' @param pattern Regex pattern for the classified filename
#' @param method_names Optional character vector of Assay names
#' @param recursive Search subdirectories recursively? (default: TRUE)
#'
#' @return data.frame with columns from classified_df plus Assay
#' @keywords internal
#' @noRd
import_classified_results <- function(results_dir,
                                      pattern      = "benchmark_classified\\.tsv$",
                                      method_names = NULL,
                                      recursive    = TRUE) {
  .bm_import(results_dir, pattern, method_names, recursive,
             validate = NULL, label = "import_classified_results")
}


#' Import benchmark metrics from multiple benchmark result folders
#'
#' Reads \code{benchmark_metrics.tsv} files from subdirectories,
#' adding an \code{Assay} column derived from the folder name.
#'
#' @param results_dir Parent directory containing benchmark subfolders
#' @param pattern Regex pattern for the metrics filename
#' @param method_names Optional character vector of Assay names
#' @param recursive Search subdirectories recursively? (default: TRUE)
#'
#' @return data.frame with columns:
#'   Assay, Comparison, TP, FP, TN, FN, Sensitivity, Specificity,
#'   Precision, NPV, Accuracy, F1, MCC, AUC, Performance
#' @keywords internal
#' @noRd
import_benchmark_metrics <- function(results_dir,
                                     pattern      = "benchmark_metrics\\.tsv$",
                                     method_names = NULL,
                                     recursive    = TRUE) {
  .bm_import(results_dir, pattern, method_names, recursive,
             validate = .bm_validate_bench_metrics, label = "import_benchmark_metrics")
}


# =============================================================================
# SECTION 3: RANKING COMPUTATION
# =============================================================================

#' Compute OpDEA ranking across multiple methods
#'
#' For each metric, aggregates across comparisons (mean and median),
#' then ranks methods in descending order (rank 1 = best).
#' Final rank = average of individual metric ranks.
#'
#' @param opdea_combined data.frame with columns: Assay, Comparison,
#'   nMCC, G_mean, pAUC_001, pAUC_005, pAUC_010
#' @param metrics Character vector of metrics to include in ranking
#'   (default: all 5 OpDEA metrics)
#'
#' @return Named list with:
#'   \describe{
#'     \item{mean_aggregated}{data.frame of mean metric values per Assay}
#'     \item{median_aggregated}{data.frame of median metric values per Assay}
#'     \item{mean_ranking}{data.frame of ranks based on mean aggregation}
#'     \item{median_ranking}{data.frame of ranks based on median aggregation}
#'     \item{opdea_combined}{Input combined data.frame}
#'     \item{n_methods}{Number of methods}
#'     \item{n_comparisons}{Number of comparisons}
#'     \item{metrics_used}{Metrics included in ranking}
#'   }
#' @keywords internal
#' @noRd
bm_compute_ranking <- function(opdea_combined,
                               metrics = .BM_METRICS) {
  .bm_validate_opdea(opdea_combined)

  # Validate requested metrics exist
  avail <- intersect(metrics, colnames(opdea_combined))
  if (length(avail) == 0)
    stop("None of the requested metrics found in opdea_combined.")
  if (length(avail) < length(metrics)) {
    missing <- setdiff(metrics, avail)
    warning("Metrics not found (skipped): ", paste(missing, collapse = ", "))
    metrics <- avail
  }

  # Aggregate
  mean_agg   <- .bm_aggregate_metrics(opdea_combined, metrics, mean)
  median_agg <- .bm_aggregate_metrics(opdea_combined, metrics, median)

  # Rank
  mean_rank   <- .bm_rank_methods(mean_agg, metrics)
  median_rank <- .bm_rank_methods(median_agg, metrics)

  list(
    mean_aggregated   = mean_agg,
    median_aggregated = median_agg,
    mean_ranking      = mean_rank,
    median_ranking    = median_rank,
    opdea_combined    = opdea_combined,
    n_methods         = length(unique(opdea_combined$Assay)),
    n_comparisons     = length(unique(opdea_combined$Comparison)),
    metrics_used      = metrics
  )
}


#' Compute OpDEA ranking separately for each comparison
#'
#' For each comparison, ranks methods directly on their metric values
#' (no aggregation needed since there is one value per method per comparison).
#'
#' @param opdea_combined data.frame with columns: Assay, Comparison,
#'   nMCC, G_mean, pAUC_001, pAUC_005, pAUC_010
#' @param metrics Character vector of metrics to include in ranking
#'
#' @return Named list with:
#'   \describe{
#'     \item{by_comparison}{Named list where each key is a comparison and
#'       each value is a data.frame with metric values, ranks, and rank_final}
#'     \item{ranking_combined}{data.frame with all comparisons in long format
#'       (includes Comparison column)}
#'   }
#' @keywords internal
#' @noRd
bm_compute_ranking_by_comparison <- function(opdea_combined,
                                             metrics = .BM_METRICS) {
  .bm_validate_opdea(opdea_combined)

  avail <- intersect(metrics, colnames(opdea_combined))
  if (length(avail) == 0)
    stop("None of the requested metrics found in opdea_combined.")
  metrics <- avail

  comps <- sort(unique(opdea_combined$Comparison))
  by_comp <- vector("list", length(comps))
  names(by_comp) <- comps

  combined_list <- vector("list", length(comps))

  for (i in seq_along(comps)) {
    comp <- comps[i]
    sub <- opdea_combined[opdea_combined$Comparison == comp,
                          c("Assay", metrics), drop = FALSE]

    # Aggregate by Assay to ensure 1 row per method (dedup)
    sub <- .bm_aggregate_metrics(sub, metrics, mean)

    # Rank directly (1 row per method)
    rank_df <- .bm_rank_methods(sub, metrics)

    # Merge values + ranks
    merged <- merge(sub, rank_df, by = "Assay", sort = FALSE)
    merged <- merged[order(merged$rank_final), , drop = FALSE]
    rownames(merged) <- NULL

    by_comp[[comp]] <- merged

    combined_entry <- merged
    combined_entry$Comparison <- comp
    combined_list[[i]] <- combined_entry
  }

  ranking_combined <- do.call(rbind, combined_list)
  rownames(ranking_combined) <- NULL

  list(
    by_comparison    = by_comp,
    ranking_combined = ranking_combined
  )
}


# =============================================================================
# SECTION 3b: EXTENDED RANKING (11-metric combined ranking)
# =============================================================================

#' Build an extended combined data.frame from opdea + benchmark metrics
#'
#' Merges OpDEA metrics and benchmark classification metrics by Assay +
#' Comparison, converts Performance to numeric, and selects the 11 extended
#' metrics for ranking.
#'
#' @param opdea_combined data.frame from import_opdea_results()
#' @param bench_metrics_combined data.frame from import_benchmark_metrics()
#' @return data.frame with columns: Assay, Comparison, and the 11 metrics
#'   in .BM_EXTENDED_METRICS
#' @keywords internal
.bm_build_extended <- function(opdea_combined, bench_metrics_combined) {
  # Compute Performance from F1 if column is missing
  if (!"Performance" %in% colnames(bench_metrics_combined) &&
      "F1" %in% colnames(bench_metrics_combined)) {
    bench_metrics_combined$Performance <- ifelse(
      is.na(bench_metrics_combined$F1), NA_character_,
      ifelse(bench_metrics_combined$F1 >= 0.9, "Excellent",
      ifelse(bench_metrics_combined$F1 >= 0.8, "Very Good",
      ifelse(bench_metrics_combined$F1 >= 0.7, "Good",
      ifelse(bench_metrics_combined$F1 >= 0.5, "Acceptable", "Poor")))))
  }

  # Merge by Assay + Comparison
  bench_cols <- intersect(
    c("Assay", "Comparison", "Sensitivity", "Specificity",
      "Precision", "NPV", "F1", "Accuracy", "MCC", "Performance"),
    colnames(bench_metrics_combined)
  )
  merged <- merge(
    bench_metrics_combined[, bench_cols, drop = FALSE],
    opdea_combined[, c("Assay", "Comparison", "nMCC", "G_mean",
                       "pAUC_005"), drop = FALSE],
    by = c("Assay", "Comparison"), all = FALSE
  )

  # Convert Performance to numeric
  if ("Performance" %in% colnames(merged)) {
    merged$Performance <- .bm_performance_to_numeric(merged$Performance)
  }

  merged
}


#' Compute extended ranking across 11 classification + OpDEA metrics
#'
#' For each of the 11 metrics, averages across comparisons per Assay (mean),
#' then ranks methods in descending order (higher = better for all metrics).
#' Final rank = mean of the 11 individual ranks.
#'
#' Note: this extended ranking can drift away from the canonical OpDEA ranking
#' (`bm_compute_ranking`, based on nMCC/G_mean/pAUC). `Performance` is collinear
#' with F1 (so F1 effectively counts twice), and in spike-in designs
#' `Accuracy`/`NPV` are dominated by the background (TN >> ), hence they have
#' little discriminative power. Use it as a complementary view, not as a
#' substitute for the OpDEA ranking.
#'
#' @param opdea_combined data.frame from import_opdea_results()
#' @param bench_metrics_combined data.frame from import_benchmark_metrics()
#' @param metrics Character vector of metrics to include (default:
#'   all 11 extended metrics)
#'
#' @return Named list with:
#'   \describe{
#'     \item{extended_combined}{Merged data.frame (long format)}
#'     \item{mean_aggregated}{data.frame of mean metric values per Assay}
#'     \item{mean_ranking}{data.frame with rank columns and rank_final}
#'     \item{n_methods}{Number of unique methods}
#'     \item{n_comparisons}{Number of unique comparisons}
#'     \item{metrics_used}{Metrics actually used}
#'   }
#' @keywords internal
#' @noRd
bm_compute_extended_ranking <- function(opdea_combined,
                                        bench_metrics_combined,
                                        metrics = .BM_EXTENDED_METRICS) {
  .bm_validate_opdea(opdea_combined)
  .bm_validate_bench_metrics(bench_metrics_combined)

  extended <- .bm_build_extended(opdea_combined, bench_metrics_combined)

  # Check requested metrics exist
  available <- intersect(metrics, colnames(extended))
  if (length(available) == 0)
    stop("None of the requested metrics found in extended data.")
  if (length(available) < length(metrics)) {
    missing <- setdiff(metrics, available)
    warning("Metrics not found (skipping): ", paste(missing, collapse = ", "))
  }
  metrics <- available

  mean_agg  <- .bm_aggregate_metrics(extended, metrics, mean)
  mean_rank <- .bm_rank_methods(mean_agg, metrics)

  list(
    extended_combined = extended,
    mean_aggregated   = mean_agg,
    mean_ranking      = mean_rank,
    n_methods         = length(unique(extended$Assay)),
    n_comparisons     = length(unique(extended$Comparison)),
    metrics_used      = metrics
  )
}


#' Horizontal bar chart of extended combined ranking
#'
#' @param extended_ranking List from bm_compute_extended_ranking()
#' @param title Optional plot title
#' @return ggplot object
#' @keywords internal
#' @noRd
bm_plot_extended_ranking_bars <- function(extended_ranking, title = NULL) {
  rank_df   <- extended_ranking$mean_ranking
  n_methods <- extended_ranking$n_methods

  rank_df$Assay <- factor(rank_df$Assay, levels = rev(rank_df$Assay))

  title <- title %||% "Extended Combined Ranking (11 metrics)"
  subtitle <- paste0(extended_ranking$n_methods, " methods, ",
                     extended_ranking$n_comparisons, " comparisons")

  ggplot2::ggplot(rank_df,
                  ggplot2::aes(x = Assay, y = rank_final, fill = rank_final)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = round(rank_final, 2)),
                       hjust = -0.2, size = 3.5) +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = NULL, y = "Final Rank (lower = better)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
}


#' Heatmap of extended ranking (11 metrics)
#'
#' @param extended_ranking List from bm_compute_extended_ranking()
#' @param title Optional plot title
#' @return ggplot object
#' @keywords internal
#' @noRd
bm_plot_extended_ranking_heatmap <- function(extended_ranking, title = NULL) {
  rank_df <- extended_ranking$mean_ranking
  metrics <- extended_ranking$metrics_used
  n_methods <- extended_ranking$n_methods

  # Order methods by rank_final (best first)
  method_order <- rank_df$Assay
  rank_cols <- c(paste0("rank_", metrics), "rank_final")

  # Reshape to long format
  long_list <- lapply(rank_cols, function(rc) {
    metric_label <- sub("^rank_", "", rc)
    data.frame(
      Assay  = rank_df$Assay,
      Metric = metric_label,
      Rank   = rank_df[[rc]],
      stringsAsFactors = FALSE
    )
  })
  long_df <- do.call(rbind, long_list)

  # Factor levels
  metric_order <- c(metrics, "final")
  long_df$Metric <- factor(long_df$Metric, levels = metric_order)
  long_df$Assay  <- factor(long_df$Assay, levels = rev(method_order))

  title <- title %||% "Extended Ranking Heatmap (11 metrics)"

  ggplot2::ggplot(long_df,
                  ggplot2::aes(x = Metric, y = Assay, fill = Rank)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(Rank, 1)),
                       size = 3, color = "black") +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid       = ggplot2::element_blank()
    )
}


# =============================================================================
# SECTION 4: VISUALIZATIONS
# =============================================================================

#' Heatmap of metric ranks by method
#'
#' Rows = Assay (ordered by rank_final), columns = metrics.
#' Fill color = rank value (green = best, red = worst).
#'
#' @param ranking_result List returned by \code{bm_compute_ranking()}
#' @param type "mean" or "median" aggregation (default: "mean")
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @keywords internal
#' @noRd
bm_plot_ranking_heatmap <- function(ranking_result,
                                    type  = "mean",
                                    title = NULL) {
  rank_df <- if (type == "median") {
    ranking_result$median_ranking
  } else {
    ranking_result$mean_ranking
  }

  metrics <- ranking_result$metrics_used
  rank_cols <- paste0("rank_", metrics)

  # Pivot to long
  plot_data <- rank_df[, c("Assay", rank_cols, "rank_final"), drop = FALSE]

  # Order Assay by rank_final
  plot_data$Assay <- factor(plot_data$Assay,
                            levels = rev(rank_df$Assay))

  long <- tidyr::pivot_longer(plot_data,
                              cols      = c(rank_cols, "rank_final"),
                              names_to  = "Metric",
                              values_to = "Rank")

  # Clean metric names for display
  long$Metric <- sub("^rank_", "", long$Metric)
  metric_order <- c(metrics, "final")
  long$Metric <- factor(long$Metric, levels = metric_order)

  n_methods <- ranking_result$n_methods

  title <- title %||% paste0("OpDEA Ranking Heatmap (", type, "-based)")
  subtitle <- paste0(ranking_result$n_methods, " methods, ",
                     ranking_result$n_comparisons, " comparisons")

  gg <- ggplot2::ggplot(long, ggplot2::aes(x = Metric, y = Assay, fill = Rank)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(Rank, 1)),
                       size = 3.5, color = "black") +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid        = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(face = "bold"),
      legend.position   = "right"
    )

  gg
}


#' Bar chart of final ranks
#'
#' Horizontal bars of rank_final by Assay, ordered best to worst.
#'
#' @param ranking_result List returned by \code{bm_compute_ranking()}
#' @param type "mean" or "median" aggregation (default: "mean")
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @keywords internal
#' @noRd
bm_plot_ranking_bars <- function(ranking_result,
                                 type  = "mean",
                                 title = NULL) {
  rank_df <- if (type == "median") {
    ranking_result$median_ranking
  } else {
    ranking_result$mean_ranking
  }

  n_methods <- ranking_result$n_methods

  rank_df$Assay <- factor(rank_df$Assay,
                          levels = rev(rank_df$Assay))

  title <- title %||% paste0("OpDEA Final Ranking (", type, "-based)")
  subtitle <- paste0(ranking_result$n_methods, " methods, ",
                     ranking_result$n_comparisons, " comparisons")

  gg <- ggplot2::ggplot(rank_df,
                        ggplot2::aes(x = Assay, y = rank_final, fill = rank_final)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = round(rank_final, 1)),
                       hjust = -0.2, size = 3.5) +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = NULL, y = "Final Rank (lower = better)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )

  gg
}


#' Heatmap of aggregated metric values (not ranks)
#'
#' Shows the actual metric values per method, with per-column color scaling.
#'
#' @param ranking_result List returned by \code{bm_compute_ranking()}
#' @param type "mean" or "median" aggregation (default: "mean")
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @keywords internal
#' @noRd
bm_plot_metrics_heatmap <- function(ranking_result,
                                    type  = "mean",
                                    title = NULL) {
  agg_df <- if (type == "median") {
    ranking_result$median_aggregated
  } else {
    ranking_result$mean_aggregated
  }

  # Order by mean ranking
  rank_df <- if (type == "median") {
    ranking_result$median_ranking
  } else {
    ranking_result$mean_ranking
  }

  metrics <- ranking_result$metrics_used

  # Pivot to long
  long <- tidyr::pivot_longer(agg_df,
                              cols      = dplyr::all_of(metrics),
                              names_to  = "Metric",
                              values_to = "Value")

  long$Assay  <- factor(long$Assay, levels = rev(rank_df$Assay))
  long$Metric <- factor(long$Metric, levels = metrics)

  # Normalize values per metric to [0, 1] for consistent color scale
  long <- dplyr::group_by(long, Metric)
  long <- dplyr::mutate(long,
    Value_scaled = if (max(Value, na.rm = TRUE) == min(Value, na.rm = TRUE)) {
      0.5
    } else {
      (Value - min(Value, na.rm = TRUE)) /
        (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))
    }
  )
  long <- dplyr::ungroup(long)

  title <- title %||% paste0("OpDEA Metric Values (", type, "-based)")
  subtitle <- paste0(ranking_result$n_methods, " methods, ",
                     ranking_result$n_comparisons, " comparisons | ",
                     "color scaled per metric")

  gg <- ggplot2::ggplot(long,
                        ggplot2::aes(x = Metric, y = Assay, fill = Value_scaled)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4f", Value)),
                       size = 3, color = "black") +
    ggplot2::scale_fill_gradient(low = "#fee0d2", high = "#2ca02c",
                                 name = "Relative\nPerformance",
                                 labels = c("Worst", "Best"),
                                 breaks = c(0, 1)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid        = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(face = "bold"),
      legend.position   = "right"
    )

  gg
}


#' Boxplot comparison of metric distributions across methods
#'
#' Faceted boxplot showing the distribution of each metric across
#' comparisons for every method.
#'
#' @param opdea_combined data.frame with Assay, Comparison, and metric columns
#' @param metrics Character vector of metrics to plot (default: all 5)
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @keywords internal
#' @noRd
bm_plot_metrics_comparison <- function(opdea_combined,
                                       metrics = .BM_METRICS,
                                       title   = NULL) {
  avail <- intersect(metrics, colnames(opdea_combined))
  if (length(avail) == 0) stop("No valid metrics found in opdea_combined.")

  long <- tidyr::pivot_longer(opdea_combined,
                              cols      = dplyr::all_of(avail),
                              names_to  = "Metric",
                              values_to = "Value")

  long$Metric <- factor(long$Metric, levels = avail)

  n_methods <- length(unique(long$Assay))
  n_comps   <- length(unique(long$Comparison))

  title <- title %||% "OpDEA Metrics Distribution by Method"
  subtitle <- paste0(n_methods, " methods, ", n_comps, " comparisons")

  # Use distinct shapes + colors for comparisons on the jittered points
  comp_levels <- sort(unique(long$Comparison))
  long$Comparison <- factor(long$Comparison, levels = comp_levels)

  # Select a colorblind-friendly palette
  if (n_comps <= 8) {
    comp_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                     "#0072B2", "#D55E00", "#CC79A7", "#999999")[seq_len(n_comps)]
  } else {
    comp_colors <- grDevices::hcl.colors(n_comps, palette = "Dynamic")
  }
  names(comp_colors) <- comp_levels

  # Shapes: cycle through distinguishable filled shapes
  shape_pool <- c(16, 17, 15, 18, 8, 4, 3, 7, 9, 10, 12, 13, 14)
  comp_shapes <- shape_pool[((seq_len(n_comps) - 1) %% length(shape_pool)) + 1]
  names(comp_shapes) <- comp_levels

  # Pastel fill palette for boxplots (one per Assay), with matching darker borders
  assay_levels <- unique(long$Assay)
  n_assays <- length(assay_levels)
  assay_fills <- grDevices::hcl.colors(n_assays, palette = "Pastel 1")
  names(assay_fills) <- assay_levels
  # Derive darker border colors by reducing luminance
  assay_borders <- vapply(assay_fills, function(hex) {
    rgb_vals <- grDevices::col2rgb(hex)[, 1] / 255
    darker <- pmax(rgb_vals * 0.55, 0)
    grDevices::rgb(darker[1], darker[2], darker[3])
  }, character(1))
  names(assay_borders) <- assay_levels

  # Map border color per row for geom_boxplot
  long$assay_border <- assay_borders[as.character(long$Assay)]

  gg <- ggplot2::ggplot(long, ggplot2::aes(x = Assay, y = Value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = Assay),
                          alpha = 0.5, outlier.shape = NA, linewidth = 0.4,
                          show.legend = FALSE) +
    ggplot2::geom_point(
      ggplot2::aes(color = Comparison, shape = Comparison),
      position = ggplot2::position_jitter(width = 0.15, seed = 42),
      size = 2.2, alpha = 0.85
    ) +
    ggplot2::scale_fill_manual(values = assay_fills) +
    ggplot2::scale_color_manual(values = comp_colors) +
    ggplot2::scale_shape_manual(values = comp_shapes) +
    ggplot2::facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = NULL, y = "Value",
                  color = "Comparison", shape = "Comparison") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
      plot.title        = ggplot2::element_text(face = "bold"),
      legend.position   = "bottom",
      legend.title      = ggplot2::element_text(face = "bold", size = 9),
      legend.text       = ggplot2::element_text(size = 8),
      strip.text        = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 3)),
      shape = ggplot2::guide_legend(nrow = 1)
    )

  gg
}


#' Heatmap of metric ranks for a single comparison
#'
#' @param ranking_by_comp List returned by \code{bm_compute_ranking_by_comparison()}
#' @param comparison Character string: which comparison to plot
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @keywords internal
#' @noRd
bm_plot_ranking_heatmap_by_comp <- function(ranking_by_comp,
                                            comparison,
                                            title = NULL) {
  if (!comparison %in% names(ranking_by_comp$by_comparison))
    stop("Comparison '", comparison, "' not found. Available: ",
         paste(names(ranking_by_comp$by_comparison), collapse = ", "))

  comp_df <- ranking_by_comp$by_comparison[[comparison]]

  rank_cols <- grep("^rank_", colnames(comp_df), value = TRUE)
  plot_data <- comp_df[, c("Assay", rank_cols), drop = FALSE]

  plot_data$Assay <- factor(plot_data$Assay, levels = rev(unique(comp_df$Assay)))

  long <- tidyr::pivot_longer(plot_data,
                              cols      = rank_cols,
                              names_to  = "Metric",
                              values_to = "Rank")

  long$Metric <- sub("^rank_", "", long$Metric)
  metric_order <- sub("^rank_", "", rank_cols)
  long$Metric <- factor(long$Metric, levels = metric_order)

  n_methods <- nrow(comp_df)

  title <- title %||% paste0("OpDEA Ranking: ", comparison)

  gg <- ggplot2::ggplot(long, ggplot2::aes(x = Metric, y = Assay, fill = Rank)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(Rank, 1)),
                       size = 3.5, color = "black") +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x    = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid     = ggplot2::element_blank(),
      plot.title     = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  gg
}


#' Bar chart of final ranks for a single comparison
#'
#' @param ranking_by_comp List returned by \code{bm_compute_ranking_by_comparison()}
#' @param comparison Character string: which comparison to plot
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @keywords internal
#' @noRd
bm_plot_ranking_bars_by_comp <- function(ranking_by_comp,
                                         comparison,
                                         title = NULL) {
  if (!comparison %in% names(ranking_by_comp$by_comparison))
    stop("Comparison '", comparison, "' not found. Available: ",
         paste(names(ranking_by_comp$by_comparison), collapse = ", "))

  comp_df <- ranking_by_comp$by_comparison[[comparison]]
  n_methods <- nrow(comp_df)

  comp_df$Assay <- factor(comp_df$Assay, levels = rev(unique(comp_df$Assay)))

  title <- title %||% paste0("OpDEA Final Ranking: ", comparison)

  gg <- ggplot2::ggplot(comp_df,
                        ggplot2::aes(x = Assay, y = rank_final, fill = rank_final)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = round(rank_final, 1)),
                       hjust = -0.2, size = 3.5) +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = title, x = NULL, y = "Final Rank (lower = better)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )

  gg
}


#' Stacked bar chart of TP/FP/FN/TN (OpDEA Figure 5 style)
#'
#' Horizontal stacked bars showing confusion matrix counts per method.
#' If \code{comparison} is NULL, generates a faceted plot with all comparisons.
#'
#' @param confusion_combined data.frame with columns:
#'   Assay, Comparison, TP, FP, TN, FN
#' @param comparison Optional: single comparison to plot.
#'   If NULL, faceted plot with all comparisons.
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @keywords internal
#' @noRd
bm_plot_confusion_stacked <- function(confusion_combined,
                                      comparison = NULL,
                                      title      = NULL) {
  .bm_validate_confusion(confusion_combined)

  # Filter to single comparison if requested
  if (!is.null(comparison)) {
    confusion_combined <- confusion_combined[
      confusion_combined$Comparison == comparison, , drop = FALSE]
    if (nrow(confusion_combined) == 0)
      stop("Comparison '", comparison, "' not found.")
  }

  # OpDEA colors: TP=blue, TN=coral, FP=teal, FN=gold
  conf_colors <- c(TP = "#5494cc", TN = "#e18283", FP = "#0d898a", FN = "#f9cc52")
  categories <- c("TP", "TN", "FP", "FN")

  # Pivot to long
  long <- tidyr::pivot_longer(
    confusion_combined[, c("Assay", "Comparison", categories), drop = FALSE],
    cols      = categories,
    names_to  = "Category",
    values_to = "Count"
  )
  long$Category <- factor(long$Category, levels = categories)

  # Order methods by TP count (descending) within each comparison
  # Use mean TP across comparisons for consistent ordering
  tp_order <- tapply(
    confusion_combined$TP,
    confusion_combined$Assay,
    mean, na.rm = TRUE
  )
  assay_order <- names(sort(tp_order, decreasing = FALSE))
  long$Assay <- factor(long$Assay, levels = assay_order)

  # Only show label if segment is large enough (>2% of total per bar)
  total_per_bar <- tapply(long$Count, list(long$Assay, long$Comparison),
                          sum, default = 0)
  long$label_text <- vapply(seq_len(nrow(long)), function(i) {
    total <- total_per_bar[as.character(long$Assay[i]),
                           as.character(long$Comparison[i])]
    if (!is.na(long$Count[i]) && long$Count[i] > total * 0.02) {
      as.character(long$Count[i])
    } else {
      ""
    }
  }, character(1))

  title <- if (!is.null(title)) {
    title
  } else if (!is.null(comparison)) {
    paste0("Confusion Matrix: ", comparison)
  } else {
    "Confusion Matrix by Method"
  }

  gg <- ggplot2::ggplot(long,
                        ggplot2::aes(x = Assay, y = Count,
                                     fill = Category)) +
    ggplot2::geom_col(width = 0.75, color = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = label_text),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3, color = "black", fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = conf_colors,
                                name = NULL) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = "Number of proteins") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      axis.line.y     = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 9),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

  # Facet if multiple comparisons
  if (is.null(comparison)) {
    gg <- gg + ggplot2::facet_wrap(~ Comparison, scales = "free_x")
  }

  gg
}


#' ROC Curves by Method for a Single Comparison
#'
#' Generates ROC (or ROC zoom) plot where each curve represents
#' a different method (Assay) for one comparison.
#'
#' @param classified_combined data.frame with columns: Assay, Comparison,
#'   truth, and a p-value column
#' @param comparison Character. Single comparison to plot (e.g. "B-A")
#' @param p_col P-value column name (default: "adj.P.Val")
#' @param zoom Logical. If TRUE, zoom to FPR 0-10% and show pAUC in legend
#' @param palette Colour palette: `"pkg::name"` for a
#'   palette from paletteer; `"brewer:Name"`, or the bare name, for an
#'   RColorBrewer palette; or a vector of colours. A palette with
#'   fewer colours than curves is interpolated. Default `"Set2"`.
#'
#' @return ggplot2 object, or NULL (with a warning) if the comparison has no
#'   data from which to draw a ROC curve
#' @keywords internal
#' @noRd
bm_plot_roc <- function(classified_combined,
                        comparison,
                        p_col   = "adj.P.Val",
                        zoom    = FALSE,
                        palette = "Set2") {

  # Filter to the requested comparison
  df <- classified_combined[classified_combined$Comparison == comparison, , drop = FALSE]

  if (nrow(df) == 0) {
    warning("No data for comparison '", comparison, "'")
    return(NULL)
  }

  # Compute score: -log10(p-value)
  df$score <- -log10(pmax(as.numeric(df[[p_col]]), 1e-300))

  # Filter valid rows
  df <- df[!is.na(df$truth) & !is.na(df$score), , drop = FALSE]

  if (nrow(df) == 0) {
    warning("No valid data for ROC in comparison '", comparison, "'")
    return(NULL)
  }

  # Build ROC objects per Assay
  assay_order <- unique(df$Assay)
  df$Assay <- factor(df$Assay, levels = assay_order)
  assay_list <- split(df, df$Assay)

  roc_list <- lapply(assay_list, function(d) {
    if (length(unique(d$truth)) < 2 || nrow(d) < 10) return(NULL)
    tryCatch(
      # direction="<": a high score (-log10 p) = case; the direction is fixed
      # (not "auto") for consistency with the tabulated pAUC.
      pROC::roc(response = d$truth, predictor = d$score,
                direction = "<", quiet = TRUE),
      error = function(e) NULL
    )
  })
  roc_list <- Filter(Negate(is.null), roc_list)

  if (length(roc_list) == 0) {
    warning("Could not generate ROC curves for comparison '", comparison, "'")
    return(NULL)
  }

  # AUC / pAUC labels
  if (zoom) {
    # McClish-corrected pAUC over the FPR 0-5% region
    aucs <- vapply(roc_list, function(r) {
      tryCatch(
        as.numeric(pROC::auc(r,
                              partial.auc = c(1, 0.95),
                              partial.auc.correct = TRUE)),
        error = function(e) NA_real_
      )
    }, numeric(1))
    # Sort the legend by decreasing pAUC
    ord      <- order(aucs, decreasing = TRUE, na.last = TRUE)
    roc_list <- roc_list[ord]
    aucs     <- aucs[ord]
    labels <- paste0(names(roc_list), " (pAUC=", sprintf("%.3f", aucs), ")")
  } else {
    aucs <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)), numeric(1))
    labels <- paste0(names(roc_list), " (AUC=", sprintf("%.3f", aucs), ")")
  }

  # Build plot
  title_text <- if (zoom) {
    paste0("ROC Zoom - ", comparison)
  } else {
    paste0("ROC - ", comparison)
  }

  gg <- pROC::ggroc(roc_list, legacy.axes = TRUE, linewidth = 1) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "gray50", linewidth = 0.5
    ) +
    ggplot2::labs(
      title    = title_text,
      subtitle = if (zoom) {
        "Low FPR region (axis 0-10%) - legend pAUC at 5% FPR"
      } else {
        "Diagonal = random classifier"
      },
      x     = "False Positive Rate (1 - Specificity)",
      y     = "True Positive Rate (Sensitivity)",
      color = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold",
                                            size = 15, color = "#1D3557"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11,
                                            color = "#6C757D"),
      legend.position = "right",
      legend.text     = ggplot2::element_text(size = 10),
      legend.title    = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Apply palette (interpolated when there are more methods than colours)
  pal_colors <- .resolve_palette(palette, n = length(roc_list))
  gg <- gg + ggplot2::scale_color_manual(values = pal_colors, labels = labels)

  # Zoom or full view
  if (zoom) {
    gg <- gg + ggplot2::coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 1))
  } else {
    gg <- gg + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  }

  gg
}


# =============================================================================
# SECTION 5: ORCHESTRATOR
# =============================================================================

# --- Steps of benchmarking_multiple() ------------------------------------------
# The orchestrator runs these in order: gather the inputs, rank the methods,
# draw the plots, export, and put the result together. Each is kept apart so it
# can be tested and reused on its own.

#' Build a plot or a table, turning a failure into a warning
#' @param label What failed, for the warning.
#' @param f Function of no arguments that builds it.
#' @return Its value, or `NULL` if it failed.
#' @keywords internal
#' @noRd
.bm_try <- function(label, f) {
  tryCatch(f(), error = function(e) {
    warning(label, " failed: ", conditionMessage(e))
    NULL
  })
}

#' One optional input, imported from `results_dir` when not given, then checked
#' @param df The table given by the user, or `NULL`.
#' @param results_dir Folder to import it from, or `NULL`.
#' @param importer Function of the folder that imports it.
#' @param validate Function that stops when the table is not valid.
#' @param what Name of the data, for the messages.
#' @param arg Name of the argument, for the warning.
#' @param verbose Print progress messages.
#' @return The table, or `NULL` when it is absent or not valid.
#' @keywords internal
#' @noRd
.bm_optional_input <- function(df, results_dir, importer, validate, what, arg,
                               verbose) {
  if (is.null(df) && !is.null(results_dir)) {
    df <- tryCatch({
      if (verbose) message("  Importing ", what, " data ...")
      importer(results_dir)
    }, error = function(e) {
      if (verbose) message("  No ", what, " data found (skipping): ", e$message)
      NULL
    })
  }
  if (is.null(df)) return(NULL)
  tryCatch({
    validate(df)
    df
  }, error = function(e) {
    warning("Invalid ", arg, " (skipping): ", e$message)
    NULL
  })
}

#' Step 1: gather and check the inputs of benchmarking_multiple()
#'
#' The OpDEA metrics are required; the confusion matrices, the classified
#' proteins (for the ROC curves) and the benchmark metrics (for the extended
#' ranking) are optional, and an absent or invalid one is `NULL`.
#'
#' @inheritParams benchmarking_multiple
#' @return List with `opdea`, `confusion`, `classified` and `bench_metrics`.
#' @keywords internal
#' @noRd
.bm_gather_inputs <- function(opdea_combined, confusion_combined,
                              classified_combined, bench_metrics_combined,
                              results_dir, pattern, method_names, recursive,
                              strip_prefix, p_col, verbose) {
  if (!is.null(opdea_combined) && !is.null(results_dir)) {
    warning("Both opdea_combined and results_dir provided. Using opdea_combined.",
            call. = FALSE)
  }
  if (is.null(opdea_combined)) {
    if (is.null(results_dir))
      stop("Provide either opdea_combined or results_dir.", call. = FALSE)
    if (verbose) message("=== IMPORTING OPDEA RESULTS ===")
    opdea_combined <- import_opdea_results(results_dir, pattern,
                                           method_names, recursive)
  }
  .bm_validate_opdea(opdea_combined)

  from_dir <- function(importer) {
    function(dir) importer(dir, method_names = method_names, recursive = recursive)
  }
  inputs <- list(
    opdea      = opdea_combined,
    confusion  = .bm_optional_input(
      confusion_combined, results_dir, from_dir(import_confusion_results),
      .bm_validate_confusion, "confusion", "confusion_combined", verbose),
    classified = .bm_optional_input(
      classified_combined, results_dir, from_dir(import_classified_results),
      function(df) .bm_validate_classified(df, p_col),
      "classified", "classified_combined", verbose),
    bench_metrics = .bm_optional_input(
      bench_metrics_combined, results_dir, from_dir(import_benchmark_metrics),
      .bm_validate_bench_metrics, "benchmark metrics", "bench_metrics_combined",
      verbose)
  )

  # Drop the folder prefix (e.g. "benchmark_") from Assay. In results_dir mode
  # Assay is the folder name (benchmark_<method>); this makes it match the
  # in-memory mode, and is harmless if the Assay carries no such prefix.
  if (!is.null(strip_prefix) && nzchar(strip_prefix)) {
    inputs <- lapply(inputs, function(df) {
      if (!is.null(df) && "Assay" %in% names(df))
        df$Assay <- sub(paste0("^", strip_prefix), "", df$Assay)
      df
    })
  }

  if (verbose) {
    message("=== MULTIPLE-METHOD OPDEA BENCHMARKING ===")
    message("  Methods: ", length(unique(inputs$opdea$Assay)),
            " | Comparisons: ", length(unique(inputs$opdea$Comparison)),
            " | Rows: ", nrow(inputs$opdea))
    if (!is.null(inputs$confusion)) message("  Confusion data: available")
    if (!is.null(inputs$classified)) message("  Classified data: available (ROC curves)")
    if (!is.null(inputs$bench_metrics)) message("  Benchmark metrics data: available (extended ranking)")
  }
  inputs
}

#' Step 2: rank the methods, overall, by comparison and on the 11 metrics
#' @param inputs Output of `.bm_gather_inputs()`.
#' @inheritParams benchmarking_multiple
#' @return List with `ranking`, `by_comp`, `extended` (or `NULL`) and `comps`.
#' @keywords internal
#' @noRd
.bm_rankings <- function(inputs, metrics, extended_metrics, verbose) {
  if (verbose) message("  Computing ranking ...")
  ranking <- bm_compute_ranking(inputs$opdea, metrics)
  if (verbose) {
    message("  Best method (mean):   ", ranking$mean_ranking$Assay[1],
            " (rank_final = ", ranking$mean_ranking$rank_final[1], ")")
    message("  Best method (median): ", ranking$median_ranking$Assay[1],
            " (rank_final = ", ranking$median_ranking$rank_final[1], ")")
  }

  if (verbose) message("  Computing ranking by comparison ...")
  by_comp <- bm_compute_ranking_by_comparison(inputs$opdea, metrics)
  comps <- names(by_comp$by_comparison)
  if (verbose) {
    for (comp in comps) {
      best <- by_comp$by_comparison[[comp]]
      message("    ", comp, ": best = ", best$Assay[1],
              " (rank_final = ", best$rank_final[1], ")")
    }
  }

  extended <- NULL
  if (!is.null(inputs$bench_metrics)) {
    if (verbose) message("  Computing extended ranking (11 metrics) ...")
    extended <- .bm_try("bm_compute_extended_ranking()", function()
      bm_compute_extended_ranking(inputs$opdea, inputs$bench_metrics,
                                  extended_metrics))
    if (!is.null(extended) && verbose) {
      message("  Best method (extended): ", extended$mean_ranking$Assay[1],
              " (rank_final = ", extended$mean_ranking$rank_final[1], ")")
    }
  }

  list(ranking = ranking, by_comp = by_comp, extended = extended, comps = comps)
}

#' Step 3: draw the plots
#'
#' A plot that fails is a warning and a `NULL`, never an error.
#'
#' @param inputs Output of `.bm_gather_inputs()`.
#' @param ranks Output of `.bm_rankings()`.
#' @inheritParams benchmarking_multiple
#' @return List of the plots, plus `selected` (the global plots requested) and
#'   the counts reported in the messages.
#' @keywords internal
#' @noRd
.bm_draw <- function(inputs, ranks, metrics, plots, p_col, verbose) {
  ranking  <- ranks$ranking
  extended <- ranks$extended
  plot_fns <- list(
    ranking_heatmap_mean   = function() bm_plot_ranking_heatmap(ranking, "mean"),
    ranking_heatmap_median = function() bm_plot_ranking_heatmap(ranking, "median"),
    ranking_bars_mean      = function() bm_plot_ranking_bars(ranking, "mean"),
    ranking_bars_median    = function() bm_plot_ranking_bars(ranking, "median"),
    metrics_heatmap_mean   = function() bm_plot_metrics_heatmap(ranking, "mean"),
    metrics_heatmap_median = function() bm_plot_metrics_heatmap(ranking, "median"),
    metrics_comparison     = function() bm_plot_metrics_comparison(inputs$opdea, metrics)
  )
  if (!is.null(extended)) {
    plot_fns$extended_ranking_bars    <- function() bm_plot_extended_ranking_bars(extended)
    plot_fns$extended_ranking_heatmap <- function() bm_plot_extended_ranking_heatmap(extended)
  }

  if (identical(plots, "all")) {
    selected <- names(plot_fns)
  } else {
    unknown <- setdiff(plots, names(plot_fns))
    if (length(unknown) > 0)
      warning("Unknown plot name(s) ignored: ", paste(unknown, collapse = ", "),
              call. = FALSE)
    selected <- intersect(plots, names(plot_fns))
    if (length(selected) == 0)
      stop("No valid plot names provided.", call. = FALSE)
  }

  n_drawn <- function(x) sum(!vapply(x, is.null, logical(1)))
  out <- list(selected = selected)

  # --- Global plots ---
  out$global <- lapply(stats::setNames(nm = selected), function(nm) {
    if (verbose) message("  Generating plot: ", nm, " ...")
    .bm_try(paste0("bm_plot_", nm, "()"), plot_fns[[nm]])
  })
  out$n_ok <- n_drawn(out$global)
  if (verbose) {
    n_fail <- length(selected) - out$n_ok
    if (n_fail > 0)
      message("  ", out$n_ok, " global plot(s) generated, ", n_fail, " failed")
    else
      message("  ", out$n_ok, " global plot(s) generated.")
  }

  # --- Per-comparison ranking plots ---
  per_comp <- lapply(stats::setNames(nm = ranks$comps), function(comp) {
    if (verbose) message("  Generating plots for comparison: ", comp, " ...")
    list(
      heatmap = .bm_try(paste0("bm_plot_ranking_heatmap_by_comp(", comp, ")"),
                        function() bm_plot_ranking_heatmap_by_comp(ranks$by_comp, comp)),
      bars    = .bm_try(paste0("bm_plot_ranking_bars_by_comp(", comp, ")"),
                        function() bm_plot_ranking_bars_by_comp(ranks$by_comp, comp)))
  })
  out$heatmap_by_comp <- lapply(per_comp, `[[`, "heatmap")
  out$bars_by_comp    <- lapply(per_comp, `[[`, "bars")
  out$n_comp_plots <- n_drawn(out$heatmap_by_comp) + n_drawn(out$bars_by_comp)
  if (verbose) message("  ", out$n_comp_plots, " per-comparison plot(s) generated.")

  # --- Confusion stacked bars: all comparisons, then one per comparison ---
  confusion <- inputs$confusion
  if (!is.null(confusion)) {
    if (verbose) message("  Generating confusion stacked bar plots ...")
    out$confusion_stacked <- .bm_try("bm_plot_confusion_stacked()",
                                     function() bm_plot_confusion_stacked(confusion))
    out$confusion_by_comp <- lapply(
      stats::setNames(nm = unique(confusion$Comparison)), function(comp)
        .bm_try(paste0("bm_plot_confusion_stacked(", comp, ")"),
                function() bm_plot_confusion_stacked(confusion, comparison = comp)))
    if (verbose)
      message("  ", (!is.null(out$confusion_stacked)) + n_drawn(out$confusion_by_comp),
              " confusion plot(s) generated.")
  }

  # --- ROC curves per comparison, one curve per method ---
  classified <- inputs$classified
  if (!is.null(classified)) {
    if (verbose) message("  Generating ROC curves per comparison ...")
    roc <- lapply(stats::setNames(nm = unique(classified$Comparison)), function(comp)
      list(
        full = .bm_try(paste0("bm_plot_roc(", comp, ")"), function()
          bm_plot_roc(classified, comp, p_col = p_col, zoom = FALSE)),
        zoom = .bm_try(paste0("bm_plot_roc(", comp, ", zoom)"), function()
          bm_plot_roc(classified, comp, p_col = p_col, zoom = TRUE))))
    out$roc_by_comp      <- lapply(roc, `[[`, "full")
    out$roc_zoom_by_comp <- lapply(roc, `[[`, "zoom")
    if (verbose)
      message("  ", n_drawn(out$roc_by_comp) + n_drawn(out$roc_zoom_by_comp),
              " ROC plot(s) generated.")
  }
  out
}

#' Step 4: write the tables and plots to `output_dir`
#' @param inputs Output of `.bm_gather_inputs()`.
#' @param ranks Output of `.bm_rankings()`.
#' @param figs Output of `.bm_draw()`.
#' @inheritParams benchmarking_multiple
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.bm_export_all <- function(output_dir, inputs, ranks, figs, export_tables,
                           export_plots, plot_width, plot_height, plot_dpi,
                           verbose) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (verbose) message("  Exporting to: ", output_dir)
  comp_dir <- file.path(output_dir, "by_comparison")
  if ((export_tables || export_plots) && !dir.exists(comp_dir))
    dir.create(comp_dir, recursive = TRUE)

  safe <- function(comp) gsub("[^A-Za-z0-9_-]", "_", comp)
  put_tsv <- function(df, file, dir = output_dir)
    .bm_export_data(df, file.path(dir, file))
  put_png <- function(p, file, dir = comp_dir)
    if (!is.null(p))
      .bm_export_gg_plot(p, file.path(dir, file), plot_width, plot_height, plot_dpi)
  ranking  <- ranks$ranking
  extended <- ranks$extended

  if (export_tables) {
    put_tsv(inputs$opdea, "bm_multiple_opdea_combined.tsv")
    for (s in c("mean_aggregated", "median_aggregated",
                "mean_ranking", "median_ranking"))
      put_tsv(ranking[[s]], paste0("bm_multiple_", s, ".tsv"))
    put_tsv(ranks$by_comp$ranking_combined, "bm_multiple_ranking_by_comparison.tsv")
    for (comp in ranks$comps)
      put_tsv(ranks$by_comp$by_comparison[[comp]],
              paste0("bm_ranking_", safe(comp), ".tsv"), comp_dir)
    if (!is.null(extended)) {
      put_tsv(extended$extended_combined, "bm_multiple_extended_combined.tsv")
      put_tsv(extended$mean_aggregated, "bm_multiple_extended_mean_aggregated.tsv")
      put_tsv(extended$mean_ranking, "bm_multiple_extended_ranking.tsv")
    }
    if (verbose)
      message("  Exported ", 6 + length(ranks$comps) + 3 * !is.null(extended),
              " TSV files.")
  }

  if (export_plots) {
    for (nm in figs$selected)
      put_png(figs$global[[nm]], paste0("bm_multiple_", nm, ".png"), output_dir)
    for (comp in ranks$comps) {
      put_png(figs$heatmap_by_comp[[comp]], paste0("bm_ranking_heatmap_", safe(comp), ".png"))
      put_png(figs$bars_by_comp[[comp]], paste0("bm_ranking_bars_", safe(comp), ".png"))
    }
    if (verbose) message("  Exported ", figs$n_ok, " global + ",
                         figs$n_comp_plots, " per-comparison PNG files.")
  }

  if (!is.null(inputs$confusion)) {
    if (export_tables)
      put_tsv(inputs$confusion, "bm_multiple_confusion_combined.tsv")
    if (export_plots) {
      put_png(figs$confusion_stacked, "bm_multiple_confusion_stacked.png", output_dir)
      for (comp in names(figs$confusion_by_comp))
        put_png(figs$confusion_by_comp[[comp]],
                paste0("bm_confusion_stacked_", safe(comp), ".png"))
    }
    if (verbose) message("  Exported confusion data and plots.")
  }

  if (!is.null(inputs$classified)) {
    if (export_tables)
      put_tsv(inputs$classified, "bm_multiple_classified_combined.tsv")
    if (export_plots) {
      for (comp in names(figs$roc_by_comp)) {
        put_png(figs$roc_by_comp[[comp]], paste0("bm_roc_", safe(comp), ".png"))
        put_png(figs$roc_zoom_by_comp[[comp]], paste0("bm_roc_zoom_", safe(comp), ".png"))
      }
    }
    if (verbose) message("  Exported classified data and ROC plots.")
  }
  invisible(NULL)
}

#' Step 5: the list benchmarking_multiple() returns
#' @param inputs Output of `.bm_gather_inputs()`.
#' @param ranks Output of `.bm_rankings()`.
#' @param figs Output of `.bm_draw()`.
#' @param metrics The OpDEA metrics used.
#' @return See the Value section of [benchmarking_multiple()].
#' @keywords internal
#' @noRd
.bm_result <- function(inputs, ranks, figs, metrics) {
  ranking  <- ranks$ranking
  extended <- ranks$extended
  result <- list(
    opdea_combined    = inputs$opdea,
    mean_aggregated   = ranking$mean_aggregated,
    median_aggregated = ranking$median_aggregated,
    mean_ranking      = ranking$mean_ranking,
    median_ranking    = ranking$median_ranking,
    ranking_by_comparison          = ranks$by_comp$by_comparison,
    ranking_by_comparison_combined = ranks$by_comp$ranking_combined
  )

  # Global plots with the gg_ prefix, then the per-comparison ones
  for (nm in figs$selected) {
    result[[paste0("gg_", nm)]] <- figs$global[[nm]]
  }
  result$gg_ranking_heatmap_by_comp <- figs$heatmap_by_comp
  result$gg_ranking_bars_by_comp    <- figs$bars_by_comp

  if (!is.null(inputs$confusion)) {
    result$confusion_combined           <- inputs$confusion
    result$gg_confusion_stacked         <- figs$confusion_stacked
    result$gg_confusion_stacked_by_comp <- figs$confusion_by_comp
  }
  if (!is.null(inputs$classified)) {
    result$classified_combined <- inputs$classified
    result$gg_roc_by_comp      <- figs$roc_by_comp
    result$gg_roc_zoom_by_comp <- figs$roc_zoom_by_comp
  }
  if (!is.null(extended)) {
    result$extended_combined        <- extended$extended_combined
    result$extended_mean_aggregated <- extended$mean_aggregated
    result$extended_ranking         <- extended$mean_ranking
    result$bench_metrics_combined   <- inputs$bench_metrics
  }

  result$parameters <- list(
    metrics          = metrics,
    extended_metrics = if (!is.null(extended)) extended$metrics_used else NULL,
    n_methods        = length(unique(inputs$opdea$Assay)),
    n_comparisons    = length(unique(inputs$opdea$Comparison)),
    comparisons      = ranks$comps,
    plots            = figs$selected
  )
  result
}


#' Multiple-method OpDEA Benchmarking
#'
#' Orchestrator that imports, ranks, visualizes, and exports results from
#' multiple normalization/imputation method combinations.
#'
#' @param opdea_combined Optional pre-built data.frame with columns:
#'   Assay, Comparison, nMCC, G_mean, pAUC_001, pAUC_005, pAUC_010.
#'   If NULL, results are imported from \code{results_dir}.
#' @param results_dir Parent directory containing benchmark subfolders
#'   (each with \code{benchmark_opdea_metrics.tsv}). Used only if
#'   \code{opdea_combined} is NULL.
#' @param pattern Regex for opdea metrics filename (default:
#'   \code{"benchmark_opdea_metrics\\.tsv$"})
#' @param method_names Optional Assay names for file import
#' @param recursive Search subdirectories? (default: TRUE)
#' @param strip_prefix Prefix removed from each \code{Assay} name (default
#'   "benchmark_"), so folder-derived names in \code{results_dir} mode match the
#'   in-memory ones (e.g. "benchmark_Rlr_none_BERT" -> "Rlr_none_BERT"). Set to
#'   NULL or "" to keep the full folder name. Harmless when the prefix is absent.
#' @param bench_metrics_combined Optional pre-built data.frame, the
#'   \code{metrics_table} of several \code{benchmarking_proteomics()} runs
#'   stacked with an Assay column. Columns: Assay, Comparison,
#'   Sensitivity, Specificity, Precision, NPV, F1, Accuracy, MCC, Performance.
#'   If NULL and \code{results_dir} is provided, imported from
#'   \code{benchmark_metrics.tsv} files. Used for extended 11-metric ranking.
#' @param metrics Character vector of metrics for OpDEA ranking
#'   (default: all 5 OpDEA metrics)
#' @param extended_metrics Character vector of metrics for extended ranking
#'   (default: all 11 classification + OpDEA metrics)
#' @param confusion_combined Optional pre-built data.frame with columns:
#'   Assay, Comparison, TP, FP, TN, FN. If NULL and \code{results_dir}
#'   is provided, imported from \code{benchmark_confusion_overall.tsv} files.
#' @param classified_combined Optional pre-built data.frame with columns:
#'   Assay, Comparison, truth, and a p-value column. If NULL and
#'   \code{results_dir} is provided, imported from
#'   \code{benchmark_classified.tsv} files. Used for multi-method ROC curves.
#' @param p_col P-value column name for ROC curves (default: "adj.P.Val")
#' @param plots Which plots to generate: "all" or character vector of names.
#'   Valid names: "ranking_heatmap_mean", "ranking_heatmap_median",
#'   "ranking_bars_mean", "ranking_bars_median",
#'   "metrics_heatmap_mean", "metrics_heatmap_median",
#'   "metrics_comparison", "extended_ranking_bars",
#'   "extended_ranking_heatmap"
#' @param verbose Print progress messages (default: TRUE)
#' @param output_dir Directory for exporting results (NULL = no export)
#' @param export_plots Export plots as PNG (default: TRUE)
#' @param export_tables Export tables as TSV (default: TRUE)
#' @param plot_width Plot width in inches (default: 12)
#' @param plot_height Plot height in inches (default: 8)
#' @param plot_dpi Plot resolution (default: 150)
#'
#' @return Named list with:
#'   \describe{
#'     \item{opdea_combined}{Combined input data.frame}
#'     \item{mean_aggregated}{Mean of metrics per Assay}
#'     \item{median_aggregated}{Median of metrics per Assay}
#'     \item{mean_ranking}{Ranking table (mean-based)}
#'     \item{median_ranking}{Ranking table (median-based)}
#'     \item{gg_ranking_heatmap_mean}{Heatmap of ranks (mean)}
#'     \item{gg_ranking_heatmap_median}{Heatmap of ranks (median)}
#'     \item{gg_ranking_bars_mean}{Bar chart of final rank (mean)}
#'     \item{gg_ranking_bars_median}{Bar chart of final rank (median)}
#'     \item{gg_metrics_heatmap_mean}{Heatmap of metric values (mean)}
#'     \item{gg_metrics_heatmap_median}{Heatmap of metric values (median)}
#'     \item{gg_metrics_comparison}{Boxplot distributions}
#'     \item{ranking_by_comparison}{Named list of ranking data.frames per comparison}
#'     \item{ranking_by_comparison_combined}{data.frame with all per-comparison rankings}
#'     \item{gg_ranking_heatmap_by_comp}{Named list of heatmaps per comparison}
#'     \item{gg_ranking_bars_by_comp}{Named list of bar charts per comparison}
#'     \item{confusion_combined}{Combined confusion data.frame (if available)}
#'     \item{gg_confusion_stacked}{Faceted confusion stacked bars (if available)}
#'     \item{gg_confusion_stacked_by_comp}{Named list of confusion plots per comparison}
#'     \item{classified_combined}{Combined classified data.frame (if available)}
#'     \item{gg_roc_by_comp}{Named list of ROC plots per comparison (if available)}
#'     \item{gg_roc_zoom_by_comp}{Named list of ROC zoom plots per comparison (if available)}
#'     \item{extended_combined}{Merged classification + OpDEA data (if bench_metrics available)}
#'     \item{extended_mean_aggregated}{Mean-aggregated extended metrics per Assay}
#'     \item{extended_ranking}{Extended ranking table (11 metrics, mean-based)}
#'     \item{gg_extended_ranking_bars}{Bar chart of extended ranking}
#'     \item{gg_extended_ranking_heatmap}{Heatmap of extended ranking}
#'     \item{parameters}{List of parameters used}
#'   }
#'
#' @examples
#' data(nadia_dia)
#' sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
#'                                     package = "NADIA"))
#' ev <- data.frame(Comparison = rep(c("B-A", "D-A"), each = 2),
#'                  Species = c("ECOLI", "YEAST"),
#'                  expected_logFC = c(1, -0.58, 2, -3.3))
#' root <- file.path(tempdir(), "nadia_bench")
#' for (nm in c("cycloess", "quantile")) {
#'   de <- process_proteomics(nadia_dia, norm_method = nm,
#'                            verbose = FALSE)$DEPs_results
#'   de$Species <- sp$PG.OrganismId[match(de$Protein.IDs, sp$PG.ProteinGroups)]
#'   benchmarking_proteomics(de, ev, output_dir = file.path(root, nm),
#'                           verbose = FALSE)
#' }
#' result <- benchmarking_multiple(results_dir = root, verbose = FALSE)
#' result$mean_ranking
#' result$gg_ranking_bars_mean
#' @export
benchmarking_multiple <- function(opdea_combined            = NULL,
                                  confusion_combined        = NULL,
                                  classified_combined       = NULL,
                                  bench_metrics_combined    = NULL,
                                  results_dir               = NULL,
                                  pattern                   = "benchmark_opdea_metrics\\.tsv$",
                                  method_names              = NULL,
                                  recursive                 = TRUE,
                                  strip_prefix              = "benchmark_",
                                  metrics                   = .BM_METRICS,
                                  extended_metrics          = .BM_EXTENDED_METRICS,
                                  p_col                     = "adj.P.Val",
                                  plots                     = "all",
                                  verbose        = TRUE,
                                  output_dir     = NULL,
                                  export_plots   = TRUE,
                                  export_tables  = TRUE,
                                  plot_width     = 12,
                                  plot_height    = 8,
                                  plot_dpi       = 150) {

  # Each step is a helper below the section header: gather the inputs, rank
  # the methods, draw the plots, export them, and put the result together.
  inputs <- .bm_gather_inputs(opdea_combined, confusion_combined,
                              classified_combined, bench_metrics_combined,
                              results_dir, pattern, method_names, recursive,
                              strip_prefix, p_col, verbose)
  ranks <- .bm_rankings(inputs, metrics, extended_metrics, verbose)
  figs  <- .bm_draw(inputs, ranks, metrics, plots, p_col, verbose)
  if (!is.null(output_dir)) {
    .bm_export_all(output_dir, inputs, ranks, figs, export_tables, export_plots,
                   plot_width, plot_height, plot_dpi, verbose)
  }
  result <- .bm_result(inputs, ranks, figs, metrics)

  if (verbose) message("=== DONE ===")
  result
}


# =============================================================================
# SECTION 6: EXAMPLE WORKFLOW
# =============================================================================

# --- Option A: In-memory loop ---
#
#
# library(dplyr)
#
# # Define method combinations to test
# combos <- list(
#   list(norm = "cycloess", imp = "combo", mar = "Impseqrob", mnar = "min"),
#   list(norm = "quantile", imp = "combo", mar = "knn",       mnar = "min"),
#   list(norm = "cycloess", imp = "none"),
#   list(norm = "log2Norm", imp = "combo", mar = "Impseqrob", mnar = "MinDet")
# )
#
# # Helper: build assay name from combo
# .make_assay_name <- function(combo) {
#   if (combo$imp == "none") {
#     paste0(combo$norm, "_none")
#   } else {
#     paste(combo$norm, combo$mar, combo$mnar, sep = "_")
#   }
# }
#
# # Run processing + benchmarking for each combo
# opdea_list       <- list()
# confusion_list   <- list()
# classified_list  <- list()
# bench_met_list   <- list()
# for (combo in combos) {
#   # Build process_proteomics args (omit mar/mnar when imp = "none")
#   proc_args <- list(
#     preprocessing = preprocessing,
#     norm_method   = combo$norm,
#     imp_method    = combo$imp
#   )
#   if (combo$imp != "none") {
#     proc_args$mar_method  <- combo$mar
#     proc_args$mnar_method <- combo$mnar
#   }
#   result <- do.call(process_proteomics, proc_args)
#
#   assay_name <- .make_assay_name(combo)
#
#   bench <- benchmarking_proteomics(
#     de_res          = result$DEPs_results,
#     species_df      = species_df,
#     expected_values = expected,
#     alpha           = 0.05,
#     output_dir      = file.path(out_dir, paste0("benchmark_", assay_name))
#   )
#
#   opdea <- bench$opdea_metrics
#   opdea$Assay <- assay_name
#   opdea_list[[length(opdea_list) + 1]] <- opdea
#
#   conf <- bench$confusion_overall
#   conf$Assay <- assay_name
#   confusion_list[[length(confusion_list) + 1]] <- conf
#
#   cls <- bench$classified_df
#   cls$Assay <- assay_name
#   classified_list[[length(classified_list) + 1]] <- cls
#
#   bmet <- bench$metrics_table
#   bmet$Assay <- assay_name
#   bench_met_list[[length(bench_met_list) + 1]] <- bmet
# }
#
# opdea_all      <- do.call(rbind, opdea_list)
# confusion_all  <- do.call(rbind, confusion_list)
# classified_all <- do.call(rbind, classified_list)
# bench_met_all  <- do.call(rbind, bench_met_list)
#
# bm_result <- benchmarking_multiple(
#   opdea_combined         = opdea_all,
#   confusion_combined     = confusion_all,
#   classified_combined    = classified_all,
#   bench_metrics_combined = bench_met_all,
#   output_dir             = file.path(out_dir, "bm_multiple"),
#   verbose                = TRUE
# )
#
# bm_result$mean_ranking                    # OpDEA ranking (5 metrics)
# bm_result$extended_ranking                # Extended ranking (11 metrics)
# bm_result$gg_extended_ranking_bars        # Bar chart of extended ranking
# bm_result$gg_extended_ranking_heatmap     # Heatmap of extended ranking
# bm_result$gg_ranking_heatmap_mean
# bm_result$gg_roc_by_comp[["B-A"]]
# bm_result$gg_roc_zoom_by_comp[["B-A"]]
#
#
# --- Option B: Import from files ---
#
#
# # Assumes folders like:
# #   <out_dir>/benchmark_cycloess_Impseqrob_min/benchmark_opdea_metrics.tsv
# #   <out_dir>/benchmark_cycloess_Impseqrob_min/benchmark_metrics.tsv        (for extended ranking)
# #   <out_dir>/benchmark_cycloess_Impseqrob_min/benchmark_confusion_overall.tsv
# #   <out_dir>/benchmark_cycloess_Impseqrob_min/benchmark_classified.tsv
# #   <out_dir>/benchmark_quantile_knn_min/...
# #   <out_dir>/benchmark_log2Norm_Impseqrob_MinDet/...
#
# bm_result <- benchmarking_multiple(
#   results_dir = out_dir,
#   output_dir  = file.path(out_dir, "bm_multiple"),
#   verbose     = TRUE
# )
#
# bm_result$mean_ranking                    # OpDEA ranking (5 metrics)
# bm_result$extended_ranking                # Extended ranking (11 metrics)
# bm_result$gg_extended_ranking_bars        # Bar chart of extended ranking
# bm_result$gg_extended_ranking_heatmap     # Heatmap of extended ranking
# bm_result$gg_ranking_bars_mean
# bm_result$gg_roc_by_comp[["B-A"]]
# bm_result$gg_roc_zoom_by_comp[["B-A"]]
