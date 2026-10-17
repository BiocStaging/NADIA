# =============================================================================
# tidyHeatmap heatmaps for proteomics data
# =============================================================================



# -----------------------------------------------------------------------------
# Wrapper and print method for heatmaps carrying a custom title
# -----------------------------------------------------------------------------

#' Wrap a Heatmap Together With Its Title
#'
#' @param hm tidyHeatmap object (InputHeatmap)
#' @param title Heatmap title
#' @param title_size Font size of the title
#' @param title_face Font face of the title
#'
#' @return proteomics_heatmap object (S3 list)
wrap_heatmap_with_title <- function(hm, title, title_size = 14, title_face = "bold") {

  structure(
    list(
      heatmap = hm,
      title = title,
      title_size = title_size,
      title_face = title_face
    ),
    class = "proteomics_heatmap"
  )
}

#' Draw a proteomics heatmap
#'
#' Renders the object built by [proteomics_heatmap()] on the current graphics
#' device, adding the title stored in the object. Printing is what actually draws
#' the heatmap: building it does not.
#'
#' @param x A `proteomics_heatmap` object, as returned by [proteomics_heatmap()].
#' @param ... Further arguments passed to `ComplexHeatmap::draw()`.
#' @return The drawn `HeatmapList`, invisibly. Called for its side effect.
#'
#' @examples
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("tidyHeatmap", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'   hm <- proteomics_heatmap(res$PCA_Input, mode = "any")
#'
#'   # Drawing needs a graphics device; send it to a temporary file
#'   f <- tempfile(fileext = ".png")
#'   grDevices::png(f)
#'   print(hm)
#'   grDevices::dev.off()
#'   file.remove(f)
#' }
#'
#' @export
print.proteomics_heatmap <- function(x, ...) {
  .check_suggests(c("ComplexHeatmap", "tidyHeatmap"), "Drawing a heatmap")
  hm <- x$heatmap

  old_ht_message <- ComplexHeatmap::ht_opt("message")
  ComplexHeatmap::ht_opt(message = FALSE)
  on.exit(ComplexHeatmap::ht_opt(message = old_ht_message), add = TRUE)

  # Convert the InputHeatmap into a ComplexHeatmap using tidyHeatmap's method
  if (inherits(hm, "InputHeatmap")) {
    # Either use as.list to pull out the components and rebuild the object,
    # or simply convert it with tidyHeatmap's own internal method
    ht <- tryCatch({
      # Try the slot directly (older versions)
      methods::slot(hm, "ht")
    }, error = function(e) {
      tryCatch({
        # Try input_heatmap (newer versions)
        methods::slot(hm, "input_heatmap")
      }, error = function(e2) {
        # Fall back on tidyHeatmap's conversion to ComplexHeatmap
        tidyHeatmap::as_ComplexHeatmap(hm)
      })
    })
  } else {
    ht <- hm
  }

  # Draw with the title
  ComplexHeatmap::draw(
    ht,
    column_title = x$title,
    column_title_gp = grid::gpar(fontsize = x$title_size, fontface = x$title_face),
    ...
  )

  invisible(x)
}



# -----------------------------------------------------------------------------
# Helper that builds the colour palette for the heatmap
# -----------------------------------------------------------------------------

#' Build the Colour Palette for the Heatmap Values
#'
#' @param palette Palette for the values, used as a colour ramp: `NULL` for
#'   the default (blue-white-red); `"pkg::name"` for a palette from paletteer;
#'   `"brewer:Name"`, or the bare name, for an RColorBrewer palette (up to `n`
#'   colours); or a vector of colours.
#' @param n Number of colours to generate
#' @param reverse Reverse the palette (default: FALSE)
#'
#' @return Vector of colours, or a colorRamp2 function
get_heatmap_palette <- function(palette = NULL,
                                n = 11,
                                reverse = FALSE) {
  # Default palette (diverging blue-white-red), used as a colour ramp
  default <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
  .resolve_palette(palette, n = n, default = default, type = "continuous",
                   reverse = reverse)
}


#' Build the Colour Palette for Categorical Annotations
#'
#' @param levels Vector of levels (categories)
#' @param palette Colour palette: `NULL` for the default; `"pkg::name"` for a
#'   palette from paletteer; `"brewer:Name"`, or the bare name, for an
#'   RColorBrewer palette; or a vector of colours, unnamed or named after the groups. A palette with
#'   fewer colours than levels is interpolated.
#'
#' @return Named vector of colours
get_annotation_palette <- function(levels, palette = NULL) {
  .resolve_palette(palette, n = length(levels), levels = levels,
                   default = function(n) grDevices::hcl.colors(n, "Dark 3"))
}


# -----------------------------------------------------------------------------
# Main function: prepare the data for the heatmap
# -----------------------------------------------------------------------------

#' Prepare Long-Format Data for tidyHeatmap
#'
#' @param data Data frame in long format with columns:
#'   - SampleID: sample identifier
#'   - FeatureID: protein/feature identifier
#'   - Intensity: intensity value (log2)
#'   - Condition: experimental condition
#'   - Replicate: replicate number
#'   - sig_any: logical flagging significance (for mode="any")
#'   - adjP_*: adjusted p-value columns (for mode="target")
#' @param mode Filtering mode: "all", "any", or "target"
#' @param alpha Significance threshold for mode "target" (default: 0.05)
#' @param comparison Name of the comparison for mode "target" (e.g. "B-A")
#' @param feature_ids Vector of specific FeatureIDs (when supplied, mode/alpha are ignored)
#' @param scale_data Scaling type: "none", "row", "column" (default: "row")
#' @param sample_order Sample order: "clustering", "condition", or a custom vector
#' @param condition_order Condition order when sample_order = "condition"
#'
#' @return Data frame in long format, ready for tidyHeatmap
prepare_heatmap_data <- function(data,
                                 mode = c("all", "any", "target"),
                                 alpha = 0.05,
                                 comparison = NULL,
                                 feature_ids = NULL,
                                 scale_data = c("row", "none", "column"),
                                 sample_order = c("clustering", "condition", "custom"),
                                 condition_order = NULL) {

  mode <- match.arg(mode)
  scale_data <- match.arg(scale_data)

  if (is.character(sample_order) && length(sample_order) == 1) {
    sample_order <- match.arg(sample_order, c("clustering", "condition", "custom"))
  }

  # Validate the required columns
  required_cols <- c("SampleID", "FeatureID", "Intensity", "Condition", "Replicate")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  data <- as.data.frame(data)

  # For mode = "target", keep only the samples belonging to the compared conditions
  if (mode == "target" && !is.null(comparison)) {
    # Extract the conditions from the comparison (e.g. "B-A" -> c("B", "A"))
    conds <- unique(trimws(strsplit(comparison, "[-|:]")[[1]]))
    if (length(conds) >= 2) {
      data <- data[data$Condition %in% conds, , drop = FALSE]
      # Update condition_order so it only lists the relevant conditions
      if (!is.null(condition_order)) {
        condition_order <- condition_order[condition_order %in% conds]
      } else {
        condition_order <- conds
      }
    }
  }

  # Get the feature IDs according to the mode, or use the ones supplied
  if (!is.null(feature_ids) && length(feature_ids) > 0) {
    # Use the supplied IDs directly
    available_ids <- unique(data$FeatureID)
    ids <- feature_ids[feature_ids %in% available_ids]
    if (length(ids) < length(feature_ids)) {
      missing <- setdiff(feature_ids, available_ids)
      warning("FeatureIDs not found in the data (ignored): ",
              paste(head(missing, 5), collapse = ", "),
              if (length(missing) > 5) paste0(" ... and ", length(missing) - 5, " more"))
    }
  } else {
    # Use the mode-based filtering
    ids <- .get_feature_ids(data, mode = mode, alpha = alpha, comparison = comparison)
  }

  if (length(ids) < 2) {
    stop("Subset does not have enough proteins (minimum 2). Found: ", length(ids))
  }

  # Filter the data
  dt <- data[data$FeatureID %in% ids, , drop = FALSE]
  dt <- dt[is.finite(dt$Intensity) & !is.na(dt$Intensity), , drop = FALSE]


  # Build the long-format data frame for the heatmap
  # Include adjP when mode = "target", for the row annotation
  if (mode == "target" && !is.null(comparison)) {
    adjp_colname <- .adjp_col(comparison)
    if (adjp_colname %in% names(dt)) {
      hm_data <- dt %>%
        select(SampleID, FeatureID, Intensity, Condition, Replicate, all_of(adjp_colname)) %>%
        rename(adjP = all_of(adjp_colname)) %>%
        distinct()
    } else {
      hm_data <- dt %>%
        select(SampleID, FeatureID, Intensity, Condition, Replicate) %>%
        distinct()
    }
  } else {
    hm_data <- dt %>%
      select(SampleID, FeatureID, Intensity, Condition, Replicate) %>%
      distinct()
  }

  # Order the samples according to the requested criterion
  if (is.character(sample_order) && length(sample_order) == 1 && sample_order == "condition") {
    # Sort by Condition and then by Replicate
    if (!is.null(condition_order)) {
      hm_data <- hm_data %>%
        mutate(Condition = factor(Condition, levels = condition_order))
    }
    hm_data <- hm_data %>%
      arrange(Condition, Replicate) %>%
      mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
  } else if (is.character(sample_order) && length(sample_order) > 1) {
    # Custom order (vector of SampleIDs)
    # Keep only the specified samples
    samples_in_data <- unique(hm_data$SampleID)
    valid_samples <- sample_order[sample_order %in% samples_in_data]

    if (length(valid_samples) == 0) {
      stop("None of the samples given in sample_order is present in the data")
    }

    if (!all(sample_order %in% samples_in_data)) {
      missing <- setdiff(sample_order, samples_in_data)
      warning("Samples not found in the data (ignored): ", paste(missing, collapse = ", "))
    }

    # Subset the data down to the specified samples
    hm_data <- hm_data %>%
      filter(SampleID %in% valid_samples) %>%
      mutate(SampleID = factor(SampleID, levels = valid_samples))
  }
  # When sample_order == "clustering" we let tidyHeatmap do the clustering

  # Apply the scaling if requested
  if (scale_data != "none") {
    # Pivot to a matrix in order to scale (use mean to collapse possible duplicates)
    mat_wide <- hm_data %>%
      select(SampleID, FeatureID, Intensity) %>%
      pivot_wider(names_from = SampleID, values_from = Intensity, values_fn = mean) %>%
      as.data.frame()

    rownames_feat <- mat_wide$FeatureID
    mat <- as.matrix(mat_wide[, -1])
    rownames(mat) <- rownames_feat

    # NA-aware scaling: base scale() does NOT accept na.rm, so a single NA cell
    # would turn the whole row/column into NA. Centre/scale by hand instead.
    if (scale_data == "row") {
      ctr <- rowMeans(mat, na.rm = TRUE)
      sdv <- apply(mat, 1, sd, na.rm = TRUE)
      # Constant features (sd 0/NA) -> undefined z-score, and they would break
      # the row clustering (all-NA row). Drop them before scaling.
      keep <- is.finite(sdv) & sdv > 0
      if (any(!keep)) {
        message(sprintf(
          "Heatmap: %d features with 0/NA variance dropped before the z-score.",
          sum(!keep)))
        mat <- mat[keep, , drop = FALSE]; ctr <- ctr[keep]; sdv <- sdv[keep]
      }
      mat <- sweep(mat, 1, ctr, "-")
      mat <- sweep(mat, 1, sdv, "/")
    } else if (scale_data == "column") {
      ctr <- colMeans(mat, na.rm = TRUE)
      sdv <- apply(mat, 2, sd, na.rm = TRUE)
      sdv[!is.finite(sdv) | sdv == 0] <- NA
      mat <- sweep(mat, 2, ctr, "-")
      mat <- sweep(mat, 2, sdv, "/")
    }

    # Replace NaN with NA
    mat[is.nan(mat)] <- NA

    # Back to long format
    mat_df <- as.data.frame(mat)
    mat_df$FeatureID <- rownames(mat)

    scaled_long <- mat_df %>%
      pivot_longer(cols = -FeatureID, names_to = "SampleID", values_to = "Intensity")

    # Merge the metadata back in (including adjP when present)
    if ("adjP" %in% names(hm_data)) {
      metadata <- hm_data %>%
        select(SampleID, FeatureID, Condition, Replicate, adjP) %>%
        distinct()

      hm_data <- scaled_long %>%
        left_join(metadata, by = c("SampleID", "FeatureID"))
    } else {
      metadata <- hm_data %>%
        select(SampleID, Condition, Replicate) %>%
        distinct()

      hm_data <- scaled_long %>%
        left_join(metadata, by = "SampleID")
    }

    # Restore the factor ordering where applicable
    if (is.character(sample_order) && length(sample_order) == 1 && sample_order == "condition") {
      if (!is.null(condition_order)) {
        hm_data <- hm_data %>%
          mutate(Condition = factor(Condition, levels = condition_order))
      }
      hm_data <- hm_data %>%
        arrange(Condition, as.numeric(Replicate)) %>%
        mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
    } else if (is.character(sample_order) && length(sample_order) > 1) {
      # Restore the custom order
      hm_data <- hm_data %>%
        mutate(SampleID = factor(SampleID, levels = sample_order))
    }
  }

  hm_data
}


# -----------------------------------------------------------------------------
# Steps of proteomics_heatmap()
#
# The main function prepares the data, then calls these in order: export the
# matrix, read the row annotations, order and split the rows, build the colour
# scale, add the annotation tiles, and export the figure. Each is kept apart so
# it can be tested and reused on its own.
# -----------------------------------------------------------------------------

#' Write the matrix a heatmap draws, one row per protein
#' @param hm_data Output of `prepare_heatmap_data()`.
#' @param export_path Path of the TSV file.
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.hm_export_matrix <- function(hm_data, export_path) {
  # Build a wide matrix with FeatureID as rows and SampleID as columns
  export_wide <- hm_data %>%
    select(FeatureID, SampleID, Intensity) %>%
    tidyr::pivot_wider(
      names_from = SampleID,
      values_from = Intensity,
      values_fn = mean
    )

  # Add adjP when present
  if ("adjP" %in% names(hm_data)) {
    adjp_data <- hm_data %>%
      select(FeatureID, adjP) %>%
      distinct()
    export_wide <- export_wide %>%
      left_join(adjp_data, by = "FeatureID")
  }

  .write_tsv(export_wide, export_path)
  message("Data exported to: ", export_path)
  invisible(NULL)
}

#' Read, check and colour the row annotations of a heatmap
#'
#' Empty values become the level "NA", drawn white and placed last.
#'
#' @param row_annotation,row_annotation_cols,row_annotation_palette See
#'   [proteomics_heatmap()].
#' @param feature_ids The proteins the heatmap shows.
#' @return `NULL` when there is no usable annotation; otherwise a list with
#'   `data` (one factor column per annotation), `cols` and `colors` (one named
#'   palette per annotation, in the order of the factor levels).
#' @keywords internal
#' @noRd
.hm_row_annotation <- function(row_annotation, row_annotation_cols,
                               row_annotation_palette, feature_ids) {
  if (is.null(row_annotation)) return(NULL)

  # Load the annotations when a file path is given
  if (is.character(row_annotation) && length(row_annotation) == 1 && file.exists(row_annotation)) {
    annot <- .read_delim_file(row_annotation)
  } else if (is.data.frame(row_annotation)) {
    annot <- as.data.frame(row_annotation)
  } else {
    warning("row_annotation must be a path to a TSV file or a data.frame",
            call. = FALSE)
    return(NULL)
  }

  if (!("FeatureID" %in% names(annot))) {
    stop("The annotation file must have a 'FeatureID' column", call. = FALSE)
  }

  # Work out which columns to use
  all_annot_cols <- setdiff(names(annot), "FeatureID")
  if (is.null(row_annotation_cols)) {
    row_annotation_cols <- all_annot_cols
  } else {
    missing_cols <- setdiff(row_annotation_cols, all_annot_cols)
    if (length(missing_cols) > 0) {
      warning("Annotation columns not found: ", paste(missing_cols, collapse = ", "),
              call. = FALSE)
      row_annotation_cols <- intersect(row_annotation_cols, all_annot_cols)
    }
  }

  # Turn empty values into the "NA" string (not NA_character_) so they can be
  # coloured white
  for (col in row_annotation_cols) {
    values <- annot[[col]]
    if (is.character(values)) {
      annot[[col]] <- ifelse(is.na(values) | trimws(values) == "", "NA", values)
    }
  }

  # Keep only the FeatureIDs present in the data
  annot <- annot[annot$FeatureID %in% feature_ids, , drop = FALSE]

  # One palette per annotation, in the EXACT order of the factor levels (real
  # levels sorted, then "NA" in white): this is what makes tidyHeatmap assign
  # the colours correctly
  colors <- list()
  for (col in row_annotation_cols) {
    all_levels <- unique(annot[[col]])
    real_levels <- sort(all_levels[all_levels != "NA"])
    ordered_levels <- if ("NA" %in% all_levels) c(real_levels, "NA") else real_levels
    annot[[col]] <- factor(annot[[col]], levels = ordered_levels)

    col_palette <- if (!is.null(row_annotation_palette) &&
                       col %in% names(row_annotation_palette))
      row_annotation_palette[[col]] else NULL
    real_colors <- get_annotation_palette(real_levels, col_palette)
    final_palette <- stats::setNames(character(length(ordered_levels)), ordered_levels)
    for (lvl in ordered_levels) {
      final_palette[lvl] <- if (lvl == "NA") "#FFFFFF" else real_colors[lvl]
    }
    colors[[col]] <- final_palette
  }

  list(data = annot, cols = row_annotation_cols, colors = colors)
}

#' Join the row annotations, then order and split the rows
#' @param hm_data Output of `prepare_heatmap_data()`.
#' @param annot Output of `.hm_row_annotation()`, or `NULL`.
#' @param row_order_by,split_rows_by See [proteomics_heatmap()].
#' @return List with `hm_data` (FeatureID becomes a factor in the requested
#'   order when there is one) and `row_split` (`NULL` when not requested).
#' @keywords internal
#' @noRd
.hm_arrange_rows <- function(hm_data, annot, row_order_by, split_rows_by) {
  row_split <- NULL
  ordering <- !is.null(row_order_by) && !identical(row_order_by, "clustering")
  set_order <- function(hm_data, feature_order) {
    hm_data %>%
      mutate(FeatureID = factor(FeatureID, levels = feature_order)) %>%
      arrange(FeatureID)
  }

  if (!is.null(annot)) {
    annot_data <- annot$data
    cols <- annot$cols

    hm_data <- hm_data %>%
      left_join(annot_data[, c("FeatureID", cols), drop = FALSE], by = "FeatureID")
    # Make sure the factors survive the join
    for (col in cols) {
      if (col %in% names(hm_data) && col %in% names(annot_data)) {
        hm_data[[col]] <- factor(hm_data[[col]], levels = levels(annot_data[[col]]))
      }
    }

    # Order by annotation(s), and by adjP when it is present
    if (ordering) {
      valid_order_cols <- cols
      if ("adjP" %in% names(hm_data)) valid_order_cols <- c(valid_order_cols, "adjP")
      order_cols <- row_order_by[row_order_by %in% valid_order_cols]

      if (length(order_cols) > 0) {
        order_data <- annot_data
        if ("adjP" %in% order_cols && "adjP" %in% names(hm_data)) {
          adjp_data <- hm_data %>%
            select(FeatureID, adjP) %>%
            distinct()
          order_data <- order_data %>%
            left_join(adjp_data, by = "FeatureID")
        }
        # Sort by several columns sequentially
        order_data_sorted <- order_data %>%
          arrange(!!!rlang::syms(order_cols))
        hm_data <- set_order(hm_data, unique(order_data_sorted$FeatureID))
      }
    }

    # Split the rows (AFTER the ordering), keeping the order of the groups
    if (!is.null(split_rows_by) && split_rows_by %in% cols) {
      if (!is.null(row_order_by)) {
        ordered_features <- levels(hm_data$FeatureID)
      } else {
        ordered_features <- unique(as.character(hm_data$FeatureID))
      }
      feature_to_split <- annot_data %>%
        select(FeatureID, all_of(split_rows_by)) %>%
        distinct()
      feature_to_split <- stats::setNames(
        feature_to_split[[split_rows_by]],
        feature_to_split$FeatureID
      )
      split_values <- feature_to_split[ordered_features]
      split_levels <- unique(split_values)
      split_levels <- split_levels[!is.na(split_levels)]
      row_split <- factor(split_values, levels = split_levels)
    }
  }

  # Order by adjP when there are no row annotations (with annotations it is
  # handled above)
  if (ordering && (is.null(annot) || length(annot$cols) == 0)) {
    if ("adjP" %in% row_order_by && "adjP" %in% names(hm_data)) {
      adjp_order <- hm_data %>%
        select(FeatureID, adjP) %>%
        distinct() %>%
        arrange(adjP)
      hm_data <- set_order(hm_data, adjp_order$FeatureID)
    }
  }

  list(hm_data = hm_data, row_split = row_split)
}

#' The colour scale of the heatmap values
#'
#' Symmetric around zero for scaled data, over the observed range otherwise.
#'
#' @param intensity The values drawn.
#' @param palette_value,reverse_palette,scale_data See [proteomics_heatmap()].
#' @return A `circlize::colorRamp2()` function.
#' @keywords internal
#' @noRd
.hm_value_colours <- function(intensity, palette_value, reverse_palette, scale_data) {
  value_colors <- get_heatmap_palette(palette = palette_value, n = 11,
                                      reverse = reverse_palette)
  val_range <- range(intensity, na.rm = TRUE)
  if (scale_data != "none") {
    abs_max <- max(abs(val_range), na.rm = TRUE)
    val_breaks <- seq(-abs_max, abs_max, length.out = length(value_colors))
  } else {
    val_breaks <- seq(val_range[1], val_range[2], length.out = length(value_colors))
  }
  circlize::colorRamp2(val_breaks, value_colors)
}

#' The heatmap data as a tibble of plain base R columns
#'
#' tidyHeatmap tests `class(x) %in% ...`, which fails when a column carries
#' several classes (as arrow/Parquet columns do).
#'
#' @param hm_data The heatmap data.
#' @return A tibble.
#' @keywords internal
#' @noRd
.hm_tidy_input <- function(hm_data) {
  hm_data <- as.data.frame(lapply(hm_data, function(col) {
    if (is.factor(col)) return(factor(as.character(col), levels = levels(col)))
    if (is.logical(col)) return(as.logical(col))
    if (is.numeric(col)) return(as.numeric(col))
    if (is.character(col)) return(as.character(col))
    col
  }), stringsAsFactors = FALSE)
  tibble::as_tibble(hm_data)
}

#' Split the columns by condition, in the order the samples appear
#' @param hm_data The heatmap data.
#' @return A factor with one value per sample.
#' @keywords internal
#' @noRd
.hm_condition_split <- function(hm_data) {
  sample_info <- hm_data %>%
    select(SampleID, Condition) %>%
    distinct()
  sample_order_vec <- unique(as.character(hm_data$SampleID))
  sample_info <- sample_info[match(sample_order_vec, as.character(sample_info$SampleID)), ]
  factor(sample_info$Condition, levels = unique(sample_info$Condition))
}

#' Add the annotation tiles: condition, row annotations and adjusted p-value
#' @param hm A tidyHeatmap heatmap.
#' @param hm_data The heatmap data.
#' @param condition_colors Palette of the conditions, or `NULL` for no tile.
#' @param annot Output of `.hm_row_annotation()`, or `NULL`.
#' @param adjp_colors Three colours for the adjusted p-value, or `NULL` for no
#'   tile.
#' @param show_legend,row_annotation_size,row_annotation_name_size See
#'   [proteomics_heatmap()].
#' @return The heatmap.
#' @keywords internal
#' @noRd
.hm_add_annotations <- function(hm, hm_data, condition_colors, annot,
                                adjp_colors, show_legend, row_annotation_size,
                                row_annotation_name_size) {
  if (!is.null(condition_colors)) {
    hm <- hm %>%
      tidyHeatmap::annotation_tile(Condition, palette = condition_colors,
                                   show_legend = show_legend)
  }

  # Custom row annotations
  if (!is.null(annot) && length(annot$cols) > 0) {
    annot_size <- row_annotation_size
    if (is.numeric(annot_size)) annot_size <- grid::unit(annot_size, "cm")
    annot_name_gp <- grid::gpar(fontsize = row_annotation_name_size)

    valid_cols <- annot$cols[annot$cols %in% names(hm_data)]
    for (i in seq_along(valid_cols)) {
      col <- valid_cols[i]
      hm <- hm %>%
        tidyHeatmap::annotation_tile(
          !!rlang::sym(col),
          palette = annot$colors[[col]],
          # Only the last annotation takes a size (avoids a tidyHeatmap warning)
          size = if (i == length(valid_cols)) annot_size else NULL,
          annotation_name_gp = annot_name_gp,
          show_legend = show_legend
        )
    }
  }

  # Adjusted p-value (low values = more significant)
  if (!is.null(adjp_colors)) {
    hm <- hm %>%
      tidyHeatmap::annotation_tile(
        adjP,
        palette = circlize::colorRamp2(c(0, 0.01, 0.05), adjp_colors),
        show_legend = show_legend
      )
  }
  hm
}

#' Draw a heatmap of any of the types proteomics_heatmap() returns
#' @param hm_obj The heatmap.
#' @return Called for its side effect.
#' @keywords internal
#' @noRd
.hm_draw <- function(hm_obj) {
  if (inherits(hm_obj, "proteomics_heatmap")) {
    # Object carrying a title - the print method already handles everything
    print(hm_obj)
  } else if (inherits(hm_obj, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(hm_obj)
  } else if (inherits(hm_obj, "InputHeatmap")) {
    ComplexHeatmap::draw(tidyHeatmap::as_ComplexHeatmap(hm_obj))
  } else {
    print(hm_obj)
  }
}

#' Export a heatmap to PNG, SVG or PDF, chosen by the file extension
#'
#' A failure is a warning, and the graphics device is closed.
#'
#' @param hm The heatmap.
#' @param export_file,plot_width,plot_height,export_dpi See
#'   [proteomics_heatmap()].
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.hm_export_figure <- function(hm, export_file, plot_width, plot_height, export_dpi) {
  file_ext <- tolower(tools::file_ext(export_file))

  export_dir <- dirname(export_file)
  if (!dir.exists(export_dir) && export_dir != ".") {
    dir.create(export_dir, recursive = TRUE)
  }

  # Pixel size for PNG output (inches x dpi)
  width_pixels <- round(plot_width * export_dpi)
  height_pixels <- round(plot_height * export_dpi)

  tryCatch({
    if (file_ext == "png") {
      grDevices::png(filename = export_file, width = width_pixels,
                     height = height_pixels, res = export_dpi)
      .hm_draw(hm)
      grDevices::dev.off()
      message("Heatmap exported to PNG: ", export_file,
              " (", width_pixels, "x", height_pixels, "px)")
    } else if (file_ext == "svg") {
      grDevices::svg(filename = export_file, width = plot_width, height = plot_height)
      .hm_draw(hm)
      grDevices::dev.off()
      message("Heatmap exported to SVG: ", export_file,
              " (", plot_width, "x", plot_height, "in)")
    } else if (file_ext == "pdf") {
      grDevices::pdf(file = export_file, width = plot_width, height = plot_height)
      .hm_draw(hm)
      grDevices::dev.off()
      message("Heatmap exported to PDF: ", export_file,
              " (", plot_width, "x", plot_height, "in)")
    } else {
      warning("Unsupported format: ", file_ext, ". Use .png, .svg or .pdf",
              call. = FALSE)
    }
  }, error = function(e) {
    # Make sure the graphics device is closed even if something fails
    try(grDevices::dev.off(), silent = TRUE)
    warning("Could not export the heatmap: ", e$message, call. = FALSE)
  })
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Main function: heatmap built with tidyHeatmap
# -----------------------------------------------------------------------------

#' Create a tidyHeatmap Heatmap for Proteomics
#'
#' @param data Data frame in long format (see prepare_heatmap_data for the structure)
#' @param mode Protein filtering mode: "all", "any", or "target"
#' @param alpha Significance threshold for DEP proteins (default: 0.05)
#' @param comparison Name of the comparison for mode "target" (e.g. "B-A")
#' @param feature_ids Vector of specific FeatureIDs to show (default: NULL, uses the mode-based filtering).
#'   When supplied, only these proteins are shown and the mode/alpha filtering is ignored.
#' @param row_annotation Row annotations (FeatureIDs). It can be:
#'   - Path to a TSV file with a "FeatureID" column plus additional categorical columns
#'   - Data frame with the same structure
#'   - NULL: no row annotations (default)
#' @param row_annotation_cols Vector of column names to display as annotations.
#'   Default: NULL (uses every column except FeatureID)
#' @param row_annotation_palette Named list of palettes, one per annotation.
#'   Example: list(Pathway = "brewer:Set1", Function = c("red", "blue", "green"))
#' @param row_annotation_size Width of the row annotation bars. It can be:
#'   - Number: interpreted as centimetres (e.g. 0.3 = 0.3cm)
#'   - A unit object: grid::unit(0.3, "cm")
#'   - NULL: use tidyHeatmap's default value
#' @param row_annotation_name_size Font size of the row annotation names (default: 8)
#' @param row_order_by Row order. It can be:
#'   - NULL: no particular ordering (default)
#'   - "clustering": order by hierarchical clustering
#'   - A column name: order by that annotation (e.g. "Specie")
#'   - A vector of columns: order sequentially (e.g. c("Specie", "Process", "Function"))
#'   - Include "adjP" when mode = "target" to order by adjusted p-value
#' @param split_rows_by Name of the annotation column used to split the rows into groups (default: NULL)
#' @param scale_data Scaling type: "none", "row", "column" (default: "row")
#' @param sample_order Sample order: "clustering", "condition", or a custom vector
#' @param condition_order Condition order when sample_order = "condition"
#' @param cluster_rows Cluster the rows (default: FALSE)
#' @param cluster_columns Cluster the columns (default: TRUE)
#' @param show_row_names Show the row names (default: TRUE when <= 50 proteins)
#' @param show_column_names Show the column names (default: TRUE)
#' @param palette_value Palette for the heatmap values (see get_heatmap_palette)
#' @param palette_annotation Palette for the Condition annotation
#' @param reverse_palette Reverse the value palette (default: FALSE)
#' @param row_title Title for the rows
#' @param column_title Title for the columns
#' @param row_title_size Font size of the row title (default: 10)
#' @param column_title_size Font size of the column title (default: 10)
#' @param show_row_title Show the row title (default: TRUE)
#' @param show_column_title Show the column title (default: TRUE)
#' @param show_annotation Show the Condition annotation (default: TRUE)
#' @param split_by_condition Split the heatmap by condition (default: FALSE)
#' @param show_adjp_annotation Show the adjP row annotation when mode="target" (default: TRUE)
#' @param palette_adjp Palette for the adjP annotation, three colours for
#'   adjP = 0, 0.01 and 0.05: `NULL` for the default (red-orange-white),
#'   `"pkg::name"` (paletteer), `"brewer:Name"` (RColorBrewer) or a vector of
#'   colours.
#' @param row_names_size Font size of the row names (default: 7)
#' @param column_names_size Font size of the column names (default: 9)
#' @param column_names_rotation Rotation of the column names in degrees (default: 45)
#' @param column_dend_height Height of the column dendrogram. It can be:
#'   - Number: interpreted as millimetres (e.g. 30 = 30mm)
#'   - A unit object: grid::unit(2, "cm")
#'   - NULL: use ComplexHeatmap's default value
#' @param row_dend_width Width of the row dendrogram. Same format as column_dend_height
#' @param show_heatmap_legend Show the heatmap legend (default: TRUE)
#' @param show_annotation_legend Show the annotation legends (default: TRUE)
#' @param border_color Border colour of the heatmap cells. It can be:
#'   - NULL or FALSE: no border (default)
#'   - TRUE: black border
#'   - A colour string: that specific colour (e.g. "black", "grey", "#CCCCCC")
#' @param heatmap_title Main heatmap title (default: NULL, no title)
#' @param heatmap_title_size Font size of the main title (default: 14)
#' @param heatmap_title_face Font face of the title: "plain", "bold", "italic", "bold.italic" (default: "bold")
#' @param export_path Path used to export the heatmap data to TSV (default: NULL, nothing is exported).
#'   The file contains FeatureID, the per-sample intensity values, and the metadata (adjP where applicable).
#' @param export_file Path used to export the plot. The format is taken from the extension:
#'   - .png: PNG image (raster)
#'   - .svg: SVG image (vector)
#'   - .pdf: PDF document (vector)
#' @param plot_width Plot width in inches (default: 10).
#'   To get pixels: pixels = inches x dpi (e.g. 10" x 300dpi = 3000px)
#' @param plot_height Plot height in inches (default: 8).
#'   To get pixels: pixels = inches x dpi (e.g. 8" x 300dpi = 2400px)
#' @param export_dpi Resolution for PNG output in dots per inch (default: 300).
#'   Higher dpi = more detail. Common values: 72 (web), 150 (draft), 300 (publication)
#'
#' @return tidyHeatmap/ComplexHeatmap object
#'
#' @examples
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("tidyHeatmap", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#'   # Proteins significant in the B-A comparison, samples grouped by condition
#'   hm <- proteomics_heatmap(
#'     res$PCA_Input,
#'     mode            = "target",
#'     comparison      = "B-A",
#'     scale_data      = "row",
#'     sample_order    = "condition",
#'     condition_order = c("A", "B", "D")
#'   )
#'   print(class(hm))   # print(hm) itself draws it on the current device
#' }
#'
#' @export
proteomics_heatmap <- function(data,
                               mode = c("all", "any", "target"),
                               alpha = 0.05,
                               comparison = NULL,
                               feature_ids = NULL,
                               row_annotation = NULL,
                               row_annotation_cols = NULL,
                               row_annotation_palette = NULL,
                               row_annotation_size = NULL,
                               row_annotation_name_size = 8,
                               row_order_by = NULL,
                               split_rows_by = NULL,
                               scale_data = c("row", "none", "column"),
                               sample_order = "clustering",
                               condition_order = NULL,
                               cluster_rows = FALSE,
                               cluster_columns = TRUE,
                               show_row_names = NULL,
                               show_column_names = TRUE,
                               palette_value = NULL,
                               palette_annotation = NULL,
                               reverse_palette = FALSE,
                               row_title = "Proteins",
                               column_title = "Samples",
                               row_title_size = 10,
                               column_title_size = 10,
                               show_row_title = TRUE,
                               show_column_title = TRUE,
                               show_annotation = TRUE,
                               split_by_condition = FALSE,
                               show_adjp_annotation = TRUE,
                               palette_adjp = NULL,
                               row_names_size = 7,
                               column_names_size = 9,
                               column_names_rotation = 45,
                               column_dend_height = NULL,
                               row_dend_width = NULL,
                               show_heatmap_legend = TRUE,
                               show_annotation_legend = TRUE,
                               border_color = NULL,
                               heatmap_title = NULL,
                               heatmap_title_size = 14,
                               heatmap_title_face = c("bold", "plain", "italic", "bold.italic"),
                               export_path = NULL,
                               export_file = NULL,
                               plot_width = 10,
                               plot_height = 8,
                               export_dpi = 300) {
  .check_suggests(c("ComplexHeatmap", "tidyHeatmap", "circlize"),
                  "proteomics_heatmap()")

  mode <- match.arg(mode)
  scale_data <- match.arg(scale_data)
  heatmap_title_face <- match.arg(heatmap_title_face)

  # Silence ComplexHeatmap's messages (use_raster, magick) while this call draws,
  # and give the user back the value they had. The messages appear when drawing,
  # so print.proteomics_heatmap() silences them the same way.
  old_ht_message <- ComplexHeatmap::ht_opt("message")
  ComplexHeatmap::ht_opt(message = FALSE)
  on.exit(ComplexHeatmap::ht_opt(message = old_ht_message), add = TRUE)

  # Resolve border_color
  rect_gp <- NULL
  if (!is.null(border_color) && !isFALSE(border_color)) {
    rect_gp <- grid::gpar(col = if (isTRUE(border_color)) "black" else border_color)
  }

  # Hide the titles when show_*_title is FALSE
  if (!show_row_title) row_title <- NULL
  if (!show_column_title) column_title <- NULL

  # Dendrogram sizes given as numbers are millimetres
  if (is.numeric(column_dend_height)) column_dend_height <- grid::unit(column_dend_height, "mm")
  if (is.numeric(row_dend_width)) row_dend_width <- grid::unit(row_dend_width, "mm")

  # ---------------------------------------------------------------------------
  # 1) Data: select the proteins, scale them, order the samples; export them
  # ---------------------------------------------------------------------------
  hm_data <- prepare_heatmap_data(
    data = data,
    mode = mode,
    alpha = alpha,
    comparison = comparison,
    feature_ids = feature_ids,
    scale_data = scale_data,
    sample_order = sample_order,
    condition_order = condition_order
  )
  if (!is.null(export_path) && nzchar(export_path)) {
    .hm_export_matrix(hm_data, export_path)
  }

  # ---------------------------------------------------------------------------
  # 2) Rows: annotations, order and split
  # ---------------------------------------------------------------------------
  annot <- .hm_row_annotation(row_annotation, row_annotation_cols,
                              row_annotation_palette, unique(hm_data$FeatureID))
  rows <- .hm_arrange_rows(hm_data, annot, row_order_by, split_rows_by)
  hm_data <- rows$hm_data

  # Show the row names only when they can be read
  if (is.null(show_row_names)) {
    show_row_names <- length(unique(hm_data$FeatureID)) <= 50
  }

  # ---------------------------------------------------------------------------
  # 3) Colours and clustering
  # ---------------------------------------------------------------------------
  palette_func <- .hm_value_colours(hm_data$Intensity, palette_value,
                                    reverse_palette, scale_data)
  condition_levels <- unique(hm_data$Condition)
  if (!is.null(condition_order)) {
    condition_levels <- condition_order[condition_order %in% condition_levels]
  }
  annotation_colors <- get_annotation_palette(condition_levels, palette_annotation)

  # A sample order other than "clustering" turns column clustering off
  cluster_cols_final <- cluster_columns
  if (is.character(sample_order) && (length(sample_order) > 1 || sample_order != "clustering")) {
    cluster_cols_final <- FALSE
  }
  # row_order_by = "clustering" turns row clustering on, any other order off;
  # split_rows_by is compatible with clustering (it clusters within each split)
  cluster_rows_final <- cluster_rows
  if (!is.null(row_order_by)) {
    cluster_rows_final <- identical(row_order_by, "clustering")
  }

  # ---------------------------------------------------------------------------
  # 4) The heatmap and its annotations
  # ---------------------------------------------------------------------------
  hm_data <- .hm_tidy_input(hm_data)

  extra_args <- list(
    column_dend_height = column_dend_height,
    row_dend_width     = row_dend_width,
    rect_gp            = rect_gp,
    row_split          = rows$row_split
  )
  extra_args <- extra_args[!vapply(extra_args, is.null, logical(1))]

  hm_args <- c(
    list(
      .data = hm_data,
      .row = rlang::sym("FeatureID"),
      .column = rlang::sym("SampleID"),
      .value = rlang::sym("Intensity"),
      scale = "none",
      cluster_rows = cluster_rows_final,
      cluster_columns = cluster_cols_final,
      palette_value = palette_func,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_gp = grid::gpar(fontsize = row_names_size),
      column_names_gp = grid::gpar(fontsize = column_names_size),
      column_names_rot = column_names_rotation,
      row_title = row_title,
      column_title = column_title,
      row_title_gp = grid::gpar(fontsize = row_title_size),
      column_title_gp = grid::gpar(fontsize = column_title_size),
      column_split = if (split_by_condition) .hm_condition_split(hm_data),
      show_heatmap_legend = show_heatmap_legend
    ),
    extra_args
  )
  hm <- do.call(tidyHeatmap::heatmap, hm_args)

  show_adjp <- show_adjp_annotation && mode == "target" && "adjP" %in% names(hm_data)
  hm <- .hm_add_annotations(
    hm, hm_data,
    condition_colors = if (show_annotation) annotation_colors,
    annot = annot,
    adjp_colors = if (show_adjp) {
      # Default: dark red -> orange -> white
      .resolve_palette(palette_adjp, n = 3,
                       default = c("#67001F", "#F4A582", "#F7F7F7"))
    },
    show_legend = show_annotation_legend,
    row_annotation_size = row_annotation_size,
    row_annotation_name_size = row_annotation_name_size
  )

  # ---------------------------------------------------------------------------
  # 5) Title and export
  # ---------------------------------------------------------------------------
  if (!is.null(heatmap_title) && nzchar(heatmap_title)) {
    hm <- wrap_heatmap_with_title(
      hm = hm,
      title = heatmap_title,
      title_size = heatmap_title_size,
      title_face = heatmap_title_face
    )
  }
  if (!is.null(export_file) && nzchar(export_file)) {
    .hm_export_figure(hm, export_file, plot_width, plot_height, export_dpi)
  }

  hm
}


# -----------------------------------------------------------------------------
# Wrapper function: build a list of heatmaps
# -----------------------------------------------------------------------------

#' Build a List of Heatmaps for Several Subsets
#'
#' Automatically builds heatmaps for "all", "any" and/or specific comparisons.
#'
#' @param data Data frame in long format (see prepare_heatmap_data for the structure)
#' @param modes Vector of modes to generate: "all", "any", and/or comparison names
#'   (default: c("all", "any"))
#' @param alpha Significance threshold for DEP proteins (default: 0.05)
#' @param feature_ids Vector of specific FeatureIDs to show (default: NULL)
#' @param row_annotation Row annotations (TSV path or data.frame)
#' @param row_annotation_cols Columns to use as row annotations
#' @param row_annotation_palette List of palettes for the row annotations
#' @param row_annotation_size Width of the row annotation bars (a number in cm, or a unit)
#' @param row_annotation_name_size Font size of the row annotation names (default: 8)
#' @param row_order_by Row order: "clustering", annotation column(s), or include "adjP" for mode="target"
#' @param split_rows_by Annotation column used to split the rows into groups
#' @param scale_data Scaling type: "none", "row", "column" (default: "row")
#' @param sample_order Sample order: "clustering", "condition", or a custom vector
#' @param condition_order Condition order when sample_order = "condition"
#' @param cluster_rows Cluster the rows (default: FALSE)
#' @param cluster_columns Cluster the columns (default: TRUE)
#' @param show_row_names Show the row names (default: auto)
#' @param show_column_names Show the column names (default: TRUE)
#' @param palette_value Palette for the heatmap values
#' @param palette_annotation Palette for the Condition annotation
#' @param reverse_palette Reverse the value palette (default: FALSE)
#' @param row_title_size Font size of the row title (default: 10)
#' @param column_title_size Font size of the column title (default: 10)
#' @param show_row_title Show the row title (default: TRUE)
#' @param show_column_title Show the column title (default: TRUE)
#' @param show_annotation Show the Condition annotation (default: TRUE)
#' @param split_by_condition Split the heatmap by condition (default: FALSE)
#' @param show_adjp_annotation Show the adjP row annotation for comparisons (default: TRUE)
#' @param palette_adjp Palette for the adjP annotation (see proteomics_heatmap)
#' @param row_names_size Font size of the row names (default: 7)
#' @param column_names_size Font size of the column names (default: 9)
#' @param column_names_rotation Rotation of the column names (default: 45)
#' @param column_dend_height Height of the column dendrogram (a number in mm, or a unit)
#' @param row_dend_width Width of the row dendrogram (a number in mm, or a unit)
#' @param show_heatmap_legend Show the heatmap legend (default: TRUE)
#' @param show_annotation_legend Show the annotation legends (default: TRUE)
#' @param border_color Border colour of the cells (NULL, TRUE, or a colour)
#' @param heatmap_title Main heatmap title (default: NULL, no title).
#'   "\{mode\}" can be used as a placeholder, and is replaced by the mode name
#' @param heatmap_title_size Font size of the main title (default: 14)
#' @param heatmap_title_face Font face of the title (default: "bold")
#' @param export_path Base path used to export the data to TSV (default: NULL).
#'   The mode name is appended to the file name (e.g. "export_all.tsv", "export_B-A.tsv")
#' @param export_modes Vector of modes to export (default: NULL, exports all of them).
#'   Only relevant when export_path is set. Example: c("all", "B-A")
#' @param export_file Base path used to export the plots. The mode is appended to the name.
#'   The format is taken from the extension (.png, .svg, .pdf)
#' @param plot_width Plot width in inches (default: 10)
#' @param plot_height Plot height in inches (default: 8)
#' @param export_dpi Resolution for PNG output in dots per inch (default: 300)
#'
#' @return Named list of tidyHeatmap objects
#'
#' @examples
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("tidyHeatmap", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#'   # "any" plus one heatmap per named comparison
#'   hm_list <- proteomics_heatmap_list(
#'     res$PCA_Input,
#'     modes           = c("any", "B-A"),
#'     sample_order    = "condition",
#'     condition_order = c("A", "B", "D"),
#'     palette_value   = "brewer:RdYlBu"
#'   )
#'   print(names(hm_list))   # print(hm_list[["B-A"]]) draws one of them
#' }
#'
#' @export
proteomics_heatmap_list <- function(data,
                                    modes = c("all", "any"),
                                    alpha = 0.05,
                                    feature_ids = NULL,
                                    row_annotation = NULL,
                                    row_annotation_cols = NULL,
                                    row_annotation_palette = NULL,
                                    row_annotation_size = NULL,
                                    row_annotation_name_size = 8,
                                    row_order_by = NULL,
                                    split_rows_by = NULL,
                                    scale_data = c("row", "none", "column"),
                                    sample_order = "clustering",
                                    condition_order = NULL,
                                    cluster_rows = FALSE,
                                    cluster_columns = TRUE,
                                    show_row_names = NULL,
                                    show_column_names = TRUE,
                                    palette_value = NULL,
                                    palette_annotation = NULL,
                                    reverse_palette = FALSE,
                                    row_title_size = 10,
                                    column_title_size = 10,
                                    show_row_title = TRUE,
                                    show_column_title = TRUE,
                                    show_annotation = TRUE,
                                    split_by_condition = FALSE,
                                    show_adjp_annotation = TRUE,
                                    palette_adjp = NULL,
                                    row_names_size = 7,
                                    column_names_size = 9,
                                    column_names_rotation = 45,
                                    column_dend_height = NULL,
                                    row_dend_width = NULL,
                                    show_heatmap_legend = TRUE,
                                    show_annotation_legend = TRUE,
                                    border_color = NULL,
                                    heatmap_title = NULL,
                                    heatmap_title_size = 14,
                                    heatmap_title_face = "bold",
                                    export_path = NULL,
                                    export_modes = NULL,
                                    export_file = NULL,
                                    plot_width = 10,
                                    plot_height = 8,
                                    export_dpi = 300) {
  .check_suggests(c("ComplexHeatmap", "tidyHeatmap", "circlize"),
                  "proteomics_heatmap_list()")

  scale_data <- match.arg(scale_data)

  # ---------------------------------------------------------------------------
  # 1) Input validation
  # ---------------------------------------------------------------------------

  required_cols <- c("SampleID", "FeatureID", "Intensity", "Condition", "Replicate")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # 2) Detect the available comparisons (adjP_* columns)
  # ---------------------------------------------------------------------------

  adjp_cols <- grep("^adjP_", names(data), value = TRUE)
  available_comparisons <- sub("^adjP_", "", adjp_cols)

  # ---------------------------------------------------------------------------
  # 3) Build a heatmap for each mode
  # ---------------------------------------------------------------------------

  hm_list <- list()

  for (m in modes) {

    # Resolve the internal mode and the comparison (where applicable)
    if (m == "all") {
      internal_mode <- "all"
      comparison <- NULL
      row_title <- "All Proteins"
    } else if (m == "any") {
      internal_mode <- "any"
      comparison <- NULL
      row_title <- "DEPs (any comparison)"
    } else {
      # It is a specific comparison (mode "target")
      if (!(m %in% available_comparisons)) {
        warning("Comparison '", m, "' not found. It will be skipped.")
        next
      }
      internal_mode <- "target"
      comparison <- m
      row_title <- paste0("DEPs (", m, ")")
    }

    # Resolve the title (substitute {mode} when present)
    current_title <- NULL
    if (!is.null(heatmap_title)) {
      current_title <- gsub("\\{mode\\}", m, heatmap_title)
    }

    # Resolve export_path (append the mode to the file name)
    current_export_path <- NULL
    if (!is.null(export_path) && nzchar(export_path)) {
      # Check whether this mode should be exported
      should_export <- is.null(export_modes) || m %in% export_modes
      if (should_export) {
        # Split off the directory, the base name and the extension
        dir_path <- dirname(export_path)
        base_name <- tools::file_path_sans_ext(basename(export_path))
        ext <- tools::file_ext(export_path)
        if (nzchar(ext)) ext <- paste0(".", ext) else ext <- ".tsv"
        current_export_path <- file.path(dir_path, paste0(base_name, "_", m, ext))
      }
    }

    # Resolve export_file (append the mode to the plot file name)
    current_export_file <- NULL
    if (!is.null(export_file) && nzchar(export_file)) {
      # Check whether this mode should be exported
      should_export <- is.null(export_modes) || m %in% export_modes
      if (should_export) {
        # Split off the directory, the base name and the extension
        dir_path <- dirname(export_file)
        base_name <- tools::file_path_sans_ext(basename(export_file))
        ext <- tools::file_ext(export_file)
        if (nzchar(ext)) ext <- paste0(".", ext) else ext <- ".png"
        current_export_file <- file.path(dir_path, paste0(base_name, "_", m, ext))
      }
    }

    # Try to build the heatmap (it can fail when there are not enough proteins)
    hm <- tryCatch({
      proteomics_heatmap(
        data = data,
        mode = internal_mode,
        alpha = alpha,
        comparison = comparison,
        feature_ids = feature_ids,
        row_annotation = row_annotation,
        row_annotation_cols = row_annotation_cols,
        row_annotation_palette = row_annotation_palette,
        row_annotation_size = row_annotation_size,
        row_annotation_name_size = row_annotation_name_size,
        row_order_by = row_order_by,
        split_rows_by = split_rows_by,
        scale_data = scale_data,
        sample_order = sample_order,
        condition_order = condition_order,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        palette_value = palette_value,
        palette_annotation = palette_annotation,
        reverse_palette = reverse_palette,
        row_title = row_title,
        column_title = "Samples",
        row_title_size = row_title_size,
        column_title_size = column_title_size,
        show_row_title = show_row_title,
        show_column_title = show_column_title,
        show_annotation = show_annotation,
        split_by_condition = split_by_condition,
        show_adjp_annotation = show_adjp_annotation,
        palette_adjp = palette_adjp,
        row_names_size = row_names_size,
        column_names_size = column_names_size,
        column_names_rotation = column_names_rotation,
        column_dend_height = column_dend_height,
        row_dend_width = row_dend_width,
        show_heatmap_legend = show_heatmap_legend,
        show_annotation_legend = show_annotation_legend,
        border_color = border_color,
        heatmap_title = current_title,
        heatmap_title_size = heatmap_title_size,
        heatmap_title_face = heatmap_title_face,
        export_path = current_export_path,
        export_file = current_export_file,
        plot_width = plot_width,
        plot_height = plot_height,
        export_dpi = export_dpi
      )
    }, error = function(e) {
      warning("Could not build the heatmap for '", m, "': ", e$message)
      return(NULL)
    })

    if (!is.null(hm)) {
      hm_list[[m]] <- hm
    }
  }

  hm_list
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# --- Load the data ---
# hm_input <- arrow::read_parquet("PCA_Input.parquet")
# hm_input <- read.delim("PCA_Input.tsv", check.names = FALSE)

# --- Basic example: heatmap of every protein ---
# hm_all <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   sample_order = "condition",
#   condition_order = c("A", "B", "C", "D")
# )
# hm_all

# --- Heatmap with column clustering ---
# hm_cluster <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   sample_order = "clustering",
#   cluster_columns = TRUE,
#   cluster_rows = FALSE
# )
# hm_cluster

# --- Heatmap of the proteins significant in any comparison ---
# hm_any <- proteomics_heatmap(
#   data = hm_input,
#   mode = "any",
#   scale_data = "row",
#   sample_order = "condition",
#   condition_order = c("A", "B", "C", "D")
# )
# hm_any

# --- Heatmap of the proteins significant in one specific comparison ---
# hm_target <- proteomics_heatmap(
#   data = hm_input,
#   mode = "target",
#   comparison = "B-A",
#   alpha = 0.05,
#   scale_data = "row"
# )
# hm_target

# --- With a custom paletteer palette ---
# hm_viridis <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   palette_value = "viridis::viridis",
#   reverse_palette = TRUE
# )
# hm_viridis

# --- With an RColorBrewer palette ---
# hm_brewer <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   palette_value = "brewer:RdYlBu",
#   palette_annotation = "brewer:Set1"
# )
# hm_brewer

# --- Build several heatmaps at once ---
# hm_list <- proteomics_heatmap_list(
#   data = hm_input,
#   modes = c("all", "any", "B-A", "C-A", "D-A"),
#   scale_data = "row",
#   sample_order = "condition",
#   condition_order = c("A", "B", "C", "D")
# )
# hm_list[["all"]]
# hm_list[["any"]]
# hm_list[["B-A"]]

# --- Custom sample order ---
# custom_order <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2")
# hm_custom <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   sample_order = custom_order,
#   cluster_columns = FALSE
# )
# hm_custom

# --- Heatmap with a custom title ---
# hm_title <- proteomics_heatmap(
#   data = hm_input,
#   mode = "any",
#   scale_data = "row",
#   heatmap_title = "Differential Expression Analysis",
#   heatmap_title_size = 16,
#   heatmap_title_face = "bold"
# )
# hm_title

# --- List of heatmaps with dynamic titles ---
# hm_list <- proteomics_heatmap_list(
#   data = hm_input,
#   modes = c("all", "any", "B-A"),
#   heatmap_title = "Proteomics Heatmap: {mode}",
#   heatmap_title_size = 14
# )
# # The titles will be: "Proteomics Heatmap: all", "Proteomics Heatmap: any", "Proteomics Heatmap: B-A"

# --- Available palettes ---
# Diverging (good for scaled data):
#   - NULL (default RdBu-like)
#   - "brewer:RdBu", "brewer:RdYlBu", "brewer:PiYG", "brewer:BrBG"
#   - "viridis::plasma", "viridis::inferno"
#
# Sequential (good for unscaled data):
#   - "viridis::viridis", "viridis::magma"
#   - "brewer:Blues", "brewer:Reds", "brewer:YlOrRd"
#
# For annotations:
#   - "brewer:Set1", "brewer:Set2", "brewer:Dark2"
#   - "ggsci::category10_d3", "ggsci::nrc_npg"
