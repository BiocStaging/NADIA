# =============================================================================
# Shared internal utilities
# =============================================================================
#
# Definitions common to every NADIA module. Each module used to redefine `%||%`
# on its own (18 copies with 3 different semantics), so the last definition
# loaded won and silently changed the behaviour of the modules already loaded.
# There is now a single canonical definition.
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================


#' Null-coalescing operator
#'
#' Returns `a` unless it is `NULL`, in which case it returns `b`. Same semantics
#' as `base::\%||\%` (available since R 4.4) and as `rlang::\%||\%`.
#'
#' It is defined here rather than imported from `base` so that the package keeps
#' working under the declared `Depends: R (>= 4.4)` without relying on the
#' operator being exported in that particular version. The definition is
#' identical.
#'
#' @param a Value to check.
#' @param b Fallback value if `a` is `NULL`.
#' @return `a` if it is not `NULL`; otherwise `b`.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


# =============================================================================
# Colour helpers
# =============================================================================
#
# These were duplicated in Boxplot_Highcharts_Final.R, PCA_Highcharts_Final.R,
# Pattern_Profiler_Highcharts.R and Heatmap_tidyHeatmap.R (8 definitions of 3
# functions). With every module living in the global environment, the active
# copy was the one from the last module loaded; in a package the winner would
# be the one from the last file in alphabetical order. They are unified here.
#
# Verified before unifying: for hexadecimal inputs the three copies of
# `hex_to_rgba` and the two of `darken_hex` produced identical output for the
# same alpha/factor. They only differed in the default value, and all 7 call
# sites always pass it explicitly. That is why NO default is declared here: a
# call without the argument must fail visibly instead of taking some arbitrary
# value.


#' Normalise a hexadecimal colour
#'
#' Strips the leading `#`, discards the alpha channel if the colour comes in
#' `RRGGBBAA` format, and returns `#RRGGBB` in upper case.
#'
#' @param hex Colour in hexadecimal format, with or without `#`.
#' @return `#RRGGBB` string in upper case.
#' @keywords internal
#' @noRd
.normalize_hex <- function(hex) {
  hex <- gsub("^#", "", hex)
  if (nchar(hex) == 8) {
    hex <- substr(hex, 1, 6)
  }
  paste0("#", toupper(hex))
}


#' Convert a hexadecimal colour to an `rgba()` string
#'
#' @param hex Hexadecimal colour.
#' @param alpha Opacity between 0 and 1. Deliberately without a default value:
#'   every call site specifies it.
#' @return `"rgba(r, g, b, a)"` string suitable for Highcharts.
#' @keywords internal
#' @noRd
.hex_to_rgba <- function(hex, alpha) {
  hex <- .normalize_hex(hex)
  rgb_vals <- grDevices::col2rgb(hex)
  sprintf("rgba(%d, %d, %d, %.2f)",
          rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha)
}


#' Darken a hexadecimal colour
#'
#' @param hex Hexadecimal colour.
#' @param factor Darkening fraction between 0 and 1. Deliberately without a
#'   default value: every call site specifies it.
#' @return Darkened `#RRGGBB` colour.
#' @keywords internal
#' @noRd
.darken_hex <- function(hex, factor) {
  hex <- .normalize_hex(hex)
  rgb_vals <- grDevices::col2rgb(hex)
  rgb_dark <- pmax(0, rgb_vals * (1 - factor))
  sprintf("#%02X%02X%02X",
          round(rgb_dark[1]), round(rgb_dark[2]), round(rgb_dark[3]))
}


# =============================================================================
# Palettes
# =============================================================================
#
# Every plotting function accepts the same palette specifications and resolves
# them here. Each module used to parse them on its own, with different results:
# some recycled a short palette and others interpolated it, some stopped on an
# unknown name and others fell back to the default with a warning, and the
# volcano plot accepted paletteer names only (Bioconductor review, #186).

#' Resolve a palette specification into colours
#'
#' @param palette `NULL`, for `default`; `"pkg::name"`, a discrete or continuous
#'   palette from paletteer; `"brewer:Name"` or a bare RColorBrewer name; a
#'   single colour; or a vector of colours, unnamed or named after `levels`.
#' @param n Number of colours wanted.
#' @param levels Optional names for the colours. A named palette must cover
#'   all of them.
#' @param default Colours used when `palette` is `NULL`: a vector, or a
#'   function of `n`.
#' @param type `"discrete"` returns exactly `n` colours, interpolating when the
#'   palette has fewer: recycling would give two groups the same colour.
#'   `"continuous"` returns the colours of a colour ramp as the palette gives
#'   them (RColorBrewer: up to `n`).
#' @param reverse Reverse the order of the colours.
#' @return Character vector of colours, named after `levels` when given.
#' @keywords internal
#' @noRd
.resolve_palette <- function(palette, n, levels = NULL, default = NULL,
                             type = c("discrete", "continuous"),
                             reverse = FALSE) {
  type <- match.arg(type)
  spec <- if (is.character(palette) && length(palette) == 1 &&
              is.null(names(palette))) palette else NULL
  brewer_names <- rownames(RColorBrewer::brewer.pal.info)
  is_colour <- function(x) !inherits(tryCatch(grDevices::col2rgb(x),
                                              error = function(e) e), "error")

  if (is.null(palette)) {
    cols <- if (is.function(default)) default(n) else default
  } else if (!is.null(names(palette)) && !is.null(levels)) {
    miss <- setdiff(levels, names(palette))
    if (length(miss) > 0)
      stop("Missing colours for levels: ", paste(miss, collapse = ", "),
           call. = FALSE)
    cols <- palette[levels]
    if (reverse) cols <- stats::setNames(rev(cols), levels)
    return(cols)
  } else if (!is.null(spec) && grepl("::", spec, fixed = TRUE)) {
    pkg <- sub("::.*", "", spec)
    .check_suggests(c("paletteer", pkg), sprintf("The palette '%s'", spec))
    cols <- tryCatch(
      as.character(paletteer::paletteer_d(spec)),
      error = function(e) tryCatch(
        as.character(paletteer::paletteer_c(spec, n = n)),
        error = function(e2) stop("Unknown palette '", spec, "': ",
                                  conditionMessage(e2), call. = FALSE)))
    cols <- vapply(cols, .normalize_hex, character(1), USE.NAMES = FALSE)
  } else if (!is.null(spec) && (startsWith(spec, "brewer:") || spec %in% brewer_names)) {
    nm <- sub("^brewer:", "", spec)
    if (!nm %in% brewer_names)
      stop("Unknown RColorBrewer palette '", nm, "'. See ",
           "rownames(RColorBrewer::brewer.pal.info) for the available ones.",
           call. = FALSE)
    maxc <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
    cols <- RColorBrewer::brewer.pal(max(3, min(n, maxc)), nm)
  } else if (is.character(palette) && all(vapply(palette, is_colour, logical(1)))) {
    cols <- unname(palette)
  } else {
    stop("Unknown palette ", paste(deparse(palette), collapse = ""),
         ". Use \"pkg::name\" (paletteer), \"brewer:Name\" (RColorBrewer), ",
         "or a vector of colours.", call. = FALSE)
  }

  if (type == "discrete") {
    if (length(cols) > n) {
      cols <- cols[seq_len(n)]
    } else if (length(cols) < n) {
      cols <- grDevices::colorRampPalette(cols)(n)
    }
  }
  if (reverse) cols <- rev(cols)
  if (!is.null(levels)) names(cols) <- levels
  cols
}


#' Qualitative colours for the method comparison plots, PRONE style
#'
#' The colours of every RColorBrewer qualitative palette, concatenated in
#' reverse, as the PRONE package draws its normalisation comparisons. Used by
#' the normalisation and imputation metric plots, which each had a copy.
#'
#' @param n Number of colours.
#' @return Character vector of `n` colours.
#' @keywords internal
#' @noRd
.prone_colors <- function(n) {
  pal_info   <- RColorBrewer::brewer.pal.info
  qual_pals  <- pal_info[pal_info$category == "qual", ]
  col_vector <- rev(unlist(mapply(RColorBrewer::brewer.pal,
                                  qual_pals$maxcolors,
                                  rownames(qual_pals))))
  rep_len(col_vector, n)
}


# =============================================================================
# Feature filtering helpers
# =============================================================================

#' Build the adjusted p-value column name for a comparison
#'
#' @param comparison Name of the comparison, e.g. `"B-A"`.
#' @return Column name, e.g. `"adjP_B-A"`.
#' @keywords internal
#' @noRd
.adjp_col <- function(comparison) {
  paste0("adjP_", comparison)
}


#' Get the feature IDs according to the filtering mode
#'
#' Unifies the two copies that lived in `Heatmap_tidyHeatmap.R` and
#' `PCA_Highcharts_Final.R`. They differed in two respects:
#'
#' * The name of the third mode: `"target"` in the heatmap copy and `"specific"`
#'   in the PCA one. **Both** are accepted here as synonyms, because the public
#'   functions of each module propagate their own vocabulary and either one
#'   would have broken the other when merged into a single namespace.
#' * The `mode = "any"` filtering: the PCA copy used `sig_any == TRUE`, which
#'   lets `NA` `FeatureID`s slip through when `sig_any` contains `NA`. The
#'   `which()` from the heatmap copy is kept, which does not.
#'
#' @param data Data frame in long format with columns `FeatureID`, `sig_any` and
#'   `adjP_*`.
#' @param mode `"all"` (everything), `"any"` (significant in at least one
#'   comparison) or `"target"`/`"specific"` (significant in `comparison`).
#' @param alpha Significance threshold for the targeted mode.
#' @param comparison Comparison to use in the targeted mode.
#' @return Vector of `FeatureID`s meeting the criterion.
#' @keywords internal
#' @noRd
.get_feature_ids <- function(data,
                             mode = c("all", "any", "target", "specific"),
                             alpha = 0.05,
                             comparison = NULL) {

  mode <- match.arg(mode)
  data <- as.data.frame(data)

  feat <- data[!duplicated(data$FeatureID), , drop = FALSE]

  if (mode == "all") {
    return(feat$FeatureID)
  }

  if (mode == "any") {
    if (!("sig_any" %in% names(feat))) {
      stop("Column 'sig_any' is required for mode = 'any'")
    }
    # which() prevents NA FeatureIDs from slipping through when sig_any contains
    # NA (unlike feat$sig_any == TRUE, which would return NA rows).
    return(feat$FeatureID[which(feat$sig_any)])
  }

  # targeted mode: "target" and "specific" are synonyms
  if (is.null(comparison)) {
    stop("Argument 'comparison' is required for mode = '", mode, "'")
  }

  col <- .adjp_col(comparison)
  if (!(col %in% names(feat))) {
    stop("Column not found: ", col)
  }

  feat$FeatureID[which(feat[[col]] <= alpha)]
}


# =============================================================================
# Reading and writing delimited files
# =============================================================================
#
# NADIA used readr for these until 0.99.2, with base R fallbacks that never ran
# because readr was in Imports. These two helpers replace both, keeping what
# readr did that the old fallbacks did not: column names are left as they are
# (`B-A`, `Abundance: A_1`), empty fields are read as NA, and every number is
# read as a double. One thing is deliberately different: a field containing a
# tab, a newline or a quote is now quoted, whereas write_tsv() wrote it as it
# was and broke the table.

#' Read a delimited text file into a data.frame
#'
#' @param path File path. Gzip-compressed files are read transparently.
#' @param sep Field separator: `"\t"` (default) or `","`.
#' @return data.frame with the column names exactly as in the file and every
#'   numeric column as double, as readr returned them.
#' @keywords internal
#' @noRd
.read_delim_file <- function(path, sep = "\t") {
  df <- utils::read.table(path, header = TRUE, sep = sep, quote = "\"",
                          na.strings = c("", "NA"), check.names = FALSE,
                          stringsAsFactors = FALSE, comment.char = "",
                          fill = TRUE, encoding = "UTF-8")
  is_int <- vapply(df, is.integer, logical(1))
  df[is_int] <- lapply(df[is_int], as.double)
  df
}

#' Format doubles with the fewest significant digits that read back as the
#' same number
#'
#' `write.table()` rounds to 15 significant digits, which loses the last digit
#' of most measured values. 17 always suffice but turn 0.1 into
#' 0.10000000000000001, so each value gets 15, 16 or 17, whichever is the
#' shortest that parses back to it. `NaN` is written as missing, as readr did.
#'
#' @param x Double vector.
#' @return Character vector, `NA` where `x` is `NA` or `NaN`.
#' @keywords internal
#' @noRd
.format_double <- function(x) {
  out <- sprintf("%.15g", x)
  fin <- which(is.finite(x))
  for (d in 16:17) {
    bad <- fin[as.numeric(out[fin]) != x[fin]]
    if (length(bad) == 0) break
    out[bad] <- sprintf(paste0("%.", d, "g"), x[bad])
    fin <- bad
  }
  out[is.na(x)] <- NA_character_
  out
}

#' Write a data.frame as a tab-separated file
#'
#' @param x data.frame (or tibble) to write.
#' @param path Output path.
#' @param na String written for missing values.
#' @return `path`, invisibly.
#' @keywords internal
#' @noRd
.write_tsv <- function(x, path, na = "NA") {
  quote_if_needed <- function(v) {
    need <- !is.na(v) & grepl("[\t\n\r\"]", v)
    v[need] <- paste0("\"", gsub("\"", "\"\"", v[need], fixed = TRUE), "\"")
    v
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE, optional = TRUE)
  is_chr <- vapply(x, function(col) is.character(col) || is.factor(col),
                   logical(1))
  x[is_chr] <- lapply(x[is_chr], function(col) quote_if_needed(as.character(col)))
  # Plain doubles only: dates and times are doubles too, with a class.
  is_dbl <- vapply(x, function(col) is.double(col) && !is.object(col),
                   logical(1))
  x[is_dbl] <- lapply(x[is_dbl], .format_double)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = quote_if_needed(names(x)), na = na,
                     fileEncoding = "UTF-8")
  invisible(path)
}


# =============================================================================
# Packages in Suggests
# =============================================================================

#' Stop with an informative error if a suggested package is not installed
#'
#' Every package in `Suggests` is checked with this helper at the entry of the
#' function that needs it, so that a missing package says what to install
#' instead of failing deep inside the code with "there is no package called".
#'
#' @param pkgs Names of the packages required.
#' @param what Short description of what needs them, used as the subject of
#'   the error message (e.g. `"volcano_highchart_list()"`).
#' @return `TRUE`, invisibly, if every package is installed.
#' @keywords internal
#' @noRd
.check_suggests <- function(pkgs, what) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1),
                               quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(what, " requires ",
         paste(sprintf("'%s'", missing_pkgs), collapse = ", "),
         ". Install ", if (length(missing_pkgs) > 1) "them" else "it",
         " with BiocManager::install(",
         if (length(missing_pkgs) > 1) "c(" else "",
         paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
         if (length(missing_pkgs) > 1) ")" else "", ").",
         call. = FALSE)
  }
  invisible(TRUE)
}


# =============================================================================
# Mfuzz dependencies that must be attached
# =============================================================================
#
# `Mfuzz` declares `Depends: Biobase, e1071` and calls `exprs()` and `cmeans()`
# unqualified, so its functions only resolve those names if both packages are on
# the search path. `requireNamespace()` is not enough: it loads the namespace but
# does not attach it. This used to be taken care of by the `library(Mfuzz)` call
# in the module header, which pulled in its Depends; a package cannot do that.
#
# The solution is to attach them only while the Pattern Profiler is running and
# to leave the search path as it was, so that the user's session is untouched.

#' Attach the dependencies that Mfuzz needs on the search path
#'
#' @return Vector with the packages that this call has attached (possibly empty
#'   if they were already attached). It must be passed to
#'   `.mfuzz_deps_detach()`.
#' @keywords internal
#' @noRd
.mfuzz_deps_attach <- function() {
  .check_suggests(c("Mfuzz", "Biobase", "e1071"), "The Pattern Profiler")
  # Only Biobase and e1071 need to be attached; Mfuzz is used via Mfuzz::
  attachable <- c("Biobase", "e1071")
  already_on <- paste0("package:", attachable) %in% search()
  # The startup banner is suppressed: attaching these is an implementation
  # detail of the clustering, not something the caller asked for, and it would
  # otherwise print even under verbose = FALSE.
  for (p in attachable[!already_on]) {
    suppressPackageStartupMessages(attachNamespace(asNamespace(p)))
  }
  attachable[!already_on]
}

#' Undo what `.mfuzz_deps_attach()` did
#'
#' @param pkgs Vector returned by `.mfuzz_deps_attach()`.
#' @return `NULL`, invisibly. Called for its side effect.
#' @keywords internal
#' @noRd
.mfuzz_deps_detach <- function(pkgs) {
  for (p in pkgs) {
    entry <- paste0("package:", p)
    if (entry %in% search()) {
      detach(entry, character.only = TRUE, unload = FALSE)
    }
  }
  invisible(NULL)
}


# =============================================================================
# Shared preprocessing helpers
# =============================================================================

#' Detect the `Abundance:` columns and derive Coding / Condition / Replicate
#'
#' Two of the wide formats declare their design in the column names themselves:
#' Proteome Discoverer writes `Abundance: <Condition>_<Replicate>`, and a DIA-NN
#' matrix can be renamed to the same convention so that it needs no annotation
#' file either. The parser is therefore shared by `preprocess_tmt()` and
#' `preprocess_diann()`.
#'
#' @param df Data frame read from the export.
#' @return A data frame with `abundance_col`, `Coding`, `R.Condition` and
#'   `R.Replicate`, one row per column found.
#' @keywords internal
#' @noRd
.parse_abundance_columns <- function(df) {
  abund_cols <- grep("^Abundance:\\s*", names(df), value = TRUE)
  if (length(abund_cols) == 0) {
    stop("No 'Abundance:' columns found in the file. ",
         "Check that it is a Proteome Discoverer export.")
  }

  coding <- sub("^Abundance:\\s*", "", abund_cols)
  m <- stringr::str_match(coding, "^(.+)_(\\d+)$")
  if (any(is.na(m[, 1]))) {
    bad <- coding[is.na(m[, 1])]
    stop(
      "Could not parse the Abundance suffixes (expected <Condition>_<Replicate>):\n  - ",
      paste(bad, collapse = "\n  - ")
    )
  }

  data.frame(
    abundance_col = abund_cols,
    Coding        = coding,
    R.Condition   = m[, 2],
    R.Replicate   = as.integer(m[, 3]),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Shared by the four preprocessing functions
# =============================================================================

#' Drop conditions that were asked for but are not in the data
#'
#' @description
#' Every `preprocess_*()` function takes `condition_order`, uses it to select the
#' samples to keep, and then makes it the levels of the `R.Condition` factor.
#' Those two roles disagree when a condition is listed but does not occur in the
#' file: the selection finds nothing to keep, while the factor gains a level with
#' no samples behind it. An empty level is not inert -- it travels into
#' `colData()` and from there into legends and axes.
#'
#' Rather than fail (the other conditions are perfectly usable), this warns and
#' narrows `condition_order` to what is actually there. `setdiff()` preserves the
#' order of its first argument, so the order that was asked for survives.
#'
#' @param condition_order Character vector as supplied by the caller.
#' @param present The conditions occurring in the data, in any form.
#' @param source Word naming where they were looked for, for the message.
#' @return `condition_order` without the absent conditions.
#' @keywords internal
#' @noRd
.drop_absent_conditions <- function(condition_order, present, source = "report") {
  absent <- setdiff(condition_order, unique(as.character(present)))
  if (length(absent) > 0) {
    warning("Conditions listed in condition_order but absent from the ", source,
            ", and dropped: ", paste(absent, collapse = ", "), call. = FALSE)
    condition_order <- setdiff(condition_order, absent)
  }
  condition_order
}

#' Read a sample sheet and turn it into a design table
#'
#' @description
#' The escape hatch for reports whose intensity columns are not named
#' `<Condition>_<Replicate>`. It is a file rather than a vector of labels on
#' purpose: a vector is positional, and column order is not something an export
#' guarantees, so a re-export that reorders the columns would relabel every run
#' without any way to notice. Matching by name either matches or fails.
#'
#' Replicate numbers come from the identifier when it ends in `_<digits>`, which
#' is what the suffix route does, so a sheet stating the obvious gives exactly
#' the same answer as no sheet at all. Identifiers that do not carry a number --
#' raw-file names, TMT tags -- are numbered sequentially within their condition.
#' Either way the number does not depend on the order the sheet lists them in,
#' so a reordered sheet cannot relabel anything.
#'
#' @param annot_path Path to a TSV with `Column` and `Condition`. Any other
#'   column is ignored: batch and covariate structure belongs in `covariate_df`
#'   of [process_proteomics()], not here.
#' @param ids Sample identifiers as they appear in the report, one per intensity
#'   column, in column order.
#' @param cols The intensity column names themselves, parallel to `ids`.
#' @param verbose Report the columns the sheet leaves out.
#' @return A data frame with `abundance_col`, `Coding`, `R.Condition` and
#'   `R.Replicate` -- the same shape `.parse_abundance_columns()` returns, so
#'   that everything downstream is blind to which route was taken.
#' @keywords internal
#' @noRd
.design_from_annotation <- function(annot_path, ids, cols, verbose = TRUE) {
  if (!file.exists(annot_path)) stop("Annotation file not found: ", annot_path)

  annot <- utils::read.delim(annot_path, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE)

  missing_cols <- setdiff(c("Column", "Condition"), names(annot))
  if (length(missing_cols) > 0) {
    stop("Required columns missing from the annotation file:\n  - ",
         paste(missing_cols, collapse = "\n  - "),
         "\nIt must have at least 'Column' and 'Condition'.")
  }
  annot$Column    <- as.character(annot$Column)
  annot$Condition <- as.character(annot$Condition)

  if (anyDuplicated(annot$Column) > 0) {
    stop("The annotation file has duplicated values in 'Column':\n  - ",
         paste(unique(annot$Column[duplicated(annot$Column)]),
               collapse = "\n  - "))
  }

  hit <- match(annot$Column, ids)
  if (anyNA(hit)) {
    stop("These samples of the annotation file have no column in the report:",
         "\n  - ", paste(annot$Column[is.na(hit)], collapse = "\n  - "),
         "\nThe report has:\n  - ", paste(ids, collapse = "\n  - "))
  }

  # Columns the sheet leaves out are dropped, which is how a pool or a blank is
  # excluded. Worth saying out loud, since it silently shrinks the experiment.
  extra <- setdiff(ids, annot$Column)
  if (length(extra) > 0 && verbose) {
    message("Note: columns ignored (not listed in the annotation file): ",
            paste(extra, collapse = ", "))
  }

  # Take the replicate from the identifier when it has one, so that a sheet
  # spelling out what the suffixes already say lands on the same answer. Fall
  # back to a per-condition counter for identifiers that carry no number, and
  # take that fallback for the whole sheet rather than per row, so the numbering
  # cannot be a mixture of two schemes.
  m <- stringr::str_match(annot$Column, "^(.+)_(\\d+)$")
  replicate <- if (anyNA(m[, 1])) {
    as.integer(stats::ave(seq_len(nrow(annot)), annot$Condition, FUN = seq_along))
  } else {
    as.integer(m[, 3])
  }

  data.frame(
    abundance_col = cols[hit],
    Coding        = annot$Column,
    R.Condition   = annot$Condition,
    R.Replicate   = replicate,
    stringsAsFactors = FALSE
  )
}

#' Resolve the design of a Proteome Discoverer export
#'
#' @description Suffixes if they are there, sample sheet if they are not. Shared
#'   by `preprocess_tmt()` and `preprocess_lfq()`, whose exports both name their
#'   intensity columns `Abundance: <sample>`.
#' @keywords internal
#' @noRd
.pd_resolve_design <- function(df, annot_path = NULL, verbose = TRUE) {
  if (is.null(annot_path)) return(.parse_abundance_columns(df))

  cols <- grep("^Abundance:\\s*", names(df), value = TRUE)
  if (length(cols) == 0) {
    stop("No 'Abundance:' columns found in the file. ",
         "Check that it is a Proteome Discoverer export.")
  }
  .design_from_annotation(annot_path, sub("^Abundance:\\s*", "", cols), cols,
                          verbose = verbose)
}

# =============================================================================
# Proteome Discoverer helpers
# =============================================================================
#
# These were duplicated between Preprocessing_LFQ.R and Preprocessing_TMT.R. The
# two copies of `.parse_gene_from_description` were identical; those of
# `.validate_pd_columns` differed only in the text of the error message.

#' Extract the Gene Name from the Description column (embedded UniProt format)
#'
#' Looks for the `GN=<gene>` pattern usual in Proteome Discoverer exports.
#'
#' @param x Vector of descriptions.
#' @return Vector of gene names, `NA` where there is no match.
#' @keywords internal
#' @noRd
.parse_gene_from_description <- function(x) {
  m <- stringr::str_match(x, "GN=([^ ]+)")
  m[, 2]
}

#' Neutral names for the protein score, with the report's own name kept
#'
#' NADIA's tables follow Spectronaut's column names, and the readers of the
#' other formats used to put their score under Spectronaut's `PG.Cscore`, so
#' a Proteome Discoverer `Sum PEP Score` read as a C-score (Bioconductor
#' review, #186). The score columns are named `PG.Score` (one value per
#' protein) and `PG.Score.RunWise_<sample>` in every format, and the metric of
#' the report they hold is kept in the `nadia_labels` attribute of the table,
#' which the interactive tables use as column label.
#'
#' @param df `protein_id` or `protein_quant` table.
#' @param global,run_wise Name of the metric behind `PG.Score` and behind
#'   `PG.Score.RunWise_*` in the report, or `NULL` when there is none.
#' @return `df` as a data frame, with the columns renamed and the attribute set.
#' @keywords internal
#' @noRd
.score_columns <- function(df, global = NULL, run_wise = NULL) {
  df <- .score_legacy_names(as.data.frame(df))
  present <- c(PG.Score = "PG.Score" %in% names(df),
               PG.Score.RunWise = any(startsWith(names(df), "PG.Score.RunWise_")))
  labels <- c(PG.Score = global, PG.Score.RunWise = run_wise)
  labels <- labels[names(labels) %in% names(present)[present]]
  attr(df, "nadia_labels") <- if (length(labels)) labels else NULL
  df
}

#' Rename the score columns of a table written before NADIA 0.99.2
#'
#' The score used to be stored under Spectronaut's `PG.Cscore` names in every
#' format. A table exported or a `.nadia` file written before keeps loading;
#' its label is the neutral one, since its source was not recorded.
#' @param df A Protein_ID or Protein_QUANT table.
#' @return `df`, renamed.
#' @keywords internal
#' @noRd
.score_legacy_names <- function(df) {
  old <- grepl("^PG\\.Cscore", names(df))
  if (any(old)) names(df)[old] <- sub("^PG\\.Cscore", "PG.Score", names(df)[old])
  df
}

#' Label of a score column, as the report calls it
#' @param df A table, possibly carrying the `nadia_labels` attribute.
#' @param column `"PG.Score"` or `"PG.Score.RunWise"`.
#' @return The report's name for the metric, or `"Score"` when unknown (a
#'   table read back from a TSV, for instance).
#' @keywords internal
#' @noRd
.score_label <- function(df, column) {
  lab <- attr(df, "nadia_labels")[column]
  if (length(lab) == 1 && !is.na(lab) && nzchar(lab)) unname(lab) else "Score"
}

#' Record where an object came from
#'
#' Attaches the call that produced an object and a fingerprint of the files it
#' was read from, as the attributes `nadia_call` and `nadia_source`.
#'
#' They are attributes rather than list elements on purpose: every preprocessing
#' function documents a return value of exactly three elements (`metadata`,
#' `protein_id`, `protein_quant`), and adding a fourth would change a contract
#' that code and documentation already depend on. An attribute travels with the
#' object without appearing in `names()` or in `str()`.
#'
#' The fingerprint is what lets someone check, years later, whether a result
#' corresponds to the raw file in front of them. The MD5 is skipped for files
#' over 500 MB, where hashing costs more than the answer is worth.
#'
#' @param x     Object to stamp.
#' @param call  The call to record, normally `match.call()` from the caller.
#' @param files Named character vector of paths; the names become the `role`
#'   column (`report`, `annotation`).
#' @return `x`, with the two attributes attached.
#' @keywords internal
#' @noRd
.nadia_stamp <- function(x, call, files = character(0)) {
  attr(x, "nadia_call") <- call

  files <- files[!vapply(files, is.null, logical(1))]
  files <- unlist(files)

  if (length(files) > 0) {
    info <- file.info(files)
    attr(x, "nadia_source") <- data.frame(
      role       = names(files),
      path       = normalizePath(files, mustWork = FALSE),
      size_bytes = as.numeric(info$size),
      mtime      = format(info$mtime, tz = "UTC", usetz = FALSE),
      md5        = vapply(seq_along(files), function(i) {
        if (is.na(info$size[i]) || info$size[i] > 500e6) {
          NA_character_
        } else {
          unname(tools::md5sum(files[i]))
        }
      }, character(1)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  x
}


#' Validate the minimum columns of a Proteome Discoverer export
#'
#' @param df Data frame read from the export.
#' @return `TRUE` invisibly; aborts with an error if any required column is
#'   missing.
#' @keywords internal
#' @noRd
.validate_pd_columns <- function(df) {
  required <- c("Accession", "Description")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "Required columns missing from the data file:\n  - ",
      paste(missing, collapse = "\n  - "),
      "\nCheck that the file is a Proteome Discoverer protein export."
    )
  }
  invisible(TRUE)
}
