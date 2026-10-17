# =============================================================================
# Preprocessing_SE.R -- a SummarizedExperiment as input
# =============================================================================
#
# The other four readers parse a search-engine report. This one takes protein
# intensities that are already in a SummarizedExperiment, built by any other
# tool (QFeatures::readQFeatures(), a limpa or MSnbase import converted to an
# SE, or by hand), and returns the same `proteomics_data` object, so the rest
# of the pipeline, the tables and the .nadia container use it unchanged
# (Bioconductor review, #186).
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================

#' Import Protein Intensities from a SummarizedExperiment
#'
#' Converts a `SummarizedExperiment` of protein intensities into the
#' `proteomics_data` object that every NADIA reader returns. It is the way in
#' for data imported with other Bioconductor tools, for instance one of the
#' assays of a `QFeatures` object. [process_proteomics()] also accepts the
#' `SummarizedExperiment` directly and calls this function; call it yourself to
#' keep the object, for example to store it with [write_nadia()].
#'
#' What is read from the object:
#' - **Intensities**: the assay `assay_name`, proteins in rows and samples in
#'   columns. With `log2_input = TRUE` the values are taken back to the linear
#'   scale, because NADIA applies the log2 transformation when it normalises.
#'   Zeros are treated as missing values, as for the other readers.
#' - **Proteins**: the row names, which must be present and unique. The gene
#'   names come from `rowData()$Genes` and the descriptions from
#'   `rowData()$Description`, when those columns exist.
#' - **Design**: the condition from `colData()[[condition_column]]` and the
#'   replicate from `colData()[[replicate_column]]`; without a replicate
#'   column the samples are numbered within each condition, in column order.
#'
#' The samples are renamed `<condition>_<replicate>`, the convention the rest
#' of NADIA relies on, and the original column names are kept in
#' `metadata$R.FileName`. Identification metrics (precursor or peptide counts,
#' scores) are not read from `rowData()`, so those columns of the tables are
#' left empty.
#'
#' @param se A `SummarizedExperiment` with protein intensities.
#' @param condition_order Optional character vector with the conditions to
#'   keep, in the order they should appear. `NULL` (default) keeps them all,
#'   in the order of the factor levels or, for text, of first appearance.
#' @param assay_name Name or index of the assay with the intensities. `NULL`
#'   (default) takes the first one.
#' @param log2_input `TRUE` if the intensities are on the log2 scale.
#'   Default `FALSE`, linear intensities.
#' @param condition_column Column of `colData(se)` with the condition of each
#'   sample. Default `"Condition"`.
#' @param replicate_column Column of `colData(se)` with the replicate number.
#'   Default `"Replicate"`; if the column is absent, the replicates are
#'   numbered.
#' @param verbose Print progress messages. Default `TRUE`.
#'
#' @return An object of class `c("se_data", "proteomics_data", "list")` with
#'   the elements `metadata`, `protein_id` and `protein_quant`, like the other
#'   readers.
#'
#' @seealso [process_proteomics()], which accepts the `SummarizedExperiment`
#'   directly.
#'
#' @examples
#' data(nadia_dia)
#' # A SummarizedExperiment built from the example data
#' quant <- as.matrix(nadia_dia$protein_quant[, grep("^PG\\.Quantity_",
#'                                                   names(nadia_dia$protein_quant))])
#' rownames(quant) <- nadia_dia$protein_quant$PG.ProteinGroups
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(intensity = quant),
#'   colData = data.frame(Condition = nadia_dia$metadata$R.Condition,
#'                        Replicate = nadia_dia$metadata$R.Replicate),
#'   rowData = data.frame(Genes = nadia_dia$protein_quant$PG.Genes))
#'
#' pre <- preprocess_se(se, verbose = FALSE)
#' pre
#'
#' # The same values on the log2 scale
#' SummarizedExperiment::assay(se, "intensity") <- log2(quant)
#' pre_log2 <- preprocess_se(se, log2_input = TRUE, verbose = FALSE)
#'
#' @export
preprocess_se <- function(
    se,
    condition_order = NULL,
    assay_name = NULL,
    log2_input = FALSE,
    condition_column = "Condition",
    replicate_column = "Replicate",
    verbose = TRUE
) {
  # ==========================================================================
  # Validation
  # ==========================================================================
  if (inherits(se, "QFeatures")) {
    stop("A QFeatures object holds several assays: pass the protein-level ",
         "one, for example qf[[\"proteins\"]].", call. = FALSE)
  }
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  if (!is.logical(log2_input) || length(log2_input) != 1 || is.na(log2_input)) {
    stop("`log2_input` must be TRUE or FALSE.", call. = FALSE)
  }
  if (ncol(se) == 0 || nrow(se) == 0) {
    stop("The SummarizedExperiment is empty.", call. = FALSE)
  }

  assays <- SummarizedExperiment::assayNames(se)
  if (is.null(assay_name)) assay_name <- 1L
  if (is.character(assay_name) && !assay_name %in% assays) {
    stop("Assay '", assay_name, "' not found. Available: ",
         paste(assays, collapse = ", "), call. = FALSE)
  }
  x <- as.matrix(SummarizedExperiment::assay(se, assay_name, withDimnames = TRUE))
  if (!is.numeric(x)) stop("The assay must be numeric.", call. = FALSE)

  ids <- rownames(se)
  if (is.null(ids)) {
    stop("The SummarizedExperiment has no row names: they are the protein ",
         "identifiers.", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop("Duplicated row names (protein identifiers): ",
         paste(utils::head(unique(ids[duplicated(ids)]), 5), collapse = ", "),
         call. = FALSE)
  }

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  if (!condition_column %in% names(cd)) {
    stop("colData(se) has no '", condition_column, "' column. Name the ",
         "column holding the condition of each sample in `condition_column`.",
         call. = FALSE)
  }
  cond <- cd[[condition_column]]
  if (anyNA(cond) || any(!nzchar(as.character(cond)))) {
    stop("The condition is missing for some samples.", call. = FALSE)
  }
  levels_found <- if (is.factor(cond)) levels(droplevels(cond)) else unique(as.character(cond))
  cond <- as.character(cond)
  # The conditions enter limma's contrasts, so they must be syntactic names.
  # An underscore is fine, as for the other readers: the replicate is read
  # from the trailing digits of <condition>_<replicate>.
  bad <- levels_found[make.names(levels_found) != levels_found]
  if (length(bad) > 0) {
    stop("Condition names must be syntactically valid R names (they become ",
         "limma contrasts such as B-A): ", paste(bad, collapse = ", "),
         call. = FALSE)
  }

  if (replicate_column %in% names(cd)) {
    rep_num <- cd[[replicate_column]]
    if (anyNA(rep_num) || !is.numeric(rep_num) || any(rep_num != round(rep_num))) {
      stop("colData(se)$", replicate_column, " must hold whole numbers.",
           call. = FALSE)
    }
    rep_num <- as.integer(rep_num)
  } else {
    rep_num <- stats::ave(seq_along(cond), cond, FUN = seq_along)
  }
  coding <- paste(cond, rep_num, sep = "_")
  if (anyDuplicated(coding)) {
    stop("Two samples share condition and replicate: ",
         paste(unique(coding[duplicated(coding)]), collapse = ", "),
         call. = FALSE)
  }

  # --- Scale: NADIA log2-transforms when it normalises ---
  finite <- x[is.finite(x)]
  if (log2_input) {
    if (length(finite) && max(finite) > 100) {
      warning("log2_input = TRUE, but the values reach ", signif(max(finite), 3),
              ": they look like linear intensities.", call. = FALSE)
    }
    x <- 2^x
  } else if (length(finite) && max(finite) < 40) {
    warning("The intensities do not exceed ", signif(max(finite), 3),
            ": they look log2-transformed. If they are, set log2_input = TRUE.",
            call. = FALSE)
  }
  if (length(finite) && any(x[is.finite(x)] < 0)) {
    stop("Negative linear intensities: check the assay and `log2_input`.",
         call. = FALSE)
  }

  # --- Conditions to keep, and their order ---
  if (is.null(condition_order)) {
    condition_order <- levels_found
  } else {
    if (!is.character(condition_order) || length(condition_order) == 0) {
      stop("condition_order must be a non-empty character vector.", call. = FALSE)
    }
    if (!any(cond %in% condition_order)) {
      stop("No sample matches condition_order = c(",
           paste0("'", condition_order, "'", collapse = ", "), ").\n",
           "  Conditions in colData: ", paste(levels_found, collapse = ", "),
           call. = FALSE)
    }
    condition_order <- .drop_absent_conditions(condition_order, cond)
  }
  keep <- which(cond %in% condition_order)
  ord <- keep[order(factor(cond[keep], levels = condition_order), rep_num[keep])]
  coding <- coding[ord]

  if (verbose) {
    message("Reading ", length(ids), " proteins x ", length(ord), " samples ",
            "from assay '", if (is.character(assay_name)) assay_name else assays[assay_name],
            "'", if (log2_input) " (log2, taken back to linear)" else "", "...")
  }

  # ==========================================================================
  # Tables, in the shape the other readers produce
  # ==========================================================================
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  text_col <- function(name) {
    if (name %in% names(rd)) as.character(rd[[name]]) else rep(NA_character_, length(ids))
  }
  base_info <- data.frame(
    PG.ProteinGroups       = ids,
    PG.ProteinDescriptions = text_col("Description"),
    PG.Genes               = text_col("Genes"),
    PG.MolecularWeight     = NA_real_,
    stringsAsFactors       = FALSE
  )
  # No identification metrics are read from the object. The precursor
  # families are created empty, as for DIA-NN, because the protein tables need
  # one of them to find the samples.
  empty_family <- function(prefix) {
    m <- matrix(NA_real_, nrow = length(ids), ncol = length(coding),
                dimnames = list(NULL, paste0(prefix, "_", coding)))
    as.data.frame(m)
  }
  quant <- as.data.frame(x[, ord, drop = FALSE])
  names(quant) <- paste0("PG.Quantity_", coding)
  rownames(quant) <- NULL

  protein_ID <- cbind(base_info, empty_family("PG.NrOfPrecursorsIdentified"))
  protein_QUANT <- cbind(base_info,
                         empty_family("PG.NrOfPrecursorsUsedForQuantification"),
                         quant)

  run_summary <- data.frame(
    R.FileName  = colnames(se)[ord] %||% coding,
    R.Condition = factor(cond[ord], levels = condition_order, ordered = TRUE),
    R.Replicate = rep_num[ord],
    Coding      = coding,
    R.ProteinGroupsIdentified = as.numeric(colSums(!is.na(x[, ord, drop = FALSE]) &
                                                     x[, ord, drop = FALSE] > 0)),
    stringsAsFactors = FALSE
  )
  rownames(run_summary) <- run_summary$Coding

  if (verbose) {
    message("Processing complete:\n",
            "  - Samples: ", nrow(run_summary), "\n",
            "  - Conditions: ", paste(condition_order, collapse = ", "), "\n",
            "  - Proteins: ", nrow(protein_QUANT))
  }

  result <- list(
    metadata      = run_summary,
    protein_id    = protein_ID,
    protein_quant = protein_QUANT
  )
  class(result) <- c("se_data", "proteomics_data", "list")
  .nadia_stamp(result, match.call())
}

# =============================================================================
# Methods for the se_data class
# =============================================================================

#' Print a Summary of Data Imported from a SummarizedExperiment
#'
#' @param x An `se_data` object, as returned by [preprocess_se()].
#' @param ... Ignored, present for compatibility with the `print` generic.
#' @return `x`, invisibly. Called for the summary it prints.
#' @examples
#' data(nadia_dia)
#' quant <- as.matrix(nadia_dia$protein_quant[, grep("^PG\\.Quantity_",
#'                                                   names(nadia_dia$protein_quant))])
#' rownames(quant) <- nadia_dia$protein_quant$PG.ProteinGroups
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(intensity = quant),
#'   colData = data.frame(Condition = nadia_dia$metadata$R.Condition))
#' print(preprocess_se(se, verbose = FALSE))
#' @export
print.se_data <- function(x, ...) {
  cat("Protein intensities imported from a SummarizedExperiment\n")
  cat("--------------------------------------------------------\n")
  cat("Samples (metadata):", nrow(x$metadata), "\n")
  cat("Proteins:", nrow(x$protein_quant), "\n")
  cat("\nConditions:",
      paste(levels(x$metadata$R.Condition) %||%
              unique(x$metadata$R.Condition), collapse = ", "),
      "\n")
  invisible(x)
}
