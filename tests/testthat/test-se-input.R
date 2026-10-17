# =============================================================================
# A SummarizedExperiment as input: preprocess_se() and process_proteomics()
#
# The way in for data imported with other tools (Bioconductor review, #186).
# The example data rebuilt as a SummarizedExperiment must give the same
# analysis as the Spectronaut report they came from.
# =============================================================================

# nadia_dia as another tool would hand it over: sample columns named after the
# raw files, condition and replicate in colData, genes in rowData
.se_from_dia <- function(log2 = FALSE) {
    data(nadia_dia, package = "NADIA", envir = environment())
    pq <- nadia_dia$protein_quant
    qc <- grep("^PG.Quantity_", names(pq), value = TRUE)
    md <- nadia_dia$metadata[match(sub("^PG.Quantity_", "", qc), nadia_dia$metadata$Coding), ]
    x <- as.matrix(pq[, qc])
    dimnames(x) <- list(pq$PG.ProteinGroups, md$R.FileName)
    if (log2) x <- log2(x)
    SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = x),
        colData = data.frame(Condition = md$R.Condition, Replicate = md$R.Replicate),
        rowData = data.frame(Genes = pq$PG.Genes))
}

test_that("a SummarizedExperiment gives the same analysis as the report", {
    data(nadia_dia, package = "NADIA")
    ref <- process_proteomics(nadia_dia, verbose = FALSE)
    res <- process_proteomics(.se_from_dia(), verbose = FALSE)
    expect_identical(res$DEPs_results, ref$DEPs_results)
    expect_identical(SummarizedExperiment::assays(res$se),
                     SummarizedExperiment::assays(ref$se))
    # Calling the converter first is the same thing
    pre <- preprocess_se(.se_from_dia(), verbose = FALSE)
    expect_identical(process_proteomics(pre, verbose = FALSE)$DEPs_results,
                     ref$DEPs_results)
})

test_that("log2 intensities give the same analysis once declared", {
    lin <- preprocess_se(.se_from_dia(), verbose = FALSE)
    lg  <- preprocess_se(.se_from_dia(log2 = TRUE), log2_input = TRUE, verbose = FALSE)
    expect_equal(lg$protein_quant, lin$protein_quant)
    # Undeclared, they are caught, and so are linear values declared as log2
    expect_warning(preprocess_se(.se_from_dia(log2 = TRUE), verbose = FALSE),
                   "look log2-transformed")
    expect_warning(preprocess_se(.se_from_dia(), log2_input = TRUE, verbose = FALSE),
                   "look like linear")
})

test_that("preprocess_se() returns the object every reader returns", {
    pre <- preprocess_se(.se_from_dia(), verbose = FALSE)
    expect_s3_class(pre, c("se_data", "proteomics_data"))
    expect_named(pre, c("metadata", "protein_id", "protein_quant"))
    expect_identical(pre$metadata$Coding,
                     paste(pre$metadata$R.Condition, pre$metadata$R.Replicate, sep = "_"))
    expect_identical(grep("^PG.Quantity_", names(pre$protein_quant), value = TRUE),
                     paste0("PG.Quantity_", pre$metadata$Coding))
    # The original sample names are kept
    expect_true(all(grepl("CursoProtQ", pre$metadata$R.FileName)))
    expect_false(all(is.na(pre$protein_quant$PG.Genes)))
    expect_output(print(pre), "SummarizedExperiment")
})

test_that("the design is read from colData, and can be subset", {
    se <- .se_from_dia()
    se$Replicate <- NULL
    se <- se[, order(se$Condition, decreasing = TRUE)]
    pre <- preprocess_se(se, condition_order = c("B", "A"), verbose = FALSE)
    # Conditions kept and ordered as asked, replicates numbered within each
    expect_identical(as.character(unique(pre$metadata$R.Condition)), c("B", "A"))
    expect_identical(pre$metadata$R.Replicate, rep(1:4, 2))
    expect_false(any(grepl("_D_|^PG.Quantity_D", names(pre$protein_quant))))
})

test_that("malformed objects are rejected with a reason", {
    se <- .se_from_dia()
    no_cond <- se; no_cond$Condition <- NULL
    expect_error(preprocess_se(no_cond, verbose = FALSE), "no 'Condition' column")
    bad <- se; bad$Condition <- paste0(bad$Condition, "_x")
    expect_error(preprocess_se(bad, verbose = FALSE), "without '_'")
    dup <- se; rownames(dup)[2] <- rownames(dup)[1]
    expect_error(preprocess_se(dup, verbose = FALSE), "Duplicated row names")
    twice <- se; twice$Replicate <- 1
    expect_error(preprocess_se(twice, verbose = FALSE), "share condition and replicate")
    expect_error(preprocess_se(se, assay_name = "nope", verbose = FALSE), "not found")
    expect_error(preprocess_se(data.frame(a = 1), verbose = FALSE), "must be a SummarizedExperiment")
    # The SE-only arguments are not silently ignored for a report object
    data(nadia_dia, package = "NADIA")
    expect_warning(
        expect_error(process_proteomics(nadia_dia, log2_input = TRUE, comparisons = "nope",
                                        verbose = FALSE)),
        "apply only to a SummarizedExperiment")
})

test_that("an imported SummarizedExperiment can be stored in a .nadia file", {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")
    pre <- preprocess_se(.se_from_dia(), verbose = FALSE)
    f <- withr::local_tempfile(fileext = ".nadia")
    write_nadia(f, pre, verbose = FALSE)
    back <- nadia_preprocessing(read_nadia(f))
    expect_s3_class(back, "se_data")
    expect_identical(back$protein_quant, pre$protein_quant)
    expect_identical(back$metadata, pre$metadata)
})
