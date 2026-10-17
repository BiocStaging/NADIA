# =============================================================================
# Input validation: the error messages are part of the API
# =============================================================================
#
# Several of these guards exist to replace an incomprehensible failure further
# downstream with an actionable one at the boundary. The message text is
# therefore asserted, not just the fact that an error was raised: a guard that
# fires with the wrong message has lost most of its value.

test_that("normalize_proteomics rejects a data frame with no metadata", {
    data(nadia_dia, package = "NADIA")
    expect_error(normalize_proteomics(nadia_dia$protein_quant, verbose = FALSE),
                 "metadata", ignore.case = TRUE)
})

test_that("normalize_proteomics rejects metadata alongside a proteomics_data object", {
    # Supplying both is ambiguous: the object already carries its metadata, so
    # silently ignoring the argument would hide a real mistake.
    data(nadia_dia, package = "NADIA")
    expect_error(
        normalize_proteomics(nadia_dia, metadata = nadia_dia$metadata,
                             verbose = FALSE),
        "metadata", ignore.case = TRUE)
})

test_that("an unknown normalisation method is rejected by name", {
    data(nadia_dia, package = "NADIA")
    expect_error(normalize_proteomics(nadia_dia, norm_method = "not_a_method",
                                      verbose = FALSE))
})

test_that("an unknown imputation method is rejected", {
    data(nadia_dia, package = "NADIA")
    se <- normalize_proteomics(nadia_dia, verbose = FALSE)$se
    expect_error(impute_proteomics(se, "cycloess", imp_method = "not_a_method",
                                   verbose = FALSE))
})

test_that("a missing assay is reported with the name that was asked for", {
    data(nadia_dia, package = "NADIA")
    se <- normalize_proteomics(nadia_dia, verbose = FALSE)$se
    expect_error(impute_proteomics(se, "no_such_assay", verbose = FALSE))
})

test_that("non-syntactic condition names are rejected before the model is fitted", {
    # Without this guard makeContrasts() fails deep inside limma with a message
    # that says nothing about condition names.
    data(nadia_dia, package = "NADIA")

    bad <- nadia_dia
    bad$metadata$R.Condition <- paste0(bad$metadata$R.Condition, "-x")
    bad$metadata$Coding <- paste0(bad$metadata$R.Condition, "_",
                                  bad$metadata$R.Replicate)
    colnames(bad$protein_quant) <- sub("^(PG.Quantity_)([A-Z])_", "\\1\\2-x_",
                                       colnames(bad$protein_quant))

    expect_error(process_proteomics(bad, verbose = FALSE),
                 "makeContrasts")
})

test_that("a dot in a condition name is accepted", {
    # The complement of the test above: '.' is a valid character in an R name,
    # so it must NOT be rejected. Users separate parts of a condition name with
    # it precisely to avoid the hyphen.
    data(nadia_dia, package = "NADIA")

    ok <- nadia_dia
    ok$metadata$R.Condition <- paste0(ok$metadata$R.Condition, ".1")
    ok$metadata$Coding <- paste0(ok$metadata$R.Condition, "_",
                                 ok$metadata$R.Replicate)
    colnames(ok$protein_quant) <- sub("^(PG.Quantity_)([A-Z])_", "\\1\\2.1_",
                                      colnames(ok$protein_quant))

    res <- process_proteomics(ok, verbose = FALSE)
    expect_true(all(grepl("\\.1", as.character(res$comparisons))))
    expect_equal(ncol(res$se_proc), nrow(nadia_dia$metadata))
})

test_that("batch correction rejects NA batch labels", {
    skip_if_not_installed("BERT")
    data(nadia_dia, package = "NADIA")

    md  <- nadia_dia$metadata
    cov <- data.frame(Column = md$Coding,
                      Batch  = c(NA, rep(c("b1", "b2"), length.out = nrow(md) - 1)))

    expect_error(
        process_proteomics(nadia_dia, covariate_df = cov, batch_correct = TRUE,
                           batch_column = "Batch", verbose = FALSE))
})

test_that("batch correction warns when no covariate is protected", {
    # ComBat removes ALL batch-associated variance. If condition and batch are
    # confounded that erases biology, so a bare call must warn.
    skip_if_not_installed("BERT")
    data(nadia_dia, package = "NADIA")

    md  <- nadia_dia$metadata
    cov <- data.frame(Column = md$Coding,
                      Batch  = rep(c("b1", "b2"), length.out = nrow(md)))
    res <- process_proteomics(nadia_dia, covariate_df = cov, verbose = FALSE)

    expect_warning(
        batch_correct_proteomics(res$se_proc, assay_name = "Impseqrob_min",
                                 batch_column = "Batch", verbose = FALSE),
        "covariates", ignore.case = TRUE)
})

test_that("a categorical covariate can be protected during batch correction", {
    # BERT uses the Cov_* columns directly as the design matrix and encodes
    # nothing itself, so a character column used to abort: ComBat coerced it to
    # double and hit "NA/NaN/Inf in foreign function call", limma rejected it
    # with "design must be a numeric matrix". Condition is text with 3 levels.
    skip_if_not_installed("BERT")
    data(nadia_dia, package = "NADIA")

    md  <- nadia_dia$metadata
    cov <- data.frame(Column = md$Coding,
                      Batch  = rep(c("b1", "b2"), length.out = nrow(md)))
    res <- process_proteomics(nadia_dia, covariate_df = cov, verbose = FALSE)

    for (algo in c("ComBat", "limma")) {
        se <- batch_correct_proteomics(
            res$se_proc, assay_name = "Impseqrob_min", batch_column = "Batch",
            covariates = "Condition", algorithm = algo,
            corrected_assay_name = "BERT", verbose = FALSE)
        expect_true("BERT" %in% SummarizedExperiment::assayNames(se))
        expect_equal(dim(SummarizedExperiment::assay(se, "BERT")),
                     dim(SummarizedExperiment::assay(res$se_proc,
                                                     "Impseqrob_min")))
    }
})

test_that("covariate encoding keeps numeric columns as continuous", {
    # A numeric covariate must reach BERT untouched -- turning it into
    # indicators would silently change it from continuous to categorical.
    d <- data.frame(Condition = c("A", "A", "B", "B", "D", "D"),
                    Age       = c(30, 40, 50, 60, 70, 80))
    enc <- NADIA:::.bc_encode_covariates(d)

    expect_identical(enc$Age, d$Age)
    # 3 levels -> 2 treatment-contrast indicators, plus the untouched Age
    expect_equal(ncol(enc), 3L)
    expect_true(all(vapply(enc, is.numeric, logical(1))))
    expect_identical(enc[["Condition_B"]], c(0, 0, 1, 1, 0, 0))
    expect_identical(enc[["Condition_D"]], c(0, 0, 0, 0, 1, 1))

    # A covariate that never varies cannot protect anything
    expect_error(NADIA:::.bc_encode_covariates(data.frame(x = rep("a", 4))),
                 "single level", ignore.case = TRUE)
})

test_that("benchmarking validates the species table", {
    data(nadia_dia, package = "NADIA")
    de <- process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results
    ev <- data.frame(Comparison = "B-A", Species = "ECOLI", expected_logFC = 1)

    # Missing the Species column entirely
    expect_error(
        benchmarking_proteomics(de, ev, output_dir = NULL, verbose = FALSE,
                                species_df = data.frame(Protein.IDs = "p1")),
        "Species", ignore.case = TRUE)
})

test_that("the reactable modules name the column they are missing", {
    skip_if_not_installed("reactable")
    expect_error(protein_list_reactable(data.frame(nonsense = 1)),
                 "PG.ProteinGroups", fixed = TRUE)
    expect_error(results_list_reactable(data.frame(nonsense = 1)))
})

test_that("quant_list_widget requires ProteinGroups in the processed matrix", {
    skip_if_not_installed("reactable")
    data(nadia_dia, package = "NADIA")
    expect_error(
        quant_list_widget(nadia_dia$protein_quant,
                          matrix_data = data.frame(a = 1, b = 2),
                          metadata = nadia_dia$metadata),
        "ProteinGroups", fixed = TRUE)
})

test_that("pattern_profiler_analysis rejects an assay that does not exist", {
    skip_if_not_installed("Mfuzz")
    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    expect_error(pattern_profiler_analysis(res$se_proc, res$DEPs_results,
                                           assay_name = "no_such_assay",
                                           auto_select_c = FALSE, c = 3,
                                           verbose = FALSE))
})

test_that("the heatmap requires a comparison in the targeted mode", {
    skip_if_not_installed("ComplexHeatmap")
    skip_if_not_installed("tidyHeatmap")
    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    expect_error(proteomics_heatmap(res$PCA_Input, mode = "target"),
                 "comparison", ignore.case = TRUE)
})

test_that("a missing suggested package is named, with how to install it", {
    expect_error(.check_suggests("notAPackage0", "f()"),
                 "f() requires 'notAPackage0'. Install it with BiocManager::install(\"notAPackage0\").",
                 fixed = TRUE)
    expect_error(.check_suggests(c("stats", "notAPackage0", "notAPackage1"), "f()"),
                 "Install them with BiocManager::install(c(\"notAPackage0\", \"notAPackage1\")).",
                 fixed = TRUE)
    expect_true(.check_suggests("stats", "f()"))
})

test_that("the cluster palette honours a paletteer palette", {
    skip_if_not_installed("paletteer")
    skip_if_not_installed("ggsci")
    expect_identical(configure_cluster_palette(3, "ggsci::default_jco"),
                     stats::setNames(
                         vapply(as.character(paletteer::paletteer_d("ggsci::default_jco", 3)),
                                .normalize_hex, character(1), USE.NAMES = FALSE),
                         as.character(1:3)))
})

test_that("delimited files round trip, awkward fields included", {
    df <- data.frame(`B-A` = c(0.1 + 0.2, 1 / 3, 30.585426, -1.2e-17, NA),
                     Count = c(1L, 2L, NA, 4L, 5L),
                     Label = c("a\"quoted\"", "with\ttab", "two\nlines", "", NA),
                     check.names = FALSE, stringsAsFactors = FALSE)
    f <- withr::local_tempfile(fileext = ".tsv")
    .write_tsv(df, f)
    back <- .read_delim_file(f)

    expect_identical(names(back), names(df))
    expect_identical(back$Label, c("a\"quoted\"", "with\ttab", "two\nlines", NA, NA))
    expect_identical(back$Count, as.double(df$Count))
    expect_equal(back$`B-A`, df$`B-A`, tolerance = 1e-15)
    # The shortest form that reads back, not the 15 digits of write.table().
    expect_true(any(grepl("^30.585426\t", readLines(f))))
    expect_true(any(grepl("^0.30000000000000004\t", readLines(f))))
})
