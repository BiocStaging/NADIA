# =============================================================================
# process_proteomics(): the branches the end-to-end tests do not reach
#
# The other test files call process_proteomics() some fifty times, quietly and
# with the default imputation. What was left was the argument checks, the
# batch-corrected export, imputation "none", the progress messages and the
# print method.
# =============================================================================

test_that("an input that is neither a NADIA object nor an SE is refused", {
    expect_error(process_proteomics(list(a = 1), verbose = FALSE),
                 "must be a 'proteomics_data' object")
    expect_error(process_proteomics(data.frame(x = 1:3), verbose = FALSE),
                 "SummarizedExperiment")
})

test_that("a Parquet export without arrow fails before anything is computed", {
    data(nadia_dia, package = "NADIA")
    out <- withr::local_tempdir()
    # Whether arrow is installed or not, the check must run first
    local_mocked_bindings(.check_suggests = function(pkgs, what) stop("needs ", pkgs))
    expect_error(process_proteomics(nadia_dia, export_dir = file.path(out, "res"),
                                    export_format = "parquet", verbose = FALSE),
                 "needs arrow")
    expect_false(dir.exists(file.path(out, "res")))
})

test_that("batch correction without the batch column says which one is missing", {
    data(nadia_dia, package = "NADIA")
    expect_error(process_proteomics(nadia_dia, batch_correct = TRUE, verbose = FALSE),
                 "batch column 'Batch' not found")
})

test_that("the batch-corrected matrix is exported next to the normalised one", {
    skip_if_not_installed("BERT")
    data(nadia_dia, package = "NADIA")
    md  <- nadia_dia$metadata
    cov <- data.frame(Column = md$Coding,
                      Batch  = ifelse(md$R.Replicate <= 2, "b1", "b2"))
    out <- withr::local_tempdir()
    res <- suppressWarnings(process_proteomics(
        nadia_dia, covariate_df = cov, batch_correct = TRUE, export_dir = out,
        export_volcano = FALSE, export_boxplot = FALSE, export_pca = FALSE,
        verbose = FALSE))
    f <- file.path(out, "matrix_log2_cycloess_BERT.tsv")
    expect_true(file.exists(f))
    bert <- utils::read.delim(f, check.names = FALSE)
    corrected <- SummarizedExperiment::assay(res$se_proc, "BERT")
    expect_identical(bert$ProteinGroups, rownames(corrected))
    expect_equal(unname(as.matrix(bert[, -1])), unname(corrected), tolerance = 1e-12)
})

test_that("imp_method = 'none' keeps the missing values, and says so", {
    data(nadia_dia, package = "NADIA")
    out <- withr::local_tempdir()
    msgs <- character()
    res <- suppressWarnings(withCallingHandlers(
        process_proteomics(nadia_dia, imp_method = "none", export_dir = out,
                           verbose = TRUE),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        }))
    # The "none" assay is the normalised one with its gaps; limma analyses it
    expect_identical(unique(res$DEPs_results$Assay), "none")
    imp  <- SummarizedExperiment::assay(res$se_proc, "none")
    norm <- SummarizedExperiment::assay(res$se_proc, "cycloess")
    expect_true(anyNA(imp))
    expect_identical(unname(is.na(imp)), unname(is.na(norm)))
    expect_true(all(file.exists(file.path(out, c(
        "matrix_log2_cycloess.tsv", "matrix_log2_cycloess_none.tsv",
        "VolcanoPlot_Input_cycloess_none.tsv", "BoxPlot_Input_cycloess_none.tsv",
        "PCA_Input_cycloess_none.tsv")))))
    # verbose = TRUE reports every stage
    expect_true(any(grepl("PREPARING DATA", msgs)))
    expect_true(any(grepl("EXPORTING FILES FOR VISUALIZATION", msgs)))
    expect_true(any(grepl("PROCESSING COMPLETED", msgs)))

    expect_output(print(res), "Proteomics Processing Result")
    expect_output(print(res), "Imputation: none")
})
