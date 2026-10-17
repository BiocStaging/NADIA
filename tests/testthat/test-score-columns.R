# =============================================================================
# The protein score: neutral column names, the report's own label
#
# The readers used to put every format's score under Spectronaut's PG.Cscore,
# so a Proteome Discoverer Sum PEP Score read as a C-score (Bioconductor
# review, #186). The columns are now PG.Score and PG.Score.RunWise_*, and the
# metric of the report travels with the table as the nadia_labels attribute.
# =============================================================================

.ext <- function(pattern) {
    f <- list.files(system.file("extdata", package = "NADIA"), pattern, full.names = TRUE)
    f[!grepl("[Aa]nnot", f)][1]
}

test_that("every reader names its score PG.Score and says what it is", {
    data(nadia_dia, package = "NADIA")
    tmt <- preprocess_tmt(.ext("^nadia_tmt"), condition_order = c("A", "B", "D"),
                          verbose = FALSE)
    lfq <- preprocess_lfq(.ext("^nadia_lfq"), verbose = FALSE)
    expected <- list(
        dia = c(PG.Score.RunWise = "Cscore (run-wise)", PG.Score = "Cscore"),
        tmt = c(PG.Score.RunWise = "Sum PEP Score",     PG.Score = "Sum PEP Score"),
        lfq = c(PG.Score.RunWise = "Score Mascot",      PG.Score = "Sum PEP Score"))
    objs <- list(dia = nadia_dia, tmt = tmt, lfq = lfq)
    for (nm in names(objs)) {
        x <- objs[[nm]]
        both <- c(names(x$protein_id), names(x$protein_quant))
        expect_false(any(grepl("Cscore", both)), info = nm)
        expect_true(any(startsWith(names(x$protein_id), "PG.Score.RunWise_")), info = nm)
        expect_true("PG.Score" %in% names(x$protein_quant), info = nm)
        expect_identical(.score_label(x$protein_id, "PG.Score.RunWise"),
                         unname(expected[[nm]]["PG.Score.RunWise"]), info = nm)
        expect_identical(.score_label(x$protein_quant, "PG.Score"),
                         unname(expected[[nm]]["PG.Score"]), info = nm)
    }
})

test_that("the interactive tables label the score as the report does", {
    skip_if_not_installed("reactable")
    lfq <- preprocess_lfq(.ext("^nadia_lfq"), verbose = FALSE)
    html <- as.character(htmltools::renderTags(
        protein_list_widget(lfq$protein_id, metadata = lfq$metadata))$html)
    expect_match(html, "Score Mascot", fixed = TRUE)
    expect_no_match(html, "Cscore")
    # Without the attribute, as from an exported TSV, the label is neutral
    plain <- lfq$protein_id
    attr(plain, "nadia_labels") <- NULL
    html <- as.character(htmltools::renderTags(protein_list_widget(plain))$html)
    expect_no_match(html, "Score Mascot|Cscore")
})

test_that("a table written before the rename still loads", {
    data(nadia_dia, package = "NADIA")
    old <- nadia_dia$protein_id
    names(old) <- sub("^PG\\.Score", "PG.Cscore", names(old))
    attr(old, "nadia_labels") <- NULL
    df <- .pl_load_data(old)
    expect_true(any(startsWith(names(df), "PG.Score.RunWise_")))
    expect_identical(.score_label(df, "PG.Score.RunWise"), "Score")
})

test_that("a .nadia file keeps the labels", {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")
    tmt <- preprocess_tmt(.ext("^nadia_tmt"), condition_order = c("A", "B", "D"),
                          verbose = FALSE)
    f <- withr::local_tempfile(fileext = ".nadia")
    write_nadia(f, tmt, verbose = FALSE)
    back <- nadia_preprocessing(read_nadia(f))
    expect_identical(attr(back$protein_quant, "nadia_labels"),
                     attr(tmt$protein_quant, "nadia_labels"))
    expect_identical(back$protein_id, tmt$protein_id)
})
