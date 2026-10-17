# =============================================================================
# Batch correction and its diagnostics, on data with a known structure
#
# 300 proteins x 12 samples, two conditions crossed with two batches (three
# samples of each condition in each batch). Batch b2 shifts every protein by
# about 3 log2 units; condition B raises the first 60 proteins by 2. About 5%
# of the values are missing. Each test checks a property this design implies.
# =============================================================================

.bc_simulated <- function(seed = 1) {
    withr::local_seed(seed)
    n <- 300
    cd <- S4Vectors::DataFrame(
        Column    = paste0("s", 1:12),
        Condition = rep(c("A", "B"), each = 6),
        Batch     = rep(rep(c("b1", "b2"), each = 3), times = 2),
        row.names = paste0("s", 1:12))
    base <- stats::rnorm(n, 20, 2)
    batch_shift <- stats::rnorm(n, 3, 0.5)
    x <- outer(base, rep(1, 12)) +
        outer(c(rep(2, 60), rep(0, n - 60)), as.numeric(cd$Condition == "B")) +
        outer(batch_shift, as.numeric(cd$Batch == "b2")) +
        matrix(stats::rnorm(n * 12, 0, 0.3), n)
    x[sample(length(x), round(0.05 * length(x)))] <- NA
    dimnames(x) <- list(paste0("P", seq_len(n)), cd$Column)
    SummarizedExperiment::SummarizedExperiment(assays = list(log2 = x),
                                               colData = cd)
}

# Mean difference b2 - b1 within each condition, averaged over proteins
.bc_batch_gap <- function(x, cd) {
    gaps <- vapply(c("A", "B"), function(cond) {
        in_c <- cd$Condition == cond
        rowMeans(x[, in_c & cd$Batch == "b2"], na.rm = TRUE) -
            rowMeans(x[, in_c & cd$Batch == "b1"], na.rm = TRUE)
    }, numeric(nrow(x)))
    mean(abs(gaps), na.rm = TRUE)
}

test_that("PVCA attributes the variance to the factor that carries it", {
    skip_if_not_installed("lme4")
    se <- .bc_simulated()
    pv <- suppressMessages(pvca_compute(se, assay_name = "log2",
                                        factors = c("Batch", "Condition"),
                                        verbose = FALSE))
    expect_true(all(c("label", "weights") %in% names(pv)))
    expect_equal(sum(pv$weights), 1, tolerance = 1e-6)
    expect_true(all(pv$weights >= 0))
    w <- stats::setNames(pv$weights, pv$label)
    # The batch shift reaches every protein; the condition only 60 of 300
    expect_gt(w[["Batch"]], w[["Condition"]])
    expect_gt(w[["Batch"]], 0.5)

    gg <- pvca_plot(pv, technical_factors = "Batch", biological_factors = "Condition")
    expect_s3_class(gg, "ggplot")
    expect_no_error(ggplot2::ggplot_build(gg))
})

test_that("pvca_analysis writes its table and plot only where it is told to", {
    skip_if_not_installed("lme4")
    se <- .bc_simulated()
    out_dir <- withr::local_tempdir()
    res <- suppressMessages(pvca_analysis(
        se, assay_name = "log2", technical_factors = "Batch",
        biological_factors = "Condition", output_dir = out_dir,
        plot_width = 4, plot_height = 3, plot_dpi = 50, verbose = FALSE))
    files <- list.files(out_dir)
    expect_true("pvca_variance_components.tsv" %in% files)
    expect_true(any(grepl("[.]png$", files)))
    back <- utils::read.delim(file.path(out_dir, "pvca_variance_components.tsv"))
    expect_equal(sum(back$weights), 1, tolerance = 1e-6)
})

test_that("ComBat removes the batch shift and keeps the protected condition", {
    skip_if_not_installed("BERT")
    se <- .bc_simulated()
    cd <- SummarizedExperiment::colData(se)
    before <- SummarizedExperiment::assay(se, "log2")
    expect_gt(.bc_batch_gap(before, cd), 2.5)

    se2 <- suppressMessages(suppressWarnings(batch_correct_proteomics(
        se, assay_name = "log2", batch_column = "Batch",
        covariates = "Condition", verbose = FALSE)))
    after <- SummarizedExperiment::assay(se2, "BERT")

    expect_identical(dimnames(after), dimnames(before))
    # No measured value is lost, and none is filled in
    expect_identical(is.na(after), is.na(before))
    expect_lt(.bc_batch_gap(after, cd), 0.3)

    # The condition effect on the first 60 proteins survives the correction
    effect <- function(x) mean(rowMeans(x[1:60, cd$Condition == "B"], na.rm = TRUE) -
                               rowMeans(x[1:60, cd$Condition == "A"], na.rm = TRUE))
    expect_equal(effect(after), effect(before), tolerance = 0.15)
})

test_that("a protein that cannot protect the covariate is corrected for the batch only", {
    skip_if_not_installed("BERT")
    se <- .bc_simulated()
    cd <- SummarizedExperiment::colData(se)
    x <- SummarizedExperiment::assay(se, "log2")
    # Complete data, then a single value left in batch b2 x condition B for
    # P5 and P6: BERT needs two per batch-by-covariate cell to protect the
    # condition
    x[is.na(x)] <- 20
    in_cell <- cd$Batch == "b2" & cd$Condition == "B"
    x["P5", in_cell] <- c(25, NA, NA)
    x["P6", in_cell] <- c(NA, 18, NA)
    SummarizedExperiment::assay(se, "log2") <- x

    se2 <- suppressMessages(suppressWarnings(batch_correct_proteomics(
        se, assay_name = "log2", batch_column = "Batch",
        covariates = "Condition", verbose = FALSE)))
    after <- SummarizedExperiment::assay(se2, "BERT")
    meta <- S4Vectors::metadata(se2)$batch_correction

    expect_identical(meta$batch_only, c("P5", "P6"))
    expect_identical(meta$unadjusted, character(0))
    # Their measured values are all kept (BERT alone would set them to NA) ...
    expect_identical(is.na(after[c("P5", "P6"), ]), is.na(x[c("P5", "P6"), ]))
    # ... and corrected: the batch gap within condition A at least halves. Not
    # more: with a single, atypical value in one cell, the batch effect of such
    # a protein is estimated from unbalanced data, which is the cost of
    # correcting it without protecting the condition.
    a <- cd$Condition == "A"
    gap <- function(v) mean(v[a & cd$Batch == "b2"]) - mean(v[a & cd$Batch == "b1"])
    for (p in c("P5", "P6"))
        expect_lt(abs(gap(after[p, ])), abs(gap(x[p, ])) / 2, label = p)
    expect_message(suppressWarnings(batch_correct_proteomics(
        se, assay_name = "log2", batch_column = "Batch", covariates = "Condition")),
        "corrected for the batch only")

    # BERT needs two features per call, so a lone such protein is kept
    # unadjusted rather than lost
    x["P6", ] <- SummarizedExperiment::assay(.bc_simulated(), "log2")["P6", ]
    x["P6", is.na(x["P6", ])] <- 20
    SummarizedExperiment::assay(se, "log2") <- x
    se3 <- suppressMessages(suppressWarnings(batch_correct_proteomics(
        se, assay_name = "log2", batch_column = "Batch",
        covariates = "Condition", verbose = FALSE)))
    meta3 <- S4Vectors::metadata(se3)$batch_correction
    expect_identical(meta3$unadjusted, "P5")
    expect_identical(SummarizedExperiment::assay(se3, "BERT")["P5", ], x["P5", ])
})

test_that("a protein ComBat cannot fit is kept, unadjusted and in its place", {
    skip_if_not_installed("BERT")
    se <- .bc_simulated()
    cd <- SummarizedExperiment::colData(se)
    x <- SummarizedExperiment::assay(se, "log2")
    x["P100", cd$Batch == "b1"] <- 21        # no variance within batch b1
    SummarizedExperiment::assay(se, "log2") <- x

    unfit <- .bc_combat_unfittable(x, cd$Batch)
    expect_identical(names(which(unfit)), "P100")

    se2 <- suppressMessages(suppressWarnings(batch_correct_proteomics(
        se, assay_name = "log2", batch_column = "Batch",
        covariates = "Condition", verbose = FALSE)))
    after <- SummarizedExperiment::assay(se2, "BERT")
    expect_identical(rownames(after), rownames(x))
    expect_identical(after["P100", ], x["P100", ])
    expect_true("P100" %in% S4Vectors::metadata(se2)$batch_correction$unadjusted)
    expect_false(isTRUE(all.equal(after["P1", ], x["P1", ])))
})

test_that("an incomplete batch annotation is rejected before any correction", {
    se <- .bc_simulated()
    se$Batch[2] <- NA
    expect_error(batch_correct_proteomics(se, assay_name = "log2",
                                          batch_column = "Batch", verbose = FALSE),
                 "NA", ignore.case = TRUE)
    expect_error(batch_correct_proteomics(.bc_simulated(), assay_name = "log2",
                                          batch_column = "Lot", verbose = FALSE),
                 "Lot")
})

test_that("the PCA coloured by covariate builds, and separates the batches", {
    se <- .bc_simulated()
    p <- suppressMessages(pca_covariates_plot(se, assay_name = "log2",
                                              covariates = c("Batch", "Condition"),
                                              verbose = FALSE))
    plots <- if (inherits(p, "ggplot")) list(p) else Filter(function(x) inherits(x, "ggplot"),
                                                             unlist(list(p), recursive = FALSE))
    expect_gte(length(plots), 1)
    for (g in plots) expect_no_error(ggplot2::ggplot_build(g))

    # With a 3-unit shift on every protein, PC1 is the batch
    scores <- .pca_compute_scores(se, assay_name = "log2", verbose = FALSE)
    sc <- scores$scores
    expect_gt(abs(stats::cor(sc$PC1, as.numeric(sc$Batch == "b2"))), 0.95)
})
