# =============================================================================
# Side effects: the filesystem, the RNG and the search path
# =============================================================================
#
# A function that writes files, reseeds the RNG or attaches packages as a side
# effect of being called is a hazard in a script and unacceptable in a package.
# These tests assert that none of that happens.

# --- The filesystem ----------------------------------------------------------

#' Run `expr` inside an empty directory and report what it created
.files_created_by <- function(expr) {
    dir <- file.path(tempdir(), paste0("nadia_probe_", as.integer(runif(1, 1e6, 9e6))))
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)

    old <- setwd(dir)
    on.exit(setwd(old), add = TRUE)

    force(expr)
    list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
}

test_that("process_proteomics writes nothing without export_dir", {
    data(nadia_dia, package = "NADIA")
    created <- .files_created_by(process_proteomics(nadia_dia, verbose = FALSE))
    expect_identical(created, character(0))
})

test_that("the preprocessors write nothing without export_dir", {
    created <- .files_created_by(
        preprocess_spectronaut(
            system.file("extdata", "nadia_dia_report.tsv.gz", package = "NADIA"),
            condition_order = c("A", "B", "D"), verbose = FALSE))
    expect_identical(created, character(0))
})

test_that("the metrics orchestrators write nothing when output_dir = NULL", {
    # These three default to export_tables = TRUE, so output_dir = NULL is the
    # only thing that stops them writing. It is the trap most likely to bite a
    # user running them inside a loop.
    data(nadia_dia, package = "NADIA")

    se <- nm_run_normalizations(nm_prepare_se(nadia_dia, verbose = FALSE),
                                methods = c("log2Norm", "cycloess"),
                                verbose = FALSE)

    created <- .files_created_by(
        normalization_metrics(se, output_dir = NULL, verbose = FALSE))
    expect_identical(created, character(0))

    se_i <- im_prepare_se(nadia_dia, norm_method = "cycloess", verbose = FALSE)
    created <- .files_created_by(
        imputation_metrics(se_i, assay_name = "cycloess",
                           methods = c("min", "MinDet"),
                           output_dir = NULL, verbose = FALSE))
    expect_identical(created, character(0))
})

test_that("process_proteomics writes to export_dir when asked", {
    # The other half of the contract: it must actually export when told to.
    data(nadia_dia, package = "NADIA")
    out <- file.path(tempdir(), "nadia_export_probe")
    on.exit(unlink(out, recursive = TRUE), add = TRUE)

    invisible(process_proteomics(nadia_dia, export_dir = out,
                                 export_format = "tsv", verbose = FALSE))
    expect_gt(length(list.files(out)), 0)
})


# --- The random number generator ---------------------------------------------
#
# NADIA never seeds the RNG itself (Bioconductor asks packages not to call
# set.seed()). The functions that need random numbers draw them from the
# caller's RNG: the same seed must give the same result, and a different seed a
# different one -- which is what shows that no seed of their own overrides it.

test_that("imputation_metrics follows the caller's seed", {
    data(nadia_dia, package = "NADIA")
    se_i <- im_prepare_se(nadia_dia, norm_method = "cycloess", verbose = FALSE)
    run <- function(seed) {
        set.seed(seed)
        imputation_metrics(se_i, assay_name = "cycloess", methods = "MinDet",
                           output_dir = NULL, verbose = FALSE)$metrics_table
    }
    a <- run(7)
    expect_equal(a, run(7))
    expect_false(isTRUE(all.equal(a, run(8))))
})

test_that("the stochastic imputation methods follow the caller's seed", {
    set.seed(1)
    x <- matrix(rnorm(60 * 6, mean = 20), nrow = 60)
    x[sample(length(x), 40)] <- NA

    # MLE is left out: norm's generator keeps state that rngseed() does not
    # reset, so a second call in the same session differs even with the same
    # seed. It is reproducible only across fresh sessions.
    methods <- "PI"
    if (requireNamespace("mice", quietly = TRUE)) methods <- c(methods, "mice")

    for (m in methods) {
        run <- function(seed) {
            set.seed(seed)
            suppressWarnings(NADIA:::.dispatch_imputation(x, m, list()))
        }
        a <- run(3)
        expect_identical(a, run(3), info = m)
        expect_false(identical(a, run(4)), info = m)
    }
})


# --- The search path ---------------------------------------------------------

test_that(".mfuzz_deps_attach/.detach leave the search path as they found it", {
    skip_if_not_installed("Mfuzz")
    skip_if_not_installed("Biobase")
    skip_if_not_installed("e1071")

    before <- search()
    attached <- NADIA:::.mfuzz_deps_attach()
    NADIA:::.mfuzz_deps_detach(attached)
    expect_identical(search(), before)
})

test_that(".mfuzz_deps_detach does not detach what the user had attached", {
    # If Biobase was already on the search path when the analysis started, it
    # must still be there afterwards: the helper returns only what IT attached,
    # and only that is detached.
    skip_if_not_installed("Mfuzz")
    skip_if_not_installed("Biobase")
    skip_if_not_installed("e1071")

    already <- "package:Biobase" %in% search()
    if (!already) {
        attachNamespace(asNamespace("Biobase"))
        on.exit(detach("package:Biobase", character.only = TRUE), add = TRUE)
    }

    attached <- NADIA:::.mfuzz_deps_attach()
    expect_false("Biobase" %in% attached)   # it did not attach it

    NADIA:::.mfuzz_deps_detach(attached)
    expect_true("package:Biobase" %in% search())
})

test_that("pattern_profiler_analysis restores the search path, even on failure", {
    skip_if_not_installed("Mfuzz")
    skip_if_not_installed("Biobase")
    skip_if_not_installed("e1071")

    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    before <- search()
    invisible(pattern_profiler_analysis(res$se_proc, res$DEPs_results,
                                        assay_name = "Impseqrob_min",
                                        auto_select_c = FALSE, c = 3,
                                        verbose = FALSE))
    expect_identical(search(), before)

    # The on.exit() must also fire when the function aborts.
    before <- search()
    expect_error(pattern_profiler_analysis(res$se_proc, res$DEPs_results,
                                           assay_name = "does_not_exist",
                                           auto_select_c = FALSE, c = 3,
                                           verbose = FALSE))
    expect_identical(search(), before)
})


# --- Other packages' options -------------------------------------------------

test_that("proteomics_heatmap and its print method leave ht_opt as they found it", {
    # ComplexHeatmap emits its messages when the heatmap is drawn, not when it
    # is built, so both calls silence them. Neither may leave them silenced.
    skip_if_not_installed("ComplexHeatmap")
    skip_if_not_installed("tidyHeatmap")

    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    old <- ComplexHeatmap::ht_opt("message")
    on.exit(ComplexHeatmap::ht_opt(message = old), add = TRUE)

    for (user_value in c(TRUE, FALSE)) {
        ComplexHeatmap::ht_opt(message = user_value)
        hm <- suppressWarnings(proteomics_heatmap(res$PCA_Input, mode = "any"))
        expect_identical(ComplexHeatmap::ht_opt("message"), user_value)

        f <- tempfile(fileext = ".png")
        grDevices::png(f)
        suppressWarnings(print(hm))
        grDevices::dev.off()
        unlink(f)
        expect_identical(ComplexHeatmap::ht_opt("message"), user_value)
    }
})
