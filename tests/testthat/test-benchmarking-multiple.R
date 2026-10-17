# =============================================================================
# benchmarking_multiple(): the OpDEA ranking of several pipelines
#
# A synthetic table with a known answer checks the ranking arithmetic, and two
# real pipelines on the spike-in example check that the in-memory and the
# file-based routes, which users mix freely, give the same ranking.
# =============================================================================

# Three methods over two comparisons. "good" beats "mid", which beats "bad",
# on every metric; the values differ between comparisons so that the mean and
# the median are not trivially equal.
.bm_synthetic <- function() {
    data.frame(
        Assay      = rep(c("good", "mid", "bad"), each = 2),
        Comparison = rep(c("B-A", "D-A"), times = 3),
        nMCC     = c(0.90, 0.80, 0.70, 0.60, 0.50, 0.40),
        G_mean   = c(0.95, 0.85, 0.75, 0.65, 0.55, 0.45),
        pAUC_001 = c(0.60, 0.50, 0.40, 0.30, 0.20, 0.10),
        pAUC_005 = c(0.80, 0.70, 0.60, 0.50, 0.40, 0.30),
        pAUC_010 = c(0.90, 0.85, 0.70, 0.65, 0.50, 0.45),
        stringsAsFactors = FALSE)
}

test_that("a method that wins every metric ranks first everywhere", {
    rk <- bm_compute_ranking(.bm_synthetic())

    for (r in list(rk$mean_ranking, rk$median_ranking)) {
        expect_identical(r$Assay, c("good", "mid", "bad"))
        rank_cols <- grep("^rank_", names(r), value = TRUE)
        for (rc in setdiff(rank_cols, "rank_final"))
            expect_equal(sort(r[[rc]]), c(1, 2, 3))
        expect_equal(r$rank_final, c(1, 2, 3))
    }
    expect_identical(rk$n_methods, 3L)
    expect_identical(rk$n_comparisons, 2L)
})

test_that("the aggregation and the final rank are what they say they are", {
    df <- .bm_synthetic()
    rk <- bm_compute_ranking(df)

    for (agg in c("mean", "median")) {
        f <- match.fun(agg)
        tab <- rk[[paste0(agg, "_aggregated")]]
        for (m in .BM_METRICS) {
            manual <- tapply(df[[m]], df$Assay, f)
            expect_equal(tab[[m]], round(as.vector(manual[tab$Assay]), 4),
                         info = paste(agg, m))
        }
        # rank_final is the mean of the per-metric ranks
        r <- rk[[paste0(agg, "_ranking")]]
        expect_equal(r$rank_final,
                     round(rowMeans(r[, paste0("rank_", .BM_METRICS)]), 2))
    }
})

test_that("ties share the average rank, and a missing value ranks last", {
    df <- .bm_synthetic()
    df[df$Assay == "mid", .BM_METRICS] <- df[df$Assay == "good", .BM_METRICS]
    r <- bm_compute_ranking(df)$mean_ranking
    expect_equal(r$rank_nMCC[r$Assay %in% c("good", "mid")], c(1.5, 1.5))

    df <- .bm_synthetic()
    df$nMCC[df$Assay == "good"] <- NA
    r <- bm_compute_ranking(df)$mean_ranking
    expect_equal(r$rank_nMCC[r$Assay == "good"], 3)
})

test_that("the ranking by comparison ranks each comparison on its own", {
    df <- .bm_synthetic()
    df$nMCC[df$Assay == "bad" & df$Comparison == "D-A"] <- 0.99  # bad wins D-A
    by_comp <- bm_compute_ranking_by_comparison(df)

    expect_setequal(names(by_comp$by_comparison), c("B-A", "D-A"))
    for (comp in c("B-A", "D-A")) {
        alone <- bm_compute_ranking(df[df$Comparison == comp, ])$mean_ranking
        got <- by_comp$by_comparison[[comp]]
        expect_identical(got$Assay, alone$Assay, info = comp)
        expect_equal(got$rank_nMCC, alone$rank_nMCC, info = comp)
    }
    da <- by_comp$by_comparison[["D-A"]]
    expect_equal(da$rank_nMCC[da$Assay == "bad"], 1)
})

test_that("the input is validated with informative errors", {
    df <- .bm_synthetic()
    expect_error(bm_compute_ranking(df[, setdiff(names(df), "Comparison")]),
                 "Missing required columns")
    expect_error(bm_compute_ranking(df[0, ]), "0 rows")
    expect_error(bm_compute_ranking(df[, setdiff(names(df), "pAUC_001")]),
                 "pAUC_001")
    # Extra metrics that are absent are skipped with a warning
    expect_warning(r <- bm_compute_ranking(df, metrics = c(.BM_METRICS, "F1")),
                   "F1")
    expect_false("rank_F1" %in% names(r$mean_ranking))
    expect_error(benchmarking_multiple(verbose = FALSE),
                 "opdea_combined or results_dir")
})

test_that("the extended ranking converts Performance to a number", {
    # The labels are those benchmarking_proteomics() assigns from F1
    expect_equal(.bm_performance_to_numeric(
        c("Excellent", "Very Good", "Good", "Acceptable", "Poor")), c(5, 4, 3, 2, 1))
})


# --- Two real pipelines on the spike-in example ------------------------------

.bm_real <- local({
    cache <- NULL
    function() {
        if (!is.null(cache)) return(cache)
        data(nadia_dia, package = "NADIA", envir = environment())
        sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
                                            package = "NADIA"))
        species_df <- data.frame(Protein.IDs = sp$PG.ProteinGroups,
                                 Species = sp$PG.OrganismId)
        expected <- data.frame(
            Comparison = rep(c("B-A", "D-A", "D-B"), each = 2),
            Species = rep(c("ECOLI", "YEAST"), times = 3),
            expected_logFC = c(1, -0.58, 2, -3.3, 1, -2.72))
        pipelines <- list(cycloess_combo = list(imp_method = "combo"),
                          cycloess_min = list(imp_method = "min"))
        bench <- lapply(pipelines, function(args) {
            res <- do.call(process_proteomics,
                           c(list(nadia_dia, verbose = FALSE), args))
            benchmarking_proteomics(res$DEPs_results, expected,
                                    species_df = species_df,
                                    output_dir = NULL, verbose = FALSE)
        })
        cache <<- list(bench = bench, expected = expected,
                       species_df = species_df, nadia_dia = nadia_dia,
                       pipelines = pipelines)
        cache
    }
})

# Stack one element of every pipeline's benchmark, naming the rows after the
# pipeline (the classified table carries the assay name; it is replaced).
.bm_stack <- function(bench, element) {
    do.call(rbind, lapply(names(bench), function(p) {
        x <- as.data.frame(bench[[p]][[element]])
        x$Assay <- NULL
        cbind(Assay = p, x, stringsAsFactors = FALSE)
    }))
}

test_that("benchmarking_multiple ranks real pipelines and draws every plot", {
    skip_if_not_installed("pROC")
    fx <- .bm_real()
    out <- suppressWarnings(benchmarking_multiple(
        opdea_combined         = .bm_stack(fx$bench, "opdea_metrics"),
        confusion_combined     = .bm_stack(fx$bench, "confusion_overall"),
        classified_combined    = .bm_stack(fx$bench, "classified_df"),
        bench_metrics_combined = .bm_stack(fx$bench, "metrics_table"),
        verbose = FALSE))

    expect_setequal(out$mean_ranking$Assay, names(fx$bench))
    # The ranking agrees with the metric values it was built from
    agg <- out$mean_aggregated
    best_nmcc <- agg$Assay[which.max(agg$nMCC)]
    expect_equal(out$mean_ranking$rank_nMCC[out$mean_ranking$Assay == best_nmcc], 1)
    expect_setequal(names(out$ranking_by_comparison), c("B-A", "D-A", "D-B"))

    # Every plot it returns is a ggplot that builds
    plots <- c(out[grep("^gg_", names(out))],
               unlist(out[grep("^gg_.*_by_comp$", names(out))], recursive = FALSE))
    plots <- Filter(function(p) inherits(p, "ggplot"), plots)
    expect_gte(length(plots), 10)
    for (p in plots) expect_no_error(ggplot2::ggplot_build(p))
})

test_that("the file-based route gives the same ranking as the in-memory one", {
    skip_if_not_installed("pROC")
    fx <- .bm_real()
    results_dir <- withr::local_tempdir()
    output_dir <- withr::local_tempdir()

    # Each pipeline is benchmarked into its own folder, as the vignette shows
    for (p in names(fx$pipelines)) {
        res <- do.call(process_proteomics,
                       c(list(fx$nadia_dia, verbose = FALSE), fx$pipelines[[p]]))
        suppressWarnings(benchmarking_proteomics(
            res$DEPs_results, fx$expected, species_df = fx$species_df,
            output_dir = file.path(results_dir, paste0("benchmark_", p)),
            verbose = FALSE))
    }
    from_files <- suppressWarnings(benchmarking_multiple(
        results_dir = results_dir, output_dir = output_dir,
        plots = "ranking_bars_mean", plot_width = 4, plot_height = 3,
        plot_dpi = 50, verbose = FALSE))
    in_memory <- suppressWarnings(benchmarking_multiple(
        opdea_combined = .bm_stack(fx$bench, "opdea_metrics"),
        plots = "ranking_bars_mean", verbose = FALSE))

    expect_equal(from_files$mean_ranking, in_memory$mean_ranking)
    expect_equal(from_files$median_ranking, in_memory$median_ranking)

    # The exported tables read back as the returned ones
    tsv <- list.files(output_dir, "[.]tsv$", recursive = TRUE)
    expect_true("bm_multiple_mean_ranking.tsv" %in% tsv)
    back <- utils::read.delim(file.path(output_dir, "bm_multiple_mean_ranking.tsv"),
                              check.names = FALSE)
    expect_equal(back$Assay, from_files$mean_ranking$Assay)
    expect_equal(back$rank_final, from_files$mean_ranking$rank_final)
    expect_true(any(grepl("[.]png$", list.files(output_dir, recursive = TRUE))))
})
