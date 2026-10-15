# =============================================================================
# Guards for the conventions established in phase 6
#
# These are source-level checks rather than behavioural ones: they are what
# stops the package from drifting back to `cat()`, to untranslated strings, or
# to default arguments pointing at files no user has. Each one corresponds to a
# fix that had already been made once and was silently undone or missed before.
# =============================================================================

.nadia_r_files <- function() {
    dir <- system.file("R", package = "NADIA")
    # Under devtools::load_all() the sources are still on disk; under an
    # installed package they are not, and the check is skipped.
    src <- file.path(dirname(system.file("DESCRIPTION", package = "NADIA")), "R")
    if (!dir.exists(src)) return(character(0))
    list.files(src, pattern = "[.]R$", full.names = TRUE)
}

# The five print.* show methods are the only place cat() belongs: there the
# output IS the return value of the method, and message() would send it to
# stderr where print() output does not go.
.NADIA_PRINT_METHODS <- c("print.spectronaut_data", "print.diann_data",
                          "print.tmt_data", "print.lfq_data",
                          "print.proteomics_result", "print.nadia_db")

test_that("cat() survives only inside the print.* show methods", {
    files <- .nadia_r_files()
    skip_if(length(files) == 0, "sources not available (installed package)")

    offenders <- list()
    for (f in files) {
        lines <- readLines(f, warn = FALSE)
        # Drop the bodies of the print.* methods: from their definition to the
        # closing brace in column 1.
        in_print <- rep(FALSE, length(lines))
        starts <- grep(paste0("^(", paste(.NADIA_PRINT_METHODS, collapse = "|"),
                              ")\\s*<-\\s*function"), lines)
        for (s in starts) {
            e <- which(lines == "}" & seq_along(lines) > s)[1]
            if (!is.na(e)) in_print[s:e] <- TRUE
        }
        # `.concat(` in the embedded JavaScript is not a cat() call, and nor
        # is a comment or roxygen line that mentions one.
        hits <- grep("(^|[^a-zA-Z._])cat\\(", lines)
        hits <- hits[!in_print[hits]]
        hits <- hits[!grepl("^\\s*#", lines[hits])]
        hits <- hits[!grepl("concat\\(", lines[hits])]
        if (length(hits))
            offenders[[basename(f)]] <- paste0(basename(f), ":", hits)
    }
    expect_equal(as.character(unlist(offenders, use.names = FALSE)),
                 character(0))
})

test_that("no Spanish is left in the user-facing strings", {
    files <- .nadia_r_files()
    skip_if(length(files) == 0, "sources not available (installed package)")

    # Words that are unambiguously Spanish and would not appear in English
    # prose. Phase 4 missed thirteen strings in Benchmarking_Single.R precisely
    # because none of them carried an accent, so the non-ASCII check passed.
    spanish <- paste0("\\b(clasificacion|metricas|significativas|barras|",
                      "curvas|volcanos|generados|generadas|creado|creada|",
                      "completado|completada|instalar|directorio|columnas|",
                      "requeridas|faltantes|proteinas|muestras|analisis|",
                      "fichero|archivo)\\b")

    offenders <- character(0)
    for (f in files) {
        lines <- readLines(f, warn = FALSE)
        hits <- grep(spanish, lines, ignore.case = TRUE, perl = TRUE)
        if (length(hits))
            offenders <- c(offenders, paste0(basename(f), ":", hits, "  ",
                                             trimws(lines[hits])))
    }
    expect_equal(offenders, character(0))
})

test_that("no exported default argument points at data-raw/ or results/", {
    exported <- getNamespaceExports("NADIA")
    offenders <- character(0)
    for (nm in exported) {
        f <- get(nm, envir = asNamespace("NADIA"))
        if (!is.function(f)) next
        defaults <- formals(f)
        for (arg in names(defaults)) {
            d <- defaults[[arg]]
            # Arguments without a default are the empty symbol, which cannot be
            # forced; skip them rather than evaluating.
            if (missing(d) || !is.character(d)) next
            if (length(d) == 1 &&
                grepl("^(\\./)?(data-raw|results)/", d))
                offenders <- c(offenders, paste0(nm, "(", arg, " = \"", d, "\")"))
        }
    }
    expect_equal(offenders, character(0))
})

test_that("the widgets that lost their defaults now say so", {
    skip_if_not_installed("reactable")
    expect_error(quant_list_widget(), "'data' and 'matrix_data' are required")
    expect_error(summary_list_widget(), "'data' is required")
})

# -----------------------------------------------------------------------------
# The RNG belongs to the caller
# -----------------------------------------------------------------------------
# Bioconductor asks packages not to call set.seed(). The functions that need
# random numbers draw them from the caller's RNG, so set.seed() before the call
# reproduces it, and a different seed gives a different answer.

test_that("no function takes a seed any more", {
    for (f in c("nm_compute_metrics", "nm_plot_metrics", "normalization_metrics",
                "im_compute_metrics", "imputation_metrics",
                "pattern_profiler_analysis", "select_optimal_clusters",
                "evaluate_cluster_range", "run_mfuzz_clustering")) {
        args <- names(formals(get(f, envir = asNamespace("NADIA"))))
        expect_false(any(c("seed", "seeds") %in% args), info = f)
    }
    # The repeated starts of the cluster-number selection are now a count.
    for (f in c("pattern_profiler_analysis", "select_optimal_clusters",
                "evaluate_cluster_range")) {
        expect_true("n_starts" %in% names(formals(get(f, envir = asNamespace("NADIA")))),
                    info = f)
    }
})

test_that("nm_compute_metrics() draws its random columns from the caller's RNG", {
    se <- .nadia_toy_se(n_prot = 40, na_prop = 0)

    set.seed(1)
    a <- nm_compute_metrics(se)$Hopkins
    set.seed(1)
    b <- nm_compute_metrics(se)$Hopkins
    set.seed(99)
    d <- nm_compute_metrics(se)$Hopkins

    expect_equal(a, b)                       # same seed, same answer
    expect_false(isTRUE(all.equal(a, d)))    # no seed of its own overrides it
})

test_that("no package in Imports is checked with requireNamespace()", {
    files <- .nadia_r_files()
    skip_if(length(files) == 0, "sources not available (installed package)")

    # A package in Imports is installed with NADIA, so the check can never
    # fail and its fallback can never run (Bioconductor review, #186).
    desc <- read.dcf(file.path(dirname(dirname(files[1])), "DESCRIPTION"),
                     fields = "Imports")[1, ]
    imports <- trimws(sub("[(].*", "", strsplit(desc, ",")[[1]]))
    imports <- imports[nzchar(imports)]
    pattern <- sprintf('requireNamespace[(]["\'](%s)["\']',
                       paste(imports, collapse = "|"))
    offenders <- unlist(lapply(files, function(f) {
        lines <- readLines(f, warn = FALSE)
        hits <- grep(pattern, lines)
        if (length(hits)) paste0(basename(f), ":", hits)
    }))
    expect_equal(as.character(offenders), character(0))

    # The literal check above misses a loop such as
    # `for (pkg in c("ggplot2", ...)) requireNamespace(pkg)`, so a package name
    # held in a variable is only allowed inside the helpers that check Suggests.
    allowed <- c(".check_suggests", ".db_require", ".mfuzz_deps_attach")
    offenders <- unlist(lapply(files, function(f) {
        exprs <- parse(f, keep.source = TRUE)
        unlist(lapply(seq_along(exprs), function(k) {
            e <- exprs[[k]]
            is_fun <- is.call(e) && is.name(e[[1]]) &&
                as.character(e[[1]]) %in% c("<-", "=") &&
                is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
            if (!is_fun || as.character(e[[2]]) %in% allowed) return(NULL)
            src <- as.character(attr(exprs, "srcref")[[k]])
            src <- src[!grepl("^\\s*#", src)]
            hits <- grep("requireNamespace[(][^\"']|vapply[(][^)]*requireNamespace",
                         src)
            if (length(hits)) paste0(basename(f), ": ", as.character(e[[2]]))
        }))
    }))
    expect_equal(as.character(offenders), character(0))
})
