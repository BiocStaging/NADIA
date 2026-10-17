# =============================================================================
# The four *_widget() tables and inst/js/nadia-tables/nadia-tables.js
#
# The buttons, chips and Excel export live in one script; each widget wraps
# its controls in an element marked with the table id and registers its
# settings. What can break without anything failing in R is that wiring:
# the JSON the script receives, or two tables on a page sharing ids.
# =============================================================================

# The JSON passed to NadiaTables.register() for the table `id`
.registered <- function(widget, id) {
    html <- as.character(htmltools::renderTags(widget)$html)
    m <- regmatches(html, regexpr(
        sprintf('NadiaTables\\.register\\("%s", .*?\\);', id), html))
    expect_length(m, 1)
    jsonlite::fromJSON(sub('^[^,]+, (.*)\\);$', "\\1", m), simplifyVector = FALSE)
}

test_that("every widget wraps its controls and registers with the shared script", {
    skip_if_not_installed("reactable")
    skip_if_not_installed("crosstalk")
    data(nadia_dia, package = "NADIA")
    meta <- nadia_dia$metadata
    meta$Coding <- paste(meta$R.Condition, meta$R.Replicate, sep = "_")
    res <- process_proteomics(nadia_dia, verbose = FALSE)
    widgets <- list(
        deps_table    = results_list_widget(head(res$DEPs_results, 50)),
        summary_table = summary_list_widget(meta))
    for (id in names(widgets)) {
        html <- as.character(htmltools::renderTags(widgets[[id]])$html)
        expect_match(html, sprintf('data-nadia-table="%s"', id), fixed = TRUE)
        # (htmltools writes the quotes of an attribute as &#39;)
        expect_match(html, sprintf("NadiaTables\\.toggleFilters\\((&#39;|')%s(&#39;|')\\)", id))
        # Nothing left of the per-widget scripts
        expect_no_match(html, "(rl|pl|ql|sl)(ToggleFilters|ExportExcel)")
        expect_type(.registered(widgets[[id]], id)$excel$layout, "character")
    }
})

test_that("one-element vectors reach the script as arrays", {
    skip_if_not_installed("reactable")
    # One sample column per condition and a single annotation column: every
    # vector the script iterates over has length one
    df <- data.frame(PG.ProteinGroups = c("P1", "P2"),
                     PG.NrOfPrecursorsIdentified_A_1 = c(3, 5),
                     PG.NrOfPrecursorsIdentified_B_1 = c(2, NA))
    cfg <- .registered(protein_list_widget(df, element_id = "p1"), "p1")
    expect_identical(cfg$conditionMode, "columns")
    expect_identical(cfg$conditions$A, list("PG.NrOfPrecursorsIdentified_A_1"))
    expect_identical(cfg$excel$groups[[1]]$columns, list("PG.ProteinGroups"))
    expect_identical(cfg$excel$groups[[1]]$headers, list("Protein Groups"))
    expect_identical(unlist(cfg$excel$widths), c(22L, 8L, 8L))
    # The chip label counts the columns a condition hides
    html <- as.character(htmltools::renderTags(protein_list_widget(df))$html)
    expect_match(html, "click to hide/show 1 cols", fixed = TRUE)
})

test_that("two tables of one kind on a page do not share ids", {
    skip_if_not_installed("reactable")
    skip_if_not_installed("crosstalk")
    data(nadia_dia, package = "NADIA")
    deps <- head(process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results, 50)
    html <- as.character(htmltools::renderTags(htmltools::tagList(
        results_list_widget(deps, element_id = "r1"),
        results_list_widget(deps, element_id = "r2")))$html)
    ids <- regmatches(html, gregexpr('id="[^"]+"', html))[[1]]
    expect_identical(anyDuplicated(ids), 0L)
    expect_match(html, 'id="r2_filter_comparison"', fixed = TRUE)
})

test_that("the protein tables declare their annotation columns as text for Excel", {
    skip_if_not_installed("reactable")
    df <- data.frame(PG.ProteinGroups = c("P1", "P2"),
                     PG.ProteinDescriptions = c("60S ribosomal protein L1", "Actin"),
                     PG.Genes = c("RPL10A", "ACTB"),
                     PG.NrOfPrecursorsIdentified_A_1 = c(3, 5))
    cfg <- .registered(protein_list_widget(df, element_id = "p1"), "p1")
    expect_setequal(unlist(cfg$excel$text),
                    c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes"))
})

test_that("the Excel export keeps text that starts with digits", {
    # A description such as "60S ribosomal protein L1" used to be written as the
    # number 60. The export runs in JavaScript; node runs it with stand-ins for
    # ExcelJS and PapaParse.
    node <- Sys.which("node")
    skip_if(!nzchar(node), "node is not available")
    script <- system.file("js", "nadia-tables", "nadia-tables.js", package = "NADIA")
    harness <- test_path("excel-harness.js")
    skip_if(!nzchar(script) || !file.exists(harness))

    settings <- withr::local_tempfile(fileext = ".json")
    table <- withr::local_tempfile(fileext = ".tsv")
    writeLines(jsonlite::toJSON(list(conditionMode = "columns", excel = list(
        layout = "grouped", sheet = "s", file = "f", headerFill = "FF000000",
        groups = list(list(name = "A", columns = I(c("PG.ProteinGroups",
                                                     "PG.ProteinDescriptions",
                                                     "PG.NrOfPrecursorsIdentified_A_1")),
                           headers = I(c("Protein Groups", "Descriptions", "A_1")))),
        text = I(c("PG.ProteinGroups", "PG.ProteinDescriptions")),
        widths = I(c(22, 42, 8)), palette = list(A = "#123456"))),
        auto_unbox = TRUE), settings)
    writeLines(c("PG.ProteinGroups\tPG.ProteinDescriptions\tPG.NrOfPrecursorsIdentified_A_1",
                 "123ABC\t60S ribosomal protein L1\t3",
                 "00042\t7SK snRNA methylphosphate capping enzyme\t1.5e2"), table)

    out <- system2(node, c(shQuote(harness), shQuote(script), shQuote(settings),
                           shQuote(table)), stdout = TRUE, stderr = TRUE)
    rows <- jsonlite::fromJSON(paste(out, collapse = ""), simplifyVector = FALSE)
    expect_identical(rows[["3"]][["1"]], "123ABC")
    expect_identical(rows[["3"]][["2"]], "60S ribosomal protein L1")
    expect_identical(rows[["4"]][["1"]], "00042")          # leading zeros kept
    expect_identical(rows[["4"]][["2"]], "7SK snRNA methylphosphate capping enzyme")
    # The numeric columns are still numbers
    expect_identical(rows[["3"]][["3"]], 3L)
    expect_identical(rows[["4"]][["3"]], 150L)
})
