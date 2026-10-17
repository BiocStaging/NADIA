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
