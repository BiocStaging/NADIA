# =============================================================================
# The agent skill shipped in inst/skill
#
# nadia_skill_path() only locates it; what can go wrong is the skill itself
# drifting from the package: a renamed file it points to, or a function it
# tells an assistant to call that is no longer exported.
# =============================================================================

test_that("nadia_skill_path() finds the skill and its entry point", {
    dir <- nadia_skill_path()
    expect_true(dir.exists(dir))
    expect_true(file.exists(nadia_skill_path("SKILL.md")))
    expect_identical(nadia_skill_path("SKILL.md"), file.path(dir, "SKILL.md"))
    # A missing file gives "", as system.file() does
    expect_identical(nadia_skill_path("not_a_file.md"), "")
})

test_that("SKILL.md has the front matter a Claude skill needs", {
    lines <- readLines(nadia_skill_path("SKILL.md"), warn = FALSE)
    delim <- which(lines == "---")
    expect_identical(delim[1], 1L)
    expect_gte(length(delim), 2)
    front <- lines[2:(delim[2] - 1)]
    expect_true(any(grepl("^name: nadia$", front)))
    description <- sub("^description: ", "", grep("^description: ", front, value = TRUE))
    expect_length(description, 1)
    expect_gt(nchar(description), 50)
})

test_that("every file the skill points to exists", {
    dir <- nadia_skill_path()
    md <- list.files(dir, "[.]md$", recursive = TRUE, full.names = TRUE)
    refs <- unique(unlist(lapply(md, function(f) {
        txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
        regmatches(txt, gregexpr("(reference|scripts)/[A-Za-z0-9_.-]+[.](md|R)", txt))[[1]]
    })))
    expect_gt(length(refs), 0)
    missing <- refs[!file.exists(file.path(dir, refs))]
    expect_identical(missing, character(0))
})

test_that("the example script is valid R", {
    for (f in list.files(nadia_skill_path("scripts"), "[.]R$", full.names = TRUE))
        expect_no_error(parse(f))
})

test_that("every NADIA function the skill tells an assistant to call is exported", {
    dir <- nadia_skill_path()
    files <- list.files(dir, "[.](md|R)$", recursive = TRUE, full.names = TRUE)
    txt <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
    called <- unique(sub("[(]$", "", regmatches(
        txt, gregexpr("(?<![A-Za-z0-9_.:$])[A-Za-z][A-Za-z0-9_.]*[(]", txt, perl = TRUE))[[1]]))
    ns <- asNamespace("NADIA")
    ours <- called[vapply(called, function(f) exists(f, envir = ns, inherits = FALSE),
                          logical(1))]
    expect_gt(length(ours), 5)
    expect_identical(setdiff(ours, getNamespaceExports("NADIA")), character(0))
})
