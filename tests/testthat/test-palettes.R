# =============================================================================
# .resolve_palette(): the one parser of palette specifications
#
# Every plotting function goes through it (Bioconductor review, #186: each
# module used to parse palettes on its own, with different results).
# =============================================================================

test_that("NULL gives the default, a vector or a function of n", {
    expect_identical(.resolve_palette(NULL, 2, default = c("#111111", "#222222")),
                     c("#111111", "#222222"))
    expect_identical(.resolve_palette(NULL, 3, default = function(n) grDevices::hcl.colors(n)),
                     grDevices::hcl.colors(3))
})

test_that("RColorBrewer palettes are accepted with and without the prefix", {
    expect_identical(.resolve_palette("brewer:Set1", 4),
                     RColorBrewer::brewer.pal(4, "Set1"))
    expect_identical(.resolve_palette("Set1", 4), .resolve_palette("brewer:Set1", 4))
    # Two colours are the first two of the smallest palette (three colours)
    expect_identical(.resolve_palette("brewer:Dark2", 2),
                     RColorBrewer::brewer.pal(3, "Dark2")[1:2])
})

test_that("colours are named after the levels, and a named palette must cover them", {
    pal <- .resolve_palette("brewer:Set1", 3, levels = c("A", "B", "D"))
    expect_identical(names(pal), c("A", "B", "D"))
    named <- c(D = "#333333", A = "#111111", B = "#222222")
    expect_identical(.resolve_palette(named, 3, levels = c("A", "B", "D")),
                     named[c("A", "B", "D")])
    expect_error(.resolve_palette(c(A = "#111111"), 2, levels = c("A", "B")),
                 "Missing colours for levels: B")
})

test_that("a palette with fewer colours than groups is interpolated, not recycled", {
    pal <- .resolve_palette(c("#000000", "#FFFFFF"), 3)
    expect_length(pal, 3)
    expect_length(unique(pal), 3)
    expect_identical(toupper(pal[c(1, 3)]), c("#000000", "#FFFFFF"))
    # A single colour fills every group
    expect_length(unique(.resolve_palette("red", 4)), 1)
})

test_that("a continuous palette keeps its own length", {
    ramp <- c("#000000", "#777777", "#FFFFFF")
    expect_identical(.resolve_palette(ramp, 11, type = "continuous"), ramp)
    expect_identical(.resolve_palette("brewer:RdBu", 7, type = "continuous"),
                     RColorBrewer::brewer.pal(7, "RdBu"))
    expect_identical(.resolve_palette(ramp, 11, type = "continuous", reverse = TRUE),
                     rev(ramp))
})

test_that("an unknown palette is an error, not a silent fallback", {
    expect_error(.resolve_palette("brewer:NoSuch", 3), "Unknown RColorBrewer palette")
    expect_error(.resolve_palette("NoSuchPalette", 3), "Unknown palette")
    expect_error(.resolve_palette(c("#111111", "notacolour"), 2), "Unknown palette")
})

test_that("paletteer palettes are resolved, discrete or continuous", {
    skip_if_not_installed("paletteer")
    skip_if_not_installed("ggsci")
    expect_identical(.resolve_palette("ggsci::default_jco", 3),
                     vapply(as.character(paletteer::paletteer_d("ggsci::default_jco"))[1:3],
                            .normalize_hex, character(1), USE.NAMES = FALSE))
    # A continuous palette gives as many colours as asked
    expect_length(.resolve_palette("grDevices::terrain.colors", 5, type = "continuous"), 5)
    expect_error(.resolve_palette("ggsci::no_such_palette", 3), "Unknown palette")
})

test_that("the plotting modules share the same behaviour", {
    # The cluster palette used to fail on a vector of colours
    expect_identical(configure_cluster_palette(3, c("#123456", "#654321", "#ABCDEF")),
                     c("1" = "#123456", "2" = "#654321", "3" = "#ABCDEF"))
    # and to fall back to its default on an unknown name, with a warning
    expect_error(configure_cluster_palette(3, "brewer:NoSuch"), "Unknown RColorBrewer")
    expect_error(get_heatmap_palette("NoSuchPalette"), "Unknown palette")
    expect_identical(unname(get_annotation_palette(c("A", "B"), "brewer:Set2")),
                     RColorBrewer::brewer.pal(3, "Set2")[1:2])
})
