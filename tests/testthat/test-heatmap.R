# =============================================================================
# proteomics_heatmap(): which proteins it shows, how it scales them, and what
# it writes
#
# The exported TSV is the matrix the heatmap draws, after selection and
# scaling, so the checks are made on it rather than on the graphical object.
# =============================================================================

.hm_res <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            data(nadia_dia, package = "NADIA", envir = environment())
            cache <<- process_proteomics(nadia_dia, verbose = FALSE)
        }
        cache
    }
})

.hm_skip <- function() {
    skip_if_not_installed("ComplexHeatmap")
    skip_if_not_installed("tidyHeatmap")
    skip_if_not_installed("circlize")
}

# Build a heatmap, export its matrix and read it back
.hm_matrix <- function(data, ...) {
    f <- withr::local_tempfile(fileext = ".tsv", .local_envir = parent.frame())
    hm <- suppressMessages(proteomics_heatmap(data, export_path = f, ...))
    tab <- utils::read.delim(f, check.names = FALSE)
    list(hm = hm, tab = tab,
         mat = as.matrix(tab[, setdiff(names(tab), c("FeatureID", "adjP"))]))
}

test_that("each mode shows the proteins it is meant to", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    feat <- pca[!duplicated(pca$FeatureID), ]

    all <- .hm_matrix(pca, mode = "all", scale_data = "none")
    expect_setequal(all$tab$FeatureID, unique(pca$FeatureID))
    expect_identical(ncol(all$mat), length(unique(pca$SampleID)))

    any <- .hm_matrix(pca, mode = "any")
    expect_setequal(any$tab$FeatureID, feat$FeatureID[which(feat$sig_any)])

    tgt <- .hm_matrix(pca, mode = "target", comparison = "B-A", alpha = 0.01)
    expected <- feat$FeatureID[which(feat[["adjP_B-A"]] <= 0.01)]
    expect_setequal(tgt$tab$FeatureID, expected)
    expect_lt(length(expected), nrow(any$tab))
    # The targeted heatmap carries the adjusted p-value of its comparison
    expect_true(all(tgt$tab$adjP <= 0.01))
})

test_that("an explicit set of proteins overrides the mode", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    ids <- unique(pca$FeatureID)[c(3, 10, 25, 40)]
    m <- .hm_matrix(pca, mode = "any", feature_ids = ids)
    expect_setequal(m$tab$FeatureID, ids)
    expect_warning(
        suppressMessages(proteomics_heatmap(pca, feature_ids = c(ids, "not_a_protein"))),
        "not_a_protein")
})

test_that("the three scalings do what they say", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    ids <- unique(pca$FeatureID)[1:60]

    raw <- .hm_matrix(pca, feature_ids = ids, scale_data = "none")
    # Unscaled values are the intensities of the input table
    one <- pca[pca$FeatureID == raw$tab$FeatureID[1], ]
    expect_equal(unname(raw$mat[1, one$SampleID]), one$Intensity)

    row <- .hm_matrix(pca, feature_ids = ids, scale_data = "row")
    expect_equal(unname(rowMeans(row$mat, na.rm = TRUE)), rep(0, nrow(row$mat)),
                 tolerance = 1e-12)
    expect_equal(unname(apply(row$mat, 1, sd, na.rm = TRUE)), rep(1, nrow(row$mat)),
                 tolerance = 1e-12)

    col <- .hm_matrix(pca, feature_ids = ids, scale_data = "column")
    expect_equal(unname(colMeans(col$mat, na.rm = TRUE)), rep(0, ncol(col$mat)),
                 tolerance = 1e-12)
})

test_that("a constant protein is dropped before the row z-score, with a message", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    ids <- unique(pca$FeatureID)[1:20]
    sub <- pca[pca$FeatureID %in% ids, ]
    sub$Intensity[sub$FeatureID == ids[1]] <- 20   # no variance: z-score undefined
    expect_message(hm <- proteomics_heatmap(sub, scale_data = "row"),
                   "1 features with 0/NA variance dropped")
    m <- .hm_matrix(sub, scale_data = "row")
    expect_false(ids[1] %in% m$tab$FeatureID)
    expect_identical(nrow(m$tab), length(ids) - 1L)
})

test_that("the heatmap draws, with annotations, palettes and a title", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    ids <- unique(pca$FeatureID)[1:30]
    annot <- data.frame(FeatureID = ids,
                        Group = rep(c("g1", "g2", "g3"), length.out = 30))
    annot_file <- withr::local_tempfile(fileext = ".tsv")
    utils::write.table(annot, annot_file, sep = "\t", quote = FALSE,
                       row.names = FALSE)
    png_file <- withr::local_tempfile(fileext = ".png")

    for (ra in list(annot, annot_file)) {
        hm <- suppressMessages(proteomics_heatmap(
            pca, mode = "target", comparison = "B-A", feature_ids = ids,
            row_annotation = ra, split_rows_by = "Group",
            sample_order = "condition", condition_order = c("A", "B", "D"),
            palette_value = "brewer:RdBu", palette_annotation = c("#1b9e77", "#d95f02", "#7570b3"),
            heatmap_title = "B-A"))
        expect_s3_class(hm, "proteomics_heatmap")
        grDevices::png(png_file)
        expect_no_error(print(hm))
        grDevices::dev.off()
    }
    expect_gt(file.size(png_file), 0)
})

test_that("the heatmap is written to PNG, PDF and SVG", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    ids <- unique(pca$FeatureID)[1:30]
    exts <- c("png", "pdf", if (capabilities("cairo")) "svg")
    for (ext in exts) {
        f <- withr::local_tempfile(fileext = paste0(".", ext))
        suppressMessages(proteomics_heatmap(pca, feature_ids = ids, export_file = f,
                                            plot_width = 5, plot_height = 4,
                                            export_dpi = 72))
        expect_gt(file.size(f), 0, label = ext)
    }
})

test_that("proteomics_heatmap_list builds one heatmap per mode", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    hms <- suppressMessages(proteomics_heatmap_list(pca, modes = c("any", "B-A", "D-A")))
    expect_identical(names(hms), c("any", "B-A", "D-A"))
    expect_true(all(vapply(hms, function(h) !is.null(h), logical(1))))
})

test_that("the palettes resolve every format the documentation accepts", {
    expect_length(get_heatmap_palette(NULL, n = 9), 9)
    expect_identical(get_heatmap_palette("brewer:RdBu", n = 9),
                     RColorBrewer::brewer.pal(9, "RdBu"))
    expect_identical(get_heatmap_palette("brewer:RdBu", n = 9, reverse = TRUE),
                     rev(RColorBrewer::brewer.pal(9, "RdBu")))
    expect_identical(get_heatmap_palette(c("#000000", "#FFFFFF")),
                     c("#000000", "#FFFFFF"))

    ann <- get_annotation_palette(c("A", "B", "D"), "brewer:Dark2")
    expect_identical(names(ann), c("A", "B", "D"))
    expect_identical(unname(ann), RColorBrewer::brewer.pal(3, "Dark2"))
})

test_that("the input is validated", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    expect_error(proteomics_heatmap(pca[, setdiff(names(pca), "Intensity")]),
                 "Intensity")
    expect_error(proteomics_heatmap(pca, mode = "target"), "comparison",
                 ignore.case = TRUE)
})

test_that("a row annotation does not have to mention every protein", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    ids <- unique(pca$FeatureID)
    shown <- unique(pca$FeatureID[pca$sig_any])
    # Half the proteins shown are annotated; one value is empty
    part <- data.frame(FeatureID = shown[seq(1, length(shown), 2)],
                       Type = c("", rep_len(c("x", "y"), ceiling(length(shown) / 2) - 1)))
    hm <- suppressMessages(proteomics_heatmap(pca, mode = "any", row_annotation = part,
                                              row_order_by = "Type", split_rows_by = "Type"))
    d <- hm@data
    expect_false(anyNA(d$FeatureID))
    expect_setequal(unique(as.character(d$FeatureID)), shown)
    # Unannotated and empty both become "NA", drawn white and placed last
    expect_identical(levels(d$Type), c("x", "y", "NA"))
    unannotated <- setdiff(shown, part$FeatureID)
    expect_true(all(d$Type[d$FeatureID %in% unannotated] == "NA"))
    expect_true(all(tail(levels(d$FeatureID), length(unannotated) + 1) %in%
                    c(unannotated, part$FeatureID[1])))
    expect_no_error(suppressMessages(tidyHeatmap::as_ComplexHeatmap(hm)))
})

test_that("each condition keeps its colour in a targeted heatmap", {
    .hm_skip()
    pca <- .hm_res()$PCA_Input
    pal <- c(A = "#FF0000", B = "#0000FF", D = "#00FF00")
    hm <- suppressMessages(proteomics_heatmap(pca, mode = "target", comparison = "D-A",
                                              alpha = 0.01, palette_annotation = pal))
    expect_identical(levels(hm@data$Condition), c("A", "D"))
    # tidyHeatmap assigns colours by level position, so the D samples used to
    # get B's colour, or none
    cm <- suppressMessages(tidyHeatmap::as_ComplexHeatmap(hm))@top_annotation@
        anno_list$Condition@color_mapping@colors
    expect_identical(toupper(substr(cm[c("A", "D")], 1, 7)), pal[c("A", "D")])
})
