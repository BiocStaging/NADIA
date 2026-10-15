# =============================================================================
# Highcharts in the vignettes
# =============================================================================
#
# Sourced from the hidden setup chunk of the vignettes that render Highcharts
# widgets (NADIA, benchmarking, pattern-profiler and visualization), so that
# they share one definition. It is not part of the package.
#
# highcharter 0.9.5 declares every Highcharts module for every widget. Some of
# those modules contain literal HTML closing tags, which can be exposed as text
# when Pandoc embeds the dependencies in a self-contained BiocStyle vignette,
# and together they add about 0.7 MB to each vignette. The examples need only
# the core and highcharts-more modules. Exporting remains available in ordinary
# NADIA use and is documented in vignette("visualization").
# =============================================================================

have_highcharter <- requireNamespace("highcharter", quietly = TRUE)

prepare_vignette_highchart <- function(widget) {
  widget <- highcharter::hc_exporting(widget, enabled = FALSE)

  highcharts_dir <- system.file(
    "htmlwidgets/lib/highcharts",
    package = "highcharter"
  )

  slim_dependency <- htmltools::htmlDependency(
    name = "highcharts",
    version = "9.3.1.9000",
    src = c(file = highcharts_dir),
    script = c("highcharts.js", "highcharts-more.js", "custom/reset.js"),
    all_files = FALSE
  )

  widget$dependencies <- c(widget$dependencies, list(slim_dependency))
  widget
}
