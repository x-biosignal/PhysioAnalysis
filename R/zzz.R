#' Package on-load hook
#'
#' @param libname Library path.
#' @param pkg Package name.
#' @keywords internal
.onLoad <- function(libname, pkg) {
  invisible(NULL)
}

# Silence R CMD check notes for ggplot2 NSE symbols used in aes().
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "time", "amplitude", "epoch", "lower", "upper",
    "freq", "power", "frequency", "value", "label",
    "channel", "y",
    "col_name", "row_name", "connectivity",
    "xend", "yend", "weight", "size", "name",
    "node", "metric"
  ))
}
