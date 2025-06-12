# In file: R/utils.R

# A curated palette of 30 distinguishable colors
DISTINGUISHABLE_COLORS_30 <- c(
  "#CBAFD4", "#835CA6", "#C2C2C2", "#2077B5", "#90B5E0", "#D62928", "#D87AB1",
  "#C49C94", "#F57E20", "#FCBA78", "#9DD089", "#F7B6D1", "#2DA048", "#62C6C0",
  "#8C574C", "#F69696", "#E5E878", "#984EA3", "#FFFF33", "#A65628",
  "#F781BF", "#999999", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#666666"
)

#' Select the First N Ordered Colors
#' This is an internal helper function, not exported to the user.
#' @param n The number of colors to select.
#' @param color_palette The palette to choose from.
#' @return A character vector of 'n' hex color codes.
select_n_ordered_colors <- function(n, color_palette = DISTINGUISHABLE_COLORS_30) {
  total_colors <- length(color_palette)
  if (!is.numeric(n) || n < 1 || n > total_colors) {
    stop(paste0("'n' must be an integer between 1 and ", total_colors, "."))
  }
  return(color_palette[1:n])
}
