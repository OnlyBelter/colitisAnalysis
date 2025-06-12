# In file: R/plotting.R

#' Create a Sorted Topic Structure Plot
#' @description Generates a structure plot with samples sorted by main topic proportion.
#' @param fit The `fastTopics` fit object.
#' @param samples_metadata The metadata data frame, aligned with the fit object.
#' @param grouping_cols Character vector of column names for grouping.
#' @param colitis_model_name The name of the model, for special conditional formatting.
#' @param k The number of topics.
#' @return A `ggplot` object.
#' @importFrom fastTopics structure_plot
#' @importFrom ggplot2 labs
#' @export
create_sorted_structure_plot <- function(fit, samples_metadata, grouping_cols, colitis_model_name, k) {
  if (length(grouping_cols) > 1) {
    grouping_values <- interaction(samples_metadata[, grouping_cols], sep = "_")
  } else {
    grouping_values <- samples_metadata[[grouping_cols[1]]]
  }
  grouping_values <- factor(grouping_values)

  if (colitis_model_name == 'Merged_dataset') {
    # Special re-leveling logic
    organs <- c("BM", "BR", "CE", "CO", "HE", "iLN", "KI", "LI", "LU", "mesLN", "SI", "SK", "SP", "TH")
    abx_cdif_levels <- as.vector(rbind(paste0(organs, "_AbxTreated"), paste0(organs, "_CDifficile")))
    other_levels <- levels(grouping_values)[!(levels(grouping_values) %in% abx_cdif_levels)]
    grouping_values <- factor(grouping_values, levels = c(abx_cdif_levels, other_levels))
  }

  L_matrix <- fit$L
  sample_names <- rownames(L_matrix)
  order_by_vec <- numeric(nrow(L_matrix))
  names(order_by_vec) <- sample_names

  for (group_name in levels(grouping_values)) {
    samples_in_group <- sample_names[grouping_values == group_name & !is.na(grouping_values)]
    if (length(samples_in_group) == 0) next
    L_subset <- L_matrix[samples_in_group, , drop = FALSE]
    mean_props <- if (nrow(L_subset) == 1) L_subset[1,] else colMeans(L_subset)
    order_by_vec[samples_in_group] <- -L_subset[, which.max(mean_props)]
  }

  p <- fastTopics::structure_plot(fit, colors = select_n_ordered_colors(k), grouping = grouping_values, order_by = order_by_vec, gap = 2) +
    ggplot2::labs(x = paste(colitis_model_name, "- k =", k, "- Sorted Structure Plot"))

  return(p)
}
