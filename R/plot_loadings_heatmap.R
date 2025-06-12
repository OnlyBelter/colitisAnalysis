# In file: R/plot_loadings_heatmap.R

#' Plot Topic Loadings as a Faceted Heatmap
#'
#' @description
#' This function takes a fastTopics fit object and creates a complex, faceted plot
#' showing the topic proportions (loadings) for each sample. It includes extensive
#' data wrangling to structure the plot by sample metadata.
#'
#' @details
#' The function assumes that the sample names (i.e., `rownames(fit$L)`) follow a
#' specific format: `"Tissue_Day_Mouse_Treatment"`, for example, `"BM_d04_1_DSS"`.
#' It uses this structure to facet the plot.
#'
#' @param fit_object A `fastTopics` fit object. The number of topics 'k' is
#'   inferred from this object.
#' @param color_palette An optional named vector of colors. The names should
#'   correspond to the topics (e.g., "k1", "k2"). If NULL, a default palette
#'   is used.
#' @param filter_by_treatment A character string. If not NULL, the plot will be
#'   filtered to only include Tissues that are present in this specified
#'   Treatment group. For example, setting this to "DSS" will only show tissues
#'   that were part of the DSS experiment. Defaults to NULL (shows all tissues).
#'
#' @return A `ggplot` object.
#'
#' @importFrom magrittr `%>%` set_names
#' @importFrom dplyr as_tibble mutate group_by ungroup arrange filter pull dense_rank
#' @importFrom tidyr pivot_longer separate
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 ggplot aes geom_col facet_grid theme_minimal theme element_text unit element_blank element_rect labs scale_fill_manual
#' @importFrom rlang .data
#' @export
plot_loadings_heatmap <- function(fit_object, color_palette = NULL, filter_by_treatment = NULL) {

  # Get k from the fit object
  k <- fit_object$k

  # --- 1. Define Color Palette ---
  if (is.null(color_palette)) {
    message(paste("Color palette not provided. Using default palette for k =", k))
    coloring <- select_n_ordered_colors(k) %>%
      magrittr::set_names(paste0("k", seq(1, k, 1)))
  } else {
    coloring <- color_palette
  }

  # --- 2. Data Wrangling ---
  # This section transforms the topic loadings matrix into a long-format
  # data frame suitable for ggplot.
  df <- dplyr::as_tibble(fit_object$L, rownames = "Sample") %>%
    tidyr::pivot_longer(
      cols = -.data$Sample,
      names_to = "K",
      values_to = "Loading"
    ) %>%
    tidyr::separate(.data$Sample, into = c("Tissue", "Day", "Mouse", "Treatment"), sep = "_", remove = FALSE) %>%
    dplyr::mutate(
      Tissue = factor(.data$Tissue, levels = sort(unique(.data$Tissue))),
      Day = factor(.data$Day, levels = sort(unique(.data$Day))),
      Mouse = factor(.data$Mouse, levels = sort(unique(.data$Mouse))),
      Treatment = factor(.data$Treatment, levels = c(
        "AbxTreated", "DSS", "Citro", "CDifficile"
      )),
      Time_numerical = as.numeric(sub("^d", "", .data$Day))
    ) %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::mutate(K = forcats::fct_reorder(.data$K, .data$Loading, .desc = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Treatment) %>%
    dplyr::arrange(.data$Time_numerical) %>%
    dplyr::mutate(dense_time = dplyr::dense_rank(.data$Time_numerical)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$Treatment, .data$Day) %>%
    dplyr::arrange(as.numeric(as.character(.data$Mouse))) %>%
    dplyr::mutate(replicate = dplyr::dense_rank(as.numeric(as.character(.data$Mouse)))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Time = paste0(.data$Day, "_", .data$replicate)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$Time_numerical, .data$dense_time) %>%
    dplyr::mutate(Time = factor(.data$Time, levels = sort(unique(.data$Time))))

  # --- 3. Optional Filtering ---
  df_to_plot <- df
  if (!is.null(filter_by_treatment)) {
    message(paste("Filtering plot to only include tissues present in the '", filter_by_treatment, "' treatment group."))
    tissues_to_keep <- df %>%
      dplyr::filter(.data$Treatment == filter_by_treatment) %>%
      dplyr::pull(.data$Tissue) %>%
      unique()

    df_to_plot <- df %>% dplyr::filter(.data$Tissue %in% tissues_to_keep)
  }

  # --- 4. Create the Plot ---
  p <- ggplot2::ggplot(df_to_plot, ggplot2::aes(x = .data$Time, y = .data$Loading, fill = .data$K)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::facet_grid(.data$Tissue ~ .data$Treatment, scales = "free_x", space = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      strip.text = ggplot2::element_text(size = 7),
      panel.spacing = ggplot2::unit(0.5, "lines"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    ) +
    ggplot2::labs(
      title = "Topic Proportions (Loadings) per Sample",
      x = "Timepoint and Replicate",
      y = "Topic Proportion"
    ) +
    ggplot2::scale_fill_manual(values = coloring,
                               breaks = names(coloring))

  return(p)
}
