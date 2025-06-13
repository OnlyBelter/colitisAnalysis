# In file: R/plot_loadings_heatmap.R

#' Plot Topic Loadings as a Faceted Heatmap
#'
#' @description
#' This function takes a fastTopics fit object and creates a complex, faceted plot
#' showing the topic proportions (loadings) for each sample.
#'
#' @details
#' The function assumes that the sample names (i.e., `rownames(fit$L)`) follow a
#' specific format: `"Tissue_Day_Mouse_Treatment"`, for example, `"BM_d04_1_DSS"`.
#' It uses this structure to facet the plot.
#'
#' @param fit_object A `fastTopics` fit object.
#' @param treatment_levels A character vector specifying which "Treatments" to
#'   include in the plot and in what order they should appear as columns. If NULL,
#'   all treatments are included in alphabetical order.
#' @param tissue_levels A character vector specifying which "Tissues" to include
#'   in the plot and in what order they should appear as rows. If NULL, all
#'   tissues are included in alphabetical order.
#' @param color_palette An optional named vector of colors. The names should
#'   correspond to the topics (e.g., "k1", "k2"). If NULL, a default palette
#'   is used.
#' @param save_to_file An optional character string specifying the full path to
#'   save the plot file (e.g., "output/plot.png"). If NULL (the default), the
#'   ggplot object is returned instead of being saved.
#' @param fig_width The width of the saved plot in inches. Defaults to 16.
#' @param fig_height The height of the saved plot in inches. Defaults to 9.
#'
#' @return If `save_to_file` is NULL, returns a `ggplot` object. If a path is provided,
#'   the function saves the plot and invisibly returns the file path.
#'
#' @importFrom magrittr `%>%` set_names
#' @importFrom dplyr as_tibble mutate group_by ungroup arrange desc
#' @importFrom tidyr pivot_longer separate
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 ggplot aes geom_col facet_grid theme_minimal theme element_text unit element_blank element_rect labs scale_fill_manual guides guide_legend ggsave
#' @importFrom rlang .data
#' @export
plot_loadings_heatmap <- function(fit_object,
                                  treatment_levels = NULL,
                                  tissue_levels = NULL,
                                  color_palette = NULL,
                                  save_to_file = NULL,
                                  fig_width = 16,
                                  fig_height = 9) {

  # --- 1. Get k and Define Color Palette ---
  k <- ncol(fit_object$L)

  if (is.null(color_palette)) {
    message(paste("Color palette not provided. Using default palette for k =", k))
    coloring <- select_n_ordered_colors(k) %>%
      magrittr::set_names(paste0("k", seq(1, k, 1)))
  } else {
    coloring <- color_palette
  }

  # --- 2. Data Wrangling (Corrected Logic) ---
  df <- dplyr::as_tibble(fit_object$L, rownames = "Sample") %>%
    tidyr::pivot_longer(
      cols = -.data$Sample,
      names_to = "K",
      values_to = "Loading"
    ) %>%
    tidyr::separate(.data$Sample, into = c("Tissue", "Day", "Mouse", "Treatment"), sep = "_", remove = FALSE) %>%
    # **CHANGE**: Create a stable factor for K. This helps the legend.
    dplyr::mutate(
      K = factor(.data$K, levels = paste0("k", 1:k)),
      Time_numerical = as.numeric(sub("^d", "", .data$Day))
    ) %>%
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

  # --- 3. Filter data and set factor levels based on user input ---
  df_to_plot <- df

  if (!is.null(treatment_levels)) {
    df_to_plot <- df_to_plot %>% dplyr::filter(.data$Treatment %in% treatment_levels)
    df_to_plot$Treatment <- factor(df_to_plot$Treatment, levels = treatment_levels)
  } else {
    df_to_plot$Treatment <- factor(df_to_plot$Treatment)
  }

  if (!is.null(tissue_levels)) {
    df_to_plot <- df_to_plot %>% dplyr::filter(.data$Tissue %in% tissue_levels)
    df_to_plot$Tissue <- factor(df_to_plot$Tissue, levels = tissue_levels)
  } else {
    df_to_plot$Tissue <- factor(df_to_plot$Tissue)
  }

  # **CHANGE**: To control stacking order (largest bar at the bottom),
  # we arrange the data frame itself before plotting.
  df_to_plot <- df_to_plot %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::arrange(dplyr::desc(.data$Loading), .by_group = TRUE) %>%
    dplyr::ungroup()


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
    ggplot2::scale_fill_manual(
      name = "Topic",
      values = coloring,
      breaks = names(coloring)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, title.position = "top"))

  # --- 5. Save or Return the Plot ---
  if (!is.null(save_to_file)) {
    message(paste("Saving plot to:", save_to_file))
    ggplot2::ggsave(
      filename = save_to_file,
      plot = p,
      width = fig_width,
      height = fig_height,
      dpi = 300,
      limitsize = FALSE
    )
    return(invisible(save_to_file))
  } else {
    return(p)
  }
}
