Run the code

```R
# Install the package
# They will need devtools installed: install.packages("devtools")
devtools::install_github("OnlyBelter/colitisAnalysis")

# Load your custom-built package
library(colitisAnalysis)

setwd('/Users/belter/github/UChicago/from_Madison/')

# Define all your inputs
IN_FILES <- c(
  "./raw_data/dataset1.txt",
  "./raw_data/dataset2.txt",
  "./raw_data/dataset3.txt"
)
MERGED_COUNTS_FILE <- "./raw_data/merged_colitis_data.csv"
SAMPLE_INFO_FILE <- "./Sample_info.csv"
OUTPUT_DIR_BASE <- "./analysis_output"
MODEL_NAME <- "Merged_dataset"
K_VALUES_TO_TEST <- c(3, 14)
GROUPING_COLUMNS <- c("Organs", "Sub_group")  # Columns in Sample_info.csv

# Merge three data sets if merged file not exists
if (!file.exists(MERGED_COUNTS_FILE)) {
  merge_multi_dataset(input_file_paths = IN_FILES,
                      output_csv_path = MERGED_COUNTS_FILE,
                      transpose_before_saving = TRUE)
}


# --- Run Workflow ---
# 1. Prepare data
prepared_data <- prepare_analysis_data(
  counts_file = MERGED_COUNTS_FILE,
  sample_info_file = SAMPLE_INFO_FILE,
  colitis_model_name = MODEL_NAME
)

# 2. Loop through k values
for (k in K_VALUES_TO_TEST) {
  message(paste0("\n--- Running workflow for k = ", k, " ---"))

  # Get model
  fit_object <- get_topic_model(
    counts_matrix = prepared_data$counts_matrix,
    k = k,
    model_dir = file.path(OUTPUT_DIR_BASE, "models"),
    model_prefix = MODEL_NAME
  )

  if (is.null(fit_object)) next

  # Create and save structure plot
  plot_object <- create_sorted_structure_plot(
    fit = fit_object,
    samples_metadata = prepared_data$samples_metadata,
    grouping_cols = GROUPING_COLUMNS,
    colitis_model_name = MODEL_NAME,
    k = k
  )
  plot_dir <- file.path(OUTPUT_DIR_BASE, "plots")
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  ggplot2::ggsave(
    file.path(plot_dir, paste0("structure_plot_", MODEL_NAME, "_k", k, ".png")),
    plot_object, width = 16.7, height = 6.7, dpi = 300, limitsize = FALSE
  )

  # Create and save structure plot by another way (Provided by Maggie)
  plot_object2 <- plot_loadings_heatmap(
    fit_object,
    treatment_levels = c("DSS", "Citro", "AbxTreated", "CDifficile"),
    tissue_levels = c("BM", "BR", "CE", "CO"),
    fig_width = 12,
    fig_height = 8,
    save_to_file = file.path(plot_dir, paste0("new_structure_plot_", MODEL_NAME, "_k", k, ".png"))
  )

  # Run DE analysis
  run_de_analysis(
    fit = fit_object,
    counts_matrix = prepared_data$counts_matrix,
    k = k,
    de_dir = file.path(OUTPUT_DIR_BASE, "de_analysis"),
    model_prefix = MODEL_NAME
  )
}
```
