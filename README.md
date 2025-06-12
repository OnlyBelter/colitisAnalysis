Run the code

```R
# Install the package
# They will need devtools installed: install.packages("devtools")
devtools::install_github("OnlyBelter/colitisAnalysis")

# Load your custom-built package
library(colitisAnalysis)

# Define all your inputs
COUNTS_FILE <- "path/to/merged_data.csv"
SAMPLE_INFO_FILE <- "path/to/Sample_info.csv"
OUTPUT_DIR_BASE <- "path/to/analysis_output"
MODEL_NAME <- "Merged_dataset"
K_VALUES_TO_TEST <- c(3, 14)
GROUPING_COLUMNS <- c("Organs", "Sub_group")

# --- Run Workflow ---
# 1. Prepare data
prepared_data <- prepare_analysis_data(
  counts_file = COUNTS_FILE,
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

  # Create and save plot
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
