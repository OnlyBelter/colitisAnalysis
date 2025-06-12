# In file: R/prepare_data.R

#' Prepare Analysis Data
#'
#' @description Loads counts and metadata, transposes counts if needed, finds common
#'   samples, and filters both datasets to ensure alignment.
#'
#' @param counts_file Path to the counts data file.
#' @param sample_info_file Path to the sample metadata file.
#' @param colitis_model_name The name of the model being run (e.g., "Merged_dataset").
#' @param filter_out_controls If TRUE, removes samples where Sub_group is "Control".
#'
#' @return A list containing the filtered `counts_matrix` and `samples_metadata` data frame.
#'
#' @importFrom utils read.csv
#' @export
prepare_analysis_data <- function(counts_file, sample_info_file, colitis_model_name, filter_out_controls = TRUE) {
  message("Starting data preparation...")

  counts_orig <- utils::read.csv(counts_file, row.names = 1)
  if (nrow(counts_orig) > ncol(counts_orig)) {
    message("Transposing counts matrix to be sample x gene...")
    counts_orig <- t(counts_orig)
  }

  sample_info_orig <- utils::read.csv(sample_info_file, row.names = 1)

  current_samples_meta <- sample_info_orig
  if (filter_out_controls) {
    current_samples_meta <- current_samples_meta[current_samples_meta$Sub_group != 'Control', ]
  }

  if (colitis_model_name != 'Merged_dataset') {
    current_samples_meta <- current_samples_meta[!is.na(current_samples_meta$Colitis_model) &
                                                   current_samples_meta$Colitis_model == colitis_model_name, ]
  }

  common_sample_ids <- intersect(rownames(counts_orig), rownames(current_samples_meta))
  if (length(common_sample_ids) == 0) stop("No common samples found between counts and metadata.")

  message(paste("Found", length(common_sample_ids), "common samples for analysis."))

  counts_filtered <- as.matrix(counts_orig[common_sample_ids, ])
  samples_metadata_filtered <- current_samples_meta[common_sample_ids, ]

  message(paste("Final filtered counts matrix dimensions:", paste(dim(counts_filtered), collapse = " x ")))

  return(list(counts_matrix = counts_filtered, samples_metadata = samples_metadata_filtered))
}
