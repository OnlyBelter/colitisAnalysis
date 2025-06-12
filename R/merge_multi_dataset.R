# In file: R/merge_multi_dataset.R

#' Merge Multiple Counts Datasets and Save to CSV
#'
#' @description
#' This function reads multiple count data files (e.g., from different experiments),
#' finds the intersection of genes (columns) across all files, merges the data
#' for these common genes, and saves the final data frame to a specified CSV file.
#'
#' @details
#' The function assumes input files have samples as rows and genes as columns.
#' Row names from the original files are preserved as row names in the merged data
#' frame, provided they are unique across all input files.
#'
#' @param input_file_paths A character vector of full paths to the count data files.
#' @param output_csv_path The full path for the output CSV file that will be created.
#' @param transpose_before_saving Logical. If TRUE (the default), the merged data frame
#'   (samples x genes) will be transposed to (genes x samples) before being saved to the CSV.
#'
#' @return Invisibly returns the merged data frame (samples x genes). The primary
#'   purpose is the side effect of saving the data to a file.
#'
#' @importFrom utils read.csv write.csv head str
#'
#' @export
#' @examples
#' \dontrun{
#'   # Define input and output files
#'   in_files <- c("path/to/dataset1.txt", "path/to/dataset2.txt")
#'   out_file <- "path/to/merged_data.csv"
#'
#'   # Run the merge and save function
#'   merged_data <- merge_multi_dataset(
#'     input_file_paths = in_files,
#'     output_csv_path = out_file,
#'     transpose_before_saving = TRUE
#'   )
#' }
merge_multi_dataset <- function(input_file_paths, output_csv_path, transpose_before_saving = TRUE) {

  # --- 1. Load Datasets ---
  message("Loading datasets...")
  datasets_list <- lapply(input_file_paths, utils::read.csv, row.names = 1)

  # --- 2. Identify Common Genes ---
  message("Identifying common genes across datasets...")
  list_of_gene_names <- lapply(datasets_list, colnames)
  common_genes <- Reduce(intersect, list_of_gene_names)
  common_genes <- sort(common_genes)

  if (length(common_genes) == 0) {
    stop("No common genes found across all datasets. Cannot merge.")
  } else {
    message(paste("Found", length(common_genes), "common genes."))
  }

  # --- 3. Filter Datasets to Common Genes ---
  # Using an unnamed list to preserve original row names during rbind
  filtered_datasets_list <- vector("list", length(datasets_list))
  for (i in seq_along(datasets_list)) {
    current_df <- datasets_list[[i]]
    # drop = FALSE ensures it remains a data frame even if only one common gene
    filtered_datasets_list[[i]] <- current_df[, common_genes, drop = FALSE]
  }

  # --- 4. Merge Datasets ---
  message("Merging datasets row-wise...")
  merged_df <- do.call(rbind, filtered_datasets_list)

  # --- 5. Display Summary and Save Output ---
  message("\nStructure of the final merged dataset (before potential transposition):")
  utils::str(merged_df)
  message("\nFirst few rows of the merged dataset:")
  print(utils::head(merged_df))

  # Prepare data for saving
  data_to_save <- if (transpose_before_saving) {
    message("Transposing data frame for saving...")
    t(merged_df)
  } else {
    merged_df
  }

  # Save to CSV
  message(paste("Saving data to:", output_csv_path))
  utils::write.csv(data_to_save, file = output_csv_path, row.names = TRUE)

  # Return the merged (but not transposed) data frame invisibly
  # so it can be assigned to a variable if desired.
  return(invisible(merged_df))
}
