# In file: R/fit_topic_model.R

#' Fit or Load a Topic Model
#' @description Fits a new topic model or loads a pre-existing one from an .rds file.
#' @param counts_matrix The sample x gene counts matrix.
#' @param k The number of topics to fit.
#' @param model_dir Directory to save/load the model.
#' @param model_prefix Prefix for the model filename.
#' @return A `fastTopics` fit object, or NULL on failure.
#' @importFrom fastTopics fit_topic_model
#' @export
get_topic_model <- function(counts_matrix, k, model_dir, model_prefix) {
  if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)

  saved_fit_filename <- file.path(model_dir, paste0(model_prefix, "_k", k, ".rds"))
  fit <- NULL

  if (file.exists(saved_fit_filename)) {
    message("Loading existing model...")
    fit <- tryCatch(readRDS(saved_fit_filename), error = function(e) {
      warning(paste("Could not load model:", e$message, "\nWill re-fit."))
      return(NULL)
    })
  }

  if (is.null(fit)) {
    message("Fitting new model...")
    fit <- tryCatch(fastTopics::fit_topic_model(counts_matrix, k = k), error = function(e) {
      warning(paste("Error during model fitting:", e$message))
      return(NULL)
    })
    if (!is.null(fit)) saveRDS(fit, file = saved_fit_filename)
  }

  return(fit)
}
