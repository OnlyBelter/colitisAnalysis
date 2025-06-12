# In file: R/de_analysis.R

#' Run DE Analysis and Extract Top Genes
#' @description Runs DE analysis or loads cached results, then extracts and saves top genes.
#' @param fit The `fastTopics` fit object.
#' @param counts_matrix The counts matrix used for the fit.
#' @param k The number of topics.
#' @param de_dir Directory to save/load DE results.
#' @param model_prefix Prefix for the DE result filename.
#' @return A data frame of top differentially expressed genes across all topics.
#' @importFrom fastTopics de_analysis
#' @importFrom utils head write.csv
#' @export
run_de_analysis <- function(fit, counts_matrix, k, de_dir, model_prefix) {
  if (!dir.exists(de_dir)) dir.create(de_dir, recursive = TRUE)
  de_file <- file.path(de_dir, paste0("de_result_", model_prefix, "_k", k, ".rds"))

  de <- NULL
  if (file.exists(de_file)) {
    message("Loading cached DE analysis result...")
    de <- readRDS(de_file)
  } else {
    message("Running DE analysis...")
    de <- tryCatch(
      fastTopics::de_analysis(fit, counts_matrix, pseudocount = 0.1, control = list(ns = 1e4, nc = 5)),
      error = function(e) {
        warning(paste("DE analysis failed for k =", k, ":", e$message))
        return(NULL)
      }
    )
    if (!is.null(de)) saveRDS(de, de_file)
  }

  if(is.null(de)) return(data.frame()) # Return empty frame if DE fails

  all_top_genes <- data.frame()
  for (topic_num in 1:k) {
    dat <- data.frame(postmean = de$postmean[, topic_num], z = de$z[, topic_num], lfsr = de$lfsr[, topic_num])
    dat <- dat[which(dat$lfsr < 0.01), ]
    if (nrow(dat) == 0) next
    dat <- dat[order(dat$postmean, decreasing = TRUE), ]
    top10 <- utils::head(dat, n = 10)
    if (nrow(top10) > 0) {
      top10$gene <- rownames(top10)
      top10$topic <- topic_num
      all_top_genes <- rbind(all_top_genes, top10)
    }
  }

  if (nrow(all_top_genes) > 0) {
    all_top_genes <- all_top_genes[, c("topic", "gene", "postmean", "z", "lfsr")]
    de_csv_file <- file.path(de_dir, paste0("top_de_genes_", model_prefix, "_k", k, ".csv"))
    utils::write.csv(all_top_genes, file = de_csv_file, row.names = FALSE)
    message(paste("Saved top DE genes to:", de_csv_file))
  }

  return(all_top_genes)
}
