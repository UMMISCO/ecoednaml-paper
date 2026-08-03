# =============================================================================
# utils.R — shared utility functions for the EcoEdnaNet pipeline
# =============================================================================

get_sample_by_prevalence <- function(matrix, prevalence_rate) {
  message(paste("You set the prevalence rate to:", paste0(prevalence_rate, "%")))
  prevalence_rate <- prevalence_rate / 100
  if (prevalence_rate == 0)
    message("You are using the overall data without any filtering.")
  if (prevalence_rate == 1)
    message("You are using species that appear in all samples.")
  if (prevalence_rate < 0 || prevalence_rate > 1)
    stop("Prevalence rate must be between 0 and 100.")
  num_samples     <- nrow(matrix)
  min_occurrences <- ceiling(prevalence_rate * num_samples)
  valid_features  <- colSums(matrix > 0) >= min_occurrences
  sampled_matrix  <- matrix[, valid_features, drop = FALSE]
  sampled_matrix[rowSums(sampled_matrix > 0) > 0, , drop = FALSE]
}
