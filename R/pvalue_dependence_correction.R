#' @title p-value correction
#' @description pvalue correction using local simes like approach
#' @param p_to_correct pvalue to correct (from the current window)
#' @param dependent_p_values pvalues from neighboring windows
#' @keywords internal
pvalue_dependence_correction <- function(p_to_correct, dependent_p_values){
  # author: Tom
  # p.to.correct is a numeric
  # dependent.p.values is a vector of p.values
  if(p_to_correct > 1| p_to_correct < 0 | any(dependent_p_values > 1) | any(dependent_p_values < 0)) {
    stop('P-values must fall in the range 0 <= p <= 1')
  }
  nr_of_p_values <- length(dependent_p_values) + length(p_to_correct)
  inv_ranks <- ((nr_of_p_values + 1) - rank(c(p_to_correct, dependent_p_values)))
  inv_rank_of_p_to_correct <- inv_ranks[1:length(p_to_correct)]
  return (min(1,(p_to_correct * nr_of_p_values / inv_rank_of_p_to_correct))) # to correct for p-values going above 1
}
