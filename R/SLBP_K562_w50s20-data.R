#' ENCODE eCLIP data for SLBP in K562, low count filtered
#'
#' This is ENCODE eCLIP data which was quantified by htseq-clip
#' in sliding-windows of max. length 50nt, the step size was 20.
#' This is not ideal data for DEWSeq since it is lacking replicates,
#' however was small enough for the inclusion of the package.
#'
#' @docType data
#'
#' @usage data(SLBP_K562_w50s20)
#'
#' @format An object of class \code{"DESeq"};
#'
#' @examples
#'
#' data(SLBP_K562_w50s20)
#' SLBP_K562_w50s20
"SLBP_K562_w50s20"
