#' @export
#' @title annotate DESeq2 results
#' @description annotate DESeq2 results using annotation file from htseq-clip
#' @param results DESeq2 results data.frame, from \code{\link{DESeq}}, \code{\link{results}}
#' @param annotationFile sliding window annotation file, can be plain either text or .gz file,
#' the file MUST be TAB separated, and MUST have the following columns:\cr
#' \code{chromosome} chromosome name \cr
#' \code{unique_id} unique id of the window \cr
#' \code{begin} window start co-ordinate \cr
#' \code{end} window end co-ordinate \cr
#' \code{strand} strand \cr
#' \code{gene_id} gene id \cr
#' \code{gene_name} gene name \cr
#' \code{gene_type} gene type annotation \cr
#' \code{gene_region} gene region \cr
#' \code{Nr_of_region} number of the current region \cr
#' \code{Total_nr_of_region} total number of regions \cr
#' \code{window_number} window number \cr
#' @return data.frame
annotateResults <- function(results,annotationFile){
    results <- na.omit(as.data.frame(results))
    if ('row' %in% colnames(results)){
        rownames(results) <- results$row
    }
    neededCols <- c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
    missingCols <- setdiff(neededCols,colnames(results))
    if(length(missingCols)>0){
        stop('Input results table doesnot have all DESeq2 columns. Missing columns:',paste(missingCols,collapse=", "),'')
    }
    resGrange <- .readAnnotation(fname=annotationFile,uniqIds = rownames(results),asGRange=FALSE)
   return(cbind(results,as.data.frame(resGrange[rownames(results),])))
}
