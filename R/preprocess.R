
#' @export
#' @importFrom methods as is
#' @importFrom BiocGenerics sort
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @title create DESeq data object
#' @description create DESeq data object from sliding window counts,
#' phenotype data and annotation data
#'
#' @param countData sliding window count data
#' @param colData phenotype data
#' @param annotObj can either be a data.frame or a file name, see details
#' @param design design of the experiment, see \code{\link[DESeq2:DESeqDataSet]{DESeqDataSet}}
#' @param tidy If TRUE, first column is of countData is treated as rownames (defalt: FALSE), see \code{\link[DESeq2:DESeqDataSet]{DESeqDataSet}}
#' @param ignoreRank ignore rank, see \code{\link[DESeq2:DESeqDataSet]{DESeqDataSet}}
#' @param start0based TRUE (default) or FALSE. If TRUE, then the start positions in \code{annotObj} is  considered to be 0-based
#'
#' @details
#' If \code{annotObj} is a file name, the input file MUST be <TAB> separated, and supports reading in .gz files. If \code{annotObj} is a data.frame,
#' \code{colnames(annotObj)} MUST not be empty.\cr
#' This function checks for the following columns after reading in the file or on data.frame:\cr
#'  chromosome: chromosome name \cr
#'  unique_id: unique id of the window \cr
#'  begin: window start co-ordinate \cr
#'  end: window end co-ordinate \cr
#'  strand: strand \cr
#'  gene_id: gene id \cr
#'  gene_name: gene name \cr
#'  gene_type: gene type annotation \cr
#'  gene_region: gene region \cr
#'  Nr_of_region: number of the current region \cr
#'  Total_nr_of_region: total number of regions \cr
#'  window_number: window number \cr
#'
#' @return DESeq object
DESeqDataSetFromSlidingWindows <- function(countData, colData, annotObj, design, tidy=FALSE,ignoreRank=FALSE, start0based=TRUE){

  stopifnot(ncol(countData) > 1)
  if(is(countData,'data.table')){
    countData <- data.frame(countData)
    warning('countData is a data.table object, converting it to data.frame. First column will be used as rownames')
    if(!tidy){
      tidy <- TRUE
    }
  }
  if (tidy) {
    # this code is copied from DESeq2 source, credits to authors
    rownms <- as.character(countData[,1])
    countData <- countData[,-1,drop=FALSE]
    rownames(countData) <- rownms
  }else if(is.null(rownames(countData))){
    stop('rownames of countData cannot be empty')
  }
  if(is(annotObj,"character")){
    annotData <- .readAnnotation(fname=annotObj,asGRange=FALSE,start0based=start0based)
  }else if(is(annotObj,"data.frame")){
    neededCols <- c('chromosome','unique_id','begin','end','strand','gene_id','gene_name','gene_type','gene_region','Nr_of_region',
                    'Total_nr_of_region','window_number')
    missingCols <- setdiff(neededCols,colnames(annotObj))
    if(length(missingCols)>0){
      stop('Input annotation file is missing required columns, needed columns:
         chromosome: chromosome name
         unique_id: unique id of the window
         begin: window start co-ordinate
         end: window end co-ordinate
         strand: strand
         gene_id: gene id
         gene_name: gene name
         gene_type: gene type annotation
         gene_region: gene region
         Nr_of_region: number of the current region
         Total_nr_of_region: total number of regions
         window_number: window number
         Missing columns:
         ',paste(missingCols,collapse=", "),'')
    }
    if(length(intersect(as.character(annotObj$unique_id),rownames(countData)))==0){
      stop('There are no common unique ids between the countData and annotObj. Please check your data sets!')
    }
    annotData <- annotObj
    rownames(annotData) <- annotData$unique_id
  }else{
    stop('annotObj MUST be a data.frame or character')
  }
  if(length(setdiff(rownames(countData),rownames(annotData)))>0){
    warning('Cannot find chromosomal positions for all entries in countData.
            countData rows with missing annotation will be removed !')
  }
  commonIds <- intersect(rownames(annotData),rownames(countData))
  countData <- countData[commonIds,]
  annotData <- annotData[commonIds,]
  # The code below is copied from DESeq2 source, credits to authors
  if (is(colData,"data.frame")){
    colData <- as(colData, "DataFrame")
  }
  # check if the rownames of colData are simply in different order
  # than the colnames of the countData, if so throw an error
  # as the user probably should investigate what's wrong
  if (!is.null(rownames(colData)) & !is.null(colnames(countData))) {
    if (all(sort(rownames(colData)) == sort(colnames(countData)))) {
      if (!all(rownames(colData) == colnames(countData))) {
        stop(paste("rownames of the colData:
  ",paste(rownames(colData),collapse=","),"
  are not in the same order as the colnames of the countData:
  ",paste(colnames(countData),collapse=",")))
      }
    }
  }
  if (is.null(rownames(colData)) & !is.null(colnames(countData))) {
    rownames(colData) <- colnames(countData)
  }
  # assemble summarized experiment
  gr <- makeGRangesFromDataFrame(df=annotData,seqnames.field='chromosome',start.field='begin',
                                       end.field='end',strand.field='strand',
                                       keep.extra.columns=TRUE,starts.in.df.are.0based=start0based)
  gr <- sortSeqlevels(gr)
  gr <- sort(gr)
  countData <- as.matrix(countData[rownames(mcols(gr)),])
  se <- SummarizedExperiment(assays = SimpleList(counts=countData),colData = colData,rowRanges=gr)
  return(DESeqDataSet(se, design = design, ignoreRank))
}
