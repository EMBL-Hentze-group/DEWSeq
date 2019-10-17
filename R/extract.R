# R function(s) to extract regions from significant windows
# Author: Sudeep Sahadevan, sudeep.sahadevan@embl.de

#' @export
#'
#' @importFrom BiocGenerics sort
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
#' @importFrom S4Vectors na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @title extract significant regions
#' @description extract significant windows from output of
#'     \code{\link{resultsDEWSeq}} using the supplied padj and
#'      log2FoldChange cut-offs and merge these significant windows
#'       to regions and create the following columns for each
#'        significant region:
#' \itemize{
#' \item \code{padj_min}: min. padj value in the region
#' \item \code{padj_mean}: average padj value in the region
#' \item \code{padj_max}: max. padj value in the region
#' \item \code{log2FoldChange_min}: min. log 2 fold change in the region
#' \item \code{log2FoldChange_mean}: average log 2 fold change in the region
#' \item \code{log2FoldChange_max}: max. log 2 fold change in the region
#' }
#'
#' @details
#' The output data.frame from this function will have the following columns:
#' \itemize{
#'  \item \code{chromosome}: chromosome name
#'  \item \code{regionStartId}: \code{unique_id} of the left most window,
#'  where an enriched region begins
#'  \item \code{region_begin}: starting position of the enriched region
#'  \item \code{region_end}: ending position of the enriched region
#'  \item \code{strand}: strand info
#'  \item \code{windows_in_region}: total number of windows
#'   that make up the enriched region
#'  \item \code{region_length}: length of the enrched region
#'  \item \code{gene_id}: gene id
#'  \item \code{gene_name}: gene name
#'  \item \code{gene_type}: gene type annotation
#'  \item \code{gene_region}: gene region
#'  \item \code{Nr_of_region}: number of the current region
#'  \item \code{Total_nr_of_region}: total number of regions
#'  \item \code{window_number}: window number
#'  \item \code{padj_min}: min. padj value in the region
#'  \item \code{padj_mean}: average padj value in the region
#'  \item \code{padj_max}: max. padj value in the region
#'  \item \code{log2FoldChange_min}: min. log 2 fold change in the region
#'  \item \code{log2FoldChange_max}: max. log 2 fold change in the region
#'  \item \code{log2FoldChange_mean}: average log 2 fold change in the region
#' }

#'
#' @param windowRes \code{data.frame}, output from \code{\link{resultsDEWSeq}}
#' @param padjCol \code{character}, name of the adjusted pvalue column
#' (default: padj)
#' @param padjThresh \code{numeric}, threshold for p-adjusted value
#' (default: 0.05)
#' @param log2FoldChangeCol \code{character}, name of the log2foldchange column
#' (default: log2FoldChange)
#' @param log2FoldChangeThresh \code{numeric}, threshold for log2foldchange value
#' (default:1)
#' @param start0based \code{logical}, TRUE (default) or FALSE.
#' If TRUE, then the start positions in \code{windowRes}
#' is  considered to be 0-based
#'
#' @examples
#'
#' data("slbpWindows")
#' # using default cut-off thresholds,
#' # 'pSlidingWindows.adj' padj value columns
#' slbpRegions <- extractRegions(slbpWindows,
#' padjCol = 'pSlidingWindows.adj')
#'
#' @return data.frame
extractRegions <- function(windowRes,
                           padjCol              = 'padj',
                           padjThresh           = 0.05,
                           log2FoldChangeCol    = 'log2FoldChange',
                           log2FoldChangeThresh = 1,
                           start0based          = TRUE){
  requiredCols <- c('chromosome',
                    'unique_id',
                    'begin',
                    'end',
                    'strand',
                    'gene_id',
                    'gene_name',
                    'gene_type',
                    'gene_region',
                    'Nr_of_region',
                    'Total_nr_of_region',
                    'window_number',
                    padjCol,
                    log2FoldChangeCol)

  missingCols <- setdiff(requiredCols,colnames(windowRes))
  if(length(missingCols)>0){
    stop('Input data.frame is missing required columns, needed columns:
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
     ',padjCol,': p-adjusted value column
     ',log2FoldChangeCol,': log2foldchange column.
      Missing columns: ',paste(missingCols,collapse=", "),'')
  }
  windowRes <- na.omit(windowRes[,requiredCols])
  sigDat <- windowRes[ windowRes[,padjCol]<=padjThresh &
                         windowRes[,log2FoldChangeCol]>=log2FoldChangeThresh, ]
  if(nrow(sigDat)==0){

    warning('There are no significant windows/regions under the current threshold!
    Please lower your significance cut-off thresholds and manually check if there are any significant windows under the threshold')
    return(NULL)
  }
  geneRange <- makeGRangesFromDataFrame(sigDat,seqnames.field = 'gene_id',
                                        start.field = 'begin',end.field = 'end',
                                        strand.field = 'strand',ignore.strand=FALSE,
                                        starts.in.df.are.0based=start0based,
                                        keep.extra.columns = TRUE)
  geneRange <-  sortSeqlevels(geneRange)
  geneRange <- sort(geneRange)
  geneReduce <- reduce(geneRange,drop.empty.ranges=TRUE,
                       with.revmap=TRUE,min.gapwidth=1)
  mcols(geneReduce)$regionStartId <- 'Undefined'
  mcols(geneReduce)$gene_name <- 'Undefined'
  mcols(geneReduce)$gene_type <- 'Undefined'
  mcols(geneReduce)$gene_region  <- 'Undefined'
  mcols(geneReduce)$chromosome <- 'Undefined'
  mcols(geneReduce)$unique_ids  <- 'Undefined'
  mcols(geneReduce)$windows_in_region <- 1
  mcols(geneReduce)$Nr_of_region <-  1
  mcols(geneReduce)$Total_nr_of_region <- 1
  mcols(geneReduce)$window_number <- 1
  mcols(geneReduce)$padj_min <- 1.0
  mcols(geneReduce)$padj_max <- 1.0
  mcols(geneReduce)$padj_mean <- 1.0
  mcols(geneReduce)$log2FoldChange_min <- 0.0
  mcols(geneReduce)$log2FoldChange_max <- 0.0
  mcols(geneReduce)$log2FoldChange_mean <- 0.0
  pb <- txtProgressBar(min=0,max=length(geneReduce),style = 3)
  for(i in seq_len(length(geneReduce))){
    # values
    mergeInd <- unlist(mcols(geneReduce)[i,1])
    mcols(geneReduce)[i,'windows_in_region'] <- length(mergeInd)
    mcols(geneReduce)[i,'padj_min'] <- min(mcols(geneRange)[mergeInd,padjCol])
    mcols(geneReduce)[i,'padj_max'] <- max(mcols(geneRange)[mergeInd,padjCol])
    mcols(geneReduce)[i,'padj_mean'] <- mean(mcols(geneRange)[mergeInd,padjCol])
    mcols(geneReduce)[i,'log2FoldChange_min'] <- min(mcols(geneRange)[mergeInd,log2FoldChangeCol])
    mcols(geneReduce)[i,'log2FoldChange_max'] <- max(mcols(geneRange)[mergeInd,log2FoldChangeCol])
    mcols(geneReduce)[i,'log2FoldChange_mean'] <- mean(mcols(geneRange)[mergeInd,log2FoldChangeCol])
    # annotations
    mcols(geneReduce)[i,'regionStartId'] <- mcols(geneRange)[min(mergeInd),'unique_id']
    mcols(geneReduce)[i,'unique_ids'] <- paste(mcols(geneRange)[mergeInd,'unique_id'],collapse = ', ')
    mcols(geneReduce)[i,'chromosome'] <- mcols(geneRange)[min(mergeInd),'chromosome']
    mcols(geneReduce)[i,'gene_name'] <- mcols(geneRange)[min(mergeInd),'gene_name']
    mcols(geneReduce)[i,'gene_type'] <- mcols(geneRange)[min(mergeInd),'gene_type']
    mcols(geneReduce)[i,'gene_region'] <- mcols(geneRange)[min(mergeInd),'gene_region']
    mcols(geneReduce)[i,'Nr_of_region'] <- mcols(geneRange)[min(mergeInd),'Nr_of_region']
    mcols(geneReduce)[i,'Total_nr_of_region'] <- mcols(geneRange)[min(mergeInd),'Total_nr_of_region']
    mcols(geneReduce)[i,'window_number'] <- mcols(geneRange)[min(mergeInd),'window_number']
    # progress
    setTxtProgressBar(pb=pb,value=i)
  }
  close(pb)
  regionRes <- as.data.frame(geneReduce)
  regionRes <- regionRes[,c('chromosome','start','end','strand','windows_in_region','width',
                            'padj_min','padj_mean','padj_max',
                            'log2FoldChange_min','log2FoldChange_mean','log2FoldChange_max',
                            'regionStartId','seqnames','gene_name','gene_type','gene_region',
                            'Nr_of_region','Total_nr_of_region','window_number','unique_ids')]
  colnames(regionRes)[2] <- 'region_begin'
  colnames(regionRes)[3] <- 'region_end'
  colnames(regionRes)[6] <- 'region_length'
  colnames(regionRes)[14] <- 'gene_id'
  if(start0based){
    regionRes$region_begin <- pmax(regionRes$region_begin-1,0)
  }
  rownames(regionRes) <- NULL
  return(regionRes[order(regionRes$chromosome,regionRes$region_begin),])
}
