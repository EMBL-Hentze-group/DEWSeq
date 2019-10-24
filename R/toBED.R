# This file contains helper functions
# Author: Sudeep Sahadevan, sudeep.sahadevan@embl.de

#' @export
#' @title windows/regions to BED
#' @description given output of \code{\link{extractRegions}}, \code{\link{resultsDEWSeq}} and significance thresholds,
#' extract significant windows, create regions by merging adjacent significant windows.
#' Finally, write the output as a BED file for visualization.
#'
#' @param windowRes \code{data.frame}, output from \code{\link{resultsDEWSeq}}
#' @param regionRes \code{data.frame}, output from \code{\link{extractRegions}}
#' @param fileName \code{character}, filename to save BED output
#' @param padjCol \code{character}, name of the adjusted pvalue column (default: padj)
#' @param padjThresh \code{numeric}, threshold for p-adjusted value (default: 0.05)
#' @param log2FoldChangeCol \code{character}, name of the log2foldchange column (default: log2FoldChange)
#' @param log2FoldChangeThresh \code{numeric}, threshold for log2foldchange value (default:1)
#' @param trackName \code{character}, name of this track, for visualization
#' @param description \code{character}, description of this track, for visualization
#'
#' @return write to file
#'
#' @examples
#'
#' data(slbpRegions)
#' data(slbpWindows)
#' outFile <- tempfile('SLBP_visualization.bed')
#' # the results are written to a temp file in this example
#' toBED(slbpWindows,slbpRegions,outFile,padjCol='pSlidingWindows.adj')
#'

toBED <- function(windowRes,regionRes,fileName,
                  padjCol='padj',padjThresh=0.05,
                  log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1,
                  trackName='sliding windows',description='sliding windows'){
    wRequiredCols <- c('unique_id','chromosome','begin','end','strand',
                       padjCol,log2FoldChangeCol)
    missingCols <- setdiff(wRequiredCols,colnames(windowRes))
    if(length(missingCols)>0){
    stop('Input "windowRes" data.frame is missing required columns, needed columns:
        chromosome: chromosome name
        unique_id: unique id of the window
        strand: strand
        begin: window start co-ordinate
        end: window end co-ordinate
        strand: strand
     ',padjCol,': p-adjusted value column
     ',log2FoldChangeCol,': log2foldchange column.
      Missing columns: ',paste(missingCols,collapse=", "),'')
    }
  rRequiredCols <- c('chromosome','regionStartId','region_begin','region_end',
                     'strand','padj_mean')
  missingCols <- setdiff(rRequiredCols,colnames(regionRes))
  if(length(missingCols)>0){
    stop('Input "regionRes" data.frame is missing required columns, needed columns:
        chromosome: chromosome name
        regionStartId: unique_id of the left most overlapping window
        region_begin: region start co-ordinate
        region_end: region end co-ordinate
        strand: strand
        windows_in_region: total number of windows in the given region
        padj_mean: avg. padj value in the region.
      Missing columns: ',paste(missingCols,collapse=", "),'')
  }
  windowRes <- as.data.frame(windowRes)
  windowRes$unique_id <- as.character(windowRes$unique_id)
  regionRes$regionStartId <- paste(as.character(regionRes$regionStartId),'region',sep='@')
  windowRes <- windowRes[ windowRes[,padjCol]<=padjThresh &
                            windowRes[,log2FoldChangeCol]>=log2FoldChangeThresh, ]
  if(nrow(windowRes)==0){
      stop('There are no significant windows/regions under the current threshold!')
  }
  windowRes <- windowRes[order(windowRes$chromosome,windowRes$begin),]
  windowRes$strand <- as.character(windowRes$strand)
  # RGB color space for a red color, to show windows/regions as up regulated
  itemRGB <- '237,28,36' # #ed1c24
  windowRes$RGB <- itemRGB
  windowRes <- windowRes[,c('chromosome','begin','end','unique_id',
                            padjCol,'strand','begin','end','RGB')]
  windowsInRegion <- sum(regionRes$windows_in_region>1)
  if(windowsInRegion>0){
    regionRes <- regionRes[ regionRes$windows_in_region>1, ]
    regionRes <- regionRes[order(regionRes$chromosome,regionRes$region_begin),]
    regionRes$strand <- as.character(regionRes$strand)
    regionRes$RGB <- itemRGB
    regionRes <- regionRes[,c('chromosome','region_begin','region_end','regionStartId',
                              'padj_mean','strand','region_begin','region_end','RGB')]
  }

  # convert pvalues to bed scores with range 100-1000
  bedrange <- c(100,1000) # min max values for bed score
  allMin <- 1-padjThresh
  allMax <- 1
  # windows
  windowRes[,padjCol] <- 1-windowRes[,padjCol]
  windowRes[,padjCol] <- round(((windowRes[,padjCol]-allMin)/(allMax-allMin))*(bedrange[2]-bedrange[1])+bedrange[1])
  write(sprintf('track name="%s" description="%s" visibility=2 itemRgb="On" useScore=1',trackName,description),file=fileName)
  # regions
  if(windowsInRegion>0){
    rmin <- floor(min(regionRes$padj_mean))
    regionRes$padj_mean <- 1-regionRes$padj_mean
    regionRes$padj_mean <- round(((regionRes$padj_mean-allMin)/(allMax-allMin))*(bedrange[2]-bedrange[1])+bedrange[1])
    # now start writing
    for(i in seq_len(nrow(regionRes))){
      write(paste(regionRes[i,],collapse="\t"),file=fileName,append=TRUE)
    }
  }
  for(i in seq_len(nrow(windowRes))){
      write(paste(windowRes[i,],collapse="\t"),file=fileName,append=TRUE)
  }
}

