# This file contains helper functions
# Author: Sudeep Sahadevan, sudeep.sahadevan@embl.de

#' @export
#' @title windows/regions to BED
#' @description given output of \code{\link{extractRegions}}, \code{\link{resultsDEWSeq}} and significance thresholds,
#' extract significant windows, create regions by merging adjacent significant windows.
#' Finally, write the output as a BED file for visualization.
#'
#' @param windowRes output data.frame from \code{\link{resultsDEWSeq}}
#' @param regionRes output data.frame from \code{\link{extractRegions}}
#' @param fileName filename to save BED output
#' @param padjCol name of the adjusted pvalue column (default: padj)
#' @param padjThresh threshold for p-adjusted value (default: 0.05)
#' @param log2FoldChangeCol name of the log2foldchange column (default: log2FoldChange)
#' @param log2FoldChangeThresh threshold for log2foldchange value (default:1)
#' @param trackName name of this track, for visualization
#' @param description description of this track, for visualization
#'
#' @examples
#' # need specific examples
#' \dontrun{
#' 'toBED(windowRes=windowRes,regionRes=regionRes,fileName="enrichedWindowsRegions.bed")'
#' }
#'
#' @return NULL
toBED <- function(windowRes,regionRes,fileName,padjCol='padj',padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1,trackName='sliding windows',
                                                                                                                            description='sliding windows'){
    wRequiredCols <- c('unique_id','chromosome','begin','end','strand',padjCol,log2FoldChangeCol)
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
  rRequiredCols <- c('chromosome','regionStartId','region_begin','region_end','strand','padj_avg')
  missingCols <- setdiff(rRequiredCols,colnames(regionRes))
  if(length(missingCols)>0){
    stop('Input "regionRes" data.frame is missing required columns, needed columns:
        chromosome: chromosome name
        regionStartId: unique_id of the left most overlapping window
        region_begin: region start co-ordinate
        region_end: region end co-ordinate
        strand: strand
        windows_in_region: total number of windows in the given region
        padj_avg: avg. padj value in the region.
      Missing columns: ',paste(missingCols,collapse=", "),'')
  }
  windowRes$unique_id <- as.character(windowRes$unique_id)
  regionRes$regionStartId <- paste(as.character(regionRes$regionStartId),'region',sep='@')
  windowRes <- windowRes[ windowRes[,padjCol]<=padjThresh & windowRes[,log2FoldChangeCol]>=log2FoldChangeThresh, ]
  if(nrow(windowRes)==0){
      stop('There are no significant windows/regions under the current threshold!')
  }
  regionRes <- regionRes[ regionRes$windows_in_region>1, ]
  windowRes <- windowRes[with(windowRes,order(chromosome,begin)),]
  regionRes <- regionRes[with(regionRes,order(chromosome,region_begin)),]
  # RGB color space for a red color, to show windows/regions as up regulated
  itemRGB <- '237,28,36' # #ed1c24
  windowRes$RGB <- itemRGB
  regionRes$RGB <- itemRGB
  windowRes <- windowRes[,c('chromosome','begin','end','unique_id',padjCol,'strand','begin','end','RGB')]
  regionRes <- regionRes[,c('chromosome','region_begin','region_end','regionStartId','padj_avg','strand','region_begin','region_end','RGB')]
  # convert pvalues to bed scores with range 100-1000
  bedrange <- c(100,1000) # min max values for bed score
  allMin <- 1-padjThresh
  allMax <- 1
  # windows
  windowRes[,padjCol] <- 1-windowRes[,padjCol]
  windowRes[,padjCol] <- round(((windowRes[,padjCol]-allMin)/(allMax-allMin))*(bedrange[2]-bedrange[1])+bedrange[1])
  # regions
  rmin <- floor(min(regionRes$padj_avg))
  regionRes$padj_avg <- 1-regionRes$padj_avg
  regionRes$padj_avg <- round(((regionRes$padj_avg-allMin)/(allMax-allMin))*(bedrange[2]-bedrange[1])+bedrange[1])
  # now start writing
  write(sprintf('track name="%s" description="%s" visibility=2 itemRgb="On" useScore=1',trackName,description),file=fileName)
  for(i in seq_len(nrow(regionRes))){
    write(paste(regionRes[i,],collapse="\t"),file=fileName,append=TRUE)
  }
  for(i in seq_len(nrow(windowRes))){
      write(paste(windowRes[i,],collapse="\t"),file=fileName,append=TRUE)
  }
}

