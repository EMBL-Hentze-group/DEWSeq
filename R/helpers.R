# This file contains helper functions

#' @title windows/regions to BED
#' @description given output of \code{\link{createRegions}} and significance thresholds, extract significant windows, 
#' regions (merged neighboring significant windows) and create a BED file for visualization
#' @param windowRes output data.frame from \code{\link{createRegions}}
#' @param padjCol name of the adjusted pvalue column (default: padj)
#' @param padjThresh threshold for p-adjusted value (default: 0.05)
#' @param log2FoldChangeCol name of the log2foldchange column (default: log2FoldChange)
#' @param log2FoldChangeThresh threshold for log2foldchange value (default:1)
#' @param trackName name of this track, for visualization
#' @param description description of this track, for visualization
#' @return NULL
toBED <- function(windowRes,fileName,padjCol='padj',padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1,trackName='sliding windows', 
                                                                                                                            description='sliding windows'){
    requiredCols <- c('unique_id','chromosome','begin','end','strand','regionStartId','region_begin','region_end',padjCol,log2FoldChangeCol)
    missingCols <- setdiff(requiredCols,colnames(windowRes))
    if(length(missingCols)>0){
    stop('Input data.frame is missing required columns, needed columns:
        chromosome: chromosome name
        unique_id: unique id of the window
        strand: strand
        begin: window start co-ordinate
        end: window end co-ordinate
        region_begin: region start co-ordinate
        region_end: region end co-ordinate
        strand: strand
        gene_id: gene id
        regionStartId: unique_id of the left most overlapping window
     ',padjCol,': p-adjusted value column
     ',log2FoldChangeCol,': log2foldchange column.
      Missing columns: ',paste(missingCols,collapse=", "),'')
  }
  windowRes$unique_id <- as.character(windowRes$unique_id)
  windowRes$regionStartId <- as.character(windowRes$regionStartId)
  sigDat <- windowRes[ windowRes[,padjCol]<=padjThresh & windowRes[,log2FoldChangeCol]>=log2FoldChangeThresh, ]
  if(nrow(sigDat)==0){
      stop('There are no significant windows/regions under the current threshold!')
  }
  sigWindows <- sigDat[,c('chromosome','begin','end','unique_id',padjCol,'strand')]
  rownames(sigWindows) <- sigWindows$unique_id
  sigRegionIds <- unique(sigDat[,'regionStartId'])
  sigRegionCols <- c('chromosome','begin','end','unique_id',padjCol,'strand')
  sigRegions <- data.frame(matrix(0,nrow=length(sigRegionIds),ncol=length(sigRegionCols),dimnames=list(sigRegionIds,sigRegionCols)),stringsAsFactors=FALSE)
  for(sr in sigRegionIds){
      sigRegions[sr,'chromosome'] <- unique(sigDat[sigDat[,'regionStartId'] == sr,'chromosome'])[1]
      sigRegions[sr,'unique_id'] <- unique(sigDat[sigDat[,'regionStartId'] == sr,'regionStartId'])[1]
      sigRegions[sr,'strand'] <- unique(sigDat[sigDat[,'regionStartId'] == sr,'strand'])[1]
      sigRegions[sr,'begin'] <- min(sigDat[sigDat[,'regionStartId'] == sr,'begin'])
      sigRegions[sr,'end'] <- max(sigDat[sigDat[,'regionStartId'] == sr,'end'])
      sigRegions[sr,padjCol] <- mean(sigDat[sigDat[,'regionStartId'] == sr,padjCol])
  }
  rownames(sigRegions) <- paste(sigRegions[,'unique_id'],'region',sep='@')
  sigRegions <- unique(rbind(sigWindows,sigRegions))
  regions <- grep(pattern = '\\@region',x = rownames(sigRegions))
  sigRegions[regions,'unique_id'] <- paste(sigRegions[regions,'unique_id'],'region',sep='@')
  sigRegions <- sigRegions[with(sigRegions,order(chromosome,begin)),]
  # convert pvalues to bed scores with range 100-1000
  tprobs <- 1-sigRegions[,padjCol]
  tmin <- 1-padjThresh
  bedrange <- c(100,1000) # min max values for bed score
  sigRegions[,padjCol] <- round(((tprobs-tmin)/(max(tprobs)-tmin))*(bedrange[2]-bedrange[1])+bedrange[1])
  # RGB color space for a red color, to show windows/regions as up regulated
  itemRGB <- '237,28,36' # #ed1c24
  # convert to BED 9 format
  sigRegions <- data.frame(sigRegions,col7=sigRegions[,'begin'],col8=sigRegions[,'end'],col9=rep(itemRGB,nrow(sigRegions)),stringsAsFactors=FALSE)
  # now start writing
  write(sprintf('track name="%s" description="%s" visibility=2 itemRgb="On" useScore=1',trackName,description),file=fileName)
  for(i in c(1:nrow(sigRegions))){
      write(paste(sigRegions[i,],collapse="\t"),file=fileName,append=TRUE)
  }
}