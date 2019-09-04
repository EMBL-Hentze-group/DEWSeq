# R functions to create region from significant windows
# Author: Sudeep Sahadevan, sudeep.sahadevan@embl.de


#' @title: helper function, create region
#' @description: merge significant neighboring windows into regions
#' @param mergeDat a gene specific result data.frame
#' @param padjCol name of the adjusted pvalue column
#' @param padjThresh threshold for p-adjusted value
#' @param log2FoldChangeCol name of the log2foldchange column
#' @param log2FoldChangeThresh threshold for log2foldchange value
#' @param ncores number of cores to uses
#' @keywords internal
createRegionsHelper <- function(mergeDat,padjCol,padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1){
  mergeDat <- mergeDat[with(mergeDat,order(begin,end)),]
  rownames(mergeDat) <- NULL
  mergeDat[,'regionStartId'] <- as.character(mergeDat[,'regionStartId'])
  mergeDat[,'unique_id'] <- as.character(mergeDat[,'unique_id'])
  mergeRes <- data.frame(regionStartId=mergeDat[,'unique_id'],region_begin=mergeDat[,'begin'],region_end=mergeDat[,'end'],padjMin=mergeDat[,padjCol],
            padjMax=mergeDat[,padjCol], log2FoldChangeMin=mergeDat[,log2FoldChangeCol],log2FoldChangeMax=mergeDat[,log2FoldChangeCol],stringsAsFactors = FALSE)
  mergeDat$row_numbers <- c(1:nrow(mergeDat))
  sigDat <- mergeDat[ mergeDat[,padjCol]<= padjThresh & mergeDat[,log2FoldChangeCol]>=log2FoldChangeThresh, c('row_numbers','regionStartId') ]
  if(nrow(sigDat)>0){
    regionIds <- unique(sigDat[,'regionStartId'])
    sigList <- vector('list',nrow(sigDat))
    indStart <- 1
    for(i in c(1:nrow(sigDat))){
      if(i==1){
        sigList[[indStart]] <- c(sigDat[i,'row_numbers'])
      }else if((sigDat[i,'row_numbers']-sigDat[i-1,'row_numbers']==1) & (sigDat[i,'regionStartId']==sigDat[i-1,'regionStartId']) ){
        sigList[[indStart]] <- c(sigList[[indStart]],sigDat[i,'row_numbers'])
      }else{
        indStart <- indStart+1
        sigList[[indStart]] <- c(sigDat[i,'row_numbers'])
      }
    }
    sigList <- sigList[lapply(sigList, length)>0]
    for(i in c(1:length(sigList))){
      mergeInd <- sort(sigList[[i]])
      mergeRes[mergeInd,'regionStartId'] <- mergeDat[mergeInd[1],'unique_id']
      if(length(mergeInd)==1){
        next
      }
      mergeRes[mergeInd,'region_begin'] <- min(mergeDat[mergeInd,'begin'])
      mergeRes[mergeInd,'region_end'] <- max(mergeDat[mergeInd,'end'])
      mergeRes[mergeInd,'padjMin'] <- min(mergeDat[mergeInd,padjCol])
      mergeRes[mergeInd,'padjMax'] <- max(mergeDat[mergeInd,padjCol])
      mergeRes[mergeInd,'log2FoldChangeMin'] <- round(min(mergeDat[mergeInd,log2FoldChangeCol]),3)
      mergeRes[mergeInd,'log2FoldChangeMax'] <- round(max(mergeDat[mergeInd,log2FoldChangeCol]),3)
    }
  }else{
    mergeRes[,'regionStartId'] <- mergeDat[,'unique_id']
  }
  mergeRes <- cbind(mergeDat[,-which(colnames(mergeDat) %in% c('regionStartId','row_numbers'))],mergeRes)
  return(mergeRes)
}


