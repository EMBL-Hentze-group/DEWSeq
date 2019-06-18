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
      mergeRes[mergeInd,'log2FoldChangeMin'] <- min(mergeDat[mergeInd,log2FoldChangeCol])
      mergeRes[mergeInd,'log2FoldChangeMax'] <- max(mergeDat[mergeInd,log2FoldChangeCol])
    }
  }else{
    mergeRes[,'regionStartId'] <- mergeDat[,'unique_id']
  }
  mergeRes <- cbind(mergeDat[,-which(colnames(mergeDat) %in% c('regionStartId','row_numbers'))],mergeRes)
  return(mergeRes)
}
#' @TODO: write a wrapper for DESeq2 wald test

#' @title create regions from significant windows
#' @description create significant regions by merging significant windows
#' given p-adjusted value and log2 fold change columns and thresholds
#' @param mergeDat data.frame, output from \code{\link{mergeWindows}}
#' @param padjCol name of the adjusted pvalue column (default: pBonferroni.adj)
#' @param padjThresh threshold for p-adjusted value (default: 0.05)
#' @param log2FoldChangeCol name of the log2foldchange column (default: log2FoldChange)
#' @param log2FoldChangeThresh threshold for log2foldchange value (default:1)
#' @param ncores number of cores to uses
#' @return data.frame
createRegions <- function(mergeDat,padjCol='pBonferroni.adj',padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1,ncores=5){
  requiredCols <- c('unique_id','chromosome','begin','end','strand','gene_id','regionStartId',padjCol,log2FoldChangeCol)
  missingCols <- setdiff(requiredCols,colnames(mergeDat))
  if(length(missingCols)>0){
    stop('Input data.frame is missing required columns, needed columns:
        chromosome: chromosome name
        unique_id: unique id of the window
        begin: window start co-ordinate
        end: window end co-ordinate
        strand: strand
        gene_id: gene id
        regionStartId: unique_id of the left most overlapping window
     ',padjCol,': p-adjusted value column
     ',log2FoldChangeCol,': log2foldchange column.
      Missing columns: ',paste(missingCols,collapse=", "),'')
  }

  rownames(mergeDat) <- NULL
  geneIds <- unique(mergeDat[,'gene_id'])
  geneList <- vector('list',length(geneIds))
  names(geneList) <- geneIds
  for(geneId in geneIds){
    geneList[[geneId]] <- mergeDat[ mergeDat['gene_id']==geneId, ]
  }
  return(geneList)
  register(BatchtoolsParam(workers = ncores), default = TRUE)
  regionDat <- do.call(rbind,bplapply(geneList, createRegionsHelper,padjCol=padjCol,
    padjThresh=padjThresh,log2FoldChangeCol=log2FoldChangeCol,log2FoldChangeThresh=log2FoldChangeThresh))
  rownames(regionDat) <- NULL
  message('\n')
  return(regionDat)
}
