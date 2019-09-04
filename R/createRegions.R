#' @export
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
  mergeDat <- na.omit(mergeDat)
  rownames(mergeDat) <- NULL
  geneIds <- unique(mergeDat[,'gene_id'])
  geneList <- vector('list',length(geneIds))
  names(geneList) <- geneIds
  for(geneId in geneIds){
    geneList[[geneId]] <- mergeDat[ mergeDat['gene_id']==geneId, ]
  }
  register(BatchtoolsParam(workers = ncores), default = TRUE)
  regionDat <- do.call(rbind,bplapply(geneList, createRegionsHelper,padjCol=padjCol,
                                      padjThresh=padjThresh,log2FoldChangeCol=log2FoldChangeCol,log2FoldChangeThresh=log2FoldChangeThresh))
  rownames(regionDat) <- NULL
  message('\n')
  return(regionDat)
  }
