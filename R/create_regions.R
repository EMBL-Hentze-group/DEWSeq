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

#' @title extract significant windows 
#' @description extract significant windows from output of \code{\link{createRegions}} using \code{regionStartId} column
#' and create the following columns for each significant region:
#'   padj_min: min. padj value in the region
#'   padj_max: max. padj value in the region
#'   padj_avg: avg. padj value in the region
#'   log2FoldChange_min: min. log 2 fold change in the region
#'   log2FoldChange_max: max. log 2 fold change in the region
#'   log2FoldChange_avg: avg. log 2 fold change in the region
#' @param windowRes output data.frame from \code{\link{createRegions}}
#' @param padjCol name of the adjusted pvalue column (default: padj)
#' @param padjThresh threshold for p-adjusted value (default: 0.05)
#' @param log2FoldChangeCol name of the log2foldchange column (default: log2FoldChange)
#' @param log2FoldChangeThresh threshold for log2foldchange value (default:1)
#' @return data.frame
extractRegions <- function(windowRes,padjCol='padj',padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1){
  requiredCols <- c('chromosome','unique_id','begin','end','strand','gene_id','gene_name','regionStartId',
    'gene_type','gene_region','Nr_of_region','Total_nr_of_region','window_number',padjCol,log2FoldChangeCol)
  missingCols <- setdiff(requiredCols,colnames(windowRes))
  if(length(missingCols)>0){
    stop('Input data.frame is missing required columns, needed columns:
        chromosome: chromosome name
        unique_id: unique id of the window
        begin: window start co-ordinate
        end: window end co-ordinate
        strand: strand
        gene_id: gene id
        regionStartId: unique_id of the left most overlapping window
        gene_id: gene id
        gene_name: gene name
        gene_type: gene type annotation
        gene_region: gene region
        Nr_of_region: number of the current region
        Total_nr_of_region: total number of regions
        window_number: window numbe
     ',padjCol,': p-adjusted value column
     ',log2FoldChangeCol,': log2foldchange column.
      Missing columns: ',paste(missingCols,collapse=", "),'')
  }
  windowRes <- na.omit(windowRes[,requiredCols])
  sigDat <- windowRes[ windowRes[,padjCol]<=padjThresh & windowRes[,log2FoldChangeCol]>=log2FoldChangeThresh, ]
  if(nrow(sigDat)==0){
    stop('There are no significant windows/regions under the current threshold!')
  }
  regionCols <- c('chromosome','region_begin','region_end','region_length','strand','regionStartId','windows_in_region','padj_min','padj_max','padj_avg','log2FoldChange_min',
    'log2FoldChange_max','log2FoldChange_avg','gene_id','gene_name','gene_type','gene_region','Nr_of_region','Total_nr_of_region')
  regionIds <- unique(sigDat[,'regionStartId'])
  regionDat <- data.frame(matrix(0,nrow=length(regionIds),ncol=length(regionCols),dimnames=list(regionIds,regionCols)),stringsAsFactors=FALSE)
  for(ri in regionIds){
    regionDat[ri,'chromosome'] <- unique(sigDat[sigDat[,'regionStartId'] == ri,'chromosome'])[1]
    regionDat[ri,'strand'] <- unique(sigDat[sigDat[,'regionStartId'] == ri,'strand'])[1]
    regionDat[ri,'regionStartId'] <- ri
    regionDat[ri,'gene_id'] <- unique(sigDat[sigDat[,'regionStartId'] == ri,'gene_id'])[1]
    regionDat[ri,'gene_name'] <- unique(sigDat[sigDat[,'regionStartId'] == ri,'gene_name'])[1]
    regionDat[ri,'gene_type'] <- unique(sigDat[sigDat[,'regionStartId'] == ri,'gene_type'])[1]
    regionDat[ri,'gene_region'] <- unique(sigDat[sigDat[,'regionStartId'] == ri,'gene_region'])[1]
    regionDat[ri,'Nr_of_region'] <- unique(sigDat[sigDat[,'regionStartId'] == ri,'Nr_of_region'])[1]
    regionDat[ri,'Total_nr_of_region'] <- unique(sigDat[sigDat[,'regionStartId'] == ri,'Total_nr_of_region'])[1]
    # min and max begin positions,
    # number of windows in region and total length of the region
    regionDat[ri,'region_begin'] <- min(sigDat[sigDat[,'regionStartId'] == ri,'begin'])
    regionDat[ri,'windows_in_region'] <- length(unique(sigDat[sigDat[,'regionStartId'] == ri,'begin']))
    regionDat[ri,'region_end'] <- max(sigDat[sigDat[,'regionStartId'] == ri,'end'])
    regionDat[ri,'region_length'] <- regionDat[ri,'region_end'] - regionDat[ri,'region_begin']
    # min max and mean for p-value
    pvals <- sigDat[sigDat[,'regionStartId'] == ri,padjCol]
    regionDat[ri,'padj_min'] <- min(pvals)
    regionDat[ri,'padj_max'] <- max(pvals)
    regionDat[ri,'padj_avg'] <- mean(pvals)
    # min max and mean for log2 fold changes
    log2fcs <- sigDat[sigDat[,'regionStartId'] == ri,log2FoldChangeCol]
    regionDat[ri,'log2FoldChange_min'] <- round(min(log2fcs),3)
    regionDat[ri,'log2FoldChange_max'] <- round(max(log2fcs),3)
    regionDat[ri,'log2FoldChange_avg'] <- round(mean(log2fcs),3)
  }
  return(regionDat)
}