#' @export
#' @title extract significant regions
#' @description extract significant windows from output of \code{\link{createRegions}} using \code{regionStartId} column,
#' merge these significant windows to regions and create the following columns for each significant region: \cr
#'   \code{padj_min}: min. padj value in the region \cr
#'   \code{padj_max}: max. padj value in the region \cr
#'   \code{padj_avg}: avg. padj value in the region \cr
#'   \code{log2FoldChange_min}: min. log 2 fold change in the region \cr
#'   \code{log2FoldChange_max}: max. log 2 fold change in the region \cr
#'   \code{log2FoldChange_avg}: avg. log 2 fold change in the region \cr
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
