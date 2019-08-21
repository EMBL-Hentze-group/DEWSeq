# R function(s) to extract regions from significant windows
# Author: Sudeep Sahadevan, sudeep.sahadevan@embl.de

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
  requiredCols <- c('chromosome','unique_id','begin','end','strand','gene_id','gene_name',
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
  sigDat <- windowRes[ windowRes[,padjCol]<=padjThresh & windowRes[,log2FoldChangeCol]>=log2FoldChangeThresh, ]
  if(nrow(sigDat)==0){
    message('There are no significant windows/regions under the current threshold!\nPlease lower your significance cut-off thresholds and manually check if there are any significant windows under the threshold')
    return(NULL)
  }
  geneRange <- GenomicRanges::makeGRangesFromDataFrame(sigDat,seqnames.field = 'gene_id',start.field = 'begin',end.field = 'end',strand.field = 'strand',
                                                      ignore.strand=FALSE,starts.in.df.are.0based=FALSE,keep.extra.columns = TRUE)
  geneRange <-  GenomeInfoDb::sortSeqlevels(geneRange)
  geneRange <- BiocGenerics::sort(geneRange)
  geneReduce <- GenomicRanges::reduce(geneRange,drop.empty.ranges=TRUE,with.revmap=TRUE,min.gapwidth=1)
  mcols(geneReduce)$regionStartId <- mcols(geneReduce)$gene_name <- mcols(geneReduce)$gene_type <- mcols(geneReduce)$gene_region  <- 'Undefined'
  mcols(geneReduce)$chromosome <- mcols(geneReduce)$unique_ids  <- 'Undefined'
  mcols(geneReduce)$windows_in_region <- mcols(geneReduce)$Nr_of_region <-  mcols(geneReduce)$Total_nr_of_region <- mcols(geneReduce)$window_number <- 1
  mcols(geneReduce)$padj_min <- mcols(geneReduce)$padj_max <- mcols(geneReduce)$padj_avg <- 1.0
  mcols(geneReduce)$log2FoldChange_min <- mcols(geneReduce)$log2FoldChange_max <- mcols(geneReduce)$log2FoldChange_avg <- 0.0
  pb <- utils::txtProgressBar(min=0,max=length(geneReduce),style = 3)
  for(i in c(1:length(geneReduce))){
    # values
    mergeInd <- unlist(mcols(geneReduce)[i,1])
    mcols(geneReduce)[i,'windows_in_region'] <- length(mergeInd)
    mcols(geneReduce)[i,'padj_min'] <- min(mcols(geneRange)[mergeInd,padjCol])
    mcols(geneReduce)[i,'padj_max'] <- max(mcols(geneRange)[mergeInd,padjCol])
    mcols(geneReduce)[i,'padj_avg'] <- mean(mcols(geneRange)[mergeInd,padjCol])
    mcols(geneReduce)[i,'log2FoldChange_min'] <- min(mcols(geneRange)[mergeInd,log2FoldChangeCol])
    mcols(geneReduce)[i,'log2FoldChange_max'] <- max(mcols(geneRange)[mergeInd,log2FoldChangeCol])
    mcols(geneReduce)[i,'log2FoldChange_avg'] <- mean(mcols(geneRange)[mergeInd,log2FoldChangeCol])
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
    utils::setTxtProgressBar(pb=pb,value=i)
  }
  close(pb)
  regionRes <- as.data.frame(geneReduce)
  regionRes <- regionRes[,c('chromosome','start','end','strand','windows_in_region','width','padj_min','padj_avg','padj_max',
                            'log2FoldChange_min','log2FoldChange_avg','log2FoldChange_max','regionStartId','seqnames','gene_name','gene_type','gene_region',
                            'Nr_of_region','Total_nr_of_region','window_number','unique_ids')]
  colnames(regionRes)[2] <- 'region_begin'
  colnames(regionRes)[3] <- 'region_end'
  colnames(regionRes)[6] <- 'region_length'
  colnames(regionRes)[14] <- 'gene_id'
  rownames(regionRes) <- NULL
  return(regionRes[with(regionRes,order(chromosome,region_begin)),])
}


#' @export
#' @title normalized counts per region
countsPerRegion <- function(windowRes,padjCol='padj',padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1,
                                normalizedCounts,selectCol='log2FoldChange',type='max'){
  requiredCols <- c('chromosome','unique_id','begin','end','strand','gene_id','gene_name',
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
  sigDat <- windowRes[ windowRes[,padjCol]<=padjThresh & windowRes[,log2FoldChangeCol]>=log2FoldChangeThresh, ]
  rownames(sigDat) <- sigDat$unique_id
  if(nrow(sigDat)==0){
    stop('There are no significant windows/regions under the current threshold!\nPlease lower your significance cut-off thresholds and manually check if there are any significant windows under the threshold')
    # return(NULL)
  }
  if(! selectCol %in% colnames(windowRes)){
    stop('Cannot find "',selectCol,'" in windowRes column names. MUST be one of ',paste(colnames(windowRes),collapse=", "),'')
  }
  if(! is.numeric(windowRes[,selectCol])){
    stop('"',selectCol,'" values must be numeric!')
  }
  ops <- match.arg(type,c('min','max'),several.ok = FALSE)
  callFn <- ifelse(ops=='min',which.min,which.max)
  normalizedCounts <- as.data.frame(na.omit(normalizedCounts))
  commonids <- intersect(rownames(normalizedCounts),rownames(windowRes))
  if(length(commonids)==0){
    stop('There are no common elements between significant windows/regions in under the current threshold and windows in normalizedCounts data.frame!\nPlease lower your significance cut-off thresholds and manually check if there are any significant windows under the threshold')
  }
  sigDat <- cbind(sigDat[commonids,],normalizedCounts[commonids,])
  sigRange <- GenomicRanges::makeGRangesFromDataFrame(sigDat,seqnames.field = 'gene_id',start.field = 'begin',end.field = 'end',strand.field = 'strand',
                                                       ignore.strand=FALSE,starts.in.df.are.0based=FALSE,keep.extra.columns = TRUE)
  sigRange <-  GenomeInfoDb::sortSeqlevels(sigRange)
  sigRange <- BiocGenerics::sort(sigRange)
  sigReduce <- GenomicRanges::reduce(sigRange,drop.empty.ranges=TRUE,with.revmap=TRUE,min.gapwidth=1)
  if(nrow(mcols(sigReduce))<1){
    stop('Cannot find overlapping windows (regions) in input results!')
  }
  outCols <- c('gene_id','region_begin','region_end','width','strand','regionStartId','unique_id','chromosome','begin','end','gene_name','gene_region',
               'Nr_of_region','Total_nr_of_region','window_number')+ colnames(normalizedCounts)
  outDat <- cbind(as.data.frame(sigReduce)[,c(1:5)],data.frame(matrix(data=NA,nrow=nrow(mcols(sigReduce)),ncol = length(outCols)-5)))
  colnames(outDat) <- outCols
  for(i in c(1:length(sigReduce))){
    mergeInd <- sigReduce$revmap[[i]]
    outDat[i,'regionStartId'] <- mcols(sigRange)[min(mergeInd),'unique_id']
    outDat[i,'unique_id'] <- mcols(sigRange)[min(mergeInd),'unique_id']
    outDat[i,'chromosome'] <- mcols(sigRange)[min(mergeInd),'chromosome']
    outDat[i,'gene_name'] <- mcols(sigRange)[min(mergeInd),'gene_name']
    outDat[i,'gene_type'] <- mcols(sigRange)[min(mergeInd),'gene_type']
    outDat[i,'gene_region'] <- mcols(sigRange)[min(mergeInd),'gene_region']
    outDat[i,'Nr_of_region'] <- mcols(sigRange)[min(mergeInd),'Nr_of_region']
    outDat[i,'Total_nr_of_region'] <- mcols(sigRange)[min(mergeInd),'Total_nr_of_region']
    #outDat[i,'window_number'] <- mcols(sigRange)[min(mergeInd),'window_number']
    outDat[i,colnames(normalizedCounts)] <- unlist(mcols(sigRange)[callFn(mcols(sigRange)[mergeInd,selectCol]),colnames(normalizedCounts)])
  }
  return(outDat)
}

