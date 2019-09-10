
#' @export
#'
#' @importFrom BiocGenerics sort
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
#' @importFrom S4Vectors na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @title stats for the top windows in each region
#' @description given window resutls and normalized counts, combine significant overlapping windows into regions and
#' for each region, pick two candidate winodws:
#' \enumerate{
#'  \item with highest log2FoldChange and
#'  \item with highest normalized mean in treatment samples (see parameter \code{treatmentCols})
#' }
#' Return a data.frame with region information and stats, and for the selected windows, the following information:
#' \itemize{
#'  \item \code{unique_id} of the window
#'  \item start and end co-ordinates
#'  \item log2FoldChange
#'  \item normalized mean expression in treatment and control samples and
#'  \item individual normalized expression in replicates
#' }
#'
#' @details
#' The output data.frame of this function has the following columns:
#' \itemize{
#'  \item \code{chromosome}: chromosome name
#'  \item \code{gene_id}: gene id
#'  \item \code{gene_name}: gene name
#'  \item \code{gene_region}: gene region
#'  \item \code{gene_type}: gene type annotation
#'  \item \code{regionStartId}: \code{unique_id} of the left most window, where a enriched region begins
#'  \item \code{region_begin}: start position of the enriched region
#'  \item \code{region_end}: end position of the enriched region
#'  \item \code{region_length}: length of the enrched region
#'  \item \code{strand}: strand info
#'  \item \code{Nr_of_region}: number of the current region
#'  \item \code{Total_nr_of_region}: total number of regions
#'  \item \code{log2FoldChange_min}: min. log 2 fold change in the region
#'  \item \code{log2FoldChange_mean}: average log 2 fold change in the region
#'  \item \code{log2FoldChange_max}: max. log 2 fold change in the region
#'  \item \code{unique_id.log2FCWindow}: unique_id of the window with largest log2FoldChange
#'  \item \code{begin.log2FCWindow}: start position of the window with largest log2FoldChange
#'  \item \code{end.log2FCWindow}: end of the window with largest log2FoldChange
#'  \item \code{log2FoldChange.log2FCWindow}: log2FoldChange  of the window with largest log2FoldChange
#'  \item \code{treatmentName.mean.log2FCWindow}: mean of the normalized expression of the treatment samples for log2FCWindow, names in \code{treatmentCols} are used to calculate mean and treatmentName is from the parameter \code{treatmentName}
#'  \item \code{controlName.mean.log2FCWindow}: mean of the normalized expression of the control samples for log2FCWindow, \code{colnames(normalizedCounts)} not found in \code{treatmentCols} are used to calculate mean and controlName is from the parameter \code{controlName}
#'  \item the next columns will be normalized expression values  of the log2FCWindow from individual treatment and control samples.
#'  \item \code{unique_id.meanWindow}:  unique_id of the window with largest mean in all treatment samples from \code{treatmentCols}
#'  \item \code{begin.meanWindow}: start position of the mean window
#'  \item \code{end.meanWindow}:  end position of the mean window
#'  \item \code{log2FoldChange.meanWindow}:log2FoldChange of the mean window
#'  \item \code{treatmentName.mean.meanWindow}: mean of the normalized expression of the treatment samples for meanWindow, names in \code{treatmentCols} are used to calculate mean and treatmentName is from the parameter \code{treatmentName}
#'  \item \code{controlName.mean.meanWindow}: mean of the normalized expression of the control samples for log2FCWindow, \code{colnames(normalizedCounts)} not found in \code{treatmentCols} are used to calculate mean and controlName is from the parameter \code{controlName}
#'  \item the next columns will be normalized expression values  of the meanWindow from individual treatment and control samples
#' }
#'
#' @param windowRes output data.frame from \code{\link{resultsDEWSeq}}
#' @param padjCol name of the adjusted pvalue column (default: padj)
#' @param padjThresh threshold for p-adjusted value (default: 0.05)
#' @param log2FoldChangeCol name of the log2foldchange column (default: log2FoldChange)
#' @param log2FoldChangeThresh threshold for log2foldchange value (default:1)
#' @param begin0based TRUE (default) or FALSE. If TRUE, then the start positions in \code{windowRes} is  considered to be 0-based
#' @param normalizedCounts data.frame with normalized read counts per window. \code{rownames(normalizedCounts)} and \code{unique_id} column from \code{windoeRes} must match
#' @param treatmentCols column names in \code{normalizedCounts} for treatment/case samples. The remaining columns in the data.frame will be considered control samples
#' @param treatmentName treatment name, see Details  (default: treatment)
#' @param controlName control name, see Details (default: control)
#' @param op can be one of \code{max} (default) or \code{min}. \code{max} returns windows with maximum log2FoldChange and mean normalized expression in the \code{treatmentCols} columns,
#' \code{min} returns windows with minimum log2FoldChange and mean normalized expression
#'
#' @return data.frame
topWindowStats <- function(windowRes,padjCol='padj',padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1,begin0based=TRUE, normalizedCounts,
                           treatmentCols,treatmentName='treatment',controlName='control', op='max'){
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
  if (length(unique(treatmentCols))<length(treatmentCols)){
    stop('There are repeating names in treatmentCols')
  }
  if(length(setdiff(treatmentCols,colnames(normalizedCounts)))>0){
    stop('mismatch treatmentCols and column names of normalizedCounts!')
  }
  controlCols <- setdiff(colnames(normalizedCounts),treatmentCols)
  op <- match.arg(op,c('min','max'),several.ok = FALSE)
  callFn <- which.max
  if(op=='min'){
    callFn <- which.min
  }
  normalizedCounts <- as.data.frame(na.omit(normalizedCounts))
  commonids <- intersect(rownames(normalizedCounts),rownames(sigDat))
  if(length(commonids)==0){
    stop('There are no common elements between significant windows/regions in under the current threshold and windows in normalizedCounts data.frame!\nPlease lower your significance cut-off thresholds and manually check if there are any significant windows under the threshold')
  }
  sigDat <- cbind(sigDat[commonids,],normalizedCounts[commonids,])
  sigRange <- makeGRangesFromDataFrame(sigDat,seqnames.field = 'gene_id',start.field = 'begin',end.field = 'end',strand.field = 'strand',
                                                      ignore.strand=FALSE,starts.in.df.are.0based=begin0based,keep.extra.columns = TRUE)
  sigRange <-  sortSeqlevels(sigRange)
  sigRange <- sort(sigRange)
  sigReduce <- reduce(sigRange,drop.empty.ranges=TRUE,with.revmap=TRUE,min.gapwidth=1)
  if(nrow(mcols(sigReduce))<1){
    stop('Cannot find overlapping windows (regions) in input results!')
  }
  log2FCWindowCols <- paste(colnames(normalizedCounts),'log2FCWindow',sep='.')
  log2FCMean <- paste(c(treatmentName,controlName),'mean.log2FCWindow',sep='.')
  meanWindowCols <- paste(colnames(normalizedCounts),'meanWindow',sep='.')
  meanMean <- paste(c(treatmentName,controlName),'mean.meanWindow',sep='.')
  outCols <- c(c('gene_id','region_begin','region_end','region_length','strand','regionStartId','chromosome','gene_name','gene_region','gene_type',
                 'Nr_of_region','Total_nr_of_region',paste0(log2FoldChangeCol,'_min'),paste0(log2FoldChangeCol,'_mean'),paste0(log2FoldChangeCol,'_max'),
                 'unique_id.log2FCWindow','begin.log2FCWindow','end.log2FCWindow',paste0(log2FoldChangeCol,'.log2FCWindow')),log2FCMean, log2FCWindowCols,
               c('unique_id.meanWindow','begin.meanWindow','end.meanWindow',paste0(log2FoldChangeCol,'.meanWindow')),meanMean,meanWindowCols)
  outDat <- cbind(as.data.frame(sigReduce)[,c(1:5)],data.frame(matrix(data=NA,nrow=nrow(mcols(sigReduce)),ncol = length(outCols)-5)))
  colnames(outDat) <- outCols
  # now fill out 0s for the mean window columns
  outDat$unique_id.meanWindow <- 'same window'
  outDat[,c('begin.meanWindow','end.meanWindow',meanMean,meanWindowCols)] <- 0
  # mean, log2foldchange comparison data.frame
  pb <- txtProgressBar(min=0,max=length(sigReduce),style = 3)
  for(i in seq_len(length(sigReduce))){
    mergeInd <- sigReduce$revmap[[i]]
    outDat[i,'regionStartId'] <- mcols(sigRange)[min(mergeInd),'unique_id']
    outDat[i,'regionStartId'] <- mcols(sigRange)[min(mergeInd),'unique_id']
    outDat[i,'chromosome'] <- mcols(sigRange)[min(mergeInd),'chromosome']
    outDat[i,'gene_name'] <- mcols(sigRange)[min(mergeInd),'gene_name']
    outDat[i,'gene_type'] <- mcols(sigRange)[min(mergeInd),'gene_type']
    outDat[i,'gene_region'] <- mcols(sigRange)[min(mergeInd),'gene_region']
    outDat[i,'gene_type'] <- mcols(sigRange)[min(mergeInd),'gene_type']
    outDat[i,'Nr_of_region'] <- mcols(sigRange)[min(mergeInd),'Nr_of_region']
    outDat[i,'Total_nr_of_region'] <- mcols(sigRange)[min(mergeInd),'Total_nr_of_region']
    #get the window with minimum/maximum log2 Foldchange
    # fill log2FoldChange first
    outDat[i,paste0(log2FoldChangeCol,'_min')] <- min(unlist(mcols(sigRange)[mergeInd,log2FoldChangeCol]))
    outDat[i,paste0(log2FoldChangeCol,'_mean')] <- mean(unlist(mcols(sigRange)[mergeInd,log2FoldChangeCol]))
    outDat[i,paste0(log2FoldChangeCol,'_max')] <- max(unlist(mcols(sigRange)[mergeInd,log2FoldChangeCol]))
    log2FCInd <- callFn(mcols(sigRange)[mergeInd,log2FoldChangeCol])
    outDat[i,'unique_id.log2FCWindow'] <- mcols(sigRange)[mergeInd[log2FCInd],'unique_id']
    outDat[i,'begin.log2FCWindow'] <- start(sigRange[mergeInd[log2FCInd]])
    outDat[i,'end.log2FCWindow'] <- end(sigRange[mergeInd[log2FCInd]])
    outDat[i,paste0(log2FoldChangeCol,'.log2FCWindow')] <- mcols(sigRange)[mergeInd[log2FCInd],log2FoldChangeCol]
    outDat[i,paste0(treatmentName,'.mean.log2FCWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[log2FCInd],treatmentCols]))
    outDat[i,paste0(controlName,'.mean.log2FCWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[log2FCInd],controlCols]))
    outDat[i,log2FCWindowCols] <- unlist(mcols(sigRange)[mergeInd[log2FCInd],colnames(normalizedCounts)])
    # fill out log2FC window data
    # get the window with minimum/maximum normalized count
    meanInd <- callFn( rowMeans( as.data.frame(mcols(sigRange)[mergeInd,treatmentCols]) ) )
    if(log2FCInd!=meanInd){
      outDat[i,'unique_id.meanWindow'] <- mcols(sigRange)[mergeInd[meanInd],'unique_id']
    }
    outDat[i,'begin.meanWindow'] <- start(sigRange[mergeInd[meanInd]])
    outDat[i,'end.meanWindow'] <- end(sigRange[mergeInd[meanInd]])
    outDat[i,paste0(log2FoldChangeCol,'.meanWindow')] <- mcols(sigRange)[mergeInd[meanInd],log2FoldChangeCol]
    outDat[i,paste0(treatmentName,'.mean.meanWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[meanInd],treatmentCols]))
    outDat[i,paste0(controlName,'.mean.meanWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[meanInd],controlCols]))
    outDat[i,meanWindowCols] <- unlist(mcols(sigRange)[mergeInd[meanInd],colnames(normalizedCounts)])
    # fill out mean window data
    # progress
    setTxtProgressBar(pb=pb,value=i)
  }
  close(pb)
  if(begin0based){ # return 0 based start positions
    outDat$region_begin <- outDat$region_begin-1
    outDat$begin.log2FCWindow <- outDat$begin.log2FCWindow-1
    outDat$begin.meanWindow <- pmax(outDat$begin.meanWindow-1,0)
  }
  colOrder <- c(c('chromosome','gene_id','gene_name','gene_region','gene_type','regionStartId','region_begin','region_end','width','strand','Nr_of_region','Total_nr_of_region',
                  paste0('min.',log2FoldChangeCol),paste0('mean.',log2FoldChangeCol),paste0('max.',log2FoldChangeCol),'unique_id.log2FCWindow','begin.log2FCWindow','end.log2FCWindow',
                  paste0(log2FoldChangeCol,'.log2FCWindow')),log2FCMean,log2FCWindowCols,c('unique_id.meanWindow','begin.meanWindow','end.meanWindow',paste0(log2FoldChangeCol,'.meanWindow')),
                meanMean,meanWindowCols)
  # sort main results
  outDat <- outDat[order(outDat$chromosome,outDat$region_begin),] # sort as a separate step to fix
  rownames(outDat) <- NULL # fix jumbled up rownames
  return(outDat[,colOrder])
}
