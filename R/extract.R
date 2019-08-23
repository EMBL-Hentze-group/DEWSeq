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
#' @description given window resutls and normalized counts, combine significant overlapping windows into regions and
#' for each region, pick two candidate winodws:\cr (i) with highest log2FoldChange and\cr (ii) with highest normalized mean in
#' treatment samples.\cr Return a data.frame with region information and stats, and for the selected windows stats:\cr
#' normalized mean expression in treatment and control samples and\cr individual expression in replicates
#' @TODO comment me!
countsPerRegion <- function(windowRes,padjCol='padj',padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1,begin0based=TRUE, normalizedCounts,
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
  # if (length(unique(controlCols))<length(controlCols)){
  #   stop('There are repeating names in controlCols')
  # }
  if (length(unique(treatmentCols))<length(treatmentCols)){
    stop('There are repeating names in treatmentCols')
  }
  if(length(setdiff(treatmentCols,colnames(normalizedCounts)))>0){
    stop('mismatch treatmentCols and column names of normalizedCounts!')
  }
  # if(length(setdiff(controlCols,colnames(normalizedCounts)))>0){
  #   stop('mismatch controlCols and column names of normalizedCounts!')
  # }
  # if (length(intersect(treatmentCols,controlCols))>0){
  #   stop('treatmentCols and controlCols should not have shared names!')
  # }
  controlCols <- setdiff(colnames(normalizedCounts),treatmentCols)
  op <- match.arg(op,c('min','max'),several.ok = FALSE)
  callFn <- which.max
  if(op=='min'){
    callFn <- which.min
  }
  # callFn <- ifelse(op=='min',which.min,which.max)
  # normalizedCounts <- as.data.frame(na.omit(normalizedCounts[,c(treatmentCols,controlCols)]))
  normalizedCounts <- as.data.frame(na.omit(normalizedCounts))
  commonids <- intersect(rownames(normalizedCounts),rownames(sigDat))
  if(length(commonids)==0){
    stop('There are no common elements between significant windows/regions in under the current threshold and windows in normalizedCounts data.frame!\nPlease lower your significance cut-off thresholds and manually check if there are any significant windows under the threshold')
  }
  sigDat <- cbind(sigDat[commonids,],normalizedCounts[commonids,])
  sigRange <- GenomicRanges::makeGRangesFromDataFrame(sigDat,seqnames.field = 'gene_id',start.field = 'begin',end.field = 'end',strand.field = 'strand',
                                                       ignore.strand=FALSE,starts.in.df.are.0based=begin0based,keep.extra.columns = TRUE)
  sigRange <-  GenomeInfoDb::sortSeqlevels(sigRange)
  sigRange <- BiocGenerics::sort(sigRange)
  sigReduce <- GenomicRanges::reduce(sigRange,drop.empty.ranges=TRUE,with.revmap=TRUE,min.gapwidth=1)
  if(nrow(mcols(sigReduce))<1){
    stop('Cannot find overlapping windows (regions) in input results!')
  }
  log2FCWindowCols <- paste(colnames(normalizedCounts),'log2FCWindow',sep='.')
  log2FCMean <- paste(c(treatmentName,controlName),'mean.log2FCWindow',sep='.')
  meanWindowCols <- paste(colnames(normalizedCounts),'meanWindow',sep='.')
  meanMean <- paste(c(treatmentName,controlName),'mean.meanWindow',sep='.')
  outCols <- c(c('gene_id','region_begin','region_end','width','strand','regionStartId','chromosome','gene_name','gene_region','gene_type',
               'Nr_of_region','Total_nr_of_region',paste0('min.',log2FoldChangeCol),paste0('mean.',log2FoldChangeCol),paste0('max.',log2FoldChangeCol),
               'unique_id.log2FCWindow','begin.log2FCWindow','end.log2FCWindow',paste0(log2FoldChangeCol,'.log2FCWindow')),log2FCMean, log2FCWindowCols, 
               c('unique_id.meanWindow','begin.meanWindow','end.meanWindow',paste0(log2FoldChangeCol,'.meanWindow')),meanMean,meanWindowCols)
  outDat <- cbind(as.data.frame(sigReduce)[,c(1:5)],data.frame(matrix(data=NA,nrow=nrow(mcols(sigReduce)),ncol = length(outCols)-5)))
  colnames(outDat) <- outCols
  # now fill out 0s for the mean window columns
  outDat$unique_id.meanWindow <- 'same window'
  outDat[,c('begin.meanWindow','end.meanWindow',meanMean,meanWindowCols)] <- 0
  # mean, log2foldchange comparison data.frame
  # compDat <- data.frame(matrix(NA,nrow=nrow(mcols(sigReduce)),ncol=8))
  # colnames(compDat) <- c('regionStartId',log2FoldChangeCol,'log2FoldChange.log2FCWindow',paste0(treatmentName,'.mean.log2FCWindow'),paste0(controlName,'.mean.log2FCWindow'),
  #   'log2FoldChange.meanWindow',paste0(treatmentName,'.mean.meanWindow'),paste0(controlName,'.mean.meanWindow'))
  pb <- utils::txtProgressBar(min=0,max=length(sigReduce),style = 3)
  for(i in c(1:length(sigReduce))){
    mergeInd <- sigReduce$revmap[[i]]
    outDat[i,'regionStartId'] <- mcols(sigRange)[min(mergeInd),'unique_id']
    compDat[i,'regionStartId'] <- mcols(sigRange)[min(mergeInd),'unique_id']
    outDat[i,'chromosome'] <- mcols(sigRange)[min(mergeInd),'chromosome']
    outDat[i,'gene_name'] <- mcols(sigRange)[min(mergeInd),'gene_name']
    outDat[i,'gene_type'] <- mcols(sigRange)[min(mergeInd),'gene_type']
    outDat[i,'gene_region'] <- mcols(sigRange)[min(mergeInd),'gene_region']
    outDat[i,'gene_type'] <- mcols(sigRange)[min(mergeInd),'gene_type']
    outDat[i,'Nr_of_region'] <- mcols(sigRange)[min(mergeInd),'Nr_of_region']
    outDat[i,'Total_nr_of_region'] <- mcols(sigRange)[min(mergeInd),'Total_nr_of_region']
    #get the window with minimum/maximum log2 Foldchange
    # fill log2FoldChange first
    # compDat[i,log2FoldChangeCol] <-  mean(unlist(mcols(sigRange)[mergeInd,log2FoldChangeCol]))
    outDat[i,paste0('min.',log2FoldChangeCol)] <- min(unlist(mcols(sigRange)[mergeInd,log2FoldChangeCol]))
    outDat[i,paste0('mean.',log2FoldChangeCol)] <- mean(unlist(mcols(sigRange)[mergeInd,log2FoldChangeCol]))
    outDat[i,paste0('max.',log2FoldChangeCol)] <- max(unlist(mcols(sigRange)[mergeInd,log2FoldChangeCol]))
    log2FCInd <- callFn(mcols(sigRange)[mergeInd,log2FoldChangeCol])
    outDat[i,'unique_id.log2FCWindow'] <- mcols(sigRange)[mergeInd[log2FCInd],'unique_id']
    outDat[i,'begin.log2FCWindow'] <- start(sigRange[mergeInd[log2FCInd]])
    outDat[i,'end.log2FCWindow'] <- end(sigRange[mergeInd[log2FCInd]])
    outDat[i,paste0(log2FoldChangeCol,'.log2FCWindow')] <- mcols(sigRange)[mergeInd[log2FCInd],log2FoldChangeCol]
    outDat[i,paste0(treatmentName,'.mean.log2FCWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[log2FCInd],treatmentCols]))
    outDat[i,paste0(controlName,'.mean.log2FCWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[log2FCInd],controlCols]))
    outDat[i,log2FCWindowCols] <- unlist(mcols(sigRange)[mergeInd[log2FCInd],colnames(normalizedCounts)])
    # fill out log2FC window data
    # compDat[i,'log2FoldChange.log2FCWindow'] <- mcols(sigRange)[mergeInd[log2FCInd],log2FoldChangeCol]
    # compDat[i,paste0(treatmentName,'.mean.log2FCWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[log2FCInd],treatmentCols]))
    # compDat[i,paste0(controlName,'.mean.log2FCWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[log2FCInd],controlCols]))
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
    # compDat[i,'log2FoldChange.meanWindow'] <- mcols(sigRange)[mergeInd[log2FCInd],log2FoldChangeCol]
    # compDat[i,paste0(treatmentName,'.mean.meanWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[log2FCInd],treatmentCols]))
    # compDat[i,paste0(controlName,'.mean.meanWindow')] <- mean(unlist(mcols(sigRange)[mergeInd[log2FCInd],controlCols]))
    # progress
    utils::setTxtProgressBar(pb=pb,value=i)
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

