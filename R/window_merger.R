
# @TODO: fix repeating code!

#' @title  local simes correction
#' @description Helper function, merge overlapping windows and
#' correct p-values for these overlapping windows using local simes like method
#' @param annRes annotated results (data.frame)
#' @param minDist minimum distance between non overlapping windows,
#' default: 0 to account for bed formatted window regions
#' @keywords internal
mergeWindowsHelper_local <- function(annRes,minDist=0){
  annRes <- annRes[with(annRes,order(begin,end)),]
  annRes <- data.frame(annRes,regionStartId=annRes[,'unique_id'],pLocal=rep(1,nrow(annRes)),stringsAsFactors = FALSE)
  endW <- c(annRes[1,'begin']+1,annRes[,'end']) # end of all previous windows, except for first entry
  endW <- endW[-length(endW)] # remove last
  beDiff <- annRes[,'begin'] - endW
  # all indices where the current window do not overlap the previous one
  posInd <- which(beDiff>=minDist) # The operator MUST BE '>='
  prevGeneId <-  c(annRes[1,'gene_id'],annRes[,'gene_id'])
  prevGeneId <- prevGeneId[-length(prevGeneId)]
  # all indices where the current gene id is different from the previous one
  diffGeneInd <- which(annRes[,'gene_id']!=prevGeneId)
  # all indices where sliding window should be split
  splitInd <- sort(posInd)
  prevInd <- 1
  for(p in splitInd){
    mergeInd <- c(prevInd:(p-1))
    if(length(mergeInd)==1){
      annRes[mergeInd[1],'pLocal'] <- annRes[mergeInd[1],'pvalue']
    }else{
      tmpRes <- annRes[mergeInd,]
      rownames(tmpRes) <- NULL
      for(m in c(1:length(mergeInd))){
        upstreamInd <- which(tmpRes[,'begin']< tmpRes[m,'begin'] & tmpRes[,'end']> tmpRes[m,'begin']) # all indices where current begin is in between start and stop
        downstremInd <- which(tmpRes[,'begin']< tmpRes[m,'end'] & tmpRes[,'end']> tmpRes[m,'end'])# all indices where current end is in between start and stop
        annRes[ mergeInd[m], 'pLocal'] <- pvalue_dependence_correction(tmpRes[m,'pvalue'],tmpRes[union(upstreamInd,downstremInd),'pvalue' ])
      }
      annRes[mergeInd,'regionStartId'] <- annRes[mergeInd[1],'unique_id']
    }
    prevInd <- p
  }
  if(prevInd<=nrow(annRes)){
    mergeInd <- c(prevInd:nrow(annRes))
    if(length(mergeInd)==1){
      annRes[mergeInd[1],'pLocal'] <- annRes[mergeInd[1],'pvalue']
    }else{
      tmpRes <- annRes[mergeInd,]
      rownames(tmpRes) <- NULL
      for(m in c(1:length(mergeInd))){
        upstreamInd <- which(tmpRes[,'begin']< tmpRes[m,'begin'] & tmpRes[,'end']> tmpRes[m,'begin']) # all indices where current begin position is in between start and stop
        downstremInd <- which(tmpRes[,'begin']< tmpRes[m,'end'] & tmpRes[,'end']> tmpRes[m,'end'])# all indices where current end is in between start and stop
        annRes[ mergeInd[m], 'pLocal'] <- pvalue_dependence_correction(tmpRes[m,'pvalue'],tmpRes[union(upstreamInd,downstremInd),'pvalue' ])
      }
      annRes[mergeInd,'regionStartId'] <- annRes[mergeInd[1],'unique_id']
    }
  }
  return(annRes)
}

# @TODO: fix repeating code!

#' @title  bonferroni correction
#' @description Helper function merge overlapping windows, correct p-values for these overlapping windows using bonferroni FWER
#' @param annRes annotated results (data.frame)
#' @param minDist minimum distance between non overlapping windows,
#' default: 0 to account for bed formatted window regions
#' @keywords internal
mergeWindowsHelper_bonferroni <- function(annRes,minDist=0){
  annRes <- annRes[with(annRes,order(begin,end)),]
  annRes <- data.frame(annRes,regionStartId=annRes[,'unique_id'],pBonferroni=rep(1,nrow(annRes)),stringsAsFactors = FALSE)
  endW <- c(annRes[1,'begin']+1,annRes[,'end']) # end of all previous windows, except for first entry
  endW <- endW[-length(endW)] # remove last
  beDiff <- annRes[,'begin'] - endW
  # all indices where the current window do not overlap the previous one
  posInd <- which(beDiff>=minDist) # The operator MUST BE '>='
  prevGeneId <-  c(annRes[1,'gene_id'],annRes[,'gene_id'])
  prevGeneId <- prevGeneId[-length(prevGeneId)]
  # all indices where the current gene id is different from the previous one
  diffGeneInd <- which(annRes[,'gene_id']!=prevGeneId)
  # all indices where sliding window should be split
  splitInd <- sort(posInd)
  prevInd <- 1
  for(p in splitInd){
    mergeInd <- c(prevInd:(p-1))
    if(length(mergeInd)==1){
      annRes[mergeInd[1],'pBonferroni'] <- annRes[mergeInd[1],'pvalue']
    }else{
      tmpRes <- annRes[mergeInd,]
      rownames(tmpRes) <- NULL
      for(m in c(1:length(mergeInd))){
        upstreamInd <- which(tmpRes[,'begin']< tmpRes[m,'begin'] & tmpRes[,'end']> tmpRes[m,'begin']) # all indices where current begin is in between start and stop
        downstremInd <- which(tmpRes[,'begin']< tmpRes[m,'end'] & tmpRes[,'end']> tmpRes[m,'end'])# all indices where current end is in between start and stop
        annRes[ mergeInd[m], 'pBonferroni'] <- p.adjust(p=tmpRes[m,'pvalue'],method='bonferroni',n=length(union(upstreamInd,downstremInd))+1)
      }
      annRes[mergeInd,'regionStartId'] <- annRes[mergeInd[1],'unique_id']
    }
    prevInd <- p
  }
  if(prevInd<=nrow(annRes)){
    mergeInd <- c(prevInd:nrow(annRes))
    if(length(mergeInd)==1){
      annRes[mergeInd[1],'pBonferroni'] <- annRes[mergeInd[1],'pvalue']
    }else{
      tmpRes <- annRes[mergeInd,]
      rownames(tmpRes) <- NULL
      for(m in c(1:length(mergeInd))){
        upstreamInd <- which(tmpRes[,'begin']< tmpRes[m,'begin'] & tmpRes[,'end']> tmpRes[m,'begin']) # all indices where current begin position is in between start and stop
        downstremInd <- which(tmpRes[,'begin']< tmpRes[m,'end'] & tmpRes[,'end']> tmpRes[m,'end'])# all indices where current end is in between start and stop
        annRes[ mergeInd[m], 'pBonferroni'] <- p.adjust(p=tmpRes[m,'pvalue'],method='bonferroni',n=length(union(upstreamInd,downstremInd))+1)
      }
      annRes[mergeInd,'regionStartId'] <- annRes[mergeInd[1],'unique_id']
    }
  }
  return(annRes)
}

#' @title p-value correction
#' @description pvalue correction using local simes like approach
#' @param p_to_correct pvalue to correct (from the current window)
#' @param dependent_p_values pvalues from neighboring windows
#' @keywords internal
pvalue_dependence_correction <- function(p_to_correct, dependent_p_values){
  # author: Tom
  # p.to.correct is a numeric
  # dependent.p.values is a vector of p.values
  if(p_to_correct > 1| p_to_correct < 0 | any(dependent_p_values > 1) | any(dependent_p_values < 0)) {
    stop('P-values must fall in the range 0 <= p <= 1')
  }
  nr_of_p_values <- length(dependent_p_values) + length(p_to_correct)
  inv_ranks <- ((nr_of_p_values + 1) - rank(c(p_to_correct, dependent_p_values)))
  inv_rank_of_p_to_correct <- inv_ranks[1:length(p_to_correct)]
  return (min(1,(p_to_correct * nr_of_p_values / inv_rank_of_p_to_correct))) # to correct for p-values going above 1
}

#' @title merge windows
#' @description merge annotated sliding windows
#' @param annRes annotated DESEq2 table
#' @param minDist minimum distance between the sliding windows
#' default 0 to account for bed formatted window regions
#' @param padjWindow p-value adjustment method for window, can be one of 'local' or 'bonferroni'
#' @param padjMethod p-value adjustment method (multiple testing correction) for the whole table, can be any input for base function p.adjust
#' @param ncores number of cores to use
mergeWindows <- function(annRes,minDist=0,padjWindow='bonferroni',padjMethod='BH',ncores=5){
  annRes <- na.omit(annRes)
  padjWindow <- match.arg(padjWindow,choices=c('local','bonferroni'),several.ok=FALSE)
  neededCols <- c('chromosome','unique_id','begin','end','strand','baseMean','log2FoldChange','lfcSE','stat','pvalue','gene_id','gene_name',
                  'gene_type','gene_region','Nr_of_region','Total_nr_of_region','window_number')
  missingCols <- setdiff(neededCols,colnames(annRes))
  if(length(missingCols)>0){
    stop('Input data.frame is missing required columns, needed columns:
         chromosome: chromosome name
         unique_id: unique id of the window
         begin: window start co-ordinate
         end: window end co-ordinate
         strand: strand
         baseMean: baseMean column from DESeq2 results
         log2FoldChange: log2FoldChange column from DESeq2 results
         lfcSE: lfcSE column from DESeq2 results
         stat: stat column from DESeq2 results
         pvalue: pvalue column from DESeq2 results
         gene_id: gene id
         gene_name: gene name
         gene_type: gene type annotation
         gene_region: gene region
         Nr_of_region: number of the current region
         Total_nr_of_region: total number of regions
         window_number: window number
         Missing columns:
         ',paste(missingCols,collapse=", "),'')
  }
  register(BatchtoolsParam(workers = ncores), default = TRUE)
  rownames(annRes) <- NULL
  geneIds <- unique(annRes[,'gene_id'])
  geneList <- vector('list',length(geneIds))
  names(geneList) <- geneIds
  for(geneId in geneIds){
    geneList[[geneId]] <- annRes[ annRes['gene_id']==geneId, ]
  }
  mergeDat <- NULL
  if (padjWindow=='local'){
        mergeDat <- do.call(rbind, bplapply(geneList, mergeWindowsHelper_local))
        mergeDat$pLocal.adj <- p.adjust(mergeDat[,'pLocal'], method = padjMethod)
    } else if (padjWindow == 'bonferroni') {
        mergeDat <- do.call(rbind, bplapply(geneList, mergeWindowsHelper_bonferroni))
        mergeDat$pBonferroni.adj <- p.adjust(mergeDat[,'pBonferroni'], method = padjMethod)
    }
  rownames(mergeDat) <- NULL
  message('\n')
  return(mergeDat)
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
#' @returns data.frame
createRegions <- function(mergeDat,padjCol='pBonferroni.padj',padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1,ncores=5){
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
  register(BatchtoolsParam(workers = ncores), default = TRUE)
  rownames(mergeDat) <- NULL
  geneIds <- unique(mergeDat[,'gene_id'])
  geneList <- vector('list',length(geneIds))
  names(geneList) <- geneIds
  for(geneId in geneIds){
    geneList[[geneId]] <- mergeDat[ mergeDat['gene_id']==geneId, ]
  }
  regionDat <- do.call(rbind, bplapply(geneList, createRegionsHelper))
  return(regionDat)
}

#' @title: helper function, create region
#' @description: merge significant neighboring windows into regions
#' @param mergeDat a gene specific result data.frame
#' @param padjCol name of the adjusted pvalue column (default: pBonferroni.adj)
#' @param padjThresh threshold for p-adjusted value (default: 0.05)
#' @param log2FoldChangeCol name of the log2foldchange column (default: log2FoldChange)
#' @param log2FoldChangeThresh threshold for log2foldchange value (default:1)
#' @param ncores number of cores to uses
#' @keywords internal
createRegionsHelper <- function(mergeDat,padjCol='pBonferroni.padj',padjThresh=0.05,log2FoldChangeCol='log2FoldChange',log2FoldChangeThresh=1){
  mergeDat <- mergeDat[with(mergeDat,order(begin,end)),]
  rownames(mergeDat) <- NULL
  mergeDat[,'regionStartId'] <- as.character(mergeDat[,'regionStartId'])
  mergeDat[,'unique_id'] <- as.character(mergeDat[,'unique_id'])
  prevRegionId <-  c(mergeDat[1,'regionStartId'],mergeDat[,'regionStartId'])
  prevRegionId <- prevRegionId[-length(prevRegionId)]
  # all indices where the current region start id is different from the previous one
  diffRegionInd <- which(mergeDat[,'regionStartId']!=prevRegionId)
  # all indices where padj value is > padjThresh OR log2FoldChange is < log2foldchangeThresh
  nonSigInd <- which(mergeDat[,padjCol]>padjThresh | mergeDat[,log2FoldChangeCol]<log2FoldChangeThresh)
  mergeRes <- data.frame(regionStartId=mergeDat[,'unique_id'],region_begin=mergeDat[,'begin'],region_end=mergeDat[,'end'],padjMin=mergeDat[,padjCol],padjMax=mergeDat[,padjCol],
            log2FoldChangeMin=mergeDat[,log2FoldChangeCol],log2FoldChangeMax=mergeDat[,log2FoldChangeCol],stringsAsFactors = FALSE)
  # if there are any significant windows, get all significant windows
  sigInd <- setdiff(c(1:nrow(mergeDat)),sort(union(diffRegionInd,nonSigInd)))
  if(length(sigInd)>0){
    sigList <- vector('list',length(sigInd))
    indStart <- 1
    for(i in c(1:length(sigInd))){
      if(i==1){
        sigList[[indStart]] <- c(sigInd[i])
      }else if(sigInd[i]-sigInd[i-1]==1){
        sigList[[indStart]] <- c(sigList[[indStart]],sigInd[i])
      }else{
        indStart <- indStart+1
        sigList[[indStart]] <- c(sigInd[i])
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
  return(cbind(mergeDat[,-which(colnames(mergeDat)=='regionStartId')],mergeRes))
}
#' @TODO: write a wrapper for DESeq2 wald test
