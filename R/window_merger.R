
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