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
  #prevGeneId <-  c(annRes[1,'gene_id'],annRes[,'gene_id'])
  #prevGeneId <- prevGeneId[-length(prevGeneId)]
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
