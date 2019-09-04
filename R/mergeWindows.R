#' @export
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
  BiocParallel::register(BatchtoolsParam(workers = ncores), default = TRUE)
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
