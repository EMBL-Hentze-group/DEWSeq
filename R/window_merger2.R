#' @export
#' @title merge windows
#' @description merge annotated sliding windows
#' @param annRes annotated DESEq2 table
#' @param minOverlap minimum overlap between sliding windows
mergeWindows <- function(annRes,minOverlap=5,padjMethod='BH'){
  annRes <- na.omit(annRes)
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
  rownames(annRes) <- NULL
  chrPos <- annRes[,c('chromosome','begin','end','strand')]
  colnames(chrPos) <- c('chr','start','end','strand','id')
  chrRange <- GenomicRanges::makeGRangesFromDataFrame(chrPos,keep.extra.columns = TRUE)
  chrOverlap <- GenomicRanges::findOverlaps(chrRange,ignore.strand=FALSE,minoverlap = minOverlap,drop.redundant=FALSE,drop.self=FALSE)
  exactOverlap <- GenomicRanges::findOverlaps(chrRange,ignore.strand=FALSE,drop.redundant=FALSE,drop.self=FALSE,type='equal')
  nOverlaps <- as.data.frame(table(queryHits(chrOverlap)))[,2]
  eOverlaps <- as.data.frame(table(queryHits(exactOverlap)))[,2]
  nOverlaps <- nOverlaps - (eOverlaps-1)
  rm(chrPos)
}
