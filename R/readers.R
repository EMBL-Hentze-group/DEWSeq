# Function(s) to read in large files
# Author: Sudeep Sahadevan, sudeep.sahadevan@embl.de

#'
#' @keywords  internal
#' @importFrom BiocGenerics sort
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
#' @importFrom data.table fread getDTthreads setDTthreads
#' @importFrom utils read.table
#'
#' @title  read annotation data
#' @description read annotation data for windows
#' This is an unexported lowlevel function to read annotations from a file.
#' The file MUST be tab separated and MUST have the following columns:\cr
#'  chromosome: chromosome name \cr
#'  unique_id: unique id of the window \cr
#'  begin: window start co-ordinate \cr
#'  end: window end co-ordinate \cr
#'  strand: strand \cr
#'  gene_id: gene id \cr
#'  gene_name: gene name \cr
#'  gene_type: gene type annotation \cr
#'  gene_region: gene region \cr
#'  Nr_of_region: number of the current region \cr
#'  Total_nr_of_region: total number of regions \cr
#'  window_number: window number \cr
#' @param fname file name/path
#' @param uniqIds filter stable and keep annotation for these unique ids
#' @param asGRange boolean, whether to reaturn a GRange object or a data.frame
#' object
#' @param checkWindowNumber check window number
#' @param start0based TRUE (default) or FALSE. If TRUE, then the start
#' @param threads number of threads for fread (default: 10)
#' positions are considered to be 0-based
#'
#'
#' @return data.frame or GRanges object
.readAnnotation <- function(fname,uniqIds=NULL,asGRange=TRUE,checkWindowNumber=TRUE,
                            start0based=TRUE,threads=10){
  neededCols <- c('chromosome','unique_id','begin','end','strand',
                  'gene_id','gene_name','gene_type','gene_region',
                  'Nr_of_region', 'Total_nr_of_region')
  if(checkWindowNumber){
    neededCols <- c(neededCols,'window_number')
  }
  # check number of available cores for fread
  allCores <- getDTthreads()
  setDTthreads(ifelse(threads>=allCores,allCores-1,threads))
  annTable <- fread(fname,sep="\t",stringsAsFactors = FALSE,header=TRUE)
  missingCols <- setdiff(neededCols,colnames(annTable))
  if(length(missingCols)>0 & checkWindowNumber){
    stop('Input annotation file is missing required columns, needed columns:
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
         Missing columns:
         ',paste(missingCols,collapse=", "),'')
  }else if(length(missingCols)>0 &!checkWindowNumber){
    stop('Input annotation file is missing required columns, needed columns:
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
         Missing columns:
         ',paste(missingCols,collapse=", "),'')
  }
  annTable <- as.data.frame(annTable)
  if(length(annTable$unique_id)==length(unique(annTable$unique_id))){
    rownames(annTable) <- annTable$unique_id
  }
  if(is.null(uniqIds)){
    if(asGRange){
      gr <- makeGRangesFromDataFrame(annTable,seqnames.field='chromosome',start.field='begin',end.field='end',
                                     strand.field='strand', ignore.strand=FALSE,
                                     keep.extra.columns=TRUE,
                                     starts.in.df.are.0based=start0based)
      gr <- sortSeqlevels(gr)
      gr <- sort(gr)
      rm(annTable)
      gc()
      return(gr)
    }else{
      return(annTable)
    }
  }else{
    commonIds <- intersect(as.character(uniqIds),rownames(annTable))
    if(length(commonIds)==0){
      stop('There are no common unique ids between the input file and uniqueIds. Please check your data sets!')
    }
    if(asGRange){
      gr <- makeGRangesFromDataFrame(annTable[commonIds,], seqnames.field='chromosome',
                                     start.field='begin',end.field='end',
                                     strand.field='strand', ignore.strand=FALSE,
                                     keep.extra.columns=TRUE,
                                     starts.in.df.are.0based=start0based)
      gr <- sortSeqlevels(gr)
      gr <- sort(gr)
      rm(annTable)
      gc()
      return(gr)
    }else{
      return(annTable[commonIds,])
    }
  }
}
