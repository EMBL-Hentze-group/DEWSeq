#' @title  read annotation data
#' @description read annotation data for windows
#' The file MUST be tab separated and MUST have the following columns:
#' chromosome: chromosome name
#' unique_id: unique id of the window
#' begin: window start co-ordinate
#' end: window end co-ordinate
#' strand: strand
#' gene_id: gene id
#' gene_name: gene name
#' gene_type: gene type annotation
#' gene_region: gene region
#' Nr_of_region: number of the current region
#' Total_nr_of_region: total number of regions
#' window_number: window number
#' @param fname file name/path
#' @param uniqIds filter stable and keep annotation for these unique ids
.readAnnotation <- function(fname,uniqIds=NULL){
  neededCols <- c('chromosome','unique_id','begin','end','strand','gene_id','gene_name','gene_type','gene_region','Nr_of_region',
                 'Total_nr_of_region','window_number')
  gzlen <- grep(pattern = '\\.gz',ignore.case = TRUE,x=fname)
  platform <- Sys.info()[['sysname']]
  if(gzlen>0 && platform=='Windows'){
    annTable <- read.table(gzfile(fname),sep="\t",stringsAsFactors=FALSE,header=TRUE)
  }else if(gzlen>0){ # assuming that zcat binary is nstalled in Linux and Mac distributions
    annTable <- data.table::fread(input=paste('zcat',fname),sep="\t",stringsAsFactors = FALSE,header=TRUE)
  }else if (gzlen==0){
    annTable <- data.table::fread(fname,sep="\t",stringsAsFactors = FALSE,header=TRUE)
  }
  missingCols <- setdiff(neededCols,colnames(annTable))
  if(length(missingCols)>0){
    stop('Input annotatino file is missing required columns, needed columns:
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
  }
  annTable <- as.data.frame(annTable)
  rownames(annTable) <- annTable$unique_id
  if(is.null(uniqIds)){
    gr <- GenomicRanges::makeGRangesFromDataFrame(annTable,seqnames.field='chromosome',start.field='begin',end.field='end',strand.field='strand',
                                                  ignore.strand=FALSE,keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)
    gr <- GenomeInfoDb::sortSeqlevels(gr)
    gr <- BiocGenerics::sort(gr)
    return(gr)
  }else{
    commonIds <- intersect(as.character(uniqIds),rownames(annTable))
    if(length(commonIds)==0){
      stop('There are no common unique ids between the input file and uniqueIds. Please check your data sets!')
    }
    gr <- GenomicRanges::makeGRangesFromDataFrame(annTable[commonIds,], seqnames.field='chromosome',start.field='begin',end.field='end',
                                                   strand.field='strand', ignore.strand=FALSE,keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)
    gr <- GenomeInfoDb::sortSeqlevels(gr)
    gr <- BiocGenerics::sort(gr)
    return(gr)
  }

}
