#' @export
#' @description
#' for the given widow description file and annotation file,
#' find overlaps between the two in terms of chromosomal locations,
#' return the number of regions (3'UTR, 5'UTR, exon, intron, CDS ...) as a
#' as a count table
countGeneRegion <- function(regionRes,annotationFile,begin0based=TRUE,minOverlap=0.5){
  requiredCols <- c('chromosome','regionStartId','region_begin','region_end','strand','gene_id','gene_name')
  missingCols <- setdiff(requiredCols,colnames(regionRes))
  if(length(missingCols)>0){
    stop('Input data.frame is missing required columns, needed columns:
        chromosome: chromosome name
        regionStartId: unique id of the left most window of the region
        begin: region start co-ordinate
        end: region end co-ordinate
        strand: strand
        gene_id: gene id\n',
         'Missing columns: ',paste(missingCols,collapse=", "),'')
  }
  regDat <- GenomicRanges::makeGRangesFromDataFrame(regionRes[,requiredCols],seqnames.field = 'gene_id',start.field = 'begin',end.field = 'end',strand.field = 'strand',
                                                    ignore.strand=FALSE,starts.in.df.are.0based=begin0based,keep.extra.columns = TRUE)
  regDat <-  GenomeInfoDb::sortSeqlevels(regDat)
  regDat <- BiocGenerics::sort(regDat)
  annDat <- .readAnnotation(fname=annotationFile,uniqIds =NULL, asGRange=FALSE,checkWindowNumber=FALSE)
  annDat <- GenomicRanges::makeGRangesFromDataFrame(annDat,seqnames.field = 'gene_id',start.field = 'begin',end.field = 'end',strand.field = 'strand',
                                                    ignore.strand=FALSE,starts.in.df.are.0based=begin0based,keep.extra.columns = TRUE)
  annDat <-  GenomeInfoDb::sortSeqlevels(annDat)
  annDat <- BiocGenerics::sort(annDat)
  resAnnOv <- GenomicRanges::findOverlaps(regDat,annDat)
  ovRegion <- GenomicRanges::pintersect(regDat[queryHits(resAnnOv)],annDat[subjectHits(resAnnOv)],drop.nohit.ranges=TRUE)
  selectedOV <- BiocGenerics::width(ovRegion)/BiocGenerics::width(regDat[queryHits(resAnnOv)])>=minOverlap
  if(!any(selectedOV)){
    stop('Cannot find overlaps between regionRes and annotationFile. Please lower minOverlap threshold and try again')
  }
  regionStartId <- mcols(regDat[queryHits(resAnnOv)[selectedOV]])[,'regionStartId']
  geneId <- gsub('\\:.*$','',regionStartId)
  outDat <- data.frame(mcols(annDat[subjectHits(resAnnOv)[selectedOV]])[,c('unique_id','gene_type','gene_name','gene_region')],regiontartId=regionStartId,gene_id=geneId)
  return(outDat)
}
