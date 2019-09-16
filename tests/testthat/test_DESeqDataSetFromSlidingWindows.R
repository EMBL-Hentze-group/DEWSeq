context("DESeqDataSetFromSlidingWindows")
test_that("DESeqDataSetFromSlidingWindows throws proper errors during object construction", {
  set.seed(1)
  library(data.table)
  phenoDat <- DataFrame(condition=factor(c('test','test','control','control')),row.names=c('S1','S2','S3','S4'))
  testCols <- c('chromosome','unique_id','begin','end','strand','gene_id',
                'gene_name','gene_type','gene_region','Nr_of_region',
                'Total_nr_of_region','window_number')
  annotDf <- data.frame(chromosome = paste('chr',c(1,1,3,3,5,10,15,20,'X','Y'),sep = ''),
                        unique_id = paste('uid',c(1:10),sep='_'), begin = c(100,124,212,215,986,763,3221,2311,1212,3322),
                        end=c(150,174,262,265,1036,813,3271,2361,1262,3372), strand=c('-', '+', '+', '+', '-', '-', '+', '+', '-', '-'),
                        gene_id=paste('gene',c(1,1,3,3,5,10,15,20,'X','Y'),sep = ''),gene_name=paste('gene',c(1,1,3,3,5,10,15,20,'X','Y'),sep = ''),
                        gene_type=c(rep('protein_coding',3),'rRNA',' miRNA',rep('lncRNA',2),'pseudogene','snoRNA','protein_coding'),
                        gene_region=c(rep('exon',5),rep('intron',5)),Nr_of_region=c(1,2,4,5,11,10,5,7,1,33),
                        Total_nr_of_region=c(6,7,9,10,16,15,10,12,6,38),window_number=c(287,614,145,329,487,855,851,630,498,858))
  expect_error(DESeqDataSetFromSlidingWindows(matrix(0,10,4),phenoDat,annotDf),'rownames of countData cannot be empty')
  expect_error(DESeqDataSetFromSlidingWindows(data.frame(matrix(0,10,4)),phenoDat,annotObj=matrix(0,10,4)),'annotObj MUST be a data.frame or character')
  errMsg <- 'Input annotation data.frame is missing required columns, needed columns:
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
         window_number'
  expect_error(DESeqDataSetFromSlidingWindows(data.frame(matrix(0,10,4)),phenoDat,annotDf[,c(1:11)]),errMsg)
  expect_error(DESeqDataSetFromSlidingWindows(data.frame(matrix(0,10,4)),phenoDat,annotDf),
  'There are no common unique ids between the countData and annotObj. Please check your data sets!')
  sampleData <- data.frame(matrix(rnbinom(48,10,0.8),12,4,dimnames=list(c(paste('uid',c(1:12),sep='_')),c('S1','S2','S3','S4'))))
  warnMsg <- 'Cannot find chromosomal positions for all entries in countData.
            countData rows with missing annotation will be removed !'
  expect_warning(DESeqDataSetFromSlidingWindows(sampleData,phenoDat,annotDf,design = ~condition),warnMsg)
  expect_warning(DESeqDataSetFromSlidingWindows(data.table(data.frame(uid=rownames(sampleData),sampleData))[c(1:10),],phenoDat,annotDf,design = ~condition),
                 'countData is a data.table object, converting it to data.frame. First column will be used as rownames')
})
