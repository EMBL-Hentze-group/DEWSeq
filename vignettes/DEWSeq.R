## ----setup, echo=FALSE, results="hide"-------------------------------------
knitr::opts_chunk$set(tidy    = FALSE,
                      cache   = FALSE,
                      dev     = "png",
                      message = FALSE,
                      error   = FALSE,
                      warning = TRUE)
suppressPackageStartupMessages(require(BiocStyle))

## ----citaction, eval = FALSE, echo = FALSE---------------------------------
#  #**Note:** if you use DEWSeq in published research, please cite:
#  #
#  #> Authors. (Year)
#  #> Title
#  #> *Genome Biology*, **15**:550.
#  #> [10.1186/s13059-014-0550-8](http://dx.doi.org/10.1186/s13059-014-0550-8)

## ---- fig.cap = "Crosslink site trunction at reverse transcription", echo = FALSE----
knitr::include_graphics("truncation.png")

## ---- fig.cap = "Truncation sites (referred to as crosslink sites) are the neighboring position of the aligned read", fig.small = TRUE, echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("truncationsite.png")

## ---- fig.cap = "Different binding modes of RNA-binding proteins", fig.wide = TRUE, echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("binding_modes.png")

## ---- fig.cap = "Chance of crosslinking", fig.small = TRUE, echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("crosslinking_chance.png")

## ---- fig.cap = "Different RNase concentrations will result in different fragment sizes.", fig.small=TRUE, echo = FALSE----
knitr::include_graphics("digestionpatterns.png")

## ---- fig.cap = "Read-throughs", fig.small = T, echo = FALSE---------------
knitr::include_graphics("readthrough.png")

## ---- fig.cap = "Early truncation events", echo = FALSE--------------------
knitr::include_graphics("earlytruncation.png")

## ---- fig.cap = "Sliding window approach with single-nucleotide position count data", echo = FALSE----
knitr::include_graphics("dewseqoverview.png")

## ---- fig.cap = "Combining significant windows", echo = FALSE, fig.small = TRUE----
# fig.width=200, out.width="200px"
knitr::include_graphics("overlapping_sliding_windows.png")

## ---- fig.cap = "htseq-clip workflow", echo = FALSE------------------------
knitr::include_graphics("htseq-clip.png")

## ----install, eval = FALSE-------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("DEWSeq")

## ----load library, eval = TRUE---------------------------------------------
require(DEWSeq)
require(IHW)
require(tidyverse)

## ----filenames, eval = TRUE------------------------------------------------
countFile <- file.path(system.file('extdata',package='DEWSeq'),'SLBP_K562_w50s20_counts.txt.gz')
annotationFile <- file.path(system.file('extdata',package='DEWSeq'),'SLBP_K562_w50s20_annotation.txt.gz')

## ----loading tidyverse, eval = TRUE, results='hide'------------------------
library(tidyverse)
library(data.table)

## ----read in count matrix, eval = TRUE-------------------------------------
countData <- fread(countFile, sep = "\t")

## ----read count matrix tidyr, eval = TRUE----------------------------------
countData <- read_tsv(countFile)

## ----read in annotation, eval = TRUE---------------------------------------
annotationData <- fread(annotationFile, sep = "\t")

## ----read annotation tidyr, eval = FALSE-----------------------------------
#  annotationData <- read_tsv(annotationFile)

## ----datasummary-----------------------------------------------------------
head(countData)
dim(countData)
head(annotationData)
dim(annotationData)

## ----create colData--------------------------------------------------------
colData <- data.frame(
  row.names = colnames(countData)[-1], # since the first column is unique_id
  type      = factor(
                c(rep("IP", 2),    ##  change this accordingly
                  rep("SMI", 1)),  ##
                levels = c("IP", "SMI"))
)

## ----example import--------------------------------------------------------
ddw <- DESeqDataSetFromSlidingWindows(countData  = countData,
                                      colData    = colData,
                                      annotObj   = annotationData,
                                      tidy       = TRUE,
                                      design     = ~type)
ddw

## ----row filtering, eval = TRUE--------------------------------------------
keep <- rowSums(counts(ddw)) >= 10
ddw <- ddw[keep,]

## ----estimate size factors, eval = TRUE------------------------------------
ddw <- estimateSizeFactors(ddw)
sizeFactors(ddw)

## ----filter for mRNAs, eval = TRUE-----------------------------------------
ddw_mRNAs <- ddw[ rowData(ddw)[,"gene_type"] == "protein_coding", ]
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)
sizeFactors(ddw)

## ----tmp significant windows-----------------------------------------------
ddw_tmp <- ddw
ddw_tmp <- estimateDispersions(ddw_tmp, fitType = "local", quiet = TRUE)
ddw_tmp <- nbinomWaldTest(ddw_tmp)

tmp_significant_windows <- 
                results(ddw_tmp,
                    contrast = c("type", "IP", "SMI"),
                    tidy = TRUE,
                    filterFun = ihw) %>% 
                dplyr::filter(padj < 0.05) %>% 
                .[["row"]]
rm("ddw_tmp")

## --------------------------------------------------------------------------
ddw_mRNAs <- ddw_mRNAs[ !rownames(ddw_mRNAs) %in% tmp_significant_windows, ]
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)

rm( list = c("tmp_significant_windows", "ddw_mRNAs"))

sizeFactors(ddw)

## ----estimate dispersions and wald test, eval = TRUE-----------------------
ddw <- estimateDispersions(ddw, fitType = "local", quiet = TRUE)
ddw <- nbinomWaldTest(ddw)

## ---- fig.cap = "Dispersion Estimates", fig.wide = TRUE--------------------
plotDispEsts(ddw)

## ----DEWSeq results, eval = TRUE-------------------------------------------
resultWindows <- resultsDEWSeq(ddw,
                              contrast = c("type", "IP", "SMI"),
                              tidy = TRUE) %>% as_tibble

resultWindows

## ----IHW multiple hypothesis correction, eval = TRUE, warning = FALSE------
resultWindows[,"p_adj_IHW"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, 
                     data = resultWindows,
                     alpha = 0.05,
                     nfolds = 10))

## ----windows tables, eval = TRUE-------------------------------------------
resultWindows <- resultWindows %>% 
                      mutate(significant = resultWindows$p_adj_IHW < 0.01)

sum(resultWindows$significant)

## ----how genes form windows------------------------------------------------
resultWindows %>%
   filter(significant) %>% 
   arrange(desc(log2FoldChange)) %>% 
   .[["gene_name"]] %>% 
   unique %>% 
   head(20)

## ----extractRegions, eval = TRUE, results = "hide"-------------------------
resultRegions <- extractRegions(windowRes  = resultWindows,
                                padjCol    = "p_adj_IHW",
                                padjThresh = 0.01, 
                                log2FoldChangeThresh = 0.5) %>% as_tibble

## ----resultRegions output, eval = TRUE-------------------------------------
resultRegions

## ----toBED output, eval = FALSE--------------------------------------------
#  toBED(windowRes = resultWindows,
#        regionRes = resultRegions,
#        fileName  = "enrichedWindowsRegions.bed")

## --------------------------------------------------------------------------
sessionInfo()

