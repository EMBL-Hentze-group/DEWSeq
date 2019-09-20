## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(tidy    = FALSE,
                      cache   = FALSE,
                      dev     = "png",
                      message = FALSE,
                      error   = FALSE,
                      warning = TRUE)
suppressPackageStartupMessages(require(BiocStyle))

## ----citaction, eval = FALSE, echo = FALSE-------------------------------
#  #**Note:** if you use DEWSeq in published research, please cite:
#  #
#  #> Authors. (Year)
#  #> Title
#  #> *Genome Biology*, **15**:550.
#  #> [10.1186/s13059-014-0550-8](http://dx.doi.org/10.1186/s13059-014-0550-8)

## ---- fig.cap = "Crosslink site trunction at reverse transcription", echo = FALSE----
knitr::include_graphics("truncation.png")

## ---- fig.cap = "Truncation sites (often referred to as crosslink sites) are the neighboring position of the aligned read.", fig.small = TRUE, echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("truncationsite.png")

## ----libs, message=FALSE-------------------------------------------------
library("IHW")
library("tidyverse")

## ---- fig.cap = "Different binding modes of RNA-binding proteins", fig.wide = TRUE, echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("binding_modes.png")

## ---- fig.cap = "Chance of crosslinking", fig.small = TRUE, echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("crosslinking_chance.png")

## ---- fig.cap = "Different RNase concentrations will result in different fragment sizes.", fig.small=TRUE, echo = FALSE----
knitr::include_graphics("digestionpatterns.png")

## ---- fig.cap = "Read-throughs", fig.small = T, echo = FALSE-------------
knitr::include_graphics("readthrough.png")

## ----early truncation, fig.cap = "Early truncations.", echo = FALSE------
knitr::include_graphics("earlytruncation.png")

## ----deseq overview image, fig.cap = "Sliding window approach with single-nucleotide position count data", echo = FALSE----
knitr::include_graphics("dewseqoverview.png")

## ---- fig.cap = "Sliding window approach with single-nucleotide position count data", echo = FALSE, fig.small = TRUE----
# fig.width=200, out.width="200px"
knitr::include_graphics("overlapping_sliding_windows.png")

## ---- echo = FALSE, eval = FALSE-----------------------------------------
#  # Analysis strategies for CLIP data

## ----load library, eval = TRUE-------------------------------------------
require(DEWSeq)

## ----load test data, eval = TRUE-----------------------------------------
data(SLBP_K562_w50s20)

## ----head of test data, eval = TRUE--------------------------------------
ddw <- SLBP_K562_w50s20
rm(SLBP_K562_w50s20)
ddw

## ----loading tidyverse, eval = TRUE, results='hide'----------------------
library(tidyverse)
library(data.table)

## ----read in count matrix, eval = FALSE----------------------------------
#  countData <- fread("path/swcounts/count_matrix.txt.gz", sep = "\t")

## ----read count matrix tidyr, eval = FALSE-------------------------------
#  countData <- read_tsv("path/swcounts/count_matrix.txt.gz")

## ----read in annotation, eval = FALSE------------------------------------
#  annotationData <- fread("path/annotation/annotation.txt.gz", sep = "\t")

## ----read annotation tidyr, eval = FALSE---------------------------------
#  annotationData <- read_tsv("path/annotation/annotation.txt.gz")

## ----create colData, eval = FALSE----------------------------------------
#  colData <- data.frame(
#    row.names = colnames(countData),
#    type      = factor(
#                  c(rep("IP", 3),    ##  change this accordingly
#                    rep("SMI", 3)),  ##
#                  levels = c("IP", "SMI"))
#  )

## ----example import, eval = FALSE----------------------------------------
#  ddw <- DESeqDataSetFromSlidingWindows(countData  = countData,
#                                        colData    = colData,
#                                        annotation = annotationData,
#                                        tidy       = TRUE,
#                                        design     = ~type)

## ----row filtering, eval = TRUE------------------------------------------
keep <- rowSums(counts(ddw)) >= 10
ddw <- ddw[keep,]

## ----estimate size factors, eval = TRUE----------------------------------
ddw <- estimateSizeFactors(ddw)
sizeFactors(ddw)

## ----filter for mRNAs, eval = TRUE---------------------------------------
ddw_mRNAs <- ddw[ rowData(ddw)[,"gene_type"] == "protein_coding", ]
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)
sizeFactors(ddw)

## ----tmp significant windows---------------------------------------------
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

## ------------------------------------------------------------------------
ddw_mRNAs <- ddw_mRNAs[ !rownames(ddw_mRNAs) %in% tmp_significant_windows, ]
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)

rm( list = c("tmp_significant_windows", "ddw_mRNAs"))

sizeFactors(ddw)

## ----estimate dispersions and wald test, eval = TRUE---------------------
ddw <- estimateDispersions(ddw, fitType = "local", quiet = TRUE)
ddw <- nbinomWaldTest(ddw)

## ---- fig.cap = "Dispersion Estimates", fig.wide = TRUE------------------
plotDispEsts(ddw)

## ----DEWSeq results, eval = TRUE-----------------------------------------
resultWindows <- resultsDEWSeq(ddw,
                              contrast = c("type", "IP", "SMI"),
                              tidy = TRUE) %>% as_tibble

resultWindows

## ----IHW multiple hypothesis correction, eval = TRUE, warning = FALSE----
resultWindows[,"p_adj_IHW"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, 
                     data = resultWindows,
                     alpha = 0.05,
                     nfolds = 10))

## ----windows tables, eval = TRUE-----------------------------------------
resultWindows <- resultWindows %>% 
                      mutate(significant = resultWindows$p_adj_IHW < 0.01)

sum(resultWindows$significant)

## ----how genes form windows----------------------------------------------
resultWindows %>%
   filter(significant) %>% 
   arrange(desc(log2FoldChange)) %>% 
   .[["gene_name"]] %>% 
   unique %>% 
   head(20)

## ----extractRegions, eval = TRUE, results = "hide"-----------------------
resultRegions <- extractRegions(windowRes  = resultWindows,
                                padjCol    = "p_adj_IHW",
                                padjThresh = 0.01, 
                                log2FoldChangeThresh = 0.5) %>% as_tibble

## ----resultRegions output, eval = TRUE-----------------------------------
resultRegions

## ----toBED output, eval = FALSE------------------------------------------
#  toBED(windowRes = resultWindows,
#        regionRes = resultRegions,
#        fileName  = "enrichedWindowsRegions.bed")

## ----notes for further dev, eval = FALSE, echo = FALSE-------------------
#  ## Plotting
#  
#  # Simple Enrichment
#  
#  ## Summing up IP
#  
#  
#  ## Plotting enrichment IP over SMI
#  
#  
#  
#  # Region-testing with DESeq2
#  
#  ## Loading data
#  
#  ### Loading Testdata
#  
#  #```{r loading region test data, eval = FALSE}
#  ## eval = TRUE
#  #data(YBX3eCLIPregion)
#  #```
#  
#  # FAQ

## ------------------------------------------------------------------------
sessionInfo()

