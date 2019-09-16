## ----setup, echo=FALSE, results="hide"-------------------------------------
knitr::opts_chunk$set(tidy    = FALSE,
                      cache   = FALSE,
                      dev     = "png",
                      message = FALSE,
                      error   = FALSE,
                      warning = TRUE)

## ---- eval = FALSE, echo = FALSE-------------------------------------------
#  #**Note:** if you use DEWSeq in published research, please cite:
#  #
#  #> Authors. (Year)
#  #> Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.
#  #> *Genome Biology*, **15**:550.
#  #> [10.1186/s13059-014-0550-8](http://dx.doi.org/10.1186/s13059-014-0550-8)

## ---- fig.cap = "Crosslink site trunction at reverse transcription", echo =FALSE----
knitr::include_graphics("truncation.png")

## ---- fig.cap = "Truncation sites (often referred to as crosslink sites) are the neighboring position of the aligned read.", fig.small = TRUE, echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("truncationsite.png")

## ----libs, message=FALSE---------------------------------------------------
library("IHW")
library("tidyverse")

## ---- fig.cap = "Different binding modes of RNA-binding proteins", fig.wide = TRUE, echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("binding_modes.png")

## ---- fig.cap = "Chance of crosslinking", fig.small = TRUE, echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("crosslinking_chance.png")

## ---- fig.cap = "Sliding window approach with single-nucleotide position count data", echo = FALSE----
# fig.width=200, out.width="200px"
knitr::include_graphics("dewseqoverview.png")

## ---- fig.cap = "Sliding window approach with single-nucleotide position count data", echo = FALSE, fig.small = TRUE----
# fig.width=200, out.width="200px"
knitr::include_graphics("overlapping_sliding_windows.png")

## ---- echo = FALSE, eval = FALSE-------------------------------------------
#  # Analysis strategies for CLIP data

## ----load library, eval = TRUE---------------------------------------------
require(DEWSeq)

## ----load test data, eval = TRUE-------------------------------------------
data(SLBP_K562_w50s20)

## ----head of test data, eval = TRUE----------------------------------------
ddw <- SLBP_K562_w50s20
ddw

## ----loading tidyverse, eval = TRUE----------------------------------------
library(tidyverse)
library(data.table)

## ----read in count matrix, eval = FALSE------------------------------------
#  countData <- fread("path/swcounts/count_matrix.txt.gz", sep = "\t")
#  
#  # or alternatively with tidyr:
#  # countData <- read_tsv("path/swcounts/count_matrix.txt.gz")

## ----read in annotation, eval = FALSE--------------------------------------
#  annotationData <- fread("path/annotation/annotation.txt.gz", sep = "\t")
#  
#  # or alternatively with tidyr:
#  #annotationData <- read_tsv("path/annotation/annotation.txt.gz")

## ----create colData, eval = FALSE------------------------------------------
#  colData <- data.frame(
#    row.names = colnames(countData),
#    type      = factor(
#                  c(rep("IP", 3),    ##  change this accordingly
#                    rep("SMI", 3)),  ##
#                  levels = c("IP", "SMI"))
#  )

## ----example import, eval = FALSE------------------------------------------
#  ddw <- DESeqDataSetFromSlidingWindows(countData  = countData,
#                                        colData    = colData,
#                                        annotation = annotationData,
#                                        tidy       = TRUE,
#                                        design     = ~type)

## ----row filtering, eval = TRUE--------------------------------------------
keep <- rowSums(counts(ddw)) >= 10
ddw <- ddw[keep,]

## ----estimate size factors, eval = FALSE-----------------------------------
#  ddw <- estimateSizeFactors(ddw)

## ----filter for mRNAs, eval = TRUE-----------------------------------------
ddw_mRNAs <- ddw[ rowData(ddw)[,"gene_type"] == "protein_coding", ]
ddw <- estimateSizeFactors(ddw_mRNAs)
rm(ddw_mRNAs)

## ----estimate dispersions and wald test, eval = TRUE-----------------------
ddw <- estimateDispersions(ddw)
ddw <- nbinomWaldTest(ddw)

## ----DEWSeq results, eval = TRUE-------------------------------------------
resultWindows <- resultsDEWSeq(ddw,
                              contrast = c("type", "IP", "SMI"),
                              tidy = TRUE) %>% as_tibble

resultWindows

## ----IHW multiple hypothesis correction, eval = TRUE, warning = FALSE------
resultWindows[,"p_adj_IHW"] <- adj_pvalues(ihw(pBonferroni ~ baseMean, 
                     data = resultWindows,
                     alpha = 0.05,
                     nfolds = 10))

## ----windows tables, eval = TRUE-------------------------------------------
resultWindows <- resultWindows %>% mutate(significant = resultWindows$p_adj_IHW < 0.05)

sum(resultWindows$significant)

## ----how genes form windows------------------------------------------------
resultWindows %>% filter(significant) %>% arrange(desc(log2FoldChange))  %>% .[["gene_name"]] %>% unique %>% head(20)

## ----extractRegions, eval = TRUE, results = "hide"-------------------------
resultRegions <- extractRegions(windowRes  = resultWindows,
                                padjCol    = "p_adj_IHW",
                                padjThresh = 0.05, 
                                log2FoldChangeThresh = 0.5) %>% as_tibble

## ----resultRegions output, eval = TRUE-------------------------------------
resultRegions

## ----toBED output, eval = FALSE--------------------------------------------
#  toBED(windowRes = resultWindows,
#        regionRes = resultRegions,
#        fileName  = "enrichedWindowsRegions.bed")

## ----eval = FALSE, echo = FALSE--------------------------------------------
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

## --------------------------------------------------------------------------
sessionInfo()

