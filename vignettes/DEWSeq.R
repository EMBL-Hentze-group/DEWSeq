## ----setup, echo=FALSE, results="hide"-------------------------------------
knitr::opts_chunk$set(tidy    = FALSE,
                      cache   = FALSE,
                      dev     = "png",
                      message = FALSE,
                      error   = FALSE,
                      warning = TRUE)

## ---- eval = F, echo = F---------------------------------------------------
#  #**Note:** if you use DEWSeq in published research, please cite:
#  #
#  #> Authors. (Year)
#  #> Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.
#  #> *Genome Biology*, **15**:550.
#  #> [10.1186/s13059-014-0550-8](http://dx.doi.org/10.1186/s13059-014-0550-8)

## ---- fig.cap = "Crosslink site trunction at reverse transcription", echo =F----
knitr::include_graphics("truncation.png")

## ---- fig.cap = "Truncation sites (often referred to as crosslink sites) are the neighboring position of the aligned read.", fig.small = T, echo = F----
# fig.width=200, out.width="200px"
knitr::include_graphics("truncationsite.png")

## ---- fig.cap = "Different binding modes of RNA-binding proteins", fig.wide = T, echo = F----
# fig.width=200, out.width="200px"
knitr::include_graphics("binding_modes.png")

## ---- fig.cap = "Chance of crosslinking", fig.small = T, echo = F----------
# fig.width=200, out.width="200px"
knitr::include_graphics("crosslinking_chance.png")

## ---- fig.cap = "Sliding window approach with single-nucleotide position count data", echo = F----
# fig.width=200, out.width="200px"
knitr::include_graphics("dewseqoverview.png")

## ---- fig.cap = "Sliding window approach with single-nucleotide position count data", echo = F, fig.small = T----
# fig.width=200, out.width="200px"
knitr::include_graphics("overlapping_sliding_windows.png")

## ---- echo = F, eval = F---------------------------------------------------
#  # Analysis strategies for CLIP data

## ----load library, eval = F------------------------------------------------
#  #eval = T
#  require(DEWSeq)

## ----load test data, eval = F----------------------------------------------
#  # eval = T
#  data(YBX3eCLIPChr1)

## ----head of test data, eval = F-------------------------------------------
#  #eval = T
#  YBX3eCLIPChr1

## ----loading tidyverse, eval = F-------------------------------------------
#  library(tidyverse)

## ----read in count matrix, eval = F----------------------------------------
#  countData <- fread("path/swcounts/count_matrix.txt.gz", sep = "\t")
#  
#  # or alternatively with tidyr:
#  # countData <- read_tsv("path/swcounts/count_matrix.txt.gz")

## ----read in annotation, eval = F------------------------------------------
#  annotationData <- fread("path/annotation/annotation.txt.gz", sep = "\t")
#  
#  # or alternatively with tidyr:
#  #annotationData <- read_tsv("path/annotation/annotation.txt.gz")

## ----create colData, eval = F----------------------------------------------
#  colData <- data.frame(
#    row.names = colnames(countData),
#    type      = factor(
#                  c(rep("IP", 3),    ##  change this accordingly
#                    rep("SMI", 3)),  ##
#                  levels = c("IP", "SMI"))
#  )

## ----example import, eval = F----------------------------------------------
#  ddw <- DESeqDataSetFromSlidingWindows(countData  = countData,
#                                        colData    = colData,
#                                        annotation = annotationData,
#                                        tidy       = T,
#                                        design     = ~type)

## ----row filtering, eval = F-----------------------------------------------
#  # eval = T
#  keep <- rowSums(counts(ddw)) >= 20
#  ddw <- ddw[keep,]

## ----estimate size factors, eval = F---------------------------------------
#  ddw <- estimateSizeFactors(ddw)

## ----filter for mRNAs, eval = F--------------------------------------------
#  # eval T
#  ddw_mRNAs <- dds[ mcol(srowRanges(dds))[,"gene_type"] == "protein_coding", ]
#  ddw <- estimateSizeFactors(ddw_mRNAs)
#  rm(ddw_mRNAs)

## ---- eval = F-------------------------------------------------------------
#  # eval T
#  ddw <- estimateDispersions(ddw)
#  ddw <- nbinomWaldTest(ddw)

## ----results, eval = F-----------------------------------------------------
#  # eval = T
#  resultWindows <- resultsDEWSeq(dds,
#                                contrast = c("type", "IP", "SMI"),
#                                tidy = T)

## ----IHW, eval = F---------------------------------------------------------
#  # eval = T
#  
#  suppressPackageStartupMessages(require(IHW))
#  
#  results[,"p_adj_IHW"] <- ihw(pBonferroni ~ baseMean,
#                       data = resultWindows,
#                       alpha = 0.05,
#                       nfolds = 10)

## ----extractRegions, eval = F----------------------------------------------
#  # eval = T
#  resultRegions <- extractRegions(windowRes = resultWindows,
#                                  padjCol   = "p_adj_IHW",
#                                  padjThresh = 0.05,
#                                  log2FoldChangeThresh = 0.5)

## ----resultRegions output, eval = F----------------------------------------
#  # eval = T
#  head(resultRegions)

## ----eval = F, echo = F----------------------------------------------------
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
#  #```{r loading region test data, eval = F}
#  ## eval = T
#  #data(YBX3eCLIPregion)
#  #```
#  
#  # Discussion
#  # FAQ

## --------------------------------------------------------------------------
sessionInfo()

