## ----setup, echo=FALSE, results="hide"-----------------------------------
#  DEWSeq package version: `r packageVersion("DESeq2")` <- this goes into abstract

knitr::opts_chunk$set(tidy=FALSE, cache=FALSE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)

## ---- eval = F, echo = F-------------------------------------------------
#  #**Note:** if you use DEWSeq in published research, please cite:
#  #
#  #> Authors. (Year)
#  #> Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.
#  #> *Genome Biology*, **15**:550.
#  #> [10.1186/s13059-014-0550-8](http://dx.doi.org/10.1186/s13059-014-0550-8)

