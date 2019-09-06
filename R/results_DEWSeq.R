#' @export
#' @title extract DEWseq results
#' @description This is a modified version of the \code{\link[DESeq2:results]{results}}  function.
#' This documentation is based on DESeq2 results function documentation\cr
#' Extract results from a DESeq analysis object\cr
#'
#' For further details, please refer documentation for \code{\link[DESeq2:results]{results}} function in DESeq2 package
#'
#' @details
#' The output data.frame from this function will have the following columns: \cr
#'  \code{chromosome} : chromosome name\cr
#'  \code{unique_id} : unique id of the window\cr
#'  \code{begin} : window start co-ordinate, see parameter \code{begin0based}\cr
#'  \code{end} : window end co-ordinate\cr
#'  \code{strand} : strand\cr
#'  \code{gene_id} : gene id\cr
#'  \code{gene_name} : gene name\cr
#'  \code{gene_type} : gene type annotation\cr
#'  \code{gene_region} : gene region\cr
#'  \code{Nr_of_region} : number of the current region\cr
#'  \code{Total_nr_of_region} : total number of regions\cr
#'  \code{window_number} : window number\cr
#'
#' The columns listed below are from DESeq2 \code{\link[DESeq2:results]{results}}, please consult DESeq2 vignettes for an explanation of these columns:\cr
#' \code{baseMean}, \code{log2FoldChange}, \code{lfcSE}, \code{stat}, \code{pvalue}, \code{padj}
#'
#' @param object a DESeqDataSet, on which one of the following functions has already been called:
#' \code{\link[DESeq2:nbinomWaldTest]{nbinomWaldTest}}\cr
#' \code{\link[DESeq2:nbinomLRT]{nbinomLRT}} is NOT supported in this version
#' @param annotationFile sliding window annotation file, can be plain either text or .gz file,
#' the file MUST be <TAB> separated, and MUST have the following columns:\cr
#' \code{chromosome}: chromosome name \cr
#' \code{unique_id}: unique id of the window, \code{rownames(object)} must match this column \cr
#' \code{begin}: window start co-ordinate, see parameter \code{begin0based} \cr
#' \code{end}: window end co-ordinate \cr
#' \code{strand}: strand \cr
#' \code{gene_id}: gene id \cr
#' \code{gene_name}: gene name \cr
#' \code{gene_type}: gene type annotation \cr
#' \code{gene_region}: gene region \cr
#' \code{Nr_of_region}: number of the current region \cr
#' \code{Total_nr_of_region}: total number of regions \cr
#' \code{window_number}: window number \cr
#' @param contrast this argument specifies what comparison to extract from the \code{object} to build a results table.
#' @param name the name of the individual effect (coefficient) for building a results table.
#' \code{name} argument is ignored if \code{contrast} is specified
#' @param cooksCutoff theshold on Cook's distance
#' @param test this is automatically detected internally if not provided.
#' @param addMLE if \code{betaPrior=TRUE} was used
#' @param tidy whether to output the results table with rownames as a first column 'row'.
#' The table will also be coerced to \code{data.frame}
#' @param parallel if FALSE, no parallelization. if TRUE, parallel
#' execution using \code{BiocParallel}, see next argument \code{BPPARAM}
#' @param BPPARAM an optional parameter object passed internally
#' to \code{\link{bplapply}} when \code{parallel=TRUE}.
#' If not specified, the parameters last registered with
#' \code{\link{register}} will be used.
#' @param minmu lower bound on the estimated count (used when calculating contrasts)
#' @param begin0based TRUE (default) or FALSE. If TRUE, then the start positions in \code{annotationFile} are  considered to be 0-based
#'
#' @return data.frame
results_DEWSeq <- function(object, annotationFile,contrast, name, listValues=c(1,-1), cooksCutoff, test,
                           addMLE=FALSE, tidy=FALSE, parallel=FALSE, BPPARAM=bpparam(), minmu=0.5,begin0based=TRUE) {
  stopifnot(is(object, "DESeqDataSet"))
  stopifnot(length(listValues)==2 & is.numeric(listValues))
  stopifnot(listValues[1] > 0 & listValues[2] < 0)
  if (!"results" %in% mcols(mcols(object))$type) {
    stop("couldn't find results. you should first run DESeq()")
  }
  if (missing(test)) {
    test <- attr(object, "test")
  }
  if (test == "Wald" & attr(object, "test") == "LRT") {
    # initially test was LRT, now need to add Wald statistics and p-values
    object <- DESeq2:::makeWaldTest(object)
  }
  if (test == "LRT") {
    stop("this function do not support Likelihood ratio test!")
  }

  if (addMLE) {
    if (!attr(object,"betaPrior")) {
      stop("addMLE=TRUE is only for when a beta prior was used. Otherwise, the log2 fold changes are already MLE")
    }
    if (!missing(name) & missing(contrast)) {
      stop("addMLE=TRUE should be used by providing character vector of length 3 to 'contrast' instead of using 'name'")
    }
  }

  if (!missing(contrast)) {
    if (attr(object,"modelMatrixType") == "user-supplied" & is.character(contrast)) {
      stop("only list- and numeric-type contrasts are supported for user-supplied model matrices")
    }
  }

  if (is(design(object), "formula")) {
    hasIntercept <- attr(terms(design(object)),"intercept") == 1
    isExpanded <- attr(object, "modelMatrixType") == "expanded"
    termsOrder <- attr(terms.formula(design(object)),"order")
    # if no intercept was used or an expanded model matrix was used,
    # and neither 'contrast' nor 'name' were specified,
    # and no interactions...
    # then we create the result table: last / first level for last variable
    if ((test == "Wald") & (isExpanded | !hasIntercept) & missing(contrast) & missing(name) & all(termsOrder < 2)) {
      designVars <- all.vars(design(object))
      lastVarName <- designVars[length(designVars)]
      lastVar <- colData(object)[[lastVarName]]
      if (is.factor(lastVar)) {
        nlvls <- nlevels(lastVar)
        contrast <- c(lastVarName, levels(lastVar)[nlvls], levels(lastVar)[1])
      }
    }
  }

  if (missing(name)) {
    name <- DESeq2:::lastCoefName(object)
  } else {
    if (length(name) != 1 | !is.character(name)) {
      stop("the argument 'name' should be a character vector of length 1")
    }
  }

  WaldResults <- paste0("WaldPvalue_",name) %in% names(mcols(object))

  # this will be used in cleanContrast, and in the lfcThreshold chunks below
  useT <- "tDegreesFreedom" %in% names(mcols(object))

  # if performing a contrast call the function cleanContrast()
  if (!missing(contrast)) {
    resNames <- resultsNames(object)
    # do some arg checking/cleaning
    contrast <- DESeq2:::checkContrast(contrast, resNames)

    ### cleanContrast call ###
    # need to go back to C++ code in order to build the beta covariance matrix
    # then this is multiplied by the numeric contrast to get the Wald statistic.
    # with 100s of samples, this can get slow, so offer parallelization
    if (!parallel) {
      res <- DESeq2:::cleanContrast(object, contrast, expanded=isExpanded, listValues=listValues,
                                    test=test, useT=useT, minmu=minmu)
    } else if (parallel) {
      # parallel execution
      nworkers <- getNworkers(BPPARAM)
      idx <- factor(sort(rep(seq_len(nworkers),length.out=nrow(object))))
      res <- do.call(rbind, bplapply(levels(idx), function(l) {
        DESeq2:::cleanContrast(object[idx == l,,drop=FALSE], contrast,
                               expanded=isExpanded, listValues=listValues,
                               test=test, useT=useT, minmu=minmu)
      }, BPPARAM=BPPARAM))
    }

  } else {
    # if not performing a contrast
    # pull relevant columns from mcols(object)
    log2FoldChange <- DESeq2:::getCoef(object, name)
    lfcSE <- DESeq2:::getCoefSE(object, name)
    stat <- DESeq2:::getStat(object, test, name)
    pvalue <- DESeq2:::getPvalue(object, test, name)
    res <- cbind(mcols(object)["baseMean"],log2FoldChange,lfcSE,stat,pvalue)
    names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
  }

  rownames(res) <- rownames(object)
  if(any(mcols(object)$allZero)){
    res <- res[-which(mcols(object)$allZero), ]
  }
  # add unshrunken MLE coefficients to the results table
  if (addMLE) {
    if (is.numeric(contrast)) stop("addMLE only implemented for: contrast=c('condition','B','A')")
    if (is.list(contrast)) stop("addMLE only implemented for: contrast=c('condition','B','A')")
    res <- cbind(res, DESeq2:::mleContrast(object, contrast))
    res <- res[,c("baseMean","log2FoldChange","lfcMLE","lfcSE","stat","pvalue")]
    # if an all zero contrast, also zero out the lfcMLE
    res$lfcMLE[ which(res$log2FoldChange == 0 & res$stat == 0) ] <- 0
  }
  # read the annotation table and keep it for later
  resGrange <- .readAnnotation(fname=annotationFile,uniqIds = rownames(res),begin0based=begin0based)
  gc()
  # prune res object for all regions/windows with stat>=0
  res <- res[res$stat>=0,]
  if(nrow(res)==0){
    stop("Cannot find any windows/regions with log2FoldChange>=0")
  }

  # recalculate p-values, in this case, only right sided test is enough
  if(useT){
    # keep degrees of freedom for the regions/windows with stat>=0
    df <- mcols(object)$tDegreesFreedom
    names(df) <- rownames(object)
    df <- df[rownames(res)]
    res$pvalue <- pt(res$stat,df=df,lower.tail = FALSE)
  }else{
    res$pvalue <- pnorm(res$stat,lower.tail = FALSE)
  }

  # calculate Cook's cutoff
  m <- nrow(attr(object,"dispModelMatrix"))
  p <- ncol(attr(object,"dispModelMatrix"))

  defaultCutoff <- qf(.99, p, m - p)
  if (missing(cooksCutoff)) {
    cooksCutoff <- defaultCutoff
  }
  stopifnot(length(cooksCutoff)==1)
  if (is.logical(cooksCutoff) & cooksCutoff) {
    cooksCutoff <- defaultCutoff
  }

  # apply cutoff based on maximum Cook's distance
  performCooksCutoff <- (is.numeric(cooksCutoff) | cooksCutoff)
  if (performCooksCutoff) {
    cooksOutlier <- mcols(object)$maxCooks > cooksCutoff

    ### BEGIN heuristic to avoid filtering genes with low count outliers
    # as according to Cook's cutoff. only for two group designs.
    # dont filter if three or more counts are larger
    if (any(cooksOutlier,na.rm=TRUE) & is(design(object), "formula")) {
      designVars <- all.vars(design(object))
      if (length(designVars) == 1) {
        var <- colData(object)[[designVars]]
        if (is(var, "factor") && nlevels(var) == 2) {
          dontFilter <- logical(sum(cooksOutlier,na.rm=TRUE))
          for (i in seq_along(dontFilter)) {
            # index along rows of object
            ii <- which(cooksOutlier)[i]
            # count for the outlier with max cooks
            outCount <- counts(object)[ii,which.max(assays(object)[["cooks"]][ii,])]
            # if three or more counts larger than the outlier
            if (sum(counts(object)[ii,] > outCount) >= 3) {
              # don't filter out the p-value for that gene
              dontFilter[i] <- TRUE
            }
          }
          # reset the outlier status for these genes
          cooksOutlier[which(cooksOutlier)][dontFilter] <- FALSE
        }
      }
    } ### END heuristic ###
    res <- res[ setdiff(rownames(res),rownames(object)[which(cooksOutlier)]),  ]

  }
  # calculate the number of windows that overlaps with each other by atleast 2 bp,
  # this way, any intron/exon junctions will not be counted, but all ovelapping windows will be counted
  resOvs <- GenomicRanges::findOverlaps(resGrange,minoverlap=1,drop.redundant=FALSE,drop.self=FALSE,ignore.strand=FALSE)
  nOvWindows <- as.data.frame(table(queryHits(resOvs)))[,2]
  names(nOvWindows) <- rownames(mcols(resGrange))
  nOvWindows <- nOvWindows[rownames(res)]
  # adjusted p-values for overlapping windows using Bonferroni correction
  res$pBonferroni <- pmin(res$pvalue*nOvWindows,1)
  # if original baseMean was positive, but now zero due to replaced counts, fill in results
  if ( sum(mcols(object)$replace, na.rm=TRUE) > 0) {
    nowZeroIds <- intersect(rownames(res),rownames(object)[which(mcols(object)$replace & mcols(object)$baseMean == 0)])
    if(length(nowZeroIds)>0){
      res[nowZeroIds,"log2FoldChange"] <- 0
      if (addMLE) { res[nowZeroIds,"lfcMLE"] <- 0 }
      res[nowZeroIds,"lfcSE"] <- 0
      res[nowZeroIds,"stat"] <- 0
      res[nowZeroIds,"pvalue"] <- 1
      res[nowZeroIds,"pBonferroni"] <- 1
    }
  }
  # correct the window p-values for FDR on the genome level
  res$pBonferroni.adj <- p.adjust(res[,'pBonferroni'], method = 'BH')
  # add prior information
  deseq2.version <- packageVersion("DESeq2")
  if (!attr(object,"betaPrior")) {
    priorInfo <- list(type="none", package="DESeq2", version=deseq2.version)
  } else {
    betaPriorVar <- attr(object, "betaPriorVar")
    priorInfo <- list(type="normal", package="DESeq2", version=deseq2.version,
                      betaPriorVar=betaPriorVar)
  }
  # make results object
  deseqRes <- DESeqResults(cbind(res,as.data.frame(resGrange[rownames(res),])), priorInfo=priorInfo)
  colnames(deseqRes) <- c('baseMean', 'log2FoldChange', 'lfcSE','stat', 'pvalue', 'pBonferroni','pBonferroni.adj', 'chromosome', 'begin','end',
                          'width',  'strand','unique_id', 'gene_id',  'gene_name','gene_type', 'gene_region', 'Nr_of_region','Total_nr_of_region',
                          'window_number')
  # 'seqnames' from Granges is always a factor!
  deseqRes$chromosome <- as.character(deseqRes$chromosome)
  if(begin0based){
    deseqRes$begin <- pmax(deseqRes$begin-1,0)
  }
  # may this helps to improve the memory usage
  rm(resOvs,nOvWindows,resGrange)
  gc()
  if (tidy) {
    colnms <- colnames(deseqRes)
    mcols(deseqRes,use.names=TRUE)["unique_id","type"] <- "results"
    mcols(deseqRes,use.names=TRUE)["unique_id","description"] <- "unique id for the window/region(row names)"
    deseqRes <- deseqRes[,c("unique_id",setdiff(colnms,'unique_id'))]
    rownames(deseqRes) <- NULL
    deseqRes <- as.data.frame(deseqRes)
  }
  return(deseqRes)
}
