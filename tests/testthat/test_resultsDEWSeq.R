# context("resultsDEWSeq")
# test_that("resultsDEWSeq throws expected errors", {
#   set.seed(1)
#   expect_error(resultsDEWSeq(matrix(0,10,4)),
#                'object MUST be of class DESeqDataSet!')
#   data("slbpDds")
#   slbpDds <- estimateSizeFactors(slbpDds)
#   slbpDds <- estimateDispersions(slbpDds)
#   slbpDds <- nbinomLRT(slbpDds,reduced = ~1)
#   expect_error(resultsDEWSeq(slbpDds),
#               'this function do not support likelihood ratio test!')
# })
