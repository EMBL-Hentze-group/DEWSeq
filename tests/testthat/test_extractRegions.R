context("extractRegions")
test_that("extractRegions throws expected errors and warnings",{
  set.seed(1)
  data("slbpWindows")
  errorMsg <- 'Input data.frame is missing required columns, needed columns:
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
     padj: p-adjusted value column
     log2FoldChange: log2foldchange column.
      Missing columns: padj'
  expect_error(extractRegions(windowRes = slbpWindows),errorMsg)
  warnMsg <- 'There are no significant windows/regions under the current threshold!
    Please lower your significance cut-off thresholds and manually check if there are any significant windows under the threshold'
  expect_warning(extractRegions(windowRes = slbpWindows,padjCol = 'pBonferroni.adj',
                                log2FoldChangeThresh = 20),warnMsg)
})
