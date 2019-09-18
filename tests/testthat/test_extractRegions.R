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
  expect_warning(extractRegions(windowRes = slbpWindows,padjCol = 'pSlidingWindows.adj',
                                log2FoldChangeThresh = 20),warnMsg)
  expect_null(suppressWarnings(extractRegions(windowRes = slbpWindows,padjCol = 'pSlidingWindows.adj',
                                              log2FoldChangeThresh = 20)))
  expectedResults <- c("ENSG00000184678.10:exon0001W00087", "ENSG00000184260.5:exon0001W00017", "ENSG00000196890.4:exon0001W00020",
                       "ENSG00000188486.3:exon0001W00053", "ENSG00000278705.1:exon0001W00004", "ENSG00000197061.4:exon0001W00020",
                       "ENSG00000180596.7:exon0001W00430", "ENSG00000180573.9:exon0001W00028", "ENSG00000168298.6:exon0001W00029",
                       "ENSG00000168298.6:exon0001W00034", "ENSG00000158373.8:exon0001W00018", "ENSG00000282988.1:exon0002W00061",
                       "ENSG00000282988.1:intron0001W00159", "ENSG00000277224.2:exon0001W00020", "ENSG00000276966.2:exon0001W00021",
                       "ENSG00000273802.2:exon0001W00053", "ENSG00000277075.2:exon0001W00022", "ENSG00000275713.2:exon0001W00018",
                       "ENSG00000158406.4:exon0001W00373", "ENSG00000124635.8:intron0001W00321", "ENSG00000196787.3:exon0001W00022",
                       "ENSG00000276180.1:exon0001W00052", "ENSG00000185130.5:exon0001W00002", "ENSG00000196747.4:exon0001W00022",
                       "ENSG00000278828.1:exon0001W00021", "ENSG00000197238.4:exon0001W00011", "ENSG00000273542.1:exon0001W00001",
                       "ENSG00000233822.4:exon0001W00065", "ENSG00000233822.4:exon0004W01400", "ENSG00000197153.4:intron0001W00115",
                       "ENSG00000274641.1:exon0001W00018", "ENSG00000169398.19:intron0010W12981")
  expect_equal(extractRegions(windowRes = slbpWindows,padjCol = 'pSlidingWindows.adj',log2FoldChangeThresh = 10)$regionStartId,expectedResults)
})
