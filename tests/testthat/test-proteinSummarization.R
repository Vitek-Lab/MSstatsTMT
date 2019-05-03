context("proteinSummarization")

test_that("proteinSummarization works", {
  
  output<-proteinSummarization(MSstatsTMT::input.pd, 
                                method = "msstats", normalization = TRUE)
  expect_equal(output, quant.pd.msstats, , tolerance=1e-4)
  
})

test_that("proteinSummarization can handle missing data", {
  #Peptide Sequence and Charge did not throw error
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-1],
                                      method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-2], 
                                      method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-3],
                                      method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-4],
                                      method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-5],
                                      method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-6], 
                                     method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-7], 
                                     method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-8], 
                                     method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-9], 
                                     method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-10], 
                                     method = "msstats",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-c(1,2)], 
                                     method = "msstats",normalization = TRUE))
  
})

test_that("proteinSummarization can handle missing input", {
  
  expect_error(proteinSummarization(MSstatsTMT::input.pd, 
                                     method = "missing",normalization = TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd, 
                                     method = "Median",normalization = "abc"))
  expect_error(proteinSummarization(MSstatsTMT::input.pd, 
                                     method = "msstats",MBimpute = "abc"))
  expect_error(proteinSummarization(MSstatsTMT::input.pd, 
                                     method = "msstats",maxQuantileforCensored = "abc"))
  
})
  