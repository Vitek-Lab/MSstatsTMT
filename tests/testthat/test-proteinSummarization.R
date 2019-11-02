context("proteinSummarization")

test_that("proteinSummarization works", {
  
  output<-proteinSummarization(MSstatsTMT::input.pd, 
                               method = "msstats", 
                               global_norm=TRUE,
                               reference_norm=TRUE)
  expect_equal(output, quant.pd.msstats, tolerance=1e-4)
  
})

test_that("proteinSummarization can handle missing data", {
  #Peptide Sequence and Charge did not throw error
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-1],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-2],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-3],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-4],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-5],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-6],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-7],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-8],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-9],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-10],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd[,-c(1,2)],
                                    method = "msstats",global_norm=TRUE, 
                                    reference_norm=TRUE))
  
})

test_that("proteinSummarization can handle missing input", {
  
  expect_error(proteinSummarization(MSstatsTMT::input.pd, 
                                    method = "missing",
                                    global_norm=TRUE, 
                                    reference_norm=TRUE))
  expect_error(proteinSummarization(MSstatsTMT::input.pd, 
                                    method = "Median",
                                    global_norm = "abc", 
                                    reference_norm="true"))
  expect_error(proteinSummarization(MSstatsTMT::input.pd,
                                    method = "msstats",MBimpute = "abc"))
  expect_error(proteinSummarization(MSstatsTMT::input.pd, 
                                    method = "msstats",maxQuantileforCensored = "abc"))
  
})
  