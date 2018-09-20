context("protein.summarization")

test_that("protein.summarization works", {
  
  output<-protein.summarization(MSstatsTMT::required.input, 
                                method = "MedianPolish", normalization = TRUE)
  expect_equal(output, MedianPolish.quant)
  
})


test_that("protein.summarization can handle missing data", {
  #Peptide Sequence and Charge did not throw error
  expect_error(protein.summarization(MSstatsTMT::required.input[,-1],
                                      method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-2], 
                                      method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-3],
                                      method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-4],
                                      method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-5],
                                      method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-6], 
                                     method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-7], 
                                     method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-8], 
                                     method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-9], 
                                     method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-10], 
                                     method = "MedianPolish",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input[,-c(1,2)], 
                                     method = "MedianPolish",normalization = TRUE))
  
})

test_that("protein.summarization can handle missing input", {
  
  expect_error(protein.summarization(MSstatsTMT::required.input, 
                                     method = "missing",normalization = TRUE))
  expect_error(protein.summarization(MSstatsTMT::required.input, 
                                     method = "Median",normalization = "abc"))
  expect_error(protein.summarization(MSstatsTMT::required.input, 
                                     method = "MSstats",MBimpute = "abc"))
  expect_error(protein.summarization(MSstatsTMT::required.input, 
                                     method = "MSstats",maxQuantileforCensored = "abc"))
  
})
  