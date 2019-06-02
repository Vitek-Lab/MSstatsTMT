context("groupComparision")

test_that("groupComparision works", {
  
  output<-groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats, moderated = TRUE)
  expect_equal(output, MSstatsTMT::test.pairwise, tolerance=1e-5)
  
})

test_that("groupComparision handle missing column", {
  
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats[,-1]))
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats[,-2]))
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats[,-3]))
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats[,-4]))
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats[,-5]))
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats[,-6]))
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats[,-7]))
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats[,-c(1,2)]))

})

test_that("groupComparision handle wrong input", {
  
  expect_error(groupComparisonTMT())
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats, 
                                   contrast.matrix = matrix(c(1,2,3,4),nrow = 2)))
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats,remove_norm_channel = "abc"))
  expect_error(groupComparisonTMT(data = MSstatsTMT::quant.pd.msstats, 
                                   moderated = "abc"))
  
})

