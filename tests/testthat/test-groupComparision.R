context("groupComparision")

test_that("groupComparision works", {
  
  output<-groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats)
  expect_equal(output,MSstatsTMT::test.pairwise)
  
})

test_that("groupComparision handle missing column", {
  
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats[,-1]))
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats[,-2]))
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats[,-3]))
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats[,-4]))
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats[,-5]))
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats[,-6]))
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats[,-7]))
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats[,-c(1,2)]))

})

test_that("groupComparision handle wrong input", {
  
  expect_error(groupComparison.TMT())
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats, 
                                   contrast.matrix = matrix(c(1,2,3,4),nrow = 2)))
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats,remove_norm_channel = "abc"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::quant.pd.msstats, 
                                   moderated = "abc"))
  
})

