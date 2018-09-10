context("groupComparision")

test_that("groupComparision works",{
  
  output<-groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant,model = "proposed")
  expect_equal(output,MSstatsTMT::groupComparision.results)
  
})
