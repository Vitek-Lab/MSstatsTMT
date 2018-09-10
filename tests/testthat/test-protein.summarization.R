context("protein.summarization")

test_that("protein.summarization works",{
  output<-protein.summarization(required.input,method = "MedianPolish",normalization = TRUE)
  expect_equal(output,MedianPolish.quant)
})