context("groupComparision")

test_that("groupComparision works", {

  output<-groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant,
                              model = "proposed")
  expect_equal(output,MSstatsTMT::groupComparision.results)

})

test_that("groupComparision handle missing column", {

  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant[, -1],
                                    model = "proposed"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant[, -2],
                                    model = "proposed"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant[, -3],
                                    model = "proposed"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant[, -4],
                                    model = "proposed"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant[, -5],
                                    model = "proposed"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant[, -6],
                                    model = "proposed"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant[, -7],
                                    model = "proposed"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant[, -c(1, 2)],
                                   model = "proposed"))

})

test_that("groupComparision handle wrong input", {

  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant,
                                   model = "Proposed"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant,
                                   contrast.matrix = matrix(c(1, 2, 3, 4), nrow = 2)))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant,
                                   model = "proposed", remove_norm_channel = "abc"))
  expect_error(groupComparison.TMT(data = MSstatsTMT::MedianPolish.quant,
                                   moderated = "abc"))

})

