context("PDtoMSstatsTMTFormat")
# context("MaxQtoMSstatsTMTFormat")
# context("SpectroMinetoMSstatsTMTFormat")

test_that("PDtoMSstatsTMTFormatn works", {
  
  output <- PDtoMSstatsTMTFormat(input = MSstatsTMT::raw.pd, annotation = MSstatsTMT::annotation.pd)
  expect_equal(output,MSstatsTMT::input.pd)
  
})

# test_that("MaxQtoMSstatsTMTFormat", {
#   
#   output <- MaxQtoMSstatsTMTFormat(evidence = MSstatsTMT::evidence, proteinGroups = MSstatsTMT::proteinGroups, annotation = MSstatsTMT::annotation.mq)
#   expect_equal(output,MSstatsTMT::input.pd)
#   
# })

# test_that("SpectroMinetoMSstatsTMTFormat", {
#   
#   input.mine <- SpectroMinetoMSstatsTMTFormat(input = raw.mine, annotation.mine)
#   expect_equal(output,MSstatsTMT::input.pd)
#   
# })
