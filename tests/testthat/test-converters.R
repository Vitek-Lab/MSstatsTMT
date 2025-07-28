context("PDtoMSstatsTMTFormat")
context("MaxQtoMSstatsTMTFormat")
context("SpectroMinetoMSstatsTMTFormat")

test_that("PDtoMSstatsTMTFormat works", {

    expect_error(PDtoMSstatsTMTFormat(input = MSstatsTMT::raw.pd[, !colnames(MSstatsTMT::raw.pd) == "Protein.Accessions"], # missing columns in input
                                    annotation = MSstatsTMT::annotation.pd))

    expect_error(PDtoMSstatsTMTFormat(input = MSstatsTMT::raw.pd,
                                      annotation = MSstatsTMT::annotation.pd[, !colnames(MSstatsTMT::annotation.pd) == "Condition"])) # missing columns in annotation

})

test_that("MaxQtoMSstatsTMTFormat works", {

    expect_error(MaxQtoMSstatsTMTFormat(evidence = MSstatsTMT::evidence, proteinGroups = MSstatsTMT::proteinGroups,
                                        annotation = MSstatsTMT::annotation.pd)) # wrong annotation file

})

test_that("SpectroMinetoMSstatsTMTFormat works", {

    expect_error(SpectroMinetoMSstatsTMTFormat(input = MSstatsTMT::raw.mine,
                                        annotation = MSstatsTMT::annotation.mine,
                                        summaryforMultipleRows = average)) # wrong argument value

})

test_that("PhilosophertoMSstatsTMTFormat works", {
  
  input_file_path = system.file("raw_data/Philosopher/msstats.csv",
                       package = "MSstatsTMT")
  annotation_file_path = system.file("raw_data/Philosopher/MSstatsTMT_annotation.csv",
                      package = "MSstatsTMT")
  input = data.table::fread(input_file_path)
  annotation = data.table::fread(annotation_file_path)
  msstats_format = PhilosophertoMSstatsTMTFormat(input, annotation)
  expect_equal(nrow(msstats_format), 550)
  
})
