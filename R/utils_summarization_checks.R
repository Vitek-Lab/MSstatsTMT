#' Check validity of parameters to proteinSummarization function
#' @inheritParams proteinSummarization
#' @return TRUE invisibly if all parameters are valid
#' @keywords internal
.checkSummarizationParams = function(data, method, global_norm, reference_norm, 
                                     remove_norm_channel, remove_empty_channel,
                                     MBimpute, maxQuantileforCensored) {
  checkmate::assertChoice(method, c("LogSum", "Median", 
                                    "MedianPolish", "msstats"))
  checkmate::assertLogical(global_norm, any.missing = FALSE)
  checkmate::assertLogical(reference_norm, any.missing = FALSE)
  checkmate::assertLogical(remove_norm_channel, any.missing = FALSE)
  checkmate::assertLogical(remove_empty_channel, any.missing = FALSE)
  checkmate::assertLogical(MBimpute, any.missing = FALSE)
  
  cols = c("ProteinName", "PeptideSequence", "Charge", "PSM",
           "Mixture", "TechRepMixture", "Run",
           "Channel", "Condition", "BioReplicate", "Intensity")
  column_diff = setdiff(cols, colnames(data))
  if (length(column_diff) > 0) {
    stop(
      paste("** Please check the required input. Columns:",
            paste(column_diff, sep = ", ", collapse = ", "), "are missing.")
    )
  }
  invisible(TRUE)
}


#' Log parameters for proteinSummarization function
#' @inheritParams proteinSummarization
#' @return TRUE invisibly after logging successfully
#' @keywords internal
.logSummarizationParams = function(method, global_norm, reference_norm,
                                   remove_norm_channel, remove_empty_channel) {
  
  msg = paste(
    "Selected options:",
    paste("Method for protein summarization :", method),
    paste("Constant median normalization between channels :", global_norm),
    paste("Reference-channel based normalization between MS runs :", 
          reference_norm),
    paste("Remove 'Norm' channels before inference:", remove_norm_channel),
    paste("Remove 'Empty' channels before inference:", remove_empty_channel),
    sep = "\n"
  )
  getOption("MSstatsTMTLog")("INFO", msg)
  invisible(TRUE)
}

#' Remove empty and normalization channels
#' @inheritParams proteinSummarization
#' @param input data.table processed by the protein summarization function
#' @return data.table
#' @keywords internal
.removeRedundantChannels = function(input, remove_empty_channel, 
                                    remove_norm_channel) {
  Condition = NULL
  
  if (remove_empty_channel & is.element("Empty", unique(input$Condition))) {
    input = input[Condition != "Empty"]
  }
  
  if (remove_norm_channel & is.element("Norm", unique(input$Condition))) {
    input = input[Condition != "Norm"]
  }
  input$Condition = factor(input$Condition)
  input
}


#' Prepare TMT data for protein-level summarization
#' @param input data.table
#' @return data.table with required column types
#' @keywords internal
.prepareForSummarization = function(input) {
  ProteinName = Intensity = Run = Channel = NULL
  log2Intensity = RunChannel = PSM = NULL
  
  input[, log2Intensity := log(Intensity, 2)]
  input[, ProteinName := as.character(ProteinName)]
  input[, PSM := as.character(PSM)]
  input[, RunChannel := paste(Run, Channel, sep = "_")]
  
  if (any(!is.na(input$Intensity) & input$Intensity < 1)) {
    input[, log2Intensity := ifelse(!is.na(Intensity) & Intensity < 1,
                                    NA, log2Intensity)]
    msg = "** Negative log2 intensities were replaced with NA."
    getOption("MSstatsTMTMsg")("INFO", msg)
  }
  input
}
