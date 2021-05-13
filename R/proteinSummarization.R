#' Summarizing peptide level quantification to protein level quantification
#'
#' We assume missing values are censored and then impute the missing values. Protein-level summarization from peptide level quantification are performed.
#' After all, global median normalization on peptide level data and normalization between MS runs using reference channels will be implemented.
#'
#' @param data Name of the output of PDtoMSstatsTMTFormat function or peptide-level quantified data from other tools. 
#' It should have columns ProteinName, PeptideSequence, Charge, PSM, Mixture, TechRepMixture, Run, Channel, Condition, BioReplicate, Intensity
#' @param method Four different summarization methods to protein-level can be performed : "msstats"(default), "MedianPolish", "Median", "LogSum".
#' @param global_norm  Global median normalization on peptide level data (equalizing the medians across all the channels and MS runs). Default is TRUE. 
#' It will be performed before protein-level summarization.
#' @param reference_norm Reference channel based normalization between MS runs on protein level data. 
#' TRUE(default) needs at least one reference channel in each MS run, annotated by 'Norm' in Condtion column. 
#' It will be performed after protein-level summarization. FALSE will not perform this normalization step. 
#' If data only has one run, then reference_norm=FALSE.
#' @param MBimpute only for method="msstats". TRUE (default) imputes missing values by Accelated failure model. 
#' FALSE uses minimum value to impute the missing value for each peptide precursor ion.
#' @param maxQuantileforCensored We assume missing values are censored. 
#' maxQuantileforCensored is Maximum quantile for deciding censored missing value, for instance, 0.999. Default is Null.
#' @param remove_norm_channel TRUE(default) removes 'Norm' channels from protein level data.
#' @param remove_empty_channel TRUE(default) removes 'Empty' channels from protein level data.
#' @param msstats_log_path path to a MSstats log file
#' @inheritParams .documentFunction
#'
#' @importFrom MSstatsConvert MSstatsSaveSessionInfo
#' @importFrom stats lm median model.matrix na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar
#'   
#' @return list that consists of two data.frames with feature-level (FeatureLevelData) and protein-level data (ProteinLevelData)
#' 
#' @export
#' 
#' @examples
#' data(input.pd)
#' quant.pd.msstats <- proteinSummarization(input.pd,
#'                                          method = "msstats",
#'                                          global_norm = TRUE,
#'                                          reference_norm = TRUE)
#' head(quant.pd.msstats$ProteinLevelData)
#' 
proteinSummarization = function(
  data, method = 'msstats', global_norm = TRUE, reference_norm = TRUE,
  remove_norm_channel = TRUE, remove_empty_channel = TRUE, MBimpute = TRUE,
  maxQuantileforCensored = NULL, use_log_file = TRUE, append = FALSE, 
  verbose = TRUE, log_file_path = NULL, msstats_log_path = NULL
){
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path, 
                                      base = "MSstatsTMT_summarization_log_",
                                      pkg_name = "MSstatsTMT")
  input = MSstatsPrepareForSummarizationTMT(
    data, method, global_norm, reference_norm,remove_norm_channel, 
    remove_empty_channel, MBimpute, maxQuantileforCensored
  )
  input = MSstatsNormalizeTMT(input, "peptides", global_norm)
  summarized = MSstatsSummarizeTMT(input,
                                   method,
                                   MBimpute,
                                   maxQuantileforCensored,
                                   msstats_log_path)
  processed = getProcessedTMT(summarized, input)
  summarized = getSummarizedTMT(summarized)
  summarized = MSstatsNormalizeTMT(summarized, "proteins", reference_norm)
  output = MSstatsSummarizationOutputTMT(summarized, processed,
                                         remove_empty_channel,
                                         remove_norm_channel)
  output
}


#' Get processed feature-level data
#' 
#' @param summarized output of the MSstatsSummarizeTMT function
#' @param input output of MSstatsNormalizeTMT function
#' 
#' @return data.table
#' 
#' @keywords internal
#' 
getProcessedTMT = function(summarized, input) {
  if (is.list(summarized) & !is.data.table(summarized)) {
    processed = summarized[[2]]
  } else {
    processed = input[, `:=`(censored = FALSE, 
                             predicted = NA, 
                             Intensity = NULL, 
                             RunChannel = NULL)]
  }
  processed
}


#' Get protein-level data from MSstatsSummarizeTMT output
#' 
#' @param summarized output of the MSstatsSummarizeTMT function
#' 
#' @return data.table
#' 
#' @keywords internal
#' 
getSummarizedTMT = function(summarized) {
  if (is.list(summarized) & !is.data.table(summarized)) {
    summarized = summarized[[1]]
  }
  summarized = .makeFactorColumnsTMT(summarized)
  summarized[, c("Run", "Protein", "Abundance", "Channel", "BioReplicate",
                 "Condition", "TechRepMixture", "Mixture"),
             with = FALSE]
}


#' Prepare output of MSstatsTMT converters for protein-level summarization
#' 
#' @inheritParams proteinSummarization
#' 
#' @return data.table
#' 
#' @keywords internal
#' 
MSstatsPrepareForSummarizationTMT = function(
  data, method, global_norm, reference_norm,remove_norm_channel, 
  remove_empty_channel, MBimpute, maxQuantileforCensored
) {
  current_msstats_log = getOption("MSstatsLog")
  current_msstats_msg = getOption("MSstatsMsg")
  MSstatsConvert::MSstatsLogsSettings(FALSE, FALSE, FALSE, NULL)
  getOption("MSstatsTMTLog")("INFO", 
                             "** MSstatsTMT - proteinSummarization function")
  getOption("MSstatsTMTMsg")("INFO", 
                             "** MSstatsTMT - proteinSummarization function")
  .checkSummarizationParams(data, method, global_norm, reference_norm, 
                            remove_norm_channel, remove_empty_channel,
                            MBimpute, maxQuantileforCensored)
  .logSummarizationParams(method, global_norm, reference_norm,
                          remove_norm_channel, remove_empty_channel)
  is_validated = inherits(data, "MSstatsValidated")
  input = data.table::as.data.table(unclass(data))
  if (!is_validated) {
    data.table::setnames(input, "Charge", "PrecursorCharge")
    input = MSstatsConvert::MSstatsBalancedDesign(
      input, c("PeptideSequence", "PrecursorCharge"),
      TRUE, FALSE, "zero_to_na"
    )
    input = data.table::as.data.table(unclass(input))
    data.table::setnames(input, "PrecursorCharge", "Charge")
  }
  options(MSstatsLog = current_msstats_log,
          MSstatsMsg = current_msstats_msg)
  input = .prepareForSummarization(input)
  input
}


#' Combine feature-level and protein-level data into single output
#'
#' @param summarized output of the getSummarizedTMT function
#' @param processed output of the getProcessedTMT function
#' @inheritParams proteinSummarization
#' 
#' @return list that consists of two data.frames with feature-level and protein-level data
#' 
#' @keywords internal
#' 
MSstatsSummarizationOutputTMT = function(summarized, processed,
                                         remove_empty_channel,
                                         remove_norm_channel) {
  summarized = .removeRedundantChannels(summarized, remove_empty_channel,
                                        remove_norm_channel)
  list(FeatureLevelData = as.data.frame(processed),
       ProteinLevelData = as.data.frame(summarized))
}
