#' Summarizing peptide level quantification to protein level quantification
#'
#' We assume missing values are censored and then impute the missing values. Protein-level summarization from peptide level quantification are performed.
#' After all, global median normalization on peptide level data and normalization between MS runs using reference channels will be implemented.
#'
#' @param data Name of the output of PDtoMSstatsTMTFormat function or peptide-level quantified data from other tools. It should have columns ProteinName, PeptideSequence, Charge, PSM, Mixture, TechRepMixture, Run, Channel, Condition, BioReplicate, Intensity
#' @param method Four different summarization methods to protein-level can be performed : "msstats"(default), "MedianPolish", "Median", "LogSum".
#' @param global_norm  Global median normalization on peptide level data (equalizing the medians across all the channels and MS runs). Default is TRUE. It will be performed before protein-level summarization.
#' @param reference_norm Reference channel based normalization between MS runs on protein level data. TRUE(default) needs at least one reference channel in each MS run, annotated by 'Norm' in Condtion column. It will be performed after protein-level summarization. FALSE will not perform this normalization step. If data only has one run, then reference_norm=FALSE.
#' @param MBimpute only for method="msstats". TRUE (default) imputes missing values by Accelated failure model. FALSE uses minimum value to impute the missing value for each peptide precursor ion.
#' @param maxQuantileforCensored We assume missing values are censored. maxQuantileforCensored is Maximum quantile for deciding censored missing value, for instance, 0.999. Default is Null.
#' @param remove_norm_channel TRUE(default) removes 'Norm' channels from protein level data.
#' @param remove_empty_channel TRUE(default) removes 'Empty' channels from protein level data.
#'   
#' @return data.frame with protein-level summarization for each run and channel
#' 
#' @export
#' 
#' @examples
#' data(input.pd)
#' quant.pd.msstats <- proteinSummarization(input.pd,
#'                                          method = "msstats",
#'                                          global_norm = TRUE,
#'                                          reference_norm = TRUE)
#' head(quant.pd.msstats)
#' 
proteinSummarization = function(
  data, method = 'msstats', global_norm = TRUE, reference_norm = TRUE,
  remove_norm_channel = TRUE, remove_empty_channel = TRUE, MBimpute = TRUE,
  maxQuantileforCensored = NULL
){

  getOption("MSstatsLog")("INFO", "** MSstatsTMT - proteinSummarization function")
  getOption("MSstatsLog")("INFO", "** MSstatsTMT - proteinSummarization function")
  
  .checkSummarizationParams(data, method, global_norm, reference_norm, 
                            remove_norm_channel, remove_empty_channel,
                            MBimpute, maxQuantileforCensored)
  .logSummarizationParams(method, global_norm, reference_norm,
                          remove_norm_channel, remove_empty_channel)
  
  input = MSstatsNormalizeTMT(data, "peptides", global_norm)
  summarized = MSstatsSummarizeTMT(input,
                                   method,
                                   global_norm,
                                   reference_norm,
                                   MBimpute,
                                   maxQuantileforCensored)
  n_runs = data.table::uniqueN(summarized$Run, na.rm = TRUE)
  summarized = MSstatsNormalizeTMT(summarized, "proteins", 
                                   reference_norm & n_runs > 1)
  summarized = .removeRedundantChannels(summarized, remove_empty_channel, 
                                        remove_norm_channel)
  as.data.frame(summarized)
}
