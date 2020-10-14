#' Protein summarization for TMT data
#' @param input data.table with TM quant data
#' @param method Four different summarization methods to protein-level 
#' can be performed : "msstats"(default), "MedianPolish", "Median", "LogSum".
#' @inheritParams .summarizeMSstats
#' 
#' @return data.table
#' 
#' @export
#' 
MSstatsSummarizeTMT = function(input, method, impute, fill_incomplete,
                               max_quantile_censored = NULL
) {
  log2Intensity = NULL
  
  annotation = unique(input[!is.na(log2Intensity),
                            c("Run", "Channel", "BioReplicate", "Condition",
                              "Mixture", "TechRepMixture", "RunChannel"),
                            with = FALSE])
  
  summarized = .summarizeTMT(input, method, annotation, impute, 
                             fill_incomplete, max_quantile_censored)
  .makeFactorColumnsTMT(summarized)
  summarized[, c("Run", "Protein", "Abundance", "Channel", "BioReplicate",
                 "Condition", "TechRepMixture", "Mixture"),
             with = FALSE]
}


#' Converts required columns to factor in summarization output
#' @param input data.table
#' @return NULL 
#' @keywords internal
.makeFactorColumnsTMT = function(input) {
  Run = Channel = Condition = TechRepMixture = Mixture = NULL
  
  input[, Run := as.factor(Run)]
  input[, Channel := as.factor(Channel)]
  input[, Condition := as.factor(Condition)]
  input[, TechRepMixture := as.factor(TechRepMixture)]
  input[, Mixture := as.factor(Mixture)]
}


#' Performs summarization for TMT data
#' @param input data.table
#' @param method "mstats"/"MedianPolish"/"LogSum"/"Median"
#' @inheritParams .summarizeMSstats
#' @return data.table
#' @keywords internal
.summarizeTMT = function(input, method, annotation, impute, fill_incomplete,
                         max_quantile_censored
) {
  if (method == "msstats") {
    summarized = .summarizeMSstats(input, annotation, impute, fill_incomplete,
                                   max_quantile_censored)
    method_msg = "MSstats"
  } else if (method == "MedianPolish") {
    summarized = .summarizeTMP(input, annotation)
    method_msg = "median polish"
  } else if (method == "LogSum") {
    summarized = .summarizeSimpleStat(input, annotation, .logSum)
    method_msg = "log(sum of intensities)"
  } else if (method == "Median") {
    summarized = .summarizeSimpleStat(input, annotation, median)
    method_msg = "median"
  }
  msg = paste0("** Protein-level summarization done by ", method_msg, ".")
  getOption("MSstatsLog")("INFO", msg)
  getOption("MSstatsMsg")("INFO", msg)
  summarized
}


#' Summarization based on MSstats
#' @param input data.table
#' @param annotation data.table with run and channel annotation
#' @param impute only for method="msstats". TRUE (default) imputes missing 
#' values by Accelated failure model. FALSE uses minimum value to impute the 
#' missing value for each peptide precursor ion.
#' @param fill_incomplete if TRUE, missing rows will be added with Intensity=NA.
#' @param max_quantile_censored We assume missing values are censored. 
#' maxQuantileforCensored is Maximum quantile for deciding censored missing 
#' value, for instance, 0.999. Default is Null.
#' @importFrom MSstatsdev dataProcess
#' @return data.table
#' @keywords internal
.summarizeMSstats = function(input, annotation, impute, fill_incomplete,
                             max_quantile_censored = NULL) {
  MSRun = NULL
  
  runs = na.omit(unique(annotation$Run))
  num_runs = length(runs)
  
  data.table::setnames(input, c("Run", "RunChannel", "Charge"),
                       c("MSRun", "Run", "PrecursorCharge"))  
  input[, FragmentIon := NA]
  input[, ProductCharge := NA]
  input[, IsotopeLabelType := "L"]
  
  summarized_results = vector("list", num_runs)
  
  for (i in seq_len(num_runs)) {
    msstats_summary = MSstatsdev::dataProcess(
      as.data.frame(input[MSRun == runs[i],
            list(ProteinName, PeptideSequence, PrecursorCharge,
                 FragmentIon, ProductCharge, Run, Condition,
                 BioReplicate, Intensity, IsotopeLabelType,
                 Fraction = 1)]),
      normalization = FALSE,
      summaryMethod = "TMP",
      censoredInt = "NA",
      MBimpute = impute,
      fillIncompleteRows = fill_incomplete,
      maxQuantileforCensored = max_quantile_censored
    )
    msstats_summary = msstats_summary$RunlevelData
    msstats_summary = msstats_summary[, c("Protein", "LogIntensities",
                                          "originalRUN")]
    summarized_results[[i]] = msstats_summary
  }

  summarized_results = data.table::rbindlist(summarized_results)
  data.table::setnames(summarized_results,
                       c("LogIntensities", "originalRUN"),
                       c("Abundance", "RunChannel"))
  summarized_results = merge(summarized_results, annotation,
                             by = "RunChannel", all.x = TRUE)

  data.table::setnames(input, c("MSRun", "Run", "PrecursorCharge"),
                       c("Run", "RunChannel", "Charge"))
  input[, FragmentIon := NULL]
  input[, ProductCharge := NULL]
  input[, IsotopeLabelType := NULL]

  summarized_results[, colnames(summarized_results) != "RunChannel",
                     with = FALSE]
  summarized_results
}


#' Summarize TMT data with median polish
#' @param input data.table
#' @param annotation data.table with run and channel annotation
.summarizeTMP = function(input, annotation) {
  log2Intensity = Run = Channel = ProteinName = RunChannel = NULL
  
  channel_len = data.table::uniqueN(annotation$Channel, na.rm = TRUE)

  input = input[order(ProteinName, Run), ]
  new_annotation = unique(input[, list(ProteinName, Run, 
                                       RunChannel = paste(Run, Channel, sep = "_"))])
  summarized = input[!is.na(log2Intensity),
                     list(MedianPolish = .medianPolish(log2Intensity, 
                                                       channel_len)),
                     by = c("Run", "ProteinName")]
  data.table::setnames(summarized, colnames(summarized),
                       c("Run", "Protein", "Abundance"))  
  summarized[, RunChannel := new_annotation$RunChannel]
  summarized = merge(summarized, 
                     annotation[, colnames(annotation) != "Run", with = FALSE], 
                     by = "RunChannel", all.x = TRUE)  
  summarized
}


#' Utility function: compute log of sum of 2^x
#' @param x numeric
#' @return numeric
#' @keywords internal
.logSum = function(x) {
  log(sum(2 ^ x, na.rm = TRUE), 2)
}


#' Summarize TMT data with a simple aggregate of log-intensities
#' @param input data.table
#' @param annotation data.table with run and channel annotation
#' @param stat_aggregate function that will be used to compute protein-level
#' summary
#' @return data.table
#' @keywords internal
.summarizeSimpleStat = function(input, annotation, stat_aggregate) {
  log2Intensity = NULL
  
  summarized = input[!is.na(log2Intensity), 
                     list(Median = stat_aggregate(log2Intensity)),
                     by = c("Run", "ProteinName", "RunChannel")]
  data.table::setnames(summarized, 
                       colnames(summarized),
                       c("Run", "Protein", "RunChannel", "Abundance"))
  summarized = merge(summarized, annotation[, colnames(annotation) != "Run",
                                            with = FALSE],
                     by = "RunChannel", all.x = TRUE)
  summarized
}


#' Tukey median polish
#' @param intensities vector of log-intensities per protein and run
#' @param num_channels number of channels
#' @importFrom stats medpolish
#' @return 
#' @keywords internal
.medianPolish <- function(intensities, num_channels){
  num_prot_runs = length(intensities)
  wide = matrix(intensities)
  dim(wide) = c(num_prot_runs / num_channels, num_channels)
  tmp_fit = stats::medpolish(wide, na.rm = TRUE, trace.iter = FALSE)
  tmp_fit$overall + tmp_fit$col
}
