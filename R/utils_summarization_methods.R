#' Protein summarization for TMT data
#' @param input data.table with TM quant data
#' @param method Four different summarization methods to protein-level 
#' can be performed : "msstats"(default), "MedianPolish", "Median", "LogSum".
#' @inheritParams .summarizeMSstats
#' 
#' @return data.table
#' 
#' @keywords internal
#' 
MSstatsSummarizeTMT = function(input, method, impute,
                               max_quantile_censored = NULL,
                               log_file_path = NULL) {
  log2Intensity = NULL
  
  annotation = unique(input[!is.na(log2Intensity),
                            c("Run", "Channel", "BioReplicate", "Condition",
                              "Mixture", "TechRepMixture", "RunChannel"),
                            with = FALSE])
  summarized = .summarizeTMT(input, method, annotation, impute, 
                             max_quantile_censored, log_file_path)
  summarized
}


#' Converts required columns to factor in summarization output
#' @param input data.table
#' @return a data table with factored columns
#' @keywords internal
.makeFactorColumnsTMT = function(input) {
  Run = Channel = Condition = TechRepMixture = Mixture = NULL
  
  input[, Run := as.factor(Run)]
  input[, Channel := as.factor(Channel)]
  input[, Condition := as.factor(Condition)]
  input[, TechRepMixture := as.factor(TechRepMixture)]
  input[, Mixture := as.factor(Mixture)]
  input
}


#' Performs summarization for TMT data
#' @param input data.table
#' @param method "mstats"/"MedianPolish"/"LogSum"/"Median"
#' @inheritParams .summarizeMSstats
#' @importFrom methods new
#' @return data.table
#' @keywords internal
.summarizeTMT = function(input, method, annotation, impute,
                         max_quantile_censored, log_file_path
) {
  
  FragmentIon <- ProductCharge <- IsotopeLabelType <- ProteinName <- 
    PeptideSequence <- PrecursorCharge <- Run <- Condition <- BioReplicate <- 
    Intensity <- PSM <- RunChannel <- NULL
  
  if (method == "msstats") {
    summarized = .summarizeMSstats(input, annotation, impute,
                                   max_quantile_censored, log_file_path)
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
  getOption("MSstatsTMTLog")("INFO", msg)
  getOption("MSstatsTMTMsg")("INFO", msg)
  summarized
}


#' Summarization based on MSstats
#' @param input data.table
#' @param annotation data.table with run and channel annotation
#' @param impute only for method="msstats". TRUE (default) imputes missing 
#' values by Accelated failure model. FALSE uses minimum value to impute the 
#' missing value for each peptide precursor ion.
#' @param max_quantile_censored We assume missing values are censored. 
#' maxQuantileforCensored is Maximum quantile for deciding censored missing 
#' value, for instance, 0.999. Default is Null.
#' @param log_file_path path to a MSstats log file
#' @importFrom MSstats dataProcess
#' @importFrom methods new
#' 
#' @return data.table
#' @keywords internal
.summarizeMSstats = function(input, annotation, impute, 
                             max_quantile_censored = NULL,
                             log_file_path = NULL) {
  MSRun = FragmentIon = ProductCharge = IsotopeLabelType = ProteinName = 
    PeptideSequence = PrecursorCharge = Run = Condition = BioReplicate =
    Intensity = PSM = RunChannel = NULL
  
  current_msstats_log = options("MSstatsLog")
  current_msstats_msg = options("MSstatsMsg")
  if (is.null(log_file_path)) {
    now = gsub("[ \\:\\-]", "_", Sys.time())
    log_file_path = paste0("MSstatsTMT_summarization_MSstats_", now, ".log")
  }
  MSstatsConvert::MSstatsLogsSettings(TRUE, TRUE, FALSE, log_file_path)
  
  runs = na.omit(unique(annotation$Run))
  num_runs = length(runs)
  
  data.table::setnames(input, c("Run", "RunChannel", "Charge"),
                       c("MSRun", "Run", "PrecursorCharge"))  
  input[, FragmentIon := NA]
  input[, ProductCharge := NA]
  input[, IsotopeLabelType := "L"]
  
  processed_data = vector("list", num_runs)
  summarized_results = vector("list", num_runs)
  for (i in seq_len(num_runs)) {
    ## For each run, use msstats dataprocess
    msg = paste("Summarizing for Run :", runs[i] ,
                "(", i, " of ", num_runs, ")")
    getOption("MSstatsTMTLog")("INFO", msg)
    getOption("MSstatsTMTMsg")("INFO", msg)
    
    single_run = input[MSRun == runs[i],
            list(ProteinName, PeptideSequence, PrecursorCharge,
                 FragmentIon, ProductCharge, Run, Condition,
                 BioReplicate, Intensity, IsotopeLabelType,
                 Fraction = 1)]
    single_run = new("MSstatsValidated", single_run)
    msstats_summary = MSstats::dataProcess(
      single_run,
      normalization = FALSE,
      summaryMethod = "TMP",
      censoredInt = "NA",
      MBimpute = impute,
      maxQuantileforCensored = max_quantile_censored,
      use_log_file = TRUE, append = TRUE, verbose = FALSE,
      log_file_path = log_file_path
    )
    feature_level_data = msstats_summary$FeatureLevelData 
    msstats_cols = c("PROTEIN", "PEPTIDE", "originalRUN", "censored",
                     "predicted", "newABUNDANCE")
    msstats_cols = intersect(msstats_cols, colnames(feature_level_data))
    feature_level_data = feature_level_data[, msstats_cols]
    processed_data[[i]] = feature_level_data
    
    protein_level_data = msstats_summary$ProteinLevelData
    protein_level_data = protein_level_data[, c("Protein", "LogIntensities",
                                          "originalRUN")]
    summarized_results[[i]] = protein_level_data
  }
  options(MSstatsLog = current_msstats_log,
          MSstatsMsg = current_msstats_msg)

  processed = data.table::rbindlist(processed_data)
  summarized_results = data.table::rbindlist(summarized_results)
  
  data.table::setnames(summarized_results,
                       c("LogIntensities", "originalRUN"),
                       c("Abundance", "RunChannel"))
  summarized_results = merge(summarized_results, annotation,
                             by = "RunChannel", all.x = TRUE)
  summarized_results = summarized_results[, colnames(summarized_results) != "RunChannel",
                                          with = FALSE]
  data.table::setnames(processed, 
                       c("PROTEIN", "PEPTIDE",
                         "originalRUN", "newABUNDANCE"),
                       c("ProteinName", "PSM", 
                         "RunChannel", "log2Intensity"))
  processed = merge(processed, annotation,
                    by = "RunChannel", all.x = TRUE)
  processed[, c("PeptideSequence", "Charge") := tstrsplit(PSM, "_", fixed=TRUE)]
  processed[, RunChannel := NULL]
  list(summarized_results, processed)
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


#' Summarize TMT data with median polish
#' @param input data.table
#' @param annotation data.table with run and channel annotation
#' @keywords internal
#' @return data.table with summaried protein intensities
.summarizeTMP = function(input, annotation) {
  log2Intensity = Run = Channel = ProteinName = RunChannel = PSM = NULL
  channel_len = data.table::uniqueN(annotation$Channel, na.rm = TRUE)
  input = input[order(Run, ProteinName, PSM, Channel), ]
  new_annotation = unique(input[, list(Run, ProteinName, Channel,
                                       RunChannel = paste(Run, Channel, sep = "_"))])
  summarized = input[,
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


#' Tukey median polish
#' @param intensities vector of log-intensities per protein and run
#' @param num_channels number of channels
#' @importFrom stats medpolish
#' @return numeric vector with length `num_channels`
#' @keywords internal
.medianPolish <- function(intensities, num_channels){
  wide = matrix(intensities, byrow = TRUE, ncol = num_channels)
  tmp_fit = stats::medpolish(wide, na.rm = TRUE, trace.iter = FALSE)
  tmp_fit$overall + tmp_fit$col
}
