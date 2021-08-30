#' A dummy function to store shared documentation items.
#' 
#' @import data.table
#' @importFrom MSstatsConvert MSstatsImport MSstatsClean MSstatsPreprocess 
#' MSstatsBalancedDesign MSstatsMakeAnnotation MSstatsSaveSessionInfo
#' MSstatsLogsSettings
#' 
#' @param fewMeasurements 'remove'(default) will remove the features that have 1 or 2 measurements across runs.
#' @param useUniquePeptide TRUE (default) removes peptides that are assigned for more than one proteins. 
#' We assume to use unique peptide for each protein.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have only 1 feature, which is the combination of peptide, precursor charge, fragment and charge. FALSE is default.
#' @param removeProtein_with1Peptide TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.
#' @param removeOxidationMpeptides TRUE will remove the peptides including 'oxidation (M)' in modification. FALSE is default.
#' @param removeMpeptides TRUE will remove the peptides including 'M' sequence. FALSE is default.
#' @param use_log_file logical. If TRUE, information about data processing
#' will be saved to a file.
#' @param append logical. If TRUE, information about data processing will be added
#' to an existing log file.
#' @param verbose logical. If TRUE, information about data processing wil be printed
#' to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. 
#' If not provided, such a file will be created automatically.
#' If `append = TRUE`, has to be a valid path to a file.
#' 
#' @return NULL.
#' @keywords internal
.documentFunction = function(fewMeasurements, 
                             useUniquePeptide,
                             summaryforMultipleRows, 
                             removeProtein_with1Feature,
                             removeProtein_with1Protein,
                             removeOxidationMpeptides,
                             removeMpeptides) {
  NULL
}


#' Generate MSstatsTMT required input format from MaxQuant output
#' 
#' @param evidence name of 'evidence.txt' data, which includes feature-level data.
#' @param proteinGroups name of 'proteinGroups.txt' data.
#' @param annotation data frame which contains column Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition. Refer to the example 'annotation.mq' for the meaning of each column.
#' @param which.proteinid Use 'Proteins' (default) column for protein name. 'Leading.proteins' or 'Leading.razor.proteins' or 'Gene.names' can be used instead to get the protein ID with single protein. However, those can potentially have the shared peptides.
#' @param rmProt_Only.identified.by.site TRUE will remove proteins with '+' in 'Only.identified.by.site' column from proteinGroups.txt, which was identified only by a modification site. FALSE is the default.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .documentFunction
#' 
#' @return data.frame of class "MSstatsTMT"
#' 
#' @export
#' 
#' @examples
#' head(evidence)
#' head(proteinGroups)
#' head(annotation.mq)
#' input.mq <- MaxQtoMSstatsTMTFormat(evidence, proteinGroups, annotation.mq)
#' head(input.mq)
#' 
MaxQtoMSstatsTMTFormat = function(
  evidence, proteinGroups, annotation, which.proteinid = 'Proteins',
  rmProt_Only.identified.by.site = FALSE, useUniquePeptide = TRUE,
  rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE, 
  summaryforMultipleRows = sum, 
  use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
  ...
) {
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path, 
                                      base = "MSstatsTMT_converter_log_")

  input = MSstatsConvert::MSstatsImport(list(evidence = evidence,
                                             protein_groups = proteinGroups), 
                                        "MSstatsTMT", "MaxQuant", ...)
  input = MSstatsConvert::MSstatsClean(
    input,
    protein_id_col = which.proteinid, 
    remove_by_site = rmProt_Only.identified.by.site,
    channel_columns = "Reporterintensitycorrected")
  annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
  
  feature_columns = c("PeptideSequence", "PrecursorCharge")
  input = MSstatsConvert::MSstatsPreprocess(
    input, 
    annotation, 
    feature_columns,
    remove_shared_peptides = useUniquePeptide,
    remove_single_feature_proteins = rmProtein_with1Feature,
    feature_cleaning = list(remove_features_with_few_measurements = rmPSM_withfewMea_withinRun,
                            summarize_multiple_psms = summaryforMultipleRows))
  input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                fix_missing = "zero_to_na")
  data.table::setnames(input, "PrecursorCharge", "Charge", skip_absent = TRUE)
  
  msg_final = paste("** Finished preprocessing. The dataset is ready",
                    "to be processed by the proteinSummarization function.")
  getOption("MSstatsLog")("INFO", msg_final)
  getOption("MSstatsMsg")("INFO", msg_final)
  getOption("MSstatsLog")("INFO", "\n")
  input
}


#' Generate MSstatsTMT required input format for OpenMS output
#' @param input MSstatsTMT report from OpenMS
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultiplePSMs sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .documentFunction
#'  
#' @return `data.frame` of class `MSstatsTMT`.
#' 
#' @export
#' 
#' @examples
#' head(raw.om)
#' input.om <- OpenMStoMSstatsTMTFormat(raw.om)
#' head(input.om)
#' 
OpenMStoMSstatsTMTFormat = function(
  input, useUniquePeptide = TRUE, rmPSM_withfewMea_withinRun = TRUE, 
  rmProtein_with1Feature = FALSE, summaryforMultiplePSMs = sum, 
  use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL, ...
) {
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path, 
                                      base = "MSstatsTMT_converter_log_")

  input = MSstatsConvert::MSstatsImport(list(input = input), 
                                        "MSstatsTMT", "OpenMS", ...)
  input = MSstatsConvert::MSstatsClean(input)
  
  feature_columns = c("PeptideSequence", "PrecursorCharge")
  input = MSstatsConvert::MSstatsPreprocess(
    input, 
    NULL, 
    feature_columns,
    remove_shared_peptides = useUniquePeptide,
    remove_single_feature_proteins = rmProtein_with1Feature,
    feature_cleaning = list(remove_features_with_few_measurements = rmPSM_withfewMea_withinRun,
                            summarize_multiple_psms = summaryforMultiplePSMs)
  )
  input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                fix_missing = "zero_to_na")
  
  data.table::setnames(input, "PrecursorCharge", "Charge", skip_absent = TRUE)
  
  msg_final = paste("** Finished preprocessing. The dataset is ready",
                    "to be processed by the proteinSummarization function.")
  getOption("MSstatsLog")("INFO", msg_final)
  getOption("MSstatsMsg")("INFO", msg_final)
  getOption("MSstatsLog")("INFO", "\n")
  input
}


#' Convert Proteome Discoverer output to MSstatsTMT format.
#' 
#' @param input PD report or a path to it.
#' @param annotation annotation with Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition columns or a path to file. Refer to the example 'annotation' for the meaning of each column.
#' @param which.proteinid Use 'Protein.Accessions'(default) column for protein name. 'Master.Protein.Accessions' can be used instead to get the protein name with single protein.
#' @param useNumProteinsColumn logical, TURE(default) remove shared peptides by information of # Proteins column in PSM sheet.
#' @param useUniquePeptide logical, if TRUE (default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum (default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .documentFunction
#' 
#' @return `data.frame` of class `MSstatsTMT`
#' 
#' @export
#' 
#' @examples
#' head(raw.pd)
#' head(annotation.pd)
#' input.pd <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd)
#' head(input.pd)
#' 
PDtoMSstatsTMTFormat <- function(
  input, annotation, which.proteinid = 'Protein.Accessions', 
  useNumProteinsColumn = TRUE, useUniquePeptide = TRUE, 
  rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE, 
  summaryforMultipleRows = sum, 
  use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL, ...
) {
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path, 
                                      base = "MSstatsTMT_converter_log_")

  input = MSstatsConvert::MSstatsImport(list(input = input),
                                        "MSstatsTMT", "ProteomeDiscoverer", ...)
  input = MSstatsConvert::MSstatsClean(
    input, 
    protein_id_column = which.proteinid,
    remove_shared = useNumProteinsColumn,
    remove_protein_groups = useNumProteinsColumn)
  annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
  
  feature_columns = c("PeptideSequence", "PrecursorCharge")
  input = MSstatsConvert::MSstatsPreprocess(
    input,
    annotation, 
    feature_columns = c("PeptideSequence", "PrecursorCharge"),
    remove_shared_peptides = useUniquePeptide,
    remove_single_feature_proteins = rmProtein_with1Feature,
    feature_cleaning = list(remove_features_with_few_measurements = rmPSM_withfewMea_withinRun,
                            summarize_multiple_psms = summaryforMultipleRows)
  )
  input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                fix_missing = "zero_to_na")
  data.table::setnames(input, "PrecursorCharge", "Charge", skip_absent = TRUE)
  
  msg_final = paste("** Finished preprocessing. The dataset is ready",
                    "to be processed by the proteinSummarization function.")
  getOption("MSstatsLog")("INFO", msg_final)
  getOption("MSstatsMsg")("INFO", msg_final)
  getOption("MSstatsLog")("INFO", "\n")
  input
}


#' Import data from SpectroMine
#'
#' @param input data name of SpectroMine PSM output. Read PSM sheet.
#' @param annotation data frame which contains column Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition. Refer to the example 'annotation.mine' for the meaning of each column.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with NA and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .documentFunction
#'  
#' @return `data.frame` of class `MSstatsTMT`
#' 
#' @export 
#' @examples
#' head(raw.mine)
#' head(annotation.mine)
#' input.mine <- SpectroMinetoMSstatsTMTFormat(raw.mine, annotation.mine)
#' head(input.mine)
#' 
SpectroMinetoMSstatsTMTFormat <- function(
  input, annotation, filter_with_Qvalue = TRUE, qvalue_cutoff = 0.01,
  useUniquePeptide = TRUE, rmPSM_withfewMea_withinRun = TRUE, 
  rmProtein_with1Feature = FALSE, summaryforMultipleRows = sum,
  use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL, ...
) {
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path, 
                                      base = "MSstatsTMT_converter_log_")
  
  input = MSstatsConvert::MSstatsImport(list(input = input), 
                                        "MSstatsTMT", "SpectroMine", ...)
  input = MSstatsConvert::MSstatsClean(input)
  annotation = MSstatsMakeAnnotation(input, annotation)
  
  pq_filter = list(score_column = "PGQValue", 
                   score_threshold = 0.01, 
                   direction = "smaller",
                   behavior = "fill", 
                   handle_na = "keep", 
                   fill_value = NA,
                   filter = TRUE, 
                   drop_column = TRUE)
  
  qval_filter = list(score_column = "Qvalue", 
                     score_threshold = qvalue_cutoff, 
                     direction = "smaller",
                     behavior = "fill", 
                     handle_na = "keep", 
                     fill_value = NA,
                     filter = filter_with_Qvalue, 
                     drop_column = TRUE)
  
  feature_columns = c("PeptideSequence", "PrecursorCharge") 
  input = MSstatsConvert::MSstatsPreprocess(
    input, 
    annotation, 
    feature_columns,
    remove_shared_peptides = useUniquePeptide,
    remove_single_feature_proteins = rmProtein_with1Feature,
    score_filtering = list(pgq = pq_filter, psm_q = qval_filter),
    feature_cleaning = list(remove_features_with_few_measurements = rmPSM_withfewMea_withinRun,
                            summarize_multiple_psms = summaryforMultipleRows)
  )
  input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                fix_missing = "zero_to_na")
  data.table::setnames(input, "PrecursorCharge", "Charge", skip_absent = TRUE)
  
  msg_final = paste("** Finished preprocessing. The dataset is ready",
                    "to be processed by the proteinSummarization function.")
  getOption("MSstatsLog")("INFO", msg_final)
  getOption("MSstatsMsg")("INFO", msg_final)
  getOption("MSstatsLog")("INFO", "\n")
  input
}
