#' Normalization for TMT data
#' @param input data.table
#' @param type "peptides" for peptide normalization between channel and run,
#' "proteins" for protein normalization
#' @param normalize logical, if TRUE, data will be normalized
#' 
#' @return data.table
#' 
#' @export
#' 
MSstatsNormalizeTMT = function(input, type, normalize) {
  if (type == "peptides") {
    .normalizePeptides(input, normalize)
  } else {
    .normalizeProteins(input, normalize)
  }
}


#' Normalization between channels (before summarization)
#' @param input data.table
#' @param normalize logical, if TRUE, `input` data will be normalized
#' @return data.table
#' @keywords internal
.normalizePeptides = function(input, normalize) {
  log2Intensity = Intensity = Run = Channel = NULL
  
  if (normalize) {
    input[, MedianLog2Int := median(log2Intensity, na.rm = TRUE),
          by = c("Run", "Channel")]
    median_baseline = median(
      unique(input[, list(Run, Channel, MedianLog2Int)])[, MedianLog2Int],
      na.rm = TRUE
    )
    input[, log2Intensity := log2Intensity + median_baseline - MedianLog2Int]
    input[, Intensity := 2 ^ log2Intensity]
  }
  input[, !(colnames(input) == "MedianLog2Int"), with = FALSE]
}


#' Normalization between MS runs (after protein summarization)
#' @param input data.table
#' @param normalize logical, if TRUE, data will be normalized
#' @return data.table
#' @keywords internal
.normalizeProteins = function(input, normalize) {
  Abundance = NormalizationAbundance = Run = Condition = MedianNormalized = NULL
  Mixture = TechRepMixture = Channel = Protein = BioReplicate = Condition = NULL
  
  if (normalize) {
    group_info = unique(input$Condition)
    if (is.element("Norm", group_info)) {
      input[!is.na(Abundance),
            NumRuns := data.table::uniqueN(Run, na.rm = TRUE),
            by = "Protein"]
      input[!is.na(Abundance),
            NumRunsWithNorm := .countRunsWithNorm(Run, Condition),
            by = "Protein"]
      input[!is.na(Abundance),
            NormalizationAbundance := .getNormalizationAbundance(Abundance, 
                                                                 Condition),
            by = c("Protein", "Run")]
      input[!is.na(Abundance), 
            MedianNormalized := .getRunsMedian(.SD),
            by = "Protein", .SDcols = c("Run", "NormalizationAbundance")]
      input[!is.na(Abundance),
            NormalizedAbundance := Abundance + MedianNormalized -
              NormalizationAbundance]
      input[,
            Abundance := ifelse(NumRuns > 1 & NumRunsWithNorm > 1,
                                NormalizedAbundance, Abundance)]
    } else {
      msg = paste("** 'Norm' information in Condition",
                  "is required for normalization.",
                  "Please check it. At this moment,", 
                  "normalization is not performed.")
      getOption("MSstatsLog")("INFO", msg)
      getOption("MSstatsMsg")("INFO", msg)
    }
  }
  input[, list(Mixture, TechRepMixture, Run, Channel, 
               Protein, Abundance, BioReplicate, Condition)]
}


#' Utility function: count runs with "Norm" channel
#' @param run vector of run labels
#' @param condition vector of condition labels
#' @return integer
#' @keywords internal
.countRunsWithNorm = function(run, condition) {
  data.table::uniqueN(run[condition == "Norm"], na.rm = TRUE)
}


#' Utility function: get mean abundance for "Norm" channels
#' @param abundance vector of abundances
#' @param condition vector of condition labels
#' @return numeric
#' @keywords internal
.getNormalizationAbundance = function(abundance, condition) {
  data = abundance[condition == "Norm"]
  if (all(is.na(data))) {
    return(NA)
  } else {
    return(mean(data, na.rm = TRUE))
  }
}


#' Utility function: get median from unique values per run
#' @param input data.table / list
#' @return numeric
#' @keywords internal
.getRunsMedian = function(input) {
  median(unique(input)$NormalizationAbundance, na.rm = TRUE)
}
