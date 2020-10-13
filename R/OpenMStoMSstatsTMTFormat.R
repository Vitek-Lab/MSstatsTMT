#' Generate MSstatsTMT required input format for OpenMS output
#'
#' Convert OpenMS MSstatsTMT report into the required input format for MSstatsTMT.
#'
#' @export
#' @import tidyr
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom data.table as.data.table setkey rbindlist
#' @param input MSstatsTMT report from OpenMS
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withMissing_withinRun TRUE will remove PSM with any missing value within each Run. Defaut is FALSE.
#' @param rmPSM_withfewMea_withinRun only for rmPSM_withMissing_withinRun = FALSE. TRUE(default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultiplePSMs sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @return input for \code{\link{proteinSummarization}} function
#' @examples
#' head(raw.om)
#' input.om <- OpenMStoMSstatsTMTFormat(raw.om)
#' head(input.om)

OpenMStoMSstatsTMTFormat <- function(input,
                                     useUniquePeptide = TRUE,
                                     rmPSM_withMissing_withinRun = FALSE,
                                     rmPSM_withfewMea_withinRun = TRUE,
                                     rmProtein_with1Feature = FALSE,
                                     summaryforMultiplePSMs = sum){
  
  RetentionTime <- Reference <- ProteinName <- PeptideSequence <- Charge <- 
    Channel <- Run <- Intensity <- NULL
  
  Feature <- total_intensity <- nmea <- totalmea <- ID <- PSM <- NULL
  
  ################################################
  ## 0. check whether input has all the requried columns
  ################################################
  .check.input(input)
  
  ################################################
  ## 1. remove the psm with all zero intensities across all channels
  ## which means by row
  ################################################
  channels <- as.character(unique(input$Channel)) # record the channels
  
  # replace zero with NA
  input[input$Intensity == 0, "Intensity"] <- NA
  
  # replace NAs
  input <- input %>% dplyr::filter(!is.na(Intensity))
  
  # remove the psm with all NA intensities in one run
  tmp <- input %>% dplyr::group_by(ProteinName, PeptideSequence, Charge, Reference, RetentionTime, Run) %>% 
    summarise(total_intensity = sum(Intensity, na.rm = TRUE)) 
  
  input <- input %>% left_join(tmp) %>% 
    filter(total_intensity != 0) %>% 
    select(-total_intensity)
  
  message('** PSMs, that have all zero intensities across channels in each run, are removed.')
  
  ################################################
  ## 2. remove peptides which are used in more than one protein
  ## we assume to use unique peptide
  ################################################
  if (useUniquePeptide) {
    
    ## double check
    pepcount <- unique(input[, c("ProteinName", "PeptideSequence")])
    pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)
    
    ## count how many proteins are assigned for each peptide
    structure <- pepcount %>% group_by(PeptideSequence) %>% summarise(npro = n_distinct(ProteinName))
    remove_peptide <- structure[structure$npro > 1, ]
    
    ## remove the peptides which are used in more than one protein
    if (nrow(remove_peptide) != 0) {
      input <- input[-which(input$PeptideSequence %in% remove_peptide$PeptideSequence), ]
      
      message('** Peptides, that are used in more than one proteins, are removed.')
    }
    rm(pepcount)
    rm(structure)
    rm(remove_peptide)
  }
  
  ##############################
  ## 3. remove features which has missing measurements within each run
  ##############################
  if (rmPSM_withMissing_withinRun) {
    
    # remove the psm with all zero intensities in one run
    tmp <- input %>% dplyr::group_by(ProteinName, PeptideSequence, Charge, Reference, RetentionTime,  Run) %>% 
      summarise(nmea = sum(!is.na(Intensity))) 
    
    input <- input %>% left_join(tmp) %>% 
      filter(nmea == length(channels)) %>% 
      select(-nmea)
    
    message('** Rows which has any missing value within a run were removed from that run.')
  }
  
  ##############################
  ##  4. remove features which has 1 or 2 measurements across runs
  ##############################
  if (rmPSM_withfewMea_withinRun){
    
    # remove the psm with all zero intensities in one run
    tmp <- input %>% dplyr::group_by(ProteinName, PeptideSequence, Charge, Reference, RetentionTime, Run) %>% 
      summarise(nmea = sum(!is.na(Intensity))) 
    
    input <- input %>% left_join(tmp) %>% 
      filter(nmea > 2) %>% 
      select(-nmea)
    
    message(paste0('** ', sum(tmp$nmea <= 2), ' features have 1 or 2 intensities across runs are removed.'))
  }
  
  ##############################
  ## 5. remove multiple measurements per feature and run
  ##############################
  input <- .AggregatePSMstoPeptideIons(data = input, summaryforMultiplePSMs = summaryforMultiplePSMs)
  message('** PSMs have been aggregated to peptide ions.')
  
  ##############################
  ## 6. remove proteins with only one peptide and charge per protein
  ##############################
  if (rmProtein_with1Feature) {
    
    ## remove protein which has only one peptide
    tmp <- unique(input[, c("ProteinName", 'PSM')])
    tmp$ProteinName <- factor(tmp$ProteinName)
    count <- xtabs( ~ ProteinName, data = tmp)
    lengthtotalprotein <- length(count)
    
    removepro <- names(count[count <= 1])
    
    if (length(removepro) > 0) {
      input <- input[-which(input$ProteinName %in% removepro), ]
      message(paste0("** ", length(removepro),
                     ' proteins, which have only one feature in a protein, are removed among ',
                     lengthtotalprotein, ' proteins.'))
    }
  }
  
  ##############################
  ## 7. combine fractios within each mixture
  ##############################
  fractions <- unique(input$Fraction) # check the number of fractions in the input data
  if (length(fractions) > 1) { # combine fractions
    input <- .combine.fractions(input)
    ## change data.table to data.frame, in order to make the same class for input, without fraction
    input <- as.data.frame(input)
    message('** Fractions belonging to same mixture have been combined.')
  }
  
  input <- input[,c("ProteinName", "PeptideSequence", "Charge", "PSM",
                    "Mixture", "TechRepMixture", "Run",
                    "Channel", "Condition", "BioReplicate", "Intensity")]
  
  return(input)
}

## Check whether the input have all the requried columns
.check.input <- function(data){
  
  required.columns <- c("ProteinName", "PeptideSequence", "Charge", "Run", "TechRepMixture", "Fraction", 
                        "Mixture", "Channel", "Condition", "BioReplicate", "Intensity", "Reference",
                        "RetentionTime")
  
  if (!all(required.columns %in% colnames(data))) {
    
    missedColumns <- which(!(required.columns %in% colnames(data)))
    stop(paste("Please check the required column in the input file. ** columns :",
               paste(required.columns[missedColumns], collapse = ", "), " are missed."))
    
  }
}


## combine PSMs to peptide ions 
.AggregatePSMstoPeptideIons <- function(data, summaryforMultiplePSMs =  sum){
  
  RetentionTime <- Reference <- ProteinName <- PeptideSequence <- Charge <- 
    Channel <- Run <- Intensity <- NULL
  
  Feature <- total_intensity <- nmea <- totalmea <- ID <- PSM <- NULL
  
  if(!(length(summaryforMultiplePSMs) == 1 & (identical(summaryforMultiplePSMs, sum) | identical(summaryforMultiplePSMs, max)))){
    stop("summaryforMultiplePSMs can only be sum or max! ")
  }
  
  subject_info <- unique(data[,c("Mixture", "TechRepMixture", "Fraction", "Run",
                                 "Channel", "BioReplicate", "Condition")])
  subject_info$Channel <- as.character(subject_info$Channel)
  channels <- unique(subject_info$Channel)
  
  input <- data %>% dplyr::select(RetentionTime, Reference, ProteinName, 
                                  PeptideSequence, Charge, 
                                  Channel, Run, Intensity) %>% 
    tidyr::spread(Channel, Intensity)
  
  input$PSM <- paste(input$PeptideSequence, input$Charge, sep="_")
  input$Feature <- paste(input$ProteinName, input$PSM, input$Run, sep="_")
  input$ID <- 1:nrow(input)
  tmt_fea_summary <- input %>% 
    dplyr::group_by(Feature) %>%
    dplyr::summarise(nmea = n()) %>%
    dplyr::ungroup() 
  
  # peptide has single PSM
  tmt_single_PSM <- tmt_fea_summary %>%
    dplyr::filter(nmea == 1)
  keep_features <- as.vector(input[input$Feature %in% tmt_single_PSM$Feature, "ID"])
  
  # peptide has multiple PSMs
  tmt_multiple_PSMs <- tmt_fea_summary %>%
    dplyr::filter(nmea > 1)
  ## summarize multiple psms to single peptide ion
  if (nrow(tmt_multiple_PSMs) > 0) { ## if there is any feature issued.
    
    ## keep selected rows among issued rows
    for (i in 1:nrow(tmt_multiple_PSMs)) {
      if(i%%500 == 0){
        message("\t Summary for ", round(i/nrow(tmt_multiple_PSMs)*100, 2), "% peptides done.")
      }
      
      sub <- input[input$Feature == tmt_multiple_PSMs$Feature[i], ]
      sub <- unique(sub)
      
      if (nrow(sub) < 2) {
        keep_features <- c(keep_features, sub$ID)
        next()
      }
      
      ## decision 1 : first use the psm which has most number of measurement
      ## count the number of measurement per row
      sub$nmea <- apply(sub[, channels], 1, function(x) sum(!is.na(x)))
      sub2 <- sub[sub$nmea == max(sub$nmea), ] ## which.max choose only one row
      sub2 <- sub2 %>% dplyr::select(-nmea) # remove the added columns
      
      if (nrow(sub2) < 2) {
        keep_features <- c(keep_features, sub2$ID)
        next()
      } else {
        ## decision 2 : ## maximum abundances among intensities for identical features within one run
        sub2$totalmea <- apply(sub2[, channels], 1, function(x) summaryforMultiplePSMs(x, na.rm = TRUE))
        sub3 <- sub2[sub2$totalmea == max(sub2$totalmea), ]
        sub3 <- sub3 %>% dplyr::select(-totalmea)
        
        keep_features <- c(keep_features, sub3$ID)
        rm(sub3)
      }
      rm(sub2)
      rm(sub)
    }
  }
  # keep the best psm for each peptide
  input <- input %>% dplyr::filter(ID %in% keep_features) 
  # keep the required columns
  input <- input[,c("ProteinName", "PeptideSequence", "Charge", "PSM", "Run", channels)]
  # make long-format data
  tmt_data_by_peptide <- input %>% 
    tidyr::gather(Channel, Intensity, -ProteinName, -PeptideSequence, -Charge, -PSM, -Run) %>%
    dplyr::group_by(Run, Channel, ProteinName, PeptideSequence, Charge, PSM) %>%
    dplyr::summarize(Intensity = mean(Intensity, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% 
    filter(Intensity > 0) %>%
    dplyr::left_join(subject_info)
  
  return(tmt_data_by_peptide)
}
