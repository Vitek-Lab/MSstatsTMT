#' Generate MSstatsTMT required input format from Proteome discoverer output
#'
#' Convert Proteome discoverer output into the required input format for MSstatsTMT.
#'
#' @export
#' @import tidyr
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom data.table as.data.table setkey rbindlist
#' @param input data name of Proteome discover PSM output.
#' @param annotation data frame which contains column Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition. Refer to the example 'annotation.pd' for the meaning of each column.
#' @param which.proteinid Use 'Protein.Accessions'(default) column for protein name. 'Master.Protein.Accessions' can be used instead to get the protein name with single protein.
#' @param useNumProteinsColumn TURE(default) remove shared peptides by information of # Proteins column in PSM sheet.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withMissing_withinRun TRUE will remove PSM with any missing value within each Run. Defaut is FALSE.
#' @param rmPSM_withfewMea_withinRun only for rmPSM_withMissing_withinRun = FALSE. TRUE(default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @return input for \code{\link{proteinSummarization}} function
#' @examples
#' head(raw.pd)
#' head(annotation.pd)
#' input.pd <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd)
#' head(input.pd)

PDtoMSstatsTMTFormat <- function(input,
                                 annotation,
                                 which.proteinid = 'Protein.Accessions',
                                 useNumProteinsColumn = TRUE,
                                 useUniquePeptide = TRUE,
                                 rmPSM_withMissing_withinRun = FALSE,
                                 rmPSM_withfewMea_withinRun = TRUE,
                                 rmProtein_with1Feature = FALSE,
                                 summaryforMultipleRows = sum){

    ################################################
    ## 0. check input for annotation
    ################################################
    .check.annotation(annotation)

    if (!all(unique(annotation$Run) %in% unique(input$Spectrum.File))) {

        stop("Please check the annotation file. 'Run' must be matched with 'Spectrum.File'. ")
    }

    ################################################
    ## todo. check design of experiments
    ################################################

    ################################################
    ## 1. which protein id : Protein Accessions vs Master Protein Accesisions
    ################################################
    ## default : Protein Accessions
    which.pro <- NULL
    which.NumProteins <- NULL

    if (which.proteinid == 'Protein.Accessions') {
        which.pro <- 'Protein.Accessions'
    } else if (which.proteinid == 'Master.Protein.Accessions'){
        which.pro <- 'Master.Protein.Accessions'
    }

    if (is.null(which.pro)) {
        stop('** Please select which columns should be used for protein ids,
             among two options (Protein.Accessions, Master.Protein.Accessions).')
    }

    if (which.pro == 'Protein.Accessions' & !is.element('Protein.Accessions', colnames(input))) {

        which.pro <- 'Master.Protein.Accessions'
        message('** Use Master.Protein.Accessions instead of Protein.Accessions.')
    }

    if (which.pro == 'Master.Protein.Accessions' & !is.element('Master.Protein.Accessions', colnames(input))) {

        which.pro <- 'Protein.Accessions'
        message('** Use Protein.Accessions instead of Master.Protein.Accessions.')
    }

    if (!is.element(which.pro, colnames(input))) {
        stop('** Please select which columns should be used for protein ids,
             among two options (Protein.Accessions, Master.Protein.Accessions).')
    }

    # Find the corresponding number of proteins or protein groups for each peptide ions
    if (which.pro == 'Protein.Accessions') {
        which.NumProteins <- c('X..Proteins', '#.Proteins', 'Number.of.Proteins')
    } else if ( which.pro == 'Master.Protein.Accessions') {
        which.NumProteins <- c('X..Protein.Groups', '#.Protein.Groups', 'Number.of.Protein.Groups')
    }

    ################################################
    ## 2. get subset of columns
    ################################################
    # make sure the input is data frame format
    input <- as.data.frame(input)

    # For PD 2.2. It needs to update for new version of PD.
    check.column.start <- 'Abundance'
    check.column.end <- ''
    if(sum(startsWith(colnames(input), check.column.start) &
           endsWith(colnames(input), check.column.end)) == 0) {

        stop("There is no channel intensity column in the input data, which should start with 'Abundance'.")
    }

    # extract the columns for channels
    pattern <- paste(check.column.start, ".*", check.column.end, sep = "")
    channels <- grep(pattern, colnames(input), value = TRUE)

    input <- input[, which(colnames(input) %in% c(which.pro, which.NumProteins,
                                                  'Annotated.Sequence', 'Charge',
                                                  'Ions.Score', 'Spectrum.File', 'Quan.Info', 'Isolation.Interference....',
                                                  channels))]

    colnames(input)[colnames(input) == 'Master.Protein.Accessions'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'Protein.Accessions'] <- 'ProteinName'

    colnames(input)[colnames(input) %in% which.NumProteins] <- 'numProtein'
    colnames(input)[colnames(input) %in% which.NumProteins] <- 'numProtein'

    colnames(input)[colnames(input) == 'Annotated.Sequence'] <- 'PeptideSequence'
    colnames(input)[colnames(input) == 'Spectrum.File'] <- 'Run'

    # remove the rows which has missing values
    tmp <- input[, channels]
    missings <- apply(tmp, 1, function(x) sum(is.na(x)))
    input <- input[missings != length(channels), ]

    # remove the rows whose protein ID is empty
    input <- input[(input$ProteinName != "") & (!is.na(input$ProteinName)), ]

    ################################################
    ## 3. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################
    if (useNumProteinsColumn) {

      
        ## remove rows with #proteins is not 1
        input <- input[input$numProtein == '1', ]

        message('** Shared PSMs (assigned in multiple proteins) are removed.')

    }

    if (useUniquePeptide) {

        # make sure Quan.Info has 'unique' value
        if('Unique' %in% unique(input$Quan.Info)){
            input <- input[input$Quan.Info == 'Unique', ]
        }

        ## double check
        pepcount <- unique(input[, c("ProteinName", "PeptideSequence")])
        pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)

        ## count how many proteins are assigned for each peptide
        structure <- aggregate(ProteinName ~., data = pepcount, length)
        remove_peptide <- structure[structure$ProteinName > 1, ]

        ## remove the peptides which are used in more than one protein
        if (sum(remove_peptide$ProteinName != 1) != 0) {
            input <- input[-which(input$PeptideSequence %in% remove_peptide$PeptideSequence), ]

            message('** Peptides, that are used in more than one proteins, are removed.')
        }
        rm(pepcount)
        rm(structure)
        rm(remove_peptide)
    }

    ##############################
    ## 4. remove features which has missing measurements within each run
    ##############################
    if (rmPSM_withMissing_withinRun) {

        tmp <- input[, channels]
        nmea <- apply(tmp, 1, function(x) sum(!is.na(x)))
        input <- input[nmea == length(channels), ] # no missing values within run
        message('** Rows which has any missing value within a run were removed from that run.')
    }

    ##############################
    ##  5. remove features which has 1 or 2 measurements across runs
    ##############################
    if (rmPSM_withfewMea_withinRun){

        tmp <- input[, channels]
        nmea <- apply(tmp, 1, function(x) sum(!is.na(x)))
        input <- input[nmea > 2, ]
        message(paste0('** ', sum(nmea <= 2), ' features have 1 or 2 intensities across runs and are removed.'))
    }

    ##############################
    ## 6. remove multiple measurements per feature and run
    ##############################
    input$fea <- paste(input$PeptideSequence, input$Charge, sep="_")

    ## check multiple measurements
    input$fea2 <- paste(input$fea, input$ProteinName, sep="_")
    input$fea2 <- factor(input$fea2)

    count <- xtabs(~ fea2 + Run, input)
    ## there are multiple measurements
    count2 <- as.data.frame(count)
    fea.multimeas <- count2[count2$Freq > 1, ]

    ## separate input by multiple measurements vs one measurement
    if (nrow(fea.multimeas) > 0) { ## if there is any feature issued.
        fea.multimeas$issue <- paste(fea.multimeas$fea2, fea.multimeas$Run, sep="_")
        input$issue <- paste(input$fea2, input$Run, sep="_")

        ## keep rows with no issue
        input.no <- input[-which(input$issue %in% unique(fea.multimeas$issue)), ]

        ## keep selected rows among issued rows
        keepinfo.select <- NULL
        for (i in seq_along(unique(fea.multimeas$issue))) {
            # message("Row ", i)
            sub <- input[input$issue == unique(fea.multimeas$issue)[i], ]
            sub <- unique(sub)
            #subfea <- fea.multimeas[fea.multimeas$issue == unique(fea.multimeas$issue)[i], ]

            if (nrow(sub) < 2) {
                keepinfo.select <- rbind(keepinfo.select, sub)
                next()
            }

            ## decision 1 : first use the rows which has most number of measurement
            ## count the number of measurement per row
            sub$nmea <- apply(sub[, channels], 1, function(x) sum(!is.na(x)))
            sub2 <- sub[sub$nmea == max(sub$nmea), ] ## which.max choose only one row
            sub2 <- sub2[, -which(colnames(sub2) %in% c('nmea'))] # remove the added columns

            if (nrow(sub2) < 2) {
                keepinfo.select <- rbind(keepinfo.select, sub2)
                next()
            } else {
                ## decision 2 : keep the row with lowest isolation interference
                if("Isolation.Interference...." %in% names(sub2) & (sum(is.na(sub2$Isolation.Interference....)) == 0)){ # make sure Isolation.Interference.... is available
                    sub3 <- sub2[sub2$Isolation.Interference.... == min(sub2$Isolation.Interference....), ]
                } else {
                    sub3 <- sub2
                }

                if (nrow(sub3) < 2) {
                    keepinfo.select <- rbind(keepinfo.select, sub3)
                    next()
                } else {
                    ## decision 3 : keep the row with higher identification score
                    if("Ions.Score" %in% names(sub3) & (sum(is.na(sub3$Ions.Score)) == 0)){ # make sure Ions.Score is available
                        sub4 <- sub3[sub3$Ions.Score == max(sub3$Ions.Score), ] ## which.max choose only one row
                    } else {
                        sub4 <- sub3
                    }

                    if (nrow(sub4) < 2) {
                        keepinfo.select <- rbind(keepinfo.select, sub4)
                        next()
                    } else {
                        ## decision 4 : ## maximum or sum up abundances among intensities for identical features within one run
                        if(!(length(summaryforMultipleRows) == 1 & (identical(summaryforMultipleRows, sum) | identical(summaryforMultipleRows, max)))){
                            stop("summaryforMultipleRows can only be sum or max! ")
                        }

                        sub4$totalmea <- apply(sub4[, channels], 1, function(x) summaryforMultipleRows(x, na.rm = TRUE))
                        sub5 <- sub4[sub4$totalmea == max(sub4$totalmea), ]
                        sub5 <- sub5[, which(colnames(sub4) != "totalmea")]

                        if (nrow(sub5) < 2) {
                            keepinfo.select <- rbind(keepinfo.select, sub5)
                            next()
                        } else {
                            # sum up or maximum abundances among intensities for identical features within one run
                            if(identical(summaryforMultipleRows, sum)){
                                summaryforMultipleRows <- max
                            } else {
                                summaryforMultipleRows <- sum
                            }

                            sub4$totalmea <- apply(sub4[, channels], 1, function(x) summaryforMultipleRows(x, na.rm = TRUE))
                            sub5 <- sub4[sub4$totalmea == max(sub4$totalmea), ]
                            sub5 <- sub5[, which(colnames(sub4) != "totalmea")]
                            keepinfo.select <- rbind(keepinfo.select, sub5)
                            
                        }
                        rm(sub5)
                    }
                    rm(sub4)
                }
                rm(sub3)
            }
            rm(sub2)
            rm(sub)
        }

        input.new <- rbind(input.no, keepinfo.select)
        input.new <- input.new[, -which(colnames(input.new) %in%
                                            c('Quan.Info', 'numProtein', 'Ions.Score', 'fea', 'fea2',
                                              'issue', 'Isolation.Interference....'))]

        message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')

    } else {
        input.new <- input[, -which(colnames(input) %in% c('Quan.Info', 'numProtein', 'Ions.Score', 'fea', 'fea2',
                                                           'issue', 'Isolation.Interference....'))]
    }

    # make long format
    input.long <- melt(input.new, id = c('ProteinName',
                                       'PeptideSequence',
                                       'Charge',
                                       'Run'),
                       variable.name = "Channel",
                       value.name = "Intensity")

    # make sure no dupliate rows
    # input.long <- input.long[!is.na(input.long$Intensity), ]
    input <- input.long
    rm(input.long)
    rm(input.new)

    ##############################
    ## 10. add annotation
    ##############################
    # match the channels from input with that in annotation file
    input$Channel <- gsub(check.column.start, "", input$Channel)
    input$Channel <- gsub("\\.", "", input$Channel)
    input$Channel <- gsub(check.column.end, "", input$Channel)

    if (!all(unique(annotation$Channel) %in% unique(input$Channel))) {

        stop("Please check the annotation file. The channel name must be matched with that in input data ")
    }

    input <- merge(input, annotation, by = c("Run", "Channel"), all.x = TRUE)

    ## check whether there is any missing 'Condition'
    noruninfo <- unique(input[is.na(input$Condition) & !is.na(input$Intensity), c("Run", "Channel")])

    if (nrow(noruninfo) > 0) {
        for(i in seq_len(nrow(noruninfo))){
            message( paste0('** Annotation for Run : ', noruninfo[i, "Run"],
                            ", Channel : ", noruninfo[i, "Channel"], " are missed.") )
        }
        stop('** Please add them to annotation file. If the channal doesn\'t have sample, please add \'Empty\'.')
    }

    input.final <- data.frame("ProteinName" = input$ProteinName,
                              "PeptideSequence" = input$PeptideSequence,
                              "Charge" = input$Charge,
                              "PSM" = paste(input$PeptideSequence, input$Charge, sep = "_"),
                              "Channel" = as.factor(input$Channel),
                              "Condition" = as.factor(input$Condition),
                              "BioReplicate" = as.factor(input$BioReplicate),
                              "Mixture" = as.factor(input$Mixture),
                              "TechRepMixture" = as.factor(input$TechRepMixture),
                              "Fraction" = as.factor(input$Fraction),
                              "Run" = as.factor(input$Run),
                              "Intensity" = input$Intensity)

    input <- input.final
    rm(input.final)

    ##############################
    ## 10. remove proteins with only one peptide and charge per protein
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
    ## 11. combine fractios within each mixture
    ##############################
    fractions <- unique(annotation$Fraction) # check the number of fractions in the input data
    if (length(fractions) > 1) { # combine fractions
        input <- .combine.fractions(input)
        ## change data.table to data.frame, in order to make the same class for input, without fraction
        input <- as.data.frame(input)
        message('** Fractions belonging to same mixture have been combined.')
    }

    input <- input[,c("ProteinName", "PeptideSequence", "Charge", "PSM",
                      "Mixture", "TechRepMixture", "Run",
                      "Channel", "Condition", "BioReplicate", "Intensity")]
    
    ## finally, make sure no duplicate rows
    input <- unique(input)
    
    return(input)
}

## Remove the peptide ions overlapped among multiple fractions of same biological mixture
## data: peptide level data, which has columns Protein, PSM, Channel, Condition, BioReplicate,
## Mixture, TechRepMixture, Fraction, Run, Intensity
.combine.fractions <- function(data){

    Mixture <- techrun <- ProteinName <- PeptideSequence <- Charge <- PSM <- TechRepMixture <- Channel <-
    Condition <- BioReplicate <- Intensity <- fea <-  Run <- . <- NULL

    # combine fractions for each technical replcate and each mixture
    data$techrun <- paste(data$Mixture, data$TechRepMixture, sep = "_")
    techruns <- unique(data$techrun)
    data <- as.data.table(data) # make sure data format
    data$Run <- as.character(data$Run)

    all.data <- list()
    for (i in seq_along(techruns)) { # for each technical replicate
        sub_data <- data[techrun == techruns[i]]
        sub_data <- sub_data[!is.na(Intensity)] # remove the NAs
        sub_data$fea <- paste(sub_data$PSM, sub_data$ProteinName, sep = "_")
        sub_data$fea <- factor(sub_data$fea) # make the feature
        sub_data$id <- paste(sub_data$fea, sub_data$Run, sep = "_")

        ## count how many fractions are assigned for each peptide ion
        structure <- aggregate(Run ~ . , data = unique(sub_data[, .(fea, Run)]), length)
        ## 1. first, keep features which are measured in one fraction
        remove_peptide_ion <- structure[structure$Run > 1, ]

        ## 2. second, if features are measured in multiple fractionations,
        ## use the fractions with maximum average.
        ## remove_peptide_ion : features that are measured in multiple fractions
        if (nrow(remove_peptide_ion) > 0) {
            # select the rows for the features that are measured in multiple fractions
            overlapped.feas <- sub_data %>% dplyr::filter(fea %in% remove_peptide_ion$fea)

            # keep the fractions with maximum average reporter ion abundance
            mean.frac.feature <- overlapped.feas %>% 
              dplyr::group_by(fea, id) %>% 
              dplyr::summarise(mean = mean(Intensity, na.rm = TRUE))
            
            remove.fraction <- mean.frac.feature %>% 
              dplyr::group_by(fea) %>% 
              dplyr::filter(mean != max(mean))
            
            # filter out the unused fractions
            sub_data <- sub_data %>% dplyr::filter(!id %in% remove.fraction$id)
            
            structure_2 <- aggregate(Run ~ . , 
                                     data = unique(as.data.table(sub_data)[,.(fea, Run)]), 
                                     length)
            
            ## check if there are features which have same mean intensities across runs
            remove_peptide_ion_2 <- structure_2[structure_2$Run > 1, ]
            
            if(nrow(remove_peptide_ion_2) > 0) {
              
              # select the rows for the features that are measured in multiple fractions
              overlapped.feas.2 <- sub_data %>% dplyr::filter(fea %in% remove_peptide_ion_2$fea)
              
              # keep the fractions with maximum summation reporter ion abundance
              sum.frac.feature <- overlapped.feas.2 %>% 
                dplyr::group_by(fea, id) %>% 
                dplyr::summarise(sum = sum(Intensity, na.rm = TRUE))
              
              remove.fraction.2 <- sum.frac.feature %>% 
                dplyr::group_by(fea) %>% 
                dplyr::filter(sum != max(sum))
              
              # filter out the unused fractions
              sub_data <- sub_data %>% dplyr::filter(!id %in% remove.fraction.2$id)
              
              structure_3 <- aggregate(Run ~ . , 
                                       data = unique(as.data.table(sub_data)[,.(fea, Run)]), 
                                       length)
              
              ## check if there are features which have same summed intensities across runs
              remove_peptide_ion_3 <- structure_3[structure_3$Run > 1, ]
              
              if(nrow(remove_peptide_ion_3) > 0) {
                
                # select the rows for the features that are measured in multiple fractions
                overlapped.feas.3 <- sub_data %>% dplyr::filter(fea %in% remove_peptide_ion_3$fea)
                
                # keep the fractions with maximum summation reporter ion abundance
                max.frac.feature <- overlapped.feas.3 %>% 
                  dplyr::group_by(fea, id) %>% 
                  dplyr::summarise(max = max(Intensity, na.rm = TRUE))
                
                remove.fraction.3 <- max.frac.feature %>% 
                  dplyr::group_by(fea) %>% 
                  dplyr::filter(max != max(max))
                
                # filter out the unused fractions
                sub_data <- sub_data %>% dplyr::filter(!id %in% remove.fraction.3$id)
                
                rm(max.frac.feature)
                rm(remove.fraction.3)
                rm(overlapped.feas.3)
              }
              
              rm(sum.frac.feature)
              rm(remove.fraction.2)
              rm(overlapped.feas.2)
              
            }
            
            rm(mean.frac.feature)
            rm(remove.fraction)
            rm(overlapped.feas)
            
            message('** For peptides overlapped between fractions of ', techruns[i],
                    ', use the fraction with maximal average abundance.')
        } 
        
        # Remove unnecessary columns
        sub_data <- sub_data %>%
          dplyr::select(ProteinName, PeptideSequence, Charge, PSM,
                   Mixture, TechRepMixture, Run,
                   Channel, Condition, BioReplicate, Intensity)

        all.data[[i]] <- as.data.table(sub_data)
    }

    data.shared.pep.rm <- rbindlist(all.data)
    # The fractions have been combined
    data.shared.pep.rm$Run <- paste(data.shared.pep.rm$Mixture, data.shared.pep.rm$TechRepMixture, sep = "_")
    return(data.shared.pep.rm)
}

## Check whether the annotation file matches with the input
.check.annotation <- function(annotation){

    required.annotation <- c("Run", "TechRepMixture", "Fraction", "Mixture", "Channel", "Condition", "BioReplicate")

    if (!all(required.annotation %in% colnames(annotation))) {

        missedAnnotation <- which(!(required.annotation %in% colnames(annotation)))
        stop(paste("Please check the required column in the annotation file. ** columns :",
                   paste(required.annotation[missedAnnotation], collapse = ", "), " are missed."))

    }
}
