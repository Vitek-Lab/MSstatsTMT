#' Generate MSstatsTMT required input format for SpectroMine output
#'
#' Convert SpectroMine output into the required input format for MSstatsTMT.
#'
#' @export
#' @import tidyr
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom data.table as.data.table setkey rbindlist
#' @param input data name of SpectroMine PSM output. Read PSM sheet.
#' @param annotation data frame which contains column Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition. Refer to the example 'annotation.mine' for the meaning of each column.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with NA and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withMissing_withinRun TRUE will remove PSM with any missing value within each Run. Defaut is FALSE.
#' @param rmPSM_withfewMea_withinRun only for rmPSM_withMissing_withinRun = FALSE. TRUE(default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @return input for \code{\link{proteinSummarization}} function
#' @examples
#' head(raw.mine)
#' head(annotation.mine)
#' input.mine <- SpectroMinetoMSstatsTMTFormat(raw.mine, annotation.mine)
#' head(input.mine)

SpectroMinetoMSstatsTMTFormat <- function(input,
                                          annotation,
                                          filter_with_Qvalue = TRUE,
                                          qvalue_cutoff = 0.01,
                                          useUniquePeptide = TRUE,
                                          rmPSM_withMissing_withinRun = FALSE,
                                          rmPSM_withfewMea_withinRun = TRUE,
                                          rmProtein_with1Feature = FALSE,
                                          summaryforMultipleRows = sum){

    PeptideSequence <- ProteinName <- Channel <- NULL
    ################################################
    ## 0. check input for annotation
    ################################################
    .check.annotation(annotation)

    if (!all(unique(annotation$Run) %in% unique(input$R.FileName))) {

        stop("Please check the annotation file. 'Run' must be matched with 'R.FileName'. ")
    }

    ################################################
    ## todo. check design of experiments
    ################################################

    ################################################
    ## 1. get subset of columns
    ################################################
    # make sure the input is data frame format
    input <- as.data.frame(input)

    # For SpectroMine. It needs to update for new version of SpectroMine.
    check.column.start <- 'PSM.'
    check.column.end <- '..Raw.'
    if(sum(startsWith(colnames(input), check.column.start) &
           endsWith(colnames(input), check.column.end)) == 0) {

        stop("There is no channel intensity column in the input data, which should start with 'PSM' and end with 'Raw'.")
    }

    # extract the columns for channels
    pattern <- paste(check.column.start, ".*", check.column.end, sep = "")
    channels <- grep(pattern, colnames(input), value = TRUE)

    # extract the required columns
    required.cols <- c('PG.ProteinAccessions', 'P.MoleculeID', 'PP.Charge',
                       'PG.QValue', 'PSM.Qvalue', 'R.FileName', channels)

    if (!all( required.cols %in% colnames(input))) {

        missing.col <- required.cols[!required.cols %in% colnames(input)]

        stop(paste0("** Please check the required input. The required input needs : ",
                    toString(missing.col)))
    }

    input <- input[, which(colnames(input) %in% required.cols)]
    colnames(input)[colnames(input) == 'PG.ProteinAccessions'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'PP.Charge'] <- 'Charge'
    colnames(input)[colnames(input) == 'P.MoleculeID'] <- 'PeptideSequence'
    colnames(input)[colnames(input) == 'R.FileName'] <- 'Run'
    colnames(input)[colnames(input) == 'PSM.Qvalue'] <- 'Qvalue'

    # remove the rows whose protein ID is empty
    input <- input[(input$ProteinName != "") & (!is.na(input$ProteinName)), ]

    ##############################
    ## 2. filter by Qvalue
    ##############################
    ## protein FDR
    input[!is.na(input$PG.QValue) & input$PG.QValue > 0.01, channels] <- NA

    message('** Intensities with great than 0.01 in PG.QValue are replaced with NA.')
    input <- input[, -which(colnames(input) %in% 'PG.QValue')]

    ## PSM qvalue
    if (filter_with_Qvalue) {
        ## when qvalue > qvalue_cutoff, replace with zero for intensity
        input[!is.na(input$Qvalue) & input$Qvalue > qvalue_cutoff, channels] <- NA

        message(paste0('** Intensities with great than ', qvalue_cutoff, ' in EG.Qvalue are replaced with NA.'))
        input <- input[, -which(colnames(input) %in% 'Qvalue')]
    }

    ##############################
    ## 3. remove featuares with all na or zero
    ## some rows have all zero values across all MS runs. They should be removed.
    ##############################
    tmp <- input[,channels]
    nmea <- apply(tmp, 1, function(x) sum(!is.na(x)))
    input <- input[nmea > 0, ]
    message(paste0('** ',
                   sum(nmea == 0), ' rows have all NAs are removed.'))

    ################################################
    ## 4. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################
    if (useUniquePeptide) {

        pepcount <- unique(input[, c("ProteinName", "PeptideSequence")]) ## Protein.group.IDs or Sequence
        pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)

        ## count how many proteins are assigned for each peptide
        structure <- pepcount %>% group_by(PeptideSequence) %>% summarise(length = n_distinct(ProteinName))
        remove_peptide <- structure[structure$length > 1, ]

        ## remove the peptides which are used in more than one protein
        if(nrow(remove_peptide) != 0){
            input <- input[-which(input$PeptideSequence %in% remove_peptide$PeptideSequence), ]

            message('** Peptides, that are used in more than one proteins, are removed.')
        } else {
            message('** All peptides are unique peptides in proteins.')
        }
        rm(pepcount)
        rm(structure)
        rm(remove_peptide)
    }

    ##############################
    ## 5. remove features which has missing measurements within each run
    ##############################
    if (rmPSM_withMissing_withinRun) {

        tmp <- input[,channels]
        nmea <- apply(tmp, 1, function(x) sum(!is.na(x)))
        input <- input[nmea == length(channels), ] # no missing values within run
        message('** Rows which has any missing value within a run were removed from that run.')
    }

    ##############################
    ##  6. remove features which has 1 or 2 measurements across runs
    ##############################
    if (rmPSM_withfewMea_withinRun){

        tmp <- input[,channels]
        nmea <- apply(tmp, 1, function(x) sum(!is.na(x)))
        input <- input[nmea > 2, ]
        message(paste0('** ', sum(nmea <= 2), ' features have 1 or 2 intensities across runs and are removed.'))
    }

    ##############################
    ## 7. remove multiple measurements per feature and run
    ##############################
    input$fea <- paste(input$PeptideSequence, input$Charge, sep = "_")

    ## check multiple measurements
    input$fea2 <- paste(input$fea, input$ProteinName, sep = "_")
    input$fea2 <- factor(input$fea2)

    count <- xtabs(~ fea2 + Run, input)
    ## there are multiple measurements
    count2 <- as.data.frame(count)
    fea.multimeas <- count2[count2$Freq > 1, ]

    ## separate input by multiple measurements vs one measurement
    if (nrow(fea.multimeas) > 0) { ## if there is any feature issued.
        fea.multimeas$issue <- paste(fea.multimeas$fea2, fea.multimeas$Run, sep = "_")
        input$issue <- paste(input$fea2, input$Run, sep = "_")
        keepinfo.select <- NULL # store the final data

        ## keep rows with no issue
        input.no <- input[-which(input$issue %in% unique(fea.multimeas$issue)), ]

        ## keep selected rows among issued rows
        for (i in seq_along(unique(fea.multimeas$issue))) {
            # message("Row ", i)
            sub <- input[input$issue == unique(fea.multimeas$issue)[i], ]
            sub <- unique(sub)

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

            } else {## decision 2 : ## maximum or sum up abundances among intensities for identical features within one run

                if(!(length(summaryforMultipleRows) == 1 & (identical(summaryforMultipleRows, sum) | identical(summaryforMultipleRows, max)))){
                    stop("summaryforMultipleRows can only be sum or max! ")
                }

                sub2$totalmea <- apply(sub2[, channels], 1, function(x) summaryforMultipleRows(x, na.rm = TRUE))
                sub3 <- sub2[sub2$totalmea == max(sub2$totalmea), ]
                sub3 <- sub3[, which(colnames(sub2) != "totalmea")]

                if (nrow(sub3) < 2) {
                    keepinfo.select <- rbind(keepinfo.select, sub3)

                } else {
                    # sum up or maximum abundances among intensities for identical features within one run
                    if(identical(summaryforMultipleRows, sum)){
                        temp_summaryforMultipleRows <- max
                    } else {
                        temp_summaryforMultipleRows <- sum
                    }

                    sub2$totalmea <- apply(sub2[, channels], 1, function(x) temp_summaryforMultipleRows(x, na.rm = TRUE))
                    sub3 <- sub2[sub2$totalmea == max(sub2$totalmea), ]
                    sub3 <- sub3[, which(colnames(sub2) != "totalmea")]
                    keepinfo.select <- rbind(keepinfo.select, sub3)
                }
                rm(sub3)
            }
            rm(sub2)
            rm(sub)
        }

        input.new <- rbind(input.no, keepinfo.select)
        input.new <- input.new[, -which(colnames(input.new) %in%
                                            c('fea', 'fea2', 'issue'))]
        message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')

    } else {
        input.new <- input[, -which(colnames(input) %in% c('fea', 'fea2', 'issue'))]
    }

    # make long format
    input.long <- melt(input.new, id = c('ProteinName',
                                       'PeptideSequence','Charge',
                                       'Run'),
                       variable.name = "Channel",
                       value.name = "Intensity")

    # make sure no dupliate rows
    # input.long <- input.long[!is.na(input.long$Intensity), ]
    input <- input.long
    rm(input.long)
    rm(input.new)

    ##############################
    ## 8. merge annotation
    ##############################
    # match the channels from input with that in annotation file
    input$Channel <- gsub(check.column.start, "", input$Channel)
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
    ## 9. remove proteins with only one peptide and charge per protein
    ##############################
    if (rmProtein_with1Feature) {

        ## remove protein which has only one peptide
        tmp <- unique(input[!is.na(input$Intensity), c("ProteinName", 'PSM')])
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
    ## 10. combine fractions within each mixture
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
                      "Channel", "BioReplicate", "Condition", "Intensity")]
    
    ## finally, make sure no duplicate rows
    input <- unique(input)
    
    return(input)
}
