#' Generate MSstatsTMT required input format from MaxQuant output
#'
#' Convert MaxQuant output into the required input format for MSstatsTMT.
#'
#' @export
#' @import tidyr
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom data.table as.data.table setkey rbindlist
#' @param evidence name of 'evidence.txt' data, which includes feature-level data.
#' @param proteinGroups name of 'proteinGroups.txt' data.
#' @param annotation data frame which contains column Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition. Refer to the example 'annotation.mq' for the meaning of each column.
#' @param which.proteinid Use 'Proteins'(default) column for protein name. 'Leading.proteins' or 'Leading.razor.proteins' or 'Gene.names' can be used instead to get the protein ID with single protein. However, those can potentially have the shared peptides.
#' @param rmProt_Only.identified.by.site TRUE will remove proteins with '+' in 'Only.identified.by.site' column from proteinGroups.txt, which was identified only by a modification site. FALSE is the default.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withMissing_withinRun TRUE will remove PSM with any missing value within each Run. Defaut is FALSE.
#' @param rmPSM_withfewMea_withinRun only for rmPSM_withMissing_withinRun = FALSE. TRUE(default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @return input for \code{\link{proteinSummarization}} function
#' @examples
#' head(evidence)
#' head(proteinGroups)
#' head(annotation.mq)
#' input.mq <- MaxQtoMSstatsTMTFormat(evidence, proteinGroups, annotation.mq)
#' head(input.mq)

MaxQtoMSstatsTMTFormat <- function(evidence,
                                 proteinGroups,
                                 annotation,
                                 which.proteinid = 'Proteins',
                                 rmProt_Only.identified.by.site = FALSE,
                                 useUniquePeptide = TRUE,
                                 rmPSM_withMissing_withinRun = FALSE,
                                 rmPSM_withfewMea_withinRun = TRUE,
                                 rmProtein_with1Feature = FALSE,
                                 summaryforMultipleRows = sum){

    PeptideSequence <- ProteinName <- NULL
     
    ## evidence.txt file
    input <- evidence

    ################################################
    ## 0. check input for annotation
    ################################################
    .check.annotation(annotation)

    if (!all(unique(annotation$Run) %in% unique(input$Raw.file))) {

        stop("Please check the annotation file. 'Run' must be matched with 'Raw.file'. ")
    }

    ################################################
    ## todo. check design of experiments
    ################################################

    ################################################
    ## 1. remove contaminant, reverse proteinID
    ## Contaminant, Reverse column in evidence
    ################################################
    ## remove contaminant
    if (is.element("Contaminant", colnames(input)) &
        is.element("+", unique(input$Contaminant))) {
        input <- input[-which(input$Contaminant %in% "+"), ]
    }

    if (is.element("Potential.contaminant", colnames(input)) &
        is.element("+", unique(input$Potential.contaminant))) {
        input <- input[-which(input$Potential.contaminant %in% "+"), ]
    }

    ## remove reversed protein id
    if (is.element("Reverse", colnames(input)) &
        is.element("+", unique(input$Reverse))) {
        input <- input[-which(input$Reverse %in% "+"), ]
    }

    message('** + Contaminant, + Reverse, + Only.identified.by.site, proteins are removed.')

    ## ? Only.identified.by.site column in proteinGroupID? : sometimes, it is not in evidence.txt
    if (rmProt_Only.identified.by.site) {
        if (is.element("Only.identified.by.site", colnames(input)) &
            is.element("+", unique(proteinGroups$Only.identified.by.site))) {
            input <- input[-which(input$Only.identified.by.site %in% "+"), ]

            message('** + Only.identified.by.site are removed.')
        }
    }

    ################################################
    ## 2. matching proteinGroupID protein list
    ################################################
    ## need to check proteinGroupID in evidence and proteinGroup.txt the same
    ## 'id' in proteinGroups.txt vs 'Protein.group.IDs' in input
    ## possible to have some combination in Protein.group.IDs in input,
    ## such as 64;1274;1155;1273 instead of 64, 1274.. separately.
    ## however, 'id' in proteinGroups.txt has only one protein id.
    ## Use the subset of input, which has the same 'id' from proteinGroups.
    ## (combination of some ids seems not to be used for intensity)

    ## first, remove contaminants
    if (is.element("Contaminant", colnames(proteinGroups)) &
        is.element("+", unique(proteinGroups$Contaminant))) {
        proteinGroups <- proteinGroups[-which(proteinGroups$Contaminant %in% "+"), ]
    }

    if (is.element("Potential.contaminant", colnames(proteinGroups)) &
        is.element("+", unique(proteinGroups$Potential.contaminant))) {
        proteinGroups <- proteinGroups[-which(proteinGroups$Potential.contaminant %in% "+"), ]
    }

    if (is.element("Reverse", colnames(proteinGroups)) &
        is.element("+", unique(proteinGroups$Reverse))) {
        proteinGroups <- proteinGroups[-which(proteinGroups$Reverse %in% "+"), ]
    }

    ## ? Only.identified.by.site column in proteinGroupID? : sometimes, it is not in evidence.txt
    if (rmProt_Only.identified.by.site) {
        if (is.element("Only.identified.by.site", colnames(proteinGroups)) &
            is.element("+", unique(proteinGroups$Only.identified.by.site))) {
            proteinGroups <- proteinGroups[-which(proteinGroups$Only.identified.by.site %in% "+"), ]
        }
    }

    ## then take proteins which are included
    input <- input[which(input$Protein.group.IDs %in% unique(proteinGroups$id)), ]

    ## then use 'protein.IDs' in proteinGroups.txt
    ## because if two 'proteins' in evidence.txt are used in one protein ID, need to use certain protein name in input.
    ## for example, protein.IDs in proteinGroups.txt are P05204;O00479.
    ## but, two 'proteins in evidence.txt, such as P05204;O00479, and P05204.

    tempname <- unique(proteinGroups[,c("Protein.IDs", "id")])
    colnames(tempname) <- c("uniquefromProteinGroups", "Protein.group.IDs")

    input <- merge(input, tempname, by <- "Protein.group.IDs")

    ################################################
    ## 3. which protein id :
    ## 1) Proteins
    ## 2) Leading.proteins
    ## 3) Leading.razor.protein
    ## 4) Gene.names
    ################################################
    ## default : Protein
    which.pro <- NULL

    if (which.proteinid == 'Proteins') {
        which.pro <- 'Proteins'
    } else if (which.proteinid == 'Leading.proteins') {
        which.pro <- 'Leading.proteins'
    } else if (which.proteinid == 'Leading.razor.protein') {
        which.pro <- 'Leading.razor.protein'
    } else if (which.proteinid == 'Gene.names') {
        which.pro <- 'Gene.names'
    }

    if (which.pro == 'Gene.names' & !is.element('Gene.names', colnames(input))) {

        which.pro <- 'Leading.razor.protein'
        message('** Use Leading.razor.protein instead of Gene.names.')
    }

    if (which.pro == 'Leading.razor.protein' & !is.element('Leading.razor.protein', colnames(input))) {

        which.pro <- 'Leading.proteins'
        message('** Use Leading.proteins instead of Leading.razor.protein.')
    }

    if (which.pro == 'Leading.proteins' & !is.element('Leading.proteins', colnames(input))) {

        which.pro <- 'Proteins'
        message('** Use Proteins instead of Leading.proteins.')
    }
    ## at least 'Proteins' should be in the evidence.txt

    if (!is.element(which.pro, colnames(input))) {
        stop('** Please select which columns should be used for protein ids,
             among four options (Proteins, Leading.proteins, Leading.razor.protein, Gene.names).')
    }

    ################################################
    ## 4. get subset of columns
    ################################################
    # For maxquant 1.5. It needs to update for new version of maxquant
    inputlevel <- 'Reporter.intensity.corrected'

    ## columns which include inputlevel
    channels <- colnames(input)[grep(inputlevel, colnames(input))]
    input<-as.data.frame(input)
    input <- input[, which(colnames(input) %in% c(which.pro,
                                                  'Modified.sequence', 'Charge',
                                                  'Raw.file',
                                                  'Score',
                                                  channels))]

    colnames(input)[colnames(input) == 'Proteins'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'Leading.proteins'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'Leading.razor.protein'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'Gene.names'] <- 'ProteinName'

    colnames(input)[colnames(input) == 'Modified.sequence'] <- 'PeptideSequence'
    colnames(input)[colnames(input) == 'Raw.file'] <- 'Run'

    ## remove the underbar(_) in peptide sequence
    input$PeptideSequence <- gsub('_', '', input$PeptideSequence)

    ################################################
    ## 5. remove the psm with all zero intensities across all channels
    ## which means by row
    ################################################

    tmp <- input[, channels]
    tmp.count <- apply(tmp, 1, function(x) sum(x == 0, na.rm = TRUE))
    nchannel <- ncol(tmp)

    ## remove row with all zero
    input <- input[tmp.count < nchannel, ]

    message('** PSMs, that have all zero intensities across channels in each run, are removed.')

    # replace zero with NA
    input.int <- input[, channels]
    input.int[input.int == 0] <- NA
    input[, channels] <- input.int

    rm(input.int)

    ################################################
    ## 6. remove peptides which are used in more than one protein
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
    ## 7. remove features which has missing measurements within each run
    ##############################
    if (rmPSM_withMissing_withinRun) {

        tmp <- input[,channels]
        nmea <- apply(tmp, 1, function(x) sum(!is.na(x)))
        input <- input[nmea == length(channels), ] # no missing values within run
        message('** Rows which has any missing value within a run were removed from that run.')
    }

    ##############################
    ##  8. remove features which has 1 or 2 measurements across runs
    ##############################
    if (rmPSM_withfewMea_withinRun){

        tmp <- input[,channels]
        nmea <- apply(tmp, 1, function(x) sum(!is.na(x)))
        input <- input[nmea > 2, ]
        message(paste0('** ', sum(nmea <= 2), ' features have 1 or 2 intensities across runs and are removed.'))
    }

    ##############################
    ## 9. remove multiple measurements per feature and run
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

        ## keep rows with no issue
        input.no <- input[-which(input$issue %in% unique(fea.multimeas$issue)), ]

        ## keep selected rows among issued rows
        keepinfo.select <- NULL
        for (i in 1:length(unique(fea.multimeas$issue))) {
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
            } else {
                ## decision 2 : keep the row with higher Score
                ## Score : Andromeda score for the best associated MS/MS spectrum.
                if("Score" %in% names(sub2) & (sum(is.na(sub2$Score)) == 0)){ # make sure Score is available
                    sub3 <- sub2[sub2$Score == max(sub2$Score), ]
                } else {
                    sub3 <- sub2
                }

                if (nrow(sub3) < 2) {
                    keepinfo.select <- rbind(keepinfo.select, sub3)
                    next()
                } else {
                    ## decision 3 : ## maximum or sum up abundances among intensities for identical features within one run
                    if(!(length(summaryforMultipleRows) == 1 & (identical(summaryforMultipleRows, sum) | identical(summaryforMultipleRows, max)))){
                        stop("summaryforMultipleRows can only be sum or max! ")
                    }

                    sub3$totalmea <- apply(sub3[, channels], 1, function(x) summaryforMultipleRows(x, na.rm = TRUE))
                    sub4 <- sub3[sub3$totalmea == max(sub3$totalmea), ]
                    sub4 <- sub4[, which(colnames(sub3) != "totalmea")]

                    if (nrow(sub4) < 2) {
                        keepinfo.select <- rbind(keepinfo.select, sub4)
                        next()
                    } else {
                        # sum up or maximum abundances among intensities for identical features within one run
                        if(identical(summaryforMultipleRows, sum)){
                            summaryforMultipleRows <- max
                        } else {
                            summaryforMultipleRows <- sum
                        }

                        sub3$totalmea <- apply(sub3[, channels], 1, function(x) summaryforMultipleRows(x, na.rm = TRUE))
                        sub4 <- sub3[sub3$totalmea == max(sub3$totalmea), ]
                        sub4 <- sub4[, which(colnames(sub3) != "totalmea")]
                        keepinfo.select <- rbind(keepinfo.select, sub4)
                    }
                    rm(sub4)
                }
                rm(sub3)
            }
            rm(sub2)
            rm(sub)
        }

        ## combine the featurew without any issue
        input.new <- rbind(input.no, keepinfo.select)
        input.new <- input.new[, -which(colnames(input.new) %in%
                                            c('Score', 'fea', 'fea2', 'issue'))]

        message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')

    } else {
        input.new <- input[, -which(colnames(input) %in% c('Score', 'fea', 'fea2'))]
    }

    # make long format
    input.long <- melt(input.new, id=c('ProteinName',
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
    ## channel prefix for channel
    input$Channel <- gsub(inputlevel, 'channel', input$Channel)

    input <- merge(input, annotation, by = c("Run", "Channel"), all.x = TRUE)

    ## check whether there is any missing 'Condition'
    noruninfo <- unique(input[is.na(input$Condition) & !is.na(input$Intensity), c("Run", "Channel")])

    if (nrow(noruninfo) > 0) {
        for(i in 1:nrow(noruninfo)){
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
                      "Channel", "BioReplicate", "Condition", "Intensity")]
    return(input)
}
