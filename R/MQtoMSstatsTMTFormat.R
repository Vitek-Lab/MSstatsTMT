#' Generate MSstatsTMT required input format for MaxQuant output
#'
#' Convert MaxQuant output into the required input format for MSstatsTMT.
#'
#' @export
#' @importFrom dplyr summarise n
#' @param evidence name of 'evidence.txt' data, which includes feature-level data.
#' @param proteinGroups name of 'proteinGroups.txt' data.
#' @param annotation data frame which contains column Run, Channel, Condition, BioReplicate, Mixture.
#' @param fraction indicates whether the data has fractions. If there are fractions, then overlapped peptide ions will be removed and then fractions are combined for each mixture.
#' @param remove.Only.identified.by.site TRUE will remove proteins with '+' in 'Only.identified.by.site' column from proteinGroups.txt, which was identified only by a modification site. FALSE is the default.
#' @param which.proteinid Use 'Proteins'(default) column for protein name. 'Leading.proteins' or 'Leading.razor.proteins' can be used instead. However, those can potentially have the shared peptides.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.
#' @param removePSM_withMissingValue_withinRun TRUE(default) will remove PSM with any missing value within each Run.
#' @param removeProtein_with1Feature TRUE(default) will remove the proteins which have only 1 peptide and charge.
#' @return input for protein.summarization function
#' @examples
#' \dontrun{
#' head(evidence)
#' head(proteinGroups)
#' head(annotation)
#' required.input <- MQtoMSstatsTMTFormat(evidence, proteinGroups, annotation)
#' head(required.input)
#' }

MQtoMSstatsTMTFormat <- function(evidence,
                                proteinGroups,
                                annotation,
                                fraction = FALSE,
                                remove.Only.identified.by.site = FALSE,
                                which.proteinid = 'Proteins',
                                useUniquePeptide = TRUE,
                                summaryforMultipleRows = sum,
                                removePSM_withMissingValue_withinRun = TRUE,
                                removeProtein_with1Feature = FALSE){

    PeptideSequence = fea2 = Run = NULL
    ## evidence.txt file
    input <- evidence

    ################################################
    ## 0. check input for annotation
    ################################################
    #required.annotation <- c("Run", "Channel", "Group", "BiologicalMixture", "Subject")
    required.annotation <- c("Run", "Channel", "Condition", "BioReplicate", "Mixture")

    if (!all(required.annotation %in% colnames(annotation))) {

        missedAnnotation <- which(!(required.annotation %in% colnames(annotation)))
        stop(paste("Please check the required column in the annotation file. ** columns :",
                   paste(required.annotation[missedAnnotation], collapse = ", "), " are missed."))

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
    if (remove.Only.identified.by.site) {
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
    if (remove.Only.identified.by.site) {
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

    input <- merge(input, tempname, by="Protein.group.IDs")

    ################################################
    ## 3. which protein id :
    ## 1) Proteins
    ## 2) Leading.proteins
    ## 3) Leading.razor.protein
    ################################################
    ## default : Protein Accessions
    which.pro <- NULL

    if (which.proteinid == 'Proteins') {
        which.pro <- 'Proteins'
    } else if (which.proteinid == 'Leading.proteins') {
        which.pro <- 'Leading.proteins'
    } else if (which.proteinid == 'Leading.razor.protein') {
        which.pro <- 'Leading.razor.protein'
    }

    if (is.null(which.pro)) {
        stop('** Please select which columns should be used for protein ids,
             among three options (Proteins, Leading.proteins, Leading.razor.protein).')
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
             among three options (Proteins, Leading.proteins, Leading.razor.protein).')
    }

    ################################################
    ## 4. get subset of columns
    ################################################

    inputlevel <- 'Reporter.intensity.corrected'
    ## columns which include inputlevel
    int.column <- colnames(input)[grep(inputlevel, colnames(input))]

    input <- input[, which(colnames(input) %in% c(which.pro,
                                                'Modified.sequence', 'Charge',
                                                'Raw.file',
                                                'Score',
                                                int.column))]

    colnames(input)[colnames(input) == 'Proteins'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'Leading.proteins'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'Leading.razor.protein'] <- 'ProteinName'

    colnames(input)[colnames(input) == 'Modified.sequence'] <- 'PeptideSequence'
    colnames(input)[colnames(input) == 'Raw.file'] <- 'Run'

    ## remove the underbar(_) in peptide sequence
    input$PeptideSequence <- gsub('_', '', input$PeptideSequence)

    ################################################
    ## 5. remove the psm with all zero intensities across all channels
    ## which means by row
    ################################################

    tmp <- input[, which(colnames(input) %in% int.column)]
    tmp.count <- apply(tmp, 1, function(x) sum(x == 0))
    nchannel <- ncol(tmp)

    ## remove row with all zero
    input <- input[tmp.count < nchannel, ]

    message('** PSMs, that have all zero intensities across channels in each run, are removed.')

    ################################################
    ## 6. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################

    if (useUniquePeptide) {

        ## double check
        pepcount <- unique(input[, c("ProteinName", "PeptideSequence")])
        pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)

        ## count how many proteins are assigned for each peptide
        structure <- pepcount %>% group_by(PeptideSequence) %>% summarise(npro=n())
        remove_peptide <- structure[structure$npro > 1, ]

        ## remove the peptides which are used in more than one protein
        if (nrow(remove_peptide) != 0) {
            input <- input[-which(input$PeptideSequence %in% remove_peptide$PeptideSequence), ]

            message('** Peptides, that are used in more than one proteins, are removed.')
        }
        rm(structure)
    }

    ##############################
    ## 7. remove multiple measurements per feature and run
    ##############################

    input$fea <- paste(input$PeptideSequence, input$Charge, sep="_")

    ## check multiple measurements
    input$fea2 <- paste(input$fea, input$ProteinName, sep="_")
    input$fea2 <- factor(input$fea2)

    structure <- input %>% group_by(fea2, Run) %>% summarise(nmeasure=n())
    ## nmeasurement should be 1, otherwise, there are multiple measurements
    fea.multimeas <- structure[structure$nmeasure > 1, ]

    ## if there is any feature issued.
    ## separate input by multiple measurements vs one measurement
    if (nrow(fea.multimeas) > 0) {
        ## make unique feature and run
        fea.multimeas$issue <- paste(fea.multimeas$fea2, fea.multimeas$Run, sep="_")
        input$issue <- paste(input$fea2, input$Run, sep="_")

        ## keep rows with no issue
        input.no <- input[-which(input$issue %in% unique(fea.multimeas$issue)), ]

        ## keep selected rows among issued rows
        keepinfo.select <- NULL

        ## each feature
        for (i in 1:length(unique(fea.multimeas$fea2))) {
            sub <- input[input$fea2 == unique(fea.multimeas$fea2)[i], ]
            subfea <- fea.multimeas[fea.multimeas$fea2 == unique(fea.multimeas$fea2)[i], ]

            ## each run within specific feature
            for (j in 1:length(unique(subfea$Run))) {
                subsub <- sub[sub$Run == unique(subfea$Run)[j], ]

                if (nrow(subsub) < 2) {
                    next()
                }

                ## decision 1 : first use the rows which has most number of measurement
                ## count the number of measurement per row
                ## not NA also not zero
                subsub$nmea <- apply(subsub[, int.column], 1, function(x) sum(!is.na(x) & x > 0))
                ## which.max choose only one row
                subsub2 <- subsub[subsub$nmea == max(subsub$nmea), ]

                if (nrow(subsub2) < 2) {
                    ## new column within if should be removed. otherwise can't be added in data.frame.
                    subsub2 <- subsub2[, -which(colnames(subsub2) %in% c("nmea"))]
                    keepinfo.select <- rbind(keepinfo.select,
                                            subsub2)
                } else {
                    ## decision 2 : keep the row with higher Score
                    ## Score : Andromeda score for the best associated MS/MS spectrum.
                    subsub3 <- subsub2[subsub2$Score == max(subsub2$Score), ] ## which.max choose only one row

                    if (nrow(subsub3) < 2) {
                        ## new column within if should be removed. otherwise can't be added in data.frame.
                        subsub3 <- subsub3[, -which(colnames(subsub3) %in% c("nmea"))]
                        keepinfo.select <- rbind(keepinfo.select, subsub3)
                    } else {
                        ## decision 3 : maximum or sum up abundances among intensities for identical features within one run
                        subsub3$totalmea <- apply(subsub3[, int.column], 1, function(x) summaryforMultipleRows(x, na.rm = TRUE))
                        subsub4 <- subsub3[subsub3$totalmea == max(subsub3$totalmea), ]

                        ## new column within if should be removed. otherwise can't be added in data.frame.
                        subsub4 <- subsub4[, -which(colnames(subsub4) %in% c("nmea", "totalmea"))]
                        keepinfo.select <- rbind(keepinfo.select, subsub4)
                        rm(subsub4)
                    }
                    rm(subsub3)
                }
                rm(subsub2)
            }
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

    input.long <- input.long[!is.na(input.long$Intensity), ]
    input <- input.long
    rm(input.long)
    rm(input.new)


    ## channel prefix for channel
    input$Channel <- gsub(inputlevel, 'channel', input$Channel)

    ##############################
    ## 8. add annotation
    ##############################

    input <- merge(input, annotation, by=c("Run", "Channel"), all.x =TRUE)

    ## check whether there is any missing 'Condition'
    noruninfo <- unique(input[is.na(input$Condition), c("Run", "Channel")])

    if (nrow(noruninfo) > 0) {
        for(i in 1:nrow(noruninfo)){
            message( paste0('** Annotation for Run : ', noruninfo[i, "Run"],
                            ", Channel : ", noruninfo[i, "Channel"], " are missed.") )
        }
        stop('** Please add them to annotation file.')
    }

    input.final <- data.frame("ProteinName" = input$ProteinName,
                            "PeptideSequence" = input$PeptideSequence,
                            "Charge" = input$Charge,
                            "PSM" = paste(input$PeptideSequence, input$Charge, sep="_"),
                            "Channel" = as.factor(input$Channel),
                            "Condition" = input$Condition,
                            "BioReplicate" = input$BioReplicate,
                            "Run" = input$Run,
                            "Mixture" = input$Mixture,
                            "Intensity" = input$Intensity)

    input <- input.final
    rm(input.final)

    ##############################
    ## 6. remove features which has missing measurements within each run
    ##############################
    ## number of channels in the dataset
    n_channels <- length(unique(input$Channel))

    if (removePSM_withMissingValue_withinRun) {

        ## it is the same across experiments. # measurement per feature.
        xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
        xtmp$eachRun <- paste(xtmp$PSM, xtmp$Run, sep="_")
        count_measure <- xtabs( ~ eachRun, xtmp)
        remove_feature_name <- count_measure[count_measure < n_channels]

        if (length(remove_feature_name) > 0) {
            xtmp <- xtmp[-which(xtmp$eachRun %in% names(remove_feature_name)), ]
        }
        input <- xtmp[, colnames(xtmp) != "eachRun"]
        message('** Features which has any missing value within a run were removed from that run.')
    }

    ##############################
    ## 7. remove proteins with only one peptide and charge per protein
    ##############################

    if (removeProtein_with1Feature) {

        ## remove protein which has only one peptide
        tmp <- unique(input[, c("ProteinName", 'PSM')])
        tmp$ProteinName <- factor(tmp$ProteinName)
        count <- xtabs( ~ ProteinName, data=tmp)
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
    ## 8. combine fractions within each mixture
    ##############################

    if (fraction) {
      input <- combine.fractions(input)
      ## change data.table to data.frame, in order to make the same class for input, without fraction
      input <- as.data.frame(input)
      message('** Fractions belonging to same mixture have been combined.')
    }

    return(input)
}
