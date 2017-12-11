#' Converter for Proteome discoverer output
#' @importFrom reshape2 melt
#' @export
#' @param input output from Proteome discoverer : PSM sheet
#' @param annotation data frame which contains column Run, Channel, Group, BiologicalMixture. BiologicalMixture indicates which runs are techinical runs or fractions of the same biological mixture.
#' @param Fraction Indicates whether the data has fractions. If there are fractions, then overlapped peptide ions will be removed and then fractions are combined for each biological mixture.
PDtoMSstatsTMTFormat <- function(input,
                                 annotation,
                                 Fraction = FALSE,
                                 useNumProteinsColumn=TRUE,
                                 useUniquePeptide=TRUE,
                                 summaryforMultipleRows=max,
                                 removePSM_withMissingValue_withinRun = TRUE,
                                 removeProtein_with1Peptide=FALSE,
                                 which.proteinid = 'Protein.Accessions'){


    ################################################
    ## 1. which protein id : Protein Accessions vs Master Protein Accesisions
    ################################################
    ## default : Protein Accessions
    which.pro <- NULL
    which.NumProteins <- NULL

    if(which.proteinid == 'Protein.Accessions'){
        which.pro <- 'Protein.Accessions'
    } else if (which.proteinid == 'Master.Protein.Accessions'){
        which.pro <- 'Master.Protein.Accessions'
    }


    if( is.null(which.pro) ){
        stop('** Please select which columns should be used for protein ids, among two options (Protein.Accessions, Master.Protein.Accessions).')
    }

    if(which.pro == 'Protein.Accessions' & !is.element('Protein.Accessions', colnames(input))){

        which.pro <- 'Master.Protein.Accessions'
        message('** Use Master.Protein.Accessions instead of Protein.Accessions.')
    }

    if(which.pro == 'Master.Protein.Accessions' & !is.element('Master.Protein.Accessions', colnames(input))){

        which.pro <- 'Protein.Accessions'
        message('** Use Protein.Accessions instead of Master.Protein.Accessions.')
    }

    if(!is.element(which.pro, colnames(input))){
        stop('** Please select which columns should be used for protein ids, among two options (Protein.Accessions, Master.Protein.Accessions).')
    }

    # Find the corresponding number of proteins or protein groups for each peptide ions
    if(which.pro == 'Protein.Accessions'){
        which.NumProteins <- 'X..Proteins'
    } else if (which.pro == 'Master.Protein.Accessions'){
        which.NumProteins <- 'X..Protein.Groups'
    }

    ################################################
    ## 2. get subset of columns
    ################################################

    channels <- as.character(unique(annotation$Channel))
    input <- input[, which(colnames(input) %in% c(which.pro, which.NumProteins, 'Annotated.Sequence', 'Charge',
                                                  'Ions.Score', 'Spectrum.File', 'Quan.Info',
                                                  channels))]

    colnames(input)[colnames(input) == 'Protein.Group.Accessions'] <- 'ProteinName'
    colnames(input)[colnames(input) == 'Protein.Accessions'] <- 'ProteinName'

    colnames(input)[colnames(input) == 'X..Proteins'] <- 'numProtein'
    colnames(input)[colnames(input) == 'X..Protein.Groups'] <- 'numProtein'

    colnames(input)[colnames(input) == 'Annotated.Sequence'] <- 'PeptideSequence'
    colnames(input)[colnames(input) == 'Spectrum.File'] <- 'Run'


    ################################################
    ## 3. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    ################################################

    if( useNumProteinsColumn ){

        ## remove rows with #proteins is not 1
        input <- input[input$numProtein == '1', ]

        message('** Rows with #Proteins, which are not equal to 1, are removed.')

    }

    if( useUniquePeptide ){

        input <- input[input$Quan.Info == 'Unique', ]

        ## double check
        pepcount <- unique(input[, c("ProteinName", "PeptideSequence")])
        pepcount$PeptideSequence <- factor(pepcount$PeptideSequence)

        ## count how many proteins are assigned for each peptide
        structure <- aggregate(ProteinName ~., data=pepcount, length)
        remove_peptide <- structure[structure$ProteinName != 1, ]

        ## remove the peptides which are used in more than one protein
        if( sum(remove_peptide$ProteinName != 1) != 0 ){
            input <- input[-which(input$PeptideSequence %in% remove_peptide$PeptideSequence), ]

            message('** Peptides, that are used in more than one proteins, are removed.')

        }
    }

    ##############################
    ### 4. remove multiple measurements per feature and run
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
    fea.multimeas$issue <- paste(fea.multimeas$fea2, fea.multimeas$Run, sep="_")
    input$issue <- paste(input$fea2, input$Run, sep="_")

    ## keep rows with no issue
    input.no <- input[-which(input$issue %in% unique(fea.multimeas$issue)), ]

    ## keep selected rows among issued rows
    keepinfo.select <- NULL
    for(i in 1:length(unique(fea.multimeas$fea2))){
        sub <- input[input$fea2 == unique(fea.multimeas$fea2)[i], ]
        subfea <- fea.multimeas[fea.multimeas$fea2 == unique(fea.multimeas$fea2)[i], ]

        for(j in 1:length(unique(subfea$Run))){
            subsub <- sub[sub$Run == unique(subfea$Run)[j], ]
            subsub$multi <- seq(1, nrow(subsub))

            if( nrow(subsub) < 2 ) next()
            ## decision1 : first use the rows which has most number of measurement
            ## count the number of measurement per row
            subsub$nmean <- apply(subsub[, channels], 1, function(x) sum(!is.na(x)))
            subsub2 <- subsub[subsub$nmea == max(subsub$nmea), ] ## which.max choose only one row

            if( nrow(subsub2) < 2 ){
                keepinfo.select <- rbind(keepinfo.select,
                                         subsub2)
            } else {
                ## decision2 : keep the row with higher identification score
                subsub3 <- subsub2[subsub2$Ions.Score == max(subsub2$Ions.Score), ] ## which.max choose only one row
                if(nrow(subsub3) < 2){
                    keepinfo.select <- rbind(keepinfo.select, subsub3)
                } else{
                    ## decision3 : ## maximum or sum up abundances among intensities for identical features within one run
                    subsub3$totalmea <- apply(subsub3[, channels], 1, function(x) summaryforMultipleRows(x, na.rm = T))
                    subsub4 <- subsub3[subsub3$totalmea == max(subsub3$totalmea), ]
                    subsub4 <- subsub4[, which(colnames(keepinfo.select) != "totalmea")]
                    keepinfo.select <- rbind(keepinfo.select, subsub4)
                    rm(subsub4)
                }
                rm(subsub3)

            }
            rm(subsub2)
            # print(paste('j=', j, sep=""))
        }
        print(paste('i=', i, sep=""))
    }
    keepinfo.select <- keepinfo.select[, -which(colnames(keepinfo.select) %in% c('multi', 'nmean'))]
    input.new <- rbind(input.no, keepinfo.select)

    input.new <- input.new[, -which(colnames(input.new) %in% c('Quan.Info', 'numProtein', 'Ions.Score', 'fea', 'fea2', 'issue'))]

    # make long format
    input.long <- reshape2::melt(input.new, id=c('ProteinName',
                                       'PeptideSequence','Charge',
                                       'Run'),
                       variable.name = "Channel",
                       value.name = "Intensity")

    message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')

    input <- input.long
    rm(input.long)

    ##############################
    ### 5. add annotation
    ##############################

    input <- merge(input, annotation, by=c("Run","Channel"), all.x =TRUE)
    noruninfo <- unique(input[is.na(input$Group), c("Run", "Channel")])

    if ( nrow(noruninfo) > 0 ) {
        for(i in 1:nrow(noruninfo)){
            message( paste('** Annotation for Run ', noruninfo[i, "Run"], " Channel ", noruninfo[i, "Channel"], ' are needed)' ))
        }
        stop('Please add them to annotation file.')
    }

    input.final <- data.frame(Protein = input$ProteinName,
                              PSM = paste(input$PeptideSequence, input$Charge, sep="_"),
                              log2Intensity = log2(input$Intensity),
                              Channel = as.factor(input$Channel),
                              Subject = paste(input$Run, input$Channel, sep="."),
                              Run = input$Run,
                              BiologicalMixture = input$BiologicalMixture,
                              Group = input$Group)

    #if( any(is.element(colnames(input), 'Fraction')) ) {
    #input.final <- data.frame(input.final,
    #Fraction = input$Fraction)
    #}

    input <- input.final
    rm(input.final)

    ##############################
    ###  6. remove features which has missing measurements within each run
    ##############################
    if(removePSM_withMissingValue_withinRun){

        ## it is the same across experiments. # measurement per feature.
        xtmp <- input[!is.na(input$log2Intensity), ]
        xtmp$eachRun <- paste(xtmp$PSM, xtmp$Run, sep="_")
        count_measure <- xtabs( ~eachRun, xtmp)
        remove_feature_name <- count_measure[count_measure < 10]

        if( length(remove_feature_name) > 0 ){
            xtmp <- xtmp[-which(xtmp$eachRun %in% names(remove_feature_name)), ]
        }
        input <- xtmp[, colnames(xtmp) != "eachRun"]
        message('** Features which has any missing value within a run were removed from that run.')
    }

    ##############################
    ###  7. remove proteins with only one peptide and charge per protein
    ##############################

    if(removeProtein_with1Peptide){

        ######## remove protein which has only one peptide
        tmp <- unique(input[, c("Protein", 'PSM')])
        tmp$Protein <- factor(tmp$Protein)
        count <- xtabs( ~ Protein, data=tmp)
        lengthtotalprotein <- length(count)

        removepro <- names(count[count <= 1])

        if (length(removepro) > 0) {
            input <- input[-which(input$Protein %in% removepro), ]
            message(paste("** ", length(removepro), ' proteins, which have only one feature in a protein, are removed among ', lengthtotalprotein, ' proteins.', sep=""))
        }
    }
    
    ##############################
    ###  8. combine fractions within each mixture
    ##############################
    
    if(Fraction){
      input <- combine.fractions(input)
      message('Fractions belonging to same mixture have been combined.')
    }
    
    return(input)
}


# Remove the peptide ions overlapped among multiple fractions of same biological mixture
# data: PSM level data, which has columns Protein, PSM, Subject, Run, Channel, log2Intensity, BiologicalMixture
combine.fractions <- function(data){
  mixtures <- unique(data$BiologicalMixture)
  data <- as.data.table(data)
  data$Run <- as.character(data$Run)
  all.data <- list()
  for(i in 1: length(mixtures)){
    sub_data <- data[BiologicalMixture == mixtures[i]]
    sub_data <- sub_data[!is.na(log2Intensity)]
    sub_data$fea <- paste(sub_data$PSM, sub_data$Protein, sep="_")
    sub_data$fea <- factor(sub_data$fea)
    
    ## count how many fractions are assigned for each peptide ion
    structure <- aggregate(Run ~., data=unique(sub_data[,.(fea, Run)]), length)
    remove_peptide_ion <- structure[structure$Run > 1, ]
    
    ## remove the peptide ions which are shared by multiple fractions
    if( sum(remove_peptide_ion$Run > 1) != 0 ){
      sub_data <- sub_data[!fea %in% remove_peptide_ion$fea ]
      message('** Peptides, that are shared by more than one fraction of mixture ', mixtures[i],', are removed.')
      
    }
    sub_data_shared_pep_rm <- sub_data[,  fea:= NULL]
    all.data[[i]] <- as.data.table(sub_data_shared_pep_rm)
  }
  data.shared.pep.rm <- rbindlist(all.data)
  data.shared.pep.rm$Run <- data.shared.pep.rm$BiologicalMixture
  data.shared.pep.rm$Subject <- paste(data.shared.pep.rm$Run, data.shared.pep.rm$Channel, sep=".")
  return(data.shared.pep.rm)
}