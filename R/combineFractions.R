# Combine multiple fractions of same biological mixture. Remove the overlapped petide ions among multiple fractions.
# data: PSM level data, which has columns Protein, PSM, Subject, Run, Channel, log2Intensity, BiologicalMixture
combine.fractions <- function(data){
    mixtures <- unique(data$BiologicalMixture)
    data <- as.data.table(data)
    data$Run <- as.character(data$Run)
    all.data <- list()

    # Combine fractions for each biological mixture
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

    #The run becomes the biological mixture
    data.shared.pep.rm$Run <- data.shared.pep.rm$BiologicalMixture
    data.shared.pep.rm$Subject <- paste(data.shared.pep.rm$Run, data.shared.pep.rm$Channel, sep=".")
    return(data.shared.pep.rm)
}
