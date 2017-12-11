#' @import data.table
# Do normalization after protein summarization. Normalization are based on the channels which have group 'Norm'.
# data: protein level data, which has columns Protein, Group, Subject, Run, Channel, Abundance
protein.normalization <-function(data) {
    data$Protein <- as.character(data$Protein) # make sure protein names are character
    proteins <- unique(data$Protein) # proteins
    data <- as.data.table(data) # make suree the input data is with data table format
    norm.data <- list()
    # do inference for each protein individually
    for(i in 1:length(proteins)) {
        message("Protein: ", i)
        sub_data <- data[Protein == proteins[i]] # data for protein i
        norm.channel <- sub_data[Group == "Norm"]
        if(length(channel) > 1){
            norm.channel = norm.channel[, .(Abundance = mean(Abundance, na.rm = T)), by = .(Protein, Run)]
        }
        norm.channel$diff <- median(norm.channel$Abundance, na.rm = T) - norm.channel$Abundance
        setkey(sub_data, Run)
        setkey(norm.channel, Run)
        norm.sub_data <- merge(sub_data, norm.channel[,.(Run, diff)], all.x=TRUE)
        norm.sub_data$Abundance <- norm.sub_data$Abundance + norm.sub_data$diff
        norm.sub_data[,diff:=NULL]
        norm.data[[proteins[i]]] <- norm.sub_data
    }
    norm.data <- rbindlist(norm.data)
    return(norm.data)
}


#' @import data.table
# Do normalization before protein summarization. Normalization are based on the channels which have group 'Norm'.
# data: PSM level data, which has columns Protein, PSM, Subject, Run, Channel, log2Intensity, BiologicalMixture
peptide.normalization <-function(data) {
    data$PSM <- as.character(data$PSM) # make sure protein names are character
    PSMs <- unique(data$PSM) # proteins
    data <- as.data.table(data) # make sure the input data is with data table format
    norm.data <- list()

    norm.data <- data
    # do inference for each protein individually
    for(i in 1:length(PSMs)) {
        if(i%%100==0){
            message("PSM: ", i)
        }
        sub_data <- u[PSM == PSMs[i]] # data for protein i
        norm.channel <- sub_data[Group == "Norm"]
        if(length(channel) > 1){
            norm.channel = norm.channel[, .(Abundance = mean(Abundance, na.rm = T)), by = .(Protein, Run)]
        }
        norm.channel$diff <- median(norm.channel$log2Intensity, na.rm = T) - norm.channel$log2Intensity
        setkey(sub_data, Run)
        setkey(norm.channel, Run)
        norm.sub_data <- merge(sub_data, norm.channel[,.(Run, diff)], all.x=TRUE)
        norm.sub_data$log2Intensity <- norm.sub_data$log2Intensity + norm.sub_data$diff
        norm.sub_data[,diff:=NULL]
        norm.data[[PSMs[i]]] <- norm.sub_data
    }
    norm.data <- rbindlist(norm.data)
    return(norm.data)
}
