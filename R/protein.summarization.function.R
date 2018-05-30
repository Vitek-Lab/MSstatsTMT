#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom stats medpolish
#' @importFrom matrixStats colMedians
#' @importFrom MASS huber
#' @importFrom affy generateExprVal.method.mas
#' @importFrom data.table :=
protein.summarization.function <- function(data, method, normalization){

    data <- as.data.table(data)
    ## make sure the protein ID is character
    data$ProteinName <- as.character(data$ProteinName)
    ## make new column: combination of run and channel
    data$runchannel <- paste(data$Run, data$Channel, sep = '_')

    ## makae intensity transformed
    data$log2Intensity <- log2(data$Intensity)

    ## Number of negative values : if intensity is less than 1, replace with zero
    ## then we don't need to worry about -Inf = log2(0)
    if (nrow(data[!is.na(data$Intensity) & data$Intensity < 1]) > 0){
        data[!is.na(data$Intensity) & data$Intensity < 1, 'log2Intensity'] <- 0
        message('** Negative log2 intensities were replaced with zero.')
    }

    ## Record the group information
    annotation <- unique(data[ , c('Run', 'Channel', 'BioReplicate', 'Condition', 'Mixture', 'runchannel')])
    data <- data[, c('ProteinName', 'PSM', 'log2Intensity', 'Run', 'Channel', 'BioReplicate', 'runchannel')]

    proteins <- unique(data$ProteinName)
    num.protein <- length(proteins)
    runs <- unique(data$Run)
    runchannel.id <- unique(data$runchannel)

    ## Store the estimated protein abundance
    protein.abundance <- matrix(rep(NA, length(runchannel.id)*length(proteins)), nrow = length(proteins))
    colnames(protein.abundance) <- runchannel.id

    ## For each protein and each run, do the summarization individually
    for (i in 1:length(proteins)) {
        message(paste("Summarizing for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))

        for (j in 1:length(runs)) {
            sub_data <- data %>% filter(ProteinName == proteins[i] & Run == runs[j])

            if (nrow(sub_data) != 0) {
                nfea <- length(unique(sub_data$PSM))
                ## Change the long format to wide format
                sub_data_wide <- sub_data %>%
                    dplyr::select(log2Intensity, PSM, runchannel) %>%
                    tidyr::spread(runchannel, log2Intensity)
                rownames(sub_data_wide) <- sub_data_wide[,1]
                sub_data_wide <- sub_data_wide[,-1]

                ## Number of negative values
                #index <- which(apply(sub_data_wide, 1, function(col) any(col < 0)))
                #if (length(index) != 0) {
                #    sub_data_wide[!is.na(sub_data_wide) & sub_data_wide < 0 ] <- 0
                #}

                if (nrow(sub_data_wide) != 0) {
                    if (nrow(sub_data_wide) == 1) { # Only one PSM for the protein
                        protein.abundance[i, colnames(sub_data_wide)] <- as.matrix(sub_data_wide)
                    } else {
                        if (method == "LogSum") {
                            # log2 (sum of original intensity)
                            protein.abundance[i, colnames(sub_data_wide)] <- log2(colSums(2^sub_data_wide, na.rm = TRUE))
                        }

                        if (method == "Median") {
                            #Median
                            protein.abundance[i, colnames(sub_data_wide)] <- colMedians(as.matrix(sub_data_wide, na.rm = TRUE))
                        }
                        if (method == "Biweight") {
                            #Biweight
                            protein.abundance[i, colnames(sub_data_wide)] <- log2(generateExprVal.method.mas(as.matrix(2^sub_data_wide))$exprs)
                        }
                        if (method == "MedianPolish") {
                            #median polish
                            meddata  <-  stats::medpolish(as.matrix(sub_data_wide), na.rm=TRUE, trace.iter = FALSE)
                            tmpresult <- meddata$overall + meddata$col
                            protein.abundance[i, colnames(sub_data_wide)] <- tmpresult[colnames(sub_data_wide)]
                        }
                        if (method == "Huber") {
                            #Huber
                            protein.abundance[i, colnames(sub_data_wide)] <- unlist(apply(as.matrix(sub_data_wide), 2,
                                                                                          function(x) huber(x, k = 1.345)$mu))
                        }
                    }
                }
            }
        } ## end for loop by run
    } ## end for loop by protein

    rownames(protein.abundance) <- proteins

    ## Get the group information for each subject
    ## Make the data long format and add the group information to protein level data frame
    res <- as.data.frame(protein.abundance)
    res$Protein <- rownames(res)
    res <- res %>% gather(runchannel, Abundance, -Protein) # Change to long format
    res$runchannel <- as.character(res$runchannel)
    annotation$runchannel <- as.character(annotation$runchannel)
    res <- left_join(res, annotation, by='runchannel')

    ## remove runchannel column
    res <- res[, -which(colnames(res) %in% 'runchannel')]

    if (normalization & length(runs) > 1) { # Do normalization based on group 'Norm'
      res <- protein.normalization(res)
    }

    return(res)
}



###########################################
## function for normalization between MS runs
## Do normalization after protein summarization.
## Normalization are based on the channels which have group 'Norm'.
## data: protein level data, which has columns Protein, Group, Subject, Run, Channel, Abundance

protein.normalization <- function(data) {

    ## check whethere there are 'Norm' info or not.
    group.info <- unique(data$Condition)

    ## if there is 'Norm' available in Condition column,
    if (is.element('Norm', group.info)) {

        data$Protein <- as.character(data$Protein) # make sure protein names are character
        proteins <- unique(data$Protein) # proteins
        num.protein <- length(proteins)
        data <- as.data.table(data) # make suree the input data is with data table format
        norm.data <- list()

        # do inference for each protein individually
        for (i in 1:length(proteins)) {

            message(paste("Normalization between MS runs for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))

            sub_data <- data[Protein == proteins[i]] # data for protein i
            sub_data <- na.omit(sub_data)
            norm.channel <- sub_data[Condition == "Norm"]
            norm.channel <- norm.channel[, .(Abundance = mean(Abundance, na.rm = TRUE)), by = .(Protein, Run)]
            norm.channel$diff <- median(norm.channel$Abundance, na.rm = TRUE) - norm.channel$Abundance
            setkey(sub_data, Run)
            setkey(norm.channel, Run)
            norm.sub_data <- merge(sub_data, norm.channel[, .(Run, diff)], all.x=TRUE)
            norm.sub_data$Abundance <- norm.sub_data$Abundance + norm.sub_data$diff
            norm.sub_data[, diff:=NULL]
            norm.data[[proteins[i]]] <- norm.sub_data
        }
        norm.data <- rbindlist(norm.data)

    } else {
        message("** 'Norm' information in Condition is required for normalization.
                Please check it. At this moment, normalization is not performed.")
        norm.data <- data
    }

    return(norm.data)
}
