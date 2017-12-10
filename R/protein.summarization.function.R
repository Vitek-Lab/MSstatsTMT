
#' @import data.table
#' @import affy
protein.summarization.function<- function(data, method){
    data <- as.data.table(data)
    # Record the group information
    annotation <- unique(data[ ,.(Run, Channel, Subject, Group, BiologicalMixture)])
    data <- data[,.(Protein, PSM, log2Intensity, Run, Channel, Subject)]

    proteins <- unique(data$Protein)
    subjects <- unique(data$Subject)
    runs <- unique(data$Run)

    # Store the estimated protein abundance
    protein.abundance <- matrix(rep(NA, length(subjects)*length(proteins)), ncol = length(subjects))
    colnames(protein.abundance) <- subjects
    # For each protein and each run, do the summarization individually
    for(i in 1:length(proteins)) {
      message("Protein: ", i)
      for(j in 1:length(runs)){
        sub_data <- data %>% filter(Protein == proteins[i] & Run == runs[j])
        if(nrow(sub_data) != 0){
          nfea <- length(unique(sub_data$PSM))
          # Change the long format to wide format
          sub_data_wide <- sub_data %>% dplyr::select(log2Intensity, PSM, Subject) %>% spread(Subject, log2Intensity)
          rownames(sub_data_wide) <- sub_data_wide[,1]
          sub_data_wide <- sub_data_wide[,-1]
          # Number of negative values
          index <- which(apply(sub_data_wide, 1, function(col) any(col < 0)))
          if(length(index) != 0){
            # MC - 20170808 : replace negative values with zero.
            sub_data_wide[!is.na(sub_data_wide) & sub_data_wide < 0 ] <- 0
            message('* replace negatives with zero')
            # end MC- 20170808
          }

          if(nrow(sub_data_wide) != 0){
            if(nrow(sub_data_wide) == 1){ # Only one PSM for the protein
              protein.abundance[i, colnames(sub_data_wide)] <- as.matrix(sub_data_wide)
            } else{
              if(method == "LogSum"){
                #Sum
                # MC- 20170808 : change to log2 (sum of intensity)
                protein.abundance[i, colnames(sub_data_wide)] <- log2(colSums(2^sub_data_wide, na.rm = TRUE))
              }
              if(method == "Median"){
                #Median
                protein.abundance[i, colnames(sub_data_wide)] <- colMedians(as.matrix(sub_data_wide, na.rm = TRUE))
              }
              if(method == "Biweight"){
                #Biweight
                protein.abundance[i, colnames(sub_data_wide)] <- log2(generateExprVal.method.mas(as.matrix(2^sub_data_wide))$exprs)
              }
              if(method == "MedianPolish"){
                #median polish
                meddata  <-  medpolish(as.matrix(sub_data_wide), na.rm=TRUE,trace.iter = FALSE)
                tmpresult <- meddata$overall + meddata$col
                protein.abundance[i, colnames(sub_data_wide)] <- tmpresult[colnames(sub_data_wide)]
              }
              if(method == "Huber"){
                #Huber
                protein.abundance[i, colnames(sub_data_wide)] <- unlist(apply(as.matrix(sub_data_wide), 2, function(x) huber(x, k = 1.345)$mu))
              }
            }
          }
        }
      }#end of second for loop


    }#end of first for loop
    rownames(protein.abundance) <- proteins
    # Get the group information for each subject
    # Make the data long format and add the group information to protein level data frame
    res <- as.data.frame(protein.abundance)
    res$Protein <- rownames(res)
    res <- res %>% gather(Subject, Abundance, -Protein) # Change to long format
    res$Subject <- as.character(res$Subject)
    annotation$Subject <- as.character(annotation$Subject)
    res <- left_join(res, annotation)
    return(res)
  }
