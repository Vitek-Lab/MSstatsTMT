#' Summarizing PSM level quantification to protein level quantification
#'
#'
#' We assume missing values are censored and then impute the missing values. Protein-level summarization from PSM level quantification should be performed before testing differentially abundant proteins.
#' After all, normalization between MS runs using normalization channels will be implemented.
#'
#' @export
#' @importFrom utils read.table sessionInfo write.table
#' @param data Name of the output of PDtoMSstatsTMTFormat function or PSM-level quantified data from other tools. It should have columns ProteinName, PeptideSequence, Charge, PSM, Mixture, TechRepMixture, Run, Channel, Condition, BioReplicate, Intensity
#' @param method Four different summarization methods to protein-level can be performed : "msstats"(default), "MedianPolish", "Median", "LogSum".
#' @param normalization Normalization between MS runs. TRUE(default) needs at least one normalization channel in each MS run, annotated by 'Norm' in Condtion column. It will be performed after protein-level summarization. FALSE will not perform normalization step. If data only has one run, then normalization=FALSE.
#' @param MBimpute only for method="msstats". TRUE (default) imputes missing values by Accelated failure model. FALSE uses minimum value to impute the missing value for each PSM.
#' @param maxQuantileforCensored We assume missing values are censored. maxQuantileforCensored is Maximum quantile for deciding censored missing value, for instance, 0.999. Default is Null.
#' @param remove_norm_channel TRUE(default) removes 'Norm' channels from protein level data.
#' @param remove_empty_channel TRUE(default) removes 'Empty' channels from protein level data.
#' @return data.frame with protein-level summarization for each run and channel
#' @examples
#' data(input.pd)
#'
#' quant.pd.msstats <- proteinSummarization(input.pd,
#'                                          method="msstats",
#'                                          normalization=TRUE)
#' head(quant.pd.msstats)

proteinSummarization <- function(data,
                                 method = 'msstats',
                                 normalization = TRUE,
                                 MBimpute = TRUE,
                                 maxQuantileforCensored = NULL,
                                 remove_norm_channel = TRUE,
                                 remove_empty_channel = TRUE){

    ## save process output in each step
    allfiles <- list.files()

    num <- 0
    filenaming <- "msstatstmt"
    finalfile <- "msstatstmt.log"

    while (is.element(finalfile,allfiles)) {
        num <- num + 1
        finalfile <- paste(paste(filenaming, num, sep ="-"), ".log", sep = "")
    }

    session <- sessionInfo()
    sink("sessionInfo.txt")
    print(session)
    sink()

    processout <- as.matrix(read.table("sessionInfo.txt", header = TRUE, sep = "\t"))
    write.table(processout, file = finalfile, row.names = FALSE)

    processout <- rbind(processout,
                        as.matrix(c(" ", " ", "MSstatsTMT - proteinSummarization function", " "),
                                  ncol = 1))


    ## check input
    required.info <- c("ProteinName", "PeptideSequence", "Charge", "PSM",
                       "Mixture", "TechRepMixture", "Run",
                       "Channel", "Condition", "BioReplicate", "Intensity")

    if (!all(required.info %in% colnames(data))) {

        missedAnnotation <- which(!(required.info %in% colnames(data)))
        stop(paste("Please check the required input. ** columns :",
                   required.info[missedAnnotation], ", are missed.", collapse = ", "))

    }

    ## check the option for method
    method.list <- c("LogSum", "Median", "MedianPolish", "msstats")

    if (sum(method == method.list) != 1) {
        stop(" 'method' must be one of the following : 'LogSum', 'Median',
             'MedianPolish', 'msstats' default is 'msstats'. ")
    }

    ## report which options are used.
    processout <- rbind(processout,
                        c(paste("Method for protein summarization :", method)))
    processout <- rbind(processout,
                        c(paste("Normalization between MS runs :", normalization)))

    processout <- rbind(processout, c(paste("Remove 'Norm' channels before inference:", remove_norm_channel)))
    processout <- rbind(processout, c(paste("Remove 'Empty' channels before inference:", remove_empty_channel)))
    
    write.table(processout, file = finalfile, row.names = FALSE)

    norm.protein.data <- .protein.summarization.function(data,
                                                         method,
                                                         normalization,
                                                         MBimpute,
                                                         maxQuantileforCensored)
    
    ## remove 'Empty' column : It should not used for further analysis
    if (remove_empty_channel & is.element('Empty', unique(norm.protein.data$Condition))) {
      norm.protein.data <- norm.protein.data[norm.protein.data$Condition != "Empty",]
      norm.protein.data$Condition <- factor(norm.protein.data$Condition)
    }
    
    ## remove 'Norm' column : It should not used for further analysis
    if (remove_norm_channel & is.element('Norm', unique(norm.protein.data$Condition))) {
      norm.protein.data <- norm.protein.data[norm.protein.data$Condition != "Norm",]
      norm.protein.data$Condition <- factor(norm.protein.data$Condition)
    }
    
    return(norm.protein.data)
}
