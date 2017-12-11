#' Summarize PSM level data to protein level
#'
#' @param data PSM level data, which has columns Protein, PSM, Subject, Run, Channel, log2Intensity
#' @param method summarization methods. Possible options: "LogSum", "Median", "Biweight", "MedianPolish", "Huber"
#' @param normalization If TRUE, do normalization after protein summarization. Normalization are based on the channels which have group 'Norm'.
#' @return Protein Abundance
#' @export
#' @examples
#' head(required.input)
#' str(required.input)
#' MedianPolish.abun <- protein.summarization(required.input, "MedianPolish")
#' head(MedianPolish.abun)

protein.summarization <- function(data,
                                  method = 'MedianPolish',
                                  normalization = TRUE){
    ### check input
    required.info <- c('Protein', 'PSM', 'Channel', 'Subject', 'Run', 'BiologicalMixture', 'Group', 'log2Intensity')

    if ( !all(required.info %in% colnames(data)) ) {

        missedAnnotation <- which(!(required.info %in% colnames(data)))
        stop(paste("Please check the required input. ** columns :", required.info[missedAnnotation], ", are missed.", collapse = ", "))

    }

    ### check the option for method
    method.list <- c("LogSum", "Median", "Biweight", "MedianPolish", "Huber")

    if(sum(method==method.list)!=1){
        stop(" 'method' must be one of the following : 'LogSum', 'Median', 'Biweight', 'MedianPolish', 'Huber' default is 'MedianPolish'. ")
    }

    return(protein.summarization.function(data, method, normalization))
}

