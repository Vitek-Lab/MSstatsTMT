#' Summarizing PSM level quantification to protein level quantification
#'
#' Protein-level summarization from PSM level quantification should be performed before testing differentially abundant proteins.
#' Then, normalization between MS runs using normalization channels will be implemented.
#'
#' @export
#' @param data Name of the output of PDtoMSstatsTMTFormat function or PSM-level quantified data from other tools. It should have columns named Protein, PSM, BiolobicalMixture, Run, Channel, Group, Subject, log2Intensity
#' @param method Five different summarization methods to protein-level can be performed : "MedianPolish"(default), "Huber", "LogSum", "Median", "Biweight".
#' @param normalization Normalization between MS runs. TRUE(default) needs at least normalization channel in each MS run, annotated by 'Norm' in Group column. It will be performed after protein-level summarization. FALSE will not perform normalization step.
#' @examples
#' head(required.input)
#' str(required.input)
#' quant.byprotein <- protein.summarization(required.input, method="MedianPolish", normalization=TRUE)
#' head(quant.byprotein)

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

