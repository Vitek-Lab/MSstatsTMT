#' Summarize PSM level data to protein level
#'
#' @param data PSM level data, which has columns Protein, PSM, Subject, Run, Channel, log2Intensity
#' @param method summarization methods. Possible options: "LogSum", "Median", "Biweight", "MedianPolish", "Huber"
#' @return Protein Abundance
#' @examples MedianPolish.abun <- protein.summarization(data, annotation,  "MedianPolish")
#' @export

protein.summarization <- function(data, method){
    #check input

    #################
    ## MC, 20171128 : start
    ## can you make sure to check all together and report the list of missing information?
    ## With current version, it seems that it will stop with one by one.

   # required.annotation <- c('Run', 'BioReplicate', 'Condition')

   # if ( !all(required.annotation %in% colnames(annotation)) ) {

   #     missedAnnotation <- which(!(required.annotation %in% colnames(annotation)))

   #     stop(paste("**", paste(required.annotation[missedAnnotation], collapse = ", "), "is not provided in Annotation. Please check the annotation."))
   # }
    ## MC, 20171128 : start
    #################


    #data
    if(is.null(data$Run)){
        stop("Please make sure the data has a colume called 'Run'!")
    }
    if(is.null(data$Channel)){
        stop("Please make sure input data has a colume called 'Channel'!")
    }
    if(is.null(data$Protein)){
        stop("Please make sure input data has a colume called 'Protein'!")
    }
    if(is.null(data$PSM)){
        stop("Please make sure input data has a colume called 'PSM'!")
    }
    if(is.null(data$log2Intensity)){
        stop("Please make sure input data has a colume called 'log2Intensity'!")
    }
    if(!all.equal(length(data$Run),length(data$Channel),length(data$Protein),length(data$PSM),length(data$log2Intensity),length(data$Subject))){
        stop("Please make sure all columes have same length")
    }
    if(!is.numeric(data$log2Intensity)){
        stop("Please make sure 'log2Intensity' is numeric!")
    }
    #annotation
    if(is.null(annotation$Run)){
        stop("Please make sure input annotation data has a colume called 'Run'!")
    }
    if(is.null(annotation$Channel)){
        stop("Please make sure input annotation data has a colume called 'Channel'!")
    }
    if(is.null(annotation$Group)){
        stop("Please make sure input annotation data has a colume called 'Group'!")
    }
    if(!all.equal(length(annotation$Run),length(annotation$Channel),length(annotation$Group))){
        stop("Please make sure all columes have same length")
    }
    #method
    method.list<-c("LogSum", "Median", "Biweight", "MedianPolish", "Huber")
    if(sum(method==method.list)!=1){
        stop(" 'Method' must be one of the following, 'LogSum', 'Median', 'Biweight', 'MedianPolish', 'Huber' default is 'LogSum' ")
    }
    return(protein.summarization.function(data, method))
}

