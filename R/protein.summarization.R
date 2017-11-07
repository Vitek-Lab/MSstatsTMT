#' Summarize PSM level data to protein level
#'
#' @param data PSM level data, which has columns Protein, PSM, Subject, Run, Channel, IonIntensity
#' @param method summarization methods. Possible options: "LogSum", "Median", "Biweight", "MedianPolish", "Huber"
#' @return The sum of \code{x} and \code{y}.//TODO
#' @examples LogSum.abun <- MSstatsTMT::protein.summarization(MSstatsTMT::test.data, MSstatsTMT::annotation.data,  "LogSum")
#' print(head(LogSum.abun))
#' @import tidyr
#' @import data.table
#' @import MASS
#' @import matrixStats
#' @import dplyr
#' @import tidyr
#' @import affy
#' @export

protein.summarization <- function(data, annotation, method){
  #check input
    #data
    if(is.null(data$Run)){
        message("Please make sure the data has a colume called 'Run'!")
    }
    if(is.null(data$Channel)){
        message("Please make sure input data has a colume called 'Channel'!")
    }
    if(is.null(data$Protein)){
        message("Please make sure input data has a colume called 'Protein'!")
    }
    if(is.null(data$PSM)){
        message("Please make sure input data has a colume called 'PSM'!")
    }
    if(is.null(data$IonIntensity)){
        message("Please make sure input data has a colume called 'IonIntensity'!")
    }
    if(is.null(data$Subject)){
        message("Please make sure input data has a colume called 'Subject'!")
    }
    if(!all.equal(length(data$Run),length(data$Channel),length(data$Protein),length(data$PSM),length(data$IonIntensity),length(data$Subject))){
        message("Please make sure all columes have same length")
    }
    if(!is.numeric(test.data$IonIntensity)){
        message("Please make sure 'IonIntensity' is numeric!")
    }
    #annotation
    if(is.null(annotation$Run)){
        message("Please make sure input annotation data has a colume called 'Run'!")
    }
    if(is.null(annotation$Channel)){
        message("Please make sure input annotation data has a colume called 'Channel'!")
    }
    if(is.null(annotation$Group)){
        message("Please make sure input annotation data has a colume called 'Group'!")
    }
    if(!all.equal(length(annotation$Run),length(annotation$Channel),length(annotation$Group))){
        message("Please make sure all columes have same length")
    }
    #method
    method.list<-c("LogSum", "Median", "Biweight", "MedianPolish", "Huber")
    if(sum(method==method.list)!=1){
        message(" 'Method' must be one of the following, 'LogSum', 'Median', 'Biweight', 'MedianPolish', 'Huber' default is 'LogSum' ")
    }
  return(protein.summarization.function(data,annotation,method))
}

