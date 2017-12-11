#' Finding differentially abundant proteins across conditions in TMT experiment
#'
#' Tests for significant changes in protein abundance across conditions based on a family of linear mixed-effects models in TMT experiment.
#' Experimental design of case-control study (patients are not repeatedly measured) or time course study (patients are repeatedly measured) is automatically determined based on proper statistical model.
#'
#' @export
#' @param data Name of the output of protein.summarization function. It should have columns named Protein, BiologicalMixture, Run, Channel, Group, Subject, Abundance.
#' @param contrast.matrix Comparison between conditions of interests. 1) default is 'pairwise', which compare all possible pairs between two conditions. 2) Otherwise, users can specify the comparisons of interest. Based on the levels of conditions, specify 1 or -1 to the conditions of interests and 0 otherwise. The levels of conditions are sorted alphabetically.
#' @param remove_norm_channel TRUE(default) removes 'Norm' channels for inference step.
#' @param model Three different statistical approached can be performed : "proposed", "limma", "ttest". "proposed" is the default.
#' @examples
#' quant.byprotein <- protein.summarization(required.input, method = "MedianPolish", normalization=TRUE)
#' test.byproposed <- groupComparison.TMT(quant.byprotein, model = "proposed")


groupComparison.TMT <- function(data,
                                contrast.matrix = 'pairwise',
                                remove_norm_channel = TRUE,
                                model = 'proposed'){

    ## save process output in each step
    allfiles <- list.files()
    filenaming <- "msstatstmt"

    if (length(grep(filenaming, allfiles)) == 0) {

        finalfile <- "msstatstmt.log"
        processout <- NULL

    } else {

        num <- 0
        finalfile <- "msstatstmt.log"

        while(is.element(finalfile, allfiles)) {
            num <- num + 1
            lastfilename <- finalfile ## in order to rea
            finalfile <- paste(paste(filenaming, num, sep="-"), ".log", sep="")
        }

        finalfile <- lastfilename
        processout <- as.matrix(read.table(finalfile, header=TRUE, sep="\t"))
    }

    processout <- rbind(processout, as.matrix(c(" ", " ", "MSstatsTMT - groupComparison.TMT function"," "), ncol=1))


    ### check input data
    required.info <- c('Protein', 'Subject', 'Abundance', 'Run', 'Channel', 'Group', 'BiologicalMixture')

    if ( !all(required.info %in% colnames(data)) ) {

        missedAnnotation <- which(!(required.info %in% colnames(data)))
        stop(paste("Please check the required input. ** columns :", required.info[missedAnnotation], ", are missed.", collapse = ", "))

    }

    ### check the option for model
    model.list <- c("ttest", "limma", "proposed")

    if( sum(model == model.list) != 1 ){
        stop(" 'model' must be one of the following : 'proposed', 'limma', 'ttest'. Default is 'proposed'. ")
    }

    ### report which options are used.
    processout <- rbind(processout, c(paste("Remove 'Norm' channels before inference:", remove_norm_channel)))
    processout <- rbind(processout, c(paste("Model for inference :", model)))

    write.table(processout, file=finalfile, row.names=FALSE)


    ### remove 'Norm' column
    if( remove_norm_channel & is.element('Norm', unique(data$Group)) ){
        data <- data %>% filter(Group != "Norm")
        data$Group <- factor(data$Group)
    }

    ### Inference
    if(model == "proposed"){
        result <- MSstatsTMT::proposed.model(data, contrast.matrix)
    } else if(model == "ttest"){
        if( is.matrix(contrast.matrix) ){ ## maybe better way later
            message("** For t-test, all pairwise comparisons will be reported.")
        }
        result <- MSstatsTMT::protein.ttest(data)
    } else if(model == "limma"){
        result <- MSstatsTMT::ebayes.limma(data, contrast.matrix)
    }

    ### check column name in order to use groupComparisonPlot from MSstats
    colnames(result)[colnames(result) == 'Comparison'] <- 'Label'
    colnames(result)[colnames(result) == 'adjusted.pvalue'] <- 'adj.pvalue'

    return(result)
}

# Proposed inference model
# data: protein level data, which has columns Protein, Group, Subject, Run, Channel, log2Intensity
# adj.method: adjusted method for multiple comparison

# Limma inference model
# data: protein level data matrix, whose columns are subjects and rows are proteins.
# label: vector with group information, whose columns are subjects and rows are proteins.
# adj.method: adjusted method for multiple comparison

# t test
# data: protein level data, which has columns Protein, Group, Subject, Run, Channel, log2Intensity
# adj.method: adjusted method for multiple comparison

