#' Finding differentially abundant proteins across conditions in TMT experiment
#'
#' Tests for significant changes in protein abundance across conditions based on a family of linear mixed-effects models in TMT experiment.
#' Experimental design of case-control study (patients are not repeatedly measured) or time course study (patients are repeatedly measured) is automatically determined based on proper statistical model.
#'
#' @export
#' @param data Name of the output of protein.summarization function. It should have columns named Protein, Mixture, Run, Channel, Condition, BioReplicate, Abundance.
#' @param contrast.matrix Comparison between conditions of interests. 1) default is 'pairwise', which compare all possible pairs between two conditions. 2) Otherwise, users can specify the comparisons of interest. Based on the levels of conditions, specify 1 or -1 to the conditions of interests and 0 otherwise. The levels of conditions are sorted alphabetically.
#' @param remove_norm_channel TRUE(default) removes 'Norm' channels for inference step.
#' @param model Three different statistical approached can be performed : "proposed", "limma", "ttest". "proposed" is the default.
#' @param moderated Only for model = 'proposed'. If moderated = TRUE, then moderated t statistic will be calculated; otherwise, ordinary t statistic will be used.
#' @param adj.method adjusted method for multiple comparison. "BH" is default.
#' @return data.frame with result of inference
#' @examples
#' quant.byprotein <- protein.summarization(required.input,
#'                                          method="msstats",
#'                                          normalization=TRUE,)
#'
#' test.byproposed <- groupComparison.TMT(quant.byprotein)
#'
#' # Only compare condition 0.125 and 1
#' levels(quant.byprotein$Condition)
#' # 'Norm' should be not considered in the contrast
#' comparison<-matrix(c(-1,0,0,1),nrow=1)
#' # Set the names of each row
#' row.names(comparison)<-"1-0.125"
#' # Set the column names
#' colnames(comparison)<- c("0.125", "0.5", "0.667", "1")
#' test.contrast <- groupComparison.TMT(data = quant.byprotein, contrast.matrix = comparison)

groupComparison.TMT <- function(data,
                                contrast.matrix = 'pairwise',
                                remove_norm_channel = TRUE,
                                model = 'proposed',
                                moderated = TRUE,
                                adj.method = "BH"){

    ## save process output in each step
    allfiles <- list.files()
    filenaming <- "msstatstmt"

    if (length(grep(filenaming, allfiles)) == 0) {

        finalfile <- "msstatstmt.log"
        processout <- NULL

    } else {

        num <- 0
        finalfile <- "msstatstmt.log"

        while (is.element(finalfile, allfiles)) {
            num <- num + 1
            lastfilename <- finalfile ## in order to rea
            finalfile <- paste(paste(filenaming, num, sep="-"), ".log", sep="")
        }

        finalfile <- lastfilename
        processout <- as.matrix(read.table(finalfile, header=TRUE, sep="\t"))
    }

    processout <- rbind(processout, as.matrix(c(" ", " ", "MSstatsTMT - groupComparison.TMT function"," "), ncol=1))


    ## check input data
    required.info <- c('Protein', 'BioReplicate', 'Abundance', 'Run', 'Channel', 'Condition', 'Mixture')

    if (!all(required.info %in% colnames(data))) {

        missedAnnotation <- which(!(required.info %in% colnames(data)))
        stop(paste("Please check the required input. ** columns :",
                   required.info[missedAnnotation],
                   ", are missed.", collapse = ", "))

    }

    ## change some column names as used in group comparison function
    colnames(data)[colnames(data) == 'BioReplicate'] <- 'Subject'
    colnames(data)[colnames(data) == 'Condition'] <- 'Group'

    ## check the option for model
    model.list <- c("ttest", "limma", "proposed")

    if (sum(model == model.list) != 1) {
        stop(" 'model' must be one of the following : 'proposed', 'limma', 'ttest'. Default is 'proposed'. ")
    }

    ## report which options are used.
    processout <- rbind(processout, c(paste("Remove 'Norm' channels before inference:", remove_norm_channel)))
    processout <- rbind(processout, c(paste("Model for inference :", model)))
    processout <- rbind(processout, c(paste("Moderated t-stat :", moderated)))

    write.table(processout, file=finalfile, row.names=FALSE)

    ## remove 'Norm' column : It should not used for inference
    if (remove_norm_channel & is.element('Norm', unique(data$Group))) {
        data <- data[data$Group != "Norm",]
        data$Group <- factor(data$Group)
    }

    ## Inference
    if (model == "proposed") {
        result <- proposed.model(data, moderated, contrast.matrix, adj.method)
    } else if (model == "ttest") {
        #if( is.matrix(contrast.matrix) ){ ## maybe better way later
        #message("** For t-test, all pairwise comparisons will be reported.")
        #}
        result <- protein.ttest(data, contrast.matrix, adj.method)
    } else if (model == "limma") {
        result <- ebayes.limma(data, contrast.matrix, adj.method)
    }

    ### check column name in order to use groupComparisonPlot from MSstats
    colnames(result)[colnames(result) == 'Comparison'] <- 'Label'
    colnames(result)[colnames(result) == 'adjusted.pvalue'] <- 'adj.pvalue'

    return(result)
}

# Proposed inference model
# data: protein level data, which has columns Protein, Group, Subject, Run, Channel, log2Intensity
# adj.method: adjusted method for multiple comparison


