#' Group comparison on protein-level data
#'
#' @param data data: protein level data, which has columns Protein, Group, Subject, Run, Channel, log2Intensity
#' @param model Possible options: "proposed", "limma", "t"
#' @return Result of groupComparison
#' @export
#' @examle
#' quant.byprotein <- protein.summarization(required.input, "MedianPolish")
#' test.byproposed <- groupComparison.TMT(quant.byprotein, model = "proposed")

groupComparison.TMT <- function(data,
                                model = "proposed"){

    ### check the option for model
    model.list <- c("ttest", "limma", "proposed")

    if( sum(model == model.list) != 1 ){
        stop(" 'model' must be one of the following : 'proposed', 'limma', 'ttest'. Default is 'proposed'. ")
    }

    if(model == "proposed"){
        result <- MSstatsTMT::proposed.model(data)
    } else if(model == "ttest"){
        result <- MSstatsTMT::protein.ttest(data)
    } else if(model == "limma"){
        result <- MSstatsTMT::ebayes.limma(data)
    }

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

