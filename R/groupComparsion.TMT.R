#' Group comparison on protein-level data
#'
#' @param data data: protein level data, which has columns Protein, Group, Subject, Run, Channel, log2Intensity
#' @param model Possible options: "proposed", "limma", "t"
#' @return Result of groupComparison
#' @export
#' @examle
#' MedianPolish.abun <- protein.summarization(input.data, "MedianPolish")
#' proposed.res <- groupComparison.TMT(MedianPolish.abun,model = "proposed")

groupComparison.TMT<-function(data,model = "proposed"){

  if(model == "proposed"){
    result<-MSstatsTMT::proposed.model(data)
  } else if(model == "t"){
    result<-MSstatsTMT::protein.ttest(data)
  } else if(model == "limma"){
    result<-MSstatsTMT::ebayes.limma(data)
  } else {
      stop("Please enter a model name, options are 'proposed','t' and 'limma' (default is 'proposed')")
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

