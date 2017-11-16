#' Group comparison on protein-level data
#'
#' @param data data: protein level data, which has columns Protein, Group, Subject, Run, Channel, IonIntensity
#' @param model Possible options: "proposed", "limma", "t"
#' @return The sum of \code{x} and \code{y}.//TODO
#' @examples TODO
#' @export
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
# data: protein level data, which has columns Protein, Group, Subject, Run, Channel, IonIntensity
# adj.method: adjusted method for multiple comparison

# Limma inference model
# data: protein level data matrix, whose columns are subjects and rows are proteins.
# label: vector with group information, whose columns are subjects and rows are proteins.
# adj.method: adjusted method for multiple comparison

# t test
# data: protein level data, which has columns Protein, Group, Subject, Run, Channel, IonIntensity
# adj.method: adjusted method for multiple comparison

