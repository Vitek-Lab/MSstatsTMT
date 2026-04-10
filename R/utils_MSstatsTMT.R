#' MSstatsTMT: A package for protein significance analysis in shotgun mass spectrometry-based proteomic experiments with tandem mass tag (TMT) labeling
#'
#' A set of tools for detecting differentially abundant peptides and proteins in shotgun mass spectrometry-based proteomic experiments with tandem mass tag (TMT) labeling.
#'
#' @section functions :
#' \itemize{
#'   \item \code{\link{PDtoMSstatsTMTFormat}} : generates MSstatsTMT required input format for Proteome discoverer output.
#'   \item \code{\link{MaxQtoMSstatsTMTFormat}} : generates MSstatsTMT required input format for MaxQuant output.
#'   \item \code{\link{SpectroMinetoMSstatsTMTFormat}} : generates MSstatsTMT required input format for SpectroMine output.
#'   \item \code{\link{OpenMStoMSstatsTMTFormat}} : generates MSstatsTMT required input format for OpenMS output.
#'   \item \code{\link{proteinSummarization}} : summarizes PSM level quantification to protein level quantification.
#'   \item \code{\link{dataProcessPlotsTMT}} : visualizes for explanatory data analysis.
#'   \item \code{\link{groupComparisonTMT}} : tests for significant changes in protein abundance across conditions.
#' }
#'
#' @name MSstatsTMT
#' @keywords internal
"_PACKAGE"

#' Example of output from PDtoMSstatsTMTFormat function
#'
#' It is made from \code{\link{raw.pd}} and \code{\link{annotation.pd}},
#' which is the output of PDtoMSstatsTMTFormat function.
#' It should include the required columns as below.
#'
#' \itemize{
#'   \item ProteinName : Protein ID
#'   \item PeptideSequence : peptide sequence
#'   \item Charge : peptide charge
#'   \item PSM : peptide ion and spectra match
#'   \item Channel : Labeling information (126, ... 131)
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item BioReplicate : Unique ID for biological subject.
#'   \item Run : MS run ID
#'   \item Mixture : Unique ID for TMT mixture.
#'   \item TechRepMixture : Unique ID for technical replicate of one TMT mixture.
#'   \item Intensity: Protein Abundance
#' }
#'
#' @format A data frame with 20110 rows and 11 variables.
#' @examples
#' head(input.pd)
#'
"input.pd"

#' Example of output from proteinSummarizaiton function
#'
#' It is made from \code{\link{input.pd}}.
#' It is the output of proteinSummarization function.
#' It is a list that consists of two data.frames with 
#' feature-level (FeatureLevelData) and protein-level data (ProteinLevelData).
#' ProteinLevelData should include the required columns as below.
#'
#' \itemize{
#'   \item Run : MS run ID
#'   \item Protein : Protein ID
#'   \item Abundance: Protein-level summarized abundance
#'   \item Channel : Labeling information (126, ... 131)
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item BioReplicate : Unique ID for biological subject.
#'   \item TechRepMixture : Unique ID for technical replicate of one TMT mixture.
#'   \item Mixture : Unique ID for TMT mixture.
#' }
#'
#' @format A data frame with 100 rows and 8 variables.
#' @examples
#' head(quant.pd.msstats$ProteinLevelData)
#'
#' @keywords internal
#' 
"quant.pd.msstats"


#' Example of output from groupComparisonTMT function
#'
#' It is the output of groupComparisonTMT function,
#' which is made from \code{\link{quant.pd.msstats}}.
#' It is a list that consists of the following elements: 
#' (1) ComparisonResult: statistical testing results; 
#' (2) FittedModel: the fitted linear models
#' ComparisonResult should include the columns as below.
#'
#' \itemize{
#'   \item Protein : Protein ID
#'   \item Label: Label of the pairwise comparision or contrast
#'   \item log2FC: Log2 fold change
#'   \item SE: Standard error of the comparsion of contrast results
#'   \item DF: Degree of freedom
#'   \item pvalue: Value of p statistic of the test
#'   \item adj.pvalue: adjusted p value
#'   \item issue: used for indicating the reason why a comparison is not testable. NA means the comparison is testable. 
#'   'oneConditionMissing' means the protein has no measurements in one conndition of the comparison.
#'   Furtherone, when 'issue = oneConditionMissing', 'log2FC = Inf' means the negative condition 
#'   (with coefficient -1 in the Label column)  is missing and 'log2FC = -Inf' means 
#'   the positive condition (with coefficient 1 in the Label column)  is missing.
#'   completeMissing' means the protein has no measurements in all the connditions of the comparison.
#'   unfittableModel' means there is no enough measurements to fit the linear model. 
#'   In other words, each condition has only one measurement.
#' }
#'
#' @format A data frame with 60 rows and 7 variables.
#' @examples
#' head(test.pairwise$ComparisonResult)
#'
#' @keywords internal
#' 
"test.pairwise"