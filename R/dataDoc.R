
#' Example of output from Proteome Discoverer for TMT10 experiments.
#'
#' Example of Proteome discover PSM sheet. It is the input for PDtoMSstatsTMTFormat function, with annotation file.
#' It includes peak intensities for 10 proteins among 15 MS runs with TMT10.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Master.Protein.Accessions
#'   \item Protein.Accessions
#'   \item Annotated.Sequence
#'   \item Charge
#'   \item Ions.Score
#'   \item Spectrum.File
#'   \item Quan.Info
#'   \item Channels : X126, ... X131
#' }
#'
#' @format A data frame with 1928 rows and 51 variables.
#' @examples
#' head(raw.input)
#'
"raw.input"

#' Example of annotation file for raw.input
#'
#' Annotation of example data, raw.input, in this package. It should be prepared by users.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Run : MS run ID
#'   \item Channel : Labeling ID
#'   \item Group : Condition (ex. Healthy, Cancer, Time0)
#'   \item BiologicalMixture : Mixture ID
#'   \item Subject : Unique ID for biological subject
#' }
#'
#' @format A data frame with 30 rows and 5 variables.
#' @examples
#' head(annotation)
#'
"annotation"


#' Example of output from PDtoMSstatsTMTFormat function
#'
#' It is formated from raw.input and annotation.
#' It is the output of PDtoMSstatsTMTFormat function and the input for protein.summarization function.
#' It should includes the required columns as below.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Protein : Protein ID
#'   \item PSM : combination of Peptide sequence and charge
#'   \item log2Intensity : log 2 transformed intensity
#'   \item Channel : Labeling ID
#'   \item Subject : Unique ID for biological subject.
#'   \item Run : MS run ID
#'   \item BiologicalMixture : Mixture ID
#'   \item Group : Condition (ex. Healthy, Cancer, Time0)
#' }
#'
#' @format A data frame with 14740 rows and 8 variables.
#' @examples
#' head(required.input)
#'
"required.input"
