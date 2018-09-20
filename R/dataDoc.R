
#' Example of output from Proteome Discoverer 2.2 for TMT10 experiments.
#'
#' Example of Proteome discover PSM sheet.
#' It is the input for PDtoMSstatsTMTFormat function, with annotation file.
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
#'   \item Channels : 126, ..., 131
#' }
#'
#' @format A data frame with 2858 rows and 50 variables.
#' @examples
#' head(raw.pd)
#'
"raw.pd"

#' Example of annotation file for raw.pd
#'
#' Annotation of example data, raw.pd, in this package.
#' It should be prepared by users.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Run : MS run ID. It should be the same as Spectrum.File info
#'   in raw.pd.
#'   \item Channel : Labeling information (126, ... 131). It should be
#'   consistent with the channel columns in raw.pd.
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item Mixture : TMT Mixture ID. It is used to indicate the technical
#'   replicates or fractions. For example, if 'Run' 1, 2, 3 are the fractions
#'   of same TMT mixture, then they should have same 'Mixture'.
#'   \item BioReplicate : Unique ID for biological subject
#' }
#'
#' @format A data frame with 150 rows and 5 variables.
#' @examples
#' head(annotation.pd)
#'
"annotation.pd"

#' Example of output from MaxQuant for TMT10 experiments.
#'
#' Example of MaxQuant peptide evidence file.
#' It is the input for MaxQtoMSstatsTMTFormat function, with proteinGroups.txt
#'  and annotation file.
#' It includes peak intensities for 10 proteins among 15 MS runs with TMT10.
#' The important variables are as follows:
#'
#' \itemize{
#'   \item Proteins
#'   \item Protein.group.IDs
#'   \item Modified.sequence
#'   \item Charge
#'   \item Raw.file
#'   \item Score
#'   \item Potential.contaminant
#'   \item Reverse
#'   \item Channels : Reporter.intensity.corrected.0, ...,
#'   Reporter.intensity.corrected.9
#' }
#'
#' @format A data frame with 1075 rows and 105 variables.
#' @examples
#' head(evidence)
#'
"evidence"

#' Example of proteinGroups file from MaxQuant for TMT10 experiments.
#'
#' Example of MaxQuant identified protein group information file.
#' It is the input for MaxQtoMSstatsTMTFormat function, with evidence.txt
#'  and annotation file.
#' It includes identified protein groups for 10 proteins among 15 MS runs with TMT10.
#' The important variables are as follows:
#'
#' \itemize{
#'   \item id
#'   \item Protein.IDs
#'   \item Only.identified.by.site
#'   \item Potential.contaminant
#'   \item Reverse
#' }
#'
#' @format A data frame with 1075 rows and 105 variables.
#' @examples
#' head(proteinGroups)
#'
"proteinGroups"

#' Example of annotation file for raw.mq
#'
#' Annotation of example data, raw.mq, in this package.
#' It should be prepared by users.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Run : MS run ID. It should be the same as Raw.file info
#'   in raw.mq
#'   \item Channel : Labeling information (channel.0, ..., channel.9).
#'   The channel index should be consistent with the channel columns in raw.mq.
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item Mixture : TMT Mixture ID. It is used to indicate the technical
#'   replicates or fractions. For example, if 'Run' 1, 2, 3 are the fractions
#'   of same TMT mixture, then they should have same 'Mixture'.
#'   \item BioReplicate : Unique ID for biological subject
#' }
#'
#' @format A data frame with 150 rows and 5 variables.
#' @examples
#' head(annotation.mq)
#'
"annotation.mq"

#' Example of output from SpectroMine for TMT6 experiments.
#'
#' Example of SpectroMine PSM sheet.
#' It is the input for SpectroMinetoMSstatsTMTFormat function, with annotation file.
#' It includes peak intensities for 10 proteins among 12 MS runs with TMT6.
#' The important variables are as follows:
#'
#' \itemize{
#'   \item PG.ProteinAccessions
#'   \item P.MoleculeID
#'   \item PP.Charge
#'   \item R.FileName
#'   \item PG.QValue
#'   \item PSM.Qvalue
#'   \item Channels : PSM.TMT6_126..Raw., ..., PSM.TMT6_131..Raw.
#' }
#'
#' @format A data frame with 170 rows and 28 variables.
#' @examples
#' head(raw.mine)
#'
"raw.mine"

#' Example of annotation file for raw.mine
#'
#' Annotation of example data, raw.mine, in this package.
#' It should be prepared by users.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Run : MS run ID. It should be the same as R.FileName info
#'   in raw.mine
#'   \item Channel : Labeling information (TMT6_126, ..., TMT6_131).
#'   The channels should be consistent with the channel columns in raw.mine.
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item Mixture : TMT Mixture ID. It is used to indicate the technical
#'   replicates or fractions. For example, if 'Run' 1, 2, 3 are the fractions
#'   of same TMT mixture, then they should have same 'Mixture'.
#'   \item BioReplicate : Unique ID for biological subject
#' }
#'
#' @format A data frame with 72 rows and 5 variables.
#' @examples
#' head(annotation.mine)
#'
"annotation.mine"

#' Example of output from PDtoMSstatsTMTFormat function
#'
#' It is calculated from raw.pd and annotation.pd
#' It is the output of PDtoMSstatsTMTFormat function
#' It should includes the required columns as below.
#' The variables are as follows:
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
#'   \item Intensity: Protein Abundance
#' }
#'
#' @format A data frame with 20110 rows and 10 variables.
#' @examples
#' head(input.pd)
#'
"input.pd"

#' Example of output from protein.summarizaiton function
#'
#' It is calculated from input.pd.
#' It is the output of protein.summarization function
#' It should includes the required columns as below.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Run : MS run ID
#'   \item Protein : Protein ID
#'   \item Abundance: Protein Abundance
#'   \item Channel : Labeling information (126, ... 131)
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item BioReplicate : Unique ID for biological subject.
#'   \item Mixture : Unique ID for TMT mixture.
#' }
#'
#' @format A data frame with 100 rows and 7 variables.
#' @examples
#' head(quant.msstats)
#'
"quant.msstats"

#' Example of output from groupComparison.TMT function
#'
#' It is calculated from the result of protein.summarizaiton function with contrast.matrix = 'pairwise'.
#' It is the output of groupComparison.TMT function
#' It should includes the required columns as below.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Protein : Protein ID
#'   \item Label: Label of the pairwise comparision or contrast
#'   \item log2FC: Log2 fold change
#'   \item SE: Standard error of the comparsion of contrast results
#'   \item DF: Degree of freedom
#'   \item pvalue: Value of p statistic of the test
#'   \item adj.pvalue: adjusted p value
#' }
#'
#' @format A data frame with 60 rows and 7 variables.
#' @examples
#' head(test.pairwise)
#'
"test.pairwise"
