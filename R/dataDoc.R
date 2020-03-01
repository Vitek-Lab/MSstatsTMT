
#' Example of output from Proteome Discoverer 2.2 for TMT-10plex experiments.
#'
#' Example of Proteome discover PSM sheet.
#' It is the input for PDtoMSstatsTMTFormat function, with annotation file.
#' Annotation file should be made by users.
#' It includes peak intensities for 10 proteins
#' among 15 MS runs with TMT-10plex.
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

#' Example of annotation file for raw.pd,
#' which is the PSM output of Proteome Discoverer
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
#'   \item Mixture : Mixture of samples labeled with different TMT reagents, which can be analyzed in
#'   a single mass spectrometry experiment. If the channal doesn't have sample, please add `Empty' under Condition.
#'   \item TechRepMixture : Technical replicate of one mixture. One mixture may have multiple technical replicates.
#'   For example, if `TechRepMixture' = 1, 2 are the two technical replicates of one mixture, then they should match
#'   with same `Mixture' value.
#'   \item Fraction : Fraction ID. One technical replicate of one mixture may be fractionated into multiple fractions to increase the analytical depth.
#'   Then one technical replicate of one mixture should correspond to multuple fractions.
#'   For example, if `Fraction' = 1, 2, 3 are three fractions of the first technical replicate of one TMT mixture of biological subjects,
#'   then they should have same `TechRepMixture'  and `Mixture' value.
#'   \item BioReplicate : Unique ID for biological subject. If the channal doesn't have sample, please add `Empty' under BioReplicate.
#' }
#'
#' @format A data frame with 150 rows and 7 variables.
#' @examples
#' head(annotation.pd)
#'
"annotation.pd"

#' Example of output from MaxQuant for TMT-10plex experiments.
#'
#' Example of evidence.txt from MaxQuant.
#' It is the input for MaxQtoMSstatsTMTFormat function, with proteinGroups.txt
#' and annotation file. Annotation file should be made by users.
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

#' Example of proteinGroups file from MaxQuant for TMT-10plex experiments.
#'
#' Example of proteinGroup.txt file from MaxQuant,
#' which is identified protein group information file.
#' It is the input for MaxQtoMSstatsTMTFormat function, with evidence.txt
#' and annotation file.
#' It includes identified protein groups for 10 proteins
#' among 15 MS runs with TMT10.
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

#' Example of annotation file for evidence, which is the output of MaxQuant.
#'
#' Annotation of example data, evidence, in this package.
#' It should be prepared by users.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Run : MS run ID. It should be the same as Raw.file info
#'   in raw.mq
#'   \item Channel : Labeling information (channel.0, ..., channel.9).
#'   The channel index should be consistent with the channel columns in raw.mq.
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item Mixture : Mixture of samples labeled with different TMT reagents, which can be analyzed in
#'   a single mass spectrometry experiment. If the channal doesn't have sample, please add `Empty' under Condition.
#'   \item TechRepMixture : Technical replicate of one mixture. One mixture may have multiple technical replicates.
#'   For example, if `TechRepMixture' = 1, 2 are the two technical replicates of one mixture, then they should match
#'   with same `Mixture' value.
#'   \item Fraction : Fraction ID. One technical replicate of one mixture may be fractionated into multiple fractions to increase the analytical depth.
#'   Then one technical replicate of one mixture should correspond to multuple fractions.
#'   For example, if `Fraction' = 1, 2, 3 are three fractions of the first technical replicate of one TMT mixture of biological subjects,
#'   then they should have same `TechRepMixture'  and `Mixture' value.
#'   \item BioReplicate : Unique ID for biological subject. If the channal doesn't have sample, please add `Empty' under BioReplicate.
#' }
#'
#' @format A data frame with 150 rows and 7 variables.
#' @examples
#' head(annotation.mq)
#'
"annotation.mq"

#' Example of output from SpectroMine for TMT-6plex experiments.
#'
#' Example of SpectroMine PSM sheet.
#' It is the output of SpectroMine and the input for SpectroMinetoMSstatsTMTFormat function,
#' with annotation file.
#' Annotation file should be made by users.
#' It includes peak intensities for 10 proteins among 12 MS runs with TMT-6plex.
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

#' Example of annotation file for raw.mine, which is the output of SpectroMine.
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
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0). If the channal doesn't have sample, please add 'Empty' under Condition.
#'   \item Mixture : Mixture of samples labeled with different TMT reagents, which can be analyzed in
#'   a single mass spectrometry experiment.
#'   \item TechRepMixture : Technical replicate of one mixture. One mixture may have multiple technical replicates.
#'   For example, if 'TechRepMixture' = 1, 2 are the two technical replicates of one mixture, then they should match
#'   with same 'Mixture' value.
#'   \item Fraction : Fraction ID. One technical replicate of one mixture may be fractionated into multiple fractions to increase the analytical depth.
#'   Then one technical replicate of one mixture should correspond to multuple fractions.
#'   For example, if 'Fraction' = 1, 2, 3 are three fractions of the first technical replicate of one TMT mixture of biological subjects,
#'   then they should have same 'TechRepMixture'  and 'Mixture' value.
#'   \item BioReplicate : Unique ID for biological subject. If the channal doesn't have sample, please add 'Empty' under BioReplicate
#' }
#'
#' @format A data frame with 72 rows and 7 variables.
#' @examples
#' head(annotation.mine)
#'
"annotation.mine"


#' Example of MSstatsTMT report from OpenMS for TMT-10plex experiments.
#'
#' Example of MSstatsTMT PSM sheet from MaxQuant.
#' It is the input for OpenMStoMSstatsTMTFormat function.
#' It includes peak intensities for 10 proteins among 27 MS runs from three TMT10 mixtures.
#' The important variables are as follows:
#'
#' \itemize{
#'   \item RetentionTime
#'   \item ProteinName
#'   \item PeptideSequence
#'   \item Charge
#'   \item Channel
#'   \item Condition
#'   \item BioReplicate
#'   \item Run
#'   \item Mixture
#'   \item TechRepMixture
#'   \item Fraction
#'   \item Intensity
#'   \item Reference
#' }
#'
#' @format A data frame with 860 rows and 13 variables.
#' @examples
#' head(raw.om)
#'
"raw.om"

#' Example of output from PDtoMSstatsTMTFormat function
#'
#' It is made from \code{\link{raw.pd}} and \code{\link{annotation.pd}},
#' which is the output of PDtoMSstatsTMTFormat function.
#' It should include the required columns as below.
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
#' It should include the required columns as below.
#' The variables are as follows:
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
#' head(quant.pd.msstats)
#'
"quant.pd.msstats"

#' Example of output from groupComparisonTMT function
#'
#' It is the output of groupComparisonTMT function,
#' which is the result of group comparions with the output of proteinSummarization function.
#' It should include the columns as below.
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
#' head(test.pairwise)
#'
"test.pairwise"
