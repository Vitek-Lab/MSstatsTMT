\name{MSstatsTMTnews}
\title{MSstatsTMT News}
\encoding{UTF-8}

\section{Changes in version 2.14.1 (2024-10-30) }{\itemize{
  \item Added interactive plotting functions using Plotly.
  \item Integrated support for ProteinProspector data conversion via MSstatsConvert.
}}

\section{Changes in version 2.6.1 (2023-02-26) }{\itemize{
  \item Minor change: update PhilosophertoMSstatsTMTFormat function
  normalization
}}

\section{Changes in version 2.4.1 (2022-09-12) }{\itemize{
  \item Minor change: fix the bug of data table in performing local
  normalization
}}

\section{Changes in version 2.2.7 (2022-02-18)}{\itemize{
  \item Minor change: extend PhilosophertoMSstatsTMTFormat function to have multiple types of input
}}

\section{Changes in version 2.2.6 (2022-02-14)}{\itemize{
  \item Major change: add PhilosophertoMSstatsTMTFormat function as converter for outputs from Philosopher
}}

\section{Changes in version 2.2.5 (2021-10-25)}{\itemize{
  \item Minor change: add different point shape to dataProcessPlotsTMT as indicator of imputed values
}}

\section{Changes in version 2.2.3 (2021-10-06)}{\itemize{
  \item Minor change: fix the bug when df.prior is infinite
}}

\section{Changes in version 2.2.0 (2021-09-21)}{\itemize{
  \item Major change: extend groupComparisonTMT() function to cover repeated measures design
  \item Allow flexible order of condition in dataProcessPlotsTMT.
  \item Fix bug in Condition label in dataProcessPlotsTMT.
  \item Improve MSstatsTestSingleProteinTMT() by directly reading lmerTest output. This may make statistics slightly different due to different numeric accuracy
  \item fix bug when condition name contains 'group'
  \item change the x-axis order in profile plot
}}

\section{Changes in version 2.0.1 (2021-06-14)}{\itemize{
  \item update comments of PD converter function
  \item fix bug in proteinSummarization() function when MBimpute = F}}

\section{Changes in version 2.0.0 (2021-05-14)}{\itemize{
  \item Refactor the pacakge to make it modulized
  \item parameter rmPSM_withMissing_withinRun was retired
  \item change converter function output from data.frame to MSstatsValidated
  \item change proteinSummarization() function output from data.frame to a list with two data.frames: FeatureLevelData and ProteinLevelData
  \item change dataProcessPlotsTMT() function input from two data.frames to a list with two data.frames, which is the output of function proteinSummarization()
  \item change groupComparisonTMT() function input from a data.frame to a list with two data.frames, which is the output of function proteinSummarization()
  \item change groupComparisonTMT() function output from data.frame to a list with three items, ComparisonResult, FittedModel and ModelQC
}}

\section{Changes in version 1.8.2 (2020-12-17)}{\itemize{
  \item Update progress bar
  \item Update groupComparisonTMT to avoid reusing the local function copied from lmer pacakge}}

\section{Changes in version 1.8.1 (2020-12-10)}{\itemize{
  \item Add citation of the MSstatsTMT paper
  \item Fix the bug in groupComparisonTMT() due to the update of dependent pacakge
  \item Fix the bug in MedianPolish summarization
  \item proteinSummarization(): replace the zero values with NA before and after peptide normalization}}

\section{Changes in version 1.6.6 (2020-10-13)}{\itemize{
  \item Fix the bug in converters due to fractions with same mean, sum and max values
  \item Fix the bug in converters due to summaryforMultipleRows
  \item Fix the bug in OpemMS converter due to duplicated rows}}

\section{Changes in version 1.6.3 (2020-06-05)}{\itemize{
  \item Allow NA in the annotation file}}

\section{Changes in version 1.6.2 (2020-06-02)}{\itemize{
  \item Fix the bug in proteinSummarization() function and make sure the input to dataProcess is data.frame}}

\section{Changes in version 1.6.1 (2020-05-10)}{\itemize{
  \item Update groupComparisonTMT() to make predictions for every protein}}

\section{Changes in version 1.4.7 (2020-04-24)}{\itemize{
  \item Update the NEWS file}}

\section{Changes in version 1.4.6 (2020-04-14)}{\itemize{
  \item Fix bug in groupComparison() for unbalanced design
  \item Use df approximation from lmerTest to perform group comparison}}

\section{Changes in version 1.4.5 (2020-03-01)}{\itemize{
  \item Add new function OpenMStoMSstatsTMTFormat()}}

\section{Changes in version 1.4.4 (2020-02-01)}{\itemize{
  \item Fix bug in PDtoMSstatsTMTFormat(): remove redundant rows when combining multiple fractions
}}

\section{Changes in version 1.4.3 (2019-12-28)}{\itemize{
  \item Fix bug in groupComparisonTMT(): very few measurements case previously doesn't work
}}
\section{Changes in version 1.4.2 (2019-12-20)}{\itemize{
\item Add the column 'issue' to the output of groupComparisonTMT()
}}
\section{Changes in version 1.4.1 (2019-10-31)}{\itemize{
\item Fix the bug in the PDtoMSstatsTMTFormat() due to different PD version
}}
\section{Changes in version 1.2.7 (2019-08-22)}{\itemize{
\item Update normalization options in proteinSummarization() function to 
include global protein normalization and local peptide normalization with respect to reference channel
\item Fix the bug in the contrast comparison of groupComparisonTMT() function.
}}
\section{Changes in version 1.2.2 (2019-06-03)}{\itemize{
\item Update the linear modeling in groupComparisonTMT() function. Implement 5 linear models.
\item Update the format of annotation file. If the channel has no sample, add 'Empty' under 
condition and BioReplicate columns.
}}
\section{Changes in version 1.2.0 (2019-05-03)}{\itemize{
\item Fix bugs in groupComparisonTMT() when using lm() function
\item Update the format of annotation file to include fraction column
\item Remove the 'fraction' option from coverter functions
\item Update the linear modeling in groupComparisonTMT() function
}}
\section{Changes in version 1.0.0 (2018-09-21)}{\itemize{
\item Submitted to Bioconductor
}}
