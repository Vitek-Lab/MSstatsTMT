\name{MSstatsTMTnews}
\title{MSstatsTMT News}
\encoding{UTF-8}

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
