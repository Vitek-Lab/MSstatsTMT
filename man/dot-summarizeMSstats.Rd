% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_summarization_methods.R
\name{.summarizeMSstats}
\alias{.summarizeMSstats}
\title{Summarization based on MSstats}
\usage{
.summarizeMSstats(
  input,
  annotation,
  impute,
  max_quantile_censored = NULL,
  log_file_path = NULL
)
}
\arguments{
\item{input}{data.table}

\item{annotation}{data.table with run and channel annotation}

\item{impute}{only for method="msstats". TRUE (default) imputes missing 
values by Accelated failure model. FALSE uses minimum value to impute the 
missing value for each peptide precursor ion.}

\item{max_quantile_censored}{We assume missing values are censored. 
maxQuantileforCensored is Maximum quantile for deciding censored missing 
value, for instance, 0.999. Default is Null.}

\item{log_file_path}{path to a MSstats log file}
}
\value{
data.table
}
\description{
Summarization based on MSstats
}
\keyword{internal}
