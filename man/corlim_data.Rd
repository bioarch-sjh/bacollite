% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corlim_data.R
\name{corlim_data}
\alias{corlim_data}
\title{Take the output of ms_fit and convert to hits for a range of correlation thresholds}
\usage{
corlim_data(fn_hits, fn_laglim, fn_corlim = seq(0, 1, 0.05))
}
\arguments{
\item{fn_hits}{the dataframe returned by ms_fit}

\item{fn_laglim}{the max permissable lag}

\item{fn_corlim}{the set of correlation thresholds}
}
\value{
a list of dataframes, one for each replicate, each containing the following fields
  \item{cl}{The correlation limit}
  \item{nh}{The number of hits above the limit}
  \item{sc}{The sum of ion counts for hits}
}
\description{
Take the output of ms_fit and convert to hits for a range of correlation thresholds
}
