% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotQTL.R
\name{plotQTL}
\alias{plotQTL}
\title{Generate plots of expression values for given eQTL and expressed heterozygote variant}
\usage{
plotQTL(modelOut, aseVar, qtlVar, otherAll = FALSE)
}
\arguments{
\item{modelOut}{The output of the Run.model function}

\item{aseVar}{The id of the allele specific expression variant.}

\item{qtlVar}{The id of the putative regulatory variant.}
}
\description{
This function takes the output of Run.model and produces various plots for the specified putative regulatory variant.
}
\examples{
#' #Plot the output of Run.model.R stored in modelOut
plotQTL(aseDat,'rs116197803', 'rs175183')
}
