% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{HapTab}
\alias{HapTab}
\title{encapsulate phaser output text tables}
\usage{
HapTab(datfr, source)
}
\arguments{
\item{datfr}{data.frame instance importing 'haplotypes.txt'}

\item{source}{character(1) tag characterizing origin}
}
\description{
encapsulate phaser output text tables
}
\examples{
HapTab(NA06986_rnahaps, "phaser NA06986")
}
