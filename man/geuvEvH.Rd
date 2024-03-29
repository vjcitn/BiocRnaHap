% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{geuvEvH}
\alias{geuvEvH}
\title{compute a data.frame with GEUVADIS expression for a selected gene
and associated haplotypes formed using entries in snps}
\usage{
geuvEvH(sym = "GSDMB", snps = c("rs2305480", "rs2305479"),
  vcf = NULL, ...)
}
\arguments{
\item{sym}{a character(1) gene symbol}

\item{snps}{a vector of SNP identifiers}

\item{vcf}{a CompressedVcf instance from VariantAnnotation}

\item{\dots}{not currently used}
}
\description{
compute a data.frame with GEUVADIS expression for a selected gene
and associated haplotypes formed using entries in snps
}
\note{
it is assumed that the entries in `snps` are elements of
a 'long-range' RNA-based haplotype as
inferred by phaser
}
\examples{
requireNamespace("geuvPack")
if (!exists("geuFPKM")) {
  require(geuvPack)
  data(geuFPKM)
  }
tab = geuvEvH(vcf=BiocRnaHap::gsd_vcf)
head(tab)
require(ggplot2)
ggplot(tab, aes(x=haplotypes, y=exprs, colour=popcode)) + 
  geom_boxplot(aes(group=haplotypes), outlier.size=0) + 
  geom_jitter() 
}
