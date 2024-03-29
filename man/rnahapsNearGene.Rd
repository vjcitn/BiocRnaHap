% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rnahapsNearGene}
\alias{rnahapsNearGene}
\title{produce a GRanges with RNA haplotypes near a specific gene}
\usage{
rnahapsNearGene(sym, haptab = NA06986_rnahaps,
  genes = ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75),
  radius = 10000, min_total_reads = 15, min_variants = 2,
  varnum_to_use = 2:6)
}
\arguments{
\item{sym}{character(1) gene symbol}

\item{haptab}{data.frame corresponding to phaser haplotypes.txt table}

\item{genes}{a GRanges as produced by ensembldb::genes}

\item{radius}{numeric(1) interval upstream and downstream to be included with gene region for selection of RNA haplotypes}

\item{min_total_reads}{numeric(1) lower bound on number of reads needed to retain a haplotype}

\item{min_variants}{lower bound on number of variants in the haplotype to warrant retention}

\item{varnum_to_use}{the values of the 'variants' field for which records are retained}
}
\description{
produce a GRanges with RNA haplotypes near a specific gene
}
\note{
grep() is used with fixed=TRUE to find sym in genes()$gene_name
}
\examples{
rnahapsNearGene("GSDMB")
}
