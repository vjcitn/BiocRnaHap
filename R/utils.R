#library(shiny)
#library(geuvPack)
#if (!exists("geuFPKM")) data(geuFPKM)
#library(ggplot2)
#library(BiocRnaHap)
#library(magrittr)
#library(dplyr)
#
#genegr = ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
#allg = sort(genegr$gene_name)

# varnum_to_use = 2:6
# haptab = NA06986_rnahaps

#' produce a filtered RNA-haplotype table as GRanges, limiting
#' to SNP-sets of specific cardinality and providing a set
#' of gene ranges 'nearest' to the initial SNP
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @param haptab data.frame corresponding to phaser haplotypes.txt table
#' @param genegr a GRanges with gene ranges
#' @param varnum_to_use the values of the 'variants' field for which records are retained
#' @examples
#' head(filter_haptab())
#' @export
filter_haptab = function(haptab=BiocRnaHap::NA06986_rnahaps,
   genegr = ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75),
   varnum_to_use=2:6) {
 haptab_gr = GenomicRanges::GRanges(haptab$contig, 
              IRanges::IRanges(haptab$start, haptab$stop))
 mcols(haptab_gr) = haptab[,-c(1:4)]
 haptab_gr_filt = haptab_gr[haptab_gr$variants %in% varnum_to_use]
 haptab_gr_filt = haptab_gr_filt[ grep("^rs", haptab_gr_filt$variant_ids)]
 gn = genegr[nearest(haptab_gr_filt, genegr)]
 gn = gn[which(nchar(gn$gene_name)>0)]
 list(haptab_filt = haptab_gr_filt, genes_near=gn)
}

#' produce a GRanges with RNA haplotypes near a specific gene
#' @param sym character(1) gene symbol
#' @param haptab data.frame corresponding to phaser haplotypes.txt table
#' @param genes a GRanges as produced by ensembldb::genes
#' @param radius numeric(1) interval upstream and downstream to be included with gene region for selection of RNA haplotypes
#' @param min_total_reads numeric(1) lower bound on number of reads needed to retain a haplotype
#' @param min_variants lower bound on number of variants in the haplotype to warrant retention
#' @param varnum_to_use the values of the 'variants' field for which records are retained
#' @examples
#' rnahapsNearGene("GSDMB")
#' @export
rnahapsNearGene = function(sym, haptab=NA06986_rnahaps,
       genes=ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75),
       radius=10000, min_total_reads=15, min_variants=2, 
       varnum_to_use = 2:6) {
 ind = grep(sym,genes$gene_name)
 if (length(ind)==0) stop("sym not found")
 if (length(ind)>1) {
    warning("multiple representations of symbol, using first")
    ind = ind[1]
    }
 haptab_gr_filt = filter_haptab(haptab)$haptab_filt
# haptab_gr = GRanges(haptab$contig, IRanges(haptab$start, haptab$start))
# mcols(haptab_gr) = haptab[,-c(1:4)]
# haptab_gr_filt = haptab_gr[haptab_gr$variants %in% varnum_to_use]
# haptab_gr_filt = haptab_gr_filt[ grep("^rs", haptab_gr_filt$variant_ids)]
 chk = subsetByOverlaps(haptab_gr_filt, genes[ind]+radius)
 chk[which(chk$reads_total>=min_total_reads & chk$variants >= min_variants)]
}
