#' query 1000 genomes genotypes at pairs of SNP
#' @importFrom VariantAnnotation ScanVcfParam readVcf
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom GenomeInfoDb seqlevels seqnames genome<- seqlevels<-
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics path
#' @importFrom DelayedArray rowRanges
#' @param rsids character vector
#' @param ch chain object from rtracklayer
#' @param src character(1) sprintf template for 1000 genomes vcf, chr[nn] to be supplied for pct-s
#' @param rsonly logical(1), exclude 'esv' records
#' @examples
#' lk2 = look1kg(c("rs71503074","rs71503075","rs12084944"))
#' table(geno(lk2)$GT[1,], geno(lk2)$GT[2,], geno(lk2)$GT[3,])
#' @export
look1kg = function(rsids=c("rs71503074","rs71503075","rs12084944"),
    ch = BiocRnaHap::ch38to19,
    src=
    "http://1000genomes.s3.amazonaws.com/release/20130502/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    rsonly=TRUE) {
  suppressMessages({
    lkups = lapply(rsids, function(x) ensRsLoc(x, ch=ch))
  })
  mygr = unlist(GRangesList(lkups))
  sn = unique(seqnames(mygr))
  if (length(sn)>1) stop("SNPs found to lie on more than 1 chromosome")
  src = sprintf(src, as.character(sn))
  seqlevels(mygr) = gsub("chr", "", seqlevels(mygr))
  p = ScanVcfParam(which=mygr)
  ans = readVcf(src, param=p)
  if (rsonly) {
    gr = rowRanges(ans)
    ok = grep("^rs", names(gr))
    ans = ans[ok,]
  }
  ans
}
