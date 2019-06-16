look1kg = function (rsids = c("rs71503074", "rs71503075", "rs12084944"), 
    ch = BiocRnaHap::ch38to19, src = "http://1000genomes.s3.amazonaws.com/release/20130502/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", 
    rsonly = TRUE, snvonly = TRUE) 
{
    isrs = grep("^rs", rsids[1])
    if (length(isrs) == 1) {
        suppressMessages({
            lkups = lapply(rsids, function(x) BiocRnaHap::ensRsLoc(x, ch = ch))
            mygr = unlist(GRangesList(lkups))
        })
    }
    else {
        ldat = strsplit(rsids, "_")
        chr = sapply(ldat, "[", 2)
        loc = sapply(ldat, "[", 3)
        mygr = GRanges(chr, IRanges(loc, width = 1))
    }
    sn = unique(seqnames(mygr))
    if (length(sn) > 1) 
        stop("SNPs found to lie on more than 1 chromosome")
    src = sprintf(src, as.character(sn))
    seqlevels(mygr) = gsub("chr", "", seqlevels(mygr))
    p = ScanVcfParam(which = mygr)
    ans = readVcf(src, param = p)
    gr = rowRanges(ans)
    if (rsonly) {
        ok = grep("^rs", names(gr))
        ans = ans[ok, ]
    }
    if (snvonly) {
        indels = which(lengths(gr$REF) != 1)
        if (length(indels) > 0) 
            ans = ans[-indels, ]
    }
    ans
}
