#' use ensembl REST API to acquire SNP location for an rsid, lifting over as necessary
#' @importFrom rtracklayer liftOver
#' @importFrom httr GET content
#' @importFrom rjson fromJSON
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @param rsid character(1)
#' @param ch a liftOver chain object as imported with rtracklayer
#' @param chk38 logical(1) to check that the query was resolved with assembly_name 'GRCh38'
#' @param target_genome character(1) just used for tag on result
#' @examples
#' ensRsLoc()
#' @export
ensRsLoc = function(rsid="rs6060535", ch=NULL, chk38=FALSE,
   target_genome="hg19") {
 stopifnot(length(rsid)==1)
 req1 = sprintf("https://rest.ensembl.org/variation/human/%s?content-type=application/json", rsid)
 info = fromJSON(httr::content(GET(req1), type="text/plain"))
 rng = GRanges(paste0("chr", info$mappings[[1]]$seq_region_name),
   IRanges(info$mappings[[1]]$start, width=1))
 assy = info$mappings[[1]]$assembly_name
 if (chk38) stopifnot(assy == "GRCh38")
 if (!is.null(ch)) rng = liftOver(rng, ch)[[1]]
 genome(rng) = target_genome
 rng
}
 
