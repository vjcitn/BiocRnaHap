# ad hoc code here to 'speed up' app performance

library(shiny)
load("geuFPKM.rda")
library(ggplot2)
library(GenomicRanges)
library(VariantAnnotation)
source("look1kg.R")
library(magrittr)
library(dplyr)


load("genegr.rda")
allg = genegr$gene_name
load("genes_avail.rda")

# in this section we filter down the phASER input table

varnum_to_use = 2:6
haptab = BiocRnaHap::NA06986_rnahaps
haptab_gr = GenomicRanges::GRanges(haptab$contig, 
              IRanges::IRanges(haptab$start, haptab$stop))
mcols(haptab_gr) = haptab[,-c(1:4)]
haptab_gr_filt = haptab_gr[haptab_gr$variants %in% varnum_to_use]
haptab_gr_filt = haptab_gr_filt[ grep("^rs", haptab_gr_filt$variant_ids)]
    
server = function(input, output, session) {
#
# to consistently populate the selectInput for gene identifier,
# we need to use an updateSelectInput, otherwise the initial
# selection is dropped
#
# here we build from genes() and haptab the collection
# of gene regions nearest to RNA haplotypes -- other
# methods of filtering genes would be reasonable -- we
# just don't want to have to deal with them all
#
#observe({
# showNotification("setting up gene list", id="setupg")
# filter_method = "overlaps"
## filter_method = "nearest"
# filter_radius = 10000
# 
# genegr = ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
# allg = sort(genegr$gene_name)
# varnum_to_use = 2:6
# haptab = NA06986_rnahaps
# haptab_gr = GRanges(haptab$contig, IRanges(haptab$start, haptab$stop))
# mcols(haptab_gr) = haptab[,-c(1:4)]
# haptab_gr_filt = haptab_gr[haptab_gr$variants %in% varnum_to_use]
# haptab_gr_filt = haptab_gr_filt[ grep("^rs", haptab_gr_filt$variant_ids)]
# if (filter_method == "nearest")
#    gn = genegr[nearest(haptab_gr_filt, genegr)]
# else gn = subsetByOverlaps(genegr, haptab_gr_filt+filter_radius)
# #print(length(gn))
# chs = sort(gn$gene_name)
# removeNotification(id="setupg")
# updateSelectInput(session, "gene", "gene", choices=chs, selected="GSDMB")
# })
 
 goodtab = reactive({
  validate(need(input$gene %in% genes_avail, "enter a valid gene symbol"))
  ind = grep(input$gene, genegr$gene_name, fixed=TRUE)
  if (length(ind)>1) {
     warning("multiple instances of symbol found, using first")
     showNotification("multiple instances of symbol found, using first")
     }
  nrtab_gr = subsetByOverlaps(haptab_gr_filt, genegr[ind]+input$radius)
 })
#
 output$tab1 = DT::renderDataTable({
  nrtab_gr = goodtab()
  validate(need(length(nrtab_gr)>0, "increase radius"))
  basic = cbind(gene=input$gene, radius=input$radius, as.data.frame(nrtab_gr))
  basic[,c("gene", "radius", "variant_ids", "variant_alleles", "reads_total",
    "reads_hap_a", "reads_hap_b", "seqnames", "start", "end", "strand")]
  })
 
 output$vcfGrab = renderUI({
   validate(need(input$tab1_rows_selected, "Click on a row"))
   tabPanel("VCF", verbatimTextOutput("snps"))
   })
 output$snps = renderPrint({
   validate(need(input$tab1_rows_selected, "Click on a row in haps tab"))
   showNotification("interrogating VCF", id="chkvcf")
   tab = goodtab()
   tab = as.data.frame(tab)
   pick = as.character(tab[ input$tab1_rows_selected, "variant_ids" ])
   snvec = strsplit(pick, ",")
   removeNotification(id="chkvcf")
   look1kg(snvec[[1]])
  })
 output$selnum = renderText({
  validate(need(input$tab1_rows_selected, "Click on a row"))
  input$tab1_rows_selected
 })


 observeEvent(input$tab1_rows_selected, {
#   tab = BiocRnaHap::rnahapsNearGene(input$gene, radius=input$radius,
#    haptab=BiocRnaHap::NA06986_rnahaps)
   tab = goodtab()
   tab = as.data.frame(tab)
#print(tab)
   pick = as.character(tab[ input$tab1_rows_selected, "variant_ids" ])
#print(pick)
   snvec = strsplit(pick, ",")
#print(snvec)
   myvcf = look1kg(snvec[[1]])
   oksn = names(rowRanges(myvcf))
print(oksn)
   glkg = geno(myvcf)$GT
   mm = geuFPKM[, intersect(colnames(glkg), colnames(geuFPKM))]
   elem = apply(glkg,2,paste, collapse=":")
   gind = grep(input$gene, rowData(geuFPKM)$gene_name, fixed=TRUE)[1]
   validate(need(length(gind)>0, 
      "please pick a gene assayed in GEUVADIS"))
   quant = as.numeric(log(assay(mm[gind,])+1))
   newdf = data.frame(quant=quant, hap=elem[colnames(mm)], pop=mm$popcode,
        stringsAsFactors=FALSE)
   squant = split(quant, elem[colnames(mm)])
   mor = sapply(squant, median)
   hc = newdf %>% dplyr::select(hap) %>% group_by(hap) %>% 
        summarise(n=n())
   okh = hc[hc$n > 2, 1]
   newdf = newdf[newdf$hap %in% okh$hap,]
#   print(head(newdf))
#save(newdf, file="newdf.rda")
   gg = ggplot(newdf, aes(x=factor(hap), y=quant)) + 
     geom_boxplot(data=newdf, aes(x=factor(hap)), outlier.size=0) + geom_jitter(aes(colour=pop)) + 
     ggtitle(paste(input$gene, "expression in GEUVADIS")) +
     xlab(paste0(snvec, collapse=", ")) + 
       theme(axis.text.x=element_text(angle=45, hjust=1))
   output$box = renderPlot(gg)
   })
 output$sessinf = renderPrint({
   sessionInfo()
   })

# handle stop button
     observeEvent(input$stopBtn, {
       ans = NULL
       if (length(input$gene)>0) {
           nrtab_gr = goodtab()
           ans = cbind(gene=input$gene, 
              radius=input$radius, as.data.frame(nrtab_gr))
        }   
       stopApp(returnValue=ans)
       })  
  output$demoplot = renderPlot({  plotDemo() })



 }

#runApp(list(ui=ui, server=server))
  
