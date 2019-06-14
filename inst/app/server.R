library(shiny)
library(geuvPack)
if (!exists("geuFPKM")) data(geuFPKM)
library(ggplot2)
library(BiocRnaHap)
library(magrittr)
library(dplyr)
    
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
observe({
 filter_method = "overlaps"
# filter_method = "nearest"
 filter_radius = 10000
 
 genegr = ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
 allg = sort(genegr$gene_name)
 varnum_to_use = 2:6
 haptab = NA06986_rnahaps
 haptab_gr = GRanges(haptab$contig, IRanges(haptab$start, haptab$stop))
 mcols(haptab_gr) = haptab[,-c(1:4)]
 haptab_gr_filt = haptab_gr[haptab_gr$variants %in% varnum_to_use]
 haptab_gr_filt = haptab_gr_filt[ grep("^rs", haptab_gr_filt$variant_ids)]
 if (filter_method == "nearest")
    gn = genegr[nearest(haptab_gr_filt, genegr)]
 else gn = subsetByOverlaps(genegr, haptab_gr_filt+filter_radius)
 print(length(gn))
 chs = sort(gn$gene_name)
 updateSelectInput(session, "gene", "gene", choices=chs, selected="GSDMB")
 })
 
 goodtab = reactive({
  nrtab_gr = rnahapsNearGene(input$gene, radius=input$radius)
 })
#
 output$tab1 = DT::renderDataTable({
  nrtab_gr = goodtab()
  cbind(gene=input$gene, radius=input$radius, as.data.frame(nrtab_gr))
  })
 output$snps = renderPrint({
   tab = goodtab()
   tab = as.data.frame(tab)
   pick = as.character(tab[ input$tab1_rows_selected, "variant_ids" ])
   snvec = strsplit(pick, ",")
   look1kg(snvec[[1]])
  })
 output$selnum = renderText({
  validate(need(input$tab1_rows_selected, "pick a row"))
  input$tab1_rows_selected
 })


 observeEvent(input$tab1_rows_selected, {
   tab = rnahapsNearGene(input$gene, radius=input$radius)
   tab = as.data.frame(tab)
#print(tab)
   pick = as.character(tab[ input$tab1_rows_selected, "variant_ids" ])
#print(pick)
   snvec = strsplit(pick, ",")
#print(snvec)
   myvcf = look1kg(snvec[[1]])
   glkg = geno(myvcf)$GT
   mm = geuFPKM[, intersect(colnames(glkg), colnames(geuFPKM))]
   elem = apply(glkg,2,paste, collapse=":")
   gind = grep(input$gene, rowData(geuFPKM)$gene_name)[1]
   validate(need(length(gind)>0, 
      "please pick a gene assayed in GEUVADIS"))
   quant = as.numeric(log(assay(mm[gind,])+1))
   newdf = data.frame(quant=quant, hap=elem[colnames(mm)], 
        stringsAsFactors=FALSE)
   squant = split(quant, elem[colnames(mm)])
   mor = sapply(squant, median)
   hc = newdf %>% dplyr::select(hap) %>% group_by(hap) %>% 
        summarise(n=n())
   okh = hc[hc$n > 2, 1]
   newdf = newdf[newdf$hap %in% okh$hap,]
#   print(head(newdf))
#save(newdf, file="newdf.rda")
   gg = ggplot(newdf, aes(x=factor(hap), y=quant, colour=factor(hap))) + 
     geom_boxplot() + geom_jitter() + 
     ggtitle(paste(input$gene, "expression")) +
     xlab(paste0(snvec, collapse=", "))
#
   pdf("abc.pdf")
   gg
   dev.off()
   output$box = renderPlot(gg)
   })


 }

#runApp(list(ui=ui, server=server))
  
