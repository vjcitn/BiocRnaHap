load("genes_avail.rda")
library(shiny)

# the comment* strings here set up the 'about' tab
# with embedded hyperlinks

comment = "This app is a prototype of a general approach to
exploring RNA-seq-derived haplotypes.  We use the output
of the"
comment000=a(href="https://www.ncbi.nlm.nih.gov/pubmed/27605262", "phASER")
comment0= a(href="https://github.com/secastel/phaser/tree/master/phaser", "algorithm")
comment00= "as a transcriptome-wide collection
of long-range haplotypes.  A subset of these are presented to the user
in a gene-centric fashion: the user picks a gene symbol and a
radius measured in bp that defines a region around the gene body
to be checked for events of multiple SNV in phase.  A table
listing these events, if any, is presented in tab 'haps'."

comment1a = "Our objective is creation of statistics related
to expression-haplotype configurations like the following:"

comment1b = "This display depends on a certain approach to
a) determining haplotypes (SNP configurations) of interest,
b) stratifying a population using these
haplotypes, c) representing gene expression distributions
for these strata.  In this case, the haplotypes are derived
from a single individual, NA06986, and the expression and
genotype data are obtained for a subset of the cell lines
studied in GEUVADIS."

comment3 = "The haplotypes presented are those with
2-6 SNPs, and total read count supporting haplotype at least 15.
These filtering parameters are controlled via BiocRnaHaps::rnahapsNearGene."

comment4 = "At this time the haplotypes are limited to those
derived in the example dataset described in the phASER tutorial,
for HapMap participant NA06986."

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("BiocRnaHap app for exploring", a(
     href='https://github.com/secastel/phaser/tree/master/phaser',"'phASER' results")),
   helpText("See the 'about' tab for all details.  Patience is required as numerous annotation and experimental resources are assembled to combine
inferred haplotypes with expression data.  Specifically, gene symbol
processing takes time, wait until the selector box below has a
white background."),
#   uiOutput("selector"),
#   selectInput("gene", "gene", choices=genes_avail, selected="ORMDL3"),
   textInput("gene", "gene", value="ORMDL3"),
   numericInput("radius", "radius", 100000, 
     min=0, max=500000, step=50000),
   textOutput("selnum"),
   helpText("This app was developed in conjunction with the",
a(href='https://github.com/NCBI-Hackathons/Computational_Medicine_1',
"NHGRI Computational Medicine Hackathon"),"led by Ben Busby"),
   actionButton("stopBtn", "Stop app."), width=3
   ),
  mainPanel(
   tabsetPanel(
    tabPanel("haps", 
         helpText("scroll table to right for genomic coordinates"),
         DT::dataTableOutput("tab1")
       ),
#    uiOutput("VCFgrab"),
    tabPanel("VCF", verbatimTextOutput("snps")),
    tabPanel("ExprsBox", plotOutput("box")),
    tabPanel("about", helpText(comment, comment000, comment0, comment00), 
             helpText(comment1a), 
             plotOutput("demoplot"),
             helpText(comment1b), 
             helpText(comment3), 
             helpText(comment4), 
             verbatimTextOutput("sessinf"))
    )
   )
  )
 )

