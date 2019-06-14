ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("BiocRnaHap app for exploring 'phaser' results"),
   uiOutput("selector"),
   selectInput("gene", "gene", choices="GSDMB"),
   numericInput("radius", "radius", 10000, 
     min=0, max=500000, step=10000),
   textOutput("selnum")
   ),
  mainPanel(
   tabsetPanel(
    tabPanel("haps", DT::dataTableOutput("tab1")),
    tabPanel("picks", verbatimTextOutput("snps")),
    tabPanel("box", plotOutput("box"))
    )
   )
  )
 )

