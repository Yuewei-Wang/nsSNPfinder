library(shiny)
library(ggplot2)
library(r3dmol)
# This example is adapted from
# RStudio Inc. (2013). Widgets. Shiny Gallery. URL:https://shiny.rstudio.com/gallery/widgets.html
# RStudio Inc. (2013). Sliders. Shiny Gallery. URL:https://shiny.rstudio.com/gallery/sliders.html
# RDocumentation. (n.d.) r3dmol (version 0.1.2). URL: https://www.rdocumentation.org/packages/r3dmol/versions/0.1.2/topics/r3dmol-shiny
# Reference for the sequence length:
# Pagès, H. (2017). SNPlocs.Hsapiens.dbSNP144.GRCh38: SNP locations for Homo sapiens (dbSNPBuild 144). R package version 0.99.20.

ui <- fluidPage(
  # input and guide information
  titlePanel("Finding nsSNP Info for Gene of Interest"),

  sidebarLayout(

    sidebarPanel(
      tags$h4("nsSNPfinder Overview:"),

      tags$p("The app initially requires the input of chromosome and the coordinates of interest, the retrieved results from bioMart leading to preliminary check of information from genome level to protein level."),

      br(),

      tags$p("Three tabs demonstrate the query results:"),

      tags$p("1. \"Summary\" gives the table information of gene that contains the input range, list the count information of nsSNPs."),

      tags$p("2. \"nsSNPplot\" displays the nsSNPs distribution plot for input region. The pointer on nsSNPplot could guide the input for \"PDB\" tab"),

      tags$p("3. \"PDB\" displays the fisrt retrieved 3D protein structure of input gene from protein data bank. The input nsSNP-involved residue will be label in red."),

      hr(),

      tags$h4("Chromosome Sequence Lengths (bp) Reference"),

      tags$body(
        p("Chromosome 1: 248956422"),
        p("Chromosome 2: 242193529"),
        p("Chromosome 3: 198295559"),
        p("Chromosome 4: 190214555"),
        p("Chromosome 5: 181538259"),
        p("Chromosome 6: 170805979"),
        p("Chromosome 7: 159345973"),
        p("Chromosome 8: 145138636"),
        p("Chromosome 9: 138394717"),
        p("Chromosome 10: 133797422"),
        p("Chromosome 11: 135086622"),
        p("Chromosome 12: 133275309"),
        p("Chromosome 13: 114364328"),
        p("Chromosome 14: 107043718"),
        p("Chromosome 15: 101991189"),
        p("Chromosome 16: 90338345"),
        p("Chromosome 17: 83257441"),
        p("Chromosome 18: 80373285"),
        p("Chromosome 19: 58617616"),
        p("Chromosome 20: 64444167"),
        p("Chromosome 21: 46709983"),
        p("Chromosome 22: 50818468"),
        p("Chromosome X: 156040895"),
        p("Chromosome Y: 57227415"))
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Summary",
                           fluidRow(
                             column(width = 12,
                                    h4("\"Summary\" tab returns the retrived nsSNP results from bioMart"),
                                    selectInput("chr", "Choose a chromosome:", c(c(1:22),"X", "Y"), selected = 3),

                                    textInput("start", "Enter the start coordinate:", 49359145),

                                    textInput("end", "Enter the end coordinate:", 49411645),

                                    actionButton(inputId = "update1",
                                                 label = "Run Summary")
                             )
                           ),
                           hr(),
                           tableOutput("summary")
                  ),
                  tabPanel("nsSNPplot",
                           fluidRow(
                             column(width = 12,
                                    h4("\"nsSNPplot\" tab displays the nsSNPs distribution of input region"),
                                    br(),
                                    helpText("NOTE: for \"nsSNPplot\" function, please limit the input range within 200-length for better resolution."),
                                    selectInput("chr2", "Choose a chromosome:", c(c(1:22),"X", "Y"), selected = 3),

                                    textInput("start2", "Enter the start coordinate:", 49395500),

                                    textInput("end2", "Enter the end coordinate:", 49395600),

                                    actionButton(inputId = "update2",
                                                 label = "Run nsSNPplot")
                             )
                           ),
                           hr(),
                           plotOutput("nssnp",
                                      click = clickOpts(id = "click")),
                           fluidRow(
                             column(width = 12,
                                    h4("nsSNP Near the Click: please click the center of label"),
                                    verbatimTextOutput("pointInfo")
                             )
                           )
                  ),
                  tabPanel("PDB",
                           fluidRow(
                             column(width = 12,
                                    h4("\"PDB\" tab displays the fisrt retrieved 3D protein structure"),

                                    selectInput("pdbChr", "Choose a chromosome:", c(c(1:22),"X", "Y"), selected = 3),

                                    textInput("gene", "Enter gene name:", "RHOA"),

                                    textInput("pos", "Enter the nsSNP location:", 49395565),

                                    actionButton(inputId = "update3",
                                                 label = "Run PDB (default choice)")
                             )
                           ),
                           hr(),
                           r3dmolOutput("pdb"),
                           fluidRow(
                             column(width = 12,
                                    h4("PDB Model List:"),
                                    tags$strong("The list of all PDB IDs of involved models:"),
                                    tags$p("PDB choice could be custmized by runing *displayPDB*."),
                                    verbatimTextOutput("listInfo"))
                           )
                  )
      )
    )
  )
)


server <- function(input, output, session) {

  snp <- eventReactive(input$update1,{
    nsSNPCalculatebyRange(as.numeric(input$chr),
                          as.numeric(input$start),
                          as.numeric(input$end))
  })

  p <- eventReactive(input$update2,{
    nsSNPFreqPlot(as.numeric(input$chr2),
                  as.numeric(input$start2),
                  as.numeric(input$end2))
  })

  defaultPDB <- eventReactive(input$update3,{
    displayPDB(as.integer(input$pdbChr),
               input$gene,
               as.integer(input$pos))
  })

  output$ref <- renderTable({
    chrVars <- c(c(1:22),"X", "Y")
    coordVars <- c("248956422", "242193529", "198295559", "190214555", "181538259",
                   "170805979", "159345973", "145138636", "138394717", "133797422",
                   "135086622", "133275309", "114364328", "107043718", "101991189",
                   "90338345", "83257441", "80373285", "58617616", "64444167",
                   "46709983", "50818468", "156040895", "57227415")
    data<-data.frame(chromosome = chrVars, length = coordVars)
    data
  })

  output$summary <- renderTable ({
    if (! is.null(snp))
      snp()
  })

  output$nssnp <- renderPlot({
    if (! is.null(p))
      p()
  })

  output$pointInfo <- renderPrint({
    vars <- data.frame(loc = c("loc"), numVars = c("numVars"))
    nearPoints(vars, input$click)
    cat("nsSNP location:\n")
    str(input$click$x)
  })

  output$listInfo <- renderPrint({
    if (! is.null(defaultPDB))
      capture.output(data<-defaultPDB())

  })
  output$pdb <- renderR3dmol({
    if (! is.null(defaultPDB))
      defaultPDB()
  })

}
shiny::shinyApp(ui = ui, server = server)
# [END]
