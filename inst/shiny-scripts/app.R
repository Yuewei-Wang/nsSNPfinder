library(shiny)
library(ggplot2)
# This example is adapted from
# RStudio Inc. (2013). Widgets. Shiny Gallery. URL:https://shiny.rstudio.com/gallery/widgets.html
# RStudio Inc. (2013). Sliders. Shiny Gallery. URL:https://shiny.rstudio.com/gallery/sliders.html
# Reference for the sequence length:
# Pagès, H. (2017). SNPlocs.Hsapiens.dbSNP144.GRCh38: SNP locations for Homo sapiens (dbSNPBuild 144). R package version 0.99.20.

ui <- fluidPage(
  # input and guide information
  titlePanel("Input Gene of Interest"),

  sidebarLayout(

    sidebarPanel(

      tags$p("The app requires the input of chromosome name and the coordinates of gene of interest."),

      tags$p("The distribution of nsSNP locations within the input range will be returned."),

      helpText("NOTE: for \"nsSNPplot\" function, please limit the input range within 200-length for better resolution."),

      br(),

      selectInput("chr", "Choose a chromosome:", c(c(1:22),"X", "Y"), selected = 3),

      textInput("start", "Enter the start coordinate:", 49359145),

      textInput("end", "Enter the end coordinate:", 49411645),

      actionButton(inputId = "update",
                   label = "Update coordinate"),

      hr(),

      tags$strong('Sequency Lengths(bp) Reference'),

      br(),

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
          p("Chromosome Y: 57227415")
        ),
      ),

    mainPanel(
      #outputs
      tabsetPanel(type = "tabs",
                  tabPanel("Summary", tableOutput("summary")),
                  tabPanel("nsSNPplot", plotOutput("nssnp"))
      )
    )
  )
)

server <- function(input, output) {
  snp <- eventReactive(input$update,{
    nsSNPCalculatebyRange(as.numeric(input$chr),
                  as.numeric(input$start),
                  as.numeric(input$end))
  })

  p <- eventReactive(input$update,{
    nsSNPFreqPlot(as.numeric(input$chr),
                          as.numeric(input$start),
                          as.numeric(input$end))
  })

  output$summary <- renderTable ({
    if (! is.null(snp))
      snp()
  })

  output$nssnp <- renderPlot({
    if (! is.null(p))
      p()
    })
}
shiny::shinyApp(ui = ui, server = server)
# [END]
