# This app is adapted from
# Grolemund, G. (2015). Learn Shiny - Video Tutorials.
# URL:https://shiny.rstudio.com/tutorial/

# The app code is adapted from
# RStudio Inc. (2013). Tabsets. Shiny Gallery.
# URL:https://shiny.rstudio.com/gallery/tabsets.html

# The app code is adapted from
# RStudio Inc. (2013). Shiny theme selector. Shiny Gallery.
# URL:https://shiny.posit.co/r/gallery/application-layout/shiny-theme-selector/

library(shiny)

# Define UI for random distribution app ----
ui <- fluidPage(
  navbarPage(

    "GEAnaly",
    tabPanel("Introduction",
             sidebarPanel(
               tags$p("This is a the Shiny App that is part of the
                      TestingPackage in R. Its purpose is to illustrate
                      the functionality of a simple Shiny App.")
               )
             ),
    tabPanel("Differential Expression Analysis & Enrichment Analysis",
             sidebarPanel(
               tags$a("Example Gene expression matrix Data and Sample information data
                      for Differential Expression Analysis & Enrichment Analysis",
                      href = "https://drive.google.com/drive/folders/1q4c-enivo-aubvL1sXCtsy5QY5U8aZPN?usp=share_link"),
               fileInput("expressionDataDiff",
                         "Gene expression matrix input: A .rda file containing the
                         matrix of un-normalized read counts of gene expression data,
                         normally obtained from RNA-seq or other sequencing experiments.
                         The count numbers in the matrix should be non-negative integers.
                         The row names of the matrix should be the genes, and the column
                         names should be the sample IDs",
                         accept = c(".csv")),
               fileInput("sampleInDiff",
                         "Sample information input: A .rda file containing the
                         comparison table between sample IDs and sample group,
                         e.g.control group or treated group. The first column
                         is the IDs or names of the samples, the second column
                         is the group that the sample belongs to. First column
                         name shoud be 'sample', second column name should be
                         'condition'",
                         accept = c(".csv")),

               numericInput('pValue',
                            'p-Value for Differential Expression Analysis:',
                            0.05, min = 0, max = 1, step = 0.05),
               numericInput('fdchange',
                            'Fold Change for Differential Expression Analysis:',
                            2),
               numericInput('pValueEn',
                            'p-Value for Enrichment Analysis:',
                            0.05, min = 0, max = 1, step = 0.05),
               selectInput('correctMethod',
                           'Correction Method for Enrichment Analysis:',
                           c("g_SCS", "fdr", "bonferroni"))
             ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Input Content",
                            tags$p("Show the content of the input data"),
                            DT::dataTableOutput("expressionDiffData"),
                            DT::dataTableOutput("sampleDiff")
                            ),
                   tabPanel("Differential Expression Analysis",
                            plotOutput('volcanoPlot'),
                            DT::dataTableOutput("labelResult"),
                            DT::dataTableOutput("sigGenes")
                            ),
                   tabPanel("Enrichment Analysis"
                            )
                 )
               )
      ),
    tabPanel("Gene Correlation Analysis",
             sidebarPanel(
               tags$a("Example Gene expression matrix Data for Correlation Analysis",
                      href = "https://drive.google.com/file/d/13qhrDyt-0yMI32HBgNOz9QuLFlTBXt9W/view?usp=share_link"),
               fileInput("corExpressionData",
                         "Gene expression matrix input: A .rda file containing
                         the matrix of un-normalized read counts of gene
                         expression data, normally obtained from RNA-seq or
                         other sequencing experiments. The count numbers in the
                         matrix should be non-negative integers. The row names
                         of the matrix should be the genes, and the column names
                         should be the sample IDs",
                         accept = c(".rda")),
               selectInput('corco',
                           'Correlation Coefficient:',
                           c("pearson", "kendall", "spearman"))
               ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Input Content"
                 ),
                 tabPanel("Gene Correlation Analysis"
                 )
               )
             )
    )
    )
)


server <- function(input, output) {

  geneExpressionDataDiff <- reactive({

    if (! is.null(input$expressionDataDiff)) {
      read.csv(input$expressionDataDiff$datapath, row.names = 1)
    }
  })

  sampleInDiff <- reactive({
    if (! is.null(input$sampleInDiff)) {
      read.csv(input$sampleInDiff$datapath)
    }
  })

  corExpressionData <- reactive({
    if (! is.null(input$corExpressionData)) {
      read.csv(input$corExpressionData$datapath, row.names = 1)
    }
  })

  output$expressionDiffData <- DT::renderDataTable(
    DT::datatable({as.data.frame(geneExpressionDataDiff())}))
  output$sampleDiff <- DT::renderDataTable(
    DT::datatable({as.data.frame(sampleInDiff())
      }
    )
  )

  deResult <- reactive(
    {GEAnaly::diffExpressionAnalysis(geneExpressionDataDiff(),
                                                 sampleInDiff())}
  )



  output$labelResult <- DT::renderDataTable(
    DT::datatable({
      labelGenes <- labelGenes(deResult(),
                               save = FALSE,
                               pValue = as.double(input$pValue),
                               foldChange = as.double(input$fdchange))
      as.data.frame(labelGenes)}
      )
    )

  output$sigGenes <- DT::renderDataTable(
    DT::datatable({
      sigGenes <- extractSignificantGene(deResult(),
                                         save = FALSE,
                                         pValue = as.double(input$pValue),
                                         foldChange = as.double(input$fdchange))
      as.data.frame(sigGenes)}
    )
  )

  output$volcanoPlot <- renderPlot(
    {
      labelGenes <- labelGenes(deResult(),
                              save = FALSE,
                              pValue = as.double(input$pValue),
                              foldChange = as.double(input$fdchange))
      visDeAnaly(labelGenes,
                 logFCT = as.double(input$fdchange),
                 padjT = as.double(input$pValue)
                 )}
  )



  # output$plot <- renderPlot({
  #
  #
  #   # Filter out the genes that have significantly different expression levels
  #   sigGenes <- extractSignificantGene(deResult, filePath = filePath)
  #
  #   # Label genes with "UP", "DOWN" and "NOCHANGE"
  #   labelGenes <- labelGenes(deResult, filePath = filePath)
  #
  #   # Visualize the differential expression analysis
  #   visDeAnaly(labelGenes, filePath = filePath)
  #
  #   # Perform enrichment analysis
  #   enrichOutputListE <- enrichAnalysis(sigGenes, filePath = filePath)
  #
  #   # Visualize enrichment analysis result as Manhattan plot and Lollipop plot
  #   visEnrichAnaly(enrichOutputListE, filePath = filePath)
  #   visEnrichAnalyLollipop(enrichOutputListE, filePath = filePath)
  # })
}

shinyApp(ui, server)


