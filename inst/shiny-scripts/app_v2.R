# This app is adapted from
# Grolemund, G. (2015). Learn Shiny - Video Tutorials.
# URL:https://shiny.rstudio.com/tutorial/

# The app code is adapted from
# RStudio Inc. (2013). Tabsets. Shiny Gallery.
# URL:https://shiny.rstudio.com/gallery/tabsets.html

# The app code is adapted from
# RStudio Inc. (2013). Shiny theme selector. Shiny Gallery.
# URL:https://shiny.posit.co/r/gallery/application-layout/shiny-theme-selector/

# https://stackoverflow.com/questions/46559251/how-to-add-multiple-line-breaks-conveniently-in-shiny

library(shiny)

linebreaks <- function(n){
  HTML(strrep(br(), n))
}

# Define UI for GEAnaly app ----
ui <- fluidPage(
  navbarPage(

    "GEAnaly",
    tabPanel("Introduction",
             sidebarPanel(
               h4("Brief Introduction to GEAnaly"),
               tags$p("Measuring gene expression levels in different samples
                      and linking the expression of genes of interest to a
                      biological pathway helps scientists understand gene
                      functions, biological processes and the genes that play
                      important roles in regulations. It is common to perform
                      the gene expression analysis for scientific purposes."),

               tags$p("`GEAnaly` is an R package that implements a pipeline,
                      which integrates gene differential expression analysis,
                      gene enrichment analysis and gene correlation analysis
                      together for the analysis of gene expression matrix data
                      and visualization of the results of each analysis in the
                      pipeline. The whole analysis pipeline can be executed in
                      a single command efficiently."),

               tags$p("The pipeline of the package is basically: Gene expression
                      matrix data → Gene differential expression analysis →
                      Select genes that are differentially expressed →
                      Enrichment analysis for those differentially expressed
                      genes, correlation analysis would be applied on Gene
                      expression matrix data."),

               tags$p("GEAnaly aims to make the analysis of gene expression data
                      as easy as possible without adjusting file formats between
                      different analyses.This is a shiny app page with
                      user-friendly UI that assists users without programming
                      experience to perform the analyses easily and
                      efficiently."),

               tags$p("Reference for GEAnaly:"),
               tags$p("Li, G. (2023) GEAnaly: An R Package For Analysis Of Gene
                      Expression Matrix. Unpublished.
                      URL:https://github.com/ReaganGen/GEAnaly"),

               tags$p("References:"),
               tags$a("Please see the References of README",
                      href = "https://github.com/ReaganGen/GEAnaly/tree/master")
               )
             ),
    tabPanel("Differential Expression Analysis & Enrichment Analysis",
             sidebarPanel(
               tags$a("Example Gene expression matrix Data and Sample information data
                      for Differential Expression Analysis & Enrichment Analysis",
                      href = "https://drive.google.com/drive/folders/1q4c-enivo-aubvL1sXCtsy5QY5U8aZPN?usp=share_link"),
               fileInput("expressionDataDiff",
                         "Gene expression matrix input: A .csv file containing the
                         matrix of un-normalized read counts of gene expression data,
                         normally obtained from RNA-seq or other sequencing experiments.
                         The count numbers in the matrix should be non-negative integers.
                         The row names of the matrix should be the genes, and the column
                         names should be the sample IDs",
                         accept = c(".csv")),
               fileInput("sampleInDiff",
                         "Sample information input: A .csv file containing the
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
                            DT::dataTableOutput("labelResult")
                            ),
                   tabPanel("Enrichment Analysis",
                            htmlOutput('manPlot'),
                            plotOutput('loPlot'),
                            DT::dataTableOutput("sigGenes")
                            )
                 )
               )
      ),
    tabPanel("Gene Correlation Analysis",
             sidebarPanel(
               tags$a("Example Gene expression matrix Data for Correlation Analysis",
                      href = "https://drive.google.com/file/d/13qhrDyt-0yMI32HBgNOz9QuLFlTBXt9W/view?usp=share_link"),
               fileInput("corExpressionData",
                         "Gene expression matrix input: A .csv file containing
                         the matrix of un-normalized read counts of gene
                         expression data, normally obtained from RNA-seq or
                         other sequencing experiments. The count numbers in the
                         matrix should be non-negative integers. The row names
                         of the matrix should be the genes, and the column names
                         should be the sample IDs",
                         accept = c(".csv")),
               selectInput('corco',
                           'Correlation Coefficient:',
                           c("pearson", "kendall", "spearman"))
               ),
             mainPanel(
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Input Content",
                          tags$p("Show the content of the input data"),
                          DT::dataTableOutput("expressionCorData")
                          ),
                 tabPanel("Gene Correlation Analysis",
                          plotOutput('heatMap'),
                          linebreaks(26),
                          DT::dataTableOutput("corMatrix")
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

  output$manPlot <- renderPlot(
    {
      sigGenes <- extractSignificantGene(deResult(),
                                         save = FALSE,
                                         pValue = as.double(input$pValue),
                                         foldChange = as.double(input$fdchange))

      enrichOutputListE <- enrichAnalysis(sigGenes,
                                          pvalueCutoff = input$pValueEn,
                                          correctionMethod = input$correctMethod)

      visEnrichAnaly(enrichOutputListE, interactive = FALSE)

      }
  )

  output$loPlot <- renderPlot(
    {
      sigGenes <- extractSignificantGene(deResult(),
                                         save = FALSE,
                                         pValue = as.double(input$pValue),
                                         foldChange = as.double(input$fdchange))

      enrichOutputListE <- enrichAnalysis(sigGenes,
                                          pvalueCutoff = input$pValueEn,
                                          correctionMethod = input$correctMethod)

      visEnrichAnalyLollipop(enrichOutputListE)

    }
  )

  output$expressionCorData <- DT::renderDataTable(
    DT::datatable({
      as.data.frame(corExpressionData())}
    )
  )

  output$corMatrix <- DT::renderDataTable(
    DT::datatable({
      corMatrix <- corrAnalysis(corExpressionData(),
                                method = input$corco)
      as.data.frame(corMatrix)}
    )
  )

  output$heatMap <- renderPlot(
    {corMatrix <- corrAnalysis(corExpressionData(),
                                method = input$corco)
      visCorrelationAnaly(corMatrix)},
    width = 900,
    height = 900
  )

}

shinyApp(ui, server)


