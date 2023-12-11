# This app is adapted from
# Grolemund, G. (2015). Learn Shiny - Video Tutorials.
# URL:https://shiny.rstudio.com/tutorial/

# The app code is adapted from
# RStudio Inc. (2023). Tabsets. Shiny Gallery.
# URL:https://shiny.rstudio.com/gallery/tabsets.html

# The app code is adapted from
# RStudio Inc. (2023). Shiny theme selector. Shiny Gallery.
# URL:https://shiny.posit.co/r/gallery/application-layout/shiny-theme-selector/

# The code for linebreaks is adapted from
# Humpelstielzchen, May 28, 2020 (19:10), comment on rankthefirst,
# "how to add multiple line breaks conveniently in shiny?, Stack Overflow,
# Oct 04, 2017 (7:52)
# URL: https://stackoverflow.com/questions/46559251/how-to-add-multiple-line-breaks-conveniently-in-shiny

library(shiny)
library(GEAnaly)

lineBreaks <- function(n) {
  HTML(strrep(br(), n))
}

# Define UI for GEAnaly app
ui <- fluidPage(
  navbarPage(

    "GEAnaly",
    tabPanel("Introduction",
             # Content for introduction
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
                      different analyses. This is a shiny app page with
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

             # Content for the data input area for Diff expression & enrichment
             sidebarPanel(
               h4("Inputs Required:"),
               tags$p("Gene expression matrix and sample information table are
                      two inputs required for differential expression analysis
                      and enrichment analysis, details about the inputs can be
                      found below. In the example data link below,
                      geneCountsDiffExpression.csv is the gene expression
                      matrix, sampleInforDiffExpression.csv is the sample
                      information table."),

               tags$a("Example Gene expression matrix Data and Sample
                      information data for Differential Expression Analysis &
                      Enrichment Analysis",
                      href = "https://drive.google.com/drive/folders/1q4c-enivo-aubvL1sXCtsy5QY5U8aZPN?usp=share_link"),

               fileInput("expressionDataDiff",
                         "Gene expression matrix input: a .csv file containing
                         the matrix of un-normalized read counts of gene
                         expression data, normally obtained from RNA-seq or
                         other sequencing experiments. The count numbers in the
                         matrix should be non-negative integers. The row names
                         of the matrix should be the genes, and the column
                         names should be the sample IDs",
                         accept = c(".csv")),

               fileInput("sampleInDiff",
                         "Sample information input: a .csv file containing the
                         comparison table between sample IDs and sample group,
                         e.g.control group or treated group. The first column
                         is the IDs or names of the samples, the second column
                         is the group that the sample belongs to. First column
                         name shoud be 'sample', second column name should be
                         'condition'",
                         accept = c(".csv")),

               numericInput('pValue',
                            'p-Value for Differential Expression Analysis: the
                            p-value is the threshold used to decide if the
                            difference in expression is significant. Usually,
                            if the p value for genes is smaller than the value
                            you choose, the gene would be regarded to be
                            significantly differentialy expressed. Default as
                            0.05, the value is between 0 and 1',
                            0.05, min = 0, max = 1, step = 0.05),
               numericInput('fdchange',
                            'Fold Change for Differential Expression Analysis:
                            Fold change is a measure describing how much the
                            expression levels of genes change between different
                            samples. Only genes that has larger absolute values
                            of fold change than the threshold would be regarded
                            as genes express differently. Default value is 2.',
                            2),
               numericInput('pValueEn',
                            'p-Value for Enrichment Analysis: threshold used to
                            decide if the enrichment of a pathway is significant.
                            Only pathways with smaller p-value are tagged as
                            significant. Default value is 0.05. The pvalue
                            should be a value between 0 and 1.',
                            0.05, min = 0, max = 1, step = 0.05),
               selectInput('correctMethod',
                           'Correction Method for Enrichment Analysis: the
                           algorithm used for multiple testing correction,
                           one of "g_SCS" (default), "fdr", "bonferroni".
                           If you are not sure choose which one, just use g_SCS',
                           c("g_SCS", "fdr", "bonferroni"))
             ),
               mainPanel(
                 # Subtabs for the differential expression analysis and
                 # enrichment analysis
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Input Content",
                            h3("The content of the input data:"),
                            h4("The content of Gene Expression Matrix Data
                               (will appear after uploading input files):"),
                            DT::dataTableOutput("expressionDiffData"),
                            h4("The content of Sample Information for samples
                               in the gene expression data (will appear after
                               uploading input files):"),
                            DT::dataTableOutput("sampleDiff")
                            ),
                   tabPanel("Differential Expression Analysis",
                            h3("The result for Differential Expression Analysis
                               (Note: The result requires 30s to 1 min to be
                               prepared after you open the subtab):"),
                            h4("Visualization for differentialy expressed genes:"),
                            plotOutput('volcanoPlot'),
                            h4("Data used for the visualization:"),
                            DT::dataTableOutput("labelResult")
                            ),
                   tabPanel("Enrichment Analysis",
                            h3("The result for Enrichment Analysis on
                               significantly differentialy expressed genes
                               (Note: The result requires 30s to 1 min to be
                               prepared after you open the subtab):"),
                            h4("Visualization for Enrichment Analysis:"),
                            plotOutput('manPlot'),
                            plotOutput('loPlot'),
                            h4("Significant Genes used for Enrichment Analysis:"),
                            DT::dataTableOutput("sigGenes")
                            )
                   )
                 )
             ),
    tabPanel("Gene Correlation Analysis",

             # Big tab for Correlation Analysis
             sidebarPanel(
               h4("Inputs Required:"),
               tags$a("Example Gene expression matrix Data for Correlation Analysis",
                      href = "https://drive.google.com/drive/folders/1P9_PX4zsRlh7J-0tV_S-jYm2SpshhBPX?usp=share_link"),
               tags$head(tags$style('h6 {color:red;}')),
               tags$h6("Please open the example data link in a new tab/window."),
               fileInput("corExpressionData",
                         "Gene expression matrix input: a .csv file containing
                         the matrix of un-normalized read counts of gene
                         expression data, normally obtained from RNA-seq or
                         other sequencing experiments. The count numbers in the
                         matrix should be non-negative integers. The row names
                         of the matrix should be the genes, and the column names
                         should be the sample IDs",
                         accept = c(".csv")),
               selectInput('corco',
                           'Correlation Coefficient: the correlation coefficient
                           used for the correlation analysis',
                           c("pearson", "kendall", "spearman"))
               ),
             mainPanel(
               # Result subtab and input subtab for correlation analysis
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Input Content",
                          h3("The content of the input data:"),
                          h4("The content of Gene Expression Matrix Data
                               (will appear after uploading input files):"),
                          DT::dataTableOutput("expressionCorData")
                          ),
                 tabPanel("Gene Correlation Analysis",
                          h3("The result for Correlation Analysis on
                               gene expression matrix:"),
                          h4("Visualization for Correlation Analysis Result:"),
                          plotOutput('heatMap'),
                          lineBreaks(26),
                          h4("Data used for the Heatmap:"),
                          DT::dataTableOutput("corMatrix")
                 )
               )
              )
             )
            )
          )

server <- function(input, output) {

  geneExpressionDataDiff <- reactive({
    # read in file when there is something uploaded
    if (! is.null(input$expressionDataDiff)) {
      read.csv(input$expressionDataDiff$datapath, row.names = 1)
    }
  })

  sampleInDiff <- reactive({
    # read in file when there is something uploaded
    if (! is.null(input$sampleInDiff)) {
      read.csv(input$sampleInDiff$datapath)
    }
  })

  corExpressionData <- reactive({
    # read in file when there is something uploaded
    if (! is.null(input$corExpressionData)) {
      read.csv(input$corExpressionData$datapath, row.names = 1)
    }
  })

  # show the content of uploaded files
  output$expressionDiffData <- DT::renderDataTable(
    DT::datatable({as.data.frame(geneExpressionDataDiff())}))

  output$sampleDiff <- DT::renderDataTable(
    DT::datatable({as.data.frame(sampleInDiff())
    }
    )
  )

  deResult <- reactive(
    # Perform differential expression analysis
    {GEAnaly::diffExpressionAnalysis(geneExpressionDataDiff(),
                                     sampleInDiff())}
  )


  # label genes with group labels
  output$labelResult <- DT::renderDataTable(
    DT::datatable({
      labelGenes <- GEAnaly::labelGenes(deResult(),
                               save = FALSE,
                               pValue = as.double(input$pValue),
                               foldChange = as.double(input$fdchange))
      as.data.frame(labelGenes)}
    )
  )

  # Filter to obtain the significant genes
  output$sigGenes <- DT::renderDataTable(
    DT::datatable({
      sigGenes <- GEAnaly::extractSignificantGene(deResult(),
                                         save = FALSE,
                                         pValue = as.double(input$pValue),
                                         foldChange = as.double(input$fdchange))
      as.data.frame(sigGenes)}
    )
  )

  # show the visualization of differential expression analysis
  output$volcanoPlot <- renderPlot(
    {
      labelGenes <- GEAnaly::labelGenes(deResult(),
                               save = FALSE,
                               pValue = as.double(input$pValue),
                               foldChange = as.double(input$fdchange))
      GEAnaly::visDeAnaly(labelGenes,
                 logFCT = as.double(input$fdchange),
                 padjT = as.double(input$pValue)
      )}
  )

  # show the visualization for enrichment analysis
  output$manPlot <- renderPlot(
    {
      sigGenes <- GEAnaly::extractSignificantGene(deResult(),
                                         save = FALSE,
                                         pValue = as.double(input$pValue),
                                         foldChange = as.double(input$fdchange))

      enrichOutputListE <- GEAnaly::enrichAnalysis(sigGenes,
                                          pvalueCutoff = input$pValueEn,
                                          correctionMethod = input$correctMethod)

      GEAnaly::visEnrichAnaly(enrichOutputListE, interactive = FALSE)

    }
  )

  output$loPlot <- renderPlot(
    {
      sigGenes <- GEAnaly::extractSignificantGene(deResult(),
                                         save = FALSE,
                                         pValue = as.double(input$pValue),
                                         foldChange = as.double(input$fdchange))

      enrichOutputListE <- GEAnaly::enrichAnalysis(sigGenes,
                                          pvalueCutoff = input$pValueEn,
                                          correctionMethod = input$correctMethod)

      GEAnaly::visEnrichAnalyLollipop(enrichOutputListE)

    }
  )

  # show data used for correlation analysis
  output$expressionCorData <- DT::renderDataTable(
    DT::datatable({
      as.data.frame(corExpressionData())}
    )
  )

  # show result of correlation analysis
  output$corMatrix <- DT::renderDataTable(
    DT::datatable({
      corMatrix <- GEAnaly::corrAnalysis(corExpressionData(),
                                method = input$corco)
      as.data.frame(corMatrix)}
    )
  )

  output$heatMap <- renderPlot(
    {corMatrix <- GEAnaly::corrAnalysis(corExpressionData(),
                               method = input$corco)
    GEAnaly::visCorrelationAnaly(corMatrix)},
    width = 900,
    height = 900
  )
}

# Create Shiny app
shinyApp(ui, server)

# [END]
