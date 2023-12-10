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
    # theme = "cerulean",  # <--- To use a theme, uncomment this
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
               fileInput("file",
                         "Gene expression matrix input: A csv file containing the
                         matrix of un-normalized read counts of gene expression data,
                         normally obtained from RNA-seq or other sequencing experiments.
                         The count numbers in the matrix should be non-negative integers.
                         The row names of the matrix should be the genes, and the column
                         names should be the sample IDs",
                         accept = c(".csv")),
               fileInput("file",
                         "Sample information input: A csv file containing the
                         comparison table between sample IDs and sample group,
                         e.g.control group or treated group. The first column
                         is the IDs or names of the samples, the second column
                         is the group that the sample belongs to. First column
                         name shoud be 'sample', second column name should be
                         'condition'",
                         accept = c(".csv")),

               numericInput('pValue',
                            'p-Value for Differential Expression Analysis',
                            0.05, min = 0, max = 1),
               numericInput('fdchange',
                            'Fold Change for Differential Expression Analysis',
                            2, min = 0, max = 1),
               numericInput('pValueEn',
                            'p-Value for Enrichment Analysis',
                            0.05, min = 1, max = 9),
               selectInput('correctMethod',
                           'Correction Method:',
                           c("g_SCS", "fdr", "bonferroni"))
             ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Input Content",
                            DT::dataTableOutput("diffdata")
                            ),
                   tabPanel("Differential Expression Analysis",
                            plotOutput("plot"),
                            DT::dataTableOutput("data")
                            ),
                   tabPanel("Enrichment Analysis",
                            plotOutput("plot"),
                            DT::dataTableOutput("data")
                            )
                 )
               )
      ),
    tabPanel("Gene Correlation Analysis",
             sidebarPanel(
               tags$a("Example Gene expression matrix Data for Correlation Analysis",
                      href = "https://drive.google.com/file/d/13qhrDyt-0yMI32HBgNOz9QuLFlTBXt9W/view?usp=share_link"),
               fileInput("file",
                         "Gene expression matrix input: A csv file containing
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
                 tabPanel("Input Content",
                          DT::dataTableOutput("cordata")
                 ),
                 tabPanel("Differential Expression Analysis",
                          plotOutput("plot"),
                          DT::dataTableOutput("data")
                 ),
                 tabPanel("Enrichment Analysis",
                          plotOutput("plot"),
                          DT::dataTableOutput("data")
                 )
               )
             )
    )
    )
)


server <- function(input, output) {
  x <- rnorm(100)
  y <- rnorm(100)

  output$plot <- renderPlot({
    if (input$plotType == "scatter") {
      plot(x, y)
    } else {
      breaks <- input$breaks
      if (breaks == "custom") {
        breaks <- input$breakCount
      }

      hist(x, breaks = breaks)
    }
  })
}

shinyApp(ui, server)
