#' Generate a volcano plot for the differential expression analysis result
#'
#' This is a function that can visualize the differential expression analysis
#' result as a volcano plot, whose x-axis is "Log2(Fold change)" and y-axis is
#' "-Log10(P.adj)". Basically, the larger the value of -Log10(P.adj) and the
#' absolute value of Log2(Fold change) are, the more significant the gene
#' expression differences between samples are.
#'
#' @param genes A data frame, which has a column called "group",
#'    indicating if the gene is up expressed or down expressed or no change for
#'    its expression. The output of function "labelGenes" is a good example.
#' @param logFCT A positive integer, the log fold change threshold used to
#'    select genes with significant large differences in expression levels.
#'    Default value is 2.
#' @param padjT  A positive double between 0 and 1, used to select genes with
#'    significant differences in expression levels. Default value is 0.05.
#' @param maxOverlap The max number of labels of genes in the figure that can
#'    overlap. Exclude text labels that overlap too many things. Defaults to 8.
#' @param colorBlindF A boolean indicating whther to use colorblind-friendly
#'    colors. Default is FALSE.
#' @param filePath A character string path to the directory that the user want
#'    to store the visualization of the differential expression analysis result.
#'    It is recommended to create a new directory to store the output .png file.
#'    Default value is NULL. Should in the format: "/Path/to/the/directory".
#' @param save A boolean value to indicate if store the result into a file
#'    at the file path. If set to TRUE, the output would be stored at the
#'    filePath. Default value is FALSE.
#'
#' @return A ggplot object, which is the visualization and save the output
#'    visualization plot to "filePath"
#'
#' @examples
#' # Example 1:
#' # Using labelled genes (labelledGenes) available with the package
#' \dontrun{
#' dim(labelledGenes) # 18309 rows, 7 columns
#'
#' # Plot a volcano plot for the differential expression analysis results
#' # with colorblind-friendly colors
#' volcanoplot <- visDeAnaly(labelledGenes,
#'                           colorBlindF = TRUE)
#'
#' # You should see visualizations saved in the current working directory,
#' # which is the path from getwd()
#' }
#'
#' @references
#' R Core Team (2023). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/
#'
#' Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. https://CRAN.R-project.org/package=ggplot2
#'
#' @export
#' @import ggplot2

visDeAnaly <- function(genes = NULL,
                       logFCT = 2,
                       padjT = 0.05,
                       maxOverlap = 8,
                       colorBlindF = FALSE,
                       filePath = NULL,
                       save = FALSE) {

  # Performing checks of user input
  if (save == TRUE & is.null(filePath) == TRUE) {
    message("Please input a file path to store the output files.
         Use ?visDeAnaly to check arguments.")
    return(invisible(NULL))
  } else {
    ;
  }

  if (is.null(genes) == TRUE) {
    message("Please input the output from function 'labelGenes'.
         Use ?visDeAnaly to check arguments.")
    return(invisible(NULL))
  } else {
    ;
  }

  # Select genes that have significantly different expression for labelling
  message("Plotting and saving Differential Expression analysis result to ", filePath)

  if (colorBlindF == TRUE){
    colors <- c("#56B4E9", "#999999", "#E69F00")
  } else {
    colors <- c('steelblue','gray','brown')
  }
  # Create the volcano plot for differential expression analysis
  volcanoPlot <- ggplot2::ggplot(genes,
                                 ggplot2::aes(x = log2FoldChange,
                                              y = -log10(padj),
                                              color = group)) +
    ggplot2::geom_point(alpha = 0.4, size = 2, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::xlab("Log2(Fold change)") +
    ggplot2::ylab("-Log10(P.adj)") +
    ggplot2::theme(plot.title = element_text(size = 15, hjust = 0.5)) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::geom_hline(yintercept = -log10(padjT), linetype = "dotdash") +
    ggplot2::geom_vline(xintercept = c(-logFCT, logFCT), linetype = "dotdash") +
    ggplot2::labs(title = "Volcano Plot For Differential Expression Analysis")

  # Save the volcano plot
  if (save == TRUE) {
    ggplot2::ggsave(file.path(filePath, "DiffExpresion_volcano_plot.png"),
                    plot = volcanoPlot,
                    dpi = 300)
  } else {
    ;
  }

  return(volcanoPlot)
}

#' Generate a heatmap for the gene correlation analysis result
#'
#' This is a function that can visualize the gene correlation analysis result
#' as a heatmap, whose x-axis and y-axis are gene names, each cell in the
#' heatmap indicates the correlation coefficient of corresponding genes.
#'
#' @param corMatrix A correlation coefficient matrix that has n * n dimensions
#'    (n is the number of genes). Each cell contains the correlation coefficient
#'    between corresponding genes. Output of function "enrichAnalysis" is a good
#'    example.
#' @param colors A vector that has three colors used. Default as
#'    c("steelblue","white","brown"). For the default value, negative values are
#'    marked with relatively "steelblue" color, positive values are marked with
#'    relatively "brown" color. 0 is marked as "white".
#' @param filePath A character string path to the directory that the user want
#'    to store the visualization of the correlation analysis result.
#'    It is recommended to create a new directory to store the output .png file.
#'    Default value is NULL. Should in the format: "/Path/to/the/directory".
#' @param save A boolean value to indicate if store the result into a file
#'    at the file path. If set to TRUE, the output would be stored at the
#'    filePath. Default value is FALSE.
#'
#' @return A ggplot object, and save the output visualization plot to "filePath"
#'
#' @examples
#' # Example 1:
#' # Using correlation matrix (geneCorResult) available with the package
#' \dontrun{
#' dim(geneCorResult) # 103 rows, 103 columns
#'
#' # Plot a heatmap for the correlation analysis results
#' visCorrelationAnaly(geneCorResult)
#'
#' # You should see visualizations saved in the current working directory,
#' # which is the path from getwd()
#' }
#'
#' @references
#' Kolde R. (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12,
#' https://CRAN.R-project.org/package=pheatmap
#'
#' R Core Team (2023). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/
#'
#' Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4
#' https://CRAN.R-project.org/package=ggplot2
#'
#' Yu G. (2023). ggplotify: Convert Plot to 'grob' or 'ggplot' Object.
#' R package version 0.1.2, <https://CRAN.R-project.org/package=ggplotify>.
#'
#' @export
#' @import ggplot2
#' @import grDevices
#' @import pheatmap
#' @import ggplotify

visCorrelationAnaly <- function(corMatrix = NULL,
                                colors = c("steelblue","white","brown"),
                                filePath = NULL,
                                save =FALSE) {

  # Performing checks of user input
  if (save == TRUE & is.null(filePath) == TRUE) {
    message("Please input a file path to store the output files.
         Use ?visCorrelationAnaly to check arguments.")
    return(invisible(NULL))
  } else {
    ;
  }

  if (is.null(corMatrix) == TRUE) {
    message("Please input the correlation matrix for visualization. Output of
         function corrAnalysis could be the input.
         Use ?visCorrelationAnaly to check arguments.")
    return(invisible(NULL))
  } else {
    ;
  }

  message("Plotting and saving gene correlation analysis result to ", filePath)

  # Generate the saving path
  if (save == TRUE) {
    savePath <- file.path(filePath, "correlation_analysis_vis.png")

    # Create the color palette and plot the heatmap
    colorPalette <- grDevices::colorRampPalette(colors)(100)
    heatmap <- pheatmap::pheatmap(corMatrix,
                                  cluster_rows = FALSE,
                                  fontsize_row = 5,
                                  fontsize_col = 5,
                                  color = colorPalette,
                                  border_color = NA,
                                  filename = savePath)
  } else {
    # Create the color palette and plot the heatmap
    colorPalette <- grDevices::colorRampPalette(colors)(100)
    heatmap <- pheatmap::pheatmap(corMatrix,
                                  cluster_rows = FALSE,
                                  fontsize_row = 5,
                                  fontsize_col = 5,
                                  color = colorPalette,
                                  border_color = NA)
  }

  ggheatMap <- ggplotify::as.ggplot(heatmap)
  return(ggheatMap)
}

#' Generate a Manhattan plot for the enrichment analysis result
#'
#' This is a function that can visualize the enrichment analysis
#' result as a Manhattan plot, whose x-axis contains different databases and
#' y-axis is "-Log10(P.adj)". Basically, the larger the value of -Log10(P.adj),
#' the more the pathway is enriched in genes from the query gene list.
#'
#' @param enrichOutputList A named list containing the named list of enrichment
#'    analysis results from gprofiler2 and a processed data frame based on the
#'    results from gprofiler2 for visualization. The output of function
#'    "enrichAnalysis" is a good example.
#' @param interactive A boolean indicating whether to output an interactive plot,
#'    Default as TRUE.
#' @param filePath A character string path to the directory that the user want
#'    to store the visualization of the enrichment analysis result.
#'    It is recommended to create a new directory to store the output .png/.html
#'    file. Default value is NULL. Should in the format: "/Path/to/the/directory".
#' @param save A boolean value to indicate if store the result into a file
#'    at the file path. If set to TRUE, the output would be stored at the
#'    filePath. Default value is FALSE.
#'
#' @return Either a plotly object (if interactive = TRUE) or a ggplot object
#'    (if interactive = FALSE), save the output visualization plot to "filePath"
#'
#' @examples
#' # Example 1:
#' # Using enrichment analysis output (enrichOutputList) available with the
#' # package
#' \dontrun{
#' enrichOutputList$gProfilerResult
#' enrichOutputList$enrichmentVis
#'
#' # Plot a Manhattan plot for the enrichment analysis results
#' visEnrichAnaly(enrichOutputList)
#'
#' # You should see visualizations saved in the current working directory,
#' # which is the path from getwd()
#' }
#'
#' @references
#' Kolberg L, Raudvere U, Kuzmin I, Vilo J, Peterson H (2020). “gprofiler2-an R package for
#' gene list functional enrichment analysis and namespace conversion toolset g:Profiler.”
#' F1000Research, 9 (ELIXIR)(709). R package version 0.2.2.
#'
#' R Core Team (2023). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/
#'
#' Vaidyanathan R, Xie Y, Allaire J, Cheng J, Sievert C, Russell K (2023). htmlwidgets: HTML
#' Widgets for R. R package version 1.6.2,
#' https://CRAN.R-project.org/package=htmlwidgets
#'
#' Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4
#' https://CRAN.R-project.org/package=ggplot2
#'
#'
#' @export
#' @import ggplot2
#' @import htmlwidgets
#' @import gprofiler2

visEnrichAnaly <- function(enrichOutputList = NULL,
                           interactive = TRUE,
                           filePath = NULL,
                           save = FALSE) {

  # Performing checks of user input
  if (save == TRUE & is.null(filePath) == TRUE) {
    message("Please input a file path to store the output files.
         Use ?visEnrichAnaly to check arguments.")
    return(invisible(NULL))
  } else {
    ;
  }

  if (is.null(enrichOutputList) == TRUE) {
    message("Please input the correlation matrix for visualization. Output of
         function corrAnalysis could be the input.
         Use ?visEnrichAnaly to check arguments.")
    return(invisible(NULL))
  } else {
    ;
  }

  message("Plotting and saving enrichment analysis result to ", filePath)

  # Get the gprofiler2 result and plot the enrichment analysis result
  gProfilerResult <- enrichOutputList$gProfilerResult
  plot <- gprofiler2::gostplot(gProfilerResult, interactive = interactive)

  # Decide if the plot is interactive or not
  if (save == TRUE) {
    if (interactive == TRUE) {
      htmlwidgets::saveWidget(plot,
                              file = file.path(filePath,
                                               "enrich_analysis_vis.html"))
    } else {
      ggplot2::ggsave(file.path(filePath, "enrich_analysis_vis.png"),
                      plot = plot,
                      dpi = 300)
    }
  } else {
    ;
  }
  return(plot)
}

#' Generate a Lollipop plot for the enrichment analysis result
#'
#' This is a function that can visualize the enrichment analysis
#' result as a Lollipop plot, whose x-axis is the gene ratio and
#' y-axis is the name of the pathways. Basically, the larger the gene ratio is,
#' the more genes enrich in the pathway. The plot shows the top 10 pathways with
#' top10 highest gene ratio.
#'
#' @param enrichOutputList A named list containing the named list of enrichment
#'    analysis results from gprofiler2 and a processed data frame based on the
#'    results from gprofiler2 for visualization. The output of function
#'    "enrichAnalysis" is a good example. Default value is NULL.
#' @param filePath A character string path to the directory that the user want
#'    to store the visualization of the enrichment analysis result.
#'    It is recommended to create a new directory to store the output .png file.
#'    Default value is NULL. Should in the format: "/Path/to/the/directory/".
#' @param save A boolean value to indicate if store the result into a file
#'    at the file path. If set to TRUE, the output would be stored at the
#'    filePath. Default value is FALSE.
#'
#' @return A ggplot object, and save the output visualization plot to "filePath"
#'
#' @examples
#' # Example 1:
#' # Using enrichment analysis output (enrichOutputList) available with the
#' # package
#' \dontrun{
#' enrichOutputList$gProfilerResult
#' enrichOutputList$enrichmentVis
#'
#' # Plot a Lollipop plot for the enrichment analysis results
#' visEnrichAnalyLollipop(enrichOutputList)
#'
#' # You should see visualizations saved in the current working directory,
#' # which is the path from getwd()
#' }
#'
#' @references
#' R Core Team (2023). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/
#'
#' Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4
#' https://CRAN.R-project.org/package=ggplot2
#'
#' @export
#' @import ggplot2

visEnrichAnalyLollipop <- function(enrichOutputList = NULL,
                                   filePath = NULL,
                                   save = FALSE) {

  # Performing checks of user input
  if (save == TRUE & is.null(filePath) == TRUE) {
    message("Please input a file path to store the output files.
         Use ?visEnrichAnalyLollipop to check arguments.")
    return(invisible(NULL))
  } else {
    ;
  }

  if (is.null(enrichOutputList) == TRUE) {
    message("Please input the correlation matrix for visualization. Output of
         function corrAnalysis could be the input.
         Use ?visEnrichAnalyLollipop to check arguments.")
    return(invisible(NULL))
  } else {
    ;
  }

  message("Plotting and saving enrichment analysis result to ", filePath)

  # Get the data for visualization and convert to a data frame
  enrichVis <- enrichOutputList$enrichmentVis
  enrichDf <- as.data.frame(enrichVis)

  # Order the pathways by values of geneRatio
  enrichDf2 <- enrichDf[order(enrichDf$geneRatio), ]
  rownames(enrichDf2) <- 1:nrow(enrichDf2)

  # Plot and save the top 10 pathways that genes enrich in as a lollipop plot
  lollipopPlot <- ggplot2::ggplot(enrichDf2[1:10, ],
                                  ggplot2::aes(x = geneRatio,
                                               y = c(1:10))) +
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = c(1:10))) +
    ggplot2::geom_point(ggplot2::aes(color = p.Val, size = intersectSize)) +
    ggplot2::scale_color_viridis_b() +
    ggplot2::scale_size_continuous(range = c(2, 10)) +
    ggplot2::theme_bw() +
    ggplot2::ylab(NULL) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(color="black",
                                                       size=14)) +
    ggplot2::scale_y_continuous(breaks = seq(1, 10),
                                label = enrichDf2[1:10, ]$Description)

  if (save == TRUE) {
    ggplot2::ggsave(file.path(filePath, "enrich_analysis_vis_Lollipop.png"),
                    plot = lollipopPlot,
                    dpi = 300,
                    width=14,
                    height=8,)
  } else {
    ;
  }

  return(lollipopPlot)
}

# [END]
