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
#'    overlap. Exclude text labels that overlap too many things. Defaults to 15.
#' @param filePath A character string path to the directory that the user want
#'    to store the visualization of the differential expression analysis result.
#'    It is recommended to create a new directory to store the output .png file.
#'    Default value is NULL. Should in the format: "/Path/to/the/directory/".
#'
#' @return A NULL, save the output visualization plot to "filePath"
#'
#' @examples
#' # Example 1:
#'
#' @references
#' R Core Team (2023). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' \href{https://www.R-project.org/}{link}.
#'
#' Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4
#' \href{ https://CRAN.R-project.org/package=ggplot2}{Link}.
#'
#' Slowikowski K (2023). _ggrepel: Automatically Position Non-Overlapping Text Labels with
#' ggplot2'_. R package version 0.9.4, \href{https://CRAN.R-project.org/package=ggrepel}{link}.
#'
#' @export
#' @import ggplot2
#' @import ggrepel

visDeAnaly <- function(genes = NULL,
                       logFCT = 2,
                       padjT = 0.05,
                       maxOverlap = 15,
                       filePath = NULL) {

  # Performing checks of user input
  if (is.null(filePath) == TRUE) {
    stop("Please input a file path to store the output files.
         Use ?visDeAnaly to check arguments.")
  } else {
    ;
  }

  if (is.null(genes) == TRUE) {
    stop("Please input the output from function 'labelGenes'.
         Use ?visDeAnaly to check arguments.")
  } else {
    ;
  }

  # Select genes that have significantly different expression for labelling
  message("Plotting and saving Differential Expression analysis result to ", filePath)

  geneToLabel <- genes[genes$padj < padjT
                            & abs(genes$log2FoldChange) >= logFCT, , drop = FALSE]

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
    ggplot2::scale_colour_manual(values = c('steelblue','gray','brown')) +
    ggplot2::geom_hline(yintercept = -log10(padjT), linetype = "dotdash") +
    ggplot2::geom_vline(xintercept = c(-logFCT, logFCT), linetype = "dotdash") +
    ggplot2::labs(title = "Volcano Plot For Differential Expression Analysis") +
    suppressWarnings(ggrepel::geom_label_repel(data = geneToLabel,
                              ggplot2::aes(label = row.names(geneToLabel)),
                              size = 3,
                              segment.color = "black",
                              show.legend = FALSE,
                              max.overlaps = maxOverlap,
                              na.rm = TRUE))

  # Save the volcano plot
  ggplot2::ggsave(file.path(filePath, "DiffExpresion_volcano_plot.png"),
                  plot = volcanoPlot,
                  dpi = 300)
  return(invisible(NULL))
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
#'    marked with relatively "steelblue" color, positive value are marked with
#'    relatively "brown" color. 0 is marked as "white".
#' @param filePath A character string path to the directory that the user want
#'    to store the visualization of the correlation analysis result.
#'    It is recommended to create a new directory to store the output .png file.
#'    Default value is NULL. Should in the format: "/Path/to/the/directory/".
#'
#' @return A NULL, save the output visualization plot to "filePath"
#'
#' @examples
#' # Example 1:
#'
#' @references
#' R Core Team (2023). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' \href{https://www.R-project.org/}{link}.
#'
#' Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4
#' \href{ https://CRAN.R-project.org/package=ggplot2}{Link}.
#'
#' Slowikowski K. (2023). ggrepel: Automatically Position Non-Overlapping Text Labels with
#' ggplot2'. R package version 0.9.4, \href{https://CRAN.R-project.org/package=ggrepel}{link}.
#'
#' Kolde R. (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12,
#' \href{https://CRAN.R-project.org/package=pheatmap}{link}.
#'
#' @export
#' @import ggplot2
#' @import ggrepel
#' @import grDevices
#' @import pheatmap

visCorrelationAnaly <- function(corMatrix = NULL,
                                colors = c("steelblue","white","brown"),
                                filePath = NULL) {

  # Performing checks of user input
  if (is.null(filePath) == TRUE) {
    stop("Please input a file path to store the output files.
         Use ?visCorrelationAnaly to check arguments.")
  } else {
    ;
  }

  if (is.null(corMatrix) == TRUE) {
    stop("Please input the correlation matrix for visualization. Output of
         function corrAnalysis could be the input.
         Use ?visCorrelationAnaly to check arguments.")
  } else {
    ;
  }

  message("Plotting and saving gene correlation analysis result to ", filePath)

  # Generate the saving path
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

  return(invisible(NULL))
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
#'    file. Default value is NULL. Should in the format: "/Path/to/the/directory/".
#'
#' @return A NULL, save the output visualization plot to "filePath"
#'
#' @examples
#' # Example 1:
#'
#' @references
#' R Core Team (2023). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' \href{https://www.R-project.org/}{link}.
#'
#' Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4
#' \href{ https://CRAN.R-project.org/package=ggplot2}{Link}.
#'
#' Vaidyanathan R, Xie Y, Allaire J, Cheng J, Sievert C, Russell K (2023). htmlwidgets: HTML
#' Widgets for R. R package version 1.6.2, \href{https://CRAN.R-project.org/package=htmlwidgets}{link}.
#'
#' Kolberg L, Raudvere U, Kuzmin I, Vilo J, Peterson H (2020). “gprofiler2-an R package for
#' gene list functional enrichment analysis and namespace conversion toolset g:Profiler.”
#' _F1000Research_, *9 (ELIXIR)*(709). R package version 0.2.2.
#'
#' @export
#' @import ggplot2
#' @import htmlwidgets
#' @import gprofiler2

visEnrichAnaly <- function(enrichOutputList = NULL,
                           interactive = TRUE,
                           filePath = NULL) {

  # Performing checks of user input
  if (is.null(filePath) == TRUE) {
    stop("Please input a file path to store the output files.
         Use ?visEnrichAnaly to check arguments.")
  } else {
    ;
  }

  if (is.null(enrichOutputList) == TRUE) {
    stop("Please input the correlation matrix for visualization. Output of
         function corrAnalysis could be the input.
         Use ?visEnrichAnaly to check arguments.")
  } else {
    ;
  }

  message("Plotting and saving enrichment analysis result to ", filePath)

  # Get the gprofiler2 result and plot the enrichment analysis result
  gProfilerResult <- enrichOutputList$gProfilerResult
  plot <- gprofiler2::gostplot(gProfilerResult, interactive = interactive)

  # Decide if the plot is interactive or not
  if (interactive == TRUE) {
    htmlwidgets::saveWidget(plot,
                            file = file.path(filePath,
                                             "enrich_analysis_vis.html"))
  } else {
    ggplot2::ggsave(file.path(filePath, "enrich_analysis_vis.png"),
                    plot = plot,
                    dpi = 300)
  }

  return(invisible(NULL))
}



#' Generate a Lollipop plot for the enrichment analysis result
#'
#' This is a function that can visualize the enrichment analysis
#' result as a Lollipop plot, whose x-axis is the gene ratio and
#' y-axis is the name of the pathways. Basically, the larger the gene ratio is,
#' the more genes enrich in the pathway.
#'
#' @param enrichOutputList A named list containing the named list of enrichment
#'    analysis results from gprofiler2 and a processed data frame based on the
#'    results from gprofiler2 for visualization. The output of function
#'    "enrichAnalysis" is a good example. Default value is NULL.
#' @param filePath A character string path to the directory that the user want
#'    to store the visualization of the enrichment analysis result.
#'    It is recommended to create a new directory to store the output .png file.
#'    Default value is NULL. Should in the format: "/Path/to/the/directory/".
#'
#' @return A NULL, save the output visualization plot to "filePath"
#'
#' @examples
#' # Example 1:
#'
#' @references
#' R Core Team (2023). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' \href{https://www.R-project.org/}{link}.
#'
#' Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4
#' \href{ https://CRAN.R-project.org/package=ggplot2}{Link}.
#'
#' @export
#' @import ggplot2

visEnrichAnalyLollipop <- function(enrichOutputList = NULL,
                                   filePath = NULL) {

  # Performing checks of user input
  if (is.null(filePath) == TRUE) {
    stop("Please input a file path to store the output files.
         Use ?visEnrichAnalyLollipop to check arguments.")
  } else {
    ;
  }

  if (is.null(enrichOutputList) == TRUE) {
    stop("Please input the correlation matrix for visualization. Output of
         function corrAnalysis could be the input.
         Use ?visEnrichAnalyLollipop to check arguments.")
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

  ggplot2::ggsave(file.path(filePath, "enrich_analysis_vis_Lollipop.png"),
                  plot = lollipopPlot,
                  dpi = 300)

  return(invisible(NULL))
}

# [END]
