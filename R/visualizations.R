library(ggplot2)
library(ggrepel)
library(pheatmap)
library(grDevices)
library(gprofiler2)
library(htmlwidgets)
library(forcats)
#要记得这里的padj和logFC要和前面filterDiffGene保持一致！！！

visDeAnaly <- function(labelGenes,
                       logFCT = 2,
                       padjT = 0.05,
                       maxOverlap = 15,
                       filePath = ".") {

  geneToLabel <- labelGenes[labelGenes$padj < padjT
                            & abs(labelGenes$log2FoldChange) >= logFCT, , drop = FALSE]
  volcanoPlot <- ggplot2::ggplot(labelGenes,
                                 ggplot2::aes(x = log2FoldChange,
                                              y = -log10(padj),
                                              color = group)) +
    ggplot2::geom_point(alpha = 0.4, size = 2, na.rm = TRUE) +  # 设置点的透明度和大小
    ggplot2::theme_bw(base_size = 12) +  #设置一个主题背景
    ggplot2::xlab("Log2(Fold change)") + # x轴名字
    ggplot2::ylab("-Log10(P.adj)") + # y轴名字
    ggplot2::theme(plot.title = element_text(size = 15, hjust = 0.5)) +
    ggplot2::scale_colour_manual(values = c('steelblue','gray','brown')) + # 各自的颜色
    ggplot2::geom_hline(yintercept = -log10(padjT), linetype = "dotdash") + #
    ggplot2::geom_vline(xintercept = c(-logFCT, logFCT), linetype = "dotdash")+ #定义差异倍数和线形
    ggplot2::labs(title = "Volcano Plot For Differential Expression Analysis") + #加上题目
    ggrepel::geom_label_repel(data = geneToLabel,
                              ggplot2::aes(label = row.names(geneToLabel)),
                              size = 3,
                              segment.color = "black",
                              show.legend = FALSE,
                              max.overlaps = maxOverlap,
                              na.rm = TRUE)
  ggplot2::ggsave(file.path(filePath, "DiffExpresion_volcano_plot.png"),
                  plot = volcanoPlot,
                  dpi = 300)
  return(invisible(NULL))
}

visCorrelationAnaly <- function(corMatrix,
                                filePath = ".",
                                colors = c("steelblue","white","brown")) {

  savePath <- file.path(filePath, "correlation_analysis_vis.png")
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


visEnrichAnaly <- function(enrichOutputList,
                           interactive = TRUE,
                           filePath = ".") {

  gProfilerResult <- enrichOutputList$gProfilerResult

  plot <- gprofiler2::gostplot(gProfilerResult, interactive = interactive)
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


visEnrichAnalyLollipop <- function(enrichOutputList,
                                   filePath = filePath) {

  enrichMap <- enrichOutputList$enrichmentMap
  enrichDf <- as.data.frame(enrichMap)
  enrichDf2 <- enrichDf[order(enrichDf$geneRatio), ]
  rownames(enrichDf2) <- 1:nrow(enrichDf2)
  lollipopPlot <- ggplot2::ggplot(enrichDf2[1:10, ],
                                  ggplot2::aes(x = geneRatio,
                                               y = c(1:10))) +
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = c(1:10))) +
    ggplot2::geom_point(ggplot2::aes(color = p.Val, size = intersectSize)) +
    ggplot2::scale_color_viridis_b() +
    ggplot2::scale_size_continuous(range = c(2, 10)) +
    ggplot2::theme_bw() +
    ylab(NULL) +
    theme(axis.text.y = element_text(color="black",
                                     size=14)) +
    ggplot2::scale_y_continuous(breaks = seq(1, 10),
                                label = enrichDf2[1:10, ]$Description)

  ggplot2::ggsave(file.path(filePath, "enrich_analysis_vis_Lollipop.png"),
                  plot = lollipopPlot,
                  dpi = 300)
}





