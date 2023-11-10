library(ggplot2)
library(ggrepel)
library(pheatmap)
library(grDevices)
#要记得这里的padj和logFC要和前面filterDiffGene保持一致！！！

visDeAnaly <- function(labelGenes,
                       logFCT = 2,
                       padjT = 0.05,
                       maxOverlap = 15,
                       filePath = "./DE_volvano_plot.png") {

  geneToLabel <- subset(labelGenes, padj < padjT & abs(labelGenes$log2FoldChange) >= logFCT)
  volcanoPlot <- ggplot2::ggplot(labelGenes,
                                 ggplot2::aes(x = log2FoldChange, y = -log10(padj),
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
  ggplot2::ggsave(filePath,
         plot = volcanoPlot,
         dpi = 300)
}

visCorrelationAnaly <- function(corMatrix,
                                filePath = "./Correlation_heatmap.png",
                                colors = c("steelblue","white","brown")) {
  colorPalette <- grDevices::colorRampPalette(colors)(100)
  heatmap <- pheatmap::pheatmap(corMatrix,
                                cluster_rows = FALSE,
                                fontsize_row = 5,
                                fontsize_col = 5,
                                color = colorPalette,
                                border_color = NA,
                                filename = filePath)

}


visEnrichAnaly <- function() {

}
