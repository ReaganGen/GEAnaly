library(gprofiler2)
enrichAnalysis <- function(significantGenes,
                           pvalueCutoff = 0.05,
                           correctionMethod = "g_SCS",
                           filePath = ".") {

  geneList <- row.names(significantGenes)
  enrichOutput <- gprofiler2::gost(query = geneList,
                                   organism = "hsapiens",
                                   user_threshold = pvalueCutoff,
                                   correction_method = correctionMethod,
                                   highlight = TRUE,
                                   evcodes = TRUE)


  gem <- enrichOutput$result[,c("term_id",
                                "p_value",
                                "term_size",
                                "query_size",
                                "intersection_size",
                                "term_name")]

  colnames(gem) <- c("GO.ID",
                     "p.Val",
                     "termSize",
                     "querySize",
                     "intersectSize",
                     "Description")

  gem$FDR <- gem$p.Val
  gem$Phenotype <- "+1"
  gem$geneRatio <- gem$intersectSize / gem$querySize

  filePathCombined <- file.path(filePath, "enrich_analysis_result.tsv")
  write.table(gem, file = filePathCombined, sep="\t", quote = F, row.names=FALSE)
  enrichOutputList <- list(gProfilerResult = enrichOutput,
                           enrichmentMap = gem)
  return(enrichOutputList)
}
