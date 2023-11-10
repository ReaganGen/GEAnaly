# DESeq2 必须要求input的sample information和count matrix里面的gene顺序一样，别忘记把这个事写代码check
# DESeq2 要求sample information读取进来必须是factor， coldata$condition <- factor(coldata$condition)需要有这一步而且要写代码check
# 检查coldata行数和countdata列数是否相同， 检查输入的数据里是不是factor,检查输入的数据里是否都要是integer
# 检查有没有行的expression数和是0的，有的话去掉并且告诉user

library("DESeq2")
DiffExpressionAnalysis <- function(geneCounts,
                                   sampleInfor) {

  # Check if gene count table has the same sample order as the sample
  # information table and check if the samples are the same
  if (identical(colnames(geneCounts), sampleInfor$sample) != TRUE) {
    stop("The content and order of the sample names in Gene Counts table
         and Sample Infotmation table should be the same.")
  } else {
    ;
  }

  # Check if all counts in the gene count table are integers
  if (all((geneCounts - floor(geneCounts)) == 0) == FALSE) {
    stop("The count in your Gene Count table should be integers.")
  } else {
    ;
  }

  sampleInfor$sample <- factor(sampleInfor$sample)
  sampleInfor$condition <- factor(sampleInfor$condition)

  # remove genes that has no expression in all samples
  geneCounts <- geneCounts[rowSums(geneCounts) != 0, ]

  # perform analysis
  deseqDataset <- DESeq2::DESeqDataSetFromMatrix(countData = geneCounts,
                                         colData = sampleInfor,
                                         design = ~condition)

  deseqDataset <- DESeq2::DESeq(deseqDataset)
  result <- DESeq2::results(deseqDataset)
  return(result)
}

