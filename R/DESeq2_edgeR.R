edgeRRawCountAnalysis <- function(path.prefix, independent.variable, control.group, experiment.group) {
  cat(paste0("\n************** edgeR analysis **************\n"))
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR"))
  }
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images"))
  }
  # likelihood ratio test and quasi-likelihood F-test
  # Read pheno_data
  pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
  pheno_data.group <- pheno_data[independent.variable][[1]]
  # gene count table
  gene.count.table <- read.csv(paste0(path.prefix, "gene_data/reads_count_matrix/gene_count_matrix.csv"))
  # transcript count table
  # transcript.count.table <- read.csv(paste0(path.prefix, "gene_data/reads_count_matrix/transcript_count_matrix.csv"))
  # gene name list
  gene.data.frame <- data.frame(gene.id = gene.count.table$gene_id)
  # gene count table without first row
  gene.count.table <- gene.count.table[-1]
  rownames(gene.count.table) <- gene.data.frame$gene.id
  # transcript name list
  # transcript.data.frame <- data.frame(gene.id = transcript.count.table$transcript_id)
  # gene count table without first row
  # transcript.count.table <- transcript.count.table[-1]
  # rownames(transcript.count.table) <- transcript.data.frame$gene.id
  # gene matrix
  gene.count.table <- as.matrix(gene.count.table)
  # set up condition
  sample.table <- as.data.frame(table(pheno_data[independent.variable]))
  control.group.size <- sample.table[sample.table$Var1 == control.group,]$Freq
  experiment.group.size <- sample.table[sample.table$Var1 == experiment.group,]$Freq

  # create DGEList object (edgeR)
  cat("\u25CF Creating 'DGEList' object from count matrix ... \n")
  deglist.object <- edgeR::DGEList(counts=gene.count.table, group = pheno_data.group, genes = gene.data.frame)

  # plotMDS(mds, col=col, labels=group)

  # nc = edgeR::cpm(rawdata[,2:9], normalized.lib.sizes=TRUE)
  # rownames(nc) <- raw.data[,1]

  # Normalization with TMM
  cat("     \u25CF Normalizing DGEList object ... \n")
  deglist.object <- edgeR::calcNormFactors(deglist.object, method="TMM")

  cat("\u25CF Plotting edgeR MDS plot ... \n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/MDS_plot.png"))
  my_colors=c(rgb(255, 47, 35,maxColorValue = 255),
              rgb(50, 147, 255,maxColorValue = 255))
  limma::plotMDS(deglist.object, top = 1000, labels = NULL, col = my_colors[as.numeric(deglist.object$samples$group)],
                 pch = 20, cex = 2)
  par(xpd=TRUE)
  legend("bottomright",inset=c(0,1), horiz=TRUE, bty="n", legend=levels(deglist.object$samples$group) , col=my_colors, pch=20 )
  dev.off()

  # de = edgeR::exactTest(y, pair=c(control.group, experiment.group))

  # countsPerMillion <- edgeR::cpm(deglist.object)
  # logcountsPerMillion <- edgeR::cpm(deglist.object, log=TRUE)
  # # countCheck <- countsPerMillion > 1

  # estimating Dispersions
  independent.variable <- as.character(pheno_data.group)
  designMat <- model.matrix(~independent.variable)

  dgList <- edgeR::estimateGLMCommonDisp(deglist.object, design=designMat)
  dgList <- edgeR::estimateGLMTrendedDisp(dgList, design=designMat)
  dgList <- edgeR::estimateGLMTagwiseDisp(dgList, design=designMat)

  cat("\u25CF Plotting edgeR MeanVar plot ... \n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/MeanVar_plot.png"))
  edgeR::plotMeanVar(dgList, show.tagwise.vars=TRUE, NBline=TRUE)
  dev.off()

  cat("\u25CF Plotting edgeR BCV plot ...\n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/BCV_plot.png"))
  edgeR::plotBCV(dgList)
  dev.off()

  fit <- edgeR::glmFit(dgList, designMat)
  lrt <- edgeR::glmLRT(fit, coef=2)
  edgeR_result <- edgeR::topTags(lrt)
  deGenes <- edgeR::decideTestsDGE(lrt, p=0.001)
  deGenes <- rownames(lrt)[as.logical(deGenes)]
  cat("\u25CF Plotting edgeR Smear plot ... \n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/Smear_plot.png"))
  edgeR::plotSmear(lrt, de.tags=deGenes)
  abline(h=c(-1, 1), col=2)
  dev.off()
}

#' @import DESeq2
DESeq2RawCountAnalysis <- function(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC) {
  cat(paste0("\n************** DESeq2 analysis **************\n"))
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2"))
  }
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/images"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/images"))
  }
  # Read pheno_data
  pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
  # gene count table
  gene.count.table <- read.csv(paste0(path.prefix, "gene_data/reads_count_matrix/gene_count_matrix.csv"))
  # transcript count table
  # transcript.count.table <- read.csv(paste0(path.prefix, "gene_data/reads_count_matrix/transcript_count_matrix.csv"))
  # gene name list
  gene.data.frame <- data.frame(gene.id = gene.count.table$gene_id)
  # gene count table without first row
  gene.count.table <- gene.count.table[-1]
  rownames(gene.count.table) <- gene.data.frame$gene.id
  # transcript name list
  # transcript.data.frame <- data.frame(gene.id = transcript.count.table$transcript_id)
  # gene count table without first row
  # transcript.count.table <- transcript.count.table[-1]
  # rownames(transcript.count.table) <- transcript.data.frame$gene.id
  # gene matrix
  gene.count.table <- as.matrix(gene.count.table)
  # set up condition
  sample.table <- as.data.frame(table(pheno_data[independent.variable]))
  control.group.size <- sample.table[sample.table$Var1 == control.group,]$Freq
  experiment.group.size <- sample.table[sample.table$Var1 == experiment.group,]$Freq
  # condition <- c(rep(control.group, control.group.size), rep(experiment.group, experiment.group.size))

  # Deseq data
  colData <- data.frame("independent.variable" = as.character(pheno_data[independent.variable][[1]]))

  rownames(colData) <- as.character(pheno_data$ids)

  # creat DESeqDataSet
  # Rows of colData correspond to columns of countData
  cat("\u25CF Creating 'DESeqDataSet' object from count matrix ... \n")
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = gene.count.table,
                                          colData = colData,
                                          design =  ~independent.variable)
  # filter out rowSums bigger than 0 !!
  cat("     \u25CF Filtering out row sum of matrix that is equal to 0 \n")
  dds <- dds[rowSums(counts(dds))>0, ]
  # performs a default analysis
  # estimation of size factors, estimation of dispersion, Negative Binomial GLM fitting and Wald statistics
  dds <- DESeq2::DESeq(dds)

  # result function
  #Set to Inf or FALSE to disable the resetting of p-values to NA.
  # cooksCutoff : this test excludes the Cook's distance of samples belonging to experimental groups with only 2 samples.
  # independentFiltering : whether independent filtering should be applied automatically
  res <- DESeq2::results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
  cat(paste0("\n\u25CF Writing '", path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2_", control.group, "_vs_", experiment.group, "'\n"))
  write.csv(res, file = paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/DESeq2_", control.group, "_vs_", experiment.group))
  cat(paste0("\u25CF Plotting DESeq2 MA plot\n"))
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/images/DESeq2_MA_plot.png"))
  DESeq2::plotMA(dds,main="MAplot")
  dev.off()

  # reorder the result by padj !!
  res.sort.padj <- res[order(res$padj),]

  # filter out res.sort.padj (padj not null, padj < value, log2FoldChange >= 1)
  sig <-res.sort.padj[(!is.na(res.sort.padj$padj)) && (res.sort.padj$padj < DESeq2.padj) && (abs(res.sort.padj$log2FoldChange) >= DESeq2.log2FC)]

  # plot plotDispEsts
  # plotDispEsts
  cat(paste0("\u25CF Plotting DESeq2 Dispersion plot\n"))
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/images/DESeq2_Dispersion_plot.png"))
  DESeq2::plotDispEsts(dds, main="Dispersion plot")
  dev.off()
}

DESeq2edgeRRawCountAnalysis <- function(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC) {
  cat(paste0("\n************** Reads count matrix analysis **************\n"))
  if(!file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/gene_count_matrix.csv"))){
    cat("'gene_count_matrix.csv' is not exist!\n")
    stop("'gene_count_matrix.csv' not exist ERROR")
  }
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/"))
  }
  DESeq2RawCountAnalysis(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC)
  edgeRRawCountAnalysis(path.prefix, independent.variable, control.group, experiment.group)
}


