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

  # Filtering
  # Self defined low abundance condition (a CPM of 1 corresponds to a count of 6-7 in the smallest sample)
  cat("     \u25CF Filtering DGEList object (cpm > 1 and rowSums >= 2) ... \n")
  keep <- rowSums(edgeR::cpm(deglist.object)>1) >= 2
  deglist.object <- deglist.object[keep, , keep.lib.sizes=FALSE]

  # Normalization with TMM (trimmed mean of M-values )
  cat("     \u25CF Normalizing DGEList object (TMM) ... \n")
  deglist.object <- edgeR::calcNormFactors(deglist.object, method="TMM")

  #  run the cpm function on a DGEList object which contains TMM normalisation factors ==> get TMM normalized counts !!

  # lib.size <- deglist.object$samples$lib.size
  # lib.size <- lib.size*deglist.object$samples$norm.factors
  # cpm.default(deglist.object$counts,lib.size=lib.size,log=log,prior.count=prior.count)
  #

  cat("\u25CF Plotting edgeR MDS plot ... \n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/MDS_plot.png"))
  my_colors=c(rgb(255, 47, 35,maxColorValue = 255),
              rgb(50, 147, 255,maxColorValue = 255))
  limma::plotMDS(deglist.object, top = 1000, labels = NULL, col = my_colors[as.numeric(deglist.object$samples$group)],
                 pch = 20, cex = 2)
  par(xpd=TRUE)
  legend("bottomright",inset=c(0,1), horiz=TRUE, bty="n", legend=levels(deglist.object$samples$group) , col=my_colors, pch=20 )
  dev.off()

  # estimating Dispersions
  # qCML method. Given a DGEList object y, we estimate the dispersions using the following commands.
  dgList <- estimateCommonDisp(deglist.object)
  dgList <- estimateTagwiseDisp(dgList)


  cat("\u25CF Plotting edgeR MeanVar plot ... \n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/MeanVar_plot.png"))
  edgeR::plotMeanVar(dgList, show.tagwise.vars=TRUE, NBline=TRUE)
  dev.off()

  cat("\u25CF Plotting edgeR BCV plot ...\n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/BCV_plot.png"))
  edgeR::plotBCV(dgList)
  dev.off()

  # Testing for DE genes
  seperate.group.data.frame <- FindControlExperiment(path.prefix, independent.variable, control.group, experiment.group)
  control.group.data.frame <- seperate.group.data.frame$control.group
  experiment.group.data.frame <- seperate.group.data.frame$experiment.group

  de <- edgeR::exactTest(dgList)
  nc = edgeR::cpm(deglist.object, normalized.lib.sizes=TRUE)
  control.cpm.data.frame <- nc[,colnames(nc) %in% as.character(control.group.data.frame$ids)]
  colnames(control.cpm.data.frame) <- paste0(as.character(control.group.data.frame$ids), ".", control.group)
  experiment.cpm.data.frame <- nc[,colnames(nc) %in% as.character(experiment.group.data.frame$ids)]
  colnames(experiment.cpm.data.frame) <- paste0(as.character(experiment.group.data.frame$ids), ".", experiment.group)
  gene.id.data.frame <- data.frame(rownames(de$table))
  colnames(gene.id.data.frame) <- "gene.id"
  edgeR_result <- cbind(gene.id.data.frame, control.cpm.data.frame, experiment.cpm.data.frame, de$table)
  write.csv(edgeR_result, file = paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/edgeR_", control.group, "_vs_", experiment.group, ".csv"))



  edgeR::plotSmear(de, de.tags = de$genes)



  fit <- edgeR::glmFit(dgList)
  lrt <- edgeR::glmLRT(fit, coef=2)
  edgeR_result <- edgeR::topTags(lrt)
  deGenes <- edgeR::decideTestsDGE(lrt, p=0.01)
  deGenes <- rownames(lrt)[as.logical(deGenes)]
  cat("\u25CF Plotting edgeR Smear plot ... \n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/Smear_plot.png"))
  edgeR::plotSmear(lrt, de.tags=deGenes)
  abline(h=c(-1, 1), col=2)
  dev.off()
  cat("\n")
}

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
  control.experiment.data.frame <- FindControlExperiment(path.prefix, independent.variable, control.group, experiment.group)
  control.group.size <- length(row.names(control.experiment.data.frame$control.group))
  experiment.group.size <- length(row.names(control.experiment.data.frame$experiment.group))
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
  # Pre-filter out rowSums bigger than 0 !!
  cat("     \u25CF Filtering out row sum of matrix that is equal to 0 \n")
  dds <- dds[rowSums(counts(dds))>0, ]

  # performs a default analysis
  # estimation of size factors, estimation of dispersion, Negative Binomial GLM fitting and Wald statistics
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("independent.variable", control.group, experiment.group))

  DESeq2::resultsNames(dds)

  dds <- estimateSizeFactors(dds)
  counts(dds, normalized=TRUE)

  # result function
  #Set to Inf or FALSE to disable the resetting of p-values to NA.
  # cooksCutoff : this test excludes the Cook's distance of samples belonging to experimental groups with only 2 samples.
  # independentFiltering : whether independent filtering should be applied automatically
  res <- DESeq2::results(dds, cooksCutoff=FALSE, independentFiltering=FALSE, contrast = c("independent.variable", control.group, experiment.group))
  cat(paste0("\n\u25CF Writing '", path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2_", control.group, "_vs_", experiment.group, "'\n"))
  write.csv(res, file = paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/DESeq2_", control.group, "_vs_", experiment.group, ".csv"))
  cat(paste0("\u25CF Plotting DESeq2 MA plot\n"))
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/images/DESeq2_MA_plot.png"))
  DESeq2::plotMA(dds,main="MAplot")
  dev.off()

  resOrdered <- res[order(res$pvalue),]
  summary(res)

  # reorder the result by padj !!
  res.sort.padj <- res[order(res$padj),]

  DESeq2::plotCounts(dds, gene=which.min(res$padj), intgroup="independent.variable")


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


