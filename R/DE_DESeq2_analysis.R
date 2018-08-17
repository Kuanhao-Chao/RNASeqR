DESeq2RawCountAnalysis <- function(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC) {
  cat(paste0("\n************** DESeq2 analysis **************\n"))
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis"))){
    dir.create(paste0(path.prefix, "RNAseq_results/DESeq2_analysis"))
  }
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic"))){
    dir.create(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic"))
  }
  #############################################
  ## Creating "DESeq2_normalized_result.csv" ##
  #############################################
  pre.de.pheno.data <- RawCountPreData(path.prefix, independent.variable, control.group, experiment.group)
  # create DGEList object (edgeR)
  cat("\u25CF Creating 'DGEList' object from count matrix ... \n")
  # creatin gene name data frame
  gene.data.frame <- data.frame(gene.id = row.names(pre.de.pheno.data$gene.count.matrix))
  # create design data.frame (independent.variable)
  colData <- data.frame("independent.variable" = as.character(pre.de.pheno.data$pheno_data[independent.variable][[1]]))
  rownames(colData) <- as.character(pre.de.pheno.data$pheno_data$ids)
  # creat DESeqDataSet
  # Rows of colData correspond to columns of countData
  cat("\u25CF Creating 'DESeqDataSet' object from count matrix ... \n")
  dds.from.matrix <- DESeq2::DESeqDataSetFromMatrix(countData = pre.de.pheno.data$gene.count.matrix,
                                                    colData = colData,
                                                    design =  ~independent.variable)
  # Pre-filter out rowSums bigger than 0 !!
  cat("     \u25CF Filtering out row sum of matrix that is equal to 0 \n")
  dds.from.matrix <- dds.from.matrix[rowSums(counts(dds.from.matrix))>0, ]

  # Run DESeq() function
  # 1.estimation of size factors, 2.estimation of dispersion, 3.gene-wise dispersion estimates, 4.mean-dispersion relationship, 5.final dispersion estimates, 6.fitting model and testing
  # Negative Binomial GLM fitting and Wald statistics
  dds_de <- DESeq2::DESeq(dds.from.matrix)
  # creating statistic result
  statistic.res <- DESeq2::results(dds_de, contrast = c("independent.variable", control.group, experiment.group))
  write.csv(statistic.res, file = paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic/statistic.csv"), row.names=FALSE)
  # Normalization method of DESeq2 is Median Ratio Normalization (MRN)
  normalized.count.table <- counts(dds_de, normalized=TRUE)
  # For control group
  control.mrn.data.frame <- data.frame(normalized.count.table[,colnames(normalized.count.table) %in% as.character(pre.de.pheno.data$control.group.data.frame$ids)])
  colnames(control.mrn.data.frame) <- paste0(as.character(pre.de.pheno.data$control.group.data.frame$ids), ".", control.group)
  write.csv(control.mrn.data.frame, file = paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic/MRN_control.csv"), row.names=FALSE)
  # For experiment group
  experiment.mrn.data.frame <- data.frame(normalized.count.table[,colnames(normalized.count.table) %in% as.character(pre.de.pheno.data$experiment.group.data.frame$ids)])
  colnames(experiment.mrn.data.frame) <- paste0(as.character(pre.de.pheno.data$experiment.group.data.frame$ids), ".", experiment.group)
  write.csv(experiment.mrn.data.frame, file = paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic/MRN_experiment.csv"), row.names=FALSE)
  # create whole data.frame
  gene.id.data.frame <- data.frame("gene_id" = row.names(control.mrn.data.frame))
  total.data.frame <- cbind(gene.id.data.frame, control.mrn.data.frame, experiment.mrn.data.frame)
  total.data.frame[paste0(control.group, ".average")] <- rowMeans(control.mrn.data.frame)
  total.data.frame[paste0(experiment.group, ".average")] <- rowMeans(experiment.mrn.data.frame)
  total.data.frame[paste0(control.group, "+", experiment.group, ".average")]<- rowMeans(total.data.frame[-1])
  DESeq2.result <- cbind(total.data.frame, statistic.res)
  # Write result into file (csv)
  write.csv(DESeq2.result, file = paste0(path.prefix, "RNAseq_results/DESeq2_analysis/DESeq2_normalized_result.csv"), row.names=FALSE)

  cat(paste0("     \u25CF Selecting differential expressed genes(DESeq2), padj-value : ", DESeq2.padj, "  log2(Fold Change) : ", DESeq2.log2FC, " ..."))
  DESeq2.result.DE <- DESeq2.result[(DESeq2.result$log2FoldChange>DESeq2.log2FC) & (DESeq2.result$padj<DESeq2.padj), ]
  DEGList.length <- length(row.names(DESeq2.result.DE))
  cat("          \u25CF ", DEGList.length, " DEG have been found !!")
  write.csv(DESeq2.result.DE, file = paste0(path.prefix, "RNAseq_results/edgeR_analysis/edgeR_normalized_DE_result.csv"), row.names=FALSE)

  #########################
  ## DESeq2 visulization ##
  #########################

  # PreDE
  if(file.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/DESeq2_normalized_result.csv")) &&
     file.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic/MRN_control.csv")) &&
     file.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic/MRN_experiment.csv")) &&
     file.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic/statistic.csv"))){
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/"))
    }
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/preDE/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/preDE/"))
    }
    # Frequency
    FrequencyPlot("DESeq2_analysis", "MRN", path.prefix, independent.variable, control.group, experiment.group)
    # Bax and Violin
    BoxViolinPlot("DESeq2_analysis", "MRN", path.prefix, independent.variable, control.group, experiment.group)
    # PCA
    PCAPlot("DESeq2_analysis", "MRN", path.prefix, independent.variable, control.group, experiment.group)
    #Correlation
    CorrelationPlot("DESeq2_analysis", "MRN", path.prefix, independent.variable, control.group, experiment.group)
    # Volcano
    VolcanoPlot("DESeq2_analysis", "MRN", path.prefix, independent.variable, control.group, experiment.group, DESeq2.padj, DESeq2.log2FC)
    # MA
    # MAPlot("DESeq2_analysis", "MRN", path.prefix, independent.variable, control.group, experiment.group, ballgown.qval)
  } else {
    stop("(\u2718) file missing ERROR.\n\n")
  }





  # # result function
  # #Set to Inf or FALSE to disable the resetting of p-values to NA.
  # # cooksCutoff : this test excludes the Cook's distance of samples belonging to experimental groups with only 2 samples.
  # # independentFiltering : whether independent filtering should be applied automatically
  # res <- DESeq2::results(dds, cooksCutoff=FALSE, independentFiltering=FALSE, contrast = c("independent.variable", control.group, experiment.group))
  # cat(paste0("\n\u25CF Writing '", path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2_", control.group, "_vs_", experiment.group, "'\n"))
  # write.csv(res, file = paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/DESeq2_", control.group, "_vs_", experiment.group, ".csv"), row.names=FALSE)
  # cat(paste0("\u25CF Plotting DESeq2 MA plot\n"))
  # png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/images/DESeq2_MA_plot.png"))
  # DESeq2::plotMA(dds,main="MAplot")
  # dev.off()
  #
  # resOrdered <- res[order(res$pvalue),]
  # summary(res)
  #
  # # reorder the result by padj !!
  # res.sort.padj <- res[order(res$padj),]
  #
  # DESeq2::plotCounts(dds, gene=which.min(res$padj), intgroup="independent.variable")
  #
  #
  # # filter out res.sort.padj (padj not null, padj < value, log2FoldChange >= 1)
  # sig <-res.sort.padj[(!is.na(res.sort.padj$padj)) && (res.sort.padj$padj < DESeq2.padj) && (abs(res.sort.padj$log2FoldChange) >= DESeq2.log2FC)]
  #
  # # plot plotDispEsts
  # # plotDispEsts
  # cat(paste0("\u25CF Plotting DESeq2 Dispersion plot\n"))
  # png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/DESeq2/images/DESeq2_Dispersion_plot.png"))
  # DESeq2::plotDispEsts(dds, main="Dispersion plot")
  # dev.off()
}
