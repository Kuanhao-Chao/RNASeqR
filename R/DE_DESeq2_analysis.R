DESeq2RawCountAnalysis <- function(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC) {
  cat("\u2618\u2618 DESeq2 analysis ...\n")
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
  raw.count <- pre.de.pheno.data$gene.count.matrix
  raw.counts.result <- RawCountGeneNameChange(raw.count, path.prefix)
  # Convert gene id to gene name
  gene.name.list <- raw.counts.result$raw.counts.name
  gene.raw.count <- raw.counts.result$raw.counts

  cat("\u25CF Creating 'DGEList' object from count matrix ... \n")
  # creatin gene name data frame
  gene.data.frame <- data.frame(gene.name = gene.name.list)
  # create design data.frame (independent.variable)
  colData <- data.frame("independent.variable" = as.character(pre.de.pheno.data$pheno_data[independent.variable][[1]]))
  rownames(colData) <- as.character(pre.de.pheno.data$pheno_data$ids)
  # creat DESeqDataSet
  # Rows of colData correspond to columns of countData
  cat("\u25CF Creating 'DESeqDataSet' object from count matrix ... \n")
  dds.from.matrix <- DESeq2::DESeqDataSetFromMatrix(countData = gene.raw.count,
                                                    colData = colData,
                                                    design =  ~independent.variable)
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
  total.data.frame <- cbind(gene.data.frame, control.mrn.data.frame, experiment.mrn.data.frame)
  total.data.frame[paste0(control.group, ".average")] <- rowMeans(control.mrn.data.frame)
  total.data.frame[paste0(experiment.group, ".average")] <- rowMeans(experiment.mrn.data.frame)
  total.data.frame[paste0(control.group, "+", experiment.group, ".average")]<- rowMeans(total.data.frame[-1])
  DESeq2.result <- cbind(total.data.frame, statistic.res)
  DESeq2.result.no.NA <- DESeq2.result[((!is.na(DESeq2.result$pvalue)) & (!is.na(DESeq2.result$padj))), ]
  # Write result into file (csv)
  write.csv(DESeq2.result.no.NA, file = paste0(path.prefix, "RNAseq_results/DESeq2_analysis/DESeq2_normalized_result.csv"), row.names=FALSE)

  cat(paste0("     \u25CF Selecting differential expressed genes(DESeq2) ==> padj-value : ", DESeq2.padj, "  log2(Fold Change) : ", DESeq2.log2FC, " ...\n"))
  DESeq2.result.DE <- DESeq2.result.no.NA[(DESeq2.result.no.NA$log2FoldChange>DESeq2.log2FC) & (DESeq2.result.no.NA$padj<DESeq2.padj), ]
  cat("          \u25CF Total '", length(row.names(DESeq2.result.DE)), "' DEG have been found !!")
  write.csv(DESeq2.result.DE, file = paste0(path.prefix, "RNAseq_results/DESeq2_analysis/DESeq2_normalized_DE_result.csv"), row.names=FALSE)

  #########################
  ## DESeq2 visulization ##
  #########################

  if(file.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/DESeq2_normalized_result.csv")) &&
     file.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic/MRN_control.csv")) &&
     file.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic/MRN_experiment.csv")) &&
     file.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/normalized_&_statistic/statistic.csv"))){
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/"))
    }
    ###############
    #### PreDE ####
    ###############
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

    ############
    #### DE ####
    ############
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/DE/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/DE/"))
    }
    # Volcano
    VolcanoPlot("DESeq2_analysis", "MRN", path.prefix, independent.variable, control.group, experiment.group, DESeq2.padj, DESeq2.log2FC)

    # MA plot
    cat(paste0("\u25CF Plotting  MA plot\n"))
    png(paste0(path.prefix, paste0("RNAseq_results/DESeq2_analysis/images/DE/MA_Plot.png")))
    p1 <- DESeq2::plotMA(dds_de, main = "MA (MD) Plot")
    print(p1)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/DE/MA_Plot.png"), "' has been created. \n\n"))

    # PCA plot
    DEPCAPlot("DESeq2_analysis", "MRN", path.prefix, independent.variable, control.group, experiment.group)

    # Heatmap
    DEHeatmap("DESeq2_analysis", "MRN", path.prefix, independent.variable, control.group, experiment.group)

    # counts plot
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/DE/Counts_Gene"))){
      dir.create(paste0(path.prefix, "RNAseq_results/DESeq2_analysis/images/DE/Counts_Gene"))
    }
    cat(paste0("\u25CF Plotting  Counts Gene plot\n"))
    for( i in DESeq2.result.DE$gene_id ) {
      png(paste0(path.prefix, paste0("RNAseq_results/DESeq2_analysis/images/DE/Counts_Gene/Counts_", i , "_Plot.png")))
      p <- DESeq2::plotCounts(dds_de, gene = i, intgroup = "independent.variable")
      print(p)
      dev.off()
      cat(paste0("(\u2714) : '", paste0("RNAseq_results/DESeq2_analysis/images/DE/Counts_Gene/Counts_", i , "_Plot.png"), "' has been created. \n\n"))
    }

    # dispersion plot
    cat(paste0("\u25CF Plotting  Dispersion plot\n"))
    png(paste0(path.prefix, paste0("RNAseq_results/DESeq2_analysis/images/DE/Dispersion_Plot.png")))
    plotDispEsts(dds_de, main="Dispersion plot")
    dev.off()
    cat(paste0("(\u2714) : '", paste0("RNAseq_results/DESeq2_analysis/images/DE/Dispersion_Plot.png"), "' has been created. \n\n"))
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
