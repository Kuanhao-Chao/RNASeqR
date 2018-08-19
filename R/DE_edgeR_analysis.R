edgeRRawCountAnalysis <- function(path.prefix, independent.variable, control.group, experiment.group, edgeR.pval, edgeR.log2FC) {
  cat("\u2618\u2618 edgeR analysis ...\n")
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis"))){
    dir.create(paste0(path.prefix, "RNAseq_results/edgeR_analysis"))
  }
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis/normalized_&_statistic"))){
    dir.create(paste0(path.prefix, "RNAseq_results/edgeR_analysis/normalized_&_statistic"))
  }
  #############################################
  ## Creating "edgeR_normalized_result.csv" ##
  ############################################
  pre.de.pheno.data <- RawCountPreData(path.prefix, independent.variable, control.group, experiment.group)
  raw.count <- pre.de.pheno.data$gene.count.matrix
  raw.counts.result <- RawCountGeneNameChange(raw.count, path.prefix)
  # Convert gene id to gene name
  gene.name.list <- raw.counts.result$raw.counts.name
  gene.raw.count <- raw.counts.result$raw.counts
  # create DGEList object (edgeR)
  cat("\u25CF Creating 'DGEList' object from count matrix ... \n")
  gene.data.frame <- data.frame(gene.name = gene.name.list)
  deglist.object <- edgeR::DGEList(counts=gene.raw.count, group = pre.de.pheno.data$pheno_data[independent.variable][[1]], genes = gene.name.list)
  # Normalization with TMM (trimmed mean of M-values )
  cat("     \u25CF Normalizing DGEList object (TMM) ... \n")
  deglist.object <- edgeR::calcNormFactors(deglist.object, method="TMM")
  # estimating Dispersions
  # quantile-adjusted conditional maximum likelihood (qCML) method for experiments with single factor.
  dgList <- estimateCommonDisp(deglist.object)
  dgList <- estimateTagwiseDisp(dgList)
  # Testing for DE genes
  de.statistic.result <- edgeR::exactTest(dgList)
  write.csv(de.statistic.result$table, file = paste0(path.prefix, "RNAseq_results/edgeR_analysis/normalized_&_statistic/statistic.csv"), row.names=FALSE)
  #  run the cpm function on a DGEList object which is normalisation by TMM factors ==> get TMM normalized counts !!
  # counts were first normalized with TMM and then be presented as cpm !
  normalized.count.table <- edgeR::cpm(dgList, normalized.lib.sizes=TRUE)
  # For control group
  control.cpm.data.frame <- data.frame(normalized.count.table[,colnames(normalized.count.table) %in% as.character(pre.de.pheno.data$control.group.data.frame$ids)])
  colnames(control.cpm.data.frame) <- paste0(as.character(pre.de.pheno.data$control.group.data.frame$ids), ".", control.group)
  write.csv(control.cpm.data.frame, file = paste0(path.prefix, "RNAseq_results/edgeR_analysis/normalized_&_statistic/TMM&CPM_control.csv"), row.names=FALSE)
  # For experiment group
  experiment.cpm.data.frame <- data.frame(normalized.count.table[,colnames(normalized.count.table) %in% as.character(pre.de.pheno.data$experiment.group.data.frame$ids)])
  colnames(experiment.cpm.data.frame) <- paste0(as.character(pre.de.pheno.data$experiment.group.data.frame$ids), ".", experiment.group)
  write.csv(experiment.cpm.data.frame, file = paste0(path.prefix, "RNAseq_results/edgeR_analysis/normalized_&_statistic/TMM&CPM_experiment.csv"), row.names=FALSE)
  # create whole data.frame
  total.data.frame <- cbind(gene.data.frame, control.cpm.data.frame, experiment.cpm.data.frame)
  total.data.frame[paste0(control.group, ".average")] <- rowMeans(control.cpm.data.frame)
  total.data.frame[paste0(experiment.group, ".average")] <- rowMeans(experiment.cpm.data.frame)
  total.data.frame[paste0(control.group, "+", experiment.group, ".average")]<- rowMeans(total.data.frame[-1])
  edgeR.result <- cbind(total.data.frame, de.statistic.result$table)
  # Write result into file (csv)
  write.csv(edgeR.result, file = paste0(path.prefix, "RNAseq_results/edgeR_analysis/edgeR_normalized_result.csv"), row.names=FALSE)

  # Slecet DE genes (condition)
  cat(paste0("     \u25CF Selecting differential expressed genes(edgeR) ==> p-value : ", edgeR.pval, "  log2(Fold Change) : ", edgeR.log2FC, " ...\n"))
  edgeR.result.DE <- edgeR.result[(edgeR.result$logFC>edgeR.log2FC) & (edgeR.result$PValue<edgeR.pval), ]
  cat("          \u25CF Total '", length(row.names(edgeR.result.DE)), "' DEG have been found !!")
  write.csv(edgeR.result.DE, file = paste0(path.prefix, "RNAseq_results/edgeR_analysis/edgeR_normalized_DE_result.csv"), row.names=FALSE)

  ########################
  ## edgeR visulization ##
  ########################

  # PreDE
  if(file.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis/edgeR_normalized_result.csv")) &&
     file.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis/normalized_&_statistic/TMM&CPM_control.csv")) &&
     file.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis/normalized_&_statistic/TMM&CPM_experiment.csv")) &&
     file.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis/normalized_&_statistic/statistic.csv"))){
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images"))){
      dir.create(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images"))
    }

    ###############
    #### PreDE ####
    ###############
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/preDE/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/preDE/"))
    }
    # Frequency
    FrequencyPlot("edgeR_analysis", "TMM&CPM", path.prefix, independent.variable, control.group, experiment.group)
    # Bax and Violin
    BoxViolinPlot("edgeR_analysis", "TMM&CPM", path.prefix, independent.variable, control.group, experiment.group)
    # PCA
    PCAPlot("edgeR_analysis", "TMM&CPM", path.prefix, independent.variable, control.group, experiment.group)
    #Correlation
    CorrelationPlot("edgeR_analysis", "TMM&CPM", path.prefix, independent.variable, control.group, experiment.group)

    ############
    #### DE ####
    ############
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/DE/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/DE/"))
    }
    # Volcano
    VolcanoPlot("edgeR_analysis", "TMM&CPM", path.prefix, independent.variable, control.group, experiment.group, edgeR.pval, edgeR.log2FC)

    # PCA plot
    DEPCAPlot("edgeR_analysis", "TMM&CPM", path.prefix, independent.variable, control.group, experiment.group)

    # Heatmap
    DEHeatmap("edgeR_analysis", "TMM&CPM", path.prefix, independent.variable, control.group, experiment.group)

    # MDS plot
    cat("\u25CF Plotting MDS plot ... \n")
    png(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/preDE/MDS_Plot.png"))
    my_colors=c(rgb(50, 147, 255,maxColorValue = 255),
                rgb(255, 47, 35,maxColorValue = 255))
    limma::plotMDS(deglist.object, top = 1000, labels = NULL, col = my_colors[as.numeric(deglist.object$samples$group)],
                   pch = 20, cex = 2)
    par(xpd=TRUE)
    legend("bottomright",inset=c(0,1), horiz=TRUE, bty="n", legend=levels(deglist.object$samples$group) , col=my_colors, pch=20 )
    title("MDS Plot")
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/preDE/MDS_Plot.png"), "' has been created. \n\n"))

    # MeanVar plot
    cat("\u25CF Plotting MeanVar plot ... \n")
    png(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/preDE/MeanVar_Plot.png"))
    edgeR::plotMeanVar(dgList, show.tagwise.vars=TRUE, NBline=TRUE)
    title("Mean-Variance Plot")
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/preDE/MeanVar_Plot.png"), "' has been created. \n\n"))

    # BCV plot
    cat("\u25CF Plotting BCV plot ...\n")
    png(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/preDE/BCV_Plot.png"))
    edgeR::plotBCV(dgList)
    title("BCV Plot")
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/preDE/BCV_Plot.png"), "' has been created. \n\n"))
  } else {
    stop("(\u2718) file missing ERROR.\n\n")
  }

#   edgeR::plotSmear(de, de.tags = de$genes)
#
#   # Fit a negative binomial generalized log-linear model
#   fit <- edgeR::glmFit(dgList)
#   # conducts likelihood ratio tests for one or more coefficients in the linear model
#   lrt <- edgeR::glmLRT(fit, coef=2)
#   topTags(lrt)
#   # MA-plot is a plot of log-intensity ratios (M-values) versus log-intensity averages (A-values)
#   # mean-difference plot (aka MA plot)
#   cat("\u25CF Plotting MD plot ...\n")
#   png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/MD_plot.png"))
#   plotMD(lrt, main = "MD (MA) Plot")
#   abline(h=c(-1, 1), col="blue")
#   dev.off()
}