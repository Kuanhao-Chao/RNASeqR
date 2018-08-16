edgeRRawCountAnalysis <- function(path.prefix, independent.variable, control.group, experiment.group) {
  cat(paste0("\n************** edgeR analysis **************\n"))
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
  # create DGEList object (edgeR)
  cat("\u25CF Creating 'DGEList' object from count matrix ... \n")
  deglist.object <- edgeR::DGEList(counts=pre.de.pheno.data$gene.count.matrix, group = pre.de.pheno.data$pheno_data[independent.variable][[1]], genes = row.names(pre.de.pheno.data$gene.count.matrix))
  # Filtering
  # Self defined low abundance condition (a CPM of 1 corresponds to a count of 6-7 in the smallest sample)
  cat("     \u25CF Filtering DGEList object (raw counts row sum bigger than 0) ... \n")
  keep <- rowSums(deglist.object$counts) > 0
  deglist.object <- deglist.object[keep, , keep.lib.sizes=FALSE]
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
  gene.id.data.frame <- data.frame("gene_id" = row.names(control.cpm.data.frame))
  total.data.frame <- cbind(gene.id.data.frame, control.cpm.data.frame, experiment.cpm.data.frame)
  total.data.frame[paste0(control.group, ".average")] <- rowMeans(control.cpm.data.frame)
  total.data.frame[paste0(experiment.group, ".average")] <- rowMeans(experiment.cpm.data.frame)
  total.data.frame[paste0(control.group, "+", experiment.group, ".average")]<- rowMeans(total.data.frame[-1])
  edgeR_result <- cbind(total.data.frame, de.statistic.result$table)
  # Write result into file (csv)
  write.csv(edgeR_result, file = paste0(path.prefix, "RNAseq_results/edgeR_analysis/edgeR_normalized_result.csv"), row.names=FALSE)

  ########################
  ## edgeR visulization ##
  ########################
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images"))){
    dir.create(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images"))
  }

  # Inside package
  cat("\u25CF Plotting edgeR MDS plot ... \n")
  png(paste0(path.prefix, "RNAseq_results/edgeR_analysis/images/MDS_plot.png"))
  my_colors=c(rgb(255, 47, 35,maxColorValue = 255),
              rgb(50, 147, 255,maxColorValue = 255))
  limma::plotMDS(deglist.object, top = 1000, labels = NULL, col = my_colors[as.numeric(deglist.object$samples$group)],
                 pch = 20, cex = 2)
  par(xpd=TRUE)
  legend("bottomright",inset=c(0,1), horiz=TRUE, bty="n", legend=levels(deglist.object$samples$group) , col=my_colors, pch=20 )
  title("MDS Plot")
  dev.off()

  cat("\u25CF Plotting edgeR MeanVar plot ... \n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/MeanVar_plot.png"))
  edgeR::plotMeanVar(dgList, show.tagwise.vars=TRUE, NBline=TRUE)
  title("Mean-Variance Plot")
  dev.off()

  cat("\u25CF Plotting edgeR BCV plot ...\n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/BCV_plot.png"))
  edgeR::plotBCV(dgList)
  title("BCV Plot")
  dev.off()



  edgeR::plotSmear(de, de.tags = de$genes)

  # Fit a negative binomial generalized log-linear model
  fit <- edgeR::glmFit(dgList)
  # conducts likelihood ratio tests for one or more coefficients in the linear model
  lrt <- edgeR::glmLRT(fit, coef=2)
  topTags(lrt)
  # MA-plot is a plot of log-intensity ratios (M-values) versus log-intensity averages (A-values)
  # mean-difference plot (aka MA plot)
  cat("\u25CF Plotting MD plot ...\n")
  png(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/edgeR/images/MD_plot.png"))
  plotMD(lrt, main = "MD (MA) Plot")
  abline(h=c(-1, 1), col="blue")
  dev.off()

  # self
}
