# TPM normalization
TPMNormalizationAnalysis <- function(path.prefix, genome.name, sample.pattern, independent.variable, control.group, experiment.group, TPM.pval, TPM.log2FC) {
  cat("\n\u2618\u2618 TPM analysis ...\n")
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/TPM_analysis/"))
  }
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/normalized_&_statistic"))){
    dir.create(paste0(path.prefix, "RNASeq_results/TPM_analysis/normalized_&_statistic"))
  }
  ###############################################
  ## Creating "TPM_normalized_result.csv" ##
  ##############################################
  control.FPKM <- read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/FPKM_control.csv"))
  experiment.FPKM <- read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/FPKM_experiment.csv"))
  statistic.FPKM <- read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/statistic.csv"))
  gene.name <- read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/gene_name.csv"))

  control.TPM <- (control.FPKM/colSums(control.FPKM))*(10**6)
  write.csv(control.TPM, file = paste0(path.prefix, "RNASeq_results/TPM_analysis/normalized_&_statistic/TPM_control.csv"), row.names=FALSE)

  experiment.TPM <- (experiment.FPKM/colSums(control.FPKM))*(10**6)
  write.csv(experiment.TPM, file = paste0(path.prefix, "RNASeq_results/TPM_analysis/normalized_&_statistic/TPM_experiment.csv"), row.names=FALSE)

  gene.id.data.frame <- data.frame(read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/gene_name.csv")))

  p.value <- unlist(lapply(seq_len(nrow(control.TPM)), function(x) { stats::t.test(control.TPM[x,], experiment.TPM[x,])$p.value }))
  fold.change <- unlist(lapply(seq_len(nrow(control.TPM)), function(x) { mean(unlist(experiment.TPM[x,])) / mean(unlist(control.TPM[x,])) }))
  statistic.T.test <- data.frame("pval" = p.value, "fc" = fold.change, "log2FC" = log2(fold.change))
  write.csv(statistic.T.test, file = paste0(path.prefix, "RNASeq_results/TPM_analysis/normalized_&_statistic/statistic.csv"), row.names=FALSE)

  total.data.frame <- cbind(gene.id.data.frame, control.TPM, experiment.TPM)
  total.data.frame[paste0(control.group, ".average")] <- rowMeans(control.TPM)
  total.data.frame[paste0(experiment.group, ".average")] <- rowMeans(experiment.TPM)
  total.data.frame[paste0(control.group, ".", experiment.group, ".average")]<- rowMeans(total.data.frame[-1])
  TPM_Ttest.result <- cbind(total.data.frame, statistic.T.test)
  TPM_Ttest.result <- rbind(TPM_Ttest.result[TPM_Ttest.result$gene.name != ".",], TPM_Ttest.result[TPM_Ttest.result$gene.name == ".",])
  # Filter out pval is NaN and qval is NaN
  write.csv(TPM_Ttest.result, file = paste0(path.prefix, "RNASeq_results/TPM_analysis/TPM_normalized_result.csv"), row.names=FALSE)

  cat(paste0("     \u25CF Selecting differential expressed genes() ==> p-value : ", TPM.pval, "  log2(Fold Change) : ", TPM.log2FC, " ...\n"))
  TPM_Ttest.result.DE <- TPM_Ttest.result[(TPM_Ttest.result$log2FC>TPM.log2FC) & (TPM_Ttest.result$pval<TPM.pval), ]
  cat(paste0("          \u25CF Total '", length(row.names(TPM_Ttest.result.DE)), "' DEG have been found !!!\n"))
  write.csv(TPM_Ttest.result.DE, file = paste0(path.prefix, "RNASeq_results/TPM_analysis/TPM_normalized_DE_result.csv"), row.names=FALSE)

  ###########################
  ## TPM&Ttest visulization ##
  ###########################
  if(file.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/TPM_normalized_result.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/normalized_&_statistic/TPM_control.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/normalized_&_statistic/TPM_experiment.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/normalized_&_statistic/statistic.csv"))){
    # Transcript Related
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/images/"))){
      dir.create(paste0(path.prefix, "RNASeq_results/TPM_analysis/images/"))
    }


    ###############
    #### PreDE ####
    ###############
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/images/preDE/"))){
      dir.create(paste0(path.prefix, "RNASeq_results/TPM_analysis/images/preDE/"))
    }
    # Frequency
    FrequencyPlot("TPM_analysis", "TPM", path.prefix, independent.variable, control.group, experiment.group)
    # Bax and Violin
    BoxViolinPlot("TPM_analysis", "TPM", path.prefix, independent.variable, control.group, experiment.group)
    # PCA
    PCAPlot("TPM_analysis", "TPM", path.prefix, independent.variable, control.group, experiment.group)
    #Correlation
    CorrelationPlot("TPM_analysis", "TPM", path.prefix, independent.variable, control.group, experiment.group)

    ############
    #### DE ####
    ############
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/images/DE/"))){
      dir.create(paste0(path.prefix, "RNASeq_results/TPM_analysis/images/DE/"))
    }
    # Volcano
    VolcanoPlot("TPM_analysis", "TPM", path.prefix, independent.variable, control.group, experiment.group, TPM.pval, TPM.log2FC)
    # MA
    MAPlot("TPM_analysis", "TPM", path.prefix, independent.variable, control.group, experiment.group, TPM.pval)
    # DE PCA plot
    DEPCAPlot("TPM_analysis", "TPM", path.prefix, independent.variable, control.group, experiment.group)
    # Heatmap
    DEHeatmap("TPM_analysis", "TPM", path.prefix, independent.variable, control.group, experiment.group)
  } else {
    stop("(\u2718) file missing ERROR.\n\n")
  }


}
