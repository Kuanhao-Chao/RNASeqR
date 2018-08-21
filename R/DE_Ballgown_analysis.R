# Run ballgown analysi
BallgownAnalysis <- function(path.prefix, genome.name, sample.pattern, independent.variable, control.group, experiment.group, ballgown.qval, ballgown.log2FC) {
  cat("\n\u2618\u2618 ballgown analysis ...\n")
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/ballgown_analysis/"))
  }
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic"))){
    dir.create(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic"))
  }
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/"))
  }
  ###############################################
  ## Creating "ballgown_normalized_result.csv" ##
  ##############################################
  pre.de.pheno.data <- RawCountPreData(path.prefix, independent.variable, control.group, experiment.group)
  # make ballgown object
  ballgown.object <- ballgown::ballgown(dataDir = paste0(path.prefix, "gene_data/ballgown"), samplePattern = sample.pattern, pData = pre.de.pheno.data$pheno_data, meas = 'all')
  # save ballgown object
  save(ballgown.object, file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/ballgown.rda"))
  # ballgown object do statistic test
  de.statistic.result <- ballgown::stattest(ballgown.object, feature="gene",covariate = independent.variable, getFC=TRUE, meas="FPKM")
  # Save this data for DESeq2 and edgeR gene name conversion !! (and transciprt plot)
  ballgown.texpr <- ballgown::texpr(ballgown.object, 'all')
  write.csv(ballgown.texpr, file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/texpr.csv"), row.names=FALSE)
  ballgown.t2g <- ballgown::indexes(ballgown.object)$t2g
  write.csv(ballgown.t2g, file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/t2g.csv"), row.names=FALSE)
  # Convert gene id to gene name
  indices <- match(de.statistic.result$id, ballgown.texpr$gene_id)
  gene_names_for_result <- ballgown.texpr$gene_name[indices]
  de.statistic.result <- data.frame(geneNames=gene_names_for_result, de.statistic.result)
  # Filter statistic
  de.statistic.result <- de.statistic.result[(de.statistic.result$fc != 1 & !is.na(de.statistic.result$pval) & !is.na(de.statistic.result$qval)), ]
  gene.id.data.frame <- data.frame("gene.name" = de.statistic.result$geneNames)
  write.csv(gene.id.data.frame, file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/gene_name.csv"), row.names=FALSE)
  de.statistic.result$feature <- NULL; de.statistic.result$id <- NULL; de.statistic.result$geneNames <- NULL
  de.statistic.result["log2FC"] <- log2(de.statistic.result$fc)
  write.csv(de.statistic.result, file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/statistic.csv"), row.names=FALSE)
  # Creating ballgown object
  FPKM.data.frame <- ballgown::gexpr(ballgown.object)
  # Filter FPKM row sum == 0
  FPKM.data.frame <- FPKM.data.frame[rowSums(FPKM.data.frame)>0, ]
  colnames(FPKM.data.frame) <- gsub("FPKM.", "", colnames(FPKM.data.frame))
  # For control group
  control.FPKM.data.frame <- data.frame(FPKM.data.frame[,colnames(FPKM.data.frame) %in% as.character(pre.de.pheno.data$control.group.data.frame$ids)])
  colnames(control.FPKM.data.frame) <- paste0(as.character(pre.de.pheno.data$control.group.data.frame$ids), ".", control.group)
  write.csv(control.FPKM.data.frame, file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/FPKM_control.csv"), row.names=FALSE)
  # For experiment group
  experiment.FPKM.data.frame <- data.frame(FPKM.data.frame[,colnames(FPKM.data.frame) %in% as.character(pre.de.pheno.data$experiment.group.data.frame$ids)])
  colnames(experiment.FPKM.data.frame) <- paste0(as.character(pre.de.pheno.data$experiment.group.data.frame$ids), ".", experiment.group)
  write.csv(experiment.FPKM.data.frame, file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/FPKM_experiment.csv"), row.names=FALSE)
  total.data.frame <- cbind(gene.id.data.frame, control.FPKM.data.frame, experiment.FPKM.data.frame)
  total.data.frame[paste0(control.group, ".average")] <- rowMeans(control.FPKM.data.frame)
  total.data.frame[paste0(experiment.group, ".average")] <- rowMeans(experiment.FPKM.data.frame)
  total.data.frame[paste0(control.group, ".", experiment.group, ".average")]<- rowMeans(total.data.frame[-1])
  ballgown.result <- cbind(total.data.frame, de.statistic.result)
  ballgown.result <- rbind(ballgown.result[ballgown.result$gene.name != ".",], ballgown.result[ballgown.result$gene.name == ".",])
  # Filter out pval is NaN and qval is NaN
  write.csv(ballgown.result, file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_normalized_result.csv"), row.names=FALSE)

  cat(paste0("     \u25CF Selecting differential expressed genes(ballgown) ==> q-value : ", ballgown.qval, "  log2(Fold Change) : ", ballgown.log2FC, " ...\n"))
  ballgown.result.DE <- ballgown.result[(ballgown.result$log2FC>ballgown.log2FC) & (ballgown.result$qval<ballgown.qval), ]
  cat(paste0("          \u25CF Total '", length(row.names(ballgown.result.DE)), "' DEG have been found !!!\n"))
  write.csv(ballgown.result.DE, file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_normalized_DE_result.csv"), row.names=FALSE)

  ###########################
  ## ballgown visulization ##
  ###########################
  if(file.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_normalized_result.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/FPKM_control.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/FPKM_experiment.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/statistic.csv"))){
    # Transcript Related
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/"))){
      dir.create(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/"))
    }
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/transcript_related/"))){
      dir.create(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/transcript_related/"))
    }

    ###############
    #### PreDE ####
    ###############
    BallgownTranscriptRelatedPlot(path.prefix)
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/preDE/"))){
      dir.create(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/preDE/"))
    }
    # Frequency
    FrequencyPlot("ballgown_analysis", "FPKM", path.prefix, independent.variable, control.group, experiment.group)
    # Bax and Violin
    BoxViolinPlot("ballgown_analysis", "FPKM", path.prefix, independent.variable, control.group, experiment.group)
    # PCA
    PCAPlot("ballgown_analysis", "FPKM", path.prefix, independent.variable, control.group, experiment.group)
    #Correlation
    CorrelationPlot("ballgown_analysis", "FPKM", path.prefix, independent.variable, control.group, experiment.group)

    ############
    #### DE ####
    ############
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/DE/"))){
      dir.create(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/DE/"))
    }
    # Volcano
    VolcanoPlot("ballgown_analysis", "FPKM", path.prefix, independent.variable, control.group, experiment.group, ballgown.qval, ballgown.log2FC)
    # MA
    MAPlot("ballgown_analysis", "FPKM", path.prefix, independent.variable, control.group, experiment.group, ballgown.qval)
    # DE PCA plot
    DEPCAPlot("ballgown_analysis", "FPKM", path.prefix, independent.variable, control.group, experiment.group)
    # Heatmap
    DEHeatmap("ballgown_analysis", "FPKM", path.prefix, independent.variable, control.group, experiment.group)
  } else {
    stop("(\u2718) file missing ERROR.\n\n")
  }
}

# Transcription related plot
BallgownTranscriptRelatedPlot <- function(path.prefix){
  # draw for distribution of transcript count per gene
  if(file.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_normalized_result.csv"))){
    cat(paste0("\u25CF Plotting Transcript related plot\n"))
    texpr.read.csv <- read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/texpr.csv"))
    t2g.read.csv <- read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/ballgown_R_object/t2g.csv"))
    transcript_gene_table <- t2g.read.csv
    counts=table(transcript_gene_table[,"g_id"])
    c_one = length(which(counts == 1))
    c_more_than_one = length(which(counts > 1))
    c_max = max(counts)
    png(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/transcript_related/Distribution_Transcript_Count_per_Gene_Plot.png"))
    hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of Transcript Count per Gene")
    legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
    legend("topright", legend_text, lty=NULL)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/Transcript_Related/Distribution_Transcript_Count_per_Gene_Plot.png"), "' has been created. \n"))

    # draw the distribution of transcript length
    full_table <- texpr.read.csv
    t.mini.length = min(full_table$length[full_table$length > 0])
    t.max.length = max(full_table$length[full_table$length > 0])
    png(paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/transcript_related/Distribution_Transcript_Length_Plot.png"))
    hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of Transcript Lengths", col="steelblue")
    legend_text = c(paste("Minimum transcript length =", t.mini.length), paste("Maximum transcript length =", t.max.length))
    legend("topright", legend_text, lty=NULL)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/ballgown_analysis/images/Transcript_Related/Distribution_Transcript_Length_Plot.png"), "' has been created. \n\n"))
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}


