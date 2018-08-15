# Run ballgown analysi
BallgownPreprocess <- function(path.prefix, genome.name, sample.pattern, independent.variable, control.group, experiment.group, ballgown.log2FC, ballgown.qval) {
  # this ballgown function is only for two group
  CheckOperatingSystem(FALSE)
  results <- ProgressGenesFiles(path.prefix = path.prefix, genome.name = genome.name, sample.pattern = sample.pattern, print = FALSE)
  if (isTRUE(results$phenodata.file.df) && results$ballgown.dirs.number.df != 0){
    # sorting 'pheno_data'
    cat(paste0("************** Ballgown data preprocessing **************\n"))
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/ballgown_analysis/"))
    }
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/normalized_&_statistic"))){
      dir.create(paste0(path.prefix, "RNAseq_results/ballgown_analysis/normalized_&_statistic"))
    }
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/ballgown_R_object/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/ballgown_analysis/ballgown_R_object/"))
    }
    #############################################
    ## Creating "edgeR_normalized_result.csv" ##
    ############################################
    pre.de.pheno.data <- RawCountPreData(path.prefix, independent.variable, control.group, experiment.group)
    # make ballgown object
    cat("\u25CF 3. Making ballgown object : \n")
    ballgown.object <- ballgown::ballgown(dataDir = paste0(path.prefix, "gene_data/ballgown"), samplePattern = sample.pattern, pData = pre.de.pheno.data$pheno_data, meas = 'all')
    # ballgown object do statistic test
    de.statistic.result <- ballgown::stattest(ballgown.object, feature="gene",covariate = independent.variable, getFC=TRUE, meas="FPKM")
    # Filter statistic
    de.statistic.result <- de.statistic.result[(de.statistic.result$fc != 1 & !is.na(de.statistic.result$pval) & !is.na(de.statistic.result$qval)), ]
    gene.id.data.frame <- data.frame("gene_id" = de.statistic.result$id)
    de.statistic.result$feature <- NULL; de.statistic.result$id <- NULL
    de.statistic.result["log2FC"] <- log2(de.statistic.result$fc)
    write.csv(de.statistic.result, file = paste0(path.prefix, "RNAseq_results/ballgown_analysis/normalized_&_statistic/statistic.csv"), row.names=FALSE)
    # save ballgown object
    bg <- ballgown.object
    save(bg, file = paste0(path.prefix, "RNAseq_results/ballgown_analysis/ballgown_R_object/ballgown.rda"))
    # Creating ballgown object
    FPKM.data.frame <- ballgown::gexpr(ballgown.object)
    # Filter FPKM row sum == 0
    FPKM.data.frame <- FPKM.data.frame[rowSums(FPKM.data.frame)>0, ]
    colnames(FPKM.data.frame) <- gsub("FPKM.", "", colnames(FPKM.data.frame))
    # For control group
    control.FPKM.data.frame <- data.frame(FPKM.data.frame[,colnames(FPKM.data.frame) %in% as.character(pre.de.pheno.data$control.group.data.frame$ids)])
    colnames(control.FPKM.data.frame) <- paste0(as.character(pre.de.pheno.data$control.group.data.frame$ids), ".", control.group)
    write.csv(control.FPKM.data.frame, file = paste0(path.prefix, "RNAseq_results/ballgown_analysis/normalized_&_statistic/FPKM_control.csv"), row.names=FALSE)
    # For experiment group
    experiment.FPKM.data.frame <- data.frame(FPKM.data.frame[,colnames(FPKM.data.frame) %in% as.character(pre.de.pheno.data$experiment.group.data.frame$ids)])
    colnames(experiment.FPKM.data.frame) <- paste0(as.character(pre.de.pheno.data$experiment.group.data.frame$ids), ".", experiment.group)
    write.csv(experiment.FPKM.data.frame, file = paste0(path.prefix, "RNAseq_results/ballgown_analysis/normalized_&_statistic/FPKM_experiment.csv"), row.names=FALSE)
    total.data.frame <- cbind(gene.id.data.frame, control.FPKM.data.frame, experiment.FPKM.data.frame)
    total.data.frame[paste0(control.group, ".average")] <- rowMeans(control.FPKM.data.frame)
    total.data.frame[paste0(experiment.group, ".average")] <- rowMeans(experiment.FPKM.data.frame)
    total.data.frame[paste0(control.group, "+", experiment.group, ".average")]<- rowMeans(total.data.frame[-1])
    ballgown_result <- cbind(total.data.frame, de.statistic.result)
    # Filter out pval is NaN and qval is NaN
    write.csv(ballgown_result, file = paste0(path.prefix, "RNAseq_results/ballgown_analysis/ballgown_normalized_result.csv"), row.names=FALSE)

    ###########################
    ## ballgown visulization ##
    ###########################
    if(file.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/ballgown_normalized_result.csv")) &&
       file.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/normalized_&_statistic/FPKM_control.csv")) &&
       file.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/normalized_&_statistic/FPKM_experiment.csv")) &&
       file.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/normalized_&_statistic/statistic.csv"))){
      # Transcript Related
      if(!dir.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/"))){
        dir.create(paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/"))
      }
      if(!dir.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/transcript_related/"))){
        dir.create(paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/transcript_related/"))
      }
      BallgownTranscriptRelatedPlot

      # PreDE
      if(!dir.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/preDE/"))){
        dir.create(paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/preDE/"))
      }
      # Frequency

      # Bax and Violin
      # PCA
      #Correlation
      #Volcano
      # MA
    } else {
      stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
    }
  }
}

# Transcription related plot
BallgownTranscriptRelatedPlot <- function(path.prefix){
  # draw for distribution of transcript count per gene
  if(file.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/ballgown_normalized_result.csv"))){
    # load gene name for further usage
    if (is.null(pkg.ballgown.data$bg_chrX)) {
      LoadBallgownObject(path.prefix)
    } else {
      cat(paste0("\u25CF Plotting Transcript related plot\n"))
      transcript_gene_table <- ballgown::indexes(pkg.ballgown.data$bg_chrX)$t2g
      counts=table(transcript_gene_table[,"g_id"])
      c_one = length(which(counts == 1))
      c_more_than_one = length(which(counts > 1))
      c_max = max(counts)
      png(paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/transcript_related/Distribution_Transcript_Count_per_Gene_Plot.png"))
      hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of Transcript Count per Gene")
      legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
      legend("topright", legend_text, lty=NULL)
      dev.off()
      cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/Transcript_Related/Distribution_transcript_count_per_gene_plot.png"), "' has been created. \n"))

      # draw the distribution of transcript length
      full_table <- ballgown::texpr(pkg.ballgown.data$bg_chrX, 'all')
      t.mini.length = min(full_table$length[full_table$length > 0])
      t.max.length = max(full_table$length[full_table$length > 0])
      png(paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/transcript_related/Distribution_Transcript_Length_Plot.png"))
      hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of Transcript Lengths", col="steelblue")
      legend_text = c(paste("Minimum transcript length =", t.mini.length), paste("Maximum transcript length =", t.max.length))
      legend("topright", legend_text, lty=NULL)
      dev.off()
      cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/ballgown_analysis/images/Transcript_Related/Distribution_transcript_length_plot.png"), "' has been created. \n\n"))
    }
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

LoadBallgownObject <- function(path.prefix) {
  if(isTRUE(file.exists(paste0(path.prefix, "RNAseq_results/ballgown_analysis/ballgown_R_object/ballgown.rda")))) {
    load(paste0(path.prefix, "RNAseq_results/ballgown_analysis/ballgown_R_object/ballgown.rda"))
    pkg.ballgown.data$bg_chrX <- bg
  } else {
    stop(paste0("(\u2718) '", paste0(path.prefix, "gene_data/ballgown/ballgown.rda"), "' haven't created yet. Please run \"BallgownPreprocess(genome.name, sample.pattern, independent.variable)\" first.\n\n"))
  }
}

#' Check ballgown object
#' @export
#' @param path.prefix path prefix of 'gene_data/', 'RNAseq_bin/', 'RNAseq_results/', 'Rscript/' and 'Rscript_out/' directories
#'
#' @return None
CheckBallgownObject <- function(path.prefix) {
  print(pkg.ballgown.data$bg_chrX)
  print(pkg.ballgown.data$bg_chrX_filt)
}


