# Run ballgown analysis
BallgownAnalysis <- function(path.prefix,
                             genome.name,
                             sample.pattern,
                             independent.variable,
                             case.group,
                             control.group,
                             ballgown.pval,
                             ballgown.log2FC,
                             phenoData.result) {
  message("\n\u2618\u2618 ballgown analysis ...\n")
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/ballgown_analysis/"))
  }
  if(!dir.exists(paste0(path.prefix,
                        "RNASeq_results/ballgown_analysis/",
                        "normalized_&_statistic"))){
    dir.create(paste0(path.prefix,
                      "RNASeq_results/ballgown_analysis/",
                      "normalized_&_statistic"))
  }
  if(!dir.exists(paste0(path.prefix,
                        "RNASeq_results/ballgown_analysis/",
                        "ballgown_R_object/"))){
    dir.create(paste0(path.prefix,
                      "RNASeq_results/ballgown_analysis/",
                      "ballgown_R_object/"))
  }
  ###############################################
  ## Creating "ballgown_normalized_result.csv" ##
  ##############################################
  phenoData.result<- phenoDataWrap(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group)
  pheno_data <- phenoData.result$pheno_data
  # Order in ID!!
  pheno_data <- pheno_data[order(pheno_data$ids),]
  # Refactor to make sure 1 is case and 2 is control
  pheno_data[independent.variable][[1]] <-
    factor(as.character(pheno_data[independent.variable][[1]]),
           levels = c(case.group, control.group))
  # make ballgown object
  ballgown.object <- ballgown::ballgown(dataDir = paste0(path.prefix,
                                                         "gene_data/ballgown"),
                                        samplePattern = sample.pattern,
                                        pData = pheno_data, meas = 'all')
  # save ballgown object"
  save(ballgown.object,
       file = paste0(path.prefix,
                     "RNASeq_results/ballgown_analysis/",
                     "ballgown_R_object/ballgown.rda"))
  # Save this data for DESeq2 and edgeR gene name conversion !!
  ballgown.texpr <- ballgown::texpr(ballgown.object, 'all')
  write.csv(ballgown.texpr,
            file = paste0(path.prefix,
                          "RNASeq_results/ballgown_analysis/",
                          "ballgown_R_object/texpr.csv"),
            row.names=FALSE)
  ballgown.t2g <- ballgown::indexes(ballgown.object)$t2g
  write.csv(ballgown.t2g,
            file = paste0(path.prefix,
                          "RNASeq_results/ballgown_analysis/",
                          "ballgown_R_object/t2g.csv"),
            row.names=FALSE)
  # ballgown object do statistic test
  de.statistic.result <- ballgown::stattest(ballgown.object,
                                            feature="gene",
                                            covariate = independent.variable,
                                            getFC=TRUE,
                                            meas="FPKM")
  de.statistic.result$feature <- NULL
  de.statistic.result$id <- NULL
  de.statistic.result$gene.name <- NULL
  de.statistic.result["log2FC"] <- log2(de.statistic.result$fc)
  FPKM.data.frame <- ballgown::gexpr(ballgown.object)
  # Change name
  colnames(FPKM.data.frame) <- gsub("FPKM.", "", colnames(FPKM.data.frame))
  ballgown.result <- cbind(FPKM.data.frame, de.statistic.result)

  # save case and control data frame
  case.group.data.frame <- phenoData.result$case.group.data.frame
  control.group.data.frame <- phenoData.result$control.group.data.frame
  # save case and control group size
  case.group.size <- phenoData.result$case.group.size
  control.group.size <- phenoData.result$control.group.size

  # Select FPKM sum not 0!!
  tmp <- seq_len(case.group.size+control.group.size)
  ballgown.result <- ballgown.result[rowSums(ballgown.result[,tmp]) > 0, ]
  # Map gene id to gene name
  indices <- match(row.names(ballgown.result), ballgown.texpr$gene_id)
  gene_names_for_result <- ballgown.texpr$gene_name[indices]
  ballgown.result <- cbind("gene.name" = gene_names_for_result, ballgown.result)
  # Store with seperation fo noval gene and gene with name
  ballgown.result <- rbind(ballgown.result[ballgown.result$gene.name != ".",],
                           ballgown.result[ballgown.result$gene.name == ".",])
  # Filter out gene (p-value = Na) and
  #                 (log2FC = Inf or log2FC = Na or log2FC = -Inf)
  ballgown.result <- ballgown.result[!is.na(ballgown.result$pval) &
                                       (ballgown.result$log2FC != Inf) &
                                       !is.na(ballgown.result$log2FC) &
                                       (ballgown.result$log2FC != -Inf), ]
  # For case group
  case.ballgown.result <-
    data.frame(ballgown.result[,colnames(ballgown.result) %in%
                                 as.character(case.group.data.frame$ids)])
  colnames(case.ballgown.result) <-
    paste0(as.character(case.group.data.frame$ids),
           ".", case.group)
  # For control group
  control.ballgown.result <-
    data.frame(ballgown.result[,colnames(ballgown.result) %in%
                                 as.character(control.group.data.frame$ids)])
  colnames(control.ballgown.result) <-
    paste0(as.character(control.group.data.frame$ids),
           ".", control.group)
  # Storing whole ballgown report
  total.data.frame <- cbind("gene.name" = ballgown.result$gene.name,
                            case.ballgown.result,
                            control.ballgown.result)
  total.data.frame[paste0(case.group, ".average")] <-
    rowMeans(case.ballgown.result)
  total.data.frame[paste0(control.group, ".average")] <-
    rowMeans(control.ballgown.result)
  total.data.frame[paste0(case.group, ".", control.group, ".average")] <-
    rowMeans(cbind(case.ballgown.result, control.ballgown.result))
  ballgown.result <- cbind(total.data.frame,
                           "fc" = ballgown.result$fc,
                           "log2FC" = ballgown.result$log2FC,
                           "pval" = ballgown.result$pval,
                           "qval" = ballgown.result$qval)
  ballgown.result <- rbind(ballgown.result[ballgown.result$gene.name != ".", ],
                           ballgown.result[ballgown.result$gene.name == ".", ])
  write.csv(data.frame("gene.name" = ballgown.result[,1]),
            file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/",
                          "ballgown_R_object/gene_name.csv"),
            row.names=FALSE)
  write.csv(ballgown.result[,c(2:(case.group.size+1))],
            file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/",
                          "normalized_&_statistic/FPKM_case.csv"),
            row.names=FALSE)
  write.csv(ballgown.result[,c((2+case.group.size):
                                 (case.group.size+control.group.size+1))],
            file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/",
                          "normalized_&_statistic/FPKM_control.csv"),
            row.names=FALSE)
  write.csv(ballgown.result[,c((2+3+case.group.size+control.group.size):
                                 (length(total.data.frame)))],
            file = paste0(path.prefix, "RNASeq_results/ballgown_analysis/",
                          "normalized_&_statistic/statistic.csv"),
            row.names=FALSE)
  write.csv(ballgown.result,
            file = paste0(path.prefix,
                          "RNASeq_results/ballgown_analysis/",
                          "ballgown_normalized_result.csv"),
            row.names=FALSE)

  message("     \u25CF Selecting differential expressed genes",
          "(ballgown) ==> p-value : ", ballgown.pval,
          "  log2(Fold Change) : ", ballgown.log2FC, " ...\n")
  ballgown.result.DE <-
    ballgown.result[((ballgown.result$log2FC>ballgown.log2FC) |
                       (ballgown.result$log2FC<(-ballgown.log2FC))) &
                      (ballgown.result$pval<ballgown.pval), ]
  message("          \u25CF Total '", nrow(ballgown.result.DE),
          "' DEG have been found !!!\n")
  write.csv(ballgown.result.DE,
            file = paste0(path.prefix,
                          "RNASeq_results/ballgown_analysis/",
                          "ballgown_normalized_DE_result.csv"),
            row.names=FALSE)

  # Check ballgown.result.DE before visulization!!
  ###########################
  ## ballgown visulization ##
  ###########################
  if(file.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/",
                        "ballgown_normalized_result.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/",
                        "normalized_&_statistic/FPKM_case.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/",
                        "normalized_&_statistic/FPKM_control.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/ballgown_analysis/",
                        "normalized_&_statistic/statistic.csv"))){
    # Transcript Related
    if(!dir.exists(paste0(path.prefix,
                          "RNASeq_results/ballgown_analysis/images/"))){
      dir.create(paste0(path.prefix,
                        "RNASeq_results/ballgown_analysis/images/"))
    }
    if(!dir.exists(paste0(path.prefix,
                          "RNASeq_results/",
                          "ballgown_analysis/images/","transcript_related/"))){
      dir.create(paste0(path.prefix,
                        "RNASeq_results/",
                        "ballgown_analysis/images/transcript_related/"))
    }
    BallgownTranscriptRelatedPlot(path.prefix)
    ###############
    #### PreDE ####
    ###############
    if (nrow(ballgown.result) > 1) {
      if(!dir.exists(paste0(path.prefix,
                            "RNASeq_results/ballgown_analysis/images/preDE/"))){
        dir.create(paste0(path.prefix,
                          "RNASeq_results/ballgown_analysis/images/preDE/"))
      }
      csv.results <- ParseResultCSV("ballgown_analysis",
                                    "FPKM",
                                    path.prefix,
                                    independent.variable,
                                    case.group,
                                    control.group)
      case.normalized <- csv.results$case
      control.normalized <- csv.results$control
      independent.variable.data.frame <- cbind(case.normalized,
                                               control.normalized)
      normalized_dataset <- read.csv(paste0(path.prefix, "RNASeq_results/",
                                            "ballgown_analysis", "/",
                                            strsplit("ballgown_analysis", "_")[[1]][1],
                                            "_normalized_result.csv"))
      # Frequency
      FrequencyPlot("ballgown_analysis",
                    "FPKM",
                    path.prefix,
                    independent.variable,
                    case.group,
                    control.group,
                    independent.variable.data.frame)
      # Bax and Violin
      BoxViolinPlot("ballgown_analysis",
                    "FPKM",
                    path.prefix,
                    independent.variable,
                    case.group,
                    control.group,
                    independent.variable.data.frame,
                    phenoData.result)
      # PCA
      PCAPlot("ballgown_analysis",
              "FPKM",
              path.prefix,
              independent.variable,
              case.group,
              control.group,
              independent.variable.data.frame,
              phenoData.result)

      #Correlation
      CorrelationPlot("ballgown_analysis",
                      "FPKM",
                      path.prefix,
                      independent.variable,
                      case.group,
                      control.group,
                      independent.variable.data.frame,
                      phenoData.result)

      ############
      #### DE ####
      ############
      if(!dir.exists(paste0(path.prefix,
                            "RNASeq_results/ballgown_analysis/images/DE/"))){
        dir.create(paste0(path.prefix,
                          "RNASeq_results/ballgown_analysis/images/DE/"))
      }

      # Volcano
      VolcanoPlot("ballgown_analysis",
                  "FPKM",
                  path.prefix,
                  independent.variable,
                  case.group,
                  control.group,
                  ballgown.pval,
                  ballgown.log2FC,
                  normalized_dataset)

      # MA
      MAPlot("ballgown_analysis",
             "FPKM",
             path.prefix,
             independent.variable,
             case.group,
             control.group,
             ballgown.pval,
             csv.results,
             normalized_dataset)

      if (nrow(ballgown.result.DE) > 1) {
        # DE PCA plot
        DEPCAPlot("ballgown_analysis",
                  "FPKM",
                  path.prefix,
                  independent.variable,
                  case.group,
                  control.group,
                  normalized_dataset,
                  phenoData.result)
        # Heatmap
        DEHeatmap("ballgown_analysis",
                  "FPKM",
                  path.prefix,
                  independent.variable,
                  case.group,
                  control.group,
                  normalized_dataset,
                  phenoData.result)
      } else {
        message("(\u26A0) Less than one differential expressed gene terms ",
                "found !!! Skip DE_PCA and DE_Heatmap visualization !!! \n\n")
      }
    } else {
      message("(\u26A0) Less than one gene terms are found !!! ",
           "Skip visualization step !!! \n\n")
    }
  } else {
    message("(\u2718) necessary file is missing!! Something ERROR happend ",
            "during ballgown analysis!! Skip visualization!!\n\n")
  }
}

# Transcription related plot
BallgownTranscriptRelatedPlot <- function(path.prefix){
  # draw for distribution of transcript count per gene
  message("\u25CF Plotting Transcript related plot\n")
  texpr.read.csv <- read.csv(paste0(path.prefix,
                                    "RNASeq_results/ballgown_analysis/",
                                    "ballgown_R_object/texpr.csv"))
  t2g.read.csv <- read.csv(paste0(path.prefix,
                                  "RNASeq_results/ballgown_analysis/",
                                  "ballgown_R_object/t2g.csv"))
  transcript_gene_table <- t2g.read.csv
  count=table(transcript_gene_table[,"g_id"])
  c_one = length(which(count == 1))
  c_more_than_one = length(which(count > 1))
  c_max = max(count)
  png(paste0(path.prefix,
             "RNASeq_results/ballgown_analysis/images/transcript_related/",
             "Distribution_Transcript_Count_per_Gene_Plot.png"),
      width=5,
      height=5,
      units="in",
      res=300)
  cex.before <- par("cex")
  par(cex = 0.7)
  hist(count,
       breaks=50,
       col="bisque4",
       xlab="Transcripts per gene",
       main="Distribution of Transcript Count per Gene",
       cex.main=1.14)
  legend_text = c(paste("Genes with one transcript =", c_one),
                  paste("Genes with more than one transcript =",
                        c_more_than_one),
                  paste("Max transcripts for single gene = ", c_max))
  legend("topright", legend_text, lty=NULL)
  par(cex = cex.before)
  dev.off()
  message("(\u2714) : '", path.prefix,
          "RNASeq_results/ballgown_analysis/images/Transcript_Related/",
          "Distribution_Transcript_Count_per_Gene_Plot.png",
          "' has been created. \n")

  # draw the distribution of transcript length
  full_table <- texpr.read.csv
  t.mini.length = min(full_table$length[full_table$length > 0])
  t.max.length = max(full_table$length[full_table$length > 0])
  png(paste0(path.prefix,
             "RNASeq_results/ballgown_analysis/images/",
             "transcript_related/Distribution_Transcript_Length_Plot.png"),
      width=5,
      height=5,
      units="in",
      res=300)
  cex.before <- par("cex")
  par(cex = 0.7)
  hist(full_table$length, breaks=50, xlab="Transcript length (bp)",
       main="Distribution of Transcript Lengths",
       col="steelblue",
       cex.main=1.14)
  legend_text = c(paste("Minimum transcript length =", t.mini.length),
                  paste("Maximum transcript length =", t.max.length))
  legend("topright", legend_text, lty=NULL)
  par(cex = cex.before)
  dev.off()
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/ballgown_analysis/images/",
          "Transcript_Related/Distribution_Transcript_Length_Plot.png",
          "' has been created. \n\n")
}


