DESeq2RawCountAnalysis <- function(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group,
                                   DESeq2.pval,
                                   DESeq2.log2FC,
                                   phenoData.result) {
  message("\n\u2618\u2618 DESeq2 analysis ...\n")
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/DESeq2_analysis"))){
    dir.create(paste0(path.prefix, "RNASeq_results/DESeq2_analysis"))
  }
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/DESeq2_analysis/",
                        "normalized_&_statistic"))){
    dir.create(paste0(path.prefix,
                      "RNASeq_results/DESeq2_analysis/normalized_&_statistic"))
  }

  #############################################
  ## Creating "DESeq2_normalized_result.csv" ##
  #############################################
  phenoData.result<- phenoDataWrap(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group)
  pheno.data <- phenoData.result$pheno_data
  rawCount.result <- RawCountWrap(path.prefix)
  raw.count <- rawCount.result$gene.count.matrix
  raw.count.gene.name <- rawCount.result$gene.count.name
  # creatin gene name data frame
  gene.data.frame <- data.frame(gene.name = raw.count.gene.name)
  # create design data.frame (independent.variable)

  # Order in ID!!
  pheno.data <- pheno.data[order(pheno.data$ids),]
  colData <- pheno.data[independent.variable]
  row.names(colData) <- pheno.data$ids
  colnames(colData) <- "independent.variable"
  colData$independent.variable <- factor(as.character(colData$
                                                        independent.variable),
                                         levels = c(case.group, control.group))
  # creat DESeqDataSet
  # Rows of colData correspond to columns of countData
  message("\u25CF Creating 'DESeqDataSet' object from count matrix ... \n")
  dds.mat <- DESeq2::DESeqDataSetFromMatrix(countData = raw.count,
                                            colData = colData,
                                            design =  ~independent.variable)
  # Run DESeq() function
  # 1.estimation of size factors, 2.estimation of dispersion,
  # 3.gene-wise dispersion estimates, 4.mean-dispersion relationship,
  # 5.final dispersion estimates, 6.fitting model and testing
  # Negative Binomial GLM fitting and Wald statistics
  dds_de <- DESeq2::DESeq(dds.mat)
  # creating statistic result
  statistic.res <- DESeq2::results(dds_de,
                                   contrast = c("independent.variable",
                                                case.group,
                                                control.group))
  colnames(statistic.res)[2] <- "log2FC"
  colnames(statistic.res)[5] <- "pval"
  # Normalization method of DESeq2 is Median Ratio Normalization (MRN)
  normalized.count.table <- DESeq2::counts(dds_de, normalized=TRUE)

  # For case group
  case.id <- as.character(phenoData.result$case.group.data.frame$ids)
  case.mrn.data.frame <-
    data.frame(normalized.count.table[,colnames(normalized.count.table) %in%
                                        case.id])
  colnames(case.mrn.data.frame) <- paste0(case.id, ".", case.group)

  # For control group
  control.id <- as.character(phenoData.result$control.group.data.frame$ids)
  control.mrn.data.frame <-
    data.frame(normalized.count.table[,colnames(normalized.count.table) %in%
                                        control.id])
  colnames(control.mrn.data.frame) <- paste0(control.id, ".", control.group)

  # create whole data.frame
  total.data.frame <- cbind(gene.data.frame,
                            case.mrn.data.frame,
                            control.mrn.data.frame)
  total.data.frame[paste0(case.group, ".average")] <-
    rowMeans(case.mrn.data.frame)
  total.data.frame[paste0(control.group, ".average")] <-
    rowMeans(control.mrn.data.frame)
  total.data.frame[paste0(case.group, ".", control.group, ".average")]<-
    rowMeans(cbind(case.mrn.data.frame, control.mrn.data.frame))
  DESeq2.result <- cbind(total.data.frame, statistic.res)
  DESeq2.result <- rbind(DESeq2.result[DESeq2.result$gene.name != ".", ],
                         DESeq2.result[DESeq2.result$gene.name == ".", ])

  # Filter out gene (p-value = Na) and
  #                 (q-value = Na) and
  #                 (log2FC = Inf or log2FC = Na or log2FC = -Inf)
  DESeq2.result <- DESeq2.result[!is.na(DESeq2.result$pval) &
                                   !is.na(DESeq2.result$padj) &
                                   (DESeq2.result$log2FC != Inf) &
                                   !is.na(DESeq2.result$log2FC) &
                                   (DESeq2.result$log2FC != -Inf), ]
  # Write result into file (csv)
  case.group.size <- phenoData.result$case.group.size
  control.group.size <- phenoData.result$control.group.size
  write.csv(DESeq2.result[,c(2:(case.group.size+1))],
            file = paste0(path.prefix,
                          "RNASeq_results/DESeq2_analysis/",
                          "normalized_&_statistic/MRN_case.csv"),
            row.names=FALSE)
  write.csv(DESeq2.result[,c((2+case.group.size):
                               (case.group.size+control.group.size+1))],
            file = paste0(path.prefix, "RNASeq_results/DESeq2_analysis/",
                          "normalized_&_statistic/MRN_control.csv"),
            row.names=FALSE)
  write.csv(DESeq2.result[,c((2+3+case.group.size+control.group.size):
                               (length(DESeq2.result)))],
            file = paste0(path.prefix, "RNASeq_results/DESeq2_analysis/",
                          "normalized_&_statistic/statistic.csv"),
            row.names=FALSE)
  write.csv(DESeq2.result,
            file = paste0(path.prefix, "RNASeq_results/DESeq2_analysis/",
                          "DESeq2_normalized_result.csv"),
            row.names=FALSE)

  message("     \u25CF Selecting differential expressed genes",
          "(DESeq2) ==> padj-value : ",
          DESeq2.pval, "  log2(Fold Change) : ",
          DESeq2.log2FC, " ...\n")
  DESeq2.result.DE <- DESeq2.result[((DESeq2.result$log2FC>DESeq2.log2FC) |
                                       (DESeq2.result$log2FC<(-DESeq2.log2FC)))&
                                      (DESeq2.result$pval<DESeq2.pval), ]
  message("          \u25CF Total '",
          length(row.names(DESeq2.result.DE)),
          "' DEG have been found !!")
  write.csv(DESeq2.result.DE,
            file = paste0(path.prefix, "RNASeq_results/DESeq2_analysis/",
                          "DESeq2_normalized_DE_result.csv"),
            row.names=FALSE)


  #########################
  ## DESeq2 visulization ##
  #########################
  if(file.exists(paste0(path.prefix, "RNASeq_results/DESeq2_analysis/",
                        "DESeq2_normalized_result.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/DESeq2_analysis/",
                        "normalized_&_statistic/MRN_case.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/DESeq2_analysis/",
                        "normalized_&_statistic/MRN_control.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/DESeq2_analysis/",
                        "normalized_&_statistic/statistic.csv"))){
    if(!dir.exists(paste0(path.prefix,
                          "RNASeq_results/DESeq2_analysis/images/"))){
      dir.create(paste0(path.prefix, "RNASeq_results/DESeq2_analysis/images/"))
    }
    if (nrow(DESeq2.result) > 1) {
      ###############
      #### PreDE ####
      ###############
      if(!dir.exists(paste0(path.prefix,
                            "RNASeq_results/DESeq2_analysis/images/preDE/"))){
        dir.create(paste0(path.prefix,
                          "RNASeq_results/DESeq2_analysis/images/preDE/"))
      }

      csv.results <- ParseResultCSV("DESeq2_analysis",
                                    "MRN",
                                    path.prefix,
                                    independent.variable,
                                    case.group,
                                    control.group)
      case.normalized <- csv.results$case
      control.normalized <- csv.results$control
      independent.variable.data.frame <- cbind(case.normalized,
                                               control.normalized)
      normalized_dataset <- read.csv(paste0(path.prefix, "RNASeq_results/",
                                            "DESeq2_analysis", "/",
                                            strsplit("DESeq2_analysis", "_")[[1]][1],
                                            "_normalized_result.csv"))

      # Frequency
      FrequencyPlot("DESeq2_analysis",
                    "MRN",
                    path.prefix,
                    independent.variable,
                    case.group,
                    control.group,
                    independent.variable.data.frame)
      # Bax and Violin
      BoxViolinPlot("DESeq2_analysis",
                    "MRN",
                    path.prefix,
                    independent.variable,
                    case.group,
                    control.group,
                    independent.variable.data.frame,
                    phenoData.result)
      # PCA
      PCAPlot("DESeq2_analysis",
              "MRN",
              path.prefix,
              independent.variable,
              case.group,
              control.group,
              independent.variable.data.frame,
              phenoData.result)

      #Correlation
      CorrelationPlot("DESeq2_analysis",
                      "MRN",
                      path.prefix,
                      independent.variable,
                      case.group,
                      control.group,
                      independent.variable.data.frame,
                      phenoData.result)
      # dispersion plot
      message("\u25CF Plotting  Dispersion plot\n")
      png(paste0(path.prefix,
                 "RNASeq_results/DESeq2_analysis/",
                 "images/preDE/Dispersion_Plot_DESeq2.png"),
          width=5,
          height=5,
          units="in",
          res=300)
      cex.before <- par("cex")
      par(xpd=TRUE, cex = 0.7)
      plotDispEsts(dds_de)
      par(cex = 0.8)
      title("Dispersion plot")
      par(cex = cex.before)
      dev.off()
      message("(\u2714) : 'RNASeq_results/DESeq2_analysis/images/DE/",
              "Dispersion_Plot_DESeq2.png' has been created. \n\n")

      ############
      #### DE ####
      ############
      if(!dir.exists(paste0(path.prefix,
                            "RNASeq_results/DESeq2_analysis/images/DE/"))){
        dir.create(paste0(path.prefix,
                          "RNASeq_results/DESeq2_analysis/images/DE/"))
      }
      # Volcano
      VolcanoPlot("DESeq2_analysis",
                  "MRN",
                  path.prefix,
                  independent.variable,
                  case.group,
                  control.group,
                  DESeq2.pval,
                  DESeq2.log2FC,
                  normalized_dataset)

      # MA plot
      message("\u25CF Plotting  MA plot\n")
      png(paste0(path.prefix,"RNASeq_results/DESeq2_analysis/",
                 "images/DE/MA_Plot_DESeq2.png"),
          width=5,
          height=5,
          units="in",
          res=300)
      cex.before <- par("cex")
      par(xpd=TRUE, cex = 0.7)
      DESeq2::plotMA(dds_de)
      par(cex = 0.8)
      title("MA (MD) Plot")
      par(cex = cex.before)
      dev.off()
      message("(\u2714) : '",path.prefix, "RNASeq_results/",
              "DESeq2_analysis/images/DE/MA_Plot_DESeq2.png",
              "' has been created. \n\n")

      # Check DESeq2.result.DE before visulization!!
      if (nrow(DESeq2.result.DE) > 1) {
        # PCA plot
        DEPCAPlot("DESeq2_analysis",
                  "MRN",
                  path.prefix,
                  independent.variable,
                  case.group,
                  control.group,
                  normalized_dataset,
                  phenoData.result)

        # Heatmap
        DEHeatmap("DESeq2_analysis",
                  "MRN",
                  path.prefix,
                  independent.variable,
                  case.group,
                  control.group,
                  normalized_dataset,
                  phenoData.result)

      } else {
        message("(\u26A0) Less than one differential expressed gene term ",
                "found !!! Skip DE_PCA and DE_Heatmap visualization !!! \n\n")
      }
    } else {
      message("(\u26A0) Less than one gene terms found !!! ",
              "Skip visualization step !!! \n\n")
    }
  } else {
    message("(\u2718) necessary file is missing!! Something ERROR ",
            "happend during DESeq2 analysis!! Skip visualization!!\n\n")
  }
}
