edgeRRawCountAnalysis <- function(path.prefix,
                                  independent.variable,
                                  case.group,
                                  control.group,
                                  edgeR.pval,
                                  edgeR.log2FC) {
  message("\n\u2618\u2618 edgeR analysis ...\n")
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/edgeR_analysis"))){
    dir.create(paste0(path.prefix, "RNASeq_results/edgeR_analysis"))
  }
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                        "normalized_&_statistic"))){
    dir.create(paste0(path.prefix,
                      "RNASeq_results/edgeR_analysis/normalized_&_statistic"))
  }
  #############################################
  ## Creating "edgeR_normalized_result.csv" ##
  ############################################
  pre.de.pheno.data <- RawCountPreData(path.prefix,
                                       independent.variable,
                                       case.group,
                                       control.group)
  raw.count <- pre.de.pheno.data$gene.count.matrix
  raw.count.gene.name <- pre.de.pheno.data$gene.count.name

  # create DGEList object (edgeR)
  message("\u25CF Creating 'DGEList' object from count matrix ... \n")
  gene.data.frame <- data.frame(gene.name = raw.count.gene.name)
  pheno.data <- pre.de.pheno.data$pheno_data
  pheno.data <- pheno.data[order(pheno.data$ids),]
  pheno.data[independent.variable][[1]] <-
    factor(as.character(pheno.data[independent.variable][[1]]),
           levels = c(case.group, control.group))

  deglist.object <-
    edgeR::DGEList(counts=raw.count,
                   group = pheno.data[independent.variable][[1]],
                   genes = raw.count.gene.name)
  # Normalization with TMM (trimmed mean of M-values )
  message("     \u25CF Normalizing DGEList object (TMM) ... \n")
  deglist.object <- edgeR::calcNormFactors(deglist.object, method="TMM")
  deglist.object$samples$normalized.lib.size <-
    deglist.object$samples$lib.size*deglist.object$samples$norm.factors
  write.csv(deglist.object$samples,
            file = paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                          "normalized_&_statistic/TMM.csv"))
  # estimating Dispersions
  # quantile-adjusted conditional maximum likelihood (qCML)
  # method for experiments with single factor.
  dgList <- edgeR::estimateCommonDisp(deglist.object)
  dgList <- edgeR::estimateTagwiseDisp(dgList)
  # Testing for DE genes
  de.statistic.result <- edgeR::exactTest(dgList)
  colnames(de.statistic.result$table)[1] <- "log2FC"
  colnames(de.statistic.result$table)[3] <- "pval"
  # run the cpm function on a DGEList object which is normalisation by
  # TMM factors ==> get TMM normalized count !!
  # count were first normalized with TMM and then be presented as cpm !
  normalized.count.table <- edgeR::cpm(dgList, normalized.lib.sizes=TRUE)

  # For case group
  case.id <- as.character(pre.de.pheno.data$case.group.data.frame$ids)
  case.cpm.data.frame <-
    data.frame(normalized.count.table[,colnames(normalized.count.table) %in%
                                        case.id])
  colnames(case.cpm.data.frame) <- paste0(case.id, ".", case.group)

  # For control group
  control.id <- as.character(pre.de.pheno.data$control.group.data.frame$ids)
  control.cpm.data.frame <-
    data.frame(normalized.count.table[,colnames(normalized.count.table) %in%
                                        control.id])
  colnames(control.cpm.data.frame) <- paste0(control.id, ".", control.group)
  # create whole data.frame
  total.data.frame <- cbind(gene.data.frame,
                            case.cpm.data.frame,
                            control.cpm.data.frame)
  total.data.frame[paste0(case.group, ".average")] <-
    rowMeans(case.cpm.data.frame)
  total.data.frame[paste0(control.group, ".average")] <-
    rowMeans(control.cpm.data.frame)
  total.data.frame[paste0(case.group, ".", control.group, ".average")]<-
    rowMeans(cbind(case.cpm.data.frame, control.cpm.data.frame))
  edgeR.result <- cbind(total.data.frame, de.statistic.result$table)
  edgeR.result <- rbind(edgeR.result[edgeR.result$gene.name != ".", ],
                        edgeR.result[edgeR.result$gene.name == ".", ])
  # Filter out gene (p-value = Na) and
  #                 (log2FC = Inf or log2FC = Na or log2FC = -Inf)
  edgeR.result <- edgeR.result[!is.na(edgeR.result$pval) &
                                 (edgeR.result$log2FC != Inf) &
                                 !is.na(edgeR.result$log2FC) &
                                 (edgeR.result$log2FC != -Inf), ]
  case.group.size <- pre.de.pheno.data$case.group.size
  control.group.size <- pre.de.pheno.data$control.group.size
  # Write out csv file
  write.csv(edgeR.result[,c(2:(case.group.size+1))],
            file = paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                          "normalized_&_statistic/TMM&CPM_case.csv"),
            row.names=FALSE)
  write.csv(edgeR.result[,c((2+case.group.size):
                              (case.group.size+control.group.size+1))],
            file = paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                          "normalized_&_statistic/TMM&CPM_control.csv"),
            row.names=FALSE)
  write.csv(edgeR.result[,c((2+3+case.group.size+control.group.size):
                              (length(edgeR.result)))],
            file = paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                          "normalized_&_statistic/statistic.csv"),
            row.names=FALSE)
  # Write result into file (csv)
  write.csv(edgeR.result,
            file = paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                          "edgeR_normalized_result.csv"),
            row.names=FALSE)

  # Slecet DE genes (condition)
  message(paste0("     \u25CF Selecting differential ",
                 "expressed genes(edgeR) ==> p-value : ",
                 edgeR.pval, "  log2(Fold Change) : ", edgeR.log2FC, " ...\n"))
  edgeR.result.DE <- edgeR.result[((edgeR.result$log2FC>edgeR.log2FC) |
                                     (edgeR.result$log2FC<(-edgeR.log2FC))) &
                                    (edgeR.result$pval<edgeR.pval), ]
  message("          \u25CF Total '",
          length(row.names(edgeR.result.DE)), "' DEG have been found !!")
  write.csv(edgeR.result.DE,
            file = paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                          "edgeR_normalized_DE_result.csv"),
            row.names=FALSE)

  # Check edgeR.result.DE before visulization!!

    ########################
    ## edgeR visulization ##
    ########################
  # PreDE
  if(file.exists(paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                        "edgeR_normalized_result.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                        "normalized_&_statistic/TMM&CPM_case.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                        "normalized_&_statistic/TMM&CPM_control.csv")) &&
     file.exists(paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                        "normalized_&_statistic/statistic.csv"))){
    if(!dir.exists(paste0(path.prefix,
                          "RNASeq_results/edgeR_analysis/images"))){
      dir.create(paste0(path.prefix, "RNASeq_results/edgeR_analysis/images"))
    }
    if (nrow(edgeR.result) > 1) {
      ###############
      #### PreDE ####
      ###############
      if(!dir.exists(paste0(path.prefix,
                            "RNASeq_results/edgeR_analysis/images/preDE/"))){
        dir.create(paste0(path.prefix,
                          "RNASeq_results/edgeR_analysis/images/preDE/"))
      }
      # Frequency
      FrequencyPlot("edgeR_analysis",
                    "TMM&CPM",
                    path.prefix,
                    independent.variable,
                    case.group,
                    control.group)
      # Bax and Violin
      BoxViolinPlot("edgeR_analysis",
                    "TMM&CPM",
                    path.prefix,
                    independent.variable,
                    case.group,
                    control.group)
      # PCA
      PCAPlot("edgeR_analysis",
              "TMM&CPM",
              path.prefix,
              independent.variable,
              case.group,
              control.group)
      #Correlation
      CorrelationPlot("edgeR_analysis",
                      "TMM&CPM",
                      path.prefix,
                      independent.variable,
                      case.group,
                      control.group)

      # MDS plot
      message("\u25CF Plotting MDS plot ... \n")
      png(paste0(path.prefix,
                 "RNASeq_results/edgeR_analysis/",
                 "images/preDE/MDS_Plot_edgeR.png"),
          width=5,
          height=5,
          units="in",
          res=300)
      my_colors=c(rgb(50, 147, 255,maxColorValue = 255),
                  rgb(255, 47, 35,maxColorValue = 255))
      cex.before <- par("cex")
      par(xpd=TRUE, cex = 0.7)
      edgeR::plotMDS.DGEList(deglist.object,
                             top = 1000,
                             labels = NULL,
                             col = my_colors[as.numeric(deglist.object$
                                                          samples$group)],
                             pch = 20)
      legend("bottomright", inset=c(0,1), horiz=TRUE, bty = "n",
             legend = levels(deglist.object$samples$group) ,
             col = my_colors, pch=20)
      par(cex = 0.8)
      title("MDS Plot")
      par(cex = cex.before)
      dev.off()
      message(paste0("(\u2714) : '",
                     paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                            "images/preDE/MDS_Plot_edgeR.png"),
                     "' has been created. \n\n"))

      # MeanVar plot
      message("\u25CF Plotting MeanVar plot ... \n")
      png(paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                 "images/preDE/MeanVar_Plot_edgeR.png"),
          width=5,
          height=5,
          units="in",
          res=300)
      cex.before <- par("cex")
      par(xpd=TRUE, cex = 0.7)
      edgeR::plotMeanVar(dgList, show.tagwise.vars=TRUE, NBline=TRUE)
      par(cex = 0.8)
      graphics::title("Mean-Variance Plot")
      par(cex = cex.before)
      dev.off()
      message(paste0("(\u2714) : '",
                     paste0(path.prefix, "RNASeq_results/edgeR_analysis/",
                            "images/preDE/MeanVar_Plot_edgeR.png"),
                     "' has been created. \n\n"))
      # BCV plot
      message("\u25CF Plotting BCV plot ...\n")
      png(paste0(path.prefix, "RNASeq_results/edgeR_analysis/images/preDE/",
                 "BCV_Plot_edgeR.png"),
          width=5,
          height=5,
          units="in",
          res=300)
      cex.before <- par("cex")
      par(xpd=TRUE, cex = 0.7)
      edgeR::plotBCV(dgList)
      par(cex = 0.8)
      graphics::title("BCV Plot")
      par(cex = cex.before)
      dev.off()
      message(paste0("(\u2714) : '",
                     paste0(path.prefix, "RNASeq_results/",
                            "edgeR_analysis/images/preDE/BCV_Plot_edgeR.png"),
                     "' has been created. \n\n"))
      ############
      #### DE ####
      ############
      if(!dir.exists(paste0(path.prefix,
                            "RNASeq_results/edgeR_analysis/images/DE/"))){
        dir.create(paste0(path.prefix,
                          "RNASeq_results/edgeR_analysis/images/DE/"))
      }
      # Volcano
      VolcanoPlot("edgeR_analysis",
                  "TMM&CPM",
                  path.prefix,
                  independent.variable,
                  case.group,
                  control.group,
                  edgeR.pval,
                  edgeR.log2FC)

      # Plot Smear plot !!
      fit <- edgeR::glmFit(dgList)
      lrt <- edgeR::glmLRT(fit, coef=2)
      png(paste0(path.prefix, "RNASeq_results/edgeR_analysis/images/DE/",
                 "Smear_Plot_edgeR.png"),
          width=5,
          height=5,
          units="in",
          res=300)
      cex.before <- par("cex")
      par(xpd=TRUE, cex = 0.7)
      # rownames for genes identified as being differentially expressed;
      # use exactTest or glmLRT to identify DE genes. Note that `tag' and `gene'
      # are synonymous here.
      edgeR::plotSmear(lrt, de.tags = dgList$genes)
      par(cex = 0.8)
      graphics::title(paste0(case.group, " vs ", control.group))
      par(cex = cex.before)
      dev.off()
      message(paste0("(\u2714) : '",
                     paste0(path.prefix, "RNASeq_results/",
                            "edgeR_analysis/images/preDE/Smear_Plot_edgeR.png"),
                     "' has been created. \n\n"))

      if (nrow(edgeR.result.DE) > 1) {
        # PCA plot
        DEPCAPlot("edgeR_analysis",
                  "TMM&CPM",
                  path.prefix,
                  independent.variable,
                  case.group,
                  control.group)

        # Heatmap
        DEHeatmap("edgeR_analysis",
                  "TMM&CPM",
                  path.prefix,
                  independent.variable,
                  case.group,
                  control.group)

      } else {
        cat ("(\u26A0) Less than one differential expressed gene term found !!",
             " Skip DE_PCA and DE_Heatmap visualization !!! \n\n")
      }
    } else {
      cat ("(\u26A0) Less than one gene terms found !!! ",
           "Skip visualization step !!! \n\n")
    }
  } else {
    message("(\u2718) necessary file is missing!! Something ERROR happend ",
            "during edgeR analysis!! Skip visualization!!\n\n")
  }
}
