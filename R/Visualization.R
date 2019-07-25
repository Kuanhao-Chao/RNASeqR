#################################
#### Alignment visualization ####
#################################
AlignmentPlot <- function(path.prefix,
                          independent.variable,
                          case.group,
                          control.group) {
  Alignment_report_reads <- read.csv(paste0(path.prefix,
                                            "RNASeq_results/Alignment_Report/",
                                            "Alignment_report_reads.csv"),
                                     row.names = 1, header = TRUE)
  Overall.map.rates <- read.csv(paste0(path.prefix,
                                       "RNASeq_results/Alignment_Report/",
                                       "Overall_Mapping_rate.csv"),
                                row.names = 1, header = TRUE)
  Alignment_report_reads$samples <- row.names(Alignment_report_reads)
  melted <- reshape2::melt(Alignment_report_reads, id.vars = c("samples"))
  melted$value <- as.numeric(melted$value)
  ggplot(data=melted,
         aes(x = melted$samples,
             y = melted$value,
             fill = melted$variable)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw() +
    geom_text(aes(label=value,
                  angle = 90),
              size = 3,
              position = position_dodge(width=1),
              check_overlap = TRUE) +
    xlab("Samples") + ylab("Mapping Reads") +
    ggtitle("Bar Plot (Each condition reads)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position="top",
          legend.text = element_text(size = 6)) +
    scale_fill_discrete(name = "Conditions")
  ggsave(paste0(path.prefix,
                "RNASeq_results/",
                "Alignment_Report/Alignment_Result_ggplot2.png"),
         dpi = 300,
         width = 7,
         height = 7)

  phenoData.result<- phenoDataWrap(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group)
  my_colors=c("#00AFBB", "#E7B800")
  color.group <- c(rep(1, phenoData.result$case.group.size),
                   rep(2, phenoData.result$case.group.size))
  Overall.map.rates$samples <- row.names(Overall.map.rates)
  colnames(Overall.map.rates) <- c("rates", "samples")
  ggplot(data=Overall.map.rates,
         aes(x = Overall.map.rates$samples,
             y = Overall.map.rates$rates)) +
    geom_bar(stat = "identity", fill=my_colors[as.numeric(color.group)]) +
    theme_bw() + xlab("Samples") + ylab("Mapping Rates") +
    ggtitle("Bar Plot (Over all Mapping Rates)") +
    geom_text(aes(label = Overall.map.rates$rates,
                  angle = 90),
              size = 4,
              position = position_dodge(width=1),
              check_overlap = TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  ggsave(paste0(path.prefix,
                "RNASeq_results/",
                "Alignment_Report/Overall_Mapping_rate_ggplot2.png"),
         dpi = 300,
         width = 7,
         height = 7)
}

###################################################
#### Pre-differential expression visualization ####
###################################################
# Frequency Plot
FrequencyPlot <- function(which.analysis,
                          which.count.normalization,
                          path.prefix,
                          independent.variable,
                          case.group,
                          control.group) {
  message("\u25CF Plotting  Frequency plot\n")
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                        which.analysis, "/images/preDE/Frequency"))){
    dir.create(paste0(path.prefix, "RNASeq_results/",
                      which.analysis, "/images/preDE/Frequency"))
  }
  csv.results <- ParseResultCSV(which.analysis,
                                which.count.normalization,
                                path.prefix,
                                independent.variable,
                                case.group,
                                control.group)
  case.normalized <- csv.results$case
  control.normalized <- csv.results$control
  independent.variable.data.frame <- cbind(case.normalized, control.normalized)
  rafalib::mypar(1, 1)
  sample.size <- length(independent.variable.data.frame)
  # , ylim = (-5, )
  # Reorder 'independent.variable.data.frame' by column mean (big to small)
  mns <- colMeans(independent.variable.data.frame, na.rm=TRUE)
  # order(mns, decreasing = FALSE)
  independent.variable.data.frame <-
    independent.variable.data.frame[,order(mns, decreasing = TRUE)]

  melted.data.normal <- reshape2::melt(independent.variable.data.frame)
  x.range.normal <- stats::quantile(melted.data.normal$value,probs=c(0,0.9))
  ggplot(aes(x = melted.data.normal$value,
             colour = melted.data.normal$variable),
         data = melted.data.normal) +
    xlim(x.range.normal[1]-20, x.range.normal[2]+20) +
    geom_density() + theme_bw() +
    xlab(which.count.normalization) + ylab("Frequency") +
    ggtitle("Frequency Plot (ggplot2)") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = "none")
  ggsave(paste0(path.prefix,
                paste0("RNASeq_results/", which.analysis,
                       "/images/preDE/Frequency/",
                       "Frequency_Plot_normalized_count_ggplot2.png")),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '", path.prefix, "RNASeq_results/", which.analysis,
          "/images/Frequency_Plot_normalized_count_ggplot2.png",
          "' has been created. \n\n")

  melted.data <- reshape2::melt(log2(independent.variable.data.frame+1))
  x.range <- stats::quantile(melted.data$value,probs=c(0,.99))
  ggplot(aes(x=melted.data$value,
             colour=melted.data$variable),
         data = melted.data) +
    xlim(x.range[1]-5, x.range[2]+5) +
    geom_density() + theme_bw() +
    xlab(bquote(~Log[2](.(which.count.normalization)+1))) + ylab("Frequency") +
    ggtitle("Frequency Plot (ggplot2)") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = "none")
  ggsave(paste0(path.prefix,
                paste0("RNASeq_results/", which.analysis,
                       "/images/preDE/Frequency/",
                       "Frequency_Plot_log_normalized_count_ggplot2.png")),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/", which.analysis, "/images/",
          "Frequency_Plot_log_normalized_count_ggplot2.png",
          "' has been created. \n\n")
}

# Box plot and violin plot
BoxViolinPlot <- function(which.analysis,
                          which.count.normalization,
                          path.prefix,
                          independent.variable,
                          case.group,
                          control.group) {
  # load gene name for further usage
  csv.results <- ParseResultCSV(which.analysis,
                                which.count.normalization,
                                path.prefix,
                                independent.variable,
                                case.group,
                                control.group)
  case.normalized <- csv.results$case
  control.normalized <- csv.results$control
  phenoData.result<- phenoDataWrap(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group)
  my_colors=c("#00AFBB", "#E7B800")
  color.group <- c(rep(1, phenoData.result$case.group.size),
                   rep(2, phenoData.result$case.group.size))
  independent.variable.data.frame <- cbind(case.normalized, control.normalized)
  log2.normalized.value = log2(independent.variable.data.frame+1)
  log2.normalized.value <- reshape2::melt(log2.normalized.value)
  colnames(log2.normalized.value) <- c("samples", which.count.normalization)
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                        which.analysis, "/images/preDE/Distribution/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/",
                      which.analysis, "/images/preDE/Distribution/"))
  }
  # Box plot
  message("\u25CF Plotting Box plot\n")
  ggplot(data = log2.normalized.value,
         aes(x=log2.normalized.value$samples,
             y=log2.normalized.value[which.count.normalization][[1]]),
         las = 2) +
    geom_boxplot(fill=my_colors[as.numeric(color.group)]) +
    theme_bw() +
    xlab("Samples") + ylab(bquote(~Log[2](.(which.count.normalization)+1))) +
    ggtitle("Box Plot (ggplot2)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  ggsave(paste0(path.prefix,
                paste0("RNASeq_results/", which.analysis,
                       "/images/preDE/Distribution/Box_Plot_ggplot2.png")),
         dpi = 300,
         width = 7,
         height = 7)
  # dev.off()
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/",which.analysis, "/images/preDE/",
          "Distribution/Box_Plot_ggplot2.png' has been created. \n\n")
  # Violin plot
  message("\u25CF Plotting Violin plot\n")
  ggplot(data = log2.normalized.value,
         aes(x=log2.normalized.value$samples,
             y=log2.normalized.value[which.count.normalization][[1]],
             color=log2.normalized.value$samples),
         las = 2) +
    geom_violin() +
    scale_color_manual(values=my_colors[as.numeric(color.group)]) +
    stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
    theme_bw() +
    xlab("Samples") + ylab(bquote(~Log[2](.(which.count.normalization)+1))) +
    ggtitle("Violin Plot (ggplot2)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = "none")
  ggsave(paste0(path.prefix,
                paste0("RNASeq_results/", which.analysis,
                       "/images/preDE/Distribution/Violin_Plot_ggplot2.png")),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/", which.analysis,
          "/images/preDE/Distribution/Violin_Plot_ggplot2.png",
          "' has been created. \n\n")
}

# PCA plot
PCAPlot <- function(which.analysis,
                    which.count.normalization,
                    path.prefix,
                    independent.variable,
                    case.group,
                    control.group){
  # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-
  #practical-guide/112-pca-principal-component-analysis-essentials/
  # load gene name for further usage
  message("\u25CF Plotting PCA related plot\n")
  csv.results <- ParseResultCSV(which.analysis,
                                which.count.normalization,
                                path.prefix,
                                independent.variable,
                                case.group,
                                control.group)
  case.normalized <- csv.results$case
  control.normalized <- csv.results$control
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                        which.analysis, "/images/preDE/PCA/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/",
                      which.analysis, "/images/preDE/PCA/"))
  }
  phenoData.result<- phenoDataWrap(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group)
  # The independent.variable group
  grp = factor(c(rep(case.group, phenoData.result$case.group.size),
                 rep(control.group,  phenoData.result$case.group.size)),
               levels = c(case.group, control.group))
  # color group
  color.group <- c(rep(1, phenoData.result$case.group.size),
                   rep(2, phenoData.result$case.group.size))
  independent.variable.data.frame <- cbind(case.normalized, control.normalized)
  normalized.trans <- data.frame(t(independent.variable.data.frame))
  normalized.trans$attribute <- grp
  pca <- FactoMineR::PCA(normalized.trans,
                         ncp=2,
                         quali.sup=length(normalized.trans),
                         graph = FALSE)
  eig.val <- factoextra::get_eigenvalue(pca)
  # Write out general PCA
  write.csv(data.frame(eig.val),
            file = paste0(path.prefix, "RNASeq_results/", which.analysis,
                          "/images/preDE/PCA/PCA_dimension_factoextra.csv"),
            row.names=TRUE)
  factoextra::fviz_eig(pca, addlabels = TRUE,
                       ylim = c(0, 50),
                       main = "PCA Dimensions") +
    labs(title ="PCA Dimensions Plot (factoextra)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  ggsave(paste0(path.prefix, paste0("RNASeq_results/", which.analysis,
                                    "/images/preDE/PCA/",
                                    "Dimension_PCA_Plot_factoextra.png")),
         dpi = 300,
         width = 7,
         height = 7)

  message("(\u2714) : '", path.prefix, "RNASeq_results/", which.analysis,
          "/images/preDE/PCA/Dimension_PCA_Plot_factoextra.png",
          "' has been created. \n")
  #var$coord: coordinates of variables to create a scatter plot
  #var$cos2: represents the quality of representation for variables on the
  #          factor map. It’s calculated as the squared coordinates:
  #          var.cos2 = var.coord * var.coord.
  #var$contrib: contains the contributions (in percentage) of the variables to
  #             the principal components. The contribution of a variable (var)
  #             to a given principal component is (in percentage) :
  #             (var.cos2 * 100) / (total cos2 of the component).
  #var <- get_pca_var(res.pca)
  factoextra::fviz_pca_ind(pca,
                           xlab = paste0("PC1(",
                                         round(data.frame(eig.val)$
                                                 variance.percent[1], 2), "%)"),
                           ylab = paste0("PC2(",
                                         round(data.frame(eig.val)$
                                                 variance.percent[2],2), "%)"),
                           legend.title = "Treatment variable",
                           legend.position = "top",
                           pointshape = 21, pointsize = 3.5, geom.ind = "point",
                           fill.ind = normalized.trans$attribute,
                           palette = c("#00AFBB", "#E7B800"),
                           habillage = normalized.trans$attribute,
                           addEllipses = TRUE) +
    labs(title ="PCA Plot (factoextra)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                "/images/preDE/PCA/PCA_Plot_factoextra.png"),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '", path.prefix, "RNASeq_results/", which.analysis,
          "/images/preDE/PCA/PCA_Plot_factoextra.png' has been created. \n")

  normalized.res.PCA <- FactoMineR::PCA(normalized.trans,
                                        scale.unit=TRUE, ncp=2,
                                        quali.sup=length(normalized.trans),
                                        graph = FALSE)
  groups <- c(case.group, control.group)
  PCA.data.frame <- data.frame("PC1" = normalized.res.PCA$ind$coord[,1],
                               "PC2" = normalized.res.PCA$ind$coord[,2])
  ggplot(PCA.data.frame,
         aes(x=normalized.res.PCA$ind$coord[,1] ,
             y=normalized.res.PCA$ind$coord[,2])) +
    geom_point(aes(color = factor(groups[as.numeric(color.group)],
                                  levels = c(case.group, control.group))),
               size = 3.5) +
    scale_color_manual(values = c("#00AFBB", "#E7B800")) +
    theme_bw() +
    xlab(paste0("PC1(", round(normalized.res.PCA$eig[,2][1], 2), "%)")) +
    ylab(paste0("PC2(", round(normalized.res.PCA$eig[,2][2], 2), "%)")) +
    ggtitle("PCA Plot (ggplot2)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position="top") +
    labs(color = independent.variable) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black")
  ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                "/images/preDE/PCA/PCA_Plot_ggplot2.png"),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/", which.analysis,
          "/images/preDE/PCA/PCA_Plot_ggplot2.png' has been created. \n\n")
}

#Correlation plot
CorrelationPlot <- function(which.analysis,
                            which.count.normalization,
                            path.prefix,
                            independent.variable,
                            case.group,
                            control.group){
  # load gene name for further usage
  message("\u25CF Plotting Correlation plot\n")
  csv.results <- ParseResultCSV(which.analysis,
                                which.count.normalization,
                                path.prefix,
                                independent.variable,
                                case.group,
                                control.group)
  case.normalized <- csv.results$case
  control.normalized <- csv.results$control
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                        which.analysis, "/images/preDE/Correlation/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/",
                      which.analysis, "/images/preDE/Correlation/"))
  }
  phenoData.result<- phenoDataWrap(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group)
  independent.variable.data.frame <- cbind(case.normalized, control.normalized)
  res <- round(stats::cor(independent.variable.data.frame,
                          method = c("pearson", "kendall", "spearman")), 3)
  max.value <- max(res)
  min.value <- min(res)
  # Correlation_dot_plot.png
  col <- colorRampPalette(c("#4477AA", "#FFFFFF", "#BB4444"))
  png(paste0(path.prefix, "RNASeq_results/", which.analysis,
             "/images/preDE/Correlation/Correlation_Dot_Plot_corrplot.png"),
      width=5,
      height=5,
      units="in",
      res=300)
  cex.before <- par("cex")
  par(cex = 0.5)
  corrplot::corrplot(res, col=col(200), type = "upper",
                     tl.col = "black", tl.srt = 45, addCoef.col = "black",
                     cl.cex = 1, mar=c(0,0,4,0),
                     cl.lim = c(min.value, max.value), is.corr = FALSE)
  mtext(expression(bold("Correlation Dot Plot (corrplot)")))
  par(cex = cex.before)
  dev.off()
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/", which.analysis,
          "/images/preDE/Correlation/Correlation_Dot_Plot_corrplot.png",
          "' has been created. \n")

  # Correlation_plot.png
  png(paste0(path.prefix, "RNASeq_results/", which.analysis,
             "/images/preDE/Correlation/",
             "Correlation_Bar_Plot_PerformanceAnalytics.png"),
      width=5,
      height=5,
      units="in",
      res=300)
  PerformanceAnalytics::chart.Correlation(res, histogram=TRUE, pch=19)
  dev.off()
  message("(\u2714) : '", path.prefix, "RNASeq_results/", which.analysis,
          "/images/preDE/Correlation/",
          "Correlation_Bar_Plot_PerformanceAnalytics.png' has been created. \n")

  # Correlation_heat_plot.png
  melted_res <- reshape2::melt(res)
  max.value <- max(melted_res$value)
  min.value <- min(melted_res$value)
  colnames(melted_res) <- c("Var1", "Var2", "value")
  ggplot(melted_res,
         aes(melted_res$Var1,
             melted_res$Var2,
             fill = melted_res$value))+
    xlab("Samples") + ylab("Samples") +
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", mid ="white", high = "red",
                         midpoint = (max.value + min.value)/2,
                         limit = c(min.value,max.value),
                         space = "Lab", name="Correlation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 30, vjust = 1,
                                     size = 10, hjust = 1))+
    labs(title ="Correlation Heat Plot (ggplot2)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          axis.title.x = element_text(size = 10 ),
          axis.title.y = element_text(size = 10)) +
    coord_fixed()
  # Print the heatmap
  ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                "/images/preDE/Correlation/Correlation_Heat_Plot_ggplot2.png"),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/", which.analysis,
          "/images/preDE/Correlation/Correlation_Heat_Plot_ggplot2.png",
          "' has been created. \n\n")
}

###################################################
#### Differential expression visualization ####
###################################################
# Volcano plot
VolcanoPlot <- function(which.analysis,
                        which.count.normalization,
                        path.prefix,
                        independent.variable,
                        case.group,
                        control.group,
                        condition.pval,
                        condition.log2FC) {
  # load gene name for further usage
  message("\u25CF Plotting Volcano plot\n")
  normalized_dataset <- read.csv(paste0(path.prefix, "RNASeq_results/",
                                        which.analysis, "/",
                                        strsplit(which.analysis, "_")[[1]][1],
                                        "_normalized_result.csv"))
  ## Volcano plot
  # Make a basic volcano plot
  log2FC.pval <- data.frame("log2FC" = normalized_dataset$log2FC,
                            "pval" = normalized_dataset$pval)
  down.regulated.gene <- log2FC.pval[((log2FC.pval$log2FC < -condition.log2FC) &
                                        (log2FC.pval$pval < condition.pval)),]
  up.regulated.gene <- log2FC.pval[((log2FC.pval$log2FC > condition.log2FC) &
                                      (log2FC.pval$pval < condition.pval)),]

  all.x.value <- c(down.regulated.gene$log2FC, up.regulated.gene$log2FC)
  x.range <- stats::quantile(all.x.value, probs=c(0.05,0.95))
  x.limit <- max(abs(x.range)) + 5

  all.y.value <- c(-log10(down.regulated.gene$pval),
                   -log10(up.regulated.gene$pval))
  y.range <- stats::quantile(all.y.value, probs=c(0.05,0.95))
  y.limit <- max(abs(y.range)) + 5
  ggplot(log2FC.pval,
         aes(x = log2FC.pval$log2FC, y=-log10(log2FC.pval$pval))) +
    geom_point(size = 0.8) +
    xlim(-x.limit, x.limit) + ylim(0, y.limit) +
    theme_bw() +
    xlab(bquote(~Log[2](fold~change))) + ylab(bquote(~-Log[10](p-value))) +
    ggtitle("Volcano (ggplot2)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position="top") +
    geom_point(data = down.regulated.gene,
               aes(x = down.regulated.gene$log2FC,
                   y=-log10(down.regulated.gene$pval)),
               colour = "#00cc00",
               size = 0.8) +
    geom_point(data = up.regulated.gene,
               aes(x = up.regulated.gene$log2FC,
                   y=-log10(up.regulated.gene$pval)),
               colour = "red",
               size = 0.8) +
    geom_hline(yintercept = -log10(condition.pval),
               linetype="dashed",
               color = "black") +
    geom_vline(xintercept = 1,
               linetype="dashed",
               color = "black") +
    geom_vline(xintercept = -1,
               linetype="dashed",
               color = "black")
  ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                "/images/DE/Volcano_Plot_graphics.png"),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '", path.prefix, "RNASeq_results/", which.analysis,
          "/images/DE/Volcano_Plot_ggplot2.png' has been created. \n\n")
}

# MA plot
MAPlot <- function(which.analysis,
                   which.count.normalization,
                   path.prefix,
                   independent.variable,
                   case.group,
                   control.group,
                   condition.pval) {
  # load gene name for further usage
  csv.results <- ParseResultCSV(which.analysis,
                                which.count.normalization,
                                path.prefix,
                                independent.variable,
                                case.group,
                                control.group)
  case.normalized <- csv.results$case
  case.size <- length(case.normalized)
  control.normalized <- csv.results$control
  control.size <- length(control.normalized)
  message("\u25CF Plotting MA plot\n")
  normalized_dataset <- read.csv(paste0(path.prefix, "RNASeq_results/",
                                        which.analysis, "/",
                                        strsplit(which.analysis, "_")[[1]][1],
                                        "_normalized_result.csv"))
  ## Ma plot
  ggplot(normalized_dataset,
         aes(x = log2(normalized_dataset[,1+case.size+control.size+3]),
             y = normalized_dataset$log2FC,
             colour = normalized_dataset$pval<condition.pval)) +
    xlab(bquote(~Log[2](.(which.count.normalization)))) +
    ylab(bquote(~Log[2](fold~change))) +
    theme_bw() +
    ggtitle("MA Plot (ggplot2)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10)) +
    scale_color_manual(values=c("#999999", "#FF0000")) +
    geom_point(size = 0.8) +
    geom_hline(yintercept=0, color="blue") +
    ylim(-6, 6) +
    theme(legend.position="top") +
    labs(color=paste0("p-value < ", condition.pval))
  ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                "/images/DE/MA_Plot_ggplot2.png"),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '", path.prefix, "RNASeq_results/", which.analysis,
          "/images/preDE/MA_Plot_ggplot2.png' has been created. \n\n")
}

# PCA plot
DEPCAPlot <- function(which.analysis,
                      which.count.normalization,
                      path.prefix,
                      independent.variable,
                      case.group,
                      control.group){
  # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-
  #practical-guide/112-pca-principal-component-analysis-essentials/
  # load gene name for further usage
  message("\u25CF Plotting PCA related plot\n")
  DE.csv.results <- read.csv(paste0(path.prefix, "RNASeq_results/",
                                    which.analysis, "/",
                                    strsplit(which.analysis,
                                             split = "_")[[1]][1],
                                    "_normalized_DE_result.csv"))
  phenoData.result<- phenoDataWrap(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group)
  # This is to get the result of normalized count
  DE.csv.normalized.count.only <-
    DE.csv.results[,2:(phenoData.result$case.group.size +
                         phenoData.result$case.group.size + 1)]
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                        which.analysis, "/images/DE/PCA/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/",
                      which.analysis, "/images/DE/PCA/"))
  }
  my_colors=c("#00AFBB", "#E7B800")
  # The independent.variable group
  grp = factor(c(rep(case.group, phenoData.result$case.group.size),
                 rep(control.group,  phenoData.result$case.group.size)),
               levels = c(case.group, control.group))
  # color group
  color.group <- c(rep(1, phenoData.result$case.group.size),
                   rep(2, phenoData.result$case.group.size))
  normalized.trans <- data.frame(t(DE.csv.normalized.count.only))
  normalized.trans$attribute <- grp
  pca = FactoMineR::PCA(normalized.trans, ncp=2,
                        quali.sup=length(normalized.trans),
                        graph = FALSE)
  eig.val <- factoextra::get_eigenvalue(pca)
  # Write out DE PCA
  write.csv(data.frame(eig.val),
            file = paste0(path.prefix, "RNASeq_results/", which.analysis,
                          "/images/DE/PCA/PCA_dimension_factoextra.csv"),
            row.names=TRUE)
  factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50)) +
    labs(title ="PCA Dimensions Plot (factoextra)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  ggsave(paste0(path.prefix,
                paste0("RNASeq_results/", which.analysis,
                       "/images/DE/PCA/Dimension_PCA_Plot_factoextra.png")),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/", which.analysis,
          "/images/DE/PCA/Dimension_PCA_Plot_factoextra.png",
          "' has been created. \n")
  #var$coord: coordinates of variables to create a scatter plot
  #var$cos2: represents the quality of representation for variables on the
  #          factor map. It’s calculated as the squared coordinates:
  #          var.cos2 = var.coord * var.coord.
  #var$contrib: contains the contributions (in percentage) of the variables
  #             to the principal components. The contribution of a variable
  #             (var) to a given principal component is (in percentage) :
  #             (var.cos2 * 100) / (total cos2 of the component).
  #var <- get_pca_var(res.pca)
  factoextra::fviz_pca_ind(pca,
                           xlab = paste0("PC1(",
                                         round(data.frame(eig.val)$
                                                 variance.percent[1], 2), "%)"),
                           ylab = paste0("PC2(",
                                         round(data.frame(eig.val)$
                                                 variance.percent[2],2), "%)"),
                           legend.title = "Treatment variable",
                           legend.position = "top",
                           pointshape = 21, pointsize = 3.5, geom.ind = "point",
                           fill.ind = normalized.trans$attribute,
                           palette = c("#00AFBB", "#E7B800"),
                           habillage = normalized.trans$attribute,
                           addEllipses = TRUE
  ) +
    labs(title ="PCA Plot (factoextra)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                "/images/DE/PCA/PCA_Plot_factoextra.png"),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/", which.analysis,
          "/images/DE/PCA/PCA_Plot_factoextra.png' has been created. \n")

  normalized.res.PCA <- FactoMineR::PCA(normalized.trans,
                                        scale.unit=TRUE, ncp=2,
                                        quali.sup=length(normalized.trans),
                                        graph = FALSE)
  groups <- c(case.group, control.group)
  PCA.data.frame <- data.frame("PC1" = normalized.res.PCA$ind$coord[,1],
                               "PC2" = normalized.res.PCA$ind$coord[,2])
  ggplot(PCA.data.frame, aes(x=normalized.res.PCA$ind$coord[,1] ,
                             y=normalized.res.PCA$ind$coord[,2])) +
    geom_point(aes(color = factor(groups[as.numeric(color.group)],
                                  levels = c(case.group, control.group))),
               size = 3.5) +
    scale_color_manual(values = c("#00AFBB", "#E7B800")) +
    theme_bw() +
    xlab(paste0("PC1(", round(normalized.res.PCA$eig[,2][1], 2), "%)")) +
    ylab(paste0("PC2(", round(normalized.res.PCA$eig[,2][2], 2), "%)")) +
    ggtitle("PCA Plot (ggplot2)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10)) +
    theme(legend.position="top") +
    labs(color = independent.variable) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black")
  ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                "/images/DE/PCA/PCA_Plot_ggplot2.png"),
         dpi = 300,
         width = 7,
         height = 7)
  message("(\u2714) : '",
          path.prefix, "RNASeq_results/", which.analysis,
          "/images/DE/PCA/PCA_Plot_ggplot2.png' has been created. \n\n")
}

DEHeatmap <- function(which.analysis,
                      which.count.normalization,
                      path.prefix,
                      independent.variable,
                      case.group,
                      control.group) {
  # load gene name for further usage
  message("\u25CF Plotting Differential Expressed Heatmap related plot\n")
  DE.csv.results <- read.csv(paste0(path.prefix, "RNASeq_results/",
                                    which.analysis, "/",
                                    strsplit(which.analysis,
                                             split = "_")[[1]][1],
                                    "_normalized_DE_result.csv"))
  phenoData.result<- phenoDataWrap(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group)
  ## Maybe change !!!! temp !!
  DE.csv.results <- DE.csv.results[DE.csv.results$gene.name != ".",]
  message("     \u25CF Checking found differential express transcript term.\n")
  if (nrow(DE.csv.results) == 0) {
    message("          \u25CF (\u26A0) No term were found.\n\n")
  } else {
    DE.csv.results <- DE.csv.results[order(abs(DE.csv.results$log2FC),
                                           decreasing = TRUE),]
    DE.csv.results <- DE.csv.results[!duplicated(DE.csv.results$gene.name),]

    # Decide sample name whether to show
    if (phenoData.result$case.group.size + phenoData.result$case.group.size > 16) {
      message("          \u25CF Total sample number: ",
              phenoData.result$case.group.size +
                phenoData.result$case.group.size,
              " (sample name will not be shown).\n")
      show.sample.name <- FALSE
    } else {
      message("          \u25CF Total sample number: ",
              phenoData.result$case.group.size +
                phenoData.result$case.group.size,
              " (sample name will be shown).\n")
      show.sample.name <- TRUE
    }
    # Decide gene name whether to show
    DE.csv.normalized.count.only <-
      DE.csv.results[,2:(phenoData.result$case.group.size +
                           phenoData.result$case.group.size + 1)]
    row.names(DE.csv.normalized.count.only) <- DE.csv.results$gene.name
    if (nrow(DE.csv.normalized.count.only) > 60) {
      message("          \u25CF Found ",
              nrow(DE.csv.normalized.count.only),
              " terms. More than 60 terms (gene names will not be shown).\n")
      show.gene.name <- FALSE
    } else {
      message("          \u25CF Found ", nrow(DE.csv.normalized.count.only),
              " terms (gene names will be shown).\n")
      show.gene.name <- TRUE
    }
    message("     \u25CF Calculating log2(",which.count.normalization, "+1).\n")
    log.data.frame <- log2(DE.csv.normalized.count.only+1)
    # Getting log control mean
    control.log.average <-
      rowMeans(log.data.frame[seq_len(phenoData.result$case.group.size)])
    message("     \u25CF Each log2(", which.count.normalization,
            "+1) minus average of control.\n")
    log.data.frame.minus <- log.data.frame - control.log.average
    df.new <- scale(log.data.frame.minus)
    phenoData.result<- phenoDataWrap(path.prefix,
                                     independent.variable,
                                     case.group,
                                     control.group)
    # The independent.variable group
    # Do for annotation ! (grouping in pheatmap)
    annotation_list = factor(c(rep(case.group,
                                   phenoData.result$case.group.size),
                               rep(control.group,
                                   phenoData.result$case.group.size)),
                             levels = c(case.group, control.group))
    annotation <- data.frame(Var1 = annotation_list)
    colnames(annotation) <- independent.variable
    # check out the row names of annotation
    rownames(annotation) <- colnames(df.new)

    my_colors_list <- c("#00AFBB","#E7B800")
    names(my_colors_list) <- c(case.group, control.group)
    anno_colors <- list(Var1 = my_colors_list)
    names(anno_colors) <- independent.variable

    # Check Na(list) or Infinite(numeric)
    if (any(is.na((df.new)) | is.infinite((df.new)))) {
      message("(\u26A0) : There are invalid value after ",
              "scaling DEG dataframe. Heatmap can't be drawn !\n\n")
    } else {
      redgreen <- c("blue", "white", "red")
      pal <- colorRampPalette(redgreen)(100)
      ## Not change distance , highlight case and control
      pheatmap::pheatmap(df.new, scale = "row",
                         xlab = "samples", ylab = "transcript names",
                         margins = c(10,8), col = pal,
                         main = "Heatmap Plot (pheatmap)",
                         cluster_rows = TRUE, cluster_cols = FALSE,
                         show_rownames = show.gene.name,
                         show_colnames = show.sample.name,
                         annotation_col = annotation,
                         annotation_colors = anno_colors,
                         filename = paste0(path.prefix, "RNASeq_results/",
                                           which.analysis,
                                           "/images/DE/",
                                           "Heatmap_Plot_pheatmap.png"),
                         fontsize = 7)
      message("(\u2714) : '", path.prefix, "RNASeq_results/",
              which.analysis, "/images/DE/Heatmap_Plot_pheatmap.png",
              "' has been created. \n\n")
    }
  }
}

