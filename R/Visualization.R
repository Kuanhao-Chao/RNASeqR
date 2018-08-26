###################################################
#### Pre-differential expression visualization ####
###################################################
# Frequency Plot
FrequencyPlot <- function(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group) {
  cat(paste0("\u25CF Plotting  Frequency plot\n"))
  csv.results <- ParseResultCSV(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group)
  control.normalized <- csv.results$control
  experiment.normalized <- csv.results$experiment
  independent.variable.data.frame <- cbind(control.normalized, experiment.normalized)
  png(paste0(path.prefix, paste0("RNASeq_results/", which.analysis, "/images/preDE/Frequency_Plot_rafalib.png")))
  rafalib::mypar(1, 1)
  sample.size <- length(independent.variable.data.frame)
  rafalib::shist(log2(independent.variable.data.frame[, 1]), unit = 0.5, type = "n", xlab = paste0("log2(", which.count.normalization, ")"),
                 main = "Frequency Plot (rafalib)", cex.main = 4, xlim = c(-10, 20))
  for (i in seq_len(sample.size)){
    rafalib::shist(log2(independent.variable.data.frame[, i]), unit = 0.5, col = i, add = TRUE, lwd = 2, lty = i, xlim = c(-10, 20))
  }
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/Frequency_Plot_rafalib.png"), "' has been created. \n\n"))
}

# Box plot and violin plot
BoxViolinPlot <- function(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group) {
  # load gene name for further usage
  csv.results <- ParseResultCSV(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group)
  control.normalized <- csv.results$control
  experiment.normalized <- csv.results$experiment
  pre.pheno_data <- RawCountPreData(path.prefix, independent.variable, control.group, experiment.group)
  my_colors=c(rgb(50, 147, 255,maxColorValue = 255),
              rgb(255, 47, 35,maxColorValue = 255))
  color.group <- c(rep(1, pre.pheno_data$control.group.size), rep(2, pre.pheno_data$experiment.group.size))
  independent.variable.data.frame <- cbind(control.normalized, experiment.normalized)
  log2.normalized.value = log2(independent.variable.data.frame+1)
  log2.normalized.value <- reshape2::melt(log2.normalized.value)
  colnames(log2.normalized.value) <- c("samples", which.count.normalization)
  # Box plot
  cat(paste0("\u25CF Plotting Box plot\n"))
  png(paste0(path.prefix, paste0("RNASeq_results/", which.analysis, "/images/preDE/Box_Plot_ggplot2.png")))
  p1 <- ggplot(data = log2.normalized.value,  aes(x=log2.normalized.value$samples, y=log2.normalized.value[which.count.normalization][[1]]), las = 2) + geom_boxplot(fill=my_colors[as.numeric(color.group)]) +
    xlab("Samples") + ylab(paste0("Log2(", which.count.normalization, "+1)")) + ggtitle("Box Plot (ggplot2)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
    theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
  print(p1)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/Box_Plot_ggplot2.png"), "' has been created. \n\n"))
  # Violin plot
  cat(paste0("\u25CF Plotting Violin plot\n"))
  png(paste0(path.prefix, paste0("RNASeq_results/", which.analysis, "/images/preDE/Violin_Plot_ggplot2.png")))
  p2 <- ggplot(data = log2.normalized.value,  aes(x=log2.normalized.value$samples, y=log2.normalized.value[which.count.normalization][[1]], color=log2.normalized.value$samples), las = 2) + geom_violin() +
    scale_color_manual(values=my_colors[as.numeric(color.group)]) + stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
    xlab("Samples") + ylab(paste0("Log2(", which.count.normalization, "+1)")) + ggtitle("Violin Plot (ggplot2)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
    theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
    theme(legend.position = "none")
  print(p2)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/Violin_Plot_ggplot2.png"), "' has been created. \n\n"))
}

# PCA plot
PCAPlot <- function(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group){
  # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
  # load gene name for further usage
  cat(paste0("\u25CF Plotting PCA related plot\n"))
  csv.results <- ParseResultCSV(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group)
  control.normalized <- csv.results$control
  experiment.normalized <- csv.results$experiment
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/PCA/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/PCA/"))
  }
  pre.pheno_data <- RawCountPreData(path.prefix, independent.variable, control.group, experiment.group)
  # The independent.variable group
  grp = factor(c(rep(control.group, pre.pheno_data$control.group.size), rep(experiment.group,  pre.pheno_data$experiment.group.size)))
  # color group
  color.group <- c(rep(1, pre.pheno_data$control.group.size), rep(2, pre.pheno_data$experiment.group.size))
  independent.variable.data.frame <- cbind(control.normalized, experiment.normalized)
  normalized.trans <- data.frame(t(independent.variable.data.frame))
  normalized.trans$attribute <- grp
  pca = FactoMineR::PCA(normalized.trans, ncp=2, quali.sup=length(normalized.trans), graph = FALSE)
  eig.val <- factoextra::get_eigenvalue(pca)
  png(paste0(path.prefix, paste0("RNASeq_results/", which.analysis, "/images/preDE/PCA/Dimension_PCA_Plot_factoextra.png")))
  p1 <- factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50), main = "PCA Dimensions") +
    labs(title ="PCA Dimensions Plot (factoextra)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
  print(p1)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/PCA/Dimension_PCA_Plot_factoextra.png"), "' has been created. \n"))
  #var$coord: coordinates of variables to create a scatter plot
  #var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
  #var$contrib: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
  #var <- get_pca_var(res.pca)
  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/PCA/PCA_Plot_factoextra.png"))
  p2 <- factoextra::fviz_pca_ind(pca, xlab = paste0("PC1(", round(data.frame(eig.val)$variance.percent[1], 2), "%)"),
                                 ylab = paste0("PC2(", round(data.frame(eig.val)$variance.percent[2],2), "%)"),
                                 legend.title = "Treatment variable", legend.position = "top",
                                 pointshape = 21, pointsize = 2.5, geom.ind = "point", # show points only (nbut not "text")
                                 habillage = normalized.trans$attribute, fill.ind = normalized.trans$attribute,
                                 col.ind = normalized.trans$attribute, # color by groups
                                 addEllipses=TRUE
                                 ) +
    labs(title ="PCA Plot (factoextra)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
  print(p2)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/PCA/PCA_Plot_factoextra.png"), "' has been created. \n"))

  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/PCA/PCA_Plot_graphics.png"))
  normalized.res.PCA = FactoMineR::PCA(normalized.trans, scale.unit=TRUE, ncp=2, quali.sup=length(normalized.trans), graph = FALSE)
  my_colors=c(rgb(50, 147, 255,maxColorValue = 255),
              rgb(255, 47, 35,maxColorValue = 255))
  plot(normalized.res.PCA$ind$coord[,1] , normalized.res.PCA$ind$coord[,2], main = "PCA Plot (graphics)",
       xlab=paste0("PC1(", round(normalized.res.PCA$eig[,2][1], 2), "%)") ,
       ylab=paste0("PC2(", round(normalized.res.PCA$eig[,2][2], 2), "%)") ,
       pch=20 , cex=3 , col=my_colors[as.numeric(color.group)] )
  #my_colors[as.numeric(res.PCA$call$quali.sup$quali.sup[,1])]
  abline(h=0 , v=0, lty= 2)
  par(xpd=TRUE)
  legend("bottomright",inset=c(0,1), horiz=TRUE, bty="n", legend=c(control.group, experiment.group) , col=my_colors, pch=20 )
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/PCA/PCA_Plot_graphics.png"), "' has been created. \n\n"))
}

#Correlation plot
CorrelationPlot <- function(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group){
  # load gene name for further usage
  cat(paste0("\u25CF Plotting Correlation plot\n"))
  csv.results <- ParseResultCSV(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group)
  control.normalized <- csv.results$control
  experiment.normalized <- csv.results$experiment
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/Correlation/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/Correlation/"))
  }
  pre.pheno_data <- RawCountPreData(path.prefix, independent.variable, control.group, experiment.group)
  independent.variable.data.frame <- cbind(control.normalized, experiment.normalized)
  res <- round(cor(independent.variable.data.frame, method = c("pearson", "kendall", "spearman")), 3)
  # Correlation_dot_plot.png
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/Correlation/Correlation_Dot_Plot_corrplot.png"))
  cex.before <- par("cex")
  par(cex = 0.7)
  corrplot::corrplot(res, col=col(200), type = "upper", tl.col = "black", tl.srt = 45, addCoef.col = "black", cl.cex = 1/par("cex"), mar=c(0,0,1,0))
  mtext("Correlation Dot Plot (corrplot)", at=7, line=-0.5, cex=1.5)
  par(cex = cex.before)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/Correlation/Correlation_Dot_Plot_corrplot.png"), "' has been created. \n"))

  # Correlation_plot.png
  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/Correlation/Correlation_Bar_Plot_PerformanceAnalytics.png"))
  p2 <- PerformanceAnalytics::chart.Correlation(res, histogram=TRUE, pch=19)
  title(main = "Correlation Bar Plot (PerformanceAnalytics)")
  print(p2)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/Correlation/Correlation_Bar_Plot_PerformanceAnalytics.png"), "' has been created. \n"))

  # Correlation_heat_plot.png
  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/Correlation/Correlation_Heat_Plot_ggplot2.png"))
  melted_res <- reshape2::melt(res)
  colnames(melted_res) <- c("Var1", "Var2", "value")
  ggheatmap <- ggplot(melted_res, aes(melted_res$Var1, melted_res$Var2, fill = melted_res$value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue",mid ="white"  ,high = "red",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Correlation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 10, hjust = 1))+
    labs(title ="Correlation Heat Plot (ggplot2)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10 ), axis.title.y = element_text(size = 10)) +
    coord_fixed()
  # Print the heatmap
  print(ggheatmap)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/Correlation/Correlation_Heat_Plot_ggplot2.png"), "' has been created. \n\n"))
}

###################################################
#### Differential expression visualization ####
###################################################
# Volcano plot
VolcanoPlot <- function(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group, condition.pval, condition.log2FC) {
  # load gene name for further usage
  cat(paste0("\u25CF Plotting Volcano plot\n"))
  normalized_dataset <- read.csv(paste0(path.prefix, "RNASeq_results/", which.analysis, "/", strsplit(which.analysis, "_")[[1]][1], "_normalized_result.csv"))
  ## Volcano plot
  # Make a basic volcano plot
  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/Volcano_Plot_graphics.png"))
  par(mar=c(5,7,5,5), cex=0.6, cex.main=2, cex.axis=1.5, cex.lab=1.5)
  topT <- as.data.frame(normalized_dataset)
  with(topT, plot(topT$log2FC, -log10(topT$pval), pch=20, main="Volcano Plot (graphics)", xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value), xlim=c(-15,15), ylim = c(0,12)))
  subset.result.red <- subset(topT, topT$pval<condition.pval & topT$log2FC>=condition.log2FC)
  with(subset.result.red, points(subset.result.red$log2FC, -log10(subset.result.red$pval), pch=20, cex=1, col="red"))
  subset.result.green <- subset(topT, topT$pval<condition.pval & topT$log2FC<=-1*condition.log2FC)
  with(subset.result.green, points(subset.result.green$log2FC, -log10(subset.result.green$pval), pch=20, cex=1, col="green"))
  # hight = -log10(pavl) = height
  abline(v=c(-1*condition.log2FC,condition.log2FC), h=-1*log10(condition.pval), col="black", lty='dashed')
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/Volcano_Plot_graphics.png"), "' has been created. \n\n"))
}

# MA plot
MAPlot <- function(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group, condition.pval) {
  # load gene name for further usage
  csv.results <- ParseResultCSV(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group)
  control.normalized <- csv.results$control
  experiment.normalized <- csv.results$experiment
  cat(paste0("\u25CF Plotting MA plot\n"))
  normalized_dataset <- read.csv(paste0(path.prefix, "RNASeq_results/", which.analysis, "/", strsplit(which.analysis, "_")[[1]][1], "_normalized_result.csv"))
  ## Ma plot
  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/MA_Plot_ggplot2.png"))
  p <- ggplot(normalized_dataset, aes(x = log2(normalized_dataset[paste0(control.group, ".", experiment.group, ".average")][[1]]), y = normalized_dataset$log2FC, colour = normalized_dataset$pval<condition.pval)) +
    xlab(paste0("Log2(", which.count.normalization, ".mean)")) +
    ylab("Log2FC") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
    labs(title = "MA Plot (ggplot2)") +
    scale_color_manual(values=c("#999999", "#FF0000")) +
    geom_point(size = 0.8) +
    geom_hline(yintercept=0, color="blue") +
    ylim(-6, 6) +
    theme(legend.position="top") +
    labs(color=paste0("p-value < ", condition.pval))
  print(p)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/preDE/MA_Plot_ggplot2.png"), "' has been created. \n\n"))
}

# PCA plot
DEPCAPlot <- function(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group){
  # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
  # load gene name for further usage
  cat(paste0("\u25CF Plotting PCA related plot\n"))
  DE.csv.results <- read.csv(paste0(path.prefix, "RNASeq_results/", which.analysis, "/", strsplit(which.analysis, split = "_")[[1]][1], "_normalized_DE_result.csv"))
  pre.de.pheno.data <- RawCountPreData(path.prefix, independent.variable, control.group, experiment.group)
  # This is to get the result of normalized counts
  DE.csv.normalized.counts.only <- DE.csv.results[,2:(pre.de.pheno.data$control.group.size + pre.de.pheno.data$experiment.group.size + 1)]
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/PCA/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/PCA/"))
  }
  # The independent.variable group
  grp = factor(c(rep(control.group, pre.de.pheno.data$control.group.size), rep(experiment.group,  pre.de.pheno.data$experiment.group.size)))
  # color group
  color.group <- c(rep(1, pre.de.pheno.data$control.group.size), rep(2, pre.de.pheno.data$experiment.group.size))
  normalized.trans <- data.frame(t(DE.csv.normalized.counts.only))
  normalized.trans$attribute <- grp
  pca = FactoMineR::PCA(normalized.trans, ncp=2, quali.sup=length(normalized.trans), graph = FALSE)
  eig.val <- factoextra::get_eigenvalue(pca)
  png(paste0(path.prefix, paste0("RNASeq_results/", which.analysis, "/images/DE/PCA/Dimension_PCA_Plot_factoextra.png")))
  p1 <- factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50), main = "PCA Dimensions") +
    labs(title ="PCA Dimensions Plot (factoextra)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
  print(p1)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/PCA/Dimension_PCA_Plot_factoextra.png"), "' has been created. \n"))
  #var$coord: coordinates of variables to create a scatter plot
  #var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
  #var$contrib: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
  #var <- get_pca_var(res.pca)
  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/PCA/PCA_Plot_factoextra.png"))
  p2 <- factoextra::fviz_pca_ind(pca, xlab = paste0("PC1(", round(data.frame(eig.val)$variance.percent[1], 2), "%)"),
                                 ylab = paste0("PC2(", round(data.frame(eig.val)$variance.percent[2],2), "%)"),
                                 legend.title = "Treatment variable", legend.position = "top",
                                 pointshape = 21, pointsize = 2.5, geom.ind = "point", # show points only (nbut not "text")
                                 habillage = normalized.trans$attribute, fill.ind = normalized.trans$attribute,
                                 col.ind = normalized.trans$attribute, # color by groups
                                 addEllipses=TRUE
  ) +
    labs(title ="PCA Plot (factoextra)") +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
  print(p2)
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/PCA/PCA_Plot_factoextra.png"), "' has been created. \n"))

  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/PCA/PCA_Plot_graphics.png"))
  normalized.res.PCA = FactoMineR::PCA(normalized.trans, scale.unit=TRUE, ncp=2, quali.sup=length(normalized.trans), graph = FALSE)
  my_colors=c(rgb(50, 147, 255,maxColorValue = 255),
              rgb(255, 47, 35,maxColorValue = 255))
  plot(normalized.res.PCA$ind$coord[,1] , normalized.res.PCA$ind$coord[,2], main = "PCA Plot (graphics)",
       xlab=paste0("PC1(", round(normalized.res.PCA$eig[,2][1], 2), "%)") ,
       ylab=paste0("PC2(", round(normalized.res.PCA$eig[,2][2], 2), "%)") ,
       pch=20 , cex=3 , col=my_colors[as.numeric(color.group)] )
  #my_colors[as.numeric(res.PCA$call$quali.sup$quali.sup[,1])]
  abline(h=0 , v=0, lty= 2)
  par(xpd=TRUE)
  legend("bottomright",inset=c(0,1), horiz=TRUE, bty="n", legend=c(control.group, experiment.group) , col=my_colors, pch=20 )
  dev.off()
  cat(paste0("(\u2714) : '", paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/PCA/PCA_Plot_graphics.png"), "' has been created. \n\n"))
}

DEHeatmap <- function(which.analysis, which.count.normalization, path.prefix, independent.variable, control.group, experiment.group) {
  # load gene name for further usage
  cat(paste0("\u25CF Plotting Differential Expressed Heatmap related plot\n"))
  DE.csv.results <- read.csv(paste0(path.prefix, "RNASeq_results/", which.analysis, "/", strsplit(which.analysis, split = "_")[[1]][1], "_normalized_DE_result.csv"))
  pre.de.pheno.data <- RawCountPreData(path.prefix, independent.variable, control.group, experiment.group)
  DE.csv.results <- DE.csv.results[DE.csv.results$gene.name != ".",]
  DE.csv.normalized.counts.only <- DE.csv.results[,2:(pre.de.pheno.data$control.group.size + pre.de.pheno.data$experiment.group.size + 1)]
  row.names(DE.csv.normalized.counts.only) <- DE.csv.results$gene.name
  cat(paste0("     \u25CF Checking found differential express transcript term.\n"))
  if (nrow(DE.csv.normalized.counts.only) == 0) {
    cat(paste0("          \u25CF (\u26A0) No term were found.\n"))
  } else {
    if (length(row.names(DE.csv.normalized.counts.only)) > 50) {
      cat(paste0("          \u25CF Found ", length(row.names(DE.csv.normalized.counts.only)), " terms. More than 50 terms (Only plot top 50 smallest p value).\n"))
      DE.csv.normalized.counts.only <- DE.csv.normalized.counts.only[seq_len(50),]
    } else {
      cat(paste0("          \u25CF Found ", length(row.names(DE.csv.normalized.counts.only)), " terms.\n"))
    }
  }
  cat(paste0("     \u25CF Calculating log2(", which.count.normalization, "+1).\n"))
  log.data.frame <- log2(DE.csv.normalized.counts.only+1)
  # Getting log control mean
  control.log.average <- rowMeans(log.data.frame[seq_len(pre.de.pheno.data$control.group.size)])
  cat(paste0("     \u25CF Each log2(", which.count.normalization, "+1) minus average of control.\n"))
  log.data.frame.minus <- log.data.frame - control.log.average
  df.new <- scale(log.data.frame.minus)
  png(paste0(path.prefix, "RNASeq_results/", which.analysis, "/images/DE/Heatmap_Plot_pheatmap.png"), width = 1000, height = 1000)
  redgreen <- c("blue", "white", "red")
  pal <- colorRampPalette(redgreen)(100)
  pheatmap::pheatmap(df.new, scale = "row", xlab = "samples", ylab = "transcript names",cexRow=1, cexCol = 1, margins = c(10,8), col = pal, main = "Heatmap Plot (pheatmap)")
  # theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
  dev.off()
  cat(paste0("(\u2714) : '", path.prefix, "RNASeq_results/", which.analysis, "/images/DE/Heatmap_Plot_pheatmap.png"), "' has been created. \n")
}

