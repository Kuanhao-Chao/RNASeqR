#' Run ballgown analysis
#'
#' @import ballgown
#' @import genefilter
#' @importFrom dplyr arrange
#' @export
BallgownPreprocess <- function(path.prefix, gene.name, sample.pattern, experiment.type, main.variable, additional.variable) {
  # this ballgown function is only for two group
  results <- ProgressGenesFiles(path.prefix = path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print = FALSE)
  if (isTRUE(results$phenodata.file.df) && results$ballgown.dirs.number.df != 0){
    # sorting 'pheno_data'
    cat(paste0("************** Ballgown data preprocessing **************\n"))
    cat("\u25CF 1. Printing origin phenodata.csv : \n")
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    print(pheno_data)
    cat('\n')
    sample.table <- as.data.frame(table(pheno_data[main.variable]))
    if (length(row.names(sample.table)) == 2) {
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/"))
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Ballgown_object/"))
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/"))
      cat("\u25CF 2. Sorting phenodata.csv : \n")
      pheno_data.arrange <- dplyr::arrange(pheno_data, unlist(pheno_data[main.variable]))
      print(pheno_data.arrange)
      cat('\n')
      # for adding FPKM column!
      sample.names <- as.character(pheno_data.arrange$ids)
      sample.names.with.main.variable <- paste0(pheno_data.arrange$ids, ".", pheno_data.arrange[main.variable][,1])
      sample.number <- length(sample.names)
      # make ballgown object

      cat("\u25CF 3. Making ballgown object : \n")
      pkg.ballgown.data$bg_chrX <- ballgown(dataDir = paste0(path.prefix, "gene_data/ballgown"), samplePattern = sample.pattern, pData = pheno_data, meas = 'all')
      bg <- pkg.ballgown.data$bg_chrX
      save(bg, file = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Ballgown_object/ballgown.rda"))
      cat("     \u25CF writing data.frame into 'ballgown.rda' in \n")
      cat('\n')
      cat("\u25CF 4. Filtering ballgown object (variance less than 1): \n")
      pkg.ballgown.data$bg_chrX_filt <- ballgown::subset(pkg.ballgown.data$bg_chrX, cond = 'genefilter::rowVars(texpr(pkg.ballgown.data$bg_chrX)) >1', genomesubset=TRUE)
      bg_filter <- pkg.ballgown.data$bg_chrX_filt
      save(bg_filter, file = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Ballgown_object/ballgown_filter.rda"))
      cat("     \u25CF writing data.frame into 'ballgown_filter.rda' in \n")
      cat('\n')

      # differential expression
      cat("\u25CF 5. Differential Transcript expression preprocessing : \n")
      cat("     \u25CF creating 'Ballgown Transcript FPKM data.frame ......\n")
      cat(c("         \u25CF  main.variable :", main.variable, "\n"))
      cat(c("         \u25CF  experiment.type :", experiment.type, "\n"))

      results_transcripts <- stattest(pkg.ballgown.data$bg_chrX_filt, feature="transcript",covariate = main.variable, getFC=TRUE, meas="FPKM")
      write.csv(results_transcripts, paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"), row.names=FALSE)
      results_transcripts$feature <- NULL
      results_transcripts.FC <- results_transcripts$fc
      results_transcripts.log2FC <- log2(results_transcripts$fc)
      results_transcripts.pval <- results_transcripts$pval
      results_transcripts.qval <- results_transcripts$qval
      results_transcripts$fc <- NULL; results_transcripts$pval <- NULL; results_transcripts$qval <- NULL; colnames(results_transcripts)[1] <- "transcriptIDs"
      results_transcripts <- data.frame(geneNames=ballgown::geneNames(pkg.ballgown.data$bg_chrX_filt), geneIDs=ballgown::geneIDs(pkg.ballgown.data$bg_chrX_filt), transcriptNames=transcriptNames(pkg.ballgown.data$bg_chrX_filt), results_transcripts)
      # adding fpkm
      # cov : average per-base read coverage
      cat("     \u25CF merging each FPKM column and calculating average FPKM ......\n")
      fpkm <- data.frame(texpr(pkg.ballgown.data$bg_chrX_filt,meas="FPKM"))
      all.mean <- c()
      for(i in 1:length(row.names(sample.table))) {
        columns.to.mean <- c()
        current.sum <- 0
        if (i-1 == 0 ) current.sum = 0
        else {
          for(z in 1:(i-1)) {
            current.sum <- current.sum + sample.table$Freq[z]
          }
        }
        for(j in 1:sample.table$Freq[i]){
          a <- paste0("FPKM.", sample.names[current.sum+j])
          results_transcripts[[sample.names.with.main.variable[current.sum+j]]] <- fpkm[[a]]
          columns.to.mean <- append(columns.to.mean, sample.names.with.main.variable[current.sum+j])
        }
        results_transcripts[[paste0(as.character(sample.table$Var1)[i], ".mean")]] <- rowMeans(results_transcripts[columns.to.mean])
        all.mean <- append(all.mean, paste0(as.character(sample.table$Var1)[i], ".mean"))
        columns.to.mean <- c()
      }
      results_transcripts[["FPKM.all.mean"]] <- rowMeans(results_transcripts[all.mean])
      results_transcripts[["FC"]] <- results_transcripts.FC
      results_transcripts[["log2FC"]] <-results_transcripts.log2FC
      results_transcripts[["pval"]] <- results_transcripts.pval
      results_transcripts[["qval"]] <- results_transcripts.qval
      cat("     \u25CF writing data.frame into 'ballgown_FPKM_result.csv' ......\n\n")
      cat("     \u25CF writing data.frame into 'ballgown_FPKM_result.png' ......\n\n")

      write.csv(results_transcripts, paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"), row.names=FALSE)
      cat("\u25CF 6. Printing DEG dataset : \n")
      print(head(results_transcripts))
      cat("\n")
      png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.png"), width = sample.number*200 + 200, height = 40*8)
      p <- gridExtra::grid.table(head(results_transcripts))
      print(p)
      dev.off()
      cat("\n")
    } else {
      stop(paste0("(\u2718) 'experiment.type' is 'two.group'. It is only available for 2-group comparisons. However ",length(row.names(sample.table)), "-group is detected.\n\n"))
    }
  }
}

#' Frequency plot
#'
#' @importFrom rafalib mypar shist
BallgownFrequencyPlot <- function() {
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))){
    # load gene name for further usage
    cat(paste0("\u25CF Plotting  Frequency plot\n"))
    FPKM_dataset <- read.csv(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))
    # frequency plot
    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Frequency_plot.png"))
    pms <- FPKM_dataset
    rafalib::mypar(1, 1)
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    sample.table <- as.data.frame(table(pheno_data[2]))
    rafalib::shist(log2(pms[, 5]), unit = 0.1, type = "n", xlab = "log (base 2) FPKM",
                   main = "Frequency Plot", cex.main = 4, xlim = c(-5, 15))
    for(i in 1:length(row.names(sample.table))){
      current.sum <- 0
      if (i-1 == 0 ) current.sum = 0
      else {
        for(z in 1:(i-1)) {
          current.sum <- current.sum + sample.table$Freq[z]
        }
      }
      for(j in 1:sample.table$Freq[i]){
        plot.column.number <- 4+j+current.sum + i -1
        rafalib::shist(log2(pms[, plot.column.number]), unit = 0.1, col = plot.column.number, add = TRUE, lwd = 2, lty = plot.column.number)
      }
    }
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Frequency_plot.png"), "' has been created. \n\n"))
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

#' DEGTranscriptRelatedPlot
#'
BallgownTranscriptRelatedPlot <- function(){
  # draw for distribution of transcript count per gene
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))){
    # load gene name for further usage
    if (is.null(pkg.ballgown.data$bg_chrX)) {
      LoadBallgownObject()
    } else {
      cat(paste0("\u25CF Plotting Transcript related plot\n"))
      if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Transcript_Related/"))){
        dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Transcript_Related/"))
      }
      transcript_gene_table <- indexes(pkg.ballgown.data$bg_chrX)$t2g
      counts=table(transcript_gene_table[,"g_id"])
      c_one = length(which(counts == 1))
      c_more_than_one = length(which(counts > 1))
      c_max = max(counts)
      png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Transcript_Related/Distribution_transcript_count_per_gene_plot.png"))
      hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of Transcript Count per Gene")
      legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
      legend("topright", legend_text, lty=NULL)
      dev.off()
      cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Transcript_Related/Distribution_transcript_count_per_gene_plot.png"), "' has been created. \n"))

      # draw the distribution of transcript length
      full_table <- texpr(pkg.ballgown.data$bg_chrX, 'all')
      t.mini.length = min(full_table$length[full_table$length > 0])
      t.max.length = max(full_table$length[full_table$length > 0])
      png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Transcript_Related/Distribution_transcript_length_plot.png"))
      hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of Transcript Lengths", col="steelblue")
      legend_text = c(paste("Minimum transcript length =", t.mini.length), paste("Maximum transcript length =", t.max.length))
      legend("topright", legend_text, lty=NULL)
      dev.off()
      cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Transcript_Related/Distribution_transcript_length_plot.png"), "' has been created. \n\n"))
    }
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

#'
BallgownBoxViolinPlot <- function() {
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))){
    # load gene name for further usage
    if (is.null(pkg.ballgown.data$bg_chrX)) {
      LoadBallgownObject()
    } else {
      cat(paste0("\u25CF Plotting FPKM Box plot\n"))
      pheno.data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
      # frequency plot
      # tropical <- c('darkorange', 'dodgerblue')
      # palette(tropical)
      my_colors=c(rgb(255, 47, 35,maxColorValue = 255),
                  rgb(50, 147, 255,maxColorValue = 255))
      fpkm = data.frame(texpr(pkg.ballgown.data$bg_chrX,meas="FPKM"))
      fpkm = log2(fpkm+1)
      fpkm <- reshape2::melt(fpkm)
      colnames(fpkm) <- c("samples", "FPKM")
      png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Box_plot.png"))
      p1 <- ggplot(data = fpkm,  aes(x=samples, y=FPKM), las = 2) + geom_boxplot(fill=my_colors[as.numeric(pheno.data[,2])]) +
        xlab("Samples") + ylab("Log2(FPKM+1)") + ggtitle("Transcript FPKM Box Plot") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
        theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
      print(p1)
      dev.off()
      cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Box_plot.png"), "' has been created. \n\n"))

      cat(paste0("\u25CF Plotting FPKM Violin plot\n"))
      png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Violin_plot.png"))
      p2 <- ggplot(data = fpkm,  aes(x=samples, y=FPKM, color=samples), las = 2) + geom_violin() +
        scale_color_manual(values=my_colors[as.numeric(pheno.data[,2])]) + stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
        xlab("Samples") + ylab("Log2(FPKM+1)") + ggtitle("Transcript FPKM Violin Plot") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
        theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
      print(p2)
      dev.off()
      cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Violin_plot.png"), "' has been created. \n\n"))
    }
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

#' DEGPCAPlot
#'
#' @importFrom FactoMineR PCA
#' @importFrom factoextra get_eigenvalue fviz_eig fviz_pca_ind
BallgownPCAPlot <- function(){
  # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))){
    # load gene name for further usage
    if (is.null(pkg.ballgown.data$bg_chrX)) {
      LoadBallgownObject()
    } else {
      cat(paste0("\u25CF Plotting PCA related plot\n"))
      if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/PCA/"))){
        dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/PCA/"))
      }
      fpkm <- data.frame(texpr(pkg.ballgown.data$bg_chrX,meas="FPKM"))
      fpkm.trans <- data.frame(t(fpkm))
      fpkm.trans.row.names <- row.names(fpkm.trans)
      fpkm.trans.row.names.clean <- gsub("FPKM.", "", fpkm.trans.row.names)
      row.names(fpkm.trans) <- fpkm.trans.row.names.clean

      fpkm.trans.sort <- fpkm.trans[ order(row.names(fpkm.trans)), ]
      pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
      pheno_data.sort <- pheno_data[order(pheno_data$id),]
      fpkm.trans.sort$attribute <- pheno_data.sort[,2]
      # scale.unit = TRUE ==> the data are scaled to unit variance before the analysis.
      # fpkm.pca <- PCA(fpkm.trans.sort[,-ncol(fpkm.trans.sort)], scale.unit = TRUE, ncp = 2, graph = FALSE)
      fpkm.pca = FactoMineR::PCA(fpkm.trans.sort, ncp=2, quali.sup=length(fpkm.trans.sort), graph = FALSE)
      eig.val <- factoextra::get_eigenvalue(fpkm.pca)
      png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/PCA/Dimension_pca_plot.png"))
      p1 <- factoextra::fviz_eig(fpkm.pca, addlabels = TRUE, ylim = c(0, 50), main = "PCA Dimensions") +
        labs(title ="PCA Dimensions") +
        theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
      print(p1)
      dev.off()
      cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/PCA/Dimension_pca_plot.png"), "' has been created. \n"))
      #var$coord: coordinates of variables to create a scatter plot
      #var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
      #var$contrib: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
      #var <- get_pca_var(res.pca)
      png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/PCA/PCA_plot_factoextra.png"))
      p2 <- factoextra::fviz_pca_ind(fpkm.pca,
                                     xlab = paste0("PC1(", round(data.frame(eig.val)$variance.percent[1], 2), "%)"), ylab = paste0("PC2(", round(data.frame(eig.val)$variance.percent[2],2), "%)"),
                                     legend.title = "Treatment variable", legend.position = "top",
                                     pointshape = 21,
                                     pointsize = 2.5,
                                     geom.ind = "point", # show points only (nbut not "text")
                                     habillage = fpkm.trans.sort$attribute,
                                     fill.ind = fpkm.trans.sort$attribute,
                                     col.ind = fpkm.trans.sort$attribute, # color by groups
                                     addEllipses=TRUE
                                     #palette = c("#00AFBB", "#E7B800"),
                                     #                 addEllipses = TRUE, # Concentration ellipses
                                     ) +
        labs(title ="Principal Component Analysis") +
        theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
      print(p2)
      dev.off()
      cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/PCA/PCA_plot_factoextra.png"), "' has been created. \n"))

      png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/PCA/PCA_plot_self.png"))
      # fpkm.trans.sort$attribute <- factor(fpkm.trans.sort$attribute)
      #length(fpkm.trans.sort)
      FPKM.res.PCA = FactoMineR::PCA(fpkm.trans.sort, scale.unit=TRUE, ncp=2, quali.sup=length(fpkm.trans.sort), graph = FALSE)

      my_colors=c(rgb(255, 47, 35,maxColorValue = 255),
                  rgb(50, 147, 255,maxColorValue = 255))

      #par(mfrow=c(1,1))
      # fpkm.trans.sort <- fpkm.trans[ order(row.names(fpkm.trans)), ]
      #FPKM.res.PCA$ind$coord <- FPKM.res.PCA$ind$coord[ order(row.names(FPKM.res.PCA$ind$coord)), ]
      plot(FPKM.res.PCA$ind$coord[,1] , FPKM.res.PCA$ind$coord[,2], xlab=paste0("PC1(", round(FPKM.res.PCA$eig[,2][1], 2), "%)") , ylab=paste0("PC2(", round(FPKM.res.PCA$eig[,2][2], 2), "%)") , pch=20 , cex=3 ,
           col=my_colors[as.numeric(FPKM.res.PCA$call$quali.sup$quali.sup[,1])] )
      #my_colors[as.numeric(res.PCA$call$quali.sup$quali.sup[,1])]
      abline(h=0 , v=0, lty= 2)
      par(xpd=TRUE)
      legend("bottomright",inset=c(0,1), horiz=TRUE, bty="n", legend=levels(FPKM.res.PCA$call$quali.sup$quali.sup[,1] ) , col=my_colors, pch=20 )
      dev.off()
      cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/PCA/PCA_plot.png"), "' has been created. \n\n"))
    }
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

#' Plot correlation plo
#'
#' @importFrom corrplot corrplot
#' @importFrom PerformanceAnalytics chart.Correlation
#' @importFrom reshape2 melt
BallgownCorrelationPlot <- function(){
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))){
    # load gene name for further usage
    cat(paste0("\u25CF Plotting Correlation plot\n"))
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Correlation/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Correlation/"))
    }
    FPKM_dataset <- read.csv(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))
    pheno.data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    pheno.data.table <- as.data.frame(table(pheno.data[2]))
    count <- 0
    select.column <- c()
    for(i in 1:length(row.names(pheno.data.table))){
      for(j in 1:pheno.data.table$Freq[i]){
        a <- j + 4 + count + (i-1)
        select.column <- c(select.column, a)
      }
      count <- count + pheno.data.table$Freq[i]
    }
    res <- round(cor(FPKM_dataset[select.column], method = c("pearson", "kendall", "spearman")), 3)
    # Correlation_dot_plot.png
    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Correlation/Correlation_dot_plot.png"))
    cex.before <- par("cex")
    par(cex = 0.7)
    corrplot::corrplot(res, col=col(200), type = "upper", tl.col = "black", tl.srt = 45, addCoef.col = "black", cl.cex = 1/par("cex"), mar=c(0,0,1,0))
    mtext("Correlation Dot Plot", at=7, line=-0.5, cex=1.5)
    par(cex = cex.before)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Correlation/Correlation_dot_plot.png"), "' has been created. \n"))

    # Correlation_plot.png
    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Correlation/Correlation_plot.png"))
    p2 <- PerformanceAnalytics::chart.Correlation(res, histogram=TRUE, pch=19)
    print(p2)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Correlation/Correlation_plot.png"), "' has been created. \n"))

    # Correlation_heat_plot.png
    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Correlation/Correlation_heat_plot.png"))
    melted_res <- reshape2::melt(res)
    colnames(melted_res) <- c("Var1", "Var2", "value")
    ggheatmap <- ggplot(melted_res, aes(Var1, Var2, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue",mid ="white"  ,high = "red",
                           midpoint = 0, limit = c(-1,1), space = "Lab",
                           name="Correlation") +
      theme_minimal()+ # minimal theme
      theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 10, hjust = 1))+
      labs(title ="Correlation Heat Plot") +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10 ), axis.title.y = element_text(size = 10)) +
      coord_fixed()
    # Print the heatmap
    print(ggheatmap)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Correlation/Correlation_heat_plot.png"), "' has been created. \n\n"))
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}


#' inner function : DEG volcanplot
BallgownVolcanoPlot <- function(ballgown.pval, ballgown.log2FC) {
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))){
    # load gene name for further usage
    cat(paste0("\u25CF Plotting Volcano plot\n"))
    FPKM_dataset <- read.csv(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))
    ## Volcano plot
    # Make a basic volcano plot
    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Volcano_plot.png"))
    par(mar=c(5,7,5,5), cex=0.6, cex.main=2, cex.axis=1.5, cex.lab=1.5)
    topT <- as.data.frame(FPKM_dataset)
    with(topT, plot(log2FC, -log10(pval), pch=20, main="Volcano Plot", xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value), xlim=c(-15,15), ylim = c(0,12)))
    # user4 input pvalue log2FC
    # pval to qvalue
    with(subset(topT, pval<ballgown.pval & abs(log2FC)>=ballgown.log2FC), points(log2FC, -log10(pval), pch=20, cex=1, col="red"))
    with(subset(topT, pval<ballgown.pval & log2FC<= -1*ballgown.log2FC), points(log2FC, -log10(pval), pch=20, cex=1, col="green"))
    # hight = -log10(pavl) = height
    abline(v=c(-1*ballgown.log2FC,ballgown.log2FC), h=-1*log10(ballgown.pval), col="black", lty='dashed')
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/Volcano_plot.png"), "' has been created. \n\n"))
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

#' DEGMAPlot
#'
#' @import ggplot2
#' @export
BallgownMAPlot <- function(ballgown.qval) {
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))){
    # load gene name for further usage
    cat(paste0("************** Plotting MA plot **************\n"))
    FPKM_dataset <- read.csv(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"))
    ## Ma plot
    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/MA_plot.png"))
    p <- ggplot(FPKM_dataset, aes(x = log2(FPKM.all.mean), y = log2FC, colour = qval<ballgown.qval)) +
      xlab("Log2(FPKM.all.mean)") +
      ylab("Log2FC") +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
      labs(title = "MA Plot") +
      scale_color_manual(values=c("#999999", "#FF0000")) +
      geom_point() +
      geom_hline(yintercept=0, color="blue") +
      ylim(-6, 6)
    print(p)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/ballgown_FPKM_result.csv"), "' has been created. \n\n"))
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

#'
#' @export
BallgownPlotAll <- function(ballgown.log2FC, ballgown.pval, ballgown.qval) {
  cat(paste0("************** Ballgown result visualization **************\n"))
  BallgownFrequencyPlot()
  BallgownTranscriptRelatedPlot()
  BallgownBoxViolinPlot()
  BallgownPCAPlot()
  BallgownCorrelationPlot()
  BallgownVolcanoPlot(ballgown.pval, ballgown.log2FC)
  BallgownMAPlot(ballgown.qval)
}


#' Check ballgown object
#' @export
CheckBallgownObject <- function() {
  print(pkg.ballgown.data$bg_chrX)
  print(pkg.ballgown.data$bg_chrX_filt)
}


#' load ballgown object
LoadBallgownObject <- function() {
  if(isTRUE(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Ballgown_object/ballgown.rda")))) {
    load(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Ballgown_object/ballgown.rda"))
    pkg.ballgown.data$bg_chrX <- bg
  } else {
    stop(paste0("(\u2718) '", paste0(path.prefix, "gene_data/ballgown/ballgown.rda"), "' haven't created yet. Please run \"BallgownPreprocess(gene.name, sample.pattern, main.variable)\" first.\n\n"))
  }
}
