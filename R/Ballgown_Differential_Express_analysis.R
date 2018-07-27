find_ENTREZID_ID_DE <- function(symbol.id, gene.df.DE) {
  ENTREZID_IDs.DE <- c(ENTREZID_IDs.DE, gene.df.DE[gene.df.DE$SYMBOL == symbol.id, ]["ENTREZID"][[1]][1])
  return(ENTREZID_IDs.DE)
}

find_ENTREZID_ID_Univ <- function(symbol.id, gene.df.Univ) {
  ENTREZID_IDs_Univ <- c(ENTREZID_IDs_Univ, gene.df.Univ[gene.df.Univ$SYMBOL == symbol.id, ]["ENTREZID"][[1]][1])
  return(ENTREZID_IDs_Univ)
}

#' BallgownPCAPlot
DEBallgownPCAPlot <- function(path.prefix, independent.variable){
  # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv"))){
    # load gene name for further usage
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/"))
    }
    cat(paste0("\u25CF Plotting Differential Expressed PCA related plot\n"))
    file.path <- paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv")
    data.read.csv <- read.csv(file = file.path)
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    sample.table <- as.data.frame(table(pheno_data[independent.variable]))
    new.data.frame.index <- c()
    for(i in 1:length(row.names(sample.table))){
      current.sum <- 0
      if (i-1 == 0 ) current.sum = 0
      else {
        for(z in 1:(i-1)) {
          current.sum <- current.sum + sample.table$Freq[z]
        }
      }
      for(j in 1:sample.table$Freq[i]){
        column.number <- 4+j+current.sum + i -1
        vector <- c(vector, column.number)
        new.data.frame.index <- c(new.data.frame.index, column.number)
      }
    }
    sub.data.frame <- data.read.csv[new.data.frame.index]
    fpkm.trans <- data.frame(t(sub.data.frame))
    fpkm.trans.split <- strsplit(row.names(fpkm.trans), split = "[.]")
    fpkm.trans.independent.variable <- c()
    for( i in 1:length(fpkm.trans.split)){
      fpkm.trans.independent.variable <- c(fpkm.trans.independent.variable, fpkm.trans.split[[i]][2])
    }
    fpkm.trans$attribute <- factor(fpkm.trans.independent.variable)
    fpkm.pca = FactoMineR::PCA(fpkm.trans, ncp=2, quali.sup=length(fpkm.trans), graph = FALSE)
    eig.val <- factoextra::get_eigenvalue(fpkm.pca)
    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/Dimension_pca_plot.png"))
    p1 <- factoextra::fviz_eig(fpkm.pca, addlabels = TRUE, ylim = c(0, 50), main = "PCA Dimensions") +
      labs(title ="PCA Dimensions") +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
    print(p1)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/Dimension_pca_plot.png"), "' has been created. \n"))
    #var$coord: coordinates of variables to create a scatter plot
    #var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
    #var$contrib: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
    #var <- get_pca_var(res.pca)



    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/PCA_plot_factoextra.png"))
    p2 <- factoextra::fviz_pca_ind(fpkm.pca,
                                   xlab = paste0("PC1(", round(data.frame(eig.val)$variance.percent[1], 2), "%)"), ylab = paste0("PC2(", round(data.frame(eig.val)$variance.percent[2],2), "%)"),
                                   legend.title = "Treatment variable", legend.position = "top",
                                   pointshape = 21,
                                   pointsize = 2.5,
                                   geom.ind = "point", # show points only (nbut not "text")
                                   habillage = fpkm.trans$attribute,
                                   fill.ind = fpkm.trans$attribute,
                                   col.ind = fpkm.trans$attribute, # color by groups
                                   addEllipses=TRUE
                                   #palette = c("#00AFBB", "#E7B800"),
                                   #                 addEllipses = TRUE, # Concentration ellipses
    ) +
      labs(title ="Principal Component Analysis") +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
    print(p2)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/PCA_plot_factoextra.png"), "' has been created. \n"))

    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/PCA_plot_self.png"))
    FPKM.res.PCA = FactoMineR::PCA(fpkm.trans, scale.unit=TRUE, ncp=2, quali.sup=length(fpkm.trans), graph = FALSE)

    my_colors=c(rgb(255, 47, 35,maxColorValue = 255),
                rgb(50, 147, 255,maxColorValue = 255))

    #par(mfrow=c(1,1))
    # fpkm.trans.sort <- fpkm.trans[ order(row.names(fpkm.trans)), ]
    #FPKM.res.PCA$ind$coord <- FPKM.res.PCA$ind$coord[ order(row.names(FPKM.res.PCA$ind$coord)), ]
    plot(FPKM.res.PCA$ind$coord[,1] , FPKM.res.PCA$ind$coord[,2], main = "PCA  Plot", xlab=paste0("PC1(", round(FPKM.res.PCA$eig[,2][1], 2), "%)") , ylab=paste0("PC2(", round(FPKM.res.PCA$eig[,2][2], 2), "%)") , pch=20 , cex=3 ,
         col=my_colors[as.numeric(FPKM.res.PCA$call$quali.sup$quali.sup[,1])] )
    #my_colors[as.numeric(res.PCA$call$quali.sup$quali.sup[,1])]
    abline(h=0 , v=0, lty= 2)
    par(xpd=TRUE)
    legend("bottomright",inset=c(0,1), horiz=TRUE, bty="n", legend=levels(FPKM.res.PCA$call$quali.sup$quali.sup[,1] ) , col=my_colors, pch=20 )
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/PCA_plot.png"), "' has been created. \n\n"))
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

DEHeatmap <- function(path.prefix) {
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv"))){
    # load gene name for further usage
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/Heatmap/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/Heatmap/"))
    }
    cat(paste0("\u25CF Plotting Differential Expressed Heatmap related plot\n"))
    file.path <- paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv")
    data.read.csv <- read.csv(file = file.path)
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    sample.table <- as.data.frame(table(pheno_data[2]))
    new.data.frame.index <- c()
    for(i in 1:length(row.names(sample.table))){
      current.sum <- 0
      if (i-1 == 0 ) current.sum = 0
      else {
        for(z in 1:(i-1)) {
          current.sum <- current.sum + sample.table$Freq[z]
        }
      }
      for(j in 1:sample.table$Freq[i]){
        column.number <- 4+j+current.sum + i -1
        vector <- c(vector, column.number)
        new.data.frame.index <- c(new.data.frame.index, column.number)
      }
    }
    sub.data.frame <- data.read.csv[new.data.frame.index]
    row.names(sub.data.frame) <- data.read.csv$transcriptNames
    df.new <- scale(sub.data.frame)
    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/Heatmap/Heatmap_plot.png"), width = 1000, height = 1000)
    # par(mar=c(10,4,4,2))
    p1 <- heatmap(df.new, scale = "row", xlab = "samples", ylab = "transcript names",cexRow=1, cexCol = 1, margins = c(10,8))
      # theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
    print(p1)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/Heatmap/Heatmap_plot.png"), "' has been created. \n"))
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}


#'
#' @export
DEGOFunctionalAnalysis <- function(path.prefix) {
  ballgown.FPKM.path <- paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv")
  ballgown.FPKM.path.univ <- paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv")
  DE.csv <- read.csv(ballgown.FPKM.path)
  Univ.csv <- read.csv(ballgown.FPKM.path.univ)
  # DE gene
  gene_list_SYMBOL <- DE.csv[DE.csv$geneNames != ".",]$FC
  gene_list_ENTREZID <- DE.csv[DE.csv$geneNames != ".",]$FC
  gene_name <- as.character(DE.csv[DE.csv$geneNames != ".",]$geneNames)
  # all ballgown gene
  gene_list_SYMBOL_univ <- Univ.csv[Univ.csv$geneNames != ".",]$FC
  gene_list_ENTREZID_univ <- Univ.csv[Univ.csv$geneNames != ".",]$FC
  gene_name_univ <- as.character(Univ.csv[Univ.csv$geneNames != ".",]$geneNames)
  # Rename gene_list_SYMBOL
  names(gene_list_SYMBOL) <-gene_name
  names(gene_list_SYMBOL_univ) <-gene_name_univ
  # Sort gene_list_SYMBOL
  gene_list_SYMBOL = sort(gene_list_SYMBOL, decreasing = TRUE)
  gene_list_SYMBOL_univ = sort(gene_list_SYMBOL_univ, decreasing = TRUE)

  # GO classification
  # GO classification : groupGO designed for gene classification based on GO distribution.
  # IDs conversion
  gene.df.DE <- clusterProfiler::bitr(gene_name, fromType = "SYMBOL",
                                      toType = c("ENTREZID", "ENSEMBL"),
                                      OrgDb = org.Hs.eg.db)
  gene.df.Univ <- clusterProfiler::bitr(gene_name_univ, fromType = "SYMBOL",
                                        toType = c("ENTREZID", "ENSEMBL"),
                                        OrgDb = org.Hs.eg.db)
  # Get ENTREZID for DE and Univ
  # DE
  ENTREZID_IDs.DE <- c()
  ENTREZID_IDs.DE <- lapply(gene_name, find_ENTREZID_ID_DE, gene.df.DE = gene.df.DE)
  names(gene_list_ENTREZID) <- ENTREZID_IDs.DE
  gene_list_ENTREZID = sort(gene_list_ENTREZID, decreasing = TRUE)
  # Univ
  ENTREZID_IDs_Univ <- c()
  ENTREZID_IDs_Univ <- lapply(gene_name_univ, find_ENTREZID_ID_Univ, gene.df.Univ = gene.df.Univ)
  names(gene_list_ENTREZID_univ) <- ENTREZID_IDs_Univ
  gene_list_ENTREZID_univ = sort(gene_list_ENTREZID_univ, decreasing = TRUE)

  # designed for gene classification on GO distribution at a specific level. "MF", "BP", "CC"
  ggo <- clusterProfiler::groupGO(gene     = names(gene_list_ENTREZID),
                                  OrgDb    = org.Hs.eg.db,   # variable
                                  ont      = "CC",           # variable
                                  level    = 3,              # Not sure
                                  readable = TRUE)
  ggo.data.frame <- data.frame(ggo)
  write.csv(ggo.data.frame, file = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO/go_classification.csv"))
  png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO/go_classification.png"))
  p1 <- barplot(ggo, drop=TRUE, showCategory=12)
  print(p1)
  dev.off()

  # GO over-representeation test
  ego <- clusterProfiler::enrichGO(gene          = names(gene_list_ENTREZID),
                                   universe      = names(gene_list_ENTREZID_univ),
                                   OrgDb         = org.Hs.eg.db,                   # variable
                                   ont           = "CC",                           # variable
                                   pAdjustMethod = "BH",                           # variable
                                   pvalueCutoff  = 0.01,                           # variable : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                                   qvalueCutoff  = 0.05,                           # variable
                                   readable      = TRUE)
  ego.data.frame <- data.frame(ego)
  write.csv(ego.data.frame, file = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO/go_overrepresentation_test.csv"))
  # bar plot
  png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO/go_overrepresentation_bar_plot.png"))
  p1 <- barplot(ego, showCategory=8)
  print(p1)
  dev.off()
  # dot plot
  png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO/go_overrepresentation_dot_plot.png"))
  p1 <- clusterProfiler::dotplot(ego)
  print(p1)
  dev.off()

  clusterProfiler::emapplot(ego)
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  clusterProfiler::cnetplot(ego, categorySize="pvalue", foldChange=geneList)
  clusterProfiler::goplot(ego)


  data(geneList, package="DOSE")
  gene <- names(geneList)[abs(geneList) > 2]
  gene.df <- bitr(gene, fromType = "ENTREZID",
                  toType = c("ENSEMBL", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  head(gene.df)

  ego <- clusterProfiler::enrichGO(gene          = gene,
                  universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  head(ego)

  ego2 <- clusterProfiler::enrichGO(gene         = gene.df$ENSEMBL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

  barplot(ego2, showCategory=8)
  clusterProfiler::dotplot(ego2)
  clusterProfiler::emapplot(ego2)
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  clusterProfiler::cnetplot(ego2, categorySize="pvalue", foldChange=geneList)
  clusterProfiler::goplot(ego2)



  # GO Gene Set Enrichment Analysis
  ego3 <- clusterProfiler::gseGO(geneList     = gene_list_ENTREZID_univ,
                                 OrgDb        = org.Hs.eg.db,
                                 ont          = "CC",
                                 nPerm        = 1000,
                                 minGSSize    = 100,
                                 maxGSSize    = 500,
                                 pvalueCutoff = 0.1,
                                 verbose      = FALSE)
  ego3 <- clusterProfiler::gseGO(geneList     = geneList,
                                 OrgDb        = org.Hs.eg.db,
                                 ont          = "CC",
                                 nPerm        = 1000,
                                 minGSSize    = 100,
                                 maxGSSize    = 500,
                                 pvalueCutoff = 0.1,
                                 verbose      = FALSE)
  head(ego3)
}




#'
#' @export
DEBallgownPlotAll <- function(path.prefix, independent.variable) {
  cat(paste0("************** Ballgown result visualization **************\n"))
  DEBallgownPCAPlot(path.prefix, independent.variable)
  DEHeatmap(path.prefix)
}

