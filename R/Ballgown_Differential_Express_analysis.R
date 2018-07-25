#' BallgownPCAPlot
DEBallgownPCAPlot <- function(path.prefix){
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
    sub.data.frame<- data.read.csv[new.data.frame.index]
    fpkm.trans <- data.frame(t(sub.data.frame))
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
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/PCA_plot_factoextra.png"), "' has been created. \n"))

    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/PCA/PCA_plot_self.png"))
    # fpkm.trans.sort$attribute <- factor(fpkm.trans.sort$attribute)
    #length(fpkm.trans.sort)
    FPKM.res.PCA = FactoMineR::PCA(fpkm.trans.sort, scale.unit=TRUE, ncp=2, quali.sup=length(fpkm.trans.sort), graph = FALSE)

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

DEGOFunctionalAnalysis <- function(path.prefix) {

  # GO classification
  # GO classification : groupGO designed for gene classification based on GO distribution.
  data(geneList, package="DOSE")
  gene <- names(geneList)[abs(geneList) > 2]
  gene_2 <- names(geneList)
  gene.df <- clusterProfiler::bitr(gene, fromType = "ENTREZID",
                                   toType = c("ENSEMBL", "SYMBOL"),
                                   OrgDb = org.Hs.eg.db)
  head(gene.df)

  ggo <- clusterProfiler::groupGO(gene     = gene,
                                  OrgDb    = org.Hs.eg.db,
                                  ont      = "CC",
                                  level    = 3,
                                  readable = TRUE)

  head(ggo)

  # GO over-representeation test
  ego <- clusterProfiler::enrichGO(gene          = gene,
                                   universe      = names(geneList),
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "CC",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
  head(ego)

  ego2 <- clusterProfiler::enrichGO(gene          = gene.df$ENSEMBL,
                                    OrgDb         = org.Hs.eg.db,
                                    keyType       = 'ENSEMBL',
                                    ont           = "CC",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.05)
  head(ego2)

  # GO Gene Set Enrichment Analysis
  ego3 <- clusterProfiler::gseGO(geneList     = geneList,
                                 OrgDb        = org.Hs.eg.db,
                                 ont          = "CC",
                                 nPerm        = 1000,
                                 minGSSize    = 100,
                                 maxGSSize    = 500,
                                 pvalueCutoff = 0.05,
                                 verbose      = FALSE)
  head(ego3)
}


own_gene_list <- paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv")
d <- read.csv(own_gene_list)
gene_list_SYMBOL <- d[d[,1] != ".",][,20]
gene_list_ENTREZID <- d[d[,1] != ".",][,20]
gene_name <- as.character(d[d[,1] != ".",][,1])

names(gene_list_SYMBOL) <-gene_name

gene_list_SYMBOL = sort(gene_list_SYMBOL, decreasing = TRUE)

DE_with_valid_name <- d[d[,1] != ".", ]

ENTREZID_IDs <- c()
ENTREZID_IDs <- lapply(gene_name, find_ENTREZID_ID)
find_ENTREZID_ID <- function(sample.id) {
  ENTREZID_IDs <- c(ENTREZID_IDs, eg[eg$SYMBOL == sample.id, ]["ENTREZID"][[1]][1])
  return(ENTREZID_IDs)
}

names(gene_list_ENTREZID) <- ENTREZID_IDs

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
length(gene)

# Cellular Component (CC)
ggo <- clusterProfiler::groupGO(gene     = names(gene_list_ENTREZID),
                                OrgDb    = org.Hs.eg.db,
                                ont      = "CC",
                                level    = 3,
                                readable = TRUE)

#Molecular Function (MF)
ggo <- clusterProfiler::groupGO(gene     = names(gene_list_ENTREZID),
                                OrgDb    = org.Hs.eg.db,
                                ont      = "MF",
                                level    = 3,
                                readable = TRUE)

# Biological Process (BP)
ggo <- clusterProfiler::groupGO(gene     = names(gene_list_ENTREZID),
                                OrgDb    = org.Hs.eg.db,
                                ont      = "BP",
                                level    = 3,
                                readable = TRUE)

# ‘kegg’, ‘ncbi-geneid’, ‘ncbi-proteinid’ or ‘uniprot’.
# ‘kegg’ ID is entrezgene ID for eukaryote species and Locus ID for prokaryotes.
eg2np <- clusterProfiler::bitr_kegg(eg$ENTREZID, fromType='kegg', toType='ncbi-proteinid', organism='hsa')
head(eg2np)

keytypes(org.Hs.eg.db)

clusterProfiler::bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')

# GO analyses (groupGO(), enrichGO() and gseGO()) support organisms that have an OrgDb object available.

gene
eg$ENTREZID
data(geneList, package="DOSE")
eg$ENTREZID <- names(gene_FC)
gene.df <- clusterProfiler::bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)



head(ggo)

hub <- AnnotationHub()
q <- query(hub, "Cricetulus")
id <- q$ah_id[length(q)]
Cgriseus <- hub[[id]]





# GO classification
# GO classification : groupGO designed for gene classification based on GO distribution.
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene_2 <- names(geneList)
gene.df <- clusterProfiler::bitr(gene, fromType = "ENTREZID",
                                 toType = c("ENSEMBL", "SYMBOL"),
                                 OrgDb = org.Hs.eg.db)
head(gene.df)

ggo <- clusterProfiler::groupGO(gene     = gene,
                                OrgDb    = org.Hs.eg.db,
                                ont      = "CC",
                                level    = 3,
                                readable = TRUE)

head(ggo)

# GO over-representeation test
ego <- clusterProfiler::enrichGO(gene          = gene,
                                 universe      = names(geneList),
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "CC",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)
head(ego)

ego2 <- clusterProfiler::enrichGO(gene          = gene.df$ENSEMBL,
                                  OrgDb         = org.Hs.eg.db,
                                  keyType       = 'ENSEMBL',
                                  ont           = "CC",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05)
head(ego2)

# GO Gene Set Enrichment Analysis
ego3 <- clusterProfiler::gseGO(geneList     = geneList,
                               OrgDb        = org.Hs.eg.db,
                               ont          = "CC",
                               nPerm        = 1000,
                               minGSSize    = 100,
                               maxGSSize    = 500,
                               pvalueCutoff = 0.05,
                               verbose      = FALSE)
head(ego3)


