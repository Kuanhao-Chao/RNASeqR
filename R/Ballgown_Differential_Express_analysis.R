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


#' This package will autamatically install org.Hs.eg.db, org.Rn.eg.db, org.Mm.eg.db. If you want to use different OrgDb annotation species, please install that annotation package and attach to your session.
#' @export
DEGOAnalysis <- function(path.prefix, GO.OrgDb.species) {
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/"))
  }
  # get return value
  DEUniv_results <- DEUnivGeneList(path.prefix = path.prefix, GO.OrgDb.species = GO.OrgDb.species)
  gene_list_SYMBOL = DEUniv_results$gene_list_SYMBOL_rt
  gene_list_SYMBOL_univ = DEUniv_results$gene_list_SYMBOL_univ_rt
  gene_list_ENTREZID = DEUniv_results$gene_list_ENTREZID_rt
  gene_list_ENTREZID_univ = DEUniv_results$gene_list_ENTREZID_univ_rt

  GO.Ontology.list <- c("MF", "BP", "CC")
  for ( i in GO.Ontology.list) {
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i))
    }
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i, "/images"))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i, "/images"))
    }
    # designed for gene classification on GO distribution at a specific level. "MF", "BP", "CC"
    ggo <- clusterProfiler::groupGO(gene     = names(gene_list_ENTREZID),
                                    OrgDb    = GO.OrgDb.species,   # variable
                                    ont      = i,           # variable
                                    level    = 3,              # Not sure
                                    readable = TRUE)
    ggo.data.frame <- data.frame(ggo)
    print(ggo.data.frame)

    # Condition 1 for GO classification ! Row number have to bigger than 1 !
    if (length(row.names(ggo.data.frame))!=0) {
      write.csv(ggo.data.frame, file = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i , "/GO_", i, "_Classification.csv"))
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i, "/images/GO_", i, "_Classification_Bar_plot.png"))
      p1 <- barplot(ggo, drop=TRUE, showCategory=12)
      print(p1)
      dev.off()
    } else {
      cat(paste0("Invalid! Row size of 'GO_", i,"_Classification.csv' is 0.\n"))
    }


    # GO over-representeation test
    ego <- clusterProfiler::enrichGO(gene          = names(gene_list_ENTREZID),
                                     # universe      = geneList,
                                     OrgDb         = GO.OrgDb.species,                   # variable
                                     ont           = i,
                                     pAdjustMethod = "BH",                            # variable : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                                     readable      = TRUE)
    ego.data.frame <- data.frame(ego)
    print(ego.data.frame)

    # Condition 2 for GO Enrichment analysis ! Row numebr have to bigger or equal to 2 !
    if (length(row.names(ego.data.frame)) >= 2) {
      write.csv(ego.data.frame, file = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i, "/GO_Enrichment.csv"))
      # bar plot
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i, "/images/GO_Enrichment_Bar_plot.png"))
      p2 <- barplot(ego, showCategory=12)
      print(p2)
      dev.off()

      # dot plot
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i, "/images/GO_Enrichment_Dot_plot.png"))
      p3 <- clusterProfiler::dotplot(ego)
      print(p3)
      dev.off()

      # have to check before run
      # no enriched term found
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i, "/images/GO_Enrichment_Map_plot.png"))
      p4 <- clusterProfiler::emapplot(ego)
      print(p4)
      dev.off()

      ## categorySize can be scaled by 'pvalue' or 'geneNum'
      # the data frame should contain at least two columns
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i, "/images/GO_Enrichment_Complex_plot.png"))
      p5 <- clusterProfiler::cnetplot(ego, categorySize="pvalue", foldChange = gene_list_ENTREZID)
      print(p5)
      dev.off()

      # keys must be supplied in a character vector with no NAs
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/GO_analysis/", i, "/images/GO_Enrichment_Induced_plot.png"))
      p6 <- clusterProfiler::goplot(ego)
      print(p6)
      dev.off()
    } else {
      cat(paste0("Invalid! Row size of 'GO_", i,"_Enrichment.csv' is smaller than 2. Can't draw.\n"))
    }
  }
}

DEKEGGAnalysis <- function(path.prefix = path.prefix, KEGG.organism) {
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/KEGG_analysis/"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/KEGG_analysis/"))
  }
  DEUniv_results <- DEUnivGeneList(path.prefix = path.prefix, GO.OrgDb.species = GO.OrgDb.species)
  gene_list_SYMBOL = DEUniv_results$gene_list_SYMBOL_rt
  gene_list_SYMBOL_univ = DEUniv_results$gene_list_SYMBOL_univ_rt
  gene_list_ENTREZID = DEUniv_results$gene_list_ENTREZID_rt
  gene_list_ENTREZID_univ = DEUniv_results$gene_list_ENTREZID_univ_rt

  # to search available KEGG database
  # clusterProfiler::search_kegg_organism('ece', by='kegg_code')

  ecoli <- clusterProfiler::search_kegg_organism('Escherichia coli', by='scientific_name')
  dim(ecoli)
  head(ecoli)

  # KEGG Enrichment test
  # organism : species supported at 'http://www.genome.jp/kegg/catalog/org_list.html'
  kk <- clusterProfiler::enrichKEGG(gene         = names(gene_list_ENTREZID),
                                    organism     = KEGG.organism)                     # variable
  head(kk)
  kk.data.frame <- data.frame(kk)
  KEGGUrl <- GetKEGGUrl(kk, 'hsa05202')
  # write the url into a file
  download.file(KEGGUrl, destfile = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/KEGG_analysis/URL.html"))




  hsa04110 <- pathview(gene.data  = names(gene_list_ENTREZID),
                       pathway.id = "hsa05202",
                       species    = "hsa")




  # KEGG Gene Set Enrichment Analysis
  kk2 <- clusterProfiler:: gseKEGG(geneList     = geneList,
                 organism     = 'hsa',
                 nPerm        = 1000,
                 minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  head(kk2)

  # KEGG Module over-representation test
  mkk <- clusterProfiler::enrichMKEGG(gene = gene,
                     organism = 'hsa')
  # KEGG Module Gene Set Enrichment Analysis
  mkk2 <- clusterProfiler::gseMKEGG(geneList = geneList,
                   species = 'hsa')

  # DAVID functional analysis
  david <- enrichDAVID(gene = gene,
                       idType = "ENTREZ_GENE_ID",
                       listType = "Gene",
                       annotation = "KEGG_PATHWAY",
                       david.user = "clusterProfiler@hku.hk")

}

DEUnivGeneList <- function(path.prefix, GO.OrgDb.species) {
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
  # GO.OrgDb.species!!
  gene.df.DE <- clusterProfiler::bitr(gene_name, fromType = "SYMBOL",
                                      toType = c("ENTREZID", "ENSEMBL"),
                                      OrgDb = GO.OrgDb.species)
  gene.df.Univ <- clusterProfiler::bitr(gene_name_univ, fromType = "SYMBOL",
                                        toType = c("ENTREZID", "ENSEMBL"),
                                        OrgDb = GO.OrgDb.species)
  # Get ENTREZID for DE and Univ
  # DE
  ENTREZID_IDs.DE <- c()
  ENTREZID_IDs.DE <- lapply(gene_name, find_ENTREZID_ID_DE, gene.df.DE = gene.df.DE, ENTREZID_IDs.DE = ENTREZID_IDs.DE)
  names(gene_list_ENTREZID) <- ENTREZID_IDs.DE
  gene_list_ENTREZID = sort(gene_list_ENTREZID, decreasing = TRUE)
  # Univ
  ENTREZID_IDs_Univ <- c()
  ENTREZID_IDs_Univ <- lapply(gene_name_univ, find_ENTREZID_ID_Univ, gene.df.Univ = gene.df.Univ, ENTREZID_IDs_Univ = ENTREZID_IDs_Univ )
  names(gene_list_ENTREZID_univ) <- ENTREZID_IDs_Univ
  gene_list_ENTREZID_univ = sort(gene_list_ENTREZID_univ, decreasing = TRUE)
  return(list(gene_list_SYMBOL_rt = gene_list_SYMBOL,
              gene_list_SYMBOL_univ_rt = gene_list_SYMBOL_univ,
              gene_list_ENTREZID_rt = gene_list_ENTREZID,
              gene_list_ENTREZID_univ_rt = gene_list_ENTREZID_univ))
}

find_ENTREZID_ID_DE <- function(symbol.id, gene.df.DE, ENTREZID_IDs.DE) {
  ENTREZID_IDs.DE <- c(ENTREZID_IDs.DE, gene.df.DE[gene.df.DE$SYMBOL == symbol.id, ]["ENTREZID"][[1]][1])
  return(ENTREZID_IDs.DE)
}

find_ENTREZID_ID_Univ <- function(symbol.id, gene.df.Univ, ENTREZID_IDs_Univ) {
  ENTREZID_IDs_Univ <- c(ENTREZID_IDs_Univ, gene.df.Univ[gene.df.Univ$SYMBOL == symbol.id, ]["ENTREZID"][[1]][1])
  return(ENTREZID_IDs_Univ)
}

GetKEGGUrl <- function(x, pathID) {
  url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", pathID, '/', x[pathID, "geneID"])
  return(url)
}

#'
#' @export
DEBallgownPlotAll <- function(path.prefix, independent.variable) {
  cat(paste0("************** Ballgown result visualization **************\n"))
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Differential_Expression/images/"))
  }
  DEBallgownPCAPlot(path.prefix, independent.variable)
  DEHeatmap(path.prefix)
}

