# BallgownPCAPlot
DEBallgownPCAPlot <- function(path.prefix, independent.variable){
  # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv"))){
    # load gene name for further usage
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/PCA/"))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/PCA/"))
    }
    cat(paste0("\u25CF Plotting Differential Expressed PCA related plot\n"))
    # read pheno data
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    # input sample FPKM
    return.sample.data.frame <- ParseFPKMBallgownResult(path.prefix, independent.variable, control.group, experiment.group, paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv"))
    # get the control and experiment sample size
    sample.table <- as.data.frame(table(pheno_data[independent.variable]))
    control.group.size <- sample.table[sample.table$Var1 == control.group,]$Freq
    experiment.group.size <- sample.table[sample.table$Var1 == experiment.group,]$Freq
    grp = factor(c(rep(control.group, control.group.size), rep(experiment.group, experiment.group.size)))
    fpkm.trans <- data.frame(t(return.sample.data.frame))
    fpkm.trans$attribute <- grp
    fpkm.pca = FactoMineR::PCA(fpkm.trans, ncp=2, quali.sup=length(fpkm.trans), graph = FALSE)
    eig.val <- factoextra::get_eigenvalue(fpkm.pca)
    png(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/PCA/Dimension_pca_plot.png"))
    p1 <- factoextra::fviz_eig(fpkm.pca, addlabels = TRUE, ylim = c(0, 50), main = "PCA Dimensions") +
      labs(title ="PCA Dimensions") +
      theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
    print(p1)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/PCA/Dimension_pca_plot.png"), "' has been created. \n"))
    #var$coord: coordinates of variables to create a scatter plot
    #var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
    #var$contrib: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
    #var <- get_pca_var(res.pca)
    png(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/PCA/PCA_plot_factoextra.png"))
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
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/PCA/PCA_plot_factoextra.png"), "' has been created. \n"))

    png(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/PCA/PCA_plot_self.png"))
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
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/PCA/PCA_plot.png"), "' has been created. \n\n"))
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

DEHeatmap <- function(path.prefix, independent.variable, control.group, experiment.group) {
  if(file.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv"))){
    # load gene name for further usage
    # if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/Heatmap/"))){
    #   dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/Heatmap/"))
    # }
    cat(paste0("\u25CF Plotting Differential Expressed Heatmap related plot\n"))
    # Retuen only independent variable FPKM
    return.data.frame <- ParseFPKMBallgownResult(path.prefix, independent.variable, control.group, experiment.group, paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv"))
    row.names(return.data.frame) <- data.read.csv$transcriptNames
    cat(paste0("     \u25CF Filtering out transcript with no name.\n"))
    return.data.frame.all.valid.name <- subset(return.data.frame, rownames(return.data.frame) != ".")
    cat(paste0("     \u25CF Checking found differential express transcript term.\n"))
    if (length(row.names(return.data.frame.all.valid.name)) == 0) {
      cat(paste0("          \u25CF (\u26A0) No term were found.\n"))
    } else {
      if (length(row.names(return.data.frame.all.valid.name)) > 50) {
        cat(paste0("          \u25CF Found ", length(row.names(return.data.frame.all.valid.name)), " terms. More than 50 terms (Only plot top 50 smallest p value).\n"))
        return.data.frame.all.valid.name <- return.data.frame.all.valid.name[seq_len(50),]
      } else {
        cat(paste0("          \u25CF Found ", length(row.names(return.data.frame.all.valid.name)), " terms.\n"))
      }
    }
    cat(paste0("     \u25CF Calculating log2(FPKM+1).\n"))
    log.data.frame <- log(return.data.frame+1)
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    sample.table <- as.data.frame(table(pheno_data[independent.variable]))
    control.group.size <- sample.table[sample.table$Var1 == control.group,]$Freq
    control.average <- rowMeans(log.data.frame[seq_len(control.group.size)])
    cat(paste0("     \u25CF Each log2(FPKM+1) minus average of control.\n"))
    log.data.frame.minus <- log.data.frame - control.average
    df.new <- scale(log.data.frame.minus)
    png(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/Heatmap_plot.png"), width = 1000, height = 1000)
    redgreen <- c("red", "white", "blue")
    pal <- colorRampPalette(redgreen)(100)
    p1 <- heatmap(df.new, scale = "row", xlab = "samples", ylab = "transcript names",cexRow=1, cexCol = 1, margins = c(10,8), col = pal)
      # theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
    print(p1)
    dev.off()
    cat(paste0("(\u2714) : '", paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/Heatmap/Heatmap_plot.png"), "' has been created. \n"))
  } else {
    stop("(\u2718) 'ballgown_FPKM_result.csv' haven't created yet.\n\n")
  }
}

#
DEBallgownPlotAll <- function(path.prefix, independent.variable) {
  cat(paste0("************** Ballgown result visualization **************\n"))
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/images/"))
  }
  DEBallgownPCAPlot(path.prefix, independent.variable)
  DEHeatmap(path.prefix)
}

# This package will autamatically install org.Hs.eg.db, org.Rn.eg.db, org.Mm.eg.db. If you want to use different OrgDb annotation species, please install that annotation package and attach to your session.
DEGOAnalysis <- function(path.prefix, OrgDb.species) {
  cat(paste0("\n************** Gene Ontology Analysis **************\n"))
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/"))){
    cat("\u25CF Creating directory : 'RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/'\n")
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/"))
  }
  # get return value
  DEUniv_results <- DEUnivGeneList(path.prefix = path.prefix, OrgDb.species = OrgDb.species)
  gene_list_SYMBOL = DEUniv_results$gene_list_SYMBOL_rt
  gene_list_SYMBOL_univ = DEUniv_results$gene_list_SYMBOL_univ_rt
  gene_list_ENTREZID = DEUniv_results$gene_list_ENTREZID_rt
  gene_list_ENTREZID_univ = DEUniv_results$gene_list_ENTREZID_univ_rt

  # DETranscript.limit <- 400
  GO.Ontology.list <- c("MF", "BP", "CC")
  for ( i in GO.Ontology.list) {
    cat("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
    cat(paste0("\u25CF Gene Ontology Analysis : '", i, "' group ... \n"))
    cat("\u25CF Checking differential expression gene number ... \n")
    # if (length(gene_list_SYMBOL) < DETranscript.limit && length(gene_list_ENTREZID) < DETranscript.limit) {
    cat(paste0("     \u25CF Differential expression gene number : ", length(gene_list_SYMBOL), "\n\n"))
    # Do GO Gene Set Enrichment Analysis
    dir_name <- paste0("GO_Gene_Set_Enrichment")
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name))
    }
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i))
    }
    cat("\u25CF GO Gene Set Enrichment Analysis ... \n")
    gse <- clusterProfiler::gseGO(geneList     = gene_list_ENTREZID_univ,
                                  OrgDb        = OrgDb.species,
                                  ont          = i,
                                  verbose      = FALSE)
    gse.data.frame <- data.frame(gse)
    if (length(row.names(gse.data.frame)) > 0) {
      # cat("     \u25CF Printing GO Gene Set Enrichment Analysis result \n")
      # print(head(gse.data.frame))
      if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images"))){
        dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images"))
      }
      # Result column must bigger than 0
      cat(paste0("     \u25CF (\u2714) GO Gene Set Enrichment test (", i,") enriched term found! \n"))
      cat(paste0("     \u25CF Writing 'GO_", i, "_Gene_Set_Enrichment.csv' \n"))
      write.csv(gse.data.frame, file = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/GO_", i, "_Gene_Set_Enrichment.csv"))
      cat(paste0("     \u25CF Checking 'GO_", i, "_Gene_Set_Enrichment.csv' result row number \n"))
      if (length(row.names(gse.data.frame)) < 5) {
        cat(paste0("          \u25CF 'GO_", i, "_Gene_Set_Enrichment.csv' result row number : ", length(row.names(gse.data.frame)), " (less than 5)\n"))
        for( GO.ID in gse.data.frame$ID) {
          cat(paste0("               \u25CF GO ID : ", GO.ID, "\n"))
          png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images/", GO.ID, "_gseGO.png"))
          cat(paste0("               \u25CF Plotting '", GO.ID, "_gseGO.png'\n"))
          p <- clusterProfiler::gseaplot(gse, geneSetID = GO.ID)
          print(p)
          dev.off()
        }
      } else {
        cat(paste0("          \u25CF 'GO_", i, "_Gene_Set_Enrichment.csv' result row number : ", length(row.names(gse.data.frame)), " (more than 5)\n"))
        for(GO.ID in gse.data.frame$ID[seq_len(5)]) {
          cat(paste0("               \u25CF GO ID : ", GO.ID, "\n"))
          png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images/", GO.ID, "_gseGO.png"))
          cat(paste0("               \u25CF Plotting '", GO.ID, "_gseGO.png'\n"))
          p <- clusterProfiler::gseaplot(gse, geneSetID = GO.ID)
          print(p)
          dev.off()
        }
      }
        cat("\n")
    } else {
      cat("     \u25CF (\u26A0) No enriched term is found.\n\n")
      file.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/GO_GSE_NO_TERM"))
    }
    # Do GO classification and GO over-representation test
    dir_name <- paste0("GO_DE_Classification_Erichment")
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name))
    }
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i))
    }
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images"))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images"))
    }
    cat("\u25CF GO Classification ... \n")
    # designed for gene classification on GO distribution at a specific level. "MF", "BP", "CC"
    ggo <- clusterProfiler::groupGO(gene     = names(gene_list_ENTREZID),
                                    OrgDb    = OrgDb.species,   # variable
                                    ont      = i,           # variable
                                    level    = 3,              # Not sure
                                    readable = TRUE)
    ggo.data.frame <- data.frame(ggo)
    # Condition 1 for GO classification ! Row number have to bigger than 1 !
    if (length(row.names(ggo.data.frame)) > 0) {
      # cat("     \u25CF Printing GO Classification result \n")
      # print(head(ggo.data.frame))
      cat(paste0("     \u25CF (\u2714) GO Classification (", i,") result found! \n"))
      cat(paste0("     \u25CF Writing 'GO_", i, "_Classification.csv' \n"))
      write.csv(ggo.data.frame, file = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i , "/GO_", i, "_Classification.csv"))
      cat(paste0("     \u25CF Plotting 'GO_", i, "_Classification_Bar_plot.png' \n\n"))
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images/GO_", i, "_Classification_Bar_plot.png"))
      p1 <- barplot(ggo, drop=TRUE, showCategory=12)
      print(p1)
      dev.off()
    }  else {
      cat("     \u25CF (\u26A0) No term is found.\n\n")
      file.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/GO_CLASSIFICATION_NO_TERM"))
    }

    cat("\u25CF GO Enrichment Test ... \n")
    # GO over-representation test
    ego <- clusterProfiler::enrichGO(gene          = names(gene_list_ENTREZID),
                                     # universe      = geneList,
                                     OrgDb         = OrgDb.species,                   # variable
                                     ont           = i,
                                     pAdjustMethod = "BH",                            # variable : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                                     readable      = TRUE)
    ego.data.frame <- data.frame(ego)
    if (length(row.names(ego.data.frame)) >= 0) {
      # cat("     \u25CF Printing GO Enrichment result \n")
      # print(head(ggo.data.frame))
      cat(paste0("     \u25CF (\u2714) GO Enrichment test (", i,") enriched term found! \n"))
      cat(paste0("     \u25CF Writing 'GO_", i, "_Enrichment.csv' \n"))
      write.csv(ego.data.frame, file = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/GO_", i, "_Enrichment.csv"))
    } else {
      cat(paste0("     \u25CF (\u26A0) No enriched term is found.\n"))
      file.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/GO_ENRICHMENT_NO_TERM"))
    }
    # Condition 2 for GO Enrichment analysis ! Row numebr have to bigger or equal to 2 !
    if (length(row.names(ego.data.frame)) >= 2) {
      # bar plot
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images/GO_", i, "_Enrichment_Bar_plot.png"))
      cat(paste0("     \u25CF Plotting 'GO_", i, "_Enrichment_Bar_plot.png' \n"))
      p2 <- barplot(ego, showCategory=12)
      print(p2)
      dev.off()

      # dot plot
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images/GO_", i, "_Enrichment_Dot_plot.png"))
      cat(paste0("     \u25CF Plotting 'GO_", i, "_Enrichment_Dot_plot.png' \n"))
      p3 <- clusterProfiler::dotplot(ego)
      print(p3)
      dev.off()

      # have to check before run
      # no enriched term found
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images/GO_", i, "_Enrichment_Map_plot.png"))
      cat(paste0("     \u25CF Plotting 'GO_", i, "_Enrichment_Map_plot.png' \n"))
      p4 <- clusterProfiler::emapplot(ego)
      print(p4)
      dev.off()

      ## categorySize can be scaled by 'pvalue' or 'geneNum'
      # the data frame should contain at least two columns
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images/GO_", i, "_Enrichment_Complex_plot.png"))
      cat(paste0("     \u25CF Plotting 'GO_", i, "_Enrichment_Complex_plot.png' \n"))
      p5 <- clusterProfiler::cnetplot(ego, categorySize="pvalue", foldChange = gene_list_ENTREZID)
      print(p5)
      dev.off()

      # keys must be supplied in a character vector with no NAs
      png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/images/GO_", i, "_Enrichment_Induced_plot.png"))
      cat(paste0("     \u25CF Plotting 'GO_", i, "_Enrichment_Induced_plot.png' \n"))
      p6 <- clusterProfiler::goplot(ego)
      print(p6)
      dev.off()
      cat("\n")
    } else {
      cat(paste0("     \u25CF Row size of 'GO_", i,"_Enrichment.csv' is smaller than 2. Can't draw.\n"))
      file.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/GO_analysis/", dir_name, "/", i, "/GO_ENRICHMENT_LESS_THAN_2"))
    }
  }
}

DEKEGGAnalysis <- function(path.prefix, OrgDb.species, KEGG.organism) {
  cat(paste0("\n************** Kyoto Encyclopedia of Genes and Genomes Analysis **************\n"))
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/"))){
    cat("\u25CF Creating directory : 'RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/'\n")
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/"))
  }
  DEUniv_results <- DEUnivGeneList(path.prefix = path.prefix, OrgDb.species = OrgDb.species)
  gene_list_SYMBOL = DEUniv_results$gene_list_SYMBOL_rt
  gene_list_SYMBOL_univ = DEUniv_results$gene_list_SYMBOL_univ_rt
  gene_list_ENTREZID = DEUniv_results$gene_list_ENTREZID_rt
  gene_list_ENTREZID_univ = DEUniv_results$gene_list_ENTREZID_univ_rt
  cat("\u25CF Checking differential expression gene number ... \n")
  cat(paste0("     \u25CF Differential expression gene number : ", length(gene_list_SYMBOL), "\n\n"))
  # DETranscript.limit <- 400
  # Do KEGG Gene Set Enrichment Analysis
  dir_name <- paste0("KEGG_Gene_Set_Enrichment")
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name))){
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name))
  }
  cat("\u25CF KEGG Gene Set Enrichment Analysis ... \n")
  # DO KEGG Gene Set Enrichment Analysis
  kk.gse <- clusterProfiler::gseKEGG(geneList     = gene_list_ENTREZID_univ,
                                     organism     = KEGG.organism,
                                     verbose      = FALSE)
  kk.gse.frame <- data.frame(kk.gse)
  if (length(row.names(kk.gse.frame)) > 0) {
    # cat("     \u25CF Printing KEGG Gene Set Enrichment Analysis result \n")
    # print(head(kk.gse.frame))
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/images"))){
      dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/images"))
    }
    # Result column must bigger than 0
    cat(paste0("     \u25CF (\u2714) KEGG Gene Set Enrichment test enriched term found! \n"))
    cat(paste0("     \u25CF Writing 'KEGG_Gene_Set_Enrichment.csv' \n"))
    write.csv(kk.gse.frame, file = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/KEGG_Gene_Set_Enrichment.csv"))
    cat(paste0("     \u25CF Checking 'KEGG_Gene_Set_Enrichment.csv' result row number \n"))
    if (length(row.names(kk.gse.frame)) < 5) {
      cat(paste0("          \u25CF 'KEGG_Gene_Set_Enrichment.csv' result row number : ", length(row.names(kk.gse.frame)), " (less than 5)\n"))
      for( KEGG.ID in kk.gse.frame$ID) {
        print(KEGG.ID)
        png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/images/", KEGG.ID, "_gseGO.png"))
        p <- clusterProfiler::gseaplot(kk.gse, geneSetID = KEGG.ID)
        print(p)
        dev.off()
      }
    } else {
      cat(paste0("          \u25CF 'KEGG_Gene_Set_Enrichment.csv' result row number : ", length(row.names(kk.gse.frame)), " (more than 5)\n"))
      for(KEGG.ID in kk.gse.frame$ID[seq_len(5)]) {
        cat(paste0("               \u25CF KEGG ID : ", KEGG.ID, "\n"))
        png(filename = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/images/", KEGG.ID, "_gseGO.png"))
        cat(paste0("               \u25CF Plotting '", KEGG.ID, "_gseKEGG.png'\n"))
        p <- clusterProfiler::gseaplot(kk.gse, geneSetID = KEGG.ID)
        print(p)
        dev.off()
      }
    }
  } else {
    cat("     \u25CF (\u26A0) No enriched term is found.\n\n")
    file.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/KEGG_GSE_NO_TERM"))
  }
  dir_name <- paste0("KEGG_DE_Erichment")
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name))){
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name))
  }
  # Do KEGG over-representation test
  # organism : species supported at 'http://www.genome.jp/kegg/catalog/org_list.html'
  cat("\u25CF KEGG Enrichment Test ... \n")
  # KEGG Enrichment test
  kk <- clusterProfiler::enrichKEGG(gene         = names(gene_list_ENTREZID),
                                    organism     = KEGG.organism)                     # variable
  kk.data.frame <- data.frame(kk)
  # Row size have to bigger than 0!
  if (length(row.names(kk.data.frame)) > 0) {
    # cat("     \u25CF Printing KEGG Enrichment result \n")
    # print(head(kk.data.frame))
    cat(paste0("     \u25CF (\u2714) KEGG Enrichment test enriched term found! \n"))
    cat(paste0("     \u25CF Writing 'KEGG_Enrichment.csv' \n"))
    write.csv(kk.data.frame, file = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/KEGG_Enrichment.csv"))
    for ( i in kk.data.frame$ID) {
      if(!dir.exists(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/", i))){
        print(paste0("     \u25CF Creating directory 'RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/", i, "\n"))
        dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/", i))
      }
      current.path <- getwd()
      # get the url from KEGG result
      cat(paste0("     \u25CF Finding '", i, "' KEGG URL ... \n"))
      KEGGUrl <- GetKEGGUrl(kk, i)
      cat(paste0("     \u25CF Writting 'URL_", i, "_Pathway.txt' \n"))
      write(KEGGUrl, file = paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/", i, "/URL_", i, "_Pathway.txt"))
      # drawing pathway picture with 'pathway' package
      pathway.dir <- paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name, "/", i, "/pathview_result/")
      if(!dir.exists(pathway.dir)){
        dir.create(pathway.dir)
      }
      setwd(pathway.dir)
      cat(paste0("     \u25CF Plotting '", i, "' pathway by package \"pathview\" \n"))
      pathview::pathview(gene.data  = names(gene_list_ENTREZID),
                         pathway.id = i,
                         species    = KEGG.organism,
                         kegg.dir   = pathway.dir)
      on.exit(setwd(current.path))
    }
  } else {
    cat(paste0("     \u25CF (\u26A0) No enriched term is found.\n"))
    file.create(paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/KEGG_analysis/", dir_name,"/KEGG_ENRICHMENT_NO_TERM"))
  }
}

DEUnivGeneList <- function(path.prefix, OrgDb.species) {
  ballgown.FPKM.path <- paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/Differential_Expression/ballgown_FPKM_DE_result.csv")
  ballgown.FPKM.path.univ <- paste0(path.prefix, "RNAseq_results/Ballgown_FPKM_analysis/ballgown_FPKM_result.csv")
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
  # OrgDb.species!!
  gene.df.DE <- clusterProfiler::bitr(gene_name, fromType = "SYMBOL",
                                      toType = c("ENTREZID", "ENSEMBL"),
                                      OrgDb = OrgDb.species)
  gene.df.Univ <- clusterProfiler::bitr(gene_name_univ, fromType = "SYMBOL",
                                        toType = c("ENTREZID", "ENSEMBL"),
                                        OrgDb = OrgDb.species)
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



