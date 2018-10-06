# This package will autamatically install org.Hs.eg.db,
#                                         org.Rn.eg.db,
#                                         org.Mm.eg.db,
#                                         org.Sc.sgd.db,
# If you want to use different OrgDb annotation species,
# please install that annotation package and attach to your session.
# How this function works :
# 1. If there is no Univ terms ==> NA
# 2. If there are Univ terms but no DE terms ==> Do Gene_set_analysis but not
#    doing classification and over-representation
# 3. If there are Univ and DE terms ==> Do Gene_set_analysis
#    and classification and over-representation
GOAnalysis <- function(which.analysis,
                       path.prefix,
                       OrgDb.species,
                       go.level,
                       input.TYPE.ID) {
  CheckGoLevel(go.level)
  message(paste0("\u2694\u2694 Gene Ontology Analysis \n"))
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                        which.analysis, "/GO_analysis/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/",
                      which.analysis, "/GO_analysis/"))
  }
  # get return value
  DE_results <- DEGeneList(which.analysis,
                           path.prefix,
                           OrgDb.species,
                           input.TYPE.ID)
  if (length(DE_results) == 1) {
    message("!! Less than 2 differential expressed valid kegg ID is found, ",
            "can not do KEGG GO analysis!!\n\n")
  } else {
    GO.Ontology.list <- c("MF", "BP", "CC")
    for ( i in GO.Ontology.list) {
      message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
              "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
              "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
      message(paste0("\u25CF '", i, "' : \n"))
      #####################
      #### Checking DE ####
      #####################
      message("     \u25CF Checking differential expression gene number ... \n")
      message("          \u25CF Differential expression gene number : ",
              length(DE_results), "\n")
      message("     Found ENTREZID gene ID :
              ", paste(head(DE_results), " "), "...\n")
      # Do GO classification
      dir_name <- "GO_DE_Classification"
      if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                            which.analysis, "/GO_analysis/", dir_name))){
        dir.create(paste0(path.prefix, "RNASeq_results/",
                          which.analysis, "/GO_analysis/", dir_name))
      }
      if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                            which.analysis, "/GO_analysis/", dir_name, "/", i))){
        dir.create(paste0(path.prefix, "RNASeq_results/",
                          which.analysis, "/GO_analysis/", dir_name, "/", i))
      }
      ###########################
      #### GO Classification ####
      ###########################
      message("     \u25CF GO Classification ... \n")
      # designed for gene classification on GO distribution
      # at a specific level. "MF", "BP", "CC"
      ggo <- clusterProfiler::groupGO(gene     = DE_results,
                                      keyType  = "ENTREZID",
                                      OrgDb    = OrgDb.species,   # variable
                                      ont      = i,           # variable
                                      level    = go.level)
      ggo.data.frame <- data.frame(ggo)
      ggo.data.frame <- ggo.data.frame[order(as.numeric(ggo.data.frame$Count),
                                             decreasing = TRUE),]
      "@"(ggo, result) <- ggo.data.frame
      # Condition 1 for GO classification ! Row number have to bigger than 1 !
      if (length(row.names(ggo.data.frame)) > 0) {
        if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis,
                              "/GO_analysis/", dir_name, "/", i, "/images"))){
          dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis,
                            "/GO_analysis/", dir_name, "/", i, "/images"))
        }
        message("          \u25CF (\u2714) GO Classification ",
                "(", i,") result found! \n")
        message("               \u25CF Writing 'GO_", i,
                "_Classification.csv' \n")
        write.csv(ggo.data.frame,
                  file = paste0(path.prefix, "RNASeq_results/", which.analysis,
                                "/GO_analysis/", dir_name, "/", i ,
                                "/GO_", i, "_Classification.csv"))

        # Visualization
        message(paste0("               \u25CF Plotting 'GO_", i,
                       "_Classification_Bar_Plot_clusterProfiler.png' \n"))
        barplot(ggo, drop=TRUE, showCategory=15, font.size = 7,
                title = "GO Classification Bar Plot") +
          theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
                axis.title.x = element_text(size = 10),
                axis.title.y = element_text(size = 10))
        ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                      "/GO_analysis/", dir_name, "/", i, "/images/",
                      "GO_", i, "_Classification_Bar_",
                      "Plot_clusterProfiler.png"),
               dpi = 300,
               width = 7,
               height = 7)
      }  else {
        message("          \u25CF (\u26A0) No term is found.\n")
        file.create(paste0(path.prefix, "RNASeq_results/",
                           which.analysis, "/GO_analysis/", dir_name, "/", i,
                           "/GO_CLASSIFICATION_NO_TERM"))
      }

      # Do GO over-representation test
      dir_name <- "GO_DE_Overrepresentation"
      if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                            which.analysis, "/GO_analysis/", dir_name))){
        dir.create(paste0(path.prefix, "RNASeq_results/",
                          which.analysis, "/GO_analysis/", dir_name))
      }
      if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                            which.analysis, "/GO_analysis/", dir_name, "/", i))){
        dir.create(paste0(path.prefix, "RNASeq_results/",
                          which.analysis, "/GO_analysis/", dir_name, "/", i))
      }
      ################################
      #### GO Over-representation ####
      ################################
      # GO over-representation test !
      message("     \u25CF GO Over-representation Test ... \n")
      # GO over-representation test
      ego <- clusterProfiler::enrichGO(gene          = DE_results,
                                       keyType       = "ENTREZID",
                                       OrgDb         = OrgDb.species,   # variable
                                       ont           = i,
                                       pAdjustMethod = "BH")
      # variable : "holm", "hochberg", "hommel",
      #            "bonferroni", "BH", "BY", "fdr", "none"
      ego.data.frame <- data.frame(ego)
      message(paste0("          \u25CF (\u2714) GO over-representation test (",
                     i,") enriched term found : ",
                     length(row.names(ego.data.frame)), "\n"))
      if (length(row.names(ego.data.frame)) > 0) {
        if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis,
                              "/GO_analysis/", dir_name, "/", i, "/images"))){
          dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis,
                            "/GO_analysis/", dir_name, "/", i, "/images"))
        }
        message(paste0("               \u25CF Writing 'GO_", i,
                       "_Overrepresentation.csv' \n"))
        write.csv(ego.data.frame,
                  file = paste0(path.prefix, "RNASeq_results/", which.analysis,
                                "/GO_analysis/", dir_name, "/", i,
                                "/GO_", i, "_Overrepresentation.csv"))
        # visualization
        # bar plot
        message(paste0("               \u25CF Plotting 'GO_", i,
                       "_Overrepresentation_Bar_Plot_clusterProfiler.png' \n"))
        barplot(ego, showCategory=12, font.size= 7,
                title = "Over-representation Bar Plot") +
          theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
                axis.title.x = element_text(size = 10),
                axis.title.y = element_text(size = 10))
        ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                      "/GO_analysis/", dir_name, "/", i, "/images/GO_", i,
                      "_Overrepresentation_Bar_Plot_clusterProfiler.png"),
               dpi = 300,
               width = 7,
               height = 7)

        # dot plot
        message(paste0("               \u25CF Plotting 'GO_", i,
                       "_Overrepresentation_Dot_Plot_clusterProfiler.png' \n"))
        clusterProfiler::dotplot(ego, font.size = 7,
                                 title = "Over-representation Dot Plot")+
          theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
                axis.title.x = element_text(size = 10),
                axis.title.y = element_text(size = 10))
        ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                      "/GO_analysis/", dir_name, "/", i, "/images/GO_", i,
                      "_Overrepresentation_Dot_Plot_clusterProfiler.png"),
               dpi = 300,
               width = 7,
               height = 7)
      } else {
        message(paste0("          \u25CF (\u26A0) No enriched term is found.\n"))
        file.create(paste0(path.prefix, "RNASeq_results/", which.analysis,
                           "/GO_analysis/", dir_name, "/", i,
                           "/GO_OVERREPRESENTATION_NO_TERM"))
      }
    }
  }
}

KEGGAnalysis <- function(which.analysis,
                         path.prefix,
                         OrgDb.species,
                         input.TYPE.ID,
                         KEGG.organism) {
  message(paste0("\u2694\u2694 Kyoto Encyclopedia of ",
                 "Genes and Genomes Analysis \n"))
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
                        which.analysis, "/KEGG_analysis/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/",
                      which.analysis, "/KEGG_analysis/"))
  }
  DE_results_kegg <- KEGGDEGeneList(which.analysis,
                                    path.prefix,
                                    OrgDb.species,
                                    input.TYPE.ID,
                                    KEGG.organism)
  if (length(DE_results_kegg$return.kegg.id) == 1) {
    message("!! Less than 2 differential expressed valid kegg ID is found, ",
            "can not do KEGG over-representation analysis!!\n\n")
  } else {
    #####################
    #### Checking DE ####
    #####################
    message("\u25CF Checking differential expression gene number ... \n")
    message(paste0("     \u25CF Differential expression gene number : ",
                   length(DE_results_kegg$return.kegg.id), "\n"))
    dir_name <- paste0("KEGG_DE_Overrepresentation")
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis,
                          "/KEGG_analysis/", dir_name))){
      dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis,
                        "/KEGG_analysis/", dir_name))
    }
    ##################################
    #### KEGG Over-representation ####
    ##################################
    # Do KEGG over-representation test
    # organism : species supported at
    #   'http://www.genome.jp/kegg/catalog/org_list.html'
    message("\u25CF KEGG Over-representation Test ... \n")
    # KEGG Over-representation test
    message("     Found kegg gene ID : ",
            paste(head(DE_results_kegg$return.kegg.id), " "), "...\n")
    kk <- clusterProfiler::enrichKEGG(gene         = DE_results_kegg$return.kegg.id,
                                      keyType      = "kegg",
                                      organism     = KEGG.organism,
                                      pvalueCutoff = 0.05)
    kk.data.frame <- data.frame(kk)
    # Row size have to bigger than 0!
    if (length(row.names(kk.data.frame)) > 0) {
      message(paste0("     \u25CF (\u2714) KEGG ",
                     "over-representation test enriched term found! ",
                     length(row.names(kk.data.frame)), "\n"))
      message(paste0("          \u25CF Writing ",
                     "'KEGG_Overrepresentation.csv' \n"))
      # All ID pathway will be plotted!
      write.csv(kk.data.frame,
                file = paste0(path.prefix, "RNASeq_results/",
                              which.analysis, "/KEGG_analysis/",
                              dir_name, "/KEGG_Overrepresentation.csv"))
      # Create directory for visualization
      if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis,
                            "/KEGG_analysis/", dir_name, "/images/"))){
        dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis,
                          "/KEGG_analysis/", dir_name, "/images/"))
      }
      # Do visualization!!
      # bar plot
      message(paste0("               \u25CF Plotting 'KEGG",
                     "_Overrepresentation_Bar_Plot_clusterProfiler.png' \n"))
      barplot(kk, showCategory=12, font.size= 7,
              title = "Over-representation Bar Plot") +
        theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_text(size = 10))
      ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                    "/KEGG_analysis/", dir_name, "/images/KEGG",
                    "_Overrepresentation_Bar_Plot_clusterProfiler.png"),
             dpi = 300,
             width = 7,
             height = 7)

      # dot plot
      message(paste0("               \u25CF Plotting 'KEGG",
                     "_Overrepresentation_Dot_Plot_clusterProfiler.png' \n"))
      clusterProfiler::dotplot(kk, font.size = 7,
                               title = "Over-representation Dot Plot")+
        theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_text(size = 10))
      ggsave(paste0(path.prefix, "RNASeq_results/", which.analysis,
                    "/KEGG_analysis/", dir_name, "/images/KEGG",
                    "_Overrepresentation_Dot_Plot_clusterProfiler.png"),
             dpi = 300,
             width = 7,
             height = 7)
      # Visualize and create url for only top 5 !!
      if (nrow(kk.data.frame) > 5) {
        kk.data.frame <- kk.data.frame[1:5, ]
      }
      for ( i in kk.data.frame$ID) {
        if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis,
                              "/KEGG_analysis/", dir_name, "/", i))){
          dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis,
                            "/KEGG_analysis/", dir_name, "/", i))
        }
        current.path <- getwd()
        # get the url from KEGG result
        message(paste0("          \u25CF Finding '", i, "' KEGG URL ... \n"))
        KEGGUrl <- GetKEGGUrl(kk, i)
        message(paste0("               \u25CF Writting 'URL_",
                       i, "_Pathway.txt' \n"))
        write(KEGGUrl,
              file = paste0(path.prefix, "RNASeq_results/", which.analysis,
                            "/KEGG_analysis/", dir_name, "/", i,
                            "/URL_", i, "_Pathway.txt"))
        # drawing pathway picture with 'pathway' package
        pathway.dir <- paste0(path.prefix, "RNASeq_results/", which.analysis,
                              "/KEGG_analysis/", dir_name, "/", i,
                              "/pathview_result/")
        if(!dir.exists(pathway.dir)){
          dir.create(pathway.dir)
        }
        setwd(pathway.dir)
        message("               \u25CF Plotting '", i,
                "' pathway by package \"pathview\" \n")
        pathview::pathview(gene.data  = DE_results_kegg$return.ENTREZID,
                           pathway.id = i,
                           species    = KEGG.organism,
                           kegg.dir   = pathway.dir)
        on.exit(setwd(current.path))
      }
    } else {
      message("     \u25CF (\u26A0) No over-representation term is found.\n")
      file.create(paste0(path.prefix, "RNASeq_results/", which.analysis,
                         "/KEGG_analysis/", dir_name,
                         "/KEGG_OVERREPRESENTATION_NO_TERM"))
    }
  }
  message("===============================================================\n",
          "===============================================================\n\n")
}

DEGeneList <- function(which.analysis,
                       path.prefix,
                       OrgDb.species,
                       input.TYPE.ID) {
  DE.path.csv <- paste0(path.prefix, "RNASeq_results/", which.analysis,
                        "/", strsplit(which.analysis, "_")[[1]][1],
                        "_normalized_DE_result.csv")
  DE.csv <- read.csv(DE.path.csv)
  # DE gene
  # First filter out "." gene name
  DE.csv <- DE.csv[DE.csv$gene.name != ".",]
  gene_name <- as.character(DE.csv$gene.name)

  ############
  #### DE ####
  ############
  if (length(gene_name) == 0) {
    message(paste0("No annotated gene terms are found in '",
                   path.prefix, "RNASeq_results/", which.analysis, "/",
                   strsplit(which.analysis, "_")[[1]][1],
                   "_normalized_DE_result.csv'\n\n"))
    return(return.id = NA)
  } else {
    gene.df.DE <- clusterProfiler::bitr(gene_name,
                                        fromType = input.TYPE.ID,
                                        toType = c("ENTREZID"),
                                        OrgDb = OrgDb.species)
    return.id <- gene.df.DE$ENTREZID
    return.id <- return.id[!is.na(return.id)]
    if (length(return.id) <= 1) {
      return(return.id = NA)
    } else {
      return(return.id = return.id)
    }
  }
}

KEGGDEGeneList <- function(which.analysis,
                           path.prefix,
                           OrgDb.species,
                           input.TYPE.ID,
                           KEGG.organism) {
  DE.path.csv <- paste0(path.prefix, "RNASeq_results/", which.analysis,
                        "/", strsplit(which.analysis, "_")[[1]][1],
                        "_normalized_DE_result.csv")
  DE.csv <- read.csv(DE.path.csv)
  # DE gene
  # First filter out "." gene name
  DE.csv <- DE.csv[DE.csv$gene.name != ".",]
  ## ==> So every gene in DE would be distinct (with the larget log2FC)
  gene_name <- as.character(DE.csv$gene.name)

  ############
  #### DE ####
  ############
  if (length(gene_name) == 0) {
    message(paste0("No annotated gene terms are found in '",
                   path.prefix, "RNASeq_results/", which.analysis, "/",
                   strsplit(which.analysis, "_")[[1]][1],
                   "_normalized_DE_result.csv'\n\n"))
    return(list(return.kegg.id = NA,
                return.ENTREZID = NA))
  } else {
    # IDs conversion
    gene.df.DE <- clusterProfiler::bitr(gene_name,
                                        fromType = input.TYPE.ID,
                                        toType = c("UNIPROT", "ENTREZID"),
                                        OrgDb = OrgDb.species)
    return.ENTREZID <- gene.df.DE$ENTREZID
    return.ENTREZID <- return.ENTREZID[!is.na(return.ENTREZID)]
    # It is ok if gene.df.DE is NA~~
    gene.df.DE.KEGG <- clusterProfiler::bitr_kegg(gene.df.DE$UNIPROT,
                                                  fromType = "uniprot",
                                                  toType = "kegg",
                                                  organism = KEGG.organism)
    return.kegg.id <- gene.df.DE.KEGG$kegg
    return.kegg.id <- return.kegg.id[!is.na(return.kegg.id)]
    if (length(return.kegg.id) <= 1) {
      return(list(return.kegg.id = NA,
                  return.ENTREZID = NA))
    } else {
      return(list(return.kegg.id = return.kegg.id,
                  return.ENTREZID = return.ENTREZID))
    }
  }
}

GetKEGGUrl <- function(x, pathID) {
  url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?",
                pathID, '/', x[pathID, "geneID"])
  return(url)
}

CheckGoLevel <- function(go.level) {
  whether.integer <- go.level%%1 == 0
  message(paste0("\u25CF Checking 'go.level' value : ", go.level, "\n"))
  if (whether.integer) {
    message("(\u2714) 'go.level' is integer! Valid\n\n")
  } else {
    message("(\u2718) 'go.level' is not integer!\n\n")
    stop("'go.level' not integer ERROR")
  }
}


