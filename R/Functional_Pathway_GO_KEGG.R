# This package will autamatically install org.Hs.eg.db, org.Rn.eg.db, org.Mm.eg.db. If you want to use different OrgDb annotation species, please install that annotation package and attach to your session.
# How this function works :
# 1. If there is no Univ terms ==> NA
# 2. If there are Univ terms but no DE terms ==> Do Gene_set_analysis but not doing classification and over-representation
# 3. If there are Univ and DE terms ==> Do Gene_set_analysis and classification and over-representation
GOAnalysis <- function(which.analysis, path.prefix, OrgDb.species, go.level, input.TYPE.ID) {
  CheckGoLevel(go.level)
  message(paste0("\u2694\u2694 Gene Ontology Analysis \n"))
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/"))
  }
  # get return value
  DE_results <- DEGeneList(which.analysis, path.prefix, OrgDb.species, input.TYPE.ID)
  gene_list_input_TYP <- DE_results$gene_list_input_TYP_rt
  # gene_list_input_TYP <- gene_list_input_TYP[(gene_list_input_TYP != "Inf") & (gene_list_input_TYP != "-Inf")]
  gene_list_ENTREZID <- DE_results$gene_list_ENTREZID_rt
  # gene_list_ENTREZID <- gene_list_ENTREZID[(gene_list_ENTREZID != "Inf") & (gene_list_ENTREZID != "-Inf")]

  GO.Ontology.list <- c("MF", "BP", "CC")
  for ( i in GO.Ontology.list) {
    message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
    message(paste0("\u25CF '", i, "' : \n"))


    #####################
    #### Checking DE ####
    #####################
    message("     \u25CF Checking differential expression gene number ... \n")
    message(paste0("          \u25CF Differential expression gene number : ", length(gene_list_input_TYP), "\n"))
    # Do GO classification and GO over-representation test

    dir_name <- paste0("GO_DE_Classification_Overrepresentation")
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name))){
      dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name))
    }
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i))){
      dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i))
    }
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i, "/images"))){
      dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i, "/images"))
    }
    ###########################
    #### GO Classification ####
    ###########################
    message("     \u25CF GO Classification ... \n")
    # designed for gene classification on GO distribution at a specific level. "MF", "BP", "CC"
    ggo <- clusterProfiler::groupGO(gene     = names(gene_list_ENTREZID),
                                    keyType  = "ENTREZID",
                                    OrgDb    = OrgDb.species,   # variable
                                    ont      = i,           # variable
                                    level    = go.level)
    ggo.data.frame <- data.frame(ggo)
    ggo.data.frame <- ggo.data.frame[order(as.numeric(ggo.data.frame$GeneRatio), decreasing = TRUE),]

    ggo@result <- ggo.data.frame
    # Condition 1 for GO classification ! Row number have to bigger than 1 !
    if (length(row.names(ggo.data.frame)) > 0) {
      message(paste0("          \u25CF (\u2714) GO Classification (", i,") result found! \n"))
      message(paste0("               \u25CF Writing 'GO_", i, "_Classification.csv' \n"))
      write.csv(ggo.data.frame, file = paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i , "/GO_", i, "_Classification.csv"))
      message(paste0("               \u25CF Plotting 'GO_", i, "_Classification_Bar_Plot_clusterProfiler.png' \n"))
      png(filename = paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i, "/images/GO_", i, "_Classification_Bar_Plot_clusterProfiler.png"), width=5, height=5, units="in", res=300)
      p1 <- barplot(ggo, drop=TRUE, showCategory=12)
      print(p1)
      dev.off()
    }  else {
      message("          \u25CF (\u26A0) No term is found.\n")
      file.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i, "/GO_CLASSIFICATION_NO_TERM"))
    }

    ################################
    #### GO Over-representation ####
    ################################
    # GO over-representation test !
    message("     \u25CF GO Over-representation Test ... \n")
    # GO over-representation test
    ego <- clusterProfiler::enrichGO(gene          = names(gene_list_ENTREZID),
                                     keyType       = "ENTREZID",
                                     OrgDb         = OrgDb.species,                   # variable
                                     ont           = i,
                                     pAdjustMethod = "BH")                            # variable : "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
    ego.data.frame <- data.frame(ego)
    message(paste0("          \u25CF (\u2714) GO over-representation test (", i,") enriched term found : ", length(row.names(ego.data.frame)), "\n"))
    if (length(row.names(ego.data.frame)) > 0) {
      message(paste0("               \u25CF Writing 'GO_", i, "_Overrepresentation.csv' \n"))
      write.csv(ego.data.frame, file = paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i, "/GO_", i, "_Overrepresentation.csv"))
      # visualization
      # bar plot
      png(filename = paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i, "/images/GO_", i, "_Overrepresentation_Bar_Plot_clusterProfiler.png"), width=5, height=5, units="in", res=300)
      message(paste0("               \u25CF Plotting 'GO_", i, "_Overrepresentation_Bar_Plot_clusterProfiler.png' \n"))
      p2 <- barplot(ego, showCategory=12)
      print(p2)
      dev.off()

      # dot plot
      png(filename = paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i, "/images/GO_", i, "_Overrepresentation_Dot_Plot_clusterProfiler.png"), width=5, height=5, units="in", res=300)
      message(paste0("               \u25CF Plotting 'GO_", i, "_Overrepresentation_Dot_Plot_clusterProfiler.png' \n"))
      p3 <- clusterProfiler::dotplot(ego)
      print(p3)
      dev.off()
    } else {
      message(paste0("          \u25CF (\u26A0) No enriched term is found.\n"))
      file.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/GO_analysis/", dir_name, "/", i, "/GO_OVERREPRESENTATION_NO_TERM"))
    }
  }
}

KEGGAnalysis <- function(which.analysis, path.prefix, OrgDb.species, input.TYPE.ID, KEGG.organism) {
  message(paste0("\u2694\u2694 Kyoto Encyclopedia of Genes and Genomes Analysis \n"))
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/"))
  }
  DE_results <- DEGeneList(which.analysis, path.prefix, OrgDb.species, input.TYPE.ID)
  gene_list_input_TYP = DE_results$gene_list_input_TYP_rt
  gene_list_ENTREZID = DE_results$gene_list_ENTREZID_rt

  # Do KEGG Gene Set Enrichment Analysis
  dir_name <- paste0("KEGG_Gene_Set_Enrichment")
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name))){
    dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name))
  }

  ######################################
  #### Gene Set Enrichment Analysis ####
  ######################################
  message("\u25CF KEGG Gene Set Enrichment Analysis ... \n")
  # DO KEGG Gene Set Enrichment Analysis


  gene_list_input_TYP
  names(gene_list_input_TYP)

  eg = clusterProfiler::bitr(names(gene_list_input_TYP), fromType="GENENAME", toType="UNIPROT", OrgDb="org.Sc.sgd.db")
  indices <- match(names(gene_list_input_TYP), eg$GENENAME)
  gene_list_input_TYP <- gene_list_input_TYP[indices]
  gene_list_input_TYP <- gene_list_input_TYP[!is.na(names(gene_list_input_TYP))]
  names(gene_list_input_TYP) <- eg$UNIPROT

  kk.gse <- clusterProfiler::gseKEGG(geneList     = gene_list_input_TYP,
                                     keyType      = "uniprot",
                                     organism     = "sce")

  clusterProfiler::search_kegg_organism('sce', by='kegg_code')

  kk.gse.frame <- data.frame(kk.gse)
  if (length(row.names(kk.gse.frame)) > 0) {
    # Result column must bigger than 0
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/images"))){
      dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/images"))
    }
    message(paste0("     \u25CF (\u2714) KEGG Gene Set Enrichment test enriched term found! \n"))
    message(paste0("          \u25CF Writing 'KEGG_Gene_Set_Enrichment.csv' \n"))
    write.csv(kk.gse.frame, file = paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/KEGG_Gene_Set_Enrichment.csv"))
    message(paste0("     \u25CF Checking 'KEGG_Gene_Set_Enrichment.csv' result row number \n"))
    if (length(row.names(kk.gse.frame)) < 5) {
      message(paste0("          \u25CF 'KEGG_Gene_Set_Enrichment.csv' result row number : ", length(row.names(kk.gse.frame)), " (less than 5)\n"))
      iterate.terms <- kk.gse.frame$ID
    } else {
      message(paste0("          \u25CF 'KEGG_Gene_Set_Enrichment.csv' result row number : ", length(row.names(kk.gse.frame)), " (more than 5)\n"))
      iterate.terms <- kk.gse.frame$ID[seq_len(5)]
    }
    for(KEGG.ID in iterate.terms) {
      message(paste0("               \u25CF KEGG ID : ", KEGG.ID, "\n"))
      png(filename = paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/images/", KEGG.ID, "_gseKEGG_Plot_clusterProfiler.png"), res = 300)
      message(paste0("                    \u25CF Plotting '", KEGG.ID, "_gseKEGG_Plot_clusterProfiler.png'\n"))
      p <- clusterProfiler::gseaplot(kk.gse, geneSetID = KEGG.ID)
      print(p)
      dev.off()
    }
  } else {
    message("     \u25CF (\u26A0) No enriched term is found.\n")
    file.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/KEGG_GSE_NO_TERM"))
  }







  #####################
  #### Checking DE ####
  #####################
  # Do GO classification and GO over-representation test
  if(is.na(gene_list_input_TYP) && is.na(gene_list_ENTREZID)) {
    # Dot doing anything! Printing log will be reported by previous check!!
  } else {
    message("\u25CF Checking differential expression gene number ... \n")
    message(paste0("     \u25CF Differential expression gene number : ", length(gene_list_input_TYP), "\n"))
    dir_name <- paste0("KEGG_DE_Overrepresentation")
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name))){
      dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name))
    }
    ##################################
    #### KEGG Over-representation ####
    ##################################
    # Do KEGG over-representation test
    # organism : species supported at 'http://www.genome.jp/kegg/catalog/org_list.html'
    message("\u25CF KEGG Over-representation Test ... \n")
    # KEGG Over-representation test
    kk <- clusterProfiler::enrichKEGG(gene     = names(gene_list_input_TYP),
                                      organism = KEGG.organism)                     # variable


    kk.data.frame <- data.frame(kk)
    # Row size have to bigger than 0!
    if (length(row.names(kk.data.frame)) > 0) {
      message(paste0("     \u25CF (\u2714) KEGG over-representation test enriched term found! \n"))
      message(paste0("          \u25CF Writing 'KEGG_Overrepresentation.csv' \n"))
      write.csv(kk.data.frame, file = paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/KEGG_Overrepresentation.csv"))
      for ( i in kk.data.frame$ID) {
        if(!dir.exists(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/", i))){
          dir.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/", i))
        }
        current.path <- getwd()
        # get the url from KEGG result
        message(paste0("          \u25CF Finding '", i, "' KEGG URL ... \n"))
        KEGGUrl <- GetKEGGUrl(kk, i)
        message(paste0("               \u25CF Writting 'URL_", i, "_Pathway.txt' \n"))
        write(KEGGUrl, file = paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/", i, "/URL_", i, "_Pathway.txt"))
        # drawing pathway picture with 'pathway' package
        pathway.dir <- paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name, "/", i, "/pathview_result/")
        if(!dir.exists(pathway.dir)){
          dir.create(pathway.dir)
        }
        setwd(pathway.dir)
        message(paste0("               \u25CF Plotting '", i, "' pathway by package \"pathview\" \n"))
        pathview::pathview(gene.data  = names(gene_list_ENTREZID),
                           pathway.id = i,
                           species    = KEGG.organism,
                           kegg.dir   = pathway.dir)
        on.exit(setwd(current.path))
      }
    } else {
      message(paste0("     \u25CF (\u26A0) No over-representation term is found.\n"))
      file.create(paste0(path.prefix, "RNASeq_results/", which.analysis, "/KEGG_analysis/", dir_name,"/KEGG_OVERREPRESENTATION_NO_TERM"))
    }
  }
  message("================================================================\n================================================================\n\n")
}

DEGeneList <- function(which.analysis, path.prefix, OrgDb.species, input.TYPE.ID) {
  DE.path.csv <- paste0(path.prefix, "RNASeq_results/", which.analysis, "/", strsplit(which.analysis, "_")[[1]][1], "_normalized_DE_result.csv")
  DE.csv <- read.csv(DE.path.csv)
  # DE gene
  # First filter out "." gene name
  DE.csv <- DE.csv[DE.csv$gene.name != ".",]
  # Second sort row by the values of log2FC (from big to small)
  DE.csv <- DE.csv[order(abs(DE.csv$log2FC), decreasing = TRUE),]
  # If gene name is repeat, choose the first one!!
  DE.csv <- DE.csv[!duplicated(DE.csv$gene.name),]
  ## ==> So every gene in DE would be distinct (with the larget log2FC)
  gene_name <- as.character(DE.csv$gene.name)
  gene_list_input_TYP <- DE.csv$log2FC
  gene_list_ENTREZID <- DE.csv$log2FC


  ############
  #### DE ####
  ############
  if (length(gene_list_input_TYP) == 0 && length(gene_list_ENTREZID) == 0) {
    message(paste0("No annotated gene terms are found in '", path.prefix, "RNASeq_results/", which.analysis, "/", strsplit(which.analysis, "_")[[1]][1], "_normalized_DE_result.csv'\n\n"))
    return(list(gene_list_input_TYP_rt = NA,
                gene_list_ENTREZID_rt = NA))
  } else {
    # Rename gene_list_input_TYP
    names(gene_list_input_TYP) <-gene_name
    # Sort gene_list_input_TYP
    gene_list_input_TYP = sort(gene_list_input_TYP, decreasing = TRUE)
    # IDs conversion
    gene.df.DE <- clusterProfiler::bitr(gene_name,
                                        fromType = input.TYPE.ID,
                                        toType = c("ENTREZID"),
                                        OrgDb = OrgDb.species)

    # input.TYPE.ID
    DE.indices <- match(gene_name, gene.df.DE[input.TYPE.ID][[1]])
    selected.gene.df.DE <- gene.df.DE[DE.indices, ]
    names(gene_list_ENTREZID) <- selected.gene.df.DE$ENTREZID
    gene_list_ENTREZID <- gene_list_ENTREZID[!is.na(names(gene_list_ENTREZID))]
    gene_list_ENTREZID = sort(gene_list_ENTREZID, decreasing = TRUE)
    return(list(gene_list_input_TYP_rt = gene_list_input_TYP,
                gene_list_ENTREZID_rt = gene_list_ENTREZID))
  }
}


GetKEGGUrl <- function(x, pathID) {
  url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", pathID, '/', x[pathID, "geneID"])
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


