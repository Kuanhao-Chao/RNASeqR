#' converting stringtie ballogwn preprocessed data to count table
#'
#' @importFrom reticulate py_config py_available
#' @export
PreDECountTable <- function(path.prefix, sample.pattern, print=TRUE) {
  # ftp server : ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-source.zip
  if (!dir.exists(paste0(path.prefix, "gene_data/ballgown/raw_count/"))){
    dir.create(paste0(path.prefix, "gene_data/ballgown/raw_count/"))
  }
  cat("************** Installing prepDE.py ************\n")
  cat(paste0(path.prefix, "gene_data/ballgown/raw_count\n"))
  current.path <- getwd()
  setwd(paste0(path.prefix, "gene_data/ballgown/raw_count/"))
  system2(command = 'curl', args = c('https://ccb.jhu.edu/software/stringtie/dl/prepDE.py', '--output', paste0(path.prefix, "gene_data/ballgown/raw_count/prepDE.py")), stdout = "", wait = TRUE)
  cat(paste0("'", path.prefix, "gene_data/ballgown/raw_count/prepDE.py' has been installed.\n\n"))
  cat("************** Creating 'sample_lst.txt' file ************\n")
  sample.files <- list.files(paste0(path.prefix, "gene_data/ballgown/"), pattern = sample.pattern)
  write.content <- print(paste0(sample.files[1], " ", path.prefix, "gene_data/ballgown/", sample.files[1] ,"/", sample.files[1], ".gtf"))
  for(i in 2:length(sample.files)){
    #print(paste0(sample.files[i], " ", path.prefix, "gene_data/ballgown/", sample.files[i],"/", sample.files[i], ".gtf"))
    write.content <- c(write.content, paste0(sample.files[i], " ", path.prefix, "gene_data/ballgown/", sample.files[i],"/", sample.files[i], ".gtf"))
  }
  write.file<-file(paste0(path.prefix, "gene_data/ballgown/raw_count/sample_lst.txt"))
  writeLines(write.content, write.file)
  close(write.file)
  cat(paste0("'", path.prefix, "gene_data/ballgown/raw_count/sample_lst.txt' has been created\n\n"))
  cat("************** Creating gene and transcript raw count file ************\n")
  #print(paste0(path.prefix, "gene_data/ballgown/raw_count/prepDE.py -i ",  path.prefix, "gene_data/ballgown/raw_count/sample_lst.txt"))
  # have to check python !!!
  if(reticulate::py_available(initialize = "TRUE")){
    cat("(\u2714) : Python is available on your device!\n")
    python.version <- as.numeric(reticulate::py_config()$version)
    cat(paste0("       Python version : ", reticulate::py_config()$version, "\n"))
    if(python.version >= 3) {
      cat("(\u270D) : Converting 'prepDE.py' from python2 to python3 \n\n")
      system2(command = '2to3', arg = paste0("-w ", path.prefix, "gene_data/ballgown/raw_count/prepDE.py"))
    } else if (python.version < 3 && python.version >= 2 ){
    }
    system2(command = 'python', args = paste0(path.prefix, "gene_data/ballgown/raw_count/prepDE.py -i ",  path.prefix, "gene_data/ballgown/raw_count/sample_lst.txt"))
    cat(paste0("'", path.prefix, "gene_data/ballgown/raw_count/gene_count_matrix.csv' has been created\n"))
    cat(paste0("'", path.prefix, "gene_data/ballgown/raw_count/transcript_count_matrix.csv' has been created\n\n"))
    on.exit(setwd(current.path))
    return(TRUE)
  } else {
    on.exit(setwd(current.path))
    stop("(\u2718)  Python is not available on this device. Please install python to run python script 'prepDE.py'\n\n")
  }
}

#' DEG analysis with edgeR
#'
#' @import edgeR
#' @export
DEGedgeRPlot <- function(path.prefix) {
  if(file.exists(paste0(path.prefix, "gene_data/ballgown/raw_count/gene_count_matrix.csv"))){
    # load gene name for further usage
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/DEG_results/edgeR"))){
      dir.create(paste0(path.prefix, "RNAseq_results/DEG_results/edgeR"))
    }
    cat(paste0("************** Plotting MDS plot (edgeR) **************\n"))
    # likelihood ratio test and quasi-likelihood F-test
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    count.table <- read.csv(paste0(path.prefix, "gene_data/ballgown/raw_count/gene_count_matrix.csv"))
    group <- pheno_data$sex
    gene.data.frame <- data.frame(gene.id=count.table$gene_id)

    # create DGEList object (edgeR)
    deglist.object <- edgeR::DGEList(counts=count.table[-1], group = group, genes = gene.data.frame)
    # Normalization
    deglist.object <- edgeR::calcNormFactors(deglist.object, method="TMM")
    png(paste0(path.prefix, "RNAseq_results/DEG_results/edgeR/MDS_plot.png"))
    my_colors=c(rgb(255, 47, 35,maxColorValue = 255),
                rgb(50, 147, 255,maxColorValue = 255))

    plotMDS(deglist.object, top = 1000, labels = NULL, col = my_colors[as.numeric(deglist.object$samples$group)],
            pch = 20, cex = 2)
    par(xpd=TRUE)
    legend("bottomright",inset=c(0,1), horiz=TRUE, bty="n", legend=levels(deglist.object$samples$group) , col=my_colors, pch=20 )
    dev.off()

    countsPerMillion <- edgeR::cpm(deglist.object)
    logcountsPerMillion <- edgeR::cpm(deglist.object, log=TRUE)
    # countCheck <- countsPerMillion > 1

    # estimating Dispersions
    sampleType<- as.character(group)
    # sampleType[grep("female", sampleType)] <- "F"
    # sampleType[grep("male", sampleType)] <- "M"
    # set up model
    designMat <- model.matrix(~sampleType)

    dgList <- edgeR::estimateGLMCommonDisp(deglist.object, design=designMat)
    dgList <- edgeR::estimateGLMTrendedDisp(dgList, design=designMat)
    dgList <- edgeR::estimateGLMTagwiseDisp(dgList, design=designMat)
    cat(paste0("************** Plotting BCV (Biological Coefficient Of Variation) plot (edgeR) **************\n"))
    png(paste0(path.prefix, "RNAseq_results/DEG_results/edgeR/BCV_plot.png"))
    p <- edgeR::plotBCV(dgList)
    print(p)
    dev.off()

    # plot smear plot
    fit <- edgeR::glmFit(dgList, designMat)
    lrt <- edgeR::glmLRT(fit, coef=2)
    edgeR_result <- edgeR::topTags(lrt)
    deGenes <- edgeR::decideTestsDGE(lrt, p=0.001)
    deGenes <- rownames(lrt)[as.logical(deGenes)]
    cat(paste0("************** Plotting smear plot (edgeR) **************\n"))
    png(paste0(path.prefix, "RNAseq_results/DEG_results/edgeR/Smear_plot.png"))
    p <- edgeR::plotSmear(lrt, de.tags=deGenes)
    print(p)
    abline(h=c(-1, 1), col=2)
    dev.off()
  } else {
    stop("'gene_count_matrix.csv' not exist ERROR")
  }
}

#' DEG analysis with DESeq2
#'
#' @import DESeq2
#' @importFrom pheatmap pheatmap
#' @export
DEDESeq2Plot <- function() {
  if(file.exists(paste0(path.prefix, "gene_data/ballgown/raw_count/gene_count_matrix.csv"))){
    # load gene name for further usage
    if(!dir.exists(paste0(path.prefix, "RNAseq_results/DEG_results/DESeq2"))){
      dir.create(paste0(path.prefix, "RNAseq_results/DEG_results/DESeq2"))
    }
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    count.table <- read.csv(paste0(path.prefix, "gene_data/ballgown/raw_count/gene_count_matrix.csv"))
    group <- pheno_data$sex
    gene.data.frame <- data.frame(gene.id=count.table$gene_id)

    # Deseq data
    colData <- as.data.frame(as.character(group))
    colnames(colData) = c("covariate")
    ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = count.table[-1],
                                             colData = colData,
                                             design = ~covariate)

    # transformation
    vsd <- DESeq2::varianceStabilizingTransformation(ddsMat)
    DESeq2::plotPCA(vsd, "covariate")

    # plot pca by ggplot2
    # data <- plotPCA(vsd, intgroup = c( "covariate"), returnData=TRUE)
    # percentVar <- round(100 * attr(data, "percentVar"))
    # ggplot(data, aes(PC1, PC2, color=covariate)) + geom_point(size=3) +
    # xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    # ylab(paste0("PC2: ",percentVar[2],"% variance"))

    ## differential
    ddsMat <- DESeq2::DESeq(ddsMat)
    # building result table
    res <- DESeq2::results(ddsMat)
    mcols(res, use.names=TRUE)

    res.05 <- DESeq2::results(ddsMat, alpha=.05)
    table(res.05$padj < .05)

    resLFC1 <- DESeq2::results(ddsMat, lfcThreshold=1)
    table(resLFC1$padj < 0.1)
    # MA plot with DESeq2
    cat(paste0("************** Plotting MA plot (DESeq2) **************\n"))
    png(paste0(path.prefix, "RNAseq_results/DEG_results/DESeq2/MA_plot.png"))
    p <- DESeq2::plotMA(res, ylim=c(-5,5))
    print(p)
    dev.off()

    mat <- assay(vsd)[ head(order(res$padj),30), ]
    mat <- mat - rowMeans(mat)
    df <- as.data.frame(colData(vsd)[,c("covariate")])
    rownames(df) <- as.character(pheno_data$ids)
    cat(paste0("************** Plotting heatmap plot (DESeq2) **************\n"))
    png(paste0(path.prefix, "RNAseq_results/DEG_results/DESeq2/Heatmap_plot.png"))
    pheatmap::pheatmap(mat, annotation_col=df)
    dev.off()
    #
    #     table(DESeq2=res$padj < 0.1, edgeR=tt.all$table$FDR < 0.1)
    #
    #     treatres <- glmTreat(fit, coef = ncol(designMat), lfc = 1)
    #     tt.treat <- topTags(treatres, n = nrow(y), sort.by = "none")
    #     table(DESeq2 = resLFC1$padj < 0.1, edgeR = tt.treat$table$FDR < 0.1)
  }
}
