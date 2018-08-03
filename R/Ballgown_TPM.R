Ballgown_TPM <- function() {
  cat("\u25CF 1. Printing origin phenodata.csv : \n")
  pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
  print(pheno_data)
  cat('\n')
  sample.table <- as.data.frame(table(pheno_data[independent.variable]))
  if (length(row.names(sample.table)) == 2) {
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/"))
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Ballgown_object/"))
    dir.create(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/images/"))

    cat("\u25CF 2. Sorting phenodata.csv : \n")
    pheno_data.arrange <- dplyr::arrange(pheno_data, unlist(pheno_data[independent.variable]))
    print(pheno_data.arrange)
    cat('\n')
    # for adding FPKM column!
    sample.names <- as.character(pheno_data.arrange$ids)
    sample.names.with.independent.variable <- paste0(pheno_data.arrange$ids, ".", pheno_data.arrange[independent.variable][,1])
    sample.number <- length(sample.names)


    # creat trancript gtf
    for (i in sample.names) {
      input.file.path <- paste0(path.prefix, "gene_data/ballgown/", i, "/", i, ".gtf")
      output.file.path <- paste0(path.prefix, "gene_data/ballgown/", i, "/", i, "_transcript.gtf")
      command <- paste0("'{if ($3==\"transcript\") print}' ", input.file.path, " | cut -f 1,4,9 > ", output.file.path)
      print(command)
      main.command <- "awk"
      command.result <- system2(command = main.command, args = c(command))
      if (command.result != 0 ) {
        cat(paste0("(\u2718) '", main.command, "' is failed !!"))
        stop(paste0("'", main.command, "' ERROR"))
      }
    }


    my_file <- "ERR204916_transcript.gtf"
    my_obj <- rtracklayer::import(my_file)

    input.file.path <- paste0(path.prefix, "gene_data/ballgown/", "ERR204916", "/")
    setwd(input.file.path)

    read.files <- read.delim(my_file)
    ens <- refGenome::ensemblGenome()
    gtf <- refGenome::read.gtf(ens, "ERR204916_transcript.gtf")

    refGenome::tableSeqids(ens)
    refGenome::extractSeqids(ens, 'chrX')
    head(refGenome::tableFeatures(ens))

    ballgown::structure(bg)

    ballgown::texpr(bg, 'FPKM')

    # make ballgown object
    cat("\u25CF 3. Making ballgown object : \n")
    pkg.ballgown.data$bg_chrX <- ballgown::ballgown(dataDir = paste0(path.prefix, "gene_data/ballgown"), samplePattern = sample.pattern, pData = pheno_data, meas = 'all')
    bg <- pkg.ballgown.data$bg_chrX
    save(bg, file = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Ballgown_object/ballgown.rda"))
    cat("     \u25CF writing data.frame into 'ballgown.rda' in \n")
    cat('\n')
    cat("\u25CF 4. Filtering ballgown object (variance less than 1): \n")
    pkg.ballgown.data$bg_chrX_filt <- ballgown::subset(pkg.ballgown.data$bg_chrX, cond = 'genefilter::rowVars(ballgown::texpr(pkg.ballgown.data$bg_chrX)) >1', genomesubset=TRUE)
    bg_filter <- pkg.ballgown.data$bg_chrX_filt
    save(bg_filter, file = paste0(path.prefix, "RNAseq_results/Ballgown_analysis/Ballgown_object/ballgown_filter.rda"))
    cat("     \u25CF writing data.frame into 'ballgown_filter.rda' in \n")
    cat('\n')
    # differential expression
    cat("\u25CF 5. Differential Transcript expression preprocessing : \n")
    cat("     \u25CF creating 'Ballgown Transcript FPKM data.frame ......\n")
    cat(c("         \u25CF  independent.variable :", independent.variable, "\n"))

    results_transcripts <- ballgown::stattest(pkg.ballgown.data$bg_chrX_filt, feature="transcript",covariate = independent.variable, getFC=TRUE, meas="FPKM")
    results_transcripts_2 <- ballgown::stattest(pkg.ballgown.data$bg_chrX, feature="transcript",covariate = independent.variable, getFC=TRUE, meas="FPKM")
    # write.csv(results_transcripts, paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.csv"), row.names=FALSE)
    results_transcripts$feature <- NULL
    results_transcripts.FC <- results_transcripts$fc
    results_transcripts.log2FC <- log2(results_transcripts$fc)
    results_transcripts.pval <- results_transcripts$pval
    results_transcripts.qval <- results_transcripts$qval
    results_transcripts$fc <- NULL; results_transcripts$pval <- NULL; results_transcripts$qval <- NULL; colnames(results_transcripts)[1] <- "transcriptIDs"
    results_transcripts <- data.frame(geneNames=ballgown::geneNames(pkg.ballgown.data$bg_chrX_filt), geneIDs=ballgown::geneIDs(pkg.ballgown.data$bg_chrX_filt), transcriptNames=ballgown::transcriptNames(pkg.ballgown.data$bg_chrX_filt), results_transcripts)
    # adding fpkm
    # cov : average per-base read coverage
    cat("     \u25CF merging each FPKM column and calculating average FPKM ......\n")
    fpkm <- data.frame(ballgown::texpr(pkg.ballgown.data$bg_chrX_filt,meas="FPKM"))
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
        results_transcripts[[sample.names.with.independent.variable[current.sum+j]]] <- fpkm[[a]]
        columns.to.mean <- append(columns.to.mean, sample.names.with.independent.variable[current.sum+j])
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
    cat("\u25CF 6. Printing Ballgown FPKM dataset : \n")
    print(head(results_transcripts))
    cat("\n")
    cat("\u25CF 7. Creating Ballgown FPKM dataset png : \n")
    png(paste0(path.prefix, "RNAseq_results/Ballgown_analysis/ballgown_FPKM_result.png"), width = sample.number*200 + 200, height = 40*8)
    p <- gridExtra::grid.table(head(results_transcripts))
    print(p)
    dev.off()
  }
}
