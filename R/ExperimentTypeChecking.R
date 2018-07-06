#'
#' Check the experiment type and seperate into three function
#'
#' @export
ExperimentTypeChecking <- function(path.prefix, gene.name, sample.pattern, experiment.type, main.variable, additional.variable) {
  # Check the how many group !
  results <- ProgressGenesFiles(path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print = FALSE)
  if (isTRUE(results$phenodata.file.df)){
    # sorting 'pheno_data'
    cat(paste0("************** Experiment group checking **************\n"))
    cat("\u25CF 1. Printing origin phenodata.csv : \n")
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    print(pheno_data)
    cat('\n')
    sample.table <- as.data.frame(table(pheno_data[main.variable]))
    groups.number <- length(row.names(sample.table))
    if (groups.number == 1) {
      stop("Group number ERROR")
    }
    if (experiment.type == "two.group") {
      if (groups.number != 2) {
        stop("Group number ERROR")
      }
      TwoGroupAnalysis(path.prefix, gene.name, sample.pattern, main.variable, additional.variable)
    } else if (experiment.type == "multi.group.pairs") {
      if (groups.number == 2) {
        stop("Group number ERROR")
      }
      MultiGroupPairsAnalysis(path.prefix, gene.name, sample.pattern, main.variable, additional.variable)
    } else if (experiment.type == "multi.group.anova") {
      if (groups.number == 2) {
        stop("Group number ERROR")
      }
      MultiGroupAnovaAnalysis(path.prefix, gene.name, sample.pattern, main.variable, additional.variable)
    }
  }
}

#' @export
TwoGroupAnalysis <- function(path.prefix, gene.name, sample.pattern, main.variable, additional.variable) {
  print(c("path.prefix : ", path.prefix))
  print(c("gene.name : ", gene.name))
  print(c("sample.pattern : ", sample.pattern))
  print(c("main.variable : ", main.variable))
  print(c("additional.variable : ", additional.variable))
  cat(paste0("Experiment type : two.group"))
  # seperate into three part analysis : FPKM, TPM, raw reads count

}

#' @export
MultiGroupPairsAnalysis <- function(path.prefix, gene.name, sample.pattern, main.variable, additional.variable) {
  print(c("path.prefix : ", path.prefix))
  print(c("gene.name : ", gene.name))
  print(c("sample.pattern : ", sample.pattern))
  print(c("main.variable : ", main.variable))
  print(c("additional.variable : ", additional.variable))
}

#' @export
MultiGroupAnovaAnalysis <- function(path.prefix, gene.name, sample.pattern, main.variable, additional.variable) {
  print(c("path.prefix : ", path.prefix))
  print(c("gene.name : ", gene.name))
  print(c("sample.pattern : ", sample.pattern))
  print(c("main.variable : ", main.variable))
  print(c("additional.variable : ", additional.variable))
}






#' Run ballgown analysis
#'
#' @import ballgown
#' @import genefilter
#' @importFrom dplyr arrange
#' @export
BallgownPreprocess <- function(path.prefix, gene.name, sample.pattern, experiment.type, main.variable, additional.variable) {
  results <- ProgressGenesFiles(path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print = FALSE)
  if (isTRUE(results$phenodata.file.df) && results$ballgown.dirs.number.df != 0){
    # sorting 'pheno_data'
    cat(paste0("************** Ballgown data preprocessing **************\n"))
    cat("\u25CF 1. Printing origin phenodata.csv : \n")
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    print(pheno_data)
    cat('\n')
    sample.table <- as.data.frame(table(pheno_data[main.variable]))
    if (length(row.names(sample.table)) == 2) {
      dir.create(paste0(path.prefix, "RNAseq_results/DEG_results/"))
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
      bg_chrX <- ballgown(dataDir = paste0(path.prefix, "gene_data/ballgown"), samplePattern = sample.pattern, pData = pheno_data, meas = 'all')
      save(bg_chrX, file = paste0(path.prefix, "gene_data/ballgown/ballgown.rda"))
      cat('\n')
      bg_chrX_filt <- ballgown::subset(bg_chrX, cond = 'rowVars(ballgown::texpr(pkg.ballgown.data$bg_chrX)) >1', genomesubset=TRUE)

      # differential expression
      cat("\u25CF 4. Differential expression preprocessing : \n")
      cat("     \u25CF creating 'Differential Expression Gene FPKM data.frame ......'\n")
      cat(c("         \u25CF  main.variable :", main.variable, "\n"))
      if (length(adjustvars) != 0) {
        cat(c("         \u25CF adjustvars :", adjustvars, "\n"))
        results_transcripts <- stattest(pkg.ballgown.data$bg_chrX_filt, feature="transcript",main.variable=main.variable, adjustvars = adjustvars, getFC=TRUE, meas="FPKM")
      } else {
        results_transcripts <- stattest(pkg.ballgown.data$bg_chrX_filt, feature="transcript",main.variable=main.variable, getFC=TRUE, meas="FPKM")
      }
      results_transcripts$feature <- NULL
      results_transcripts.FC <- results_transcripts$fc
      results_transcripts.log2FC <- log2(results_transcripts$fc)
      results_transcripts.pval <- results_transcripts$pval
      results_transcripts.qval <- results_transcripts$qval
      results_transcripts$fc <- NULL; results_transcripts$pval <- NULL; results_transcripts$qval <- NULL; colnames(results_transcripts)[1] <- "transcriptIDs"
      results_transcripts <- data.frame(geneNames=ballgown::geneNames(pkg.ballgown.data$bg_chrX_filt), geneIDs=ballgown::geneIDs(pkg.ballgown.data$bg_chrX_filt), transcriptNames=transcriptNames(pkg.ballgown.data$bg_chrX_filt), results_transcripts)
      # adding fpkm
      # cov : average per-base read coverage
      cat("     \u25CF merging each FPKM column and calculating average FPKM ......'\n")
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
      cat("     \u25CF writing data.frame into 'FPKM_DEG_result.csv' ......'\n\n")
      write.csv(results_transcripts, paste0(path.prefix, "RNAseq_results/DEG_results/FPKM_DEG_result.csv"), row.names=FALSE)
      cat("\u25CF 5. Printing DEG dataset : \n")
      print(head(results_transcripts))
      cat("\n")
    } else {
      stop(paste0("(\u2718) This pipline is only available for 2-group comparisons.\n", "      ",length(row.names(sample.table)), "-group is detected.\n\n"))
    }
  }
}

#' Check ballgown object
#' @export
CheckBallgownObject <- function() {
  print(pkg.ballgown.data$bg_chrX)
  print(pkg.ballgown.data$bg_chrX_filt)
}

#' load ballgown object
#' @export
LoadBallgownObject <- function() {
  if(isTRUE(file.exists(paste0(path.prefix, "gene_data/ballgown/ballgown.rda")))) {
    load(paste0(path.prefix, "gene_data/ballgown/ballgown.rda"))
    pkg.ballgown.data$bg_chrX <- bg
  } else {
    stop(paste0("(\u2718) '", paste0(path.prefix, "gene_data/ballgown/ballgown.rda"), "' haven't created yet. Please run \"BallgownPreprocess(gene.name, sample.pattern, main.variable)\" first.\n\n"))
  }
}
