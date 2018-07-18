#' Run RNA seq raw read processing (alignment, assembly, expressive expression) in background.
#'
#' To process RNA seq raw read in background, you have to run 'RNASeqWorkFlowParam()' to create S4 object beforehands.
#' If you want to process RNA seq raw read in R shell, please see 'RNAseqRawReadProcess()' function.
#' This function do 4 things : 1. Hisat2 : raw reads align to reference genome. 2. Stringtie : Alignments assembly into transcript. 3. Stringtie : Expressive gene and transcript quantification
#' 4. Compare aligned genome with reference genome
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related parameters
#'
#' @export
#' @example
#' exp <- RNASeqWorkFlowParam(path.prefix = "/home/rnaseq", input.path.prefix = "/home", gene.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNAseqEnvironmentSet_CMD(RNASeqWorkFlowParam <- exp)
RNAseqRawReadProcess_CMD <- function(RNASeqWorkFlowParam, num.parallel.threads = 8) {
  CheckS4Object(RNASeqWorkFlowParam)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  gene.name <- RNASeqWorkFlowParam@gene.name
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  python.variable <- RNASeqWorkFlowParam@python.variable
  indexes.optional <- RNASeqWorkFlowParam@indexes.optional
  # not print but if the prefix is invalid, then 'Prefix path '", path.prefix, "' is invalid. Please try another one.' will be printed.
  results <- ProgressGenesFiles(path.prefix = path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print=TRUE)
  if (isTRUE(results$gtf.file.logic.df) && isTRUE(results$fa.file.logic.df) && results$fastq.gz.files.number.df != 0 && isTRUE(results$phenodata.file.df)) {
    # If precheck doesn't have .ht2 files is fine
    ExportPath(path.prefix = path.prefix)
    if (isTRUE(CheckToolAll(print=TRUE))) {
      cat("(\u2714) : Successful in RNAseq-pipeline precheck. \n\n")
      fileConn<-file(paste0(path.prefix, "Rscript/Raw_Read_Process.R"))
      first <- "library(RNASeqWorkflow)"
      second <- paste0('RNAseqRawReadProcess(path.prefix = "', path.prefix, '", input.path.prefix = "', input.path.prefix, '", gene.name = "', gene.name, '", sample.pattern = "', sample.pattern, '", python.variable.answer = ', python.variable.answer, ', python.variable.version = ', python.variable.version, ', num.parallel.threads = ', num.parallel.threads, ', indexes.optional = ', indexes.optional, ')')
      writeLines(c(first, second), fileConn)
      close(fileConn)
      cat(paste0("\u2605 '", path.prefix, "Rscript/Raw_Read_Process.R' has been created.\n"))
      system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Raw_Read_Process.R ", path.prefix, "Rscript_out/Raw_Read_Process.Rout"), stdout = "", wait = FALSE)
      cat(paste0("\u2605 RNAseq alignment, assembly, quantification, mergence, comparison, reads process are doing in the background. Check current progress in '", path.prefix, "Rscript_out/Raw_Read_Process.Rout'\n\n"))
    }
  }
}


#' rna seq pipline
#' @export
RNAseqRawReadProcess <- function(path.prefix, input.path.prefix, gene.name, sample.pattern, python.variable.answer, python.variable.version, num.parallel.threads = 8, indexes.optional) {
  Check.tools.result <- CheckToolAll()
  # check condition need to fix
  if (Check.tools.result) {
    ExportPath(path.prefix)
    check.results <- ProgressGenesFiles(path.prefix = path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print=FALSE)
    if (isTRUE(check.results$gtf.file.logic.df) && isTRUE(check.results$fa.file.logic.df) && (check.results$fastq.gz.files.number.df != 0)) {
      if (check.results$ht2.files.number.df == 0 && indexes.optional) {
        CreateHisat2Index(gene.name, sample.pattern)
      }
      Hisat2AlignmentDefault(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
      SamtoolsToBam(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
      StringTieAssemble(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
      StringTieMergeTrans(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
      GffcompareRefSample(path.prefix, gene.name, sample.pattern)
      StringTieToBallgown(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
      finals <- ProgressGenesFiles(path.prefix, gene.name, sample.pattern, print=TRUE)
      PreDECountTable(path.prefix= path.prefix, sample.pattern = sample.pattern, python.variable.answer = python.variable.answer, python.variable.version = python.variable.version, print=TRUE)
      file.prepDE.py <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/prepDE.py"))
      file.sample.lst.txt <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/sample_lst.txt"))
      file.gene_count_matrix <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/gene_count_matrix.csv"))
      file.transcript_count_matrix <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/transcript_count_matrix.csv"))
      if (isTRUE(finals$gtf.file.logic.df) && isTRUE(finals$fa.file.logic.df) &&
          finals$fastq.gz.files.number.df != 0 &&
          isTRUE(finals$phenodata.file.df) &&
          finals$ht2.files.number.df != 0 &&
          finals$sam.files.number.df != 0 &&
          finals$bam.files.number.df != 0 &&
          finals$gtf.files.number.df != 0 &&
          isTRUE(finals$stringtie_merged.gtf.file.df) &&
          finals$gffcompare.related.dirs.number.df != 0 &&
          finals$ballgown.dirs.number.df != 0 &&
          file.prepDE.py && file.sample.lst.txt && file.gene_count_matrix && file.transcript_count_matrix) {
        cat("\n")
        cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
        cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
        cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
      } else {
        stop(paste0("(\u2718) Output files ERROR \n"))
      }
    }
    else{
      stop(paste0("(\u2718) Necessary files are lost.\n     Please check whether 'ref_genes/", gene.name, ".gtf' , 'ref_genome/", gene.name, ".fa' , 'samples_.fastq.gz/XXX_", gene.name, "_*.fastq.gz' are exit.\n\n" ))
    }
  } else {
    stop("RNAseqRawReadProcess enviroment ERROR")
  }
}
