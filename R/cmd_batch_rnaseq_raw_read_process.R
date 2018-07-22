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
RNAseqRawReadProcess_CMD <- function(RNASeqWorkFlowParam, num.parallel.threads = 8, run = TRUE, check.s4.print = TRUE) {
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  gene.name <- RNASeqWorkFlowParam@gene.name
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  python.variable <- RNASeqWorkFlowParam@python.variable
  python.variable.answer <- python.variable$check.answer
  python.variable.version <- python.variable$python.version
  indexes.optional <- RNASeqWorkFlowParam@indexes.optional
  # not print but if the prefix is invalid, then 'Prefix path '", path.prefix, "' is invalid. Please try another one.' will be printed.
  # If precheck doesn't have .ht2 files is fine
  # ExportPath(path.prefix = path.prefix)
  fileConn<-file(paste0(path.prefix, "Rscript/Raw_Read_Process.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0('RNAseqRawReadProcess(path.prefix = "', path.prefix, '", input.path.prefix = "', input.path.prefix, '", gene.name = "', gene.name, '", sample.pattern = "', sample.pattern, '", python.variable.answer = ', python.variable.answer, ', python.variable.version = ', python.variable.version, ', num.parallel.threads = ', num.parallel.threads, ', indexes.optional = ', indexes.optional, ')')
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Raw_Read_Process.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Raw_Read_Process.R ", path.prefix, "Rscript_out/Raw_Read_Process.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 RNAseq alignment, assembly, quantification, mergence, comparison, reads process are doing in the background. Check current progress in '", path.prefix, "Rscript_out/Raw_Read_Process.Rout'\n\n"))
  }
}


#' rna seq pipline
#' @export
RNAseqRawReadProcess <- function(path.prefix, input.path.prefix, gene.name, sample.pattern, python.variable.answer, python.variable.version, num.parallel.threads = 8, indexes.optional) {
  ExportPath(path.prefix)
  PreRNAseqRawReadProcess(path.prefix = path.prefix, gene.name = gene.name, sample.pattern = sample.pattern)
  check.results <- ProgressGenesFiles(path.prefix = path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print=FALSE)
  if (check.results$ht2.files.number.df == 0 && indexes.optional) {
    CreateHisat2Index(gene.name, sample.pattern)
  }
  Hisat2AlignmentDefault(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  SamtoolsToBam(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  StringTieAssemble(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  StringTieMergeTrans(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  GffcompareRefSample(path.prefix, gene.name, sample.pattern)
  StringTieToBallgown(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  # StringTieReEstimate(path.prefix, gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  finals <- ProgressGenesFiles(path.prefix, gene.name, sample.pattern, print=TRUE)
  PreDECountTable(path.prefix= path.prefix, sample.pattern = sample.pattern, python.variable.answer = python.variable.answer, python.variable.version = python.variable.version, print=TRUE)
  file.prepDE.py <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/prepDE.py"))
  file.sample.lst.txt <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/sample_lst.txt"))
  file.gene_count_matrix <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/gene_count_matrix.csv"))
  file.transcript_count_matrix <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/transcript_count_matrix.csv"))
  PostRNAseqRawReadProcess(path.prefix = path.prefix, sample.pattern = sample.pattern)
}

PreRNAseqRawReadProcess <- function(path.prefix, gene.name, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNAseqRawReadProcess()' environment pre-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  chrX.gtf <- file.exists(paste0(path.prefix, "gene_data/ref_genes/chrX.gtf"))
  chrX.fa <- file.exists(paste0(path.prefix, "gene_data/ref_genome/chrX.fa"))
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  check.tool.result <- CheckToolAll()
  check.results <- ProgressGenesFiles(path.prefix = path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print=FALSE)
  check.progress.results.bool <- check.results$gtf.file.logic.df && check.results$fa.file.logic.df && (check.results$fastq.gz.files.number.df != 0)
  validity <- phenodata.csv && chrX.gtf && chrX.fa && check.tool.result && (length(raw.fastq) != 0) && check.progress.results.bool
  if (!isTRUE(validity)) {
    stop("RNAseqRawReadProcess() environment ERROR")
  }
  cat("     (\u2714) : RNAseqRawReadProcess() pre-check is valid\n\n")
}

PostRNAseqRawReadProcess <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNAseqRawReadProcess()' environment post-check ...\n")
  # Still need to add condition
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNAseqRawReadProcess() post-check ERROR")
  }
  cat("     (\u2714) : RNAseqRawReadProcess() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}

