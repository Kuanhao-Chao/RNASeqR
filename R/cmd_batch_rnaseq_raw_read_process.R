#' exp <- RNASeqWorkFlowParam(path.prefix = "/home/rnaseq", input.path.prefix = "/home", genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNAseqEnvironmentSet_CMD(RNASeqWorkFlowParam <- exp)
#'
#' @title Raw reads process (alignment, assembly, expressive quantification) of RNA-Seq in background.
#'
#' @description Process raw reads for RNA-Seq workflow in background.
#' This function do 5 things :
#' 1. 'Hisat2' : aligns raw reads to reference genome. If \code{indexes.optional} in \code{RNASeqWorkFlowParam} is \code{FALSE}, Hisat2 indexes will be created.
#' 2. 'Samtools' : converts '.sam' files to '.bam' files.
#' 3. 'Stringtie' : assembles alignments into transcript.
#' 4. 'Gffcompare' : examines how transcripts compare with the reference annotation.
#' 5. 'Stringtie' : creates input files for ballgown, edgeR and DESeq2.
#' Before running this function, \code{RNAseqEnvironmentSet_CMD()} (or\code{RNAseqEnvironmentSet()}) must be executed successfully.
#' If you want to process raw reads for the following RNA-Seq workflow in R shell, please see \code{RNAseqRawReadProcess()} function.
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related parameters
#' @param num.parallel.threads Specify the number of processing threads (CPUs) to use for transcript assembly. The default is 1.
#' @param run Default value is \code{TRUE}. If \code{TRUE}, 'Rscript/Environment_Set.R' will be created and executed. The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#' If \code{False}, 'Rscript/Environment_Set.R' will be created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of checking \code{RNASeqWorkFlowParam} will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking \code{RNASeqWorkFlowParam} will not be in 'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' ## Before run this function, make sure \code{RNAseqEnvironmentSet_CMD()} (or\code{RNAseqEnvironmentSet()}) is executed successfully.
#' RNAseqRawReadProcess_CMD(RNASeqWorkFlowParam = exp)}
RNAseqRawReadProcess_CMD <- function(RNASeqWorkFlowParam, num.parallel.threads = 1, run = TRUE, check.s4.print = TRUE) {
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  genome.name <- RNASeqWorkFlowParam@genome.name
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
  second <- paste0('RNAseqRawReadProcess(path.prefix = "', path.prefix, '", input.path.prefix = "', input.path.prefix, '", genome.name = "', genome.name, '", sample.pattern = "', sample.pattern, '", python.variable.answer = ', python.variable.answer, ', python.variable.version = ', python.variable.version, ', num.parallel.threads = ', num.parallel.threads, ', indexes.optional = ', indexes.optional, ')')
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Raw_Read_Process.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Raw_Read_Process.R ", path.prefix, "Rscript_out/Raw_Read_Process.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 RNAseq alignment, assembly, quantification, mergence, comparison, reads process are doing in the background. Check current progress in '", path.prefix, "Rscript_out/Raw_Read_Process.Rout'\n\n"))
  }
}

#' exp <- RNASeqWorkFlowParam(path.prefix = "/home/rnaseq", input.path.prefix = "/home", genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNAseqEnvironmentSet_CMD(RNASeqWorkFlowParam <- exp)
#'
#' @title Raw reads process (alignment, assembly, expressive quantification) of RNA-Seq in background.
#'
#' @description Process raw reads for RNA-Seq workflow in R shell.
#' It is strongly advised to run \code{RNAseqRawReadProcess_CMD()} directly. Running this function directly is not recommended.
#' This function do 5 things :
#' 1. 'Hisat2' : aligns raw reads to reference genome. If \code{indexes.optional} in \code{RNASeqWorkFlowParam} is \code{FALSE}, Hisat2 indexes will be created.
#' 2. 'Samtools' : converts '.sam' files to '.bam' files.
#' 3. 'Stringtie' : assembles alignments into transcript.
#' 4. 'Gffcompare' : examines how transcripts compare with the reference annotation.
#' 5. 'Stringtie' : creates input files for ballgown, edgeR and DESeq2.
#' Before running this function, \code{RNAseqEnvironmentSet_CMD()} (or\code{RNAseqEnvironmentSet()}) must be executed successfully.
#' If you want to process raw reads for the following RNA-Seq workflow in background, please see \code{RNAseqRawReadProcess_CMD()} function.
#'
#' @param path.prefix path prefix of 'gene_data/', 'RNAseq_bin/', 'RNAseq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @param input.path.prefix path prefix of 'input_files/' directory
#' @param genome.name variable of genome name defined in this RNA-Seq workflow (ex. genome.name.fa, genome.name.gtf)
#' @param sample.pattern  regular expression of raw fastq.gz files under 'input_files/raw_fastq.gz'
#' @param python.variable.answer logical value whether python is available on the device
#' @param python.variable.version python version of the device
#' @param indexes.optional logical value whether 'indexes/' is exit in 'input_files/'
#' @param num.parallel.threads Specify the number of processing threads (CPUs) to use for transcript assembly. The default is 1.
#' @param run Default value is \code{TRUE}. If \code{TRUE}, 'Rscript/Environment_Set.R' will be created and executed. The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#' If \code{False}, 'Rscript/Environment_Set.R' will be created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of checking \code{RNASeqWorkFlowParam} will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking \code{RNASeqWorkFlowParam} will not be in 'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' ## Before run this function, make sure \code{RNAseqEnvironmentSet_CMD()} (or\code{RNAseqEnvironmentSet()}) is executed successfully.
#' RNAseqRawReadProcess(path.prefix = exp@@path.prefix, input.path.prefix = exp@@input.path.prefix,
#'                      genome.name = exp@@genome.name, sample.pattern = exp@@sample.pattern, python.variable.answer = exp@@python.variable[0],
#'                      python.variable.version = exp@@python.variable[1], indexes.optional = exp@@indexes.optional)}
RNAseqRawReadProcess <- function(path.prefix, input.path.prefix, genome.name, sample.pattern, python.variable.answer, python.variable.version, num.parallel.threads = 1, indexes.optional) {
  CheckOperatingSystem(FALSE)
  ExportPath(path.prefix)
  PreRNAseqRawReadProcess(path.prefix = path.prefix, genome.name = genome.name, sample.pattern = sample.pattern)
  check.results <- ProgressGenesFiles(path.prefix = path.prefix, genome.name = genome.name, sample.pattern = sample.pattern, print=FALSE)
  if (check.results$ht2.files.number.df == 0 && indexes.optional) {
    CreateHisat2Index(genome.name, sample.pattern)
  }
  Hisat2AlignmentDefault(path.prefix, genome.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  SamtoolsToBam(path.prefix, genome.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  StringTieAssemble(path.prefix, genome.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  StringTieMergeTrans(path.prefix, genome.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  GffcompareRefSample(path.prefix, genome.name, sample.pattern)
  StringTieToBallgown(path.prefix, genome.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  # StringTieReEstimate(path.prefix, genome.name, sample.pattern, num.parallel.threads = num.parallel.threads)
  finals <- ProgressGenesFiles(path.prefix, genome.name, sample.pattern, print=TRUE)
  PreDECountTable(path.prefix= path.prefix, sample.pattern = sample.pattern, python.variable.answer = python.variable.answer, python.variable.version = python.variable.version, print=TRUE)
  file.prepDE.py <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/prepDE.py"))
  file.sample.lst.txt <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/sample_lst.txt"))
  file.gene_count_matrix <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/gene_count_matrix.csv"))
  file.transcript_count_matrix <- file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/transcript_count_matrix.csv"))
  PostRNAseqRawReadProcess(path.prefix = path.prefix, genome.name = genome.name, sample.pattern = sample.pattern)
}

PreRNAseqRawReadProcess <- function(path.prefix, genome.name, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNAseqRawReadProcess()' environment pre-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  chrX.gtf <- file.exists(paste0(path.prefix, "gene_data/ref_genes/chrX.gtf"))
  chrX.fa <- file.exists(paste0(path.prefix, "gene_data/ref_genome/chrX.fa"))
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  check.tool.result <- CheckToolAll()
  check.results <- ProgressGenesFiles(path.prefix = path.prefix, genome.name = genome.name, sample.pattern = sample.pattern, print=FALSE)
  check.progress.results.bool <- check.results$gtf.file.logic.df && check.results$fa.file.logic.df && (check.results$fastq.gz.files.number.df != 0)
  validity <- phenodata.csv && chrX.gtf && chrX.fa && check.tool.result && (length(raw.fastq) != 0) && check.progress.results.bool
  if (!isTRUE(validity)) {
    stop("RNAseqRawReadProcess() environment ERROR")
  }
  cat("(\u2714) : RNAseqRawReadProcess() pre-check is valid\n\n")
}

PostRNAseqRawReadProcess <- function(path.prefix, genome.name, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNAseqRawReadProcess()' environment post-check ...\n")
  # Still need to add condition
  gene_abundance <- dir.exists(paste0(path.prefix, "gene_data/gene_abundance/"))
  check.results <- ProgressGenesFiles(path.prefix = path.prefix, genome.name = genome.name, sample.pattern = sample.pattern, print=FALSE)
  ht2.bool <- (check.results$ht2.files.number.df) != 0
  sam.bool <- (check.results$sam.files.number.df) != 0
  bam.bool <- (check.results$bam.files.number.df) != 0
  gtf.bool <- (check.results$gtf.files.number.df) != 0
  merged.bool <- check.results$stringtie_merged.gtf.file.df
  gffcompare.bool <- (check.results$gffcompare.related.dirs.number.df) != 0
  ballgown.bool <- (check.results$ballgown.dirs.number.df) != 0
  validity <- gene_abundance && ht2.bool && sam.bool && bam.bool && gtf.bool && merged.bool && gffcompare.bool && ballgown.bool
  if (!isTRUE(validity)) {
    stop("RNAseqRawReadProcess() post-check ERROR")
  }
  cat("(\u2714) : RNAseqRawReadProcess() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}

