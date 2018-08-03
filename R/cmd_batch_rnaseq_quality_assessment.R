#' @title Quality assessment of '.fastq.gz' files for RNA-Seq workflow in background.
#'
#' @description Assess the quality of '.fastq.gz' files for RNA-Seq workflow in background. This step is optional in the whole RNA-Seq workflow.
#' This function reports the quality assessment result in two packages : \code{Rqc} and \code{systemPipeR}
#' For \code{Rqc}, 'RNAseq_results/QA_results/Rqc/Rqc_report.html' will be created
#' For \code{systemPipeR}, 'RNAseq_results/QA_results/Rqc/systemPipeR/fastqReport.pdf' will be created
#' If you want to assess the quality of '.fastq.gz' files for the following RNA-Seq workflow in R shell, please see \code{RNAseqQualityAssessment()} function.
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related parameters
#' @param run Default value is \code{TRUE}. If \code{TRUE}, 'Rscript/Environment_Set.R' will be created and executed. The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#' If \code{False}, 'Rscript/Environment_Set.R' will be created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of checking \code{RNASeqWorkFlowParam} will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking \code{RNASeqWorkFlowParam} will not be in 'Rscript_out/Environment_Set.Rout'
#'
#' @export
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNAseqQualityAssessment_CMD(RNASeqWorkFlowParam = exp)}
RNAseqQualityAssessment_CMD <- function(RNASeqWorkFlowParam, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Assessment.R"))
  first <- "library(RNASeqWorkflow)"
  second <- "library(ggplot2)"
  third <- "library(ShortRead)"
  fourth <- paste0("RNAseqQualityAssessment(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', sample.pattern = '", sample.pattern, "')")
  writeLines(c(first, second, third, fourth), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Quality_Assessment.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Assessment.R ", path.prefix, "Rscript_out/Quality_Assessment.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Quality_Assessment.Rout'\n\n"))
  }
}

#' @title Quality assessment of '.fastq.gz' files for RNA-Seq workflow in R shell
#'
#' @description Assess the quality of '.fastq.gz' files for RNA-Seq workflow in R shell. This step is optional in the whole RNA-Seq workflow.
#' It is strongly advised to run \code{RNAseqQualityAssessment_CMD()} directly. Running this function directly is not recommended.
#' This function reports the quality assessment result in two packages : \code{Rqc} and \code{systemPipeR}
#' For \code{Rqc}, 'RNAseq_results/QA_results/Rqc/Rqc_report.html' will be created
#' For \code{systemPipeR}, 'RNAseq_results/QA_results/Rqc/systemPipeR/fastqReport.pdf' will be created
#' If you want to assess the quality of '.fastq.gz' files for the following RNA-Seq workflow in background, please see \code{RNAseqQualityAssessment_CMD()} function.
#'
#' @param path.prefix path prefix of 'gene_data/', 'RNAseq_bin/', 'RNAseq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @param input.path.prefix path prefix of 'input_files/' directory
#' @param sample.pattern  regular expression of raw fastq.gz files under 'input_files/raw_fastq.gz'
#'
#' @export
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNAseqEnvironmentSet(path.prefix = exp@@path.prefix, input.path.prefix = exp@@input.path.prefix,
#'                      sample.pattern = exp@@sample.pattern)}
RNAseqQualityAssessment <- function(path.prefix, input.path.prefix, sample.pattern) {
  CheckOperatingSystem(FALSE)
  PreCheckRNAseqQualityAssessment(path.prefix = path.prefix, sample.pattern = sample.pattern)

  if(!dir.exists(paste0(path.prefix, "RNAseq_results/QA_results/"))){
    dir.create(file.path(paste0(path.prefix, 'RNAseq_results/QA_results/')), showWarnings = FALSE)
  }
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/QA_results/Rqc/"))){
    dir.create(file.path(paste0(path.prefix, 'RNAseq_results/QA_results/Rqc/')), showWarnings = FALSE)
  }
  trimmed.raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  cat(paste0("************** Quality Assessment **************\n"))
  folder <- paste0(path.prefix, "gene_data/raw_fastq.gz/")
  files <- list.files(folder, sample.pattern, full.names=TRUE)
  cat(paste0("\u25CF 1. R package \"Rqc\" quality assessment\n"))
  cat(paste0("     \u25CF  Running 'rqcQA()' ...  Please wait \u231B\u231B\u231B\n"))
  qa <- Rqc::rqcQA(files)
  cat(paste0("     \u25CF  Creating 'rqc_report.html' ...  Please wait \u231B\u231B\u231B\n"))
  reportFile <- Rqc::rqcReport(qa)
  file.rename(from = reportFile, to = paste0(path.prefix, "RNAseq_results/QA_results/Rqc/Rqc_report.html"))
  cat(paste0("     (\u2714) : Rqc assessment success ~~\n\n"))

  if(!dir.exists(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/"))){
    dir.create(file.path(paste0(path.prefix, 'RNAseq_results/QA_results/systemPipeR/')), showWarnings = FALSE)
  }
  cat(paste0("\u25CF 2. R package \"systemPipeR\" quality assessment\n"))
  current.path <- getwd()
  setwd(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR"))
  # create 'rnaseq' directory
  systemPipeRdata::genWorkenvir(workflow="rnaseq")
  # create targets.txt
  cat(paste0("     \u25CF  Writing \"data.list.txt\""))
  raw.fastq.data.frame <- data.frame("FileName" = files, "SampleName" = 1:length(files), "SampleLong" = 1:length(files), "Experiment" = 1:length(files), "Date" = 1:length(files))
  write.table(raw.fastq.data.frame, "data.list.txt", sep="\t", row.names = FALSE, quote=FALSE)
  args <- systemPipeR::systemArgs(sysma="rnaseq/param/trim.param", mytargets="data.list.txt")
  cat(paste0("     \u25CF  Running 'seeFastq()' ...  Please wait \u231B\u231B\u231B\n"))
  fqlist <- systemPipeR::seeFastq(fastq=systemPipeR::infile1(args), batchsize=10000, klength=8)
  cat(paste0("     \u25CF  Creating 'fastqReport.pdf' ...  Please wait \u231B\u231B\u231B\n"))
  pdf(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/fastqReport.pdf"), height=18, width=4*length(fqlist))
  systemPipeR::seeFastqPlot(fqlist)
  dev.off()
  on.exit(setwd(current.path))
  cat(paste0("     \u25CF  Removing 'rnaseq' directory...  Please wait \u231B\u231B\u231B\n"))
  unlink(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/rnaseq"), recursive = TRUE)
  cat(paste0("     (\u2714) : systemPipeR assessment success ~~\n\n"))

  if(!dir.exists(paste0(path.prefix, "RNAseq_results/QA_results/ShortRead/"))){
    dir.create(file.path(paste0(path.prefix, 'RNAseq_results/QA_results/ShortRead/')), showWarnings = FALSE)
  }
  cat(paste0("\u25CF 3. R package \"ShortRead\" quality assessment\n"))
  files <- list.files(folder, sample.pattern, full.names=TRUE)

  cat(paste0("     \u25CF  Running 'qa()' ...  Please wait \u231B\u231B\u231B\n"))
  qaSummary <- ShortRead::qa(files, type="fastq")
  cat(paste0("     \u25CF  Creating 'ShortRead_report.html' ...  Please wait \u231B\u231B\u231B\n"))
  resultFile <- ShortRead::report(qaSummary)
  file.rename(from = resultFile, to = paste0(path.prefix, "RNAseq_results/QA_results/ShortRead/ShortRead_report.html"))
  cat(paste0("     (\u2714) : ShortRead assessment success ~~\n\n"))
  PostCheckRNAseqQualityAssessment(path.prefix = path.prefix)
}

PreCheckRNAseqQualityAssessment <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C 'RNAseqQualityAssessment()' environment pre-check ...\n")
  # 1. must have .fastq.gz files
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  validity <- length(raw.fastq) != 0
  if (!isTRUE(validity)) {
    stop("RNAseqQualityAssessment environment() ERROR")
  }
  cat("     (\u2714) : RNAseqQualityAssessment() pre-check is valid\n\n")
}

PostCheckRNAseqQualityAssessment <- function(path.prefix) {
  cat("\u269C\u265C\u265C\u265C 'RNAseqQualityAssessment()' environment post-check ...\n")
  # Assessment results exist
  file.rqc.result <- file.exists(paste0(path.prefix, "RNAseq_results/QA_results/Rqc/Rqc_report.html"))
  file.systemPipeR.data <- file.exists(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/data.list.txt"))
  file.systemPipeR.result <- file.exists(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/fastqReport.pdf"))
  file.ShortRead.result <- file.exists(paste0(path.prefix, "RNAseq_results/QA_results/ShortRead/ShortRead_report.html"))
  validity <- file.rqc.result && file.systemPipeR.data && file.systemPipeR.result && file.ShortRead.result
  if (!isTRUE(validity)) {
    if (!file.rqc.result) {
      cat(paste0("'", path.prefix, "RNAseq_results/QA_results/Rqc/Rqc_report.html' is missing!\n"))
    }
    if (!file.systemPipeR.data) {
      cat(paste0("'", path.prefix, "RNAseq_results/QA_results/systemPipeR/data.list.txt' is missing!\n"))
    }
    if (!file.systemPipeR.result){
      cat(paste0("'", path.prefix, "RNAseq_results/QA_results/systemPipeR/fastqReport.pdf' is missing!\n"))
    }
    if (!file.ShortRead.result) {
      cat(paste0("'", path.prefix, "RNAseq_results/QA_results/ShortRead/ShortRead_report.html' is missing!\n"))
    }
    stop("RNAseqQualityAssessment() post-check ERROR")
  } else {
    cat("     (\u2714) : RNAseqQualityAssessment() post-check is valid\n\n")
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  }
}

