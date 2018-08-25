#' @title Quality assessment of '.fastq.gz' files for RNA-Seq workflow in background.
#'
#' @description Assess the quality of '.fastq.gz' files for RNA-Seq workflow in background. This step is optional in the whole RNA-Seq workflow.
#' This function reports the quality assessment result in packages \code{systemPipeR}
#' For \code{systemPipeR}, 'RNASeq_results/QA_results/Rqc/systemPipeR/fastqReport.pdf' will be created
#' If you want to assess the quality of '.fastq.gz' files for the following RNA-Seq workflow in R shell, please see \code{RNASeqQualityAssessment()} function.
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related parameters
#' @param run Default value is \code{TRUE}. If \code{TRUE}, 'Rscript/Environment_Set.R' will be created and executed. The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#' If \code{False}, 'Rscript/Environment_Set.R' will be created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of checking \code{RNASeqWorkFlowParam} will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking \code{RNASeqWorkFlowParam} will not be in 'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNASeqQualityAssessment_CMD(RNASeqWorkFlowParam = exp)}
RNASeqQualityAssessment_CMD <- function(RNASeqWorkFlowParam, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Assessment.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNASeqQualityAssessment(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', sample.pattern = '", sample.pattern, "')")
  writeLines(c(first, second), fileConn)
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
#' It is strongly advised to run \code{RNASeqQualityAssessment_CMD()} directly. Running \code{RNASeqQualityAssessment_CMD()} will create 'Quality_Assessment.Rout' file in 'Rscript_out/' directory.
#' This function reports the quality assessment result in one packages \code{systemPipeR}
#' For \code{systemPipeR}, 'RNASeq_results/QA_results/Rqc/systemPipeR/fastqReport.pdf' will be created
#' If you want to assess the quality of '.fastq.gz' files for the following RNA-Seq workflow in background, please see \code{RNASeqQualityAssessment_CMD()} function.
#'
#' @param path.prefix path prefix of 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @param input.path.prefix path prefix of 'input_files/' directory
#' @param sample.pattern  sample.pattern  Regular expression of paired-end fastq.gz files under 'input_files/raw_fastq.gz'. Expression not includes \code{_[1,2].fastq.gz}.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNASeqEnvironmentSet(path.prefix = exp@@path.prefix, input.path.prefix = exp@@input.path.prefix,
#'                      sample.pattern = exp@@sample.pattern)}
RNASeqQualityAssessment <- function(path.prefix, input.path.prefix, sample.pattern) {
  CheckOperatingSystem(FALSE)
  PreCheckRNASeqQualityAssessment(path.prefix = path.prefix, sample.pattern = sample.pattern)
  QA_results_subfiles <- list.files(paste0(path.prefix, "RNASeq_results/QA_results"), pattern = "QA_[0-9]*")
  QA.count <- length(QA_results_subfiles) + 1
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/QA_results/"))){
    dir.create(file.path(paste0(path.prefix, 'RNASeq_results/QA_results/')), showWarnings = FALSE)
  }
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count))){
    dir.create(file.path(paste0(path.prefix, 'RNASeq_results/QA_results/QA_', QA.count)), showWarnings = FALSE)
  }
  # if(!dir.exists(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/Rqc/"))){
  #   dir.create(file.path(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/Rqc/")), showWarnings = FALSE)
  # }
  trimmed.raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  cat(paste0("************** Quality Assessment **************\n"))
  folder <- paste0(path.prefix, "gene_data/raw_fastq.gz")
  files <- list.files(folder, sample.pattern, full.names=TRUE)
  # cat(paste0("\u25CF 1. R package \"Rqc\" quality assessment\n"))
  # cat(paste0("     \u25CF  Running 'rqcQA()' ...  Please wait \u231B\u231B\u231B\n"))
  # qa <- Rqc::rqcQA(files)
  # cat(paste0("     \u25CF  Creating 'rqc_report.html' ...  Please wait \u231B\u231B\u231B\n"))
  # reportFile <- Rqc::rqcReport(qa)
  # file.rename(from = reportFile, to = paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/Rqc/Rqc_report.html"))
  # cat(paste0("     (\u2714) : Rqc assessment success ~~\n\n"))

  if(!dir.exists(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/systemPipeR/"))){
    dir.create(file.path(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/systemPipeR/")), showWarnings = FALSE)
  }
  cat(paste0("\u25CF 2. R package \"systemPipeR\" quality assessment\n"))
  current.path <- getwd()
  setwd(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/systemPipeR/"))
  # create 'RNASeq' directory
  systemPipeRdata::genWorkenvir(workflow="rnaseq")
  # create targets.txt
  cat(paste0("     \u25CF  Writing \"data.list.txt\"\n"))
  raw.fastq.data.frame <- data.frame("FileName" = files, "SampleName" = gsub(".fastq.gz", "", basename(files)), "SampleLong" = seq_len(length(files)))
  write.table(raw.fastq.data.frame, "data.list.txt", sep="\t", row.names = FALSE, quote=FALSE)
  args <- systemPipeR::systemArgs(sysma="rnaseq/param/trim.param", mytargets="data.list.txt")
  cat(paste0("     \u25CF  Running 'seeFastq()' ...  Please wait \u231B\u231B\u231B\n"))
  fqlist <- systemPipeR::seeFastq(fastq=systemPipeR::infile1(args), batchsize=10000, klength=8)
  cat(paste0("     \u25CF  Creating 'fastqReport.pdf' ...  Please wait \u231B\u231B\u231B\n"))
  pdf(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/systemPipeR/fastqReport.pdf"), height=18, width=4*length(fqlist))
  systemPipeR::seeFastqPlot(fqlist)
  dev.off()
  on.exit(setwd(current.path))
  cat(paste0("     \u25CF  Removing 'RNASeq' directory...  Please wait \u231B\u231B\u231B\n"))
  unlink(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/systemPipeR/rnaseq"), recursive = TRUE)
  cat(paste0("     (\u2714) : systemPipeR assessment success ~~\n\n"))
  # if(!dir.exists(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/ShortRead/"))){
  #   dir.create(file.path(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/ShortRead/")), showWarnings = FALSE)
  # }
  # cat(paste0("\u25CF 3. R package \"ShortRead\" quality assessment\n"))
  # files <- list.files(folder, sample.pattern, full.names=TRUE)
  # cat(paste0("     \u25CF  Running 'qa()' ...  Please wait \u231B\u231B\u231B\n"))
  # qaSummary <- ShortRead::qa(files, type="fastq")
  # cat(paste0("     \u25CF  Creating 'ShortRead_report.html' ...  Please wait \u231B\u231B\u231B\n"))
  # resultFile <- ShortRead::report(qaSummary)
  # file.rename(from = resultFile, to = paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/ShortRead/ShortRead_report.html"))
  # cat(paste0("     (\u2714) : ShortRead assessment success ~~\n\n"))
  PostCheckRNASeqQualityAssessment(path.prefix = path.prefix)
}

PreCheckRNASeqQualityAssessment <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C 'RNASeqQualityAssessment()' environment pre-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  chrX.gtf <- file.exists(paste0(path.prefix, "gene_data/ref_genes/chrX.gtf"))
  chrX.fa <- file.exists(paste0(path.prefix, "gene_data/ref_genome/chrX.fa"))
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  ExportPath(path.prefix)
  check.tool.result <- CheckToolAll()
  validity <- phenodata.csv && chrX.gtf && chrX.fa && check.tool.result && (length(raw.fastq) != 0)
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  validity <- length(raw.fastq) != 0
  if (!isTRUE(validity)) {
    stop("RNASeqQualityAssessment environment() ERROR")
  }
  cat("(\u2714) : RNASeqQualityAssessment() pre-check is valid\n\n")
}

PostCheckRNASeqQualityAssessment <- function(path.prefix) {
  cat("\u269C\u265C\u265C\u265C 'RNASeqQualityAssessment()' environment post-check ...\n")
  # Assessment results exist
  QA_results_subfiles <- list.files(paste0(path.prefix, "RNASeq_results/QA_results"), pattern = "QA_[0-9]*")
  # Don't need to plus one
  QA.count <- length(QA_results_subfiles)
  # file.rqc.result <- file.exists(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/Rqc/Rqc_report.html"))
  file.systemPipeR.data <- file.exists(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/systemPipeR/data.list.txt"))
  file.systemPipeR.result <- file.exists(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/systemPipeR/fastqReport.pdf"))
  validity <- file.systemPipeR.data && file.systemPipeR.result
  if (!isTRUE(validity)) {
    # if (!file.rqc.result) {
    #   cat(paste0("'", path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/Rqc/Rqc_report.html' is missing!\n"))
    # }
    if (!file.systemPipeR.data) {
      cat(paste0("'", path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/systemPipeR/data.list.txt' is missing!\n"))
    }
    if (!file.systemPipeR.result){
      cat(paste0("'", path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/systemPipeR/fastqReport.pdf' is missing!\n"))
    }
    # if (!file.ShortRead.result) {
    #   cat(paste0("'", path.prefix, "RNASeq_results/QA_results/QA_", QA.count, "/ShortRead/ShortRead_report.html' is missing!\n"))
    # }
    stop("RNASeqQualityAssessment() post-check ERROR")
  } else {
    cat("(\u2714) : RNASeqQualityAssessment() post-check is valid\n\n")
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  }
}

