#'
#' @export
RNAseqQualityAssessment_CMD <- function(RNASeqWorkFlowParam) {
  # check input param
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Control.R"))
  first <- "library(RNASeqWorkflow)"
  second <- "library(ggplot2)"
  third <- paste0("QualityControlRqc(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', sample.pattern = '", sample.pattern, "')")
  writeLines(c(first, second, third), fileConn)
  close(fileConn)
  system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Control.R ", path.prefix, "Rscript_out/Quality_Control.Rout"), stdout = "", wait = FALSE)
  cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Quality_Control.Rout'\n\n"))
}

#' Quality control
#' @import Rqc
#' @import ggplot2
#' @import systemPipeR
#' @import systemPipeRdata
#' @export
QualityControlRqc <- function(path.prefix, input.path.prefix, sample.pattern) {
  cat(paste0("\n************** Quality Assessment **************\n"))
  folder <- paste0(path.prefix, "gene_data/raw_fastq.gz/")

  files <- list.files(folder, sample.pattern, full.names=TRUE)
  cat(paste0("     \u25CF  R package \"Rqc\" quality assessment\n"))
  cat(paste0("          \u25CF  Running 'rqcQA()' ...  Please wait \u231B\u231B\u231B\n"))
  qa <- rqcQA(files, workers=1)
  cat(paste0("          \u25CF  Creating 'rqc_report.html' ...  Please wait \u231B\u231B\u231B\n"))
  reportFile <- rqcReport(qa)
  file.rename(from = reportFile, to = paste0(path.prefix, "RNAseq_results/QA_results/Rqc/rqc_report.html"))
  cat(paste0("          (\u2714) : Rqc assessment success ~~\n\n"))

  cat(paste0("     \u25CF  R package \"systemPipeR\" quality assessment\n"))
  current.path <- getwd()
  setwd(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR"))
  # create 'rnaseq' directory
  genWorkenvir(workflow="rnaseq")
  # create targets.txt
  cat(paste0("          \u25CF  Writing \"data.list.txt\"\n"))
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz'), pattern = sample.pattern, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE)
  raw.fastq.data.frame <- data.frame("FileName" = raw.fastq, "SampleName" = 1:length(raw.fastq), "SampleLong" = 1:length(raw.fastq), "Experiment" = 1:length(raw.fastq), "Date" = 1:length(raw.fastq))
  write.table(raw.fastq.data.frame, "data.list.txt", sep="\t", row.names = FALSE, quote=FALSE)
  args <- systemArgs(sysma="rnaseq/param/trim.param", mytargets="data.list.txt")
  cat(paste0("          \u25CF  Running 'seeFastq()' ...  Please wait \u231B\u231B\u231B\n"))
  fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
  cat(paste0("          \u25CF  Creating 'fastqReport.pdf' ...  Please wait \u231B\u231B\u231B\n"))
  pdf(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/fastqReport.pdf"), height=18, width=4*length(fqlist))
  seeFastqPlot(fqlist)
  dev.off()
  on.exit()
  cat(paste0("          \u25CF  Removing 'rnaseq' directory...  Please wait \u231B\u231B\u231B\n"))
  unlink(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/rnaseq"), recursive = TRUE)
  cat(paste0("          (\u2714) : systemPipeR assessment success ~~\n\n"))
}
