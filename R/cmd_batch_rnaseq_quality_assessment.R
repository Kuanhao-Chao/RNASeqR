#'
#' @export
RNAseqQualityAssessment_CMD <- function(RNASeqWorkFlowParam) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Assessment.R"))
  first <- "library(RNASeqWorkflow)"
  second <- "library(ggplot2)"
  third <- paste0("RNAseqQualityAssessment(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', sample.pattern = '", sample.pattern, "')")
  writeLines(c(first, second, third), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Quality_Assessment.R' has been created.\n"))
  system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Assessment.R ", path.prefix, "Rscript_out/Quality_Assessment.Rout"), stdout = "", wait = FALSE)
  cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Quality_Assessment.Rout'\n\n"))
}

#' Quality control
#' @export
RNAseqQualityAssessment <- function(path.prefix, input.path.prefix, sample.pattern) {
  trimmed.raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/trimmed_fastq.gz'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  if (dir.exists(paste0(path.prefix, "gene_data/raw_fastq.gz/trimmed_fastq.gz")) && (length(trimmed.raw.fastq) != 0)) {
    cat(paste0("\n************** Quality Assessment **************\n"))
    folder <- paste0(path.prefix, "gene_data/raw_fastq.gz/")
    files <- list.files(folder, sample.pattern, full.names=TRUE)
    cat(paste0("     \u25CF  R package \"Rqc\" quality assessment\n"))
    cat(paste0("          \u25CF  Running 'rqcQA()' ...  Please wait \u231B\u231B\u231B\n"))
    qa <- Rqc::rqcQA(files, workers=1)
    cat(paste0("          \u25CF  Creating 'rqc_report.html' ...  Please wait \u231B\u231B\u231B\n"))
    reportFile <- Rqc::rqcReport(qa)
    file.rename(from = reportFile, to = paste0(path.prefix, "RNAseq_results/QA_results/Rqc/Rqc_report.html"))
    cat(paste0("          (\u2714) : Rqc assessment success ~~\n\n"))

    cat(paste0("     \u25CF  R package \"systemPipeR\" quality assessment\n"))
    current.path <- getwd()
    setwd(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR"))
    # create 'rnaseq' directory
    systemPipeRdata::genWorkenvir(workflow="rnaseq")
    # create targets.txt
    cat(paste0("          \u25CF  Writing \"data.list.txt\""))
    raw.fastq.data.frame <- data.frame("FileName" = files, "SampleName" = 1:length(files), "SampleLong" = 1:length(files), "Experiment" = 1:length(files), "Date" = 1:length(files))
    write.table(raw.fastq.data.frame, "data.list.txt", sep="\t", row.names = FALSE, quote=FALSE)
    args <- systemPipeR::systemArgs(sysma="rnaseq/param/trim.param", mytargets="data.list.txt")
    cat(paste0("          \u25CF  Running 'seeFastq()' ...  Please wait \u231B\u231B\u231B\n"))
    fqlist <- systemPipeR::seeFastq(fastq=systemPipeR::infile1(args), batchsize=10000, klength=8)
    cat(paste0("          \u25CF  Creating 'fastqReport.pdf' ...  Please wait \u231B\u231B\u231B\n"))
    pdf(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/fastqReport.pdf"), height=18, width=4*length(fqlist))
    systemPipeR::seeFastqPlot(fqlist)
    dev.off()
    on.exit(setwd(current.path))
    cat(paste0("          \u25CF  Removing 'rnaseq' directory...  Please wait \u231B\u231B\u231B\n"))
    unlink(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/rnaseq"), recursive = TRUE)
    cat(paste0("          (\u2714) : systemPipeR assessment success ~~\n\n"))

    cat(paste0("     \u25CF  R package \"ShortRead\" quality assessment\n"))
    files <- list.files(folder, sample.pattern, full.names=TRUE)

    cat(paste0("          \u25CF  Running 'qa()' ...  Please wait \u231B\u231B\u231B\n"))
    qaSummary <- ShortRead::qa(files, type="fastq")
    cat(paste0("          \u25CF  Creating 'ShortRead_report.html' ...  Please wait \u231B\u231B\u231B\n"))
    resultFile <- ShortRead::report(qaSummary)
    file.rename(from = resultFile, to = paste0(path.prefix, "RNAseq_results/QA_results/ShortRead/ShortRead_report.html"))
    cat(paste0("          (\u2714) : ShortRead assessment success ~~\n\n"))
    CheckQAFiles(path.prefix)
  } else {
    stop("RNAseqQualityAssessment environment ERROR")
  }
}

CheckQAFiles <- function(path.prefix) {
  file.rqc.result <- file.exists(paste0(path.prefix, "RNAseq_results/QA_results/Rqc/Rqc_report.html"))
  file.systemPipeR.data <- file.exists(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/data.list.txt"))
  file.systemPipeR.result <- file.exists(paste0(path.prefix, "RNAseq_results/QA_results/systemPipeR/fastqReport.pdf"))
  file.ShortRead.result <- file.exists(paste0(path.prefix, "RNAseq_results/QA_results/ShortRead/ShortRead_report.html"))
  if (file.rqc.result && file.systemPipeR.data && file.systemPipeR.result && file.ShortRead.result) {
    cat("\n")
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  } else {
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
    stop("QA file missing ERROR")
  }
}
