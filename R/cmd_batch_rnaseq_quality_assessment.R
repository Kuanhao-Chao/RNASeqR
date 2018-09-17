#' @title RNASeqQualityAssessment_CMD
#'
#' @description
#'   Assess the quality of '.fastq.gz' files for RNA-Seq workflow in background.
#'   This step is optional in the whole RNA-Seq workflow. \cr
#'   This function reports the quality assessment result in packages
#'   \code{systemPipeR}
#'   For \code{systemPipeR},
#'   'RNASeq_results/QA_results/Rqc/systemPipeR/fastqReport.pdf'
#'   will be created. \cr If you want to assess the quality of '.fastq.gz'
#'   files for the following RNA-Seq workflow in R shell,
#'   please see \code{RNASeqQualityAssessment()} function.
#'
#' @param RNASeqWorkFlowParam S4 object instance of
#'   experiment-related parameters
#' @param run Default value is \code{TRUE}. If \code{TRUE},
#'   'Rscript/Environment_Set.R' will be created and executed.
#'   The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{False}, 'Rscript/Environment_Set.R' will be
#'   created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of
#'   checking \code{RNASeqWorkFlowParam} will be reported in
#'   'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqWorkFlowParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqQualityAssessment_CMD(RNASeqWorkFlowParam = yeast)}
RNASeqQualityAssessment_CMD <- function(RNASeqWorkFlowParam,
                                        run                 = TRUE,
                                        check.s4.print      = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- "@"(RNASeqWorkFlowParam, path.prefix)
  input.path.prefix <- "@"(RNASeqWorkFlowParam, input.path.prefix)
  sample.pattern <- "@"(RNASeqWorkFlowParam, sample.pattern)
  fileConn <- file(paste0(path.prefix, "Rscript/Quality_Assessment.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNASeqQualityAssessment(path.prefix = '", path.prefix,
                   "', input.path.prefix = '", input.path.prefix,
                   "', sample.pattern = '", sample.pattern, "')")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message(paste0("\u2605 '", path.prefix,
                 "Rscript/Quality_Assessment.R' has been created.\n"))
  if (run) {
    system2(command = "nohup",
            args = paste0("R CMD BATCH ", path.prefix,
                                             "Rscript/Quality_Assessment.R ",
                          path.prefix, "Rscript_out/Quality_Assessment.Rout"),
            stdout = "",
            wait = FALSE)
    message(paste0("\u2605 Tools are installing in the background. ",
                   "Check current progress in '", path.prefix,
                   "Rscript_out/Quality_Assessment.Rout'\n\n"))
  }
}

#' @title RNASeqQualityAssessment
#'
#' @description
#'   Assess the quality of '.fastq.gz' files for RNA-Seq workflow in R shell.
#'   This step is optional in the whole RNA-Seq workflow. \cr
#'   This function reports the quality assessment result in packages
#'   \code{systemPipeR}
#'   For \code{systemPipeR},
#'   'RNASeq_results/QA_results/Rqc/systemPipeR/fastqReport.pdf'
#'   will be created. \cr If you want to assess the quality of '.fastq.gz'
#'   files for the following RNA-Seq workflow in background,
#'   please see \code{RNASeqQualityAssessment_CMD()} function.
#'
#' @param path.prefix path prefix of 'gene_data/', 'RNASeq_bin/',
#'   'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @param input.path.prefix path prefix of 'input_files/' directory
#' @param sample.pattern  sample.pattern  Regular expression of paired-end
#'   fastq.gz files under 'input_files/raw_fastq.gz'.
#'   Expression not includes \code{_[1,2].fastq.gz}.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqQualityAssessment(path.prefix       = yeast@@path.prefix,
#'                         input.path.prefix = yeast@@input.path.prefix,
#'                         sample.pattern    = yeast@@sample.pattern)}
RNASeqQualityAssessment <- function(path.prefix,
                                    input.path.prefix,
                                    sample.pattern) {
  CheckOperatingSystem(FALSE)
  PreCheckRNASeqQualityAssessment(path.prefix, sample.pattern)
  QA_results_subfiles <- list.files(paste0(path.prefix,
                                           "RNASeq_results/QA_results"),
                                    pattern = "QA_[0-9]*")
  QA.count <- length(QA_results_subfiles) + 1
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/QA_results/"))){
    dir.create(file.path(paste0(path.prefix, "RNASeq_results/QA_results/")),
               showWarnings = FALSE)
  }
  if(!dir.exists(paste0(path.prefix,
                        "RNASeq_results/QA_results/QA_",
                        QA.count))){
    dir.create(file.path(paste0(path.prefix,
                                'RNASeq_results/QA_results/QA_', QA.count)),
               showWarnings = FALSE)
  }
  # if(!dir.exists(paste0(path.prefix, "RNASeq_results/QA_results/QA_",
  #                       QA.count, "/Rqc/"))){
  #   dir.create(file.path(paste0(path.prefix,
  #                               "RNASeq_results/QA_results/QA_",
  #                               QA.count,
  #                               "/Rqc/")),
  #              showWarnings = FALSE)
  # }
  trimmed.raw.fastq <- list.files(path = paste0(path.prefix,
                                                "gene_data/raw_fastq.gz/"),
                                  pattern = sample.pattern,
                                  all.files = FALSE,
                                  full.names = FALSE,
                                  recursive = FALSE,
                                  ignore.case = FALSE)
  message(paste0("************** Quality Assessment **************\n"))
  folder <- paste0(path.prefix, "gene_data/raw_fastq.gz")
  files <- list.files(folder, sample.pattern, full.names = TRUE)
  # message(paste0("\u25CF 1. R package \"Rqc\" quality assessment\n"))
  # message(paste0("     \u25CF  Running 'rqcQA()' ...  ",
  #                "Please wait \u231B\u231B\u231B\n"))
  # qa <- Rqc::rqcQA(files)
  # message(paste0("     \u25CF  Creating 'rqc_report.html' ...  ",
  #                "Please wait \u231B\u231B\u231B\n"))
  # reportFile <- Rqc::rqcReport(qa)
  # file.rename(from = reportFile,
  #             to = paste0(path.prefix,
  #                         "RNASeq_results/QA_results/QA_",
  #                         QA.count,
  #                         "/Rqc/Rqc_report.html"))
  # message(paste0("     (\u2714) : Rqc assessment success ~~\n\n"))
  if(!dir.exists(paste0(path.prefix,
                        "RNASeq_results/QA_results/QA_",
                        QA.count,
                        "/systemPipeR/"))){
    dir.create(file.path(paste0(path.prefix,
                                "RNASeq_results/QA_results/QA_",
                                QA.count, "/systemPipeR/")),
               showWarnings = FALSE)
  }
  message(paste0("\u25CF 2. R package \"systemPipeR\" quality assessment\n"))
  current.path <- getwd()
  setwd(paste0(path.prefix,
               "RNASeq_results/QA_results/QA_",
               QA.count,
               "/systemPipeR/"))
  # create 'RNASeq' directory
  systemPipeRdata::genWorkenvir(workflow="rnaseq")
  # create targets.txt
  message(paste0("     \u25CF  Writing \"data.list.txt\"\n"))
  raw.fastq.data.frame <- data.frame("FileName" = files,
                                     "SampleName" = gsub(".fastq.gz",
                                                         "",
                                                         basename(files)),
                                     "SampleLong" = seq_len(length(files)))
  write.table(raw.fastq.data.frame,
              "data.list.txt",
              sep="\t",
              row.names = FALSE,
              quote = FALSE)
  args <- systemPipeR::systemArgs(sysma="rnaseq/param/trim.param",
                                  mytargets="data.list.txt")
  message(paste0("     \u25CF  Running 'seeFastq()' ...  ",
                 "Please wait \u231B\u231B\u231B\n"))
  fqlist <- systemPipeR::seeFastq(fastq=systemPipeR::infile1(args),
                                  batchsize=10000, klength=8)
  message(paste0("     \u25CF  Creating 'fastqReport.pdf' ...  ",
                 "Please wait \u231B\u231B\u231B\n"))
  pdf(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count,
             "/systemPipeR/fastqReport.pdf"),
      height=18, width=4*length(fqlist))
  systemPipeR::seeFastqPlot(fqlist)
  dev.off()
  on.exit(setwd(current.path))
  message(paste0("     \u25CF  Removing 'RNASeq' directory...  ",
                 "Please wait \u231B\u231B\u231B\n"))
  unlink(paste0(path.prefix, "RNASeq_results/QA_results/QA_",
                QA.count, "/systemPipeR/rnaseq"),
         recursive = TRUE)
  message(paste0("     (\u2714) : systemPipeR assessment success ~~\n\n"))
  # if(!dir.exists(paste0(path.prefix,
  #                       "RNASeq_results/QA_results/QA_",
  #                       QA.count,
  #                       "/ShortRead/"))){
  #   dir.create(file.path(paste0(path.prefix,
  #                               "RNASeq_results/QA_results/QA_",
  #                               QA.count,
  #                               "/ShortRead/")),
  #              showWarnings = FALSE)
  # }
  # message(paste0("\u25CF 3. R package \"ShortRead\" quality assessment\n"))
  # files <- list.files(folder, sample.pattern, full.names=TRUE)
  # message(paste0("     \u25CF  Running 'qa()' ...  ",
  #                "Please wait \u231B\u231B\u231B\n"))
  # qaSummary <- ShortRead::qa(files, type="fastq")
  # message(paste0("     \u25CF  Creating 'ShortRead_report.html' ...  ",
  #                "Please wait \u231B\u231B\u231B\n"))
  # resultFile <- ShortRead::report(qaSummary)
  # file.rename(from = resultFile,
  #             to = paste0(path.prefix,
  #                         "RNASeq_results/QA_results/QA_",
  #                         QA.count,
  #                         "/ShortRead/ShortRead_report.html"))
  # message(paste0("     (\u2714) : ShortRead assessment success ~~\n\n"))
  PostCheckRNASeqQualityAssessment(path.prefix = path.prefix)
}

PreCheckRNASeqQualityAssessment <- function(path.prefix, sample.pattern) {
  message("\u269C\u265C\u265C\u265C ",
          "'RNASeqQualityAssessment()' ",
          "environment pre-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  ref.gtf <- file.exists(paste0(path.prefix, "gene_data/ref_genes/ref.gtf"))
  ref.fa <- file.exists(paste0(path.prefix, "gene_data/ref_genome/ref.fa"))
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'),
                          pattern = sample.pattern,
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
  check.tool.result <- CheckToolAll(path.prefix)
  validity <- phenodata.csv && ref.gtf && ref.fa &&
    check.tool.result && (length(raw.fastq) != 0)
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'),
                          pattern = sample.pattern,
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
  validity <- length(raw.fastq) != 0
  if (!isTRUE(validity)) {
    stop("RNASeqQualityAssessment environment() ERROR")
  }
  message("(\u2714) : RNASeqQualityAssessment() pre-check is valid\n\n")
}

PostCheckRNASeqQualityAssessment <- function(path.prefix) {
  message("\u269C\u265C\u265C\u265C ",
          "'RNASeqQualityAssessment()' ",
          "environment post-check ...\n")
  # Assessment results exist
  QA_results_subfiles <- list.files(paste0(path.prefix,
                                           "RNASeq_results/QA_results"),
                                    pattern = "QA_[0-9]*")
  # Don't need to plus one
  QA.count <- length(QA_results_subfiles)
  # file.rqc.result <- file.exists(paste0(path.prefix,
  #                                       "RNASeq_results/QA_results/QA_",
  #                                       QA.count, "/Rqc/Rqc_report.html"))
  file.systemPipeR.data <- file.exists(paste0(path.prefix,
                                              "RNASeq_results/QA_results/QA_",
                                              QA.count,
                                              "/systemPipeR/data.list.txt"))
  file.systemPipeR.result <- file.exists(paste0(path.prefix,
                                                "RNASeq_results/QA_results/QA_",
                                                QA.count,
                                                "/systemPipeR/fastqReport.pdf"))
  validity <- file.systemPipeR.data && file.systemPipeR.result
  if (!isTRUE(validity)) {
    # if (!file.rqc.result) {
      # message(paste0("'", path.prefix,
      #                "RNASeq_results/QA_results/QA_",
      #                QA.count, "/Rqc/Rqc_report.html' is missing!\n"))
    # }
    if (!file.systemPipeR.data) {
      message(paste0("'", path.prefix, "RNASeq_results/QA_results/QA_",
                     QA.count,
                     "/systemPipeR/data.list.txt' is missing!\n"))
    }
    if (!file.systemPipeR.result){
      message(paste0("'", path.prefix,
                     "RNASeq_results/QA_results/QA_",
                     QA.count,
                     "/systemPipeR/fastqReport.pdf' is missing!\n"))
    }
    # if (!file.ShortRead.result) {
      # message(paste0("'", path.prefix,
      #                "RNASeq_results/QA_results/QA_",
      #                QA.count,
      #                "/ShortRead/ShortRead_report.html' is missing!\n"))
    # }
    stop("RNASeqQualityAssessment() post-check ERROR")
  } else {
    message("(\u2714) : RNASeqQualityAssessment() post-check is valid\n\n")
    message(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
                   "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
                   "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
                   "\u2605\u2605\u2605\u2605\n"))
    message(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
                   "\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605",
                   "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
    message(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
                   "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
                   "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
                   "\u2605\u2605\u2605\u2605\n"))
  }
}

