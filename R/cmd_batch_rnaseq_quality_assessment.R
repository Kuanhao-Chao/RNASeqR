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
#' @param RNASeqRParam S4 object instance of
#'   experiment-related parameters
#' @param run Default value is \code{TRUE}. If \code{TRUE},
#'   'Rscript/Environment_Set.R' will be created and executed.
#'   The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{False}, 'Rscript/Environment_Set.R' will be
#'   created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of
#'   checking \code{RNASeqRParam} will be reported in
#'   'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqQualityAssessment_CMD(RNASeqRParam = yeast)}
RNASeqQualityAssessment_CMD <- function(RNASeqRParam,
                                        run            = TRUE,
                                        check.s4.print = TRUE) {
  # check input param
  which.s4.object <- CheckS4Object_All(RNASeqRParam, check.s4.print)
  if (which.s4.object == "RNASeqRParam_Sam") {
    stop("'RNASeqQualityAssessment_CMD' must use RNASeqRParam S4 object ")
  } else if (which.s4.object == "RNASeqRParam") {
    CheckOperatingSystem(FALSE)
    path.prefix <- "@"(RNASeqRParam, path.prefix)
    INSIDE.path.prefix <- "@"(RNASeqRParam, path.prefix)
    saveRDS(RNASeqRParam,
            file = paste0(INSIDE.path.prefix,
                          "gene_data/RNASeqRParam.rds"))
    fileConn <- file(paste0(path.prefix, "Rscript/Quality_Assessment.R"))
    first <- "library(RNASeqR)"
    second <- paste0("RNASeqQualityAssessment(RNASeqRParam = 'INSIDE'",
                     ", which.trigger = 'INSIDE'",
                     ", INSIDE.path.prefix = '", INSIDE.path.prefix,"')")
    writeLines(c(first, second), fileConn)
    close(fileConn)
    message("\u2605 '", path.prefix,
            "Rscript/Quality_Assessment.R' has been created.\n")
    if (run) {
      R.home.lib <- R.home()
      R.home.bin <- gsub("/lib/R", "/bin/R", R.home.lib)
      system2(command = "nohup",
              args = paste0(R.home.bin, " CMD BATCH ",
                            path.prefix,
                            "Rscript/Quality_Assessment.R ",
                            path.prefix, "Rscript_out/Quality_Assessment.Rout"),
              stdout = "",
              wait = FALSE)
      message("\u2605 Tools are installing in the background. ",
              "Check current progress in '", path.prefix,
              "Rscript_out/Quality_Assessment.Rout'\n\n")
    }
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
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param which.trigger Default value is \code{OUTSIDE}. User should not change
#'   this value.
#' @param INSIDE.path.prefix Default value is \code{NA}. User should not change
#'   this value.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of
#'   checking \code{RNASeqRParam} will be reported in
#'   'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqQualityAssessment(RNASeqRParam = yeast)}
RNASeqQualityAssessment <- function(RNASeqRParam,
                                    which.trigger      = "OUTSIDE",
                                    INSIDE.path.prefix = NA,
                                    check.s4.print     = TRUE) {
  CheckOperatingSystem(FALSE)
  # If `which.trigger` is OUTSIDE, then directory must be built
  # If `which.trigger` is INSIDE, then directory must not be
  #  built here(will created in CMD)
  if (isS4(RNASeqRParam) &
      which.trigger == "OUTSIDE" &
      is.na(INSIDE.path.prefix)) {
    # This is an external call!!
    # Check the S4 object(user input)
  } else if (RNASeqRParam == "INSIDE" &
             which.trigger == "INSIDE" &
             !is.na(INSIDE.path.prefix)) {
    # This is an internal call!!
    # Load the S4 object that saved in CMD process
    RNASeqRParam <- readRDS(paste0(INSIDE.path.prefix,
                                   "gene_data/RNASeqRParam.rds"))
  }
  which.s4.object <- CheckS4Object_All(RNASeqRParam, check.s4.print)
  if (which.s4.object == "RNASeqRParam_Sam") {
    stop("'RNASeqQualityAssessment_CMD' must use RNASeqRParam S4 object ")
  } else if (which.s4.object == "RNASeqRParam") {
    path.prefix <- "@"(RNASeqRParam, path.prefix)
    input.path.prefix <- "@"(RNASeqRParam, input.path.prefix)
    sample.pattern <- "@"(RNASeqRParam, sample.pattern)
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
    message("************** Quality Assessment **************\n")
    folder <- paste0(path.prefix, "gene_data/raw_fastq.gz")
    files <- list.files(folder, sample.pattern, full.names = TRUE)
    # Other package reports
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
    message("\u25CF 2. R package \"systemPipeR\" quality assessment\n")
    current.path <- getwd()
    setwd(paste0(path.prefix,
                 "RNASeq_results/QA_results/QA_",
                 QA.count,
                 "/systemPipeR/"))
    # create 'RNASeq' directory
    systemPipeRdata::genWorkenvir(workflow="rnaseq")
    # create targets.txt
    message("     \u25CF  Writing \"data.list.txt\"\n")
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
    message("     \u25CF  Running 'seeFastq()' ...  ",
            "Please wait \u231B\u231B\u231B\n")
    fqlist <- systemPipeR::seeFastq(fastq=systemPipeR::infile1(args),
                                    batchsize=10000, klength=8)
    message("     \u25CF  Creating 'fastqReport.pdf' ...  ",
            "Please wait \u231B\u231B\u231B\n")
    pdf(paste0(path.prefix, "RNASeq_results/QA_results/QA_", QA.count,
               "/systemPipeR/fastqReport.pdf"),
        height=18, width=4*length(fqlist))
    systemPipeR::seeFastqPlot(fqlist)
    dev.off()
    on.exit(setwd(current.path))
    message("     \u25CF  Removing 'RNASeq' directory...  ",
            "Please wait \u231B\u231B\u231B\n")
    unlink(paste0(path.prefix, "RNASeq_results/QA_results/QA_",
                  QA.count, "/systemPipeR/rnaseq"),
           recursive = TRUE)
    message("     (\u2714) : systemPipeR assessment success ~~\n\n")
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
      message("'", path.prefix, "RNASeq_results/QA_results/QA_",
              QA.count,
              "/systemPipeR/data.list.txt' is missing!\n")
    }
    if (!file.systemPipeR.result){
      message("'", path.prefix,
              "RNASeq_results/QA_results/QA_",
              QA.count, "/systemPipeR/fastqReport.pdf' is missing!\n")
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
    message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
            "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
            "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
            "\u2605\u2605\u2605\u2605\n")
    message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
            "\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605",
            "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
    message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
            "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
            "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
            "\u2605\u2605\u2605\u2605\n")
  }
}

