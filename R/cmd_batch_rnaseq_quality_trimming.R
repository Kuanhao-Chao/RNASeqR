#' @title Quality trimming of '.fastq.gz' files for RNA-Seq workflow in background.
#'
#' @description Trim '.fastq.gz' files for RNA-Seq workflow in background. This step is optional in the whole RNA-Seq workflow.
#' The trimming method is implemented by R package \code{ShortRead}
#' If you want to trim '.fastq.gz' files for the RNA-Seq workflow in R shell, please see \code{RNASeqQualityTrimming()} function.
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related parameters
#' @param cum.error Default \code{1}. Cut of threshold of cumulative probability of error per base.
#' @param reads.length.limit Default \code{36}. The shortest base pair length of short reads
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
#' RNASeqQualityTrimming_CMD(RNASeqWorkFlowParam = exp)}
RNASeqQualityTrimming_CMD <- function(RNASeqWorkFlowParam, cum.error = 1, reads.length.limit = 36, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Trimming.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNASeqQualityTrimming(path.prefix = '", path.prefix, "', sample.pattern = '", sample.pattern, "', cum.error = ", cum.error, ", reads.length.limit = ", reads.length.limit, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Quality_Trimming.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Trimming.R ", path.prefix, "Rscript_out/Quality_Trimming.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Quality_Trimming.Rout'\n\n"))
  }
}


#' @title Quality trimming of '.fastq.gz' files for RNA-Seq workflow in R shell
#'
#' @description Trim '.fastq.gz' files for RNA-Seq workflow in R shell. This step is optional in the whole RNA-Seq workflow.
#' It is strongly advised to run \code{RNASeqQualityTrimming_CMD()} directly. Running \code{RNASeqQualityTrimming_CMD()} will create 'Quality_Trimming.Rout' file in 'Rscript_out/' directory.
#' The trimming method is implemented by R package \code{ShortReads}
#' If you want to trim '.fastq.gz' files for the RNA-Seq workflow in background, please see \code{RNASeqQualityTrimming_CMD()} function.
#'
#' @param path.prefix Path prefix of 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @param sample.pattern  sample.pattern  Regular expression of paired-end fastq.gz files under 'input_files/raw_fastq.gz'. Expression not includes \code{_[1,2].fastq.gz}.
#' @param cum.error Default \code{1}. Cut of threshold of cumulative probability of error per base.
#' @param reads.length.limit Default \code{36}. The shortest base pair length of short reads
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNASeqEnvironmentSet(path.prefix = exp@@path.prefix, sample.pattern = exp@@sample.pattern)}
RNASeqQualityTrimming <- function(path.prefix, sample.pattern, cum.error = 1, reads.length.limit = 36) {
  CheckOperatingSystem(FALSE)
  PreCheckRNASeqQualityTrimming(path.prefix = path.prefix, sample.pattern = sample.pattern)
  cat(paste0("************** Quality Trimming **************\n"))
  if(!dir.exists(paste0(path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/"))){
    dir.create(file.path(paste0(path.prefix, 'gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/')), showWarnings = FALSE)
  }
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  raw.fastq.unique <- unique(gsub("[1-2]*.fastq.gz$", replacement = "", raw.fastq))
  lapply(raw.fastq.unique, myFilterAndTrim, path.prefix = path.prefix, cum.error = cum.error, reads.length.limit = reads.length.limit)
  cat("\n")
  PostCheckRNASeqQualityTrimming(path.prefix, sample.pattern)
}


myFilterAndTrim <- function(fl.name, path.prefix, cum.error, reads.length.limit) {
  # adding print log
  # file1 and file2 is original fastq.gz without trimmed
  cat(paste0("\u25CF \"", gsub("_", "", fl.name), "\" quality trimming\n"))
  file1 <- paste0(path.prefix, "gene_data/raw_fastq.gz/", fl.name, "1.fastq.gz")
  file2 <- paste0(path.prefix, "gene_data/raw_fastq.gz/", fl.name, "2.fastq.gz")
  if (file.exists(file1) && file.exists(file2)) {
    # file1.output and file2.output are the new original fastq.gz file name
    file1.untrimmed <- paste0(path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/", fl.name, "1.fastq.gz")
    file2.untrimmed <- paste0(path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/", fl.name, "2.fastq.gz")
    cat(paste0("     \u25CF Moving \"", basename(file1), "\" to \"", path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/\"\n"))
    cat(paste0("     \u25CF Moving \"", basename(file2), "\" to \"", path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/\"\n"))
    file.rename(from = file1, to = file1.untrimmed)
    file.rename(from = file2, to = file2.untrimmed)
    #Sequence complexity (H) is calculated based on the dinucleotide composition using the formula (Shannon entropy):
    cat(paste0("     \u25CF Start trimming ...\n"))
    file1.read <- ShortRead::readFastq(file1.untrimmed)
    file2.read <- ShortRead::readFastq(file2.untrimmed)

    cat(paste0("          \u25CF Getting quality score list as PhredQuality ...\n"))
    # get quality score list as PhredQuality
    qual1 <- as(Biostrings::quality(file1.read), "matrix")
    qual2 <- as(Biostrings::quality(file2.read), "matrix")

    # Calculate probability error per base (through column) ==> Q = -10log10(P)   or  P = 10^(-Q/10)
    cat(paste0("          \u25CF Calculating probability error per base ...\n"))
    pe1 <- apply(qual1, MARGIN = 2, function(x){10^(-(x/10))})
    pe2 <- apply(qual2, MARGIN = 2, function(x){10^(-(x/10))})
    # Calculate cpm of error
    cat(paste0("          \u25CF Calculating cumulative distribution probability of error per base ...\n"))
    cum.pe1 <- apply(pe1, MARGIN = 1, cumsum)
    cum.pe2 <- apply(pe2, MARGIN = 1, cumsum)

    # Get the trimming position of each file
    cat(paste0("          \u25CF Filtering out cumulative distribution probability of error per base < 1 ...\n"))
    trimPos1 <- apply(cum.pe1, 2, function(x) { min(min(which(x > cum.error)), length(x)) } )
    trimPos2 <- apply(cum.pe2, 2, function(x) { min(min(which(x > cum.error)), length(x)) } )

    # Get the trimPos for pair-end files
    cat(paste0("          \u25CF Finding trimming position for paired-end ...\n"))
    trimPos.together <- mapply(function(list1, list2) {min(list1, list2)}, list1 = trimPos1, list2 = trimPos2)

    trimmed.file1 <- ShortRead::narrow(x = file1.read, start = 1, end = trimPos.together)
    trimmed.file2 <- ShortRead::narrow(x = file2.read, start = 1, end = trimPos.together)

    ## drop reads that are less than 36nt
    cat(paste0("     \u25CF Removing reads that are less than ", reads.length.limit, " base pairs ...\n"))
    trimmed.file1 <- trimmed.file1[ShortRead::width(trimmed.file1) >= reads.length.limit]
    trimmed.file2 <- trimmed.file2[ShortRead::width(trimmed.file2) >= reads.length.limit]

    # write new fastaq files
    cat(paste0("     \u25CF Creating trimmed pair-end files ...\n"))
    ShortRead::writeFastq(trimmed.file1, file1, "w")
    ShortRead::writeFastq(trimmed.file2, file2, "w")

    cat(paste0("     \u25CF \"", file1, "\" has been created.\n"))
    cat(paste0("     \u25CF \"", file2, "\" has been created.\n\n"))
  } else {
    stop("paired-end file ERROR")
  }
}

PreCheckRNASeqQualityTrimming <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C 'RNASeqQualityTrimming()' environment pre-check ...\n")
  # have fastq.gz files
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  validity <- length(raw.fastq) != 0
  if (!isTRUE(validity)) {
    stop("CheckRNASeqQualityTrimming() pre-check ERROR")
  }
  cat("(\u2714) : RNASeqQualityTrimming() pre-check is valid\n\n")
}

PostCheckRNASeqQualityTrimming <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C 'RNASeqQualityTrimming()' environment post-check ...\n")
  # have fastq.gz and trimmed fastq.gz files
  trimmed.raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/original_untrimmed_fastq.gz'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  validity <- (length(trimmed.raw.fastq) != 0) && (length(raw.fastq) != 0)
  if (!isTRUE(validity)) {
    stop("RNASeqQualityTrimming() post-check ERROR")
  }
  cat("(\u2714) : RNASeqQualityTrimming() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}
