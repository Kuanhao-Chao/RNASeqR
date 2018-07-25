#'
#' @export
RNAseqQualityTrimming_CMD <- function(RNASeqWorkFlowParam, truncateStartBases = 0, truncateEndBases = 0, complexity = NULL, minLength = 50, nBases = 2, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Trimming.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqQualityTrimming(path.prefix = '", path.prefix, "', sample.pattern = '", sample.pattern, "', truncateStartBases = ", truncateStartBases, ", truncateEndBases = ", truncateEndBases, ", complexity = ", complexity, ", minLength = ", minLength, ", nBases = ", nBases,  ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Quality_Trimming.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Trimming.R ", path.prefix, "Rscript_out/Quality_Trimming.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Quality_Trimming.Rout'\n\n"))
  }
}


#' @export
RNAseqQualityTrimming <- function(path.prefix, sample.pattern, truncateStartBases = 0, truncateEndBases = 0, complexity = NULL, minLength = 50, nBases = 2) {
  CheckOperatingSystem(FALSE)
  PreCheckRNAseqQualityTrimming(path.prefix = path.prefix, sample.pattern = sample.pattern)
  cat(paste0("************** Quality Trimming **************\n"))
  samples.fastq.untrim.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.fastq.untrim.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/'.\n     Please check whether the directory is already exit.\n"))
  }
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  raw.fastq.unique <- unique(gsub("[1-2]*.fastq.gz$", replace = "", raw.fastq))
  lapply(raw.fastq.unique, myFilterAndTrim, path.prefix = path.prefix, truncateStartBases = truncateStartBases,
         truncateEndBases = truncateEndBases, complexity = complexity, minLength = minLength, nBases = nBases)
  cat("\n")
  PostCheckRNAseqQualityTrimming(path.prefix = path.prefix, sample.pattern = sample.pattern)
}


myFilterAndTrim <- function(fl.name, path.prefix, truncateStartBases, truncateEndBases, complexity, minLength, nBases) {
  # adding print log
  # file1 and file2 is original fastq.gz without trimmed
  cat(paste0("\u25CF \"", gsub("_", "", fl.name), "\" quality trimming\n"))
  file1 <- paste0(path.prefix, "gene_data/raw_fastq.gz/", fl.name, "1.fastq.gz")
  file2 <- paste0(path.prefix, "gene_data/raw_fastq.gz/", fl.name, "2.fastq.gz")
  if (file.exists(file1) && file.exists(file2)) {
    # file1.output and file2.output are the new original fastq.gz file name
    file1.output <- paste0(path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/", fl.name, "1.fastq.gz")
    file2.output <- paste0(path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/", fl.name, "2.fastq.gz")
    cat(paste0("     \u25CF Moving \"", basename(file1), "\" to \"", path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/\"\n"))
    cat(paste0("     \u25CF Moving \"", basename(file2), "\" to \"", path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/\"\n"))
    file.rename(from = file1, to = file1.output)
    file.rename(from = file2, to = file2.output)
    #Sequence complexity (H) is calculated based on the dinucleotide composition using the formula (Shannon entropy):
    cat(paste0("     \u25CF Start trimming ...\n"))
    QuasR::preprocessReads(filename = file1.output, outputFilename = file1, filenameMate = file2.output, outputFilenameMate = file2,
                           truncateStartBases = 0, truncateEndBases = 0, complexity = NULL, minLength = 50, nBases = 2)
    cat(paste0("     \u25CF \"", file1, "\" has been created.\n"))
    cat(paste0("     \u25CF \"", file2, "\" has been created.\n"))
  } else {
    stop("paired-end file ERROR")
  }
}

PreCheckRNAseqQualityTrimming <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C 'RNAseqQualityTrimming()' environment pre-check ...\n")
  # have fastq.gz files
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  validity <- length(raw.fastq) != 0
  if (!isTRUE(validity)) {
    stop("CheckRNAseqQualityTrimming() pre-check ERROR")
  }
  cat("     (\u2714) : RNAseqQualityTrimming() pre-check is valid\n\n")
}

PostCheckRNAseqQualityTrimming <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C 'RNAseqQualityTrimming()' environment post-check ...\n")
  # have fastq.gz and trimmed fastq.gz files
  trimmed.raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/original_untrimmed_fastq.gz'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  validity <- (length(trimmed.raw.fastq) != 0) && (length(raw.fastq) != 0)
  if (!isTRUE(validity)) {
    stop("RNAseqQualityTrimming() post-check ERROR")
  }
  cat("     (\u2714) : RNAseqQualityTrimming() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}
