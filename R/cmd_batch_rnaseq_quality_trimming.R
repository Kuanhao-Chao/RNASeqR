#'
#' @export
RNAseqQualityTrimming_CMD <- function(RNASeqWorkFlowParam, trimming.score = 30, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Trimming.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqQualityTrimming(path.prefix = '", path.prefix, "', sample.pattern = '", sample.pattern, "', trimming.score = ", trimming.score, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Quality_Trimming.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Trimming.R ", path.prefix, "Rscript_out/Quality_Trimming.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Quality_Trimming.Rout'\n\n"))
  }
}


#' @export
RNAseqQualityTrimming <- function(path.prefix, sample.pattern, trimming.score = 30) {
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  if (dir.exists(paste0(path.prefix, "gene_data/raw_fastq.gz/trimmed_fastq.gz/")) && (length(raw.fastq) != 0)) {
    # raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz'), pattern = sample.pattern, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE)
    raw.fastq.unique <- unique(gsub("[1-2]*.fastq.gz$", replace = "", raw.fastq))
    lapply(raw.fastq, myFilterAndTrim, path.prefix = path.prefix, trimming.score = trimming.score)
  } else {
    stop("RNAseqQualityTrimming environment ERROR")
  }
}


myFilterAndTrim <- function(fl.name, path.prefix, minLength, nBases) {
  # adding print log
  file1 <- paste0(path.prefix, "gene_data/raw_fastq.gz/", fl.name, "1.fastq.gz")
  file2 <- paste0(path.prefix, "gene_data/raw_fastq.gz/", fl.name, "2.fastq.gz")
  if (file.exists(file1) && file.exists(file2)) {
    file1.output <- paste0(path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/", fl.name, "1.fastq.gz")
    file2.output <- paste0(path.prefix, "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/", fl.name, "2.fastq.gz")
    file.rename(from = file1, to = file1.output)
    file.rename(from = file2, to = file2.output)
    #Sequence complexity (H) is calculated based on the dinucleotide composition using the formula (Shannon entropy):
    QuasR::preprocessReads(filename = file1.output, outputFilename = file1, filenameMate = file2.output, outputFilenameMate = file2,
                           complexity = NULL, minLength = 50, nBases = 2)
  } else {
    stop("paired-end file ERROR")
  }
}
