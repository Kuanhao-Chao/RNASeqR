#'
#' @export
RNAseqQualityTrimming_CMD <- function(RNASeqWorkFlowParam, trimming.score = 30) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Trimming.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqQualityTrimming(path.prefix = '", path.prefix, "', sample.pattern = '", sample.pattern, "', trimming.score = ", trimming.score, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Quality_Trimming.R' has been created.\n"))
  system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Trimming.R ", path.prefix, "Rscript_out/Quality_Trimming.Rout"), stdout = "", wait = FALSE)
  cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Quality_Trimming.Rout'\n\n"))
}


RNAseqQualityTrimming <- function(path.prefix, sample.pattern, trimming.score = 30) {
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  if (dir.exists(paste0(path.prefix, "gene_data/raw_fastq.gz/trimmed_fastq.gz/")) && (length(raw.fastq) != 0)) {
    raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE)
    lapply(raw.fastq, myFilterAndTrim, path.prefix = path.prefix, trimming.score = trimming.score)
  } else {
    stop("RNAseqQualityTrimming environment ERROR")
  }
}


myFilterAndTrim <- function(fl, path.prefix, trimming.score) {
  # adding print log
  ## open input stream
  rfq <- ShortRead::readFastq(fl)
  Biostrings::quality(rfq)
  ShortRead::alphabet(Biostrings::quality(rfq))[1]
  trimming.alph <- names(encoding(quality(rfq))[trimming.score+1])
  ## trim and filter, e.g., reads cannot contain 'N'...
  rfq.Na <- rfq[ShortRead::nFilter()(rfq)]  # see ?srFilter for pre-defined filters
  ## trim as soon as 2 of 5 nucleotides has quality encoding less
  ## than "4" (phred score 20)
  rfq.trim <- ShortRead::trimTailw(rfq.Na, 2, trimming.alph, 2)
  ## drop reads that are less than 50nt
  rfq.less <- rfq.trim[BiocGenerics::width(rfq.trim) >= 50]
  ## append to destination
  trimmed.file.name <- basename(fl)
  destination <- paste0(path.prefix,"/gene_data/raw_fastq.gz/trimmed_fastq.gz/", trimmed.file.name)
  ShortRead::writeFastq(rfq.less, destination, "a")
}
