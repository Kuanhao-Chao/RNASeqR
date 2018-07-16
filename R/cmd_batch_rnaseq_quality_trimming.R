#'
#' @export
RNAseqQualityTrimming_CMD <- function(RNASeqWorkFlowParam) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_trimming.R"))
  first <- "library(RNASeqWorkflow)"
  second <- "library(ggplot2)"
  third <- paste0("QualityTrimming(path.prefix = '", path.prefix, "')")
  writeLines(c(first, second, third), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Quality_trimming.R' has been created.\n"))
  system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_trimming.R ", path.prefix, "Rscript_out/Quality_trimming.Rout"), stdout = "", wait = FALSE)
  cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Quality_trimming.Rout'\n\n"))
}


QualityTrimming <- function(path.prefix) {
  raw.fastq <- list.files(path = paste0(input.path.prefix, 'input_files/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE)
  lapply(raw.fastq, myFilterAndTrim)
}


myFilterAndTrim <- function(fl, path.prefix) {
  # adding print log
  ## open input stream
  rfq <- ShortRead::readFastq(fl)
  ShortRead::quality(rfq)
  ShortRead::alphabet(quality(rfq))
  ## trim and filter, e.g., reads cannot contain 'N'...
  rfq.Na <- rfq[ShortRead::nFilter()(rfq)]  # see ?srFilter for pre-defined filters
  ## trim as soon as 2 of 5 nucleotides has quality encoding less
  ## than "4" (phred score 20)
  rfq.trim <- ShortRead::trimTailw(rfq.Na, 2, "4", 2)
  ## drop reads that are less than 36nt
  rfq.less <- ShortRead::rfq.trim[width(rfq.trim) >= 36]
  ## append to destination
  trimmed.file.name <- basename(fl)
  destination_2 <- paste0(path.prefix,"/gene_data/raw_fastq.gz/trimmed_fastq.gz/", trimmed.file.name)
  ShortRead::writeFastq(rfq.less, destination, "a")
  # rfq1 = trimEnds(rfq, "4")
  # quality(rfq1)
  # repeat {
  #   ## input chunk
  #   fq <- yield(stream)
  #   if (length(fq) == 0)
  #     break
  #   ## trim and filter, e.g., reads cannot contain 'N'...
  #   fq <- fq[nFilter()(fq)]  # see ?srFilter for pre-defined filters
  #   ## trim as soon as 2 of 5 nucleotides has quality encoding less
  #   ## than "4" (phred score 20)
  #   fq <- trimTailw(fq, 2, "4", 2)
  #   ## drop reads that are less than 36nt
  #   fq <- fq[width(fq) >= 36]
  #   ## append to destination
  #   writeFastq(fq, destination, "a")
  # }
}
