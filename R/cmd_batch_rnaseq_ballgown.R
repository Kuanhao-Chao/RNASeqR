#' @export
RNAseqBallgownProcess_CMD <- function(RNASeqWorkFlowParam, ballgown.log2FC = 1, ballgown.qval = 0.05, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  gene.name <- RNASeqWorkFlowParam@gene.name
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  independent.variable <- RNASeqWorkFlowParam@independent.variable
  fileConn<-file(paste0(path.prefix, "Rscript/Ballgown_Process.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqBallgownProcess(path.prefix = '", path.prefix, "', gene.name = '", gene.name, "', sample.pattern = '", sample.pattern, "', independent.variable = '",independent.variable, "', ballgown.log2FC = ", ballgown.log2FC, ", ballgown.qval = ", ballgown.qval, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Ballgown_Process.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Ballgown_Process.R ", path.prefix, "Rscript_out/Ballgown_Process.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Ballgown_Process.Rout'\n\n"))
  }
}

#' @export
RNAseqBallgownProcess <- function(path.prefix, gene.name, sample.pattern, independent.variable, ballgown.log2FC = 1, ballgown.qval = 0.05) {
  CheckOperatingSystem(FALSE)
  PreRNAseqBallgownProcess(path.prefix = path.prefix, sample.pattern = sample.pattern)
  if (file.exists(paste0(path.prefix, "Rscript_out/Raw_Read_Process.Rout"))) {
    Hisat2ReportAssemble(path.prefix, gene.name, sample.pattern)
  }
  BallgownPreprocess(path.prefix, gene.name, sample.pattern, independent.variable, ballgown.log2FC, ballgown.qval)
  BallgownPlotAll(path.prefix, independent.variable, ballgown.log2FC, ballgown.qval)
  PostRNAseqBallgownProcess(path.prefix = path.prefix, sample.pattern = sample.pattern)
}


PreRNAseqBallgownProcess <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNAseqBallgownProcess()' environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNAseqBallgownProcess() environment ERROR")
  }
  cat("     (\u2714) : RNAseqBallgownProcess() pre-check is valid\n\n")
}

PostRNAseqBallgownProcess <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNAseqBallgownProcess()' environment post-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNAseqBallgownProcess() post-check ERROR")
  }
  cat("     (\u2714) : RNAseqBallgownProcess() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}
