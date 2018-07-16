#' @export
RNAseqBallgownProcess_CMD <- function(RNASeqWorkFlowParam, ballgown.log2FC = 1, ballgown.pval = 0.05, ballgown.qval = 0.05) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  gene.name <- RNASeqWorkFlowParam@gene.name
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  independent.variable <- RNASeqWorkFlowParam@independent.variable

  fileConn<-file(paste0(path.prefix, "Rscript/Ballgown_Process.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqBallgownProcess(path.prefix = '", path.prefix, "', gene.name = '", gene.name, "', sample.pattern = '", sample.pattern, "', independent.variable = '",independent.variable, "', ballgown.log2FC = ", ballgown.log2FC, ", ballgown.pval = ", ballgown.pval, ", ballgown.qval = ", ballgown.qval, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Ballgown_Process.R' has been created.\n"))
  system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Ballgown_Process.R ", path.prefix, "Rscript_out/Ballgown_Process.Rout"), stdout = "", wait = FALSE)
  cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Ballgown_Process.Rout'\n\n"))
}

#'
#' @export
RNAseqBallgownProcess <- function(path.prefix, gene.name, sample.pattern, independent.variable, ballgown.log2FC, ballgown.pval, ballgown.qval) {
  Hisat2ReportAssemble(path.prefix, gene.name, sample.pattern)
  BallgownPreprocess(path.prefix, gene.name, sample.pattern, independent.variable)
  BallgownPlotAll(path.prefix, ballgown.log2FC, ballgown.pval, ballgown.qval)
  CheckBallgownResult()
}

#' inner function
CheckBallgownResult <- function() {
  if (TRUE) {
    cat("\n")
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
    cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  }
}
