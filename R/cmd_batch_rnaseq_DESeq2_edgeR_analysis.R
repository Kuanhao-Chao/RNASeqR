#'
#' @export
RNASeqDESeq2edgeR_CMD <- function(RNASeqWorkFlowParam, DESeq2.padj = 0.01, DESeq2.log2FC = 1, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix = RNASeqWorkFlowParam@path.prefix
  independent.variable = RNASeqWorkFlowParam@independent.variable
  control.group = RNASeqWorkFlowParam@control.group
  experiment.group = RNASeqWorkFlowParam@experiment.group
  fileConn<-file(paste0(path.prefix, "Rscript/DESeq2_edgeR_Process.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNASeqDESeq2edgeR(path.prefix = '", path.prefix, "', independent.variable = '", independent.variable, "', control.group = '", control.group, "', experiment.group = '",experiment.group,  "', DESeq2.padj = ",DESeq2.padj, ", DESeq2.log2FC = ",DESeq2.log2FC, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/DESeq2_edgeR_Process.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/DESeq2_edgeR_Process.R ", path.prefix, "Rscript_out/DESeq2_edgeR_Process.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/DESeq2_edgeR_Process.Rout'\n\n"))
  }
}

#'
#' @export
RNASeqDESeq2edgeR <- function(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC) {
  CheckOperatingSystem(FALSE)
  PreRNASeqDESeq2edgeR()
  DESeq2edgeRRawCountAnalysis(path.prefix = path.prefix, independent.variable = independent.variable, control.group = control.group, experiment.group = experiment.group, DESeq2.padj = DESeq2.padj, DESeq2.log2FC = DESeq2.log2FC)
  PostRNASeqDESeq2edgeR()
}

PreRNASeqDESeq2edgeR <- function() {
  cat("\u269C\u265C\u265C\u265C RNASeqBallgownProcess()' environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqBallgownProcess() environment ERROR")
  }
  cat("(\u2714) : RNASeqBallgownProcess() pre-check is valid\n\n")
}

PostRNASeqDESeq2edgeR <- function() {
  cat("\u269C\u265C\u265C\u265C RNASeqBallgownProcess()' environment post-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqBallgownProcess() post-check ERROR")
  }
  cat("(\u2714) : RNASeqBallgownProcess() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}
