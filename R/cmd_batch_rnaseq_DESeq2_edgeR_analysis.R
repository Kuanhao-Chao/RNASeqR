#'
#' @export
RNAseqDESeq2edgeR_CMD <- function(RNASeqWorkFlowParam, DESeq2.padj = 0.01, DESeq2.log2FC = 1, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix = RNASeqWorkFlowParam@path.prefix
  independent.variable = RNASeqWorkFlowParam@independent.variable
  control.group = RNASeqWorkFlowParam@control.group
  experiment.group = RNASeqWorkFlowParam@experiment.group
  fileConn<-file(paste0(path.prefix, "Rscript/DESeq2_edgeR_Process.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqDESeq2edgeR(path.prefix = '", path.prefix, "', independent.variable = '", independent.variable, "', control.group = '", control.group, "', experiment.group = '",experiment.group,  "', DESeq2.padj = ",DESeq2.padj, ", DESeq2.log2FC = ",DESeq2.log2FC, ")")
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
RNAseqDESeq2edgeR <- function(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC) {
  CheckOperatingSystem(FALSE)
  PreRNAseqDESeq2edgeR()
  DESeq2edgeRRawCountAnalysis(path.prefix = path.prefix, independent.variable = independent.variable, control.group = control.group, experiment.group = experiment.group, DESeq2.padj = DESeq2.padj, DESeq2.log2FC = DESeq2.log2FC)
  PostRNAseqDESeq2edgeR()
}

PreRNAseqDESeq2edgeR <- function() {
  cat("\u269C\u265C\u265C\u265C RNAseqBallgownProcess()' environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNAseqBallgownProcess() environment ERROR")
  }
  cat("(\u2714) : RNAseqBallgownProcess() pre-check is valid\n\n")
}

PostRNAseqDESeq2edgeR <- function() {
  cat("\u269C\u265C\u265C\u265C RNAseqBallgownProcess()' environment post-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNAseqBallgownProcess() post-check ERROR")
  }
  cat("(\u2714) : RNAseqBallgownProcess() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}
