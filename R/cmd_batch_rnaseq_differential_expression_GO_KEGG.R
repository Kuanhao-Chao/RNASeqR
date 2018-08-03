#'
#' @export
RNAseqDeGoKegg_CMD <- function(RNASeqWorkFlowParam, OrgDb.species, KEGG.organism, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  independent.variable <- RNASeqWorkFlowParam@independent.variable
  fileConn<-file(paste0(path.prefix, "Rscript/Differential_Expression_GO_KEGG.R"))
  first <- "library(RNASeqWorkflow)"
  second <- "library(pathview)"
  third <- paste0("RNAseqDeGoKegg(path.prefix = '", path.prefix, "', independent.variable = '", independent.variable, "', OrgDb.species = '", OrgDb.species, "', KEGG.organism = '",KEGG.organism, "')")
  writeLines(c(first, second, third), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Differential_Expression_GO_KEGG.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Differential_Expression_GO_KEGG.R ", path.prefix, "Rscript_out/Differential_Expression_GO_KEGG.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Differential_Expression_GO_KEGG.Rout'\n\n"))
  }
}

#'
#' @export
RNAseqDeGoKegg <- function(path.prefix, independent.variable, OrgDb.species, KEGG.organism) {
  CheckOperatingSystem(FALSE)
  DEBallgownPlotAll(path.prefix, independent.variable)
  DEGOAnalysis(path.prefix, OrgDb.species)
  DEKEGGAnalysis(path.prefix, OrgDb.species, KEGG.organism)
  PostRNAseqDeGoKegg()
}

PreRNAseqDeGoKegg <- function() {
  cat("\u269C\u265C\u265C\u265C RNAseqBallgownProcess()' environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNAseqBallgownProcess() environment ERROR")
  }
  cat("(\u2714) : RNAseqBallgownProcess() pre-check is valid\n\n")
}

PostRNAseqDeGoKegg <- function() {
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
