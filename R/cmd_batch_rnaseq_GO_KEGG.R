#'
#' @export
RNAseqGoKegg_CMD <- function(RNASeqWorkFlowParam, OrgDb.species, KEGG.organism, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  independent.variable <- RNASeqWorkFlowParam@independent.variable
  fileConn<-file(paste0(path.prefix, "Rscript/GO_KEGG_Analysis.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqGoKegg(path.prefix = '", path.prefix, "', independent.variable = '", independent.variable, "', OrgDb.species = '", OrgDb.species, "', KEGG.organism = '",KEGG.organism, "')")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/GO_KEGG_Analysis.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/GO_KEGG_Analysis.R ", path.prefix, "Rscript_out/GO_KEGG_Analysis.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/GO_KEGG_Analysis.Rout'\n\n"))
  }
}

#'
#' @export
RNAseqGoKegg <- function(path.prefix, independent.variable, OrgDb.species, KEGG.organism) {
  CheckOperatingSystem(FALSE)
  PreRNAseqGoKegg()
  which.analyses <- c("ballgown_analysis", "DESeq2_analysis", "edgeR_analysis")
  cat("\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618  Start 'Gene Ontology', 'Kyoto Encyclopedia of Genes and Genomes' analyses  \u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618\n")
  for(which.analysis in which.analyses) {
    cat("\u2618\u2618 ", strsplit(which.analysis, "_")[[1]][1] , " analysis ...\n")
    GOAnalysis(which.analysis, path.prefix, OrgDb.species)
    KEGGAnalysis(which.analysis, path.prefix, OrgDb.species, KEGG.organism)
  }
  PostRNAseqGoKegg()
}

PreRNAseqGoKegg <- function() {
  cat("\u269C\u265C\u265C\u265C RNAseqGoKegg()' environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNAseqGoKegg() environment ERROR")
  }
  cat("(\u2714) : RNAseqGoKegg() pre-check is valid\n\n")
}

PostRNAseqGoKegg <- function() {
  cat("\u269C\u265C\u265C\u265C RNAseqGoKegg()' environment post-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNAseqGoKegg() post-check ERROR")
  }
  cat("(\u2714) : RNAseqGoKegg() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}
