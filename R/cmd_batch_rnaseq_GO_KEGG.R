#' @title GO functional analysis and KEGG pathway analysis for RNA-Seq workflow in background
#'
#' @description Run Gene Ontology(GO) and Kyoto Encyclopedia of Genes and Genomes(KEGG) analysis in background.
#' This function do Gene Ontology(GO) and Kyoto Encyclopedia of Genes and Genomes(KEGG) analysis :
#' 1. Gene Ontology(GO) :
#'    All gene set :
#'        1. Do GO gene set enrichment analysis.
#'    Differential expressed gene set :
#'        1. Do GO function classification analysis.
#'        2. Do GO function enrichment analysis.
#'        3. Visualization : bar plot, dot plot etc.
#' 2. Kyoto Encyclopedia of Genes and Genomes(KEGG) :
#'    All gene set :
#'        1. Do KEGG gene set enrichment analysis.
#'    Differential expressed gene set :
#'        1. Do KEGG pathway gene set enrichment analysis
#'        2. Do KEGG pathway enrichment analysis
#'        3. Pathway visulization with \code{pathview} package. KEGG webpage pathway url will also be created
#' If you want to do GO functional analysis and KEGG pathway analysis for the following RNA-Seq workflow in R shell, please see \code{RNASeqGoKegg()} function.
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related parameters
#' @param OrgDb.species the genome wide annotation packages of species on Bioconductor. Currently, there are 19 supported genome wide annotation packages of species.
#' @param KEGG.organism the species that are supported for KEGG analysis. Currently, there are more than 5000 supported species genome. Check the valid species terms on https://www.genome.jp/kegg/catalog/org_list.html
#' @param run Default value is \code{TRUE}. If \code{TRUE}, 'Rscript/Environment_Set.R' will be created and executed. The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#' If \code{False}, 'Rscript/Environment_Set.R' will be created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of checking \code{RNASeqWorkFlowParam} will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking \code{RNASeqWorkFlowParam} will not be in 'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNASeqEnvironmentSet_CMD(RNASeqWorkFlowParam = exp)}
RNASeqGoKegg_CMD <- function(RNASeqWorkFlowParam, OrgDb.species, go.level = 3, input.TYPE.ID, KEGG.organism, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  independent.variable <- RNASeqWorkFlowParam@independent.variable
  fileConn<-file(paste0(path.prefix, "Rscript/GO_KEGG_Analysis.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNASeqGoKegg(path.prefix = '", path.prefix, "', independent.variable = '", independent.variable, "', OrgDb.species = '", OrgDb.species, "', go.level = ", go.level, ", input.TYPE.ID = '", input.TYPE.ID, "', KEGG.organism = '",KEGG.organism, "')")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message(paste0("\u2605 '", path.prefix, "Rscript/GO_KEGG_Analysis.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/GO_KEGG_Analysis.R ", path.prefix, "Rscript_out/GO_KEGG_Analysis.Rout"), stdout = "", wait = FALSE)
    message(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/GO_KEGG_Analysis.Rout'\n\n"))
  }
}

#' @title GO functional analysis and KEGG pathway analysis for RNA-Seq workflow in R shell
#'
#' @description Run Gene Ontology(GO) and Kyoto Encyclopedia of Genes and Genomes(KEGG) analysis in background.
#' It is strongly advised to run \code{RNASeqGoKegg_CMD()} directly. Running \code{RNASeqGoKegg_CMD()} will create 'GO_KEGG_Analysis.Rout' file in 'Rscript_out/' directory.
#' This function do Gene Ontology(GO) and Kyoto Encyclopedia of Genes and Genomes(KEGG) analysis :
#' 1. Gene Ontology(GO) :
#'    All gene set :
#'        1. Do GO gene set enrichment analysis.
#'    Differential expressed gene set :
#'        1. Do GO function classification analysis.
#'        2. Do GO function enrichment analysis.
#'        3. Visualization : bar plot, dot plot etc.
#' 2. Kyoto Encyclopedia of Genes and Genomes(KEGG) :
#'    All gene set :
#'        1. Do KEGG gene set enrichment analysis.
#'    Differential expressed gene set :
#'        1. Do KEGG pathway gene set enrichment analysis
#'        2. Do KEGG pathway enrichment analysis
#'        3. Pathway visulization with \code{pathview} package. KEGG webpage pathway url will also be created
#' If you want to do GO functional analysis and KEGG pathway analysis for the following RNA-Seq workflow in background, please see \code{RNASeqGoKegg_CMD()} function.
#' @param path.prefix path prefix of 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @param independent.variable independent variable for the biological experiment design of two-group RNA-Seq workflow
#' @param OrgDb.species the genome wide annotation packages of species on Bioconductor. Currently, there are 19 supported genome wide annotation packages of species.
#' @param KEGG.organism the species that are supported for KEGG analysis. Currently, there are more than 5000 supported species genome. Check the valid species terms on https://www.genome.jp/kegg/catalog/org_list.html
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNASeqEnvironmentSet_CMD(RNASeqWorkFlowParam = exp)}
RNASeqGoKegg <- function(path.prefix, independent.variable, OrgDb.species, go.level = 3, input.TYPE.ID, KEGG.organism) {
  CheckOperatingSystem(FALSE)
  PreRNASeqGoKegg()
  raw.read.avail <- RawReadCountAvailability(path.prefix)
  if (raw.read.avail) {
    which.analyses <- c("ballgown_analysis", "TPM_analysis", "DESeq2_analysis", "edgeR_analysis")
  } else {
    which.analyses <- c("ballgown_analysis", "TPM_analysis")
  }
  message("\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618  Start 'Gene Ontology', 'Kyoto Encyclopedia of Genes and Genomes' analyses  \u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618\n")
  for(which.analysis in which.analyses) {
    message("\u2618\u2618 ", strsplit(which.analysis, "_")[[1]][1] , " analysis ...\n")
    GOAnalysis(which.analysis, path.prefix, OrgDb.species, go.level, input.TYPE.ID)
    KEGGAnalysis(which.analysis, path.prefix, OrgDb.species, input.TYPE.ID, KEGG.organism)
  }
  PostRNASeqGoKegg()
}

PreRNASeqGoKegg <- function() {
  message("\u269C\u265C\u265C\u265C RNASeqGoKegg()' environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqGoKegg() environment ERROR")
  }
  message("(\u2714) : RNASeqGoKegg() pre-check is valid\n\n")
}

PostRNASeqGoKegg <- function() {
  message("\u269C\u265C\u265C\u265C RNASeqGoKegg()' environment post-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqGoKegg() post-check ERROR")
  }
  message("(\u2714) : RNASeqGoKegg() post-check is valid\n\n")
  message(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  message(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  message(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}
