#' @title Ballgown analysis for RNA-Seq workflow in background
#'
#' @description Use ballgown R package for statistical analysis of assembled transcriptomes, including flexible differential expression analysis, and sample FPKM visualization as well as pre differential analysis visualization for the following RNA-Seq workflow in background.
#' This function do 3 things :
#' 1. Create ballgown object, and write "ballgown_FPKM_result.csv" file.
#' 2. Find the differential expressed genes, and write "ballgown_FPKM_DE_result.csv"
#' 3. Sample data visualization.
#'    For transcript-related : Creating transcript_count_per_gene_plot, transcript_length_plot.
#'    For sample profile : Creating Box plot, Violin plot, Frequency plot, PCA related plots, Correlation plots.
#'    For pre differential expressed gene analysis : Creating MA plot, Volcano plot.
#'    For differential express gene analysis : Creating Heatmap plot, PCA related plots.
#' If you want to run ballgown analysis for the following RNA-Seq workflow in R shell, please see \code{RNASeqBallgownDESeq2EdgeRProcess()} function.
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related parameters
#' @param ballgown.log2FC Default \code{1}. Set the threshold of log2 fold change to filter out differential expressed gene.
#' @param ballgown.qval Default \code{0.05}. Set the threshold of q-value to filter out differential expressed gene.
#' @param run Default value is \code{TRUE}. If \code{TRUE}, 'Rscript/Environment_Set.R' will be created and executed. The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#' If \code{False}, 'Rscript/Environment_Set.R' will be created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of checking \code{RNASeqWorkFlowParam} will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking \code{RNASeqWorkFlowParam} will not be in 'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNASeqEnvironmentSet_CMD(RNASeqWorkFlowParam = exp)}
RNASeqBallgownDESeq2EdgeRProcess_CMD <- function(RNASeqWorkFlowParam, ballgown.qval = 0.05, ballgown.log2FC = 1, DESeq2.padj = 0.1, DESeq2.log2FC = 1, edgeR.pval = 0.05, edgeR.log2FC = 1, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  genome.name <- RNASeqWorkFlowParam@genome.name
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  independent.variable <- RNASeqWorkFlowParam@independent.variable
  control.group <- RNASeqWorkFlowParam@control.group
  experiment.group <- RNASeqWorkFlowParam@experiment.group
  fileConn<-file(paste0(path.prefix, "Rscript/Ballgown_DESeq2_edgeR_Process.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNASeqBallgownDESeq2EdgeRProcess(path.prefix = '", path.prefix, "', genome.name = '", genome.name, "', sample.pattern = '", sample.pattern, "', independent.variable = '",independent.variable,  "', control.group = '",control.group,  "', experiment.group = '",experiment.group, "', ballgown.log2FC = ", ballgown.log2FC, ", ballgown.qval = ", ballgown.qval, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Ballgown_DESeq2_edgeR_Process.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Ballgown_DESeq2_edgeR_Process.R ", path.prefix, "Rscript_out/Ballgown_DESeq2_edgeR_Process.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Ballgown_DESeq2_edgeR_Process.Rout'\n\n"))
  }
}

#' @title Ballgown analysis for RNA-Seq workflow in R shell
#'
#' @description Use ballgown R package for statistical analysis of assembled transcriptomes, including flexible differential expression analysis, and sample FPKM visualization as well as pre differential analysis visualization for the following RNA-Seq workflow in background.
#' It is strongly advised to run \code{RNASeqBallgownDESeq2EdgeRProcess_CMD()} directly. Running this function directly is not recommended.
#' This function do 3 things :
#' 1. Create ballgown object, and write "ballgown_FPKM_result.csv" file.
#' 2. Find the differential expressed genes, and write "ballgown_FPKM_DE_result.csv"
#' 3. Sample data visualization.
#'    For transcript-related : Creating transcript_count_per_gene_plot, transcript_length_plot.
#'    For sample profile : Creating Box plot, Violin plot, Frequency plot, PCA related plots, Correlation plots.
#'    For pre differential expressed gene analysis : Creating MA plot, Volcano plot.
#'    For differential express gene analysis : Creating Heatmap plot, PCA related plots.
#' If you want to run ballgown analysis for the following RNA-Seq workflow in background, please see \code{RNASeqBallgownDESeq2EdgeRProcess_CMD()} function.
#'
#' @param path.prefix path prefix of 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @param genome.name variable of genome name defined in this RNA-Seq workflow (ex. genome.name.fa, genome.name.gtf)
#' @param sample.pattern  regular expression of raw fastq.gz files under 'input_files/raw_fastq.gz'
#' @param independent.variable independent variable for the biological experiment design of two-group RNA-Seq workflow
#' @param control.group group name of the control group
#' @param experiment.group group name of the experiment group
#' @param ballgown.log2FC Default \code{1}. Set the threshold of log2 fold change to filter out differential expressed gene.
#' @param ballgown.qval Default \code{0.05}. Set the threshold of q-value to filter out differential expressed gene.
#'
#' @return None
#' @export
#' @examples
#' \dontrun{
#' input_file_dir <- system.file(package = "RNASeqWorkflow", "exdata")
#' exp <- RNASeqWorkFlowParam(path.prefix = "/tmp/", input.path.prefix = input_file_dir, genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNASeqEnvironmentSet_CMD(RNASeqWorkFlowParam = exp)}
RNASeqBallgownDESeq2EdgeRProcess <- function(path.prefix, genome.name, sample.pattern, independent.variable, control.group, experiment.group, ballgown.log2FC = 1, ballgown.qval = 0.05, DESeq2.padj = 0.1, DESeq2.log2FC = 1, edgeR.pval = 0.05, edgeR.log2FC = 1) {
  CheckOperatingSystem(FALSE)
  PreRNASeqBallgownDESeq2EdgeRProcess(path.prefix = path.prefix, sample.pattern = sample.pattern)
  if (file.exists(paste0(path.prefix, "Rscript_out/Raw_Read_Process.Rout"))) {
    Hisat2ReportAssemble(path.prefix, genome.name, sample.pattern)
  }
  cat("\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618  Start 'ballgown', 'DESeq2' 'edgeR' analyses  \u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618\n")
  BallgownAnalysis(path.prefix, genome.name, sample.pattern, independent.variable, control.group, experiment.group, ballgown.qval, ballgown.log2FC)
  DESeq2RawCountAnalysis(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC)
  edgeRRawCountAnalysis(path.prefix, independent.variable, control.group, experiment.group, edgeR.pval, edgeR.log2FC)
  PostRNASeqBallgownDESeq2EdgeRProcess(path.prefix = path.prefix, sample.pattern = sample.pattern)
}


PreRNASeqBallgownDESeq2EdgeRProcess <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNASeqBallgownDESeq2EdgeRProcess()' environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqBallgownDESeq2EdgeRProcess() environment ERROR")
  }
  cat("(\u2714) : RNASeqBallgownDESeq2EdgeRProcess() pre-check is valid\n\n")
}

PostRNASeqBallgownDESeq2EdgeRProcess <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNASeqBallgownDESeq2EdgeRProcess()' environment post-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqBallgownDESeq2EdgeRProcess() post-check ERROR")
  }
  cat("(\u2714) : RNASeqBallgownDESeq2EdgeRProcess() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}
