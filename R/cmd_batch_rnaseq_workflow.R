#'
#'@export
RNAseqWorkFlow <- function(RNASeqWorkFlowParam, num.parallel.threads = 8, trimming.score = 30, ballgown.log2FC = 1, ballgown.pval = 0.05, ballgown.qval = 0.05) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam)
  os.type <- RNASeqWorkFlowParam@os.type
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  gene.name <- RNASeqWorkFlowParam@gene.name
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  indexes.optional <- RNASeqWorkFlowParam@indexes.optional
  independent.variable <- RNASeqWorkFlowParam@independent.variable
  python.variable <- RNASeqWorkFlowParam@python.variable
  python.variable.answer <- python.variable$check.answer
  python.variable.version <- python.variable$python.version
  MkdirAll(path.prefix)
  fileConn<-file(paste0(path.prefix, "Rscript/Environment_Set.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqEnvironmentSet(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', gene.name = '", gene.name, "', sample.pattern = '", sample.pattern, "', indexes.optional = ",indexes.optional, ", os.type = '", os.type, "')")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Environment_Set.R' has been created.\n"))
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Trimming.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqQualityTrimming(path.prefix = '", path.prefix, "', sample.pattern = '", sample.pattern, "', trimming.score = ", trimming.score, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Quality_Trimming.R' has been created.\n"))
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Assessmnet.R"))
  first <- "library(RNASeqWorkflow)"
  second <- "library(ggplot2)"
  third <- paste0("RNASeqQualityAssessment(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', sample.pattern = '", sample.pattern, "')")
  writeLines(c(first, second, third), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Quality_Assessmnet.R' has been created.\n"))
  # If precheck doesn't have .ht2 files is fine
  fileConn<-file(paste0(path.prefix, "Rscript/Raw_Read_Process.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0('RNAseqRawReadProcess(path.prefix = "', path.prefix, '", input.path.prefix = "', input.path.prefix, '", gene.name = "', gene.name, '", sample.pattern = "', sample.pattern, '", python.variable.answer = ', python.variable.answer, ', python.variable.version = ', python.variable.version, ', num.parallel.threads = ', num.parallel.threads, ', indexes.optional = ', indexes.optional, ')')
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Raw_Read_Process.R' has been created.\n"))
  fileConn<-file(paste0(path.prefix, "Rscript/Ballgown_Process.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqBallgownProcess(path.prefix = '", path.prefix, "', gene.name = '", gene.name, "', sample.pattern = '", sample.pattern, "', independent.variable = '",independent.variable, "', ballgown.log2FC = ", ballgown.log2FC, ", ballgown.pval = ", ballgown.pval, ", ballgown.qval = ", ballgown.qval, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Ballgown_Process.R' has been created.\n"))
  fileConn<-file(paste0(path.prefix, "Rscript/RNASEQ_WORKFLOW.R"))
  first <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Environment_Set.R ", path.prefix, "Rscript_out/Environment_Set.Rout"), "\", stdout = \"\", wait = TRUE )")
  second <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Control.R ", path.prefix, "Rscript_out/Quality_Control.Rout"), "\", stdout = \"\", wait = TRUE )")
  third <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Raw_Read_Process.R ", path.prefix, "Rscript_out/Raw_Read_Process.Rout"), "\", stdout = \"\", wait = TRUE )")
  fourth <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Ballgown_Process.R ", path.prefix, "Rscript_out/Ballgown_Process.Rout"), "\", stdout = \"\", wait = TRUE )")
  writeLines(c(first, second, third, fourth), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/RNASEQ_WORKFLOW.R' has been created.\n"))
  system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/RNASEQ_WORKFLOW.R ", path.prefix, "Rscript_out/RNASEQ_WORKFLOW.Rout"), stdout = "", wait = FALSE )
  cat(paste0("\u2605 RNAseq workflow is running in the background. Check current progress in '", path.prefix, "Rscript_out/'\n\n"))
}
