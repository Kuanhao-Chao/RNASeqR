#'
#'@export
RNAseqWorkFlow <- function(RNASeqWorkFlowParam, num.parallel.threads = 8) {
  # check input param
  os.type <- RNASeqWorkFlowParam@os.type
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  gene.name <- RNASeqWorkFlowParam@gene.name
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  indexes.optional <- RNASeqWorkFlowParam@indexes.optional
  MkdirAll(path.prefix)
  r_script.dir <- dir.create(file.path(paste0(path.prefix, 'Rscript/')), showWarnings = FALSE) == 0
  r_script.out.dir <- dir.create(file.path(paste0(path.prefix, 'Rscript_out/')), showWarnings = FALSE) == 0
  fileConn<-file(paste0(path.prefix, "Rscript/RNAseqEnvironmentSet.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqEnvironmentSet(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', gene.name = '", gene.name, "', sample.pattern = '", sample.pattern, "', indexes.optional = ",indexes.optional, ", os.type = '", os.type, "')")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/RNAseqEnvironmentSet.Rout'\n\n"))
  fileConn<-file(paste0(path.prefix, "Rscript/Quality_Control.R"))
  first <- "library(RNASeqWorkflow)"
  second <- "library(ggplot2)"
  third <- paste0("QualityControlRqc(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', sample.pattern = '", sample.pattern, "')")
  writeLines(c(first, second, third), fileConn)
  close(fileConn)
  cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Quality_Control.Rout'\n\n"))
  # If precheck doesn't have .ht2 files is fine
  fileConn<-file(paste0(path.prefix, "Rscript/RNASEQ_RAW_READ_PROCESS.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0('RNAseqRawReadProcess(path.prefix = "', path.prefix, '", input.path.prefix = "', input.path.prefix, '", gene.name = "', gene.name, '", sample.pattern = "', sample.pattern, '", num.parallel.threads = ', num.parallel.threads, ', indexes.optional = ', indexes.optional, ')')
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 RNAseq alignment, assembly, mergence, comparison, reads preprocess are doing in the background. Check current progress in '", path.prefix, "Rscript_out/RNASEQ_PIPELINE.Rout'\n\n"))
  fileConn<-file(paste0(path.prefix, "Rscript/RNASEQ_WORKFLOW.R"))
  first <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/RNAseqEnvironmentSet.R ", path.prefix, "Rscript_out/RNAseqEnvironmentSet.Rout"), "\", stdout = \"\", wait = TRUE )")
  second <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Control.R ", path.prefix, "Rscript_out/Quality_Control.Rout"), "\", stdout = \"\", wait = TRUE )")
  third <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/RNASEQ_RAW_READ_PROCESS.R ", path.prefix, "Rscript_out/RNASEQ_RAW_READ_PROCESS.Rout"), "\", stdout = \"\", wait = TRUE )")
  writeLines(c(first, second, third), fileConn)
  close(fileConn)
  system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/RNASEQ_WORKFLOW.R ", path.prefix, "Rscript_out/RNASEQ_WORKFLOW.Rout"), stdout = "", wait = FALSE )
}
