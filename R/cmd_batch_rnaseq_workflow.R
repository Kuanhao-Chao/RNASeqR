#' #'
#' #'@export
#' RNASeqWorkFlow <- function(RNASeqWorkFlowParam, num.parallel.threads = 1, trimming.score = 30, ballgown.log2FC = 1, ballgown.qval = 0.05, run = TRUE, check.s4.print = TRUE) {
#'   # check input param
#'   CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
#'   CheckOperatingSystem(FALSE)
#'   RNASeqEnvironmentSet_CMD(RNASeqWorkFlowParam, run = FALSE, check.s4.print = FALSE)
#'   RNASeqQualityAssessment_CMD(RNASeqWorkFlowParam, run = FALSE, check.s4.print = FALSE)
#'   RNASeqQualityTrimming_CMD(RNASeqWorkFlowParam, run = FALSE, check.s4.print = FALSE)
#'   RNASeqReadProcess_CMD(RNASeqWorkFlowParam, num.parallel.threads = 1, run = FALSE, check.s4.print = FALSE)
#'   RNASeqBallgownProcess_CMD(RNASeqWorkFlowParam, ballgown.log2FC = 1, ballgown.qval = 0.05, run = FALSE, check.s4.print = FALSE)
#'   fileConn<-file(paste0(path.prefix, "Rscript/RNASeq_WORKFLOW.R"))
#'   first <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Environment_Set.R ", path.prefix, "Rscript_out/Environment_Set.Rout"), "\", stdout = \"\", wait = TRUE )")
#'   second <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Trimming.R ", path.prefix, "Rscript_out/Quality_Trimming.Rout"), "\", stdout = \"\", wait = TRUE )")
#'   third <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Quality_Assessment.R ", path.prefix, "Rscript_out/Quality_Assessment.Rout"), "\", stdout = \"\", wait = TRUE )")
#'   fourth <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Read_Process.R ", path.prefix, "Rscript_out/Read_Process.Rout"), "\", stdout = \"\", wait = TRUE )")
#'   fifth <- paste0("system2(command = 'nohup', args = \"", paste0("R CMD BATCH ", path.prefix, "Rscript/Ballgown_Process.R ", path.prefix, "Rscript_out/Ballgown_Process.Rout"), "\", stdout = \"\", wait = TRUE )")
#'   writeLines(c(first, second, third, fourth, fifth), fileConn)
#'   close(fileConn)
#'   cat(paste0("\u2605 '", path.prefix, "Rscript/RNASeq_WORKFLOW.R' has been created.\n"))
#'   if (run) {
#'     system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/RNASeq_WORKFLOW.R ", path.prefix, "Rscript_out/RNASeq_WORKFLOW.Rout"), stdout = "", wait = FALSE )
#'     cat(paste0("\u2605 RNASeq workflow is running in the background. Check current progress in '", path.prefix, "Rscript_out/'\n\n"))
#'   }
#' }
