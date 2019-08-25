#' @title All_Steps_Interface_CMD
#'
#' @description
#'   A functios to run all the steps with in one function. This function execute
#'   in the background:\cr
#'   \enumerate{
#'     \item Create file directories.\cr
#'     \item Install necessary tools. \cr
#'     \item Export 'RNASeq_bin/' to the R environment. \cr
#'     \item Check command of tools. \cr
#'   }
#'   First it will create 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/',
#'   'Rscript/', 'Rscript_out/' directories. \cr Afterwards, 'Hisat2',
#'   'Stringtie', 'Gffcompare' will be installed under
#'   'RNASeq_bin/Download/' and be unpacked under 'RNASeq_bin/Unpacked/'. \cr
#'   'RNASeq_bin/' will be added to the R environment and
#'   validity of tools will be checked.\cr
#'   Any ERROR occurs will be reported and the program will be terminated.\cr
#'   If you want to set up the environment for the following RNA-Seq workflow
#'   in R shell, please see \code{RNASeqEnvironmentSet()} function.
#'
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param RNASeqQualityAssessment.RUN Default value is \code{TRUE}. Set \code{FALSE}
#' to skip "Quality Assessment" step.
#' @param RNASeqReadProcess.RUN Default value is \code{TRUE}. Set \code{FALSE}
#' to skip "RNASeq Read Process" step.
#' @param RNASeqDifferentialAnalysis.RUN Default value is \code{TRUE}. Set \code{FALSE}
#' to skip "RNASeq Differential Analysis" step.
#' @param RNASeqGoKegg.RUN Default value is \code{TRUE}. Set \code{FALSE}
#' to skip "RNASeq Go & Kegg" step.
#' @param run Default value is \code{TRUE}. If \code{TRUE},
#'   'Rscript/Environment_Set.R' will be created and executed.
#'   The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{False}, 'Rscript/Environment_Set.R' will be
#'   created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE},
#'   the result of checking \code{RNASeqRParam} will be reported in
#'   'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqEnvironmentSet_CMD(yeast)}
All_Steps_Interface_CMD <- function(RNASeqRParam,
                                    RNASeqQualityAssessment.RUN    = TRUE,
                                    RNASeqReadProcess.RUN          = TRUE,
                                    RNASeqDifferentialAnalysis.RUN = TRUE,
                                    RNASeqGoKegg.RUN               = TRUE,
                                    OrgDb.species,
                                    go.level = 3,
                                    input.TYPE.ID,
                                    KEGG.organism,
                                    run                = TRUE,
                                    check.s4.print     = TRUE) {
  CheckS4Object_All(RNASeqRParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  INSIDE.path.prefix <- "@"(RNASeqRParam, path.prefix)
  saveRDS(RNASeqRParam,
          file = paste0(INSIDE.path.prefix,
                        "gene_data/RNASeqRParam.rds"))
  fileConn <- file(paste0(path.prefix, "Rscript/All_Step_Analysis.R"))
  first <- "library(RNASeqR)"
# This should be changed
  second <- paste0("All_Steps_Interface(RNASeqRParam = 'INSIDE'",
                   ", which.trigger = 'INSIDE'",
                   ", INSIDE.path.prefix = '", INSIDE.path.prefix,
                   "', RNASeqQualityAssessment.RUN = ", RNASeqQualityAssessment.RUN,
                   ", RNASeqReadProcess.RUN = ", RNASeqReadProcess.RUN,
                   ", RNASeqDifferentialAnalysis.RUN = ", RNASeqDifferentialAnalysis.RUN,
                   ", RNASeqGoKegg.RUN = ", RNASeqGoKegg.RUN,
                   ", OrgDb.species = ", OrgDb.species,
                   ", go.level = ", go.level,
                   ", input.TYPE.ID = ", input.TYPE.ID,
                   ", KEGG.organism = ", KEGG.organism,")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message("\u2605 '", path.prefix,
          "Rscript/All_Step_Analysis.R' has been created.\n")
  if (run) {
    R.home.lib <- R.home()
    R.home.bin <- gsub("/lib/R", "/bin/R", R.home.lib)
    system2(command = "nohup",
            args = paste0(R.home.bin, " CMD BATCH ",
                          path.prefix,
                          "Rscript/All_Step_Analysis.R ",
                          path.prefix,
                          "Rscript_out/All_Step_Analysis.Rout"),
            stdout = "", wait = FALSE)
    message("\u2605 Tools are installing in the background. ",
            "Check current progress in '",
            path.prefix, "Rscript_out/Differential_Analysis.Rout'\n\n")
  }
}

#' @title RNASeqEnvironmentSet
#'
#' @description
#'   Set up the environment for the following RNA-Seq workflow in R shell\cr
#'   This function do 4 things :\cr
#'   \enumerate{
#'     \item Create file directories.\cr
#'     \item Install necessary tools. \cr
#'     \item Export 'RNASeq_bin/' to the R environment. \cr
#'     \item Check command of tools. \cr
#'   }
#'   First it will create 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/',
#'   'Rscript/', 'Rscript_out/' directories. \cr Afterwards, 'Hisat2',
#'   'Stringtie', 'Gffcompare' will be installed under
#'   'RNASeq_bin/Download/' and be unpacked under 'RNASeq_bin/Unpacked/'. \cr
#'   'RNASeq_bin/' will be added to the R environment and
#'   validity of tools will be checked.\cr
#'   Any ERROR occurs will be reported and the program will be terminated.\cr
#'   If you want to set up the environment for the following RNA-Seq workflow
#'   in background, please see \code{RNASeqEnvironmentSet_CMD()} function.
#'
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param which.trigger Default value is \code{OUTSIDE}. User should not change
#'   this value.
#' @param INSIDE.path.prefix Default value is \code{NA}. User should not change
#'   this value.
#' @param RNASeqQualityAssessment.RUN Default value is \code{TRUE}. Set \code{FALSE}
#' to skip "Quality Assessment" step.
#' @param RNASeqReadProcess.RUN Default value is \code{TRUE}. Set \code{FALSE}
#' to skip "RNASeq Read Process" step.
#' @param RNASeqDifferentialAnalysis.RUN Default value is \code{TRUE}. Set \code{FALSE}
#' to skip "RNASeq Differential Analysis" step.
#' @param RNASeqGoKegg.RUN Default value is \code{TRUE}. Set \code{FALSE}
#' to skip "RNASeq Go & Kegg" step.
#' @param OrgDb.species the genome wide annotation packages of species on
#'   Bioconductor. Currently, there are 19 supported genome wide annotation
#'   packages of species.
#' @param go.level the depth of acyclic graph in GO analysis
#' @param input.TYPE.ID The gene name type in OrgDb.species annotation packahge.
#' @param KEGG.organism the species that are supported for KEGG analysis.
#'   Currently, there are more than 5000 supported species genome.
#'   Check the valid species terms on
#'   https://www.genome.jp/kegg/catalog/org_list.html
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE},
#'   the result of checking \code{RNASeqRParam} will be reported in
#'   'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqEnvironmentSet(RNASeqRParam = yeast)}
All_Steps_Interface <- function(RNASeqRParam,
                                which.trigger      = "OUTSIDE",
                                INSIDE.path.prefix = NA,
                                RNASeqQualityAssessment.RUN    = TRUE,
                                RNASeqReadProcess.RUN          = TRUE,
                                RNASeqDifferentialAnalysis.RUN = TRUE,
                                RNASeqGoKegg.RUN               = TRUE,
                                OrgDb.species,
                                go.level = 3,
                                input.TYPE.ID,
                                KEGG.organism,
                                check.s4.print     = TRUE) {
  CheckOperatingSystem(FALSE)
  # If `which.trigger` is OUTSIDE, then directory must be built
  # If `which.trigger` is INSIDE, then directory must not be
  #  built here(will created in CMD)
  if (isS4(RNASeqRParam) &
      which.trigger == "OUTSIDE" &
      is.na(INSIDE.path.prefix)) {
    # This is an external call!!
    # Check the S4 object(user input)
    CheckS4Object_All(RNASeqRParam, check.s4.print)
  } else if (RNASeqRParam == "INSIDE" &
             which.trigger == "INSIDE" &
             !is.na(INSIDE.path.prefix)) {
    # This is an internal call!!
    # Load the S4 object that saved in CMD process
    RNASeqRParam <- readRDS(paste0(INSIDE.path.prefix,
                                   "gene_data/RNASeqRParam.rds"))
  }
  if (RNASeqQualityAssessment.RUN) {
    RNASeqQualityAssessment(RNASeqRParam)
  }
  if (RNASeqReadProcess.RUN) {
    RNASeqReadProcess(RNASeqRParam)
  }
  if (RNASeqDifferentialAnalysis.RUN) {
    RNASeqDifferentialAnalysis(RNASeqRParam)
  }
  if (RNASeqGoKegg.RUN) {
    RNASeqGoKegg(RNASeqRParam,
                 OrgDb.species = OrgDb.species,
                 go.level = 3,
                 input.TYPE.ID = input.TYPE.ID,
                 KEGG.organism = KEGG.organism)
  }
}
