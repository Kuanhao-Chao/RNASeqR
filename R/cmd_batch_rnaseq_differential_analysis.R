#' @title RNASeqDifferentialAnalysis_CMD
#'
#' @description
#'   This function will run differential analysis on ballgown,
#'   DESeq2 and edgeR in R shell. \cr
#'   This function do following things : \cr
#'   \enumerate{
#'     \item ballgown analysis \cr
#'     Raw reads are normalized into FPKM values \cr
#'     The main statistic test in ballgown is paramatic F-test comparing nested
#'     linear models \cr
#'     \item DESeq2 analysis  \cr
#'      Median of rations normalization(MRN) is used in DESeq2 for raw reads
#'       count normalization.  \cr
#'      Sequencing depth and RNA composition is taken into consideration is this
#'       normalization method.  \cr
#'      The main statistic test in DESeq2 is negative binomial distribution. \cr
#'     \item edgeR analysis  \cr
#'      Raw reads are normalized by TMM and library size.
#'       (run \code{calcNormFactors()} to get a DGEList,
#'      and then run \code{cpm()} on that DGEList)  \cr
#'      The main statistic test in edgeR is trimmed mean of M-values(TMM).\cr
#'   }
#'   If you want to run differential analysis on ballgown,
#'   DESeq2, edgeR for the following RNA-Seq workflow in R shell,
#'   please see \code{RNASeqDifferentialAnalysis()} function.
#'
#' @param RNASeqRParam S4 object instance of
#'   experiment-related parameters
#' @param which.trigger Default value is \code{OUTSIDE}. User should not change
#'   this value.
#' @param INSIDE.path.prefix Default value is \code{NA}. User should not change
#'   this value.
#' @param Pre_DE.visualization Default \code{TRUE}. Whether to visualize pre-DE
#' analysis results.
#' @param Post_DE.visualization Default \code{TRUE}. Whether to visualize
#' post-DE analysis results.
#' @param ballgown.run Default \code{TRUE}. Logical value whether to run
#'   ballgown differential analysis.
#' @param ballgown.pval Default \code{0.05}. Set the threshold of ballgown
#'   p-value to filter out differential expressed gene.
#' @param ballgown.log2FC Default \code{1}. Set the threshold of ballgown
#'   log2 fold change to filter out differential expressed gene.
#' @param DESeq2.run Default \code{TRUE}. Logical value whether to run
#'   DESeq2 differential analysis.
#' @param DESeq2.pval Default \code{0.05}. Set the threshold of DESeq2 p-value
#'   to filter out differential expressed gene.
#' @param DESeq2.log2FC Default \code{1}. Set the threshold of DESeq2 log2
#'   fold change to filter out differential expressed gene.
#' @param edgeR.run Default \code{TRUE}. Logical value whether to run
#'   edgeR differential analysis.
#' @param edgeR.pval Default \code{0.05}. Set the threshold of edgeR p-value
#'   to filter out differential expressed gene.
#' @param edgeR.log2FC Default \code{1}. Set the threshold of edgeR log2
#'   fold change to filter out differential expressed gene.
#' @param run Default value is \code{TRUE}. If \code{TRUE},
#'   'Rscript/Environment_Set.R' will be created and executed. The output log
#'   will be stored in 'Rscript_out/Environment_Set.Rout'. If \code{False},
#'   'Rscript/Environment_Set.R' will be created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of
#'   checking \code{RNASeqRParam} will be reported in
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
#' RNASeqDifferentialAnalysis_CMD(RNASeqRParam = yeast)}
RNASeqDifferentialAnalysis_CMD <- function(RNASeqRParam,
                                           which.trigger      = "OUTSIDE",
                                           INSIDE.path.prefix = NA,
                                           Pre_DE.visualization = TRUE,
                                           Post_DE.visualization = TRUE,
                                           ballgown.run    = TRUE,
                                           ballgown.pval   = 0.05,
                                           ballgown.log2FC = 1,
                                           DESeq2.run      = TRUE,
                                           DESeq2.pval     = 0.1,
                                           DESeq2.log2FC   = 1,
                                           edgeR.run       = TRUE,
                                           edgeR.pval      = 0.05,
                                           edgeR.log2FC    = 1,
                                           run             = TRUE,
                                           check.s4.print  = TRUE) {
  # check input param
  CheckS4Object(RNASeqRParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  INSIDE.path.prefix <- "@"(RNASeqRParam, path.prefix)
  saveRDS(RNASeqRParam,
          file = paste0(INSIDE.path.prefix,
                        "gene_data/RNASeqRParam.rds"))
  fileConn <- file(paste0(path.prefix, "Rscript/Differential_Analysis.R"))
  first <- "library(RNASeqR)"
  second <- paste0("RNASeqDifferentialAnalysis(RNASeqRParam = 'INSIDE'",
                   ", which.trigger = 'INSIDE'",
                   ", INSIDE.path.prefix = '", INSIDE.path.prefix,
                   "', ballgown.run = ", ballgown.run,
                   ", ballgown.pval = ", ballgown.pval,
                   ", ballgown.log2FC = ", ballgown.log2FC,
                   ", DESeq2.run = ", DESeq2.run,
                   ", DESeq2.pval = ", DESeq2.pval,
                   ", DESeq2.log2FC = ", DESeq2.log2FC,
                   ", edgeR.run = ", edgeR.run,
                   ", edgeR.pval = ", edgeR.pval,
                   ", edgeR.log2FC = ", edgeR.log2FC,")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message("\u2605 '", path.prefix,
          "Rscript/Differential_Analysis.R' has been created.\n")
  if (run) {
    R.home.lib <- R.home()
    R.home.bin <- gsub("/lib/R", "/bin/R", R.home.lib)
    system2(command = "nohup",
            args = paste0(R.home.bin, " CMD BATCH ",
                          path.prefix,
                          "Rscript/Differential_Analysis.R ",
                          path.prefix,
                          "Rscript_out/Differential_Analysis.Rout"),
            stdout = "", wait = FALSE)
    message("\u2605 Tools are installing in the background. ",
            "Check current progress in '",
            path.prefix, "Rscript_out/Differential_Analysis.Rout'\n\n")
  }
}

#' @title RNASeqDifferentialAnalysis
#'
#' @description
#'   This function will run differential analysis on ballgown,
#'   DESeq2 and edgeR in background. \cr
#'   This function do following things : \cr
#'   \enumerate{
#'     \item ballgown analysis \cr
#'     Raw reads are normalized into FPKM values \cr
#'     The main statistic test in ballgown is paramatic F-test comparing nested
#'     linear models \cr
#'     \item DESeq2 analysis  \cr
#'      Median of rations normalization(MRN) is used in DESeq2 for raw reads
#'       count normalization.  \cr
#'      Sequencing depth and RNA composition is taken into consideration is this
#'       normalization method.  \cr
#'      The main statistic test in DESeq2 is negative binomial distribution. \cr
#'     \item edgeR analysis  \cr
#'      Raw reads are normalized by TMM and library size.
#'       (run \code{calcNormFactors()} to get a DGEList,
#'      and then run \code{cpm()} on that DGEList)  \cr
#'      The main statistic test in edgeR is trimmed mean of M-values(TMM).\cr
#'   }
#'   If you want to run differential analysis on ballgown,
#'   DESeq2, edgeR for the following RNA-Seq workflow in background,
#'   please see \code{RNASeqDifferentialAnalysis()} function.
#'
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param which.trigger Default value is \code{OUTSIDE}. User should not change
#'   this value.
#' @param INSIDE.path.prefix Default value is \code{NA}. User should not change
#'   this value.
#' @param ballgown.run Default \code{TRUE}. Logical value whether to run
#'   ballgown differential analysis.
#' @param ballgown.pval Default \code{0.05}. Set the threshold of ballgown
#'   p-value to filter out differential expressed gene.
#' @param ballgown.log2FC Default \code{1}. Set the threshold of ballgown
#'   log2 fold change to filter out differential expressed gene.
#' @param DESeq2.run Default \code{TRUE}. Logical value whether to run
#'   DESeq2 differential analysis.
#' @param DESeq2.pval Default \code{0.05}. Set the threshold of DESeq2 p-value
#'   to filter out differential expressed gene.
#' @param DESeq2.log2FC Default \code{1}. Set the threshold of DESeq2 log2
#'   fold change to filter out differential expressed gene.
#' @param edgeR.run Default \code{TRUE}. Logical value whether to run
#'   edgeR differential analysis.
#' @param edgeR.pval Default \code{0.05}. Set the threshold of edgeR p-value
#'   to filter out differential expressed gene.
#' @param edgeR.log2FC Default \code{1}. Set the threshold of edgeR log2
#'   fold change to filter out differential expressed gene.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of
#'   checking \code{RNASeqRParam} will be reported in
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
#' RNASeqDifferentialAnalysis(RNASeqRParam = yeast)}
RNASeqDifferentialAnalysis <- function(RNASeqRParam,
                                       which.trigger      = "OUTSIDE",
                                       INSIDE.path.prefix = NA,
                                       Pre_DE.visualization = TRUE,
                                       Post_DE.visualization = TRUE,
                                       ballgown.run    = TRUE,
                                       ballgown.pval   = 0.05,
                                       ballgown.log2FC = 1,
                                       DESeq2.run      = TRUE,
                                       DESeq2.pval     = 0.1,
                                       DESeq2.log2FC   = 1,
                                       edgeR.run       = TRUE,
                                       edgeR.pval      = 0.05,
                                       edgeR.log2FC    = 1,
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
    CheckS4Object(RNASeqRParam, check.s4.print)
  } else if (RNASeqRParam == "INSIDE" &
             which.trigger == "INSIDE" &
             !is.na(INSIDE.path.prefix)) {
    # This is an internal call!!
    # Load the S4 object that saved in CMD process
    RNASeqRParam <- readRDS(paste0(INSIDE.path.prefix,
                                   "gene_data/RNASeqRParam.rds"))
  }
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  genome.name <- "@"(RNASeqRParam, genome.name)
  sample.pattern <- "@"(RNASeqRParam, sample.pattern)
  independent.variable <- "@"(RNASeqRParam, independent.variable)
  case.group <- "@"(RNASeqRParam, case.group)
  control.group <- "@"(RNASeqRParam, control.group)


  # 1. Pre-DE assessment visualization
  phenoData.result<- phenoDataWrap(path.prefix,
                                   independent.variable,
                                   case.group,
                                   control.group)
  if (Pre_DE.visualization) {
    PreRNASeqDifferentialAnalysis(path.prefix = path.prefix,
                                  sample.pattern = sample.pattern)
  }
  if (file.exists(paste0(path.prefix,
                         "RNASeq_results/Alignment_Report/",
                         "Alignment_report_reads.csv")) &
      file.exists(paste0(path.prefix,
                         "RNASeq_results/Alignment_Report/",
                         "Overall_Mapping_rate.csv"))) {
    AlignmentPlot(path.prefix,
                  independent.variable,
                  case.group,
                  control.group,
                  phenoData.result)
  }
  message("\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618  ",
          "Start Differential Expression Analysis  ",
          "\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618\n")
  # 2. DE analysis
  if (ballgown.run) {
    BallgownAnalysis(path.prefix,
                     genome.name,
                     sample.pattern,
                     independent.variable,
                     case.group,
                     control.group,
                     ballgown.pval,
                     ballgown.log2FC,
                     phenoData.result)
  }
  raw.read.avail <- RawReadCountAvailability(path.prefix)
  if (raw.read.avail) {
    if (DESeq2.run) {
      DESeq2RawCountAnalysis(path.prefix,
                             independent.variable,
                             case.group,
                             control.group,
                             DESeq2.pval,
                             DESeq2.log2FC,
                             phenoData.result)
    }
    if (edgeR.run) {
      edgeRRawCountAnalysis(path.prefix,
                            independent.variable,
                            case.group,
                            control.group,
                            edgeR.pval,
                            edgeR.log2FC,
                            phenoData.result)
    }
  }
  if (Post_DE.visualization) {
    # 3. Post-DE assessment visualization
    PostRNASeqDifferentialAnalysis(path.prefix = path.prefix,
                                   sample.pattern = sample.pattern)
  }
}


PreRNASeqDifferentialAnalysis <- function(path.prefix, sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqDifferentialAnalysis()' ",
          "environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqDifferentialAnalysis() environment ERROR")
  }
  message("(\u2714) : RNASeqDifferentialAnalysis() pre-check is valid\n\n")
}

PostRNASeqDifferentialAnalysis <- function(path.prefix, sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqDifferentialAnalysis()' ",
          "environment post-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqDifferentialAnalysis() post-check ERROR")
  }
  message("(\u2714) : RNASeqDifferentialAnalysis() post-check is valid\n\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\n")
}
