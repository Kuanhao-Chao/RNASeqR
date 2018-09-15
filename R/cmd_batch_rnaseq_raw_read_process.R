#' @title RNASeqReadProcess_CMD
#'
#' @description
#'   Process raw reads for RNA-Seq workflow in background. \cr
#'   This function do 5 things : \cr
#'   \enumerate{
#'     \item 'Hisat2' : aligns raw reads to reference genome.
#'       If \code{indices.optional} in \code{RNASeqWorkFlowParam} is
#'       \code{FALSE}, Hisat2 indices will be created.\cr
#'     \item 'Samtools': converts '.sam' files to '.bam' files.\cr
#'     \item 'Stringtie': assembles alignments into transcript.\cr
#'     \item 'Gffcompare': examines how transcripts compare with the
#'       reference annotation. \cr
#'     \item 'Stringtie': creates input files for ballgown, edgeR and DESeq2.\cr
#'     \item raw reads count: create raw reads count for DESeq2 and edgeR \cr
#'   }\cr
#'   Before running this function, \code{RNASeqEnvironmentSet_CMD()} or
#'   \code{RNASeqEnvironmentSet()} must be executed successfully. \cr
#'   If you want to process raw reads for the following RNA-Seq workflow
#'   in R shell, please see \code{RNASeqReadProcess()} function.\cr
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related
#'   parameters
#' @param num.parallel.threads Specify the number of processing threads (CPUs)
#'   to use for transcript assembly. The default is 1.
#' @param run Default value is \code{TRUE}. If \code{TRUE},
#'   'Rscript/Environment_Set.R' will be created and executed.
#'   The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{False}, 'Rscript/Environment_Set.R' will be
#'   created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE},
#'   the result of checking \code{RNASeqWorkFlowParam}
#'   will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE},
#'   the result of checking \code{RNASeqWorkFlowParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'.
#' @param Hisat2.Index.run Whether to run 'HISAT2 index' step in this function
#'   step. Default value is \code{TRUE}. Set \code{FALSE} to skip
#'   'HISAT2 index' step.
#' @param Hisat2.Alignment.run Whether to run 'HISAT2 alignment' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2 alignment' step.
#' @param Samtools.Bam.run Whether to run 'SAMTools SAM to BAM' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'SAMTools SAM to BAM' step.
#' @param StringTie.Assemble.run Whether to run 'StringTie assembly' step in
#'   this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie assembly' step.
#' @param StringTie.Merge.Trans.run Whether to run 'StringTie GTF merging' step
#'   in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie GTF merging' step.
#' @param Gffcompare.Ref.Sample.run Whether to run 'Gffcompare comparison' step
#'   in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'Gffcompare comparison' step.
#' @param StringTie.Ballgown.run Whether to run 'StringTie ballgown creation'
#'   step in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie ballgown creation' step.
#' @param PreDECountTable.run Whether to run 'gene raw reads count creation'
#'   step in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'gene raw reads count creation' step.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' ## Before run this function, make sure \code{RNASeqEnvironmentSet_CMD()}
#' ## (or\code{RNASeqEnvironmentSet()}) is executed successfully.
#' RNASeqReadProcess_CMD(RNASeqWorkFlowParam = exp,
#'                       num.parallel.threads = 10)}
RNASeqReadProcess_CMD <- function(RNASeqWorkFlowParam,
                                  num.parallel.threads      = 1,
                                  run                       = TRUE,
                                  check.s4.print            = TRUE,
                                  Hisat2.Index.run          = TRUE,
                                  Hisat2.Alignment.run      = TRUE,
                                  Samtools.Bam.run          = TRUE,
                                  StringTie.Assemble.run    = TRUE,
                                  StringTie.Merge.Trans.run = TRUE,
                                  Gffcompare.Ref.Sample.run = TRUE,
                                  StringTie.Ballgown.run    = TRUE,
                                  PreDECountTable.run       = TRUE) {
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- "@"(RNASeqWorkFlowParam, path.prefix)
  input.path.prefix <- "@"(RNASeqWorkFlowParam, input.path.prefix)
  genome.name <- "@"(RNASeqWorkFlowParam, genome.name)
  sample.pattern <- "@"(RNASeqWorkFlowParam, sample.pattern)
  python.variable <- "@"(RNASeqWorkFlowParam, python.variable)
  python.variable.answer <- python.variable$check.answer
  python.variable.version <- python.variable$python.version
  python.2to3 <- "@"(RNASeqWorkFlowParam, python.2to3)
  indices.optional <- "@"(RNASeqWorkFlowParam, indices.optional)
  fileConn<-file(paste0(path.prefix, "Rscript/Read_Process.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0('RNASeqReadProcess(path.prefix = "', path.prefix,
                   '", input.path.prefix = "', input.path.prefix,
                   '", genome.name = "', genome.name,
                   '", sample.pattern = "', sample.pattern,
                   '", python.variable.answer = ', python.variable.answer,
                   ', python.variable.version = ', python.variable.version,
                   ', python.2to3 = ', python.2to3,
                   ', num.parallel.threads = ', num.parallel.threads,
                   ', indices.optional = ', indices.optional,
                   ', Hisat2.Index.run = ', Hisat2.Index.run,
                   ', Hisat2.Alignment.run = ', Hisat2.Alignment.run,
                   ', Samtools.Bam.run = ', Samtools.Bam.run,
                   ', StringTie.Assemble.run = ', StringTie.Assemble.run,
                   ', StringTie.Merge.Trans.run = ', StringTie.Merge.Trans.run,
                   ', Gffcompare.Ref.Sample.run = ', Gffcompare.Ref.Sample.run,
                   ', StringTie.Ballgown.run = ', StringTie.Ballgown.run,
                   ', PreDECountTable.run = ', PreDECountTable.run, ')')
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message(paste0("\u2605 '", path.prefix,
                 "Rscript/Read_Process.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup',
            args = paste0("R CMD BATCH ", path.prefix,
                          "Rscript/Read_Process.R ", path.prefix,
                          "Rscript_out/Read_Process.Rout"),
            stdout = "",
            wait = FALSE)
    message(paste0("\u2605 RNASeq alignment, assembly, quantification, ",
                   "mergence, comparison, reads process are doing in the ",
                   "background. Check current progress in '", path.prefix,
                   "Rscript_out/Read_Process.Rout'\n\n"))
  }
}


#' @title RNASeqReadProcess
#'
#' @description
#'   Process raw reads for RNA-Seq workflow in R shell \cr
#'   This function do 5 things : \cr
#'   \enumerate{
#'     \item 'Hisat2' : aligns raw reads to reference genome.
#'       If \code{indices.optional} in \code{RNASeqWorkFlowParam} is
#'       \code{FALSE}, Hisat2 indices will be created.\cr
#'     \item 'Samtools': converts '.sam' files to '.bam' files.\cr
#'     \item 'Stringtie': assembles alignments into transcript.\cr
#'     \item 'Gffcompare': examines how transcripts compare with the
#'       reference annotation. \cr
#'     \item 'Stringtie': creates input files for ballgown, edgeR and DESeq2.\cr
#'     \item raw reads count: create raw reads count for DESeq2 and edgeR \cr
#'   }\cr
#'   Before running this function, \code{RNASeqEnvironmentSet_CMD()} or
#'   \code{RNASeqEnvironmentSet()} must be executed successfully.
#'   If you want to process raw reads for the following RNA-Seq workflow in
#'   background, please see \code{RNASeqReadProcess_CMD()} function.
#'
#' @param path.prefix path prefix of 'gene_data/', 'RNASeq_bin/',
#'   'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @param input.path.prefix path prefix of 'input_files/' directory
#' @param genome.name Variable of genome name defined in this RNA-Seq workflow
#'   (ex. \code{genome.name}.fa, \code{genome.name}.gtf)
#' @param sample.pattern  Regular expression of paired-end fastq.gz files under
#'   'input_files/raw_fastq.gz'. Expression not includes \code{_[1,2].fastq.gz}.
#' @param python.variable.answer logical value whether python is available
#'   on the device
#' @param python.variable.version python version of the device
#' @param python.2to3 logical value whether \code{2to3} is available
#'   on the device
#' @param num.parallel.threads Specify the number of processing threads (CPUs)
#'   to use for transcript assembly. The default is 1.
#' @param indices.optional logical value whether 'indices/' is
#'   exit in 'input_files/'
#' @param Hisat2.Index.run Whether to run 'HISAT2 index' step in this function
#'   step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2 index' step.
#' @param Hisat2.Alignment.run Whether to run 'HISAT2 alignment' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2 alignment' step.
#' @param Samtools.Bam.run Whether to run 'SAMTools SAM to BAM' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'SAMTools SAM to BAM' step.
#' @param StringTie.Assemble.run Whether to run 'StringTie assembly' step in
#'   this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie assembly' step.
#' @param StringTie.Merge.Trans.run Whether to run 'StringTie GTF merging' step
#'   in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie GTF merging' step.
#' @param Gffcompare.Ref.Sample.run Whether to run 'Gffcompare comparison' step
#'   in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'Gffcompare comparison' step.
#' @param StringTie.Ballgown.run Whether to run 'StringTie ballgown creation'
#'   step in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie ballgown creation' step.
#' @param PreDECountTable.run Whether to run 'gene raw reads count creation'
#'   step in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'gene raw reads count creation' step.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' ## Before run this function, make sure \code{RNASeqEnvironmentSet_CMD()}
#' ##(or\code{RNASeqEnvironmentSet()}) is executed successfully.
#' RNASeqReadProcess(path.prefix             = yeast@@path.prefix,
#'                   input.path.prefix       = yeast@@input.path.prefix,
#'                   genome.name             = yeast@@genome.name,
#'                   sample.pattern          = yeast@@sample.pattern,
#'                   python.variable.answer  = yeast@@python.variable[0],
#'                   python.variable.version = yeast@@python.variable[1],
#'                   indices.optional        = yeast@@indices.optional,
#'                   num.parallel.threads    = 10)}
RNASeqReadProcess <- function(path.prefix,
                              input.path.prefix,
                              genome.name,
                              sample.pattern,
                              python.variable.answer,
                              python.variable.version,
                              python.2to3,
                              num.parallel.threads      = 1,
                              indices.optional,
                              Hisat2.Index.run          = TRUE,
                              Hisat2.Alignment.run      = TRUE,
                              Samtools.Bam.run          = TRUE,
                              StringTie.Assemble.run    = TRUE,
                              StringTie.Merge.Trans.run = TRUE,
                              Gffcompare.Ref.Sample.run = TRUE,
                              StringTie.Ballgown.run    = TRUE,
                              PreDECountTable.run       = TRUE) {
  CheckOperatingSystem(FALSE)
  ExportPath(path.prefix)
  PreRNASeqReadProcess(path.prefix, genome.name, sample.pattern)
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)
  if (check.results$ht2.files.number.df == 0 &&
      !indices.optional & Hisat2.Index.run) {
    CreateHisat2Index(path.prefix, genome.name, sample.pattern)
  }
  if (Hisat2.Alignment.run) {
    Hisat2AlignmentDefault(path.prefix,
                           genome.name,
                           sample.pattern,
                           num.parallel.threads)
  }
  if (Samtools.Bam.run) {
    SamtoolsToBam(path.prefix,
                  genome.name,
                  sample.pattern,
                  num.parallel.threads)
  }
  if (StringTie.Assemble.run) {
    StringTieAssemble(path.prefix,
                      genome.name,
                      sample.pattern,
                      num.parallel.threads)
  }
  if (StringTie.Merge.Trans.run) {
    StringTieMergeTrans(path.prefix,
                        genome.name,
                        sample.pattern,
                        num.parallel.threads)
  }
  if (Gffcompare.Ref.Sample.run) {
    GffcompareRefSample(path.prefix,
                        genome.name,
                        sample.pattern)
  }
  if (StringTie.Ballgown.run) {
    StringTieToBallgown(path.prefix,
                        genome.name,
                        sample.pattern,
                        num.parallel.threads)
  }
  finals <- ProgressGenesFiles(path.prefix,
                               genome.name,
                               sample.pattern,
                               print=TRUE)
  if (PreDECountTable.run) {
    PreDECountTable(path.prefix,
                    sample.pattern,
                    python.variable.answer,
                    python.variable.version,
                    python.2to3,
                    print=TRUE)
  }
  PostRNASeqReadProcess(path.prefix,
                        genome.name,
                        sample.pattern)
}

PreRNASeqReadProcess <- function(path.prefix, genome.name, sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqReadProcess()' ",
          "environment pre-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  ref.gtf <- file.exists(paste0(path.prefix,
                                "gene_data/ref_genes/", genome.name, ".gtf"))
  ref.fa <- file.exists(paste0(path.prefix,
                               "gene_data/ref_genome/", genome.name, ".fa"))
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'),
                          pattern = sample.pattern,
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
  check.tool.result <- CheckToolAll(path.prefix)
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)
  check.progress.results.bool <- check.results$gtf.file.logic.df &&
    check.results$fa.file.logic.df &&
    (check.results$fastq.gz.files.number.df != 0)
  validity <- phenodata.csv && ref.gtf && ref.fa && check.tool.result &&
    (length(raw.fastq) != 0) && check.progress.results.bool
  if (!isTRUE(validity)) {
    stop("RNASeqReadProcess() environment ERROR")
  }
  message("(\u2714) : RNASeqReadProcess() pre-check is valid\n\n")
}

PostRNASeqReadProcess <- function(path.prefix, genome.name, sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqReadProcess()' ",
          "environment post-check ...\n")
  # Still need to add condition
  gene_abundance <- dir.exists(paste0(path.prefix, "gene_data/gene_abundance/"))
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)
  ht2.bool <- (check.results$ht2.files.number.df) != 0
  sam.bool <- (check.results$sam.files.number.df) != 0
  bam.bool <- (check.results$bam.files.number.df) != 0
  gtf.bool <- (check.results$gtf.files.number.df) != 0
  merged.bool <- check.results$stringtie_merged.gtf.file.df
  gffcompare.bool <- (check.results$gffcompare.related.dirs.number.df) != 0
  ballgown.bool <- (check.results$ballgown.dirs.number.df) != 0
  validity <- gene_abundance && ht2.bool && sam.bool && bam.bool && gtf.bool &&
    merged.bool && gffcompare.bool && ballgown.bool
  if (!isTRUE(validity)) {
    stop("RNASeqReadProcess() post-check ERROR")
  }
  message("(\u2714) : RNASeqReadProcess() post-check is valid\n\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
}

