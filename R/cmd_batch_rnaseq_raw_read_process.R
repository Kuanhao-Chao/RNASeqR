#' @title RNASeqReadProcess_CMD
#'
#' @description
#'   Process raw reads for RNA-Seq workflow in background. \cr
#'   This function do 5 things : \cr
#'   \enumerate{
#'     \item 'Hisat2' : aligns raw reads to reference genome.
#'       If \code{indices.optional} in \code{RNASeqRParam} is
#'       \code{FALSE}, Hisat2 indices will be created.\cr
#'     \item 'Rsamtools': converts '.sam' files to '.bam' files.\cr
#'     \item 'Stringtie': assembles alignments into transcript.\cr
#'     \item 'Gffcompare': examines how transcripts compare with the
#'       reference annotation. \cr
#'     \item 'Stringtie': creates input files for ballgown, edgeR and DESeq2.\cr
#'     \item raw reads count: create raw reads count for DESeq2 and edgeR \cr
#'   }
#'   Before running this function, \code{RNASeqEnvironmentSet_CMD()} or
#'   \code{RNASeqEnvironmentSet()} must be executed successfully. \cr
#'   If you want to process raw reads for the following RNA-Seq workflow
#'   in R shell, please see \code{RNASeqReadProcess()} function.\cr
#'
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param SAMtools.or.Rsamtools Default value is \code{Rsamtools}. User can set
#'   to \code{SAMtools} to use command-line-based 'samtools' instead.
#' @param num.parallel.threads Specify the number of processing threads (CPUs)
#'   to use for each step. The default is 1.
#' @param Rsamtools.nCores The number of cores to use when running
#'   'Rsamtools' step.
#' @param Hisat2.Index.run Whether to run 'HISAT2 index' step in this function
#'   step. Default value is \code{TRUE}. Set \code{FALSE} to skip
#'   'HISAT2 index' step.
#' @param Hisat2.Alignment.run Whether to run 'HISAT2 alignment' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2 alignment' step.
#'
#'
#'
#'
#' @param STAR.Alignment.run Whether to run 'STAR alignment' step in this
#'   function step. Default value is \code{FALSE}.
#'   Set \code{TRUE} to run 'STAR alignment' step.
#'
#'
#'
#'
#'
#' @param Rsamtools.Bam.run Whether to run 'Rsamtools SAM to BAM' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'Rsamtools SAM to BAM' step.
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
#' @param run Default value is \code{TRUE}. If \code{TRUE},
#'   'Rscript/Environment_Set.R' will be created and executed.
#'   The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{False}, 'Rscript/Environment_Set.R' will be
#'   created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE},
#'   the result of checking \code{RNASeqRParam}
#'   will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE},
#'   the result of checking \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' ## Before run this function, make sure \code{RNASeqEnvironmentSet_CMD()}
#' ## (or\code{RNASeqEnvironmentSet()}) is executed successfully.
#' RNASeqReadProcess_CMD(RNASeqRParam = yeast,
#'                       num.parallel.threads = 10)}
RNASeqReadProcess_CMD <- function(RNASeqRParam,
                                  SAMtools.or.Rsamtools     = "Rsamtools",
                                  num.parallel.threads      = 1,
                                  Hisat2.Index.run          = TRUE,
                                  Hisat2.Index.num.parallel.threads = 1,
                                  Hisat2.Index.large.index = FALSE,
                                  Hisat2.Index.local.ftab.chars = 6,
                                  Hisat2.Index.local.off.rate = 3,
                                  Hisat2.Index.ftab.chars = 10,
                                  Hisat2.Index.off.rate = 4,
                                  Hisat2.Alignment.run      = TRUE,
                                  Hisat2.Alignment.num.parallel.threads = 1,
                                  Hisat2.Alignment.skip = 0,
                                  Hisat2.Alignment.qupto = "None",
                                  Hisat2.Alignment.trim5 = 0,
                                  Hisat2.Alignment.trim3 = 0,
                                  Hisat2.Alignment.phred = 33,
                                  Hisat2.Alignment.int.quals = FALSE,
                                  Hisat2.Alignment.n.ceil.1.function.type = "L",
                                  Hisat2.Alignment.n.ceil.2.constant.term = 0,
                                  Hisat2.Alignment.n.ceil.3.coefficient = 0.15,
                                  Hisat2.Alignment.mp.MX = 6,
                                  Hisat2.Alignment.mp.MN = 2,
                                  Hisat2.Alignment.sp.MX = 2,
                                  Hisat2.Alignment.sp.MN = 1,
                                  Hisat2.Alignment.np = 1,
                                  Hisat2.Alignment.rdg.1 = 5,
                                  Hisat2.Alignment.rdg.2 = 3,
                                  Hisat2.Alignment.rfg.1 = 5,
                                  Hisat2.Alignment.rfg.2 = 3,
                                  Hisat2.Alignment.score.min.1.function.type = "L",
                                  Hisat2.Alignment.score.min.2.constant.term = 0,
                                  Hisat2.Alignment.score.min.3.coefficient = -0.2,
                                  Hisat2.Alignment.pen.cansplice = 0,
                                  Hisat2.Alignment.penc.noncansplice = 12,
                                  Hisat2.Alignment.pen.canintronlen.1.function.type = "G",
                                  Hisat2.Alignment.pen.canintronlen.2.constant.term = -8,
                                  Hisat2.Alignment.pen.canintronlen.3.coefficient = 1,
                                  Hisat2.Alignment.pen.noncanintronlen.1.function.type = "G",
                                  Hisat2.Alignment.pen.noncanintronlen.2.constant.term = -8,
                                  Hisat2.Alignment.pen.noncanintronlen.3.coefficient = 1,
                                  Hisat2.Alignment.min.intronlen = 20,
                                  Hisat2.Alignment.max.intronlen = 500000,
                                  Hisat2.Alignment.rna.strandness = "None",
                                  Hisat2.Alignment.k = 5,
                                  Hisat2.Alignment.max.seeds = 5,
                                  Hisat2.Alignment.secondary = FALSE,
                                  Hisat2.Alignment.minins = 0,
                                  Hisat2.Alignment.maxins = 500,
                                  Hisat2.Alignment.seed = 0,
                                  STAR.Index.num.parallel.threads = 1,
                                  STAR.Index.sjdbOverhang.Read.length = 100,
                                  STAR.Index.genomeSAindexNbases = 14,
                                  STAR.Index.genomeChrBinNbits = 18,
                                  STAR.Index.genomeSAsparseD = 1,
                                  STAR.Alignment.run        = FALSE,
                                  STAR.Alignment.num.parallel.threads = 1,
                                  STAR.Alignment.genomeLoad = "NoSharedMemory",
                                  STAR.Alignment.readMapNumber = -1,
                                  STAR.Alignment.clip3pNbases = 0,
                                  STAR.Alignment.clip5pNbases = 0,
                                  STAR.Alignment.clip3pAdapterSeq = "-",
                                  STAR.Alignment.clip3pAdapterMMp = 0.1,
                                  STAR.Alignment.clip3pAfterAdapterNbases = 0,
                                  STAR.Alignment.limitGenomeGenerateRAM = 31000000000,
                                  STAR.Alignment.limitIObufferSize = 150000000,
                                  STAR.Alignment.limitOutSAMoneReadBytes = 100000,
                                  STAR.Alignment.limitOutSJoneRead = 1000,
                                  STAR.Alignment.limitOutSJcollapsed = 1000000,
                                  STAR.Alignment.limitBAMsortRAM = 0,
                                  STAR.Alignment.outReadsUnmapped = "None",
                                  STAR.Alignment.outQSconversionAdd = 0,
                                  STAR.Alignment.outSAMprimaryFlag = "OneBestScore",
                                  STAR.Alignment.outSAMmapqUnique = 255,
                                  STAR.Alignment.scoreGap = 0,
                                  STAR.Alignment.scoreGapNoncan = -8,
                                  STAR.Alignment.scoreGapGCAG = -4,
                                  STAR.Alignment.scoreGapATAC = -8,
                                  STAR.Alignment.scoreGenomicLengthLog2scale = -0.25,
                                  STAR.Alignment.scoreDelOpen = -2,
                                  STAR.Alignment.scoreDelBase = -2,
                                  STAR.Alignment.scoreInsOpen = -2,
                                  STAR.Alignment.scoreInsBase = -2,
                                  STAR.Alignment.scoreStitchSJshift = 1,
                                  STAR.Alignment.seedSearchStartLmax = 50,
                                  STAR.Alignment.seedSearchStartLmaxOverLread = 1.0,
                                  STAR.Alignment.seedSearchLmax = 0,
                                  STAR.Alignment.seedMultimapNmax = 10000,
                                  STAR.Alignment.seedPerReadNmax = 1000,
                                  STAR.Alignment.seedPerWindowNmax = 50,
                                  STAR.Alignment.seedNoneLociPerWindow = 10,
                                  STAR.Alignment.alignIntronMin = 21,
                                  STAR.Alignment.alignIntronMax = 0,
                                  STAR.Alignment.alignMatesGapMax = 0,
                                  STAR.Alignment.alignSJoverhangMin = 5,
                                  STAR.Alignment.alignSJDBoverhangMin = 3,
                                  STAR.Alignment.alignSplicedMateMapLmin = 0,
                                  STAR.Alignment.alignSplicedMateMapLminOverLmate = 0.66,
                                  STAR.Alignment.alignWindowsPerReadNmax = 10000,
                                  STAR.Alignment.alignTranscriptsPerWindowNmax = 100,
                                  STAR.Alignment.alignTranscriptsPerReadNmax = 10000,
                                  STAR.Alignment.alignEndsType = "Local",
                                  STAR.Alignment.winAnchorMultimapNmax = 50,
                                  STAR.Alignment.winBinNbits = 16,
                                  STAR.Alignment.winAnchorDistNbins = 9,
                                  STAR.Alignment.winFlankNbins = 4,
                                  Rsamtools.Bam.run         = TRUE,
                                  Samtools.Bam.num.parallel.threads = 1,
                                  Rsamtools.nCores          = 1,
                                  StringTie.Assemble.run    = TRUE,
                                  Stringtie.Assembly.num.parallel.threads = 1,
                                  Stringtie.Assembly.f = 0.1,
                                  Stringtie.Assembly.m = 200,
                                  Stringtie.Assembly.c = 2.5,
                                  Stringtie.Assembly.g = 50,
                                  Stringtie.Assembly.M = 0.95,
                                  StringTie.Merge.Trans.run = TRUE,
                                  Stringtie.Merge.num.parallel.threads = 1,
                                  Gffcompare.Ref.Sample.run = TRUE,
                                  StringTie.Ballgown.run    = TRUE,
                                  Stringtie.2.Ballgown.num.parallel.threads = 1,
                                  PreDECountTable.run       = TRUE,
                                  run                       = TRUE,
                                  check.s4.print            = TRUE) {
  CheckS4Object(RNASeqRParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  INSIDE.path.prefix <- "@"(RNASeqRParam, path.prefix)
  saveRDS(RNASeqRParam,
          file = paste0(INSIDE.path.prefix,
                        "gene_data/RNASeqRParam.rds"))
  fileConn<-file(paste0(path.prefix, "Rscript/Read_Process.R"))
  first <- "library(RNASeqR)"
  second <- paste0("RNASeqReadProcess(RNASeqRParam = 'INSIDE'",
                   ", which.trigger = 'INSIDE'",
                   ", INSIDE.path.prefix = '", INSIDE.path.prefix,
                   "', SAMtools.or.Rsamtools = '", SAMtools.or.Rsamtools,
                   "', num.parallel.threads = ", num.parallel.threads,
                   ", Rsamtools.nCores = ", Rsamtools.nCores,
                   ", Hisat2.Index.run = ", Hisat2.Index.run,
                   ", Hisat2.Alignment.run = ", Hisat2.Alignment.run,
                   ", Rsamtools.Bam.run = ", Rsamtools.Bam.run,
                   ", StringTie.Assemble.run = ", StringTie.Assemble.run,
                   ", StringTie.Merge.Trans.run = ", StringTie.Merge.Trans.run,
                   ", Gffcompare.Ref.Sample.run = ", Gffcompare.Ref.Sample.run,
                   ", StringTie.Ballgown.run = ", StringTie.Ballgown.run,
                   ", PreDECountTable.run = ", PreDECountTable.run, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message("\u2605 '", path.prefix,
          "Rscript/Read_Process.R' has been created.\n")
  if (run) {
    R.home.lib <- R.home()
    R.home.bin <- gsub("/lib/R", "/bin/R", R.home.lib)
    system2(command = "nohup",
            args = paste0(R.home.bin, " CMD BATCH ",
                          path.prefix,
                          "Rscript/Read_Process.R ", path.prefix,
                          "Rscript_out/Read_Process.Rout"),
            stdout = "",
            wait = FALSE)
    message("\u2605 RNASeq alignment, assembly, quantification, ",
            "mergence, comparison, reads process are doing in the ",
            "background. Check current progress in '", path.prefix,
            "Rscript_out/Read_Process.Rout'\n\n")
  }
}


#' @title RNASeqReadProcess
#'
#' @description
#'   Process raw reads for RNA-Seq workflow in R shell \cr
#'   This function do 5 things : \cr
#'   \enumerate{
#'     \item 'Hisat2' : aligns raw reads to reference genome.
#'       If \code{indices.optional} in \code{RNASeqRParam} is
#'       \code{FALSE}, Hisat2 indices will be created.\cr
#'     \item 'Rsamtools': converts '.sam' files to '.bam' files.\cr
#'     \item 'Stringtie': assembles alignments into transcript.\cr
#'     \item 'Gffcompare': examines how transcripts compare with the
#'       reference annotation. \cr
#'     \item 'Stringtie': creates input files for ballgown, edgeR and DESeq2.\cr
#'     \item raw reads count: create raw reads count for DESeq2 and edgeR \cr
#'   }
#'   Before running this function, \code{RNASeqEnvironmentSet_CMD()} or
#'   \code{RNASeqEnvironmentSet()} must be executed successfully.
#'   If you want to process raw reads for the following RNA-Seq workflow in
#'   background, please see \code{RNASeqReadProcess_CMD()} function.
#'
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param which.trigger Default value is \code{OUTSIDE}. User should not change
#'   this value.
#' @param INSIDE.path.prefix Default value is \code{NA}. User should not change
#'   this value.
#' @param SAMtools.or.Rsamtools Default value is \code{Rsamtools}. User can set
#'   to \code{SAMtools} to use command-line-based 'samtools' instead.
#' @param num.parallel.threads Specify the number of processing threads (CPUs)
#'   to use for each step. The default is 1.
#' @param Rsamtools.nCores The number of cores to use when running
#'   'Rsamtools' step.
#' @param Hisat2.Index.run Whether to run 'HISAT2 index' step in this function
#'   step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2 index' step.
#' @param Hisat2.Alignment.run Whether to run 'HISAT2 alignment' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2 alignment' step.
#'
#'
#'
#' @param STAR.Alignment.run Whether to run 'STAR alignment' step in this
#'   function step. Default value is \code{FALSE}.
#'   Set \code{TRUE} to run 'STAR alignment' step.
#'
#'
#'
#'
#' @param Rsamtools.Bam.run Whether to run 'Rsamtools SAM to BAM' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'Rsamtools SAM to BAM' step.
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
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE},
#'   the result of checking \code{RNASeqRParam}
#'   will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE},
#'   the result of checking \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' ## Before run this function, make sure \code{RNASeqEnvironmentSet_CMD()}
#' ##(or\code{RNASeqEnvironmentSet()}) is executed successfully.
#' RNASeqReadProcess(RNASeqRParam         = yeast,
#'                   num.parallel.threads = 10)}
RNASeqReadProcess <- function(RNASeqRParam,
                              which.trigger             = "OUTSIDE",
                              INSIDE.path.prefix        = NA,
                              SAMtools.or.Rsamtools     = "Rsamtools",
                              num.parallel.threads      = 1,
                              Hisat2.Index.run          = TRUE,
                              Hisat2.Index.num.parallel.threads = 1,
                              Hisat2.Index.large.index = FALSE,
                              Hisat2.Index.local.ftab.chars = 6,
                              Hisat2.Index.local.off.rate = 3,
                              Hisat2.Index.ftab.chars = 10,
                              Hisat2.Index.off.rate = 4,
                              Hisat2.Alignment.run      = TRUE,
                              Hisat2.Alignment.num.parallel.threads = 1,
                              Hisat2.Alignment.skip = 0,
                              Hisat2.Alignment.qupto = "None",
                              Hisat2.Alignment.trim5 = 0,
                              Hisat2.Alignment.trim3 = 0,
                              Hisat2.Alignment.phred = 33,
                              Hisat2.Alignment.int.quals = FALSE,
                              Hisat2.Alignment.n.ceil.1.function.type = "L",
                              Hisat2.Alignment.n.ceil.2.constant.term = 0,
                              Hisat2.Alignment.n.ceil.3.coefficient = 0.15,
                              Hisat2.Alignment.mp.MX = 6,
                              Hisat2.Alignment.mp.MN = 2,
                              Hisat2.Alignment.sp.MX = 2,
                              Hisat2.Alignment.sp.MN = 1,
                              Hisat2.Alignment.np = 1,
                              Hisat2.Alignment.rdg.1 = 5,
                              Hisat2.Alignment.rdg.2 = 3,
                              Hisat2.Alignment.rfg.1 = 5,
                              Hisat2.Alignment.rfg.2 = 3,
                              Hisat2.Alignment.score.min.1.function.type = "L",
                              Hisat2.Alignment.score.min.2.constant.term = 0,
                              Hisat2.Alignment.score.min.3.coefficient = -0.2,
                              Hisat2.Alignment.pen.cansplice = 0,
                              Hisat2.Alignment.penc.noncansplice = 12,
                              Hisat2.Alignment.pen.canintronlen.1.function.type = "G",
                              Hisat2.Alignment.pen.canintronlen.2.constant.term = -8,
                              Hisat2.Alignment.pen.canintronlen.3.coefficient = 1,
                              Hisat2.Alignment.pen.noncanintronlen.1.function.type = "G",
                              Hisat2.Alignment.pen.noncanintronlen.2.constant.term = -8,
                              Hisat2.Alignment.pen.noncanintronlen.3.coefficient = 1,
                              Hisat2.Alignment.min.intronlen = 20,
                              Hisat2.Alignment.max.intronlen = 500000,
                              Hisat2.Alignment.rna.strandness = "None",
                              Hisat2.Alignment.k = 5,
                              Hisat2.Alignment.max.seeds = 5,
                              Hisat2.Alignment.secondary = FALSE,
                              Hisat2.Alignment.minins = 0,
                              Hisat2.Alignment.maxins = 500,
                              Hisat2.Alignment.seed = 0,
                              STAR.Index.num.parallel.threads = 1,
                              STAR.Index.sjdbOverhang.Read.length = 100,
                              STAR.Index.genomeSAindexNbases = 14,
                              STAR.Index.genomeChrBinNbits = 18,
                              STAR.Index.genomeSAsparseD = 1,
                              STAR.Alignment.run        = FALSE,
                              STAR.Alignment.num.parallel.threads = 1,
                              STAR.Alignment.genomeLoad = "NoSharedMemory",
                              STAR.Alignment.readMapNumber = -1,
                              STAR.Alignment.clip3pNbases = 0,
                              STAR.Alignment.clip5pNbases = 0,
                              STAR.Alignment.clip3pAdapterSeq = "-",
                              STAR.Alignment.clip3pAdapterMMp = 0.1,
                              STAR.Alignment.clip3pAfterAdapterNbases = 0,
                              STAR.Alignment.limitGenomeGenerateRAM = 31000000000,
                              STAR.Alignment.limitIObufferSize = 150000000,
                              STAR.Alignment.limitOutSAMoneReadBytes = 100000,
                              STAR.Alignment.limitOutSJoneRead = 1000,
                              STAR.Alignment.limitOutSJcollapsed = 1000000,
                              STAR.Alignment.limitBAMsortRAM = 0,
                              STAR.Alignment.outReadsUnmapped = "None",
                              STAR.Alignment.outQSconversionAdd = 0,
                              STAR.Alignment.outSAMprimaryFlag = "OneBestScore",
                              STAR.Alignment.outSAMmapqUnique = 255,
                              STAR.Alignment.scoreGap = 0,
                              STAR.Alignment.scoreGapNoncan = -8,
                              STAR.Alignment.scoreGapGCAG = -4,
                              STAR.Alignment.scoreGapATAC = -8,
                              STAR.Alignment.scoreGenomicLengthLog2scale = -0.25,
                              STAR.Alignment.scoreDelOpen = -2,
                              STAR.Alignment.scoreDelBase = -2,
                              STAR.Alignment.scoreInsOpen = -2,
                              STAR.Alignment.scoreInsBase = -2,
                              STAR.Alignment.scoreStitchSJshift = 1,
                              STAR.Alignment.seedSearchStartLmax = 50,
                              STAR.Alignment.seedSearchStartLmaxOverLread = 1.0,
                              STAR.Alignment.seedSearchLmax = 0,
                              STAR.Alignment.seedMultimapNmax = 10000,
                              STAR.Alignment.seedPerReadNmax = 1000,
                              STAR.Alignment.seedPerWindowNmax = 50,
                              STAR.Alignment.seedNoneLociPerWindow = 10,
                              STAR.Alignment.alignIntronMin = 21,
                              STAR.Alignment.alignIntronMax = 0,
                              STAR.Alignment.alignMatesGapMax = 0,
                              STAR.Alignment.alignSJoverhangMin = 5,
                              STAR.Alignment.alignSJDBoverhangMin = 3,
                              STAR.Alignment.alignSplicedMateMapLmin = 0,
                              STAR.Alignment.alignSplicedMateMapLminOverLmate = 0.66,
                              STAR.Alignment.alignWindowsPerReadNmax = 10000,
                              STAR.Alignment.alignTranscriptsPerWindowNmax = 100,
                              STAR.Alignment.alignTranscriptsPerReadNmax = 10000,
                              STAR.Alignment.alignEndsType = "Local",
                              STAR.Alignment.winAnchorMultimapNmax = 50,
                              STAR.Alignment.winBinNbits = 16,
                              STAR.Alignment.winAnchorDistNbins = 9,
                              STAR.Alignment.winFlankNbins = 4,
                              Rsamtools.Bam.run         = TRUE,
                              Samtools.Bam.num.parallel.threads = 1,
                              Rsamtools.nCores          = 1,
                              StringTie.Assemble.run    = TRUE,
                              Stringtie.Assembly.num.parallel.threads = 1,
                              Stringtie.Assembly.f = 0.1,
                              Stringtie.Assembly.m = 200,
                              Stringtie.Assembly.c = 2.5,
                              Stringtie.Assembly.g = 50,
                              Stringtie.Assembly.M = 0.95,
                              StringTie.Merge.Trans.run = TRUE,
                              Stringtie.Merge.num.parallel.threads = 1,
                              Gffcompare.Ref.Sample.run = TRUE,
                              StringTie.Ballgown.run    = TRUE,
                              Stringtie.2.Ballgown.num.parallel.threads = 1,
                              PreDECountTable.run       = TRUE,
                              check.s4.print            = TRUE) {
  CheckOperatingSystem(FALSE)
  if (SAMtools.or.Rsamtools != "Rsamtools" &
      SAMtools.or.Rsamtools != "SAMtools") {
    stop("'SAMtools.or.Rsamtools' must be 'Rsamtools' or 'SAMtools'!!")
  }
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
  # To find 'HISAT2', 'StringTie' and 'Gffcompare'
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  input.path.prefix <- "@"(RNASeqRParam, input.path.prefix)
  genome.name <- "@"(RNASeqRParam, genome.name)
  sample.pattern <- "@"(RNASeqRParam, sample.pattern)
  independent.variable <- "@"(RNASeqRParam, independent.variable)
  case.group <- "@"(RNASeqRParam, case.group)
  control.group <- "@"(RNASeqRParam, control.group)
  python.variable <- "@"(RNASeqRParam, python.variable)
  python.variable.answer <- python.variable$check.answer
  python.variable.version <- python.variable$python.version
  python.2to3 <- "@"(RNASeqRParam, python.2to3)
  indices.optional <- "@"(RNASeqRParam, indices.optional)
  ExportPath(path.prefix)
  PreRNASeqReadProcess(path.prefix, genome.name, sample.pattern)
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)


  # Check alignemnt selection first !
  if (isTRUE(Hisat2.Alignment.run) && isTRUE(STAR.Alignment.run)) {
    stop("Hisat2 and STAR can not run at the same time !")
  } else if (isTRUE(Hisat2.Alignment.run) && !isTRUE(STAR.Alignment.run)) {
    message("Hisat2 is selected as aligner in the pipeline !")
  } else if (!isTRUE(Hisat2.Alignment.run) && isTRUE(STAR.Alignment.run)) {
    message("STAR is selected as aligner in the pipeline !")
  }


  if (check.results$ht2.files.number.df == 0 &&
      !indices.optional & Hisat2.Index.run) {
    CreateHisat2Index(path.prefix,
                      genome.name,
                      sample.pattern,
                      splice.site.info = TRUE,
                      exon.info = TRUE,
                      Hisat2.Index.num.parallel.threads,
                      Hisat2.large.index,
                      Hisat2.local.ftab.chars,
                      Hisat2.local.off.rate,
                      Hisat2.ftab.chars,
                      Hisat2.off.rate)
  }

  if (Hisat2.Alignment.run) {
    Hisat2AlignmentDefault(path.prefix,
                           genome.name,
                           sample.pattern,
                           independent.variable,
                           case.group,
                           control.group,
                           Hisat2.Alignment.num.parallel.threads,
                           Hisat2.Alignment.skip,
                           Hisat2.Alignment.qupto,
                           Hisat2.Alignment.trim5,
                           Hisat2.Alignment.trim3,
                           Hisat2.Alignment.phred,
                           Hisat2.Alignment.int.quals,
                           Hisat2.Alignment.n.ceil.1.function.type,
                           Hisat2.Alignment.n.ceil.2.constant.term,
                           Hisat2.Alignment.n.ceil.3.coefficient,
                           Hisat2.Alignment.mp.MX,
                           Hisat2.Alignment.mp.MN,
                           Hisat2.Alignment.sp.MX,
                           Hisat2.Alignment.sp.MN,
                           Hisat2.Alignment.np,
                           Hisat2.Alignment.rdg.1,
                           Hisat2.Alignment.rdg.2,
                           Hisat2.Alignment.rfg.1,
                           Hisat2.Alignment.rfg.2,
                           Hisat2.Alignment.score.min.1.function.type,
                           Hisat2.Alignment.score.min.2.constant.term,
                           Hisat2.Alignment.score.min.3.coefficient,
                           Hisat2.Alignment.pen.cansplice,
                           Hisat2.Alignment.penc.noncansplice,
                           Hisat2.Alignment.pen.canintronlen.1.function.type,
                           Hisat2.Alignment.pen.canintronlen.2.constant.term,
                           Hisat2.Alignment.pen.canintronlen.3.coefficient,
                           Hisat2.Alignment.pen.noncanintronlen.1.function.type,
                           Hisat2.Alignment.pen.noncanintronlen.2.constant.term,
                           Hisat2.Alignment.pen.noncanintronlen.3.coefficient,
                           Hisat2.Alignment.min.intronlen,
                           Hisat2.Alignment.max.intronlen,
                           Hisat2.Alignment.rna.strandness,
                           Hisat2.Alignment.k,
                           Hisat2.Alignment.max.seeds,
                           Hisat2.Alignment.secondary,
                           Hisat2.Alignment.minins,
                           Hisat2.Alignment.maxins,
                           Hisat2.Alignment.seed)
    }

  if (STAR.Alignment.run) {
    STARAlignmentDefault(path.prefix,
                         genome.name,
                         sample.pattern,
                         STAR.Index.num.parallel.threads,
                         STAR.Index.sjdbOverhang.Read.length,
                         STAR.Index.genomeSAindexNbases,
                         STAR.Index.genomeChrBinNbits,
                         STAR.Index.genomeSAsparseD)


    STARAlignmentDefault(path.prefix,
                         genome.name,
                         sample.pattern,
                         STAR.Alignment.num.parallel.threads,
                         STAR.Alignment.genomeLoad,
                         STAR.Alignment.readMapNumber,
                         STAR.Alignment.clip3pNbases,
                         STAR.Alignment.clip5pNbases,
                         STAR.Alignment.clip3pAdapterSeq,
                         STAR.Alignment.clip3pAdapterMMp,
                         STAR.Alignment.clip3pAfterAdapterNbases,
                         STAR.Alignment.limitGenomeGenerateRAM,
                         STAR.Alignment.limitIObufferSize,
                         STAR.Alignment.limitOutSAMoneReadBytes,
                         STAR.Alignment.limitOutSJoneRead,
                         STAR.Alignment.limitOutSJcollapsed,
                         STAR.Alignment.limitBAMsortRAM,
                         STAR.Alignment.outReadsUnmapped,
                         STAR.Alignment.outQSconversionAdd,
                         STAR.Alignment.outSAMprimaryFlag,
                         STAR.Alignment.outSAMmapqUnique,
                         STAR.Alignment.scoreGap,
                         STAR.Alignment.scoreGapNoncan,
                         STAR.Alignment.scoreGapGCAG,
                         STAR.Alignment.scoreGapATAC,
                         STAR.Alignment.scoreGenomicLengthLog2scale,
                         STAR.Alignment.scoreDelOpen,
                         STAR.Alignment.scoreDelBase,
                         STAR.Alignment.scoreInsOpen,
                         STAR.Alignment.scoreInsBase,
                         STAR.Alignment.scoreStitchSJshift,
                         STAR.Alignment.seedSearchStartLmax,
                         STAR.Alignment.seedSearchStartLmaxOverLread,
                         STAR.Alignment.seedSearchLmax,
                         STAR.Alignment.seedMultimapNmax,
                         STAR.Alignment.seedPerReadNmax,
                         STAR.Alignment.seedPerWindowNmax,
                         STAR.Alignment.seedNoneLociPerWindow,
                         STAR.Alignment.alignIntronMin,
                         STAR.Alignment.alignIntronMax,
                         STAR.Alignment.alignMatesGapMax,
                         STAR.Alignment.alignSJoverhangMin,
                         STAR.Alignment.alignSJDBoverhangMin,
                         STAR.Alignment.alignSplicedMateMapLmin,
                         STAR.Alignment.alignSplicedMateMapLminOverLmate,
                         STAR.Alignment.alignWindowsPerReadNmax,
                         STAR.Alignment.alignTranscriptsPerWindowNmax,
                         STAR.Alignment.alignTranscriptsPerReadNmax,
                         STAR.Alignment.alignEndsType,
                         STAR.Alignment.winAnchorMultimapNmax,
                         STAR.Alignment.winBinNbits,
                         STAR.Alignment.winAnchorDistNbins,
                         STAR.Alignment.winFlankNbins)
  }


  if (Rsamtools.Bam.run) {
    RSamtoolsToBam(SAMtools.or.Rsamtools,
                   Samtools.Bam.num.parallel.threads,
                   path.prefix,
                   genome.name,
                   sample.pattern,
                   Rsamtools.nCores)
  }
  if (StringTie.Assemble.run) {
    StringTieAssemble(path.prefix,
                      genome.name,
                      sample.pattern,
                      Stringtie.Assembly.num.parallel.threads,
                      Stringtie.Assembly.f,
                      Stringtie.Assembly.m,
                      Stringtie.Assembly.c,
                      Stringtie.Assembly.g,
                      Stringtie.Assembly.M)
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
  if (PreDECountTable.run &
      ((python.variable.answer & python.variable.version == 2) |
       python.variable.answer & python.variable.version == 3 & python.2to3)) {
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

