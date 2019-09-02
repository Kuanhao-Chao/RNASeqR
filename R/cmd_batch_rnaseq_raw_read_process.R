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
#' parameters
#' @param SAMtools.or.Rsamtools Default value is \code{Rsamtools}. User can set
#' to \code{SAMtools} to use command-line-based 'samtools' instead.
#' @param Hisat2.Index.run Whether to run 'HISAT2 index' step in this function
#' step. Default value is \code{TRUE}. Set \code{FALSE} to skip
#' 'HISAT2 index' step.
#' @param Hisat2.Index.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Hisat2 index step. The default is \code{"1"}
#' @param Hisat2.Index.large.index Hisat2 index terminal '--large-index' option.
#' Default value is \code{FALSE}
#' @param Hisat2.Index.local.ftab.chars Hisat2 index terminal '-t/--ftabchars'
#' option. Default value is \code{"6"}
#' @param Hisat2.Index.local.off.rate Hisat2 index terminal '--localoffrate'
#' option. Default value is \code{"3"}
#' @param Hisat2.Index.ftab.chars Hisat2 index terminal '--localftabchars'
#' option. Default value is \code{"10"}
#' @param Hisat2.Index.off.rate Hisat2 index terminal '--offrate' option.
#' Default value is \code{"4"}
#' @param Hisat2.Alignment.run Whether to run 'HISAT2 alignment' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2 alignment' step.
#' @param Hisat2.Alignment.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Hisat2 alignment step. The default is \code{"1"}
#' @param Hisat2.Alignment.skip Hisat2 alignment terminal '-s/--skip' option.
#' Default value is \code{"0"}
#' @param Hisat2.Alignment.trim5 Hisat2 alignment terminal '-5/--trim5' option.
#' Default value is \code{"0"}
#' @param Hisat2.Alignment.trim3 Hisat2 alignment terminal '-3/--trim3' option.
#' Default value is \code{"0"}
#' @param Hisat2.Alignment.n.ceil.1.function.type Hisat2 alignment terminal
#' '--n-ceil' option. Default value is \code{"L"}
#' @param Hisat2.Alignment.n.ceil.2.constant.term Hisat2 alignment terminal
#' '--n-ceil' option. Default value is \code{"0"}
#' @param Hisat2.Alignment.n.ceil.3.coefficient Hisat2 alignment terminal
#' '--n-ceil' option. Default value is \code{"0.15"}
#' @param Hisat2.Alignment.mp.MX Hisat2 alignment terminal '--mp MX' option.
#' Default value is \code{"6"}
#' @param Hisat2.Alignment.mp.MN Hisat2 alignment terminal '--mp MN' option.
#' Default value is \code{"2"}
#' @param Hisat2.Alignment.sp.MX Hisat2 alignment terminal '--sp MX' option.
#' Default value is \code{"2"}
#' @param Hisat2.Alignment.sp.MN Hisat2 alignment terminal '--sp MN' option.
#' Default value is \code{"1"}
#' @param Hisat2.Alignment.np Hisat2 alignment terminal '--np' option.
#' Default value is \code{"1"}
#' @param Hisat2.Alignment.rdg.1 Hisat2 alignment terminal '--rdg' first option.
#' Default value is \code{"5"}
#' @param Hisat2.Alignment.rdg.2 Hisat2 alignment terminal '--rdg' first option.
#' Default value is \code{"3"}
#' @param Hisat2.Alignment.rfg.1 Hisat2 alignment terminal '--rfg' first option.
#' Default value is \code{"5"}
#' @param Hisat2.Alignment.rfg.2 Hisat2 alignment terminal '--rfg' first option.
#' Default value is \code{"3"}
#' @param Hisat2.Alignment.score.min.1.function.type Hisat2 alignment terminal
#' '--rdg' first option. Default value is \code{"L"}
#' @param Hisat2.Alignment.score.min.2.constant.term Hisat2 alignment terminal
#' '--rdg' first option. Default value is \code{"0"}
#' @param Hisat2.Alignment.score.min.3.coefficient Hisat2 alignment terminal
#' '--rdg' first option. Default value is \code{"-0.2"}
#' @param Hisat2.Alignment.pen.cansplice Hisat2 alignment terminal
#' '--pen-cansplice' first option. Default value is \code{"-0"}
#' @param Hisat2.Alignment.penc.noncansplice Hisat2 alignment terminal
#' '--pen-noncansplice' option. Default value is \code{"12"}
#' @param Hisat2.Alignment.pen.canintronlen.1.function.type Hisat2 alignment
#' terminal '--pen-canintronlen' first option. Default value is \code{"G"}
#' @param Hisat2.Alignment.pen.canintronlen.2.constant.term Hisat2 alignment
#' terminal '--pen-canintronlen' second option. Default value is \code{"-8"}
#' @param Hisat2.Alignment.pen.canintronlen.3.coefficient Hisat2 alignment
#' terminal '--pen-canintronlen' third option. Default value is \code{"1"}
#' @param Hisat2.Alignment.pen.noncanintronlen.1.function.type Hisat2 alignment
#' terminal '--pen-noncanintronlen' first option. Default value is \code{"G"}
#' @param Hisat2.Alignment.pen.noncanintronlen.2.constant.term Hisat2 alignment
#' terminal '--pen-noncanintronlen' second option. Default value is \code{"-8"}
#' @param Hisat2.Alignment.pen.noncanintronlen.3.coefficient Hisat2 alignment
#' terminal '--pen-noncanintronlen' third option. Default value is \code{"1"}
#' @param Hisat2.Alignment.min.intronlen Hisat2 alignment terminal
#' '--min-intronlen' option. Default value is \code{"20"}
#' @param Hisat2.Alignment.max.intronlen Hisat2 alignment terminal '--max-intronlen'
#' option. Default value is \code{"20"}
#' @param Hisat2.Alignment.rna.strandness Hisat2 alignment terminal '--rna-strandness'
#' option. Default value is \code{"None"}
#' @param Hisat2.Alignment.k Hisat2 alignment terminal '-k' option.
#' Default value is \code{"5"}
#' @param Hisat2.Alignment.max.seeds Hisat2 alignment terminal '--max-seeds' option.
#' Default value is \code{"5"}
#' @param Hisat2.Alignment.secondary Hisat2 alignment terminal '--secondary' option.
#' Default value is \code{"FALSE"}
#' @param Hisat2.Alignment.minins Hisat2 alignment terminal '-I/--minins' option.
#' Default value is \code{"0"}
#' @param Hisat2.Alignment.maxins Hisat2 alignment terminal '-X/--maxins' option.
#' Default value is \code{"500"}
#' @param Hisat2.Alignment.seed Hisat2 alignment terminal '-X/--maxins' option.
#' Default value is \code{"0"}
#' @param STAR.Index.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for STAR index step. The default is \code{"1"}
#' @param STAR.Index.sjdbOverhang.Read.length STAR index terminal
#' '--sjdbOverhang' option. Default value is \code{"100"}
#' @param STAR.Index.genomeSAindexNbases STAR index terminal
#' '--genomeSAindexNbases' option. Default value is \code{"14"}
#' @param STAR.Index.genomeChrBinNbits STAR index terminal
#' '--genomeChrBinNbits' option. Default value is \code{"18"}
#' @param STAR.Index.genomeSAsparseD STAR index terminal
#' '--genomeSAsparseD' option. Default value is \code{"1"}
#' @param STAR.Alignment.run Whether to run 'STAR index' step in this function
#' step. Default value is \code{FALSE}. Set \code{TRUE} to run STAR alignment step.
#' (need to set Hisat2.Index.run to \code{FALSE})
#' @param STAR.Alignment.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for STAR alignment step. The default is \code{"1"}
#' @param STAR.Alignment.genomeLoad STAR alignment terminal '--genomeLoad'
#' option. Default value is \code{"NoSharedMemory"}
#' @param STAR.Alignment.readMapNumber STAR alignment terminal '--readMapNumber'
#' option. Default value is \code{"-1"}
#' @param STAR.Alignment.clip3pNbases STAR alignment terminal '--clip3pNbases'
#' option. Default value is \code{"0"}
#' @param STAR.Alignment.clip5pNbases STAR alignment terminal '--clip5pNbases'
#' option. Default value is \code{"0"}
#' @param STAR.Alignment.clip3pAdapterSeq STAR alignment terminal '--clip3pAdapterSeq'
#' option. Default value is \code{"-"}
#' @param STAR.Alignment.clip3pAdapterMMp STAR alignment terminal '--clip3pAdapterMMp'
#' option. Default value is \code{"0.1"}
#' @param STAR.Alignment.clip3pAfterAdapterNbases STAR alignment terminal
#' '--clip3pAfterAdapterNbases' option. Default value is \code{"0"}
#' @param STAR.Alignment.limitGenomeGenerateRAM  STAR alignment terminal
#' '--limitGenomeGenerateRAM' option. Default value is \code{"31000000000"}
#' @param STAR.Alignment.limitIObufferSize STAR alignment terminal
#' '--limitIObufferSize' option. Default value is \code{"150000000"}
#' @param STAR.Alignment.limitOutSAMoneReadBytes STAR alignment terminal
#' '--limitOutSAMoneReadBytes' option. Default value is \code{"100000"}
#' @param STAR.Alignment.limitOutSJoneRead STAR alignment terminal
#' '--limitOutSJoneRead' option. Default value is \code{"1000"}
#' @param STAR.Alignment.limitOutSJcollapsed STAR alignment terminal
#' '--limitOutSJcollapsed' option. Default value is \code{"1000000"}
#' @param STAR.Alignment.limitBAMsortRAM STAR alignment terminal
#' '--limitBAMsortRAM' option. Default value is \code{"0"}
#' @param STAR.Alignment.outReadsUnmapped STAR alignment terminal
#' '--outReadsUnmapped' option. Default value is \code{"None"}
#' @param STAR.Alignment.outQSconversionAdd STAR alignment terminal
#' '--outQSconversionAdd' option. Default value is \code{"0"}
#' @param STAR.Alignment.outSAMprimaryFlag STAR alignment terminal
#' '--outSAMprimaryFlag' option. Default value is \code{"OneBestScore"}
#' @param STAR.Alignment.outSAMmapqUnique STAR alignment terminal
#' '--outSAMmapqUnique' option. Default value is \code{"255"}
#' @param STAR.Alignment.scoreGap STAR alignment terminal
#' '--scoreGap' option. Default value is \code{"0"}
#' @param STAR.Alignment.scoreGapNoncan STAR alignment terminal
#' '--scoreGapNoncan' option. Default value is \code{"-8"}
#' @param STAR.Alignment.scoreGapGCAG STAR alignment terminal
#' '--scoreGapGCAG' option. Default value is \code{"-4"}
#' @param STAR.Alignment.scoreGapATAC STAR alignment terminal
#' '--scoreGapATAC' option. Default value is \code{"-8"}
#' @param STAR.Alignment.scoreGenomicLengthLog2scale STAR alignment terminal
#' '--scoreGenomicLengthLog2scale' option. Default value is \code{"-0.25"}
#' @param STAR.Alignment.scoreDelOpen STAR alignment terminal
#' '--scoreDelOpen' option. Default value is \code{"-2"}
#' @param STAR.Alignment.scoreDelBase STAR alignment terminal
#' '--scoreDelBase' option. Default value is \code{"-2"}
#' @param STAR.Alignment.scoreInsOpen STAR alignment terminal
#' '--scoreInsOpen' option. Default value is \code{"-2"}
#' @param STAR.Alignment.scoreInsBase STAR alignment terminal
#' '--scoreInsBase' option. Default value is \code{"-2"}
#' @param STAR.Alignment.scoreStitchSJshift STAR alignment terminal
#' '--scoreStitchSJshift' option. Default value is \code{"1"}
#' @param STAR.Alignment.seedSearchStartLmax STAR alignment terminal
#' '--scoreStitchSJshift' option. Default value is \code{"50"}
#' @param STAR.Alignment.seedSearchStartLmaxOverLread STAR alignment terminal
#' '--seedSearchStartLmaxOverLread' option. Default value is \code{"1.0"}
#' @param STAR.Alignment.seedSearchLmax STAR alignment terminal
#' '--seedSearchLmax' option. Default value is \code{"0"}
#' @param STAR.Alignment.seedMultimapNmax STAR alignment terminal
#' '--seedMultimapNmax' option. Default value is \code{"10000"}
#' @param STAR.Alignment.seedPerReadNmax STAR alignment terminal
#' '--seedPerReadNmax' option. Default value is \code{"1000"}
#' @param STAR.Alignment.seedPerWindowNmax STAR alignment terminal
#' '--seedPerWindowNmax' option. Default value is \code{"50"}
#' @param STAR.Alignment.seedNoneLociPerWindow STAR alignment terminal
#' '--seedNoneLociPerWindow' option. Default value is \code{"10"}
#' @param STAR.Alignment.alignIntronMin STAR alignment terminal
#' '--alignIntronMin' option. Default value is \code{"21"}
#' @param STAR.Alignment.alignIntronMax STAR alignment terminal
#' '--alignIntronMax' option. Default value is \code{"0"}
#' @param STAR.Alignment.alignMatesGapMax STAR alignment terminal
#' '--alignMatesGapMax' option. Default value is \code{"0"}
#' @param STAR.Alignment.alignSJoverhangMin STAR alignment terminal
#' '--alignSJoverhangMin' option. Default value is \code{"5"}
#' @param STAR.Alignment.alignSJDBoverhangMin STAR alignment terminal
#' '--alignSJDBoverhangMin' option. Default value is \code{"3"}
#' @param STAR.Alignment.alignSplicedMateMapLmin STAR alignment terminal
#' '--alignSplicedMateMapLmin' option. Default value is \code{"0"}
#' @param STAR.Alignment.alignSplicedMateMapLminOverLmate STAR alignment terminal
#' '--alignSplicedMateMapLminOverLmate' option. Default value is \code{"0.66"}
#' @param STAR.Alignment.alignWindowsPerReadNmax STAR alignment terminal
#' '--alignWindowsPerReadNmax' option. Default value is \code{"10000"}
#' @param STAR.Alignment.alignTranscriptsPerWindowNmax STAR alignment terminal
#' '--alignTranscriptsPerWindowNmax' option. Default value is \code{"100"}
#' @param STAR.Alignment.alignTranscriptsPerReadNmax STAR alignment terminal
#' '--alignTranscriptsPerReadNmax' option. Default value is \code{"10000"}
#' @param STAR.Alignment.alignEndsType STAR alignment terminal
#' '--alignEndsType' option. Default value is \code{"Local"}
#' @param STAR.Alignment.winAnchorMultimapNmax STAR alignment terminal
#' '--winAnchorMultimapNmax' option. Default value is \code{"50"}
#' @param STAR.Alignment.winBinNbits STAR alignment terminal
#' '--winBinNbits' option. Default value is \code{"16"}
#' @param STAR.Alignment.winAnchorDistNbins STAR alignment terminal
#' '--winAnchorDistNbins' option. Default value is \code{"9"}
#' @param STAR.Alignment.winFlankNbins STAR alignment terminal
#' '--winFlankNbins' option. Default value is \code{"4"}
#' @param Rsamtools.Bam.run Whether to run 'Rsamtools SAM to BAM' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'Rsamtools SAM to BAM' step.
#' @param Samtools.Bam.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Samtools sam to bam step. The default is \code{"1"}
#' @param Rsamtools.nCores The number of cores to use when running
#'  'Rsamtools' step. Default value is \code{1}
#' @param StringTie.Assemble.run Whether to run 'StringTie assembly' step in
#'   this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie assembly' step.
#' @param Stringtie.Assembly.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Stringtie assembly. The default is \code{"1"}
#' @param Stringtie.Assembly.f Stringtie assembly terminal
#' '-f' option. Default value is \code{"0.1"}
#' @param Stringtie.Assembly.m Stringtie assembly terminal
#' '-m' option. Default value is \code{"200"}
#' @param Stringtie.Assembly.c Stringtie assembly terminal
#' '-c' option. Default value is \code{"2.5"}
#' @param Stringtie.Assembly.g Stringtie assembly terminal
#' '-g' option. Default value is \code{"50"}
#' @param Stringtie.Assembly.M Stringtie assembly terminal
#' '-M' option. Default value is \code{"0.95"}
#' @param StringTie.Merge.Trans.run Whether to run 'StringTie GTF merging' step
#'   in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie GTF merging' step.
#' @param Stringtie.Merge.num.parallel.threads  Specify the number of processing
#' threads (CPUs) to use for Stringtie merge step. The default is \code{"1"}
#' @param Gffcompare.Ref.Sample.run Whether to run 'Gffcompare comparison' step
#'   in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'Gffcompare comparison' step.
#' @param StringTie.Ballgown.run Whether to run 'StringTie ballgown creation'
#'   step in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie ballgown creation' step.
#' @param Stringtie.2.Ballgown.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Stringtie to ballgown step. The default is \code{"1"}
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
#'
# Parameters: 119
RNASeqReadProcess_CMD <- function(RNASeqRParam,
                                  SAMtools.or.Rsamtools     = "Rsamtools",
                                  Hisat2.Index.run          = TRUE,
                                  Hisat2.Index.num.parallel.threads = "1",
                                  Hisat2.Index.large.index = FALSE,
                                  Hisat2.Index.local.ftab.chars = "6",
                                  Hisat2.Index.local.off.rate = "3",
                                  Hisat2.Index.ftab.chars = "10",
                                  Hisat2.Index.off.rate = "4",
                                  Hisat2.Alignment.run      = TRUE,
                                  Hisat2.Alignment.num.parallel.threads = "1",
                                  Hisat2.Alignment.skip = "0",
                                  Hisat2.Alignment.trim5 = "0",
                                  Hisat2.Alignment.trim3 = "0",
                                  Hisat2.Alignment.n.ceil.1.function.type = "L",
                                  Hisat2.Alignment.n.ceil.2.constant.term = "0",
                                  Hisat2.Alignment.n.ceil.3.coefficient = "0.15",
                                  Hisat2.Alignment.mp.MX = "6",
                                  Hisat2.Alignment.mp.MN = "2",
                                  Hisat2.Alignment.sp.MX = "2",
                                  Hisat2.Alignment.sp.MN = "1",
                                  Hisat2.Alignment.np = "1",
                                  Hisat2.Alignment.rdg.1 = "5",
                                  Hisat2.Alignment.rdg.2 = "3",
                                  Hisat2.Alignment.rfg.1 = "5",
                                  Hisat2.Alignment.rfg.2 = "3",
                                  Hisat2.Alignment.score.min.1.function.type = "L",
                                  Hisat2.Alignment.score.min.2.constant.term = "0",
                                  Hisat2.Alignment.score.min.3.coefficient = "-0.2",
                                  Hisat2.Alignment.pen.cansplice = "0",
                                  Hisat2.Alignment.penc.noncansplice = "12",
                                  Hisat2.Alignment.pen.canintronlen.1.function.type = "G",
                                  Hisat2.Alignment.pen.canintronlen.2.constant.term = "-8",
                                  Hisat2.Alignment.pen.canintronlen.3.coefficient = "1",
                                  Hisat2.Alignment.pen.noncanintronlen.1.function.type = "G",
                                  Hisat2.Alignment.pen.noncanintronlen.2.constant.term = "-8",
                                  Hisat2.Alignment.pen.noncanintronlen.3.coefficient = "1",
                                  Hisat2.Alignment.min.intronlen = "20",
                                  Hisat2.Alignment.max.intronlen = "500000",
                                  Hisat2.Alignment.rna.strandness = "None",
                                  Hisat2.Alignment.k = "5",
                                  Hisat2.Alignment.max.seeds = "5",
                                  Hisat2.Alignment.secondary = FALSE,
                                  Hisat2.Alignment.minins = "0",
                                  Hisat2.Alignment.maxins = "500",
                                  Hisat2.Alignment.seed = "0",
                                  STAR.Index.num.parallel.threads = "1",
                                  STAR.Index.sjdbOverhang.Read.length = "100",
                                  STAR.Index.genomeSAindexNbases = "14",
                                  STAR.Index.genomeChrBinNbits = "18",
                                  STAR.Index.genomeSAsparseD = "1",
                                  STAR.Alignment.run        = FALSE,
                                  STAR.Alignment.num.parallel.threads = "1",
                                  STAR.Alignment.genomeLoad = "NoSharedMemory",
                                  STAR.Alignment.readMapNumber = "-1",
                                  STAR.Alignment.clip3pNbases = "0",
                                  STAR.Alignment.clip5pNbases = "0",
                                  STAR.Alignment.clip3pAdapterSeq = "-",
                                  STAR.Alignment.clip3pAdapterMMp = "0.1",
                                  STAR.Alignment.clip3pAfterAdapterNbases = "0",
                                  STAR.Alignment.limitGenomeGenerateRAM = "31000000000",
                                  STAR.Alignment.limitIObufferSize = "150000000",
                                  STAR.Alignment.limitOutSAMoneReadBytes = "100000",
                                  STAR.Alignment.limitOutSJoneRead = "1000",
                                  STAR.Alignment.limitOutSJcollapsed = "1000000",
                                  STAR.Alignment.limitBAMsortRAM = "0",
                                  STAR.Alignment.outReadsUnmapped = "None",
                                  STAR.Alignment.outQSconversionAdd = "0",
                                  STAR.Alignment.outSAMprimaryFlag = "OneBestScore",
                                  STAR.Alignment.outSAMmapqUnique = "255",
                                  STAR.Alignment.scoreGap = "0",
                                  STAR.Alignment.scoreGapNoncan = "-8",
                                  STAR.Alignment.scoreGapGCAG = "-4",
                                  STAR.Alignment.scoreGapATAC = "-8",
                                  STAR.Alignment.scoreGenomicLengthLog2scale = "-0.25",
                                  STAR.Alignment.scoreDelOpen = "-2",
                                  STAR.Alignment.scoreDelBase = "-2",
                                  STAR.Alignment.scoreInsOpen = "-2",
                                  STAR.Alignment.scoreInsBase = "-2",
                                  STAR.Alignment.scoreStitchSJshift = "1",
                                  STAR.Alignment.seedSearchStartLmax = "50",
                                  STAR.Alignment.seedSearchStartLmaxOverLread = "1.0",
                                  STAR.Alignment.seedSearchLmax = "0",
                                  STAR.Alignment.seedMultimapNmax = "10000",
                                  STAR.Alignment.seedPerReadNmax = "1000",
                                  STAR.Alignment.seedPerWindowNmax = "50",
                                  STAR.Alignment.seedNoneLociPerWindow = "10",
                                  STAR.Alignment.alignIntronMin = "21",
                                  STAR.Alignment.alignIntronMax = "0",
                                  STAR.Alignment.alignMatesGapMax = "0",
                                  STAR.Alignment.alignSJoverhangMin = "5",
                                  STAR.Alignment.alignSJDBoverhangMin = "3",
                                  STAR.Alignment.alignSplicedMateMapLmin = "0",
                                  STAR.Alignment.alignSplicedMateMapLminOverLmate = "0.66",
                                  STAR.Alignment.alignWindowsPerReadNmax = "10000",
                                  STAR.Alignment.alignTranscriptsPerWindowNmax = "100",
                                  STAR.Alignment.alignTranscriptsPerReadNmax = "10000",
                                  STAR.Alignment.alignEndsType = "Local",
                                  STAR.Alignment.winAnchorMultimapNmax = "50",
                                  STAR.Alignment.winBinNbits = "16",
                                  STAR.Alignment.winAnchorDistNbins = "9",
                                  STAR.Alignment.winFlankNbins = "4",
                                  Rsamtools.Bam.run         = TRUE,
                                  Samtools.Bam.num.parallel.threads = "1",
                                  Rsamtools.nCores          = "1",
                                  StringTie.Assemble.run    = TRUE,
                                  Stringtie.Assembly.num.parallel.threads = "1",
                                  Stringtie.Assembly.f = "0.1",
                                  Stringtie.Assembly.m = "200",
                                  Stringtie.Assembly.c = "2.5",
                                  Stringtie.Assembly.g = "50",
                                  Stringtie.Assembly.M = "0.95",
                                  StringTie.Merge.Trans.run = TRUE,
                                  Stringtie.Merge.num.parallel.threads = "1",
                                  Gffcompare.Ref.Sample.run = TRUE,
                                  StringTie.Ballgown.run    = TRUE,
                                  Stringtie.2.Ballgown.num.parallel.threads = "1",
                                  PreDECountTable.run       = TRUE,
                                  run                       = TRUE,
                                  check.s4.print            = TRUE) {
  which.s4.object <- CheckS4Object_All(RNASeqRParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  INSIDE.path.prefix <- "@"(RNASeqRParam, path.prefix)
  saveRDS(RNASeqRParam,
          file = paste0(INSIDE.path.prefix,
                        "gene_data/RNASeqRParam.rds"))
  fileConn<-file(paste0(path.prefix, "Rscript/Read_Process.R"))
  first <- "library(RNASeqR)"
  if (which.s4.object == "RNASeqRParam") {
  } else if (which.s4.object == "RNASeqRParam_Sam") {
    Hisat2.Index.run = FALSE
    Hisat2.Alignment.run = FALSE
    STAR.Alignment.run = FALSE
  } else if (which.s4.object == "RNASeqRParam_Bam") {
    Hisat2.Index.run = FALSE
    Hisat2.Alignment.run = FALSE
    STAR.Alignment.run = FALSE
    Rsamtools.Bam.run = FALSE
  }
  second <- paste0("RNASeqReadProcess(RNASeqRParam = 'INSIDE'",
                   ", which.trigger = 'INSIDE'",
                   ", INSIDE.path.prefix = '", INSIDE.path.prefix,
                   "', SAMtools.or.Rsamtools ='", SAMtools.or.Rsamtools,
                   "', Hisat2.Index.run =", Hisat2.Index.run,
                   ", Hisat2.Index.num.parallel.threads = '", Hisat2.Index.num.parallel.threads,
                   "', Hisat2.Index.large.index = ", Hisat2.Index.large.index,
                   ", Hisat2.Index.local.ftab.chars = '", Hisat2.Index.local.ftab.chars,
                   "', Hisat2.Index.local.off.rate = '", Hisat2.Index.local.off.rate,
                   "', Hisat2.Index.ftab.chars = '", Hisat2.Index.ftab.chars,
                   "', Hisat2.Index.off.rate = '", Hisat2.Index.off.rate,
                   "', Hisat2.Alignment.run  = ", Hisat2.Alignment.run,
                   ", Hisat2.Alignment.num.parallel.threads = '", Hisat2.Alignment.num.parallel.threads,
                   "', Hisat2.Alignment.skip = '", Hisat2.Alignment.skip,
                   "', Hisat2.Alignment.trim5 = '", Hisat2.Alignment.trim5,
                   "', Hisat2.Alignment.trim3 = '", Hisat2.Alignment.trim3,
                   "', Hisat2.Alignment.n.ceil.1.function.type = '", Hisat2.Alignment.n.ceil.1.function.type,
                   "', Hisat2.Alignment.n.ceil.2.constant.term = '", Hisat2.Alignment.n.ceil.2.constant.term,
                   "', Hisat2.Alignment.n.ceil.3.coefficient = '", Hisat2.Alignment.n.ceil.3.coefficient,
                   "', Hisat2.Alignment.mp.MX = '", Hisat2.Alignment.mp.MX,
                   "', Hisat2.Alignment.mp.MN = '", Hisat2.Alignment.mp.MN,
                   "', Hisat2.Alignment.sp.MX = '", Hisat2.Alignment.sp.MX,
                   "', Hisat2.Alignment.sp.MN = '", Hisat2.Alignment.sp.MN,
                   "', Hisat2.Alignment.np = '", Hisat2.Alignment.np,
                   "', Hisat2.Alignment.rdg.1 = '", Hisat2.Alignment.rdg.1,
                   "', Hisat2.Alignment.rdg.2 = '", Hisat2.Alignment.rdg.2,
                   "', Hisat2.Alignment.rfg.1 = '", Hisat2.Alignment.rfg.1,
                   "', Hisat2.Alignment.rfg.2 = '", Hisat2.Alignment.rfg.2,
                   "', Hisat2.Alignment.score.min.1.function.type = '", Hisat2.Alignment.score.min.1.function.type,
                   "', Hisat2.Alignment.score.min.2.constant.term = '", Hisat2.Alignment.score.min.2.constant.term,
                   "', Hisat2.Alignment.score.min.3.coefficient = '", Hisat2.Alignment.score.min.3.coefficient,
                   "', Hisat2.Alignment.pen.cansplice = '", Hisat2.Alignment.pen.cansplice,
                   "', Hisat2.Alignment.penc.noncansplice = '", Hisat2.Alignment.penc.noncansplice,
                   "', Hisat2.Alignment.pen.canintronlen.1.function.type = '", Hisat2.Alignment.pen.canintronlen.1.function.type,
                   "', Hisat2.Alignment.pen.canintronlen.2.constant.term = '", Hisat2.Alignment.pen.canintronlen.2.constant.term,
                   "', Hisat2.Alignment.pen.canintronlen.3.coefficient = '", Hisat2.Alignment.pen.canintronlen.3.coefficient,
                   "', Hisat2.Alignment.pen.noncanintronlen.1.function.type = '", Hisat2.Alignment.pen.noncanintronlen.1.function.type,
                   "', Hisat2.Alignment.pen.noncanintronlen.2.constant.term = '", Hisat2.Alignment.pen.noncanintronlen.2.constant.term,
                   "', Hisat2.Alignment.pen.noncanintronlen.3.coefficient = '", Hisat2.Alignment.pen.noncanintronlen.3.coefficient,
                   "', Hisat2.Alignment.min.intronlen = '", Hisat2.Alignment.min.intronlen,
                   "', Hisat2.Alignment.max.intronlen = '", Hisat2.Alignment.max.intronlen,
                   "', Hisat2.Alignment.rna.strandness = '", Hisat2.Alignment.rna.strandness,
                   "', Hisat2.Alignment.k = '", Hisat2.Alignment.k,
                   "', Hisat2.Alignment.max.seeds = '", Hisat2.Alignment.max.seeds,
                   "', Hisat2.Alignment.secondary = ", Hisat2.Alignment.secondary,
                   ", Hisat2.Alignment.minins = '", Hisat2.Alignment.minins,
                   "', Hisat2.Alignment.maxins = '", Hisat2.Alignment.maxins,
                   "', Hisat2.Alignment.seed = '", Hisat2.Alignment.seed,
                   "', STAR.Index.num.parallel.threads = '", STAR.Index.num.parallel.threads,
                   "', STAR.Index.sjdbOverhang.Read.length = '", STAR.Index.sjdbOverhang.Read.length,
                   "', STAR.Index.genomeSAindexNbases = '", STAR.Index.genomeSAindexNbases,
                   "', STAR.Index.genomeChrBinNbits = '", STAR.Index.genomeChrBinNbits,
                   "', STAR.Index.genomeSAsparseD = '", STAR.Index.genomeSAsparseD,
                   "', STAR.Alignment.run = ", STAR.Alignment.run,
                   ", STAR.Alignment.num.parallel.threads = '", STAR.Alignment.num.parallel.threads,
                   "', STAR.Alignment.genomeLoad = '", STAR.Alignment.genomeLoad,
                   "', STAR.Alignment.readMapNumber = '", STAR.Alignment.readMapNumber,
                   "', STAR.Alignment.clip3pNbases = '", STAR.Alignment.clip3pNbases,
                   "', STAR.Alignment.clip5pNbases = '", STAR.Alignment.clip5pNbases,
                   "', STAR.Alignment.clip3pAdapterSeq = '", STAR.Alignment.clip3pAdapterSeq,
                   "', STAR.Alignment.clip3pAdapterMMp = '", STAR.Alignment.clip3pAdapterMMp,
                   "', STAR.Alignment.clip3pAfterAdapterNbases = '", STAR.Alignment.clip3pAfterAdapterNbases,
                   "', STAR.Alignment.limitGenomeGenerateRAM = '", STAR.Alignment.limitGenomeGenerateRAM,
                   "', STAR.Alignment.limitIObufferSize = '", STAR.Alignment.limitIObufferSize,
                   "', STAR.Alignment.limitOutSAMoneReadBytes = '", STAR.Alignment.limitOutSAMoneReadBytes,
                   "', STAR.Alignment.limitOutSJoneRead = '", STAR.Alignment.limitOutSJoneRead,
                   "', STAR.Alignment.limitOutSJcollapsed = '", STAR.Alignment.limitOutSJcollapsed,
                   "', STAR.Alignment.limitBAMsortRAM = '", STAR.Alignment.limitBAMsortRAM,
                   "', STAR.Alignment.outReadsUnmapped = '", STAR.Alignment.outReadsUnmapped,
                   "', STAR.Alignment.outQSconversionAdd = '", STAR.Alignment.outQSconversionAdd,
                   "', STAR.Alignment.outSAMprimaryFlag = '", STAR.Alignment.outSAMprimaryFlag,
                   "', STAR.Alignment.outSAMmapqUnique = '", STAR.Alignment.outSAMmapqUnique,
                   "', STAR.Alignment.scoreGap = '", STAR.Alignment.scoreGap,
                   "', STAR.Alignment.scoreGapNoncan = '", STAR.Alignment.scoreGapNoncan,
                   "', STAR.Alignment.scoreGapGCAG = '", STAR.Alignment.scoreGapGCAG,
                   "', STAR.Alignment.scoreGapATAC = '", STAR.Alignment.scoreGapATAC,
                   "', STAR.Alignment.scoreGenomicLengthLog2scale = '", STAR.Alignment.scoreGenomicLengthLog2scale,
                   "', STAR.Alignment.scoreDelOpen = '", STAR.Alignment.scoreDelOpen,
                   "', STAR.Alignment.scoreDelBase = '", STAR.Alignment.scoreDelBase,
                   "', STAR.Alignment.scoreInsOpen = '", STAR.Alignment.scoreInsOpen,
                   "', STAR.Alignment.scoreInsBase = '", STAR.Alignment.scoreInsBase,
                   "', STAR.Alignment.scoreStitchSJshift = '", STAR.Alignment.scoreStitchSJshift,
                   "', STAR.Alignment.seedSearchStartLmax = '", STAR.Alignment.seedSearchStartLmax,
                   "', STAR.Alignment.seedSearchStartLmaxOverLread = '", STAR.Alignment.seedSearchStartLmaxOverLread,
                   "', STAR.Alignment.seedSearchLmax = '", STAR.Alignment.seedSearchLmax,
                   "', STAR.Alignment.seedMultimapNmax = '", STAR.Alignment.seedMultimapNmax,
                   "', STAR.Alignment.seedPerReadNmax = '", STAR.Alignment.seedPerReadNmax,
                   "', STAR.Alignment.seedPerWindowNmax = '", STAR.Alignment.seedPerWindowNmax,
                   "', STAR.Alignment.seedNoneLociPerWindow = '", STAR.Alignment.seedNoneLociPerWindow,
                   "', STAR.Alignment.alignIntronMin = '", STAR.Alignment.alignIntronMin,
                   "', STAR.Alignment.alignIntronMax = '", STAR.Alignment.alignIntronMax,
                   "', STAR.Alignment.alignMatesGapMax = '", STAR.Alignment.alignMatesGapMax,
                   "', STAR.Alignment.alignSJoverhangMin = '", STAR.Alignment.alignSJoverhangMin,
                   "', STAR.Alignment.alignSJDBoverhangMin = '", STAR.Alignment.alignSJDBoverhangMin,
                   "', STAR.Alignment.alignSplicedMateMapLmin = '", STAR.Alignment.alignSplicedMateMapLmin,
                   "', STAR.Alignment.alignSplicedMateMapLminOverLmate = '", STAR.Alignment.alignSplicedMateMapLminOverLmate,
                   "', STAR.Alignment.alignWindowsPerReadNmax = '", STAR.Alignment.alignWindowsPerReadNmax,
                   "', STAR.Alignment.alignTranscriptsPerWindowNmax = '", STAR.Alignment.alignTranscriptsPerWindowNmax,
                   "', STAR.Alignment.alignTranscriptsPerReadNmax = '", STAR.Alignment.alignTranscriptsPerReadNmax,
                   "', STAR.Alignment.alignEndsType = '", STAR.Alignment.alignEndsType,
                   "', STAR.Alignment.winAnchorMultimapNmax = '", STAR.Alignment.winAnchorMultimapNmax,
                   "', STAR.Alignment.winBinNbits = '", STAR.Alignment.winBinNbits,
                   "', STAR.Alignment.winAnchorDistNbins = '", STAR.Alignment.winAnchorDistNbins,
                   "', STAR.Alignment.winFlankNbins = '", STAR.Alignment.winFlankNbins,
                   "', Rsamtools.Bam.run = ", Rsamtools.Bam.run,
                   ", Samtools.Bam.num.parallel.threads = '", Samtools.Bam.num.parallel.threads,
                   "', Rsamtools.nCores = '", Rsamtools.nCores,
                   "', StringTie.Assemble.run = ", StringTie.Assemble.run,
                   ", Stringtie.Assembly.num.parallel.threads = '", Stringtie.Assembly.num.parallel.threads,
                   "', Stringtie.Assembly.f = '", Stringtie.Assembly.f,
                   "', Stringtie.Assembly.m = '", Stringtie.Assembly.m,
                   "', Stringtie.Assembly.c = '", Stringtie.Assembly.c,
                   "', Stringtie.Assembly.g = '", Stringtie.Assembly.g,
                   "', Stringtie.Assembly.M = '", Stringtie.Assembly.M,
                   "', StringTie.Merge.Trans.run = ", StringTie.Merge.Trans.run,
                   ", Stringtie.Merge.num.parallel.threads = '", Stringtie.Merge.num.parallel.threads,
                   "', Gffcompare.Ref.Sample.run = ", Gffcompare.Ref.Sample.run,
                   ", StringTie.Ballgown.run = ", StringTie.Ballgown.run,
                   ", Stringtie.2.Ballgown.num.parallel.threads = '", Stringtie.2.Ballgown.num.parallel.threads,
                   "', PreDECountTable.run = ", PreDECountTable.run,
                   ", check.s4.print = ", check.s4.print, ")")
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
#' @param Hisat2.Index.run Whether to run 'HISAT2 index' step in this function
#' step. Default value is \code{TRUE}. Set \code{FALSE} to skip
#' 'HISAT2 index' step.
#' @param Hisat2.Index.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Hisat2 index step. The default is \code{"1"}
#' @param Hisat2.Index.large.index Hisat2 index terminal '--large-index' option.
#' Default value is \code{FALSE}
#' @param Hisat2.Index.local.ftab.chars Hisat2 index terminal '-t/--ftabchars'
#' option. Default value is \code{"6"}
#' @param Hisat2.Index.local.off.rate Hisat2 index terminal '--localoffrate'
#' option. Default value is \code{"3"}
#' @param Hisat2.Index.ftab.chars Hisat2 index terminal '--localftabchars'
#' option. Default value is \code{"10"}
#' @param Hisat2.Index.off.rate Hisat2 index terminal '--offrate' option.
#' Default value is \code{"4"}
#' @param Hisat2.Alignment.run Whether to run 'HISAT2 alignment' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2 alignment' step.
#' @param Hisat2.Alignment.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Hisat2 alignment step. The default is \code{"1"}
#' @param Hisat2.Alignment.skip Hisat2 alignment terminal '-s/--skip' option.
#' Default value is \code{"0"}
#' @param Hisat2.Alignment.trim5 Hisat2 alignment terminal '-5/--trim5' option.
#' Default value is \code{"0"}
#' @param Hisat2.Alignment.trim3 Hisat2 alignment terminal '-3/--trim3' option.
#' Default value is \code{"0"}
#' @param Hisat2.Alignment.n.ceil.1.function.type Hisat2 alignment terminal
#' '--n-ceil' option. Default value is \code{"L"}
#' @param Hisat2.Alignment.n.ceil.2.constant.term Hisat2 alignment terminal
#' '--n-ceil' option. Default value is \code{"0"}
#' @param Hisat2.Alignment.n.ceil.3.coefficient Hisat2 alignment terminal
#' '--n-ceil' option. Default value is \code{"0.15"}
#' @param Hisat2.Alignment.mp.MX Hisat2 alignment terminal '--mp MX' option.
#' Default value is \code{"6"}
#' @param Hisat2.Alignment.mp.MN Hisat2 alignment terminal '--mp MN' option.
#' Default value is \code{"2"}
#' @param Hisat2.Alignment.sp.MX Hisat2 alignment terminal '--sp MX' option.
#' Default value is \code{"2"}
#' @param Hisat2.Alignment.sp.MN Hisat2 alignment terminal '--sp MN' option.
#' Default value is \code{"1"}
#' @param Hisat2.Alignment.np Hisat2 alignment terminal '--np' option.
#' Default value is \code{"1"}
#' @param Hisat2.Alignment.rdg.1 Hisat2 alignment terminal '--rdg' first option.
#' Default value is \code{"5"}
#' @param Hisat2.Alignment.rdg.2 Hisat2 alignment terminal '--rdg' first option.
#' Default value is \code{"3"}
#' @param Hisat2.Alignment.rfg.1 Hisat2 alignment terminal '--rfg' first option.
#' Default value is \code{"5"}
#' @param Hisat2.Alignment.rfg.2 Hisat2 alignment terminal '--rfg' first option.
#' Default value is \code{"3"}
#' @param Hisat2.Alignment.score.min.1.function.type Hisat2 alignment terminal
#' '--rdg' first option. Default value is \code{"L"}
#' @param Hisat2.Alignment.score.min.2.constant.term Hisat2 alignment terminal
#' '--rdg' first option. Default value is \code{"0"}
#' @param Hisat2.Alignment.score.min.3.coefficient Hisat2 alignment terminal
#' '--rdg' first option. Default value is \code{"-0.2"}
#' @param Hisat2.Alignment.pen.cansplice Hisat2 alignment terminal
#' '--pen-cansplice' first option. Default value is \code{"-0"}
#' @param Hisat2.Alignment.penc.noncansplice Hisat2 alignment terminal
#' '--pen-noncansplice' option. Default value is \code{"12"}
#' @param Hisat2.Alignment.pen.canintronlen.1.function.type Hisat2 alignment
#' terminal '--pen-canintronlen' first option. Default value is \code{"G"}
#' @param Hisat2.Alignment.pen.canintronlen.2.constant.term Hisat2 alignment
#' terminal '--pen-canintronlen' second option. Default value is \code{"-8"}
#' @param Hisat2.Alignment.pen.canintronlen.3.coefficient Hisat2 alignment
#' terminal '--pen-canintronlen' third option. Default value is \code{"1"}
#' @param Hisat2.Alignment.pen.noncanintronlen.1.function.type Hisat2 alignment
#' terminal '--pen-noncanintronlen' first option. Default value is \code{"G"}
#' @param Hisat2.Alignment.pen.noncanintronlen.2.constant.term Hisat2 alignment
#' terminal '--pen-noncanintronlen' second option. Default value is \code{"-8"}
#' @param Hisat2.Alignment.pen.noncanintronlen.3.coefficient Hisat2 alignment
#' terminal '--pen-noncanintronlen' third option. Default value is \code{"1"}
#' @param Hisat2.Alignment.min.intronlen Hisat2 alignment terminal
#' '--min-intronlen' option. Default value is \code{"20"}
#' @param Hisat2.Alignment.max.intronlen Hisat2 alignment terminal '--max-intronlen'
#' option. Default value is \code{"20"}
#' @param Hisat2.Alignment.rna.strandness Hisat2 alignment terminal '--rna-strandness'
#' option. Default value is \code{"None"}
#' @param Hisat2.Alignment.k Hisat2 alignment terminal '-k' option.
#' Default value is \code{"5"}
#' @param Hisat2.Alignment.max.seeds Hisat2 alignment terminal '--max-seeds' option.
#' Default value is \code{"5"}
#' @param Hisat2.Alignment.secondary Hisat2 alignment terminal '--secondary' option.
#' Default value is \code{"FALSE"}
#' @param Hisat2.Alignment.minins Hisat2 alignment terminal '-I/--minins' option.
#' Default value is \code{"0"}
#' @param Hisat2.Alignment.maxins Hisat2 alignment terminal '-X/--maxins' option.
#' Default value is \code{"500"}
#' @param Hisat2.Alignment.seed Hisat2 alignment terminal '-X/--maxins' option.
#' Default value is \code{"0"}
#' @param STAR.Index.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for STAR index step. The default is \code{"1"}
#' @param STAR.Index.sjdbOverhang.Read.length STAR index terminal
#' '--sjdbOverhang' option. Default value is \code{"100"}
#' @param STAR.Index.genomeSAindexNbases STAR index terminal
#' '--genomeSAindexNbases' option. Default value is \code{"14"}
#' @param STAR.Index.genomeChrBinNbits STAR index terminal
#' '--genomeChrBinNbits' option. Default value is \code{"18"}
#' @param STAR.Index.genomeSAsparseD STAR index terminal
#' '--genomeSAsparseD' option. Default value is \code{"1"}
#' @param STAR.Alignment.run Whether to run 'STAR index' step in this function
#' step. Default value is \code{FALSE}. Set \code{TRUE} to run STAR alignment step.
#' (need to set Hisat2.Index.run to \code{FALSE})
#' @param STAR.Alignment.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for STAR alignment step. The default is \code{"1"}
#' @param STAR.Alignment.genomeLoad STAR alignment terminal '--genomeLoad'
#' option. Default value is \code{"NoSharedMemory"}
#' @param STAR.Alignment.readMapNumber STAR alignment terminal '--readMapNumber'
#' option. Default value is \code{"-1"}
#' @param STAR.Alignment.clip3pNbases STAR alignment terminal '--clip3pNbases'
#' option. Default value is \code{"0"}
#' @param STAR.Alignment.clip5pNbases STAR alignment terminal '--clip5pNbases'
#' option. Default value is \code{"0"}
#' @param STAR.Alignment.clip3pAdapterSeq STAR alignment terminal '--clip3pAdapterSeq'
#' option. Default value is \code{"-"}
#' @param STAR.Alignment.clip3pAdapterMMp STAR alignment terminal '--clip3pAdapterMMp'
#' option. Default value is \code{"0.1"}
#' @param STAR.Alignment.clip3pAfterAdapterNbases STAR alignment terminal
#' '--clip3pAfterAdapterNbases' option. Default value is \code{"0"}
#' @param STAR.Alignment.limitGenomeGenerateRAM  STAR alignment terminal
#' '--limitGenomeGenerateRAM' option. Default value is \code{"31000000000"}
#' @param STAR.Alignment.limitIObufferSize STAR alignment terminal
#' '--limitIObufferSize' option. Default value is \code{"150000000"}
#' @param STAR.Alignment.limitOutSAMoneReadBytes STAR alignment terminal
#' '--limitOutSAMoneReadBytes' option. Default value is \code{"100000"}
#' @param STAR.Alignment.limitOutSJoneRead STAR alignment terminal
#' '--limitOutSJoneRead' option. Default value is \code{"1000"}
#' @param STAR.Alignment.limitOutSJcollapsed STAR alignment terminal
#' '--limitOutSJcollapsed' option. Default value is \code{"1000000"}
#' @param STAR.Alignment.limitBAMsortRAM STAR alignment terminal
#' '--limitBAMsortRAM' option. Default value is \code{"0"}
#' @param STAR.Alignment.outReadsUnmapped STAR alignment terminal
#' '--outReadsUnmapped' option. Default value is \code{"None"}
#' @param STAR.Alignment.outQSconversionAdd STAR alignment terminal
#' '--outQSconversionAdd' option. Default value is \code{"0"}
#' @param STAR.Alignment.outSAMprimaryFlag STAR alignment terminal
#' '--outSAMprimaryFlag' option. Default value is \code{"OneBestScore"}
#' @param STAR.Alignment.outSAMmapqUnique STAR alignment terminal
#' '--outSAMmapqUnique' option. Default value is \code{"255"}
#' @param STAR.Alignment.scoreGap STAR alignment terminal
#' '--scoreGap' option. Default value is \code{"0"}
#' @param STAR.Alignment.scoreGapNoncan STAR alignment terminal
#' '--scoreGapNoncan' option. Default value is \code{"-8"}
#' @param STAR.Alignment.scoreGapGCAG STAR alignment terminal
#' '--scoreGapGCAG' option. Default value is \code{"-4"}
#' @param STAR.Alignment.scoreGapATAC STAR alignment terminal
#' '--scoreGapATAC' option. Default value is \code{"-8"}
#' @param STAR.Alignment.scoreGenomicLengthLog2scale STAR alignment terminal
#' '--scoreGenomicLengthLog2scale' option. Default value is \code{"-0.25"}
#' @param STAR.Alignment.scoreDelOpen STAR alignment terminal
#' '--scoreDelOpen' option. Default value is \code{"-2"}
#' @param STAR.Alignment.scoreDelBase STAR alignment terminal
#' '--scoreDelBase' option. Default value is \code{"-2"}
#' @param STAR.Alignment.scoreInsOpen STAR alignment terminal
#' '--scoreInsOpen' option. Default value is \code{"-2"}
#' @param STAR.Alignment.scoreInsBase STAR alignment terminal
#' '--scoreInsBase' option. Default value is \code{"-2"}
#' @param STAR.Alignment.scoreStitchSJshift STAR alignment terminal
#' '--scoreStitchSJshift' option. Default value is \code{"1"}
#' @param STAR.Alignment.seedSearchStartLmax STAR alignment terminal
#' '--scoreStitchSJshift' option. Default value is \code{"50"}
#' @param STAR.Alignment.seedSearchStartLmaxOverLread STAR alignment terminal
#' '--seedSearchStartLmaxOverLread' option. Default value is \code{"1.0"}
#' @param STAR.Alignment.seedSearchLmax STAR alignment terminal
#' '--seedSearchLmax' option. Default value is \code{"0"}
#' @param STAR.Alignment.seedMultimapNmax STAR alignment terminal
#' '--seedMultimapNmax' option. Default value is \code{"10000"}
#' @param STAR.Alignment.seedPerReadNmax STAR alignment terminal
#' '--seedPerReadNmax' option. Default value is \code{"1000"}
#' @param STAR.Alignment.seedPerWindowNmax STAR alignment terminal
#' '--seedPerWindowNmax' option. Default value is \code{"50"}
#' @param STAR.Alignment.seedNoneLociPerWindow STAR alignment terminal
#' '--seedNoneLociPerWindow' option. Default value is \code{"10"}
#' @param STAR.Alignment.alignIntronMin STAR alignment terminal
#' '--alignIntronMin' option. Default value is \code{"21"}
#' @param STAR.Alignment.alignIntronMax STAR alignment terminal
#' '--alignIntronMax' option. Default value is \code{"0"}
#' @param STAR.Alignment.alignMatesGapMax STAR alignment terminal
#' '--alignMatesGapMax' option. Default value is \code{"0"}
#' @param STAR.Alignment.alignSJoverhangMin STAR alignment terminal
#' '--alignSJoverhangMin' option. Default value is \code{"5"}
#' @param STAR.Alignment.alignSJDBoverhangMin STAR alignment terminal
#' '--alignSJDBoverhangMin' option. Default value is \code{"3"}
#' @param STAR.Alignment.alignSplicedMateMapLmin STAR alignment terminal
#' '--alignSplicedMateMapLmin' option. Default value is \code{"0"}
#' @param STAR.Alignment.alignSplicedMateMapLminOverLmate STAR alignment terminal
#' '--alignSplicedMateMapLminOverLmate' option. Default value is \code{"0.66"}
#' @param STAR.Alignment.alignWindowsPerReadNmax STAR alignment terminal
#' '--alignWindowsPerReadNmax' option. Default value is \code{"10000"}
#' @param STAR.Alignment.alignTranscriptsPerWindowNmax STAR alignment terminal
#' '--alignTranscriptsPerWindowNmax' option. Default value is \code{"100"}
#' @param STAR.Alignment.alignTranscriptsPerReadNmax STAR alignment terminal
#' '--alignTranscriptsPerReadNmax' option. Default value is \code{"10000"}
#' @param STAR.Alignment.alignEndsType STAR alignment terminal
#' '--alignEndsType' option. Default value is \code{"Local"}
#' @param STAR.Alignment.winAnchorMultimapNmax STAR alignment terminal
#' '--winAnchorMultimapNmax' option. Default value is \code{"50"}
#' @param STAR.Alignment.winBinNbits STAR alignment terminal
#' '--winBinNbits' option. Default value is \code{"16"}
#' @param STAR.Alignment.winAnchorDistNbins STAR alignment terminal
#' '--winAnchorDistNbins' option. Default value is \code{"9"}
#' @param STAR.Alignment.winFlankNbins STAR alignment terminal
#' '--winFlankNbins' option. Default value is \code{"4"}
#' @param Rsamtools.Bam.run Whether to run 'Rsamtools SAM to BAM' step in this
#'   function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'Rsamtools SAM to BAM' step.
#' @param Samtools.Bam.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Samtools sam to bam step. The default is \code{"1"}
#' @param Rsamtools.nCores The number of cores to use when running
#'  'Rsamtools' step. Default value is \code{1}
#' @param StringTie.Assemble.run Whether to run 'StringTie assembly' step in
#'   this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie assembly' step.
#' @param Stringtie.Assembly.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Stringtie assembly. The default is \code{"1"}
#' @param Stringtie.Assembly.f Stringtie assembly terminal
#' '-f' option. Default value is \code{"0.1"}
#' @param Stringtie.Assembly.m Stringtie assembly terminal
#' '-m' option. Default value is \code{"200"}
#' @param Stringtie.Assembly.c Stringtie assembly terminal
#' '-c' option. Default value is \code{"2.5"}
#' @param Stringtie.Assembly.g Stringtie assembly terminal
#' '-g' option. Default value is \code{"50"}
#' @param Stringtie.Assembly.M Stringtie assembly terminal
#' '-M' option. Default value is \code{"0.95"}
#' @param StringTie.Merge.Trans.run Whether to run 'StringTie GTF merging' step
#'   in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie GTF merging' step.
#' @param Stringtie.Merge.num.parallel.threads  Specify the number of processing
#' threads (CPUs) to use for Stringtie merge step. The default is \code{"1"}
#' @param Gffcompare.Ref.Sample.run Whether to run 'Gffcompare comparison' step
#'   in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'Gffcompare comparison' step.
#' @param StringTie.Ballgown.run Whether to run 'StringTie ballgown creation'
#'   step in this function step. Default value is \code{TRUE}.
#'   Set \code{FALSE} to skip 'StringTie ballgown creation' step.
#' @param Stringtie.2.Ballgown.num.parallel.threads Specify the number of processing
#' threads (CPUs) to use for Stringtie to ballgown step. The default is \code{"1"}
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
# Parameters: 120
RNASeqReadProcess <- function(RNASeqRParam,
                              which.trigger             = "OUTSIDE",
                              INSIDE.path.prefix        = NA,
                              SAMtools.or.Rsamtools     = "Rsamtools",
                              Hisat2.Index.run          = TRUE,
                              Hisat2.Index.num.parallel.threads = "1",
                              Hisat2.Index.large.index = FALSE,
                              Hisat2.Index.local.ftab.chars = "6",
                              Hisat2.Index.local.off.rate = "3",
                              Hisat2.Index.ftab.chars = "10",
                              Hisat2.Index.off.rate = "4",
                              Hisat2.Alignment.run      = TRUE,
                              Hisat2.Alignment.num.parallel.threads = "1",
                              Hisat2.Alignment.skip = "0",
                              Hisat2.Alignment.trim5 = "0",
                              Hisat2.Alignment.trim3 = "0",
                              Hisat2.Alignment.n.ceil.1.function.type = "L",
                              Hisat2.Alignment.n.ceil.2.constant.term = "0",
                              Hisat2.Alignment.n.ceil.3.coefficient = "0.15",
                              Hisat2.Alignment.mp.MX = "6",
                              Hisat2.Alignment.mp.MN = "2",
                              Hisat2.Alignment.sp.MX = "2",
                              Hisat2.Alignment.sp.MN = "1",
                              Hisat2.Alignment.np = "1",
                              Hisat2.Alignment.rdg.1 = "5",
                              Hisat2.Alignment.rdg.2 = "3",
                              Hisat2.Alignment.rfg.1 = "5",
                              Hisat2.Alignment.rfg.2 = "3",
                              Hisat2.Alignment.score.min.1.function.type = "L",
                              Hisat2.Alignment.score.min.2.constant.term = "0",
                              Hisat2.Alignment.score.min.3.coefficient = "-0.2",
                              Hisat2.Alignment.pen.cansplice = "0",
                              Hisat2.Alignment.penc.noncansplice = "12",
                              Hisat2.Alignment.pen.canintronlen.1.function.type = "G",
                              Hisat2.Alignment.pen.canintronlen.2.constant.term = "-8",
                              Hisat2.Alignment.pen.canintronlen.3.coefficient = "1",
                              Hisat2.Alignment.pen.noncanintronlen.1.function.type = "G",
                              Hisat2.Alignment.pen.noncanintronlen.2.constant.term = "-8",
                              Hisat2.Alignment.pen.noncanintronlen.3.coefficient = "1",
                              Hisat2.Alignment.min.intronlen = "20",
                              Hisat2.Alignment.max.intronlen = "500000",
                              Hisat2.Alignment.rna.strandness = "None",
                              Hisat2.Alignment.k = "5",
                              Hisat2.Alignment.max.seeds = "5",
                              Hisat2.Alignment.secondary = FALSE,
                              Hisat2.Alignment.minins = "0",
                              Hisat2.Alignment.maxins = "500",
                              Hisat2.Alignment.seed = "0",
                              STAR.Index.num.parallel.threads = "1",
                              STAR.Index.sjdbOverhang.Read.length = "100",
                              STAR.Index.genomeSAindexNbases = "14",
                              STAR.Index.genomeChrBinNbits = "18",
                              STAR.Index.genomeSAsparseD = "1",
                              STAR.Alignment.run        = FALSE,
                              STAR.Alignment.num.parallel.threads = "1",
                              STAR.Alignment.genomeLoad = "NoSharedMemory",
                              STAR.Alignment.readMapNumber = "-1",
                              STAR.Alignment.clip3pNbases = "0",
                              STAR.Alignment.clip5pNbases = "0",
                              STAR.Alignment.clip3pAdapterSeq = "-",
                              STAR.Alignment.clip3pAdapterMMp = "0.1",
                              STAR.Alignment.clip3pAfterAdapterNbases = "0",
                              STAR.Alignment.limitGenomeGenerateRAM = "31000000000",
                              STAR.Alignment.limitIObufferSize = "150000000",
                              STAR.Alignment.limitOutSAMoneReadBytes = "100000",
                              STAR.Alignment.limitOutSJoneRead = "1000",
                              STAR.Alignment.limitOutSJcollapsed = "1000000",
                              STAR.Alignment.limitBAMsortRAM = "0",
                              STAR.Alignment.outReadsUnmapped = "None",
                              STAR.Alignment.outQSconversionAdd = "0",
                              STAR.Alignment.outSAMprimaryFlag = "OneBestScore",
                              STAR.Alignment.outSAMmapqUnique = "255",
                              STAR.Alignment.scoreGap = "0",
                              STAR.Alignment.scoreGapNoncan = "-8",
                              STAR.Alignment.scoreGapGCAG = "-4",
                              STAR.Alignment.scoreGapATAC = "-8",
                              STAR.Alignment.scoreGenomicLengthLog2scale = "-0.25",
                              STAR.Alignment.scoreDelOpen = "-2",
                              STAR.Alignment.scoreDelBase = "-2",
                              STAR.Alignment.scoreInsOpen = "-2",
                              STAR.Alignment.scoreInsBase = "-2",
                              STAR.Alignment.scoreStitchSJshift = "1",
                              STAR.Alignment.seedSearchStartLmax = "50",
                              STAR.Alignment.seedSearchStartLmaxOverLread = "1.0",
                              STAR.Alignment.seedSearchLmax = "0",
                              STAR.Alignment.seedMultimapNmax = "10000",
                              STAR.Alignment.seedPerReadNmax = "1000",
                              STAR.Alignment.seedPerWindowNmax = "50",
                              STAR.Alignment.seedNoneLociPerWindow = "10",
                              STAR.Alignment.alignIntronMin = "21",
                              STAR.Alignment.alignIntronMax = "0",
                              STAR.Alignment.alignMatesGapMax = "0",
                              STAR.Alignment.alignSJoverhangMin = "5",
                              STAR.Alignment.alignSJDBoverhangMin = "3",
                              STAR.Alignment.alignSplicedMateMapLmin = "0",
                              STAR.Alignment.alignSplicedMateMapLminOverLmate = "0.66",
                              STAR.Alignment.alignWindowsPerReadNmax = "10000",
                              STAR.Alignment.alignTranscriptsPerWindowNmax = "100",
                              STAR.Alignment.alignTranscriptsPerReadNmax = "10000",
                              STAR.Alignment.alignEndsType = "Local",
                              STAR.Alignment.winAnchorMultimapNmax = "50",
                              STAR.Alignment.winBinNbits = "16",
                              STAR.Alignment.winAnchorDistNbins = "9",
                              STAR.Alignment.winFlankNbins = "4",
                              Rsamtools.Bam.run         = TRUE,
                              Samtools.Bam.num.parallel.threads = "1",
                              Rsamtools.nCores          = "1",
                              StringTie.Assemble.run    = TRUE,
                              Stringtie.Assembly.num.parallel.threads = "1",
                              Stringtie.Assembly.f = "0.1",
                              Stringtie.Assembly.m = "200",
                              Stringtie.Assembly.c = "2.5",
                              Stringtie.Assembly.g = "50",
                              Stringtie.Assembly.M = "0.95",
                              StringTie.Merge.Trans.run = TRUE,
                              Stringtie.Merge.num.parallel.threads = "1",
                              Gffcompare.Ref.Sample.run = TRUE,
                              StringTie.Ballgown.run    = TRUE,
                              Stringtie.2.Ballgown.num.parallel.threads = "1",
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
  } else if (RNASeqRParam == "INSIDE" &
             which.trigger == "INSIDE" &
             !is.na(INSIDE.path.prefix)) {
    # This is an internal call!!
    # Load the S4 object that saved in CMD process
    RNASeqRParam <- readRDS(paste0(INSIDE.path.prefix,
                                   "gene_data/RNASeqRParam.rds"))
  }

  which.s4.object <- CheckS4Object_All(RNASeqRParam, check.s4.print)
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
  fastq.gz.type <- "@"(RNASeqRParam, fastq.gz.type)
  ExportPath(path.prefix)
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)

  if (which.s4.object == "RNASeqRParam") {
    indices.optional <- "@"(RNASeqRParam, indices.optional)
    PreRNASeqReadProcess(path.prefix, genome.name, sample.pattern)
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
      # Parameters: 11
      CreateHisat2Index(path.prefix,
                        genome.name,
                        sample.pattern,
                        splice.site.info = TRUE,
                        exon.info = TRUE,
                        Hisat2.Index.num.parallel.threads,
                        Hisat2.Index.large.index,
                        Hisat2.Index.local.ftab.chars,
                        Hisat2.Index.local.off.rate,
                        Hisat2.Index.ftab.chars,
                        Hisat2.Index.off.rate)
    }

    if (Hisat2.Alignment.run) {
      # Parameters: 43
      Hisat2AlignmentDefault(path.prefix,
                             fastq.gz.type,
                             genome.name,
                             sample.pattern,
                             independent.variable,
                             case.group,
                             control.group,
                             Hisat2.Alignment.num.parallel.threads,
                             Hisat2.Alignment.skip,
                             Hisat2.Alignment.trim5,
                             Hisat2.Alignment.trim3,
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
      # Parameters: 8
      CreateSTARIndex(path.prefix,
                      genome.name,
                      sample.pattern,
                      STAR.Index.num.parallel.threads,
                      STAR.Index.sjdbOverhang.Read.length,
                      STAR.Index.genomeSAindexNbases,
                      STAR.Index.genomeChrBinNbits,
                      STAR.Index.genomeSAsparseD)
      # Parameters: 53
      STARAlignmentDefault(path.prefix,
                           fastq.gz.type,
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
      # Parameters: 6
      RSamtoolsToBam(SAMtools.or.Rsamtools,
                     Samtools.Bam.num.parallel.threads,
                     path.prefix,
                     genome.name,
                     sample.pattern,
                     Rsamtools.nCores)
    }
  } else if (which.s4.object == "RNASeqRParam_Sam") {
    PreRNASeqReadProcess_Sam(path.prefix, genome.name, sample.pattern)
    if (Rsamtools.Bam.run) {
      # Parameters: 6
      RSamtoolsToBam(SAMtools.or.Rsamtools,
                     Samtools.Bam.num.parallel.threads,
                     path.prefix,
                     genome.name,
                     sample.pattern,
                     Rsamtools.nCores)
    }
  } else if (which.s4.object == "RNASeqRParam_Bam") {
    PreRNASeqReadProcess_Bam(path.prefix, genome.name, sample.pattern)
  }
  if (StringTie.Assemble.run) {
    # Parameters: 9
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
    # Parameters: 4
    StringTieMergeTrans(path.prefix,
                        genome.name,
                        sample.pattern,
                        Stringtie.Merge.num.parallel.threads)
  }
  if (Gffcompare.Ref.Sample.run) {
    # Parameters: 3
    GffcompareRefSample(path.prefix,
                        genome.name,
                        sample.pattern)
  }
  if (StringTie.Ballgown.run) {
    # Parameters: 4
    StringTieToBallgown(path.prefix,
                        genome.name,
                        sample.pattern,
                        Stringtie.2.Ballgown.num.parallel.threads)
  }
  finals <- ProgressGenesFiles(path.prefix,
                               genome.name,
                               sample.pattern,
                               print=TRUE)
  if (PreDECountTable.run &
      ((python.variable.answer & python.variable.version == 2) |
       python.variable.answer & python.variable.version == 3 & python.2to3)) {
    # Parameters: 6
    PreDECountTable(path.prefix,
                    sample.pattern,
                    python.variable.answer,
                    python.variable.version,
                    python.2to3,
                    print=TRUE)
  }
  if (which.s4.object == "RNASeqRParam") {
    PostRNASeqReadProcess(path.prefix,
                          genome.name,
                          sample.pattern)
  } else if (which.s4.object == "RNASeqRParam_Sam") {
    PostRNASeqReadProcess_Sam(path.prefix,
                          genome.name,
                          sample.pattern)
  } else if (which.s4.object == "RNASeqRParam_Bam") {
    PostRNASeqReadProcess_Bam(path.prefix,
                              genome.name,
                              sample.pattern)
  }
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

PreRNASeqReadProcess_Sam <- function(path.prefix, genome.name, sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqReadProcess()' ",
          "environment pre-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  ref.gtf <- file.exists(paste0(path.prefix,
                                "gene_data/ref_genes/", genome.name, ".gtf"))
  raw.sam <- list.files(path = paste0(path.prefix, 'gene_data/raw_sam/'),
                          pattern = sample.pattern,
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
  check.tool.result <- CheckTool_Sam_Bam(path.prefix)
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)
  check.progress.results.bool <- check.results$gtf.file.logic.df &&
    (check.results$sam.files.number.df != 0)
  validity <- phenodata.csv && ref.gtf && check.tool.result &&
    (length(raw.sam) != 0) && check.progress.results.bool
  if (!isTRUE(validity)) {
    stop("RNASeqReadProcess() environment ERROR")
  }
  message("(\u2714) : RNASeqReadProcess() pre-check is valid\n\n")
}

PreRNASeqReadProcess_Bam <- function(path.prefix, genome.name, sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqReadProcess()' ",
          "environment pre-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  ref.gtf <- file.exists(paste0(path.prefix,
                                "gene_data/ref_genes/", genome.name, ".gtf"))
  raw.bam <- list.files(path = paste0(path.prefix, 'gene_data/raw_bam/'),
                        pattern = sample.pattern,
                        all.files = FALSE,
                        full.names = FALSE,
                        recursive = FALSE,
                        ignore.case = FALSE)
  check.tool.result <- CheckTool_Sam_Bam(path.prefix)
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)
  check.progress.results.bool <- check.results$gtf.file.logic.df &&
    (check.results$bam.files.number.df != 0)
  validity <- phenodata.csv && ref.gtf && check.tool.result &&
    (length(raw.bam) != 0) && check.progress.results.bool
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

PostRNASeqReadProcess_Sam <- function(path.prefix, genome.name, sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqReadProcess()' ",
          "environment post-check ...\n")
  # Still need to add condition
  gene_abundance <- dir.exists(paste0(path.prefix, "gene_data/gene_abundance/"))
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)
  sam.bool <- (check.results$sam.files.number.df) != 0
  bam.bool <- (check.results$bam.files.number.df) != 0
  gtf.bool <- (check.results$gtf.files.number.df) != 0
  merged.bool <- check.results$stringtie_merged.gtf.file.df
  gffcompare.bool <- (check.results$gffcompare.related.dirs.number.df) != 0
  ballgown.bool <- (check.results$ballgown.dirs.number.df) != 0
  validity <- gene_abundance && sam.bool && bam.bool && gtf.bool &&
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

PostRNASeqReadProcess_Bam <- function(path.prefix, genome.name, sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqReadProcess()' ",
          "environment post-check ...\n")
  # Still need to add condition
  gene_abundance <- dir.exists(paste0(path.prefix, "gene_data/gene_abundance/"))
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)
  bam.bool <- (check.results$bam.files.number.df) != 0
  gtf.bool <- (check.results$gtf.files.number.df) != 0
  merged.bool <- check.results$stringtie_merged.gtf.file.df
  gffcompare.bool <- (check.results$gffcompare.related.dirs.number.df) != 0
  ballgown.bool <- (check.results$ballgown.dirs.number.df) != 0
  validity <- gene_abundance && bam.bool && gtf.bool &&
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
