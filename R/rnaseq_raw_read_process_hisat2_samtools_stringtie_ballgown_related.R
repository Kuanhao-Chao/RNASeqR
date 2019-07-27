# Creating Hisat2 index
# Creat Hisat2 index for further use
CreateHisat2Index <- function (path.prefix,
                               genome.name,
                               sample.pattern,
                               splice.site.info = TRUE,
                               exon.info = TRUE,
                               Hisat2.Index.num.parallel.threads,
                               Hisat2.large.index,
                               Hisat2.local.ftab.chars,
                               Hisat2.local.off.rate,
                               Hisat2.ftab.chars,
                               Hisat2.off.rate) {
  if (isTRUE(CheckHisat2(print=FALSE))){
    if (!is.logical(splice.site.info) || !is.logical(exon.info)) {
      stop("(\u2718) Please make sure the type of ",
           "'splice.site.info' and 'exon.info' are logical.\n")
    } else {
      # Check file progress
      check.results <- ProgressGenesFiles(path.prefix,
                                          genome.name,
                                          sample.pattern,
                                          print=TRUE)
      message("\n\u2618\u2618\u2618 Index Creation :\n")
      message("************** Creating Hisat2 Index **************\n")
      if (check.results$gtf.file.logic.df && check.results$fa.file.logic.df){
        command.list <- c()
        command.list <- c(command.list, "* Creating Hisat2 Index : ")
        current.path <- getwd()
        setwd(paste0(path.prefix, "gene_data/indices/"))
        if (isTRUE(splice.site.info)) {
          whole.command <-  paste0(path.prefix, "gene_data/ref_genes/",
                                   genome.name, ".gtf > ", genome.name, ".ss")
          main.command <- "extract_splice_sites.py"
          message("Input command : ",
                  paste(main.command, whole.command), "\n")
          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command,
                                    args = whole.command)
          # break is return is not 0 (means success!)
          if (command.result != 0 ) {
            on.exit(setwd(current.path))
            message("(\u2718) '", main.command, "' is failed !!")
            stop("'", main.command, "' ERROR")
          }
          message("\n")
        }
        if (isTRUE(exon.info)) {
          whole.command <- paste0(path.prefix, 'gene_data/ref_genes/',
                                  genome.name, '.gtf > ', genome.name, '.exon')
          main.command <- "extract_exons.py"
          message("Input command : ",
                  paste(main.command, whole.command), "\n")
          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command,
                                    args = whole.command)
          if (command.result != 0 ) {
            on.exit(setwd(current.path))
            message("(\u2718) '", main.command, "' is failed !!")
            stop("'", main.command, "' ERROR")
          }
          message("\n")
        }
        if (isTRUE(splice.site.info) && isTRUE(exon.info)) {
          if (isTRUE(Hisat2.large.index)) {
            Hisat2.large.index.value <- "--large-index"
          } else {
            Hisat2.large.index.value <- ""
          }
          whole.command <- paste(paste(Hisat2.large.index.value,
                                       '-p', Hisat2.Index.num.parallel.threads,
                                       '--localftabchars', Hisat2.local.ftab.chars,
                                       '--localoffrate', Hisat2.local.off.rate,
                                       '--ftabchars', Hisat2.ftab.chars,
                                       '--offrate', Hisat2.off.rate,
                                       '--ss', paste0(genome.name, '.ss'),
                                       "--exon", paste0(genome.name, '.exon'),
                                       paste0(path.prefix,
                                              'gene_data/ref_genome/',
                                              genome.name, '.fa'),
                                       paste0(genome.name, '_tran')))
          main.command <- "hisat2-build"
          message("Input command : ", paste(main.command, whole.command), "\n")

          # --large-index
          #num.parallel.threads.hisat2.index  -p : 1
          # --localftabchars : 6
          # --localoffrate : 3
          # --ftabchars : 10
          # --offrate : 4

          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command,
                                    args = whole.command)
          if (command.result != 0 ) {
            on.exit(setwd(current.path))
            message("(\u2718) '", main.command, "' is failed !!")
            stop("'", main.command, "' ERROR")
          }
          message("\n")
        } else if (isTRUE(splice.site.info) && !isTRUE(exon.info)) {
          whole.command <- paste(paste('--ss', paste0(genome.name, '.ss'),
                                       paste0(path.prefix,
                                              'gene_data/ref_genome/',
                                              genome.name, '.fa'),
                                       paste0(genome.name, '_tran')))
          main.command <- "hisat2-build"
          message("Input command : ", paste(main.command, whole.command), "\n")
          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command,
                                    args = whole.command)
          if (command.result != 0 ) {
            on.exit(setwd(current.path))
            message("(\u2718) '", main.command, "' is failed !!")
            stop("'", main.command, "' ERROR")
          }
          message("\n")
        } else if (!splice.site.info && exon.info) {
          whole.command <- paste(paste('--exon', paste0(genome.name, '.exon'),
                                       paste0(path.prefix,
                                              'gene_data/ref_genome/',
                                              genome.name, '.fa'),
                                       paste0(genome.name, '_tran')))
          main.command <- "hisat2-build"
          message("Input command : ", paste(main.command, whole.command), "\n")
          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command,
                                    args = whole.command)
          if (command.result != 0 ) {
            on.exit(setwd(current.path))
            message("(\u2718) '", main.command, "' is failed !!")
            stop("'", main.command, "' ERROR")
          }
          message("\n")
        } else {
          whole.command <- paste(paste(paste0(path.prefix,
                                              'gene_data/ref_genome/',
                                              genome.name, '.fa'),
                                       paste0(genome.name, '_tran')))
          main.command <- "hisat2-build"
          message("Input command : ", paste(main.command, whole.command), "\n")
          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command,
                                    args = whole.command)
          if (command.result != 0 ) {
            on.exit(setwd(current.path))
            message("(\u2718) '", main.command, "' is failed !!")
            stop("'", main.command, "' ERROR")
          }
          message("\n")
        }
        command.list <- c(command.list, "\n")
        fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
        write(command.list, fileConn, append = TRUE)
        on.exit(setwd(current.path))
        message("'", path.prefix, "gene_data/indices/",
                genome.name, "_tran.*.ht2' has been created.\n\n")
      } else {
        on.exit(setwd(current.path))
        stop("(\u2718) '", genome.name, ".gtf' or",
             " '", genome.name, ".fa'is missing.\n\n")
      }
    }
  }
}

# hisat2 alignment default
Hisat2AlignmentDefault <- function(path.prefix,
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
                                   Hisat2.Alignment.seed) {
  if (isTRUE(CheckHisat2(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Alignment :\n")
    message("************** Hisat2 Alignment **************\n")
    if (check.results$ht2.files.number.df != 0 &&
        check.results$fastq.gz.files.number.df != 0){
      command.list <- c()
      command.list <- c(command.list, "* Hisat2 Alignment : ")
      # Map reads to each alignment
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      deleteback <- gsub("[1-2]*.fastq.gz$", replacement = "",
                         check.results$fastq.gz.files.df)
      sample.table.r.value <- gsub(paste0("[A-Z, a-z]*[0-9]*_"),
                                   replacement = "",
                                   deleteback)
      if (isTRUE(length(unique(sample.table.r.value)) != 1)){
        on.exit(setwd(current.path))
        stop("(\u2718) Inconsistent formats. Please check files are all",
             "'XXX_*.fastq.gz'\n\n")
      } else {
        sample.table <- table(gsub(paste0("_[1-2]*.fastq.gz$"),
                                   replacement = "",
                                   check.results$fastq.gz.files.df))
        iteration.num <- length(sample.table)
        sample.name <- names(sample.table)
        sample.value <- as.vector(sample.table)
        alignment.result <- data.frame(matrix(0, ncol = 0, nrow = 5))
        total.map.rates <- c()
        for( i in seq_len(iteration.num)){
          current.sub.command <- ""
          total.sub.command <- ""
          for ( j in seq_len(sample.value[i])){
            current.sub.command <- paste(paste0("-", j),
                                         paste0("raw_fastq.gz/",
                                                sample.name[i], "_",
                                                sample.table.r.value[1],
                                                j, ".fastq.gz"))
            total.sub.command <- paste(total.sub.command, current.sub.command)
          }
          if (Hisat2.Alignment.qupto == "None") {
            qupto.value = ""
          } else if (isTRUE(strtoi(Hisat2.Alignment.qupto)%%1==0)) {
            # 'qupto' must be integer !
            Hisat2.Alignment.qupto.value <- paste("--qseq", Hisat2.Alignment.qupto)
          } else {
            stop("'qupto' must be 'None' or integer")
          }
          if (isTRUE(Hisat2.Alignment.int.quals)) {
            Hisat2.Alignment.int.quals.value = ""
          } else {
            Hisat2.Alignment.int.quals.value = "--int-quals"
          }
          if (Hisat2.Alignment.rna.strandness == "FR") {
            Hisat2.Alignment.rna.strandness.value = "--rna-strandness FR"
          } else if (Hisat2.Alignment.rna.strandness == "RF") {
            Hisat2.Alignment.rna.strandness.value = "--rna-strandness RF"
          } else if (Hisat2.Alignment.rna.strandness == "None") {
            Hisat2.Alignment.rna.strandness.value = ""
          } else {
            stop("'rna.strandness' variable is out of range. It should be 'FR' or 'RF', 'None'.")
          }
          if (isTRUE(Hisat2.Alignment.secondary)) {
            Hisat2.Alignment.secondary.value = "--secondary"
          } else {
            Hisat2.Alignment.secondary.value = ""
          }
          whole.command <- paste("--skip", Hisat2.Alignment.skip, Hisat2.Alignment.qupto.value,
                                 "--trim5", Hisat2.Alignment.trim5,
                                 "--trim3", Hisat2.Alignment.trim3,
                                 paste0("--phred", Hisat2.Alignment.phred),
                                 Hisat2.Alignment.int.quals.value,
                                 "--n-ceil", Hisat2.Alignment.n.ceil.1.function.type,
                                 Hisat2.Alignment.n.ceil.2.constant.term,
                                 Hisat2.Alignment.n.ceil.3.coefficient,
                                 "--mp", Hisat2.Alignment.mp.MX, Hisat2.Alignment.mp.MN,
                                 "--sp", Hisat2.Alignment.sp.MX, Hisat2.Alignment.sp.MN,
                                 "--np", Hisat2.Alignment.np,
                                 "--rdg", Hisat2.Alignment.rdg.1, Hisat2.Alignment.rdg.2,
                                 "--rfg", Hisat2.Alignment.rfg.1, Hisat2.Alignment.rfg.2,
                                 "--score-min" ,Hisat2.Alignment.score.min.1.function.type,
                                 Hisat2.Alignment.score.min.2.constant.term,
                                 Hisat2.Alignment.score.min.3.coefficient,
                                 "--pen-cansplice", Hisat2.Alignment.pen.cansplice,
                                 "--pen-noncansplice", Hisat2.Alignment.penc.noncansplice,
                                 "--pen-canintronlen",
                                 Hisat2.Alignment.pen.canintronlen.1.function.type,
                                 Hisat2.Alignment.pen.canintronlen.2.constant.term,
                                 Hisat2.Alignment.pen.canintronlen.3.coefficient,
                                 "--pen-noncanintronlen",
                                 Hisat2.Alignment.pen.noncanintronlen.1.function.type,
                                 Hisat2.Alignment.pen.noncanintronlen.2.constant.term,
                                 Hisat2.Alignment.pen.noncanintronlen.3.coefficient,
                                 "--min-intronlen", Hisat2.Alignment.min.intronlen,
                                 "--max-intronlen", Hisat2.Alignment.max.intronlen,
                                 Hisat2.Alignment.rna.strandness.value,
                                 "-k", Hisat2.Alignment.k,
                                 "--max-seeds", Hisat2.Alignment.max.seeds,
                                 Hisat2.Alignment.secondary.value,
                                 "--minins", Hisat2.Alignment.minins,
                                 "--maxins", Hisat2.Alignment.maxins,
                                 "--seed", Hisat2.Alignment.seed, "--time",
                                 "-p", Hisat2.Alignment.num.parallel.threads,"--dta -x",
                                 paste0("indices/", genome.name, "_tran"),
                                 total.sub.command, "-S",
                                 paste0("raw_sam/", sample.name[i],".sam") )
          if (i != 1) message("\n")
          main.command <- "hisat2"
          message("Input command : ", paste(main.command, whole.command), "\n")
          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command,
                                    args = whole.command,
                                    stderr = TRUE, stdout = TRUE)
          # Total reads
          total.reads <- gsub(" reads; of these:", "", command.result[1])
          # aligned concordantly exactly 1 time
          concordantly_1_time <- as.numeric(gsub(" \\([0-9]*.[0-9]*%) aligned concordantly exactly 1 time", "", command.result[4]))
          # aligned concordantly >1 times
          concordantly_more_1_times <- as.numeric(gsub(" \\([0-9]*.[0-9]*%) aligned concordantly >1 times", "", command.result[5]))
          # aligned dicordantly 1 time
          dicordantly_1_time <- as.numeric(gsub(" \\([0-9]*.[0-9]*%) aligned discordantly 1 time", "", command.result[8]))
          # aligned 0 times concordantly or discordantly
          not_dicordantly_concordantly <- as.numeric(gsub(" pairs aligned 0 times concordantly or discordantly; of these:", "", command.result[10]))
          # Total mapping rate
          total.map.rate <- as.numeric(gsub("% overall alignment rate", "", command.result[15]))
          total.map.rates <- c(total.map.rates, total.map.rate)
          one.result <- c(total.reads, concordantly_1_time, concordantly_more_1_times, dicordantly_1_time, not_dicordantly_concordantly)
          alignment.result[[sample.name[i]]] <- one.result
          if (length(command.result) == 0) {
            on.exit(setwd(current.path))
            message("(\u2718) '", main.command, "' is failed !!")
            stop("'", main.command, "' ERROR")
          }
        }
        row.names(alignment.result) <- c("total_reads", "concordantly_1", "concordantly_more_1", "dicordantly_1", "not_both")
        alignment.result[] <- lapply(alignment.result, as.character)
        alignment.result[] <- lapply(alignment.result, as.numeric)
        if(!dir.exists(paste0(path.prefix, "RNASeq_results/Alignment_Report/"))){
          dir.create(paste0(path.prefix, "RNASeq_results/Alignment_Report/"))
        }
        trans.df  <- data.frame(t(alignment.result))
        write.csv(trans.df,
                  file = paste0(path.prefix,
                                "RNASeq_results/",
                                "Alignment_Report/Alignment_report_reads.csv"),
                  row.names = TRUE)
        trans.df.portion  <- data.frame(t(alignment.result))
        trans.df.portion$concordantly_1 <- round(trans.df.portion$concordantly_1 / trans.df.portion$total_reads, 4)
        trans.df.portion$concordantly_more_1 <- round(trans.df.portion$concordantly_more_1 / trans.df.portion$total_reads, 4)
        trans.df.portion$dicordantly_1 <- round(trans.df.portion$dicordantly_1 / trans.df.portion$total_reads, 4)
        trans.df.portion$not_both <- round(trans.df.portion$not_both / trans.df.portion$total_reads, 4)
        write.csv(trans.df.portion,
                  file = paste0(path.prefix,
                                "RNASeq_results/",
                                "Alignment_Report/Alignment_report_proportion.csv"),
                  row.names = TRUE)
        names(total.map.rates) <- sample.name
        total.map.rates <- data.frame(total.map.rates)
        write.csv(total.map.rates,
                  file = paste0(path.prefix,
                                "RNASeq_results/",
                                "Alignment_Report/Overall_Mapping_rate.csv"),
                  row.names = TRUE)
        message("\n")
        command.list <- c(command.list, "\n")
        fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
        write(command.list, fileConn, append = TRUE)
        on.exit(setwd(current.path))
      }
    } else {
      on.exit(setwd(current.path))
      stop("(\u2718) '", genome.name, "_tran.*.ht2' ",
           "or 'XXX_*.fastq.gz' is missing.\n\n")
    }
  }
}


# Creating STAR index
# Creat STAR index for further use
CreateSTARIndex <- function (path.prefix,
                             genome.name,
                             sample.pattern,
                             STAR.Index.num.parallel.threads = 1,
                             STAR.Index.sjdbOverhang.Read.length = 100,
                             STAR.Index.genomeSAindexNbases = 14,
                             STAR.Index.genomeChrBinNbits = 18,
                             STAR.Index.genomeSAsparseD = 1) {
  if (isTRUE(CheckSTAR(print=FALSE))){
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Index Creation :\n")
    message("************** Creating STAR Index **************\n")
    if (check.results$gtf.file.logic.df && check.results$fa.file.logic.df){
      command.list <- c()
      command.list <- c(command.list, "* Creating STAR Index : ")
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/indices/"))
      whole.command <- paste("--genomeSAindexNbases", genomeSAindexNbases,
                             "--genomeChrBinNbits", genomeChrBinNbits,
                             "--genomeSAsparseD", genomeSAsparseD,
                             "--runThreadN", num.parallel.threads.star.index,
                             "--runMode", "genomeGenerate",
                             "--genomeDir",
                             paste0(path.prefix, "gene_data/indices/"),
                             "--genomeFastaFiles",
                             paste0(path.prefix, 'gene_data/ref_genome/',
                                    genome.name, '.fa'),
                             "--sjdbGTFfile",
                             paste0(path.prefix, "gene_data/ref_genes/",
                                    genome.name, ".gtf"),
                             "--sjdbOverhang", sjdbOverhang.Read.length-1)
      main.command <- "STAR"
      message("Input command : ", paste(main.command, whole.command), "\n")
      command.list <- c(command.list,
                        paste("    command :", main.command, whole.command))
      command.result <- system2(command = main.command,
                                args = whole.command)
      if (command.result != 0 ) {
        on.exit(setwd(current.path))
        message("(\u2718) '", main.command, "' is failed !!")
        stop("'", main.command, "' ERROR")
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
      message("'", path.prefix, "gene_data/indices/",
              genome.name, "STAR indices has been created.\n\n")
    } else {
      on.exit(setwd(current.path))
      stop("(\u2718) '", genome.name, ".gtf' or",
           " '", genome.name, ".fa'is missing.\n\n")
    }
  }
}

# STAR alignment default
STARAlignmentDefault <- function(path.prefix,
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
                                 STAR.Alignment.winFlankNbins) {
  if (isTRUE(CheckSTAR(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Alignment :\n")
    message("************** STAR Alignment **************\n")
    samples.star.dir <- dir.create(file.path(paste0(path.prefix,
                                                    "gene_data/raw_star")),
                                   showWarnings = FALSE) == 0
    if (!isTRUE(samples.star.dir)) {
      message("     (\u2714) : Create '", path.prefix,
              "gene_data/raw_star/'.\n")
    } else {
      message("     (\u26A0) : '", path.prefix,
              "gene_data/raw_star/' has already be created.\n")
    }
    if (check.results$ht2.files.number.df != 0 &&
        check.results$fastq.gz.files.number.df != 0){
      command.list <- c()
      command.list <- c(command.list, "* STAR Alignment : ")
      # Map reads to each alignment
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      deleteback <- gsub("[1-2]*.fastq.gz$", replacement = "",
                         check.results$fastq.gz.files.df)
      sample.table.r.value <- gsub(paste0("[A-Z, a-z]*[0-9]*_"),
                                   replacement = "",
                                   deleteback)
      if (isTRUE(length(unique(sample.table.r.value)) != 1)){
        on.exit(setwd(current.path))
        stop("(\u2718) Inconsistent formats. Please check files are all",
             "'XXX_*.fastq.gz'\n\n")
      } else {
        sample.table <- table(gsub(paste0("_[1-2]*.fastq.gz$"),
                                   replacement = "",
                                   check.results$fastq.gz.files.df))
        iteration.num <- length(sample.table)
        sample.name <- names(sample.table)
        sample.value <- as.vector(sample.table)
        alignment.result <- data.frame(matrix(0, ncol = 0, nrow = 8))
        total.map.rates <- c()
        for( i in seq_len(iteration.num)){
          current.sub.command <- ""
          total.sub.command <- ""
          for ( j in seq_len(sample.value[i])){
            current.sub.command <- paste(paste0("raw_fastq.gz/",
                                                sample.name[i], "_",
                                                sample.table.r.value[1],
                                                j, ".fastq.gz"))
            total.sub.command <- paste(total.sub.command, current.sub.command)
          }

          samples.star.dir.in <- dir.create(file.path(paste0(path.prefix,
                                                             "gene_data/raw_star/",
                                                             sample.name[i])),
                                            showWarnings = FALSE) == 0
          samples.star.dir.output <- paste0(path.prefix, "gene_data/raw_star/",
                                            sample.name[i], "/")
          if (!isTRUE(samples.star.dir.in)) {
            message("     (\u2714) : Create '", samples.star.dir.output, ".\n")
          } else {
            message("     (\u26A0) : '", samples.star.dir.output,
                    " has already be created.\n")
          }
          # --genomeLoad : NoSharedMemory
          # --readMapNumber : -1
          # --clip3pNbases : 0
          # --clip5pNbases : 0
          # --clip3pAdapterSeq : -
          # --clip3pAdapterMMp : 0.1
          # --clip3pAfterAdapterNbases : 0
          # --limitGenomeGenerateRAM : 31000000000
          # --limitIObufferSize : 150000000
          # --limitOutSAMoneReadBytes : 100000
          # --limitOutSJoneRead : 1000
          # --limitOutSJcollapsed : 1000000
          # --limitBAMsortRAM : 0
          # --outReadsUnmapped : None
          # --outQSconversionAdd : 0
          # --outSAMprimaryFlag : OneBestScore
          # --outSAMmapqUnique : 255
          # --scoreGap : 0
          # --scoreGapNoncan : -8
          # --scoreGapGCAG : -4
          # --scoreGapATAC : -8
          # --scoreGenomicLengthLog2scale : -0.25
          # --scoreDelOpen : -2
          # --scoreDelBase : -2
          # --scoreInsOpen : -2
          # --scoreInsBase : -2
          # --scoreStitchSJshift : 1
          # --seedSearchStartLmax : 50
          # --seedSearchStartLmaxOverLread : 1.0
          # --seedSearchLmax : 0
          # --seedMultimapNmax : 10000
          # --seedPerReadNmax : 1000
          # --seedPerWindowNmax : 50
          # --seedNoneLociPerWindow : 10
          # --alignIntronMin : 21
          # --alignIntronMax : 0
          # --alignMatesGapMax : 0
          # --alignSJoverhangMin : 5
          # --alignSJDBoverhangMin : 3
          # --alignSplicedMateMapLmin : 0
          # --alignSplicedMateMapLminOverLmate : 0.66
          # --alignWindowsPerReadNmax : 10000
          # --alignTranscriptsPerWindowNmax : 100
          # --alignTranscriptsPerReadNmax : 10000
          # --alignEndsType : Local
          # --winAnchorMultimapNmax : 50
          # --winBinNbits : 16
          # --winAnchorDistNbins : 9
          # --winFlankNbins : 4

          whole.command <- paste("--runThreadN", STAR.Alignment.num.parallel.threads,
                                 "--runMode", "alignReads",
                                 "--genomeLoad", STAR.Alignment.genomeLoad,
                                 "--readMapNumber", STAR.Alignment.readMapNumber,
                                 "--clip3pNbases", STAR.Alignment.clip3pNbases,
                                 "--clip5pNbases", STAR.Alignment.clip5pNbases,
                                 "--clip3pAdapterSeq", STAR.Alignment.clip3pAdapterSeq,
                                 "--clip3pAdapterMMp", STAR.Alignment.clip3pAdapterMMp,
                                 "--clip3pAfterAdapterNbases", STAR.Alignment.clip3pAfterAdapterNbases,
                                 "--limitGenomeGenerateRAM", STAR.Alignment.limitGenomeGenerateRAM,
                                 "--limitIObufferSize", STAR.Alignment.limitIObufferSize,
                                 "--limitOutSAMoneReadBytes", STAR.Alignment.limitOutSAMoneReadBytes,
                                 "--limitOutSJoneRead", STAR.Alignment.limitOutSJoneRead,
                                 "--limitOutSJcollapsed", STAR.Alignment.limitOutSJcollapsed,
                                 "--limitBAMsortRAM", STAR.Alignment.limitBAMsortRAM,
                                 "--outReadsUnmapped", STAR.Alignment.outReadsUnmapped,
                                 "--outQSconversionAdd", STAR.Alignment.outQSconversionAdd,
                                 "--outSAMprimaryFlag", STAR.Alignment.outSAMprimaryFlag,
                                 "--outSAMmapqUnique", STAR.Alignment.outSAMmapqUnique,
                                 "--scoreGap", STAR.Alignment.scoreGap,
                                 "--scoreGapNoncan", STAR.Alignment.scoreGapNoncan,
                                 "--scoreGapGCAG", STAR.Alignment.scoreGapGCAG,
                                 "--scoreGapATAC", STAR.Alignment.scoreGapATAC,
                                 "--scoreGenomicLengthLog2scale", STAR.Alignment.scoreGenomicLengthLog2scale,
                                 "--scoreDelOpen", STAR.Alignment.scoreDelOpen,
                                 "--scoreDelBase", STAR.Alignment.scoreDelBase,
                                 "--scoreInsOpen", STAR.Alignment.scoreInsOpen,
                                 "--scoreInsBase", STAR.Alignment.scoreInsBase,
                                 "--scoreStitchSJshift", STAR.Alignment.scoreStitchSJshift,
                                 "--seedSearchStartLmax", STAR.Alignment.seedSearchStartLmax,
                                 "--seedSearchStartLmaxOverLread", STAR.Alignment.seedSearchStartLmaxOverLread,
                                 "--seedSearchLmax", STAR.Alignment.seedSearchLmax,
                                 "--seedMultimapNmax", STAR.Alignment.seedMultimapNmax,
                                 "--seedPerReadNmax", STAR.Alignment.seedPerReadNmax,
                                 "--seedPerWindowNmax", STAR.Alignment.seedPerWindowNmax,
                                 "--seedNoneLociPerWindow", STAR.Alignment.seedNoneLociPerWindow,
                                 "--alignIntronMin", STAR.Alignment.alignIntronMin,
                                 "--alignIntronMax", STAR.Alignment.alignIntronMax,
                                 "--alignMatesGapMax", STAR.Alignment.alignMatesGapMax,
                                 "--alignSJoverhangMin", STAR.Alignment.alignSJoverhangMin,
                                 "--alignSJDBoverhangMin", STAR.Alignment.alignSJDBoverhangMin,
                                 "--alignSplicedMateMapLmin", STAR.Alignment.alignSplicedMateMapLmin,
                                 "--alignSplicedMateMapLminOverLmate", STAR.Alignment.alignSplicedMateMapLminOverLmate,
                                 "--alignWindowsPerReadNmax", STAR.Alignment.alignWindowsPerReadNmax,
                                 "--alignTranscriptsPerWindowNmax", STAR.Alignment.alignTranscriptsPerWindowNmax,
                                 "--alignTranscriptsPerReadNmax", STAR.Alignment.alignTranscriptsPerReadNmax,
                                 "--alignEndsType", STAR.Alignment.alignEndsType,
                                 "--winAnchorMultimapNmax", STAR.Alignment.winAnchorMultimapNmax,
                                 "--winBinNbits", STAR.Alignment.winBinNbits,
                                 "--winAnchorDistNbins", STAR.Alignment.winAnchorDistNbins,
                                 "--winFlankNbins", STAR.Alignment.winFlankNbins,
                                 "--genomeDir",
                                 paste0(path.prefix, "gene_data/indices/"),
                                 "--readFilesIn", total.sub.command,
                                 "--outFileNamePrefix", samples.star.dir.output,
                                 "--readFilesCommand", "zcat")
          if (i != 1) message("\n")
          main.command <- "STAR"
          message("Input command : ", paste(main.command, whole.command), "\n")
          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command,
                                    args = whole.command,
                                    stderr = TRUE, stdout = TRUE)
          Log.final.out <-  paste0(path.prefix, "gene_data/raw_star/", sample.name[i], "/Log.final.out")
          Log.final.out.read.in <- read.delim(Log.final.out, header = FALSE, sep = "\t", dec = ".")
          Log.final.out.read.in.num <- suppressWarnings(as.numeric(levels(Log.final.out.read.in["V2"][[1]])[as.integer(Log.final.out.read.in["V2"][[1]])]))
          # Total reads
          total.reads <- Log.final.out.read.in.num[5]
          # Uniquely_mapped
          uniquely_mapping <- Log.final.out.read.in.num[8]
          # mapped to multiple loci
          multi_mapping_multiple_loci <- Log.final.out.read.in.num[23]
          # mapped to too many loci
          multi_mapping_many_loci <- Log.final.out.read.in.num[25]

          # reads unmapped: too many mismatches
          unmapped_many_mismatches <- Log.final.out.read.in.num[28]
          # reads unmapped: too short
          unmapped_short <- Log.final.out.read.in.num[30]
          # reads unmapped: other
          unmapped_other <- Log.final.out.read.in.num[32]
          # number of chimeric reads
          chimeric_reads <- Log.final.out.read.in.num[35]
          # Total mapping rate
          total.map.rate <- round(((uniquely_mapping + multi_mapping_multiple_loci + multi_mapping_many_loci) / total.reads) * 100, digits = 2)
          total.map.rates <- c(total.map.rates, total.map.rate)
          one.result <- c(total.reads, uniquely_mapping,
                          multi_mapping_multiple_loci, multi_mapping_many_loci,
                          unmapped_many_mismatches, unmapped_short,
                          unmapped_other, chimeric_reads)

          alignment.result[[sample.name[i]]] <- one.result
          if (length(command.result) == 0) {
            on.exit(setwd(current.path))
            message("(\u2718) '", main.command, "' is failed !!")
            stop("'", main.command, "' ERROR")
          }
          file.symlink(paste0(path.prefix, "gene_data/raw_star/", sample.name[i], "/Aligned.out.sam"),
                       paste0(path.prefix, "gene_data/raw_sam/", sample.name[i], ".sam"))
        }
        row.names(alignment.result) <- c("total_reads", "uniquely_mapping",
                                         "multi_mapping_multiple_loci",
                                         "multi_mapping_many_loci",
                                         "unmapped_too_many_mismatches",
                                         "unmapped_too_short",
                                         "unmapped_other",
                                         "chimeric_reads")
        alignment.result[] <- lapply(alignment.result, as.character)
        alignment.result[] <- lapply(alignment.result, as.numeric)
        if(!dir.exists(paste0(path.prefix, "RNASeq_results/Alignment_Report/"))){
          dir.create(paste0(path.prefix, "RNASeq_results/Alignment_Report/"))
        }
        trans.df  <- data.frame(t(alignment.result))
        write.csv(trans.df,
                  file = paste0(path.prefix,
                                "RNASeq_results/",
                                "Alignment_Report/Alignment_report_reads.csv"),
                  row.names = TRUE)
        trans.df.portion  <- data.frame(t(alignment.result))
        trans.df.portion$uniquely_mapping <- round(trans.df.portion$uniquely_mapping / trans.df.portion$total_reads, 4)
        trans.df.portion$multi_mapping_multiple_loci <- round(trans.df.portion$multi_mapping_multiple_loci / trans.df.portion$total_reads, 4)
        trans.df.portion$multi_mapping_many_loci <- round(trans.df.portion$multi_mapping_many_loci / trans.df.portion$total_reads, 4)
        trans.df.portion$unmapped_too_many_mismatches <- round(trans.df.portion$unmapped_too_many_mismatches / trans.df.portion$total_reads, 4)
        trans.df.portion$unmapped_too_short <- round(trans.df.portion$unmapped_too_short / trans.df.portion$total_reads, 4)
        trans.df.portion$unmapped_other <- round(trans.df.portion$unmapped_other / trans.df.portion$total_reads, 4)
        trans.df.portion$chimeric_reads <- round(trans.df.portion$chimeric_reads / trans.df.portion$total_reads, 4)
        write.csv(trans.df.portion,
                  file = paste0(path.prefix,
                                "RNASeq_results/",
                                "Alignment_Report/Alignment_report_proportion.csv"),
                  row.names = TRUE)
        names(total.map.rates) <- sample.name

        total.map.rates <- data.frame(total.map.rates)
        write.csv(total.map.rates,
                  file = paste0(path.prefix,
                                "RNASeq_results/",
                                "Alignment_Report/Overall_Mapping_rate.csv"),
                  row.names = TRUE)
        message("\n")
        command.list <- c(command.list, "\n")
        fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
        write(command.list, fileConn, append = TRUE)
        on.exit(setwd(current.path))
      }
    } else {
      on.exit(setwd(current.path))
      stop("(\u2718) 'indices/*' ",
           "or 'XXX_*.fastq.gz' is missing.\n\n")
    }
  }
}

# use 'Rsamtools' to sort and convert the SAM files to BAM
RSamtoolsToBam <- function(SAMtools.or.Rsamtools,
                           Samtools.Bam.num.parallel.threads,
                           path.prefix,
                           genome.name,
                           sample.pattern,
                           Rsamtools.nCores){
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=TRUE)
  message("\n\u2618\u2618\u2618 'SAM' to 'BAM' :\n")
  message("************** ",
          "Rsamtools converting '.sam' to '.bam' **************\n")
  command.list <- c()
  if (check.results$sam.files.number.df != 0){
    if (SAMtools.or.Rsamtools == "SAMtools") {
      check.sam <- CheckSamtools()
      if (check.sam) {
        command.list <- c(command.list,
                          "* SAMtools Converting '.sam' to '.bam' : ")
        sample.table <- table(gsub(paste0(".sam$"),
                                   replacement = "",
                                   check.results$sam.files.df))
        sample.name <- names(sample.table)
        for( i in sample.name){
          whole.command <- paste("sort -@", Samtools.Bam.num.parallel.threads,
                                 "-o", paste0(path.prefix, "gene_data/raw_bam/", i, ".bam"),
                                 paste0(path.prefix, "gene_data/raw_sam/", i, ".sam"))
          if (i != 1) message("\n")
          main.command <- "samtools"
          message("Input command :", paste(main.command, whole.command), "\n")
          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command, args = whole.command)
          if (command.result != 0 ) {
            message("(\u2718) '", main.command, "' is failed !!")
            stop("'", main.command, "' ERROR")
          }
        }
      } else {
        message("(\u2718) : 'SAMtools.or.Rsamtools' parameter can't be set to 'SAMtools'",
                " because 'samtools' command is not found in R shell through ",
                "'system2()'. Please make sure 'samtools' is available",
                " on your device or change 'SAMtools.or.Rsamtools' parameter ",
                "to 'Rsamtools'!!\n\n")
      }
    } else if (SAMtools.or.Rsamtools == "Rsamtools") {
      command.list <- c(command.list,
                        "* Rsamtools Converting '.sam' to '.bam' : ")
      sample.table <- table(gsub(paste0(".sam$"),
                                 replacement = "",
                                 check.results$sam.files.df))
      iteration.num <- length(sample.table)
      sample.name <- names(sample.table)
      sam.files <- lapply(sample.name, function(x) {paste0(path.prefix,
                                                           "gene_data/raw_sam/",
                                                           x, ".sam")})
      bam.files <- lapply(sample.name,
                          function(x) {paste0(path.prefix,
                                              "gene_data/raw_bam/",x)})
      command.list <- c(command.list,
                        paste0("     Running 'asBam' in Rsamtools in parallel:",
                               Rsamtools.nCores, "\n"))
      mcmapply(FUN = function(input.sam,output.bam) {
        message("     Input SAM file: '", input.sam, "'\n")
        message("     Creating BAM file: '", output.bam, ".bam'\n")
        # Remove tmp file
        # tmp_file <- tempfile()
        # unlink(list.files(dirname(dirname(tmp_file)),
        #                   pattern = "Rtmp*", full.names = TRUE),
        #        recursive = TRUE, force = TRUE)
        Rsamtools::asBam(file = input.sam,
                         destination = output.bam,
                         overwrite = TRUE)
        message("     Ouput BAM file: '", output.bam, ".bam' is created\n")
      }, sam.files, bam.files, mc.cores = Rsamtools.nCores)
      # command.list <- c(command.list, output.log)

      # iteration.num <- length(sample.table)
      # sample.name <- names(sample.table)
      # sample.value <- as.vector(sample.table)
      # for( i in seq_len(iteration.num)){
      #   input.sam <- paste0(path.prefix,
      #                       "gene_data/raw_sam/",
      #                       sample.name[i], ".sam")
      #   output.bam <- paste0(path.prefix,
      #                        "gene_data/raw_bam/",
      #                        sample.name[i])
      #   message("     Input SAM file: '", input.sam, "'\n")
      #   message("     Creating BAM file: '", output.bam, ".bam'\n")
      #   ba.file <- Rsamtools::asBam(file = input.sam,
      #                               destination = output.bam,
      #                               overwrite = TRUE,
      #                               maxMemory = Rsamtools.maxMemory)
      #   message("     Ouput BAM file: '", output.bam, ".bam' is created\n")
      #   command.list <- c(command.list, ba.file)
      # }
    }
    message("\n")
    command.list <- c(command.list, "\n")
    fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
    write(command.list, fileConn, append = TRUE)
  } else {
    stop(c("(\u2718) 'XXX.sam' is missing.\n\n"))
  }
}

# stringtie assemble and quantify expressed genes and transcripts
StringTieAssemble <- function(path.prefix,
                              genome.name,
                              sample.pattern,
                              Stringtie.Assembly.num.parallel.threads,
                              Stringtie.Assembly.f,
                              Stringtie.Assembly.m,
                              Stringtie.Assembly.c,
                              Stringtie.Assembly.g,
                              Stringtie.Assembly.M) {
  if (isTRUE(CheckStringTie(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Assembly :\n")
    message("************** Stringtie assembly **************\n")
    if (check.results$bam.files.number.df != 0){
      command.list <- c()
      command.list <- c(command.list, "* Stringtie assembly : ")
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      sample.name <- sort(gsub(paste0(".bam$"),
                               replacement = "",
                               check.results$bam.files.df))
      iteration.num <- length(sample.name)
      for( i in seq_len(iteration.num)){
        whole.command <- paste("-p", Stringtie.Assembly.num.parallel.threads,
                               "-f", Stringtie.Assembly.f, "-m", Stringtie.Assembly.m,
                               "-c", Stringtie.Assembly.c, "-g", Stringtie.Assembly.g,
                               "-M", Stringtie.Assembly.M,
                               "-G", paste0("ref_genes/", genome.name, ".gtf"),
                               "-o", paste0("raw_gtf/", sample.name[i], ".gtf"),
                               "-l", sample.name[i],
                               paste0("raw_bam/", sample.name[i], ".bam"))
        if (i != 1) message("\n")
        main.command <- "stringtie"
        message("Input command :", paste(main.command, whole.command), "\n")
        command.list <- c(command.list,
                          paste("    command :", main.command, whole.command))
        command.result <- system2(command = main.command, args = whole.command)
        if (command.result != 0 ) {
          on.exit(setwd(current.path))
          message("(\u2718) '", main.command, "' is failed !!")
          stop("'", main.command, "' ERROR")
        }
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
    } else {
      on.exit(setwd(current.path))
      stop("(\u2718) '", genome.name, ".gtf' or 'XXX.bam' is missing.\n\n")
    }
  }
}

# stringtie merge transcripts from all samples
StringTieMergeTrans <- function(path.prefix,
                                genome.name,
                                sample.pattern,
                                Stringtie.Merge.num.parallel.threads) {
  if (isTRUE(CheckStringTie(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Transcript Merging :\n")
    message("************** ",
            "Stringtie merging transcripts **************\n")
    if ( isTRUE(check.results$gtf.file.logic.df) &&
         check.results$gtf.files.number.df != 0){
      command.list <- c()
      command.list <- c(command.list, "* Stringtie Merging Transcripts : ")
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      if(!dir.exists(paste0(path.prefix, "gene_data/merged/"))){
        dir.create(paste0(path.prefix, "gene_data/merged/"))
      }
      sample.table <- table(check.results$gtf.files.df)
      iteration.num <- length(sample.table)
      sample.name <- names(sample.table)
      sample.value <- as.vector(sample.table)
      write.content <- paste0("raw_gtf/",sample.name[1])
      for (i in 2:iteration.num){
        write.content <- c(write.content, paste0("raw_gtf/" ,sample.name[i]))
      }
      write.file<-file("merged/mergelist.txt")
      writeLines(write.content, write.file)
      close(write.file)
      whole.command <- paste("--merge -p", num.parallel.threads.stringtie.merge,
                             "-G", paste0("ref_genes/", genome.name, ".gtf"),
                             "-o", "merged/stringtie_merged.gtf",
                             "merged/mergelist.txt")
      main.command <- "stringtie"
      message("Input command :", paste(main.command, whole.command), "\n")
      command.list <- c(command.list,
                        paste("    command :", main.command, whole.command))
      command.result <- system2(command = main.command, args = whole.command)
      if (command.result != 0 ) {
        on.exit(setwd(current.path))
        message("(\u2718) '", main.command, "' is failed !!")
        stop("'", main.command, "' ERROR")
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
    } else {
      on.exit(setwd(current.path))
      stop("(\u2718) :'", genome.name, ".gtf' or 'XXX.gtf' is missing.\n\n")
    }
  }
}

# stringtie estimate transcript abundances and create table count for Ballgown
StringTieToBallgown <- function(path.prefix,
                                genome.name,
                                sample.pattern,
                                Stringtie.2.Ballgown.num.parallel.threads) {
  if (isTRUE(CheckStringTie(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Ballgown Table count Creation :\n")
    message("************** ",
            "Stringtie creating table count for Ballgown **************\n")
    if ((check.results$bam.files.number.df != 0) &&
        isTRUE(check.results$stringtie_merged.gtf.file.df)){
      command.list <- c()
      command.list <- c(command.list,
                        "* Stringtie Creating Table count for Ballgown : ")
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      sample.table <- table(gsub(paste0(".bam$"),
                                 replacement = "",
                                 check.results$bam.files.df))
      iteration.num <- length(sample.table)
      sample.name <- names(sample.table)
      sample.value <- as.vector(sample.table)
      for( i in seq_len(iteration.num)){
        # '-e' only estimate the abundance of given reference transcripts
        # (requires -G)
        whole.command <- paste("-e -B -p", num.parallel.threads.stringtie.2.ballgown,
                               "-G", "merged/stringtie_merged.gtf",
                               "-o", paste0("ballgown/", sample.name[i],
                                            "/", sample.name[i], ".gtf"),
                               "-A", paste0("gene_abundance/", sample.name[i],
                                            "/", sample.name[i], ".tsv"),
                               paste0("raw_bam/", sample.name[i], ".bam"))
        main.command <- 'stringtie'
        if (i != 1) message("\n")
        message("Input command :", paste(main.command, whole.command), "\n")
        command.list <- c(command.list,
                          paste("    command :", main.command, whole.command))
        command.result <- system2(command = main.command, args = whole.command)
        if (command.result != 0 ) {
          on.exit(setwd(current.path))
          message("(\u2718) '", main.command, "' is failed !!")
          stop("'", main.command, "' ERROR")
        }
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
    } else {
      on.exit(setwd(current.path))
      stop("(\u2718) 'stringtie_merged.gtf' or 'XXX.bam' is missing.\n\n")
    }
  }
}

# Examine how the transcripts compare with the reference annotation
GffcompareRefSample <- function(path.prefix,
                                genome.name,
                                sample.pattern) {
  if (isTRUE(CheckGffcompare(print=FALSE))){
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Merged `GTF` & Reference Comparison :\n")
    message("************** ",
            "Gffcompare comparing transcripts between ",
            "merged and reference **************\n")
    if ( isTRUE(check.results$stringtie_merged.gtf.file.df) &&
         isTRUE(check.results$gtf.file.logic.df)){
      command.list <- c()
      command.list <- c(command.list, "* Gffcompare Comparing ",
                        "Transcripts between Merged and Reference : ")
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      whole.command <- paste("-r", paste0("ref_genes/", genome.name, ".gtf"),
                             "-G -o merged/merged merged/stringtie_merged.gtf")
      main.command <- "gffcompare"
      message("Input command : ", paste(main.command, whole.command), "\n")
      command.list <- c(command.list,
                        paste("    command :", main.command, whole.command))
      command.result <- system2(command = main.command, args = whole.command)
      if (command.result != 0 ) {
        on.exit(setwd(current.path))
        message("(\u2718) '", main.command, "' is failed !!")
        stop("'", main.command, "' ERROR")
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
    } else {
      on.exit(setwd(current.path))
      stop("(\u2718) '", genome.name,
           ".gtf' or 'stringtie_merged.gtf' is missing.\n\n")
    }
  }
}

# converting stringtie ballogwn preprocessed data to count table
PreDECountTable <- function(path.prefix,
                            sample.pattern,
                            python.variable.answer,
                            python.variable.version,
                            python.2to3,
                            print=TRUE) {
  # ftp server :
  # 'prepDE.py' is from
  #  ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-source.zip
  message("\n\u2618\u2618\u2618 Raw Reads Count Creation :\n")
  message("************** Reading prepDE.py ************\n")
  if(!dir.exists(paste0(path.prefix, "gene_data/reads_count_matrix/"))){
    dir.create(file.path(paste0(path.prefix, 'gene_data/reads_count_matrix/')),
               showWarnings = FALSE)
  }
  url <- "https://ccb.jhu.edu/software/stringtie/dl/prepDE.py"
  command.list <- c()
  command.list <- c(command.list, "* Installing 'prepDE.py' : ")
  message(path.prefix, "gene_data/reads_count_matrix\n")
  current.path <- getwd()
  setwd(paste0(path.prefix, "gene_data/reads_count_matrix/"))
  file.download <- getURL(url, download.file,
                          paste0(path.prefix,
                                 "gene_data/reads_count_matrix/prepDE.py"))

  message("Using R function : 'download.file()' is called. \n")
  command.list <- c(command.list,
                    "    Using R function : 'download.file()' is called.")
  command.list <- c(command.list, "\n")
  if (file.download != 0 ) {
    on.exit(setwd(current.path))
    message("(\u2718) '", main.command, "' is failed !!")
    stop("'", main.command, "' ERROR")
  }
  message("'", path.prefix,
          "gene_data/reads_count_matrix/prepDE.py' has been installed.\n\n")
  message("************** Creating 'sample_lst.txt' file ************\n")
  sample.files <- list.files(paste0(path.prefix, "gene_data/ballgown/"),
                             pattern = sample.pattern)
  write.content <- paste0(sample.files[1], " ", path.prefix,
                          "gene_data/ballgown/", sample.files[1] ,
                          "/", sample.files[1], ".gtf")
  for(i in 2:length(sample.files)){
    write.content <- c(write.content,
                       paste0(sample.files[i], " ", path.prefix,
                              "gene_data/ballgown/",  sample.files[i],
                              "/",sample.files[i], ".gtf"))
  }
  write.file<-file(paste0(path.prefix,
                          "gene_data/reads_count_matrix/sample_lst.txt"))
  writeLines(write.content, write.file)
  close(write.file)
  message("'", path.prefix,
          "gene_data/reads_count_matrix/sample_lst.txt' has been created\n\n")
  message("************** ",
          "Creating gene and transcript raw count file ",
          "************\n")
  # have to check python !!!
  if (python.variable.answer) {
    message("(\u2714) : Python is available on your device!\n")
    message("       Python version : ",
            reticulate::py_config()$version, "\n")
    if(python.variable.version >= 3) {
      ## If Python3 ==> check whether 2to3 variable is valid!
      if (isTRUE(python.2to3)) {
        message("(\u270D) : Converting 'prepDE.py' ",
                "from python2 to python3 \n\n")
        whole.command <- paste0("-W ", path.prefix,
                                "gene_data/reads_count_matrix/prepDE.py ",
                                "--no-diffs")
        main.command <- "2to3"
        message("Input command :",
                paste(main.command, whole.command), "\n\n")
        command.list <- c(command.list, "* Coverting prepDE.py to Python3 : ")
        command.list <- c(command.list,
                          paste("    command :", main.command, whole.command))
        command.list <- c(command.list, "\n")
        command.result <- system2(command = main.command, args = whole.command)
        if (command.result != 0 ) {
          on.exit(setwd(current.path))
          message("(\u2718) '", main.command, "' is failed !!")
          stop("'", main.command, "' ERROR")
        }
      } else {
        on.exit(setwd(current.path))
        message("(\u26A0) 2to3 command is not available on your device !\n\n'")
        return(TRUE)
      }
    }
    message("\n(\u2714) : Creating gene and transcript raw count files now !\n")
    whole.command <- paste0(path.prefix,
                            "gene_data/reads_count_matrix/prepDE.py -i ",
                            path.prefix,
                            "gene_data/reads_count_matrix/sample_lst.txt")
    main.command <- "python"
    message("Input command :", paste(main.command, whole.command), "\n\n")
    command.list <- c(command.list,
                      "* Creating Gene and Transcript Raw Count File : ")
    command.list <- c(command.list,
                      paste("    command :", main.command, whole.command))
    command.result <- system2(command = main.command,
                              args = whole.command,
                              wait = TRUE)
    if (command.result != 0 ) {
      on.exit(setwd(current.path))
      message("(\u2718) '", main.command, "' is failed !!")
      stop("'", main.command, "' ERROR")
    }
    message("'", path.prefix,
            "gene_data/reads_count_matrix/",
            "gene_count_matrix.csv' has been created\n")
    message("'", path.prefix,
            "gene_data/reads_count_matrix/",
            "transcript_count_matrix.csv' has been created\n\n")
    command.list <- c(command.list, "\n")
    fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
    write(command.list, fileConn, append = TRUE)
    on.exit(setwd(current.path))
    return(TRUE)
  } else {
    on.exit(setwd(current.path))
    message("(\u26A0) Python is not available on your device!! ",
            "Please install python to run ",
            "python script 'prepDE.py'. ",
            "Raw reads count table creation is skipped!!\n\n'")
    return(TRUE)
  }
}

#
# # Report Hisat2 assemble rate
# Hisat2ReportAssemble <- function(path.prefix,
#                                  genome.name,
#                                  sample.pattern){
#   check.results <- ProgressGenesFiles(path.prefix,
#                                       genome.name,
#                                       sample.pattern,
#                                       print=FALSE)
#   message("\n************** Reporting Hisat2 Alignment **************\n")
#   if (check.results$phenodata.file.df &&
#       check.results$bam.files.number.df != 0){
#     file.read <- paste0(path.prefix, "Rscript_out/Read_Process.Rout")
#     sample.name <- sort(gsub(paste0(".bam$"),
#                              replacement = "",
#                              check.results$bam.files.df))
#     iteration.num <- length(sample.name)
#     load.data <- readChar(file.read, file.info(file.read)$size)
#     # overall alignment rate
#     overall.alignment <- strsplit(load.data, "\n")
#     overall.alignment.with.NA <-
#       stringr::str_extract(overall.alignment[[1]],
#                            "[0-9]*.[0-9]*% overall alignment rate")
#     overall.alignment.result <-
#       overall.alignment.with.NA[!is.na(overall.alignment.with.NA)]
#     overall.alignment.result.cut <- gsub(" overall alignment rate",
#                                          " ",
#                                          overall.alignment.result)
#     # different mapping rate
#     first.split <- strsplit(load.data,
#                             paste0("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*",
#                                    "\\* Hisat2 Alignment \\*",
#                                    "\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\n"))
#     second.split <- strsplit(first.split[[1]][2],
#                              paste0("\n\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*",
#                                     "\\* Current progress of RNA-seq files in"))
#     split.lines <- strsplit(second.split[[1]][1], "\n")
#     alignment.rate.with.NA <-
#       stringr::str_extract(split.lines[[1]],
#                            "[0-9]* \\([0-9]*.[0-9]*%\\) aligned concordantly")
#     alignment.first.result <-
#       alignment.rate.with.NA[!is.na(alignment.rate.with.NA)]
#     alignment.first.result.cut1 <- gsub(") aligned concordantly",
#                                         " ",
#                                         alignment.first.result)
#     alignment.first.result.cut2 <- gsub("[0-9]* \\(",
#                                         " ",
#                                         alignment.first.result.cut1)
#     report.data.frame <- data.frame(matrix(0, ncol = 0, nrow = 3))
#     row.names(report.data.frame) <- c("Unique mapping rate",
#                                       "Multiple mapping rate",
#                                       "Overall alignment rate")
#     for( i in seq_len(iteration.num)){
#       add.column <- c()
#       for( j in (i*3-1):(i*3)){
#         add.column <- c(add.column, alignment.first.result.cut2[j])
#       }
#       add.column <- c(add.column, overall.alignment.result.cut[i])
#       report.data.frame[[(sample.name[i])]] <- add.column
#     }
#     if(!dir.exists(paste0(path.prefix, "RNASeq_results/Alignment_Report/"))){
#       dir.create(paste0(path.prefix, "RNASeq_results/Alignment_Report/"))
#     }
#     write.csv(report.data.frame,
#               file = paste0(path.prefix,
#                             "RNASeq_results/",
#                             "Alignment_Report/Alignment_report.csv"))
#     png(paste0(path.prefix,
#                "RNASeq_results/Alignment_Report/Alignment_report.png"),
#         width = iteration.num*100 + 200, height = 40*4)
#     p <- gridExtra::grid.table(report.data.frame)
#     print(p)
#     dev.off()
#     message("Results are in ",
#             paste0("'", path.prefix, "RNASeq_results/Alignment_Report/'"),
#             "\n\n")
#   }
# }

