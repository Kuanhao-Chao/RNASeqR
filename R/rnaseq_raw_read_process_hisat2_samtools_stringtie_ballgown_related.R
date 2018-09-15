# Creating Hisat2 index
# Creat Hisat2 index for further use
CreateHisat2Index <- function (path.prefix,
                               genome.name,
                               sample.pattern,
                               splice.site.info = TRUE,
                               exon.info = TRUE) {
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
      message(paste0("************** Creating Hisat2 Index **************\n"))
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
            message("(\u2718) '", main.command, "' is failed !!")
            stop(paste0("'", main.command, "' ERROR"))
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
            message(paste0("(\u2718) '", main.command, "' is failed !!"))
            stop(paste0("'", main.command, "' ERROR"))
          }
          message("\n")
        }
        if (isTRUE(splice.site.info) && isTRUE(exon.info)) {
          whole.command <- paste(paste('--ss', paste0(genome.name, '.ss'),
                                       "--exon", paste0(genome.name, '.exon'),
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
            message("(\u2718) '", main.command, "' is failed !!")
            stop(paste0("'", main.command, "' ERROR"))
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
            message("(\u2718) '", main.command, "' is failed !!")
            stop(paste0("'", main.command, "' ERROR"))
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
            message("(\u2718) '", main.command, "' is failed !!")
            stop(paste0("'", main.command, "' ERROR"))
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
            message("(\u2718) '", main.command, "' is failed !!")
            stop(paste0("'", main.command, "' ERROR"))
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
        stop(c(paste0("(\u2718) '", genome.name, ".gtf' "), "or",
               paste0(" '", genome.name, ".fa'"), "is missing.\n\n"))
      }
    }
  }
}

# hisat2 alignment default
Hisat2AlignmentDefault <- function(path.prefix,
                                   genome.name,
                                   sample.pattern,
                                   num.parallel.threads = 8) {
  if (isTRUE(CheckHisat2(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Alignment :\n")
    message(paste0("************** Hisat2 Alignment **************\n"))
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
        stop("(\u2718) Inconsistent formats. Please check files are all",
             paste0("'XXX_*.fastq.gz'"), "\n\n")
      } else {
        sample.table <- table(gsub(paste0("_[1-2]*.fastq.gz$"),
                                   replacement = "",
                                   check.results$fastq.gz.files.df))
        iteration.num <- length(sample.table)
        sample.name <- names(sample.table)
        sample.value <- as.vector(sample.table)
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
          whole.command <- paste("-p", num.parallel.threads,"--dta -x",
                                 paste0("indices/", genome.name, "_tran"),
                                 total.sub.command, "-S",
                                 paste0("raw_sam/", sample.name[i],".sam") )
          if (i != 1) message("\n")
          main.command <- "hisat2"
          message("Input command : ", paste(main.command, whole.command), "\n")
          command.list <- c(command.list,
                            paste("    command :", main.command, whole.command))
          command.result <- system2(command = main.command,
                                    args = whole.command)
          if (command.result != 0 ) {
            message(paste0("(\u2718) '", main.command, "' is failed !!"))
            stop(paste0("'", main.command, "' ERROR"))
          }
        }
        message("\n")
        command.list <- c(command.list, "\n")
        fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
        write(command.list, fileConn, append = TRUE)
        on.exit(setwd(current.path))
      }
    } else {
      stop(c(paste0("(\u2718) '", genome.name, "_tran.*.ht2' "),
             "or 'XXX_*.fastq.gz' is missing.\n\n"))
    }
  }
}

# Report Hisat2 assemble rate
Hisat2ReportAssemble <- function(path.prefix,
                                 genome.name,
                                 sample.pattern){
  check.results <- ProgressGenesFiles(path.prefix,
                                      genome.name,
                                      sample.pattern,
                                      print=FALSE)
  message("\n************** Reporting Hisat2 Alignment **************\n")
  if (check.results$phenodata.file.df &&
      check.results$bam.files.number.df != 0){
    file.read <- paste0(path.prefix, "Rscript_out/Read_Process.Rout")
    sample.name <- sort(gsub(paste0(".bam$"),
                             replacement = "",
                             check.results$bam.files.df))
    iteration.num <- length(sample.name)
    load.data <- readChar(file.read, file.info(file.read)$size)
    # overall alignment rate
    overall.alignment <- strsplit(load.data, "\n")
    overall.alignment.with.NA <-
      stringr::str_extract(overall.alignment[[1]],
                           "[0-9]*.[0-9]*% overall alignment rate")
    overall.alignment.result <-
      overall.alignment.with.NA[!is.na(overall.alignment.with.NA)]
    overall.alignment.result.cut <- gsub(" overall alignment rate",
                                         " ",
                                         overall.alignment.result)
    # different mapping rate
    first.split <- strsplit(load.data,
                            paste0("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*",
                                   "\\* Hisat2 Alignment \\*",
                                   "\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\n"))
    second.split <- strsplit(first.split[[1]][2],
                             paste0("\n\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*",
                                    "\\* Current progress of RNA-seq files in"))
    split.lines <- strsplit(second.split[[1]][1], "\n")
    alignment.rate.with.NA <-
      stringr::str_extract(split.lines[[1]],
                           "[0-9]* \\([0-9]*.[0-9]*%\\) aligned concordantly")
    alignment.first.result <-
      alignment.rate.with.NA[!is.na(alignment.rate.with.NA)]
    alignment.first.result.cut1 <- gsub(") aligned concordantly",
                                        " ",
                                        alignment.first.result)
    alignment.first.result.cut2 <- gsub("[0-9]* \\(",
                                        " ",
                                        alignment.first.result.cut1)
    report.data.frame <- data.frame(matrix(0, ncol = 0, nrow = 3))
    row.names(report.data.frame) <- c("Unique mapping rate",
                                      "Multiple mapping rate",
                                      "Overall alignment rate")
    for( i in seq_len(iteration.num)){
      add.column <- c()
      for( j in (i*3-1):(i*3)){
        add.column <- c(add.column, alignment.first.result.cut2[j])
      }
      add.column <- c(add.column, overall.alignment.result.cut[i])
      report.data.frame[[(sample.name[i])]] <- add.column
    }
    if(!dir.exists(paste0(path.prefix, "RNASeq_results/Alignment_Report/"))){
      dir.create(paste0(path.prefix, "RNASeq_results/Alignment_Report/"))
    }
    write.csv(report.data.frame,
              file = paste0(path.prefix,
                            "RNASeq_results/",
                            "Alignment_Report/Alignment_report.csv"))
    png(paste0(path.prefix,
               "RNASeq_results/Alignment_Report/Alignment_report.png"),
        width = iteration.num*100 + 200, height = 40*4)
    p <- gridExtra::grid.table(report.data.frame)
    print(p)
    dev.off()
    message("Results are in ",
            paste0("'", path.prefix, "RNASeq_results/Alignment_Report/'"),
            "\n\n")
  }
}

# use 'samtools' to sort and convert the SAM files to BAM
SamtoolsToBam <- function(path.prefix,
                          genome.name,
                          sample.pattern,
                          num.parallel.threads = 8) {
  if (isTRUE(CheckSamtools(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 'SAM' to 'BAM' :\n")
    message(paste0("************** ",
                   "Samtools converting '.sam' to '.bam' ",
                   "**************\n"))
    if (check.results$sam.files.number.df != 0){
      command.list <- c()
      command.list <- c(command.list,
                        "* Samtools Converting '.sam' to '.bam' : ")
      # Map reads to each alignment
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      sample.table <- table(gsub(paste0(".sam$"),
                                 replacement = "",
                                 check.results$sam.files.df))
      iteration.num <- length(sample.table)
      sample.name <- names(sample.table)
      sample.value <- as.vector(sample.table)
      for( i in seq_len(iteration.num)){
        whole.command <- paste("sort -@", num.parallel.threads,
                               "-o", paste0("raw_bam/", sample.name[i], ".bam"),
                               paste0("raw_sam/", sample.name[i], ".sam"))
        if (i != 1) message("\n")
        main.command <- "samtools"
        message(c("Input command :", paste(main.command, whole.command), "\n"))
        command.list <- c(command.list,
                          paste("    command :", main.command, whole.command))
        command.result <- system2(command = main.command, args = whole.command)
        if (command.result != 0 ) {
          message(paste0("(\u2718) '", main.command, "' is failed !!"))
          stop(paste0("'", main.command, "' ERROR"))
        }
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
    } else {
      stop(c("(\u2718) 'XXX.sam' is missing.\n\n"))
    }
  }
}

# stringtie assemble and quantify expressed genes and transcripts
StringTieAssemble <- function(path.prefix,
                              genome.name,
                              sample.pattern,
                              num.parallel.threads = 8) {
  if (isTRUE(CheckStringTie(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Assembly :\n")
    message(paste0("************** Stringtie assembly **************\n"))
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
        whole.command <- paste("-p", num.parallel.threads,
                               "-G",paste0("ref_genes/", genome.name, ".gtf"),
                               "-o", paste0("raw_gtf/", sample.name[i], ".gtf"),
                               "-l", sample.name[i],
                               paste0("raw_bam/", sample.name[i], ".bam"))
        if (i != 1) message("\n")
        main.command <- "stringtie"
        message(c("Input command :", paste(main.command, whole.command), "\n"))
        command.list <- c(command.list,
                          paste("    command :", main.command, whole.command))
        command.result <- system2(command = main.command, args = whole.command)
        if (command.result != 0 ) {
          message(paste0("(\u2718) '", main.command, "' is failed !!"))
          stop(paste0("'", main.command, "' ERROR"))
        }
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
    } else {
      stop(c(paste0("(\u2718) '", genome.name, ".gtf' "),
             "or 'XXX.bam' is missing.\n\n"))
    }
  }
}

# stringtie merge transcripts from all samples
StringTieMergeTrans <- function(path.prefix,
                                genome.name,
                                sample.pattern,
                                num.parallel.threads = 8) {
  if (isTRUE(CheckStringTie(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Transcript Merging :\n")
    message(paste0("************** ",
                   "Stringtie merging transcripts ",
                   "**************\n"))
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
      whole.command <- paste("--merge -p", num.parallel.threads,
                             "-G", paste0("ref_genes/", genome.name, ".gtf"),
                             "-o", "merged/stringtie_merged.gtf",
                             "merged/mergelist.txt")
      main.command <- "stringtie"
      message(c("Input command :", paste(main.command, whole.command), "\n"))
      command.list <- c(command.list,
                        paste("    command :", main.command, whole.command))
      command.result <- system2(command = main.command, args = whole.command)
      if (command.result != 0 ) {
        message(paste0("(\u2718) '", main.command, "' is failed !!"))
        stop(paste0("'", main.command, "' ERROR"))
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
    } else {
      stop(c(paste0("(\u2718) :'", genome.name, ".gtf' "),
             "or", " 'XXX.gtf' is missing.\n\n"))
    }
  }
}

# stringtie estimate transcript abundances and create table count for Ballgown
StringTieToBallgown <- function(path.prefix,
                                genome.name,
                                sample.pattern,
                                num.parallel.threads = 8) {
  if (isTRUE(CheckStringTie(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix,
                                        genome.name,
                                        sample.pattern,
                                        print=TRUE)
    message("\n\u2618\u2618\u2618 Ballgown Table count Creation :\n")
    message(paste0("************** ",
                   "Stringtie creating table count for Ballgown ",
                   "**************\n"))
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
        whole.command <- paste("-e -B -p", num.parallel.threads,
                               "-G", "merged/stringtie_merged.gtf",
                               "-o", paste0("ballgown/", sample.name[i],
                                            "/", sample.name[i], ".gtf"),
                               "-A", paste0("gene_abundance/", sample.name[i],
                                            "/", sample.name[i], ".tsv"),
                               paste0("raw_bam/", sample.name[i], ".bam"))
        main.command <- 'stringtie'
        if (i != 1) message("\n")
        message(c("Input command :", paste(main.command, whole.command), "\n"))
        command.list <- c(command.list,
                          paste("    command :", main.command, whole.command))
        command.result <- system2(command = main.command, args = whole.command)
        if (command.result != 0 ) {
          message(paste0("(\u2718) '", main.command, "' is failed !!"))
          stop(paste0("'", main.command, "' ERROR"))
        }
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
    } else {
      stop(c(paste0("(\u2718) 'stringtie_merged.gtf' "),
             "or", " 'XXX.bam' is missing.\n\n"))
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
    message(paste0("************** ",
                   "Gffcompare comparing transcripts between ",
                   "merged and reference **************\n"))
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
        message("(\u2718) '", main.command, "' is failed !!")
        stop(paste0("'", main.command, "' ERROR"))
      }
      message("\n")
      command.list <- c(command.list, "\n")
      fileConn <- paste0(path.prefix, "RNASeq_results/COMMAND.txt")
      write(command.list, fileConn, append = TRUE)
      on.exit(setwd(current.path))
    } else {
      stop(c(paste0("(\u2718) '", genome.name, ".gtf' "),
             "or", paste0(" 'stringtie_merged.gtf'"), "is missing.\n\n"))
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
  #  ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-source.zip
  message("\n\u2618\u2618\u2618 Raw Reads Count Creation :\n")
  message("************** Installing prepDE.py ************\n")
  if(!dir.exists(paste0(path.prefix, "gene_data/reads_count_matrix/"))){
    dir.create(file.path(paste0(path.prefix, 'gene_data/reads_count_matrix/')),
               showWarnings = FALSE)
  }
  command.list <- c()
  command.list <- c(command.list, "* Installing prepDE.py : ")
  message(paste0(path.prefix, "gene_data/reads_count_matrix\n"))
  current.path <- getwd()
  setwd(paste0(path.prefix, "gene_data/reads_count_matrix/"))
  whole.command <- paste('https://ccb.jhu.edu/software/stringtie/dl/prepDE.py',
                         '--output',
                         paste0(path.prefix,
                                "gene_data/reads_count_matrix/prepDE.py"))
  main.command <- "curl"
  message(c("Input command :", paste(main.command, whole.command), "\n"))
  command.list <- c(command.list,
                    paste("    command :", main.command, whole.command))
  command.list <- c(command.list, "\n")
  command.result <- system2(command = main.command,
                            args = whole.command,
                            stdout = "",
                            wait = TRUE)
  if (command.result != 0 ) {
    message(paste0("(\u2718) '", main.command, "' is failed !!"))
    stop(paste0("'", main.command, "' ERROR"))
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
    message(paste0("       Python version : ",
                   reticulate::py_config()$version, "\n"))
    if(python.variable.version >= 3) {
      ## If Python3 ==> check whether 2to3 variable is valid!
      if (isTRUE(python.2to3)) {
        message("(\u270D) : Converting 'prepDE.py' ",
                "from python2 to python3 \n\n")
        whole.command <- paste0("-W ", path.prefix,
                                "gene_data/reads_count_matrix/prepDE.py ",
                                "--no-diffs")
        main.command <- "2to3"
        message(c("Input command :",
                  paste(main.command, whole.command), "\n\n"))
        command.list <- c(command.list, "* Coverting prepDE.py to Python3 : ")
        command.list <- c(command.list,
                          paste("    command :", main.command, whole.command))
        command.list <- c(command.list, "\n")
        command.result <- system2(command = main.command, args = whole.command)
        if (command.result != 0 ) {
          message(paste0("(\u2718) '", main.command, "' is failed !!"))
          stop(paste0("'", main.command, "' ERROR"))
        }
      } else {
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
    message(c("Input command :", paste(main.command, whole.command), "\n\n"))
    command.list <- c(command.list,
                      "* Creating Gene and Transcript Raw Count File : ")
    command.list <- c(command.list,
                      paste("    command :", main.command, whole.command))
    command.result <- system2(command = main.command,
                              args = whole.command,
                              wait = TRUE)
    if (command.result != 0 ) {
      message(paste0("(\u2718) '", main.command, "' is failed !!"))
      stop(paste0("'", main.command, "' ERROR"))
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
    ## Fix !!
    on.exit(setwd(current.path))
    message(paste0("(\u26A0) Python is not available on your device!! ",
                   "Please install python to run ",
                   "python script 'prepDE.py'\n\n'"))
    return(TRUE)
  }
}
