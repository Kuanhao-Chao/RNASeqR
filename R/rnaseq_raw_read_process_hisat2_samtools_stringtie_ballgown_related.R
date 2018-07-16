#' Creating Hisat2 index
#' Creat Hisat2 index for further use
CreateHisat2Index <- function (path.prefix, gene.name, sample.pattern, splice.site.info = TRUE, exon.info = TRUE) {
  ## need to change 'Let the user can choose where to put their files'
  ## need to learn how to get the filenames under a directory
  ## need to use 'real filename' rather than the 'chrX'
  if (isTRUE(CheckHisat2(print=FALSE))){
    if (!is.logical(splice.site.info) || !is.logical(exon.info)) {
      stop("(\u2718) Please make sure the type of 'splice.site.info' and 'exon.info' are logical.\n")
    } else {
      check.results <- ProgressGenesFiles(path.prefix, path.prefix, gene.name, sample.pattern, print=TRUE)
      cat(paste0("\n************** Creating Hisat2 Index **************\n"))
      if (isTRUE(check.results$gtf.file.logic.df) && isTRUE(check.results$fa.file.logic.df)){
        current.path <- getwd()
        setwd(paste0(path.prefix, "gene_data/indexes/"))
        if (isTRUE(splice.site.info)) {
          cat(c("Input command :", paste("extract_splice_sites.py", paste0(path.prefix, 'gene_data/ref_genes/', gene.name, '.gtf'), '>', paste0(gene.name, '.ss')), "\n"))
          system2(command = 'extract_splice_sites.py', args = c(paste0(path.prefix, 'gene_data/ref_genes/', gene.name, '.gtf'), '>', paste0(gene.name, '.ss')))
          cat("\n")
        }
        if (isTRUE(exon.info)) {
          cat(c("Input command :", paste("extract_exons.py", paste0(path.prefix, 'gene_data/ref_genes/', gene.name, '.gtf'), '>', paste0(gene.name, '.exon')), "\n"))
          system2(command = 'extract_exons.py', args = c(paste0(path.prefix, 'gene_data/ref_genes/', gene.name, '.gtf'), '>', paste0(gene.name, '.exon')))
          cat("\n")
        }

        if (isTRUE(splice.site.info) && isTRUE(exon.info)) {
          cat(c("Input command :", paste("hisat2-build", paste('--ss', paste0(gene.name, '.ss'), '--exon', paste0(gene.name, '.exon'), paste0(path.prefix, 'gene_data/ref_genome/', gene.name, '.fa'), paste0(gene.name, '_tran')), "\n")))
          system2(command = 'hisat2-build', args = c('--ss', paste0(gene.name, '.ss'), '--exon', paste0(gene.name, '.exon'), paste0(path.prefix, 'gene_data/ref_genome/', gene.name, '.fa'), paste0(gene.name, '_tran')))
          cat("\n")
        } else if (isTRUE(splice.site.info) && !isTRUE(exon.info)) {
          cat(c("Input command :", paste("hisat2-build", paste('--ss', paste0(gene.name, '.ss'), paste0(path.prefix, 'gene_data/ref_genome/', gene.name, '.fa'), paste0(gene.name, '_tran')), "\n")))
          system2(command = 'hisat2-build', args = c('--ss', paste0(gene.name, '.ss'), paste0(path.prefix, 'gene_data/ref_genome/', gene.name, '.fa'), paste0(gene.name, '_tran')))
          cat("\n")
        } else if (!isTRUE(splice.site.info) && isTRUE(exon.info)) {
          cat(c("Input command :", paste("hisat2-build", paste('--exon', paste0(gene.name, '.exon'), paste0(path.prefix, 'gene_data/ref_genome/', gene.name, '.fa'), paste0(gene.name, '_tran')), "\n")))
          system2(command = 'hisat2-build', args = c('--exon', paste0(gene.name, '.exon'), paste0(path.prefix, 'gene_data/ref_genome/', gene.name, '.fa'), paste0(gene.name, '_tran')))
          cat("\n")
        } else {
          cat(c("Input command :", paste("hisat2-build", paste(paste0(path.prefix, 'gene_data/ref_genome/', gene.name, '.fa'), paste0(gene.name, '_tran')), "\n")))
          system2(command = 'hisat2-build', args = c(paste0(path.prefix, 'gene_data/ref_genome/', gene.name, '.fa'), paste0(gene.name, '_tran')))
          cat("\n")
        }
        on.exit(setwd(current.path))
        cat(paste0("'", path.prefix, "gene_data/indexes/", gene.name, "_tran.*.ht2' has been created.\n\n"))
      } else {
        stop(c(paste0("(\u2718) '", gene.name, ".gtf' "), "or", paste0(" '", gene.name, ".fa'"), "is missing.\n\n"))
      }
    }
  }
}

#' hisat2 alignment default
Hisat2AlignmentDefault <- function(path.prefix, gene.name, sample.pattern, num.parallel.threads = 8) {
  if (isTRUE(CheckHisat2(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix, gene.name, sample.pattern, print=TRUE)
    cat(paste0("\n************** Hisat2 Aligning **************\n"))
    if (check.results$ht2.files.number.df != 0 && check.results$fastq.gz.files.number.df != 0){
      # Map reads to each alignment
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      # Determine 'r'/'R'/''
      deleteback <- gsub("[1-2]*.fastq.gz$", replace = "", check.results$fastq.gz.files.df)
      sample.table.r.value <- gsub(paste0("[A-Z, a-z]*[0-9]*_"), replace = "", deleteback)
      if (isTRUE(length(unique(sample.table.r.value)) != 1)){
        stop("(\u2718) Inconsistent formats. Please check files are all",  paste0("'XXX_r*.fastq.gz'"), "OR",  paste0("'XXX_R*.fastq.gz'"), "OR",  paste0("'XXX_*.fastq.gz'"), "\n\n")
      } else {
        sample.table <- table(gsub(paste0("_[R]*[r]*[1-2]*.fastq.gz$"), replace = "", check.results$fastq.gz.files.df))
        iteration.num <- length(sample.table)
        sample.name <- names(sample.table)
        sample.value <- as.vector(sample.table)
        for( i in 1:iteration.num){
          current.sub.command <- ""
          total.sub.command <- ""
          for ( j in 1:sample.value[i]){
            current.sub.command <- paste(paste0("-", j),  paste0("raw_fastq.gz/", sample.name[i], "_", sample.table.r.value[1], j, ".fastq.gz"))
            total.sub.command <- paste(total.sub.command, current.sub.command)
          }
          whole.command <- paste("-p", num.parallel.threads,"--dta -x", paste0("indexes/", gene.name, "_tran"), total.sub.command, "-S", paste0("raw_sam/", sample.name[i],".sam") )
          if (i != 1) cat("\n")
          cat(c("Input command :", paste("hisat2", whole.command), "\n"))
          system2(command = 'hisat2', args = whole.command)
        }
        cat("\n")
        on.exit(setwd(current.path))
      }
    } else {
      stop(c(paste0("(\u2718) '", gene.name, "_tran.*.ht2' "), "or 'XXX_*.fastq.gz' is missing.\n\n"))
    }
  }
}

#' Report Hisat2 assemble rate
#'
#'
Hisat2ReportAssemble <- function(path.prefix, gene.name, sample.pattern){
  check.results <- ProgressGenesFiles(path.prefix, gene.name, sample.pattern, print=FALSE)
  cat(paste0("\n************** Reporting hisat2 alignment **************\n"))
  if (isTRUE(check.results$phenodata.file.df) && check.results$bam.files.number.df != 0){
    file.read <- paste0(path.prefix, "Rscript_out/Raw_Read_Process.Rout")
    sample.name <- sort(gsub(paste0(".bam$"), replace = "", check.results$bam.files.df))
    iteration.num <- length(sample.name)
    load.data <- readChar(file.read, file.info(file.read)$size)
    # overall alignment rate
    overall.alignment <- strsplit(load.data, "\n")
    overall.alignment.with.NA <- stringr::str_extract(overall.alignment[[1]], "[0-9]*.[0-9]*% overall alignment rate")
    overall.alignment.result <- overall.alignment.with.NA[!is.na(overall.alignment.with.NA)]
    overall.alignment.result.cut <- gsub(" overall alignment rate", " ", overall.alignment.result)
    # different mapping rate
    first.split <- strsplit(load.data, "\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* Hisat2 Aligning \\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\n")
    second.split <- strsplit(first.split[[1]][2], "\n\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\* Current progress of RNA-seq files in")
    split.lines <- strsplit(second.split[[1]][1], "\n")
    alignment.rate.with.NA <- stringr::str_extract(split.lines[[1]], "[0-9]* \\([0-9]*.[0-9]*%\\) aligned concordantly")
    alignment.first.result <- alignment.rate.with.NA[!is.na(alignment.rate.with.NA)]
    alignment.first.result.cut1 <- gsub(") aligned concordantly", " ", alignment.first.result)
    alignment.first.result.cut2 <- gsub("[0-9]* \\(", " ", alignment.first.result.cut1)
    report.data.frame <- data.frame(matrix(0, ncol = 0, nrow = 3))
    row.names(report.data.frame) <- c("Unique mapping rate", "Multiple mapping rate", "Overall alignment rate")
    for( i in 1:iteration.num){
      add.column <- c()
      for( j in (i*3-1):(i*3)){
        add.column <- c(add.column, alignment.first.result.cut2[j])
      }
      add.column <- c(add.column, overall.alignment.result.cut[i])
      report.data.frame[[(sample.name[i])]] <- add.column
    }
    dir.create(paste0(path.prefix, "RNAseq_results/Alignment_Report/"))
    write.csv(report.data.frame, file = paste0(path.prefix, "RNAseq_results/Alignment_Report/Alignment_report.csv"))
    png(paste0(path.prefix, "RNAseq_results/Alignment_Report/Alignment_report.png"), width = iteration.num*100 + 200, height = 40*4)
    p <- gridExtra::grid.table(report.data.frame)
    print(p)
    dev.off()
    cat(c("Results are in", paste0("'", path.prefix, "RNAseq_results/Alignment_Report/'"), "\n\n"))
  }
}

#' use 'samtools' to sort and convert the SAM files to BAM
SamtoolsToBam <- function(path.prefix, gene.name, sample.pattern, num.parallel.threads = 8) {
  if (isTRUE(CheckSamtools(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix, gene.name, sample.pattern, print=TRUE)
    cat(paste0("\n************** Samtools converting '.sam' to '.bam' **************\n"))
    if (check.results$sam.files.number.df != 0){
      # Map reads to each alignment
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      sample.table <- table(gsub(paste0(".sam$"), replace = "", check.results$sam.files.df))
      iteration.num <- length(sample.table)
      sample.name <- names(sample.table)
      sample.value <- as.vector(sample.table)
      for( i in 1:iteration.num){
        whole.command <- paste("sort -@", num.parallel.threads, "-o", paste0("raw_bam/", sample.name[i], ".bam"), paste0("raw_sam/", sample.name[i], ".sam"))
        if (i != 1) cat("\n")
        cat(c("Input command :", paste("samtools", whole.command), "\n"))
        system2(command = "samtools", args = whole.command)
      }
      cat("\n")
      on.exit(setwd(current.path))
    } else {
      stop(c("(\u2718) 'XXX.sam' is missing.\n\n"))
    }
  }
}

#' stringtie assemble and quantify expressed genes and transcripts
StringTieAssemble <- function(path.prefix, gene.name, sample.pattern, num.parallel.threads = 8) {
  if (isTRUE(CheckStringTie(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix, gene.name, sample.pattern, print=TRUE)
    cat(paste0("\n************** Stringtie assembling **************\n"))
    if (check.results$bam.files.number.df != 0 && isTRUE(check.results$gtf.file.logic.df)){
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      sample.name <- sort(gsub(paste0(".bam$"), replace = "", check.results$bam.files.df))
      iteration.num <- length(sample.name)
      for( i in 1:iteration.num){
        whole.command <- paste("-p", num.parallel.threads, "-G",paste0("ref_genes/", gene.name, ".gtf"), "-A", paste0("gene_abundance/", sample.name[i],"/", sample.name[i], ".tsv"), "-o", paste0("raw_gtf/", sample.name[i], ".gtf"), "-l", sample.name[i], paste0("raw_bam/", sample.name[i], ".bam"))
        if (i != 1) cat("\n")
        cat(c("Input command :", paste("stringtie", whole.command), "\n"))
        system2(command = "stringtie", args = whole.command)
      }
      cat("\n")
      on.exit(setwd(current.path))
    } else {
      stop(c(paste0("(\u2718) '", gene.name, ".gtf' "), "or 'XXX.bam' is missing.\n\n"))
    }
  }
}

#' stringtie merge transcripts from all samples
StringTieMergeTrans <- function(path.prefix, gene.name, sample.pattern, num.parallel.threads = 8) {
  if (isTRUE(CheckStringTie(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix, gene.name, sample.pattern, print=TRUE)
    cat(paste0("\n************** Stringtie merging transcripts **************\n"))
    if ( isTRUE(check.results$gtf.file.logic.df) && check.results$gtf.files.number.df != 0){
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      dir.create(file.path(paste0(path.prefix, 'gene_data/merged/')), showWarnings = FALSE)
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
      whole.command <- paste("--merge -p", num.parallel.threads, "-G", paste0("ref_genes/", gene.name, ".gtf"), "-o", "merged/stringtie_merged.gtf", "merged/mergelist.txt")
      cat(c("Input command :", paste("stringtie", whole.command), "\n"))
      system2(command = "stringtie", args = whole.command)
      cat("\n")
      on.exit(setwd(current.path))
    } else {
      stop(c(paste0("(\u2718) :'", gene.name, ".gtf' "), "or", " 'XXX.gtf' is missing.\n\n"))
    }
  }
}

#' stringtie estimate transcript abundances and create table counts for Ballgown
StringTieToBallgown <- function(path.prefix, gene.name, sample.pattern, num.parallel.threads = 8) {
  if (isTRUE(CheckStringTie(print=FALSE))) {
    check.results <- ProgressGenesFiles(path.prefix, gene.name, sample.pattern, print=TRUE)
    cat(paste0("\n************** Stringtie creating table counts for Ballgown **************\n"))
    if ((check.results$bam.files.number.df != 0) && isTRUE(check.results$stringtie_merged.gtf.file.df)){
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      sample.table <- table(gsub(paste0(".bam$"), replace = "", check.results$bam.files.df))
      iteration.num <- length(sample.table)
      sample.name <- names(sample.table)
      sample.value <- as.vector(sample.table)
      for( i in 1:iteration.num){
        # '-e' only estimate the abundance of given reference transcripts (requires -G)
        whole.command <- paste("-e -B -p", num.parallel.threads, "-G", "merged/stringtie_merged.gtf", "-o", paste0("ballgown/", sample.name[i],"/", sample.name[i], ".gtf"), paste0("raw_bam/", sample.name[i], ".bam"))
        if (i != 1) cat("\n")
        cat(c("Input command :", paste("stringtie", whole.command), "\n"))
        system2(command = "stringtie", args = whole.command)
      }
      cat("\n")
      on.exit(setwd(current.path))
    } else {
      stop(c(paste0("(\u2718) 'stringtie_merged.gtf' "), "or", " 'XXX.bam' is missing.\n\n"))
    }
  }
}

#' Examine how the transcripts compare with the reference annotation
GffcompareRefSample <- function(path.prefix, gene.name, sample.pattern) {
  if (isTRUE(CheckGffcompare(print=FALSE))){
    check.results <- ProgressGenesFiles(path.prefix, gene.name, sample.pattern, print=TRUE)
    cat(paste0("\n************** Gffcompare comparing transcripts between merged and reference **************\n"))
    if ( isTRUE(check.results$stringtie_merged.gtf.file.df) && isTRUE(check.results$gtf.file.logic.df)){
      current.path <- getwd()
      setwd(paste0(path.prefix, "gene_data/"))
      whole.command <- paste("-r", paste0("ref_genes/", gene.name, ".gtf"), "-G -o merged/merged merged/stringtie_merged.gtf")
      cat(c("Input command :", paste("gffcompare", whole.command), "\n"))
      system2(command = "gffcompare", args = whole.command)
      cat("\n")
      on.exit(setwd(current.path))
    } else {
      stop(c(paste0("(\u2718) '", gene.name, ".gtf' "), "or", paste0(" 'stringtie_merged.gtf'"), "is missing.\n\n"))
    }
  }
}
