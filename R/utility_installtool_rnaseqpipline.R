# check 'gene_data' and subdirectory files exit
ProgressGenesFiles <- function(path.prefix, genome.name, sample.pattern, print = TRUE) {
  if (print) {
    cat(paste0("************** Current progress of RNA-seq files in '", paste0(path.prefix, "gene_data/'"), " **************\n"))
  }
  # 1. Check .gtf file
  gtf.file <- file.exists(paste0(path.prefix, "gene_data", '/ref_genes/', genome.name, '.gtf'))
  if (isTRUE(gtf.file)) {
    if(print){
      cat(c("(\u2714) :", paste0("'",path.prefix, "gene_data", '/ref_genes/', genome.name, '.gtf', "'"), "is exit\n\n"))
    }
  } else {
    cat(c("(\u2718) :", paste0("'",path.prefix, "gene_data", '/ref_genes/', genome.name, '.gtf', "'"), "is not exit\n"))
    cat(c("     Put the", paste0("'",genome.name,".gtf", "'"), "file in", paste0("'",path.prefix, "gene_data", '/ref_genes/', "'"), "to fix the error.\n\n"))
  }
  # 2. Check .fa file
  fa.file <- file.exists(paste0(path.prefix, "gene_data", '/ref_genome/', genome.name, '.fa'))
  if (isTRUE(fa.file)) {
    if(print){
      cat(c("(\u2714) :",paste0("'", path.prefix, "gene_data", '/ref_genome/', genome.name, '.fa', "'"), "is exit\n\n"))
    }
  } else {
    cat(c("(\u2718) :",paste0("'", path.prefix, "gene_data", '/ref_genome/', genome.name, '.fa', "'"), "is not exit\n"))
    cat(c("     Put the", paste0("'",genome.name,".fa", "'"), "file in", paste0("'",path.prefix, "gene_data", '/ref_genome/', "'"), "to fix the error.\n\n"))
  }
  # 3. Check .fastq.gz file
  fastq.gz.files <- list.files(path = paste0(path.prefix, "gene_data", '/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  fastq.gz.files.number <- length(fastq.gz.files)
  if (fastq.gz.files.number != 0){
    if(print){
      for (i in fastq.gz.files){
        cat(c("(\u2714) :", paste0("'", path.prefix, "gene_data", '/raw_fastq.gz/', i, "'"), "is exit\n"))
      }
      cat(c("Total:", fastq.gz.files.number, "files\n\n"))
    }
  }else {
    cat(c("(\u2718) :", paste0('\'',path.prefix, "gene_data", '/raw_fastq.gz/XXX_*.fastq.gz\''), "is not exit\n"))
    cat(c("     Put the", paste0('XXX_*.fastq.gz'), "file in", paste0("'",path.prefix, "gene_data", '/raw_fastq.gz/', "'"), "to fix the error.\n\n"))
  }
  # 4. Check phenodata file
  phenodata.file <- file.exists(paste0(path.prefix, "gene_data", '/phenodata.csv'))
  if (isTRUE(phenodata.file)) {
    if(print){
      cat(c("(\u2714) :", paste0("'",path.prefix, "gene_data", '/phenodata.csv', "'"), "is exit\n\n"))
    }
  } else {
    cat(c("(\u2718) :", paste0("'",path.prefix, "gene_data", '/phenodata.csv', "'"), "is not exit\n"))
    cat(c("     Put the", paste0("'phenodata.csv'"), "file in", paste0("'",path.prefix, "gene_data", '/', "'"), "to fix the error.\n\n"))
  }
  ht2.files <- list.files(path = paste0(path.prefix, "gene_data", '/indexes/'), pattern = paste0("^", genome.name, "_tran.[0-9]*.ht2$"), all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  ht2.files.number <- length(ht2.files)
  if (ht2.files.number != 0) {
    if (print) {
      for (i in ht2.files){
        cat(c("(\u2714) :",paste0("'",path.prefix, "gene_data", '/indexes/', i, "'"), "is exit\n"))
      }
      cat(c("Total:", ht2.files.number, "files\n\n"))
    }
  } else {
    if (print) {
      cat(c("(\u231B) :", paste0('\'',path.prefix, "gene_data", '/indexes/', genome.name, '_tran.*.ht2\''), "is not exit\n"))
    }
  }
  sam.files <- list.files(path = paste0(path.prefix, "gene_data", '/raw_sam/'), pattern = paste0( "^[A-Z, a-z]*", "[0-9]*", "[A-Z, a-z]*", ".sam$"), all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  sam.files.number <- length(sam.files)
  if (sam.files.number != 0){
    if (print) {
      for (i in sam.files){
        cat(c("(\u2714) :", paste0("'", path.prefix, "gene_data", '/raw_sam/', i, "'"), "is exit\n"))
      }
      cat(c("Total:", sam.files.number, "files\n\n"))
    }
  }else {
    if (print) {
      cat(c("(\u231B) :", paste0('\'', path.prefix, "gene_data", '/raw_sam/XXX.sam\''), "is not exit\n"))
    }
  }
  bam.files <- list.files(path = paste0(path.prefix, "gene_data", '/raw_bam/'), pattern = paste0( "^[A-Z, a-z]*", "[0-9]*", "[A-Z, a-z]*", ".bam$"), all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  bam.files.number <- length(bam.files)
  if (bam.files.number != 0){
    if (print) {
      for (i in bam.files){
        cat(c("(\u2714) :", paste0("'", path.prefix, "gene_data", '/raw_bam/', i, "'"), "is exit\n"))
      }
      cat(c("Total:", bam.files.number, "files\n\n"))
    }
  }else {
    if (print) {
      cat(c("(\u231B) :", paste0('\'', path.prefix, "gene_data", '/raw_sam/XXX.bam\''), "is not exit\n"))
    }
  }
  gtf.files <- list.files(path = paste0(path.prefix, "gene_data", '/raw_gtf/'), pattern = paste0("^[A-Z, a-z]*", "[0-9]*", "[A-Z, a-z]*", ".gtf$"), all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  gtf.files.number <- length(gtf.files)
  if (gtf.files.number != 0){
    if (print) {
      for (i in gtf.files){
        cat(c("(\u2714) :", paste0("'", path.prefix, "gene_data", '/raw_gtf/', i, "'"), "is exit\n"))
      }
      cat(c("Total:", gtf.files.number, "files\n\n"))
    }
  }else {
    if (print) {
      cat(c("(\u231B) :", paste0('\'',path.prefix, "gene_data", '/raw_gtf/XXX.gtf\''), "is not exit\n"))
    }
  }
  stringtie_merged.gtf.file <- file.exists(paste0(path.prefix, "gene_data", '/merged/stringtie_merged.gtf'))
  if (isTRUE(stringtie_merged.gtf.file)) {
    if (print) {
      cat(c("(\u2714) :", paste0("'", path.prefix, "gene_data", "/merged/stringtie_merged.gtf", "'"), "is exit\n\n"))
    }
  } else {
    if (print) {
      cat(c("(\u231B) :", paste0("'", path.prefix, "gene_data", "/merged/stringtie_merged.gtf", "'"), "is not exit\n"))
    }
  }
  gffcompare.related.dirs <- list.files(path = paste0(path.prefix, "gene_data", '/merged/'), pattern = "^merged.", all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  gffcompare.related.dirs.number <- length(gffcompare.related.dirs)
  if (gffcompare.related.dirs.number != 0){
    if (print) {
      for (i in gffcompare.related.dirs){
        cat(c("(\u2714) :", paste0("'", path.prefix, "gene_data", '/merged/', i, "'"), "is exit\n"))
      }
      cat(c("Total:", gffcompare.related.dirs.number, "files\n\n"))
    }
  }else {
    if (print) {
      cat(c("(\u231B) :", paste0('\'', path.prefix, "gene_data", '/merged/', "merged.", "XXX/"), "is not exit\n"))
    }
  }
  ballgown.dirs <- list.files(path = paste0(path.prefix, "gene_data", '/ballgown/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  ballgown.dirs.number <- length(ballgown.dirs)
  if (ballgown.dirs.number != 0){
    if (print) {
      for (i in ballgown.dirs){
        cat(c("(\u2714) :", paste0("'", path.prefix, "gene_data", '/ballgown/', i, "'"), "is exit\n"))
      }
      cat(c("Total:", ballgown.dirs.number, "directories\n\n"))
    }
  }else {
    if (print) {
      cat(c("(\u231B) :", paste0('\'',path.prefix, "gene_data", '/ballgown/', gsub(".fastq.gz", replacement = "", sample.pattern), "/"), "is not exit\n"))
    }
  }
  return(list(gtf.file.logic.df = gtf.file, fa.file.logic.df = fa.file,
              fastq.gz.files.number.df = fastq.gz.files.number,
              fastq.gz.files.df = fastq.gz.files,
              phenodata.file.df = phenodata.file,
              ht2.files.number.df = ht2.files.number,
              ht2.files.df = ht2.files,
              sam.files.number.df = sam.files.number,
              sam.files.df = sam.files,
              bam.files.number.df = bam.files.number,
              bam.files.df = bam.files,
              gtf.files.number.df = gtf.files.number,
              gtf.files.df = gtf.files,
              stringtie_merged.gtf.file.df = stringtie_merged.gtf.file,
              gffcompare.related.dirs.df = gffcompare.related.dirs,
              gffcompare.related.dirs.number.df = gffcompare.related.dirs.number,
              ballgown.dirs.number.df = ballgown.dirs.number,
              ballgown.dirs.df = ballgown.dirs))
}

# Add '~/RNAseq_bin/ to R environment "PATH"
ExportPath <- function(path.prefix) {
  cat("************** Adding PATH to R environment ************\n")
  old.path <- Sys.getenv("PATH")
  Sys.setenv(
    PATH = paste(old.path, paste0(path.prefix, "RNAseq_bin"), sep = ":")
  )
  cat("\u27a4\u27a4 R environment 'PATH' : ", Sys.getenv("PATH"), "\n\n")
}

ParseFPKMBallgownResult <- function(path.prefix, independent.variable, control.group, experiment.group, file) {
  data.read.csv <- read.csv(file = file)
  pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
  sample.table <- as.data.frame(table(pheno_data[independent.variable]))
  control.group.size <- sample.table[sample.table$Var1 == control.group,]$Freq
  experiment.group.size <- sample.table[sample.table$Var1 == experiment.group,]$Freq
  # For control group
  control.group.range <- 5:(5+control.group.size-1)
  # For cexperiment group
  experiment.group.range <- (5+control.group.size+1):(5+control.group.size+experiment.group.size)
  return.data.frame.index <- c(control.group.range, experiment.group.range)
  return.data <- data.read.csv[return.data.frame.index]
  return(return.data)
}

FindControlExperiment <- function(path.prefix, independent.variable, control.group, experiment.group) {
  pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
  # gene count table
  gene.count.table <- read.csv(paste0(path.prefix, "gene_data/reads_count_matrix/gene_count_matrix.csv"))
  # transcript count table
  # transcript.count.table <- read.csv(paste0(path.prefix, "gene_data/reads_count_matrix/transcript_count_matrix.csv"))
  control.group.data.frame <- pheno_data[pheno_data[independent.variable] == control.group, ]
  experiment.group.data.frame <- pheno_data[pheno_data[independent.variable] == experiment.group, ]
  return(list("control.group" = control.group.data.frame, "experiment.group" = experiment.group.data.frame))
}
