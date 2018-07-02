#' inner function : Copy input files directory
#' @export
CopyInputDir <- function(path.prefix, input.path.prefix, gene.name, sample.pattern, optional = indexes.optional) {
  current.path <- getwd()
  setwd(paste0(path.prefix, "gene_data/"))
  cat(c("************** Copying", paste0("'", input.path.prefix, "input_files/'"), "************\n"))
  cat(c("Copying From :", paste0(input.path.prefix, "input_files/", gene.name, ".gtf"), "\n"))
  file.copy(paste0(input.path.prefix, "input_files/", gene.name, ".gtf"), paste0(getwd(), "/ref_genes/", gene.name, ".gtf"))
  cat(c("          To :"), paste0(getwd(), "/ref_genes/", gene.name, ".gtf", "\n"))
  cat(c("Copying From :", paste0(input.path.prefix, "input_files/", gene.name, ".fa"),  "\n"))
  file.copy(paste0(input.path.prefix, "input_files/", gene.name, ".fa"), paste0(getwd(), "/ref_genome/", gene.name, ".fa"))
  cat(c("          To :"), paste0(getwd(), "/ref_genome/", gene.name, ".fa", "\n"))
  cat(c("Copying From :", paste0(input.path.prefix, "input_files/", "raw_fastq.gz/"),  "\n"))
  file.copy(paste0(input.path.prefix, "input_files/", "raw_fastq.gz/"), paste0(getwd(), "/"), overwrite = TRUE, recursive = TRUE)
  cat(c("          To :", paste0(getwd(),"/raw_fastq.gz/"), "\n"))
  cat(c("Copying From :", paste0(input.path.prefix, "input_files/phenodata.csv"),  "\n"))
  file.copy(paste0(input.path.prefix, "input_files/phenodata.csv"), paste0(getwd(), "/phenodata.csv"))
  cat(c("          To :"), paste0(getwd(), "/phenodata.csv\n"))
  on.exit(setwd(current.path))
  if (isTRUE(optional)) {
    cat(c("Copying From :", paste0(input.path.prefix, "input_files/", "indexes/"),  "\n"))
    file.copy(paste0(input.path.prefix, "input_files/", "indexes/"), paste0(getwd(), "/"), overwrite = TRUE, recursive = TRUE)
    cat(c("          To :", paste0(getwd(),"/indexes/"), "\n\n"))
  } else {
    cat("\n")
  }
}

#' Create sample gene directory
MkdirGeneDir <- function(path.prefix) {
  cat("************** Creating 'gene-data/' directory ************\n")
  gene_data.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/')), showWarnings = FALSE) == 0
  if (!isTRUE(gene_data.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "gene_data/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "gene_data/'.\n     Please check whether the directory is already exit.\n"))
  }
  ref.genome.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/ref_genome/')), showWarnings = FALSE) == 0
  if (!isTRUE(ref.genome.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "gene_data/ref_genome/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "gene_data/ref_genome/'.\n     Please check whether the directory is already exit.\n"))
  }
  ref.genes.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/ref_genes/')), showWarnings = FALSE) == 0
  if (!isTRUE(ref.genes.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "gene_data/ref_genes/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "gene_data/ref_genes/'.\n     Please check whether the directory is already exit.\n"))
  }
  indexes.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/indexes/')), showWarnings = FALSE) == 0
  if (!isTRUE(indexes.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "gene_data/indexes/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "gene_data/indexes/'.\n     Please check whether the directory is already exit.\n"))
  }
  samples.fastq.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/raw_fastq.gz')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.fastq.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "gene_data/raw_fastq.gz/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "gene_data/raw_fastq.gz/'.\n     Please check whether the directory is already exit.\n"))
  }
  samples.sam.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/raw_sam')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.sam.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "gene_data/raw_sam/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "gene_data/raw_sam/'.\n     Please check whether the directory is already exit.\n"))
  }
  samples.bam.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/raw_bam')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.bam.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "gene_data/raw_bam/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "gene_data/raw_bam/'.\n     Please check whether the directory is already exit.\n"))
  }
  samples.gtf.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/raw_gtf')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.gtf.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "gene_data/raw_gtf/'.\n\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "gene_data/raw_gtf/'.\n     Please check whether the directory is already exit.\n\n"))
  }
}

#' Make RNAseq_bin/ directory
MkdirRNAseq_bin <- function(path.prefix) {
  cat("************** Creating 'RNAseq_bin/' directory ************\n")
  RNAseq_bin.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_bin/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_bin.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "RNAseq_bin/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "RNAseq_bin/'.\n     Please check whether the directory is already exit.\n"))
  }
  download.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_bin/Download/')), showWarnings = FALSE) == 0
  if (!isTRUE(download.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "RNAseq_bin/Download/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "RNAseq_bin/Download/'.\n     Please check whether the directory is already exit.\n"))
  }
  unpacked.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_bin/Unpacked/')), showWarnings = FALSE) == 0
  if (!isTRUE(unpacked.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "RNAseq_bin/Unpacked/'.\n\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "RNAseq_bin/Unpacked/'.\n     Please check whether the directory is already exit.\n\n"))
  }
}

#' Make RNAseq_results/ directory
MkdirRNAseq_results <- function(path.prefix){
  cat("************** Creating 'RNAseq_results/' directory ************\n")
  RNAseq_results.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_results/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_results.dir)) {
    cat(paste0("(\u2714) : Create '", path.prefix, "RNAseq_results/'.\n"))
  } else {
    cat(paste0("(\u2718) : Fail to create '", path.prefix, "RNAseq_results/'.\n     Please check whether the directory is already exit.\n"))
  }
}

#' Check sample gene and binary directory
#' @export
CheckDirAll <- function(path.prefix, print = TRUE) {
  if (print) {
    cat(c("************** Checking 'gene_data/' and 'RNAseq_bin/' directories ************\n"))
  }
  gene.data.logic <- dir.exists(paste0(path.prefix, "gene_data/ref_genes/")) &&
    dir.exists(paste0(path.prefix, "gene_data/ref_genome/")) &&
    dir.exists(paste0(path.prefix, "gene_data/indexes/")) &&
    dir.exists(paste0(path.prefix, "gene_data/raw_fastq.gz/")) &&
    dir.exists(paste0(path.prefix, "gene_data/raw_sam/")) &&
    dir.exists(paste0(path.prefix, "gene_data/raw_bam/")) &&
    dir.exists(paste0(path.prefix, "gene_data/raw_gtf/"))
  rnaseq.bin.logic <- dir.exists(paste0(path.prefix, "RNAseq_bin/Download/")) && dir.exists(paste0(path.prefix, "RNAseq_bin/Unpacked/"))
  if (!isTRUE(gene.data.logic)){
    cat(c("There are directories missing in", paste0(path.prefix, "gene_data/.")), "\n")
  }
  if (!isTRUE(rnaseq.bin.logic)){
    cat(c("There are directories missing in", paste0(path.prefix, "RNAseq_bin/.")), "\n")
  }
  if (isTRUE(gene.data.logic) && isTRUE(rnaseq.bin.logic)){
    if (print) {
      cat(c("(\u2714) : Directories are all correct.\n\n"))
    }
    return(TRUE)
  } else {
    cat(c("(\u2718) : Please run 'MkdirAll(path.prefix)' to add the missing directories.\n\n"))
    return(FALSE)
  }
}

#' Create sample gene and binary directory
#' @export
MkdirAll <- function(path.prefix) {
  MkdirGeneDir(path.prefix)
  MkdirRNAseq_bin(path.prefix)
  MkdirRNAseq_results(path.prefix)
}

#' check 'gene_data' and subdirectory files exit
#' @export
ProgressGenesFiles <- function(gene.name, sample.pattern, print = TRUE) {
  if (print) {
    cat(paste0("************** Current progress of RNA-seq files in '", paste0(path.prefix, "gene_data/'"), " **************\n"))
  }
  # 1. Check .gtf file
  gtf.file <- file.exists(paste0(path.prefix, "gene_data", '/ref_genes/', gene.name, '.gtf'))
  if (isTRUE(gtf.file)) {
    if(print){
      cat(c("(\u2714) :", paste0("'",path.prefix, "gene_data", '/ref_genes/', gene.name, '.gtf', "'"), "is exit\n\n"))
    }
  } else {
    cat(c("(\u2718) :", paste0("'",path.prefix, "gene_data", '/ref_genes/', gene.name, '.gtf', "'"), "is not exit\n"))
    cat(c("     Put the", paste0("'",gene.name,".gtf", "'"), "file in", paste0("'",path.prefix, "gene_data", '/ref_genes/', "'"), "to fix the error.\n\n"))
  }
  # 2. Check .fa file
  fa.file <- file.exists(paste0(path.prefix, "gene_data", '/ref_genome/', gene.name, '.fa'))
  if (isTRUE(fa.file)) {
    if(print){
      cat(c("(\u2714) :",paste0("'", path.prefix, "gene_data", '/ref_genome/', gene.name, '.fa', "'"), "is exit\n\n"))
    }
  } else {
    cat(c("(\u2718) :",paste0("'", path.prefix, "gene_data", '/ref_genome/', gene.name, '.fa', "'"), "is not exit\n"))
    cat(c("     Put the", paste0("'",gene.name,".fa", "'"), "file in", paste0("'",path.prefix, "gene_data", '/ref_genome/', "'"), "to fix the error.\n\n"))
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

  pheno_data <- read.csv(paste0(path.prefix, "gene_data", "/phenodata.csv"))
  sample.table <- as.data.frame(table(pheno_data[1]))
  compare.sample.table <- as.data.frame(table(pheno_data[2]))
  extract.fastq.gz.sample.names <- unique(gsub("_[1-2]*.fastq.gz", "", fastq.gz.files))
  if(length(extract.fastq.gz.sample.names) == length(row.names(sample.table))){
    # check the number of sample(rows) in 'phenodata.csv' is same as the number of unique samples in 'fastq.gz.files'
    if(identical(sort(extract.fastq.gz.sample.names), sort(as.character(sample.table[1][,1])))){
      # check all the name in bath 'phenodata.csv' and 'fastq.gz.files' are same
      if (length(row.names(compare.sample.table)) == 2) {
        ## only two groups are allowed
        invalid.column.number <- 0
        adjustvars <- c()
        for( i in 1:length(pheno_data) ){
          if(length(pheno_data[i][is.na(pheno_data[i]) == TRUE]) != 0) {
            invalid.column.number <- invalid.column.number + 1
            cat(paste0("(\u2718) : There are missing values in column '", names(pheno_data[i]), "'\n" ))
          }
          if(i > 2){
            adjustvars <- append(adjustvars, names(pheno_data[i]))
          }
        }
        if (invalid.column.number != 0) {
          cat(paste0("(\u2718) : ", invalid.column.number, " columns are invalid. Please fix 'phenodata.csv'\n\n" ))
          adjustvars <- c()
        }
      }
    } else {
      cat(paste0("(\u2718) : Sample names in 'phenodata.csv' and 'XXX.fastq.gz' are different. Please check 'phenodata.csv' matches the 'XXX.fastq.gz'\n\n" ))
    }
  } else {
    cat(paste0("(\u2718) : The sample's numbers in 'phenodata.csv' and 'XXX.fastq.gz' are different. Please check 'phenodata.csv' matches the 'XXX.fastq.gz'\n\n" ))
  }
  if (isTRUE(gtf.file)) {
    if(print){
      cat(c("(\u2714) :", paste0("'",path.prefix, "gene_data", '/phenodata.csv', "'"), "is exit\n\n"))
    }
  } else {
    cat(c("(\u2718) :", paste0("'",path.prefix, "gene_data", '/phenodata.csv', "'"), "is not exit\n"))
    cat(c("     Put the", paste0("'phenodata.csv'"), "file in", paste0("'",path.prefix, "gene_data", '/', "'"), "to fix the error.\n\n"))
  }
  ht2.files <- list.files(path = paste0(path.prefix, "gene_data", '/indexes/'), pattern = paste0("^", gene.name, "_tran.[0-9]*.ht2$"), all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  ht2.files.number <- length(ht2.files)
  if (ht2.files.number != 0) {
    if (print) {
      for (i in ht2.files){
        cat(c("(\u2714) :",paste0("'",path.prefix, "gene_data", '/indexes/', i, "'"), "is exit\n"))
      }
      cat(c("Total:", ht2.files.number, "files\n\n"))
    }
  } else {
    cat(c("(\u231B) :", paste0('\'',path.prefix, "gene_data", '/indexes/', gene.name, '_tran.*.ht2\''), "is not exit\n"))
    cat("       Files haven't created yet. Run 'CreateHisat2Index()' to generate '*.ht2' files or download from 'https://ccb.jhu.edu/software/hisat2/index.shtml'\n\n")
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
    cat(c("(\u231B) :", paste0('\'', path.prefix, "gene_data", '/raw_sam/XXX.sam\''), "is not exit\n"))
    cat("       Files haven't created yet. Run 'Hisat2AlignmentDefault()' to generate 'XXX.sam' files\n\n")
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
    cat(c("(\u231B) :", paste0('\'', path.prefix, "gene_data", '/raw_sam/XXX.bam\''), "is not exit\n"))
    cat("       Files haven't created yet. Run 'SamtoolsToBam()' to generate 'XXX.bam' files\n\n")
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
    cat(c("(\u231B) :", paste0('\'',path.prefix, "gene_data", '/raw_gtf/XXX.gtf\''), "is not exit\n"))
    cat("       Files haven't created yet. Run 'StringTieAssemble()' to generate 'XXX.gtf' files\n\n")
  }
  stringtie_merged.gtf.file <- file.exists(paste0(path.prefix, "gene_data", '/merged/stringtie_merged.gtf'))
  if (isTRUE(stringtie_merged.gtf.file)) {
    if (print) {
      cat(c("(\u2714) :", paste0("'", path.prefix, "gene_data", "/merged/stringtie_merged.gtf", "'"), "is exit\n\n"))
    }
  } else {
    cat(c("(\u231B) :", paste0("'", path.prefix, "gene_data", "/merged/stringtie_merged.gtf", "'"), "is not exit\n"))
    cat("       Files haven't created yet. Run 'StringTieMergeTrans()' to generate 'stringtie_merged.gtf' files\n\n")
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
    cat(c("(\u231B) :", paste0('\'', path.prefix, "gene_data", '/merged/', "merged.", "XXX/"), "is not exit\n"))
    cat("       Directories haven't created yet. Run 'GffcompareRefSample()' to generate", paste0("merged/", "merged.XXX/"), "files\n\n")
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
    cat(c("(\u231B) :", paste0('\'',path.prefix, "gene_data", '/ballgown/', gsub(".fastq.gz", replace = "", sample.pattern), "/"), "is not exit\n"))
    cat("       Directories haven't created yet. Run 'StringTieToBallgown()' to generate", paste0("ballgown/", sample.pattern, "/"), "directories\n\n")
  }
  return(list(gtf.file.logic.df = gtf.file, fa.file.logic.df = fa.file,
              fastq.gz.files.number.df = fastq.gz.files.number,
              fastq.gz.files.df = fastq.gz.files,
              phenodata.file.df = phenodata.file,
              phenodata.invalid.column.number.df = invalid.column.number,
              phenodata.adjustvars.df = adjustvars,
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

#' Add '~/RNAseq_bin/ to R environment "PATH"
#' @export
ExportPath <- function(path.prefix) {
  cat("************** Adding PATH to R environment ************\n")
  old.path <- Sys.getenv("PATH")
  Sys.setenv(
    PATH = paste(old.path, paste0(path.prefix, "RNAseq_bin"), sep = ":")
  )
  cat("\u27a4\u27a4 R environment 'PATH' : ", Sys.getenv("PATH"), "\n\n")
}

