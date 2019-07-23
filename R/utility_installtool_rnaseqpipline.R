# Querying Web Resource template
getURL <- function(URL, FUN, ..., N.TRIES=1L) {
  N.TRIES <- as.integer(N.TRIES)
  stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

  while (N.TRIES > 0L) {
    result <- tryCatch(FUN(URL, ...), error=identity)
    if (!inherits(result, "error"))
      break
    N.TRIES <- N.TRIES - 1L
  }

  if (N.TRIES == 0L) {
    stop("'getURL()' failed:",
         "\n  URL: ", URL,
         "\n  error: ", conditionMessage(result))
  }
  result
}

# check 'gene_data' and subdirectory files exit
ProgressGenesFiles <- function(path.prefix,
                               genome.name,
                               sample.pattern,
                               print = TRUE) {
  if (print) {
    message("\uD83D\uDD0D\uD83D\uDD0D\uD83D\uDD0D Checking Files ...\n")
    message("************** Current progress of RNA-seq files in '",
            path.prefix, "gene_data/' **************\n")
  }
  # 1. Check .gtf file
  gtf.file <- file.exists(paste0(path.prefix, "gene_data",
                                 '/ref_genes/', genome.name, '.gtf'))
  if (isTRUE(gtf.file)) {
    if(print){
      message("(\u2714) :'",path.prefix, "gene_data", '/ref_genes/',
              genome.name, '.gtf', "' is exit\n\n")
    }
  } else {
    message("(\u2718) :'",path.prefix, "gene_data/ref_genes/",
            genome.name, ".gtf' is not exit\n")
    message("     Put the '",genome.name, ".gtf' file in '", path.prefix,
            "gene_data/ref_genes/' to fix the error.\n\n")
  }
  # 2. Check .fa file
  fa.file <- file.exists(paste0(path.prefix, "gene_data",
                                '/ref_genome/', genome.name, '.fa'))
  if (isTRUE(fa.file)) {
    if(print){
      message("(\u2714) : '", path.prefix, "gene_data/ref_genome/",
              genome.name, ".fa' is exit\n\n")
    }
  } else {
    message("(\u2718) : '", path.prefix, "gene_data/ref_genome/",
            genome.name, ".fa' is not exit\n")
    message("     Put the '",genome.name,".fa' file in '",path.prefix,
            "gene_data/ref_genome/' to fix the error.\n\n")
  }
  # 3. Check .fastq.gz file
  fastq.gz.files <- list.files(path = paste0(path.prefix,
                                             "gene_data", "/raw_fastq.gz/"),
                               pattern = sample.pattern,
                               all.files = FALSE,
                               full.names = FALSE,
                               recursive = FALSE,
                               ignore.case = FALSE)
  fastq.gz.files.number <- length(fastq.gz.files)
  if (fastq.gz.files.number != 0){
    if(print){
      for (i in fastq.gz.files){
        message("(\u2714) : '", path.prefix,
                "gene_data/raw_fastq.gz/", i, "' is exit\n")
      }
      message("Total: ", fastq.gz.files.number, " file\n\n")
    }
  }else {
    message("(\u2718) : '", path.prefix,
            "gene_data/raw_fastq.gz/XXX_*.fastq.gz' is not exit\n")
    message("     Put the 'XXX_*.fastq.gz' file in '",
            path.prefix, "gene_data/raw_fastq.gz/' to fix the error.\n\n")
  }
  # 4. Check phenodata file
  phenodata.file <- file.exists(paste0(path.prefix,
                                       "gene_data", '/phenodata.csv'))
  if (phenodata.file) {
    if(print){
      message("(\u2714) : '", path.prefix,
              "gene_data/phenodata.csv' is exit\n\n")
    }
  } else {
    message("(\u2718) : '", path.prefix,
            "gene_data/phenodata.csv' is not exit\n")
    message("     Put the 'phenodata.csv' file in '",
            path.prefix, "gene_data/' to fix the error.\n\n")
  }
  ht2.files <- list.files(path = paste0(path.prefix, "gene_data", '/indices/'),
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
  ht2.files.number <- length(ht2.files)
  if (ht2.files.number != 0) {
    if (print) {
      for (i in ht2.files){
        message("(\u2714) : '", path.prefix,
                "gene_data/indices/", i, "' is exit\n")
      }
      message("Total: ", ht2.files.number, " file\n\n")
    }
  } else {
    if (print) {
      message("(\u231B) : '", path.prefix, "gene_data/indices/*' is not exit\n")
    }
  }
  # Target sam files !!
  target.sam.files <- unique(gsub("_[1,2].fastq.gz", ".sam", fastq.gz.files))
  # Check whether sam file is exist !
  actual.found.sam.files <- list.files(path = paste0(path.prefix,
                                                     "gene_data", '/raw_sam/'),
                                       pattern = "*.sam$",
                                       all.files = FALSE,
                                       full.names = FALSE,
                                       recursive = FALSE,
                                       ignore.case = FALSE)
  sam.files.number <- length(actual.found.sam.files)
  if (sam.files.number != 0){
    if (print) {
      sam.check <- identical(target.sam.files, actual.found.sam.files)
      if (sam.check) {
        for ( i in actual.found.sam.files ) {
          message("(\u2714) : '", path.prefix,
                  "gene_data/raw_sam/", i, "' is exit\n")
        }
        message("Total: ", sam.files.number, " file\n\n")
      } else {
        message("(\u2718) : .SAM files checking is invalid!! \n\n")
        stop("SAM files checking ERROR")
      }
    }
  } else {
    if (print) {
      message("(\u231B) : '", path.prefix,
              "gene_data/raw_sam/XXX.sam' is not exit\n")
    }
  }

  # Target bam file
  target.bam.files <- unique(gsub("_[1,2].fastq.gz", ".bam", fastq.gz.files))
  actual.found.bam.files <- list.files(path = paste0(path.prefix,
                                                     "gene_data", '/raw_bam/'),
                                       pattern = paste0("*.bam$"),
                                       all.files = FALSE,
                                       full.names = FALSE,
                                       recursive = FALSE,
                                       ignore.case = FALSE)
  bam.files.number <- length(actual.found.bam.files)
  if (bam.files.number != 0){
    if (print) {
      bam.check <- identical(target.bam.files, actual.found.bam.files)
      if (bam.check) {
        for (i in actual.found.bam.files){
          message("(\u2714) : '", path.prefix,
                  "gene_data/raw_bam/", i, "' is exit\n")
        }
        message("Total: ", bam.files.number, " files\n\n")
      } else {
        message("(\u2718) : .BAM files checking is invalid!! \n\n")
        stop("BAM files checking ERROR")
      }
    }
  }else {
    if (print) {
      message("(\u231B) : '", path.prefix,
              "gene_data/raw_bam/XXX.bam' is not exit\n")
    }
  }

  # Target gtf file
  target.gtf.files <- unique(gsub("_[1,2].fastq.gz", ".gtf", fastq.gz.files))
  actual.found.gtf.files <- list.files(path = paste0(path.prefix,
                                                     "gene_data", '/raw_gtf/'),
                                       pattern = paste0("*.gtf$"),
                                       all.files = FALSE,
                                       full.names = FALSE,
                                       recursive = FALSE,
                                       ignore.case = FALSE)
  gtf.files.number <- length(actual.found.gtf.files)
  if (gtf.files.number != 0){
    if (print) {
      gtf.check <- identical(target.gtf.files, actual.found.gtf.files)
      if (gtf.check) {
        for (i in actual.found.gtf.files){
          message("(\u2714) : '", path.prefix,
                  "gene_data/raw_gtf/", i, "' is exit\n")
        }
        message("Total: ", gtf.files.number, " files\n\n")
      } else {
        message("(\u2718) : .GTF files checking is invalid!! \n\n")
        stop("GTF files checking ERROR")
      }
    }
  }else {
    if (print) {
      message("(\u231B) : '", path.prefix,
              "gene_data/raw_gtf/XXX.gtf' is not exit\n")
    }
  }
  # Check stringtie merged gtf
  stringtie_merged.gtf.file <- file.exists(paste0(path.prefix,
                                                  "gene_data/merged/",
                                                  "stringtie_merged.gtf"))
  # Check gffcompare result
  gffcompare.related.dirs <- list.files(path = paste0(path.prefix,
                                                      "gene_data", '/merged/'),
                                        pattern = "^merged.",
                                        all.files = FALSE,
                                        full.names = FALSE,
                                        recursive = FALSE,
                                        ignore.case = FALSE)
  gffcompare.related.dirs.number <- length(gffcompare.related.dirs)
  if (gffcompare.related.dirs.number != 0){
    if (print) {
      for (i in gffcompare.related.dirs){
        message("(\u2714) : '", path.prefix,
                "gene_data/merged/", i, "' is exit\n")
      }
      message("Total: ", gffcompare.related.dirs.number, " files\n\n")
    }
  }else {
    if (print) {
      message("(\u231B) : '", path.prefix,
              "gene_data/merged/merged.XXX/' is not exit\n")
    }
  }
  ballgown.dirs <- list.files(path = paste0(path.prefix,
                                            "gene_data", '/ballgown/'),
                              pattern = sample.pattern,
                              all.files = FALSE,
                              full.names = FALSE,
                              recursive = FALSE,
                              ignore.case = FALSE)
  ballgown.dirs.number <- length(ballgown.dirs)
  if (ballgown.dirs.number != 0){
    if (print) {
      for (i in ballgown.dirs){
        message("(\u2714) : '", path.prefix,
                "gene_data/ballgown/", i, "' is exit\n")
      }
      message("Total: ", ballgown.dirs.number, " directories\n\n")
    }
  }else {
    if (print) {
      message("(\u231B) : '", path.prefix, "gene_data/ballgown/",
              gsub(".fastq.gz", replacement = "",
                   sample.pattern), "/' is not exit\n")
    }
  }
  return(list(gtf.file.logic.df = gtf.file, fa.file.logic.df = fa.file,
              fastq.gz.files.number.df = fastq.gz.files.number,
              fastq.gz.files.df = fastq.gz.files,
              phenodata.file.df = phenodata.file,
              ht2.files.number.df = ht2.files.number,
              ht2.files.df = ht2.files,
              sam.files.number.df = sam.files.number,
              sam.files.df = actual.found.sam.files,
              bam.files.number.df = bam.files.number,
              bam.files.df = actual.found.bam.files,
              gtf.files.number.df = gtf.files.number,
              gtf.files.df = actual.found.gtf.files,
              stringtie_merged.gtf.file.df = stringtie_merged.gtf.file,
              gffcompare.related.dirs.df = gffcompare.related.dirs,
              gffcompare.related.dirs.number.df =gffcompare.related.dirs.number,
              ballgown.dirs.number.df = ballgown.dirs.number,
              ballgown.dirs.df = ballgown.dirs))
}

# Add '~/RNASeq_bin/ to R environment "PATH"
ExportPath <- function(path.prefix) {
  message("************** Adding PATH to R environment ************\n")
  old.path <- Sys.getenv("PATH")
  Sys.setenv(
    PATH = paste(old.path, paste0(path.prefix, "RNASeq_bin"), sep = ":")
  )
  message("\u27a4\u27a4 R environment 'PATH' : ", Sys.getenv("PATH"), "\n\n")
}

ParseResultCSV <- function(which.analysis,
                           which.count.normalization,
                           path.prefix,
                           independent.variable,
                           case.group,
                           control.group) {
  case.normalized.csv <- paste0(path.prefix, "RNASeq_results/",
                                which.analysis, "/normalized_&_statistic/",
                                which.count.normalization, "_case.csv")
  control.normalized.csv <- paste0(path.prefix, "RNASeq_results/",
                                   which.analysis, "/normalized_&_statistic/",
                                   which.count.normalization, "_control.csv")
  statistic.csv <- paste0(path.prefix, "RNASeq_results/", which.analysis,
                          "/normalized_&_statistic/statistic.csv")
  read.case.normalized.csv <- read.csv(case.normalized.csv)
  read.control.normalized.csv <- read.csv(control.normalized.csv)
  read.statistic.csv <- read.csv(statistic.csv)
  return(list("case" = read.case.normalized.csv,
              "control" = read.control.normalized.csv,
              "statistic" = read.statistic.csv))
}



phenoDataWrap <- function(path.prefix,
                          independent.variable,
                          case.group,
                          control.group) {
  # Read pheno_data
  pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
  case.group.data.frame <-
    pheno_data[pheno_data[independent.variable] == case.group, ]
  control.group.data.frame <-
    pheno_data[pheno_data[independent.variable] == control.group, ]
  case.group.size <- length(row.names(case.group.data.frame))
  control.group.size <- length(row.names(control.group.data.frame))
  return(list("pheno_data" = pheno_data,
              "case.group.data.frame" = case.group.data.frame,
              "control.group.data.frame" = control.group.data.frame,
              "case.group.size" = case.group.size,
              "control.group.size" = control.group.size))
}

RawCountWrap <- function(path.prefix) {
  # read in gene count table
  gene.count.table.raw <- read.csv(paste0(path.prefix,
                                          "gene_data/reads_count_matrix/",
                                          "gene_count_matrix.csv"))
  # Process raw reads count!!
  gene.count.table.raw <- RawCountGeneNameChange(gene.count.table.raw,
                                                 path.prefix)
  # gene count table without first row
  gene.reads.count.table <- gene.count.table.raw$raw.count
  gene.reads.count.gene.name <-
    as.character(gene.count.table.raw$raw.count.name)
  # gene matrix
  gene.count.matrix <- as.matrix(gene.reads.count.table)
  return(list("gene.count.matrix" = gene.count.matrix,
              "gene.count.name" = gene.reads.count.gene.name))
}

RawCountGeneNameChange <- function(raw.count, path.prefix){
  # Convert gene id to gene name
  gene.id <- raw.count$gene_id
  ballgown.texpr <- read.csv(paste0(path.prefix,
                                    "RNASeq_results/ballgown_analysis/",
                                    "ballgown_R_object/texpr.csv"))
  indices <- match(gene.id, ballgown.texpr$gene_id)
  gene_names_for_result <- as.character(ballgown.texpr$gene_name[indices])
  row.names(raw.count) <- raw.count$gene_id
  raw.count$gene_id <- gene_names_for_result
  colnames(raw.count)[1] <- "gene.name"
  # Pre-filter out rowSums bigger than 0 !!
  raw.count <- raw.count[rowSums(raw.count[-1])>0, ]
  # seperate novel gene and known gene
  novel.gene.raw.count <- raw.count[raw.count$gene.name == ".", ]
  known.gene.raw.count <- raw.count[raw.count$gene.name != ".", ]
  # aggregate know gene with same name !
  novel.know.gene.raw.count <- rbind(known.gene.raw.count, novel.gene.raw.count)
  return(list("raw.count" = novel.know.gene.raw.count[-1],
              "raw.count.name" = novel.know.gene.raw.count$gene.name))
}

RawReadCountAvailability <- function(path.prefix) {
  file.prepDE.py <-
    file.exists(paste0(path.prefix,
                       "gene_data/reads_count_matrix/prepDE.py"))
  file.sample.lst.txt <-
    file.exists(paste0(path.prefix,
                       "gene_data/reads_count_matrix/sample_lst.txt"))
  file.gene_count_matrix <-
    file.exists(paste0(path.prefix,
                       "gene_data/reads_count_matrix/gene_count_matrix.csv"))
  if (file.prepDE.py & file.sample.lst.txt & file.gene_count_matrix) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

