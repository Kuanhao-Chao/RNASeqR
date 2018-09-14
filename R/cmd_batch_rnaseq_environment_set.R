#' @title Environment setting for RNA-Seq workflow in background
#'
#' @description Set up the environment for the following RNA-Seq workflow in background.
#'   This function do 4 things :
#'     1. Create file directories.
#'     2. Install necessary tools.
#'     3. Export 'RNASeq_bin/' to the R environment.
#'     4. Check command of tools.
#'   First it will create 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/', 'Rscript/',
#'   'Rscript_out/' directories. Afterwards, 'Hisat2', 'Stringtie', 'Samtools',
#'   'Gffcompare' will be installed under 'RNASeq_bin/Download/' and be unpacked under
#'   'RNASeq_bin/Unpacked/'. 'RNASeq_bin/' will be added to the R environment and validity of
#'   tools will be checked. Any ERROR occurs will be reported and the program will be terminated.
#'   If you want to set up the environment for the following RNA-Seq workflow in R shell,
#'   please see \code{RNASeqEnvironmentSet()} function.
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related parameters
#' @param run Default value is \code{TRUE}. If \code{TRUE}, 'Rscript/Environment_Set.R'
#'   will be created and executed. The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{False}, 'Rscript/Environment_Set.R' will be created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE}, the result of checking \code{RNASeqWorkFlowParam}
#'   will be reported in 'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqWorkFlowParam} will not be in 'Rscript_out/Environment_Set.Rout'
#' @param install.hisat2 Whether to install 'HISAT2' in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'HISAT2' installation.
#' @param install.stringtie Whether to install 'StringTie' in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'StringTie' installation.
#' @param install.gffcompare Whether to install 'Gffcompare' in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'Gffcompare' installation.
#' @param install.samtools Whether to install 'SAMtools' in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'SAMtools' installation.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqEnvironmentSet_CMD(yeast)}
RNASeqEnvironmentSet_CMD <- function(RNASeqWorkFlowParam,
                                     install.hisat2     = TRUE,
                                     install.stringtie  = TRUE,
                                     install.gffcompare = TRUE,
                                     install.samtools   = TRUE,
                                     run                = TRUE,
                                     check.s4.print     = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  os.type <- "@"(RNASeqWorkFlowParam, os.type)
  path.prefix <- "@"(RNASeqWorkFlowParam, path.prefix)
  input.path.prefix <- "@"(RNASeqWorkFlowParam, input.path.prefix)
  genome.name <- "@"(RNASeqWorkFlowParam, genome.name)
  sample.pattern <- "@"(RNASeqWorkFlowParam, sample.pattern)
  indices.optional <- "@"(RNASeqWorkFlowParam, indices.optional)
  MkdirAll(path.prefix)
  fileConn <- file(paste0(path.prefix, "Rscript/Environment_Set.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNASeqEnvironmentSet(path.prefix = '", path.prefix,
                   "', input.path.prefix = '", input.path.prefix,
                   "', genome.name = '", genome.name,
                   "', sample.pattern = '", sample.pattern,
                   "', indices.optional = ", indices.optional,
                   ", os.type = '", os.type,
                   "', install.hisat2 = ", install.hisat2,
                   ", install.stringtie = ", install.stringtie,
                   ", install.gffcompare = ", install.gffcompare,
                   ", install.samtools = ", install.samtools,
                   ", mkdir.bool = ", FALSE, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message("\u2605 '", path.prefix,
          "Rscript/Environment_Set.R' has been created.\n")
  if (run) {
    system2(command = "nohup",
            args = paste0("R CMD BATCH ",
                          path.prefix,
                          "Rscript/Environment_Set.R ",
                          path.prefix,
                          "Rscript_out/Environment_Set.Rout"),
            stdout = "", wait = FALSE)
    message("\u2605 Tools are installing in the background. ",
            "Check current progress in '", path.prefix,
            "Rscript_out/Environment_Set.Rout'\n\n")
  }
}

#' @title Environment setting for RNA-Seq workflow in R shell
#'
#' @description Set up the environment for the following RNA-Seq workflow in R shell.
#'   It is strongly advised to run \code{RNASeqEnvironmentSet_CMD()} directly.
#'   Running \code{RNASeqEnvironmentSet_CMD()} will create 'Environment_Set.Rout' file in 'Rscript_out/' directory.
#'   This function do 4 things :
#'     1. Create file directories.
#'     2. Install necessary tools.
#'     3. Export 'RNASeq_bin/' to the R environment.
#'     4. Check command of tools.
#'   First it will create 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/', 'Rscript/', 'Rscript_out/'
#'   directories. Afterwards, 'Hisat2', 'Stringtie', 'Samtools', 'Gffcompare' will be installed
#'   under 'RNASeq_bin/Download/' and be unpacked under 'RNASeq_bin/Unpacked/'. 'RNASeq_bin/'
#'   will be added to the R environment and validity of tools will be checked. Any ERROR occurs will be
#'   reported and the program will be terminated.
#'   If you want to set up the environment for the following RNA-Seq workflow in background,
#'   please see \code{RNASeqEnvironmentSet_CMD()} function.
#'
#' @param path.prefix path prefix of 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/',
#'   'Rscript/' and 'Rscript_out/' directories
#' @param input.path.prefix path prefix of 'input_files/' directory
#' @param genome.name genome.name variable of genome name defined in this RNA-Seq workflow
#'   (ex. \code{genome.name}.fa, \code{genome.name}.gtf)
#' @param sample.pattern  Regular expression of paired-end fastq.gz files under
#'   'input_files/raw_fastq.gz'. Expression not includes \code{_[1,2].fastq.gz}.
#' @param indices.optional logical value whether 'indices/' is exit in 'input_files/'
#' @param os.type 'linux' or 'osx'. The operating system type
#' @param install.hisat2 Whether to install 'HISAT2' in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'HISAT2' installation.
#' @param install.stringtie Whether to install 'StringTie' in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'StringTie' installation.
#' @param install.gffcompare Whether to install 'Gffcompare' in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'Gffcompare' installation.
#' @param install.samtools Whether to install 'SAMtools' in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'SAMtools' installation.
#' @param mkdir.bool Default \code{TRUE}. If \code{TRUE}, environment directories will be created.
#'   If \code{FALSE}, no directories will be created. When executing RNASeqEnvironmentSet(),
#'   'mkdir.bool' should always be \code{TRUE}
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqEnvironmentSet(path.prefix       = yeast@@path.prefix,
#'                      input.path.prefix = yeast@@input.path.prefix,
#'                      genome.name       = yeast@@genome.name,
#'                      sample.pattern    = yeast@@sample.pattern,
#'                      indices.optional  = yeast@@indices.optional,
#'                      os.type           = yeast@@os.type)}
RNASeqEnvironmentSet <- function(path.prefix,
                                 input.path.prefix,
                                 genome.name,
                                 sample.pattern,
                                 indices.optional,
                                 os.type,
                                 install.hisat2     = TRUE,
                                 install.stringtie  = TRUE,
                                 install.gffcompare = TRUE,
                                 install.samtools   = TRUE,
                                 mkdir.bool         = TRUE) {
  CheckOperatingSystem(FALSE)
  if (mkdir.bool) {
    MkdirAll(path.prefix)
  }
  PreRNASeqEnvironmentSet(path.prefix, sample.pattern)
  CopyInputDir(path.prefix,
               input.path.prefix,
               genome.name,
               sample.pattern,
               indices.optional)
  InstallAll(path.prefix,
             os.type,
             install.hisat2,
             install.stringtie,
             install.gffcompare,
             install.samtools)
  ExportPath(path.prefix)
  PostRNASeqEnvironmentSet(path.prefix,
                           genome.name,
                           sample.pattern)
}

# Create sample gene directory
MkdirGeneDir <- function(path.prefix) {
  message("\u25CF 1. Creating 'gene-data/' directory\n")
  gene_data.dir <- dir.create(file.path(paste0(path.prefix, "gene_data/")),
                              showWarnings = FALSE) == 0
  if (!isTRUE(gene_data.dir)) {
    message("     (\u2714) : Create '", path.prefix, "gene_data/'.\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "gene_data/' has already be created.\n")
  }
  ref.genome.dir <- dir.create(file.path(paste0(path.prefix,
                                                "gene_data/ref_genome/")),
                               showWarnings = FALSE) == 0
  if (!isTRUE(ref.genome.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "gene_data/ref_genome/'.\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "gene_data/ref_genome/' has already be created.\n")
  }
  ref.genes.dir <- dir.create(file.path(paste0(path.prefix,
                                               "gene_data/ref_genes/")),
                              showWarnings = FALSE) == 0
  if (!isTRUE(ref.genes.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "gene_data/ref_genes/'.\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "gene_data/ref_genes/' has already be created.\n")
  }
  indices.dir <- dir.create(file.path(paste0(path.prefix,
                                             "gene_data/indices/")),
                            showWarnings = FALSE) == 0
  if (!isTRUE(indices.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "gene_data/indices/'.\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "gene_data/indices/' has already be created.\n")
  }
  samples.fastq.dir <- dir.create(file.path(paste0(path.prefix,
                                                   "gene_data/raw_fastq.gz")),
                                  showWarnings = FALSE) == 0
  if (!isTRUE(samples.fastq.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "gene_data/raw_fastq.gz/'.\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "gene_data/raw_fastq.gz/' has already be created.\n")
  }
  samples.sam.dir <- dir.create(file.path(paste0(path.prefix,
                                                 "gene_data/raw_sam")),
                                showWarnings = FALSE) == 0
  if (!isTRUE(samples.sam.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "gene_data/raw_sam/'.\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "gene_data/raw_sam/' has already be created.\n")
  }
  samples.bam.dir <- dir.create(file.path(paste0(path.prefix,
                                                 "gene_data/raw_bam/")),
                                showWarnings = FALSE) == 0
  if (!isTRUE(samples.bam.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "gene_data/raw_bam/'.\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "gene_data/raw_bam/' has already be created.\n")
  }
  samples.gtf.dir <- dir.create(file.path(paste0(path.prefix,
                                                 "gene_data/raw_gtf/")),
                                showWarnings = FALSE) == 0
  if (!isTRUE(samples.gtf.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "gene_data/raw_gtf/'.\n\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "gene_data/raw_gtf/' has already be created.\n\n")
  }
}

# Make RNASeq_bin/ directory
MkdirRNASeq_bin <- function(path.prefix) {
  message("\u25CF 2. Creating 'RNASeq_bin/' directory\n")
  RNASeq_bin.dir <- dir.create(file.path(paste0(path.prefix, "RNASeq_bin/")),
                               showWarnings = FALSE) == 0
  if (!isTRUE(RNASeq_bin.dir)) {
    message("     (\u2714) : Create '", path.prefix, "RNASeq_bin/'.\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "RNASeq_bin/' has already be created.\n")
  }
  download.dir <- dir.create(file.path(paste0(path.prefix,
                                              "RNASeq_bin/Download/")),
                             showWarnings = FALSE) == 0
  if (!isTRUE(download.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "RNASeq_bin/Download/'.\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "RNASeq_bin/Download/' has already be created.\n")
  }
  unpacked.dir <- dir.create(file.path(paste0(path.prefix,
                                              "RNASeq_bin/Unpacked/")),
                             showWarnings = FALSE) == 0
  if (!isTRUE(unpacked.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "RNASeq_bin/Unpacked/'.\n\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "RNASeq_bin/Unpacked/' has already be created.\n\n")
  }
}

# Make RNASeq_results/ directory
MkdirRNASeq_results <- function(path.prefix){
  message("\u25CF 3. Creating 'RNASeq_results/' directory\n")
  RNASeq_results.dir <- dir.create(file.path(paste0(path.prefix,
                                                    "RNASeq_results/")),
                                   showWarnings = FALSE) == 0
  if (!isTRUE(RNASeq_results.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "RNASeq_results/'.\n\n")
  } else {
    message("     (\u26A0)) : '", path.prefix,
            "RNASeq_results/' has already be created.\n\n")
  }
}

# Make
MkdirRscript_Rscript_out <- function(path.prefix) {
  message("\u25CF 4. Creating 'Rscript/' directory\n")
  RNASeq_rscript.dir <- dir.create(file.path(paste0(path.prefix, "Rscript/")),
                                   showWarnings = FALSE) == 0
  if (!isTRUE(RNASeq_rscript.dir)) {
    message("     (\u2714) : Create '", path.prefix, "Rscript/'.\n\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
                   "Rscript/' has already be created.\n")
  }
  message("\u25CF 5. Creating 'Rscript_out/' directory\n")
  RNASeq_rscript_out.dir <- dir.create(file.path(paste0(path.prefix,
                                                        "Rscript_out/")),
                                       showWarnings = FALSE) == 0
  if (!isTRUE(RNASeq_rscript_out.dir)) {
    message("     (\u2714) : Create '", path.prefix,
            "Rscript_out/'.\n\n")
  } else {
    message("     (\u26A0) : '", path.prefix,
            "Rscript_out/' has already be created.\n")
  }
}

# Create sample gene and binary directory
MkdirAll <- function(path.prefix) {
  message("************** Creating Directories ************\n")
  MkdirGeneDir(path.prefix)
  MkdirRNASeq_bin(path.prefix)
  MkdirRNASeq_results(path.prefix)
  MkdirRscript_Rscript_out(path.prefix)
}

# inner function : Copy input files directory
CopyInputDir <- function(path.prefix,
                         input.path.prefix,
                         genome.name,
                         sample.pattern,
                         indices.optional) {
  message(c("************** Directory Copying ************\n"))
  message(c("     \u25CF Copying From :",
            paste0(input.path.prefix, "input_files/",
                   genome.name, ".gtf"), "\n"))
  file.symlink(paste0(input.path.prefix, "input_files/", genome.name, ".gtf"),
               paste0(path.prefix, "gene_data/ref_genes/", genome.name, ".gtf"))
  message(c("     \u25CF           To :"),
          paste0(path.prefix, "gene_data/ref_genes/",
                 genome.name, ".gtf", "\n"))
  message(c("     \u25CF Copying From :",
            paste0(input.path.prefix, "input_files/",
                   genome.name, ".fa"),  "\n"))
  file.symlink(paste0(input.path.prefix, "input_files/", genome.name, ".fa"),
               paste0(path.prefix, "gene_data/ref_genome/", genome.name, ".fa"))
  message(c("     \u25CF           To :"),
          paste0(path.prefix, "gene_data/ref_genome/",
                 genome.name, ".fa", "\n"))
  message(c("     \u25CF Copying From :",
            paste0(input.path.prefix, "input_files/", "raw_fastq.gz/"),  "\n"))
  raw_fastq.gz.subfiles <- list.files(path = paste0(input.path.prefix,
                                                    "input_files/",
                                                    "raw_fastq.gz"),
                                      pattern = sample.pattern,
                                      recursive = TRUE,
                                      full.names = TRUE)
  vapply(raw_fastq.gz.subfiles,
         function(x) file.symlink(x, paste0(path.prefix,
                                            "gene_data/raw_fastq.gz")),
         FUN.VALUE = TRUE)
  message(c("     \u25CF           To :",
            paste0(path.prefix, "gene_data/raw_fastq.gz/"), "\n"))
  message(c("     \u25CF Copying From :",
            paste0(input.path.prefix, "input_files/phenodata.csv"),  "\n"))
  file.symlink(paste0(input.path.prefix, "input_files/phenodata.csv")
               , paste0(path.prefix, "gene_data/phenodata.csv"))
  message(c("     \u25CF           To :"),
          paste0(path.prefix, "gene_data/phenodata.csv\n"))
  if (isTRUE(indices.optional)) {
    message(c("     \u25CF Copying From :",
              paste0(input.path.prefix, "input_files/", "indices/"),  "\n"))
    indices.subfiles <- list.files(path = paste0(input.path.prefix,
                                                 "input_files/", "indices"),
                                   pattern = paste0(genome.name,
                                                    "_tran.[0-9].ht2"),
                                   recursive = TRUE,
                                   full.names = TRUE)
    vapply(indices.subfiles,
           function(x) file.symlink(x, paste0(path.prefix,
                                              "gene_data/indices")),
           FUN.VALUE = TRUE)
    message(c("     \u25CF           To :",
              paste0(path.prefix, "gene_data/indices/"), "\n\n"))
  } else {
    message("\n")
  }
}

# Install Hisat2 binay
InstallHisat2Bianry <- function(path.prefix, os.type){
  os <- os.type
  url <- "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/"
  if (os == "linux"){
    os.file.name.zip <- "hisat2-2.1.0-Linux_x86_64.zip"
    os.file.name <- "hisat2-2.1.0"
    url <- paste0(url, os.file.name.zip)
  } else if (os == "osx"){
    os.file.name <- "hisat2-2.1.0-OSX_x86_64.zip"
    os.file.name <- "hisat2-2.1.0"
    url <- paste0(url, os.file.name)
  } else if (os == "windows"){
    stop("Hisat2 is not supporting windows.\n\n")
    return(FALSE)
  } else {
    stop("Unknow operating system.\n\n")
    return(FALSE)
  }
  message("************** Installing Hisat2 (",
          os.file.name.zip, ") ************\n")
  download.file(url, paste0(path.prefix,
                            "RNASeq_bin/Download/",
                            os.file.name.zip))
  message("\n************** Unpacking Hisat2 (",
          os.file.name.zip, ") ************\n")
  if (dir.exists(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))) {
    unlink(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name),
           recursive = TRUE)
  }
  main.command <- "unzip"
  command.result <- system2(command = main.command,
                            args = paste0(path.prefix, "RNASeq_bin/Download/",
                                          os.file.name.zip, " -d ", path.prefix,
                                          "RNASeq_bin/Unpacked/"))
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  message("hisat2 binaries\n\n")
  current.path <- getwd()
  setwd(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name, "/"))
  message("\n************** Moving Hisat2 Binary ************")
  main.command <- "cp"
  command.result <- system2(command = main.command,
                            args = c("hisat2*", "*.py",
                                     paste0(path.prefix, "RNASeq_bin/")),
                            stderr = FALSE)
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  on.exit(setwd(current.path))
  message("\n'", path.prefix, "RNASeq_bin/Download/",
          os.file.name.zip, "' has been installed.\n")
  message("Hisat2 has been unpacked. ('", path.prefix,
          "RNASeq_bin/Unpacked/", os.file.name, "/')", "\n\n")
  return(TRUE)
}

# Install stringtie binary
InstallStringTieBinary <- function(path.prefix, os.type){
  os <- os.type
  url <- "http://ccb.jhu.edu/software/stringtie/dl/"
  if (os == "linux"){
    os.file.name.zip <- "stringtie-1.3.4d.Linux_x86_64.tar.gz"
    os.file.name <- "stringtie-1.3.4d.Linux_x86_64"
    url <- paste0(url, os.file.name.zip)
  } else if (os == "osx"){
    os.file.name <- "stringtie-1.3.4d.OSX_x86_64.tar.gz"
    os.file.name <- "stringtie-1.3.4d.OSX_x86_64"
    url <- paste0(url, os.file.name)
  } else if (os == "windows"){
    stop("Stringtie is not supporting windows.\n\n")
    return(FALSE)
  } else {
    stop("Unknow operating system.\n\n")
    return(FALSE)
  }
  message("************** Installing stringtie (",
          os.file.name.zip, ") ************\n")
  download.file(url,
                paste0(path.prefix, "RNASeq_bin/Download/", os.file.name.zip))
  message("\n************** Unpacking stringtie (",
          os.file.name.zip, ") ************\n")
  if (dir.exists(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))) {
    unlink(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name),
           recursive = TRUE)
  }
  main.command <- "tar"
  command.result <- system2(command = main.command,
                            args = c("xvzf",
                                     paste0(path.prefix,
                                            "RNASeq_bin/Download/",
                                            os.file.name.zip),
                                     "-C",
                                     paste0(path.prefix,
                                            "RNASeq_bin/Unpacked/")))
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  current.path <- getwd()
  setwd(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))
  message("\n************** Moving stringtie Binary ************")
  main.command <- "cp"
  command.result <- system2(command = main.command,
                            args = c("stringtie",
                                     paste0(path.prefix, "RNASeq_bin/")),
                            stderr = FALSE)
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  on.exit(setwd(current.path))
  message("\n'", path.prefix,
          "RNASeq_bin/Download/",
          os.file.name.zip, "' has been installed.\n")
  message("StringTie has been unpacked. ('",
          path.prefix, "RNASeq_bin/Unpacked/",
          os.file.name, "')", "\n\n")
  return(TRUE)
}

# Install Gffcompare binary
InstallGffcompareBinary <- function(path.prefix, os.type){
  os <- os.type
  current.path <- getwd()
  url <- "http://ccb.jhu.edu/software/stringtie/dl/"
  if (os == "linux"){
    os.file.name.zip <- "gffcompare-0.10.4.Linux_x86_64.tar.gz"
    os.file.name <- "gffcompare-0.10.4.Linux_x86_64"
    url <- paste0(url, os.file.name.zip)
  } else if (os == "osx"){
    os.file.name <- "gffcompare-0.10.4.OSX_x86_64.tar.gz"
    os.file.name <- "gffcompare-0.10.4.OSX_x86_64"
    url <- paste0(url, os.file.name)
  } else if (os == "windows"){
    stop("Gffcompare is not supporting windows.\n\n")
    return(FALSE)
  } else {
    stop("Unknow operating system.\n\n")
    return(FALSE)
  }
  message("************** Installing gffcompare (",
          os.file.name.zip, ") ************\n")
  download.file(url,
                paste0(path.prefix, "RNASeq_bin/Download/", os.file.name.zip))
  message("\n************** Unpacking gffcompare (",
          os.file.name.zip, ") ************\n")
  if (dir.exists(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))) {
    unlink(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name),
           recursive = TRUE)
  }
  main.command <- "tar"
  command.result <- system2(command = main.command,
                            args = c("xvzf",
                                     paste0(path.prefix,
                                            "RNASeq_bin/Download/",
                                            os.file.name.zip),
                                     "-C",
                                     paste0(path.prefix,
                                            "RNASeq_bin/Unpacked/")))
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  current.path <- getwd()
  setwd(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))
  message("\n************** Moving gffcompare Binary ************")
  main.command <- "cp"
  command.result <- system2(command = main.command,
                            args = c("gffcompare",
                                     paste0(path.prefix, "RNASeq_bin/")),
                            stderr = FALSE)
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  on.exit(setwd(current.path))
  message("\n'", path.prefix, "RNASeq_bin/Download/",
          os.file.name.zip, "' has been installed.\n")
  message("Gffcompare has been unpacked. ('", path.prefix,
          "RNASeq_bin/Unpacked/", os.file.name, "')", "\n\n")
  return(TRUE)
}

# Install Samtools binary
InstallSamtoolsBinary <- function(path.prefix, os.type){
  os <- os.type
  current.path <- getwd()
  url <- paste0("https://github.com/samtools/samtools/releases/",
                "download/1.8/samtools-1.8.tar.bz2")
  if (os == "linux"){
    os.file.name.zip <- "samtools-1.8.tar.bz2"
    os.file.name <- "samtools-1.8"
  } else if (os == "osx"){
    os.file.name <- "samtools-1.8.tar.bz2"
    os.file.name <- "samtools-1.8"
  } else if (os == "windows"){
    stop("Samtools is not supporting windows.\n\n")
    return(FALSE)
  } else {
    stop("Unknow operating system.\n\n")
    return(FALSE)
  }
  message("************** Installing samtools (",
          os.file.name.zip, ") ************\n")
  main.command <- "curl"
  command.result <- system2(command = main.command,
                            args = c("-L",
                                     paste0("https://github.com/samtools/",
                                            "samtools/releases/download/1.8/",
                                            "samtools-1.8.tar.bz2"), ">",
                                     paste0(path.prefix, "RNASeq_bin/Download/",
                                            os.file.name.zip)),
                            stdout = "",
                            wait = TRUE)
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  message("\n************** Unpacking samtools (",
          os.file.name.zip, ") ************\n")
  if (dir.exists(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))) {
    unlink(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name),
           recursive = TRUE)
  }
  main.command <- "tar"
  command.result <- system2(command = main.command,
                            args = c("jxvf", paste0(path.prefix,
                                                    "RNASeq_bin/Download/",
                                                    os.file.name.zip),
                                     "-C",
                                     paste0(path.prefix,
                                            "RNASeq_bin/Unpacked/")))
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  current.path <- getwd()
  message("\n************** Making samtools (",
          os.file.name, ") ************")
  setwd(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))
  main.command <- "make"
  command.result <- system2(command = main.command,
                            args = "clean",
                            stderr = FALSE)
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  command.result <- system2(command = main.command)
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop(paste0("'", main.command, "' ERROR"))
  }
  message("\n************** Moving samtools Binary ************")
  file.copy("samtools", paste0(path.prefix, "RNASeq_bin/"))
  on.exit(setwd(current.path))
  message("\n'", path.prefix, "RNASeq_bin/Download/",
          os.file.name.zip, "' has been installed.\n")
  message("Samtools has been unpacked. ('", path.prefix,
          "RNASeq_bin/Unpacked/", os.file.name, "')", "\n\n")
  return(TRUE)
}

# Install 'Hisat2', 'StringTie', 'Gffcompare', 'Samtools'
InstallAll <- function(path.prefix,
                       os.type,
                       install.hisat2,
                       install.stringtie,
                       install.gffcompare,
                       install.samtools) {
  message("\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618  ",
          "Start installing ... ",
          "\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618\n")
  install.software <- ""
  if (install.hisat2) {
    install.software <- paste0(install.software, "\u25CF'hisat2' ")
  }
  if (install.stringtie) {
    install.software <- paste0(install.software, "\u25CF'stringtie' ")
  }
  if (install.gffcompare) {
    install.software <- paste0(install.software, "\u25CF'gffcompare' ")
  }
  if (install.samtools) {
    install.software <- paste0(install.software, "\u25CF'samtools' ")
  }
  if (!install.hisat2 & !install.stringtie &
      !install.gffcompare & !install.samtools) {
    message("   \u261E\u261E skipping installation process ... \n")
  } else {
    message("   \u261E\u261E ",
            install.software,
            "will be installed. ... \n")
    message("   \u261E\u261E  Compressed files will be in '",
            path.prefix, "RNASeq_bin/Download/'", "\n")
    message("   \u261E\u261E  Unpacked files will be in '",
            path.prefix, "RNASeq_bin/Unpacked/'", "\n")
    message("   \u261E\u261E  Binary files will be copied to '",
            path.prefix, "RNASeq_bin/'", "\n\n")
  }
  if (install.hisat2) {
    message("\u2618\u2618 Hisat2 processing ...\n")
    InstallHisat2Bianry(path.prefix, os.type)
  }
  if (install.stringtie) {
    message("\u2618\u2618 Stringtie processing ...\n")
    InstallStringTieBinary(path.prefix, os.type)
  }
  if (install.gffcompare) {
    message("\u2618\u2618 Gffcompare processing ...\n")
    InstallGffcompareBinary(path.prefix, os.type)
  }
  if (install.samtools) {
    message("\u2618\u2618 Samtools processing ...\n")
    InstallSamtoolsBinary(path.prefix, os.type)
  }
}

# Check 'hisat2'
CheckHisat2 <- function(print=TRUE){
  if (print) {
    message("\u25CF  Checking hisat2 command\n")
  }
  hisat2.installed <- system("hisat2 --version",
                             ignore.stdout = !print,
                             ignore.stderr = !print) == 0
  if (isTRUE(hisat2.installed)){
    if (isTRUE(print)){
      message("(\u2714) : 'hisat2' is installed\n\n")
    }
    return(TRUE)
  }
  else{
    message("(\u2718) : 'hisat2' command is not found on this device. ",
            "Please run 'InstallAll()' to install the necessary ",
            "programs or 'ExportPath' to update the path.\n\n")
    return(FALSE)
  }
}

# Check 'stringtie'
CheckStringTie <- function(print=TRUE){
  if (print){
    message("\u25CF  Checking stringtie command\n")
  }
  stringtie.installed <- system("stringtie --version",
                                ignore.stdout = !print,
                                ignore.stderr = !print) == 0
  if (isTRUE(stringtie.installed)){
    if (isTRUE(print)){
      message("(\u2714) : 'stringtie' is installed\n\n")
    }
    return(TRUE)
  }
  else{
    message("(\u2718) : 'stringtie' command is not found on this device. ",
            "Please run 'InstallAll()' to install the necessary programs ",
            "or 'ExportPath' to update the path.\n\n")
    return(FALSE)
  }
}

# Check 'Gffcompare'
CheckGffcompare <- function(print=TRUE) {
  if (print) {
    message("\u25CF  Checking gffcompare command\n")
  }
  gffcompare.old <- system( "gffcompare --version",
                            ignore.stdout = !print,
                            ignore.stderr = !print) == 0
  if (isTRUE(gffcompare.old)){
    if (isTRUE(print)){
      message("(\u2714) : 'gffcompare' is installed\n\n")
    }
    return(TRUE)
  }
  else{
    message("(\u2718) : \'gffcompare\' command is not found on this device. ",
            "Please run 'InstallAll()' to install the necessary programs ",
            "or 'ExportPath' to update the path.\n\n")
    return(FALSE)
  }
}

# Check 'Samtools'
CheckSamtools <- function(print=TRUE){
  if (print) {
    message("\u25CF  Checking samtools command\n")
  }
  samtools.old <- system( "samtools --version",
                          ignore.stdout = !print,
                          ignore.stderr = !print) == 0
  if (isTRUE(samtools.old)){
    if (isTRUE(print)){
      message("(\u2714) : 'samtools' is installed\n\n")
    }
    return(TRUE)
  }
  else{
    message("(\u2718) : \'samtools\' command is not found on this device. ",
            "Please run 'InstallAll()' to install the necessary programs ",
            "or 'ExportPath' to update the path.\n\n")
    return(FALSE)
  }
}

#' @title Check 'Hisat2', 'StringTie', 'Samtools' and 'Gffcompare' for this workflow
#'
#' @description Check whether 'Hisat2', 'Stringtie', 'Samtools' and 'Gffcompare' are installed on the workstation
#'
#' @param print If \code{TRUE}, detailed information will be printed. If \code{FALSE}, detailed information will not be printed.
#'
#' @return None
#' @export
#' @example
#' CheckToolAll()
CheckToolAll <- function(print=TRUE) {
  message("************** Checking Availability of Commands ************\n")
  hisat2.check <- CheckHisat2(print)
  stringtie.check <- CheckStringTie(print)
  gff.check <- CheckGffcompare(print)
  samtool.check <- CheckSamtools(print)
  if (isTRUE(hisat2.check) &&
      isTRUE(stringtie.check) &&
      isTRUE(gff.check) &&
      isTRUE(samtool.check)){
    return(TRUE)
  } else {
    stop(paste0("(\u2718) Necessary program is missing.\n     ",
                "1. Check 'INSTALL_TOOLS.Rout' whether tools are ",
                "properly installed.\n     ",
                "2. Run 'ExportPath()' to set the environment.\n\n"))
    return(FALSE)
  }
}

PreRNASeqEnvironmentSet <- function(path.prefix, sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqEnvironmentSet()' ",
          "environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqEnvironmentSet environment() ERROR")
  }
  message("(\u2714) : RNASeqEnvironmentSet() pre-check is valid\n\n")
}

PostRNASeqEnvironmentSet <- function(path.prefix,
                                     genome.name,
                                     sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqEnvironmentSet()' ",
          "environment post-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  ref.gtf <- file.exists(paste0(path.prefix,
                                "gene_data/ref_genes/",
                                genome.name, ".gtf"))
  ref.fa <- file.exists(paste0(path.prefix,
                               "gene_data/ref_genome/",
                               genome.name, ".fa"))
  raw.fastq <- list.files(path = paste0(path.prefix, "gene_data/raw_fastq.gz/"),
                          pattern = sample.pattern,
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
  check.tool.result <- CheckToolAll()
  validity <- phenodata.csv && ref.gtf && ref.fa &&
    check.tool.result && (length(raw.fastq) != 0)
  if (!isTRUE(validity)) {
    stop("RNASeqEnvironmentSet() post-check ERROR")
  }
  message("(\u2714) : RNASeqEnvironmentSet() post-check is valid\n\n")
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
