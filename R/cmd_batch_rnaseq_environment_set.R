#' Set up the environment for the following RNA seq pipeline in background.
#'
#' To set up the environment for the following RNA seq pipeline in background, create 'RNASeqWorkFlowParam' first, and then run 'RNAseqEnvironmentSet_CMD(RNASeqWorkFlowParam)'.
#' If you want to set up the environment for the following RNA seq pipeline in R shell, please see 'RNAseqEnvironmentSet()' function.
#' This function do 4 things : 1. Create file directories. 2. Install necessary tools. 3. Export 'RNAseq_bin/' to the R environment. 4. Check command of tools.
#' First it will create 'gene_data/', 'RNAseq_bin/', 'RNAseq_results/', 'Rscript/', 'Rscript_out/' directories. Afterwards, 'Hisat2', 'Stringtie', 'Samtools',
#' 'Gffcompare' will be installed under 'RNAseq_bin/'. 'RNAseq_bin/' will be added to the R environment and validity of tools will be checked. Any ERROR occurs will be reported.
#'
#' @param RNASeqWorkFlowParam S4 object instance of experiment-related parameters
#'
#' @export
#' @example
#' exp <- RNASeqWorkFlowParam(path.prefix = "/home/rnaseq", input.path.prefix = "/home", gene.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNAseqEnvironmentSet_CMD(RNASeqWorkFlowParam <- exp)
RNAseqEnvironmentSet_CMD <- function(RNASeqWorkFlowParam, run = TRUE, check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqWorkFlowParam, check.s4.print)
  os.type <- RNASeqWorkFlowParam@os.type
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  gene.name <- RNASeqWorkFlowParam@gene.name
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  indexes.optional <- RNASeqWorkFlowParam@indexes.optional
  MkdirAll(path.prefix)
  fileConn<-file(paste0(path.prefix, "Rscript/Environment_Set.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqEnvironmentSet(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', gene.name = '", gene.name, "', sample.pattern = '", sample.pattern, "', indexes.optional = ",indexes.optional, ", os.type = '", os.type, "')")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  cat(paste0("\u2605 '", path.prefix, "Rscript/Environment_Set.R' has been created.\n"))
  if (run) {
    system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/Environment_Set.R ", path.prefix, "Rscript_out/Environment_Set.Rout"), stdout = "", wait = FALSE)
    cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/Environment_Set.Rout'\n\n"))
  }
}

#' Set up the environment for the RNA seq pipeline in R shell
#'
#' 1. To set up the environment for the following RNA seq pipeline in R shell, create 'RNASeqWorkFlowParam' first, get the value from the instance of this S4 object, and then put set the corresponding value to the input of this function.
#' 2. To set up the environment for the following RNA seq pipeline in background, create 'RNASeqWorkFlowParam' first, and run 'RNAseqEnvironmentSet_CMD(RNASeqWorkFlowParam)'.
#' It is recommended to set up the environment in background.
#'
#' @param os.type 'linux' or 'osx'. The operating system type
#' @param path.prefix the directory holding installations and analysis results
#' @param input.path.prefix user has to prepared valid 'input_files/' under this directory
#' @param gene.name gene name defined in this RNA-Seq workflow (ex. gene.name.fa, gene.name.gtf)
#' @param sample.pattern  sample pattern describing the name of raw fastq.gz files
#' @param indexes.optional logical value whether indexes/ is exit in 'input_files/'
#'
#' @export
#'
#' @example
#' exp <- RNASeqWorkFlowParam(path.prefix = "/home/rnaseq", input.path.prefix = "/home", gene.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
#' RNAseqEnvironmentSet(path.prefix = exp@@path.prefix, input.path.prefix = exp@@input.path.prefix, gene.name = exp@@gene.name,
#'                      sample.pattern = exp@@sample.pattern, indexes.optional = exp@@indexes.optional, os.type = exp@@os.type)
#'
RNAseqEnvironmentSet <- function(path.prefix, input.path.prefix, gene.name, sample.pattern, indexes.optional, os.type) {
  PreRNAseqEnvironmentSet(path.prefix = path.prefix, sample.pattern = sample.pattern)
  CopyInputDir(path.prefix = path.prefix, input.path.prefix = input.path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, indexes.optional = indexes.optional)
  InstallAll(path.prefix = path.prefix, os.type = os.type)
  ExportPath(path.prefix = path.prefix)
  PostRNAseqEnvironmentSet(path.prefix = path.prefix, sample.pattern = sample.pattern)
}

#' Create sample gene directory
MkdirGeneDir <- function(path.prefix) {
  cat("\u25CF 1. Creating 'gene-data/' directory\n")
  gene_data.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/')), showWarnings = FALSE) == 0
  if (!isTRUE(gene_data.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/'.\n     Please check whether the directory is already exit.\n"))
  }
  ref.genome.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/ref_genome/')), showWarnings = FALSE) == 0
  if (!isTRUE(ref.genome.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/ref_genome/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/ref_genome/'.\n     Please check whether the directory is already exit.\n"))
  }
  ref.genes.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/ref_genes/')), showWarnings = FALSE) == 0
  if (!isTRUE(ref.genes.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/ref_genes/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/ref_genes/'.\n     Please check whether the directory is already exit.\n"))
  }
  indexes.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/indexes/')), showWarnings = FALSE) == 0
  if (!isTRUE(indexes.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/indexes/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/indexes/'.\n     Please check whether the directory is already exit.\n"))
  }
  samples.fastq.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/raw_fastq.gz')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.fastq.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/raw_fastq.gz/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/raw_fastq.gz/'.\n     Please check whether the directory is already exit.\n"))
  }
  samples.sam.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/raw_sam')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.sam.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/raw_sam/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/raw_sam/'.\n     Please check whether the directory is already exit.\n"))
  }
  samples.bam.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/raw_bam')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.bam.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/raw_bam/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/raw_bam/'.\n     Please check whether the directory is already exit.\n"))
  }
  samples.gtf.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/raw_gtf')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.gtf.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/raw_gtf/'.\n\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/raw_gtf/'.\n     Please check whether the directory is already exit.\n\n"))
  }
  samples.gene_abundance.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/gene_abundance')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.gene_abundance.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/gene_abundance/'.\n\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/gene_abundance/'.\n     Please check whether the directory is already exit.\n\n"))
  }
  samples.reads_count_matrix.dir <- dir.create(file.path(paste0(path.prefix, 'gene_data/reads_count_matrix')), showWarnings = FALSE) == 0
  if (!isTRUE(samples.reads_count_matrix.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "gene_data/reads_count_matrix/'.\n\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "gene_data/reads_count_matrix/'.\n     Please check whether the directory is already exit.\n\n"))
  }
}

#' Make RNAseq_bin/ directory
MkdirRNAseq_bin <- function(path.prefix) {
  cat("\u25CF 2. Creating 'RNAseq_bin/' directory\n")
  RNAseq_bin.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_bin/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_bin.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_bin/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_bin/'.\n     Please check whether the directory is already exit.\n"))
  }
  download.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_bin/Download/')), showWarnings = FALSE) == 0
  if (!isTRUE(download.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_bin/Download/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_bin/Download/'.\n     Please check whether the directory is already exit.\n"))
  }
  unpacked.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_bin/Unpacked/')), showWarnings = FALSE) == 0
  if (!isTRUE(unpacked.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_bin/Unpacked/'.\n\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_bin/Unpacked/'.\n     Please check whether the directory is already exit.\n\n"))
  }
}

#' Make RNAseq_results/ directory
MkdirRNAseq_results <- function(path.prefix){
  cat("\u25CF 3. Creating 'RNAseq_results/' directory\n")
  RNAseq_results.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_results/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_results.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_results/'.\n"))
  } else {
    cat(paste0("     (\u26A0)) : Fail to create '", path.prefix, "RNAseq_results/'.\n     Please check whether the directory is already exit.\n"))
  }
  RNAseq_results_Raw_reads.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_results/Raw_reads/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_results_Raw_reads.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_results/Raw_reads/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_results/Raw_reads/'.\n     Please check whether the directory is already exit.\n"))
  }
  RNAseq_results_Raw_reads_gene.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_results/Raw_reads/gene/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_results_Raw_reads_gene.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_results/Raw_reads/gene/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_results/Raw_reads/gene/'.\n     Please check whether the directory is already exit.\n"))
  }
  RNAseq_results_Raw_reads_transcript.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_results/Raw_reads/transcript/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_results_Raw_reads_transcript.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_results/Raw_reads/transcript/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_results/Raw_reads/transcript/'.\n     Please check whether the directory is already exit.\n"))
  }
  RNAseq_results_quality_control.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_results/QA_results/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_results_quality_control.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_results/QA_results/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_results/QA_results/'.\n     Please check whether the directory is already exit.\n"))
  }
  RNAseq_results_quality_control_Rqc.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_results/QA_results/Rqc/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_results_quality_control_Rqc.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_results/QA_results/Rqc/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_results/QA_results/Rqc/'.\n     Please check whether the directory is already exit.\n"))
  }
  RNAseq_results_quality_control_systemPipeR.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_results/QA_results/systemPipeR/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_results_quality_control_systemPipeR.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_results/QA_results/systemPipeR/'.\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_results/QA_results/systemPipeR/'.\n     Please check whether the directory is already exit.\n"))
  }
  RNAseq_results_quality_control_shortread.dir <- dir.create(file.path(paste0(path.prefix, 'RNAseq_results/QA_results/ShortRead/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_results_quality_control_shortread.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "RNAseq_results/QA_results/ShortRead/'.\n\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "RNAseq_results/QA_results/ShortRead/'.\n     Please check whether the directory is already exit.\n"))
  }
}

#' Make
MkdirRscript_Rscript_out <- function(path.prefix) {
  cat("\u25CF 4. Creating 'Rscript/' directory\n")
  RNAseq_rscript.dir <- dir.create(file.path(paste0(path.prefix, 'Rscript/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_rscript.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "Rscript/'.\n\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "Rscript/'.\n     Please check whether the directory is already exit.\n"))
  }
  cat("\u25CF 5. Creating 'Rscript_out/' directory\n")
  RNAseq_rscript_out.dir <- dir.create(file.path(paste0(path.prefix, 'Rscript_out/')), showWarnings = FALSE) == 0
  if (!isTRUE(RNAseq_rscript_out.dir)) {
    cat(paste0("     (\u2714) : Create '", path.prefix, "Rscript_out/'.\n\n"))
  } else {
    cat(paste0("     (\u26A0) : Fail to create '", path.prefix, "Rscript_out/'.\n     Please check whether the directory is already exit.\n"))
  }
}

#' Create sample gene and binary directory
MkdirAll <- function(path.prefix) {
  cat("************** Creating Directories ************\n")
  MkdirGeneDir(path.prefix)
  MkdirRNAseq_bin(path.prefix)
  MkdirRNAseq_results(path.prefix)
  MkdirRscript_Rscript_out(path.prefix)
}

#' inner function : Copy input files directory
CopyInputDir <- function(path.prefix, input.path.prefix, gene.name, sample.pattern, indexes.optional = indexes.optional) {
  current.path <- getwd()
  setwd(paste0(path.prefix, "gene_data/"))
  cat(c("************** Directory Copying ************\n"))
  cat(c("     \u25CF Copying From :", paste0(input.path.prefix, "input_files/", gene.name, ".gtf"), "\n"))
  file.copy(paste0(input.path.prefix, "input_files/", gene.name, ".gtf"), paste0(getwd(), "/ref_genes/", gene.name, ".gtf"))
  cat(c("     \u25CF           To :"), paste0(getwd(), "/ref_genes/", gene.name, ".gtf", "\n"))
  cat(c("     \u25CF Copying From :", paste0(input.path.prefix, "input_files/", gene.name, ".fa"),  "\n"))
  file.copy(paste0(input.path.prefix, "input_files/", gene.name, ".fa"), paste0(getwd(), "/ref_genome/", gene.name, ".fa"))
  cat(c("     \u25CF           To :"), paste0(getwd(), "/ref_genome/", gene.name, ".fa", "\n"))
  cat(c("     \u25CF Copying From :", paste0(input.path.prefix, "input_files/", "raw_fastq.gz/"),  "\n"))
  file.copy(paste0(input.path.prefix, "input_files/", "raw_fastq.gz/"), paste0(getwd(), "/"), overwrite = TRUE, recursive = TRUE)
  cat(c("     \u25CF           To :", paste0(getwd(),"/raw_fastq.gz/"), "\n"))
  cat(c("     \u25CF Copying From :", paste0(input.path.prefix, "input_files/phenodata.csv"),  "\n"))
  file.copy(paste0(input.path.prefix, "input_files/phenodata.csv"), paste0(getwd(), "/phenodata.csv"))
  cat(c("     \u25CF           To :"), paste0(getwd(), "/phenodata.csv\n"))
  on.exit(setwd(current.path))
  if (isTRUE(indexes.optional)) {
    cat(c("     \u25CF Copying From :", paste0(input.path.prefix, "input_files/", "indexes/"),  "\n"))
    file.copy(paste0(input.path.prefix, "input_files/", "indexes/"), paste0(getwd(), "/"), overwrite = TRUE, recursive = TRUE)
    cat(c("     \u25CF           To :", paste0(getwd(),"/indexes/"), "\n\n"))
  } else {
    cat("\n")
  }
}

#' Install Hisat2 binay
InstallHisat2Bianry <- function(path.prefix, os.type){
  os <- os.type
  url <- 'ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/'
  # setwd(paste0(path.prefix, "RNAseq_bin/"))
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
  cat(paste0("************** Installing Hisat2 ", "(", os.file.name.zip, ") ************\n"))
  download.file(url, paste0(path.prefix, "RNAseq_bin/Download/", os.file.name.zip))
  cat(paste0("\n************** Unpacking Hisat2 ", "(", os.file.name.zip, ") ************\n"))
  system2(command = 'unzip', args = paste0(path.prefix, "RNAseq_bin/Download/", os.file.name.zip," -d ", path.prefix, "RNAseq_bin/Unpacked/"))
  current.path <- getwd()
  setwd(paste0(path.prefix, "RNAseq_bin/Unpacked/", os.file.name, "/"))
  cat("\n************** Moving Hisat2 Binary ************")
  system2(command = 'cp', args = c("hisat2*", "*.py", paste0(path.prefix, "RNAseq_bin/")), stderr = FALSE)
  on.exit(setwd(current.path))
  cat(paste0("\n'", path.prefix, "RNAseq_bin/Download/", os.file.name.zip,"' has been installed.\n"))
  cat(paste0("Hisat2 has been unpacked. ('", path.prefix, "RNAseq_bin/Unpacked/", os.file.name, "/')"), "\n\n")
  return(TRUE)
}

#' Install stringtie binary
InstallStringTieBinary <- function(path.prefix, os.type){
  os <- os.type
  url <- 'http://ccb.jhu.edu/software/stringtie/dl/'
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
  cat(paste0("************** Installing stringtie ", "(", os.file.name.zip, ") ************\n"))
  download.file(url, paste0(path.prefix, "RNAseq_bin/Download/", os.file.name.zip))
  cat(paste0("\n************** Unpacking stringtie ", "(", os.file.name.zip, ") ************\n"))
  system2(command = 'tar', args = c("xvzf", paste0(path.prefix, "RNAseq_bin/Download/", os.file.name.zip), "-C", paste0(path.prefix, "RNAseq_bin/Unpacked/")))
  current.path <- getwd()
  setwd(paste0(path.prefix, "RNAseq_bin/Unpacked/", os.file.name))
  cat("\n************** Moving stringtie Binary ************")
  system2(command = 'cp', args = c("stringtie", paste0(path.prefix, "RNAseq_bin/")), stderr = FALSE)
  on.exit(setwd(current.path))
  cat(paste0("\n'", path.prefix, "RNAseq_bin/Download/", os.file.name.zip,"' has been installed.\n"))
  cat(paste0("StringTie has been unpacked. ('", path.prefix, "RNAseq_bin/Unpacked/", os.file.name, "')"), "\n\n")
  return(TRUE)
}

#' Install Gffcompare binary
InstallGffcompareBinary <- function(path.prefix, os.type){
  os <- os.type
  current.path <- getwd()
  url <- 'http://ccb.jhu.edu/software/stringtie/dl/'
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
  cat(paste0("************** Installing gffcompare ", "(", os.file.name.zip, ") ************\n"))
  download.file(url, paste0(path.prefix, "RNAseq_bin/Download/", os.file.name.zip))
  cat(paste0("\n************** Unpacking gffcompare ", "(", os.file.name.zip, ") ************\n"))
  system2(command = 'tar', args = c("xvzf", paste0(path.prefix, "RNAseq_bin/Download/", os.file.name.zip), "-C", paste0(path.prefix, "RNAseq_bin/Unpacked/")))
  current.path <- getwd()
  setwd(paste0(path.prefix, "RNAseq_bin/Unpacked/", os.file.name))
  cat("\n************** Moving gffcompare Binary ************")
  system2(command = 'cp', args = c("gffcompare", paste0(path.prefix, "RNAseq_bin/")), stderr = FALSE)
  on.exit(setwd(current.path))
  cat(paste0("\n'", path.prefix, "RNAseq_bin/Download/", os.file.name.zip,"' has been installed.\n"))
  cat(paste0("Gffcompare has been unpacked. ('", path.prefix, "RNAseq_bin/Unpacked/", os.file.name, "')"), "\n\n")
  return(TRUE)
}

#' Install Samtools binary
InstallSamtoolsBinary <- function(path.prefix, os.type){
  os <- os.type
  current.path <- getwd()
  url <- 'https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2'
  if (os == "linux"){
    os.file.name.zip <- "samtools-1.8.tar.bz2"
    os.file.name <- "samtools-1.8"
    #url <- paste0(url, os.file.name.zip)
  } else if (os == "osx"){
    os.file.name <- "samtools-1.8.tar.bz2"
    os.file.name <- "samtools-1.8"
    #url <- paste0(url, os.file.name)
  } else if (os == "windows"){
    stop("Samtools is not supporting windows.\n\n")
    return(FALSE)
  } else {
    stop("Unknow operating system.\n\n")
    return(FALSE)
  }
  cat(paste0("************** Installing samtools ", "(", os.file.name.zip, ") ************\n"))
  system2(command = 'curl', args = c('-L', 'https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2', '>', paste0(path.prefix, "RNAseq_bin/Download/", os.file.name.zip)), stdout = "", wait = TRUE)
  cat(paste0("\n************** Unpacking samtools ", "(", os.file.name.zip, ") ************\n"))
  system2(command = 'tar', args = c("jxvf", paste0(path.prefix, "RNAseq_bin/Download/", os.file.name.zip), "-C", paste0(path.prefix, "RNAseq_bin/Unpacked/")))
  current.path <- getwd()
  cat(paste0("\n************** Making samtools ", "(", os.file.name, ") ************"))
  setwd(paste0(path.prefix, "RNAseq_bin/Unpacked/", os.file.name))
  system2(command = 'make', args = "clean", stderr = FALSE)
  system2(command = 'make')
  cat("\n************** Moving samtools Binary ************")
  file.copy("samtools", paste0(path.prefix, "RNAseq_bin/"))
  on.exit(setwd(current.path))
  cat(paste0("\n'", path.prefix, "RNAseq_bin/Download/", os.file.name.zip,"' has been installed.\n"))
  cat(paste0("Samtools has been unpacked. ('", path.prefix, "RNAseq_bin/Unpacked/", os.file.name, "')"), "\n\n")
  return(TRUE)
}

#' Install 'Hisat2', 'StringTie', 'Gffcompare', 'Samtools'
#'
#' This function is not designed to run directly. Please run 'InstallToolsCMD()' to accomplish the whole installing process.
#'
#' @param path.prefix directory that will store all the files created throughout the pipeline
#' @param os.type the operating system type of the current workstation. Only 'linux' and 'osx'are valid
InstallAll <- function(path.prefix, os.type) {
  cat("\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618  Start installing ... \u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618\n")
  cat("   \u261E\u261E  \u25CF'hisat2', \u25CF'stringtie', \u25CF'gffcompare', \u25CF'samtools' will be installed. ... \n")
  cat(paste0("   \u261E\u261E  Compressed files will be in '", path.prefix, "RNAseq_bin/Download/'"), "\n")
  cat(paste0("   \u261E\u261E  Unpacked files will be in '", path.prefix, "RNAseq_bin/Unpacked/'"), "\n")
  cat(paste0("   \u261E\u261E  Binary files will be copied to '", path.prefix, "RNAseq_bin/'"), "\n\n")
  cat("\u2618\u2618 Hisat2 processing ...\n")
  InstallHisat2Bianry(path.prefix, os.type)
  cat("\u2618\u2618 Stringtie processing ...\n")
  InstallStringTieBinary(path.prefix, os.type)
  cat("\u2618\u2618 Gffcompare processing ...\n")
  InstallGffcompareBinary(path.prefix, os.type)
  cat("\u2618\u2618 Samtools processing ...\n")
  InstallSamtoolsBinary(path.prefix, os.type)
  #return(CheckToolAll(print=FALSE))
}


#' Check 'hisat2'
#'
#' Check whether 'hisat2' is installed on the workstation
#'
CheckHisat2 <- function(print=TRUE){
  if (print) {
    cat(paste0("     \u25CF  Checking hisat2 command\n"))
  }
  hisat2.installed <- system('hisat2 --version', ignore.stdout = !print , ignore.stderr = !print)==0
  if( isTRUE(hisat2.installed)){
    if(isTRUE(print)){
      cat("(\u2714) : 'hisat2' is installed\n\n")
    }
    return(TRUE)
  }
  else{
    cat("(\u2718) : 'hisat2' command is not found on this device. Please run 'InstallAll()' to install the necessary programs or 'ExportPath' to update the path.\n\n")
    return(FALSE)
  }
}

#' Check s'tringtie'
#'
#' Check whether 'stringtie' is installed on the workstation
#'
CheckStringTie <- function(print=TRUE){
  if (print){
    cat(paste0("     \u25CF  Checking stringtie command\n"))
  }
  stringtie.installed <- system( 'stringtie --version', ignore.stdout = !print, ignore.stderr = !print)==0
  if( isTRUE(stringtie.installed)){
    if(isTRUE(print)){
      cat("(\u2714) : 'stringtie' is installed\n\n")
    }
    return(TRUE)
  }
  else{
    cat("(\u2718) : 'stringtie' command is not found on this device. Please run 'InstallAll()' to install the necessary programs or 'ExportPath' to update the path.\n\n")
    return(FALSE)
  }
}

#' Check Gffcompare
#'
#' Check whether Gffcompare is installed on the workstation
#'
CheckGffcompare <- function(print=TRUE) {
  if(print) {
    cat(paste0("     \u25CF  Checking gffcompare command\n"))
  }
  gffcompare.old <- system( 'gffcompare --version', ignore.stdout = !print, ignore.stderr = !print)==0
  if( isTRUE(gffcompare.old)){
    if(isTRUE(print)){
      cat("(\u2714) : 'gffcompare' is installed\n\n")
    }
    return(TRUE)
  }
  else{
    cat("(\u2718) : \'gffcompare\' command is not found on this device. Please run 'InstallAll()' to install the necessary programs or 'ExportPath' to update the path.\n\n")
    return(FALSE)
  }
}

#' Check Samtools
#'
#' Check whether Samtools is installed on the workstation
#'
CheckSamtools <- function(print=TRUE){
  if (print) {
    cat(paste0("     \u25CF  Checking samtools command\n"))
  }
  samtools.old <- system( 'samtools --version', ignore.stdout = !print, ignore.stderr = !print)==0
  if( isTRUE(samtools.old)){
    if(isTRUE(print)){
      cat("(\u2714) : 'samtools' is installed\n\n")
    }
    return(TRUE)
  }
  else{
    cat("(\u2718) : \'samtools\' command is not found on this device. Please run 'InstallAll()' to install the necessary programs or 'ExportPath' to update the path.\n\n")
    return(FALSE)
  }
}

#' Check necessary programs for this pipeline
#'
#' Check whether Hisat2, Stringtie, Samtools, Gffcompare are installed on the workstation
#'
#' @export
CheckToolAll <- function(print=TRUE) {
  cat("************** Checking Availability of Commands ************\n")
  hisat2.check <- CheckHisat2(print=print)
  stringtie.check <- CheckStringTie(print=print)
  gff.check <- CheckGffcompare(print=print)
  samtool.check <- CheckSamtools(print=print)
  if (isTRUE(hisat2.check) && isTRUE(stringtie.check) && isTRUE(gff.check) && isTRUE(samtool.check)){
    return(TRUE)
  } else {
    stop("(\u2718) Necessary program is missing.\n     1. Check 'INSTALL_TOOLS.Rout' whether tools are properly installed.\n     2. Run 'ExportPath()' to set the environment.\n\n")
    return(FALSE)
  }
}


PreRNAseqEnvironmentSet <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNAseqEnvironmentSet()' environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNAseqEnvironmentSet environment ERROR")
  }
  cat("     (\u2714) : RNAseqEnvironmentSet pre-check is valid\n\n")
}

PostRNAseqEnvironmentSet <- function(path.prefix, sample.pattern) {
  cat("\u269C\u265C\u265C\u265C RNAseqEnvironmentSet()' environment post-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  chrX.gtf <- file.exists(paste0(path.prefix, "gene_data/ref_genes/chrX.gtf"))
  chrX.fa <- file.exists(paste0(path.prefix, "gene_data/ref_genome/chrX.fa"))
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  check.tool.result <- CheckToolAll()
  validity <- phenodata.csv && chrX.gtf && chrX.fa && check.tool.result && (length(raw.fastq) != 0)
  if (!isTRUE(validity)) {
    stop("RNAseqQualityTrimming() post-check ERROR")
  }
  cat("     (\u2714) : RNAseqQualityTrimming() post-check is valid\n\n")
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
  cat(paste0("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n"))
}
