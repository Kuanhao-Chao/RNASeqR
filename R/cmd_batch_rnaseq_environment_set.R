#' @title RNASeqEnvironmentSet_CMD
#'
#' @description
#'   Set up the environment for the following RNA-Seq workflow in background.\cr
#'   This function do 4 things :\cr
#'   \enumerate{
#'     \item Create file directories.\cr
#'     \item Install necessary tools. \cr
#'     \item Export 'RNASeq_bin/' to the R environment. \cr
#'     \item Check command of tools. \cr
#'   }
#'   First it will create 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/',
#'   'Rscript/', 'Rscript_out/' directories. \cr Afterwards, 'Hisat2',
#'   'Stringtie', 'Gffcompare' will be installed under
#'   'RNASeq_bin/Download/' and be unpacked under 'RNASeq_bin/Unpacked/'. \cr
#'   'RNASeq_bin/' will be added to the R environment and
#'   validity of tools will be checked.\cr
#'   Any ERROR occurs will be reported and the program will be terminated.\cr
#'   If you want to set up the environment for the following RNA-Seq workflow
#'   in R shell, please see \code{RNASeqEnvironmentSet()} function.
#'
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param install.hisat2 Whether to install 'HISAT2' in this function step.
#'   Default value is\code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2' installation.
#' @param install.STAR Whether to install 'STAR' in this function step.
#'   Default value is\code{TRUE}.
#'   Set \code{FALSE} to skip 'STAR' installation.
#' @param install.stringtie Whether to install 'StringTie'
#'   in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'StringTie' installation.
#' @param install.gffcompare Whether to install 'Gffcompare'
#'   in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'Gffcompare' installation.
#' @param run Default value is \code{TRUE}. If \code{TRUE},
#'   'Rscript/Environment_Set.R' will be created and executed.
#'   The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{False}, 'Rscript/Environment_Set.R' will be
#'   created without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE},
#'   the result of checking \code{RNASeqRParam} will be reported in
#'   'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqEnvironmentSet_CMD(yeast)}
RNASeqEnvironmentSet_CMD <- function(RNASeqRParam,
                                     install.hisat2     = TRUE,
                                     install.STAR       = TRUE,
                                     install.stringtie  = TRUE,
                                     install.gffcompare = TRUE,
                                     run                = TRUE,
                                     check.s4.print     = TRUE) {
  # check input param
  which.s4.object <- CheckS4Object_All(RNASeqRParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  # Remove all the files and restart !!
  unlink(paste0(RNASeqRParam@path.prefix, "*"), recursive = TRUE)
  # Create the main directory for RNA-Seq analysis
  MkdirAll(RNASeqRParam@path.prefix)
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  INSIDE.path.prefix <- "@"(RNASeqRParam, path.prefix)
  saveRDS(RNASeqRParam,
          file = paste0(INSIDE.path.prefix,
                        "gene_data/RNASeqRParam.rds"))
  fileConn <- file(paste0(INSIDE.path.prefix, "Rscript/Environment_Set.R"))
  first <- "library(RNASeqR)"
  if (which.s4.object == "RNASeqRParam") {
    second <- paste0("RNASeqEnvironmentSet(RNASeqRParam = 'INSIDE'",
                     ", which.trigger = 'INSIDE'",
                     ", INSIDE.path.prefix = '", INSIDE.path.prefix,
                     "', install.hisat2 = ", install.hisat2,
                     ", install.STAR = ", install.STAR,
                     ", install.stringtie = ", install.stringtie,
                     ", install.gffcompare = ", install.gffcompare,")")
  } else if (which.s4.object == "RNASeqRParam_Sam" || which.s4.object == "RNASeqRParam_Bam") {
    second <- paste0("RNASeqEnvironmentSet(RNASeqRParam = 'INSIDE'",
                     ", which.trigger = 'INSIDE'",
                     ", INSIDE.path.prefix = '", INSIDE.path.prefix,
                     "', install.hisat2 = ", "FALSE",
                     ", install.STAR = ", "FALSE",
                     ", install.stringtie = ", install.stringtie,
                     ", install.gffcompare = ", install.gffcompare,")")
  }
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message("\u2605 '", path.prefix,
          "Rscript/Environment_Set.R' has been created.\n")
  if (run) {
    R.home.lib <- R.home()
    R.home.bin <- gsub("/lib/R", "/bin/R", R.home.lib)
    system2(command = "nohup",
            args = paste0(R.home.bin, " CMD BATCH ",
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

#' @title RNASeqEnvironmentSet
#'
#' @description
#'   Set up the environment for the following RNA-Seq workflow in R shell\cr
#'   This function do 4 things :\cr
#'   \enumerate{
#'     \item Create file directories.\cr
#'     \item Install necessary tools. \cr
#'     \item Export 'RNASeq_bin/' to the R environment. \cr
#'     \item Check command of tools. \cr
#'   }
#'   First it will create 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/',
#'   'Rscript/', 'Rscript_out/' directories. \cr Afterwards, 'Hisat2',
#'   'Stringtie', 'Gffcompare' will be installed under
#'   'RNASeq_bin/Download/' and be unpacked under 'RNASeq_bin/Unpacked/'. \cr
#'   'RNASeq_bin/' will be added to the R environment and
#'   validity of tools will be checked.\cr
#'   Any ERROR occurs will be reported and the program will be terminated.\cr
#'   If you want to set up the environment for the following RNA-Seq workflow
#'   in background, please see \code{RNASeqEnvironmentSet_CMD()} function.
#'
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param which.trigger Default value is \code{OUTSIDE}. User should not change
#'   this value.
#' @param INSIDE.path.prefix Default value is \code{NA}. User should not change
#'   this value.
#' @param install.hisat2 Whether to install 'HISAT2' in this function step.
#'   Default value is\code{TRUE}.
#'   Set \code{FALSE} to skip 'HISAT2' installation.
#' @param install.STAR Whether to install 'STAR' in this function step.
#'   Default value is\code{TRUE}.
#'   Set \code{FALSE} to skip 'STAR' installation.
#' @param install.stringtie Whether to install 'StringTie'
#'   in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'StringTie' installation.
#' @param install.gffcompare Whether to install 'Gffcompare'
#'   in this function step. Default value is \code{TRUE}.
#'   Set\code{FALSE} to skip 'Gffcompare' installation.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE},
#'   the result of checking \code{RNASeqRParam} will be reported in
#'   'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'.
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqEnvironmentSet(RNASeqRParam = yeast)}
RNASeqEnvironmentSet <- function(RNASeqRParam,
                                 which.trigger      = "OUTSIDE",
                                 INSIDE.path.prefix = NA,
                                 install.hisat2     = TRUE,
                                 install.STAR       = TRUE,
                                 install.stringtie  = TRUE,
                                 install.gffcompare = TRUE,
                                 check.s4.print     = TRUE) {
  CheckOperatingSystem(FALSE)
  # If `which.trigger` is OUTSIDE, then directory must be built
  # If `which.trigger` is INSIDE, then directory must not be
  #  built here(will created in CMD)
  if (isS4(RNASeqRParam) &
      which.trigger == "OUTSIDE" &
      is.na(INSIDE.path.prefix)) {
    # Remove all the files and restart !!
    unlink(paste0(RNASeqRParam@path.prefix, "*"), recursive = TRUE)
    # This is an external call!!
    MkdirAll(RNASeqRParam@path.prefix)
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
  os.type <- "@"(RNASeqRParam, os.type)
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  input.path.prefix <- "@"(RNASeqRParam, input.path.prefix)
  genome.name <- "@"(RNASeqRParam, genome.name)
  sample.pattern <- "@"(RNASeqRParam, sample.pattern)
  PreRNASeqEnvironmentSet(path.prefix, sample.pattern)
  if (which.s4.object == "RNASeqRParam") {
    indices.optional <- "@"(RNASeqRParam, indices.optional)
    CopyInputDir(path.prefix,
                 input.path.prefix,
                 genome.name,
                 sample.pattern,
                 indices.optional)
    InstallAll(path.prefix,
               os.type,
               install.hisat2,
               install.STAR,
               install.stringtie,
               install.gffcompare)
    ExportPath(path.prefix)
    PostRNASeqEnvironmentSet(path.prefix,
                             genome.name,
                             sample.pattern)
  } else if (which.s4.object == "RNASeqRParam_Sam") {
    CopyInputDir_Sam(path.prefix,
                     input.path.prefix,
                     genome.name,
                     sample.pattern)
    InstallAll(path.prefix,
               os.type,
               FALSE,
               FALSE,
               install.stringtie,
               install.gffcompare)
    ExportPath(path.prefix)
    PostRNASeqEnvironmentSet_Sam(path.prefix,
                                 genome.name,
                                 sample.pattern)
  }  else if (which.s4.object == "RNASeqRParam_Bam") {
    CopyInputDir_Bam(path.prefix,
                     input.path.prefix,
                     genome.name,
                     sample.pattern)
    InstallAll(path.prefix,
               os.type,
               FALSE,
               FALSE,
               install.stringtie,
               install.gffcompare)
    ExportPath(path.prefix)
    PostRNASeqEnvironmentSet_Bam(path.prefix,
                                 genome.name,
                                 sample.pattern)
  }

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
  bool.path.prefix <- dir.exists(path.prefix)
  if (bool.path.prefix) {
    MkdirGeneDir(path.prefix)
    MkdirRNASeq_bin(path.prefix)
    MkdirRNASeq_results(path.prefix)
    MkdirRscript_Rscript_out(path.prefix)
  } else {
    message("'", path.prefix, "' has not been created yet!!\n\n")
  }
}

# inner function : Copy input files directory
CopyInputDir <- function(path.prefix,
                         input.path.prefix,
                         genome.name,
                         sample.pattern,
                         indices.optional) {
  message("************** Directory Copying ************\n")
  message("     \u25CF Copying From :",input.path.prefix,
          "input_files/", genome.name, ".gtf\n")
  gtf.des <- paste0(path.prefix, "gene_data/ref_genes/", genome.name, ".gtf")
  file.remove(gtf.des)
  file.symlink(paste0(input.path.prefix, "input_files/", genome.name, ".gtf"),
               gtf.des)
  message("     \u25CF           To :", gtf.des, "\n")
  fa.des <- paste0(path.prefix, "gene_data/ref_genome/", genome.name, ".fa")
  message("     \u25CF Copying From :",
          input.path.prefix, "input_files/",
          genome.name, ".fa",  "\n")
  file.remove(fa.des)
  file.symlink(paste0(input.path.prefix, "input_files/", genome.name, ".fa"),
               fa.des)
  message("     \u25CF           To :", fa.des, "\n")
  message("     \u25CF Copying From :",
          input.path.prefix, "input_files/", "raw_fastq.gz/",  "\n")
  unlink(paste0(path.prefix, "gene_data/raw_fastq.gz/"), recursive = TRUE)
  dir.create(paste0(path.prefix, "gene_data/raw_fastq.gz/"))
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
  message("     \u25CF           To :",
          path.prefix, "gene_data/raw_fastq.gz/\n")
  message("     \u25CF Copying From :",
          input.path.prefix, "input_files/phenodata.csv\n")
  pheno.des <- paste0(path.prefix, "gene_data/phenodata.csv")
  file.remove(pheno.des)
  file.symlink(paste0(input.path.prefix, "input_files/phenodata.csv"),
               pheno.des)
  message("     \u25CF           To :", pheno.des, "\n")
  if (isTRUE(indices.optional)) {
    message("     \u25CF Copying From :",
            input.path.prefix, "input_files/", "indices/\n")
    unlink(paste0(path.prefix, "gene_data/indices/"), recursive = TRUE)
    dir.create(paste0(path.prefix, "gene_data/indices/"))
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


# inner function : Copy input files directory
CopyInputDir_Sam <- function(path.prefix,
                             input.path.prefix,
                             genome.name,
                             sample.pattern) {
  message("************** Directory Copying ************\n")
  message("     \u25CF Copying From :",input.path.prefix,
          "input_files/", genome.name, ".gtf\n")
  gtf.des <- paste0(path.prefix, "gene_data/ref_genes/", genome.name, ".gtf")
  file.remove(gtf.des)
  file.symlink(paste0(input.path.prefix, "input_files/", genome.name, ".gtf"),
               gtf.des)
  message("     \u25CF           To :", gtf.des, "\n")
  message("     \u25CF Copying From :",
          input.path.prefix, "input_files/", "raw_sam/",  "\n")
  unlink(paste0(path.prefix, "gene_data/raw_sam/"), recursive = TRUE)
  dir.create(paste0(path.prefix, "gene_data/raw_sam/"))
  raw_sam.subfiles <- list.files(path = paste0(input.path.prefix,
                                               "input_files/",
                                               "raw_sam"),
                                 pattern = sample.pattern,
                                 recursive = TRUE,
                                 full.names = TRUE)
  vapply(raw_sam.subfiles,
         function(x) file.symlink(x, paste0(path.prefix,
                                            "gene_data/raw_sam")),
         FUN.VALUE = TRUE)
  message("     \u25CF           To :",
          path.prefix, "gene_data/raw_sam/\n")
  message("     \u25CF Copying From :",
          input.path.prefix, "input_files/phenodata.csv\n")
  pheno.des <- paste0(path.prefix, "gene_data/phenodata.csv")
  file.remove(pheno.des)
  file.symlink(paste0(input.path.prefix, "input_files/phenodata.csv"),
               pheno.des)
  message("     \u25CF           To :", pheno.des, "\n")
}


# inner function : Copy input files directory
CopyInputDir_Bam <- function(path.prefix,
                             input.path.prefix,
                             genome.name,
                             sample.pattern) {
  message("************** Directory Copying ************\n")
  message("     \u25CF Copying From :",input.path.prefix,
          "input_files/", genome.name, ".gtf\n")
  gtf.des <- paste0(path.prefix, "gene_data/ref_genes/", genome.name, ".gtf")
  file.remove(gtf.des)
  file.symlink(paste0(input.path.prefix, "input_files/", genome.name, ".gtf"),
               gtf.des)
  message("     \u25CF           To :", gtf.des, "\n")
  message("     \u25CF Copying From :",
          input.path.prefix, "input_files/", "raw_bam/",  "\n")
  unlink(paste0(path.prefix, "gene_data/raw_bam/"), recursive = TRUE)
  dir.create(paste0(path.prefix, "gene_data/raw_bam/"))
  raw_bam.subfiles <- list.files(path = paste0(input.path.prefix,
                                               "input_files/",
                                               "raw_bam"),
                                 pattern = sample.pattern,
                                 recursive = TRUE,
                                 full.names = TRUE)
  vapply(raw_bam.subfiles,
         function(x) file.symlink(x, paste0(path.prefix,
                                            "gene_data/raw_bam")),
         FUN.VALUE = TRUE)
  message("     \u25CF           To :",
          path.prefix, "gene_data/raw_bam/\n")
  message("     \u25CF Copying From :",
          input.path.prefix, "input_files/phenodata.csv\n")
  pheno.des <- paste0(path.prefix, "gene_data/phenodata.csv")
  file.remove(pheno.des)
  file.symlink(paste0(input.path.prefix, "input_files/phenodata.csv"),
               pheno.des)
  message("     \u25CF           To :", pheno.des, "\n")
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
    os.file.name.zip <- "hisat2-2.1.0-OSX_x86_64.zip"
    os.file.name <- "hisat2-2.1.0"
    url <- paste0(url, os.file.name.zip)
  } else if (os == "windows"){
    stop("Hisat2 is not supporting windows.\n\n")
    return(FALSE)
  } else {
    stop("Unknow operating system.\n\n")
    return(FALSE)
  }
  message("************** Installing Hisat2 (",
          os.file.name.zip, ") ************\n")
  file.download <- getURL(url, download.file, paste0(path.prefix,
                                                     "RNASeq_bin/Download/",
                                                     os.file.name.zip))
  if (file.download == 0) {
    message("HISAT2 binaries were downloaded successfully!!\n")
  }
  message("\n************** Unpacking Hisat2 (",
          os.file.name.zip, ") ************\n")
  if (dir.exists(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))) {
    unlink(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name),
           recursive = TRUE)
  }

  # utils::unzip(zipfile = paste0(path.prefix, "RNASeq_bin/Download/",
  #                               os.file.name.zip),
  #              exdir = paste0(path.prefix,
  #                             "RNASeq_bin/Unpacked/"),
  #              overwrite = TRUE)
  main.command <- "unzip"
  command.result <- system2(command = main.command,
                            args = paste0(path.prefix, "RNASeq_bin/Download/",
                                          os.file.name.zip, " -d ", path.prefix,
                                          "RNASeq_bin/Unpacked/"))
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop("'", main.command, "' ERROR")
  }
  message("hisat2 binaries\n\n")
  message("\n************** Moving Hisat2 Binary ************")

  file.target.hisat <- list.files(path = paste0(path.prefix,
                                                "RNASeq_bin/Unpacked/",
                                                os.file.name),
                                  pattern = "hisat2*",
                                  full.names = TRUE)
  file.target.py <- list.files(path = paste0(path.prefix,
                                             "RNASeq_bin/Unpacked/",
                                             os.file.name),
                               pattern = "*.py",
                               full.names = TRUE)
  file.target <- unique(c(file.target.hisat, file.target.py))
  file.move <- file.copy(from = file.target,
                         to = paste0(path.prefix, "RNASeq_bin/"),
                         overwrite = TRUE)
  if (command.result != 0 ) {
    message("(\u2718) 'file.copy()' is failed !!")
    stop("'file.copy()' ERROR")
  }
  message("\n'", path.prefix, "RNASeq_bin/Download/",
          os.file.name.zip, "' has been installed.\n")
  message("Hisat2 has been unpacked. ('", path.prefix,
          "RNASeq_bin/Unpacked/", os.file.name, "/')", "\n\n")
  return(TRUE)
}

InstallStarBianry <- function(path.prefix, os.type){
  os <- os.type
  url <- "https://github.com/alexdobin/STAR/archive/"
  if (os == "linux"){
    os.file.name.zip <- "star-2.7.1a.zip"
    os.file.name <- "start-2.7.1a"
    url <- paste0(url, "2.7.1a.zip")
  } else if (os == "osx"){
    os.file.name.zip <- "star-2.7.1a.zip"
    os.file.name <- "start-2.7.1a"
    url <- paste0(url, "2.7.1a.zip")
  } else if (os == "windows"){
    stop("STAR is not supporting windows.\n\n")
    return(FALSE)
  } else {
    stop("Unknow operating system.\n\n")
    return(FALSE)
  }
  message("************** Installing STAR (",
          os.file.name.zip, ") ************\n")
  file.download <- getURL(url, download.file, paste0(path.prefix,
                                                     "RNASeq_bin/Download/",
                                                     os.file.name.zip))
  if (file.download == 0) {
    message("STAR binaries were downloaded successfully!!\n")
  }
  message("\n************** Unpacking STAR (",
          os.file.name.zip, ") ************\n")
  if (dir.exists(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))) {
    unlink(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name),
           recursive = TRUE)
  }

  # utils::unzip(zipfile = paste0(path.prefix, "RNASeq_bin/Download/",
  #                               os.file.name.zip),
  #              exdir = paste0(path.prefix,
  #                             "RNASeq_bin/Unpacked/"),
  #              overwrite = TRUE)
  main.command <- "unzip"
  command.result <- system2(command = main.command,
                            args = paste0(path.prefix, "RNASeq_bin/Download/",
                                          os.file.name.zip, " -d ", path.prefix,
                                          "RNASeq_bin/Unpacked/"))
  if (command.result != 0 ) {
    message("(\u2718) '", main.command, "' is failed !!")
    stop("'", main.command, "' ERROR")
  }
  message("STAR binaries\n\n")
  message("\n************** Moving STAR Binary ************")
  if (os == "osx") {
    file.target.star <- list.files(path = paste0(path.prefix,
                                                 "RNASeq_bin/Unpacked/",
                                                 "STAR-2.7.1a/bin/",
                                                 "MacOSX_x86_64"),
                                   pattern = "STAR*",
                                   full.names = TRUE)
  } else if (os == "linux") {
    file.target.star <- list.files(path = paste0(path.prefix,
                                                 "RNASeq_bin/Unpacked/",
                                                 "STAR-2.7.1a/bin/",
                                                 "Linux_x86_64"),
                                   pattern = "STAR*",
                                   full.names = TRUE)
  }
  file.move <- file.copy(from = file.target.star,
                         to = paste0(path.prefix, "RNASeq_bin/"),
                         overwrite = TRUE)
  if (command.result != 0 ) {
    message("(\u2718) 'file.copy()' is failed !!")
    stop("'file.copy()' ERROR")
  }
  message("\n'", path.prefix, "RNASeq_bin/Download/",
          os.file.name.zip, "' has been installed.\n")
  message("STAR has been unpacked. ('", path.prefix,
          "RNASeq_bin/Unpacked/STAR-2.7.1a/')", "\n\n")
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
    os.file.name.zip <- "stringtie-1.3.4d.OSX_x86_64.tar.gz"
    os.file.name <- "stringtie-1.3.4d.OSX_x86_64"
    url <- paste0(url, os.file.name.zip)
  } else if (os == "windows"){
    stop("Stringtie is not supporting windows.\n\n")
    return(FALSE)
  } else {
    stop("Unknow operating system.\n\n")
    return(FALSE)
  }
  message("************** Installing stringtie (",
          os.file.name.zip, ") ************\n")
  file.download <- getURL(url, download.file, paste0(path.prefix,
                                                     "RNASeq_bin/Download/",
                                                     os.file.name.zip))
  if (file.download == 0) {
    message("StringTie binaries were downloaded successfully!!\n")
  }
  message("\n************** Unpacking stringtie (",
          os.file.name.zip, ") ************\n")
  # If file directory exists!! Remove it first!!
  if (dir.exists(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))) {
    unlink(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name),
           recursive = TRUE)
  }
  file.unpacked <- utils::untar(tarfile = paste0(path.prefix,
                                                 "RNASeq_bin/Download/",
                                                 os.file.name.zip),
                                exdir = paste0(path.prefix,
                                               "RNASeq_bin/Unpacked/"),
                                compressed = "gzip")
  if (file.unpacked != 0 ) {
    message("(\u2718) 'utils::untar()' is failed !!")
    stop("'utils::untar()' ERROR")
  }
  message("\n************** Moving stringtie Binary ************")
  file.move <- file.copy(from = paste0(path.prefix, "RNASeq_bin/Unpacked/",
                                       os.file.name, "/stringtie"),
                         to = paste0(path.prefix, "RNASeq_bin/"),
                         overwrite = TRUE)
  if (!file.move) {
    message("(\u2718) 'file.copy()' is failed !!")
    stop("'file.copy()' ERROR")
  }
  message("\n'", path.prefix,
          "RNASeq_bin/Download/",
          os.file.name.zip, "' has been installed.\n")
  message("StringTie has been unpacked. ('",
          path.prefix, "RNASeq_bin/Unpacked/",
          os.file.name, "')\n\n")
  return(TRUE)
}

# Install Gffcompare binary
InstallGffcompareBinary <- function(path.prefix, os.type){
  os <- os.type
  url <- "http://ccb.jhu.edu/software/stringtie/dl/"
  if (os == "linux"){
    os.file.name.zip <- "gffcompare-0.10.4.Linux_x86_64.tar.gz"
    os.file.name <- "gffcompare-0.10.4.Linux_x86_64"
    url <- paste0(url, os.file.name.zip)
  } else if (os == "osx"){
    os.file.name.zip <- "gffcompare-0.10.4.OSX_x86_64.tar.gz"
    os.file.name <- "gffcompare-0.10.4.OSX_x86_64"
    url <- paste0(url, os.file.name.zip)
  } else if (os == "windows"){
    stop("Gffcompare is not supporting windows.\n\n")
    return(FALSE)
  } else {
    stop("Unknow operating system.\n\n")
    return(FALSE)
  }
  message("************** Installing gffcompare (",
          os.file.name.zip, ") ************\n")
  file.download <- getURL(url, download.file, paste0(path.prefix,
                                                     "RNASeq_bin/Download/",
                                                     os.file.name.zip))
  if (file.download == 0) {
    message("Gffcompare binaries were downloaded successfully!!\n")
  }
  message("\n************** Unpacking gffcompare (",
          os.file.name.zip, ") ************\n")
  # If file directory exists!! Remove it first!!
  if (dir.exists(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name))) {
    unlink(paste0(path.prefix, "RNASeq_bin/Unpacked/", os.file.name),
           recursive = TRUE)
  }
  file.unpacked <- utils::untar(tarfile = paste0(path.prefix,
                                                 "RNASeq_bin/Download/",
                                                 os.file.name.zip),
                                exdir = paste0(path.prefix,
                                               "RNASeq_bin/Unpacked/"),
                                compressed = "gzip")
  if (file.unpacked != 0 ) {
    message("(\u2718) 'utils::untar()' is failed !!")
    stop("'utils::untar()' ERROR")
  }
  message("\n************** Moving gffcompare Binary ************")
  file.move <- file.copy(from = paste0(path.prefix, "RNASeq_bin/Unpacked/",
                                       os.file.name, "/gffcompare"),
                         to = paste0(path.prefix, "RNASeq_bin/"),
                         overwrite = TRUE)
  if (!file.move) {
    message("(\u2718) 'file.copy()' is failed !!")
    stop("'file.copy()' ERROR")
  }
  message("\n'", path.prefix, "RNASeq_bin/Download/",
          os.file.name.zip, "' has been installed.\n")
  message("Gffcompare has been unpacked. ('", path.prefix,
          "RNASeq_bin/Unpacked/", os.file.name, "')", "\n\n")
  return(TRUE)
}

# Install 'Hisat2', 'StringTie', 'Gffcompare'
InstallAll <- function(path.prefix,
                       os.type,
                       install.hisat2,
                       install.STAR,
                       install.stringtie,
                       install.gffcompare) {
  message("\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618  ",
          "Start installing ... ",
          "\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618\n")
  install.software <- ""
  if (install.hisat2) {
    install.software <- paste0(install.software, "\u25CF'hisat2' ")
  }
  if (install.STAR) {
    install.software <- paste0(install.software, "\u25CF'STAR' ")
  }
  if (install.stringtie) {
    install.software <- paste0(install.software, "\u25CF'stringtie' ")
  }
  if (install.gffcompare) {
    install.software <- paste0(install.software, "\u25CF'gffcompare' ")
  }
  if (!install.hisat2 & !install.stringtie &
      !install.gffcompare) {
    message("   \u261E\u261E skipping installation process ... \n")
  } else {
    message("   \u261E\u261E ",
            install.software,
            "will be installed. ... \n")
    message("   \u261E\u261E  Compressed files will be in '",
            path.prefix, "RNASeq_bin/Download/'\n")
    message("   \u261E\u261E  Unpacked files will be in '",
            path.prefix, "RNASeq_bin/Unpacked/'\n")
    message("   \u261E\u261E  Binary files will be copied to '",
            path.prefix, "RNASeq_bin/'", "\n\n")
  }
  if (install.hisat2) {
    message("\u2618\u2618 Hisat2 processing ...\n")
    InstallHisat2Bianry(path.prefix, os.type)
  }
  if (install.STAR) {
    message("\u2618\u2618 STAR processing ...\n")
    InstallStarBianry(path.prefix, os.type)
  }
  if (install.stringtie) {
    message("\u2618\u2618 Stringtie processing ...\n")
    InstallStringTieBinary(path.prefix, os.type)
  }
  if (install.gffcompare) {
    message("\u2618\u2618 Gffcompare processing ...\n")
    InstallGffcompareBinary(path.prefix, os.type)
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
    message("(\u2718) : 'hisat2' command is not found on this device!! \n\n")
    return(FALSE)
  }
}

# Check 'STAR'
CheckSTAR <- function(print=TRUE){
  if (print) {
    message("\u25CF  Checking STAR command\n")
  }
  STAR.installed <- system("STAR --version",
                             ignore.stdout = !print,
                             ignore.stderr = !print) == 0
  if (isTRUE(STAR.installed)){
    if (isTRUE(print)){
      message("(\u2714) : 'STAR' is installed\n\n")
    }
    return(TRUE)
  } else {
    message("(\u2718) : 'STAR' command is not found on this device!! \n\n")
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
    message("(\u2718) : 'stringtie' command is not found on this device!! \n\n")
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
    message("(\u2718) : \'gffcompare\' ",
            "command is not found on this device!!\n\n")
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
      message("(\u2714) : 'samtools' is installed\n")
    }
    return(TRUE)
  }
  else{
    message("(\u2718) : \'samtools\' command is not found on this device!!\n")
    return(FALSE)
  }
}

#' @title CheckToolAll
#'
#' @description Check whether 'Hisat2', 'Stringtie' and 'Gffcompare'
#'   are installed on the workstation
#'
#' @param print If \code{TRUE}, detailed information will be printed.
#'   If \code{FALSE}, detailed information will not be printed.
#'
#' @param path.prefix path prefix of 'gene_data/', 'RNASeq_bin/',
#'   'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories.
#'
#' @return None
#' @export
#' @examples
#' data(yeast)
#' \dontrun{
#' CheckToolAll(yeast@@path.prefix,
#'              print=TRUE)}
CheckToolAll <- function(path.prefix, print=TRUE) {
  message("************** Checking Availability of Commands ************\n")
  ExportPath(path.prefix)
  hisat2.check <- CheckHisat2(print)
  STAR.check <- CheckSTAR(print)
  stringtie.check <- CheckStringTie(print)
  gff.check <- CheckGffcompare(print)
  if (isTRUE(hisat2.check) &&
      isTRUE(stringtie.check) &&
      isTRUE(gff.check) &&
      isTRUE(STAR.check)){
    return(TRUE)
  } else {
    stop("(\u2718) Necessary program is missing.\n     ",
         "1. Check 'INSTALL_TOOLS.Rout' whether tools are ",
         "properly installed.\n     ",
         "2. Run 'ExportPath()' to set the environment.\n\n")
    return(FALSE)
  }
}


CheckTool_Sam_Bam <- function(path.prefix, print=TRUE) {
  message("************** Checking Availability of Commands ************\n")
  ExportPath(path.prefix)
  stringtie.check <- CheckStringTie(print)
  gff.check <- CheckGffcompare(print)
  if (isTRUE(gff.check) &&
      isTRUE(gff.check)){
    return(TRUE)
  } else {
    stop("(\u2718) Necessary program is missing.\n     ",
         "1. Check 'INSTALL_TOOLS.Rout' whether tools are ",
         "properly installed.\n     ",
         "2. Run 'ExportPath()' to set the environment.\n\n")
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
  check.tool.result <- CheckToolAll(path.prefix)
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

PostRNASeqEnvironmentSet_Sam <- function(path.prefix,
                                         genome.name,
                                         sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqEnvironmentSet()' ",
          "environment post-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  ref.gtf <- file.exists(paste0(path.prefix,
                                "gene_data/ref_genes/",
                                genome.name, ".gtf"))
  raw.sam <- list.files(path = paste0(path.prefix, "gene_data/raw_sam/"),
                        pattern = sample.pattern,
                        all.files = FALSE,
                        full.names = FALSE,
                        recursive = FALSE,
                        ignore.case = FALSE)
  check.tool.result <- CheckTool_Sam_Bam(path.prefix)
  validity <- phenodata.csv && ref.gtf &&
    check.tool.result && (length(raw.sam) != 0)
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

PostRNASeqEnvironmentSet_Bam <- function(path.prefix,
                                         genome.name,
                                         sample.pattern) {
  message("\u269C\u265C\u265C\u265C RNASeqEnvironmentSet()' ",
          "environment post-check ...\n")
  phenodata.csv <- file.exists(paste0(path.prefix, "gene_data/phenodata.csv"))
  ref.gtf <- file.exists(paste0(path.prefix,
                                "gene_data/ref_genes/",
                                genome.name, ".gtf"))
  raw.bam <- list.files(path = paste0(path.prefix, "gene_data/raw_bam/"),
                        pattern = sample.pattern,
                        all.files = FALSE,
                        full.names = FALSE,
                        recursive = FALSE,
                        ignore.case = FALSE)
  check.tool.result <- CheckTool_Sam_Bam(path.prefix)
  validity <- phenodata.csv && ref.gtf &&
    check.tool.result && (length(raw.bam) != 0)
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

#
CheckS4Object_All <- function(RNASeqRParam, print = TRUE) {
  if (isS4(RNASeqRParam) &&
      class(RNASeqRParam)[1] == "RNASeqRParam") {
    if (print) {
      message("************** Checking validity of S4 input ************\n")
      message("     (\u2714) : input is valid ",
              "'RNASeqRParam' instance! \n\n")
    }
    return("RNASeqRParam")
  } else if (isS4(RNASeqRParam) &&
             class(RNASeqRParam)[1] == "RNASeqRParam_Sam") {
    if (print) {
      message("************** Checking validity of S4 input ************\n")
      message("     (\u2714) : input is valid ",
              "'RNASeqRParam_Sam' instance! \n\n")
    }
    return("RNASeqRParam_Sam")
  } else if (isS4(RNASeqRParam) &&
             class(RNASeqRParam)[1] == "RNASeqRParam_Bam") {
    if (print) {
      message("************** Checking validity of S4 input ************\n")
      message("     (\u2714) : input is valid ",
              "'RNASeqRParam_Bam' instance! \n\n")
    }
    return("RNASeqRParam_Bam")
  } else {
    message("(\u2718) : input is not a valid ",
            "'RNASeqRParam' or 'RNASeqRParam_Sam' instance!.\n" )
    stop("Invalid 'RNASeqRParam' or 'RNASeqRParam_Sam' input ERROR")
  }
}
