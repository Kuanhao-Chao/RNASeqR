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
RNAseqEnvironmentSet_CMD <- function(RNASeqWorkFlowParam) {
  # check input param
  os.type <- RNASeqWorkFlowParam@os.type
  path.prefix <- RNASeqWorkFlowParam@path.prefix
  input.path.prefix <- RNASeqWorkFlowParam@input.path.prefix
  gene.name <- RNASeqWorkFlowParam@gene.name
  sample.pattern <- RNASeqWorkFlowParam@sample.pattern
  indexes.optional <- RNASeqWorkFlowParam@indexes.optional
  MkdirAll(path.prefix)
  r_script.dir <- dir.create(file.path(paste0(path.prefix, 'Rscript/')), showWarnings = FALSE) == 0
  r_script.out.dir <- dir.create(file.path(paste0(path.prefix, 'Rscript_out/')), showWarnings = FALSE) == 0
  fileConn<-file(paste0(path.prefix, "Rscript/RNAseqEnvironmentSet.R"))
  first <- "library(RNASeqWorkflow)"
  second <- paste0("RNAseqEnvironmentSet(path.prefix = '", path.prefix, "', input.path.prefix = '", input.path.prefix, "', gene.name = '", gene.name, "', sample.pattern = '", sample.pattern, "', indexes.optional = ",indexes.optional, ", os.type = '", os.type, "')")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/RNAseqEnvironmentSet.R ", path.prefix, "Rscript_out/RNAseqEnvironmentSet.Rout"), stdout = "", wait = FALSE)
  cat(paste0("\u2605 Tools are installing in the background. Check current progress in '", path.prefix, "Rscript_out/RNAseqEnvironmentSet.Rout'\n\n"))
}

#' inner function : start running in the background
RNAseqEnvironmentSet_CMD_run <- function() {
  system2(command = 'nohup', args = paste0("R CMD BATCH ", path.prefix, "Rscript/RNAseqEnvironmentSet.R ", path.prefix, "Rscript_out/RNAseqEnvironmentSet.Rout"), stdout = "", wait = FALSE)
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
  CopyInputDir(path.prefix = path.prefix, input.path.prefix = input.path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, indexes.optional = indexes.optional)
  InstallAll(path.prefix = path.prefix, os.type = os.type)
  ExportPath(path.prefix = path.prefix)
  CheckToolAll()
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
#' @export
InstallAll <- function(path.prefix, os.type) {
  cat("\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618  Start installing ... \u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618\n")
  cat("   \u261E\u261E  \u25CF'hisat2', \u25CF'stringtie', \u25CF'gffcompare', \u25CF'samtools' will be installed. ... \n")
  cat(paste0("   \u261E\u261E  Compressed files will be in '", path.prefix, "RNAseq_bin/Download/'"), "\n")
  cat(paste0("   \u261E\u261E  Unpacked files will be in '", path.prefix, "RNAseq_bin/Unpacked/'"), "\n")
  cat(paste0("   \u261E\u261E  Binary files will be copied to '", path.prefix, "RNAseq_bin/'"), "\n\n")
  InstallHisat2Bianry(path.prefix, os.type)
  InstallStringTieBinary(path.prefix, os.type)
  InstallGffcompareBinary(path.prefix, os.type)
  InstallSamtoolsBinary(path.prefix, os.type)
  #return(CheckToolAll(print=FALSE))
}


#' Check 'hisat2'
#'
#' Check whether 'hisat2' is installed on the workstation
#'
#' @export
CheckHisat2 <- function(print=TRUE){
  if (print) {
    cat("************** Checking hisat2 command ************\n")
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
#' @export
CheckStringTie <- function(print=TRUE){
  if (print){
    cat("************** Checking stringtie command ************\n")
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
#' @export
CheckGffcompare <- function(print=TRUE) {
  if(print) {
    cat("************** Checking gffcompare command ************\n")
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
#' @export
CheckSamtools <- function(print=TRUE){
  if (print) {
    cat("************** Checking samtools command ************\n")
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
