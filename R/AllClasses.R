#' RNASeq workflow
#' An S4 class for checking and storing RNA-Seq workflow parameters of this package
#' @aliases RNASeq
#'
#' @slot os.type 'linux' or 'osx'. The operating system type
#' @slot python.variable A list storing python environment. \code{(check.answer, python.version)}
#' @slot python.2to3 logical value whether \code{2to3} command is available on the workstation
#' @slot path.prefix path prefix of 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @slot input.path.prefix path prefix of 'input_files/' directory
#' @slot genome.name Variable of genome name defined in this RNA-Seq workflow (ex. \code{genome.name}.fa, \code{genome.name}.gtf)
#' @slot sample.pattern  Regular expression of paired-end fastq.gz files under 'input_files/raw_fastq.gz'. Expression not includes \code{_[1,2].fastq.gz}.
#' @slot independent.variable independent variable for the biological experiment design of two-group RNA-Seq workflow
#' @slot case.group group name of the case group
#' @slot control.group group name of the control group
#' @slot indices.optional logical value whether 'indices/' is exit in 'input_files/'
#'
#' @name RNASeqWorkFlowParam-class
#'
#' @rdname RNASeqWorkFlowParam-class
#'
#' @exportClass RNASeqWorkFlowParam
#' @author Kuan-Hao Chao
#' @examples
#' \dontrun{
#' data(workflowParam)
#' class(workflowParam) #"RNASeqWorkFlowParam"
#' workflowParam@@path.prefix
#' workflowParam@@input.path.prefix
#' workflowParam@@genome.name
#' workflowParam@@sample.pattern
#' workflowParam@@independent.variable
#' workflowParam@@case.group
#' workflowParam@@control.group}
setClass("RNASeqWorkFlowParam",
         representation(
           os.type = "character",
           python.variable = "list",
           python.2to3 = "logical",
           path.prefix = "character",
           input.path.prefix = "character",
           genome.name = "character",
           sample.pattern = "character",
           independent.variable = "character",
           case.group = "character",
           control.group = "character",
           indices.optional = "logical"
         )
)


#' Constructor function for RNASeqWorkFlowParam objects
#'
#' @name RNASeqWorkFlowParam-constructor
#'
#' @param path.prefix path prefix of 'gene_data/', 'RNASeq_bin/', 'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories
#' @param input.path.prefix path prefix of 'input_files/' directory
#' @param genome.name variable of genome name defined in this RNA-Seq workflow (ex. \code{genome.name}.fa, \code{genome.name}.gtf)
#' @param sample.pattern  Regular expression of paired-end fastq.gz files under 'input_files/raw_fastq.gz'. Expression not includes \code{_[1,2].fastq.gz}.
#' @param independent.variable independent variable for the biological experiment design of two-group RNA-Seq workflow
#' @param case.group group name of the case group
#' @param control.group group name of the control group
#'
#' @return an object of class \code{RNASeqWorkFlowParam}
#'
#' @author kuan-hao Chao
#'
#' @rdname RNASeqWorkFlowParam-constructor
#'
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' \dontrun{
#' exp <- RNASeqWorkFlowParam(path.prefix = "/home/RNASeq", input.path.prefix = "/home", genome.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            independent.variable = "two.group", case.group = "treatment", control.group = "cell")
#' }
RNASeqWorkFlowParam <- function(path.prefix = NA, input.path.prefix = NA, genome.name = NA, sample.pattern = NA,
                                independent.variable = NA, case.group = NA, control.group = NA) {
  # check input parameters
  CheckInputParamNa(path.prefix, input.path.prefix, genome.name, sample.pattern,
                    independent.variable, case.group, control.group)
  # 1. check operating system
  characters.os.type <- CheckOperatingSystem()
  # 2. check python version
  python.version.list <- CheckPython()
  bool.python.avail <- python.version.list$check.answer
  numeric.python.version <- python.version.list$python.version
  two.to.three.result <- Check2to3()
  # 3. check validity of path.prefix
  bool.prefix.path <- CheckPrefixPath(path.prefix = path.prefix)
  if (bool.prefix.path){
    # add '/' to the path.prefix
    if (substr(path.prefix, nchar(path.prefix), nchar(path.prefix)) != '/') {
      path.prefix <- paste0(path.prefix, '/')
    }
  }
  # 4. check input.path.prefix
  bool.input.path.prefix <- CheckInputPrefixPath(input.path.prefix = input.path.prefix)
  if (bool.input.path.prefix){
    # add '/' to the path.prefix
    if (substr(input.path.prefix, nchar(input.path.prefix), nchar(input.path.prefix)) != '/') {
      input.path.prefix <- paste0(input.path.prefix, '/')
    }
  }
  # check sample.pattern is valid file name !!
  # 5. check sample.pattern that can't have '.fastq.gz'
  fast.gz.extend <- tools::file_ext(sample.pattern)
  if (fast.gz.extend == "gz" || fast.gz.extend == "fastq") {
    cat("(\u2718) : 'sample.pattern' can't include file extension(.gz or .fastq or .fastq.gz)\n\n")
    stop("'sample.pattern' with extension error.")
  }
  # 6. check 'input_files/' necessary files with 'genome.name', 'sample.pattern'
  input.dir.files.list <- CheckInputDirFiles(input.path.prefix = input.path.prefix, genome.name = genome.name, sample.pattern = sample.pattern)
  bool.input.dir.files <- input.dir.files.list$check.answer
  # This determine whether to run 'CreateHisat2Index'
  bool.input.dir.indices <- input.dir.files.list$optional.indices.bool
  # below still need to fix
  # 7. check 'phenodata'
  bool.phenodata <- CheckPhenodata(input.path.prefix = input.path.prefix, genome.name = genome.name, sample.pattern = sample.pattern, independent.variable = independent.variable)
  # 8. check 'case.group' and 'control.group'
  bool.check.control.control.group <- CheckCaseControlGroup(input.path.prefix = input.path.prefix, independent.variable = independent.variable, case.group = case.group, control.group = control.group)

  if ((characters.os.type == "linux" || characters.os.type == "osx") && bool.python.avail && bool.prefix.path &&
      bool.input.path.prefix && bool.input.dir.files && bool.phenodata && bool.check.control.control.group) {
    cat(paste0("\n**************************************\n"))
    cat(paste0("************** Success! **************\n"))
    cat(paste0("**************************************\n"))
    new("RNASeqWorkFlowParam",os.type = characters.os.type, python.variable = python.version.list, python.2to3 = two.to.three.result,
        path.prefix = path.prefix, input.path.prefix = input.path.prefix, genome.name = genome.name, sample.pattern = sample.pattern,
        independent.variable = independent.variable, case.group = case.group, control.group = control.group,
        indices.optional = bool.input.dir.indices)
  }
}

# inner function : check whether input values are NA
CheckInputParamNa <- function(path.prefix, input.path.prefix, genome.name, sample.pattern,
                            independent.variable, case.group, control.group) {
  cat(c("************** Checking input parameters ************\n"))
  if (is.na(path.prefix) || is.na(input.path.prefix) || is.na(genome.name) || is.na(sample.pattern) || is.na(independent.variable) || is.na(case.group) || is.na(control.group)) {
    if (is.na(path.prefix)) {
      cat("(\u2718) : 'path.prefix' is missing.\n\n")
    }
    if (is.na(input.path.prefix)) {
      cat("(\u2718) : 'input.path.prefix' is missing.\n\n")
    }
    if (is.na(genome.name)) {
      cat("(\u2718) : 'genome.name' is missing.\n\n")
    }
    if (is.na(sample.pattern)) {
      cat("(\u2718) : 'sample.pattern' is missing.\n\n")
    }
    if (is.na(independent.variable)) {
      cat("(\u2718) : 'independent.variable' is missing.\n\n")
    }
    if (is.na(case.group)) {
      cat("(\u2718) : 'case.group' is missing.\n\n")
    }
    if (is.na(control.group)) {
      cat("(\u2718) : 'control.group' is missing.\n\n")
    }
    stop("Input parameters ERROR")
  } else {
    cat(paste0("(\u2714) : Input parameters are all not NA !! \n\n"))
  }
}

# inner function : get operating system
CheckOperatingSystem <- function(print = TRUE){
  if (print) {
    cat(c("************** Checking operating system type ************\n"))
  }
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  os <- tolower(os)
  if (os != "linux" && os != "osx") {
    cat(paste0("(\u2718) : Your system operating system is ", os, "\n"))
    cat("       This program supports only 'linux' and 'osx'\n\n")
    stop("Operating system ERROR")
  }
  if (print) {
    cat(paste0("(\u2714) : Your system operating system is '", os, "'\n\n"))
  }
  return(os)
}

# inner fucntion : check python version
CheckPython <- function() {
  cat(c("************** Checking python version ************\n"))
  # have to check python !!!
  if(reticulate::py_available(initialize = "TRUE")){
    cat("(\u2714) : Python is available on your device!\n")
    python.version <- as.numeric(reticulate::py_config()$version)
    cat(paste0("       \u25CF Python version : ", reticulate::py_config()$version, "\n\n"))
    if(python.version >= 3) {
      return.value <- list("check.answer" = TRUE, "python.version" = 3)
    } else if (python.version < 3 && python.version >= 2 ){
      return.value <- list("check.answer" = TRUE, "python.version" = 2)
    }
    return(return.value)
  } else {
    cat("(\u2718) : Python is not available on this device. Please install python.(It's fine for python2 and python3)\n\n")
    stop("Python unavaliable ERROR")
  }
}

# inner function : check 2to3 availability
Check2to3 <- function() {
  cat(c("************** Checking 2to3 command ************\n"))
  a <- system2("2to3", "--help", stdout = FALSE)
  if (a == 0) {
    cat("(\u2714) : 2to3 command is available on your device!\n\n")
    return(TRUE)
  } else {
    cat("(\u2718) : 2to3 command is not available on your device!\n\n")
    return(FALSE)
  }
}

# inner function : check prefix.path
CheckPrefixPath <- function(path.prefix) {
  # Check the prefix exist
  if (isTRUE(dir.exists(path.prefix))){
    cat(c("************** Setting prefix path ************\n"))
    if (substr(path.prefix, nchar(path.prefix), nchar(path.prefix)) != '/') {
      path.prefix <- paste0(path.prefix, '/')
    }
    cat(paste0("(\u270e) : The following files will be installed under '", path.prefix, "'\n\n"))
    return(TRUE)
  } else if (is.na(path.prefix)) {
    cat("(\u2718) : Please give value to prefix path.\n\n")
    stop("'prefix.path' NA ERROR")
  } else {
    cat(paste0("(\u2718) : Prefix path '", path.prefix, "' is invalid. Please try another one.\n\n"))
    stop("'prefix.path' invalid ERROR")
  }
}

# inner function : check prefix.path
CheckInputPrefixPath <- function(input.path.prefix) {
  # Check the prefix exist
  if (isTRUE(dir.exists(input.path.prefix))){
    if (substr(input.path.prefix, nchar(input.path.prefix), nchar(input.path.prefix)) != '/') {
      input.path.prefix <- paste0(input.path.prefix, '/')
    }
    path.prefix.input_files <- paste0(input.path.prefix, "input_files/")
    if (isTRUE(dir.exists(path.prefix.input_files))){
      cat(c("************** Setting input prefix path ************\n"))
      cat(paste0("(\u270e) : 'input_files' is found under '", input.path.prefix, "' directory.\n"))
      cat(paste0("       \u25CF Validity of 'input_files/' will be checked.\n\n"))
      return(TRUE)
    } else {
      cat(paste0("(\u2718) : 'input_files/' is not found under '", input.path.prefix, "'. Please make sure the spelling is correct or try another input.path.prefix.\n\n"))
      stop("'input.prefix' path invalid ERROR")
    }
  } else if (is.na(input.path.prefix)) {
    cat("(\u2718) : Please give value to input prefix path.\n\n")
    stop("'input.prefix' path NA ERROR")
  } else {
    cat(paste0("(\u2718) : Input prefix path '", input.path.prefix, "' is invalid. Please try another one.\n\n"))
    stop("'input.prefix' path invalid ERROR")
  }
}

# inner function : check input.path.prefix
CheckInputDirFiles <- function(input.path.prefix, genome.name, sample.pattern) {
  cat(c("************** Checking hierarchy of", paste0("'", input.path.prefix, 'input_files/\''), "************\n"))
  # only check whether exist
  gtf.file <- file.exists(paste0(input.path.prefix, "input_files/",genome.name, ".gtf"))
  # only check whether exist
  fa.file <- file.exists(paste0(input.path.prefix, "input_files/",genome.name, ".fa"))
  # check exist and rules
  raw.fastq.dir <- dir.exists(paste0(input.path.prefix, "input_files/raw_fastq.gz/"))
  # check exist and rules
  phenodata.file <- file.exists(paste0(input.path.prefix, "input_files/phenodata.csv"))
  # check exist and rules ( this is optional)
  ht2.dir <- dir.exists(paste0(input.path.prefix, "input_files/indices/"))
  optional.indices.bool <- FALSE
  # check whether sample pattern matches the file names~
  if (isTRUE(raw.fastq.dir)) {
    raw.fastq <- list.files(path = paste0(input.path.prefix, 'input_files/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
    extract.fastq.gz.sample.names <- unique(gsub("_[1-2]*.fastq.gz", "", raw.fastq))
    check.fastq.gz.1 <- vapply(extract.fastq.gz.sample.names, function(x) paste0(x, "_1.fastq.gz"), USE.NAMES=FALSE, FUN.VALUE = "a")
    check.fastq.gz.2 <- vapply(extract.fastq.gz.sample.names, function(x) paste0(x, "_2.fastq.gz"), USE.NAMES=FALSE, FUN.VALUE = "a")
    # checking the valid file naming of '.fastq.gz'
    for ( i in seq_along(check.fastq.gz.1)) {
      bool.check.1 <- check.fastq.gz.1[i] %in% raw.fastq
      if (!bool.check.1) {
        cat(paste0("(\u2718) : ", check.fastq.gz.1[i], " is not found in 'input_files/raw_fastq.gz/\n"))
        stop("Pair-end file ERROR")
      }
    }
    for ( i in seq_along(check.fastq.gz.2)) {
      bool.check.2 <- check.fastq.gz.2[i] %in% raw.fastq
      if (!bool.check.2) {
        cat(paste0("(\u2718) : ", check.fastq.gz.2[i], " is not found in 'input_files/raw_fastq.gz/\n"))
        stop("Pair-end file ERROR")
        }
    }
    if (sample.pattern == "") {
      cat(paste0("(\u2718) : The value of 'sample.pattern' is \"\" !!\n"))
      stop("'sample.pattern' ERROR")
    }
    raw.fastq.number <- length(raw.fastq)
  }
  if (isTRUE(ht2.dir)) {
    ht2.files <- list.files(path = paste0(input.path.prefix, 'input_files/indices/'), pattern = paste0("^", genome.name, "_tran.[0-9]*.ht2$"), all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
    ht2.files.number <- length(ht2.files)
  }
  if (!isTRUE(gtf.file)) {
    grab.gtf.file <- Sys.glob(file.path(path = paste0(input.path.prefix, 'input_files'), "*.gtf"))
    cat(paste0("(\u2718) : '", genome.name, ".gtf (user input)' and '", grab.gtf.file, " (find in directory)' ", " are mismatched.\n"))
  } else {
    cat(paste0("(\u2714) : '", genome.name, ".gtf'", " is in 'input_files'\n"))
  }
  if (!isTRUE(fa.file)) {
    grab.fa.file <- Sys.glob(file.path(path = paste0(input.path.prefix, 'input_files'), "*.fa"))
    cat(paste0("(\u2718) : '", genome.name, ".fa (user input)' and '", grab.fa.file, " (find in directory)' ", " are mismatched.\n"))
  } else {
    cat(paste0("(\u2714) : '", genome.name, ".fa'", " is in 'input_files'\n"))
  }
  if (!isTRUE(raw.fastq.dir)) {
    cat(c("(\u2718) : 'raw_fastq.gz/' is missing.\n"))
  } else {
    cat(paste0("(\u2714) : 'raw_fastq.gz/' is in 'input_files'\n"))
  }
  if (isTRUE(raw.fastq.dir) && raw.fastq.number == 0) {
    cat(paste0("(\u2718) : Sample pattern \"", sample.pattern ,"\" is not found in 'raw_fastq.gz/'.\n"))
  } else if (isTRUE(raw.fastq.dir) && raw.fastq.number >= 0) {
    for (i in raw.fastq) {
      cat(paste0("(\u2714) : 'raw_fastq.gz/", i, "'"), "is in 'input_files/'\n")
    }
  }
  if (!isTRUE(phenodata.file)) {
    cat(paste0("(\u2718) : '", "phenodata.csv is missing.\n"))
  } else {
    cat(paste0("(\u2714) : '", "phenodata.csv is in 'input_files'\n"))
  }
  if (!isTRUE(ht2.dir)) {
    ## not exist
    cat(c("(\u26A0) : 'indices/' is optional. You can download the corresponding 'XXX.ht2' from 'https://ccb.jhu.edu/software/hisat2/index.shtml' to speed up the process.\n"))
  } else {
    ## exist
    cat(paste0("(\u2714) : 'indices/' is in 'input_files'\n"))
  }
  if (isTRUE(ht2.dir) && ht2.files.number == 0) {
    cat(c("(\u26A0) : 'indices/' directory has been created but there are no samples in 'indices/' or files' names", paste0("\"^", genome.name, "_tran.[0-9]*.ht2$\""), "in 'indices/' are not found.\n      No files will be copied.\n      (1). Check whether files name", paste0("'", genome.name, "_tran.[0-9]*.ht2'"), "matches the files in 'indices' directory.\n      (2). If you don't have", paste0("'", genome.name, "_tran.[0-9].ht2'"), "files, remove 'indices' directory\n\n"))
    stop("'input_files/indices/' ERROR")
  } else if (isTRUE(ht2.dir) && ht2.files.number >= 0) {
    optional.indices.bool <- TRUE
    for (i in ht2.files) {
      cat(paste0("(\u2714) : 'indices/", i, "'"), "is in 'input_files/'\n")
    }
  }
  if (isTRUE(gtf.file) && isTRUE(raw.fastq.dir) && isTRUE(raw.fastq.dir) && raw.fastq.number != 0 && isTRUE(phenodata.file)) {
    cat(c(paste0("\n(\u2714) : '", input.path.prefix,"input_files/", "'"), "is valid !\n"))
    if (isTRUE(ht2.dir) && ht2.files.number != 0) {
      cat(paste0("(\u2714) : optional directory 'indices/' is valid !\n"))
      optional.indices.bool <- TRUE
    }
    cat("\n")
    return.value <- list("check.answer" = TRUE, "optional.indices.bool" = optional.indices.bool)
    return(return.value)
  } else {
    stop("'input_files/' checking ERROR")
  }
}

# inner function : check validity of phenodata
# must have column : "ids" + "2 column"
# other column : "treatment", "tissue", "cell_type", "genotype", "time", "dosage"
CheckPhenodata <- function(input.path.prefix, genome.name, sample.pattern, independent.variable) {
  # have to sort the column !! and sort them in the correct order
  cat(c("************** Checking phenodata  ************\n"))
  pheno_data <- read.csv(paste0(input.path.prefix, "input_files/phenodata.csv"))
  # Covert all column to character
  pheno_data <- data.frame(lapply(pheno_data, as.character), stringsAsFactors=FALSE)
  # Checl 'ids' is in the 'phenodata.csv'
  if (!("ids" %in% colnames(pheno_data))) {
    cat(paste0("(\u2718) : 'ids' can't find in the column of phenodata.csv.\n\n"))
    stop("'ids' invalid ERROR")
  }
  pheno_data.ids.list <- pheno_data["ids"]
  # "id" : must be distinct, same as input_files raw reads name !
  cat("\u25B6 Checking whether \"raw_fastq.gz files\" matches \"'ids' of phenodata.csv\" \n")
  raw.fastq <- list.files(path = paste0(input.path.prefix, 'input_files/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  extract.fastq.gz.sample.names <- unique(gsub("_[1-2]*.fastq.gz", "", raw.fastq))
  bool.length <- length(extract.fastq.gz.sample.names) == length(pheno_data$ids)
  bool.identical <- identical(sort(extract.fastq.gz.sample.names), sort(pheno_data$ids))
  if (!bool.length || !bool.identical) {
    cat(paste0("(\u2718) : 'ids' column doesn't match the smaple_id in 'input_files/raw_fastq.gz'. Please check the file.\n" ))
    stop("'ids' mismatch ERROR")
  }
  ids.list <- paste(sort(extract.fastq.gz.sample.names), collapse = ", ")
  cat("(\u2714) : Column 'ids' of phenodata.csv is valid. \n")
  cat(paste0("      \u25CF sample ids are : \"", ids.list,"\"\n"))
  # Check again independent.variable in the list
  if (!(independent.variable %in% colnames(pheno_data))) {
    cat(paste0("(\u2718) : 'independent.variable' : '", independent.variable, "' can't find in the column of phenodata.csv.\n\n"))
    stop("'independent.variable' invalid ERROR")
  }
  cat(paste0("(\u2714) : 'independent.variable' : '", independent.variable, "' is in the column of phenodata.csv. \n\n"))
  cat(paste0("\u25B6 Checking whether '", independent.variable, "' is a two-group 'independent.variable' ...\n"))
  length.independent.variable <- length(table(pheno_data[independent.variable]))
  if (!(length.independent.variable == 2)) {
    cat(paste0("(\u2718) : 'independent.variable' : '", independent.variable, "' is a ", length.independent.variable,"-group 'independent.variable'. Not 2-group.\n"))
    cat(paste0("          groups that found : ", paste0(names(table(pheno_data[independent.variable])), collapse = ", "), "\n"))
    stop("'independent.variable' none-two-group ERROR")
  }
  cat("(\u2714) : Column 'independent.variable' : '",independent.variable, "' of phenodata.csv is valid. \n")
  cat(paste0("      \u25CF 'independent.variable' : '", independent.variable, "'\n\n"))
  return(TRUE)
}

# inner function
CheckCaseControlGroup <- function(input.path.prefix, independent.variable, case.group, control.group) {
  cat(c("************** Checking 'case.group' & 'control.group' ************\n"))
  pheno_data <- read.csv(paste0(input.path.prefix, "/input_files/phenodata.csv"))
  # Covert all column to character
  pheno_data <- data.frame(lapply(pheno_data, as.character), stringsAsFactors=FALSE)
  cont.in <- case.group %in% as.character(data.frame(table(pheno_data[independent.variable]))$Var1)
  exp.in <- control.group %in% as.character(data.frame(table(pheno_data[independent.variable]))$Var1)
  # Check 'case.group' is on group of 'independent.variable'
  if (!cont.in) {
    cat(paste0("(\u2718) : 'case.group' : '", case.group, "' is not a group of in 'independent.variable'.\n\n"))
    stop("'case.group' invalid ERROR")
  }
  if (!exp.in) {
    cat(paste0("(\u2718) : 'control.group' : '", control.group, "' is not a group of in 'independent.variable'.\n\n"))
    stop("'control.group' invalid ERROR")
  }
  if (exp.in && cont.in) {
    cat(paste0("(\u2714) :    'case.group' : '", case.group, "' is a group of in 'independent.variable'.\n"))
    cat(paste0("(\u2714) : 'control.group' : '", control.group, "' is a group of in 'independent.variable'.\n\n"))
    if (control.group == case.group) {
      cat(paste0("(\u2718) : 'case.group' and 'control.group' are same.\n\n"))
      stop("'case.group' &'control.group' same ERROR")
    }
  return(TRUE)
  }
}

#
CheckS4Object <- function(RNASeqWorkFlowParam, print = TRUE) {
  if (isS4(RNASeqWorkFlowParam) && class(RNASeqWorkFlowParam)[1] == "RNASeqWorkFlowParam") {
    if (print) {
      cat(c("************** Checking validity of S4 input ************\n"))
      cat(paste0("     (\u2714) : input is valid 'RNASeqWorkFlowParam' instance! \n\n"))
    }
  } else {
    cat(paste0("(\u2718) : input is not a valid 'RNASeqWorkFlowParam' instance!.\n" ))
    stop("Invalid 'RNASeqWorkFlowParam' input ERROR")
  }
}
