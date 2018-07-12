#' RNASeq
#' An S4 class for storing RNA-Seq workflow parameters of this package
#' @aliases RNASeq
#'
#' @slot os.type 'linux' or 'osx'. The operating system type
#' @slot python.variable a list with (whether python is valid, python version)
#' @slot path.prefix the directory holding installations and analysis results
#' @slot input.path.prefix user has to prepared valid 'input_files/' under this directory
#' @slot gene.name gene name defined in this RNA-Seq workflow (ex. gene.name.fa, gene.name.gtf)
#' @slot sample.pattern  sample pattern describing the name of raw fastq.gz files
#' @slot experiment.type set the type of the RNA-Seq analysis workflow. Character of one of "two.group", "multi.group.pairs", "multi.group.anova"
#' @slot main.variable main sample grouping variable
#' @slot additional.variable additional sample information
#' @slot indexes.optional logical value whether indexes/ is exit in 'input_files/'
#'
#' @name RNASeqWorkFlowParam-class
#'
#' @rdname RNASeqWorkFlowParam-class
#'
#' @exportClass RNASeqWorkFlowParam
#' @author Kuan-Hao, Chao
#' @examples
#' data(workflowParam)
#' class(workflowParam) #"RNASeqWorkFlowParam"
#' workflowParam@@path.prefix
#' workflowParam@@input.path.prefix
#' workflowParam@@gene.name
#' workflowParam@@sample.pattern
#' workflowParam@@experiment.type
#' workflowParam@@main.variable
#' workflowParam@@additional.variable
setClass("RNASeqWorkFlowParam",
         representation(
           os.type = "character",
           python.variable = "list",
           path.prefix = "character",
           input.path.prefix = "character",
           gene.name = "character",
           sample.pattern = "character",
           experiment.type = "character",
           main.variable = "character",
           additional.variable = "character",
           indexes.optional = "logical"
         )
)


#' Constructor function for RNASeqWorkFlowParam objects
#'
#' @name RNASeqWorkFlowParam-constructor
#'
#' @param  path.prefix the directory holding installations and analysis results
#' @param input.path.prefix user has to prepared valid 'input_files/' under this directory
#' @param gene.name gene name defined in this RNA-Seq workflow (ex. gene.name.fa, gene.name.gtf)
#' @param sample.pattern  sample pattern describing the name of raw fastq.gz files
#' @param experiment.type set the type of the RNA-Seq analysis workflow. Character of one of "two.group", "pair-wise.group", "multi.group"
#' @param main.variable main sample grouping variable
#' @param additional.variable additional sample information
#'
#' @return an object of class \code{RNASeqWorkFlowParam}
#'
#' @author kuan-hao Chao
#'
#' @rdname RNASeqWorkFlowParam-constructor
#'
#' @importFrom tools file_ext
#' @export
#' @example
#' exp <- RNASeqWorkFlowParam(path.prefix = "/home/rnaseq", input.path.prefix = "/home", gene.name = "hg19", sample.pattern = "SRR[0-9]",
#'                            experiment.type = "two.group", main.variable = "treatment", additional.variable = "cell")
RNASeqWorkFlowParam <- function(path.prefix = NA, input.path.prefix = NA, gene.name = NA, sample.pattern = NA,
                                experiment.type = NA, main.variable = NA, additional.variable = NA) {
  # check input parameters
  CheckInputParam(path.prefix, input.path.prefix, gene.name, sample.pattern,
                  experiment.type, main.variable, additional.variable)
  # 1. check operating system
  characters.os.type <- CheckOperatingSystem(print = TRUE)
  # 2. check python version
  python.version.list <- CheckPython(print = TRUE)
  bool.python.avail <- python.version.list$check.answer
  numeric.python.version <- python.version.list$python.version
  # 3. check validity of path.prefix
  bool.prefix.path <- CheckPrefixPath(path.prefix = path.prefix, print = TRUE)
  if (bool.prefix.path){
    # add '/' to the path.prefix
    if (substr(path.prefix, nchar(path.prefix), nchar(path.prefix)) != '/') {
      path.prefix <- paste0(path.prefix, '/')
    }
  }
  # 4. check input.path.prefix
  bool.input.path.prefix <- CheckInputPrefixPath(input.path.prefix = input.path.prefix, print = TRUE)
  if (bool.input.path.prefix){
    # add '/' to the path.prefix
    if (substr(input.path.prefix, nchar(input.path.prefix), nchar(input.path.prefix)) != '/') {
      input.path.prefix <- paste0(input.path.prefix, '/')
    }
  }
  # 5. check sample.pattern that can't have '.fastq.gz'
  fast.gz.extend <- tools::file_ext(sample.pattern)
  if (fast.gz.extend == "gz" || fast.gz.extend == "fastq") {
    cat("(\u2718) : 'sample.pattern' can't include file extension(.gz or .fastq or .fastq.gz)\n\n")
    stop("sample.pattern with extension error.")
  }
  # 6. check 'input_files/' necessary files with 'gene.name', 'sample.pattern'
  input.dir.files.list <- CheckInputDirFiles(input.path.prefix = input.path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print = TRUE)
  bool.input.dir.files <- input.dir.files.list$check.answer
  bool.input.dir.indexes <- input.dir.files.list$optional.indexes.bool
  # 7. check 'experiment.type'
  bool.experiment.type <- CheckExperimentType(experiment.type = experiment.type)

  # below still need to fix
  # 8. check 'phenodata'
  phenodata.return <- CheckPhenodata(input.path.prefix = input.path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print=TRUE)
  bool.phenodata <- phenodata.return$check.answer
  # 9. check 'main variable'
  bool.check.main.var <- CheckMainVar(input.path.prefix = input.path.prefix, main.variable = main.variable, experiment.type = experiment.type,
                                      treatment.num = phenodata.return$treatment.num, tissue.num = phenodata.return$tissue.num, cell_type.num = phenodata.return$cell_type.num,
                                      genotype.num = phenodata.return$genotype.num, time.num = phenodata.return$time.num, dosage.time = phenodata.return$dosage.time,
                                      print=TRUE)
  # 10. check 'additional.variable'
  bool.check.add.var <- CheckAddVar(additional.variable = additional.variable, main.variable = main.variable)

  if ((characters.os.type == "linux" || characters.os.type == "osx") && bool.python.avail && bool.prefix.path &&
      bool.input.path.prefix && bool.input.dir.files && bool.experiment.type && bool.phenodata && bool.check.main.var && bool.check.add.var) {
    cat(paste0("\n**************************************\n"))
    cat(paste0("************** Success! **************\n"))
    cat(paste0("**************************************\n"))
    new("RNASeqWorkFlowParam",os.type = characters.os.type, python.variable = python.version.list, path.prefix = path.prefix,
        input.path.prefix = input.path.prefix, gene.name = gene.name, sample.pattern = sample.pattern,
        experiment.type = experiment.type, main.variable = main.variable, additional.variable = additional.variable,
        indexes.optional = bool.input.dir.indexes)
  }
}

#' inner function : check whether input values are NA
CheckInputParam <- function(path.prefix = NA, input.path.prefix = NA, gene.name = NA, sample.pattern = NA,
                            experiment.type = NA, main.variable = NA, additional.variable = NA) {
  cat(c("************** Checking input parameters ************\n"))
  if (is.na(path.prefix) || is.na(input.path.prefix) || is.na(gene.name) || is.na(sample.pattern) || is.na(experiment.type) || is.na(main.variable) || is.na(additional.variable)) {
    if (is.na(path.prefix)) {
      cat("(\u2718) : 'path.prefix' is missing.\n\n")
    }
    if (is.na(input.path.prefix)) {
      cat("(\u2718) : 'input.path.prefix' is missing.\n\n")
    }
    if (is.na(gene.name)) {
      cat("(\u2718) : 'gene.name' is missing.\n\n")
    }
    if (is.na(sample.pattern)) {
      cat("(\u2718) : 'sample.pattern' is missing.\n\n")
    }
    if (is.na(experiment.type)) {
      cat("(\u2718) : 'experiment.type' is missing.\n\n")
    }
    if (is.na(main.variable)) {
      cat("(\u2718) : 'main.variable' is missing.\n\n")
    }
    if (is.na(additional.variable)) {
      cat("(\u2718) : 'additional.variable' is missing.\n\n")
    }
    stop("Input parameters ERROR")
  } else {
    cat(paste0("(\u2714) : Input parameters are all valid valid !! \n\n"))
  }
}

#' inner function : check prefix.path
CheckPrefixPath <- function(path.prefix = NA_character_, print = TRUE) {
  # Check the prefix exist
  if (isTRUE(dir.exists(path.prefix))){
    if (print) {
      cat(c("************** Setting prefix path ************\n"))
      if (substr(path.prefix, nchar(path.prefix), nchar(path.prefix)) != '/') {
        path.prefix <- paste0(path.prefix, '/')
      }
      cat(paste0("(\u270e) : The following files will be installed under '", path.prefix, "'\n\n"))
    }
    return(TRUE)
  } else if (is.na(path.prefix)) {
    cat("(\u2718) : Please give value to prefix path.\n\n")
    stop("Prefix path NA ERROR")
    return(FALSE)
  } else {
    cat(paste0("(\u2718) : Prefix path '", path.prefix, "' is invalid. Please try another one.\n\n"))
    stop("Prefix path invalid ERROR")
    return(FALSE)
  }
}

#' inner function : check prefix.path
CheckInputPrefixPath <- function(input.path.prefix = NA_character_, print = TRUE) {
  # Check the prefix exist
  if (isTRUE(dir.exists(input.path.prefix))){
    if (substr(input.path.prefix, nchar(input.path.prefix), nchar(input.path.prefix)) != '/') {
      input.path.prefix <- paste0(input.path.prefix, '/')
    }
    path.prefix.input_files <- paste0(input.path.prefix, "input_files/")
    if (isTRUE(dir.exists(path.prefix.input_files))){
      if (print) {
        cat(c("************** Setting input prefix path ************\n"))
        cat(paste0("(\u270e) : 'input_files' is found under '", input.path.prefix, "' directory.\n"))
        cat(paste0("       Validity of 'input_files/' will be checked.\n\n"))
      }
      return(TRUE)
    } else {
      cat(paste0("(\u2718) : 'input_files/' is not found under '", input.path.prefix, "'. Please make sure the spelling is correct or try another input.path.prefix.\n\n"))
      stop("Input prefix path invalid ERROR")
      return(FALSE)
    }
  } else if (is.na(input.path.prefix)) {
    cat("(\u2718) : Please give value to input prefix path.\n\n")
    stop("Input prefix path NA ERROR")
    return(FALSE)
  } else {
    cat(paste0("(\u2718) : Input prefix path '", input.path.prefix, "' is invalid. Please try another one.\n\n"))
    stop("Input prefix path invalid ERROR")
    return(FALSE)
  }
}


#' inner function : check input.path.prefix
CheckInputDirFiles <- function(input.path.prefix = NA_character_, gene.name = NA_character_, sample.pattern = NA_character_, print=TRUE) {
  if (print) {
    cat(c("************** Checking hierarchy of", paste0("'", input.path.prefix, 'input_files/\''), "************\n"))
  }
  # only check whether exist
  gtf.file <- file.exists(paste0(input.path.prefix, "input_files/",gene.name, ".gtf"))
  # only check whether exist
  fa.file <- file.exists(paste0(input.path.prefix, "input_files/",gene.name, ".fa"))
  # check exist and rules
  raw.fastq.dir <- dir.exists(paste0(input.path.prefix, "input_files/raw_fastq.gz/"))
  # check exist and rules
  phenodata.file <- file.exists(paste0(input.path.prefix, "input_files/phenodata.csv"))
  # check exist and rules ( this is optional)
  ht2.dir <- dir.exists(paste0(input.path.prefix, "input_files/indexes/"))
  optional.indexes.bool <- FALSE
  # check whether sample pattern matches the file names~
  if (isTRUE(raw.fastq.dir)) {
    raw.fastq <- list.files(path = paste0(input.path.prefix, 'input_files/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
    extract.fastq.gz.sample.names <- unique(gsub("_[1-2]*.fastq.gz", "", raw.fastq))
    check.fastq.gz.1 <- sapply(extract.fastq.gz.sample.names, function(x) paste0(x, "_1.fastq.gz"),USE.NAMES=FALSE)
    check.fastq.gz.2 <- sapply(extract.fastq.gz.sample.names, function(x) paste0(x, "_2.fastq.gz"),USE.NAMES=FALSE)
    # checking the valid file naming of '.fastq.gz'
    for ( i in 1:length(check.fastq.gz.1)) {
      bool.check.1 <- check.fastq.gz.1[i] %in% raw.fastq
      if (!bool.check.1) {
        cat(paste0("(\u2718) : ", check.fastq.gz.1[i], " is not found in 'input_files/raw_fastq.gz/\n"))
        stop("Pair-end file ERROR")
      }
    }
    for ( i in 1:length(check.fastq.gz.2)) {
      bool.check.2 <- check.fastq.gz.2[i] %in% raw.fastq
      if (!bool.check.2) {
        cat(paste0("(\u2718) : ", check.fastq.gz.2[i], " is not found in 'input_files/raw_fastq.gz/\n"))
        stop("Pair-end file ERROR")
        }
    }
    raw.fastq.number <- length(raw.fastq)
  }
  if (isTRUE(ht2.dir)) {
    ht2.files <- list.files(path = paste0(input.path.prefix, 'input_files/indexes/'), pattern = paste0("^", gene.name, "_tran.[0-9]*.ht2$"), all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
    ht2.files.number <- length(ht2.files)
  }
  if (!isTRUE(gtf.file)) {
    grab.gtf.file <- Sys.glob(file.path(path = paste0(input.path.prefix, 'input_files'), "*.gtf"))
    cat(paste0("(\u2718) : '", gene.name, ".gtf (user input)' and '", grab.gtf.file, " (find in directory)' ", " are mismatched.\n"))
  } else {
    if (print) {
      cat(paste0("(\u2714) : '", gene.name, ".gtf'", " is in 'input_files'\n"))
    }
  }
  if (!isTRUE(fa.file)) {
    grab.fa.file <- Sys.glob(file.path(path = paste0(input.path.prefix, 'input_files'), "*.fa"))
    cat(paste0("(\u2718) : '", gene.name, ".fa (user input)' and '", grab.fa.file, " (find in directory)' ", " are mismatched.\n"))
  } else {
    if (print) {
      cat(paste0("(\u2714) : '", gene.name, ".fa'", " is in 'input_files'\n"))
    }
  }
  if (!isTRUE(raw.fastq.dir)) {
    cat(c("(\u2718) : 'raw_fastq.gz/' is missing.\n"))
  } else {
    if (print) {
      cat(paste0("(\u2714) : 'raw_fastq.gz/' is in 'input_files'\n"))
    }
  }
  if (isTRUE(raw.fastq.dir) && raw.fastq.number == 0) {
    cat(paste0("(\u2718) : Sample pattern \"", sample.pattern ,"\" is not found in 'raw_fastq.gz/'.\n"))
  } else if (isTRUE(raw.fastq.dir) && raw.fastq.number >= 0) {
    if (print){
      for (i in raw.fastq) {
        cat(paste0("(\u2714) : 'raw_fastq.gz/", i, "'"), "is in 'input_files/'\n")
      }
    }
  }
  if (!isTRUE(phenodata.file)) {
    cat(paste0("(\u2718) : '", "phenodata.csv is missing.\n"))
  } else {
    if (print) {
      cat(paste0("(\u2714) : '", "phenodata.csv is in 'input_files'\n"))
    }
  }
  if (!isTRUE(ht2.dir)) {
    ## not exist
    if (print) {
      cat(c("(\u26A0) : 'indexes/' is optional. You can download the corresponding 'XXX.ht2' from 'https://ccb.jhu.edu/software/hisat2/index.shtml' to speed up the process.\n"))
    }
  } else {
    if (print) {
      ## exist
      cat(paste0("(\u2714) : 'indexes/' is in 'input_files'\n"))
    }
  }
  if (isTRUE(ht2.dir) && ht2.files.number == 0) {
    cat(c("(\u26A0) : 'indexes/' directory has been created but there are no samples in 'indexes/' or files' names", paste0("\"^", gene.name, "_tran.[0-9]*.ht2$\""), "in 'indexes/' are not found.\n      No files will be copied.\n      (1). Check whether files name", paste0("'", gene.name, "_tran.[0-9]*.ht2'"), "matches the files in 'indexes' directory.\n      (2). If you don't have", paste0("'", gene.name, "_tran.[0-9].ht2'"), "files, remove 'indexes' directory\n\n"))
    stop("'input_files/indexes/' ERROR")
  } else if (isTRUE(ht2.dir) && ht2.files.number >= 0) {
    optional.indexes.bool <- TRUE
    if (print){
      for (i in ht2.files) {
        cat(paste0("(\u2714) : 'indexes/", i, "'"), "is in 'input_files/'\n")
      }
    }
  }
  if (isTRUE(gtf.file) && isTRUE(raw.fastq.dir) && isTRUE(raw.fastq.dir) && raw.fastq.number != 0 && isTRUE(phenodata.file)) {
    if (print) {
      cat(c(paste0("\n(\u2714) : '", input.path.prefix,"input_files/", "'"), "is valid !\n"))
      if (isTRUE(ht2.dir) && ht2.files.number != 0) {
        cat(paste0("(\u2714) : optional directory 'indexes/' is valid !\n\n"))
        optional.indexes.bool <- TRUE
      }
    }
    return.value <- list("check.answer" = TRUE, "optional.indexes.bool" = optional.indexes.bool)
    return(return.value)
  } else {
    stop("'input_files/' checking ERROR")
  }
}

#' inner function : check experiment type
CheckExperimentType <- function(experiment.type = NA_character_, print = TRUE) {
  if (print) {
    cat(c("************** Checking experiment.type ************\n"))
  }
  if ((experiment.type != "two.group") &&  (experiment.type != "pair-wise.group") && (experiment.type != "multi.group")) {
    cat("(\u2718) : 'experiment.type' can only be 'two.group' or 'pair-wise.group' or 'multi.group'.\n\n")
    stop("Experiment.type ERROR")
  } else {
    cat(paste0("(\u2714) : 'experiment.type' is \"", experiment.type, "\"\n\n"))
    return(TRUE)
  }
}

#' inner fucntion : check python version
#' @importFrom reticulate py_config py_available
CheckPython <- function(print = TRUE) {
  if (print) {
    cat(c("************** Checking python version ************\n"))
  }
  # have to check python !!!
  if(reticulate::py_available(initialize = "TRUE")){
    cat("(\u2714) : Python is available on your device!\n")
    python.version <- as.numeric(reticulate::py_config()$version)
    cat(paste0("       Python version : ", reticulate::py_config()$version, "\n\n"))
    if(python.version >= 3) {
      return.value <- list("check.answer" = TRUE, "python.version" = 3)
    } else if (python.version < 3 && python.version >= 2 ){
      return.value <- list("check.answer" = TRUE, "python.version" = 2)
    }
  } else {
    cat("(\u2718) : Python is not available on this device. Please install python.(It's fine for python2 and python3)\n\n")
    stop("Python unavaliable ERROR")
  }
}

#' inner function : get operating system
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
    cat(paste0("(\u2718) : Your system operating system is", os, "\n"))
    cat("       This program supports only linux and osx.\n\n")
    stop("Operating system ERROR")
  }
  cat(paste0("(\u2714) : Your system operating system is '", os, "'\n\n"))
  return(os)
}

#' inner function : check validity of phenodata
#' must have column : "ids" + "2 column"
#' other column : "treatment", "tissue", "cell_type", "genotype", "time", "dosage"
CheckPhenodata <- function(input.path.prefix = NA_character_, gene.name = NA_character_, sample.pattern = NA_character_, print=TRUE) {
  # have to sort the column !! and sort them in the correct order

  cat(c("************** Checking phenodata  ************\n"))
  bool.check.valid <- TRUE
  pheno_data <- read.csv(paste0(input.path.prefix, "/input_files/phenodata.csv"))
  # Covert all column to character
  pheno_data <- data.frame(lapply(pheno_data, as.character), stringsAsFactors=FALSE)
  # Columns that read from the file
  pheno_data_cols<- colnames(pheno_data)
  # Make all the column
  columns.any.NA <- apply(pheno_data, 2, function(x) any(is.na(x)))
  # Files must have these columns in order
  must_have_column <- c("ids", "treatment", "tissue", "cell_type", "genotype", "time", "dosage")
  # check column names are matched
  cat("     \u25CF Checking column names of 'phenodata.csv'\n")
  for ( i in 1:7 ) {
    if (pheno_data_cols[i] != must_have_column[i]) {
      bool.check.valid <- FALSE
      cat(paste0("(\u2718) : Your ", i, " column name is '",pheno_data_cols[i], "' not correpond to the expected '", must_have_column[i],"'.\n" ))
      cat("       Expected column names : \" 'ids', 'treatment', 'tissue', 'cell_type', 'genotype', 'time', 'dosage' \".\n\n")
      stop("Column names ERROR")
    }
  }
  cat("         (\u2714) : column names are valid. ('ids', 'treatment', 'tissue', 'cell_type', 'genotype', 'time', 'dosage')\n")

  # "id" : must be distinct, same as input_files raw reads name !
  cat("     \u25CF Checking whether \"raw_fastq.gz files\" matches \"'ids' of phenodata.csv\" \n")
  raw.fastq <- list.files(path = paste0(input.path.prefix, 'input_files/raw_fastq.gz/'), pattern = sample.pattern, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  extract.fastq.gz.sample.names <- unique(gsub("_[1-2]*.fastq.gz", "", raw.fastq))
  bool.length <- length(extract.fastq.gz.sample.names) == length(pheno_data$ids)
  bool.identical <- identical(sort(extract.fastq.gz.sample.names), sort(pheno_data$ids))
  if (!bool.length || !bool.identical) {
    cat(paste0("(\u2718) : 'ids' column doesn't match the smaple_id in 'input_files/raw_fastq.gz'. Please check the file.\n" ))
    stop("Ids column ERROR")
  }
  ids.list <- paste(sort(extract.fastq.gz.sample.names), collapse = " ")
  cat("         (\u2714) : Column 'ids' of phenodata.csv is valid. \n")
  cat(paste0("            sample ids are : \"", ids.list,"\"\n"))
  # first and second (ids and treatment) must can't have NA value : it should be false!
  # check each column : if there is any NA
  columns.any.NA <- apply(pheno_data, 2, function(x) any(is.na(x)))
  # check each column : if all is not NA
  columns.all.not.NA <- apply(pheno_data, 2, function(x) all(!is.na(x)))
  # check each column : if there is all NA
  columns.all.NA <- apply(pheno_data, 2, function(x) all(is.na(x)))
  cat("     \u25CF Checking 'ids' and 'treatment' of 'phenodata.csv' (can't have NA)\n")
  if (columns.any.NA[[1]] || columns.any.NA[[2]]) {
    if (columns.any.NA[[1]]) {
      cat(paste0("         (\u2718) : There are NA values in 'ids' column. Please check the files.\n" ))
    }
    if (columns.any.NA[[2]]) {
      cat(paste0("         (\u2718) : There are NA values in or in 'treatment' column. Please check the files.\n" ))
    }
    stop("Necessary column NA ERROR")
  }
  # 'tissue', 'cell_type', 'genotype', 'time', 'dosage' : if there is any one value that is not NA in each column, then that columns must can't have NA vlue
  # condition : must 'all NA' or 'all have value'
  cat("     \u25CF Checking 'tissue', 'cell_type', 'genotype', 'time', 'dosage' of 'phenodata.csv'\n")
  if ( (!columns.all.NA[[3]] && !columns.all.not.NA[[3]]) || (!columns.all.NA[[4]] && !columns.all.not.NA[[4]]) ||
       (!columns.all.NA[[5]] && !columns.all.not.NA[[5]]) || (!columns.all.NA[[6]] && !columns.all.not.NA[[6]]) ||
       (!columns.all.NA[[7]] && !columns.all.not.NA[[7]]) ) {
    if (!columns.all.NA[[3]] && !columns.all.not.NA[[3]]) {
      cat(paste0("         (\u2718) : Invalid column value in 'tissue'. Please check the files.\n" ))
    }
    if (!columns.all.NA[[4]] && !columns.all.not.NA[[4]]) {
      cat(paste0("         (\u2718) : Invalid column value in 'cell_type'. Please check the files.\n" ))
    }
    if (!columns.all.NA[[5]] && !columns.all.not.NA[[5]]) {
      cat(paste0("         (\u2718) : Invalid column value in 'genotype'. Please check the files.\n" ))
    }
    if (!columns.all.NA[[6]] && !columns.all.not.NA[[6]]) {
      cat(paste0("         (\u2718) : Invalid column value in 'time'. Please check the files.\n" ))
    }
    if (!columns.all.NA[[7]] && !columns.all.not.NA[[7]]) {
      cat(paste0("         (\u2718) : Invalid column value in 'dosage'. Please check the files.\n" ))
    }
    stop("Optional column NA ERROR")
  }
  # calculate groups number of all samples
  cat("     \u25CF Calculating groups numbers of each column in 'phenodata.csv'\n")
  ids.table <- table(pheno_data$ids)
  treatment.table <- table(pheno_data$treatment)
  tissue.table <- table(pheno_data$tissue)
  cell_type.table <- table(pheno_data$cell_type)
  genotype.table <- table(pheno_data$genotype)
  time.table <- table(pheno_data$time)
  dosage.table <- table(pheno_data$dosage)
  cat(paste0("         \u25CF There are ", length(ids.table), " sample ids\n"))
  cat(paste0("         \u25CF Groups\n"))
  cat(paste0("            \u25CF treatment : ", length(treatment.table), "\n"))
  cat(paste0("            \u25CF tissue : ", length(tissue.table), "\n"))
  cat(paste0("            \u25CF cell_type : ", length(cell_type.table), "\n"))
  cat(paste0("            \u25CF genotype : ", length(genotype.table), "\n"))
  cat(paste0("            \u25CF time : ", length(time.table), "\n"))
  cat(paste0("            \u25CF dosage : ", length(dosage.table), "\n\n"))
  return.value <- list("check.answer" = bool.check.valid, "ids.num" = length(ids.table),
                       "treatment.num" = length(treatment.table), "tissue.num" = length(tissue.table),
                       "cell_type.num" = length(cell_type.table), "genotype.num" = length(genotype.table),
                       "time.num" = length(time.table), "dosage.time" = length(dosage.table))
  return(return.value)
}

#' inner function : check main variable
CheckMainVar <- function(input.path.prefix = NA_character_, main.variable = NA_character_, experiment.type = NA_character_,
                         treatment.num = NA_character_, tissue.num = NA_character_, cell_type.num = NA_character_,
                         genotype.num = NA_character_, time.num = NA_character_, dosage.time = NA_character_, print=TRUE) {
  cat(c("************** Checking main.variable ************\n"))
  if (main.variable == "ids") {
    cat(paste0("(\u2718) : 'main.variable' can't be 'ids'.\n" ))
    stop("Main variable ERROR")
  } else if (main.variable == "treatment" || main.variable == "tissue" || main.variable == "cell_type" ||
             main.variable == "genotype" || main.variable == "time" || main.variable == "dosage") {
    # now pheno_data is valid
    pheno_data <- read.csv(paste0(input.path.prefix, "/input_files/phenodata.csv"))
    main.variable.group.num <- length(table(pheno_data[main.variable]))
    cat(paste0("     \u25CF        input 'main.variable' : \"", main.variable, "\"\n"))
    cat(paste0("     \u25CF 'main.variable' group number : ", main.variable.group.num, "\n"))
    cat(paste0("     \u25CF      input 'experiment.type' : \"", experiment.type, "\"\n"))
    if (experiment.type == "two.group") {
      if (main.variable.group.num != 2) {
        cat(paste0("(\u2718) : 'main.variable' group number must be 2 !! Not matching experiment.type. Experiment input invalid.\n" ))
        stop("experiment.type & main.variable.group.num not matching ERROR")
      }
    } else if (experiment.type == "multi.group.pairs") {
      if (main.variable.group.num <= 2) {
        cat(paste0("(\u2718) : 'main.variable' group number must be 3 or more !! Not matching experiment.type. Experiment input invalid.\n" ))
        stop("experiment.type & main.variable.group.num not matching ERROR")
      }
    } else if (experiment.type == "multi.group.anova") {
      if (main.variable.group.num <= 2) {
        cat(paste0("(\u2718) : 'main.variable' group number must be 3 or more !! Not matching experiment.type. Experiment input invalid.\n" ))
        stop("experiment.type & main.variable.group.num not matching ERROR")
      }
    }
    cat(paste0("     (\u2714) : valid 'main variable'\n\n"))
    return(TRUE)
  } else {
    cat(paste0("(\u2718) : 'main.variable' is not matching any column name of 'phenodata.csv'. Please check your input. \n" ))
    cat("      'main.variable' must be 'treatment', 'tissue', 'cell_type', 'genotype', 'time' or 'dosage'. \n")
    stop("Main variable ERROR")
  }
}

#" inner function : check additional variable
CheckAddVar <- function(additional.variable = NA_character_, main.variable = NA_character_) {
  cat(c("************** Checking additional.variable ************\n"))
  if (additional.variable == "ids") {
    cat(paste0("(\u2718) : 'additional.variable' can't be 'ids'.\n" ))
    stop("Main variable ERROR")
  } else if (additional.variable == "treatment" || additional.variable == "tissue" || additional.variable == "cell_type" ||
             additional.variable == "grenotype" || additional.variable == "time" || additional.variable == "dosage") {
    if (additional.variable != main.variable) {
      # check additional variable column : if all is not NA
      if (all(!is.na(pheno_data[additional.variable]))) {
        cat(paste0("     \u25CF  input 'additional.variable' : \"", additional.variable, "\"\n"))
        cat(paste0("     (\u2714) : valid 'additional variable'\n\n"))
        return(TRUE)
      }
      cat(paste0("(\u2718) : There are 'NA' in 'additional.variable' column.\n" ))
      stop("Additional variable NA ERROR")
    } else {
      cat(paste0("(\u2718) : 'additional.variable' can't be 'main.variable'.\n" ))
      stop("Main variable Additional variable same ERROR")
    }
  } else {
    cat(paste0("(\u2718) : 'additional.variable' is not matching any column name of 'phenodata.csv'. Please check your input. \n" ))
    cat("      'additional.variable' must be 'treatment', 'tissue', 'cell_type', 'genotype', 'time' or 'dosage'. \n")
    stop("Additional variable ERROR")
  }
}

CheckS4Object <- function(RNASeqWorkFlowParam, print = TRUE) {
  if (isS4(RNASeqWorkFlowParam) && class(RNASeqWorkFlowParam)[1] == "RNASeqWorkFlowParam") {
    cat(c("************** Checking validity of S4 input ************\n"))
    if (print) {
      cat(paste0("     (\u2714) : input is valid 'RNASeqWorkFlowParam' instance! \n\n"))
    }
  } else {
    cat(paste0("(\u2718) : input is not a valid 'RNASeqWorkFlowParam' instance!.\n" ))
    stop("Invalid 'RNASeqWorkFlowParam' input ERROR")
  }
}
