#' @title RNASeqRParam_Bam
#'
#' @description  An S4 class for checking and storing RNA-Seq workflow
#'  parameters starting with BAM files.
#'
#'
#' @slot os.type 'linux' or 'osx'. The operating system type.
#' @slot python.variable A list storing python environment.
#'    \code{(check.answer, python.version)}
#' @slot python.2to3 Logical value whether \code{2to3} command is available
#'    on the workstation.
#' @slot path.prefix Path prefix of 'gene_data/', 'RNASeq_bin/',
#'    'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories.
#' @slot input.path.prefix Path prefix of 'input_files/' directory,
#' @slot genome.name Variable of genome name defined in this RNA-Seq workflow
#'   (ex. \code{genome.name}.fa, \code{genome.name}.gtf).
#' @slot sample.pattern  Regular expression of paired-end fastq.gz files under
#'   'input_files/raw_bam'. Expression not includes \code{_[1,2].fastq.gz}.
#' @slot independent.variable Independent variable for the biological.
#' experiment design of two-group RNA-Seq workflow.
#' @slot case.group Group name of the case group.
#' @slot control.group Group name of the control group.
#'
#' @name RNASeqRParam_Bam-class
#'
#' @rdname RNASeqRParam_Bam-class
#'
#' @exportClass RNASeqRParam_Bam
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' "@@"(yeast, os.type)
#' "@@"(yeast, python.variable)
#' "@@"(yeast, python.2to3)
#' "@@"(yeast, path.prefix)
#' "@@"(yeast, input.path.prefix)
#' "@@"(yeast, genome.name)
#' "@@"(yeast, sample.pattern)
#' "@@"(yeast, independent.variable)
#' "@@"(yeast, case.group)
#' "@@"(yeast, control.group)
setClass("RNASeqRParam_Bam",
         representation(
           os.type              = "character",
           python.variable      = "list",
           python.2to3          = "logical",
           path.prefix          = "character",
           input.path.prefix    = "character",
           genome.name          = "character",
           sample.pattern       = "character",
           independent.variable = "character",
           case.group           = "character",
           control.group        = "character"
         )
)



#' @title RNASeqR_Bam
#' @description  Constructor function for RNASeqRParam_Bam objects
#'
#' @name RNASeqRParam_Bam-constructor
#'
#' @param path.prefix Path prefix of 'gene_data/', 'RNASeq_bin/',
#'   'RNASeq_results/', 'Rscript/' and 'Rscript_out/' directories.
#' @param input.path.prefix Path prefix of 'input_files/' directory.
#' @param genome.name variable of genome name defined in this RNA-Seq workflow
#'   (ex. \code{genome.name}.fa, \code{genome.name}.gtf).
#' @param sample.pattern  Regular expression of paired-end fastq.gz files under
#'   'input_files/raw_bam'. Expression not includes \code{_[1,2].fastq.gz}.
#' @param independent.variable Independent variable for the biological
#'   experiment design of two-group RNA-Seq workflow.
#' @param case.group Group name of the case group.
#' @param control.group Group name of the control group.
#'
#' @return an object of class \code{RNASeqRParam_Bam}
#'
#' @author Kuan-Hao Chao
#'
#' @rdname RNASeqRParam_Bam-constructor
#'
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' input_files.path <- system.file("extdata/", package = "RNASeqRData")
#' rnaseq_result.path <- "/Users/chaokuan-hao/TRY"
#' # rnaseq_result.path <- tempdir(check = TRUE)
#' \dontrun{
#' exp <- RNASeqRParam_Bam(path.prefix          = rnaseq_result.path,
#'                         input.path.prefix    = input_files.path,
#'                         genome.name          = "Saccharomyces_cerevisiae_XV_Ensembl",
#'                         sample.pattern       = "SRR[0-9]*_XV",
#'                         independent.variable = "state",
#'                         case.group           = "60mins_ID20_amphotericin_B",
#'                         control.group        = "60mins_ID20_control")
#' }
RNASeqRParam_Bam <- function(path.prefix          = NA,
                             input.path.prefix    = NA,
                             genome.name          = NA,
                             sample.pattern       = NA,
                             independent.variable = NA,
                             case.group           = NA,
                             control.group        = NA) {
  # check input parameters
  CheckInputParamNa(path.prefix,
                    input.path.prefix,
                    genome.name,
                    sample.pattern,
                    independent.variable,
                    case.group,
                    control.group,
                    "")
  # 1. check operating system
  characters.os.type <- CheckOperatingSystem()
  # 2. check python version
  python.version.list <- CheckPython()
  bool.python.avail <- python.version.list$check.answer
  two.to.three.result <- Check2to3()
  # 3. check validity of path.prefix
  bool.prefix.path <- CheckPrefixPath(path.prefix = path.prefix)
  if (bool.prefix.path){
    # add '/' to the path.prefix
    if (substr(path.prefix, nchar(path.prefix), nchar(path.prefix)) != "/") {
      path.prefix <- paste0(path.prefix, "/")
    }
  }
  # 4. check input.path.prefix
  bool.input.path.prefix <- CheckInputPrefixPath(input.path.prefix)
  if (bool.input.path.prefix){
    # add '/' to the path.prefix
    if (substr(input.path.prefix,
               nchar(input.path.prefix),
               nchar(input.path.prefix)) != "/") {
      input.path.prefix <- paste0(input.path.prefix, "/")
    }
  }
  # 5. check 'phenodata'
  bool.phenodata <- CheckPhenodata_Bam(input.path.prefix,
                                       genome.name,
                                       sample.pattern,
                                       independent.variable)
  # 6. check 'input_files/' necessary files with 'genome.name', 'sample.pattern'
  input.dir.files.list <- CheckInputDirFiles_Bam(input.path.prefix,
                                                 genome.name,
                                                 sample.pattern)
  bool.input.dir.files <- input.dir.files.list
  # This determine whether to run 'CreateHisat2Index'

  # 7. check 'case.group' and 'control.group'
  bool.control.control.group <- CheckCaseControlGroup(input.path.prefix,
                                                      independent.variable,
                                                      case.group, control.group)

  if ( (characters.os.type == "linux" || characters.os.type == "osx") &&
       bool.python.avail && bool.prefix.path && bool.input.path.prefix &&
       bool.input.dir.files && bool.phenodata && bool.control.control.group) {
    message("\n**************************************\n")
    message("************** Success! **************\n")
    message("**************************************\n")
    new("RNASeqRParam_Bam",
        os.type              = characters.os.type,
        python.variable      = python.version.list,
        python.2to3          = two.to.three.result,
        path.prefix          = path.prefix,
        input.path.prefix    = input.path.prefix,
        genome.name          = genome.name,
        sample.pattern       = sample.pattern,
        independent.variable = independent.variable,
        case.group           = case.group,
        control.group        = control.group)
  }
}


# inner function : check input.path.prefix
CheckInputDirFiles_Bam <- function(input.path.prefix, genome.name, sample.pattern) {
  message("************** Checking hierarchy of '",
          input.path.prefix, "input_files/' ************\n")
  # only check whether exist
  gtf.file <- file.exists(paste0(input.path.prefix,
                                 "input_files/",
                                 genome.name, ".gtf"))
  # check exist and rules
  raw.bam.dir <- dir.exists(paste0(input.path.prefix,
                                   "input_files/raw_bam/"))
  # check exist and rules
  phenodata.file <- file.exists(paste0(input.path.prefix,
                                       "input_files/phenodata.csv"))
  # check whether sample pattern matches the file names~
  if (raw.bam.dir) {
    raw.bam <- list.files(path = paste0(input.path.prefix,
                                        "input_files/raw_bam/"),
                          pattern = "*.bam$",
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
    raw.bam.number <- length(raw.bam)
    if (raw.bam.number != 0) {
      message("(\u2714) : '*.bam' is in 'input_files'\n")
    } else {
      message("(\u231B) : '", input.path.prefix,
              "input_files/raw_bam/XXX.bam' is not exit\n")
    }
  } else {
    message("(\u231B) : '", input.path.prefix,
            "input_files/raw_bam/XXX.bam' is not exit\n")
  }
  if (!isTRUE(gtf.file)) {
    grab.gtf.file <- Sys.glob(file.path(path = paste0(input.path.prefix,
                                                      "input_files"), "*.gtf"))
    message("(\u2718) : '", genome.name,
            ".gtf (user input)' and '", grab.gtf.file,
            " (find in directory)' ", " are mismatched.\n")
  } else {
    message("(\u2714) : '", genome.name,
            ".gtf'", " is in 'input_files'\n")
  }
  if (!isTRUE(phenodata.file)) {
    message("(\u2718) : '", "phenodata.csv is missing.\n")
  } else {
    message("(\u2714) : '", "phenodata.csv is in 'input_files'\n")
  }
  if (gtf.file && phenodata.file) {
    message("\n(\u2714) : '", input.path.prefix,
            "input_files/", "' is valid !\n")
    message("\n")
    return.value <- TRUE
    return(return.value)
  } else {
    stop("'input_files/' checking ERROR")
  }
}


# inner function : check validity of phenodata
CheckPhenodata_Bam <- function(input.path.prefix,
                               genome.name,
                               sample.pattern,
                               independent.variable) {
  # have to sort the column !! and sort them in the correct order
  message("************** Checking phenodata  ************\n")
  pheno_data <- read.csv(paste0(input.path.prefix, "input_files/phenodata.csv"))
  # Covert all column to character
  pheno_data <- data.frame(lapply(pheno_data, as.character),
                           stringsAsFactors = FALSE)
  # Check 'ids' is in the 'phenodata.csv'
  if (!("ids" %in% colnames(pheno_data))) {
    message("(\u2718) : 'ids' can't find in the ",
            "column of phenodata.csv.\n\n")
    stop("'ids' invalid ERROR")
  }
  # "id" : must be distinct, same as input_files raw reads name !
  message("\u25B6 Checking whether \"raw_bam files\" ",
          "matches \"'ids' of phenodata.csv\" \n")
  raw.bam <- list.files(path = paste0(input.path.prefix,
                                      "input_files/raw_bam/"),
                        pattern = sample.pattern,
                        all.files = FALSE,
                        full.names = FALSE,
                        recursive = FALSE,
                        ignore.case = FALSE)
  extract.bam.sample.names <- unique(gsub(".bam", "", raw.bam))
  bool.length <- length(extract.bam.sample.names) == length(pheno_data$ids)
  bool.identical <- identical(sort(extract.bam.sample.names),
                              sort(pheno_data$ids))
  if (!bool.length || !bool.identical) {
    message("(\u2718) : 'ids' column doesn't match the smaple_id in ",
            "'input_files/raw_bam'. Please check the file and ",
            "'input_files/raw_bam' directory.\n" )
    stop("'ids' mismatch ERROR")
  }
  ids.list <- paste(sort(extract.bam.sample.names), collapse = ", ")
  message("(\u2714) : Column 'ids' of phenodata.csv is valid. \n")
  message("      \u25CF sample ids are : \"", ids.list, "\"\n")
  # Check again independent.variable in the list
  if (!(independent.variable %in% colnames(pheno_data))) {
    message("(\u2718) : 'independent.variable' : '",
            independent.variable,
            "' can't find in the column of phenodata.csv.\n\n")
    stop("'independent.variable' invalid ERROR")
  }
  message("(\u2714) : 'independent.variable' : '",
          independent.variable,
          "' is in the column of phenodata.csv. \n\n")
  message("\u25B6 Checking whether '",
          independent.variable,
          "' is a two-group 'independent.variable' ...\n")
  length.independent.variable <- length(table(pheno_data[independent.variable]))
  if (!(length.independent.variable == 2)) {
    message("(\u2718) : 'independent.variable' : '",
            independent.variable, "' is a ",
            length.independent.variable,
            "-group 'independent.variable'. Not 2-group.\n")
    message("          groups that found : ",
            paste0(names(table(pheno_data[independent.variable])),
                   collapse = ", "), "\n")
    stop("'independent.variable' none-two-group ERROR")
  }
  message("(\u2714) : Column 'independent.variable' : '",
          independent.variable, "' of phenodata.csv is valid. \n")
  message("      \u25CF 'independent.variable' : '",
          independent.variable, "'\n\n")
  return(TRUE)
}
