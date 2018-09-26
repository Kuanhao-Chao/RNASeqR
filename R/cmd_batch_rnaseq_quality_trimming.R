#' @title RNASeqQualityTrimming_CMD
#'
#' @description
#'   Trim '.fastq.gz' files for RNA-Seq workflow in background.
#'   This step is optional in the whole RNA-Seq workflow.
#'   The trimming method is implemented by R package \code{ShortRead} \cr
#'   If you want to trim '.fastq.gz' files for the RNA-Seq workflow in R shell,
#'   please see \code{RNASeqQualityTrimming()} function.
#'
#' @param RNASeqRParam S4 object instance of
#'   experiment-related parameters
#' @param cum.error Default \code{1}.
#'   Cut of threshold of cumulative probability of error per base.
#' @param trimming.position Default \code{NA}. If value is \code{NA}, trimming
#'   points of all paired-end files will be calculated by CDF of error rate,
#'   while they will all be the value of \code{trimming.position}.
#' @param reads.length.limit Default \code{36}.
#'   The shortest base pair length of short reads
#' @param run Default value is \code{TRUE}.
#'   If \code{TRUE}, 'Rscript/Environment_Set.R' will be created and executed.
#'   The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{False}, 'Rscript/Environment_Set.R'
#'   will be created without executed.
#' @param check.s4.print Default \code{TRUE}.
#'   If \code{TRUE}, the result of checking \code{RNASeqRParam}
#'   will be reported in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqQualityTrimming_CMD(RNASeqRParam = yeast)}
RNASeqQualityTrimming_CMD <- function(RNASeqRParam,
                                      cum.error          = 1,
                                      trimming.position  = NA,
                                      reads.length.limit = 36,
                                      run                = TRUE,
                                      check.s4.print     = TRUE) {
  # check input param
  CheckS4Object(RNASeqRParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  INSIDE.path.prefix <- "@"(RNASeqRParam, path.prefix)
  saveRDS(RNASeqRParam,
          file = paste0(INSIDE.path.prefix,
                        "gene_data/RNASeqRParam.rds"))
  fileConn <- file(paste0(path.prefix, "Rscript/Quality_Trimming.R"))
  first <- "library(RNASeqR)"
  second <- paste0("RNASeqQualityTrimming(RNASeqRParam = 'INSIDE'",
                   ", which.trigger = 'INSIDE'",
                   ", INSIDE.path.prefix = '", INSIDE.path.prefix,
                   "', cum.error = ", cum.error,
                   ", trimming.position = ", trimming.position,
                   ", reads.length.limit = ", reads.length.limit, ")")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message(paste0("\u2605 '", path.prefix,
                 "Rscript/Quality_Trimming.R' has been created.\n"))
  if (run) {
    R.home.lib <- R.home()
    R.home.bin <- gsub("/lib/R", "/bin/R", R.home.lib)
    system2(command = "nohup",
            args = paste0(R.home.bin, " CMD BATCH ",
                          path.prefix,
                          "Rscript/Quality_Trimming.R ", path.prefix,
                          "Rscript_out/Quality_Trimming.Rout"),
            stdout = "",
            wait = FALSE)
    message(paste0("\u2605 Tools are installing in the background. ",
                   "Check current progress in '", path.prefix,
                   "Rscript_out/Quality_Trimming.Rout'\n\n"))
  }
}


#' @title RNASeqQualityTrimming
#'
#' @description
#'   Trim '.fastq.gz' files for RNA-Seq workflow in R shell.
#'   This step is optional in the whole RNA-Seq workflow.
#'   It is strongly advised to run \code{RNASeqQualityTrimming_CMD()} directly.
#'   Running \code{RNASeqQualityTrimming_CMD()} will create
#'   'Quality_Trimming.Rout' file in 'Rscript_out/' directory.
#'   The trimming method is implemented by R package \code{ShortReads} \cr
#'   If you want to trim '.fastq.gz' files for the RNA-Seq workflow
#'   in background, please see \code{RNASeqQualityTrimming_CMD()} function.
#'
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param which.trigger Default value is \code{OUTSIDE}. User should not change
#'   this value.
#' @param INSIDE.path.prefix Default value is \code{NA}. User should not change
#'   this value.
#' @param cum.error Default \code{1}.
#'   Cut of threshold of cumulative probability of error per base.
#' @param trimming.position Default \code{NA}. If value is \code{NA}, trimming
#'   points of all paired-end files will be calculated by CDF of error rate,
#'   while they will all be the value of \code{trimming.position}.
#' @param reads.length.limit Default \code{36}.
#'   The shortest base pair length of short reads
#' @param check.s4.print Default \code{TRUE}.
#'   If \code{TRUE}, the result of checking \code{RNASeqRParam}
#'   will be reported in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqQualityTrimming(RNASeqRParam = yeast)}
RNASeqQualityTrimming <- function(RNASeqRParam,
                                  which.trigger      = "OUTSIDE",
                                  INSIDE.path.prefix = NA,
                                  cum.error = 1,
                                  trimming.position  = NA,
                                  reads.length.limit = 36,
                                  check.s4.print     = TRUE) {
  CheckOperatingSystem(FALSE)
  # If `which.trigger` is OUTSIDE, then directory must be built
  # If `which.trigger` is INSIDE, then directory must not be
  #  built here(will created in CMD)
  if (isS4(RNASeqRParam) &
      which.trigger == "OUTSIDE" &
      is.na(INSIDE.path.prefix)) {
    # This is an external call!!
    # Check the S4 object(user input)
    CheckS4Object(RNASeqRParam, check.s4.print)
  } else if (RNASeqRParam == "INSIDE" &
             which.trigger == "INSIDE" &
             !is.na(INSIDE.path.prefix)) {
    # This is an internal call!!
    # Load the S4 object that saved in CMD process
    RNASeqRParam <- readRDS(paste0(INSIDE.path.prefix,
                                   "gene_data/RNASeqRParam.rds"))
  }
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  sample.pattern <- "@"(RNASeqRParam, sample.pattern)
  PreCheckRNASeqQualityTrimming(path.prefix = path.prefix,
                                sample.pattern = sample.pattern)
  message(paste0("************** Quality Trimming **************\n"))
  if (!dir.exists(paste0(path.prefix,
                        "gene_data/raw_fastq.gz/",
                        "original_untrimmed_fastq.gz/"))){
    dir.create(file.path(paste0(path.prefix,
                                "gene_data/raw_fastq.gz/",
                                "original_untrimmed_fastq.gz/")),
               showWarnings = FALSE)
  }
  raw.fastq <- list.files(path = paste0(path.prefix, "gene_data/raw_fastq.gz"),
                          pattern = sample.pattern,
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
  raw.fastq.unique <- unique(gsub("_[1-2]*.fastq.gz$",
                                  replacement = "", raw.fastq))
  lapply(raw.fastq.unique, myFilterAndTrim,
         path.prefix = path.prefix,
         cum.error = cum.error,
         trimming.position  = trimming.position,
         reads.length.limit = reads.length.limit)
  message("\n")
  PostCheckRNASeqQualityTrimming(path.prefix, sample.pattern)
}


myFilterAndTrim <- function(fl.name,
                            path.prefix,
                            cum.error,
                            trimming.position,
                            reads.length.limit) {
  # adding print log
  # file1 and file2 is original fastq.gz without trimmed
  message(paste0("\u25CF \"", fl.name, "\" quality trimming\n"))
  file1 <- paste0(path.prefix,
                  "gene_data/raw_fastq.gz/", fl.name, "_1.fastq.gz")
  file2 <- paste0(path.prefix,
                  "gene_data/raw_fastq.gz/", fl.name, "_2.fastq.gz")
  if (file.exists(file1) && file.exists(file2)) {
    # file1.output and file2.output are the new original fastq.gz file name
    file1.untrimmed <- paste0(path.prefix,
                              "gene_data/raw_fastq.gz/",
                              "original_untrimmed_fastq.gz/",
                              fl.name, "_1.fastq.gz")
    file2.untrimmed <- paste0(path.prefix, "gene_data/raw_fastq.gz/",
                              "original_untrimmed_fastq.gz/",
                              fl.name, "_2.fastq.gz")
    message(paste0("     \u25CF Moving \"", basename(file1), "\" to \"",
                   path.prefix,
                   "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/\"\n"))
    message(paste0("     \u25CF Moving \"", basename(file2), "\" to \"",
                   path.prefix,
                   "gene_data/raw_fastq.gz/original_untrimmed_fastq.gz/\"\n"))
    file.rename(from = file1, to = file1.untrimmed)
    file.rename(from = file2, to = file2.untrimmed)
    # Sequence complexity (H) is calculated based on the dinucleotide
    # composition using the formula (Shannon entropy):
    message(paste0("     \u25CF Reading \"", basename(file1.untrimmed),
                   "\" file in R, please wait\n"))
    message(paste0("     \u25CF Reading '", basename(file2.untrimmed),
                   "' file in R, please wait\n"))
    file1.read <- ShortRead::readFastq(file1.untrimmed)
    file2.read <- ShortRead::readFastq(file2.untrimmed)
    if (length(file1.read) == length(file2.read)) {
      total.reads.number <- length(file1.read)
    } else {
      message(paste0("'", fl.name, "_1.fastq.gz'"),
              " reads number: ", length(file1.read), " bps")
      message(paste0("'", fl.name, "_2.fastq.gz'"),
              " reads number: ", length(file2.read), " bps")
      message(paste0("'", fl.name, "_1.fastq.gz'"), " and ",
              paste0("'", fl.name, "_2.fastq.gz'"),
              " have different number of reads number!!\n")
      stop("paired-end files have different length!!")
    }
    max.width <- max(Biostrings::width(file1.read),
                     Biostrings::width(file2.read))
    message(paste0("     \u25CF Start trimming (max width: ",
                   max.width, " bps) ...\n"))
    if (is.na(trimming.position)) {
      # It means that users didn't provide trimming position !!
      message(paste0("          \u25CF Getting quality score ",
                     "list as PhredQuality ...\n"))
      # get quality score list as PhredQuality
      qual1 <- as(Biostrings::quality(file1.read), "matrix")
      qual2 <- as(Biostrings::quality(file2.read), "matrix")

      # Calculate probability error per base (through column) ==>
      #   Q = -10log10(P)  or  P = 10^(-Q/10)
      message(paste0("          \u25CF Calculating probability ",
                     "error per base ...\n"))
      pe1 <- apply(qual1, MARGIN = 2, function(x){10^(-(x/10))})
      pe2 <- apply(qual2, MARGIN = 2, function(x){10^(-(x/10))})
      # Calculate cpm of error
      message(paste0("          \u25CF Calculating cumulative distribution ",
                     "probability of error per base ...\n"))
      cum.pe1 <- apply(pe1, MARGIN = 1, cumsum)
      cum.pe2 <- apply(pe2, MARGIN = 1, cumsum)

      # Get the trimming position of each file
      message(paste0("          \u25CF Filtering out cumulative distribution ",
                     "probability of error per base < 1 ...\n"))
      trimPos1 <- apply(cum.pe1, 2, function(x) { min(min(which(x > cum.error)),
                                                      length(x)) } )
      trimPos2 <- apply(cum.pe2, 2, function(x) { min(min(which(x > cum.error)),
                                                      length(x)) } )
      # Get the trimPos for pair-end files
      message(paste0("          \u25CF Finding trimming ",
                     "position for paired-end ...\n"))
      trimPos.together <- mapply(function(list1, list2) {min(list1, list2)},
                                 list1 = trimPos1,
                                 list2 = trimPos2)
    } else {
      # trimming.position is not na!!!!
      if (trimming.position%%1 == 0 &
          trimming.position > reads.length.limit &
          trimming.position <= max.width) {
        # Valid trimming.position
        total.reads.number
        trimPos.together <- rep(trimming.position, total.reads.number)
      } else {
        message("Invalid trimming.position: ", trimming.position, " !!")
        stop("'trimming.position' ERROR")
      }
    }
    trimmed.file1 <- ShortRead::narrow(x = file1.read,
                                       start = 1,
                                       end = trimPos.together)
    trimmed.file2 <- ShortRead::narrow(x = file2.read,
                                       start = 1,
                                       end = trimPos.together)

    ## drop reads that are less than 36nt
    message(paste0("     \u25CF Removing reads that are less than ",
                   reads.length.limit, " base pairs ...\n"))
    trimmed.file1 <- trimmed.file1[ShortRead::width(trimmed.file1) >=
                                     reads.length.limit]
    trimmed.file2 <- trimmed.file2[ShortRead::width(trimmed.file2) >=
                                     reads.length.limit]

    # write new fastaq files
    message(paste0("     \u25CF Creating trimmed pair-end files ...\n"))
    ShortRead::writeFastq(trimmed.file1, file1, "w")
    ShortRead::writeFastq(trimmed.file2, file2, "w")

    message(paste0("     \u25CF \"", file1, "\" has been created.\n"))
    message(paste0("     \u25CF \"", file2, "\" has been created.\n\n"))
  } else {
    stop("paired-end file ERROR")
  }
}

PreCheckRNASeqQualityTrimming <- function(path.prefix, sample.pattern) {
  message("\u269C\u265C\u265C\u265C 'RNASeqQualityTrimming()' ",
          "environment pre-check ...\n")
  # have fastq.gz files
  raw.fastq <- list.files(path = paste0(path.prefix, 'gene_data/raw_fastq.gz/'),
                          pattern = sample.pattern,
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
  validity <- length(raw.fastq) != 0
  if (!isTRUE(validity)) {
    stop("CheckRNASeqQualityTrimming() pre-check ERROR")
  }
  message("(\u2714) : RNASeqQualityTrimming() pre-check is valid\n\n")
}

PostCheckRNASeqQualityTrimming <- function(path.prefix, sample.pattern) {
  message("\u269C\u265C\u265C\u265C 'RNASeqQualityTrimming()' ",
          "environment post-check ...\n")
  # have fastq.gz and trimmed fastq.gz files
  trimmed.raw.fastq <- list.files(path = paste0(path.prefix,
                                                'gene_data/raw_fastq.gz/'),
                                  pattern = sample.pattern,
                                  all.files = FALSE,
                                  full.names = FALSE,
                                  recursive = FALSE,
                                  ignore.case = FALSE)
  raw.fastq <- list.files(path = paste0(path.prefix,
                                        "gene_data/raw_fastq.gz/",
                                        "original_untrimmed_fastq.gz/"),
                          pattern = sample.pattern,
                          all.files = FALSE,
                          full.names = FALSE,
                          recursive = FALSE,
                          ignore.case = FALSE)
  validity <- (length(trimmed.raw.fastq) != 0) && (length(raw.fastq) != 0)
  if (!isTRUE(validity)) {
    stop("RNASeqQualityTrimming() post-check ERROR")
  }
  message("(\u2714) : RNASeqQualityTrimming() post-check is valid\n\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
}
