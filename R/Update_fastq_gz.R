#' @title Update_Fastq_gz
#'
#' @description
#'   This function let users update their trimmed fastq.gz files automatically.
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param prepared_fastq_gz absolute path to the prepared 'raw_fastq.gz' 
#'  directory.
#' @param target_samples list of samples that are going to update. Default
#'  value is \code{ALL}
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqDifferentialAnalysis(RNASeqRParam = yeast)}
Update_Fastq_gz <- function(RNASeqRParam, 
                            prepared_fastq_gz, 
                            target_samples     = "ALL") {
  path.prefix <- RNASeqRParam@path.prefix
  independent.variable <- RNASeqRParam@independent.variable
  case.group <- RNASeqRParam@case.group
  control.group <- RNASeqRParam@control.group
  bool_prepared_fastq_gz <- dir.exists(prepared_fastq_gz)
  if(bool_prepared_fastq_gz) {
    if (target_samples != "ALL") {
      Checker_result <- TRUE
      phenoData.result<- phenoDataWrap(path.prefix,
                                       independent.variable,
                                       case.group,
                                       control.group)
      all.samples <- levels(phenoData.result$pheno_data$ids)
      for (i in target_samples) {
        each_checker <- any(i == all.samples)
        if (!each_checker) {
          Checker_result <- FALSE
        }
      }
      if (Checker_result) {
        # Input samples list is valid. Now check files are inside input_files
        prepared.files <- list.files(prepared_fastq_gz)
        for (i in target_samples) {
          message("Checking '", i, "' sample validity ...")
          each_checker_1 <- any(paste0(i, "_1.fastq.gz") == prepared.files)
          each_checker_2 <- any(paste0(i, "_2.fastq.gz") == prepared.files)
          if (!each_checker_1 || !each_checker_2) {
            message(paste0("target_samples list is invalid. Please check all target ",
                        "samples are inside 'phenodata.csv'."))
            stop(paste0(i, " file does not match. Please check your target", 
                        " file naming is correct."))
          }
        }
        # Pass every check! Now update files.
        for (i in target_samples) {
          message("Moving '", i, "' sample into project ...")
          unlink(paste0(path.prefix, "gene_data/raw_fastq.gz/", i, "_1.fastq.gz"))
          unlink(paste0(path.prefix, "gene_data/raw_fastq.gz/", i, "_2.fastq.gz"))
          paste0(prepared_fastq_gz, '/', i, "_1.fastq.gz")
          file.symlink(paste0(prepared_fastq_gz, '/', i, "_1.fastq.gz"),
                       paste0(path.prefix, "gene_data/raw_fastq.gz"))
          file.symlink(paste0(prepared_fastq_gz, '/', i, "_2.fastq.gz"),
                       paste0(path.prefix, "gene_data/raw_fastq.gz"))
          message("    ", i, " fastq.gz files are updated.")
        }
        message("Fastq.gz files are updated!")
      } else {
        stop(paste0("target_samples list is invalid. Please check all target ",
                    "samples are inside 'phenodata.csv'."))
      }
    } else if (target_samples == "ALL") {
      message("Checking all samples validity ...")
      prepared.files <- list.files(prepared_fastq_gz)
      origin.files <- list.files(paste0(path.prefix, "gene_data/raw_fastq.gz"))
      if (prepared.files == origin.files) {
        different_number <- length(setdiff(prepared.files, origin.files))
        if(different_number == 0) {
          message("   Samples are valid. Start to update .fastq.gz files")
          
          
          unlink(paste0(path.prefix, "gene_data/raw_fastq.gz/"), recursive = TRUE)
          
          dir.create(paste0(path.prefix, "gene_data/raw_fastq.gz/"))
          raw_fastq.gz.subfiles <- list.files(path = prepared_fastq_gz,
                                              pattern = sample.pattern,
                                              recursive = TRUE,
                                              full.names = TRUE)
          vapply(raw_fastq.gz.subfiles,
                 function(x) file.symlink(x, paste0(path.prefix,
                                                    "gene_data/raw_fastq.gz")),
                 FUN.VALUE = TRUE)
          message("All fastq.gz files are updated.")
        } else {
          stop("   Samples do not match. Please check your fastq.gz files in",
               " 'prepared_fastq_gz':", prepared_fastq_gz)
        }
      } 
    }
  }
}
