#' Create 'RNAseqPipline.R' that user can
#' @export
RNAseqPipelineCMD <- function(path.prefix = "NOT_SET_YET", input.path.prefix = "NOT_SET_YET", gene.name = "NO_DATA", sample.pattern = "NO_DATA", num.parallel.threads = 8) {
  if (path.prefix == "NOT_SET_YET" || input.path.prefix == "NOT_SET_YET" ||gene.name == "NO_DATA" || sample.pattern == "NO_DATA"){
    if (path.prefix == "NOT_SET_YET") {
      cat("(\u2718) : 'path.prefix' is missing.\n\n")
    }
    if (input.path.prefix == "NOT_SET_YET") {
      cat("(\u2718) : 'input.path.prefix' is missing.\n\n")
    }
    if (gene.name == "NO_DATA") {
      cat("(\u2718) : 'gene.name' is missing.\n\n")
    }
    if (sample.pattern == "NO_DATA") {
      cat("(\u2718) : 'sample.pattern' is missing.\n\n")
    }
  } else {
    if (isTRUE(SetPrefixPath(path.prefix = path.prefix, print = FALSE))) {
      # not print but if the prefix is invalid, then 'Prefix path '", path.prefix, "' is invalid. Please try another one.' will be printed.
      if (isTRUE(CheckPrefixPath(path.prefix = pkg.global.path.prefix$data_path, print = TRUE))) {
        if (isTRUE(CheckDirAll(print = TRUE))) {
          results <- ProgressGenesFiles(gene.name = gene.name, sample.pattern = sample.pattern, print=TRUE)
          if (isTRUE(results$gtf.file.logic.df) && isTRUE(results$fa.file.logic.df) && results$fastq.gz.files.number.df != 0 && isTRUE(results$phenodata.file.df) && (results$phenodata.invalid.column.number.df == 0)) {
            # If precheck doesn't have .ht2 files is fine
            ExportPath()
            if (isTRUE(CheckToolAll(print=TRUE))) {
              cat("(\u2714) : Successful in RNAseq-pipeline precheck. \n\n")
              fileConn<-file(paste0(pkg.global.path.prefix$data_path, "Rscript/RNASEQ_PIPELINE.R"))
              first <- "library(RNASeq)"
              second <- paste0('RNAseqPipeline(path.prefix = "', path.prefix, '", input.path.prefix = "', input.path.prefix, '", gene.name = "', gene.name, '", sample.pattern = "', sample.pattern, '", num.parallel.threads = "', num.parallel.threads, '")')
              writeLines(c(first, second), fileConn)
              close(fileConn)
              system2(command = 'nohup', args = paste0("R CMD BATCH ", pkg.global.path.prefix$data_path, "/Rscript/RNASEQ_PIPELINE.R ", pkg.global.path.prefix$data_path, "/Rscript_out/RNASEQ_PIPELINE.Rout"), stdout = "", wait = FALSE)
              cat(paste0("\u2605 RNAseq alignment, assembly, mergence, comparison, reads preprocess are doing in the background. Check current progress in '", pkg.global.path.prefix$data_path, "/Rscript_out/RNASEQ_PIPELINE.Rout'\n\n"))
            }
          }
        }
      }
    }
  }
}


#' rna seq pipline
#' @export
RNAseqPipeline <- function(path.prefix = "NOT_SET_YET", input.path.prefix = "NOT_SET_YET", gene.name = "NO_DATA", sample.pattern = "NO_DATA", num.parallel.threads = 8) {
  if (isTRUE(SetPrefixPath(path.prefix))){
    if (isTRUE(CheckDirAll(print = TRUE))){
      if (gene.name == "NO_DATA" || sample.pattern == "NO_DATA"){
        if (gene.name == "NO_DATA") {
          cat("(\u2718) : gene.name is missing.\n\n")
        }
        if (sample.pattern == "NO_DATA") {
          cat("(\u2718) : sample.pattern is missing.\n\n")
        }
      } else {
        if (isTRUE(CheckInputDirFiles(input.path.prefix = input.path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print=TRUE))) {
          ExportPath()
          check.results <- ProgressGenesFiles(gene.name = gene.name, sample.pattern = sample.pattern, print=FALSE)
          if (isTRUE(check.results$gtf.file.logic.df) && isTRUE(check.results$fa.file.logic.df) && (check.results$fastq.gz.files.number.df != 0)){
            if (check.results$ht2.files.number.df == 0) {
              CreateHisat2Index(gene.name, sample.pattern)
            }
            Hisat2AlignmentDefault(gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
            SamtoolsToBam(gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
            StringTieAssemble(gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
            StringTieMergeTrans(gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
            GffcompareRefSample(gene.name, sample.pattern)
            StringTieToBallgown(gene.name, sample.pattern, num.parallel.threads = num.parallel.threads)
            finals <- ProgressGenesFiles(gene.name, sample.pattern, print=TRUE)

            if (isTRUE(finals$gtf.file.logic.df) && isTRUE(finals$fa.file.logic.df) &&
                finals$fastq.gz.files.number.df != 0 &&
                isTRUE(finals$phenodata.file.df) &&
                finals$phenodata.invalid.column.number.df == 0 &&
                finals$ht2.files.number.df != 0 &&
                finals$sam.files.number.df != 0 &&
                finals$bam.files.number.df != 0 &&
                finals$gtf.files.number.df != 0 &&
                isTRUE(finals$stringtie_merged.gtf.file.df) &&
                finals$gffcompare.related.dirs.number.df != 0 &&
                finals$ballgown.dirs.number.df != 0) {
              cat(paste0("\n**************************************\n"))
              cat(paste0("************** Success! **************\n"))
              cat(paste0("**************************************\n"))
            }
          }
          else{
            stop(paste0("(\u2718) Necessary files are lost.\n     Please check whether 'ref_genes/", gene.name, ".gtf' , 'ref_genome/", gene.name, ".fa' , 'samples_.fastq.gz/XXX_", gene.name, "_*.fastq.gz' are exit.\n\n" ))
          }
        }
      }
    }
  }
}
