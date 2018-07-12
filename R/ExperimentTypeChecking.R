#'
#' Check the experiment type and seperate into three function
#'
#' @export
ExperimentTypeChecking <- function(path.prefix, gene.name, sample.pattern, experiment.type, main.variable, additional.variable) {
  # Check the how many group !
  results <- ProgressGenesFiles(path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print = FALSE)
  if (isTRUE(results$phenodata.file.df)){
    # sorting 'pheno_data'
    cat(paste0("************** Experiment group checking **************\n"))
    cat("\u25CF 1. Printing origin phenodata.csv : \n")
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    print(pheno_data)
    cat('\n')
    sample.table <- as.data.frame(table(pheno_data[main.variable]))
    groups.number <- length(row.names(sample.table))
    if (groups.number == 1) {
      stop("Group number ERROR")
    }
    if (experiment.type == "two.group") {
      if (groups.number != 2) {
        stop("Group number ERROR")
      }
      TwoGroupAnalysis(path.prefix, gene.name, sample.pattern, main.variable, additional.variable)
    } else if (experiment.type == "multi.group.pairs") {
      if (groups.number == 2) {
        stop("Group number ERROR")
      }
      MultiGroupPairsAnalysis(path.prefix, gene.name, sample.pattern, main.variable, additional.variable)
    } else if (experiment.type == "multi.group.anova") {
      if (groups.number == 2) {
        stop("Group number ERROR")
      }
      MultiGroupAnovaAnalysis(path.prefix, gene.name, sample.pattern, main.variable, additional.variable)
    }
  }
}

#' @export
TwoGroupAnalysis <- function(path.prefix, gene.name, sample.pattern, main.variable, additional.variable) {
  print(c("path.prefix : ", path.prefix))
  print(c("gene.name : ", gene.name))
  print(c("sample.pattern : ", sample.pattern))
  print(c("main.variable : ", main.variable))
  print(c("additional.variable : ", additional.variable))
  cat(paste0("Experiment type : two.group"))
  # seperate into three part analysis : FPKM, TPM, raw reads count

}

#' @export
MultiGroupPairsAnalysis <- function(path.prefix, gene.name, sample.pattern, main.variable, additional.variable) {
  print(c("path.prefix : ", path.prefix))
  print(c("gene.name : ", gene.name))
  print(c("sample.pattern : ", sample.pattern))
  print(c("main.variable : ", main.variable))
  print(c("additional.variable : ", additional.variable))
}

#' @export
MultiGroupAnovaAnalysis <- function(path.prefix, gene.name, sample.pattern, main.variable, additional.variable) {
  print(c("path.prefix : ", path.prefix))
  print(c("gene.name : ", gene.name))
  print(c("sample.pattern : ", sample.pattern))
  print(c("main.variable : ", main.variable))
  print(c("additional.variable : ", additional.variable))
}
