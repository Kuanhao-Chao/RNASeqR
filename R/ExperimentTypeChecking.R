#'
#' Check the experiment type and seperate into three function
#'
#' @export
ExperimentTypeChecking <- function(path.prefix, gene.name, sample.pattern, experiment.type, main.variable, additional.variable) {
  # Check the how many group !
  results <- ProgressGenesFiles(path.prefix, gene.name = gene.name, sample.pattern = sample.pattern, print = FALSE)
  if (isTRUE(results$phenodata.file.df)){
    # sorting 'pheno_data'
    cat(paste0("************** Ballgown data preprocessing **************\n"))
    cat("\u25CF 1. Printing origin phenodata.csv : \n")
    pheno_data <- read.csv(paste0(path.prefix, "gene_data/phenodata.csv"))
    print(pheno_data)
    cat('\n')
    sample.table <- as.data.frame(table(pheno_data[main.variable]))
    groups.number <- length(row.names(sample.table))
    if (experiment.type == "two.group") {
      if (groups.number != 2) {
        stop("Group number ERROR")
      }
    } else if (experiment.type == "multi.group.pairs") {

    } else if (experiment.type == "multi.group.anova") {

    }
  }
}
