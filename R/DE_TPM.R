# TPM normalization
TPMNormalizationAnalysis <- function() {
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/TPM_Ttest_analysis/"))){
    dir.create(paste0(path.prefix, "RNASeq_results/TPM_Ttest_analysis/"))
  }
  if(!dir.exists(paste0(path.prefix, "RNASeq_results/TPM_Ttest_analysis/normalized_&_statistic"))){
    dir.create(paste0(path.prefix, "RNASeq_results/TPM_Ttest_analysis/normalized_&_statistic"))
  }
  control.FPKM <- read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/FPKM_control.csv"))
  experiment.FPKM <- read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/FPKM_experiment.csv"))
  statistic.FPKM <- read.csv(paste0(path.prefix, "RNASeq_results/ballgown_analysis/normalized_&_statistic/statistic.csv"))

  control.TPM <- (control.FPKM/colSums(control.FPKM))*(10**6)

  experiment.TPM <- (experiment.FPKM/colSums(control.FPKM))*(10**6)

  p.value <- unlist(lapply(1:nrow(control.TPM), function(x) { t.test(control.TPM[x,], experiment.TPM[x,])$p.value }))
  fold.change <- unlist(lapply(1:nrow(control.TPM), function(x) { mean(unlist(experiment.TPM[x,])) / mean(unlist(control.TPM[x,])) }))

  statistic <- data.frame("pval" = p.value, "fc" = fold.change, "log2FC" = log2(fold.change))

  data.frame(p.value)

  p1$p.value
}
