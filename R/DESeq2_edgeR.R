DESeq2edgeRRawCountAnalysis <- function(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC) {
  cat(paste0("\n************** Reads count matrix analysis **************\n"))
  if(!file.exists(paste0(path.prefix, "gene_data/reads_count_matrix/gene_count_matrix.csv"))){
    cat("'gene_count_matrix.csv' is not exist!\n")
    stop("'gene_count_matrix.csv' not exist ERROR")
  }
  if(!dir.exists(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/"))){
    dir.create(paste0(path.prefix, "RNAseq_results/Reads_Count_Matrix_analysis/"))
  }
  DESeq2RawCountAnalysis(path.prefix, independent.variable,  control.group, experiment.group, DESeq2.padj, DESeq2.log2FC)
  edgeRRawCountAnalysis(path.prefix, independent.variable, control.group, experiment.group)
}


