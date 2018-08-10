#' RNASeqWorkflow-package
#'
#' @name RNASeqWorkflow
#' @importFrom tools file_ext
#' @importFrom reticulate py_available py_config
#' @importFrom ballgown ballgown texpr subset stattest geneNames geneIDs transcriptNames indexes structure
#' @importFrom dplyr arrange count
#' @importFrom genefilter rowVars
#' @importFrom gridExtra grid.table
#' @importFrom rafalib mypar shist
#' @importFrom reshape2 melt
#' @importFrom FactoMineR PCA
#' @importFrom factoextra get_eigenvalue fviz_eig fviz_pca_ind
#' @importFrom corrplot corrplot
#' @importFrom PerformanceAnalytics chart.Correlation
#' @importFrom clusterProfiler bitr enrichGO groupGO dotplot emapplot cnetplot goplot enrichKEGG gseGO gseaplot gseKEGG
#' @importFrom pathview pathview
#' @importFrom refGenome tableSeqids extractSeqids tableFeatures ensemblGenome read.gtf
#' @importFrom genefilter rowVars
#' @importFrom rtracklayer import
#' @importFrom Rqc rqcQA rqcReport
#' @importFrom systemPipeRdata genWorkenvir
#' @importFrom systemPipeR systemArgs seeFastq infile1 seeFastqPlot
#' @importFrom QuasR preprocessReads
#' @import DESeq2
#' @import edgeR
#' @importFrom stringr str_extract
#' @import ggplot2
#' @import pheatmap pheatmap
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @importFrom graphics abline barplot hist legend mtext par
#' @importFrom grDevices colorRampPalette dev.off pdf png rgb
#' @importFrom methods new
#' @importFrom stats cor heatmap model.matrix
#' @importFrom utils data download.file head read.csv read.delim write.csv write.table
#' @importFrom clusterProfiler gseGO gseaplot
#' @importFrom rtracklayer mcols
#'
NULL

