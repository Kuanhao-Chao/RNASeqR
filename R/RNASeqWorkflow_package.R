#' RNASeqWorkflow-package
#'
#' @name RNASeqWorkflow
#' @importFrom clusterProfiler bitr enrichGO groupGO dotplot emapplot cnetplot goplot enrichKEGG gseGO gseaplot gseKEGG
#' @importFrom pathview pathview
#' @import DOSE
#' @import ggplot2
#' @import graphics
#' @import utils
#' @import grDevices
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @importFrom tools file_ext
#' @importFrom reticulate py_available py_config
#' @importFrom systemPipeRdata genWorkenvir
#' @importFrom systemPipeR systemArgs seeFastq infile1 seeFastqPlot
#' @importFrom ShortRead readFastq narrow width writeFastq
#' @importFrom Biostrings quality
#' @importFrom ballgown ballgown texpr indexes stattest gexpr
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results counts plotMA
#' @importFrom edgeR DGEList calcNormFactors estimateCommonDisp estimateTagwiseDisp exactTest cpm plotMDS.DGEList plotMeanVar plotBCV
#' @importFrom stats t.test quantile cor
#' @importFrom stringr str_extract
#' @importFrom gridExtra grid.table
#' @importFrom rafalib mypar
#' @importFrom reshape2 melt
#' @importFrom FactoMineR PCA
#' @importFrom factoextra get_eigenvalue fviz_eig fviz_pca_ind
#' @importFrom corrplot corrplot
#' @importFrom PerformanceAnalytics chart.Correlation
#' @importFrom pheatmap pheatmap
NULL

