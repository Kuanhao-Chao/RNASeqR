## ----style, echo=FALSE, results="asis", message=FALSE----------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE)

## ----install, eval=FALSE, warning=FALSE------------------------------------
#  source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script
#  biocLite("RNASeqWorkflow") # Installs RNASeqWorkflow
#  biocLite("RNASeqWorkflowData") # Installs RNASeqWorkflowData

## ----fig.width=10,  echo=FALSE---------------------------------------------
library(png)
library(grid)
img <- readPNG("figure/whole_file_structure.png")
grid.raster(img, just = "center")

## ----fig.width=10, echo=FALSE----------------------------------------------
img <- readPNG("figure/input_files_structure.png")
 grid.raster(img, just = "center")

## ---- warning=FALSE--------------------------------------------------------
library(RNASeqWorkflow)
library(RNASeqWorkflowData)

## --------------------------------------------------------------------------
input_files.path <- system.file("extdata/", package = "RNASeqWorkflowData")
rnaseq_result.path <- "/tmp/yeast_example/"
dir.create(rnaseq_result.path)
list.files(input_files.path, recursive = TRUE)

## ---- warning=FALSE--------------------------------------------------------
exp <- RNASeqWorkflowParam(path.prefix = rnaseq_result.path, 
                           input.path.prefix = input_files.path, 
                           genome.name = "Saccharomyces_cerevisiae_XV_Ensembl", 
                           sample.pattern = "SRR[0-9]*_XV",
                           independent.variable = "state", 
                           case.group = "60mins_ID20_amphotericin_B", 
                           control.group = "60mins_ID20_control")

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqEnvironmentSet_CMD(exp)

## ---- warning=FALSE--------------------------------------------------------
RNASeqEnvironmentSet(exp@path.prefix,
                     exp@input.path.prefix,
                     exp@genome.name,
                     exp@sample.pattern,
                     exp@indices.optional,
                     exp@os.type)

## ----fig.width=20, echo=FALSE----------------------------------------------
img <- readPNG("figure/fastqReport.png")
 grid.raster(img, just = "center")

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqQualityAssessment_CMD(exp)

## ---- warning=FALSE, eval=FALSE--------------------------------------------
#  RNASeqQualityAssessment(exp@path.prefix,
#                          exp@input.path.prefix,
#                          exp@sample.pattern)

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqQualityTrimming_CMD(exp)

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqQualityTrimming(exp@path.prefix,
#                        exp@sample.pattern)

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqReadProcess_CMD(exp)

## ---- warning=FALSE--------------------------------------------------------
python.variable <- "@"(exp, python.variable)
python.variable.answer <- python.variable$check.answer
python.variable.version <- python.variable$python.version
RNASeqReadProcess(exp@path.prefix,
                  exp@input.path.prefix,
                  exp@genome.name,
                  exp@sample.pattern,
                  python.variable.answer,
                  python.variable.version,
                  exp@python.2to3,
                  num.parallel.threads = 10,
                  exp@indices.optional)

## ----fig.width=10, echo=FALSE----------------------------------------------
img <- readPNG("figure/Alignment_report.png")
 grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/Frequency/Frequency_Plot_normalized_count_ggplot2.png")
 grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/Frequency/Frequency_Plot_log_normalized_count_ggplot2.png")
 grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/Distribution/Box_Plot_ggplot2.png")
 grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/Distribution/Violin_Plot_ggplot2.png")
 grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/PCA/Dimension_PCA_Plot_factoextra.png")
 grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/PCA/PCA_Plot_factoextra.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/PCA/PCA_Plot_ggplot2.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/Correlation/Correlation_Heat_Plot_ggplot2.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/Correlation/Correlation_Dot_Plot_corrplot.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/preDE/Correlation/Correlation_Bar_Plot_PerformanceAnalytics.png")
grid.raster(img, just = "center")

## ----fig.width=8, fig.height=8, echo=FALSE---------------------------------
img <- readPNG("figure/DE/Volcano_Plot_graphics.png")
 grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/DE/PCA/Dimension_PCA_Plot_factoextra.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/DE/PCA/PCA_Plot_factoextra.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/DE/PCA/PCA_Plot_ggplot2.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/DE/Heatmap_Plot_pheatmap.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/transcript_related/Distribution_Transcript_Count_per_Gene_Plot.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/transcript_related/Distribution_Transcript_Length_Plot.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/ballgown_MA_Plot_ggplot2.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/TPM_MA_Plot_ggplot2.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/DESeq2_Dispersion_Plot_DESeq2.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/DESeq2_MA_Plot_DESeq2.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/edgeR_MeanVar_Plot_edgeR.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/edgeR_BCV_Plot_edgeR.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/edgeR_MDS_Plot_edgeR.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/edgeR_Smear_Plot_edgeR.png")
grid.raster(img, just = "center")

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqDifferentialAnalysis_CMD(exp)

## --------------------------------------------------------------------------
RNASeqDifferentialAnalysis(exp@path.prefix,
                           exp@genome.name,
                           exp@sample.pattern,
                           exp@independent.variable,
                           exp@case.group,
                           exp@control.group,
                           ballgown.pval = 0.05,
                           ballgown.log2FC = 1,
                           TPM.pval = 0.05,
                           TPM.log2FC = 1,
                           DESeq2.pval = 0.1,
                           DESeq2.log2FC = 1,
                           edgeR.pval = 0.05,
                           edgeR.log2FC = 1)

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/GO_analysis/GO_CC_Classification_Bar_Plot_clusterProfiler.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/GO_analysis/GO_MF_Overrepresentation_Bar_Plot_clusterProfiler.png")
grid.raster(img, just = "center")

## ----fig.width=8, height=8, echo=FALSE-------------------------------------
img <- readPNG("figure/GO_analysis/GO_MF_Overrepresentation_Dot_Plot_clusterProfiler.png")
grid.raster(img, just = "center")

## ---- eval = FALSE---------------------------------------------------------
#  RNASeqGoKegg_CMD(exp,
#                   OrgDb.species = "org.Sc.sgd.db",
#                   go.level = 3,
#                   input.TYPE.ID = "GENENAME",
#                   KEGG.organism = "sce")

## ---- eval = FALSE---------------------------------------------------------
#  RNASeqGoKegg(exp@path.prefix,
#               exp@independent.variable,
#               OrgDb.species = "org.Sc.sgd.db",
#               go.level = 3,
#               input.TYPE.ID = "GENENAME",
#               KEGG.organism = "sce")

## --------------------------------------------------------------------------
sessionInfo()

