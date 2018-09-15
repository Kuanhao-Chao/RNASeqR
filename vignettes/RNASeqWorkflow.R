## ----style, echo=FALSE, results="asis", message=FALSE----------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE)

## ----install, eval=FALSE, warning=FALSE------------------------------------
#  source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script
#  biocLite("RNASeqWorkflow") # Installs RNASeqWorkflow
#  biocLite("RNASeqWorkflowData") # Installs RNASeqWorkflowData

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

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqQualityAssessment_CMD(exp)

## ---- warning=FALSE--------------------------------------------------------
RNASeqQualityAssessment(exp@path.prefix,
                        exp@input.path.prefix,
                        exp@sample.pattern)

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqQualityTrimming_CMD(exp)

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqQualityTrimming(exp@path.prefix,
#                        exp@sample.pattern)

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqReadProcess_CMD(exp)

## ---- warning=FALSE--------------------------------------------------------
python.variable <- exp@python.variable
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
toLatex(sessionInfo())

