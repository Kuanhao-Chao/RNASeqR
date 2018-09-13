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
                     exp@indices.optional, exp@os.type)

## ---- eval=FALSE-----------------------------------------------------------
#  RNASeqQualityAssessment_CMD(exp)

