#' @title RNASeqGoKegg_CMD
#'
#' @description
#'   Run Gene Ontology(GO) and Kyoto Encyclopedia of Genes and Genomes(KEGG)
#'   analysis in background. \cr
#'   This function do Gene Ontology(GO) and Kyoto Encyclopedia of Genes and
#'   Genomes(KEGG) analysis : \cr
#'   \enumerate{
#'     \item Gene Ontology(GO) :\cr
#'      \enumerate{
#'        \item Do GO function classification analysis. \cr
#'        \item Do GO function enrichment analysis. \cr
#'        \item Visualization : bar plot, dot plot etc. \cr
#'      }
#'     \item Kyoto Encyclopedia of Genes and Genomes(KEGG) :\cr
#'      \enumerate{
#'        \item Do KEGG pathway enrichment analysis \cr
#'        \item Pathway visulization with \code{pathview} package. KEGG webpage
#'         pathway url will also be created \cr
#'      }
#'   }
#'   If you want to do GO functional analysis and KEGG pathway analysis for the
#'   following RNA-Seq workflow in R shell,
#'   please see \code{RNASeqGoKegg()} function.
#'
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param OrgDb.species the genome wide annotation packages of species on
#'   Bioconductor. Currently, there are 19 supported genome wide annotation
#'   packages of species.
#' @param go.level the depth of acyclic graph in GO analysis
#' @param input.TYPE.ID The gene name type in OrgDb.species annotation packahge.
#' @param KEGG.organism the species that are supported for KEGG analysis.
#'   Currently, there are more than 5000 supported species genome.
#'   Check the valid species terms on
#'   https://www.genome.jp/kegg/catalog/org_list.html
#' @param run Default value is \code{TRUE}. If \code{TRUE},
#'   'Rscript/Environment_Set.R' will be created and executed.
#'   The output log will be stored in 'Rscript_out/Environment_Set.Rout'.
#'   If \code{False}, 'Rscript/Environment_Set.R' will be created
#'   without executed.
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE},
#'   the result of checking \code{RNASeqRParam} will be reported in
#'   'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqGoKegg_CMD(RNASeqRParam  = yeast,
#'                  OrgDb.species = "org.Sc.sgd.db",
#'                  go.level = 3,
#'                  input.TYPE.ID = "GENENAME",
#'                  KEGG.organism = "sce")
#' }
RNASeqGoKegg_CMD <- function(RNASeqRParam,
                             OrgDb.species,
                             go.level = 3,
                             input.TYPE.ID,
                             KEGG.organism,
                             run = TRUE,
                             check.s4.print = TRUE) {
  # check input param
  CheckS4Object(RNASeqRParam, check.s4.print)
  CheckOperatingSystem(FALSE)
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  INSIDE.path.prefix <- "@"(RNASeqRParam, path.prefix)
  saveRDS(RNASeqRParam,
          file = paste0(INSIDE.path.prefix,
                        "gene_data/RNASeqRParam.rds"))
  fileConn<-file(paste0(path.prefix, "Rscript/GO_KEGG_Analysis.R"))
  first <- "library(RNASeqR)"
  second <- paste0("RNASeqGoKegg(RNASeqRParam = 'INSIDE'",
                   ", which.trigger = 'INSIDE'",
                   ", INSIDE.path.prefix = '", INSIDE.path.prefix,
                   "', OrgDb.species = '", OrgDb.species,
                   "', go.level = ", go.level,
                   ", input.TYPE.ID = '", input.TYPE.ID,
                   "', KEGG.organism = '",KEGG.organism, "')")
  writeLines(c(first, second), fileConn)
  close(fileConn)
  message(paste0("\u2605 '", path.prefix,
                 "Rscript/GO_KEGG_Analysis.R' has been created.\n"))
  if (run) {
    R.home.lib <- R.home()
    R.home.bin <- gsub("/lib/R", "/bin/R", R.home.lib)
    system2(command = 'nohup',
            args = paste0(R.home.bin, " CMD BATCH ",
                          path.prefix,
                          "Rscript/GO_KEGG_Analysis.R ",
                          path.prefix, "Rscript_out/GO_KEGG_Analysis.Rout"),
            stdout = "", wait = FALSE)
    message("\u2605 Tools are installing in the background. ",
            "Check current progress in '", path.prefix,
            "Rscript_out/GO_KEGG_Analysis.Rout'\n\n")
  }
}

#' @title RNASeqGoKegg
#'
#' @description
#'   Run Gene Ontology(GO) and Kyoto Encyclopedia of Genes and Genomes(KEGG)
#'   analysis in R shell. \cr
#'   This function do Gene Ontology(GO) and Kyoto Encyclopedia of Genes and
#'   Genomes(KEGG) analysis : \cr
#'   \enumerate{
#'     \item Gene Ontology(GO) :\cr
#'      \enumerate{
#'        \item Do GO function classification analysis. \cr
#'        \item Do GO function enrichment analysis. \cr
#'        \item Visualization : bar plot, dot plot etc. \cr
#'      }
#'     \item Kyoto Encyclopedia of Genes and Genomes(KEGG) :\cr
#'      \enumerate{
#'        \item Do KEGG pathway enrichment analysis \cr
#'        \item Pathway visulization with \code{pathview} package. KEGG webpage
#'         pathway url will also be created \cr
#'      }
#'   }
#'   If you want to do GO functional analysis and KEGG pathway analysis
#'   for the following RNA-Seq workflow in background,
#'    please see \code{RNASeqGoKegg_CMD()} function.
#' @param RNASeqRParam S4 object instance of experiment-related
#'   parameters
#' @param which.trigger Default value is \code{OUTSIDE}. User should not change
#'   this value.
#' @param INSIDE.path.prefix Default value is \code{NA}. User should not change
#'   this value.
#' @param OrgDb.species the genome wide annotation packages of species
#'   on Bioconductor. Currently, there are 19 supported genome wide
#'   annotation packages of species.
#' @param go.level the depth of acyclic graph in GO analysis
#' @param input.TYPE.ID The gene name type in OrgDb.species annotation packahge.
#' @param KEGG.organism the species that are supported for KEGG analysis.
#'   Currently, there are more than 5000 supported species genome.
#'   Check the valid species terms on
#'   https://www.genome.jp/kegg/catalog/org_list.html
#' @param check.s4.print Default \code{TRUE}. If \code{TRUE},
#'   the result of checking \code{RNASeqRParam} will be reported in
#'   'Rscript_out/Environment_Set.Rout'. If \code{FALSE}, the result of checking
#'   \code{RNASeqRParam} will not be in
#'   'Rscript_out/Environment_Set.Rout'
#'
#' @return None
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(yeast)
#' \dontrun{
#' RNASeqGoKegg(RNASeqRParam  = yeast,
#'              OrgDb.species = "org.Sc.sgd.db",
#'              go.level = 3,
#'              input.TYPE.ID = "GENENAME",
#'              KEGG.organism = "sce")}
RNASeqGoKegg <- function(RNASeqRParam,
                         which.trigger      = "OUTSIDE",
                         INSIDE.path.prefix = NA,
                         OrgDb.species,
                         go.level = 3,
                         input.TYPE.ID,
                         KEGG.organism,
                         check.s4.print = TRUE) {
  CheckOperatingSystem(FALSE)
  # If `which.trigger` is OUTSIDE, then directory must be built
  # If `which.trigger` is INSIDE, then directory must not be
  #  built here(will created in CMD)
  if (isS4(RNASeqRParam) &
      which.trigger == "OUTSIDE" &
      is.na(INSIDE.path.prefix)) {
    # This is an external call!!
    # Check the S4 object(user input)
    CheckS4Object(RNASeqRParam, check.s4.print)
  } else if (RNASeqRParam == "INSIDE" &
             which.trigger == "INSIDE" &
             !is.na(INSIDE.path.prefix)) {
    # This is an internal call!!
    # Load the S4 object that saved in CMD process
    RNASeqRParam <- readRDS(paste0(INSIDE.path.prefix,
                                   "gene_data/RNASeqRParam.rds"))
  }
  path.prefix <- "@"(RNASeqRParam, path.prefix)
  independent.variable <- "@"(RNASeqRParam, independent.variable)
  PreRNASeqGoKegg()
  raw.read.avail <- RawReadCountAvailability(path.prefix)
  ballgown.bool <- dir.exists(paste0(path.prefix,
                                     "RNASeq_results/ballgown_analysis/"))
  DESeq2.bool <- dir.exists(paste0(path.prefix,
                                   "RNASeq_results/DESeq2_analysis/"))
  edgeR.bool <- dir.exists(paste0(path.prefix,
                                  "RNASeq_results/edgeR_analysis/"))
  which.analyses <- c()
  if (ballgown.bool) {
    which.analyses <- c(which.analyses, "ballgown_analysis")
  }
  if (raw.read.avail) {
    if (DESeq2.bool) {
      which.analyses <- c(which.analyses, "DESeq2_analysis")
    }
    if (edgeR.bool) {
      which.analyses <- c(which.analyses, "edgeR_analysis")
    }
  }
  message("\u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618  ",
          "Start 'Gene Ontology', 'Kyoto Encyclopedia of Genes and Genomes' ",
          "analyses  \u2618\u2618\u2618\u2618\u2618\u2618\u2618\u2618\n")
  for(which.analysis in which.analyses) {
    message("\u2618\u2618 ",
            strsplit(which.analysis, "_")[[1]][1] , " analysis ...\n")
    GOAnalysis(which.analysis,
               path.prefix,
               OrgDb.species,
               go.level,
               input.TYPE.ID)
    KEGGAnalysis(which.analysis,
                 path.prefix,
                 OrgDb.species,
                 input.TYPE.ID,
                 KEGG.organism)
  }
  PostRNASeqGoKegg()
}

PreRNASeqGoKegg <- function() {
  message("\u269C\u265C\u265C\u265C RNASeqGoKegg()' ",
          "environment pre-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqGoKegg() environment ERROR")
  }
  message("(\u2714) : RNASeqGoKegg() pre-check is valid\n\n")
}

PostRNASeqGoKegg <- function() {
  message("\u269C\u265C\u265C\u265C RNASeqGoKegg()' ",
          "environment post-check ...\n")
  validity <- TRUE
  if (!isTRUE(validity)) {
    stop("RNASeqGoKegg() post-check ERROR")
  }
  message("(\u2714) : RNASeqGoKegg() post-check is valid\n\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605 Success!! \u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\n")
  message("\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605",
          "\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\u2605\n")
}
