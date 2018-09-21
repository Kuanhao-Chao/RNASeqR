# # TPM normalization
# TPMNormalizationAnalysis <- function(path.prefix,
#                                      genome.name,
#                                      sample.pattern,
#                                      independent.variable,
#                                      case.group,
#                                      control.group,
#                                      TPM.pval,
#                                      TPM.log2FC) {
#   message("\n\u2618\u2618 TPM analysis ...\n")
#   if(!dir.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/"))){
#     dir.create(paste0(path.prefix, "RNASeq_results/TPM_analysis/"))
#   }
#   if(!dir.exists(paste0(path.prefix,
#                         "RNASeq_results/TPM_analysis/normalized_&_statistic"))){
#     dir.create(paste0(path.prefix,
#                       "RNASeq_results/TPM_analysis/normalized_&_statistic"))
#   }
#   ###############################################
#   ## Creating "TPM_normalized_result.csv" ##
#   ##############################################
#   case.FPKM <- read.csv(paste0(path.prefix,
#                                "RNASeq_results/ballgown_analysis/",
#                                "normalized_&_statistic/FPKM_case.csv"))
#   control.FPKM <- read.csv(paste0(path.prefix,
#                                   "RNASeq_results/ballgown_analysis/",
#                                   "normalized_&_statistic/FPKM_control.csv"))
#   statistic.FPKM <- read.csv(paste0(path.prefix,
#                                     "RNASeq_results/ballgown_analysis/",
#                                     "normalized_&_statistic/statistic.csv"))
#   gene.name <- read.csv(paste0(path.prefix,
#                                "RNASeq_results/ballgown_analysis/",
#                                "ballgown_R_object/gene_name.csv"))
#
#
#   case.TPM <- t(t(case.FPKM) / colSums(case.FPKM)) * 10**6
#   control.TPM <- t(t(control.FPKM) / colSums(control.FPKM)) * 10**6
#   gene.id.data.frame <-
#     data.frame(read.csv(paste0(path.prefix,
#                                "RNASeq_results/ballgown_analysis/",
#                                "ballgown_R_object/gene_name.csv")))
#   p.value <- unlist(lapply(seq_len(nrow(case.TPM)),
#                            function(x) {
#                              stats::t.test(case.TPM[x,],control.TPM[x,])$p.value
#                              }))
#   fold.change <- unlist(lapply(seq_len(nrow(case.TPM)),
#                                function(x) {
#                                  mean(unlist(control.TPM[x,]) + 1) /
#                                    mean(unlist(case.TPM[x,]) + 1)
#                                  }))
#   statistic.T.test <- data.frame("pval" = p.value,
#                                  "fc" = fold.change,
#                                  "log2FC" = log2(fold.change))
#   total.data.frame <- cbind(gene.id.data.frame,
#                             case.TPM,
#                             control.TPM)
#   total.data.frame[paste0(case.group, ".average")] <- rowMeans(case.TPM)
#   total.data.frame[paste0(control.group, ".average")] <- rowMeans(control.TPM)
#   total.data.frame[paste0(case.group, ".", control.group, ".average")]<-
#     rowMeans(cbind(case.TPM, control.TPM))
#   TPM.result <- cbind(total.data.frame, statistic.T.test)
#   TPM.result <- rbind(TPM.result[TPM.result$gene.name != ".",],
#                       TPM.result[TPM.result$gene.name == ".",])
#   # Filter out gene (p-value = Na) and
#   #                 (log2FC = Inf or log2FC = Na or log2FC = -Inf)
#   TPM.result <- TPM.result[!is.na(TPM.result$pval) &
#                              (TPM.result$log2FC != Inf) &
#                              !is.na(TPM.result$log2FC) &
#                              (TPM.result$log2FC != -Inf), ]
#   case.group.size <- length(case.FPKM)
#   control.group.size <- length(control.FPKM)
#   # write out csv files
#   write.csv(TPM.result[,c(2:(case.group.size+1))],
#             file = paste0(path.prefix,
#                           "RNASeq_results/TPM_analysis/",
#                           "normalized_&_statistic/TPM_case.csv"),
#             row.names=FALSE)
#   write.csv(TPM.result[,c((2+case.group.size):
#                             (case.group.size+control.group.size+1))],
#             file = paste0(path.prefix, "RNASeq_results/TPM_analysis/",
#                           "normalized_&_statistic/TPM_control.csv"),
#             row.names=FALSE)
#   write.csv(TPM.result[,c((2+3+case.group.size+control.group.size):
#                             (length(TPM.result)))],
#             file = paste0(path.prefix, "RNASeq_results/TPM_analysis/",
#                           "normalized_&_statistic/statistic.csv"),
#             row.names=FALSE)
#   write.csv(TPM.result,
#             file = paste0(path.prefix, "RNASeq_results/TPM_analysis/",
#                           "TPM_normalized_result.csv"),
#             row.names=FALSE)
#
#   message(paste0("     \u25CF Selecting differential expressed genes() ",
#                  "==> p-value : ", TPM.pval,
#                  "  log2(Fold Change) : ", TPM.log2FC, " ...\n"))
#   TPM.result.DE <- TPM.result[((TPM.result$log2FC>TPM.log2FC) |
#                                  (TPM.result$log2FC<(-TPM.log2FC))) &
#                                 (TPM.result$pval<TPM.pval), ]
#   message(paste0("          \u25CF Total '", length(row.names(TPM.result.DE)),
#                  "' DEG have been found !!!\n"))
#   write.csv(TPM.result.DE,
#             file = paste0(path.prefix, "RNASeq_results/TPM_analysis/",
#                           "TPM_normalized_DE_result.csv"),
#             row.names=FALSE)
#
#   # Check TPM.result.DE before visulization!!
#
#   ###########################
#   ## TPM&Ttest visulization ##
#   ###########################
#   if(file.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/",
#                         "TPM_normalized_result.csv")) &&
#      file.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/",
#                         "normalized_&_statistic/TPM_case.csv")) &&
#      file.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/",
#                         "normalized_&_statistic/TPM_control.csv")) &&
#      file.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/",
#                         "normalized_&_statistic/statistic.csv"))){
#     # Transcript Related
#     if(!dir.exists(paste0(path.prefix, "RNASeq_results/TPM_analysis/images/"))){
#       dir.create(paste0(path.prefix, "RNASeq_results/TPM_analysis/images/"))
#     }
#     if (nrow(TPM.result) > 1) {
#       ###############
#       #### PreDE ####
#       ###############
#       if(!dir.exists(paste0(path.prefix, "RNASeq_results/",
#                             "TPM_analysis/images/preDE/"))){
#         dir.create(paste0(path.prefix, "RNASeq_results/",
#                           "TPM_analysis/images/preDE/"))
#       }
#       # Frequency
#       FrequencyPlot("TPM_analysis",
#                     "TPM",
#                     path.prefix,
#                     independent.variable,
#                     case.group,
#                     control.group)
#       # Bax and Violin
#       BoxViolinPlot("TPM_analysis",
#                     "TPM",
#                     path.prefix,
#                     independent.variable,
#                     case.group,
#                     control.group)
#       # PCA
#       PCAPlot("TPM_analysis",
#               "TPM",
#               path.prefix,
#               independent.variable,
#               case.group,
#               control.group)
#       #Correlation
#       CorrelationPlot("TPM_analysis",
#                       "TPM",
#                       path.prefix,
#                       independent.variable,
#                       case.group,
#                       control.group)
#       ############
#       #### DE ####
#       ############
#       if(!dir.exists(paste0(path.prefix,
#                             "RNASeq_results/TPM_analysis/images/DE/"))){
#         dir.create(paste0(path.prefix,
#                           "RNASeq_results/TPM_analysis/images/DE/"))
#       }
#       # Volcano
#       VolcanoPlot("TPM_analysis",
#                   "TPM",
#                   path.prefix,
#                   independent.variable,
#                   case.group,
#                   control.group,
#                   TPM.pval,
#                   TPM.log2FC)
#       # MA
#       MAPlot("TPM_analysis",
#              "TPM",
#              path.prefix,
#              independent.variable,
#              case.group,
#              control.group,
#              TPM.pval)
#       if (nrow(TPM.result.DE) > 1) {
#         # DE PCA plot
#         DEPCAPlot("TPM_analysis",
#                   "TPM",
#                   path.prefix,
#                   independent.variable,
#                   case.group,
#                   control.group)
#         # Heatmap
#         DEHeatmap("TPM_analysis",
#                   "TPM",
#                   path.prefix,
#                   independent.variable,
#                   case.group,
#                   control.group)
#       } else {
#         cat ("(\u26A0) Less than one differential expressed gene term found !!",
#              " Skip DE_PCA and DE_Heatmap visualization !!! \n\n")
#       }
#     } else {
#       cat ("(\u26A0) Less than one gene terms found !!! ",
#            "Skip visualization step !!! \n\n")
#     }
#   } else {
#     message("(\u2718) necessary file is missing!! Something ERROR happend ",
#             "during edgeR analysis!! Skip visualization!!\n\n")
#   }
# }
