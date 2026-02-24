library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(cowplot)
library(CellChat)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(psych)
library(car)

options(future.globals.maxSize = 80 * 1024^3)

sc <- readRDS("02_Publicdata/ver6/tmp/seuratObj_add_label_latest_20251027.rds")
sc <- SetIdent(sc, value = "CellType_Class2")
sc_nodoublet <- subset(sc, idents = "Doublet", invert = TRUE)


# CellChat
## Part 0: Data input & processing and initialization of CellChat object
sc_OA4 <- subset(sc_nodoublet, subset = orig.ident == "OA_4")
sc_OA2 <- subset(sc_nodoublet, subset = orig.ident == "OA_2")

table(sc_OA2$CellType_Class2)
#   Adipocyte                     B cell                 CD4 T cell                 CD8 T cell             Dendritic cell           Endothelial cell 
#           2                         58                        249                        139                        134                         75 
# Granulocyte    Lining-layer fibroblast                 Macrophage                   Monocyte                    NK cell                   NKT cell 
#          10                        439                       1357                         95                         66                        111 
#    Pericyte                Plasma cell Sublining-layer fibroblast                       Treg                       VSMC 
#          11                         91                       1196                         23                          5 
table(sc_OA4$CellType_Class2)
#                  B cell                 CD4 T cell                 CD8 T cell             Dendritic cell           Endothelial cell                Granulocyte 
#                      34                         62                         38                        136                         21                         14 
# Lining-layer fibroblast                 Macrophage                   Monocyte                    NK cell                   NKT cell                   Pericyte 
#                    1149                       3666                        108                         17                         32                         19 
#             Plasma cell Sublining-layer fibroblast                       Treg                       VSMC 
#                      22                        608                          3                         10 

data.input_OA4 <- sc_OA4[["RNA"]]$data 
labels_OA4 <- Idents(sc_OA4)
meta_OA4 <- data.frame(labels = labels_OA4, row.names = names(labels_OA4))
cellchat_OA4 <- createCellChat(object = data.input_OA4, meta = meta_OA4, group.by = "labels")

data.input_OA2 <- sc_OA2[["RNA"]]$data 
labels_OA2 <- Idents(sc_OA2)
meta_OA2 <- data.frame(labels = labels_OA2, row.names = names(labels_OA2))
cellchat_OA2 <- createCellChat(object = data.input_OA2, meta = meta_OA2, group.by = "labels")

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_OA4@DB <- CellChatDB.use
cellchat_OA2@DB <- CellChatDB.use

cellchat_OA4 <- subsetData(cellchat_OA4)
cellchat_OA2 <- subsetData(cellchat_OA2)
cellchat_OA4 <- identifyOverExpressedGenes(cellchat_OA4)
cellchat_OA2 <- identifyOverExpressedGenes(cellchat_OA2)
cellchat_OA4 <- identifyOverExpressedInteractions(cellchat_OA4)
cellchat_OA2 <- identifyOverExpressedInteractions(cellchat_OA2)

cellchat_OA4 <- computeCommunProb(cellchat_OA4, type = "triMean")
cellchat_OA2 <- computeCommunProb(cellchat_OA2, type = "triMean")
cellchat_OA4 <- filterCommunication(cellchat_OA4, min.cells = 10)
cellchat_OA2 <- filterCommunication(cellchat_OA2, min.cells = 10)

cellchat_OA4 <- computeCommunProbPathway(cellchat_OA4)
cellchat_OA2 <- computeCommunProbPathway(cellchat_OA2)

cellchat_OA4 <- aggregateNet(cellchat_OA4)
cellchat_OA2 <- aggregateNet(cellchat_OA2)

cellchat_OA4 <- netAnalysis_computeCentrality(cellchat_OA4, slot.name = "netP")
cellchat_OA2 <- netAnalysis_computeCentrality(cellchat_OA2, slot.name = "netP")

group.new = levels(cellchat_OA2@idents)
cellchat_OA4 <- liftCellChat(cellchat_OA4, group.new)

object.list <- list(OA1 = cellchat_OA2, OA2 = cellchat_OA4)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# saveRDS(cellchat, file = "02_Publicdata/ver6/CellChat/ver2_overview/cellchat_merged.rds")


## Part I: Identify altered interactions and cell populations
gg1 <- netVisual_heatmap(cellchat, font.size = 6, font.size.title = 6)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", font.size = 6, font.size.title = 6)
png("99_Fig/ver6/CellChat/ver2_overview/1_interaction_heatmap.png", width = 7, height = 4, units = "in", res = 300)
gg1 + gg2
dev.off()
pdf("99_Fig/ver6/CellChat/ver2_overview/1_interaction_heatmap.pdf", width = 7, height = 4)
gg1 + gg2
dev.off()


# Part II: Identify altered signaling with distinct network architecture and interaction strength
# cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")
# netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

# cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural")
# netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
# netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Lining-layer fibroblast", targets.use = "Endothelial cell", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Lining-layer fibroblast", targets.use = "Endothelial cell", stacked = F, do.stat = TRUE)
png("99_Fig/ver6/CellChat/ver2_overview/2_signaling_pathway_comparison_lining_to_end.png", width = 9, height = 12, units = "in", res = 300)
gg1 + gg2
dev.off()
pdf("99_Fig/ver6/CellChat/ver2_overview/2_signaling_pathway_comparison_lining_to_end.pdf", width = 9, height = 12)
gg1 + gg2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Lining-layer fibroblast", targets.use = "Pericyte", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Lining-layer fibroblast", targets.use = "Pericyte", stacked = F, do.stat = TRUE)
png("02_Publicdata/ver6/CellChat/ver2_overview/2_signaling_pathway_comparison_lining_to_peri.png", width = 9, height = 12, units = "in", res = 300)
gg1 + gg2
dev.off()
pdf("02_Publicdata/ver6/CellChat/ver2_overview/2_signaling_pathway_comparison_lining_to_peri.pdf", width = 9, height = 12)
gg1 + gg2
dev.off()


# Part III: Identify the up-gulated and down-regulated signaling ligand-receptor pairs
levels(cellchat@meta$labels)
# [1] "Sublining-layer fibroblast" "Lining-layer fibroblast"    "Endothelial cell"           "CD8 T cell"                 "Macrophage"                 "Monocyte"                  
# [7] "Dendritic cell"             "VSMC"                       "CD4 T cell"                 "Pericyte"                   "NKT cell"                   "NK cell"                   
# [13] "Adipocyte"                  "B cell"                     "Plasma cell"                "Granulocyte"                "Treg" 

pos.dataset = "OA2"
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, idents.use = "Adipocyte", invert = TRUE, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "OA2",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "OA1",ligand.logFC = -0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

df_up <- findEnrichedSignaling(object.list[[2]], features = gene.up, pattern ="outgoing")
write_csv(df_up_down, "9.Publicdata/ver6/CellChat/ver2_overview/3_up_in_OA4_gene_info.csv")
df_down <- findEnrichedSignaling(object.list[[2]], features = gene.down, pattern ="outgoing")
write_csv(df_up_down, "9.Publicdata/ver6/CellChat/ver2_overview/3_down_in_OA4_gene_info.csv")

png("9.Publicdata/ver6/CellChat/ver2_overview/3_net_fibro_to_vessel_DEA.png", width = 30, height = 20, units = "in", res = 300)
circos.clear()
par(fig = c(0, 0.5, 0, 1), mar = c(0, 0, 0, 0), new = FALSE)
plot.new()
netVisual_chord_gene(object.list[[2]], sources.use = c("Lining-layer fibroblast","Sublining-layer fibroblast"), targets.use = c("Endothelial cell","Pericyte"), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
circos.clear()
par(fig = c(0.5, 1, 0, 1), mar = c(1, 1, 1, 1), new = TRUE)
plot.new()
netVisual_chord_gene(object.list[[1]], sources.use = c("Lining-layer fibroblast","Sublining-layer fibroblast"), targets.use = c("Endothelial cell","Pericyte"), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

pdf("9.Publicdata/ver6/CellChat/ver2_overview/3_net_fibro_to_vessel_DEA.pdf", width = 30, height = 20)
circos.clear()
par(fig = c(0, 0.5, 0, 1), mar = c(0, 0, 0, 0), new = FALSE)
plot.new()
netVisual_chord_gene(object.list[[2]], sources.use = c("Lining-layer fibroblast","Sublining-layer fibroblast"), targets.use = c("Endothelial cell","Pericyte"), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
circos.clear()
par(fig = c(0.5, 1, 0, 1), mar = c(1, 1, 1, 1), new = TRUE)
plot.new()
netVisual_chord_gene(object.list[[1]], sources.use = c("Lining-layer fibroblast","Sublining-layer fibroblast"), targets.use = c("Endothelial cell","Pericyte"), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()


sessionInfo()
# R version 4.3.3 (2024-02-29 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8    LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                    LC_TIME=Japanese_Japan.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] car_3.1-3             carData_3.0-5         psych_2.5.3           pheatmap_1.0.12       circlize_0.4.16       ComplexHeatmap_2.18.0 CellChat_2.1.1        Biobase_2.62.0        BiocGenerics_0.48.1   igraph_2.1.4          cowplot_1.1.3         patchwork_1.3.0      
# [13] lubridate_1.9.4       forcats_1.0.0         stringr_1.5.1         purrr_1.0.4           readr_2.1.5           tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1         tidyverse_2.0.0       dplyr_1.1.4           Seurat_5.1.0          SeuratObject_5.0.2   
# [25] sp_2.2-0             
# 
# loaded via a namespace (and not attached):
#  [1] RcppAnnoy_0.0.22       splines_4.3.3          later_1.4.1            polyclip_1.10-7        ggnetwork_0.5.13       fastDummies_1.7.5      lifecycle_1.0.4        rstatix_0.7.2          doParallel_1.0.17      globals_0.16.3         lattice_0.22-5        
# [12] MASS_7.3-60.0.1        backports_1.5.0        magrittr_2.0.3         plotly_4.10.4          sass_0.4.9             jquerylib_0.1.4        httpuv_1.6.15          NMF_0.28               sctransform_0.4.1      spam_2.11-1            spatstat.sparse_3.1-0 
# [23] reticulate_1.42.0      pbapply_1.7-2          RColorBrewer_1.1-3     abind_1.4-8            Rtsne_0.17             presto_1.0.0           IRanges_2.36.0         S4Vectors_0.40.2       ggrepel_0.9.6          irlba_2.3.5.1          listenv_0.9.1         
# [34] spatstat.utils_3.1-3   goftest_1.2-3          RSpectra_0.16-2        spatstat.random_3.3-3  fitdistrplus_1.2-2     parallelly_1.43.0      svglite_2.1.3          leiden_0.4.3.1         codetools_0.2-19       tidyselect_1.2.1       shape_1.4.6.1         
# [45] farver_2.1.2           matrixStats_1.5.0      stats4_4.3.3           spatstat.explore_3.4-2 jsonlite_2.0.0         GetoptLong_1.0.5       BiocNeighbors_1.20.2   Formula_1.2-5          progressr_0.15.1       ggridges_0.5.6         ggalluvial_0.12.5     
# [56] survival_3.5-8         iterators_1.0.14       systemfonts_1.2.1      foreach_1.5.2          tools_4.3.3            ragg_1.3.3             sna_2.8                ica_1.0-3              Rcpp_1.0.14            glue_1.8.0             mnormt_2.1.1          
# [67] gridExtra_2.3          withr_3.0.2            BiocManager_1.30.25    fastmap_1.2.0          digest_0.6.37          timechange_0.3.0       R6_2.6.1               mime_0.13              textshaping_1.0.0      colorspace_2.1-1       scattermore_1.2       
# [78] tensor_1.5             spatstat.data_3.1-6    generics_0.1.3         data.table_1.17.0      FNN_1.1.4.1            httr_1.4.7             htmlwidgets_1.6.4      uwot_0.2.3             pkgconfig_2.0.3        gtable_0.3.6           registry_0.5-1        
# [89] lmtest_0.9-40          htmltools_0.5.8.1      dotCall64_1.2          clue_0.3-66            scales_1.3.0           png_0.1-8              spatstat.univar_3.1-2  rstudioapi_0.17.1      tzdb_0.5.0             reshape2_1.4.4         rjson_0.2.23          
# [100] coda_0.19-4.1          statnet.common_4.11.0  nlme_3.1-164           zoo_1.8-13             cachem_1.1.0           GlobalOptions_0.1.2    KernSmooth_2.23-22     parallel_4.3.3         miniUI_0.1.1.1         pillar_1.10.1          vctrs_0.6.5           
# [111] RANN_2.6.2             promises_1.3.2         ggpubr_0.6.0           xtable_1.8-4           cluster_2.1.6          cli_3.6.4              compiler_4.3.3         rlang_1.1.5            crayon_1.5.3           rngtools_1.5.2         future.apply_1.11.3   
# [122] ggsignif_0.6.4         labeling_0.4.3         plyr_1.8.9             stringi_1.8.7          network_1.19.0         viridisLite_0.4.2      deldir_2.0-4           gridBase_0.4-7         BiocParallel_1.36.0    munsell_0.5.1          lazyeval_0.2.2        
# [133] spatstat.geom_3.3-6    Matrix_1.6-5           RcppHNSW_0.6.0         hms_1.1.3              future_1.34.0          shiny_1.10.0           ROCR_1.0-11            broom_1.0.7            bslib_0.9.0

