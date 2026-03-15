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

