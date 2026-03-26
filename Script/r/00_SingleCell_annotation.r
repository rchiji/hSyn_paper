library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

options(future.globals.maxSize = 80 * 1024^3)


# Preprocessing and coarse cell type annotation **********
count1 <- Read10X(data.dir = "./02_Publicdata/scData/nonOA_1/")
count2 <- Read10X(data.dir = "./02_Publicdata/scData/nonOA_2/")
count3 <- Read10X(data.dir = "./02_Publicdata/scData/nonOA_3/")
count4 <- Read10X(data.dir = "./02_Publicdata/scData/nonOA_4/")
count5 <- Read10X(data.dir = "./02_Publicdata/scData/OA_1/")
count6 <- Read10X(data.dir = "./02_Publicdata/scData/OA_2/")
count7 <- Read10X(data.dir = "./02_Publicdata/scData/OA_3/")
count8 <- Read10X(data.dir = "./02_Publicdata/scData/OA_4/")
count9 <- Read10X(data.dir = "./02_Publicdata/scData/OA_5/")
count10 <- Read10X(data.dir = "./02_Publicdata/scData/OA_6/")
count11 <- Read10X(data.dir = "./02_Publicdata/scData/OA_7/")
count12 <- Read10X(data.dir = "./02_Publicdata/scData/OA_8/")
count13 <- Read10X(data.dir = "./02_Publicdata/scData/OA_9/")
count14 <- Read10X(data.dir = "./02_Publicdata/scData/OA_10/")

seuratObj1 <- CreateSeuratObject(counts = count1, project = "nonOA_1", min.cells = 3)
seuratObj2 <- CreateSeuratObject(counts = count2, project = "nonOA_2", min.cells = 3)
seuratObj3 <- CreateSeuratObject(counts = count3, project = "nonOA_3", min.cells = 3)
seuratObj4 <- CreateSeuratObject(counts = count4, project = "nonOA_4", min.cells = 3)
seuratObj5 <- CreateSeuratObject(counts = count5, project = "OA_1", min.cells = 3)
seuratObj6 <- CreateSeuratObject(counts = count6, project = "OA_2", min.cells = 3)
seuratObj7 <- CreateSeuratObject(counts = count7, project = "OA_3", min.cells = 3)
seuratObj8 <- CreateSeuratObject(counts = count8, project = "OA_4", min.cells = 3)
seuratObj9 <- CreateSeuratObject(counts = count9, project = "OA_5", min.cells = 3)
seuratObj10 <- CreateSeuratObject(counts = count10, project = "OA_6", min.cells = 3)
seuratObj11 <- CreateSeuratObject(counts = count11, project = "OA_7", min.cells = 3)
seuratObj12 <- CreateSeuratObject(counts = count12, project = "OA_8", min.cells = 3)
seuratObj13 <- CreateSeuratObject(counts = count13, project = "OA_9", min.cells = 3)
seuratObj14 <- CreateSeuratObject(counts = count14, project = "OA_10", min.cells = 3)

sc <- merge(seuratObj1,c(seuratObj2,seuratObj3,seuratObj4,seuratObj5,seuratObj6,seuratObj7,seuratObj8,seuratObj9,seuratObj10,seuratObj11,seuratObj12,seuratObj13,seuratObj14),
             add.cell.ids = c("nonOA_1","nonOA_2","nonOA_3","nonOA_4","OA_1","OA_2","OA_3","OA_4","OA_5","OA_6","OA_7","OA_8","OA_9","OA_10"))

table(sc@meta.data$orig.ident)
# -> nonOA_1 nonOA_2 nonOA_3 nonOA_4    OA_1   OA_10    OA_2    OA_3    OA_4    OA_5    OA_6    OA_7    OA_8    OA_9 
#       7611   14130   15363   12018    9042    1707    5300   10901    7613   10478   12805   10539    5100    3883

# saveRDS(sc, "02_Publicdata/ver6/tmp/seuratObj_after_merge_20251022.rds")


sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(sc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
sc <- subset(sc, subset = nFeature_RNA>500 & nFeature_RNA<6000 & percent.mt<20)

table(sc@meta.data$orig.ident)
# -> nonOA_1 nonOA_2 nonOA_3 nonOA_4    OA_1   OA_10    OA_2    OA_3    OA_4    OA_5    OA_6    OA_7    OA_8    OA_9 
#       6577   11741   11507   10100    7925    1619    4154   10138    6142    9265   10983    8754    5083    3840 


sc <- NormalizeData(sc)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sc <- JoinLayers(sc)

sc <- CellCycleScoring(sc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

sc[["RNA"]] <- split(sc[["RNA"]], f = sc$orig.ident)
sc <- SCTransform(sc, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

# saveRDS(sc, "02_Publicdata/ver6/tmp/seuratObj_after_SCT_20251022.rds")


sc <- RunPCA(sc)
sc
# -> An object of class Seurat 
# 62707 features across 107828 samples within 2 assays 
# Active assay: SCT (31224 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
#  1 dimensional reduction calculated: pca

sc <- IntegrateLayers(object = sc,
                             method = RPCAIntegration,
                             normalization.method = "SCT",
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = FALSE
                             )
sc <- JoinLayers(sc, assay = "RNA")
sc
# -> An object of class Seurat 
# 62707 features across 107828 samples within 2 assays 
# Active assay: SCT (31224 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, integrated.rpca

# saveRDS(sc, "02_Publicdata/ver6/tmp/seuratObj_after_integration_20251022.rds")


sc_dim10 <- RunUMAP(sc, dims = 1:10, reduction="integrated.rpca", reduction.name = "umap.rpca")
sc_dim20 <- RunUMAP(sc, dims = 1:20, reduction="integrated.rpca", reduction.name = "umap.rpca")
sc_dim30 <- RunUMAP(sc, dims = 1:30, reduction="integrated.rpca", reduction.name = "umap.rpca")

DimPlot(sc_dim10, group.by = "orig.ident", shuffle = TRUE) + 
    DimPlot(sc_dim20, group.by = "orig.ident", shuffle = TRUE) + 
    DimPlot(sc_dim30, group.by = "orig.ident", shuffle = TRUE)

sc_dim10 <- FindNeighbors(sc_dim10, reduction="integrated.rpca", dims = 1:10)
sc_dim20 <- FindNeighbors(sc_dim20, reduction="integrated.rpca", dims = 1:20)
sc_dim30 <- FindNeighbors(sc_dim30, reduction="integrated.rpca", dims = 1:30)

sc_dim10 <- FindClusters(sc_dim10, resolution = 0.1)
sc_dim20 <- FindClusters(sc_dim20, resolution = 0.1)
sc_dim30 <- FindClusters(sc_dim30, resolution = 0.1)
sc_dim10 <- FindClusters(sc_dim10, resolution = 0.2)
sc_dim20 <- FindClusters(sc_dim20, resolution = 0.2)
sc_dim30 <- FindClusters(sc_dim30, resolution = 0.2)
sc_dim10 <- FindClusters(sc_dim10, resolution = 0.3)
sc_dim20 <- FindClusters(sc_dim20, resolution = 0.3)
sc_dim30 <- FindClusters(sc_dim30, resolution = 0.3)

DimPlot(sc_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T) + 
    DimPlot(sc_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T) + 
    DimPlot(sc_dim30, shuffle = TRUE, raster = TRUE, label = T, repel = T)
DimPlot(sc_dim10, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T) +
    DimPlot(sc_dim20, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T) +
    DimPlot(sc_dim30, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T)
FeaturePlot(sc_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","ITGAX","CD3D","MS4A1","JCHAIN"), pt.size = 2)
FeaturePlot(sc_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","ITGAX","CD3D","MS4A1","JCHAIN"), pt.size = 2)
FeaturePlot(sc_dim30, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","ITGAX","CD3D","MS4A1","JCHAIN"), pt.size = 2)
sc <- sc_dim30
DimPlot(sc, split.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2)
DimPlot(sc, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T) + 
    DimPlot(sc, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T) + 
    DimPlot(sc, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T)

# FeaturePlot(sc, raster = T, features = c("PDGFRA","ISLR","CLIC5","PRG4","THY1","CD34","PI16","CXCL12","MMP3","MMP13","CDH11","COL1A1","COL3A1","ACAN", "ADIPOQ"))
# FeaturePlot(sc, raster = T, features = c("PTPRC","CD14","CD68","MRC1","CX3CR1","TIMD4","ITGAX","CD3D","CD4","CD8A","CD19", "MS4A1", "TNFRSF17","SDC1","FCGR3A","CEACAM8"))
# FeaturePlot(sc, raster = T, features = c("PECAM1","CDH5","PDGFRB","RGS5","ACTA2","TAGLN"))
# DotPlot(sc2, features = c("PDGFRA", "PRG4", "CLIC5", "THY1", "COL1A1",
#                           "CDH5", "RGS5", "ACTA2",
#                           "PTPRC", "CD68", "TIMD4", "MRC1", "RBP4","IL1B", "ITGAX", "CD3D", "CD19", "SDC1"),
#                         cluster.idents = TRUE) + RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred")

sc <- SetIdent(sc, value = "SCT_snn_res.0.3")


sc <- PrepSCTFindMarkers(sc)

all.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc) <- "RNA"
sc[["RNA"]] <- split(sc[["RNA"]], f = sc$orig.ident)
sc <- FindVariableFeatures(sc)
sc <- ScaleData(sc, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
sc <- JoinLayers(sc, assay = "RNA")
DefaultAssay(sc) <- "SCT"

DoHeatmap(sc, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors)
DotPlot(sc, features = all.markers_top10$gene) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred")


sc <- RenameIdents(sc, "0"="Lining-layer fibroblast 1")
sc <- RenameIdents(sc, "1"="Sublining-layer fibroblast 1")
sc <- RenameIdents(sc, "2"="Myeloid cell 1")
sc <- RenameIdents(sc, "3"="Sublining-layer fibroblast 2")
sc <- RenameIdents(sc, "4"="Endothelial cell")
sc <- RenameIdents(sc, "5"="Myeloid cell 2")
sc <- RenameIdents(sc, "6"="Mural cell")
sc <- RenameIdents(sc, "7"="Myeloid cell 3")
sc <- RenameIdents(sc, "8"="Sublining-layer fibroblast 3")
sc <- RenameIdents(sc, "9"="T/NK cell")
sc <- RenameIdents(sc, "10"="Sublining-layer fibroblast 4")
sc <- RenameIdents(sc, "11"="Lining-layer fibroblast 2")
sc <- RenameIdents(sc, "12"="Lining-layer fibroblast 3")
sc <- RenameIdents(sc, "13"="Myeloid cell 4")
sc <- RenameIdents(sc, "14"="Plasma cell")
sc <- RenameIdents(sc, "15"="Granulocyte")
sc <- RenameIdents(sc, "16"="B cell")
sc <- RenameIdents(sc, "17"="Myeloid cell 5")
sc <- RenameIdents(sc, "18"="Adipocyte")

sc@meta.data$SCT_snn_res.0.3_rename <- sc@active.ident

sc <- RenameIdents(sc, "Lining-layer fibroblast 1"="Stromal cell")
sc <- RenameIdents(sc, "Lining-layer fibroblast 2"="Stromal cell")
sc <- RenameIdents(sc, "Lining-layer fibroblast 3"="Stromal cell")
sc <- RenameIdents(sc, "Sublining-layer fibroblast 1"="Stromal cell")
sc <- RenameIdents(sc, "Sublining-layer fibroblast 2"="Stromal cell")
sc <- RenameIdents(sc, "Sublining-layer fibroblast 3"="Stromal cell")
sc <- RenameIdents(sc, "Sublining-layer fibroblast 4"="Stromal cell")
sc <- RenameIdents(sc, "Adipocyte"="Stromal cell")
sc <- RenameIdents(sc, "Myeloid cell 1"="Myeloid cell")
sc <- RenameIdents(sc, "Myeloid cell 2"="Myeloid cell")
sc <- RenameIdents(sc, "Myeloid cell 3"="Myeloid cell")
sc <- RenameIdents(sc, "Myeloid cell 4"="Myeloid cell")
sc <- RenameIdents(sc, "Myeloid cell 5"="Myeloid cell")
sc <- RenameIdents(sc, "Granulocyte"="Myeloid cell")
sc <- RenameIdents(sc, "T/NK cell"="Lymphoid cell")
sc <- RenameIdents(sc, "B cell"="Lymphoid cell")
sc <- RenameIdents(sc, "Plasma cell"="Lymphoid cell")

sc@meta.data$CellType_Class1_tmp <- sc@active.ident

all.markers_rename <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_rename_top10 <- all.markers_rename %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_rename_top30 <- all.markers_rename %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DotPlot(sc, features = all.markers_rename_top10$gene) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred")

sc_stromal <- subset(sc, idents = "Stromal cell")
sc_myeloid <- subset(sc, idents = "Myeloid cell")
sc_lymphoid <- subset(sc, idents = "Lymphoid cell")
sc_endothelial <- subset(sc, idents = "Endothelial cell")
sc_mural <- subset(sc, idents = "Mural cell")

# saveRDS(all.markers, "02_Publicdata/ver6/tmp/seuratObj_after_rename_20251022.rds")
# saveRDS(all.markers_rename, "02_Publicdata/ver6/tmp/seuratObj_after_rename_20251022.rds")

# saveRDS(sc, "02_Publicdata/ver6/tmp/seuratObj_after_rename_20251022.rds")
# saveRDS(sc_stromal, "02_Publicdata/ver6/tmp/seuratObj_stromal_20251022.rds")
# saveRDS(sc_myeloid, "02_Publicdata/ver6/tmp/seuratObj_myeloid_20251022.rds")
# saveRDS(sc_lymphoid, "02_Publicdata/ver6/tmp/seuratObj_lymphoid_20251022.rds")
# saveRDS(sc_endothelial, "02_Publicdataa/ver6/tmp/seuratObj_endothelial_20251022.rds")
# saveRDS(sc_mural, "02_Publicdata/ver6/tmp/seuratObj_mural_20251022.rds")



# Fine-grained annotation of each cell type **************
# Stromal
sc_stromal <- readRDS("02_Publicdata/ver6/tmp/seuratObj_stromal_20251022.rds")

sc_stromal
# An object of class Seurat 
# 62707 features across 56867 samples within 2 assays 
# Active assay: SCT (31224 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, integrated.rpca, umap_rpca

sc_stromal <- SetIdent(sc_stromal, value = "SCT_snn_res.0.3_rename")
sc_fib <- subset(sc_stromal, idents = c("Lining-layer fibroblast 1","Lining-layer fibroblast 2","Lining-layer fibroblast 3","Sublining-layer fibroblast 1","Sublining-layer fibroblast 2","Sublining-layer fibroblast 3","Sublining-layer fibroblast 4"))

DefaultAssay(sc_fib) <- "RNA"
sc_fib <- DietSeurat(sc_fib, assays = "RNA")
sc_fib[["RNA"]] <- split(sc_fib[["RNA"]], f = sc_fib$orig.ident)

sc_fib <- SCTransform(sc_fib, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_fib <- RunPCA(sc_fib)

sc_fib <- IntegrateLayers(object = sc_fib,
                             method = RPCAIntegration,
                             normalization.method = "SCT",
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = FALSE
                             )

sc_fib <- JoinLayers(sc_fib, assay = "RNA")

sc_fib_dim10 <- RunUMAP(sc_fib, dims = 1:10, reduction="integrated.rpca")
sc_fib_dim20 <- RunUMAP(sc_fib, dims = 1:20, reduction="integrated.rpca")
sc_fib_dim30 <- RunUMAP(sc_fib, dims = 1:30, reduction="integrated.rpca")

sc_fib_dim10 <- FindNeighbors(sc_fib_dim10, dims = 1:10, reduction="integrated.rpca")
sc_fib_dim20 <- FindNeighbors(sc_fib_dim20, dims = 1:20, reduction="integrated.rpca")
sc_fib_dim30 <- FindNeighbors(sc_fib_dim30, dims = 1:30, reduction="integrated.rpca")

sc_fib_dim10 <- FindClusters(sc_fib_dim10, resolution = 0.1)
sc_fib_dim20 <- FindClusters(sc_fib_dim20, resolution = 0.1)
sc_fib_dim30 <- FindClusters(sc_fib_dim30, resolution = 0.1)
sc_fib_dim10 <- FindClusters(sc_fib_dim10, resolution = 0.2)
sc_fib_dim20 <- FindClusters(sc_fib_dim20, resolution = 0.2)
sc_fib_dim30 <- FindClusters(sc_fib_dim30, resolution = 0.2)
sc_fib_dim10 <- FindClusters(sc_fib_dim10, resolution = 0.3)
sc_fib_dim20 <- FindClusters(sc_fib_dim20, resolution = 0.3)
sc_fib_dim30 <- FindClusters(sc_fib_dim30, resolution = 0.3)

DimPlot(sc_fib_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_fib_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_fib_dim30, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save
DimPlot(sc_fib_dim10, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_fib_dim20, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_fib_dim30, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2) 

FeaturePlot(sc_fib_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_fib_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_fib_dim30, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

sc_fib <- sc_fib_dim10

DimPlot(sc_fib, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_fib, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_fib, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

sc_fib <- SetIdent(sc_fib, value = "SCT_snn_res.0.2")

DimPlot(sc_fib, split.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 3)
FeaturePlot(sc_fib, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3) # save
FeaturePlot(sc_fib, raster = T, features = c("CLIC5","CD55","PRG4","THY1","CD34","PI16","CXCL12","CD74","IL6","MMP3","LRRC15","MMP11","SPP1","ACAN"), pt.size = 3) # save

sc_fib <- PrepSCTFindMarkers(sc_fib)
all.markers <- FindAllMarkers(sc_fib, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DotPlot(sc_fib, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_fib <- RenameIdents(sc_fib, "0"="Sublining-layer fibroblast")
sc_fib <- RenameIdents(sc_fib, "1"="Sublining-layer fibroblast")
sc_fib <- RenameIdents(sc_fib, "2"="Lining-layer fibroblast")
sc_fib <- RenameIdents(sc_fib, "3"="Lining-layer fibroblast")
sc_fib <- RenameIdents(sc_fib, "4"="Sublining-layer fibroblast")
sc_fib <- RenameIdents(sc_fib, "5"="Sublining-layer fibroblast")
sc_fib <- RenameIdents(sc_fib, "6"="Lining-layer fibroblast")
sc_fib <- RenameIdents(sc_fib, "7"="Doublet")

table(sc_fib@active.ident)
#                    Doublet    Lining-layer fibroblast Sublining-layer fibroblast 
#                       1414                      21782                      33593 

DimPlot(sc_fib, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save


sc_fib_2 <- subset(sc_fib, idents = c("Lining-layer fibroblast","Sublining-layer fibroblast"))

DefaultAssay(sc_fib_2) <- "RNA"
sc_fib_2 <- DietSeurat(sc_fib_2, assays = "RNA")
sc_fib_2[["RNA"]] <- split(sc_fib_2[["RNA"]], f = sc_fib_2$orig.ident)

sc_fib_2 <- SCTransform(sc_fib_2, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_fib_2 <- RunPCA(sc_fib_2)

sc_fib_2 <- IntegrateLayers(object = sc_fib_2,
                             method = RPCAIntegration,
                             normalization.method = "SCT",
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = FALSE
                             )

sc_fib_2 <- JoinLayers(sc_fib_2, assay = "RNA")

sc_fib_2_dim10 <- RunUMAP(sc_fib_2, dims = 1:10, reduction="integrated.rpca")
sc_fib_2_dim20 <- RunUMAP(sc_fib_2, dims = 1:20, reduction="integrated.rpca")
sc_fib_2_dim30 <- RunUMAP(sc_fib_2, dims = 1:30, reduction="integrated.rpca")

sc_fib_2_dim10 <- FindNeighbors(sc_fib_2_dim10, dims = 1:10, reduction="integrated.rpca")
sc_fib_2_dim20 <- FindNeighbors(sc_fib_2_dim20, dims = 1:20, reduction="integrated.rpca")
sc_fib_2_dim30 <- FindNeighbors(sc_fib_2_dim30, dims = 1:20, reduction="integrated.rpca")

sc_fib_2_dim10 <- FindClusters(sc_fib_2_dim10, resolution = 0.1)
sc_fib_2_dim20 <- FindClusters(sc_fib_2_dim20, resolution = 0.1)
sc_fib_2_dim30 <- FindClusters(sc_fib_2_dim30, resolution = 0.1)
sc_fib_2_dim10 <- FindClusters(sc_fib_2_dim10, resolution = 0.2)
sc_fib_2_dim20 <- FindClusters(sc_fib_2_dim20, resolution = 0.2)
sc_fib_2_dim30 <- FindClusters(sc_fib_2_dim30, resolution = 0.2)
sc_fib_2_dim10 <- FindClusters(sc_fib_2_dim10, resolution = 0.3)
sc_fib_2_dim20 <- FindClusters(sc_fib_2_dim20, resolution = 0.3)
sc_fib_2_dim30 <- FindClusters(sc_fib_2_dim30, resolution = 0.3)

DimPlot(sc_fib_2_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_fib_2_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_fib_2_dim30, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

FeaturePlot(sc_fib_2_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_fib_2_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_fib_2_dim30, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

FeaturePlot(sc_fib_2_dim10, raster = T, features = c("CLIC5","CD55","PRG4","THY1","CD34","PI16","CXCL12","CD74","IL6","MMP3","LRRC15","MMP11","SPP1","ACAN"), pt.size = 3)
FeaturePlot(sc_fib_2_dim20, raster = T, features = c("CLIC5","CD55","PRG4","THY1","CD34","PI16","CXCL12","CD74","IL6","MMP3","LRRC15","MMP11","SPP1","ACAN"), pt.size = 3)
FeaturePlot(sc_fib_2_dim30, raster = T, features = c("CLIC5","CD55","PRG4","THY1","CD34","PI16","CXCL12","CD74","IL6","MMP3","LRRC15","MMP11","SPP1","ACAN"), pt.size = 3)

sc_fib_2 <- sc_fib_2_dim10

DimPlot(sc_fib_2, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_fib_2, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_fib_2, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

FeaturePlot(sc_fib_2, raster = T, features = c("CLIC5","CD55","PRG4","THY1","CD34","MFAP5","PI16","APOD","CXCL14","CXCL12","CD74","IL6","MMP3","LRRC15","SPP1","ACAN"), pt.size = 3) # save

sc_fib_2 <- SetIdent(sc_fib_2, value = "SCT_snn_res.0.2")

DotPlot(sc_fib_2, features = c("CLIC5","CD55","PRG4","LRRC15","MMP3","THY1","CD34","MFAP5","PI16","APOD","CXCL14","CXCL12","CD74","IL6","SPP1","ACAN")) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_fib_2 <- PrepSCTFindMarkers(sc_fib_2)
all.markers <- FindAllMarkers(sc_fib_2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc_fib_2) <- "RNA"
sc_fib_2[["RNA"]] <- split(sc_fib_2[["RNA"]], f = sc_fib_2$orig.ident)
sc_fib_2 <- NormalizeData(sc_fib_2)
sc_fib_2 <- FindVariableFeatures(sc_fib_2)
sc_fib_2 <- ScaleData(sc_fib_2)
DefaultAssay(sc_fib_2) <- "SCT"
sc_fib_2 <- JoinLayers(sc_fib_2, assay = "RNA")

DoHeatmap(sc_fib_2, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_fib_2, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Stromal/all.markers_Stromal_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Stromal/top10.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Stromal/top30.txt", sep = "\t")

lining_fib.markers <- FindMarkers(sc_fib_2, ident.1 = "1", ident.2 = "4", min.pct = 0.1, logfc.threshold = 0.5)

sc_fib_2 <- RenameIdents(sc_fib_2, "0"="Sublining-layer fibroblast-MFAP5+PI16+")
sc_fib_2 <- RenameIdents(sc_fib_2, "1"="Lining-layer fibroblast-2")
sc_fib_2 <- RenameIdents(sc_fib_2, "2"="Sublining-layer fibroblast-APOE+CXCL12+")
sc_fib_2 <- RenameIdents(sc_fib_2, "3"="Sublining-layer fibroblast-COMP+")
sc_fib_2 <- RenameIdents(sc_fib_2, "4"="Lining-layer fibroblast-1")
sc_fib_2 <- RenameIdents(sc_fib_2, "5"="Sublining-layer fibroblast-APOD+CXCL14+")
sc_fib_2 <- RenameIdents(sc_fib_2, "6"="Lining-layer fibroblast-MMP3+")

DimPlot(sc_fib_2, shuffle = TRUE, raster = TRUE, repel = T, pt.size = 2) +
    DimPlot(sc_fib_2, group.by = "orig.ident", shuffle = TRUE, raster = TRUE, repel = T, pt.size = 2) # save

sc_stromal$CellType_tmp <- NA_character_
sc_stromal$CellType_tmp[sc_stromal$SCT_snn_res.0.3_rename == "Adipocyte"] <- "Adipocyte"
sc_stromal$CellType_tmp[colnames(sc_fib)[Idents(sc_fib) == "Doublet"]] <- "Doublet"
fib2_map <- setNames(as.character(Idents(sc_fib_2)), colnames(sc_fib_2))
sc_stromal$CellType_tmp[names(fib2_map)] <- fib2_map

 table(sc_stromal$CellType_tmp)
#                               Adipocyte                                 Doublet               Lining-layer fibroblast-1 
#                                     78                                    1414                                    4287 
#               Lining-layer fibroblast-2           Lining-layer fibroblast-MMP3+ Sublining-layer fibroblast-APOD+CXCL14+ 
#                                   14745                                    3240                                    3790 
# Sublining-layer fibroblast-APOE+CXCL12+        Sublining-layer fibroblast-COMP+  Sublining-layer fibroblast-MFAP5+PI16+ 
#                                    8571                                    4590                                   16152

# saveRDS(sc_fib, "02_Publicdata/ver6/Stromal/seuratObj_Fibro_20251022.rds")
# saveRDS(sc_fib_2, "02_Publicdata/ver6/Stromal/seuratObj_Fibro_remove_doublet_20251022.rds")
# saveRDS(sc_stromal, "02_Publicdata/ver6/tmp/seuratObj_stromal_20251022_update.rds")


# Myeloid cell
sc_myeloid <- readRDS("02_Publicdata/ver6/tmp/seuratObj_myeloid_20251022.rds")

sc_myeloid
# An object of class Seurat 
# 62707 features across 27333 samples within 2 assays 
# Active assay: SCT (31224 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
#  3 dimensional reductions calculated: pca, integrated.rpca, umap.rpca

sc_myeloid <- SetIdent(sc_myeloid, value = "SCT_snn_res.0.3_rename")
sc_MonoMacroDC <- subset(sc_myeloid, idents = c("Myeloid cell 1","Myeloid cell 2","Myeloid cell 3","Myeloid cell 4","Myeloid cell 5"))

DefaultAssay(sc_MonoMacroDC) <- "RNA"
sc_MonoMacroDC <- DietSeurat(sc_MonoMacroDC, assays = "RNA")
sc_MonoMacroDC[["RNA"]] <- split(sc_MonoMacroDC[["RNA"]], f = sc_MonoMacroDC$orig.ident)

sc_MonoMacroDC <- SCTransform(sc_MonoMacroDC, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_MonoMacroDC <- RunPCA(sc_MonoMacroDC)

sc_MonoMacroDC <- IntegrateLayers(object = sc_MonoMacroDC,
                             method = RPCAIntegration,
                             normalization.method = "SCT",
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = FALSE
                             )

sc_MonoMacroDC <- JoinLayers(sc_MonoMacroDC, assay = "RNA")

sc_MonoMacroDC_dim10 <- RunUMAP(sc_MonoMacroDC, dims = 1:10, reduction="integrated.rpca")
sc_MonoMacroDC_dim20 <- RunUMAP(sc_MonoMacroDC, dims = 1:20, reduction="integrated.rpca")
sc_MonoMacroDC_dim30 <- RunUMAP(sc_MonoMacroDC, dims = 1:30, reduction="integrated.rpca")

sc_MonoMacroDC_dim10 <- FindNeighbors(sc_MonoMacroDC_dim10, dims = 1:10, reduction="integrated.rpca")
sc_MonoMacroDC_dim20 <- FindNeighbors(sc_MonoMacroDC_dim20, dims = 1:20, reduction="integrated.rpca")
sc_MonoMacroDC_dim30 <- FindNeighbors(sc_MonoMacroDC_dim30, dims = 1:30, reduction="integrated.rpca")

sc_MonoMacroDC_dim10 <- FindClusters(sc_MonoMacroDC_dim10, resolution = 0.1)
sc_MonoMacroDC_dim20 <- FindClusters(sc_MonoMacroDC_dim20, resolution = 0.1)
sc_MonoMacroDC_dim30 <- FindClusters(sc_MonoMacroDC_dim30, resolution = 0.1)
sc_MonoMacroDC_dim10 <- FindClusters(sc_MonoMacroDC_dim10, resolution = 0.2)
sc_MonoMacroDC_dim20 <- FindClusters(sc_MonoMacroDC_dim20, resolution = 0.2)
sc_MonoMacroDC_dim30 <- FindClusters(sc_MonoMacroDC_dim30, resolution = 0.2)
sc_MonoMacroDC_dim10 <- FindClusters(sc_MonoMacroDC_dim10, resolution = 0.3)
sc_MonoMacroDC_dim20 <- FindClusters(sc_MonoMacroDC_dim20, resolution = 0.3)
sc_MonoMacroDC_dim30 <- FindClusters(sc_MonoMacroDC_dim30, resolution = 0.3)

DimPlot(sc_MonoMacroDC_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_MonoMacroDC_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_MonoMacroDC_dim30, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save
DimPlot(sc_MonoMacroDC_dim10, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_MonoMacroDC_dim20, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_MonoMacroDC_dim30, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2)
FeaturePlot(sc_MonoMacroDC_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_MonoMacroDC_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_MonoMacroDC_dim30, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

sc_MonoMacroDC <- sc_MonoMacroDC_dim10

DimPlot(sc_MonoMacroDC, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_MonoMacroDC, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_MonoMacroDC, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

sc_MonoMacroDC <- SetIdent(sc_MonoMacroDC, value = "SCT_snn_res.0.2")

DimPlot(sc_MonoMacroDC, split.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 3)
FeaturePlot(sc_MonoMacroDC, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3) # save
FeaturePlot(sc_MonoMacroDC, raster = T, features = c("HLA-DRA","CD68","ADGRE1","MRC1","CX3CR1","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33"), pt.size = 3) # save

sc_MonoMacroDC <- PrepSCTFindMarkers(sc_MonoMacroDC)
all.markers <- FindAllMarkers(sc_MonoMacroDC, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DotPlot(sc_MonoMacroDC, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_MonoMacroDC <- RenameIdents(sc_MonoMacroDC, "0"="Macrophage")
sc_MonoMacroDC <- RenameIdents(sc_MonoMacroDC, "1"="Macrophage")
sc_MonoMacroDC <- RenameIdents(sc_MonoMacroDC, "2"="Macrophage")
sc_MonoMacroDC <- RenameIdents(sc_MonoMacroDC, "3"="Mono/DC")
sc_MonoMacroDC <- RenameIdents(sc_MonoMacroDC, "4"="Macrophage")
sc_MonoMacroDC <- RenameIdents(sc_MonoMacroDC, "5"="Macrophage")
sc_MonoMacroDC <- RenameIdents(sc_MonoMacroDC, "6"="Doublet")

table(sc_MonoMacroDC@active.ident)
#    Doublet Macrophage    Mono/DC 
#       1421      22705       2739 

DimPlot(sc_MonoMacroDC, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save


sc_MonoMacroDC_2 <- subset(sc_MonoMacroDC, idents = c("Macrophage","Mono/DC"))

DefaultAssay(sc_MonoMacroDC_2) <- "RNA"
sc_MonoMacroDC_2 <- DietSeurat(sc_MonoMacroDC_2, assays = "RNA")
sc_MonoMacroDC_2[["RNA"]] <- split(sc_MonoMacroDC_2[["RNA"]], f = sc_MonoMacroDC_2$orig.ident)

sc_MonoMacroDC_2 <- SCTransform(sc_MonoMacroDC_2, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_MonoMacroDC_2 <- RunPCA(sc_MonoMacroDC_2)

sc_MonoMacroDC_2 <- IntegrateLayers(object = sc_MonoMacroDC_2,
                             method = RPCAIntegration,
                             normalization.method = "SCT",
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = FALSE
                             )

sc_MonoMacroDC_2 <- JoinLayers(sc_MonoMacroDC_2, assay = "RNA")

sc_MonoMacroDC_2 <- RunPCA(sc_MonoMacroDC_2)

sc_MonoMacroDC_2_dim10 <- RunUMAP(sc_MonoMacroDC_2, dims = 1:10, reduction="integrated.rpca")
sc_MonoMacroDC_2_dim20 <- RunUMAP(sc_MonoMacroDC_2, dims = 1:20, reduction="integrated.rpca")

sc_MonoMacroDC_2_dim10 <- FindNeighbors(sc_MonoMacroDC_2_dim10, dims = 1:10, reduction="integrated.rpca")
sc_MonoMacroDC_2_dim20 <- FindNeighbors(sc_MonoMacroDC_2_dim20, dims = 1:20, reduction="integrated.rpca")

sc_MonoMacroDC_2_dim10 <- FindClusters(sc_MonoMacroDC_2_dim10, resolution = 0.1)
sc_MonoMacroDC_2_dim20 <- FindClusters(sc_MonoMacroDC_2_dim20, resolution = 0.1)
sc_MonoMacroDC_2_dim10 <- FindClusters(sc_MonoMacroDC_2_dim10, resolution = 0.2)
sc_MonoMacroDC_2_dim20 <- FindClusters(sc_MonoMacroDC_2_dim20, resolution = 0.2)
sc_MonoMacroDC_2_dim10 <- FindClusters(sc_MonoMacroDC_2_dim10, resolution = 0.3)
sc_MonoMacroDC_2_dim20 <- FindClusters(sc_MonoMacroDC_2_dim20, resolution = 0.3)
sc_MonoMacroDC_2_dim10 <- FindClusters(sc_MonoMacroDC_2_dim10, resolution = 0.4)
sc_MonoMacroDC_2_dim20 <- FindClusters(sc_MonoMacroDC_2_dim20, resolution = 0.4)
sc_MonoMacroDC_2_dim10 <- FindClusters(sc_MonoMacroDC_2_dim10, resolution = 0.5)
sc_MonoMacroDC_2_dim20 <- FindClusters(sc_MonoMacroDC_2_dim20, resolution = 0.5)
sc_MonoMacroDC_2_dim10 <- FindClusters(sc_MonoMacroDC_2_dim10, resolution = 0.6)
sc_MonoMacroDC_2_dim20 <- FindClusters(sc_MonoMacroDC_2_dim20, resolution = 0.6)

DimPlot(sc_MonoMacroDC_2_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_MonoMacroDC_2_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save
FeaturePlot(sc_MonoMacroDC_2_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_MonoMacroDC_2_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

FeaturePlot(sc_MonoMacroDC_2_dim10, raster = T, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33"), pt.size = 3)
FeaturePlot(sc_MonoMacroDC_2_dim20, raster = T, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33"), pt.size = 3)

sc_MonoMacroDC_2 <- sc_MonoMacroDC_2_dim20

DimPlot(sc_MonoMacroDC_2, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_MonoMacroDC_2, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_MonoMacroDC_2, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_MonoMacroDC_2, group.by = "SCT_snn_res.0.4", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_MonoMacroDC_2, group.by = "SCT_snn_res.0.5", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_MonoMacroDC_2, group.by = "SCT_snn_res.0.6", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

FeaturePlot(sc_MonoMacroDC_2, raster = T, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33"), pt.size = 3) # save

sc_MonoMacroDC_2 <- SetIdent(sc_MonoMacroDC_2, value = "SCT_snn_res.0.1")

DotPlot(sc_MonoMacroDC_2, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33")) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_MonoMacroDC_2 <- PrepSCTFindMarkers(sc_MonoMacroDC_2)
all.markers <- FindAllMarkers(sc_MonoMacroDC_2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc_MonoMacroDC_2) <- "RNA"
sc_MonoMacroDC_2[["RNA"]] <- split(sc_MonoMacroDC_2[["RNA"]], f = sc_MonoMacroDC_2$orig.ident)
sc_MonoMacroDC_2 <- NormalizeData(sc_MonoMacroDC_2)
sc_MonoMacroDC_2 <- FindVariableFeatures(sc_MonoMacroDC_2)
sc_MonoMacroDC_2 <- ScaleData(sc_MonoMacroDC_2)
DefaultAssay(sc_MonoMacroDC_2) <- "SCT"
sc_MonoMacroDC_2 <- JoinLayers(sc_MonoMacroDC_2, assay = "RNA")

DoHeatmap(sc_MonoMacroDC_2, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_MonoMacroDC_2, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Myeloid/all.markers_Myeloid_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Myeloid/top10.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Myeloid/top30.txt", sep = "\t")

sc_MonoMacroDC_2 <- RenameIdents(sc_MonoMacroDC_2, "0"="Macrophage")
sc_MonoMacroDC_2 <- RenameIdents(sc_MonoMacroDC_2, "1"="Macrophage")
sc_MonoMacroDC_2 <- RenameIdents(sc_MonoMacroDC_2, "2"="Dendritic cell")
sc_MonoMacroDC_2 <- RenameIdents(sc_MonoMacroDC_2, "3"="Macrophage")


sc_DC <- subset(sc_MonoMacroDC_2, idents = "Dendritic cell")

DefaultAssay(sc_DC) <- "RNA"
sc_DC <- DietSeurat(sc_DC, assays = "RNA")
sc_DC[["RNA"]] <- split(sc_DC[["RNA"]], f = sc_DC$orig.ident)

sc_DC <- SCTransform(sc_DC, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_DC <- RunPCA(sc_DC)

sc_DC <- IntegrateLayers(object = sc_DC,
                             method = RPCAIntegration,
                             normalization.method = "SCT",
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = FALSE,
                             k.weight = 59 # default: 100. An error will occur if the value is below the default threshold (Number of anchor cells is less than k.weight. Consider lowering k.weight to less than 59 or increase k.anchor.).
                             )

sc_DC <- JoinLayers(sc_DC, assay = "RNA")

sc_DC_dim10 <- RunUMAP(sc_DC, dims = 1:10, reduction="integrated.rpca")
sc_DC_dim20 <- RunUMAP(sc_DC, dims = 1:20, reduction="integrated.rpca")

sc_DC_dim10 <- FindNeighbors(sc_DC_dim10, dims = 1:10, reduction="integrated.rpca")
sc_DC_dim20 <- FindNeighbors(sc_DC_dim20, dims = 1:20, reduction="integrated.rpca")

sc_DC_dim10 <- FindClusters(sc_DC_dim10, resolution = 0.1)
sc_DC_dim20 <- FindClusters(sc_DC_dim20, resolution = 0.1)
sc_DC_dim10 <- FindClusters(sc_DC_dim10, resolution = 0.2)
sc_DC_dim20 <- FindClusters(sc_DC_dim20, resolution = 0.2)
sc_DC_dim10 <- FindClusters(sc_DC_dim10, resolution = 0.3)
sc_DC_dim20 <- FindClusters(sc_DC_dim20, resolution = 0.3)
sc_DC_dim10 <- FindClusters(sc_DC_dim10, resolution = 0.4)
sc_DC_dim20 <- FindClusters(sc_DC_dim20, resolution = 0.4)

DimPlot(sc_DC_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) + 
    DimPlot(sc_DC_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) # save
FeaturePlot(sc_DC_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_DC_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

FeaturePlot(sc_DC_dim10, raster = T, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33"), pt.size = 3)
FeaturePlot(sc_DC_dim20, raster = T, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33"), pt.size = 3)

sc_DC <- sc_DC_dim20

DimPlot(sc_DC, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) + 
    DimPlot(sc_DC, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) + 
    DimPlot(sc_DC, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) +
    DimPlot(sc_DC, group.by = "SCT_snn_res.0.4", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3)  # save

FeaturePlot(sc_DC, raster = T, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33"), pt.size = 3)

sc_DC <- SetIdent(sc_DC, value = "SCT_snn_res.0.4")

DotPlot(sc_DC, features = c("HLA-DRA","CD68","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33")) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_DC <- PrepSCTFindMarkers(sc_DC)
all.markers <- FindAllMarkers(sc_DC, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc_DC) <- "RNA"
sc_DC[["RNA"]] <- split(sc_DC[["RNA"]], f = sc_DC$orig.ident)
sc_DC <- NormalizeData(sc_DC)
sc_DC<- FindVariableFeatures(sc_DC)
sc_DC <- ScaleData(sc_DC)
DefaultAssay(sc_DC) <- "SCT"
sc_DC <- JoinLayers(sc_DC, assay = "RNA")

DoHeatmap(sc_DC, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_DC, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Myeloid/DC/all.markers_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Myeloid/DC/top10.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Myeloid/DC/top30.txt", sep = "\t")

sc_DC <- RenameIdents(sc_DC, "0"="cDC2")
sc_DC <- RenameIdents(sc_DC, "1"="Monocyte-Classical")
sc_DC <- RenameIdents(sc_DC, "2"="Monocyte-Int")
sc_DC <- RenameIdents(sc_DC, "3"="cDC2")
sc_DC <- RenameIdents(sc_DC, "4"="mregDC")
sc_DC <- RenameIdents(sc_DC, "5"="Monocyte-nonClassical")
sc_DC <- RenameIdents(sc_DC, "6"="cDC1")

all.markers <- FindAllMarkers(sc_DC, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DoHeatmap(sc_DC, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_DC, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Myeloid/DC/all.markers_update_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Myeloid/DC/top10_update.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Myeloid/DC/top30_update.txt", sep = "\t")


sc_Mac <- subset(sc_MonoMacroDC_2, idents = "Macrophage")

DefaultAssay(sc_Mac) <- "RNA"
sc_sc_MacDC <- DietSeurat(sc_Mac, assays = "RNA")
sc_Mac[["RNA"]] <- split(sc_Mac[["RNA"]], f = sc_Mac$orig.ident)

sc_Mac <- SCTransform(sc_Mac, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_Mac <- RunPCA(sc_Mac)

sc_Mac <- IntegrateLayers(object = sc_Mac,
                             method = RPCAIntegration,
                             normalization.method = "SCT",
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = FALSE
                             )

sc_Mac <- JoinLayers(sc_Mac, assay = "RNA")

sc_Mac_dim10 <- RunUMAP(sc_Mac, dims = 1:10, reduction="integrated.rpca")
sc_Mac_dim20 <- RunUMAP(sc_Mac, dims = 1:20, reduction="integrated.rpca")

sc_Mac_dim10 <- FindNeighbors(sc_Mac_dim10, dims = 1:10, reduction="integrated.rpca")
sc_Mac_dim20 <- FindNeighbors(sc_Mac_dim20, dims = 1:20, reduction="integrated.rpca")

sc_Mac_dim10 <- FindClusters(sc_Mac_dim10, resolution = 0.1)
sc_Mac_dim20 <- FindClusters(sc_Mac_dim20, resolution = 0.1)
sc_Mac_dim10 <- FindClusters(sc_Mac_dim10, resolution = 0.2)
sc_Mac_dim20 <- FindClusters(sc_Mac_dim20, resolution = 0.2)
sc_Mac_dim10 <- FindClusters(sc_Mac_dim10, resolution = 0.3)
sc_Mac_dim20 <- FindClusters(sc_Mac_dim20, resolution = 0.3)
#sc_Mac_dim10 <- FindClusters(sc_Mac_dim10, resolution = 0.4)
#sc_Mac_dim20 <- FindClusters(sc_Mac_dim20, resolution = 0.4)

DimPlot(sc_Mac_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_Mac_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save
FeaturePlot(sc_Mac_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_Mac_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

FeaturePlot(sc_Mac_dim10, raster = T, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33"), pt.size = 3)
FeaturePlot(sc_Mac_dim20, raster = T, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","CX3CR1","TIMD4","CD14","FCGR3A","CCR2","CD1A","CD1C","ITGAX","ITGAM","ANPEP","CD33"), pt.size = 3)

sc_Mac <- sc_Mac_dim20

DimPlot(sc_Mac, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_Mac, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_Mac, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

FeaturePlot(sc_Mac, raster = T, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","TREM2","TIMD4","LYVE1","FOLR2","ISG15","CLEC10A","S100A12","SPP1"), pt.size = 3)

sc_Mac <- SetIdent(sc_Mac, value = "SCT_snn_res.0.3")

DotPlot(sc_Mac, features = c("HLA-DRA","CD68","NOS2","MRC1","MERTK","TREM2","TIMD4","LYVE1","FOLR2","ISG15","CLEC10A","S100A12","SPP1")) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_Mac <- PrepSCTFindMarkers(sc_Mac)
all.markers <- FindAllMarkers(sc_Mac, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc_Mac) <- "RNA"
sc_Mac[["RNA"]] <- split(sc_Mac[["RNA"]], f = sc_Mac$orig.ident)
sc_Mac <- NormalizeData(sc_Mac)
sc_Mac<- FindVariableFeatures(sc_Mac)
sc_Mac <- ScaleData(sc_Mac)
DefaultAssay(sc_Mac) <- "SCT"
sc_Mac <- JoinLayers(sc_Mac, assay = "RNA")

DoHeatmap(sc_Mac, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_Mac, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Myeloid/Macrophage/all.markers_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Myeloid/Macrophage/top10.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Myeloid/Macrophage/top30.txt", sep = "\t")

sc_Mac <- RenameIdents(sc_Mac, "0"="Macrophage-TIMD4+")
sc_Mac <- RenameIdents(sc_Mac, "1"="Macrophage-Chemokine high")
sc_Mac <- RenameIdents(sc_Mac, "2"="Macrophage-Chemokine high")
sc_Mac <- RenameIdents(sc_Mac, "3"="Macrophage-LYVE1+")
sc_Mac <- RenameIdents(sc_Mac, "4"="Macrophage-TIMD4+APOE+")
sc_Mac <- RenameIdents(sc_Mac, "5"="Macrophage-CLEC10A+")
sc_Mac <- RenameIdents(sc_Mac, "6"="Macrophage-Chemokine high")
sc_Mac <- RenameIdents(sc_Mac, "7"="Macrophage-SPP1+")
sc_Mac <- RenameIdents(sc_Mac, "8"="Macrophage-ISG15+")

all.markers <- FindAllMarkers(sc_Mac, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DoHeatmap(sc_Mac, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_Mac, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Myeloid/Macrophage/all.markers_update_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Myeloid/Macrophage/top10_update.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Myeloid/Macrophage/top30_update.txt", sep = "\t")


sc_gra <- subset(sc_myeloid, idents = "Granulocyte")
DefaultAssay(sc_gra) <- "RNA"
sc_gra <- DietSeurat(sc_gra, assays = "RNA")
sc_gra [["RNA"]] <- split(sc_gra [["RNA"]], f = sc_gra $orig.ident)
sc_gra  <- SCTransform(sc_gra , vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)
sc_gra <- JoinLayers(sc_gra, assay = "RNA")
sc_gra <- RunPCA(sc_gra)
sc_gra <- RunUMAP(sc_gra, dims = 1:10)
sc_gra <- FindNeighbors(sc_gra, dims = 1:10)
sc_gra <- FindClusters(sc_gra, resolution = 0.1)
DimPlot(sc_gra) # save
FeaturePlot(sc_gra, features = c("TPSAB1","FCER1A","KIT")) # save


sc_myeloid$CellType_tmp <- NA_character_
sc_myeloid$CellType_tmp[sc_myeloid$SCT_snn_res.0.3_rename == "Granulocyte"] <- "Mast cell"
sc_myeloid$CellType_tmp[colnames(sc_MonoMacroDC)[Idents(sc_MonoMacroDC) == "Doublet"]] <- "Doublet"
DC_map <- setNames(as.character(Idents(sc_DC)), colnames(sc_DC))
Mac_map <- setNames(as.character(Idents(sc_Mac)), colnames(sc_Mac))
sc_myeloid$CellType_tmp[names(DC_map)] <- DC_map
sc_myeloid$CellType_tmp[names(Mac_map)] <- Mac_map

table(sc_myeloid$CellType_tmp)
#                      cDC1                      cDC2                   Doublet Macrophage-Chemokine high       Macrophage-CLEC10A+ 
#                        72                      1098                      1421                      9846                      1006 
#         Macrophage-ISG15+         Macrophage-LYVE1+          Macrophage-SPP1+         Macrophage-TIMD4+    Macrophage-TIMD4+APOE+ 
#                       154                      2967                       390                      6360                      2283 
#                 Mast cell        Monocyte-Classical              Monocyte-Int     Monocyte-nonClassical                    mregDC 
#                       468                       505                       384                       110                       269 

# saveRDS(sc_MonoMacroDC, "02_Publicdata/ver6/Myeloid/seuratObj_Myeloid_20251022.rds")
# saveRDS(sc_MonoMacroDC_2, "02_Publicdata/ver6/Myeloid/seuratObj_Myeloid_remove_doublet_20251022.rds")
# saveRDS(sc_DC, "02_Publicdata/ver6/Myeloid/seuratObj_Myeloid_DC_20251022.rds")
# saveRDS(sc_Mac, "02_Publicdata/ver6/Myeloid/seuratObj_Myeloid_Macrophage_20251022.rds")
# saveRDS(sc_gra, "02_Publicdata/ver6/tmp/seuratObj_Myeloid_granulocyte_20251022.rds")
# saveRDS(sc_myeloid, "02_Publicdata/ver6/tmp/seuratObj_myeloid_20251022_update.rds")


# Lymphoid cell
sc_lymphoid <- readRDS("02_Publicdata/ver6/tmp/seuratObj_lymphoid_20251022.rds")

sc_lymphoid
# An object of class Seurat 
# 62707 features across 6411 samples within 2 assays 
# Active assay: SCT (31224 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
#  3 dimensional reductions calculated: pca, integrated.rpca, umap.rpca

DefaultAssay(sc_lymphoid) <- "RNA"
sc_lymphoid <- DietSeurat(sc_lymphoid, assays = "RNA")
sc_lymphoid[["RNA"]] <- split(sc_lymphoid[["RNA"]], f = sc_lymphoid$orig.ident)

sc_lymphoid <- SCTransform(sc_lymphoid, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_lymphoid <- JoinLayers(sc_lymphoid, assay = "RNA")

sc_lymphoid <- RunPCA(sc_lymphoid)

sc_lymphoid_dim10 <- RunUMAP(sc_lymphoid, dims = 1:10)
sc_lymphoid_dim20 <- RunUMAP(sc_lymphoid, dims = 1:20)
sc_lymphoid_dim30 <- RunUMAP(sc_lymphoid, dims = 1:30)

sc_lymphoid_dim10 <- FindNeighbors(sc_lymphoid_dim10, dims = 1:10)
sc_lymphoid_dim20 <- FindNeighbors(sc_lymphoid_dim20, dims = 1:20)
sc_lymphoid_dim30 <- FindNeighbors(sc_lymphoid_dim30, dims = 1:30)

sc_lymphoid_dim10 <- FindClusters(sc_lymphoid_dim10, resolution = 0.1)
sc_lymphoid_dim20 <- FindClusters(sc_lymphoid_dim20, resolution = 0.1)
sc_lymphoid_dim30 <- FindClusters(sc_lymphoid_dim30, resolution = 0.1)
sc_lymphoid_dim10 <- FindClusters(sc_lymphoid_dim10, resolution = 0.2)
sc_lymphoid_dim20 <- FindClusters(sc_lymphoid_dim20, resolution = 0.2)
sc_lymphoid_dim30 <- FindClusters(sc_lymphoid_dim30, resolution = 0.2)
sc_lymphoid_dim10 <- FindClusters(sc_lymphoid_dim10, resolution = 0.3)
sc_lymphoid_dim20 <- FindClusters(sc_lymphoid_dim20, resolution = 0.3)
sc_lymphoid_dim30 <- FindClusters(sc_lymphoid_dim30, resolution = 0.3)

DimPlot(sc_lymphoid_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_lymphoid_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_lymphoid_dim30, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save
DimPlot(sc_lymphoid_dim10, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_lymphoid_dim20, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_lymphoid_dim30, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2)
FeaturePlot(sc_lymphoid_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_lymphoid_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_lymphoid_dim30, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

sc_lymphoid <- sc_lymphoid_dim30

DimPlot(sc_lymphoid, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_lymphoid, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_lymphoid, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

sc_lymphoid <- SetIdent(sc_lymphoid, value = "SCT_snn_res.0.3")

DimPlot(sc_lymphoid, split.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 3)
FeaturePlot(sc_lymphoid, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3) # save
FeaturePlot(sc_lymphoid, raster = T, features = c("CD3E","CD4","CD8A","FOXP3","CD19","MS4A1","SDC1","TNFRSF17","IL3RA"), pt.size = 3) # save

sc_lymphoid <- PrepSCTFindMarkers(sc_lymphoid)
all.markers <- FindAllMarkers(sc_lymphoid, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DotPlot(sc_lymphoid, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_lymphoid <- RenameIdents(sc_lymphoid, "0"="T cell-nonTreg")
sc_lymphoid <- RenameIdents(sc_lymphoid, "1"="T cell-nonTreg")
sc_lymphoid <- RenameIdents(sc_lymphoid, "2"="T cell-nonTreg")
sc_lymphoid <- RenameIdents(sc_lymphoid, "3"="Plasma")
sc_lymphoid <- RenameIdents(sc_lymphoid, "4"="B cell")
sc_lymphoid <- RenameIdents(sc_lymphoid, "5"="Doublet")
sc_lymphoid <- RenameIdents(sc_lymphoid, "6"="T cell-nonTreg")
sc_lymphoid <- RenameIdents(sc_lymphoid, "7"="Doublet")
sc_lymphoid <- RenameIdents(sc_lymphoid, "8"="Doublet")
sc_lymphoid <- RenameIdents(sc_lymphoid, "9"="pDC")
sc_lymphoid <- RenameIdents(sc_lymphoid, "10"="Treg")
sc_lymphoid <- RenameIdents(sc_lymphoid, "11"="Plasma")

table(sc_lymphoid@active.ident)
#         Plasma           Treg            pDC        Doublet T cell-nonTreg         B cell 
#            614             51             51            499           4819            377 

DimPlot(sc_lymphoid, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save


sc_T_NK <- subset(sc_lymphoid, idents = "T cell-nonTreg")

DefaultAssay(sc_T_NK) <- "RNA"
sc_T_NK <- DietSeurat(sc_T_NK, assays = "RNA")
sc_T_NK[["RNA"]] <- split(sc_T_NK[["RNA"]], f = sc_T_NK$orig.ident)

sc_T_NK <- SCTransform(sc_T_NK, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_T_NK <- JoinLayers(sc_T_NK, assay = "RNA")

sc_T_NK <- RunPCA(sc_T_NK)

sc_T_NK_dim10 <- RunUMAP(sc_T_NK, dims = 1:10)
sc_T_NK_dim20 <- RunUMAP(sc_T_NK, dims = 1:20)
sc_T_NK_dim30 <- RunUMAP(sc_T_NK, dims = 1:30)

sc_T_NK_dim10 <- FindNeighbors(sc_T_NK_dim10, dims = 1:10)
sc_T_NK_dim20 <- FindNeighbors(sc_T_NK_dim20, dims = 1:20)
sc_T_NK_dim30 <- FindNeighbors(sc_T_NK_dim30, dims = 1:30)

sc_T_NK_dim10 <- FindClusters(sc_T_NK_dim10, resolution = 0.1)
sc_T_NK_dim20 <- FindClusters(sc_T_NK_dim20, resolution = 0.1)
sc_T_NK_dim30 <- FindClusters(sc_T_NK_dim30, resolution = 0.1)
sc_T_NK_dim10 <- FindClusters(sc_T_NK_dim10, resolution = 0.2)
sc_T_NK_dim20 <- FindClusters(sc_T_NK_dim20, resolution = 0.2)
sc_T_NK_dim30 <- FindClusters(sc_T_NK_dim30, resolution = 0.2)
sc_T_NK_dim10 <- FindClusters(sc_T_NK_dim10, resolution = 0.3)
sc_T_NK_dim20 <- FindClusters(sc_T_NK_dim20, resolution = 0.3)
sc_T_NK_dim30 <- FindClusters(sc_T_NK_dim30, resolution = 0.3)
sc_T_NK_dim10 <- FindClusters(sc_T_NK_dim10, resolution = 0.4)
sc_T_NK_dim20 <- FindClusters(sc_T_NK_dim20, resolution = 0.4)
sc_T_NK_dim30 <- FindClusters(sc_T_NK_dim30, resolution = 0.4)
sc_T_NK_dim10 <- FindClusters(sc_T_NK_dim10, resolution = 0.5)
sc_T_NK_dim20 <- FindClusters(sc_T_NK_dim20, resolution = 0.5)
sc_T_NK_dim30 <- FindClusters(sc_T_NK_dim30, resolution = 0.5)
sc_T_NK_dim10 <- FindClusters(sc_T_NK_dim10, resolution = 0.6)
sc_T_NK_dim20 <- FindClusters(sc_T_NK_dim20, resolution = 0.6)
sc_T_NK_dim30 <- FindClusters(sc_T_NK_dim30, resolution = 0.6)

DimPlot(sc_T_NK_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T_NK_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_T_NK_dim30, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save
FeaturePlot(sc_T_NK_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_T_NK_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_T_NK_dim30, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

FeaturePlot(sc_T_NK_dim10, raster = T, features = c("CD3E","CD4","IFNG","IL4","IL9","IL17A","IL22","CXCR5","PDCD1","IL21","CXCL13","FOXP3","CD8A","GZMB","NCAM1","FCGR3A"), pt.size = 3)
FeaturePlot(sc_T_NK_dim20, raster = T, features = c("CD3E","CD4","IFNG","IL4","IL9","IL17A","IL22","CXCR5","PDCD1","IL21","CXCL13","FOXP3","CD8A","GZMB","NCAM1","FCGR3A"), pt.size = 3)
FeaturePlot(sc_T_NK_dim30, raster = T, features = c("CD3E","CD4","IFNG","IL4","IL9","IL17A","IL22","CXCR5","PDCD1","IL21","CXCL13","FOXP3","CD8A","GZMB","NCAM1","FCGR3A"), pt.size = 3)

sc_T_NK <- sc_T_NK_dim10

DimPlot(sc_T_NK, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T_NK, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T_NK, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_T_NK, group.by = "SCT_snn_res.0.4", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T_NK, group.by = "SCT_snn_res.0.5", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T_NK, group.by = "SCT_snn_res.0.6", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

FeaturePlot(sc_T_NK, raster = T, features = c("CD3E","CD4","IFNG","IL4","IL9","IL17A","IL22","CXCR5","PDCD1","IL21","CXCL13","FOXP3","CD8A","PRF1","GZMB","NCAM1","FCGR3A","CCR7","SELL","CD69","IL2RA","ITGAE"), pt.size = 3) # save

sc_T_NK <- SetIdent(sc_T_NK, value = "SCT_snn_res.0.5")

DotPlot(sc_T_NK, features = c("CD3E","CD4","IFNG","IL4","IL21","CXCR5","PDCD1","CXCL13","FOXP3","CD8A","PRF1","GZMB","NCAM1","FCGR3A","CCR7","SELL","CD69","IL2RA","ITGAE")) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_T_NK <- PrepSCTFindMarkers(sc_T_NK)
all.markers <- FindAllMarkers(sc_T_NK, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc_T_NK) <- "RNA"
sc_T_NK[["RNA"]] <- split(sc_T_NK[["RNA"]], f = sc_T_NK$orig.ident)
sc_T_NK <- NormalizeData(sc_T_NK)
sc_T_NK <- FindVariableFeatures(sc_T_NK)
sc_T_NK <- ScaleData(sc_T_NK)
DefaultAssay(sc_T_NK) <- "SCT"
sc_T_NK <- JoinLayers(sc_T_NK, assay = "RNA")

DoHeatmap(sc_T_NK, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_T_NK, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Lymphoid/T_NK/all.markers_T_NK_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Lymphoid/T_NK/top10.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Lymphoid/T_NK/top30.txt", sep = "\t")

sc_T_NK <- RenameIdents(sc_T_NK, "0"="CD4 T")
sc_T_NK <- RenameIdents(sc_T_NK, "1"="CD8 T")
sc_T_NK <- RenameIdents(sc_T_NK, "2"="CD8 T")
sc_T_NK <- RenameIdents(sc_T_NK, "3"="CD4 T")
sc_T_NK <- RenameIdents(sc_T_NK, "4"="CD8 T")
sc_T_NK <- RenameIdents(sc_T_NK, "5"="NK cell")
sc_T_NK <- RenameIdents(sc_T_NK, "6"="NK cell")

DimPlot(sc_T_NK)


sc_NK <- subset(sc_T_NK, idents = "NK cell")

DefaultAssay(sc_NK) <- "RNA"
sc_NK <- DietSeurat(sc_NK, assays = "RNA")

sc_NK <- SCTransform(sc_NK, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_NK <- RunPCA(sc_NK)

sc_NK_dim10 <- RunUMAP(sc_NK, dims = 1:10)
sc_NK_dim20 <- RunUMAP(sc_NK, dims = 1:20)

sc_NK_dim10 <- FindNeighbors(sc_NK_dim10, dims = 1:10)
sc_NK_dim20 <- FindNeighbors(sc_NK_dim20, dims = 1:20)

sc_NK_dim10 <- FindClusters(sc_NK_dim10, resolution = 0.1)
sc_NK_dim20 <- FindClusters(sc_NK_dim20, resolution = 0.1)
sc_NK_dim10 <- FindClusters(sc_NK_dim10, resolution = 0.2)
sc_NK_dim20 <- FindClusters(sc_NK_dim20, resolution = 0.2)
sc_NK_dim10 <- FindClusters(sc_NK_dim10, resolution = 0.3)
sc_NK_dim20 <- FindClusters(sc_NK_dim20, resolution = 0.3)

DimPlot(sc_NK_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) + 
    DimPlot(sc_NK_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) # save
FeaturePlot(sc_NK_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_NK_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_NK_dim10, raster = T, features = c("CD52","KLRC2","CCL5","IL32","VIM","GZMK","SELL","EIF3G","CD44","GAS5","XCL2","TPT1","EEF1A1","IL7R","CXCR4","IER2","KLRB1","ACTB","ACTG1","CORO1A","RNF213","NCL","NEAT1","C1orf56","TXNIP"), pt.size = 5)
FeaturePlot(sc_NK_dim20, raster = T, features = c("CD52","KLRC2","CCL5","IL32","VIM","GZMK","SELL","EIF3G","CD44","GAS5","XCL2","TPT1","EEF1A1","IL7R","CXCR4","IER2","KLRB1","ACTB","ACTG1","CORO1A","RNF213","NCL","NEAT1","C1orf56","TXNIP"), pt.size = 5)

sc_NK <- sc_NK_dim10

DimPlot(sc_NK, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) + 
    DimPlot(sc_NK, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) + 
    DimPlot(sc_NK, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 3) # save

FeaturePlot(sc_NK_dim10, raster = T, features = c("CD52","KLRC2","CCL5","IL32","VIM","GZMK","SELL","EIF3G","CD44","GAS5","XCL2","TPT1","EEF1A1","IL7R","CXCR4","IER2","KLRB1","ACTB","ACTG1","CORO1A","RNF213","NCL","NEAT1","C1orf56","TXNIP"), pt.size = 5)

sc_NK <- SetIdent(sc_NK, value = "SCT_snn_res.0.2")

DotPlot(sc_NK, features = c("CD3E","CD4","IFNG","IL21","CXCR5","PDCD1","CD8A","PRF1","GZMB","NCAM1","FCGR3A","CCR7","SELL","CD69","IL2RA","ITGAE"), cluster.idents = TRUE) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred")

DotPlot(sc_NK, features = c("CD52","KLRC2","CCL5","IL32","VIM","GZMK","SELL","EIF3G","CD44","GAS5","XCL2","TPT1","EEF1A1","IL7R","CXCR4","IER2","KLRB1","ACTB","ACTG1","CORO1A","RNF213","NCL","NEAT1","C1orf56","TXNIP")) +  # Genes reported in "https://www.nature.com/articles/s41590-024-01883-0/figures/2".
RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save


sc_NK <- PrepSCTFindMarkers(sc_NK)
all.markers <- FindAllMarkers(sc_NK, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc_NK) <- "RNA"
sc_NK <- NormalizeData(sc_NK)
sc_NK <- FindVariableFeatures(sc_NK)
sc_NK <- ScaleData(sc_NK)
DefaultAssay(sc_NK) <- "SCT"

DoHeatmap(sc_NK, features = all.markers_top10$gene, size = 3.5) + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_NK, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Lymphoid/T_NK/NK/all.markers_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Lymphoid/T_NK/NK/top10.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Lymphoid/T_NK/NK/top30.txt", sep = "\t")

sc_NK <- RenameIdents(sc_NK, "0"="NK cell-3")
sc_NK <- RenameIdents(sc_NK, "1"="NK cell-1")
sc_NK <- RenameIdents(sc_NK, "2"="NK cell-2")


sc_T <- subset(sc_T_NK, idents = c("CD4 T", "CD8 T"))

DefaultAssay(sc_T) <- "RNA"
sc_T <- DietSeurat(sc_T, assays = "RNA")
sc_T[["RNA"]] <- split(sc_T[["RNA"]], f = sc_T$orig.ident)

sc_T <- SCTransform(sc_T, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_T <- JoinLayers(sc_T, assay = "RNA")

sc_T <- RunPCA(sc_T)

sc_T_dim10 <- RunUMAP(sc_T, dims = 1:10)
sc_T_dim20 <- RunUMAP(sc_T, dims = 1:20)
sc_T_dim30 <- RunUMAP(sc_T, dims = 1:30)

sc_T_dim10 <- FindNeighbors(sc_T_dim10, dims = 1:10)
sc_T_dim20 <- FindNeighbors(sc_T_dim20, dims = 1:20)
sc_T_dim30 <- FindNeighbors(sc_T_dim30, dims = 1:30)

sc_T_dim10 <- FindClusters(sc_T_dim10, resolution = 0.1)
sc_T_dim20 <- FindClusters(sc_T_dim20, resolution = 0.1)
sc_T_dim30 <- FindClusters(sc_T_dim30, resolution = 0.1)
sc_T_dim10 <- FindClusters(sc_T_dim10, resolution = 0.2)
sc_T_dim20 <- FindClusters(sc_T_dim20, resolution = 0.2)
sc_T_dim30 <- FindClusters(sc_T_dim30, resolution = 0.2)
sc_T_dim10 <- FindClusters(sc_T_dim10, resolution = 0.3)
sc_T_dim20 <- FindClusters(sc_T_dim20, resolution = 0.3)
sc_T_dim30 <- FindClusters(sc_T_dim30, resolution = 0.3)
sc_T_dim10 <- FindClusters(sc_T_dim10, resolution = 0.4)
sc_T_dim20 <- FindClusters(sc_T_dim20, resolution = 0.4)
sc_T_dim30 <- FindClusters(sc_T_dim30, resolution = 0.4)
sc_T_dim10 <- FindClusters(sc_T_dim10, resolution = 0.5)
sc_T_dim20 <- FindClusters(sc_T_dim20, resolution = 0.5)
sc_T_dim30 <- FindClusters(sc_T_dim30, resolution = 0.5)
sc_T_dim10 <- FindClusters(sc_T_dim10, resolution = 0.6)
sc_T_dim20 <- FindClusters(sc_T_dim20, resolution = 0.6)
sc_T_dim30 <- FindClusters(sc_T_dim30, resolution = 0.6)

DimPlot(sc_T_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_T_dim30, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2)  # save
FeaturePlot(sc_T_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_T_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_T_dim30, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

FeaturePlot(sc_T_dim10, raster = T, features = c("CD4","ANXA1","GATA3","PLP2","CRIP2","TIMP1","S100A10","IL7R","AQP3","TCF7","KLRB1","RORA","TSHZ2","IFI44L","ITGB1","ANK3","ZBTB16","TXK","PCNX1","IL17A","CCL4","HLA-DRB1","CD69","CCR7"), pt.size = 3)
FeaturePlot(sc_T_dim20, raster = T, features = c("CD4","ANXA1","GATA3","PLP2","CRIP2","TIMP1","S100A10","IL7R","AQP3","TCF7","KLRB1","RORA","TSHZ2","IFI44L","ITGB1","ANK3","ZBTB16","TXK","PCNX1","IL17A","CCL4","HLA-DRB1","CD69","CCR7"), pt.size = 3)
FeaturePlot(sc_T_dim30, raster = T, features = c("CD4","ANXA1","GATA3","PLP2","CRIP2","TIMP1","S100A10","IL7R","AQP3","TCF7","KLRB1","RORA","TSHZ2","IFI44L","ITGB1","ANK3","ZBTB16","TXK","PCNX1","IL17A","CCL4","HLA-DRB1","CD69","CCR7"), pt.size = 3)

FeaturePlot(sc_T_dim10, raster = T, features = c("CD8A","CCR7","SELL","LEF1","TCF7","IL7R","ANXA1","CD69","ITGAE","GZMK","GZMB","PRF1","GNLY","CX3CR1","CCL5","TBX21","ZNF683","IKZF2","IKZF3","KLRC2"), pt.size = 3)
FeaturePlot(sc_T_dim20, raster = T, features = c("CD8A","CCR7","SELL","LEF1","TCF7","IL7R","ANXA1","CD69","ITGAE","GZMK","GZMB","PRF1","GNLY","CX3CR1","CCL5","TBX21","ZNF683","IKZF2","IKZF3","KLRC2"), pt.size = 3)
FeaturePlot(sc_T_dim30, raster = T, features = c("CD8A","CCR7","SELL","LEF1","TCF7","IL7R","ANXA1","CD69","ITGAE","GZMK","GZMB","PRF1","GNLY","CX3CR1","CCL5","TBX21","ZNF683","IKZF2","IKZF3","KLRC2"), pt.size = 3)

sc_T <- sc_T_dim20

DimPlot(sc_T, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_T, group.by = "SCT_snn_res.0.4", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T, group.by = "SCT_snn_res.0.5", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_T, group.by = "SCT_snn_res.0.6", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

FeaturePlot(sc_T, raster = T, features = c("CD4","ANXA1","GATA3","PLP2","CRIP2","TIMP1","S100A10","IL7R","AQP3","TCF7","KLRB1","RORA","TSHZ2","IFI44L","ITGB1","ANK3","ZBTB16","TXK","PCNX1","IL17A","CCL4","HLA-DRB1","CD69","CCR7"), pt.size = 3)
FeaturePlot(sc_T, raster = T, features = c(,"CCR7","SELL","LEF1","TCF7","IL7R","ANXA1","CD69","ITGAE","GZMK","GZMB","PRF1","GNLY","CX3CR1","CCL5","TBX21","ZNF683","IKZF2","IKZF3","KLRC2"), pt.size = 3)

sc_T <- SetIdent(sc_T, value = "SCT_snn_res.0.5")

DotPlot(sc_T, features = c("CD4","ANXA1","GATA3","PLP2","CRIP2","TIMP1","S100A10","IL7R","AQP3","TCF7","KLRB1","RORA","TSHZ2","IFI44L","ITGB1","ANK3","ZBTB16","TXK","PCNX1","IL17A","CCL4","HLA-DRB1","CD69","ITGAE","CCR7","SELL","LEF1","CD8A","GZMK","GZMB","PRF1","GNLY","CX3CR1","CCL5","TBX21","ZNF683","IKZF2","IKZF3","KLRC2")) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_T <- PrepSCTFindMarkers(sc_T)
all.markers <- FindAllMarkers(sc_T, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc_T) <- "RNA"
sc_T[["RNA"]] <- split(sc_T[["RNA"]], f = sc_T$orig.ident)
sc_T <- NormalizeData(sc_T)
sc_T <- FindVariableFeatures(sc_T)
sc_T <- ScaleData(sc_T)
DefaultAssay(sc_T) <- "SCT"
sc_T <- JoinLayers(sc_T, assay = "RNA")

DoHeatmap(sc_T, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_T, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Lymphoid/T_NK/T/all.markers_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Lymphoid/T_NK/T/top10.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Lymphoid/T_NK/T/top30.txt", sep = "\t")

sc_T <- RenameIdents(sc_T, "0"="CD4 T-Naive")
sc_T <- RenameIdents(sc_T, "1"="CD8 T-ResidentMemory")
sc_T <- RenameIdents(sc_T, "2"="NKT cell")
sc_T <- RenameIdents(sc_T, "3"="CD8 T-CCL4+")
sc_T <- RenameIdents(sc_T, "4"="CD4 T-Activated")
sc_T <- RenameIdents(sc_T, "5"="CD4 T-Effector")


sc_lymphoid$CellType_tmp <- NA_character_
sc_lymphoid$CellType_tmp[sc_lymphoid@active.ident == "B cell"] <- "B cell"
sc_lymphoid$CellType_tmp[sc_lymphoid@active.ident == "Plasma"] <- "Plasma cell"
sc_lymphoid$CellType_tmp[sc_lymphoid@active.ident == "pDC"] <- "pDC"
sc_lymphoid$CellType_tmp[sc_lymphoid@active.ident == "Treg"] <- "Treg"
sc_lymphoid$CellType_tmp[sc_lymphoid@active.ident == "Doublet"] <- "Doublet"
NK_map <- setNames(as.character(Idents(sc_NK)), colnames(sc_NK))
T_map <- setNames(as.character(Idents(sc_T)), colnames(sc_T))
sc_lymphoid$CellType_tmp[names(NK_map)] <- NK_map
sc_lymphoid$CellType_tmp[names(T_map)] <- T_map

table(sc_lymphoid$CellType_tmp)
#               B cell      CD4 T-Activated       CD4 T-Effector          CD4 T-Naive          CD8 T-CCL4+ CD8 T-ResidentMemory 
#                  377                  362                  312                 1346                  560                  776 
#              Doublet            NK cell-1            NK cell-2            NK cell-3             NKT cell                  pDC 
#                  499                  304                   78                  325                  756                   51 
#          Plasma cell                 Treg 
#                  614                   51 


# saveRDS(sc_lymphoid, "02_Publicdata/ver6/Lymphoid/seuratObj_Lymphoid_20251022.rds")
# saveRDS(sc_T_NK, "02_Publicdata/ver6/Lymphoid/seuratObj_T_NK_20251022.rds")
# saveRDS(sc_NK, "02_Publicdata/ver6/Lymphoid/NK/seuratObj_NK_20251022.rds")
# saveRDS(sc_T, "02_Publicdata/ver6/Lymphoid/T/seuratObj_Tcell_20251022_latest.rds")
# saveRDS(sc_lymphoid, "02_Publicdata/ver6/tmp/seuratObj_lymphoid_20251022_update.rds")


# Endothelial cell
sc_ec <- readRDS("02_Publicdata/ver6/tmp/seuratObj_endothelial_20251022.rds")

sc_ec 
# An object of class Seurat 
# 62707 features across 10176 samples within 2 assays 
# Active assay: SCT (31224 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
#  3 dimensional reductions calculated: pca, integrated.rpca, umap_rpca

DefaultAssay(sc_ec) <- "RNA"
sc_ec <- DietSeurat(sc_ec, assays = "RNA")
sc_ec[["RNA"]] <- split(sc_ec[["RNA"]], f = sc_ec$orig.ident)

sc_ec <- SCTransform(sc_ec, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_ec <- JoinLayers(sc_ec, assay = "RNA")

sc_ec <- RunPCA(sc_ec)

sc_ec_dim10 <- RunUMAP(sc_ec, dims = 1:10)
sc_ec_dim20 <- RunUMAP(sc_ec, dims = 1:20)

sc_ec_dim10 <- FindNeighbors(sc_ec_dim10, dims = 1:10)
sc_ec_dim20 <- FindNeighbors(sc_ec_dim20, dims = 1:20)

sc_ec_dim10 <- FindClusters(sc_ec_dim10, resolution = 0.1)
sc_ec_dim20 <- FindClusters(sc_ec_dim20, resolution = 0.1)
sc_ec_dim10 <- FindClusters(sc_ec_dim10, resolution = 0.2)
sc_ec_dim20 <- FindClusters(sc_ec_dim20, resolution = 0.2)
sc_ec_dim10 <- FindClusters(sc_ec_dim10, resolution = 0.3)
sc_ec_dim20 <- FindClusters(sc_ec_dim20, resolution = 0.3)

DimPlot(sc_ec_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_ec_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save
DimPlot(sc_ec_dim10, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_ec_dim20, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2)
FeaturePlot(sc_ec_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_ec_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

sc_ec <- sc_ec_dim20

DimPlot(sc_ec, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_ec, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_ec, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

sc_ec <- SetIdent(sc_ec, value = "SCT_snn_res.0.3")

DimPlot(sc_ec, split.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2)
FeaturePlot(sc_ec, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3) # save
FeaturePlot(sc_ec, raster = T, features = c("PECAM1","CDH5","VWF","SEMA3G","SOX17","RGCC","ACKR1","NR2F2","CCL21", "NOTCH4","SPARC","LIFR","ICAM1"), pt.size = 3) # save

sc_ec <- PrepSCTFindMarkers(sc_ec)
all.markers <- FindAllMarkers(sc_ec, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DotPlot(sc_ec, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_ec <- RenameIdents(sc_ec, "0"="ec")
sc_ec <- RenameIdents(sc_ec, "1"="ec")
sc_ec <- RenameIdents(sc_ec, "2"="ec")
sc_ec <- RenameIdents(sc_ec, "3"="ec")
sc_ec <- RenameIdents(sc_ec, "4"="ec")
sc_ec <- RenameIdents(sc_ec, "5"="Doublet")
sc_ec <- RenameIdents(sc_ec, "6"="ec")
sc_ec <- RenameIdents(sc_ec, "7"="ec")
sc_ec <- RenameIdents(sc_ec, "8"="Doublet")
sc_ec <- RenameIdents(sc_ec, "9"="Doublet")

table(sc_ec@active.ident)
# Doublet      ec 
#     583    9593 


sc_ec_2 <- subset(sc_ec, idents = "ec")

DefaultAssay(sc_ec_2) <- "RNA"
sc_ec_2 <- DietSeurat(sc_ec_2, assays = "RNA")
sc_ec_2[["RNA"]] <- split(sc_ec_2[["RNA"]], f = sc_ec_2$orig.ident)

sc_ec_2 <- SCTransform(sc_ec_2, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_ec_2 <- JoinLayers(sc_ec_2, assay = "RNA")

sc_ec_2 <- RunPCA(sc_ec_2)

sc_ec_2_dim10 <- RunUMAP(sc_ec_2, dims = 1:10)
sc_ec_2_dim20 <- RunUMAP(sc_ec_2, dims = 1:20)

sc_ec_2_dim10 <- FindNeighbors(sc_ec_2_dim10, dims = 1:10)
sc_ec_2_dim20 <- FindNeighbors(sc_ec_2_dim20, dims = 1:20)

sc_ec_2_dim10 <- FindClusters(sc_ec_2_dim10, resolution = 0.1)
sc_ec_2_dim20 <- FindClusters(sc_ec_2_dim20, resolution = 0.1)
sc_ec_2_dim10 <- FindClusters(sc_ec_2_dim10, resolution = 0.2)
sc_ec_2_dim20 <- FindClusters(sc_ec_2_dim20, resolution = 0.2)
sc_ec_2_dim10 <- FindClusters(sc_ec_2_dim10, resolution = 0.3)
sc_ec_2_dim20 <- FindClusters(sc_ec_2_dim20, resolution = 0.3)

DimPlot(sc_ec_2_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_ec_2_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save
FeaturePlot(sc_ec_2_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_ec_2_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

FeaturePlot(sc_ec_2_dim10, raster = T, features = c("PECAM1","CDH5","VWF","SEMA3G","SOX17","RGCC","ACKR1","NR2F2","CCL21", "NOTCH4","SPARC","LIFR","ICAM1"), pt.size = 3)
FeaturePlot(sc_ec_2_dim20, raster = T, features = c("PECAM1","CDH5","VWF","SEMA3G","SOX17","RGCC","ACKR1","NR2F2","CCL21", "NOTCH4","SPARC","LIFR","ICAM1"), pt.size = 3)

sc_ec_2 <- sc_ec_2_dim20

DimPlot(sc_ec_2, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_ec_2, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_ec_2, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

FeaturePlot(sc_ec_2, raster = T, features = c("PECAM1","CDH5","VWF","SEMA3G","SOX17","RGCC","ACKR1","NR2F2","CCL21", "NOTCH4","SPARC","LIFR","ICAM1"), pt.size = 3) # save

sc_ec_2 <- SetIdent(sc_ec_2, value = "SCT_snn_res.0.2")

DotPlot(sc_ec_2, features = c("CCL21","ACKR1","NR2F2","SEMA3G","NOTCH4","RGCC")) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_ec_2 <- PrepSCTFindMarkers(sc_ec_2)
all.markers <- FindAllMarkers(sc_ec_2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc_ec_2) <- "RNA"
sc_ec_2[["RNA"]] <- split(sc_ec_2[["RNA"]], f = sc_ec_2$orig.ident)
sc_ec_2 <- NormalizeData(sc_ec_2)
sc_ec_2 <- FindVariableFeatures(sc_ec_2)
sc_ec_2 <- ScaleData(sc_ec_2)
DefaultAssay(sc_ec_2) <- "SCT"
sc_ec_2 <- JoinLayers(sc_ec_2, assay = "RNA")

DoHeatmap(sc_ec_2, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_ec_2, features = all.markers_top10$gene) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Endothelial/all.markers_endothelial_remove_doublet_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Endothelial/top10.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Endothelial/top30.txt", sep = "\t")

sc_ec_2 <- RenameIdents(sc_ec_2, "0"="Venous EC-1")
sc_ec_2 <- RenameIdents(sc_ec_2, "1"="Capillary EC-1")
sc_ec_2 <- RenameIdents(sc_ec_2, "2"="Venous EC-2")
sc_ec_2 <- RenameIdents(sc_ec_2, "3"="Arterial EC")
sc_ec_2 <- RenameIdents(sc_ec_2, "4"="Capillary EC-2")
sc_ec_2 <- RenameIdents(sc_ec_2, "5"="Capillary EC-1")
sc_ec_2 <- RenameIdents(sc_ec_2, "6"="Lymphatic EC")

all.markers <- FindAllMarkers(sc_ec_2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DoHeatmap(sc_ec_2, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_ec_2, features = all.markers_top10$gene) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

# saveRDS(all.markers, "02_Publicdata/ver6/Endothelial/all.markers_endothelial_remove_doublet_update_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Endothelial/top10_update.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Endothelial/top30_update.txt", sep = "\t")


sc_ec$CellType_tmp <- NA_character_
sc_ec$CellType_tmp[colnames(sc_ec)[Idents(sc_ec) == "Doublet"]] <- "Doublet"
ec_map <- setNames(as.character(Idents(sc_ec_2)), colnames(sc_ec_2))
sc_ec$CellType_tmp[names(ec_map)] <- ec_map

table(sc_ec$CellType_tmp)
#    Arterial EC Capillary EC-1 Capillary EC-2        Doublet   Lymphatic EC    Venous EC-1    Venous EC-2 
#            672           1874            498            583            217           5496            836

# saveRDS(sc_ec_2, "02_Publicdata/ver6/Endothelial/seuratObj_endothelial_remove_doublet_20251022.rds")
# saveRDS(sc_ec, "02_Publicdata/ver6/tmp/seuratObj_endothelial_20251022_update.rds")


# Mural cell
sc_mural <- readRDS("02_Publicdata/ver6/tmp/seuratObj_mural_20251022.rds")

sc_mural 
# An object of class Seurat 
# 62707 features across 7041 samples within 2 assays 
# Active assay: SCT (31224 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
#  3 dimensional reductions calculated: pca, integrated.rpca, umap.rpca

DefaultAssay(sc_mural) <- "RNA"
sc_mural <- DietSeurat(sc_mural, assays = "RNA")
sc_mural[["RNA"]] <- split(sc_mural[["RNA"]], f = sc_mural$orig.ident)

sc_mural <- SCTransform(sc_mural, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_mural <- JoinLayers(sc_mural, assay = "RNA")

sc_mural <- RunPCA(sc_mural)

sc_mural_dim10 <- RunUMAP(sc_mural, dims = 1:10)
sc_mural_dim20 <- RunUMAP(sc_mural, dims = 1:20)

sc_mural_dim10 <- FindNeighbors(sc_mural_dim10, dims = 1:10)
sc_mural_dim20 <- FindNeighbors(sc_mural_dim20, dims = 1:20)

sc_mural_dim10 <- FindClusters(sc_mural_dim10, resolution = 0.1)
sc_mural_dim20 <- FindClusters(sc_mural_dim20, resolution = 0.1)
sc_mural_dim10 <- FindClusters(sc_mural_dim10, resolution = 0.2)
sc_mural_dim20 <- FindClusters(sc_mural_dim20, resolution = 0.2)
sc_mural_dim10 <- FindClusters(sc_mural_dim10, resolution = 0.3)
sc_mural_dim20 <- FindClusters(sc_mural_dim20, resolution = 0.3)

DimPlot(sc_mural_dim10, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_mural_dim20, shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save
DimPlot(sc_mural_dim10, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_mural_dim20, group.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2)
FeaturePlot(sc_mural_dim10, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)
FeaturePlot(sc_mural_dim20, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

sc_mural <- sc_mural_dim10

DimPlot(sc_mural, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) +
    DimPlot(sc_mural, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_mural, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

sc_mural <- SetIdent(sc_mural, value = "SCT_snn_res.0.3")

DimPlot(sc_mural, split.by = "orig.ident", shuffle = TRUE, raster = T, label = T, repel = T, pt.size = 2)
FeaturePlot(sc_mural, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3) # save
FeaturePlot(sc_mural, raster = T, features = c("PDGFRB","RGS5","ABCC9","KCNJ8","AGT","ACTA2","MYH11","RERGL","CASQ2","KCNAB1","HMCN2","FLNC"), pt.size = 3) # save

sc_mural <- PrepSCTFindMarkers(sc_mural)
all.markers <- FindAllMarkers(sc_mural, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DotPlot(sc_mural, features = unique(all.markers_top10$gene)) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_mural <- RenameIdents(sc_mural, "0"="mural")
sc_mural <- RenameIdents(sc_mural, "1"="mural")
sc_mural <- RenameIdents(sc_mural, "2"="mural")
sc_mural <- RenameIdents(sc_mural, "3"="Doublet")
sc_mural <- RenameIdents(sc_mural, "4"="mural")
sc_mural <- RenameIdents(sc_mural, "5"="Doublet")


sc_mural_2 <- subset(sc_mural, idents = "mural")

DefaultAssay(sc_mural_2) <- "RNA"
sc_mural_2 <- DietSeurat(sc_mural_2, assays = "RNA")
sc_mural_2[["RNA"]] <- split(sc_mural_2[["RNA"]], f = sc_mural_2$orig.ident)

sc_mural_2 <- SCTransform(sc_mural_2, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)

sc_mural_2 <- JoinLayers(sc_mural_2, assay = "RNA")

sc_mural_2 <- RunPCA(sc_mural_2)

sc_mural_2 <- RunUMAP(sc_mural_2, dims = 1:10)

sc_mural_2 <- FindNeighbors(sc_mural_2, dims = 1:10)

sc_mural_2 <- FindClusters(sc_mural_2, resolution = 0.1)
sc_mural_2 <- FindClusters(sc_mural_2, resolution = 0.2)
sc_mural_2 <- FindClusters(sc_mural_2, resolution = 0.3)

FeaturePlot(sc_mural_2, raster = T, features = c("PDGFRA","PDGFRB","PECAM1","PTPRC","CD68","CD3D","MS4A1"), pt.size = 3)

DimPlot(sc_mural_2, group.by = "SCT_snn_res.0.1", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_mural_2, group.by = "SCT_snn_res.0.2", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) + 
    DimPlot(sc_mural_2, group.by = "SCT_snn_res.0.3", shuffle = TRUE, raster = TRUE, label = T, repel = T, pt.size = 2) # save

FeaturePlot(sc_mural_2, raster = T, features = c("PDGFRB","RGS5","ABCC9","KCNJ8","AGT","ACTA2","MYH11","RERGL","CASQ2","KCNAB1","HMCN2","FLNC"), pt.size = 3) # save

sc_mural_2 <- SetIdent(sc_mural_2, value = "SCT_snn_res.0.2")

DotPlot(sc_mural_2, features = c("PDGFRB","RGS5","ABCC9","KCNJ8","AGT","ACTA2","MYH11","RERGL","CASQ2","KCNAB1","HMCN2","FLNC"), cluster.idents = TRUE) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_mural_2 <- PrepSCTFindMarkers(sc_mural_2)
all.markers <- FindAllMarkers(sc_mural_2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
all.markers_top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
all.markers_top30 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

DefaultAssay(sc_mural_2) <- "RNA"
sc_mural_2[["RNA"]] <- split(sc_mural_2[["RNA"]], f = sc_mural_2$orig.ident)
sc_mural_2 <- NormalizeData(sc_mural_2)
sc_mural_2 <- FindVariableFeatures(sc_mural_2)
sc_mural_2 <- ScaleData(sc_mural_2)
DefaultAssay(sc_mural_2) <- "SCT"
sc_mural_2 <- JoinLayers(sc_mural_2, assay = "RNA")

DoHeatmap(sc_mural_2, features = all.markers_top10$gene, size = 3.5, assay = "RNA", slot = "scale.data") + scale_fill_gradientn(colors = colors) # save
DotPlot(sc_mural_2, features = all.markers_top10$gene, cluster.idents = TRUE) + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "darkred") # save

sc_mural_2 <- RenameIdents(sc_mural_2, "0"="Pericyte")
sc_mural_2 <- RenameIdents(sc_mural_2, "1"="VSMC-1")
sc_mural_2 <- RenameIdents(sc_mural_2, "2"="VSMC-2")

# saveRDS(all.markers, "02_Publicdata/ver6/Mural/all.markers_mural_remove_doublet_20251022.rds")
# write.table(all.markers_top10, "02_Publicdata/ver6/Mural/top10.txt", sep = "\t")
# write.table(all.markers_top30, "02_Publicdata/ver6/Mural/top30.txt", sep = "\t")


sc_mural$CellType_tmp <- NA_character_
sc_mural$CellType_tmp[colnames(sc_mural)[Idents(sc_mural) == "Doublet"]] <- "Doublet"
mural_map <- setNames(as.character(Idents(sc_mural_2)), colnames(sc_mural_2))
sc_mural$CellType_tmp[names(mural_map)] <- mural_map

table(sc_mural$CellType_tmp)
#  Doublet Pericyte   VSMC-1   VSMC-2 
#     1034     2695     2359      953 

# saveRDS(sc_mural_2, "02_Publicdata/ver6/Mural/seuratObj_mural_remove_doublet_20251022.rds")
# saveRDS(sc_mural, "02_Publicdata/ver6/tmp/seuratObj_mural_20251022_update.rds")



# Integration of annotations and downstream analyses *****
sc <- readRDS("02_Publicdata/ver6/tmp/seuratObj_after_rename_20251022.rds")
sc_stromal <- readRDS("02_Publicdata/ver6/tmp/seuratObj_stromal_20251022_update.rds")
sc_myeloid <- readRDS("02_Publicdata/ver6/tmp/seuratObj_myeloid_20251022_update.rds")
sc_lymphoid <- readRDS("02_Publicdata/ver6/tmp/seuratObj_lymphoid_20251022_update.rds")
sc_ec <- readRDS("02_Publicdata/ver6//tmp/seuratObj_endothelial_20251022_update.rds")
sc_mural <- readRDS("02_Publicdata/ver6/tmp/seuratObj_mural_20251022_update.rds")

anno_map <- c(
  setNames(as.character(sc_stromal$CellType_tmp), colnames(sc_stromal)),
  setNames(as.character(sc_myeloid$CellType_tmp), colnames(sc_myeloid)),
  setNames(as.character(sc_lymphoid$CellType_tmp), colnames(sc_lymphoid)),
  setNames(as.character(sc_ec$CellType_tmp), colnames(sc_ec)),
  setNames(as.character(sc_mural$CellType_tmp), colnames(sc_mural))
)

sc$CellType_Class3_tmp <- anno_map[ match(colnames(sc), names(anno_map)) ]

map_CellType_tmp_Class1 <- c(
  "Lining-layer fibroblast-1" = "Stromal",
  "Lining-layer fibroblast-2" = "Stromal",
  "Lining-layer fibroblast-MMP3+" = "Stromal",
  "Sublining-layer fibroblast-APOE+CXCL12+" = "Stromal",
  "Sublining-layer fibroblast-APOD+CXCL14+" = "Stromal",
  "Sublining-layer fibroblast-MFAP5+PI16+" = "Stromal",
  "Sublining-layer fibroblast-COMP+" = "Stromal",
  "Adipocyte" = "Stromal",
  "Macrophage-TIMD4+" = "Myeloid cell",
  "Macrophage-TIMD4+APOE+" = "Myeloid cell",
  "Macrophage-LYVE1+" = "Myeloid cell",
  "Macrophage-CLEC10A+" = "Myeloid cell",
  "Macrophage-ISG15+" = "Myeloid cell",
  "Macrophage-SPP1+" = "Myeloid cell",
  "Macrophage-Chemokine high" = "Myeloid cell",
  "Monocyte-Classical" = "Myeloid cell",
  "Monocyte-Int" = "Myeloid cell",
  "Monocyte-nonClassical" = "Myeloid cell",
  "cDC1" = "Myeloid cell",
  "cDC2" = "Myeloid cell",
  "mregDC" = "Myeloid cell",
  "Mast cell" = "Myeloid cell",
  "CD4 T-Naive" = "Lymphoid cell",
  "CD4 T-Activated" = "Lymphoid cell",
  "CD4 T-Effector" = "Lymphoid cell",
  "CD8 T-CCL4+" = "Lymphoid cell",
  "CD8 T-ResidentMemory" = "Lymphoid cell",
  "Treg" = "Lymphoid cell",
  "NKT cell" = "Lymphoid cell",
  "NK cell-1" = "Lymphoid cell",
  "NK cell-2" = "Lymphoid cell",
  "NK cell-3" = "Lymphoid cell",
  "B cell" = "Lymphoid cell",
  "Plasma cell" = "Lymphoid cell",
  "pDC" = "Lymphoid cell",
  "Arterial EC" = "Endothelial cell",
  "Capillary EC-1" = "Endothelial cell",
  "Capillary EC-2" = "Endothelial cell",
  "Venous EC-1" = "Endothelial cell",
  "Venous EC-2" = "Endothelial cell",
  "Lymphatic EC" = "Endothelial cell",
  "VSMC-1" = "Mural cell",
  "VSMC-2" = "Mural cell",
  "Pericyte" = "Mural cell",
  "Doublet" = "Doublet"
)        

map_CellType_tmp_Class2 <- c(
  "Lining-layer fibroblast-1" = "Lining-layer fibroblast",
  "Lining-layer fibroblast-2" = "Lining-layer fibroblast",
  "Lining-layer fibroblast-MMP3+" = "Lining-layer fibroblast",
  "Sublining-layer fibroblast-APOE+CXCL12+" = "Sublining-layer fibroblast",
  "Sublining-layer fibroblast-APOD+CXCL14+" = "Sublining-layer fibroblast",
  "Sublining-layer fibroblast-MFAP5+PI16+" = "Sublining-layer fibroblast",
  "Sublining-layer fibroblast-COMP+" = "Sublining-layer fibroblast",
  "Adipocyte" = "Adipocyte",
  "Macrophage-TIMD4+" = "Macrophage",
  "Macrophage-TIMD4+APOE+" = "Macrophage",
  "Macrophage-LYVE1+" = "Macrophage",
  "Macrophage-CLEC10A+" = "Macrophage",
  "Macrophage-ISG15+" = "Macrophage",
  "Macrophage-SPP1+" = "Macrophage",
  "Macrophage-Chemokine high" = "Macrophage",
  "Monocyte-Classical" = "Monocyte",
  "Monocyte-Int" = "Monocyte",
  "Monocyte-nonClassical" = "Monocyte",
  "cDC1" = "Dendritic cell",
  "cDC2" = "Dendritic cell",
  "mregDC" = "Dendritic cell",
  "Mast cell" = "Granulocyte",
  "CD4 T-Naive" = "CD4 T cell",
  "CD4 T-Activated" = "CD4 T cell",
  "CD4 T-Effector" = "CD4 T cell",
  "CD8 T-CCL4+" = "CD8 T cell",
  "CD8 T-ResidentMemory" = "CD8 T cell",
  "Treg" = "Treg",
  "NKT cell" = "NKT cell",
  "NK cell-1" = "NK cell",
  "NK cell-2" = "NK cell",
  "NK cell-3" = "NK cell",
  "B cell" = "B cell",
  "Plasma cell" = "Plasma cell",
  "pDC" = "Dendritic cell",
  "Arterial EC" = "Endothelial cell",
  "Capillary EC-1" = "Endothelial cell",
  "Capillary EC-2" = "Endothelial cell",
  "Venous EC-1" = "Endothelial cell",
  "Venous EC-2" = "Endothelial cell",
  "Lymphatic EC" = "Endothelial cell",
  "VSMC-1" = "VSMC",
  "VSMC-2" = "VSMC",
  "Pericyte" = "Pericyte",
  "Doublet" = "Doublet"
)             

sc$CellType_Class1 <- unname(map_CellType_tmp_Class1[sc$CellType_Class3_tmp ])
sc$CellType_Class2 <- unname(map_CellType_tmp_Class2[sc$CellType_Class3_tmp ])
sc$CellType_Class3 <- sc$CellType_Class3_tmp

sc$CellType_Class1_tmp <- NULL
sc$CellType_Class3_tmp <- NULL

DimPlot(sc, group.by = "CellType_Class1", shuffle = TRUE, raster = TRUE, label = T, repel = T) + 
    DimPlot(sc, group.by = "CellType_Class2", shuffle = TRUE, raster = TRUE, label = T, repel = T) + 
    DimPlot(sc, group.by = "CellType_Class3", shuffle = TRUE, raster = TRUE, label = T, repel = T) # check

# saveRDS(sc, "02_Publicdata/ver6/tmp/seuratObj_add_label_latest_20251027.rds")


sc_nonDoublet <- sc
sc_nonDoublet <- SetIdent(sc_nonDoublet, value = "CellType_Class1")
sc_nonDoublet <- subset(sc_nonDoublet, idents = "Doublet", invert = TRUE)

sc_nonDoublet <- SetIdent(sc_nonDoublet, value = "CellType_Class2")
sc_nonDoublet$CellType_Class2 <- factor(sc_nonDoublet$CellType_Class2, levels = c(
  "Lining-layer fibroblast",
  "Sublining-layer fibroblast",
  "Adipocyte",
  "Endothelial cell",
  "VSMC",
  "Pericyte",
  "Macrophage", 
  "Monocyte",
  "Dendritic cell",
  "Granulocyte",
  "CD4 T cell",
  "CD8 T cell",
  "Treg",
  "NKT cell",
  "NK cell",  
  "B cell",
  "Plasma cell"
))

cluster_colors <- c(
  "Lining-layer fibroblast" = "#D73027",
  "Sublining-layer fibroblast" = "#FFA500",
  "Adipocyte" = "#FFD700",
  "Endothelial cell" = "#006400",
  "VSMC" = "#4C8B45",
  "Pericyte" = "#A1D99B",
  "Macrophage" = "#003366", 
  "Monocyte" = "#005B8E",
  "Dendritic cell" = "#4DA6B8",
  "Granulocyte" = "#8DD3C7",
  "CD4 T cell" = "#54278F",
  "CD8 T cell" = "#6A3D9A",
  "Treg" = "#8E6BBE",
  "NKT cell" = "#BFA6D9",
  "NK cell" = "#E6DFF2",
  "B cell" = "#AE017E",
  "Plasma cell" = "#FBB4B9"
)


# p_class1 <- DimPlot(sc_nonDoublet, group.by = "CellType_Class1", shuffle = TRUE, raster = TRUE) +
#   labs(title = "", x = "", y = "") +
#   theme_classic() +
#   theme(
#     plot.title = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.line = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x  = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     legend.title = element_text(size = 5),
#     legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
#     legend.key.size  = unit(2, "mm"),
#     legend.box.margin = margin(t = 0, r = 0, b = 0, l = -4, unit = "mm"),
#     plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
#   ) +
#   guides(colour = guide_legend(
#     override.aes = list(size = 1)
#   ))
# p_class2 <- DimPlot(sc_nonDoublet, group.by = "CellType_Class2", shuffle = TRUE, raster = TRUE) +
#   labs(title = "", x = "UMAP 1", y = "UMAP2") +
#   theme_classic() +
#   theme(
#     plot.title = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.line = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x  = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     legend.title = element_text(size = 5),
#     legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
#     legend.key.size  = unit(2, "mm"),
#     legend.box.margin = margin(t = 0, r = 0, b = -1, l = -4, unit = "mm"),
#     plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
#   ) +
#   guides(colour = guide_legend(
#     override.aes = list(size = 1)
#   ))
# p_class3 <- DimPlot(sc_nonDoublet, group.by = "CellType_Class3", shuffle = TRUE, raster = TRUE) +
#   labs(title = "", x = "UMAP 1", y = "UMAP2") +
#   theme_classic() +
#   theme(
#     plot.title = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.line = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x  = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     legend.title = element_text(size = 5),
#     legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
#     legend.key.size  = unit(2, "mm"),
#     legend.box.margin = margin(t = 0, r = 0, b = 0, l = -4, unit = "mm"),
#     plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
#   ) +
#   guides(colour = guide_legend(
#     override.aes = list(size = 1)
#   ))
# 
# p_UMAP <- p_class1 + p_class2 + p_class3 &
#   theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"))
# ggsave("02_Publicdata/ver6/All/after/Dim_Class_comparison_nonDoublet.png", plot = p_UMAP, height = 1.5, width = 9)
# ggsave("02_Publicdata/ver6/All/after/Dim_Class_comparison_nonDoublet.pdf", plot = p_UMAP, height = 1.5, width = 9)

p_class2_v2 <- DimPlot(sc_nonDoublet, cols = cluster_colors, group.by = "CellType_Class2", shuffle = TRUE, raster = TRUE) +
  labs(title = "", x = "UMAP 1", y = "UMAP2") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 5),
    legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
    legend.key.size  = unit(2, "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -1, l = -4, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  ) +
  guides(colour = guide_legend(
    override.aes = list(size = 1)
  ))
ggsave("99_Fig/sup_fig9/Dim_Class2.png", plot = p_class2_v2, height = 1.5, width = 2.5)
ggsave("99_Fig/sup_fig9/Dim_Class2.pdf", plot = p_class2_v2, height = 1.5, width = 2.5)

meta <- sc_nonDoublet@meta.data
df_stacked_bar_class1 <- meta %>% 
  count(orig.ident, CellType_Class1)
df_stacked_bar_class2 <- meta %>% 
  count(orig.ident, CellType_Class2)
df_stacked_bar_class3 <- meta %>% 
  count(orig.ident, CellType_Class3)

# p_stacked_bar_class1 <- ggplot(df_stacked_bar_class1, aes(x = orig.ident, y = n, fill = CellType_Class1)) +
#   geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.85) +
#   labs(y = "Fraction of cell type", fill = "Cell type") +
#   scale_y_continuous(labels = function(x) x * 100) +
#   theme_classic() +
#   theme(
#     plot.title = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.x = element_line(linewidth = 0.1),
#     axis.ticks.y = element_line(linewidth = 0.1),
#     axis.line = element_line(linewidth = 0.1),
#     axis.title.x = element_blank(),
#     axis.text.x  = element_text(, size = 5, angle = 45, hjust = 1),
#     axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
#     axis.text.y  = element_text(size = 5),
#     legend.title = element_text(size = 5, margin = margin(b = 2)),
#     legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
#     legend.key.size = unit(2, "mm"),
#     legend.margin = margin(t = -1, b = -1, unit = "mm"),
#     legend.box.margin = margin(t = 0, r = 0, b = 0, l = -4, unit = "mm"),
#     plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
#   )
# ggsave("02_Publicdata/ver6/All/after/Stacked_Class1_nonDoublet.png", p_stacked_bar_class1, height = 1.5, width = 2.5)
# ggsave("02_Publicdata/ver6/All/after/Stacked_Class1_nonDoublet.pdf", p_stacked_bar_class1, height = 1.5, width = 2.5)
# 
# p_stacked_bar_class2 <- ggplot(df_stacked_bar_class2, aes(x = orig.ident, y = n, fill = CellType_Class2)) +
#   geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.85) +
#   labs(y = "Fraction of cell type", fill = "Cell type") +
#   scale_y_continuous(labels = function(x) x * 100) +
#   theme_classic() +
#   theme(
#     plot.title = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.x = element_line(linewidth = 0.1),
#     axis.ticks.y = element_line(linewidth = 0.1),
#     axis.line = element_line(linewidth = 0.1),
#     axis.title.x = element_blank(),
#     axis.text.x  = element_text(, size = 5, angle = 45, hjust = 1),
#     axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
#     axis.text.y  = element_text(size = 5),
#     legend.title = element_text(size = 5, margin = margin(b = 2)),
#     legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
#     legend.key.size = unit(2, "mm"),
#     legend.margin = margin(t = -1, b = -1, unit = "mm"),
#     legend.box.margin = margin(t = 0, r = 0, b = -7.5, l = -4, unit = "mm"),
#     plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
#   )
# ggsave("02_Publicdata/ver6/All/after/Stacked_Class2_nonDoublet.png", p_stacked_bar_class2, height = 1.5, width = 2.8)
# ggsave("02_Publicdata/ver6/All/after/Stacked_Class2_nonDoublet.pdf", p_stacked_bar_class2, height = 1.5, width = 2.8)
# 
# p_stacked_bar_class3 <- ggplot(df_stacked_bar_class3, aes(x = orig.ident, y = n, fill = CellType_Class3)) +
#   geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.85) +
#   labs(y = "Fraction of cell type", fill = "Cell type") +
#   scale_y_continuous(labels = function(x) x * 100) +
#   theme_classic() +
#   theme(
#     plot.title = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.x = element_line(linewidth = 0.1),
#     axis.ticks.y = element_line(linewidth = 0.1),
#     axis.line = element_line(linewidth = 0.1),
#     axis.title.x = element_blank(),
#     axis.text.x  = element_text(, size = 5, angle = 45, hjust = 1),
#     axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
#     axis.text.y  = element_text(size = 5),
#     legend.title = element_text(size = 5, margin = margin(b = 2)),
#     legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
#     legend.key.size = unit(2, "mm"),
#     legend.margin = margin(t = -1, b = -1, unit = "mm"),
#     legend.box.margin = margin(t = 0, r = 0, b = -6, l = -4, unit = "mm"),
#     plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
#   )
# ggsave("02_Publicdata/ver6/All/after/Stacked_Class3_nonDoublet.png", p_stacked_bar_class3, height = 1.5, width = 5.55)
# ggsave("02_Publicdata/ver6/All/after/Stacked_Class3_nonDoublet.pdf", p_stacked_bar_class3, height = 1.5, width = 5.55)

p_stacked_bar_class2_v2 <- ggplot(df_stacked_bar_class2, aes(x = orig.ident, y = n, fill = CellType_Class2)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.85) +
  labs(y = "Fraction of cell type", fill = "Cell type") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(, size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
    axis.text.y  = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 2)),
    legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -7.5, l = -4, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )
ggsave("99_Fig/sup_fig9/Stacked_Class2.png", p_stacked_bar_class2_v2, height = 1.5, width = 2.8)
ggsave("99_Fig/sup_fig9/Stacked_Class2.pdf", p_stacked_bar_class2_v2, height = 1.5, width = 2.8)


all.markers <- FindAllMarkers(sc_nonDoublet, assay = "SCT", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, recorrect_umi = FALSE)
all.markers_top5 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

p_Dot <- DotPlot(sc_nonDoublet, features = unique(all.markers_top5$gene), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 5),
    axis.text.x  = element_text(size = 4, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 5),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig9/Dot_DEG_Class2.png", p_Dot, height = 2, width = 8)
ggsave("99_Fig/sup_fig9/Dot_DEG_Class2.pdf", p_Dot, height = 2, width = 8)

