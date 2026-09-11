library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

options(future.globals.maxSize = 80 * 1024^3)

sc <- readRDS("02_Publicdata/RDS/seuratObj_add_label_latest_20251027.rds")
sc <- SetIdent(sc, value = "CellType_Class1")


sc_stromal <- subset(sc, idents = "Stromal")

DefaultAssay(sc_stromal) <- "RNA"
sc_stromal <- DietSeurat(sc_stromal, assays = "RNA")
sc_stromal[["RNA"]] <- split(sc_stromal[["RNA"]], f = sc_stromal$orig.ident)
sc_stromal <- SCTransform(sc_stromal, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)
sc_stromal <- RunPCA(sc_stromal)
sc_stromal <- IntegrateLayers(object = sc_stromal,
                             method = RPCAIntegration,
                             normalization.method = "SCT",
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             verbose = FALSE
                             )
sc_stromal <- JoinLayers(sc_stromal, assay = "RNA")
sc_stromal <- RunUMAP(sc_stromal, dims = 1:10, reduction="integrated.rpca")

sc_stromal$CellType_Class3 <- factor(sc_stromal$CellType_Class3, levels = c(
  "Lining-layer fibroblast-1",
  "Lining-layer fibroblast-2",
  "Lining-layer fibroblast-MMP3+",
  "Sublining-layer fibroblast-APOE+CXCL12+",
  "Sublining-layer fibroblast-APOD+CXCL14+",
  "Sublining-layer fibroblast-MFAP5+PI16+",
  "Sublining-layer fibroblast-COMP+",
  "Adipocyte"
))

cluster_colors_stromal <- c(
  "Lining-layer fibroblast-1" = "#5E0000",
  "Lining-layer fibroblast-2" = "#C73A3A",
  "Lining-layer fibroblast-MMP3+" = "#F0A0A0",
  "Sublining-layer fibroblast-APOE+CXCL12+" = "#001A33",
  "Sublining-layer fibroblast-APOD+CXCL14+" = "#285580",
  "Sublining-layer fibroblast-MFAP5+PI16+" = "#5FA6DA",
  "Sublining-layer fibroblast-COMP+" = "#D6EAF8",
  "Adipocyte" = "grey"
)

p_class3 <- DimPlot(sc_stromal, cols = cluster_colors_stromal, group.by = "CellType_Class3", shuffle = TRUE, raster = TRUE, pt.size = 1.5) +
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
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -4, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  ) +
  guides(colour = guide_legend(
    override.aes = list(size = 1)
  ))
ggsave("99_Fig/fig6/Dim_Class3_stromal.png", plot = p_class3, height = 1.5, width = 3)
ggsave("99_Fig/fig6/Dim_Class3_stromal.pdf", plot = p_class3, height = 1.5, width = 3)


sc_stromal <- SetIdent(sc_stromal, value = "CellType_Class2")


sc_fib_lining <- subset(sc_stromal, idents = "Lining-layer fibroblast")
sc_fib_lining <- SetIdent(sc_fib_lining, value = "CellType_Class3")

meta <- sc_fib_lining@meta.data
df_stacked_bar_class3 <- meta %>% 
  count(orig.ident, CellType_Class3)

p_stacked_bar_class3 <- ggplot(df_stacked_bar_class3, aes(x = orig.ident, y = n, fill = CellType_Class3)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.85) +
  labs(y = "Fraction of cell type", fill = "Cell type") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = cluster_colors_stromal) +
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
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -4, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )
ggsave("99_Fig/fig6/Stacked_Class3_stromal_lining.png", p_stacked_bar_class3, height = 1.5, width = 3)
ggsave("99_Fig/fig6/Stacked_Class3_stromal_lining.pdf", p_stacked_bar_class3, height = 1.5, width = 3)

all.markers <- FindAllMarkers(sc_fib_lining, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, recorrect_umi = FALSE)
all.markers_top10 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p_Dot <- DotPlot(sc_fib_lining, features = unique(all.markers_top10$gene), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_DEG_stromal_lining.png", p_Dot, height = 1.5, width = 5)
ggsave("99_Fig/sup_fig10/Dot_class3_DEG_stromal_lining.pdf", p_Dot, height = 1.5, width = 5)

p_Dot2 <- DotPlot(sc_fib_lining, features = c("CLIC5","CD55","PRG4","LRRC15","MMP3"), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_stromal_lining.png", p_Dot2, height = 1.5, width = 3)
ggsave("99_Fig/sup_fig10/Dot_class3_stromal_lining.pdf", p_Dot2, height = 1.5, width = 3)

p_Dot3 <- DotPlot(sc_fib_lining, features = c("MMP3","MMP9","MMP13","ADAMTS4","CTSK","CTSS","CTSG","TPSB2","TPSAB1"), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_catabolic_stromal_lining.png", p_Dot3, height = 1.5, width = 4)
ggsave("99_Fig/sup_fig10/Dot_class3_catabolic_stromal_lining.pdf", p_Dot3, height = 1.5, width = 4)


sc_fib_sublining <- subset(sc_stromal, idents = "Sublining-layer fibroblast")
sc_fib_sublining <- SetIdent(sc_fib_sublining, value = "CellType_Class3")

meta <- sc_fib_sublining@meta.data
df_stacked_bar_class3 <- meta %>% 
  count(orig.ident, CellType_Class3)

p_stacked_bar_class3 <- ggplot(df_stacked_bar_class3, aes(x = orig.ident, y = n, fill = CellType_Class3)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.85) +
  labs(y = "Fraction of cell type", fill = "Cell type") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = cluster_colors_stromal) +
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
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -4, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )
ggsave("99_Fig/fig6/Stacked_Class3_stromal_sublining.png", p_stacked_bar_class3, height = 1.5, width = 3.4)
ggsave("99_Fig/fig6/Stacked_Class3_stromal_sublining.pdf", p_stacked_bar_class3, height = 1.5, width = 3.4)

new_order <- c(
  "Sublining-layer fibroblast-COMP+",
  "Sublining-layer fibroblast-APOE+CXCL12+",
  "Sublining-layer fibroblast-APOD+CXCL14+",
  "Sublining-layer fibroblast-MFAP5+PI16+"
)
Idents(sc_fib_sublining) <- factor(Idents(sc_fib_sublining), levels = new_order)

all.markers <- FindAllMarkers(sc_fib_sublining, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, recorrect_umi = FALSE)
all.markers_top10 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p_Dot <- DotPlot(sc_fib_sublining, features = unique(all.markers_top10$gene), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_DEG_stromal_sublining.png", p_Dot, height = 1.5, width = 6)
ggsave("99_Fig/sup_fig10/Dot_class3_DEG_stromal_sublining.pdf", p_Dot, height = 1.5, width = 6)

p_Dot2 <- DotPlot(sc_fib_sublining, features = c("COL1A1","COL3A1","LRRC15","MMP3","THY1","SPP1","ACAN","CDH11","CXCL12","CD74","IL6","APOD","CXCL14","CD34","MFAP5","PI16"), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -10, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_stromal_sublining.png", p_Dot2, height = 1.5, width = 4.5)
ggsave("99_Fig/sup_fig10/Dot_class3_stromal_sublining.pdf", p_Dot2, height = 1.5, width = 4.5)

p_Dot3 <- DotPlot(sc_fib_sublining, features = c("MMP3","MMP9","MMP13","ADAMTS4","CTSK","CTSS","CTSG","TPSB2","TPSAB1"), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -10, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_catabolic_stromal_sublining.png", p_Dot3, height = 1.5, width = 4)
ggsave("99_Fig/sup_fig10/Dot_class3_catabolic_stromal_sublining.pdf", p_Dot3, height = 1.5, width = 4)



sc_vessel <- subset(sc, idents = c("Endothelial cell","Mural cell"))

DefaultAssay(sc_vessel) <- "RNA"
sc_vessel <- DietSeurat(sc_vessel, assays = "RNA")
sc_vessel[["RNA"]] <- split(sc_vessel[["RNA"]], f = sc_vessel$orig.ident)
sc_vessel <- SCTransform(sc_vessel, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), min_cells = 3)
sc_vessel <- RunPCA(sc_vessel)
sc_vessel <- IntegrateLayers(object = sc_vessel,
                             method = RPCAIntegration,
                             normalization.method = "SCT",
                             orig.reduction = "pca",
                             new.reduction = "integrated.rpca",
                             k.weight = 45,  # An error occurs under the default settings
                             verbose = FALSE
                             )
sc_vessel <- JoinLayers(sc_vessel, assay = "RNA")
sc_vessel <- RunUMAP(sc_vessel, dims = 1:10, reduction="integrated.rpca")

sc_vessel$CellType_Class3 <- factor(sc_vessel$CellType_Class3, levels = c(
  "Arterial EC",
  "Venous EC-1",
  "Venous EC-2",
  "Capillary EC-1",
  "Capillary EC-2",
  "Lymphatic EC",
  "VSMC-1",
  "VSMC-2",
  "Pericyte"
))

cluster_colors_vessel <- c(
  "Arterial EC" = "#D73027",
  "Venous EC-1" = "#99C2E6",
  "Venous EC-2" = "#003366",
  "Capillary EC-1" = "#E6DFF2",
  "Capillary EC-2" = "#54278F",
  "Lymphatic EC" = "#AE017E",
  "VSMC-1" = "#006400", 
  "VSMC-2" = "#A1D99B",
  "Pericyte" = "#FFD700"
)

p_class3 <- DimPlot(sc_vessel, cols = cluster_colors_vessel, group.by = "CellType_Class3", shuffle = TRUE, raster = TRUE, pt.size = 2) +
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
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -6, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  ) +
  guides(colour = guide_legend(
    override.aes = list(size = 1)
  ))
ggsave("99_Fig/fig6/Dim_Class3_vessel.png", plot = p_class3, height = 1.5, width = 2.25)
ggsave("99_Fig/fig6/Dim_Class3_vessel.pdf", plot = p_class3, height = 1.5, width = 2.25)


sc_vessel <- SetIdent(sc_vessel, value = "CellType_Class1")

sc_endothelial <- subset(sc_vessel, idents = "Endothelial cell")
sc_endothelial <- SetIdent(sc_endothelial, value = "CellType_Class3")

meta <- sc_endothelial@meta.data
df_stacked_bar_class3 <- meta %>% 
  count(orig.ident, CellType_Class3)

p_stacked_bar_class3 <- ggplot(df_stacked_bar_class3, aes(x = orig.ident, y = n, fill = CellType_Class3)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.85) +
  labs(y = "Fraction of cell type", fill = "Cell type") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = cluster_colors_vessel) +
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
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -4, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )
ggsave("99_Fig/fig6/Stacked_Class3_vessel_ec.png", p_stacked_bar_class3, height = 1.5, width = 2.48)
ggsave("99_Fig/fig6/Stacked_Class3_vessel_ec.pdf", p_stacked_bar_class3, height = 1.5, width = 2.48)

new_order <- c(
  "Arterial EC",
  "Capillary EC-1",
  "Capillary EC-2",
  "Venous EC-1",
  "Venous EC-2",
  "Lymphatic EC"
)
Idents(sc_endothelial) <- factor(Idents(sc_endothelial), levels = new_order)

all.markers <- FindAllMarkers(sc_endothelial, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, recorrect_umi = FALSE)
all.markers_top10 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p_Dot <- DotPlot(sc_endothelial, features = unique(all.markers_top10$gene), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -5, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_DEG_vessel_ec.png", p_Dot, height = 1.5, width = 6.5)
ggsave("99_Fig/sup_fig10/Dot_class3_DEG_vessel_ec.pdf", p_Dot, height = 1.5, width = 6.5)

p_Dot2 <- DotPlot(sc_endothelial, features = c("SEMA3G","SOX17","NOTCH4","RGCC","SPARC","LIFR","ICAM1","ACKR1","NR2F2","CCL21"), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -5, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_vessel_ec.png", p_Dot2, height = 1.5, width = 3.5)
ggsave("99_Fig/sup_fig10/Dot_class3_vessel_ec.pdf", p_Dot2, height = 1.5, width = 3.5)

p_Dot3 <- DotPlot(sc_endothelial, features = c("MMP3","MMP9","ADAMTS4","CTSK","CTSS","CTSG","TPSB2","TPSAB1"), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -10, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_catabolic_vessel_ec.png", p_Dot3, height = 1.5, width = 3)
ggsave("99_Fig/sup_fig10/Dot_class3_catabolic_vessel_ec.pdf", p_Dot3, height = 1.5, width = 3)


sc_mural <- subset(sc_vessel, idents = "Mural cell")
sc_mural <- SetIdent(sc_mural, value = "CellType_Class3")

meta <- sc_mural@meta.data
df_stacked_bar_class3 <- meta %>% 
  count(orig.ident, CellType_Class3)

p_stacked_bar_class3 <- ggplot(df_stacked_bar_class3, aes(x = orig.ident, y = n, fill = CellType_Class3)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.85) +
  labs(y = "Fraction of cell type", fill = "Cell type") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = cluster_colors_vessel) +
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
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -4, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )
ggsave("99_Fig/fig6/Stacked_Class3_vessel_mural.png", p_stacked_bar_class3, height = 1.5, width = 2.3)
ggsave("99_Fig/fig6/Stacked_Class3_vessel_mural.pdf", p_stacked_bar_class3, height = 1.5, width = 2.3)

new_order <- c(
  "VSMC-1",
  "VSMC-2",
  "Pericyte"
)
Idents(sc_mural) <- factor(Idents(sc_mural), levels = new_order)

all.markers <- FindAllMarkers(sc_mural, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, recorrect_umi = FALSE)
all.markers_top10 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

p_Dot <- DotPlot(sc_mural, features = unique(all.markers_top10$gene), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -10, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_DEG_vessel_mural.png", p_Dot, height = 1.5, width = 4.5)
ggsave("99_Fig/sup_fig10/Dot_class3_DEG_vessel_mural.pdf", p_Dot, height = 1.5, width = 4.5)

p_Dot2 <- DotPlot(sc_mural, features = c("ACTA2","MYH11","RERGL","CASQ2","KCNAB1","HMCN2","FLNC","PDGFRB","RGS5","ABCC9","KCNJ8","AGT"), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -10, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fisup_fig10g7/Dot_class3_vessel_mural.png", p_Dot2, height = 1.5, width = 4)
ggsave("99_Fig/sup_fig10/Dot_class3_vessel_mural.pdf", p_Dot2, height = 1.5, width = 4)

p_Dot3 <- DotPlot(sc_mural, features = c("MMP3","MMP9","ADAMTS4","CTSK","CTSS","CTSG","TPSB2","TPSAB1"), dot.scale = 2) + 
  scale_color_gradient2(low = "#B3D3E8", mid = "grey90", high = "#B94E5A", midpoint = 0) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 4)),
    legend.text  = element_text(size = 5),
    legend.key.size  = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(t = 0, r = 0, b = -10, l = -2, unit = "mm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
  )
ggsave("99_Fig/sup_fig10/Dot_class3_catabolic_vessel_mural.png", p_Dot3, height = 1.5, width = 3)
ggsave("99_Fig/sup_fig10/Dot_class3_catabolic_vessel_mural.pdf", p_Dot3, height = 1.5, width = 3)

