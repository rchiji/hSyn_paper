library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

options(future.globals.maxSize = 80 * 1024^3)

sc <- readRDS("02_Publicdata/ver6/tmp/seuratObj_add_label_latest_20251027.rds")
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

p_class3 <- DimPlot(sc_stromal, cols = cluster_colors_stromal, group.by = "CellType_Class3", shuffle = TRUE, raster = TRUE) +
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
ggsave("99_Fig/sup_fig7/Dot_class3_DEG_stromal_lining.png", p_Dot, height = 1.5, width = 5)
ggsave("99_Fig/sup_fig7/Dot_class3_DEG_stromal_lining.pdf", p_Dot, height = 1.5, width = 5)

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
ggsave("99_Fig/sup_fig7/Dot_class3_stromal_lining.png", p_Dot2, height = 1.5, width = 3)
ggsave("99_Fig/sup_fig7/Dot_class3_stromal_lining.pdf", p_Dot2, height = 1.5, width = 3)

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
ggsave("99_Fig/sup_fig7/Dot_class3_catabolic_stromal_lining.png", p_Dot3, height = 1.5, width = 4)
ggsave("99_Fig/sup_fig7/Dot_class3_catabolic_stromal_lining.pdf", p_Dot3, height = 1.5, width = 4)


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
ggsave("99_Fig/sup_fig7/Dot_class3_DEG_stromal_sublining.png", p_Dot, height = 1.5, width = 6)
ggsave("99_Fig/sup_fig7/Dot_class3_DEG_stromal_sublining.pdf", p_Dot, height = 1.5, width = 6)

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
ggsave("99_Fig/sup_fig7/Dot_class3_stromal_sublining.png", p_Dot2, height = 1.5, width = 4.5)
ggsave("99_Fig/sup_fig7/Dot_class3_stromal_sublining.pdf", p_Dot2, height = 1.5, width = 4.5)

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
ggsave("99_Fig/sup_fig7/Dot_class3_catabolic_stromal_sublining.png", p_Dot3, height = 1.5, width = 4)
ggsave("99_Fig/sup_fig7/Dot_class3_catabolic_stromal_sublining.pdf", p_Dot3, height = 1.5, width = 4)



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

p_class3 <- DimPlot(sc_vessel, cols = cluster_colors_vessel, group.by = "CellType_Class3", shuffle = TRUE, raster = TRUE) +
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
ggsave("99_Fig/sup_fig7/Dot_class3_DEG_vessel_ec.png", p_Dot, height = 1.5, width = 6.5)
ggsave("99_Fig/sup_fig7/Dot_class3_DEG_vessel_ec.pdf", p_Dot, height = 1.5, width = 6.5)

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
ggsave("99_Fig/sup_fig7/Dot_class3_vessel_ec.png", p_Dot2, height = 1.5, width = 3.5)
ggsave("99_Fig/sup_fig7/Dot_class3_vessel_ec.pdf", p_Dot2, height = 1.5, width = 3.5)

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
ggsave("99_Fig/sup_fig7/Dot_class3_catabolic_vessel_ec.png", p_Dot3, height = 1.5, width = 3)
ggsave("99_Fig/sup_fig7/Dot_class3_catabolic_vessel_ec.pdf", p_Dot3, height = 1.5, width = 3)


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
ggsave("99_Fig/sup_fig7/Dot_class3_DEG_vessel_mural.png", p_Dot, height = 1.5, width = 4.5)
ggsave("99_Fig/sup_fig7/Dot_class3_DEG_vessel_mural.pdf", p_Dot, height = 1.5, width = 4.5)

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
ggsave("99_Fig/sup_fig7/Dot_class3_vessel_mural.png", p_Dot2, height = 1.5, width = 4)
ggsave("99_Fig/sup_fig7/Dot_class3_vessel_mural.pdf", p_Dot2, height = 1.5, width = 4)

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
ggsave("99_Fig/sup_fig7/Dot_class3_catabolic_vessel_mural.png", p_Dot3, height = 1.5, width = 3)
ggsave("99_Fig/sup_fig7/Dot_class3_catabolic_vessel_mural.pdf", p_Dot3, height = 1.5, width = 3)


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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] cowplot_1.1.3      patchwork_1.3.0    lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1      purrr_1.0.4        readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0    dplyr_1.1.4        Seurat_5.1.0       SeuratObject_5.0.2
# [15] sp_2.2-0          
# 
# loaded via a namespace (and not attached):
#  [1] RColorBrewer_1.1-3          rstudioapi_0.17.1           jsonlite_2.0.0              magrittr_2.0.3              spatstat.utils_3.1-3        farver_2.1.2                zlibbioc_1.48.2             ragg_1.3.3                  vctrs_0.6.5                
# [10] ROCR_1.0-11                 DelayedMatrixStats_1.24.0   spatstat.explore_3.4-2      RCurl_1.98-1.17             S4Arrays_1.2.1              htmltools_0.5.8.1           SparseArray_1.2.4           sctransform_0.4.1           parallelly_1.43.0          
# [19] KernSmooth_2.23-22          htmlwidgets_1.6.4           ica_1.0-3                   plyr_1.8.9                  plotly_4.10.4               zoo_1.8-13                  igraph_2.1.4                mime_0.13                   lifecycle_1.0.4            
# [28] pkgconfig_2.0.3             Matrix_1.6-5                R6_2.6.1                    fastmap_1.2.0               GenomeInfoDbData_1.2.11     MatrixGenerics_1.14.0       fitdistrplus_1.2-2          future_1.34.0               shiny_1.10.0               
# [37] digest_0.6.37               colorspace_2.1-1            S4Vectors_0.40.2            tensor_1.5                  RSpectra_0.16-2             irlba_2.3.5.1               GenomicRanges_1.54.1        textshaping_1.0.0           labeling_0.4.3             
# [46] progressr_0.15.1            spatstat.sparse_3.1-0       timechange_0.3.0            httr_1.4.7                  polyclip_1.10-7             abind_1.4-8                 compiler_4.3.3              withr_3.0.2                 fastDummies_1.7.5          
# [55] MASS_7.3-60.0.1             DelayedArray_0.28.0         tools_4.3.3                 lmtest_0.9-40               httpuv_1.6.15               future.apply_1.11.3         goftest_1.2-3               glmGamPoi_1.14.3            glue_1.8.0                 
# [64] nlme_3.1-164                promises_1.3.2              grid_4.3.3                  Rtsne_0.17                  cluster_2.1.6               reshape2_1.4.4              generics_0.1.3              gtable_0.3.6                spatstat.data_3.1-6        
# [73] tzdb_0.5.0                  data.table_1.17.0           hms_1.1.3                   XVector_0.42.0              BiocGenerics_0.48.1         spatstat.geom_3.3-6         RcppAnnoy_0.0.22            ggrepel_0.9.6               RANN_2.6.2                 
# [82] pillar_1.10.1               spam_2.11-1                 RcppHNSW_0.6.0              limma_3.58.1                later_1.4.1                 splines_4.3.3               lattice_0.22-5              survival_3.5-8              deldir_2.0-4               
# [91] tidyselect_1.2.1            miniUI_0.1.1.1              pbapply_1.7-2               gridExtra_2.3               IRanges_2.36.0              SummarizedExperiment_1.32.0 scattermore_1.2             stats4_4.3.3                Biobase_2.62.0             
# [100] statmod_1.5.0               matrixStats_1.5.0           stringi_1.8.7               lazyeval_0.2.2              codetools_0.2-19            cli_3.6.4                   uwot_0.2.3                  xtable_1.8-4                reticulate_1.42.0          
# [109] systemfonts_1.2.1           munsell_0.5.1               GenomeInfoDb_1.38.8         Rcpp_1.0.14                 globals_0.16.3              spatstat.random_3.3-3       png_0.1-8                   spatstat.univar_3.1-2       parallel_4.3.3             
# [118] presto_1.0.0                dotCall64_1.2               sparseMatrixStats_1.14.0    bitops_1.0-9                listenv_0.9.1               viridisLite_0.4.2           scales_1.3.0                ggridges_0.5.6              crayon_1.5.3               
# [127] leiden_0.4.3.1              rlang_1.1.5

