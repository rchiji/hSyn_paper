library(dplyr)
library(tidyverse)
library(compositions)
library(edgeR)
library(DESeq2)
library(WGCNA)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(pheatmap)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
data <- read.delim("01_formatted/gene_read_count_preprocessed.txt", sep = "\t", header = TRUE, row.names = 1)


# WGCNA
## Preprocess
Min.CPM = 1
nSamples = 80
data_filtered <- round(data,0)
data_filtered <- data_filtered[which(apply(cpm(DGEList(counts = data_filtered)), 1, function(y) {sum(y>=Min.CPM)}) >= nSamples/2), ] # 13692 genes

data_norm <- data_filtered
tmp <- rep("A",dim(data_norm)[2]); tmp[1] <- "B"
colData <- cbind(colnames(data_norm), tmp); colnames(colData) = c("sample", "groups")
dds <- DESeqDataSetFromMatrix(countData = data_norm, colData = colData, design = ~ groups)

data_norm <- vst(dds)
data_norm <- assay(data_norm)

data_norm <- data_norm[order(apply(data_norm,1,sd), decreasing = T),]
data_sd_5000 <- as.matrix(data_norm[1:5000,])

data_wgcna <- t(as.matrix(data_sd_5000))


## WGCNA
powers  <-  c(1:30)
sft <- pickSoftThreshold(data_wgcna, powerVector = powers, networkType = "signed", verbose = 5)

# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"));
#  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red");
#  abline(h=0.90,col="red")
# plot(sft$fitIndices[,1], sft$fitIndices[,5],　xlab="Soft Threshold (power)",　ylab="Mean Connectivity", type="n",　main = paste("Mean connectivity"));
#  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

softPower <- 13
adjacency <- adjacency(data_wgcna, power = softPower, type="signed")
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1-TOM

k = as.vector(apply(adjacency, 2, sum, na.rm=T))
# hist(k)

scaleFreePlot(k, main="Check scale free topology\n")
#   scaleFreeRsquared slope
# 1              0.89 -1.59

geneTree <- hclust(as.dist(dissTOM), method = "average")

minModuleSize <- 30

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 1, pamStage = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 1055  713  597  532  422  340  318  224  221  164  119  114   73   62   46 

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# black        blue       brown        cyan       green greenyellow        grey     magenta        pink      purple         red      salmon         tan   turquoise      yellow 
#   224         597         532          46         340         114        1055         164         221         119         318          62          73         713         422 

module_palette <- c(
  "black" = "#2B2B2B", "blue" = "#4C72B0", "brown" = "#8C6D5A", "cyan" = "#5FA8A8", "green" = "#5E8C61",
  "greenyellow" = "#A9C75F", "grey" = "#7F7F7F", "magenta" = "#C0509A", "pink" = "#F29AB2", "purple" = "#8A5FBF",
  "red" = "#D94C4C", "salmon" = "#F08080", "tan" = "#C7A97B", "turquoise" = "#40B8B8", "yellow" = "#D6C36B"
)
dynamicColors_hex <- module_palette[dynamicColors]

pdf("99_Fig/sup_fig4/DendroAndColors.pdf", width = 4.5, height = 3, useDingbats = FALSE)
par(bg = "white")
plotDendroAndColors(geneTree, dynamicColors_hex, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()


MEList <- moduleEigengenes(data_wgcna, colors = dynamicColors)
MEs <- MEList$eigengenes
MEs_avr <- MEList$averageExpr

gene_module <- cbind(data_sd_5000, dynamicColors)



## Enrichment
gene_module <- as.data.frame(gene_module)
gene_module$gene <- rownames(gene_module)
gene_module <- separate(gene_module, gene, c("gene_id", "gene_name"), sep="[_]", extra="merge")

module_colors <- c("black", "blue", "brown", "cyan", "green", "greenyellow", "grey", "magenta", "pink", "purple", "red", "salmon", "tan", "turquoise", "yellow")


### GO
msigdb_go <- msigdbr(species = "Homo sapiens", collection = "C5")
msigdb_go <- msigdb_go %>%
  filter(!gs_subcollection %in% c("HPO", "GO:CC", "GO:MF"))

enrichment_results_go <- list() 
for (color in module_colors) {
  gene_data <- unique(gene_module$gene_name[gene_module$dynamicColors == color])
  enrichment_results_go[[color]] <- enricher(gene = gene_data,
                                          TERM2GENE = msigdb_go[, c("gs_name", "gene_symbol")])
}


results_go <- lapply(enrichment_results_go, function(x) x@result)

results_go <- lapply(results_go, function(df) {
  df[df$qvalue < 0.05, , drop = FALSE]
})

for (color in names(results_go)) {
  results_go[[color]] <- results_go[[color]] %>%
    mutate(color = color,
           log_qvalue = -log10(qvalue))
}

df_results_go <- bind_rows(results_go)

df_results_go_top <- df_results_go %>%
  group_by(color) %>%
  arrange(qvalue) %>%
  slice_head(n = 5) %>%
  ungroup()

color_order <- df_results_go_top %>%
  group_by(color) %>%
  summarise(mean_logq = mean(log_qvalue, na.rm = TRUE)) %>%
  arrange(desc(mean_logq)) %>%
  pull(color)

df_results_go_top$color <- factor(df_results_go_top$color, levels = color_order)

df_results_go_top <- df_results_go_top %>%
  arrange(color, qvalue)

df_results_go_top$Description_color <- interaction(df_results_go_top$Description, df_results_go_top$color)
df_results_go_top$Description_color <- factor(
  df_results_go_top$Description_color,
  levels = rev(unique(df_results_go_top$Description_color))
)

df_results_go_top$plot_color <- module_palette[as.character(df_results_go_top$color)]

labels <- df_results_go_top$Description_color %>%
  gsub("\\..*$", "", .) %>%
  sub("_", ": ", .) %>%
  gsub("_", " ", .)

custom_label <- setNames(labels, df_results_go_top$Description_color)

p_go <- ggplot(df_results_go_top, aes(x = Description_color, y = log_qvalue)) +
  geom_bar(stat = "identity", aes(fill = plot_color), width = 0.9) +
  coord_flip() +
  labs(y = "-log10(q-value)") +
  scale_fill_identity() +
  scale_x_discrete(labels = custom_label) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    panel.grid.major = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 4), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 4))

ggsave("99_Fig/sup_fig4/Entropy_enrichment_go.png", plot = p_go , width = 6, height = 4.5)
ggsave("99_Fig/sup_fig4/Entropy_enrichment_go.pdf", plot = p_go , width = 6, height = 4.5)


### cell type
msigdb_celltype <- msigdbr(species = "Homo sapiens", collection = "C8")

enrichment_results_celltype <- list() 
for (color in module_colors) {
  gene_data <- unique(gene_module$gene_name[gene_module$dynamicColors == color])
  enrichment_results_celltype[[color]] <- enricher(gene = gene_data,
                                                   TERM2GENE = msigdb_celltype[, c("gs_name", "gene_symbol")])
}


results_celltype <- lapply(enrichment_results_celltype, function(x) x@result)

results_celltype <- lapply(results_celltype, function(df) {
  df[df$qvalue < 0.05, , drop = FALSE]
})

for (color in names(results_celltype)) {
  results_celltype[[color]] <- results_celltype[[color]] %>%
    mutate(color = color,
           log_qvalue = -log10(qvalue))
}

df_results_celltype <- bind_rows(results_celltype)

df_results_celltype_top <- df_results_celltype %>%
  group_by(color) %>%
  arrange(qvalue) %>%
  slice_head(n = 5) %>%
  ungroup()

color_order <- df_results_celltype_top %>%
  group_by(color) %>%
  summarise(mean_logq = mean(log_qvalue, na.rm = TRUE)) %>%
  arrange(desc(mean_logq)) %>%
  pull(color)

df_results_celltype_top$color <- factor(df_results_celltype_top$color, levels = color_order)

df_results_celltype_top <- df_results_celltype_top %>%
  arrange(color, qvalue)

df_results_celltype_top$Description_color <- interaction(df_results_celltype_top$Description, df_results_celltype_top$color)
df_results_celltype_top$Description_color <- factor(
  df_results_celltype_top$Description_color,
  levels = rev(unique(df_results_celltype_top$Description_color))
)

df_results_celltype_top$plot_color <- module_palette[as.character(df_results_celltype_top$color)]

labels <- df_results_celltype_top$Description_color %>%
  gsub("\\..*$", "", .) %>%
  gsub("_", " ", .)

custom_label <- setNames(labels, df_results_celltype_top$Description_color)

p_cell <- ggplot(df_results_celltype_top, aes(x = Description_color, y = log_qvalue)) +
  geom_bar(stat = "identity", aes(fill = plot_color), width = 0.9) +
  coord_flip() +
  labs(y = "-log10(q-value)") +
  scale_fill_identity() +
  scale_x_discrete(labels = custom_label) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    panel.grid.major = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 4), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 4))

ggsave("99_Fig/sup_fig4/Entropy_enrichment_celltype.png", plot = p_cell , width = 3.5, height = 4.5)
ggsave("99_Fig/sup_fig4/Entropy_enrichment_celltype.pdf", plot = p_cell , width = 3.5, height = 4.5)



# correlation
donor_keep <- intersect(rownames(df), rownames(MEs))

df_filtered <- df[donor_keep,]
MEs_filtered <- MEs[donor_keep,]


## CLR
res_clr <- clr(acomp(df_filtered[,6:17])) 
df_clr <- as.data.frame(res_clr)


## Correlation
custom_labels_tissue <- c(
  "Fibrous_tissue__dense_irregular" = "Fibrous tissue (dense, irregular)",
  "Fibrous_tissue__dense_regular" = "Fibrous tissue (dense, regular)",
  "Fibrous_tissue__loose" = "Fibrous tissue (loose)",
  "Micro_vessel" = "Micro vessel",
  "Large_vessel" = "Large vessel",
  "Lining_thickness" = "Lining thickness",
  "Micro_vessel_ratio" = "Micro vessel ratio",
  "Large_vessel_ratio" = "Large vessel ratio",
  "TLS_ratio" = "TLS ratio"
)

custom_labels_module <- c(
  "MEblack" = "",
  "MEblue" = "",
  "MEbrown" = "",
  "MEcyan" = "",
  "MEgreen" = "",
  "MEgreenyellow" = "",
  "MEgrey" = "",
  "MEmagenta" = ,
  "MEpink" = "",
  "MEpurple" = "",
  "MEred" = "",
  "MEsalmon" = "",
  "MEtan" = "",
  "MEturquoise" = "",
  "MEyellow" = ""
)

colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)


### composition
cor_mat <- cor(df_clr, MEs_filtered, method = "spearman")

p_mat <- matrix(NA,
                nrow = nrow(cor_mat),
                ncol = ncol(cor_mat),
                dimnames = dimnames(cor_mat))

for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    p_mat[i, j] <- cor.test(df_clr[, i],
                            t(MEs_filtered)[j, ],
                            method = "spearman",
                            exact = FALSE)$p.value
  }
}

padj_mat <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
                   nrow = nrow(p_mat),
                   ncol = ncol(p_mat),
                   dimnames = dimnames(p_mat))

sig_mat <- matrix(
  "",
  nrow = nrow(padj_mat),
  ncol = ncol(padj_mat),
  dimnames = dimnames(padj_mat)
)
sig_mat[padj_mat < 0.05]  <- "*"
sig_mat[padj_mat < 0.01]  <- "**"
sig_mat[padj_mat < 0.001] <- "***"

rownames(cor_mat) <- ifelse(rownames(cor_mat) %in% names(custom_labels_tissue),
                            custom_labels_tissue[rownames(cor_mat)], rownames(cor_mat))
rownames(sig_mat) <- ifelse(rownames(sig_mat) %in% names(custom_labels_tissue),
                            custom_labels_tissue[rownames(sig_mat)], rownames(sig_mat))

colnames(cor_mat) <- ifelse(colnames(cor_mat) %in% names(custom_labels_module),
                            custom_labels_module[colnames(cor_mat)], colnames(cor_mat))
colnames(sig_mat) <- ifelse(colnames(sig_mat) %in% names(custom_labels_module),
                            custom_labels_module[colnames(sig_mat)], colnames(sig_mat))

pdf("99_Fig/sup_fig4/Heatmap_cor_tissue_genemodule.pdf", width = 4.5, height = 3)
pheatmap(cor_mat, 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_method = "ward.D2",
         display_numbers = sig_mat,
         fontsize = 6
)
dev.off()

