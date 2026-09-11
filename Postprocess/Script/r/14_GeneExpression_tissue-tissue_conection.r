library(pheatmap)
library(dplyr)
library(tibble)
library(tidyverse)
library(edgeR)
library(DESeq2)

colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

Cluster_colors <- c("1" = "#D8C6C2", 
                    "2" = "#DD7670",
                    "3" = "#7F1212")

df <- read.delim(file = "01_formatted/annotations_full_tissue_conection_cluster.txt", sep = "\t", row.names = 1)
df_OA <- df[df$Diagnosis == "OA",]
donor_OA <- rownames(df_OA)

data <- read.delim("01_formatted/gene_read_count_preprocessed.txt", row.names = 1, header = T)
donor_select <- donor_OA[donor_OA %in% colnames(data)]
data <- data[, donor_select]


Min.CPM = 1
nSamples = 67
data_filtered <- round(data,0)
data_filtered <- data_filtered[which(apply(cpm(DGEList(counts = data_filtered)), 1, function(y) {sum(y>=Min.CPM)}) >= nSamples/2), ] # 13697 genes

data_norm <- data_filtered
tmp <- rep("A",dim(data_norm)[2]); tmp[1] <- "B"
colData <- cbind(colnames(data_norm), tmp); colnames(colData) = c("sample", "groups")
dds <- DESeqDataSetFromMatrix(countData = data_norm, colData = colData, design = ~ groups)

data_norm <- vst(dds)
data_norm <- assay(data_norm)
data_norm <- as.data.frame(data_norm)


anno <- data.frame(cluster = df[donor_select, "cluster", drop = FALSE])

identical(rownames(anno), colnames(data_norm))
# [1] TRUE


anno$cluster <- factor(anno$cluster, levels = c("1", "2", "3"))
anno <- anno[order(anno$cluster), , drop = FALSE]
data_norm <- data_norm[, rownames(anno), drop = FALSE]

data_order <- data_norm[order(apply(data_norm,1,sd), decreasing = T),]
top200_genes <- rownames(data_order[1:200,])

tmp <- data_norm %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::filter(gene %in% top200_genes) %>% 
  separate(gene, c("id", "symbol"), sep = "_") %>% 
  select(-id) %>% 
  column_to_rownames(var = "symbol")

pdf("99_Fig/fig3/Heatmap_variable_top200.pdf", width = 7, height = 12)
pheatmap(tmp, 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  
         color = colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = anno,
         annotation_colors = list(cluster = Cluster_colors),
         clustering_method = "ward.D2",
         fontsize = 5
)
dev.off()

