library(pheatmap)
library(dplyr)
library(tibble)
library(tidyverse)

colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

Cluster_colors <- c("1" = "#D8C6C2", 
                    "2" = "#F3B2A6",
                    "3" = "#C73A3A",
                    "4" = "#7F1212")

df <- read.delim(file = "01_formatted/annotations_full_tissue_conection_cluster.txt", sep = "\t", row.names = 1)

data <- read.delim("01_formatted/gene_read_count_normalized.txt", row.names = 1, header = T)

samples <- intersect(rownames(df), colnames(data))

anno <- data.frame(cluster = df[samples, "cluster", drop = FALSE])

identical(rownames(anno), colnames(data))
# [1] TRUE


anno$cluster <- factor(anno$cluster, levels = c("3", "4", "2", "1"))
anno <- anno[order(anno$cluster), , drop = FALSE]
data <- data[, rownames(anno), drop = FALSE]

data_order <- data[order(apply(data,1,sd), decreasing = T),]
top100_genes <- rownames(data_order[1:100,])

tmp <- data %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::filter(gene %in% top100_genes) %>% 
  separate(gene, c("id", "symbol"), sep = "_") %>% 
  select(-id) %>% 
  column_to_rownames(var = "symbol")

pdf("99_Fig/sup_fig5/Heatmap_variable_top100.pdf", width = 7, height = 12)
pheatmap(tmp, 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  
         color = colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = anno,
         annotation_colors = list(cluster = Cluster_colors),
         clustering_method = "ward.D2"
)
dev.off()


mat <- data
max_difference <- do.call(rbind, apply(mat, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))
max_difference$diff <- max_difference$max - max_difference$min
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

top100_genes_diff <- rownames(max_difference[1:100,])

tmp <- data %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::filter(gene %in% top100_genes_diff) %>% 
  separate(gene, c("id", "symbol"), sep = "_") %>% 
  select(-id) %>% 
  column_to_rownames(var = "symbol")

pdf("99_Fig/sup_fig5/Heatmap_variable_top100_diff.pdf", width = 7, height = 12)
pheatmap(tmp, 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  
         color = colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = anno,
         annotation_colors = list(cluster = Cluster_colors),
         clustering_method = "ward.D2"
)
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
# [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8    LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                    LC_TIME=Japanese_Japan.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     ggplot2_3.5.1   tidyverse_2.0.0 tibble_3.2.1    pheatmap_1.0.12
# 
# loaded via a namespace (and not attached):
#  [1] vctrs_0.6.5        cli_3.6.4          rlang_1.1.5        stringi_1.8.7      pkgload_1.4.0      generics_0.1.3     glue_1.8.0         colorspace_2.1-1   hms_1.1.3          scales_1.3.0       grid_4.3.3         munsell_0.5.1      tzdb_0.5.0         lifecycle_1.0.4   
# [15] compiler_4.3.3     RColorBrewer_1.1-3 timechange_0.3.0   pkgconfig_2.0.3    rstudioapi_0.17.1  farver_2.1.2       R6_2.6.1           tidyselect_1.2.1   pillar_1.10.1      magrittr_2.0.3     tools_4.3.3        withr_3.0.2        gtable_0.3.6  

