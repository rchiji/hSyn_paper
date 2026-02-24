library(NMF)
library(viridis)
library(pheatmap)

df <- read.delim("01_formatted/annotations_full_tissue_conection_cluster.txt", sep = "\t", row.names = 1)
Degree_df <- read.csv("00_src/Degree_df_241220.csv", sep = ",", header = TRUE, row.names = 1)


anno <- df[,c(1,3,4,44)]
anno$Diagnosis[anno$Diagnosis == "nonOA"] <- "ACLR"
anno$Diagnosis <- factor(anno$Diagnosis, levels = c("ACLR", "OA", "RA", "SLE"))
anno$cluster <- factor(anno$cluster, levels = c("1", "2", "3", "4"))

Sex_colors <- c("Male" = "#4E6FAE",
                "Female" = "#B94E5A")

Diagnosis_colors <- c("ACLR" = "#1f77b4", 
                      "RA" = "#2ca02c",
                      "OA" = "#ff7f0e",
                      "SLE" = "grey80")

Cluster_colors <- c("1" = "#D8C6C2", 
                    "2" = "#F3B2A6",
                    "3" = "#C73A3A",
                    "4" = "#7F1212")

Degree_df <- t(Degree_df)


# NMF
res_nmf <- nmf(Degree_df, rank = 6, seed=123, .options = "t")

W <- basis(res_nmf)
H <- coef(res_nmf)
V_reconstructed <- W %*% H

range01 <- function(x){
    x <- (x-min(x))/(max(x)-min(x))
    return(x)
}

W_scaled <- apply(W, 1, range01)
H_scaled <- apply(H, 1, range01)

conversion <- c(
  "Fibro.dense.irregular." = "Fibrous tissue (dense, irregular)", 
  "Fibro.dense.regular." = "Fibrous tissue (dense, regular)", 
  "Fibro.loose." = "Fibrous tissue (loose)", 
  "Immune.cells" = "TLS", 
  "adipose" = "Adipose", 
  "lining" = "Lining",
  "muscle" = "Muscle",
  "plasma" = "Plasma",
  "RBC" = "RBC",
  "Stroma" = "Stroma",
  "vessel" = "Micro vessel",
  "vessel.large." = "Large vessel"
)

original_names <- rownames(Degree_df)
new_names <- sapply(original_names, function(name) {
  parts <- unlist(strsplit(name, "_")) 
  converted_parts <- sapply(parts, function(part) {
    if (part %in% names(conversion)) {
      conversion[part] 
    } else {
      part 
    }
  })
  paste(converted_parts, collapse = "_")
})
rownames(Degree_df) <- new_names

pattern_cluster <- apply(W_scaled, 2, which.max)
cluster_names <- sort(unique(gsub(pattern = "_.*", replacement = "", names(pattern_cluster))))

anno_row <- data.frame(
  A = gsub(pattern = "_.*", replacement = "", rownames(Degree_df)),
  B = gsub(pattern = ".*_", replacement = "", rownames(Degree_df)),
  row.names = rownames(Degree_df)
  )
anno_row <- anno_row[,c(2,1)]

anno_row_color <- turbo(n = length(cluster_names))
names(anno_row_color) <- c("Adipose","Fibrous tissue (dense, irregular)","Fibrous tissue (dense, regular)","Fibrous tissue (loose)",
                           "TLS","Lining","Muscle","Plasma","RBC","Stroma","Micro vessel","Large vessel")

colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

anno <- anno[colnames(Degree_df), , drop = FALSE]
cluster_order <- c("2", "4", "3", "1")
order <- order(match(as.character(anno$cluster), cluster_order))
Degree_df <- Degree_df[, order, drop = FALSE]
anno <- anno[order, , drop = FALSE]

pdf(file = "99_Fig/fig4/pattern_heatmap.pdf", width = 8, height = 7.5)
pheatmap(Degree_df[rowSums(Degree_df) > 0.075,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = "row",
         breaks = seq(-1.5, 1.5, length.out=100),
         color = colors,
         annotation_row = anno_row,
         annotation_col = anno,
         annotation_colors = list(A = anno_row_color,
                                  B = anno_row_color,
                                  Sex = Sex_colors,
                                  Diagnosis = Diagnosis_colors,
                                  cluster = Cluster_colors)
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
# [1] pheatmap_1.0.12     viridis_0.6.5       viridisLite_0.4.2   NMF_0.28            Biobase_2.62.0      BiocGenerics_0.48.1 cluster_2.1.6       rngtools_1.5.2      registry_0.5-1     
# 
# loaded via a namespace (and not attached):
#  [1] gtable_0.3.6        dplyr_1.1.4         compiler_4.3.3      BiocManager_1.30.25 tidyselect_1.2.1    Rcpp_1.0.14         stringr_1.5.1       parallel_4.3.3      gridExtra_2.3       scales_1.3.0        ggplot2_3.5.1       R6_2.6.1            plyr_1.8.9         
# [14] generics_0.1.3      iterators_1.0.14    tibble_3.2.1        munsell_0.5.1       pillar_1.10.1       RColorBrewer_1.1-3  rlang_1.1.5         stringi_1.8.7       doParallel_1.0.17   cli_3.6.4           magrittr_2.0.3      digest_0.6.37       foreach_1.5.2      
# [27] grid_4.3.3          rstudioapi_0.17.1   gridBase_0.4-7      lifecycle_1.0.4     vctrs_0.6.5         glue_1.8.0          farver_2.1.2        codetools_0.2-19    colorspace_2.1-1    reshape2_1.4.4      tools_4.3.3         pkgconfig_2.0.3    

