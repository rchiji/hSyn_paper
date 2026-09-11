library(NMF)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(dunn.test)

df <- read.delim("01_formatted/annotations_full_tissue_conection_cluster.txt", sep = "\t", row.names = 1)
Degree_df <- read.csv("00_src/Degree_df_241220.csv", sep = ",", header = TRUE, row.names = 1)


anno <- df[,c(1,3,4,44)]
anno$Diagnosis[anno$Diagnosis == "nonOA"] <- "non-OA structural"
anno$Diagnosis[anno$Diagnosis %in% c("SLE", "SSc")] <- "Other autoimmune diseases"
anno$Diagnosis <- factor(anno$Diagnosis, levels = c("non-OA structural", "OA", "RA", "Other autoimmune diseases"))
anno$cluster <- factor(anno$cluster, levels = c("1", "2", "3"))

Sex_colors <- c("Male" = "#4E6FAE",
                "Female" = "#B94E5A")

Diagnosis_colors <- c("non-OA structural" = "#1f77b4", 
                      "RA" = "#2ca02c",
                      "OA" = "#ff7f0e",
                      "Other autoimmune diseases" = "grey80")

Cluster_colors <- c("1" = "#D8C6C2", 
                    "2" = "#DD7670",
                    "3" = "#7F1212")

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

anno_row_color <- c(
  "Adipose" = "#808000",
  "Fibrous tissue (dense, irregular)" = "#67da17",
  "Fibrous tissue (dense, regular)" = "#ecffb3",
  "Fibrous tissue (loose)" = "#83edaa",
  "Lining" = "#44d2e5",
  "TLS" = "#a05aa0",
  "Plasma" = "#30d383",
  "Stroma" = "#698c69",
  "RBC" = "magenta",
  "Micro vessel" = "#59220a",
  "Large vessel" = "#f68d14",
  "Muscle" = "#b299dc"
)

colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

anno <- anno[colnames(Degree_df), , drop = FALSE]
cluster_order <- c("1","2","3")
order <- order(match(as.character(anno$cluster), cluster_order))
Degree_df <- Degree_df[, order, drop = FALSE]
anno <- anno[order, , drop = FALSE]

pdf(file = "99_Fig/fig3/pattern_heatmap.pdf", width = 8, height = 4)
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

