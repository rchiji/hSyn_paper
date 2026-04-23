library(NMF)
library(viridis)
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
anno$cluster <- factor(anno$cluster, levels = c("1", "2", "3", "4"))

Sex_colors <- c("Male" = "#4E6FAE",
                "Female" = "#B94E5A")

Diagnosis_colors <- c("non-OA structural" = "#1f77b4", 
                      "RA" = "#2ca02c",
                      "OA" = "#ff7f0e",
                      "Other autoimmune diseases" = "grey80")

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


connection_select <- c("Lining_Micro vessel","Fibrous tissue (loose)_Lining","Fibrous tissue (loose)_Micro vessel",
                       "Micro vessel_Adipose",
                       "Large vessel_Adipose",
                       "Fibrous tissue (dense, irregular)_Micro vessel",
                       "TLS_TLS","Stroma_Stroma")
Degree_df_select <- Degree_df %>% 
  t() %>% 
  as.data.frame() %>% 
  select(all_of(connection_select))
Degree_df_select <- Degree_df_select[rownames(df), ]
Degree_df_select[,c(9,10)] <- df[,c(1,44)]
Degree_df_select <- Degree_df_select[Degree_df_select$Diagnosis == "OA",]

# Statistics
connections <- colnames(Degree_df_select[,1:8])

results <- list()
for (connection in connections) {
  kw_result <- kruskal.test(as.formula(paste0("`", connection, "` ~ cluster")), data = Degree_df_select)
  dunn_result <- dunn.test(Degree_df_select[[connection]], Degree_df_select$cluster, method = "bh")
  significant_idx <- which(dunn_result$P.adjusted < 0.1)
  significant_comparisons <- dunn_result$comparisons[significant_idx]
  significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
  significant_z <- dunn_result$Z[significant_idx]
  results[[connection]] <- list(
    test_type = "Kruskal-Wallis",
    main_test = kw_result,
    posthoc_test = dunn_result,
    significant_comparisons = significant_comparisons,
    significant_p_adjusted = significant_p_adjusted,
    significant_z = significant_z 
  )
}

for (connection in connections) {
  cat("\n--- Results for", connection, "---\n")
  print(results[[connection]]$main_test)
  if (length(results[[connection]]$significant_comparisons) > 0) {
    cat("\n** Significant Dunn Post-hoc Test Results for", connection, "**\n")
    cat("=============================================\n")
    for (i in seq_along(results[[connection]]$significant_comparisons)) {
      cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
                  results[[connection]]$significant_comparisons[i], 
                  results[[connection]]$significant_p_adjusted[i],
                  results[[connection]]$significant_z[i]))
    }
    cat("=============================================\n")
  }
}

# Notes: Significant differences are shown.
# ** Significant Dunn Post-hoc Test Results for Lining_Micro vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00759 - Z-value: -2.57182
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00816 - Z-value: 2.64540
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 5.10282
# >> Comparison: 1 - 4      - Adjusted P-value: 0.08936 - Z-value: 1.44331
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00007 - Z-value: 4.07791
# >> Comparison: 3 - 4      - Adjusted P-value: 0.09604 - Z-value: -1.30444
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous tissue (loose)_Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: -4.87000
# >> Comparison: 1 - 3      - Adjusted P-value: 0.06493 - Z-value: 1.60629
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.06207
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00002 - Z-value: 4.27876
# >> Comparison: 3 - 4      - Adjusted P-value: 0.03172 - Z-value: -2.03061
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous tissue (loose)_Micro vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: -4.76006
# >> Comparison: 1 - 3      - Adjusted P-value: 0.02842 - Z-value: 1.98295
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.36211
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00001 - Z-value: 4.32891
# >> Comparison: 3 - 4      - Adjusted P-value: 0.01747 - Z-value: -2.26853
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro vessel_Adipose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01621 - Z-value: -2.21121
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00001 - Z-value: -4.34519
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00666 - Z-value: -2.61674
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00000 - Z-value: -6.29286
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00001 - Z-value: -4.55386
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07711 - Z-value: -1.42478
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large vessel_Adipose **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00006 - Z-value: -4.25316
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00003 - Z-value: -4.23929
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00598 - Z-value: -2.74894
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00597 - Z-value: -2.65403
# >> Comparison: 3 - 4      - Adjusted P-value: 0.05263 - Z-value: 1.70753
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous tissue (dense, irregular)_Micro vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 5.61823
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00816 - Z-value: 2.64540
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01530 - Z-value: -2.23374
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00000 - Z-value: 5.43069
# >> Comparison: 3 - 4      - Adjusted P-value: 0.01516 - Z-value: 2.32227
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS_TLS **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.07097 - Z-value: -1.80565
# >> Comparison: 1 - 3      - Adjusted P-value: 0.07939 - Z-value: 1.50503
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00398 - Z-value: 3.20991
# >> Comparison: 1 - 4      - Adjusted P-value: 0.07921 - Z-value: -1.61825
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00460 - Z-value: -2.96085
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma_Stroma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00979 - Z-value: 2.72003
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00925 - Z-value: 2.95944
# >> Comparison: 1 - 4      - Adjusted P-value: 0.07910 - Z-value: 1.75596
# =============================================


Degree_df_long <- Degree_df_select %>%
  pivot_longer(cols = 1:8, names_to = "Connection", values_to = "Value")

Degree_df_long$Connection <- factor(Degree_df_long$Connection, levels = connection_select)
Degree_df_long$cluster <- factor(Degree_df_long$cluster)

plot_tissue <- ggplot(Degree_df_long, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_wrap(~ Connection, scales = "free", strip.position = "bottom", ncol = 3) + 
  labs(y = "Degree of connections") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#F3B2A6", "3" = "#C73A3A", "4" = "#7F1212")) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 5, margin = margin(t = 0)),
    strip.background = element_blank(),
    strip.placement  = "outside",
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, , margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin     = margin(0, 0, 0, 0))

ggsave("99_Fig/fig4/tSNE_cluster_connection.png", plot = plot_tissue, width = 4, height = 2.5)
ggsave("99_Fig/fig4/tSNE_cluster_connrction.pdf", plot = plot_tissue, width = 4, height = 2.5)

