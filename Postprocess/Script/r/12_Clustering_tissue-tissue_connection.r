library(Rtsne)
library(ggplot2)
library(dunn.test)
library(dplyr)
library(tidyverse)
library(compositions)
library(emmeans)
library(broom)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
score <- read.delim("01_formatted/annotations_full_OARSI_Krenn.txt", sep = "\t", row.names = 1)
df[,37:41]  <- score[,1:5]
df <- df[df$Joint == "Knee",]
keep <- rownames(df)

Degree_df <- read.csv("00_src/Degree_df_241220.csv", sep = ",", header = TRUE, row.names = 1)
Degree_df <- Degree_df[keep,]

entropy <- read.delim(file = "00_src/Entropy_group.txt", sep = "\t", row.names = 1)
entropy <- as.data.frame(entropy[intersect(rownames(df), rownames(entropy)),])


pca_result <- prcomp(Degree_df[,1:78], scale. = FALSE)

plot(pca_result, type = "l")

summary(pca_result)
# Importance of components:
#                            PC1     PC2     PC3     PC4    PC5     PC6     PC7     PC8      PC9     PC10     PC11    PC12     PC13     PC14
# Standard deviation     0.09322 0.07334 0.03804 0.03353 0.0201 0.01599 0.01492 0.01157 0.009451 0.007568 0.007107 0.00645 0.005015 0.004721
# Proportion of Variance 0.48176 0.29816 0.08021 0.06233 0.0224 0.01417 0.01234 0.00742 0.004950 0.003170 0.002800 0.00231 0.001390 0.001240
# Cumulative Proportion  0.48176 0.77992 0.86013 0.92246 0.9449 0.95903 0.97136 0.97878 0.983730 0.986910 
# ...

pca_data <- data.frame(pca_result$x)


# # silhouette
# library(cluster)
# set.seed(123)
# silhouette_results <- lapply(2:8, function(k) {
#   km <- kmeans(pca_data[,1:3], centers = k)
#   sil <- silhouette(km$cluster, dist(pca_data[,1:3]))
#   data.frame(
#     k = k,
#     silhouette_score = mean(sil[,"sil_width"])
#   )
# })
# 
# silhouette_summary <- bind_rows(silhouette_results)
# silhouette_summary
# #   k silhouette_score
# # 1 2        0.3520918
# # 2 3        0.3623745
# # 3 4        0.3070929
# # 4 5        0.3221676
# # 5 6        0.3306463
# # 6 7        0.3022818
# # 7 8        0.3089276
# 
# write.table(silhouette_summary, "03_res/Clustering/tissue_conection_cluster_k.txt", sep = "\t", row.names = FALSE)
# 
# plot_silhouette_k <- ggplot(silhouette_summary, aes(x = k, y = silhouette_score)) +
#   geom_line(linewidth = 0.1) +
#   geom_point(size = 0.5) +
#   scale_x_continuous(breaks = 2:8) +
#   labs(x = "Number of clusters", y = "Average silhouette score") +
#   theme_classic() +
#   theme(
#     panel.border = element_blank(),
#     axis.ticks = element_line(linewidth = 0.1),
#     axis.line = element_line(linewidth = 0.1),
#     axis.title.x = element_text(size = 6),
#     axis.text.x = element_text(size = 5),
#     axis.title.y = element_text(size = 6),
#     axis.text.y = element_text(size = 5)
#   )
# 
# ggsave("99_Fig/fig3/silhouette_cluster_number.png", plot = plot_silhouette_k, width = 2, height = 1.5)
# ggsave("99_Fig/fig3/silhouette_cluster_number.pdf", plot = plot_silhouette_k, width = 2, height = 1.5)


# library(fpc)
# set.seed(123)
# 
# cb <- clusterboot(
#   pca_data[,1:3],
#   B = 1000,
#   bootmethod = "boot",
#   clustermethod = kmeansCBI,
#   krange = 3,
#   seed = 123
# )
# 
# cb$bootmean
# [1] 0.9018762 0.8880804 0.8300214


# tSNE
set.seed(123)
tsne_res <- Rtsne(pca_data[,1:3], dims= 2, perplexity = 36)
tsne_data <- data.frame(tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2])
rownames(tsne_data) <- rownames(Degree_df)


# clustering
kmeans_result <- kmeans(pca_data[,1:3], centers = 3)

tsne_data$cluster <- as.factor(kmeans_result$cluster)
tsne_data$Diagnosis <- df$Diagnosis

tsne_data$Diagnosis_plot <- recode(
  tsne_data$Diagnosis,
  "SLE" = "Other autoimmune diseases",
  "SSc" = "Other autoimmune diseases"
)

plot_tsne <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = cluster, shape = Diagnosis_plot)) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#DD7670", "3" = "#7F1212")) +
  scale_shape_manual(values = c("nonOA" = 16, "RA" = 15, "OA" = 17, "Other autoimmune diseases" = 18),
                     labels = c("nonOA" = "non-OA structural", "RA" = "RA", "OA" = "OA", "Other autoimmune diseases" = "Other autoimmune diseases")) +
  labs(x = "tSNE1", y = "tSNE2", shape = "Diagnosis", color = "kmeans cluster") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size = 5), 
    axis.text = element_blank(),
    legend.text = element_text(size = 5, margin = margin(l = 0)),
    legend.title = element_text(size = 5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), 
    legend.key.size = unit(0.4, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -8)
  )

ggsave("99_Fig/fig3/tSNE.png", plot = plot_tsne, width = 2.25, height = 1.25)
ggsave("99_Fig/fig3/tSNE.pdf", plot = plot_tsne, width = 2.25, height = 1.25)

df[,42:44] <- tsne_data[,1:3]

# write.table(df, "01_formatted/annotations_full_tissue_conection_cluster.txt", sep = "\t", row.names = TRUE, col.names = NA)


# Proportion
df_stacked_bar <- df
df_stacked_bar$Entropy <- entropy$`entropy[intersect(rownames(df), rownames(entropy)), ]`

df_stacked_bar <- df_stacked_bar %>% 
  count(cluster, Entropy)

df_stacked_bar$Entropy <- factor(df_stacked_bar$Entropy, levels = c("Low", "Mid", "High"))

plot_stacked_bar <- ggplot(df_stacked_bar, aes(x = cluster, y = n, fill = Entropy)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.9) +
  labs(y = "Fraction of entropy type", fill = "Entropy") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = c("Low" = "#bfe6bf", "Mid" = "#66cc66", "High" = "#336633")) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 5),
    axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
    axis.text.y  = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 2)),
    legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, -3, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )

ggsave("99_Fig/fig3/cluster_entropy_proportion.png", plot = plot_stacked_bar, width = 1.25, height = 1.25)
ggsave("99_Fig/fig3/cluster_entropy_proportion.pdf", plot = plot_stacked_bar, width = 1.25, height = 1.25)


# Statistics
df_OA <- df[df$Diagnosis == "OA", ]

scores <- colnames(df_OA[,18:32])
scores_2 <- colnames(df_OA[,c(35,36,41,37:40)])

## The same applies to "score_2".
results <- list()
for (score in scores) {
  kw_result <- kruskal.test(as.formula(paste(score, "~ cluster")), data = df_OA)
  dunn_result <- dunn.test(df_OA[[score]], df_OA$cluster, method = "bh")
  significant_idx <- which(dunn_result$P.adjusted < 0.1)
  significant_comparisons <- dunn_result$comparisons[significant_idx]
  significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
  significant_z <- dunn_result$Z[significant_idx]
  results[[score]] <- list(
    test_type = "Kruskal-Wallis",
    main_test = kw_result,
    posthoc_test = dunn_result,
    significant_comparisons = significant_comparisons,
    significant_p_adjusted = significant_p_adjusted,
    significant_z = significant_z 
  )
}

for (score in scores) {
  cat("\n--- Results for", score, "---\n")
  print(results[[score]]$main_test)
    if (length(results[[score]]$significant_comparisons) > 0) {
    cat("\n** Significant Dunn Post-hoc Test Results for", score, "**\n")
    cat("=============================================\n")
    for (i in seq_along(results[[score]]$significant_comparisons)) {
      cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
                  results[[score]]$significant_comparisons[i], 
                  results[[score]]$significant_p_adjusted[i],
                  results[[score]]$significant_z[i]))
    }
    cat("=============================================\n")
  }
}

# Notes: Significant differences are shown.
# ** Significant Dunn Post-hoc Test Results for Pre__KOOS_Symptom **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.05432 - Z-value: -1.79639
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08368 - Z-value: 1.91272
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_3Month__KOOS_QOL **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01802 - Z-value: -2.51166
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_QOL **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.06798 - Z-value: 1.69202
# >> Comparison: 2 - 3      - Adjusted P-value: 0.03611 - Z-value: 2.25600
# =============================================
# 
# ** Significant Dunn Post-hoc Test Results for Synovitis_score **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.07060 - Z-value: 1.67400
# >> Comparison: 1 - 3      - Adjusted P-value: 0.07133 - Z-value: 1.98133
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 5.15059
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00502 - Z-value: 2.71155
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02670 - Z-value: -1.93162
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Stroma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.09323 - Z-value: 1.32113
# >> Comparison: 1 - 3      - Adjusted P-value: 0.05968 - Z-value: -1.75314
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00995 - Z-value: -2.71475
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Infiltrate **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00176 - Z-value: 3.24532
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00160 - Z-value: 3.07045
# =============================================
# ** Significant Dunn Post-hoc Test Results for Total **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 4.73100
# >> Comparison: 1 - 3      - Adjusted P-value: 0.02510 - Z-value: 2.12642
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01767 - Z-value: -2.10438
# =============================================


Degree_df <- t(Degree_df)

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
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 5.88658
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00467 - Z-value: 2.59962
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00586 - Z-value: -2.66025
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous tissue (loose)_Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 6.48820
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: 4.96951
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous tissue (loose)_Micro vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 6.50944
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: 4.94533
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro vessel_Adipose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00036 - Z-value: -3.49393
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00070 - Z-value: 3.19613
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 5.87415
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large vessel_Adipose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00003 - Z-value: -4.24343
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00002 - Z-value: 4.21394
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous tissue (dense, irregular)_Micro vessel **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.62647
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00001 - Z-value: -4.31866
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS_TLS **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00070 - Z-value: 3.49900
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00312 - Z-value: 2.86530
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma_Stroma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.05384 - Z-value: 1.60868
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00981 - Z-value: -2.48156
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00044 - Z-value: -3.61997
# =============================================


Degree_df_long <- Degree_df_select %>%
  pivot_longer(cols = 1:8, names_to = "Connection", values_to = "Value")

Degree_df_long$Connection <- factor(Degree_df_long$Connection, levels = connection_select)
Degree_df_long$cluster <- factor(Degree_df_long$cluster)

plot_tissue <- ggplot(Degree_df_long, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_wrap(~ Connection, scales = "free", strip.position = "bottom", ncol = 4) + 
  labs(y = "Degree of connections") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#DD7670", "3" = "#7F1212")) +
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

ggsave("99_Fig/fig3/cluster_connection.png", plot = plot_tissue, width = 4, height = 2)
ggsave("99_Fig/fig3/cluster_connrction.pdf", plot = plot_tissue, width = 4, height = 2)


all_component <- colnames(df_OA)[6:17]
major_component <- c(
  "Adipose",
  "Fibrous_tissue__dense_irregular",
  "Fibrous_tissue__dense_regular",
  "Fibrous_tissue__loose"
)
minor_component <- setdiff(all_component, major_component)

res_clr_major <- clr(acomp(df_OA[,all_component])) 

df_clr_major <- as.data.frame(res_clr_major)
df_clr_major$cluster <- df_OA$cluster


# statistics
## shapiro_test
shapiro_test_results <- df_clr_major %>% 
  summarise(across(1:12, ~ shapiro.test(.)$p.value))
shapiro_test_results
#       Adipose Fibrous_tissue__dense_irregular Fibrous_tissue__dense_regular Fibrous_tissue__loose       TLS    Plasma       Stroma       Lining      Muscle       RBC Micro_vessel Large_vessel
# 1 0.001620628                      0.01034722                     0.2144998             0.8376733 0.3980882 0.1245807 5.468267e-05 4.213218e-06 9.57414e-06 0.1283558    0.6595747   0.01884009

## kruskal.test
results <- list()
for (comp in all_component) {
  kw_result <- kruskal.test(as.formula(paste(comp, "~ cluster")), data = df_clr_major)
  dunn_result <- dunn.test(df_clr_major[[comp]], df_clr_major$cluster, method = "bh")
  significant_idx <- which(dunn_result$P.adjusted < 0.1)
  significant_comparisons <- dunn_result$comparisons[significant_idx]
  significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
  significant_z <- dunn_result$Z[significant_idx]
  results[[comp]] <- list(
    test_type = "Kruskal-Wallis",
    main_test = kw_result,
    posthoc_test = dunn_result,
    significant_comparisons = significant_comparisons,
    significant_p_adjusted = significant_p_adjusted,
    significant_z = significant_z
  )
}

for (comp in all_component) {
  cat("\n--- Results for", comp, "---\n")
  print(results[[comp]]$main_test)
  if (length(results[[comp]]$significant_comparisons) > 0) {
    cat("\n** Significant Dunn Post-hoc Test Results for", comp, "**\n")
    cat("=============================================\n")
    for (i in seq_along(results[[comp]]$significant_comparisons)) {
      cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
                  results[[comp]]$significant_comparisons[i], 
                  results[[comp]]$significant_p_adjusted[i],
                  results[[comp]]$significant_z[i]))
    }
    cat("=============================================\n")
  }
}

# Notes: Significant differences are shown.
# ** Significant Dunn Post-hoc Test Results for Adipose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: -4.56837
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00032 - Z-value: 3.41780
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.99066
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00001 - Z-value: -4.34304
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.81590
# >> Comparison: 2 - 3      - Adjusted P-value: 0.05820 - Z-value: -1.57010
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00023 - Z-value: -3.60743
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.14282
# >> Comparison: 2 - 3      - Adjusted P-value: 0.05626 - Z-value: -1.58694
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01821 - Z-value: 2.09222
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00001 - Z-value: 4.45765
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01800 - Z-value: 2.25716
# =============================================


row_sum <- rowSums(df_OA[, minor_component], na.rm = TRUE)

df_recal <- df_OA
df_recal[, minor_component] <- df_OA[, minor_component] / row_sum

res_clr_minor <- clr(acomp(df_recal[,minor_component])) 

df_clr_minor <- as.data.frame(res_clr_minor)
df_clr_minor$cluster <- df_OA$cluster


# statistics
## shapiro_test
shapiro_test_results <- df_clr_minor %>% 
  summarise(across(1:8, ~ shapiro.test(.)$p.value))
shapiro_test_results 
# TLS     Plasma      Stroma       Lining       Muscle        RBC Micro_vessel Large_vessel
# 1 0.2667346 0.05463301 0.001368096 2.547756e-05 8.975871e-05 0.04094924    0.1877178    0.4792868

## kruskal.test
results <- list()
for (comp in minor_component) {
  kw_result <- kruskal.test(as.formula(paste(comp, "~ cluster")), data = df_clr_minor)
  dunn_result <- dunn.test(df_clr_minor[[comp]], df_clr_minor$cluster, method = "bh")
  significant_idx <- which(dunn_result$P.adjusted < 0.1)
  significant_comparisons <- dunn_result$comparisons[significant_idx]
  significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
  significant_z <- dunn_result$Z[significant_idx]
  results[[comp]] <- list(
    test_type = "Kruskal-Wallis",
    main_test = kw_result,
    posthoc_test = dunn_result,
    significant_comparisons = significant_comparisons,
    significant_p_adjusted = significant_p_adjusted,
    significant_z = significant_z
  )
}

for (comp in minor_component) {
  cat("\n--- Results for", comp, "---\n")
  print(results[[comp]]$main_test)
  if (length(results[[comp]]$significant_comparisons) > 0) {
    cat("\n** Significant Dunn Post-hoc Test Results for", comp, "**\n")
    cat("=============================================\n")
    for (i in seq_along(results[[comp]]$significant_comparisons)) {
      cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
                  results[[comp]]$significant_comparisons[i], 
                  results[[comp]]$significant_p_adjusted[i],
                  results[[comp]]$significant_z[i]))
    }
    cat("=============================================\n")
  }
}

# Notes: Significant differences are shown.
# ** Significant Dunn Post-hoc Test Results for TLS **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01749 - Z-value: 2.52218
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01263 - Z-value: 2.39004
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00027 - Z-value: 3.56703
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00009 - Z-value: 4.01880
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00708 - Z-value: -2.82538
# >> Comparison: 2 - 3      - Adjusted P-value: 0.03938 - Z-value: -1.93899
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00008 - Z-value: 4.04485
# >> Comparison: 1 - 3      - Adjusted P-value: 0.03868 - Z-value: 1.94669
# >> Comparison: 2 - 3      - Adjusted P-value: 0.04623 - Z-value: -1.68255
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.02255 - Z-value: -2.16924
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00684 - Z-value: -2.83680
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00001 - Z-value: -4.60620
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00127 - Z-value: 3.13989
# =============================================


## Plot
df_OA[,major_component] <- df_OA[,major_component]
df_OA[,minor_component] <- df_OA[,minor_component]

df_long <- df_OA %>%
  pivot_longer(cols = 6:17, names_to = "Tissue", values_to = "Value")

df_long$Tissue <- factor(df_long$Tissue, levels = c("Fibrous_tissue__dense_irregular", "Fibrous_tissue__dense_regular", "Fibrous_tissue__loose", "Lining", "TLS", "Plasma", "Stroma", "RBC", "Micro_vessel", "Large_vessel", "Adipose", "Muscle"))
df_long$cluster <- factor(df_long$cluster)

custom_labels <- c(
     "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense.irregular)",
     "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense.regular)",
     "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
     "Micro_vessel" = "Micro vessel",
     "Large_vessel" = "Large vessel"
 )

df_long$Tissue <- factor(
  df_long$Tissue,
  levels = c(
    "Lining",
    "Micro_vessel",
    "RBC",
    "Fibrous_tissue__loose",
    "TLS",
    "Plasma",
    "Large_vessel",
    "Adipose",
    "Fibrous_tissue__dense_irregular",
    "Fibrous_tissue__dense_regular",
    "Stroma",
    "Muscle"
  )
)

plot_tissue <- ggplot(df_long, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels), strip.position = "bottom", ncol = 6) + 
  labs(y = "CLR(Proportion)") +
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

ggsave("99_Fig/sup_fig7/cluster_tissue_proportion_CLR.png", plot = plot_tissue, width = 6, height = 2)
ggsave("99_Fig/sup_fig7/cluster_tissue_proportion_CLR.pdf", plot = plot_tissue, width = 6, height = 2)


df_long_KOOS <- df_OA %>%
  pivot_longer(cols = 18:32, names_to = "KOOS", values_to = "Value") %>%
  mutate(KOOS_group = case_when(
    str_detect(KOOS, "Pre") ~ "Pre",
    str_detect(KOOS, "Post_12") ~ "Post12",
    TRUE ~ "Other"
  )
)

df_long_KOOS_pre <- df_long_KOOS %>% filter(KOOS_group == "Pre")
df_long_KOOS_Post12 <- df_long_KOOS %>% filter(KOOS_group == "Post12")

custom_labels_KOOS_pre <- c(
  "Pre__KOOS_Pain" = "KOOS Pain",
  "Pre__KOOS_Symptom" = "KOOS Symptom",
  "Pre__KOOS_QOL" = "KOOS QOL",
  "Pre__KOOS_ADL" = "KOOS ADL",
  "Pre__KOOS_Sports" = "KOOS Sports"
)
custom_labels_KOOS_post12 <- c(
  "Post_12Month__KOOS_Pain" = "KOOS Pain",
  "Post_12Month__KOOS_Symptom" = "KOOS Symptom",
  "Post_12Month__KOOS_QOL" = "KOOS QOL",
  "Post_12Month__KOOS_ADL" = "KOOS ADL",
  "Post_12Month__KOOS_Sports" = "KOOS Sports"
)

plot_KOOS_Pre <- ggplot(df_long_KOOS_pre, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0, 110)) +
  facet_wrap(~ KOOS, scales = "fixed", labeller = labeller(KOOS = custom_labels_KOOS_pre), strip.position = "bottom", ncol = 5) + 
  labs(x = "Pre-Operation", y = "Score") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#DD7670", "3" = "#7F1212")) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
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
    legend.margin = margin(0, 0, 0, 0))

ggsave("99_Fig/fig3/cluster_KOOS_pre.png", plot = plot_KOOS_Pre , width = 4, height = 1.25, bg = "white")
ggsave("99_Fig/fig3/cluster_KOOS_pre.pdf", plot = plot_KOOS_Pre , width = 4, height = 1.25, bg = "white")

plot_KOOS_Post12 <- ggplot(df_long_KOOS_Post12, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0, 110)) +
  facet_wrap(~ KOOS, scales = "fixed", labeller = labeller(KOOS = custom_labels_KOOS_post12), strip.position = "bottom", ncol = 5) + 
  labs(x = "Post-Operation 12-Month", y = "Score") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#DD7670", "3" = "#7F1212")) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
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
    legend.margin = margin(0, 0, 0, 0))

ggsave("99_Fig/fig3/cluster_KOOS_post12.png", plot = plot_KOOS_Post12, width = 4, height = 1.25)
ggsave("99_Fig/fig3/cluster_KOOS_post12.pdf", plot = plot_KOOS_Post12, width = 4, height = 1.25)


df_long_Krenn <- df_OA %>% 
  pivot_longer(
    cols = 37:39,
    names_to = "Krenn",
    values_to = "Value"
  )

plot_krenn_synovitis <- ggplot(df_long_Krenn, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 4.5, 1)) +
  coord_cartesian(ylim = c(0, 4.5)) +
  facet_wrap(~ Krenn, scales = "fixed", strip.position = "bottom", ncol = 3) +
  labs(x = "Krenn synovitis score",  y = "Score") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#DD7670", "3" = "#7F1212")) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    strip.text = element_text(size = 5, margin = margin(t = 0)),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, , margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0))

ggsave("99_Fig/fig3/cluster_krenn_synovitis.png", plot = plot_krenn_synovitis, width = 2.75, height = 1.25)
ggsave("99_Fig/fig3/cluster_Krenn_synovitis.pdf", plot = plot_krenn_synovitis, width = 2.75, height = 1.25)


df_KL <- df_OA %>%
  filter(!is.na(KL_grade)) %>%            
  count(cluster, KL_grade)

df_KL$KL_grade <- factor(df_KL$KL_grade, levels = c("4", "3"))

plot_KL <- ggplot(df_KL, aes(x = cluster, y = n, fill = KL_grade)) +
  geom_bar(stat = "identity", position = position_fill(reverse = FALSE), width = 0.9, color = "black", linewidth = 0.1) +
  labs(y = "KL grade\nratio %", fill = "KL_grade") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = c("3" = "white", "4" = "grey")) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 5),
    axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
    axis.text.y  = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 2)),
    legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, -3, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )

ggsave("99_Fig/fig3/cluster_KL.png", plot = plot_KL, width = 1.3, height = 1)
ggsave("99_Fig/fig3/cluster_KL.pdf", plot = plot_KL, width = 1.3, height = 1)



# # For response to reviewers
# library(Rtsne)
# library(ggplot2)
# library(dunn.test)
# library(dplyr)
# library(tidyverse)
# library(compositions)
# library(emmeans)
# library(broom)
# 
# df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
# score <- read.delim("01_formatted/annotations_full_OARSI_Krenn.txt", sep = "\t", row.names = 1)
# df[,37:41]  <- score[,1:5]
# df <- df[df$Joint == "Knee" & df$Diagnosis == "OA",]
# keep <- rownames(df)
# 
# Degree_df <- read.csv("00_src/Degree_df_241220.csv", sep = ",", header = TRUE, row.names = 1)
# Degree_df <- Degree_df[keep,]
# 
# entropy <- read.delim(file = "00_src/Entropy_group.txt", sep = "\t", row.names = 1)
# entropy <- as.data.frame(entropy[intersect(rownames(df), rownames(entropy)),])
# 
# 
# pca_result <- prcomp(Degree_df[,1:78], scale. = FALSE)
# 
# plot(pca_result, type = "l")
# 
# summary(pca_result)
# # Importance of components:
# #                            PC1     PC2     PC3    PC4     PC5     PC6     PC7      PC8      PC9     PC10     PC11     PC12     PC13     PC14
# # Standard deviation     0.08863 0.07383 0.03375 0.0266 0.01942 0.01638 0.01218 0.009805 0.008287 0.007295 0.006112 0.005035 0.004564 0.004338
# # Proportion of Variance 0.48071 0.33359 0.06972 0.0433 0.02309 0.01642 0.00907 0.005880 0.004200 0.003260 0.002290 0.001550 0.001270 0.001150
# # Cumulative Proportion  0.48071 0.81430 0.88401 0.9273 0.95041 0.96683 0.97590 0.981780 0.985990 0.989240
# # ...
# 
# pca_data <- data.frame(pca_result$x)
# 
# 
# # # silhouette
# # library(cluster)
# # set.seed(123)
# # silhouette_results <- lapply(2:8, function(k) {
# #   km <- kmeans(pca_data[,1:3], centers = k)
# #   sil <- silhouette(km$cluster, dist(pca_data[,1:3]))
# #   data.frame(
# #     k = k,
# #     silhouette_score = mean(sil[,"sil_width"])
# #   )
# # })
# # 
# # silhouette_summary <- bind_rows(silhouette_results)
# # silhouette_summary
# # k silhouette_score
# # 1 2        0.3575037
# # 2 3        0.3885601
# # 3 4        0.3496392
# # 4 5        0.3280046
# # 5 6        0.3382076
# # 6 7        0.3102634
# # 7 8        0.3140369
# # 
# # write.table(silhouette_summary, "03_res/Clustering/tissue_conection_cluster_k.txt", sep = "\t", row.names = FALSE)
# # 
# # plot_silhouette_k <- ggplot(silhouette_summary, aes(x = k, y = silhouette_score)) +
# #   geom_line(linewidth = 0.1) +
# #   geom_point(size = 0.5) +
# #   scale_x_continuous(breaks = 2:8) +
# #   labs(x = "Number of clusters", y = "Average silhouette score") +
# #   theme_classic() +
# #   theme(
# #     panel.border = element_blank(),
# #     axis.ticks = element_line(linewidth = 0.1),
# #     axis.line = element_line(linewidth = 0.1),
# #     axis.title.x = element_text(size = 6),
# #     axis.text.x = element_text(size = 5),
# #     axis.title.y = element_text(size = 6),
# #     axis.text.y = element_text(size = 5)
# #   )
# 
# 
# # tSNE
# set.seed(123)
# tsne_res <- Rtsne(pca_data[,1:3], dims= 2, perplexity = 26)
# tsne_data <- data.frame(tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2])
# rownames(tsne_data) <- rownames(Degree_df)
# 
# 
# # clustering
# kmeans_result <- kmeans(pca_data[,1:3], centers = 3)
# 
# tsne_data$cluster <- as.factor(kmeans_result$cluster)
# tsne_data$Diagnosis <- df$Diagnosis
# 
# tsne_data$Diagnosis_plot <- recode(
#   tsne_data$Diagnosis,
#   "SLE" = "Other autoimmune diseases",
#   "SSc" = "Other autoimmune diseases"
# )
# 
# plot_tsne <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = cluster, shape = Diagnosis_plot)) + 
#   geom_point(size = 0.5) +
#   scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#DD7670", "3" = "#7F1212")) +
#   scale_shape_manual(values = c("nonOA" = 16, "RA" = 15, "OA" = 17, "Other autoimmune diseases" = 18),
#                      labels = c("nonOA" = "non-OA structural", "RA" = "RA", "OA" = "OA", "Other autoimmune diseases" = "Other autoimmune diseases")) +
#   labs(x = "tSNE1", y = "tSNE2", shape = "Diagnosis", color = "kmeans cluster") +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
#     axis.ticks = element_blank(),
#     axis.line = element_blank(),
#     axis.title = element_text(size = 5), 
#     axis.text = element_blank(),
#     legend.text = element_text(size = 5, margin = margin(l = 0)),
#     legend.title = element_text(size = 5),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_blank(), 
#     legend.key.size = unit(0.4, "lines"),
#     legend.margin = margin(0, 0, 0, 0),
#     legend.box.margin = margin(0, 0, 0, -8)
#   )
# 
# ggsave("99_Fig/fig3/Revise/tSNE.png", plot = plot_tsne, width = 1.75, height = 1.25)
# ggsave("99_Fig/fig3/Revise/tSNE.pdf", plot = plot_tsne, width = 1.75, height = 1.25)
# 
# df[,42:44] <- tsne_data[,1:3]
# 
# # write.table(df, "01_formatted/Revise/annotations_full_tissue_conection_cluster.txt", sep = "\t", row.names = TRUE, col.names = NA)
# 
# 
# # Proportion
# df_stacked_bar <- df
# df_stacked_bar$Entropy <- entropy$`entropy[intersect(rownames(df), rownames(entropy)), ]`
# 
# df_stacked_bar <- df_stacked_bar %>% 
#   count(cluster, Entropy)
# 
# df_stacked_bar$Entropy <- factor(df_stacked_bar$Entropy, levels = c("Low", "Mid", "High"))
# 
# plot_stacked_bar <- ggplot(df_stacked_bar, aes(x = cluster, y = n, fill = Entropy)) +
#   geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.9) +
#   labs(y = "Fraction of entropy type", fill = "Entropy") +
#   scale_y_continuous(labels = function(x) x * 100) +
#   scale_fill_manual(values = c("Low" = "#bfe6bf", "Mid" = "#66cc66", "High" = "#336633")) +
#   theme_classic() +
#   theme(
#     plot.title = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.x = element_line(linewidth = 0.1),
#     axis.ticks.y = element_line(linewidth = 0.1),
#     axis.line = element_line(linewidth = 0.1),
#     axis.title.x = element_blank(),
#     axis.text.x  = element_text(size = 5),
#     axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
#     axis.text.y  = element_text(size = 5),
#     legend.title = element_text(size = 5, margin = margin(b = 2)),
#     legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
#     legend.key.size = unit(2, "mm"),
#     legend.margin = margin(t = -1, b = -1, unit = "mm"),
#     legend.box.margin = margin(0, 0, 0, -3, unit = "mm"),
#     plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
#   )
# 
# ggsave("99_Fig/fig3/Revise/tSNE_cluster_entropy_proportion.png", plot = plot_stacked_bar, width = 1.25, height = 1.25)
# ggsave("99_Fig/fig3/Revise/tSNE_cluster_entropy_proportion.pdf", plot = plot_stacked_bar, width = 1.25, height = 1.25)
# 
# 
# # Statistics
# df_OA <- df
# 
# scores <- colnames(df_OA[,18:32])
# scores_2 <- colnames(df_OA[,c(35,36,41,37:40)])
# 
# ## The same applies to "score_2".
# results <- list()
# for (score in scores) {
#   kw_result <- kruskal.test(as.formula(paste(score, "~ cluster")), data = df_OA)
#   dunn_result <- dunn.test(df_OA[[score]], df_OA$cluster, method = "bh")
#   significant_idx <- which(dunn_result$P.adjusted < 0.1)
#   significant_comparisons <- dunn_result$comparisons[significant_idx]
#   significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
#   significant_z <- dunn_result$Z[significant_idx]
#   results[[score]] <- list(
#     test_type = "Kruskal-Wallis",
#     main_test = kw_result,
#     posthoc_test = dunn_result,
#     significant_comparisons = significant_comparisons,
#     significant_p_adjusted = significant_p_adjusted,
#     significant_z = significant_z 
#   )
# }
# 
# for (score in scores) {
#   cat("\n--- Results for", score, "---\n")
#   print(results[[score]]$main_test)
#     if (length(results[[score]]$significant_comparisons) > 0) {
#     cat("\n** Significant Dunn Post-hoc Test Results for", score, "**\n")
#     cat("=============================================\n")
#     for (i in seq_along(results[[score]]$significant_comparisons)) {
#       cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
#                   results[[score]]$significant_comparisons[i], 
#                   results[[score]]$significant_p_adjusted[i],
#                   results[[score]]$significant_z[i]))
#     }
#     cat("=============================================\n")
#   }
# }
# 
# # Notes: Significant differences are shown.
# # ** Significant Dunn Post-hoc Test Results for Post_3Month__KOOS_QOL **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.09335 - Z-value: -1.86463
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_Symptom **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.06007 - Z-value: 1.75017
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_QOL **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.01374 - Z-value: 2.35912
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.02697 - Z-value: 2.36608
# # =============================================
# # 
# # ** Significant Dunn Post-hoc Test Results for Synovitis_score **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.07587 - Z-value: 1.95501
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Krenn_Lining **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 5.22605
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00955 - Z-value: 2.49101
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.02061 - Z-value: -2.04140
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Krenn_Stroma **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.07170 - Z-value: 1.46322
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.03798 - Z-value: -1.95447
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00323 - Z-value: -3.06826
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Krenn_Infiltrate **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00508 - Z-value: 2.93019
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00426 - Z-value: 2.76546
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Total **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00001 - Z-value: 4.61097
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.03682 - Z-value: 1.78881
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.02160 - Z-value: -2.18621
# # =============================================
# 
# 
# Degree_df <- t(Degree_df)
# 
# conversion <- c(
#   "Fibro.dense.irregular." = "Fibrous tissue (dense, irregular)", 
#   "Fibro.dense.regular." = "Fibrous tissue (dense, regular)", 
#   "Fibro.loose." = "Fibrous tissue (loose)", 
#   "Immune.cells" = "TLS", 
#   "adipose" = "Adipose", 
#   "lining" = "Lining",
#   "muscle" = "Muscle",
#   "plasma" = "Plasma",
#   "RBC" = "RBC",
#   "Stroma" = "Stroma",
#   "vessel" = "Micro vessel",
#   "vessel.large." = "Large vessel"
# )
# 
# original_names <- rownames(Degree_df)
# new_names <- sapply(original_names, function(name) {
#   parts <- unlist(strsplit(name, "_")) 
#   converted_parts <- sapply(parts, function(part) {
#     if (part %in% names(conversion)) {
#       conversion[part] 
#     } else {
#       part 
#     }
#   })
#   paste(converted_parts, collapse = "_")
# })
# rownames(Degree_df) <- new_names
# 
# connection_select <- c("Lining_Micro vessel","Fibrous tissue (loose)_Lining","Fibrous tissue (loose)_Micro vessel",
#                        "Micro vessel_Adipose",
#                        "Large vessel_Adipose",
#                        "Fibrous tissue (dense, irregular)_Micro vessel",
#                        "TLS_TLS","Stroma_Stroma")
# Degree_df_select <- Degree_df %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   select(all_of(connection_select))
# Degree_df_select <- Degree_df_select[rownames(df), ]
# Degree_df_select[,c(9,10)] <- df[,c(1,44)]
# Degree_df_select <- Degree_df_select[Degree_df_select$Diagnosis == "OA",]
# 
# # Statistics
# connections <- colnames(Degree_df_select[,1:8])
# 
# results <- list()
# for (connection in connections) {
#   kw_result <- kruskal.test(as.formula(paste0("`", connection, "` ~ cluster")), data = Degree_df_select)
#   dunn_result <- dunn.test(Degree_df_select[[connection]], Degree_df_select$cluster, method = "bh")
#   significant_idx <- which(dunn_result$P.adjusted < 0.1)
#   significant_comparisons <- dunn_result$comparisons[significant_idx]
#   significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
#   significant_z <- dunn_result$Z[significant_idx]
#   results[[connection]] <- list(
#     test_type = "Kruskal-Wallis",
#     main_test = kw_result,
#     posthoc_test = dunn_result,
#     significant_comparisons = significant_comparisons,
#     significant_p_adjusted = significant_p_adjusted,
#     significant_z = significant_z 
#   )
# }
# 
# for (connection in connections) {
#   cat("\n--- Results for", connection, "---\n")
#   print(results[[connection]]$main_test)
#   if (length(results[[connection]]$significant_comparisons) > 0) {
#     cat("\n** Significant Dunn Post-hoc Test Results for", connection, "**\n")
#     cat("=============================================\n")
#     for (i in seq_along(results[[connection]]$significant_comparisons)) {
#       cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
#                   results[[connection]]$significant_comparisons[i], 
#                   results[[connection]]$significant_p_adjusted[i],
#                   results[[connection]]$significant_z[i]))
#     }
#     cat("=============================================\n")
#   }
# }
# 
# # Notes: Significant differences are shown.
# # ** Significant Dunn Post-hoc Test Results for Lining_Micro vessel **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 6.42233
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00537 - Z-value: 2.55077
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00210 - Z-value: -2.98925
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Fibrous tissue (loose)_Lining **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 6.84748
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: 4.91107
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Fibrous tissue (loose)_Micro vessel **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 6.88249
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: 4.90279
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Micro vessel_Adipose **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00024 - Z-value: -3.60273
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00087 - Z-value: 3.13050
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 5.97130
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Large vessel_Adipose **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00003 - Z-value: -4.11887
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00002 - Z-value: 4.38405
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Fibrous tissue (dense, irregular)_Micro vessel **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.68541
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: -4.91990
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for TLS_TLS **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00399 - Z-value: 3.00425
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.01241 - Z-value: 2.39651
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Stroma_Stroma **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00384 - Z-value: -2.79971
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00035 - Z-value: -3.68028
# # =============================================
# 
# 
# Degree_df_long <- Degree_df_select %>%
#   pivot_longer(cols = 1:8, names_to = "Connection", values_to = "Value")
# 
# Degree_df_long$Connection <- factor(Degree_df_long$Connection, levels = connection_select)
# Degree_df_long$cluster <- factor(Degree_df_long$cluster)
# 
# plot_tissue <- ggplot(Degree_df_long, aes(x = cluster, y = Value, color = cluster)) +
#   geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
#   geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
#   facet_wrap(~ Connection, scales = "free", strip.position = "bottom", ncol = 4) + 
#   labs(y = "Degree of connections") +
#   scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#DD7670", "3" = "#7F1212")) +
#   theme_classic() +
#   theme(
#     panel.border = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_line(linewidth = 0.1),
#     axis.line = element_line(linewidth = 0.1),
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     strip.text = element_text(size = 5, margin = margin(t = 0)),
#     strip.background = element_blank(),
#     strip.placement  = "outside",
#     axis.title.y = element_text(size = 6),
#     axis.text.y = element_text(size = 5),
#     legend.title = element_text(size = 5),
#     legend.text = element_text(size = 5, , margin = margin(l = -5)),
#     legend.key.height = unit(0.25, "cm"),
#     legend.box.margin = margin(0, 0, 0, 0),
#     legend.margin     = margin(0, 0, 0, 0))
# 
# ggsave("99_Fig/fig3/Revise/tSNE_cluster_connection.png", plot = plot_tissue, width = 4, height = 1.75)
# ggsave("99_Fig/fig3/Revise/tSNE_cluster_connrction.pdf", plot = plot_tissue, width = 4, height = 1.75)

