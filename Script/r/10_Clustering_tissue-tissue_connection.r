library(Rtsne)
library(ggplot2)
library(dunn.test)
library(dplyr)
library(tidyverse)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
score <- read.delim("00_src/annotations_full_OARSI_Krenn.txt", sep = "\t", row.names = 1)
df[,37:41]  <- score[,1:5]
df <- df[df$Joint == "Knee",]

Degree_df <- read.csv("00_src/Degree_df_241220.csv", sep = ",", header = TRUE, row.names = 1)

entropy <- read.delim(file = "00_src/Entropy_group.txt", sep = "\t", row.names = 1)
entropy <- as.data.frame(entropy[intersect(rownames(df), rownames(entropy)),])


# tSNE
set.seed(123)
tsne_res <- Rtsne(Degree_df[,1:78], dims= 2)
tsne_data <- data.frame(tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2])
rownames(tsne_data) <- rownames(Degree_df)


# clustering
set.seed(123)
kmeans_result <- kmeans(tsne_data, centers = 4)

tsne_data$cluster <- as.factor(kmeans_result$cluster)
tsne_data$Diagnosis <- df$Diagnosis


plot_tsne <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = cluster, shape = Diagnosis)) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#F3B2A6", "3" = "#C73A3A", "4" = "#7F1212")) +
  scale_shape_manual(values = c("nonOA" = 16, "RA" = 15, "OA" = 17, "SLE" = 18),
                     labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA", "SLE" = "SLE")) +
  labs(x = "tSNE1", y = "tSNE2", shape = "Diagnosis", color = "kmeans cluster") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size = 6), 
    axis.text = element_blank(),
    legend.text = element_text(size = 5, margin = margin(l = 0)),
    legend.title = element_text(size = 5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), 
    legend.key.size = unit(0.4, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -8)
  )

ggsave("99_Fig/fig4/tSNE.png", plot = plot_tsne, width = 2, height = 1.5)
ggsave("99_Fig/fig4/tSNE.pdf", plot = plot_tsne, width = 2, height = 1.5)

df[,42:44] <- tsne_data[,1:3]


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

ggsave("99_Fig/fig4/tSNE_cluster_entropy_proportion.png", plot = plot_stacked_bar, width = 1.25, height = 1.25)
ggsave("99_Fig/fig4/tSNE_cluster_entropy_proportion.pdf", plot = plot_stacked_bar, width = 1.25, height = 1.25)


# Statistics
df_OA <- df[df$Diagnosis == "OA",]

scores <- colnames(df_OA[,18:32])
scores_2 <- colnames(df_OA[,c(35,36,41,37:40)])
tissue <- colnames(df_OA[,6:17])

## The same applies to "score_2" and "tissue".
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

cat("\n========== Kruskal-Wallis Results ==========\n")
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

# ###
# ** Significant Dunn Post-hoc Test Results for Pre__KOOS_Pain **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.06592 - Z-value: 2.29083
# =============================================
# ** Significant Dunn Post-hoc Test Results for Pre__KOOS_Symptom **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.03192 - Z-value: -2.30295
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08003 - Z-value: -1.75051
# >> Comparison: 3 - 4      - Adjusted P-value: 0.01688 - Z-value: 2.76876
# =============================================
# ** Significant Dunn Post-hoc Test Results for Pre__KOOS_ADL **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.09952 - Z-value: 2.12999
# =============================================
# ** Significant Dunn Post-hoc Test Results for Pre__KOOS_QOL **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.04230 - Z-value: 2.45468
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_3Month__KOOS_Sports **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08043 - Z-value: -2.21429
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_3Month__KOOS_QOL **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02606 - Z-value: -2.62410
# >> Comparison: 3 - 4      - Adjusted P-value: 0.09469 - Z-value: 1.85830
# =============================================
# ###
# ** Significant Dunn Post-hoc Test Results for Synovitis_score **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.03687 - Z-value: -2.24795
# >> Comparison: 2 - 3      - Adjusted P-value: 0.03884 - Z-value: 2.48525
# >> Comparison: 2 - 4      - Adjusted P-value: 0.03972 - Z-value: 2.05663
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.06834 - Z-value: -1.68953
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00462 - Z-value: 2.83254
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00002 - Z-value: 4.51048
# >> Comparison: 1 - 4      - Adjusted P-value: 0.07702 - Z-value: 1.42540
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00194 - Z-value: 3.21739
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07428 - Z-value: -1.53901
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Infiltrate **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.04106 - Z-value: -2.20613
# >> Comparison: 2 - 3      - Adjusted P-value: 0.04434 - Z-value: 2.43769
# >> Comparison: 3 - 4      - Adjusted P-value: 0.08326 - Z-value: -1.73210
# =============================================
# ** Significant Dunn Post-hoc Test Results for Total **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.03037 - Z-value: 2.32178
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00127 - Z-value: 3.52527
# >> Comparison: 2 - 4      - Adjusted P-value: 0.06341 - Z-value: 1.72492
# >> Comparison: 3 - 4      - Adjusted P-value: 0.06128 - Z-value: -1.87144
# =============================================
# ###
# ** Significant Dunn Post-hoc Test Results for Adipose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01294 - Z-value: -2.29800
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.84883
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00003 - Z-value: -4.12995
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00000 - Z-value: -6.17628
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00002 - Z-value: -4.32226
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 6.64320
# >> Comparison: 1 - 3      - Adjusted P-value: 0.07784 - Z-value: 1.51514
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00001 - Z-value: -4.34779
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00000 - Z-value: 5.62518
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00025 - Z-value: 3.58963
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 4.93057
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00170 - Z-value: -3.25521
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00611 - Z-value: 2.74202
# >> Comparison: 2 - 4      - Adjusted P-value: 0.03482 - Z-value: -1.99144
# >> Comparison: 3 - 4      - Adjusted P-value: 0.09862 - Z-value: 1.39055
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: -5.26862
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.16966
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00002 - Z-value: 4.21687
# >> Comparison: 3 - 4      - Adjusted P-value: 0.01782 - Z-value: -2.26095
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.07562 - Z-value: -1.77666
# >> Comparison: 1 - 3      - Adjusted P-value: 0.09697 - Z-value: 1.51688
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00417 - Z-value: 3.19648
# >> Comparison: 3 - 4      - Adjusted P-value: 0.02143 - Z-value: -2.45001
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00010 - Z-value: -4.14294
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00006 - Z-value: 4.11134
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00041 - Z-value: -3.45820
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00044 - Z-value: -3.51645
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00413 - Z-value: 2.99376
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00639 - Z-value: 3.07131
# >> Comparison: 1 - 4      - Adjusted P-value: 0.05568 - Z-value: 1.91352
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01516 - Z-value: -2.32228
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00905 - Z-value: 2.61032
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 4.84218
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00062 - Z-value: 3.53251
# >> Comparison: 3 - 4      - Adjusted P-value: 0.06768 - Z-value: -1.58571
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: 1 - 4      - Adjusted P-value: 0.07559 - Z-value: 1.95664
# =============================================
# ** Significant Dunn Post-hoc Test Results for RBC **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.03547 - Z-value: 2.51737
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00017 - Z-value: 3.86321
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00006 - Z-value: 4.25628
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00440 - Z-value: 2.75547
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00168 - Z-value: 3.14191
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00562 - Z-value: -2.89871
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00150 - Z-value: -3.48123
# >> Comparison: 3 - 4      - Adjusted P-value: 0.01646 - Z-value: 2.39861
# =============================================


## Plot
df_OA[, 6:17] <- log10(df_OA[, 6:17] + 1e-6)

df_long <- df_OA %>%
  pivot_longer(cols = 6:17, names_to = "Tissue", values_to = "Value")

df_long$Tissue <- factor(df_long$Tissue, levels = c("Fibrous_tissue__dense_irregular", "Fibrous_tissue__dense_regular", "Fibrous_tissue__loose", "Lining", "TLS", "Plasma", "Stroma", "RBC", "Micro_vessel", "Large_vessel", "Adipose", "Muscle"))

custom_labels <- c(
  "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense.irregular)",
  "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense.regular)",
  "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
  "Micro_vessel" = "Micro vessel",
  "Large_vessel" = "Large vessel"
)

plot_tissue <- ggplot(df_long, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels), strip.position = "bottom", ncol = 6) + 
  labs(y = "log10(Proportion)") +
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

ggsave("99_Fig/fig4/tSNE_cluster_tissue_proportion.png", plot = plot_tissue, width = 6, height = 2)
ggsave("99_Fig/fig4/tSNE_cluster_tissue_proportion.pdf", plot = plot_tissue, width = 6, height = 2)


df_long_KOOOS <- df_OA %>%
  pivot_longer(cols = 18:32, names_to = "KOOS", values_to = "Value") %>%
  mutate(KOOS_group = case_when(
    str_detect(KOOS, "Pre") ~ "Pre",
    str_detect(KOOS, "Post_12") ~ "Post12",
    TRUE ~ "Other"
  )
  )

df_long_KOOOS_pre <- df_long_KOOOS %>% filter(KOOS_group == "Pre")
df_long_KOOOS_Post12 <- df_long_KOOOS %>% filter(KOOS_group == "Post12")

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

plot_KOOS_Pre <- ggplot(df_long_KOOOS_pre, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0, 110)) +
  facet_wrap(~ KOOS, scales = "fixed", labeller = labeller(KOOS = custom_labels_KOOS_pre), strip.position = "bottom", ncol = 5) + 
  labs(x = "Pre-Operation", y = "Score") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#F3B2A6", "3" = "#C73A3A", "4" = "#7F1212")) +
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

ggsave("99_Fig/fig4/tSNE_cluster_KOOS_pre.png", plot = plot_KOOS_Pre , width = 5, height = 1.25, bg = "white")
ggsave("99_Fig/fig4/tSNE_cluster_KOOS_pre.pdf", plot = plot_KOOS_Pre , width = 5, height = 1.25, bg = "white")

plot_KOOS_Post12 <- ggplot(df_long_KOOOS_Post12, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0, 110)) +
  facet_wrap(~ KOOS, scales = "fixed", labeller = labeller(KOOS = custom_labels_KOOS_post12), strip.position = "bottom", ncol = 5) + 
  labs(x = "Post-Operation 12-Month", y = "Score") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#F3B2A6", "3" = "#C73A3A", "4" = "#7F1212")) +
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

ggsave("99_Fig/fig4/tSNE_cluster_KOOS_post12.png", plot = plot_KOOS_Post12, width = 5, height = 1.25)
ggsave("99_Fig/fig4/tSNE_cluster_KOOS_post12.pdf", plot = plot_KOOS_Post12, width = 5, height = 1.25)


df_long_score <- df_OA %>%
  pivot_longer(cols = c(35,36,40,41), names_to = "Score", values_to = "Value")

df_long_synovitis <- df_long_score %>% filter(Score == "Synovitis_score")
df_long_Krenn <- df_long_score %>% filter(Score == "Total")
df_long_OARSI <- df_long_score %>% filter(Score == "OARSI_score")

plot_synovitis <- ggplot(df_long_synovitis, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 4, 1)) +
  coord_cartesian(ylim = c(0, 4.5)) +
  labs(x = "Synovitis score",  y = "Score") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#F3B2A6", "3" = "#C73A3A", "4" = "#7F1212")) +
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
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, , margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0))

ggsave("99_Fig/fig4/tSNE_cluster_synovitis.png", plot = plot_synovitis, width = 1.5, height = 1.25)
ggsave("99_Fig/fig4/tSNE_cluster_synovitis.pdf", plot = plot_synovitis, width = 1.5, height = 1.25)

plot_krenn_synovitis <- ggplot(df_long_Krenn, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 8.5, 2)) +
  coord_cartesian(ylim = c(0, 8.5)) +
  labs(x = "Krenn synovitis score",  y = "Score") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#F3B2A6", "3" = "#C73A3A", "4" = "#7F1212")) +
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
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, , margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0))

ggsave("99_Fig/fig4/tSNE_cluster_krenn_synovitis.png", plot = plot_krenn_synovitis, width = 1.5, height = 1.25)
ggsave("99_Fig/fig4/tSNE_cluster_Krenn_synovitis.pdf", plot = plot_krenn_synovitis, width = 1.5, height = 1.25)

plot_oarsi <- ggplot(df_long_OARSI, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 7, 2)) +
  coord_cartesian(ylim = c(0, 7)) +
  labs(x = "OARSI score",  y = "Score") +
  scale_color_manual(values = c("1" = "#D8C6C2", "2" = "#F3B2A6", "3" = "#C73A3A", "4" = "#7F1212")) +
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
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, , margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0))

ggsave("99_Fig/fig4/tSNE_cluster_OARSI.png", plot = plot_oarsi, width = 1.5, height = 1.25)
ggsave("99_Fig/fig4/tSNE_cluster_OARSI.pdf", plot = plot_oarsi, width = 1.5, height = 1.25)


df_OARSI <- df_OA %>%
  filter(!is.na(OARSI_score)) %>%            
  count(cluster, OARSI_score)

df_OARSI$OARSI_score <- factor(df_OARSI$OARSI_score, levels = c("6", "5", "4"))

plot_oarsi_stacked <- ggplot(df_OARSI, aes(x = cluster, y = n, fill = OARSI_score)) +
  geom_bar(stat = "identity", position = position_fill(reverse = FALSE), width = 0.9, color = "black", linewidth = 0.1) +
  labs(y = "OARSI score\nratio %", fill = "OARSI score") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = c("4" = "white", "5" = "grey", "6" = "black")) +
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

ggsave("99_Fig/fig4/tSNE_cluster_OARSI_stackedbar.png", plot = plot_oarsi_stacked, width = 1.5, height = 1.25)
ggsave("99_Fig/fig4/tSNE_cluster_OARSI_stackedbar.pdf", plot = plot_oarsi_stacked, width = 1.5, height = 1.25)


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

ggsave("99_Fig/fig4/tSNE_cluster_KL.png", plot = plot_KL, width = 1.5, height = 1.25)
ggsave("99_Fig/fig4/tSNE_cluster_KL.pdf", plot = plot_KL, width = 1.5, height = 1.25)


write.table(df, "01_formatted/annotations_full_tissue_conection_cluster.txt", sep = "\t", row.names = TRUE, col.names = NA)


sessionInfo()
# R version 4.3.3 (2024-02-29 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
# 
# 
# locale:
#  [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8    LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                    LC_TIME=Japanese_Japan.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    tidyverse_2.0.0 dplyr_1.1.4     dunn.test_1.3.6 ggplot2_3.5.1   Rtsne_0.17     
# 
# loaded via a namespace (and not attached):
#  [1] vctrs_0.6.5       cli_3.6.4         rlang_1.1.5       stringi_1.8.7     generics_0.1.3    glue_1.8.0        colorspace_2.1-1  hms_1.1.3         scales_1.3.0      grid_4.3.3        munsell_0.5.1     tzdb_0.5.0        lifecycle_1.0.4   compiler_4.3.3    timechange_0.3.0 
# [16] Rcpp_1.0.14       pkgconfig_2.0.3   rstudioapi_0.17.1 R6_2.6.1          tidyselect_1.2.1  pillar_1.10.1     magrittr_2.0.3    tools_4.3.3       withr_3.0.2       gtable_0.3.6 

