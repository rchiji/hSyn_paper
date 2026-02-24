library(Rtsne)
library(ggplot2)
library(dunn.test)
library(dplyr)
library(tidyverse)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
score <- read.delim("00_src/annotations_full_OARSI_Krenn.txt", sep = "\t", row.names = 1)
df[,37:41]  <- score[,1:5]
df <- df[df$Joint == "Knee",]


# tSNE
set.seed(123)
tsne_res <- Rtsne(df[,6:17], dims= 2, perplexity = 36)
tsne_data <- data.frame(tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2])
rownames(tsne_data) <- rownames(df)


# clustering
set.seed(123)
kmeans_result <- kmeans(tsne_data, centers = 4)

tsne_data$cluster <- as.factor(kmeans_result$cluster)
tsne_data$Diagnosis <- df$Diagnosis

plot_tsne <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = cluster, shape = Diagnosis)) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("1" = "#C2BAB4", "2" = "#FEC089", "3" = "#F06A00", "4" = "#8C2D04")) +
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

ggsave("99_Fig/fig3/tSNE.png", plot = plot_tsne, width = 2, height = 1.5)
ggsave("99_Fig/fig3/tSNE.pdf", plot = plot_tsne, width = 2, height = 1.5)

df[,42:44] <- tsne_data[,1:3]


# Proportion
df_stacked_bar <- df %>% 
  count(cluster, Diagnosis)

df_stacked_bar$Diagnosis <- factor(df_stacked_bar$Diagnosis, levels = c("OA", "RA", "SLE", "nonOA"))

plot_stacked_bar <- ggplot(df_stacked_bar, aes(x = cluster, y = n, fill = Diagnosis)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.9) +
  labs(y = "Fraction of diagnosis type", fill = "Diagnosis") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = c("nonOA" = "#1f77b4", "RA" = "#2ca02c", "OA" = "#ff7f0e", "SLE" = "grey80"),
                    labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA", "SLE" = "SLE")) +
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

ggsave("99_Fig/fig3/tSNE_cluster_diagnosis_proportion.png", plot = plot_stacked_bar, width = 2.3, height = 1.6)
ggsave("99_Fig/fig3/tSNE_cluster_diagnosis_proportion.pdf", plot = plot_stacked_bar, width = 2.3, height = 1.6)


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

# ### KOOS *************************************************************************
# ** Significant Dunn Post-hoc Test Results for Post_3Month__KOOS_QOL **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.05870 - Z-value: -1.89041
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08454 - Z-value: 1.90822
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_Pain **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.03564 - Z-value: -2.51570
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_Symptom **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.02594 - Z-value: 2.38034
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02458 - Z-value: 2.24794
# >> Comparison: 3 - 4      - Adjusted P-value: 0.03405 - Z-value: -2.53174
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_QOL **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.08595 - Z-value: -1.71719
# =============================================
# ### Cartilage annd Synovium score ************************************************
# ** Significant Dunn Post-hoc Test Results for Synovitis_score **
# =============================================
# >> Comparison: 2 - 4      - Adjusted P-value: 0.06816 - Z-value: -2.00059
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07012 - Z-value: -2.26730
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00047 - Z-value: 3.77827
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00684 - Z-value: 2.83648
# >> Comparison: 2 - 4      - Adjusted P-value: 0.02230 - Z-value: -2.28525
# >> Comparison: 3 - 4      - Adjusted P-value: 0.08541 - Z-value: -1.58098
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Stroma **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01186 - Z-value: -2.41330
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01456 - Z-value: -2.44303
# >> Comparison: 1 - 4      - Adjusted P-value: 0.01666 - Z-value: -2.53940
# >> Comparison: 2 - 4      - Adjusted P-value: 0.02999 - Z-value: -2.57600
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Infiltrate **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.09372 - Z-value: 1.67610
# >> Comparison: 1 - 3      - Adjusted P-value: 0.02724 - Z-value: 2.60902
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07354 - Z-value: -1.96836
# =============================================
# ** Significant Dunn Post-hoc Test Results for Total **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01029 - Z-value: 2.70371
# >> Comparison: 1 - 3      - Adjusted P-value: 0.08015 - Z-value: 1.61242
# >> Comparison: 2 - 4      - Adjusted P-value: 0.01813 - Z-value: -2.74544
# >> Comparison: 3 - 4      - Adjusted P-value: 0.09105 - Z-value: -1.68989
# =============================================
# ### Tissue *************************************************************************
# ** Significant Dunn Post-hoc Test Results for Adipose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: -6.31198
# >> Comparison: 1 - 3      - Adjusted P-value: 0.06782 - Z-value: 1.49224
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.98089
# >> Comparison: 1 - 4      - Adjusted P-value: 0.01125 - Z-value: -2.35039
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00021 - Z-value: 3.70980
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00033 - Z-value: -3.51383
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01018 - Z-value: -2.46843
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.97168
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00011 - Z-value: -3.87752
# >> Comparison: 1 - 4      - Adjusted P-value: 0.03152 - Z-value: -1.93876
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00006 - Z-value: 4.10804
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01351 - Z-value: -2.36524
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.42533
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00063 - Z-value: -3.41634
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00007 - Z-value: 4.08887
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 7.20208
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.39069
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00224 - Z-value: 2.89983
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00006 - Z-value: -4.01281
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00018 - Z-value: -3.67240
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.06342 - Z-value: 1.72478
# >> Comparison: 1 - 3      - Adjusted P-value: 0.05482 - Z-value: 2.09077
# >> Comparison: 2 - 4      - Adjusted P-value: 0.07421 - Z-value: -1.78529
# >> Comparison: 3 - 4      - Adjusted P-value: 0.09689 - Z-value: -2.14071
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.07020 - Z-value: 1.67670
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00328 - Z-value: 3.06391
# >> Comparison: 2 - 3      - Adjusted P-value: 0.06142 - Z-value: 1.63352
# >> Comparison: 2 - 4      - Adjusted P-value: 0.03264 - Z-value: -2.13653
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00179 - Z-value: -3.43279
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00303 - Z-value: -3.28740
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00490 - Z-value: -2.94172
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00553 - Z-value: -2.77442
# >> Comparison: 2 - 4      - Adjusted P-value: 0.01268 - Z-value: -2.38861
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00002 - Z-value: 4.56213
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01066 - Z-value: 2.55378
# >> Comparison: 2 - 3      - Adjusted P-value: 0.09976 - Z-value: -1.38429
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00124 - Z-value: -3.34362
# >> Comparison: 3 - 4      - Adjusted P-value: 0.08746 - Z-value: -1.56916
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.06082 - Z-value: -2.04812
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08807 - Z-value: -2.17866
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 5.07856
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00061 - Z-value: -3.42838
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00030 - Z-value: -3.71682
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.02914 - Z-value: -2.58586
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02642 - Z-value: 2.37361
# >> Comparison: 2 - 4      - Adjusted P-value: 0.05171 - Z-value: 1.94557
# =============================================


# Plot
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
  scale_color_manual(values = c("1" = "#C2BAB4", "2" = "#FEC089", "3" = "#F06A00", "4" = "#8C2D04")) +
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

ggsave("99_Fig/fig3/tSNE_cluster_tissue_proportion.png", plot = plot_tissue, width = 6, height = 2)
ggsave("99_Fig/fig3/tSNE_cluster_tissue_proportion.pdf", plot = plot_tissue, width = 6, height = 2)


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
  scale_color_manual(values = c("1" = "#C2BAB4", "2" = "#FEC089", "3" = "#F06A00", "4" = "#8C2D04")) +
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

ggsave("99_Fig/fig3/tSNE_cluster_KOOS_pre.png", plot = plot_KOOS_Pre , width = 5, height = 1.25, bg = "white")
ggsave("99_Fig/fig3/tSNE_cluster_KOOS_pre.pdf", plot = plot_KOOS_Pre , width = 5, height = 1.25, bg = "white")

plot_KOOS_Post12 <- ggplot(df_long_KOOOS_Post12, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0, 110)) +
  facet_wrap(~ KOOS, scales = "fixed", labeller = labeller(KOOS = custom_labels_KOOS_post12), strip.position = "bottom", ncol = 5) + 
  labs(x = "Post-Operation 12-Month", y = "Score") +
  scale_color_manual(values = c("1" = "#C2BAB4", "2" = "#FEC089", "3" = "#F06A00", "4" = "#8C2D04")) +
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

ggsave("99_Fig/fig3/tSNE_cluster_KOOS_post12.png", plot = plot_KOOS_Post12, width = 5, height = 1.25)
ggsave("99_Fig/fig3/tSNE_cluster_KOOS_post12.pdf", plot = plot_KOOS_Post12, width = 5, height = 1.25)


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
#  [1] gtable_0.3.6      compiler_4.3.3    tidyselect_1.2.1  Rcpp_1.0.14       systemfonts_1.2.1 scales_1.3.0      textshaping_1.0.0 R6_2.6.1          labeling_0.4.3    generics_0.1.3    munsell_0.5.1     pillar_1.10.1     tzdb_0.5.0        rlang_1.1.5       stringi_1.8.7    
# [16] timechange_0.3.0  cli_3.6.4         withr_3.0.2       magrittr_2.0.3    grid_4.3.3        rstudioapi_0.17.1 hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5       glue_1.8.0        farver_2.1.2      ragg_1.3.3        colorspace_2.1-1  tools_4.3.3       pkgconfig_2.0.3  

