library(Rtsne)
library(ggplot2)
library(dunn.test)
library(dplyr)
library(tidyverse)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
score <- read.delim("01_formatted/annotations_full_OARSI_Krenn.txt", sep = "\t", row.names = 1)
df[,37:41]  <- score[,1:5]
df <- df[df$Joint == "Knee",]


# tSNE
set.seed(123)
tsne_res <- Rtsne(df[,6:17], dims= 2, perplexity = 36)
tsne_data <- data.frame(tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2])
rownames(tsne_data) <- rownames(df)


# clustering
kmeans_result <- kmeans(tsne_data, centers = 4)

tsne_data$cluster <- as.factor(kmeans_result$cluster)
tsne_data$Diagnosis <- df$Diagnosis

plot_tsne <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = cluster, shape = Diagnosis)) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("1" = "#C2BAB4", "2" = "#FEC089", "3" = "#F06A00", "4" = "#8C2D04")) +
  scale_shape_manual(values = c("nonOA" = 16, "RA" = 15, "OA" = 17, "SLE" = 18, "SSc" = 18),
                     labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA", "SLE" = "Other autoimmune diseases", "SSc" = "Other autoimmune diseases")) +
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

ggsave("99_Fig/fig3/tSNE.png", plot = plot_tsne, width = 2.5, height = 1.5)
ggsave("99_Fig/fig3/tSNE.pdf", plot = plot_tsne, width = 2.5, height = 1.5)

df[,42:44] <- tsne_data[,1:3]

# write.table(df, "01_formatted/annotations_full_tissue_proportion_cluster.txt", sep = "\t", row.names = TRUE, col.names = NA)


# Proportion
df_stacked_bar <- df %>% 
  count(cluster, Diagnosis)

df_stacked_bar$Diagnosis <- factor(df_stacked_bar$Diagnosis, levels = c("OA", "RA", "SLE", "SSc", "nonOA"))

plot_stacked_bar <- ggplot(df_stacked_bar, aes(x = cluster, y = n, fill = Diagnosis)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.9) +
  labs(y = "Fraction of diagnosis type", fill = "Diagnosis") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = c("nonOA" = "#1f77b4", "RA" = "#2ca02c", "OA" = "#ff7f0e", "SLE" = "grey80", "SSc" = "grey80"),
                    labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA", "SLE" = "Other autoimmune diseases", "SSc" = "Other autoimmune diseases")) +
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

ggsave("99_Fig/fig3/tSNE_cluster_diagnosis_proportion.png", plot = plot_stacked_bar, width = 2.5, height = 1.6)
ggsave("99_Fig/fig3/tSNE_cluster_diagnosis_proportion.pdf", plot = plot_stacked_bar, width = 2.5, height = 1.6)


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
# ** Significant Dunn Post-hoc Test Results for Post_3Month__KOOS_QOL **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.06297 - Z-value: -1.85937
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08610 - Z-value: 1.90024
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_Pain **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.04642 - Z-value: -2.42114
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_Symptom **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.04969 - Z-value: 2.39626
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02340 - Z-value: 2.26679
# >> Comparison: 3 - 4      - Adjusted P-value: 0.02935 - Z-value: -2.33461
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_QOL **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07866 - Z-value: -1.75852
# =============================================
# 
# ** Significant Dunn Post-hoc Test Results for Synovitis_score **
# =============================================
# >> Comparison: 2 - 4      - Adjusted P-value: 0.06816 - Z-value: -2.00059
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07012 - Z-value: -2.26730
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00051 - Z-value: 3.75753
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00716 - Z-value: 2.82188
# >> Comparison: 2 - 4      - Adjusted P-value: 0.01841 - Z-value: -2.35731
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07298 - Z-value: -1.65807
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Stroma **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01223 - Z-value: -2.40194
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01509 - Z-value: -2.43032
# >> Comparison: 1 - 4      - Adjusted P-value: 0.01989 - Z-value: -2.47663
# >> Comparison: 2 - 4      - Adjusted P-value: 0.03618 - Z-value: -2.51040
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Infiltrate **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.09372 - Z-value: 1.67608
# >> Comparison: 1 - 3      - Adjusted P-value: 0.02757 - Z-value: 2.60490
# >> Comparison: 3 - 4      - Adjusted P-value: 0.04991 - Z-value: -2.12880
# =============================================
# ** Significant Dunn Post-hoc Test Results for Total **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01085 - Z-value: 2.68600
# >> Comparison: 1 - 3      - Adjusted P-value: 0.08296 - Z-value: 1.59543
# >> Comparison: 2 - 4      - Adjusted P-value: 0.01285 - Z-value: -2.85646
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07176 - Z-value: -1.80062
# =============================================
# 
# ** Significant Dunn Post-hoc Test Results for Adipose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: -6.24292
# >> Comparison: 1 - 3      - Adjusted P-value: 0.06542 - Z-value: 1.51078
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.93965
# >> Comparison: 1 - 4      - Adjusted P-value: 0.01489 - Z-value: -2.24414
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00025 - Z-value: 3.66660
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00045 - Z-value: -3.42920
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01124 - Z-value: -2.43257
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.94099
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00011 - Z-value: -3.87771
# >> Comparison: 1 - 4      - Adjusted P-value: 0.03536 - Z-value: -1.88870
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00008 - Z-value: 4.05524
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01312 - Z-value: -2.37619
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.39298
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00074 - Z-value: -3.37424
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00004 - Z-value: 4.19855
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 7.14408
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.34313
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00284 - Z-value: 2.82415
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00008 - Z-value: -3.93734
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00022 - Z-value: -3.61890
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.06132 - Z-value: 1.74060
# >> Comparison: 1 - 3      - Adjusted P-value: 0.05117 - Z-value: 2.11871
# >> Comparison: 2 - 4      - Adjusted P-value: 0.03870 - Z-value: -2.06732
# >> Comparison: 3 - 4      - Adjusted P-value: 0.04926 - Z-value: -2.39943
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.07285 - Z-value: 1.65892
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00367 - Z-value: 3.02981
# >> Comparison: 2 - 3      - Adjusted P-value: 0.06384 - Z-value: 1.61458
# >> Comparison: 2 - 4      - Adjusted P-value: 0.02544 - Z-value: -2.23462
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00145 - Z-value: -3.48980
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00278 - Z-value: -3.31154
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00460 - Z-value: -2.96124
# >> Comparison: 1 - 4      - Adjusted P-value: 0.01045 - Z-value: -2.56049
# >> Comparison: 2 - 4      - Adjusted P-value: 0.02237 - Z-value: -2.17242
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00002 - Z-value: 4.53945
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01160 - Z-value: 2.52400
# >> Comparison: 2 - 3      - Adjusted P-value: 0.09788 - Z-value: -1.39463
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00096 - Z-value: -3.41337
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07508 - Z-value: -1.64433
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.05619 - Z-value: -2.08071
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08284 - Z-value: -2.20272
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 5.04974
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00069 - Z-value: -3.39480
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00020 - Z-value: -3.81841
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.02751 - Z-value: -2.60569
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02449 - Z-value: 2.40157
# >> Comparison: 2 - 4      - Adjusted P-value: 0.02648 - Z-value: 2.21914
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

