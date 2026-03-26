library(Rtsne)
library(ggplot2)
library(dunn.test)
library(dplyr)
library(tidyverse)
library(compositions)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
score <- read.delim("01_formatted/annotations_full_OARSI_Krenn.txt", sep = "\t", row.names = 1)
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
  scale_shape_manual(values = c("nonOA" = 16, "RA" = 15, "OA" = 17, "SLE" = 18, "SSc" = 18),
                     labels = c("nonOA" = "Trauma", "RA" = "RA", "OA" = "OA", "SLE" = "Other autoimmune diseases", "SSc" = "Other autoimmune diseases")) +
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

ggsave("99_Fig/fig4/tSNE.png", plot = plot_tsne, width = 2.5, height = 1.5)
ggsave("99_Fig/fig4/tSNE.pdf", plot = plot_tsne, width = 2.5, height = 1.5)

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

ggsave("99_Fig/fig4/tSNE_cluster_entropy_proportion.png", plot = plot_stacked_bar, width = 1.25, height = 1.25)
ggsave("99_Fig/fig4/tSNE_cluster_entropy_proportion.pdf", plot = plot_stacked_bar, width = 1.25, height = 1.25)


# Statistics
df_OA <- df[df$Diagnosis == "OA",]

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
# ** Significant Dunn Post-hoc Test Results for Pre__KOOS_Pain **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.08345 - Z-value: 2.19988
# =============================================
# ** Significant Dunn Post-hoc Test Results for Pre__KOOS_Symptom **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.03245 - Z-value: -2.29676
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08056 - Z-value: -1.74746
# >> Comparison: 3 - 4      - Adjusted P-value: 0.01084 - Z-value: 2.91003
# =============================================
# ** Significant Dunn Post-hoc Test Results for Pre__KOOS_QOL **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.05752 - Z-value: 2.34214
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_3Month__KOOS_Sports **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08091 - Z-value: -2.21197
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_3Month__KOOS_QOL **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02777 - Z-value: -2.60248
# =============================================
# 
# ** Significant Dunn Post-hoc Test Results for Synovitis_score **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.03687 - Z-value: -2.24795
# >> Comparison: 2 - 3      - Adjusted P-value: 0.03884 - Z-value: 2.48525
# >> Comparison: 2 - 4      - Adjusted P-value: 0.03972 - Z-value: 2.05663
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.07107 - Z-value: -1.67078
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00480 - Z-value: 2.81991
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00002 - Z-value: 4.48033
# >> Comparison: 1 - 4      - Adjusted P-value: 0.08421 - Z-value: 1.37731
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00272 - Z-value: 3.11958
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07456 - Z-value: -1.53711
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Infiltrate **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.04112 - Z-value: -2.20557
# >> Comparison: 2 - 3      - Adjusted P-value: 0.04396 - Z-value: 2.44087
# >> Comparison: 1 - 4      - Adjusted P-value: 0.07949 - Z-value: -1.61653
# >> Comparison: 3 - 4      - Adjusted P-value: 0.05675 - Z-value: -1.90523
# =============================================
# ** Significant Dunn Post-hoc Test Results for Total **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.03122 - Z-value: 2.31131
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00140 - Z-value: 3.49938
# >> Comparison: 2 - 4      - Adjusted P-value: 0.08396 - Z-value: 1.58948
# >> Comparison: 3 - 4      - Adjusted P-value: 0.05477 - Z-value: -1.92073
# =============================================


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
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00235 - Z-value: -2.88543
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.92943
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00023 - Z-value: -3.68902
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00000 - Z-value: -5.49821
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00178 - Z-value: -3.03945
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: 5.57528
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: -5.35828
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00000 - Z-value: 4.57418
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00001 - Z-value: 4.50077
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00001 - Z-value: 4.73494
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00031 - Z-value: -3.71045
# >> Comparison: 1 - 4      - Adjusted P-value: 0.01799 - Z-value: 2.36592
# >> Comparison: 2 - 4      - Adjusted P-value: 0.02465 - Z-value: -2.13378
# >> Comparison: 3 - 4      - Adjusted P-value: 0.05871 - Z-value: 1.65535
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00001 - Z-value: -4.61852
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00449 - Z-value: 2.84178
# >> Comparison: 1 - 4      - Adjusted P-value: 0.08703 - Z-value: -1.57162
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00618 - Z-value: 2.86885
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01609 - Z-value: 2.78439
# >> Comparison: 3 - 4      - Adjusted P-value: 0.03835 - Z-value: -2.23270
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00052 - Z-value: -3.47172
# >> Comparison: 1 - 3      - Adjusted P-value: 0.08892 - Z-value: 1.44592
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00001 - Z-value: 4.63982
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00260 - Z-value: -2.92266
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00007 - Z-value: -4.08879
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: 1 - 4      - Adjusted P-value: 0.03003 - Z-value: 2.57550
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00459 - Z-value: -2.83437
# >> Comparison: 1 - 3      - Adjusted P-value: 0.05974 - Z-value: 1.75272
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00003 - Z-value: 4.39350
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00024 - Z-value: 3.77536
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.08809 - Z-value: -1.70554
# >> Comparison: 1 - 4      - Adjusted P-value: 0.06721 - Z-value: 1.69747
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00142 - Z-value: 3.49635
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00127 - Z-value: -3.33784
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00026 - Z-value: -3.92512
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00442 - Z-value: 2.84634
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
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08524 - Z-value: 1.90466
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00015 - Z-value: -3.79888
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00005 - Z-value: 4.29722
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00033 - Z-value: -3.51200
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00008 - Z-value: -4.03045
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.05255 - Z-value: 2.10799
# >> Comparison: 1 - 4      - Adjusted P-value: 0.01512 - Z-value: 2.80448
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.06103 - Z-value: -1.87328
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00601 - Z-value: 2.87738
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00704 - Z-value: 3.04253
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.03899 - Z-value: 2.48390
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.01794 - Z-value: 2.74897
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00043 - Z-value: -3.62452
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00001 - Z-value: -4.60713
# >> Comparison: 2 - 4      - Adjusted P-value: 0.04480 - Z-value: -1.88278
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00679 - Z-value: 2.70717
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

plot_tissue <- ggplot(df_long, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
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

ggsave("99_Fig/sup_fig7/tSNE_cluster_tissue_proportion_CLR.png", plot = plot_tissue, width = 6, height = 2)
ggsave("99_Fig/sup_fig7/tSNE_cluster_tissue_proportion_CLR.pdf", plot = plot_tissue, width = 6, height = 2)


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


tmp <- read.delim("01_formatted/annotations_full_tissue_count_thickness.txt", sep = "\t", row.names = 1)
tmp <- tmp[tmp$Joint == "Knee",]
df <- read.delim("01_formatted/annotations_full_tissue_conection_cluster.txt", sep = "\t", row.names = 1)
df[,45:48] <- tmp[,40:43]
df <- df[df$Diagnosis == "OA",]

component <- colnames(df)[45:48]

custom_labels_ratio <- c(
  "Micro_vessel_ratio" = "Micro vessel",
  "Large_vessel_ratio" = "Large vessel",
  "TLS_ratio" = "TLS"
)

cluster_colors <- c("1" = "#D8C6C2", "2" = "#F3B2A6", "3" = "#C73A3A", "4" = "#7F1212")

df$cluster <- factor(df$cluster, levels = c("1", "2", "3", "4"))

results <- list()
for (comp in component) {
  kw_result <- kruskal.test(as.formula(paste(comp, "~ cluster")), data = df)
  dunn_result <- dunn.test(df[[comp]], df$cluster, method = "bh")
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

for (comp in component) {
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

# ** Significant Dunn Post-hoc Test Results for Lining_thickness **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.04831 - Z-value: 2.40657
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro_vessel_ratio **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00002 - Z-value: 4.24022
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 5.41393
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00059 - Z-value: 3.35782
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00001 - Z-value: 4.61240
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel_ratio **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.03242 - Z-value: -2.13925
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02083 - Z-value: -2.69948
# >> Comparison: 3 - 4      - Adjusted P-value: 0.04800 - Z-value: 2.14442
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS_ratio **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.04784 - Z-value: 2.14575
# >> Comparison: 3 - 4      - Adjusted P-value: 0.03202 - Z-value: -2.55319
# =============================================


# plot
df$cluster <- factor(df$cluster, levels = c(1, 2, 3, 4))

df_long <- df %>%
  pivot_longer(cols = all_of(component),
               names_to = "Tissue",
               values_to = "Value")

df_long_lining <- df_long %>%
  filter(Tissue %in% "Lining_thickness")

df_long_ratio_tls <- df_long %>%
  filter(Tissue %in% "TLS_ratio")

df_long_ratio_vessel <- df_long %>%
  filter(Tissue %in% c("Micro_vessel_ratio", "Large_vessel_ratio"))

df_long_ratio_vessel$Tissue <- factor(df_long_ratio_vessel$Tissue,
                                      levels = c("Micro_vessel_ratio", "Large_vessel_ratio"))
p1 <- ggplot(df_long_lining, aes(x = cluster, y = Value, color = cluster)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  labs(x = "Lining", y = "Mean thickness") +
  scale_color_manual(values = cluster_colors) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 5),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )
ggsave("99_Fig/fig4/OA_cluster_Lining_thickness.png", plot = p1, width = 1.8, height = 1.25)
ggsave("99_Fig/fig4/OA_cluster_Lining_thickness.pdf", plot = p1, width = 1.8, height = 1.25)

p2 <- ggplot(df_long_ratio_tls, aes(x = cluster, y = Value, color = cluster)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  labs(x = "TLS", y = "Number/Area") +
  scale_color_manual(values = cluster_colors) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 5),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )
ggsave("99_Fig/fig4/OA_cluster_tls_ratio.png", plot = p2, width = 1.8, height = 1.25)
ggsave("99_Fig/fig4/OA_cluster_tls_ratio.pdf", plot = p2, width = 1.8, height = 1.25)

p3 <- ggplot(df_long_ratio_vessel, aes(x = cluster, y = Value, color = cluster)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels_ratio), strip.position = "bottom", ncol = 3) +
  labs(x = "", y = "Number/Area") +
  scale_color_manual(values = cluster_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 5, margin = margin(t = 0)),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )
ggsave("99_Fig/fig4/OA_cluster_vessel_ratio.png", plot = p3, width = 3, height = 1.25)
ggsave("99_Fig/fig4/OA_cluster_vessel_ratio.pdf", plot = p3, width = 3, height = 1.25)

