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
                    labels = c("nonOA" = "Trauma", "RA" = "RA", "OA" = "OA", "SLE" = "Other autoimmune diseases", "SSc" = "Other autoimmune diseases")) +
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
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: -5.01555
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00799 - Z-value: 2.47512
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.84673
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00028 - Z-value: 3.62911
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00056 - Z-value: -3.37317
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00208 - Z-value: -3.07831
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00002 - Z-value: -4.48304
# >> Comparison: 2 - 3      - Adjusted P-value: 0.04847 - Z-value: -1.84783
# >> Comparison: 2 - 4      - Adjusted P-value: 0.04950 - Z-value: 1.73634
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00157 - Z-value: 3.27809
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00459 - Z-value: -2.83465
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -4.97152
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00804 - Z-value: -2.55162
# >> Comparison: 2 - 4      - Adjusted P-value: 0.02766 - Z-value: 1.99449
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00004 - Z-value: 4.17512
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00026 - Z-value: 3.75291
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: 5.96977
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00432 - Z-value: 2.76108
# >> Comparison: 1 - 4      - Adjusted P-value: 0.00633 - Z-value: 2.55733
# >> Comparison: 3 - 4      - Adjusted P-value: 0.00048 - Z-value: -3.49443
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS **
# =============================================
# >> Comparison: 3 - 4      - Adjusted P-value: 0.08239 - Z-value: -1.91945
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01619 - Z-value: 2.78235
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01878 - Z-value: 2.49709
# >> Comparison: 2 - 4      - Adjusted P-value: 0.04858 - Z-value: -1.97228
# >> Comparison: 3 - 4      - Adjusted P-value: 0.05183 - Z-value: -1.81770
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02262 - Z-value: -2.67205
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00003 - Z-value: 4.43927
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00144 - Z-value: 3.30147
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00466 - Z-value: -2.82977
# >> Comparison: 3 - 4      - Adjusted P-value: 0.03691 - Z-value: -1.96676
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.06883 - Z-value: -2.27438
# >> Comparison: 2 - 3      - Adjusted P-value: 0.04794 - Z-value: -1.97791
# >> Comparison: 3 - 4      - Adjusted P-value: 0.04931 - Z-value: 2.13362
# =============================================
# ** Significant Dunn Post-hoc Test Results for RBC **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.05088 - Z-value: 2.38759
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00288 - Z-value: 3.30220
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00795 - Z-value: 2.78815
# >> Comparison: 1 - 4      - Adjusted P-value: 0.05159 - Z-value: 1.94651
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00353 - Z-value: -3.04198
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00617 - Z-value: 3.08175
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00357 - Z-value: 2.91398
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
#         TLS     Plasma      Stroma       Lining       Muscle        RBC Micro_vessel Large_vessel
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
# >> Comparison: 1 - 3      - Adjusted P-value: 0.06323 - Z-value: 2.03198
# >> Comparison: 3 - 4      - Adjusted P-value: 0.09542 - Z-value: -2.14681
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.02551 - Z-value: 2.23352
# >> Comparison: 1 - 3      - Adjusted P-value: 0.02004 - Z-value: 2.71238
# >> Comparison: 2 - 4      - Adjusted P-value: 0.06181 - Z-value: -1.73685
# >> Comparison: 3 - 4      - Adjusted P-value: 0.03442 - Z-value: -2.27431
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.02431 - Z-value: -2.40425
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01156 - Z-value: -2.88980
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00246 - Z-value: 3.34601
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01452 - Z-value: 2.58708
# >> Comparison: 2 - 4      - Adjusted P-value: 0.02735 - Z-value: -2.20654
# >> Comparison: 3 - 4      - Adjusted P-value: 0.07551 - Z-value: -1.64157
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01104 - Z-value: -2.90431
# >> Comparison: 2 - 3      - Adjusted P-value: 0.07253 - Z-value: -1.79578
# >> Comparison: 3 - 4      - Adjusted P-value: 0.01670 - Z-value: 2.53848
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00127 - Z-value: -3.52421
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00149 - Z-value: 3.29274
# >> Comparison: 2 - 4      - Adjusted P-value: 0.00117 - Z-value: 3.24619
# =============================================


# Plot
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
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels), strip.position = "bottom", ncol = 4) + 
  labs(y = "CLR(Proportion)") +
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

ggsave("99_Fig/fig3/tSNE_cluster_tissue_proportion_CLR.png", plot = plot_tissue, width = 4, height = 3)
ggsave("99_Fig/fig3/tSNE_cluster_tissue_proportion_CLR.pdf", plot = plot_tissue, width = 4, height = 3)


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

