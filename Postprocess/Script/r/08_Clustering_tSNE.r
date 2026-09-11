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


# # silhouette
# library(cluster)
# 
# set.seed(123)
# silhouette_results <- lapply(2:8, function(k) {
#   km <- kmeans(df[,6:17], centers = k)
#   sil <- silhouette(km$cluster, dist(df[,6:17]))
#   data.frame(
#     k = k,
#     silhouette_score = mean(sil[,"sil_width"])
#   )
# })
# 
# silhouette_summary <- bind_rows(silhouette_results)
# silhouette_summary
#   k silhouette_score
# 1 2        0.3714283
# 2 3        0.4653234
# 3 4        0.3541900
# 4 5        0.3719614
# 5 6        0.3632589
# 6 7        0.3491665
# 7 8        0.3369623
#
# write.table(silhouette_summary, "03_res/Clustering/tissue_composition_cluster_k.txt", sep = "\t", row.names = FALSE)
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
# ggsave("99_Fig/fig2/silhouette_cluster_number.png", plot = plot_silhouette_k, width = 2, height = 1.5)
# ggsave("99_Fig/fig2/silhouette_cluster_number.pdf", plot = plot_silhouette_k, width = 2, height = 1.5)


# library(fpc)
# set.seed(123)
# 
# cb <- clusterboot(
#   df[, 6:17],
#   B = 1000,
#   bootmethod = "boot",
#   clustermethod = kmeansCBI,
#   krange = 3,
#   seed = 123
# )
# 
# cb$bootmean
# [1] 0.9527773 0.9635779 0.9570855


# tSNE
set.seed(123)
tsne_res <- Rtsne(df[,6:17], dims= 2, perplexity = 36)
tsne_data <- data.frame(tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2])
rownames(tsne_data) <- rownames(df)


# clustering
kmeans_result <- kmeans(df[,6:17], centers = 3)

tsne_data$cluster <- as.factor(kmeans_result$cluster)
tsne_data$Diagnosis <- df$Diagnosis

tsne_data$Diagnosis_plot <- recode(
  tsne_data$Diagnosis,
  "SLE" = "Other autoimmune diseases",
  "SSc" = "Other autoimmune diseases"
)

plot_tsne <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = cluster, shape = Diagnosis_plot)) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("1" = "#FEC089", "2" = "#F06A00", "3" = "#8C2D04")) +
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

ggsave("99_Fig/fig2/tSNE.png", plot = plot_tsne, width = 2.25, height = 1.25)
ggsave("99_Fig/fig2/tSNE.pdf", plot = plot_tsne, width = 2.25, height = 1.25)

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
                    labels = c("nonOA" = "non-OA structural", "RA" = "RA", "OA" = "OA", "SLE" = "Other autoimmune diseases", "SSc" = "Other autoimmune diseases")) +
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

ggsave("99_Fig/fig2/cluster_diagnosis_proportion.png", plot = plot_stacked_bar, width = 2, height = 1.25)
ggsave("99_Fig/fig2/cluster_diagnosis_proportion.pdf", plot = plot_stacked_bar, width = 2, height = 1.25)


# Statistics
df_OA <- df[df$Diagnosis =="OA",]

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
# >> Comparison: 1 - 2      - Adjusted P-value: 0.02524 - Z-value: -2.39041
# >> Comparison: 2 - 3      - Adjusted P-value: 0.07634 - Z-value: 1.63624
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_Pain **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.03093 - Z-value: 2.31486
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_Symptom **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01729 - Z-value: 2.52639
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02379 - Z-value: 2.14786
# =============================================
# ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_QOL **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.04566 - Z-value: 1.87434
# >> Comparison: 2 - 3      - Adjusted P-value: 0.08535 - Z-value: 1.90407
# =============================================
# ** Significant Dunn Post-hoc Test Results for Synovitis_score **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.07060 - Z-value: 1.67400
# >> Comparison: 1 - 3      - Adjusted P-value: 0.07133 - Z-value: 1.98133
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00032 - Z-value: 3.69897
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00491 - Z-value: 2.71934
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Stroma **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.04342 - Z-value: -1.89648
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01064 - Z-value: -2.69228
# =============================================
# ** Significant Dunn Post-hoc Test Results for Krenn_Infiltrate **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.07046 - Z-value: 1.67493
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01546 - Z-value: 2.56534
# =============================================
# ** Significant Dunn Post-hoc Test Results for Total **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00208 - Z-value: 3.19750
# >> Comparison: 1 - 3      - Adjusted P-value: 0.05991 - Z-value: 1.75136
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
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: -5.12291
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00161 - Z-value: 2.94656
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.70172
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00197 - Z-value: -3.00803
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -4.68453
# >> Comparison: 2 - 3      - Adjusted P-value: 0.02583 - Z-value: -1.94591
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00333 - Z-value: -2.84468
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -4.92569
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01088 - Z-value: -2.29456
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00266 - Z-value: 2.91597
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: 5.31085
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00477 - Z-value: 2.59200
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
# >> Comparison: 1 - 3      - Adjusted P-value: 0.06132 - Z-value: 2.04476
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.01476 - Z-value: 2.33227
# >> Comparison: 1 - 3      - Adjusted P-value: 0.01060 - Z-value: 2.69365
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.05197 - Z-value: -1.81654
# >> Comparison: 2 - 3      - Adjusted P-value: 0.01382 - Z-value: -2.60391
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00220 - Z-value: 3.18057
# >> Comparison: 1 - 3      - Adjusted P-value: 0.03407 - Z-value: 2.00070
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: 1 - 3      - Adjusted P-value: 0.00803 - Z-value: -2.78480
# >> Comparison: 2 - 3      - Adjusted P-value: 0.09016 - Z-value: -1.55390
# =============================================
# ** Significant Dunn Post-hoc Test Results for RBC **
# =============================================
# >> Comparison: 2 - 3      - Adjusted P-value: 0.09576 - Z-value: 1.85330
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: 1 - 2      - Adjusted P-value: 0.00010 - Z-value: -3.98320
# >> Comparison: 2 - 3      - Adjusted P-value: 0.00082 - Z-value: 3.26533
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
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels), strip.position = "bottom", ncol = 4) + 
  labs(y = "CLR(Proportion)") +
  scale_color_manual(values = c("1" = "#FEC089", "2" = "#F06A00", "3" = "#8C2D04")) +
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

ggsave("99_Fig/fig2/cluster_tissue_proportion_CLR.png", plot = plot_tissue, width = 3.5, height = 2.5)
ggsave("99_Fig/fig2/cluster_tissue_proportion_CLR.pdf", plot = plot_tissue, width = 3.5, height = 2.5)


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
  scale_color_manual(values = c("1" = "#FEC089", "2" = "#F06A00", "3" = "#8C2D04")) +
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

ggsave("99_Fig/sup_fig5/cluster_KOOS_pre.png", plot = plot_KOOS_Pre , width = 3.5, height = 1.25, bg = "white")
ggsave("99_Fig/sup_fig5/cluster_KOOS_pre.pdf", plot = plot_KOOS_Pre , width = 3.5, height = 1.25, bg = "white")

plot_KOOS_Post12 <- ggplot(df_long_KOOS_Post12, aes(x = cluster, y = Value, color = cluster)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0, 110)) +
  facet_wrap(~ KOOS, scales = "fixed", labeller = labeller(KOOS = custom_labels_KOOS_post12), strip.position = "bottom", ncol = 5) + 
  labs(x = "Post-Operation 12-Month", y = "Score") +
  scale_color_manual(values = c("1" = "#FEC089", "2" = "#F06A00", "3" = "#8C2D04")) +
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

ggsave("99_Fig/sup_fig5/cluster_KOOS_post12.png", plot = plot_KOOS_Post12, width = 3.5, height = 1.25)
ggsave("99_Fig/sup_fig5/cluster_KOOS_post12.pdf", plot = plot_KOOS_Post12, width = 3.5, height = 1.25)



# For response to reviewers
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
# 
# 
# # # silhouette
# # library(cluster)
# # 
# # set.seed(123)
# # silhouette_results <- lapply(2:8, function(k) {
# #   km <- kmeans(df[,6:17], centers = k)
# #   sil <- silhouette(km$cluster, dist(df[,6:17]))
# #   data.frame(
# #     k = k,
# #     silhouette_score = mean(sil[,"sil_width"])
# #   )
# # })
# # 
# # silhouette_summary <- bind_rows(silhouette_results)
# # silhouette_summary
# #   k silhouette_score
# # 1 2        0.3870862
# # 2 3        0.4871587
# # 3 4        0.4899119
# # 4 5        0.3150136
# # 5 6        0.3058462
# # 6 7        0.2923956
# # 7 8        0.2872692
# #
# # write.table(silhouette_summary, "03_res/Clustering/Revise/tissue_composition_cluster_k.txt", sep = "\t", row.names = FALSE)
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
# tsne_res <- Rtsne(df[,6:17], dims= 2, perplexity = 26)
# tsne_data <- data.frame(tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2])
# rownames(tsne_data) <- rownames(df)
# 
# 
# # clustering
# kmeans_result <- kmeans(df[,6:17], centers = 3)
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
#   scale_color_manual(values = c("1" = "#FEC089", "2" = "#F06A00", "3" = "#8C2D04")) +
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
# ggsave("99_Fig/fig2/Revise/tSNE.png", plot = plot_tsne, width = 1.75, height = 1.25)
# ggsave("99_Fig/fig2/Revise/tSNE.pdf", plot = plot_tsne, width = 1.75, height = 1.25)
# 
# df[,42:44] <- tsne_data[,1:3]
# 
# # write.table(df, "01_formatted/Revise/annotations_full_tissue_proportion_cluster.txt", sep = "\t", row.names = TRUE, col.names = NA)
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
#   if (length(results[[score]]$significant_comparisons) > 0) {
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
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.01296 - Z-value: -2.62601
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.06574 - Z-value: 1.70792
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_Pain **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.04080 - Z-value: 2.20866
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.08479 - Z-value: 1.58462
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_Symptom **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.01919 - Z-value: 2.48949
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.02030 - Z-value: 2.21059
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Post_12Month__KOOS_QOL **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.06084 - Z-value: 1.74422
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.05612 - Z-value: 2.08126
# # =============================================
# 
# 
# all_component <- colnames(df_OA)[6:17]
# major_component <- c(
#   "Adipose",
#   "Fibrous_tissue__dense_irregular",
#   "Fibrous_tissue__dense_regular",
#   "Fibrous_tissue__loose"
# )
# minor_component <- setdiff(all_component, major_component)
# 
# res_clr_major <- clr(acomp(df_OA[,all_component])) 
# 
# df_clr_major <- as.data.frame(res_clr_major)
# df_clr_major$cluster <- df_OA$cluster
# 
# 
# # statistics
# ## shapiro_test
# shapiro_test_results <- df_clr_major %>% 
#   summarise(across(1:12, ~ shapiro.test(.)$p.value))
# shapiro_test_results 
# #       Adipose Fibrous_tissue__dense_irregular Fibrous_tissue__dense_regular Fibrous_tissue__loose       TLS    Plasma       Stroma       Lining      Muscle       RBC Micro_vessel Large_vessel
# # 1 0.001620628                      0.01034722                     0.2144998             0.8376733 0.3980882 0.1245807 5.468267e-05 4.213218e-06 9.57414e-06 0.1283558    0.6595747   0.01884009
# 
# ## kruskal.test
# results <- list()
# for (comp in all_component) {
#   kw_result <- kruskal.test(as.formula(paste(comp, "~ cluster")), data = df_clr_major)
#   dunn_result <- dunn.test(df_clr_major[[comp]], df_clr_major$cluster, method = "bh")
#   significant_idx <- which(dunn_result$P.adjusted < 0.1)
#   significant_comparisons <- dunn_result$comparisons[significant_idx]
#   significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
#   significant_z <- dunn_result$Z[significant_idx]
#   results[[comp]] <- list(
#     test_type = "Kruskal-Wallis",
#     main_test = kw_result,
#     posthoc_test = dunn_result,
#     significant_comparisons = significant_comparisons,
#     significant_p_adjusted = significant_p_adjusted,
#     significant_z = significant_z
#   )
# }
# 
# for (comp in all_component) {
#   cat("\n--- Results for", comp, "---\n")
#   print(results[[comp]]$main_test)
#   if (length(results[[comp]]$significant_comparisons) > 0) {
#     cat("\n** Significant Dunn Post-hoc Test Results for", comp, "**\n")
#     cat("=============================================\n")
#     for (i in seq_along(results[[comp]]$significant_comparisons)) {
#       cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
#                   results[[comp]]$significant_comparisons[i], 
#                   results[[comp]]$significant_p_adjusted[i],
#                   results[[comp]]$significant_z[i]))
#     }
#     cat("=============================================\n")
#   }
# }
# 
# # Notes: Significant differences are shown.
# # ** Significant Dunn Post-hoc Test Results for Adipose **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00000 - Z-value: -5.12024
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00196 - Z-value: 2.88413
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00000 - Z-value: 6.67065
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00436 - Z-value: -2.75803
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00001 - Z-value: -4.61935
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.01662 - Z-value: -2.12927
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00488 - Z-value: -2.72089
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: -4.89888
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00782 - Z-value: -2.41739
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00149 - Z-value: 3.09307
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00000 - Z-value: 5.38095
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00503 - Z-value: 2.57374
# # =============================================
# 
# 
# row_sum <- rowSums(df_OA[, minor_component], na.rm = TRUE)
# 
# df_recal <- df_OA
# df_recal[, minor_component] <- df_OA[, minor_component] / row_sum
# 
# res_clr_minor <- clr(acomp(df_recal[,minor_component])) 
# 
# df_clr_minor <- as.data.frame(res_clr_minor)
# df_clr_minor$cluster <- df_OA$cluster
# 
# 
# # statistics
# ## shapiro_test
# shapiro_test_results <- df_clr_minor %>% 
#   summarise(across(1:8, ~ shapiro.test(.)$p.value))
# shapiro_test_results 
# #         TLS     Plasma      Stroma       Lining       Muscle        RBC Micro_vessel Large_vessel
# # 1 0.2667346 0.05463301 0.001368096 2.547756e-05 8.975871e-05 0.04094924    0.1877178    0.4792868
# 
# ## kruskal.test
# results <- list()
# for (comp in minor_component) {
#   kw_result <- kruskal.test(as.formula(paste(comp, "~ cluster")), data = df_clr_minor)
#   dunn_result <- dunn.test(df_clr_minor[[comp]], df_clr_minor$cluster, method = "bh")
#   significant_idx <- which(dunn_result$P.adjusted < 0.1)
#   significant_comparisons <- dunn_result$comparisons[significant_idx]
#   significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
#   significant_z <- dunn_result$Z[significant_idx]
#   results[[comp]] <- list(
#     test_type = "Kruskal-Wallis",
#     main_test = kw_result,
#     posthoc_test = dunn_result,
#     significant_comparisons = significant_comparisons,
#     significant_p_adjusted = significant_p_adjusted,
#     significant_z = significant_z
#   )
# }
# 
# for (comp in minor_component) {
#   cat("\n--- Results for", comp, "---\n")
#   print(results[[comp]]$main_test)
#   if (length(results[[comp]]$significant_comparisons) > 0) {
#     cat("\n** Significant Dunn Post-hoc Test Results for", comp, "**\n")
#     cat("=============================================\n")
#     for (i in seq_along(results[[comp]]$significant_comparisons)) {
#       cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
#                   results[[comp]]$significant_comparisons[i], 
#                   results[[comp]]$significant_p_adjusted[i],
#                   results[[comp]]$significant_z[i]))
#     }
#     cat("=============================================\n")
#   }
# }
# 
# # Notes: Significant differences are shown.
# # ** Significant Dunn Post-hoc Test Results for TLS **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.04764 - Z-value: 2.14740
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Plasma **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.01418 - Z-value: 2.34730
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00997 - Z-value: 2.71389
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Stroma **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.04511 - Z-value: -1.87976
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.01951 - Z-value: -2.48359
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Lining **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00116 - Z-value: 3.36145
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.02775 - Z-value: 2.08581
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Muscle **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.07514 - Z-value: -1.64398
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.00548 - Z-value: -2.90688
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.07906 - Z-value: -1.41142
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for RBC **
# # =============================================
# # >> Comparison: 1 - 3      - Adjusted P-value: 0.07983 - Z-value: 1.61440
# # =============================================
# # ** Significant Dunn Post-hoc Test Results for Large_vessel **
# # =============================================
# # >> Comparison: 1 - 2      - Adjusted P-value: 0.00029 - Z-value: -3.72772
# # >> Comparison: 2 - 3      - Adjusted P-value: 0.00143 - Z-value: 3.10430
# # =============================================
# 
# 
# # Plot
# df_OA[,major_component] <- df_OA[,major_component]
# df_OA[,minor_component] <- df_OA[,minor_component]
# 
# df_long <- df_OA %>%
#   pivot_longer(cols = 6:17, names_to = "Tissue", values_to = "Value")
# 
# df_long$Tissue <- factor(df_long$Tissue, levels = c("Fibrous_tissue__dense_irregular", "Fibrous_tissue__dense_regular", "Fibrous_tissue__loose", "Lining", "TLS", "Plasma", "Stroma", "RBC", "Micro_vessel", "Large_vessel", "Adipose", "Muscle"))
# df_long$cluster <- factor(df_long$cluster)
# 
# custom_labels <- c(
#   "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense.irregular)",
#   "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense.regular)",
#   "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
#   "Micro_vessel" = "Micro vessel",
#   "Large_vessel" = "Large vessel"
# )
# 
# df_long$Tissue <- factor(
#   df_long$Tissue,
#   levels = c(
#     "Lining",
#     "Micro_vessel",
#     "RBC",
#     "Fibrous_tissue__loose",
#     "TLS",
#     "Plasma",
#     "Large_vessel",
#     "Adipose",
#     "Fibrous_tissue__dense_irregular",
#     "Fibrous_tissue__dense_regular",
#     "Stroma",
#     "Muscle"
#   )
# )
# 
# plot_tissue <- ggplot(df_long, aes(x = cluster, y = Value, color = cluster)) +
#   geom_boxplot(width = 0.8, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
#   geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
#   facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels), strip.position = "bottom", ncol = 4) + 
#   labs(y = "CLR(Proportion)") +
#   scale_color_manual(values = c("1" = "#FEC089", "2" = "#F06A00", "3" = "#8C2D04")) +
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
# ggsave("99_Fig/fig2/Revise/cluster_tissue_proportion_CLR.png", plot = plot_tissue, width = 3.5, height = 2.5)
# ggsave("99_Fig/fig2/Revise/cluster_tissue_proportion_CLR.pdf", plot = plot_tissue, width = 3.5, height = 2.5)

