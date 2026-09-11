library(dplyr)
library(dunn.test)
library(tidyverse)
library(compositions)
library(edgeR)
library(DESeq2)
library(WGCNA)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
anno <- read.delim(file = "00_src/Entropy_group.txt", sep = "\t", row.names = 1)


# Tissue proportion
df$Entropy <-  anno$Entropy
df_OA <- df[df$Diagnosis == "OA",]

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
df_clr_major$Entropy <- df_OA$Entropy


## statistics
### shapiro_test
shapiro_test_results <- df_clr_major %>% 
  summarise(across(1:12, ~ shapiro.test(.)$p.value))
shapiro_test_results 
#       Adipose Fibrous_tissue__dense_irregular Fibrous_tissue__dense_regular Fibrous_tissue__loose       TLS    Plasma       Stroma       Lining      Muscle       RBC Micro_vessel Large_vessel
# 1 0.001620628                      0.01034722                     0.2144998             0.8376733 0.3980882 0.1245807 5.468267e-05 4.213218e-06 9.57414e-06 0.1283558    0.6595747   0.01884009

### kruskal.test
results <- list()
for (comp in all_component) {
  kw_result <- kruskal.test(as.formula(paste(comp, "~ Entropy")), data = df_clr_major)
  dunn_result <- dunn.test(df_clr_major[[comp]], df_clr_major$Entropy, method = "bh")
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
# >> Comparison: High - Low - Adjusted P-value: 0.01121 - Z-value: -2.67500
# >> Comparison: High - Mid - Adjusted P-value: 0.02786 - Z-value: -2.08419
# >> Comparison: Low - Mid  - Adjusted P-value: 0.07757 - Z-value: 1.42164
# =============================================
  

row_sum <- rowSums(df_OA[, minor_component], na.rm = TRUE)

df_recal <- df_OA
df_recal[, minor_component] <- df_OA[, minor_component] / row_sum

res_clr_minor <- clr(acomp(df_recal[,minor_component])) 

df_clr_minor <- as.data.frame(res_clr_minor)
df_clr_minor$Entropy <- df_OA$Entropy


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
  kw_result <- kruskal.test(as.formula(paste(comp, "~ Entropy")), data = df_clr_minor)
  dunn_result <- dunn.test(df_clr_minor[[comp]], df_clr_minor$Entropy, method = "bh")
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
# >> Comparison: High - Low - Adjusted P-value: 0.01149 - Z-value: 2.66672
# >> Comparison: High - Mid - Adjusted P-value: 0.04726 - Z-value: 1.67197
# >> Comparison: Low - Mid  - Adjusted P-value: 0.06850 - Z-value: -1.68839
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.00025 - Z-value: 3.75991
# >> Comparison: High - Mid - Adjusted P-value: 0.01130 - Z-value: 2.28007
# >> Comparison: Low - Mid  - Adjusted P-value: 0.01126 - Z-value: -2.43218
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.00351 - Z-value: -3.04296
# >> Comparison: High - Mid - Adjusted P-value: 0.08378 - Z-value: -1.38009
# >> Comparison: Low - Mid  - Adjusted P-value: 0.01699 - Z-value: 2.27928
# =============================================
# ** Significant Dunn Post-hoc Test Results for RBC **
# =============================================
# >> Comparison: High - Mid - Adjusted P-value: 0.07656 - Z-value: 1.95113
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.03721 - Z-value: -2.24435
# >> Comparison: Low - Mid  - Adjusted P-value: 0.08581 - Z-value: 1.36699
# =============================================


## Plot
df_OA[,major_component] <- df_OA[,major_component]
df_OA[,minor_component] <- df_OA[,minor_component]

df_long <- df_OA %>%
  pivot_longer(cols = 6:17, names_to = "Tissue", values_to = "Value")

df_long$Tissue <- factor(df_long$Tissue, levels = c("Fibrous_tissue__dense_irregular", "Fibrous_tissue__dense_regular", "Fibrous_tissue__loose", "Lining", "TLS", "Plasma", "Stroma", "RBC", "Micro_vessel", "Large_vessel", "Adipose", "Muscle"))
df_long$Entropy <- factor(df_long$Entropy, levels = c("Low", "Mid", "High"))

custom_labels <- c(
     "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense.irregular)",
     "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense.regular)",
     "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
     "Micro_vessel" = "Micro vessel",
     "Large_vessel" = "Large vessel"
 )

p <- ggplot(df_long, aes(x = Entropy, y = Value, color = Entropy)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.1) +
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels), strip.position = "bottom", ncol = 6) +
  labs(y = "log10(Proportion)") +
  scale_color_manual(values = c("Low" = "#bfe6bf", "Mid" = "#66cc66", "High" = "#336633")) +
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

ggsave("99_Fig/sup_fig6/Entropy_tissue_proportion_CLR.png", plot = p , width = 6, height = 2)
ggsave("99_Fig/sup_fig6/Entropy_tissue_proportion_CLR.pdf", plot = p , width = 6, height = 2)

