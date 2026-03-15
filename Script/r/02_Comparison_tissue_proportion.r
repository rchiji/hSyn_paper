library(dplyr)
library(tidyverse)
library(dunn.test)
library(ggplot2)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)

df_OA_nonOA_RA <- df[!df$Diagnosis %in% c("SLE", "SSc"), ]
df_OA_nonOA_RA <- df_OA_nonOA_RA[df_OA_nonOA_RA$Joint == "Knee",]

component <- colnames(df_OA_nonOA_RA[,6:17])


# statistics
## shapiro_test
shapiro_test_results <- df_OA_nonOA_RA %>% 
  summarise(across(6:17, ~ shapiro.test(.)$p.value))
shapiro_test_results 
#        Adipose Fibrous_tissue__dense_irregular Fibrous_tissue__dense_regular Fibrous_tissue__loose          TLS       Plasma       Stroma       Lining       Muscle          RBC Micro_vessel Large_vessel
# 1 8.335292e-07                    3.509635e-08                  1.856782e-17            0.01145201 2.528112e-21 2.191833e-19 3.852913e-22 7.103635e-13 5.921528e-23 1.052633e-16 0.0006813595 5.238785e-14

## kruskal.test
results <- list()
for (comp in component) {
  kw_result <- kruskal.test(as.formula(paste(comp, "~ Diagnosis")), data = df_OA_nonOA_RA)
  dunn_result <- dunn.test(df_OA_nonOA_RA[[comp]], df_OA_nonOA_RA$Diagnosis, method = "bh")
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

# Notes: Significant differences are shown.
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.08805 - Z-value: 1.35285
# >> Comparison: nonOA - RA - Adjusted P-value: 0.08687 - Z-value: -1.57255
# >> Comparison: OA - RA    - Adjusted P-value: 0.00384 - Z-value: -3.01580
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00028 - Z-value: -3.74079
# >> Comparison: OA - RA    - Adjusted P-value: 0.00024 - Z-value: -3.59462
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: OA - RA    - Adjusted P-value: 0.05296 - Z-value: 2.10485
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00008 - Z-value: -4.05409
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00005 - Z-value: -3.98247
# >> Comparison: OA - RA    - Adjusted P-value: 0.07704 - Z-value: -1.42526
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00054 - Z-value: -3.56943
# >> Comparison: nonOA - RA - Adjusted P-value: 0.01394 - Z-value: -2.35359
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00228 - Z-value: -2.96310
# >> Comparison: OA - RA    - Adjusted P-value: 0.00088 - Z-value: -3.43900
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.02876 - Z-value: -2.07111
# >> Comparison: OA - RA    - Adjusted P-value: 0.01327 - Z-value: 2.61787
# =============================================
# ** Significant Dunn Post-hoc Test Results for RBC **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00022 - Z-value: -3.79351
# >> Comparison: OA - RA    - Adjusted P-value: 0.04369 - Z-value: 1.89376
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.07775 - Z-value: 1.94453
# =============================================


# Plot
custom_labels <- c(
  "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense, irregular)",
  "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense, regular)",
  "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
  "Micro_vessel" = "Micro vessel",
  "Large_vessel" = "Large vessel"
)

df_log <- df_OA_nonOA_RA[,6:17]
df_log <- log10(df_log + 1e-6)

df_log$Diagnosis <- df_OA_nonOA_RA[,1]
df_log$Diagnosis <- factor(df_log$Diagnosis, levels = c("nonOA", "RA", "OA"))

df_log_long <- df_log %>%
  pivot_longer(cols = all_of(component),
               names_to = "Tissue",
               values_to = "Value")
colnames(df_log_long) <- c("Diagnosis", "Tissue", "Value")

df_log_long$Tissue <- factor(df_log_long$Tissue, levels = c("Fibrous_tissue__dense_irregular", "Fibrous_tissue__dense_regular", "Fibrous_tissue__loose", "Lining", "TLS", "Plasma", "Stroma", "RBC", "Micro_vessel", "Large_vessel", "Adipose",   "Muscle"))

p <- ggplot(df_log_long, aes(x = Diagnosis, y = Value, color = Diagnosis)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels), strip.position = "bottom", ncol = 6) +
  labs(y = "log10(Proportion)") +
  scale_color_manual(values = c("nonOA" = "#1f77b4", "RA" = "#2ca02c", "OA" = "#ff7f0e"),
                     labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA")) +
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

ggsave("99_Fig/fig2/nonOA_OA_RA_tissue_proportion.png", plot = p , width = 6, height = 2)
ggsave("99_Fig/fig2/nonOA_OA_RA_tissue_proportion.pdf", plot = p , width = 6, height = 2)

