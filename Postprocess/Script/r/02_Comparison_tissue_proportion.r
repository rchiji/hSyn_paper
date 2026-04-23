library(compositions)
library(dplyr)
library(dunn.test)
library(tidyverse)
library(ggplot2)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)


df_OA_nonOA_RA <- df[!df$Diagnosis %in% c("SLE", "SSc"), ]
df_OA_nonOA_RA <- df_OA_nonOA_RA[df_OA_nonOA_RA$Joint %in% "Knee",]

all_component <- colnames(df_OA_nonOA_RA)[6:17]
major_component <- c(
  "Adipose",
  "Fibrous_tissue__dense_irregular",
  "Fibrous_tissue__dense_regular",
  "Fibrous_tissue__loose"
)
minor_component <- setdiff(all_component, major_component)

res_clr_major <- clr(acomp(df_OA_nonOA_RA[,all_component])) 

df_clr_major <- as.data.frame(res_clr_major)
df_clr_major$Diagnosis <- df_OA_nonOA_RA$Diagnosis


# statistics
## shapiro_test
shapiro_test_results <- df_clr_major %>% 
  summarise(across(1:12, ~ shapiro.test(.)$p.value))
shapiro_test_results 
#        Adipose Fibrous_tissue__dense_irregular Fibrous_tissue__dense_regular Fibrous_tissue__loose         TLS      Plasma      Stroma      Lining       Muscle        RBC Micro_vessel Large_vessel
# 1 0.0008058735                    0.0003998907                     0.4726107             0.5232412 0.006739443 0.002289207 3.58699e-06 4.23945e-05 4.271288e-07 0.01672667  0.003029626   0.04232758

## kruskal.test
results <- list()
for (comp in all_component) {
  kw_result <- kruskal.test(as.formula(paste(comp, "~ Diagnosis")), data = df_clr_major)
  dunn_result <- dunn.test(df_clr_major[[comp]], df_clr_major$Diagnosis, method = "bh")
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
# >> Comparison: OA - RA    - Adjusted P-value: 0.08414 - Z-value: 1.58846
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: nonOA - RA - Adjusted P-value: 0.06265 - Z-value: -1.73051
# >> Comparison: OA - RA    - Adjusted P-value: 0.04351 - Z-value: -2.18338
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00055 - Z-value: -3.56121
# >> Comparison: OA - RA    - Adjusted P-value: 0.00056 - Z-value: -3.37247
# =============================================
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.05018 - Z-value: -1.83228
# >> Comparison: OA - RA    - Adjusted P-value: 0.04376 - Z-value: 2.18112
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: nonOA - RA - Adjusted P-value: 0.05845 - Z-value: 2.06454
# >> Comparison: OA - RA    - Adjusted P-value: 0.09972 - Z-value: 1.28314
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00593 - Z-value: 2.88177
# >> Comparison: OA - RA    - Adjusted P-value: 0.05416 - Z-value: -1.79779
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00790 - Z-value: -2.55793
# >> Comparison: OA - RA    - Adjusted P-value: 0.00108 - Z-value: 3.38214
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.01143 - Z-value: 2.66853
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00661 - Z-value: 2.61898
# =============================================
# ** Significant Dunn Post-hoc Test Results for RBC **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00001 - Z-value: -4.50896
# >> Comparison: OA - RA    - Adjusted P-value: 0.01699 - Z-value: 2.27928
# =============================================


row_sum <- rowSums(df_OA_nonOA_RA[, minor_component], na.rm = TRUE)

df_recal <- df_OA_nonOA_RA
df_recal[, minor_component] <- df_OA_nonOA_RA[, minor_component] / row_sum

res_clr_minor <- clr(acomp(df_recal[,minor_component])) 

df_clr_minor <- as.data.frame(res_clr_minor)
df_clr_minor$Diagnosis <- df_OA_nonOA_RA$Diagnosis


# statistics
## shapiro_test
shapiro_test_results <- df_clr_minor %>% 
  summarise(across(1:8, ~ shapiro.test(.)$p.value))
shapiro_test_results 
# TLS      Plasma       Stroma       Lining       Muscle         RBC Micro_vessel Large_vessel
# 1 0.01425803 0.001757575 5.974902e-05 9.252408e-05 1.735906e-05 0.003593411   0.01299333    0.5435163

## kruskal.test
results <- list()
for (comp in minor_component) {
  kw_result <- kruskal.test(as.formula(paste(comp, "~ Diagnosis")), data = df_clr_minor)
  dunn_result <- dunn.test(df_clr_minor[[comp]], df_clr_minor$Diagnosis, method = "bh")
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
# >> Comparison: nonOA - OA - Adjusted P-value: 0.01400 - Z-value: -2.35204
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00522 - Z-value: -2.92175
# >> Comparison: OA - RA    - Adjusted P-value: 0.05909 - Z-value: -1.56242
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00871 - Z-value: 2.52365
# >> Comparison: OA - RA    - Adjusted P-value: 0.00268 - Z-value: -3.12400
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.01744 - Z-value: -2.26934
# >> Comparison: nonOA - RA - Adjusted P-value: 0.08694 - Z-value: 1.35987
# >> Comparison: OA - RA    - Adjusted P-value: 0.00064 - Z-value: 3.52103
# =============================================
# ** Significant Dunn Post-hoc Test Results for Muscle **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.02509 - Z-value: 2.39259
# >> Comparison: nonOA - RA - Adjusted P-value: 0.03824 - Z-value: 1.95164
# =============================================
# ** Significant Dunn Post-hoc Test Results for RBC **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00001 - Z-value: -4.58236
# >> Comparison: nonOA - RA - Adjusted P-value: 0.05361 - Z-value: -1.80243
# >> Comparison: OA - RA    - Adjusted P-value: 0.05084 - Z-value: 1.63675
# =============================================


# Plot
df_OA_nonOA_RA[,major_component] <- df_clr_major[,major_component]
df_OA_nonOA_RA[,minor_component] <- df_clr_minor[,minor_component]

custom_labels <- c(
  "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense, irregular)",
  "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense, regular)",
  "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
  "Micro_vessel" = "Micro vessel",
  "Large_vessel" = "Large vessel"
)

df_OA_nonOA_RA$Diagnosis <- factor(df_OA_nonOA_RA$Diagnosis, levels = c("nonOA", "RA", "OA"))

df_long <- df_OA_nonOA_RA %>%
  pivot_longer(cols = all_of(all_component),
               names_to = "Tissue",
               values_to = "Value")

df_long$Tissue <- factor(df_long$Tissue, levels = c("Fibrous_tissue__dense_irregular", "Fibrous_tissue__dense_regular", "Fibrous_tissue__loose", "Lining", "TLS", "Plasma", "Stroma", "RBC", "Micro_vessel", "Large_vessel", "Adipose",   "Muscle"))

p <- ggplot(df_long, aes(x = Diagnosis, y = Value, color = Diagnosis)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels), strip.position = "bottom", ncol = 6) +
  labs(y = "CLR(Proportion)") +
  scale_color_manual(values = c("nonOA" = "#1f77b4", "RA" = "#2ca02c", "OA" = "#ff7f0e"),
                     labels = c("nonOA" = "non-OA structural", "RA" = "RA", "OA" = "OA")) +
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

ggsave("99_Fig/sup_fig3/nonOA_OA_RA_tissue_proportion_CLR.png", plot = p , width = 6, height = 2)
ggsave("99_Fig/sup_fig3/nonOA_OA_RA_tissue_proportion_CLR.pdf", plot = p , width = 6, height = 2)

