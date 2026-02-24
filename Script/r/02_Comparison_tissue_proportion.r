library(dplyr)
library(tidyverse)
library(dunn.test)
library(ggplot2)


df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)

df_OA_nonOA_RA <- df[df$Diagnosis != "SLE",]
df_OA_nonOA_RA <- df_OA_nonOA_RA[df_OA_nonOA_RA$Joint == "Knee",]

component <- colnames(df_OA_nonOA_RA[,6:17])


# statistics
## shapiro_test
shapiro_test_results <- df_OA_nonOA_RA %>% 
  summarise(across(6:17, ~ shapiro.test(.)$p.value))
shapiro_test_results 
#        Adipose Fibrous_tissue__dense_irregular Fibrous_tissue__dense_regular Fibrous_tissue__loose          TLS       Plasma       Stroma      Lining       Muscle          RBC Micro_vessel Large_vessel
# 1 1.003404e-06                    2.705474e-08                   1.60232e-17            0.01359999 1.987955e-21 1.693646e-19 3.202193e-22 5.46652e-13 4.694428e-23 8.222257e-17 0.0005710695 7.777662e-14

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

# --- Results for Adipose ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Adipose by Diagnosis
# Kruskal-Wallis chi-squared = 2.1599, df = 2, p-value = 0.3396
# 
# 
# --- Results for Fibrous_tissue__dense_irregular ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Fibrous_tissue__dense_irregular by Diagnosis
# Kruskal-Wallis chi-squared = 10.115, df = 2, p-value = 0.006361
# 
# 
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_irregular **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.08763 - Z-value: 1.35552
# >> Comparison: nonOA - RA - Adjusted P-value: 0.08356 - Z-value: -1.59185
# >> Comparison: OA - RA    - Adjusted P-value: 0.00352 - Z-value: -3.04230
# =============================================
# 
# --- Results for Fibrous_tissue__dense_regular ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Fibrous_tissue__dense_regular by Diagnosis
# Kruskal-Wallis chi-squared = 15.528, df = 2, p-value = 0.0004247
# 
# 
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__dense_regular **
# =============================================
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00027 - Z-value: -3.74903
# >> Comparison: OA - RA    - Adjusted P-value: 0.00027 - Z-value: -3.56907
# =============================================
# 
# --- Results for Fibrous_tissue__loose ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Fibrous_tissue__loose by Diagnosis
# Kruskal-Wallis chi-squared = 4.5066, df = 2, p-value = 0.1051
# 
# 
# ** Significant Dunn Post-hoc Test Results for Fibrous_tissue__loose **
# =============================================
# >> Comparison: OA - RA    - Adjusted P-value: 0.05100 - Z-value: 2.12010
# =============================================
# 
# --- Results for TLS ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  TLS by Diagnosis
# Kruskal-Wallis chi-squared = 20.361, df = 2, p-value = 3.789e-05
# 
# 
# ** Significant Dunn Post-hoc Test Results for TLS **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00009 - Z-value: -4.01335
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00005 - Z-value: -3.98194
# >> Comparison: OA - RA    - Adjusted P-value: 0.07171 - Z-value: -1.46315
# =============================================
# 
# --- Results for Plasma ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Plasma by Diagnosis
# Kruskal-Wallis chi-squared = 12.963, df = 2, p-value = 0.001532
# 
# 
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00052 - Z-value: -3.57777
# >> Comparison: nonOA - RA - Adjusted P-value: 0.01335 - Z-value: -2.36968
# =============================================
# 
# --- Results for Stroma ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Stroma by Diagnosis
# Kruskal-Wallis chi-squared = 11.725, df = 2, p-value = 0.002844
# 
# 
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00252 - Z-value: -2.93253
# >> Comparison: OA - RA    - Adjusted P-value: 0.00121 - Z-value: -3.35012
# =============================================
# 
# --- Results for Lining ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Lining by Diagnosis
# Kruskal-Wallis chi-squared = 9.696, df = 2, p-value = 0.007844
# 
# 
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.02945 - Z-value: -2.06135
# >> Comparison: OA - RA    - Adjusted P-value: 0.01373 - Z-value: 2.60639
# =============================================
# 
# --- Results for Muscle ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Muscle by Diagnosis
# Kruskal-Wallis chi-squared = 3.236, df = 2, p-value = 0.1983
# 
# 
# --- Results for RBC ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  RBC by Diagnosis
# Kruskal-Wallis chi-squared = 15.989, df = 2, p-value = 0.0003373
# 
# 
# ** Significant Dunn Post-hoc Test Results for RBC **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00024 - Z-value: -3.77169
# >> Comparison: OA - RA    - Adjusted P-value: 0.04746 - Z-value: 1.85722
# =============================================
# 
# --- Results for Micro_vessel ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Micro_vessel by Diagnosis
# Kruskal-Wallis chi-squared = 4.0286, df = 2, p-value = 0.1334
# 
# 
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.07353 - Z-value: 1.96843
# =============================================
# 
# --- Results for Large_vessel ---
# 
# 	Kruskal-Wallis rank sum test
# 
# data:  Large_vessel by Diagnosis
# Kruskal-Wallis chi-squared = 1.0987, df = 2, p-value = 0.5773


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
#  [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0 dunn.test_1.3.6 dplyr_1.1.4    
# 
# loaded via a namespace (and not attached):
#  [1] gtable_0.3.6      compiler_4.3.3    tidyselect_1.2.1  systemfonts_1.2.1 scales_1.3.0      textshaping_1.0.0 R6_2.6.1          labeling_0.4.3    generics_0.1.3    munsell_0.5.1     pillar_1.10.1     tzdb_0.5.0        rlang_1.1.5       stringi_1.8.7     timechange_0.3.0 
# [16] cli_3.6.4         withr_3.0.2       magrittr_2.0.3    grid_4.3.3        rstudioapi_0.17.1 hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5       glue_1.8.0        farver_2.1.2      ragg_1.3.3        colorspace_2.1-1  tools_4.3.3       pkgconfig_2.0.3  