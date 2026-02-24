library(dplyr)
library(tidyverse)
library(dunn.test)
library(ggplot2)

df <- read.delim(file = "00_src/annotations_full.txt", sep = "\t", row.names = 1)


# Aggregation
## Micro vessel
micro_vessel_paths <- list.files(path = "00_src/Tissue_number_thickness/count_vessel/", full.names = TRUE)
micro_vessel_list <- lapply(micro_vessel_paths, read.delim)

samplenames <- sapply(micro_vessel_paths, function(x) {
  x <- basename(x)
  gsub(pattern = "_HE.*",replacement = "",x)
  }
)
names(micro_vessel_list) <- samplenames

micro_vessel_list <- mapply(function(df, name) {
    df$SampleName <- name
    df
}, micro_vessel_list, samplenames, SIMPLIFY = FALSE)

df_micro_vessel <- do.call(rbind, micro_vessel_list)

## Large vessel
large_vessel_paths <- list.files(path = "00_src/Tissue_number_thickness/count_vessel_large_filtered/", full.names = TRUE)
large_vessel_list <- lapply(large_vessel_paths, read.delim)

samplenames <- sapply(large_vessel_paths, function(x) {
  x <- basename(x)
  gsub(pattern = "_HE.*",replacement = "",x)
  })
names(large_vessel_list) <- samplenames

large_vessel_list <- mapply(function(df, name) {
    df$SampleName <- name
    df
}, large_vessel_list, samplenames, SIMPLIFY = FALSE)

df_large_vessel <- do.call(rbind, large_vessel_list)

## TLS
tls_paths <- list.files(path = "00_src/Tissue_number_thickness/count_tls_filtered/", full.names = TRUE)
tls_list <- lapply(tls_paths, read.delim)

samplenames <- sapply(tls_paths, function(x) {
  x <- basename(x)
  gsub(pattern = "_HE.*",replacement = "",x)
})
names(tls_list) <- samplenames

tls_list <- mapply(function(df, name) {
  df$SampleName <- name
  df
}, tls_list, samplenames, SIMPLIFY = FALSE)

df_tls <- do.call(rbind, tls_list)

## Lining
lining_paths <- list.files(path = "00_src/Tissue_number_thickness/thickness_lining/", full.names = TRUE)
lining_list <- lapply(lining_paths, read.delim)

samplenames <- sapply(large_vessel_paths, function(x) {
  x <- basename(x)
  gsub(pattern = "_HE.*",replacement = "",x)
})
names(lining_list) <- samplenames

mean_thickness <- sapply(lining_list, function(df) {
  mean(df$Diameter.mean, na.rm = TRUE)
})

df_lining <- data.frame(
  SampleName = names(mean_thickness),
  Thickness  = mean_thickness,
  row.names  = names(mean_thickness)
)

## merge
df_micro_vessel <- df_micro_vessel %>%
  rename(Micro_vessel_number = Count)
df_large_vessel <- df_large_vessel %>%
  rename(Large_vessel_number = Count)
df_tls <- df_tls %>% 
  rename(TLS_number = Count)
df_lining <- df_lining %>% 
  rename(Lining_thickness = Thickness)

df_tissue <- list(
  df_micro_vessel,
  df_large_vessel,
  df_tls,
  df_lining
) %>% 
  Reduce(function(x, y) full_join(x, y, by = "SampleName"), .)

df_tissue <- column_to_rownames(df_tissue, var = "SampleName")


# Tissue Area
file_paths <- list.files(path = "00_src/Annotation_ratio_v2/", full.names = TRUE)
file_list <- lapply(file_paths, read.delim)

samplenames <- sapply(file_paths, function(x) {
  x <- basename(x)
  gsub(pattern = "_HE.*",replacement = "",x)
  })
names(file_list) <- samplenames

levels <- c("adipose", "Fibro(dense,irregular)", "Fibro(dense,regular)", 
"Fibro(loose)", "Immune cells", "lining", "muscle", "plasma", 
"RBC", "Stroma", "vessel", "vessel(large)")

summary_list <- lapply(file_list, function(x){
  x$Classification <- factor(x$Classification, levels = levels)
  x <- x %>% 
    group_by(Classification) %>% 
    summarise(
      Area_sum = sum(Area.µm.2),
      count = n()
    )
  return(x)
  })


# Density
total_area <- sapply(summary_list, function(df) sum(df$Area_sum, na.rm = TRUE))
df_total_area <- as.data.frame(total_area)
df_tissue$Area <- df_total_area$total_area
df_tissue <- df_tissue %>%
  mutate(
    Micro_vessel_ratio = Micro_vessel_number / Area,
    Large_vessel_ratio = Large_vessel_number / Area,
    TLS_ratio = TLS_number / Area
  )
df[,37:44] <- df_tissue[,c(1:4,6:8,5)]
# write.table(df, "01_formatted/annotations_full_tissue_count_thickness.txt", sep = "\t", row.names = TRUE, col.names = NA)

df_OA_nonOA_RA <- df[df$Diagnosis != "SLE",]
df_OA_nonOA_RA <- df_OA_nonOA_RA[df_OA_nonOA_RA$Joint == "Knee",]

component <- colnames(df_OA_nonOA_RA[,37:44])


# statistics
## shapiro_test
shapiro_test_results <- df_OA_nonOA_RA %>% summarise(across(37:43, ~ shapiro.test(.)$p.value))
shapiro_test_results 
#   Micro_vessel_number Large_vessel_number  TLS_number Lining_thickness Micro_vessel_ratio Large_vessel_ratio    TLS_ratio
# 1          3.6835e-12        8.742707e-13 1.19011e-20      0.006392078        0.009115167        2.14505e-12 2.181534e-20

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

cat("\n========== Kruskal-Wallis Results ==========\n")
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

# --- Results for Micro_vessel_number ---
# data:  Micro_vessel_number by Diagnosis
# Kruskal-Wallis chi-squared = 29.292, df = 2, p-value = 4.358e-07
# ** Significant Dunn Post-hoc Test Results for Micro_vessel_number **
#   =============================================
#   >> Comparison: nonOA - OA - Adjusted P-value: 0.00000 - Z-value: -5.34741
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00014 - Z-value: -3.73326
# =============================================
#   
# --- Results for Large_vessel_number ---
# data:  Large_vessel_number by Diagnosis
# Kruskal-Wallis chi-squared = 7.6526, df = 2, p-value = 0.02179
# ** Significant Dunn Post-hoc Test Results for Large_vessel_number **
# =============================================
#   >> Comparison: nonOA - OA - Adjusted P-value: 0.00871 - Z-value: -2.75859
# >> Comparison: nonOA - RA - Adjusted P-value: 0.06115 - Z-value: -1.74187
# =============================================
#   
# --- Results for TLS_number ---
# data:  TLS_number by Diagnosis
# Kruskal-Wallis chi-squared = 17.772, df = 2, p-value = 0.0001383
# ** Significant Dunn Post-hoc Test Results for TLS_number **
# =============================================
#   >> Comparison: nonOA - OA - Adjusted P-value: 0.00033 - Z-value: -3.51518
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00014 - Z-value: -3.91540
# >> Comparison: OA - RA    - Adjusted P-value: 0.03620 - Z-value: -1.79662
# =============================================
#   
# --- Results for Lining_thickness ---
# data:  Lining_thickness by Diagnosis
# Kruskal-Wallis chi-squared = 21.438, df = 2, p-value = 2.212e-05
# ** Significant Dunn Post-hoc Test Results for Lining_thickness **
# =============================================
#   >> Comparison: nonOA - OA - Adjusted P-value: 0.00045 - Z-value: -3.43160
# >> Comparison: OA - RA    - Adjusted P-value: 0.00053 - Z-value: 3.57025
# =============================================
#   
# --- Results for Micro_vessel_ratio ---
# data:  Micro_vessel_ratio by Diagnosis
# Kruskal-Wallis chi-squared = 0.54213, df = 2, p-value = 0.7626
# 
# --- Results for Large_vessel_ratio ---
# data:  Large_vessel_ratio by Diagnosis
# Kruskal-Wallis chi-squared = 0.2435, df = 2, p-value = 0.8854
# 
# --- Results for TLS_ratio ---
# data:  TLS_ratio by Diagnosis
# Kruskal-Wallis chi-squared = 17.49, df = 2, p-value = 0.0001593
# ** Significant Dunn Post-hoc Test Results for TLS_ratio **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00032 - Z-value: -3.52286
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00017 - Z-value: -3.85943
# >> Comparison: OA - RA    - Adjusted P-value: 0.04246 - Z-value: -1.72284
# =============================================

# Plot
custom_labels_num <- c(
  "Micro_vessel_number" = "Micro vessel",
  "Large_vessel_number" = "Large vessel",
  "TLS_number" = "TLS"
)
custom_labels_ratio <- c(
  "Micro_vessel_ratio" = "Micro vessel",
  "Large_vessel_ratio" = "Large vessel",
  "TLS_ratio" = "TLS"
)

df_OA_nonOA_RA$Diagnosis <- factor(df_OA_nonOA_RA$Diagnosis,　levels = c("nonOA", "RA", "OA"))

df_long <- df_OA_nonOA_RA %>%
    pivot_longer(cols = all_of(component),
                 names_to = "Tissue",
                 values_to = "Value")

df_long_lining <- df_long %>%
  filter(Tissue %in% "Lining_thickness")
df_long_number_tls <- df_long %>%
  filter(Tissue %in% "TLS_number")
df_long_ratio_tls <- df_long %>%
  filter(Tissue %in% "TLS_ratio")
df_long_number_vessel <- df_long %>%
  filter(Tissue %in% c("Micro_vessel_number", "Large_vessel_number"))
df_long_ratio_vessel <- df_long %>%
  filter(Tissue %in% c("Micro_vessel_ratio", "Large_vessel_ratio"))

df_long_number_vessel$Tissue <- factor(df_long_number_vessel$Tissue, levels = c("Micro_vessel_number", "Large_vessel_number"))
df_long_ratio_vessel$Tissue <- factor(df_long_ratio_vessel$Tissue, levels = c("Micro_vessel_ratio", "Large_vessel_ratio"))

p1 <- ggplot(df_long_lining, aes(x = Diagnosis, y = Value, color = Diagnosis)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  labs(x = "Lining", y = "Mean thickness") +
  scale_color_manual(values = c("nonOA" = "#1f77b4", "RA" = "#2ca02c", "OA" = "#ff7f0e"),
                     labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA")) +
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
    legend.margin     = margin(0, 0, 0, 0))
ggsave("99_Fig/fig2/nonOA_OA_RA_Lining_thickness.png", plot = p1 , width = 1.8, height = 1.25)
ggsave("99_Fig/fig2/nonOA_OA_RA_Lining_thickness.pdf", plot = p1 , width = 1.8, height = 1.25)

p2 <- ggplot(df_long_number_tls, aes(x = Diagnosis, y = Value, color = Diagnosis)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  labs(x = "TLS", y = "Number") +
  scale_color_manual(values = c("nonOA" = "#1f77b4", "RA" = "#2ca02c", "OA" = "#ff7f0e"),
                     labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA")) +
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
    legend.margin     = margin(0, 0, 0, 0))
ggsave("99_Fig/fig2/nonOA_OA_RA_tls_num.png", plot = p2 , width = 1.8, height = 1.25)
ggsave("99_Fig/fig2/nonOA_OA_RA_tls_num.pdf", plot = p2 , width = 1.8, height = 1.25)

p3 <- ggplot(df_long_ratio_tls, aes(x = Diagnosis, y = Value, color = Diagnosis)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  labs(x = "TLS", y = "Number/Area") +
  scale_color_manual(values = c("nonOA" = "#1f77b4", "RA" = "#2ca02c", "OA" = "#ff7f0e"),
                     labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA")) +
  scale_y_continuous(limits = c(0, 3e-06)) +
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
    legend.margin     = margin(0, 0, 0, 0))
ggsave("99_Fig/fig2/nonOA_OA_RA_tls_ratio.png", plot = p3 , width = 1.8, height = 1.25)
ggsave("99_Fig/fig2/nonOA_OA_RA_tls_ratio.pdf", plot = p3 , width = 1.8, height = 1.25)

p4 <- ggplot(df_long_number_vessel, aes(x = Diagnosis, y = Value, color = Diagnosis)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels_num), strip.position = "bottom", ncol = 3) +
  labs(x = "", y = "Number") +
  scale_color_manual(values = c("nonOA" = "#1f77b4", "RA" = "#2ca02c", "OA" = "#ff7f0e"),
                     labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA")) +
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
    strip.placement  = "outside",
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin     = margin(0, 0, 0, 0))
ggsave("99_Fig/fig2/nonOA_OA_RA_vessel_num.png", plot = p4, width = 3, height = 1.25)
ggsave("99_Fig/fig2/nonOA_OA_RA_vessel_num.pdf", plot = p4, width = 3, height = 1.25)

p5 <- ggplot(df_long_ratio_vessel, aes(x = Diagnosis, y = Value, color = Diagnosis)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1) +
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels_ratio), strip.position = "bottom", ncol = 3) +
  labs(x = "", y = "Number/Area") +
  scale_color_manual(values = c("nonOA" = "#1f77b4", "RA" = "#2ca02c", "OA" = "#ff7f0e"),
                     labels = c("nonOA" = "ACLR", "RA" = "RA", "OA" = "OA")) +
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
    strip.placement  = "outside",
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin     = margin(0, 0, 0, 0))
ggsave("99_Fig/fig2/nonOA_OA_RA_vessel_ratio.png", plot = p5, width = 3, height = 1.25)
ggsave("99_Fig/fig2/nonOA_OA_RA_vessel_ratio.pdf", plot = p5, width = 3, height = 1.25)


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
#  [1] dunn.test_1.3.6 lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0 dplyr_1.1.4    
# 
# loaded via a namespace (and not attached):
#  [1] gtable_0.3.6      compiler_4.3.3    tidyselect_1.2.1  systemfonts_1.2.1 scales_1.3.0      textshaping_1.0.0 R6_2.6.1          labeling_0.4.3    generics_0.1.3    munsell_0.5.1     pillar_1.10.1     tzdb_0.5.0        rlang_1.1.5       stringi_1.8.7     timechange_0.3.0 
# [16] cli_3.6.4         withr_3.0.2       magrittr_2.0.3    grid_4.3.3        rstudioapi_0.17.1 hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5       glue_1.8.0        farver_2.1.2      ragg_1.3.3        colorspace_2.1-1  tools_4.3.3       pkgconfig_2.0.3  

