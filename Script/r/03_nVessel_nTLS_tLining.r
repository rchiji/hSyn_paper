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
file_paths <- list.files(path = "00_src/Annotation_ratio/", full.names = TRUE)
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

df_OA_nonOA_RA <- df[!df$Diagnosis %in% c("SLE", "SSc"), ]
df_OA_nonOA_RA <- df_OA_nonOA_RA[df_OA_nonOA_RA$Joint == "Knee",]

component <- colnames(df_OA_nonOA_RA[,37:43])


# statistics
## shapiro_test
shapiro_test_results <- df_OA_nonOA_RA %>% summarise(across(37:43, ~ shapiro.test(.)$p.value))
shapiro_test_results 
# Micro_vessel_number Large_vessel_number   TLS_number Lining_thickness Micro_vessel_ratio Large_vessel_ratio    TLS_ratio
# 1         3.37017e-12        1.339206e-12 1.521199e-20       0.00695195         0.01054021       1.603972e-12 2.790525e-20

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
# ** Significant Dunn Post-hoc Test Results for Micro_vessel_number **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00000 - Z-value: -5.32507
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00014 - Z-value: -3.74517
# =============================================
# ** Significant Dunn Post-hoc Test Results for Large_vessel_number **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.01021 - Z-value: -2.70601
# >> Comparison: nonOA - RA - Adjusted P-value: 0.05905 - Z-value: -1.75807
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS_number **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00029 - Z-value: -3.54653
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00014 - Z-value: -3.90987
# >> Comparison: OA - RA    - Adjusted P-value: 0.03925 - Z-value: -1.75942
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining_thickness **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00051 - Z-value: -3.39952
# >> Comparison: OA - RA    - Adjusted P-value: 0.00059 - Z-value: 3.54483
# =============================================
# ** Significant Dunn Post-hoc Test Results for TLS_ratio **
# =============================================
# >> Comparison: nonOA - OA - Adjusted P-value: 0.00028 - Z-value: -3.55432
# >> Comparison: nonOA - RA - Adjusted P-value: 0.00017 - Z-value: -3.85346
# >> Comparison: OA - RA    - Adjusted P-value: 0.04599 - Z-value: -1.68508
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

