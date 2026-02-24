library(data.table)
library(dplyr)
library(tidyverse)
library(tibble)
library(ggplot2)

library(psych)
library(ggsci)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(readxl)
library(viridis)
library(reshape2)
library(patchwork)

turbo3_hex <- c("#30123B", "#FABA39", "#7A0403")
names(turbo3_hex) <- as.character(0:2)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)

file_paths <- list.files(path = "00_src/SLICTile_cluster_res0.05_vessel/", full.names = TRUE)


samplenames <- sapply(file_paths, function(x) {
  x <- basename(x)
  tools::file_path_sans_ext(x)
})

file_list <- lapply(file_paths, fread)
names(file_list) <- samplenames

file_list <- lapply(file_list, function(x){
  df <- as.data.frame(x)
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  df
})

file_list_add_name <- lapply(names(file_list), function(donor_name) {
  data <- file_list[[donor_name]]
  data$Donor <- donor_name
  return(data)
})

data_all <- do.call(rbind, file_list_add_name)


all_combinations <- expand.grid(Donor = unique(data_all$Donor),
                                Cluster = unique(data_all$Cluster))

df_summary <- all_combinations %>%
  left_join(data_all %>% group_by(Donor, Cluster) %>% summarise(Count = n(), .groups = "drop"),
            by = c("Donor", "Cluster")) %>%
  replace_na(list(Count = 0))

df_summary <- df_summary %>%
  group_by(Donor) %>%
  mutate(Ratio = Count / sum(Count))

df_wide <- df_summary %>%
  select(Donor, Cluster, Ratio) %>%
  pivot_wider(names_from = Cluster, values_from = Ratio, values_fill = list(Ratio = 0))
df_wide <- df_wide %>% 
  column_to_rownames(var = "Donor")

distance_matrix <- dist(df_wide) 
hc <- hclust(distance_matrix)

sorted_donors <- rownames(df_wide)[hc$order]

df_summary <- df_summary %>%
  mutate(Donor = factor(Donor, levels = sorted_donors)) %>%
  arrange(Donor)

df_summary$Cluster <- factor(df_summary$Cluster)

plot_cluster_ratio <- ggplot(df_summary, aes(x = Donor, y = Ratio, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Proportion of community type", fill = "Community") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = turbo3_hex) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
    axis.text.y  = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 2)),
    legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, -3, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )
ggsave("99_Fig/fig5/Vessel_res0.05_community_ratio.png", plot = plot_cluster_ratio, width = 6.5, height = 1.5)
ggsave("99_Fig/fig5/Vessel_res0.05_community_ratio.pdf", plot = plot_cluster_ratio, width = 6.5, height = 1.5)


sessionInfo()
R version 4.3.3 (2024-02-29 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
  [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8    LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                    LC_TIME=Japanese_Japan.utf8    

time zone: Asia/Tokyo
tzcode source: internal

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] lubridate_1.9.4   forcats_1.0.0     stringr_1.5.1     purrr_1.0.4       readr_2.1.5       tidyr_1.3.1       tibble_3.2.1      ggplot2_3.5.1     tidyverse_2.0.0   dplyr_1.1.4       data.table_1.17.0

loaded via a namespace (and not attached):
  [1] SummarizedExperiment_1.32.0 gtable_0.3.6                xfun_0.51                   Biobase_2.62.0              lattice_0.22-5              tzdb_0.5.0                  vctrs_0.6.5                 tools_4.3.3                 bitops_1.0-9               
[10] generics_0.1.3              stats4_4.3.3                parallel_4.3.3              pkgconfig_2.0.3             Matrix_1.6-5                S4Vectors_0.40.2            lifecycle_1.0.4             GenomeInfoDbData_1.2.11     farver_2.1.2               
[19] compiler_4.3.3              textshaping_1.0.0           statmod_1.5.0               munsell_0.5.1               DESeq2_1.42.1               codetools_0.2-19            GenomeInfoDb_1.38.8         RCurl_1.98-1.17             baySeq_2.36.2              
[28] pillar_1.10.1               crayon_1.5.3                BiocParallel_1.36.0         DelayedArray_0.28.0         limma_3.58.1                abind_1.4-8                 tidyselect_1.2.1            locfit_1.5-9.12             stringi_1.8.7              
[37] labeling_0.4.3              grid_4.3.3                  colorspace_2.1-1            cli_3.6.4                   SparseArray_1.2.4           magrittr_2.0.3              S4Arrays_1.2.1              edgeR_4.0.16                withr_3.0.2                
[46] scales_1.3.0                timechange_0.3.0            XVector_0.42.0              matrixStats_1.5.0           ROC_1.78.0                  ragg_1.3.3                  hms_1.1.3                   evaluate_1.0.3              knitr_1.50                 
[55] GenomicRanges_1.54.1        IRanges_2.36.0              rlang_1.1.5                 Rcpp_1.0.14                 glue_1.8.0                  BiocGenerics_0.48.1         TCC_1.42.0                  rstudioapi_0.17.1           R6_2.6.1                   
[64] systemfonts_1.2.1           MatrixGenerics_1.14.0       zlibbioc_1.48.2 

