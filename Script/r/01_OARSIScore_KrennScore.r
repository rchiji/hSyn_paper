library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)

data <- read.delim("00_src/annotations_full_OARSI_Krenn.txt", sep = "\t", header = T, row.names = 1)
data_OA <- data[data$Diagnosis == "OA",]

# cartilage
OARSI_data <- as.data.frame(data_OA$OARSI_score)

colnames(OARSI_data) <- "OARSI_score"
rownames(OARSI_data) <- rownames(data_OA)
OARSI_data <- na.omit(OARSI_data)
OARSI_data$label <- " "

Plot_OARSI <- ggplot(OARSI_data, aes(x = label, y = OARSI_score)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1, color ="#009E73") +
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 6, 1)) +
  labs(y = "OARSI score") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5)
  )

ggsave("99_Fig/sup_fig1/OARSI_score.png", plot = Plot_OARSI , width = 1.25, height = 1.5)
ggsave("99_Fig/sup_fig1/OARSI_score.pdf", plot = Plot_OARSI , width = 1.25, height = 1.5)

# synovium
krenn_data <- data_OA[, c("Krenn_Lining", "Krenn_Stroma", "Krenn_Infiltrate")]

colnames(krenn_data) <- c("Lining layer\nthickness", 
                          "Sub-lining layer\ncellularity", 
                          "Immune cell\ninfiltration")

krenn_data_long <- krenn_data %>% 
  pivot_longer(cols = everything(), names_to = "Krenn_Type", values_to = "Score")

Plot_Krenn <- ggplot(krenn_data_long, aes(x = Krenn_Type, y = Score)) +
  geom_violin(trim = FALSE, width = 0.5, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.1, color ="#C9657B") +
  scale_y_continuous(limits = c(-1, 5), breaks = seq(0, 5, 1)) +
  labs(y = "Krenn Synovitis score") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5)
  )

ggsave("99_Fig/sup_fig1/Krenn_score.png", plot = Plot_Krenn, width = 2.75, height = 1.5)
ggsave("99_Fig/sup_fig1/Krenn_score.pdf", plot = Plot_Krenn, width = 2.75, height = 1.5)


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
#  [1] ggbeeswarm_0.7.2 lubridate_1.9.4  forcats_1.0.0    stringr_1.5.1    purrr_1.0.4      readr_2.1.5      tidyr_1.3.1      tibble_3.2.1     ggplot2_3.5.1    tidyverse_2.0.0  dplyr_1.1.4     
# 
# loaded via a namespace (and not attached):
#  [1] vctrs_0.6.5       cli_3.6.4         rlang_1.1.5       stringi_1.8.7     generics_0.1.3    glue_1.8.0        colorspace_2.1-1  hms_1.1.3         scales_1.3.0      grid_4.3.3        munsell_0.5.1     tzdb_0.5.0        lifecycle_1.0.4   vipor_0.4.7       compiler_4.3.3   
# [16] timechange_0.3.0  pkgconfig_2.0.3   rstudioapi_0.17.1 beeswarm_0.4.0    R6_2.6.1          tidyselect_1.2.1  pillar_1.10.1     magrittr_2.0.3    tools_4.3.3       withr_3.0.2       gtable_0.3.6

