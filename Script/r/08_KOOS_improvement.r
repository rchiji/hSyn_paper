library(dplyr)
library(tidyverse)
library(ggplot2)

df <- read.delim(file = "00_src/annotations_full.txt", sep = "\t", row.names = 1)   
df_OA <- df[df$Diagnosis == "OA",]


df_pain <- df_OA[,c(18,23,28)]
df_symptom <- df_OA[,c(19,24,29)]

df_pain <- na.omit(df_pain)
df_symptom <- na.omit(df_symptom)

df_pain <- df_pain %>% 
  mutate(ImprovementValue = Post_12Month__KOOS_Pain - Pre__KOOS_Pain,
         Improvement = ifelse(ImprovementValue < 10, "< 10", ">= 10"))

df_pain_long <- df_pain %>%
  rownames_to_column("Donor") %>%
  pivot_longer(cols = c("Pre__KOOS_Pain", "Post_3Month__KOOS_Pain", "Post_12Month__KOOS_Pain"),
               names_to = "KOOS", values_to = "Value")

df_pain_long$KOOS <- factor(df_pain_long$KOOS, 
                               levels = c("Pre__KOOS_Pain", "Post_3Month__KOOS_Pain", "Post_12Month__KOOS_Pain"))

p <- ggplot(df_pain_long, aes(x = KOOS, y = Value, group = Donor, color = Improvement)) +
  geom_line(linewidth = 0.1) +
  geom_point(size = 0.4, shape = 1) +
  labs(y = "KOOS Pain") +
  scale_color_manual(name = "Improvement", values = c("< 10" = "red", ">= 10" = "black")) +
  scale_x_discrete(labels = c("Pre", "Post 3-Month", "Post 12-Month")) +
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5, margin = margin(l = 0)),
        legend.key.height = unit(0.25, "cm"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin     = margin(0, 0, 0, 0))

ggsave("99_Fig/fig3/KOOS_pain_improvement_3point.png", plot = p, width = 2, height = 2.5)
ggsave("99_Fig/fig3/KOOS_pain_improvement_3point.pdf", plot = p, width = 2, height = 2.5)


df_symptom <- df_symptom %>%
  mutate(ImprovementRate = (Post_12Month__KOOS_Symptom - Pre__KOOS_Symptom) / Pre__KOOS_Symptom,
         Improvement = ifelse(ImprovementRate < 0.25, "< 0.25", ">= 0.25"))

df_symptom <- df_symptom %>% 
    mutate(ImprovementValue = Post_12Month__KOOS_Symptom - Pre__KOOS_Symptom,
           Improvement = ifelse(ImprovementValue < 10, "< 10", ">= 10"))

df_symptom_long <- df_symptom %>%
  rownames_to_column("Donor") %>%
  pivot_longer(cols = c("Pre__KOOS_Symptom", "Post_3Month__KOOS_Symptom", "Post_12Month__KOOS_Symptom"),
               names_to = "KOOS", values_to = "Value")

df_symptom_long$KOOS <- factor(df_symptom_long$KOOS, 
                               levels = c("Pre__KOOS_Symptom", "Post_3Month__KOOS_Symptom", "Post_12Month__KOOS_Symptom"))

p <- ggplot(df_symptom_long, aes(x = KOOS, y = Value, group = Donor, color = Improvement)) +
  geom_line(linewidth = 0.1) +
  geom_point(size = 0.4, shape = 1) +
  labs(y = "KOOS Symptom") +
  scale_color_manual(name = "Improvement", values = c("< 10" = "red", ">= 10" = "black")) +
  scale_x_discrete(labels = c("Pre", "Post 3-Month", "Post 12-Month")) +
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5, margin = margin(l = 0)),
        legend.key.height = unit(0.25, "cm"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin     = margin(0, 0, 0, 0))

ggsave("99_Fig/fig3/KOOS_symptom_improvement_3point.png", plot = p, width = 2, height = 2.5)
ggsave("99_Fig/fig3/KOOS_symptom_improvement_3point.pdf", plot = p, width = 2, height = 2.5)


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
#  [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0 dplyr_1.1.4    
# 
# loaded via a namespace (and not attached):
#  [1] gtable_0.3.6      compiler_4.3.3    tidyselect_1.2.1  systemfonts_1.2.1 scales_1.3.0      textshaping_1.0.0 R6_2.6.1          labeling_0.4.3    generics_0.1.3    munsell_0.5.1     pillar_1.10.1     tzdb_0.5.0        rlang_1.1.5       stringi_1.8.7     timechange_0.3.0 
# [16] cli_3.6.4         withr_3.0.2       magrittr_2.0.3    grid_4.3.3        rstudioapi_0.17.1 hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5       glue_1.8.0        farver_2.1.2      ragg_1.3.3        colorspace_2.1-1  tools_4.3.3       pkgconfig_2.0.3
