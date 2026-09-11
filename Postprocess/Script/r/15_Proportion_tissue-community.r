library(data.table)
library(dplyr)
library(tidyverse)
library(tibble)
library(ggplot2)

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
ggsave("99_Fig/fig4/Vessel_res0.05_community_ratio.png", plot = plot_cluster_ratio, width = 6.5, height = 1.5)
ggsave("99_Fig/fig4/Vessel_res0.05_community_ratio.pdf", plot = plot_cluster_ratio, width = 6.5, height = 1.5)



df_community <- df_wide
df_number_thickness <- read.delim("01_formatted/annotations_full_tissue_count_thickness.txt", sep = "\t", row.names = 1)

donor_keep <- intersect(rownames(df_community), rownames(df_number_thickness))

df_community <- df_community[donor_keep,]
df_number_thickness <- df_number_thickness[donor_keep,]

df_all <- df_community
colnames(df_all) <- c("community1","community0","community2")
df_all$lining_thickness <- df_number_thickness$Lining_thickness

shapiro.test(df_all$community0)
# Shapiro-Wilk normality test
# data:  df_all$community0
# W = 0.93446, p-value = 4.02e-05
shapiro.test(df_all$lining_thickness)
# Shapiro-Wilk normality test
# data:  df_all$lining_thickness
# W = 0.96529, p-value = 0.00577

res <- cor.test(df_all$community0, df_all$lining_thickness, method = "spearman")
res
# Spearman's rank correlation rho
# data:  df_all$community0 and df_all$lining_thickness
# S = 125453, p-value = 2.117e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4344259

p <- ggplot(df_all, aes(x = community0, y = lining_thickness)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "#7A0403", linewidth = 0.2) +
  labs(x = "Proportion of Community 0", y = "Lining mean thickness") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    axis.text.y  = element_text(size = 5),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )
ggsave("99_Fig/fig4/cor_community0_lining.png", plot = p, width = 2, height = 2)
ggsave("99_Fig/fig4/cor_community0_lining.pdf", plot = p, width = 2, height = 2)


df_all$micro_vessel_ratio <- df_number_thickness$Micro_vessel_ratio

shapiro.test(df_all$community0)
# Shapiro-Wilk normality test
# data:  df_all$community0
# W = 0.93446, p-value = 4.02e-05
shapiro.test(df_all$micro_vessel_ratio)
# Shapiro-Wilk normality test
# data:  df_all$micro_vessel_ratio
# W = 0.96717, p-value = 0.008162

res <- cor.test(df_all$community0, df_all$micro_vessel_ratio, method = "spearman")
res
# Spearman's rank correlation rho
# data:  df_all$community0 and df_all$micro_vessel_ratio
# S = 126195, p-value = 2.585e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4310807

p <- ggplot(df_all, aes(x = community0, y = micro_vessel_ratio)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "#7A0403", linewidth = 0.2) +
  labs(x = "Proportion of Community 0", y = "Micro vessel ratio") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    axis.text.y  = element_text(size = 5),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )
ggsave("99_Fig/fig4/cor_community0_microvessel.png", plot = p, width = 2, height = 2)
ggsave("99_Fig/fig4/cor_community0_microvessel.pdf", plot = p, width = 2, height = 2)



