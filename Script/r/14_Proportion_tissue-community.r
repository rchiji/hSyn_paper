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

