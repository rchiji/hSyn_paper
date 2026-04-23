library(dplyr)
library(tidyverse)
library(ggplot2)

df <- read.delim(file = "01_formatted/annotations_full_tissue_proportion_cluster.txt", sep = "\t", row.names = 1)   
df_OA <- df[df$Diagnosis == "OA",]


df_pain <- df_OA[,c(18,23,28)]
df_symptom <- df_OA[,c(19,24,29)]

df_pain <- na.omit(df_pain)
df_symptom <- na.omit(df_symptom)

cluster <- df_OA$cluster
names(cluster) <- rownames(df_OA)


df_pain <- df_pain %>%
  mutate(ImprovementValue = Post_12Month__KOOS_Pain - Pre__KOOS_Pain,
         Improvement_strict = ifelse(ImprovementValue < 10, "< 10", ">= 10"))

df_pain_long <- df_pain %>%
  rownames_to_column("Donor") %>%
  pivot_longer(cols = c("Pre__KOOS_Pain", "Post_3Month__KOOS_Pain", "Post_12Month__KOOS_Pain"),
               names_to = "KOOS", values_to = "Value")

df_pain_long$KOOS <- factor(df_pain_long$KOOS, 
                            levels = c("Pre__KOOS_Pain", "Post_3Month__KOOS_Pain", "Post_12Month__KOOS_Pain"))

p <- ggplot(df_pain_long, aes(x = KOOS, y = Value, group = Donor, color = Improvement_strict)) +
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

ggsave("99_Fig/sup_fig5/KOOS_pain_improvement_3point.png", plot = p, width = 2, height = 2.5)
ggsave("99_Fig/sup_fig5/KOOS_pain_improvement_3point.pdf", plot = p, width = 2, height = 2.5)

df_pain <- df_pain %>%
  mutate(Improvement = ifelse(ImprovementValue < 10, "Remained stable", "Improved"))

df_pain$cluster <- cluster[rownames(df_pain)]
df_stacked_bar <- df_pain
df_stacked_bar <- df_stacked_bar %>% 
  count(Improvement, cluster)

df_stacked_bar$cluster <- factor(df_stacked_bar$cluster, levels = c("1", "2", "3","4"))

plot_stacked_bar <- ggplot(df_stacked_bar, aes(x = Improvement, y = n, fill = cluster)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.9) +
  labs(y = "Fraction of cluster", fill = "cluster") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = c("1" = "#C2BAB4", "2" = "#FEC089", "3" = "#F06A00", "4" = "#8C2D04")) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 5),
    axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
    axis.text.y  = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 2)),
    legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, -3, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )

ggsave("99_Fig/fig3/KOOS_pain_improvement_cluster_proportion.png", plot = plot_stacked_bar, width = 1, height = 1.25)
ggsave("99_Fig/fig3/KOOS_pain_improvement_cluster_proportion.pdf", plot = plot_stacked_bar, width = 1, height = 1.25)


df_symptom <- df_symptom %>% 
  mutate(ImprovementValue = Post_12Month__KOOS_Symptom - Pre__KOOS_Symptom,
         Improvement_strict = ifelse(ImprovementValue < 10, "< 10", ">= 10"))

df_symptom_long <- df_symptom %>%
  rownames_to_column("Donor") %>%
  pivot_longer(cols = c("Pre__KOOS_Symptom", "Post_3Month__KOOS_Symptom", "Post_12Month__KOOS_Symptom"),
               names_to = "KOOS", values_to = "Value")

df_symptom_long$KOOS <- factor(df_symptom_long$KOOS, 
                               levels = c("Pre__KOOS_Symptom", "Post_3Month__KOOS_Symptom", "Post_12Month__KOOS_Symptom"))

p <- ggplot(df_symptom_long, aes(x = KOOS, y = Value, group = Donor, color = Improvement_strict)) +
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

ggsave("99_Fig/sup_fig5/KOOS_symptom_improvement_3point.png", plot = p, width = 2, height = 2.5)
ggsave("99_Fig/sup_fig5/KOOS_symptom_improvement_3point.pdf", plot = p, width = 2, height = 2.5)

df_symptom <- df_symptom %>% 
  mutate(Improvement = ifelse(ImprovementValue < 10, "Remained stable", "Improved"))

df_symptom$cluster <- cluster[rownames(df_symptom)]
df_stacked_bar <- df_symptom

df_stacked_bar <- df_stacked_bar %>% 
  count(Improvement, cluster)

df_stacked_bar$cluster <- factor(df_stacked_bar$cluster, levels = c("1", "2", "3","4"))

plot_stacked_bar <- ggplot(df_stacked_bar, aes(x = Improvement, y = n, fill = cluster)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.9) +
  labs(y = "Fraction of cluster", fill = "cluster") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = c("1" = "#C2BAB4", "2" = "#FEC089", "3" = "#F06A00", "4" = "#8C2D04")) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 5),
    axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm") ),
    axis.text.y  = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 2)),
    legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, -3, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )

ggsave("99_Fig/fig3/KOOS_symptom_improvement_cluster_proportion.png", plot = plot_stacked_bar, width = 1, height = 1.25)
ggsave("99_Fig/fig3/KOOS_symptom_improvement_cluster_proportion.pdf", plot = plot_stacked_bar, width = 1, height = 1.25)

