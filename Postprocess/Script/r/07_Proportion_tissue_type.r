library(dplyr)
library(tidyverse)
library(tibble)
library(ggplot2)

df <- read.delim(file = "00_src/annotations_full.txt", sep = "\t", row.names = 1)

custom_labels <- c(
     "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense, irregular)",
     "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense, regular)",
     "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
     "Micro_vessel" = "Micro vessel",
     "Large_vessel" = "Large vessel"
 )

tissue <- c(
  "Adipose",
  "Fibrous_tissue__dense_irregular",
  "Fibrous_tissue__dense_regular",
  "Fibrous_tissue__loose",
  "Lining",
  "TLS",
  "Plasma",
  "Stroma",
  "RBC",
  "Micro_vessel",
  "Large_vessel",
  "Muscle"
)

tissue_colors <- c(
  "Adipose" = "#808000",
  "Fibrous_tissue__dense_irregular" = "#67da17",
  "Fibrous_tissue__dense_regular" = "#ecffb3",
  "Fibrous_tissue__loose" = "#83edaa",
  "Lining" = "#44d2e5",
  "TLS" = "#a05aa0",
  "Plasma" = "#30d383",
  "Stroma" = "#698c69",
  "RBC" = "magenta",
  "Micro_vessel" = "#59220a",
  "Large_vessel" = "#f68d14",
  "Muscle" = "#b299dc"
)


df_plot <- df %>%
  rownames_to_column(var = "Donor") %>%
  select(Donor, Diagnosis, all_of(tissue))

diagnosis_order <- c("OA","nonOA","RA","SLE","SSc")

reorder_diagnosis <- function(data, group, cols) {
  tmp <- data %>% filter(Diagnosis == group)
  if (nrow(tmp) == 1) {
    return(tmp$Donor)
  }
  mat <- tmp %>%
    select(Donor, all_of(cols)) %>%
    column_to_rownames(var = "Donor") %>%
    as.matrix()
  hc <- hclust(dist(mat))
  rownames(mat)[hc$order]
}

sorted_donors <- unlist(lapply(diagnosis_order, function(x) {
  reorder_diagnosis(df_plot, x, tissue)
}))

df_summary <- df_plot %>%
  pivot_longer(
    cols = all_of(tissue),
    names_to = "TissueType",
    values_to = "Ratio"
  ) %>%
  mutate(
    Donor = factor(Donor, levels = sorted_donors),
    Diagnosis = factor(Diagnosis, levels = diagnosis_order),
    TissueType = factor(TissueType, levels = tissue)
  ) %>%
  arrange(Diagnosis, Donor)

plot_tissue_ratio <- ggplot(df_summary, aes(x = Donor, y = Ratio, fill = TissueType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Proportion of tissue type", fill = "Tissue type") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = tissue_colors,
                    labels = function(x) ifelse(x %in% names(custom_labels), custom_labels[x], x)) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.y = element_text(size = 6, margin = margin(r = 0, unit = "mm")),
    axis.text.y  = element_text(size = 5),
    legend.title = element_text(size = 5, margin = margin(b = 2)),
    legend.text  = element_text(size = 5, margin = margin(l = 0.5, unit = "mm")),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, -3, unit = "mm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "mm")
  )
ggsave("99_Fig/fig3/tissue_ratio_diagnosis_ordered.png", plot = plot_tissue_ratio, width = 6.5, height = 1.5)
ggsave("99_Fig/fig3/tissue_ratio_diagnosis_ordered.pdf", plot = plot_tissue_ratio, width = 6.5, height = 1.5)

