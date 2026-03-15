library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)

data <- read.delim("01_formatted/annotations_full_OARSI_Krenn.txt", sep = "\t", header = T, row.names = 1)
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

