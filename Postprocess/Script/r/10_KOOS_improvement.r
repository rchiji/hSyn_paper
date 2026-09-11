library(dplyr)
library(tibble)
library(tidyverse)
library(ggplot2)

df <- read.delim(file = "01_formatted/annotations_full_tissue_proportion_cluster.txt", sep = "\t", row.names = 1)   


df_OA <- df %>%
  filter(Diagnosis == "OA") %>%
  mutate(
    cluster = factor(cluster, levels = c("1", "2", "3"))
  )

df_symptom <- df_OA %>%
  rownames_to_column("Donor") %>%
  filter(
    !is.na(Pre__KOOS_Symptom),
    !is.na(Post_12Month__KOOS_Symptom)
  ) %>%
  select(
    Donor,
    cluster,
    Pre__KOOS_Symptom,
    Post_3Month__KOOS_Symptom,
    Post_12Month__KOOS_Symptom
  ) %>%
  mutate(
    ImprovementValue =
      Post_12Month__KOOS_Symptom - Pre__KOOS_Symptom,
    Outcome = case_when(
      Post_12Month__KOOS_Symptom >= 80.5 & ImprovementValue >= 9 ~ "Improved", # Connelly, JW. et al. J Bone Joint Surg Am. 101, 995-1003 (2019) / Migliorini, F. et al. Knee Surg Relat Res. 36, 3 (2024).
      Post_12Month__KOOS_Symptom >= 80.5 & ImprovementValue < 9 ~ "Favorable baseline",
      Post_12Month__KOOS_Symptom < 80.5 ~ "No meaningful improvement"
    )
  )

df_symptom_long <- df_symptom %>%
  pivot_longer(
    cols = c(
      Pre__KOOS_Symptom,
      Post_12Month__KOOS_Symptom
    ),
    names_to = "Time",
    values_to = "KOOS"
  ) %>%
  mutate(
    Time = factor(
      Time,
      levels = c(
        "Pre__KOOS_Symptom",
        "Post_12Month__KOOS_Symptom"
      ),
      labels = c(
        "Pre",
        "Post 12-Month"
      )
    ),
    Outcome = factor(
      Outcome,
      levels = c(
        "Improved",
        "Favorable baseline",
        "No meaningful improvement"
      )
    )
  )

plot_symptom_trajectory <- ggplot(df_symptom_long, aes(x = Time, y = KOOS, group = Donor, color = Outcome, linetype = Outcome)) +
  geom_line(linewidth = 0.1, alpha = 0.6) +
  geom_point(size = 0.05) +
  geom_hline(yintercept = 80.5, linetype = "dotted", linewidth = 0.1, color = "grey40") +  
  facet_wrap(~ cluster, nrow = 1) +
  scale_color_manual(values = c("Improved" = "black", "Favorable baseline" = "black", "No meaningful improvement" = "red")) +
  scale_linetype_manual(values = c("Improved" = "solid", "Favorable baseline" = "dashed", "No meaningful improvement" = "solid")
  ) +
  labs(x = NULL, y = "KOOS Symptom", color = "Outcome") +
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 5),
        strip.background = element_blank(),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5, margin = margin(l = 0)),
        legend.key.height = unit(0.25, "cm"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin     = margin(0, 0, 0, 0)
      )

ggsave("99_Fig/fig2/KOOS_symptom_improvement_2point.png", plot = plot_symptom_trajectory, width = 2.75, height = 1.25)
ggsave("99_Fig/fig2/KOOS_symptom_improvement_2point.pdf", plot = plot_symptom_trajectory, width = 2.75, height = 1.25)


df_symptom_bar <- df_symptom %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = c(
        "Improved",
        "Favorable baseline",
        "No meaningful improvement"
      )
    )
  ) %>%
  count(cluster, Outcome) %>%
  group_by(cluster) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  ungroup()

plot_symptom_bar <- ggplot(df_symptom_bar, aes(x = cluster, y = proportion, fill = Outcome)) +
  geom_col(width = 0.9, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = c(
    "Improved" = "black",
    "Favorable baseline" = "grey",
    "No meaningful improvement" = "red"
  )) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Cluster", y = "Patients (%)", fill = "Outcome") +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(t = -1, b = -1, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, -2, unit = "mm")
  )

ggsave("99_Fig/fig2/KOOS_symptom_outcome_by_cluster.png", plot = plot_symptom_bar, width = 2.2, height = 1)
ggsave("99_Fig/fig2/KOOS_symptom_outcome_by_cluster.pdf", plot = plot_symptom_bar, width = 2.2, height = 1)


df_symptom <- df_symptom %>%
  mutate(
    PoorOutcome = ifelse(
      Outcome == "No meaningful improvement",
      "No meaningful improvement",
      "Other"
    )
  )
table(df_symptom$cluster, df_symptom$PoorOutcome)
#   No meaningful improvement  Other
# 1                         6     30
# 2                         5     18
# 3                         6      6
tab <- table(
  df_symptom$cluster,
  df_symptom$PoorOutcome
)
fisher.test(tab)
# Fisher's Exact Test for Count Data
# 
# data:  tab_poor
# p-value = 0.07154
# alternative hypothesis: two.sided

clusters <- levels(df_symptom$cluster)
pair_results <- data.frame()

for (i in 1:(length(clusters) - 1)) {
  for (j in (i + 1):length(clusters)) {
    dat_sub <- df_symptom %>%
      filter(cluster %in% c(clusters[i], clusters[j])) %>%
      droplevels()
    tab <- table(
      dat_sub$cluster,
      dat_sub$PoorOutcome
    )
    test <- fisher.test(tab)
    pair_results <- rbind(
      pair_results,
      data.frame(
        comparison = paste(clusters[i], "vs", clusters[j]),
        odds_ratio = unname(test$estimate),
        p_value = test$p.value
      )
    )
  }
}

pair_results$p_adjusted <- p.adjust(
  pair_results$p_value,
  method = "BH"
)
pair_results
#   comparison  odds_ratio     p_value  p_adjusted
# 1     1 vs 2   0.7241295  0.73588366   0.7358837
# 2     1 vs 3   0.2085206  0.04852673   0.1455802
# 3     2 vs 3   0.2894087  0.12971280   0.1945692

