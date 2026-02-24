library(compositions)
library(vegan)
library(ggplot2)
library(ggrepel)
library(patchwork)

df <- read.delim(file = "00_src/annotations_full.txt", sep = "\t", row.names = 1)

custom_labels <- c(
     "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense, irregular)",
     "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense, regular)",
     "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
     "Micro_vessel" = "Micro vessel",
     "Large_vessel" = "Large vessel"
 )

diagnosis_colors <- c("nonOA" = "#1f77b4",
                      "RA" = "#2ca02c",
                      "OA" = "#ff7f0e",
                      "SLE" = "grey80")

annotation_colors <- list(Diagnosis = diagnosis_colors)


df_knee <- df[df$Joint == "Knee",]

res_knee_clr <- clr(acomp(df_knee[,6:17])) 

knee_dist <- dist(res_knee_clr, method = "euclidean")

anno_knee <- df_knee[,1:5]

mds_knee_clr <- cmdscale(knee_dist, k = 2, eig = TRUE)

mdsData_knee <- as.data.frame(mds_knee_clr$points[,1:2])

mds_df_knee <- cbind(mdsData_knee, anno_knee)

# fit MDS model
env_fit_knee <- envfit(mds_knee_clr, res_knee_clr)

# correlation between dimension value and tissue ratio
res_knee_clr_df <- as.data.frame(res_knee_clr)

cor_df_knee <- data.frame(cor(res_knee_clr_df, mds_df_knee[, c(1, 2)]))

# adjust for plot scale
cor_df_knee <- cor_df_knee * ordiArrowMul(env_fit_knee)

cor_df_knee$distance <- sqrt(cor_df_knee$V1 ** 2 + cor_df_knee$V2 ** 2)

cor_df_knee_renamed <- cor_df_knee
rownames(cor_df_knee_renamed) <- gsub("Fibrous_tissue__dense_irregular", "Fibrous tissue (dense, irregular)", rownames(cor_df_knee_renamed))
rownames(cor_df_knee_renamed) <- gsub("Fibrous_tissue__dense_regular", "Fibrous tissue (dense, regular)", rownames(cor_df_knee_renamed))
rownames(cor_df_knee_renamed) <- gsub("Fibrous_tissue__loose", "Fibrous tissue (loose)", rownames(cor_df_knee_renamed))
rownames(cor_df_knee_renamed) <- gsub("Micro_vessel", "Micro vessel", rownames(cor_df_knee_renamed))
rownames(cor_df_knee_renamed) <- gsub("Large_vessel", "Large vessel", rownames(cor_df_knee_renamed))

# plot
MDSplot_knee <- ggplot(mds_df_knee, aes(x = V1, y = V2)) +
  geom_point(shape = 21, stroke = 0.1, alpha = 0.7, aes(fill = Diagnosis, size = Age, colour = Sex)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = diagnosis_colors,
                    labels = c("nonOA" = "ACLR", "OA" = "OA", "RA" = "RA", "SLE" = "SLE")) +
  scale_size_continuous(range = c(0.25, 1.5)) +
  labs(x = "MDS1", y = "MDS2") +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        axis.ticks   = element_line(linewidth = 0.1),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "lines"),
        legend.spacing.y = unit(0.05, "lines"),
        panel.grid = element_blank())
ggsave("99_Fig/fig2/betaDiversity_MDSplot.png", plot = MDSplot_knee , width = 2.5, height = 2)
ggsave("99_Fig/fig2/betaDiversity_MDSplot.pdf", plot = MDSplot_knee , width = 2.5, height = 2)

MDScorplot_knee <- ggplot(data = cor_df_knee_renamed, aes(x = 0, y = 0, xend = V1, yend = V2, colour = distance)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.2) +
  scale_color_gradient2(low = "gray",mid = "orange", high = "red") +
  geom_text_repel(data = cor_df_knee_renamed, size = 1.5, aes(x = V1, y = V2, label = rownames(cor_df_knee_renamed)),
                  color = "black", hjust = -0.2, vjust = -0.2,max.overlaps = Inf, segment.size = 0.1) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        axis.ticks   = element_line(linewidth = 0.1),
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.key.width  = unit(0.5, "lines"),
        panel.grid = element_blank())
ggsave("99_Fig/fig2/betaDiversity_MDSplot_Contribution.png", plot = MDScorplot_knee , width = 2.5, height = 2)
ggsave("99_Fig/fig2/betaDiversity_MDSplot_Contribution.pdf", plot = MDScorplot_knee , width = 2.5, height = 2)


sessionInfo()
# R version 4.3.3 (2024-02-29 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8    LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                    LC_TIME=Japanese_Japan.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.3.0    ggrepel_0.9.6      ggplot2_3.5.1      vegan_2.6-10       lattice_0.22-5     permute_0.9-7      compositions_2.0-8
# 
# loaded via a namespace (and not attached):
#   [1] Matrix_1.6-5        gtable_0.3.6        dplyr_1.1.4         compiler_4.3.3      tidyselect_1.2.1    Rcpp_1.0.14         bayesm_3.1-6        parallel_4.3.3      cluster_2.1.6       textshaping_1.0.0   systemfonts_1.2.1   splines_4.3.3       scales_1.3.0       
# [14] R6_2.6.1            labeling_0.4.3      generics_0.1.3      robustbase_0.99-4-1 MASS_7.3-60.0.1     tibble_3.2.1        munsell_0.5.1       pillar_1.10.1       rlang_1.1.5         cli_3.6.4           withr_3.0.2         magrittr_2.0.3      mgcv_1.9-1         
# [27] grid_4.3.3          rstudioapi_0.17.1   lifecycle_1.0.4     nlme_3.1-164        DEoptimR_1.1-3-1    vctrs_0.6.5         glue_1.8.0          tensorA_0.36.2.1    farver_2.1.2        ragg_1.3.3          colorspace_2.1-1    tools_4.3.3         pkgconfig_2.0.3   

