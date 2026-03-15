library(dplyr)
library(dunn.test)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(WGCNA)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
anno <- read.delim(file = "00_src/Entropy_group.txt", sep = "\t", row.names = 1)


# Tissue proportion
df$Entropy <-  anno$Entropy
df_OA <- df[df$Diagnosis == "OA",]

component <- colnames(df_OA[,6:17])

## statistics
### shapiro_test
shapiro_test_results <- df_OA %>% 
  summarise(across(6:17, ~ shapiro.test(.)$p.value))
shapiro_test_results 
#        Adipose Fibrous_tissue__dense_irregular Fibrous_tissue__dense_regular Fibrous_tissue__loose          TLS       Plasma       Stroma       Lining       Muscle          RBC Micro_vessel Large_vessel
# 1 5.564255e-05                    3.078382e-08                  8.143933e-16            0.04639187 2.087342e-18 3.910671e-16 1.158127e-19 3.342318e-10 7.259322e-20 6.129428e-14    0.1940028 2.524883e-11

### kruskal.test
results <- list()
for (comp in component) {
  kw_result <- kruskal.test(as.formula(paste(comp, "~ Entropy")), data = df_OA)
  dunn_result <- dunn.test(df_OA[[comp]], df_OA$Entropy, method = "bh")
  significant_idx <- which(dunn_result$P.adjusted < 0.1)
  significant_comparisons <- dunn_result$comparisons[significant_idx]
  significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
  significant_z <- dunn_result$Z[significant_idx]
  results[[comp]] <- list(
    test_type = "Kruskal-Wallis",
    main_test = kw_result,
    posthoc_test = dunn_result,
    significant_comparisons = significant_comparisons,
    significant_p_adjusted = significant_p_adjusted,
    significant_z = significant_z
  )
}

for (comp in component) {
  cat("\n--- Results for", comp, "---\n")
  print(results[[comp]]$main_test)
  if (length(results[[comp]]$significant_comparisons) > 0) {
    cat("\n** Significant Dunn Post-hoc Test Results for", comp, "**\n")
    cat("=============================================\n")
    for (i in seq_along(results[[comp]]$significant_comparisons)) {
      cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
                  results[[comp]]$significant_comparisons[i], 
                  results[[comp]]$significant_p_adjusted[i],
                  results[[comp]]$significant_z[i]))
    }
    cat("=============================================\n")
  }
}

# Notes: Significant differences are shown.
# ** Significant Dunn Post-hoc Test Results for TLS **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.00024 - Z-value: 3.77233
# >> Comparison: High - Mid - Adjusted P-value: 0.00637 - Z-value: 2.63168
# >> Comparison: Low - Mid  - Adjusted P-value: 0.01354 - Z-value: -2.21029
# =============================================
# ** Significant Dunn Post-hoc Test Results for Plasma **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.00808 - Z-value: 2.78300
# >> Comparison: High - Mid - Adjusted P-value: 0.07589 - Z-value: 1.43329
# >> Comparison: Low - Mid  - Adjusted P-value: 0.03661 - Z-value: -1.97022
# =============================================
# ** Significant Dunn Post-hoc Test Results for Stroma **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.05975 - Z-value: 1.75261
# >> Comparison: High - Mid - Adjusted P-value: 0.00623 - Z-value: 2.86641
# =============================================
# ** Significant Dunn Post-hoc Test Results for Lining **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.00000 - Z-value: 5.32515
# >> Comparison: High - Mid - Adjusted P-value: 0.00000 - Z-value: 4.54977
# >> Comparison: Low - Mid  - Adjusted P-value: 0.00520 - Z-value: -2.56228
# =============================================
# ** Significant Dunn Post-hoc Test Results for RBC **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.00032 - Z-value: 3.51974
# >> Comparison: High - Mid - Adjusted P-value: 0.00003 - Z-value: 4.30143
# =============================================
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.00000 - Z-value: 5.09741
# >> Comparison: High - Mid - Adjusted P-value: 0.00001 - Z-value: 4.39789
# >> Comparison: Low - Mid  - Adjusted P-value: 0.00767 - Z-value: -2.42416
# =============================================


## Plot
custom_labels <- c(
  "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense, irregular)",
  "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense, regular)",
  "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
  "Micro_vessel" = "Micro vessel",
  "Large_vessel" = "Large vessel"
)

df_log <- df_OA[,6:17]
df_log <- log10(df_log + 1e-6)

df_log$Entropy <- df_OA$Entropy
df_log$Entropy <- factor(df_log$Entropy, levels = c("Low", "Mid", "High"))

df_log_long <- df_log %>%
  pivot_longer(cols = all_of(component),
               names_to = "Tissue",
               values_to = "Value")
colnames(df_log_long) <- c("Entropy", "Tissue", "Value")

df_log_long$Tissue <- factor(df_log_long$Tissue, levels = c("Fibrous_tissue__dense_irregular", "Fibrous_tissue__dense_regular", "Fibrous_tissue__loose", "Lining", "TLS", "Plasma", "Stroma", "RBC", "Micro_vessel", "Large_vessel", "Adipose",   "Muscle"))

p <- ggplot(df_log_long, aes(x = Entropy, y = Value, color = Entropy)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.1) +
  facet_wrap(~ Tissue, scales = "free", labeller = labeller(Tissue = custom_labels), strip.position = "bottom", ncol = 6) +
  labs(y = "log10(Proportion)") +
  scale_color_manual(values = c("Low" = "#bfe6bf", "Mid" = "#66cc66", "High" = "#336633")) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 5, margin = margin(t = 0)),
    strip.background = element_blank(),
    strip.placement  = "outside",
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, , margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin     = margin(0, 0, 0, 0))

ggsave("99_Fig/sup_fig4/Entropy_tissue_proportion.png", plot = p , width = 6, height = 2)
ggsave("99_Fig/sup_fig4/Entropy_tissue_proportion.pdf", plot = p , width = 6, height = 2)



# Gene expression
data <- read.delim("01_formatted/gene_read_count_preprocessed.txt", sep = "\t", header = TRUE, row.names = 1)

Min.CPM = 1
nSamples = 80
data_filtered <- round(data,0)
data_filtered <- data_filtered[which(apply(cpm(DGEList(counts = data_filtered)), 1, function(y) {sum(y>=Min.CPM)}) >= nSamples/2), ] # 13692 genes

data_norm <- data_filtered
tmp <- rep("A",dim(data_norm)[2]); tmp[1] <- "B"
colData <- cbind(colnames(data_norm), tmp); colnames(colData) = c("sample", "groups")
dds <- DESeqDataSetFromMatrix(countData = data_norm, colData = colData, design = ~ groups)

data_norm <- vst(dds)
data_norm <- assay(data_norm)

data_norm <- data_norm[order(apply(data_norm,1,sd), decreasing = T),]
data_sd_5000 <- as.matrix(data_norm[1:5000,])

data_wgcna <- t(as.matrix(data_sd_5000))


## WGCNA
powers  <-  c(1:30)
sft <- pickSoftThreshold(data_wgcna, powerVector = powers, networkType = "signed", verbose = 5)

# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"));
#  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red");
#  abline(h=0.90,col="red")
# plot(sft$fitIndices[,1], sft$fitIndices[,5],　xlab="Soft Threshold (power)",　ylab="Mean Connectivity", type="n",　main = paste("Mean connectivity"));
#  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

softPower <- 13
adjacency <- adjacency(data_wgcna, power = softPower, type="signed")
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1-TOM

k = as.vector(apply(adjacency, 2, sum, na.rm=T))
# hist(k)

scaleFreePlot(k, main="Check scale free topology\n")
#   scaleFreeRsquared slope
# 1              0.89 -1.59

geneTree <- hclust(as.dist(dissTOM), method = "average")

minModuleSize <- 30

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 1, pamStage = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 1055  713  597  532  422  340  318  224  221  164  119  114   73   62   46 

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# black        blue       brown        cyan       green greenyellow        grey     magenta        pink      purple         red      salmon         tan   turquoise      yellow 
#   224         597         532          46         340         114        1055         164         221         119         318          62          73         713         422 

pdf("99_Fig/sup_fig4/DendroAndColors.pdf", width = 4.5, height = 3, useDingbats = FALSE)
par(bg = "white")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()


MEList <- moduleEigengenes(data_wgcna, colors = dynamicColors)
MEs <- MEList$eigengenes
MEs_avr <- MEList$averageExpr

gene_module <- cbind(data_sd_5000, dynamicColors)

# write.table(gene_module, "03_res/GeneExpression/wgcna_gene_module.txt", sep = "\t", row.names = T, col.names = NA)
# write.table(MEs, "03_res/GeneExpression/wgcna_MEs.txt", sep = "\t", row.names = T, col.names = NA)

common <- intersect(rownames(MEs), rownames(anno))
MEs <- MEs[common, ]
MEs$Entropy <- anno[intersect(rownames(MEs), rownames(anno)), "Entropy"]


## statistics
### shapiro_test
shapiro_test_results <- MEs %>% 
  summarise(across(1:15, ~ shapiro.test(.)$p.value))
shapiro_test_results 
#        MEblack       MEblue      MEbrown       MEcyan      MEgreen MEgreenyellow     MEgrey   MEmagenta    MEpink     MEpurple      MEred     MEsalmon       MEtan MEturquoise     MEyellow
# 1 9.196207e-06 7.564519e-05 2.910548e-10 3.119264e-16 6.807494e-09  6.131658e-05 0.05574899 0.000847163 0.1132316 0.0001065685 0.04963042 4.400191e-05 0.000771942   0.8518036 9.121651e-05

### kruskal.test
module <- colnames(MEs[,1:15])
results <- list()
for (mod in module) {
  kw_result <- kruskal.test(as.formula(paste(mod, "~ Entropy")), data = MEs)
  dunn_result <- dunn.test(MEs[[mod]], MEs$Entropy, method = "bh")
  significant_idx <- which(dunn_result$P.adjusted < 0.1)
  significant_comparisons <- dunn_result$comparisons[significant_idx]
  significant_p_adjusted <- dunn_result$P.adjusted[significant_idx]
  significant_z <- dunn_result$Z[significant_idx]
  results[[mod]] <- list(
    test_type = "Kruskal-Wallis",
    main_test = kw_result,
    posthoc_test = dunn_result,
    significant_comparisons = significant_comparisons,
    significant_p_adjusted = significant_p_adjusted,
    significant_z = significant_z
  )
}

for (mod in module) {
  cat("\n--- Results for", mod, "---\n")
  print(results[[mod]]$main_test)
  if (length(results[[mod]]$significant_comparisons) > 0) {
    cat("\n** Significant Dunn Post-hoc Test Results for", mod, "**\n")
    cat("=============================================\n")
    for (i in seq_along(results[[mod]]$significant_comparisons)) {
      cat(sprintf(">> Comparison: %-10s - Adjusted P-value: %.5f - Z-value: %.5f\n", 
                  results[[mod]]$significant_comparisons[i], 
                  results[[mod]]$significant_p_adjusted[i],
                  results[[mod]]$significant_z[i]))
    }
    cat("=============================================\n")
  }
}

# Notes: Significant differences are shown.
# ** Significant Dunn Post-hoc Test Results for MEgreen **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.04965 - Z-value: 2.13084
# >> Comparison: Low - Mid  - Adjusted P-value: 0.05521 - Z-value: -1.78903
# =============================================
# ** Significant Dunn Post-hoc Test Results for MEgreenyellow **
# =============================================
# >> Comparison: High - Mid - Adjusted P-value: 0.04476 - Z-value: -1.88311
# >> Comparison: Low - Mid  - Adjusted P-value: 0.07909 - Z-value: -1.93717
# =============================================
# ** Significant Dunn Post-hoc Test Results for MEgrey **
# =============================================
# >> Comparison: High - Mid - Adjusted P-value: 0.04484 - Z-value: 2.17152
# =============================================
# ** Significant Dunn Post-hoc Test Results for MEpink **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.00814 - Z-value: -2.78059
# >> Comparison: High - Mid - Adjusted P-value: 0.09829 - Z-value: -1.29138
# >> Comparison: Low - Mid  - Adjusted P-value: 0.02956 - Z-value: 2.05988
# =============================================
# ** Significant Dunn Post-hoc Test Results for MEpurple **
# =============================================
# >> Comparison: High - Mid - Adjusted P-value: 0.09147 - Z-value: 1.54662
# =============================================
# ** Significant Dunn Post-hoc Test Results for MEred **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.06556 - Z-value: 1.70924
# =============================================
# ** Significant Dunn Post-hoc Test Results for MEsalmon **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.01689 - Z-value: -2.28142
# >> Comparison: Low - Mid  - Adjusted P-value: 0.01774 - Z-value: 2.51736
# =============================================
# ** Significant Dunn Post-hoc Test Results for MEyellow **
# =============================================
# >> Comparison: High - Low - Adjusted P-value: 0.03083 - Z-value: -2.31609
# >> Comparison: High - Mid - Adjusted P-value: 0.09453 - Z-value: -1.31335
# >> Comparison: Low - Mid  - Adjusted P-value: 0.08766 - Z-value: 1.56802
# =============================================


## plot
MEs_long <- MEs %>%
  pivot_longer(cols = all_of(module),
               names_to = "Module",
               values_to = "Value")

MEs_long$Entropy <- factor(MEs_long$Entropy, levels = c("Low", "Mid", "High"))

p <- ggplot(MEs_long, aes(x = Entropy, y = Value, colour = Entropy)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.1, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.1) +
  facet_wrap(~ Module, scales = "free", strip.position = "bottom", ncol = 5) +  # Adjust number of columns
  labs(y = "Module Eigengenes") +
  scale_color_manual(values = c("Low" = "#bfe6bf", "Mid" = "#66cc66", "High" = "#336633")) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 5, margin = margin(t = 0)),
    strip.background = element_blank(),
    strip.placement  = "outside",
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5, , margin = margin(l = -5)),
    legend.key.height = unit(0.25, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin     = margin(0, 0, 0, 0))

ggsave("99_Fig/sup_fig4/Entropy_module_MEs.png", plot = p , width = 5, height = 3.5)
ggsave("99_Fig/sup_fig4/Entropy_module_MEs.pdf", plot = p , width = 5, height = 3.5)


## Enrichment
gene_module <- as.data.frame(gene_module)
gene_module$gene <- rownames(gene_module)
gene_module <- separate(gene_module, gene, c("gene_id", "gene_name"), sep="[_]", extra="merge")

module_colors <- c("black", "blue", "brown", "cyan", "green", "greenyellow", "grey", "magenta", "pink", "purple", "red", "salmon", "tan", "turquoise", "yellow")

### GO
msigdb_go <- msigdbr(species = "Homo sapiens", collection = "C5")
msigdb_go <- msigdb_go %>%
  filter(!gs_subcollection %in% c("HPO", "GO:CC", "GO:MF"))

enrichment_results_go <- list() 
for (color in module_colors) {
  gene_data <- unique(gene_module$gene_name[gene_module$dynamicColors == color])
  enrichment_results_go[[color]] <- enricher(gene = gene_data,
                                          TERM2GENE = msigdb_go[, c("gs_name", "gene_symbol")])
}

# saveRDS(enrichment_results_go, "03_res/GeneExpression/wgcna_gene_enrichment_go.rds")

results_go <- lapply(enrichment_results_go, function(x) x@result)

results_go <- lapply(results_go, function(df) {
  df[df$qvalue < 0.05, , drop = FALSE]
})

for (color in names(results_go)) {
  results_go[[color]] <- results_go[[color]] %>%
    mutate(color = color,
           log_qvalue = -log10(qvalue))
}

df_results_go <- bind_rows(results_go)

df_results_go_top <- df_results_go %>%
  group_by(color) %>%
  arrange(qvalue) %>%
  slice_head(n = 5) %>%
  ungroup()


df_results_go_top$Description_color <- interaction(df_results_go_top$Description, df_results_go_top$color)
df_results_go_top$Description_color <- factor(df_results_go_top$Description_color, levels = rev(unique(df_results_go_top$Description_color)))

custom_label <- setNames(gsub("\\..*$", "", df_results_go_top$Description_color), df_results_go_top$Description_color)

p <- ggplot(df_results_go_top, aes(x = Description_color, y = log_qvalue)) +
  geom_bar(stat = "identity", aes(fill = color), width = 0.1) +
  geom_point(aes(color = color), shape = 1, size = 2, show.legend = FALSE) +
  coord_flip() +
  labs(y = "-log10(q-value)") +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_discrete(labels = custom_label) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    panel.grid.major = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 5), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 5))

ggsave("99_Fig/sup_fig4/Entropy_enrichment_go.png", plot = p , width = 10, height = 4.5)
ggsave("99_Fig/sup_fig4/Entropy_enrichment_go.pdf", plot = p , width = 10, height = 4.5)

p <- ggplot(df_results_go_top[df_results_go_top$color == "yellow",], aes(x = Description_color, y = log_qvalue)) +
  geom_bar(stat = "identity", aes(fill = color), width = 0.05) +
  geom_point(aes(color = color), shape = 1, size = 2, show.legend = FALSE) +
  coord_flip() +
  labs(y = "-log10(q-value)") +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_discrete(labels = custom_label) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    panel.grid.major = element_line(linewidth = 0.1),
    axis.line = element_line(linewidth = 0.1),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 5), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 5))

ggsave("99_Fig/sup_fig4/Entropy_enrichment_go_yellow.png", plot = p , width = 4, height = 1)
ggsave("99_Fig/sup_fig4/Entropy_enrichment_go_yellow.pdf", plot = p , width = 4, height = 1)

