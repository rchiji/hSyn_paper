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
#        Adipose Fibrous_tissue__dense_irregular Fibrous_tissue__dense_regular Fibrous_tissue__loose          TLS       Plasma       Stroma       Lining       Muscle         RBC Micro_vessel Large_vessel
# 1 6.910185e-05                    2.351286e-08                  7.082465e-16            0.05572155 1.560488e-18 2.846042e-16 9.381429e-20 2.524819e-10 5.474948e-20 4.62055e-14    0.1670972 3.787407e-11

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

# --- Results for Adipose ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Adipose by Entropy
# Kruskal-Wallis chi-squared = 3.9371, df = 2, p-value = 0.1397
# 
# 
# ** Significant Dunn Post-hoc Test Results for Adipose **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.09854 - Z-value: -1.84053
# =============================================
#   
#   --- Results for Fibrous_tissue__dense_irregular ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Fibrous_tissue__dense_irregular by Entropy
# Kruskal-Wallis chi-squared = 1.3751, df = 2, p-value = 0.5028
# 
# 
# --- Results for Fibrous_tissue__dense_regular ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Fibrous_tissue__dense_regular by Entropy
# Kruskal-Wallis chi-squared = 0.48984, df = 2, p-value = 0.7828
# 
# 
# --- Results for Fibrous_tissue__loose ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Fibrous_tissue__loose by Entropy
# Kruskal-Wallis chi-squared = 0.32403, df = 2, p-value = 0.8504
# 
# 
# --- Results for TLS ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  TLS by Entropy
# Kruskal-Wallis chi-squared = 15.853, df = 2, p-value = 0.0003611
# 
# 
# ** Significant Dunn Post-hoc Test Results for TLS **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.00027 - Z-value: 3.74240
# >> Comparison: High - Mid - Adjusted P-value: 0.00478 - Z-value: 2.72776
# >> Comparison: Low - Mid  - Adjusted P-value: 0.01674 - Z-value: -2.12633
# =============================================
#   
#   --- Results for Plasma ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Plasma by Entropy
# Kruskal-Wallis chi-squared = 7.8818, df = 2, p-value = 0.01943
# 
# 
# ** Significant Dunn Post-hoc Test Results for Plasma **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.00811 - Z-value: 2.78158
# >> Comparison: High - Mid - Adjusted P-value: 0.07370 - Z-value: 1.44878
# >> Comparison: Low - Mid  - Adjusted P-value: 0.03696 - Z-value: -1.96620
# =============================================
#   
#   --- Results for Stroma ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Stroma by Entropy
# Kruskal-Wallis chi-squared = 7.7385, df = 2, p-value = 0.02087
# 
# 
# ** Significant Dunn Post-hoc Test Results for Stroma **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.06159 - Z-value: 1.73851
# >> Comparison: High - Mid - Adjusted P-value: 0.01098 - Z-value: 2.68199
# =============================================
#   
#   --- Results for Lining ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Lining by Entropy
# Kruskal-Wallis chi-squared = 36.001, df = 2, p-value = 1.522e-08
# 
# 
# ** Significant Dunn Post-hoc Test Results for Lining **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.00000 - Z-value: 5.35388
# >> Comparison: High - Mid - Adjusted P-value: 0.00000 - Z-value: 4.60447
# >> Comparison: Low - Mid  - Adjusted P-value: 0.00503 - Z-value: -2.57384
# =============================================
#   
#   --- Results for Muscle ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Muscle by Entropy
# Kruskal-Wallis chi-squared = 2.5574, df = 2, p-value = 0.2784
# 
# 
# --- Results for RBC ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  RBC by Entropy
# Kruskal-Wallis chi-squared = 22.93, df = 2, p-value = 1.049e-05
# 
# 
# ** Significant Dunn Post-hoc Test Results for RBC **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.00031 - Z-value: 3.53381
# >> Comparison: High - Mid - Adjusted P-value: 0.00002 - Z-value: 4.36519
# =============================================
#   
#   --- Results for Micro_vessel ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Micro_vessel by Entropy
# Kruskal-Wallis chi-squared = 33.137, df = 2, p-value = 6.374e-08
# 
# 
# ** Significant Dunn Post-hoc Test Results for Micro_vessel **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.00000 - Z-value: 5.10439
# >> Comparison: High - Mid - Adjusted P-value: 0.00001 - Z-value: 4.46219
# >> Comparison: Low - Mid  - Adjusted P-value: 0.00807 - Z-value: -2.40570
# =============================================
#   
#   --- Results for Large_vessel ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  Large_vessel by Entropy
# Kruskal-Wallis chi-squared = 0.52605, df = 2, p-value = 0.7687


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
data <- read.delim("00_src/gene_read_count.txt", sep = "\t", header = TRUE, row.names = 1)

data <- data[!is.na(data$Gene.type) & data$Gene.type == "protein_coding", ]
data$id <- rownames(data)
data <- mutate(data, gene = paste(data$id,data$gene_name, sep = "_"))
row.names(data) <- data$gene
data <- data[,2:87]
data <- na.omit(data) # 19933 genes

Min.CPM = 1
nSamples = 86
data_filtered <- round(data,0)
data_filtered <- data_filtered[which(apply(cpm(DGEList(counts = data_filtered)), 1, function(y) {sum(y>=Min.CPM)}) >= nSamples/2), ] # 13678 genes

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
# 1              0.91  -1.6

geneTree <- hclust(as.dist(dissTOM), method = "average")

minModuleSize <- 30

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 1, pamStage = FALSE, minClusterSize = minModuleSize)

table(dynamicMods)
dynamicMods_deepSplit1_pamStageFALSE
# dynamicMods
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 1089  652  580  468  446  311  307  251  233  179  174  111   88   63   48 

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# dynamicMods
# dynamicColors
# black        blue       brown        cyan       green greenyellow        grey     magenta        pink      purple         red      salmon         tan   turquoise      yellow 
#   251         580         468          48         311         111        1089         179         233         174         307          63          88         652         446


pdf("99_Fig/sup_fig4/DendroAndColors.pdf", width = 4.5, height = 3, useDingbats = FALSE)
par(bg = "white")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()


MEList <- moduleEigengenes(data_wgcna, colors = dynamicColors)
MEs <- MEList$eigengenes
MEs_avr <- MEList$averageExpr

gene_module <- cbind(data_sd_5000, dynamicColors)

# write.table(gene_module, "01_formatted/wgcna_gene_module.txt", sep = "\t", row.names = T, col.names = NA)
# write.table(MEs, "01_formatted/wgcna_MEs.txt", sep = "\t", row.names = T, col.names = NA)

common <- intersect(rownames(MEs), rownames(anno))
MEs <- MEs[common, ]
MEs$Entropy <- anno[intersect(rownames(MEs), rownames(anno)), "Entropy"]


## statistics
### shapiro_test
shapiro_test_results <- MEs %>% 
  summarise(across(1:15, ~ shapiro.test(.)$p.value))
shapiro_test_results 
# MEblack    MEblue      MEbrown       MEcyan    MEgreen MEgreenyellow   MEgrey    MEmagenta      MEpink    MEpurple       MEred    MEsalmon       MEtan  MEturquoise     MEyellow
# 1 0.06357255 0.8071947 5.547269e-05 2.555369e-16 0.08399078   0.002388821 0.828709 3.809294e-06 0.003232176 3.24445e-05 1.34634e-09 1.88453e-05 0.004754396 0.0007322564 5.549834e-10

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

# --- Results for MEblack ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEblack by Entropy
# Kruskal-Wallis chi-squared = 7.3672, df = 2, p-value = 0.02513
# 
# 
# ** Significant Dunn Post-hoc Test Results for MEblack **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.01384 - Z-value: -2.60355
# >> Comparison: High - Mid - Adjusted P-value: 0.06903 - Z-value: -1.68471
# >> Comparison: Low - Mid  - Adjusted P-value: 0.05154 - Z-value: 1.63009
# =============================================
#   
#   --- Results for MEblue ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEblue by Entropy
# Kruskal-Wallis chi-squared = 2.5346, df = 2, p-value = 0.2816
# 
# 
# --- Results for MEbrown ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEbrown by Entropy
# Kruskal-Wallis chi-squared = 5.5975, df = 2, p-value = 0.06089
# 
# 
# ** Significant Dunn Post-hoc Test Results for MEbrown **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.03421 - Z-value: -2.27663
# >> Comparison: Low - Mid  - Adjusted P-value: 0.07471 - Z-value: 1.44160
# =============================================
#   
#   --- Results for MEcyan ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEcyan by Entropy
# Kruskal-Wallis chi-squared = 6.4465, df = 2, p-value = 0.03983
# 
# 
# ** Significant Dunn Post-hoc Test Results for MEcyan **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.01943 - Z-value: -2.48510
# >> Comparison: Low - Mid  - Adjusted P-value: 0.01592 - Z-value: 2.30403
# =============================================
#   
#   --- Results for MEgreen ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEgreen by Entropy
# Kruskal-Wallis chi-squared = 3.4856, df = 2, p-value = 0.175
# 
# 
# ** Significant Dunn Post-hoc Test Results for MEgreen **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.09541 - Z-value: 1.85494
# >> Comparison: Low - Mid  - Adjusted P-value: 0.07986 - Z-value: -1.61424
# =============================================
#   
#   --- Results for MEgreenyellow ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEgreenyellow by Entropy
# Kruskal-Wallis chi-squared = 4.3598, df = 2, p-value = 0.1131
# 
# 
# ** Significant Dunn Post-hoc Test Results for MEgreenyellow **
#   =============================================
#   >> Comparison: High - Mid - Adjusted P-value: 0.06342 - Z-value: 2.03076
# =============================================
#   
#   --- Results for MEgrey ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEgrey by Entropy
# Kruskal-Wallis chi-squared = 0.68517, df = 2, p-value = 0.7099
# 
# 
# --- Results for MEmagenta ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEmagenta by Entropy
# Kruskal-Wallis chi-squared = 1.8627, df = 2, p-value = 0.394
# 
# 
# --- Results for MEpink ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEpink by Entropy
# Kruskal-Wallis chi-squared = 1.9286, df = 2, p-value = 0.3813
# 
# 
# --- Results for MEpurple ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEpurple by Entropy
# Kruskal-Wallis chi-squared = 6.64, df = 2, p-value = 0.03615
# 
# 
# ** Significant Dunn Post-hoc Test Results for MEpurple **
#   =============================================
#   >> Comparison: High - Mid - Adjusted P-value: 0.04075 - Z-value: -1.92416
# >> Comparison: Low - Mid  - Adjusted P-value: 0.04773 - Z-value: -2.14669
# =============================================
#   
#   --- Results for MEred ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEred by Entropy
# Kruskal-Wallis chi-squared = 3.7921, df = 2, p-value = 0.1502
# 
# 
# ** Significant Dunn Post-hoc Test Results for MEred **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.08639 - Z-value: 1.89877
# =============================================
#   
#   --- Results for MEsalmon ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEsalmon by Entropy
# Kruskal-Wallis chi-squared = 4.3734, df = 2, p-value = 0.1123
# 
# 
# ** Significant Dunn Post-hoc Test Results for MEsalmon **
#   =============================================
#   >> Comparison: High - Low - Adjusted P-value: 0.04559 - Z-value: -1.87508
# >> Comparison: Low - Mid  - Adjusted P-value: 0.05953 - Z-value: 2.05701
# =============================================
#   
#   --- Results for MEtan ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEtan by Entropy
# Kruskal-Wallis chi-squared = 2.3663, df = 2, p-value = 0.3063
# 
# 
# --- Results for MEturquoise ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEturquoise by Entropy
# Kruskal-Wallis chi-squared = 3.9509, df = 2, p-value = 0.1387
# 
# 
# --- Results for MEyellow ---
#   
#   Kruskal-Wallis rank sum test
# 
# data:  MEyellow by Entropy
# Kruskal-Wallis chi-squared = 2.9797, df = 2, p-value = 0.2254
  
  
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

# saveRDS(enrichment_results_go, "01_formatted/wgcna_gene_enrichment_go.rds")

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

p <- ggplot(df_results_go_top[df_results_go_top$color == "brown",], aes(x = Description_color, y = log_qvalue)) +
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

ggsave("99_Fig/sup_fig4/Entropy_enrichment_go_brown.png", plot = p , width = 4, height = 1)
ggsave("99_Fig/sup_fig4/Entropy_enrichment_go_brown.pdf", plot = p , width = 4, height = 1)


### cell type
msigdb_celltype <- msigdbr(species = "Homo sapiens", collection = "C8")

enrichment_results_celltype <- list() 
for (color in module_colors) {
  gene_data <- unique(gene_module$gene_name[gene_module$dynamicColors == color])
  enrichment_results_celltype[[color]] <- enricher(gene = gene_data,
                                             TERM2GENE = msigdb_celltype[, c("gs_name", "gene_symbol")])
}

# saveRDS(enrichment_results_celltype, "01_formatted/wgcna_gene_enrichment_celltype.rds")

results_celltype <- lapply(enrichment_results_celltype, function(x) x@result)

results_celltype <- lapply(results_celltype, function(df) {
  df[df$qvalue < 0.05, , drop = FALSE]
})

for (color in names(results_celltype)) {
  results_celltype[[color]] <- results_celltype[[color]] %>%
    mutate(color = color,
           log_qvalue = -log10(qvalue))
}

df_results_celltype <- bind_rows(results_celltype)

df_results_celltype_top <- df_results_celltype %>%
  group_by(color) %>%
  arrange(qvalue) %>%
  slice_head(n = 5) %>%
  ungroup()


df_results_celltype_top$Description_color <- interaction(df_results_celltype_top$Description, df_results_celltype_top$color)
df_results_celltype_top$Description_color <- factor(df_results_celltype_top$Description_color, levels = rev(unique(df_results_celltype_top$Description_color)))

custom_label <- setNames(gsub("\\..*$", "", df_results_celltype_top$Description_color), df_results_celltype_top$Description_color)

p <- ggplot(df_results_celltype_top, aes(x = Description_color, y = log_qvalue)) +
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

ggsave("99_Fig/sup_fig4/Entropy_enrichment_celltype.png", plot = p , width = 10, height = 6)
ggsave("99_Fig/sup_fig4/Entropy_enrichment_celltype.pdf", plot = p , width = 10, height = 6)

p <- ggplot(df_results_celltype_top[df_results_celltype_top$color == "brown",], aes(x = Description_color, y = log_qvalue)) +
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

ggsave("99_Fig/sup_fig4/Entropy_enrichment_celltype_brown.png", plot = p , width = 4, height = 1)
ggsave("99_Fig/sup_fig4/Entropy_enrichment_celltype_brown.png", plot = p , width = 4, height = 1)


sessionInfo()
# R version 4.3.3 (2024-02-29 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
# 
# 
# locale:
# [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8    LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                    LC_TIME=Japanese_Japan.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] clusterProfiler_4.10.1      org.Hs.eg.db_3.18.0         AnnotationDbi_1.64.1        msigdbr_24.1.0              WGCNA_1.73                  fastcluster_1.2.6           dynamicTreeCut_1.63-1       DESeq2_1.42.1               SummarizedExperiment_1.32.0
# [10] Biobase_2.62.0              MatrixGenerics_1.14.0       matrixStats_1.5.0           GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.48.1         edgeR_4.0.16               
# [19] limma_3.58.1                lubridate_1.9.4             forcats_1.0.0               stringr_1.5.1               purrr_1.0.4                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.1              
# [28] tidyverse_2.0.0             dunn.test_1.3.6             dplyr_1.1.4                
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      jsonlite_2.0.0          rstudioapi_0.17.1       magrittr_2.0.3          farver_2.1.2            rmarkdown_2.29          fs_1.6.5                zlibbioc_1.48.2         vctrs_0.6.5             memoise_2.0.1           RCurl_1.98-1.17        
#  [12] ggtree_3.10.1           base64enc_0.1-3         htmltools_0.5.8.1       S4Arrays_1.2.1          curl_6.2.2              gridGraphics_0.5-1      SparseArray_1.2.4       Formula_1.2-5           htmlwidgets_1.6.4       plyr_1.8.9              impute_1.76.0          
#  [23] cachem_1.1.0            igraph_2.1.4            lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3         gson_0.1.0              Matrix_1.6-5            R6_2.6.1                fastmap_1.2.0           GenomeInfoDbData_1.2.11 aplot_0.2.5            
#  [34] digest_0.6.37           enrichplot_1.22.0       colorspace_2.1-1        patchwork_1.3.0         Hmisc_5.2-3             RSQLite_2.3.9           timechange_0.3.0        polyclip_1.10-7         httr_1.4.7              abind_1.4-8             compiler_4.3.3         
#  [45] bit64_4.6.0-1           withr_3.0.2             doParallel_1.0.17       htmlTable_2.4.3         backports_1.5.0         BiocParallel_1.36.0     viridis_0.6.5           DBI_1.2.3               ggforce_0.4.2           MASS_7.3-60.0.1         DelayedArray_0.28.0    
#  [56] HDO.db_0.99.1           tools_4.3.3             foreign_0.8-86          scatterpie_0.2.4        ape_5.8-1               nnet_7.3-19             glue_1.8.0              nlme_3.1-164            GOSemSim_2.28.1         shadowtext_0.1.4        grid_4.3.3             
#  [67] checkmate_2.3.2         cluster_2.1.6           reshape2_1.4.4          fgsea_1.28.0            generics_0.1.3          gtable_0.3.6            tzdb_0.5.0              preprocessCore_1.64.0   data.table_1.17.0       hms_1.1.3               tidygraph_1.3.1        
#  [78] XVector_0.42.0          ggrepel_0.9.6           foreach_1.5.2           pillar_1.10.1           yulab.utils_0.2.0       babelgene_22.9          splines_4.3.3           tweenr_2.0.3            treeio_1.26.0           lattice_0.22-5          survival_3.5-8         
#  [89] bit_4.6.0               tidyselect_1.2.1        GO.db_3.18.0            locfit_1.5-9.12         Biostrings_2.70.3       knitr_1.50              gridExtra_2.3           xfun_0.51               graphlayouts_1.2.2      statmod_1.5.0           stringi_1.8.7          
# [100] lazyeval_0.2.2          ggfun_0.1.8             evaluate_1.0.3          codetools_0.2-19        ggraph_2.2.1            qvalue_2.34.0           ggplotify_0.1.2         cli_3.6.4               rpart_4.1.23            munsell_0.5.1           Rcpp_1.0.14            
# [111] png_0.1-8               parallel_4.3.3          assertthat_0.2.1        blob_1.2.4              DOSE_3.28.2             bitops_1.0-9            tidytree_0.4.6          viridisLite_0.4.2       scales_1.3.0            crayon_1.5.3            rlang_1.1.5            
# [122] cowplot_1.1.3           fastmatch_1.1-6         KEGGREST_1.42.0












