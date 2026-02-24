library(dplyr)
library(DESeq2)
library(edgeR)

df <- read.delim(file = "01_formatted/annotations_full_tissue_conection_cluster.txt", sep = "\t", row.names = 1)
df_OA <- df[df$Diagnosis == "OA",]

data <- read.delim("00_src/gene_read_count.txt", sep = "\t", header = TRUE, row.names = 1)
# Gene ID duplication
# Gene IDs such as “PAR_Y_XX” appear to denote homologous genes located between the X and Y chromosomes (i.e., genes in the pseudoautosomal region).
# We identified 44 such duplicated entries; however, for all of them, the count values were identical between the homologous gene pairs. Therefore, we removed the PAR_Y entries in Excel.

data <- data[, colnames(data) != "D070"]
# The entries were removed due to the absence of histological data.

data <- data[!is.na(data$Gene.type) & data$Gene.type == "protein_coding", ]
data <- data[,1:86]
data <- na.omit(data)

data$id <- rownames(data)
data <- data %>% 
  mutate(gene = paste(data$id, data$gene_name, sep = "_"))
row.names(data) <- data$gene
data <- data[,2:86]

write.table(data, "01_formatted/gene_read_count_preprocessed.txt", sep = "\t", row.names = T, col.names = NA)


donor_OA <- rownames(df_OA)
donor_select <- donor_OA[donor_OA %in% colnames(data)]
data <- data[, donor_select]

Min.CPM = 1
nSamples = 72
data_filtered <- round(data,0)
data_filtered <- data_filtered[which(apply(cpm(DGEList(counts = data_filtered)), 1, function(y) {sum(y>=Min.CPM)}) >= nSamples/2), ] # 13678 genes

data_norm <- data_filtered
tmp <- rep("A",dim(data_norm)[2]); tmp[1] <- "B"
colData <- cbind(colnames(data_norm), tmp); colnames(colData) = c("sample", "groups")
dds <- DESeqDataSetFromMatrix(countData = data_norm, colData = colData, design = ~ groups)

data_norm <- vst(dds)
data_norm <- assay(data_norm) # 13694 genes

write.table(data_norm, "01_formatted/gene_read_count_normalized.txt", sep = "\t", row.names = T, col.names = NA)


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
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] edgeR_4.0.16                limma_3.58.1                dplyr_1.1.4                 DESeq2_1.42.1               SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0       matrixStats_1.5.0           GenomicRanges_1.54.1       
# [10] GenomeInfoDb_1.38.8         IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.48.1        
# 
# loaded via a namespace (and not attached):
#  [1] Matrix_1.6-5            gtable_0.3.6            compiler_4.3.3          crayon_1.5.3            tidyselect_1.2.1        Rcpp_1.0.14             bitops_1.0-9            parallel_4.3.3          scales_1.3.0            statmod_1.5.0           BiocParallel_1.36.0    
# [12] lattice_0.22-5          ggplot2_3.5.1           R6_2.6.1                XVector_0.42.0          generics_0.1.3          S4Arrays_1.2.1          tibble_3.2.1            DelayedArray_0.28.0     munsell_0.5.1           GenomeInfoDbData_1.2.11 pillar_1.10.1          
# [23] rlang_1.1.5             SparseArray_1.2.4       cli_3.6.4               magrittr_2.0.3          zlibbioc_1.48.2         locfit_1.5-9.12         grid_4.3.3              rstudioapi_0.17.1       lifecycle_1.0.4         vctrs_0.6.5             glue_1.8.0             
# [34] codetools_0.2-19        abind_1.4-8             RCurl_1.98-1.17         colorspace_2.1-1        pkgconfig_2.0.3         tools_4.3.3    

