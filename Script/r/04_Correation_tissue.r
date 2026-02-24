library(compositions)
library(psych)
library(pheatmap)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)

custom_labels <- c(
  "Fibrous_tissue__dense_irregular" = "Fibrous tissue\n(dense.irregular)",
  "Fibrous_tissue__dense_regular" = "Fibrous tissue\n(dense.regular)",
  "Fibrous_tissue__loose" = "Fibrous tissue\n(loose)",
  "Micro_vessel" = "Micro vessel",
  "Large_vessel" = "Large vessel"
)


res <- df[,6:17]
res_clr <- clr(acomp(res))
res_clr_df <- as.data.frame(res_clr)

cor_tissue <- corr.test((res_clr_df), method = "spearman")
hoge_tissue <- cor_tissue
tmp_tissue <- matrix(data = NA, nrow = nrow(hoge_tissue$p), ncol = ncol(hoge_tissue$p))
tmp_tissue[upper.tri(tmp_tissue)] <- hoge_tissue$p.adj
p <- t(tmp_tissue)
p[upper.tri(p)] <- hoge_tissue$p.adj
p_num <- p
for(i in 1:nrow(p)){
    for(j in 1:ncol(p)){
        tmp <- as.numeric(p[i,j])
        if(is.na(tmp)){
            symbol <- ""
        } else if(tmp < 0.001){
            symbol <- "***"
        } else if(tmp < 0.01){
            symbol <- "**"
        } else if(tmp < 0.05){
            symbol <- "*"
        } else{
            symbol <- ""
        }
        p[i,j] <- symbol
    }
}

mat <- hoge_tissue$r

rownames(mat) <- ifelse(rownames(mat) %in% names(custom_labels),
                        custom_labels[rownames(mat)], rownames(mat))
colnames(mat) <- ifelse(colnames(mat) %in% names(custom_labels),
                        custom_labels[colnames(mat)], colnames(mat))

colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

png("99_Fig/fig2/corr_tissue_clr.png", width = 3.25, height = 3.25, units = "in", res = 300)
plot_tissue  <- pheatmap(mat,
                         display_numbers = p,
                         angle_col = 90,
                         fontsize = 5,
                         cellwidth = 10,
                         cellheight = 10,
                         color = colors,
                         treeheight_row = 20,
                         treeheight_col = 20,
                         )
dev.off()
pdf("99_Fig/fig2/corr_tissue_clr.pdf", width = 3.25, height = 3.25)
plot_tissue  <- pheatmap(mat,
                         display_numbers = p,
                         angle_col = 90,
                         fontsize = 5,
                         cellwidth = 10,
                         cellheight = 10,
                         color = colors,
                         treeheight_row = 20,
                         treeheight_col = 20,
)
dev.off()


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
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] pheatmap_1.0.12    psych_2.5.3        compositions_2.0-8
# 
# loaded via a namespace (and not attached):
#  [1] tensorA_0.36.2.1    R6_2.6.1            bayesm_3.1-6        RColorBrewer_1.1-3  DEoptimR_1.1-3-1    lattice_0.22-5      gtable_0.3.6        glue_1.8.0          parallel_4.3.3      lifecycle_1.0.4     cli_3.6.4           scales_1.3.0        grid_4.3.3         
# [14] robustbase_0.99-4-1 mnormt_2.1.1        compiler_4.3.3      rstudioapi_0.17.1   tools_4.3.3         nlme_3.1-164        munsell_0.5.1       colorspace_2.1-1    Rcpp_1.0.14         rlang_1.1.5         MASS_7.3-60.0.1