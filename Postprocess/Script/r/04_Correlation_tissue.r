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

max_val <- max(abs(mat), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 101)

colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

png("99_Fig/sup_fig3/corr_tissue_clr.png", width = 3.25, height = 3.25, units = "in", res = 300)
plot_tissue  <- pheatmap(mat,
                         display_numbers = p,
                         angle_col = 90,
                         fontsize = 5,
                         cellwidth = 10,
                         cellheight = 10,
                         color = colors,
                         breaks = breaks,
                         treeheight_row = 20,
                         treeheight_col = 20,
                         )
dev.off()
pdf("99_Fig/sup_fig3/corr_tissue_clr.pdf", width = 3.25, height = 3.25)
plot_tissue  <- pheatmap(mat,
                         display_numbers = p,
                         angle_col = 90,
                         fontsize = 5,
                         cellwidth = 10,
                         cellheight = 10,
                         color = colors,
                         breaks = breaks,
                         treeheight_row = 20,
                         treeheight_col = 20,
)
dev.off()

