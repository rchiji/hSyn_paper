library(compositions)
library(dplyr)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(glmnet)

df <- read.delim("00_src/annotations_full.txt", sep = "\t", row.names = 1)
data <- read.delim("01_formatted/gene_read_count_preprocessed.txt", sep = "\t", header = TRUE, row.names = 1)


# Reformatting
donor_keep <- intersect(rownames(df), colnames(data))

df_filtered <- df[donor_keep,]
data_filtered <- data[,donor_keep]


# CLR
res_clr <- clr(acomp(df_filtered[,6:17])) 
df_clr <- as.data.frame(res_clr)


# Preprocess
Min.CPM = 1
nSamples = 79
data_filtered <- round(data_filtered,0)
data_filtered <- data_filtered[which(apply(cpm(DGEList(counts = data_filtered)), 1, function(y) {sum(y>=Min.CPM)}) >= nSamples/2), ]

data_norm <- data_filtered
tmp <- rep("A",dim(data_norm)[2]); tmp[1] <- "B"
colData <- cbind(colnames(data_norm), tmp); colnames(colData) = c("sample", "groups")
dds <- DESeqDataSetFromMatrix(countData = data_norm, colData = colData, design = ~ groups)

data_norm <- vst(dds)
data_norm <- assay(data_norm)


# Correlation
cor_mat <- cor(df_clr, t(data_norm), method = "spearman")

p_mat <- matrix(NA,
                nrow = nrow(cor_mat),
                ncol = ncol(cor_mat),
                dimnames = dimnames(cor_mat))

for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    p_mat[i, j] <- cor.test(df_clr[, i],
                            data_norm[j, ],
                            method = "spearman",
                            exact = FALSE)$p.value
  }
}

padj_mat <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
                   nrow = nrow(p_mat),
                   ncol = ncol(p_mat),
                   dimnames = dimnames(p_mat))

sig_genes <- apply(padj_mat < 0.05, 2, any)
cor_sig <- cor_mat[, sig_genes]
padj_sig <- padj_mat[, sig_genes]
pheatmap(cor_sig)

sig_genes_2 <- apply(cor_sig > 0.5, 2, any)
cor_sig_2 <- cor_sig[, sig_genes_2]
padj_sig_2 <- padj_sig[, sig_genes_2]
pheatmap(cor_sig_2)

colnames(cor_sig_2) <- sub("^[^_]+_", "", colnames(cor_sig_2))
custom_labels <- c(
  "Fibrous_tissue__dense_irregular" = "Fibrous tissue (dense, irregular)",
  "Fibrous_tissue__dense_regular" = "Fibrous tissue (dense, regular)",
  "Fibrous_tissue__loose" = "Fibrous tissue (loose)",
  "Micro_vessel" = "Micro vessel",
  "Large_vessel" = "Large vessel"
)
rownames(cor_sig_2) <- ifelse(rownames(cor_sig_2) %in% names(custom_labels),
                        custom_labels[rownames(cor_sig_2)], rownames(cor_sig_2))
colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

pdf("99_Fig/fig2/Heatmap_cor_tissue_gene.pdf", width = 12, height = 3)
pheatmap(cor_sig_2, 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_method = "ward.D2"
)
dev.off()


# Gene number
sig_mat <- (padj_mat < 0.05) & (cor_mat > 0.5)

n_sig_genes <- rowSums(sig_mat)

df_gene_count <- data.frame(
  tissue = names(n_sig_genes),
  n_genes = as.numeric(n_sig_genes)
)

df_gene_count$tissue <- factor(df_gene_count$tissue, levels = df_gene_count$tissue[order(df_gene_count$n_genes)])

p1 <- ggplot(df_gene_count, aes(x = tissue, y = n_genes)) +
  geom_col(width = 0.9, linewidth = 0.1, fill = "white", color = "black") +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  scale_x_discrete(labels = custom_labels) +
  labs(y = "Number of correlated genes") +
  coord_flip() +
  theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 0.1, fill = FALSE),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_blank(),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 6)
  )
ggsave("99_Fig/fig2/Barplot_significant_gene_counts.png", plot = p1 , width = 3, height = 2)
ggsave("99_Fig/fig2/Barplot_significant_gene_counts.pdf", plot = p1 , width = 3, height = 2)


# Elastic net regression
set.seed(123)
X <- t(data_norm)

Prediction_CrossValidation <- function(y, X, nfold = 5, alpha = 0.5) {
  len_y <- length(y)
  fold_id <- sample(rep(1:nfold, length.out = len_y))
  pred <- rep(NA_real_, len_y)
  for (n in 1:nfold) {
    test_idx  <- which(fold_id == n)
    train_idx <- setdiff(seq_len(len_y), test_idx)
    cvfit <- cv.glmnet(
      x = X[train_idx, , drop = FALSE],
      y = y[train_idx],
      family = "gaussian",
      alpha = alpha,
      nfolds = 5,
      standardize = TRUE
    )
    pred[test_idx] <- as.numeric(
      predict(cvfit, newx = X[test_idx, , drop = FALSE], s = "lambda.min")
    )
  }
  data.frame(
    observed = y,
    predicted = pred
  )
}

res_list <- lapply(colnames(df_clr), function(x) {
  out <- Prediction_CrossValidation(y = df_clr[, x], X = X)
  out$tissue <- x
  out
})
names(res_list) <- colnames(df_clr)
res <- do.call(rbind, res_list)

df_res <- do.call(rbind, lapply(names(res_list), function(x) {
  tmp <- res_list[[x]]
  data.frame(
    tissue = x,
    spearman_r = cor(tmp$observed, tmp$predicted, method = "spearman"),
    r2 = 1 - sum((tmp$observed - tmp$predicted)^2) / sum((tmp$observed - mean(tmp$observed))^2),
    rmse = sqrt(mean((tmp$observed - tmp$predicted)^2))
  )
}))

df_res$tissue <- factor(df_res$tissue, levels = df_res$tissue[order(df_res$spearman_r)])

p2 <- ggplot(df_res, aes(x = tissue, y = spearman_r)) +
  geom_col(width = 0.9, linewidth = 0.1, fill = "white", color = "black") +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  scale_y_continuous(limits = c(-0.4, 1), breaks = seq(-0.4, 1, 0.2)) +
  scale_x_discrete(labels = custom_labels) +
  labs(y = "Spearman correlation") +
  coord_flip() +
  theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 0.1, fill = FALSE),
    axis.ticks.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_line(linewidth = 0.1),
    axis.line = element_blank(),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 6)
    )
ggsave("99_Fig/fig2/Barplot_prediction_performance.png", plot = p2 , width = 3, height = 2)
ggsave("99_Fig/fig2/Barplot_prediction_performance.pdf", plot = p2 , width = 3, height = 2)

