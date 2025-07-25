library(dplyr)
library(pheatmap)
library(ggplot2)
library(patchwork)
library(ggsci)

# ~~~~~~~~~~~~~~~~~~~ ----
## *****************************************************************************
## データ読み込み ----
# 241205
# QuPathから書き出したObject情報を集計
## *****************************************************************************

file_paths <- list.files(path = "Annotation_ratio_241204/", full.names = TRUE)
samplenames <- sapply(file_paths, function(x) {
  x <- basename(x)
  gsub(pattern = "_HE.*",replacement = "",x)
  })
file_list <- lapply(file_paths, read.delim)
names(file_list) <- samplenames

levels <- c("adipose", "Fibro(dense,irregular)", "Fibro(dense,regular)", 
"Fibro(loose)", "Immune cells", "lining", "muscle", "plasma", 
"RBC", "Stroma", "vessel", "vessel(large)")

## * classificationごとに面積を集計 ----
tmp <- file_list[[1]]
tmp$Classification <- factor(tmp$Classification, levels = levels)
tmp %>% 
  group_by(Classification) %>% 
  summarise(
    Area_sum = sum(Area.µm.2),
    count = n()
  )
  # # A tibble: 12 × 3
  #    Classification          Area_sum count
  #    <fct>                      <dbl> <int>
  #  1 adipose                 4522880.   186
  #  2 Fibro(dense,irregular) 30133933.  1637
  #  3 Fibro(dense,regular)     685179.    76
  #  4 Fibro(loose)           11125381.  1349
  #  5 Immune cells            7034939.   888
  #  6 lining                  4755151.  1651
  #  7 muscle                      462.    21
  #  8 plasma                    30800.   209
  #  9 RBC                      398791.  1420
  # 10 Stroma                     9460.    28
  # 11 vessel                  1145526.  2382
  # 12 vessel(large)            293762.    10


## * 面積比率計算 ----

summary_list <- lapply(file_list, function(x){
  x$Classification <- factor(x$Classification, levels = levels)
  x <- x %>% 
    group_by(Classification) %>% 
    summarise(
      Area_sum = sum(Area.µm.2),
      count = n()
    )
  return(x)
  })


tmp <- summary_list[[1]]
tmp$Area_ratio <- tmp$Area_sum / sum(tmp$Area_sum)
tmp$average_area <- tmp$Area_sum / tmp$count
tmp
  # # A tibble: 12 × 5
  #    Classification          Area_sum count Area_ratio average_area
  #    <chr>                      <dbl> <int>      <dbl>        <dbl>
  #  1 Fibro(dense,irregular) 30133933.  1637 0.501           18408. 
  #  2 Fibro(dense,regular)     685179.    76 0.0114           9016. 
  #  3 Fibro(loose)           11125381.  1349 0.185            8247. 
  #  4 Immune cells            7034939.   888 0.117            7922. 
  #  5 RBC                      398791.  1420 0.00663           281. 
  #  6 Stroma                     9460.    28 0.000157          338. 
  #  7 adipose                 4522880.   186 0.0752          24317. 
  #  8 lining                  4755151.  1651 0.0791           2880. 
  #  9 muscle                      462.    21 0.00000768         22.0
  # 10 plasma                    30800.   209 0.000512          147. 
  # 11 vessel                  1145526.  2382 0.0190            481. 
  # 12 vessel(large)            293762.    10 0.00488         29376.


summary_list <- lapply(summary_list, function(x){
  x$Area_ratio <- x$Area_sum / sum(x$Area_sum)
  x$average_area <- x$Area_sum / x$count
  return(x)
  })


value <- tmp$Area_ratio
names(value) <- tmp$Classification
value
          #      adipose Fibro(dense,irregular)   Fibro(dense,regular)           Fibro(loose) 
          # 7.521053e-02           5.010942e-01           1.139377e-02           1.850029e-01 
          # Immune cells                 lining                 muscle                 plasma 
          # 1.169833e-01           7.907294e-02           7.677537e-06           5.121726e-04 
          #          RBC                 Stroma                 vessel          vessel(large) 
          # 6.631453e-03           1.573087e-04           1.904884e-02           4.884934e-03


res <- do.call(rbind, lapply(summary_list,function(x){
  value <- x$Area_ratio
  names(value) <- x$Classification
  value <- value[levels]
  value[is.na(value)] <- 0
  return(value)
}))

pdf(file = "plots_R/002_Tissue_ratio_pheatmap.pdf", width = 15, height = 4)
pheatmap(t(res))
dev.off()

pdf(file = "plots_R/003_Tissue_ratio_pheatmap.pdf", width = 15, height = 4)
pheatmap(t(res), scale = "column")
dev.off()

pdf(file = "plots_R/004_Tissue_ratio_pheatmap.pdf", width = 15, height = 4)
pheatmap(t(res), scale = "row")
dev.off()

# ~~~~~~~~~~~~~~~~~~~ ----
## *****************************************************************************
## 患者メタデータ ----
## *****************************************************************************

library(readxl)
annotations_metadata_all_patients <- read_excel(
  "~/東京大学 整形外科/ヒト滑膜/annotations_metadata_all patients.xlsx"
  )
anno <- data.frame(annotations_metadata_all_patients)
rownames(anno) <- anno[,1]
rownames(anno) <- gsub(pattern = "R", replacement = "_R", rownames(anno))
rownames(anno) <- gsub(pattern = "L", replacement = "_L", rownames(anno))
anno <- anno[rownames(res), c("Diagnosis","Joint","Age","Sex","BMI")]



pdf(file = "plots_R/002_Tissue_ratio_pheatmap.pdf", width = 15, height = 8)
pheatmap(t(res), annotation_col = anno)
dev.off()

pdf(file = "plots_R/003_Tissue_ratio_pheatmap.pdf", width = 15, height = 8)
pheatmap(t(res), scale = "column", annotation_col = anno)
dev.off()

pdf(file = "plots_R/004_Tissue_ratio_pheatmap.pdf", width = 15, height = 8)
pheatmap(t(res), scale = "row", annotation_col = anno)
dev.off()



# ~~~~~~~~~~~~~~~~~~~ ----
## *****************************************************************************
## β diversity ----
## *****************************************************************************

library(vegan)
bc_dist <- vegdist( res, method = "bray")

# MDS
mds_beta <- cmdscale(bc_dist, k = 2, eig = TRUE)
mdsData <- as.data.frame(mds_beta$points[,1:2])

mds_df <- cbind(mdsData, res, anno)

MDSplot <- 
  ggplot(mds_df, aes(x=V1,y=V2)) +
  geom_point() +
  theme_bw() + 
  theme(aspect.ratio = 1)

MDSplot <- ggplot(mds_df, aes(x=V1,y=V2)) +
  geom_point(shape=21, aes(fill = Diagnosis,
                           size = Age, 
                           colour = Sex)) +
  scale_color_manual(values = c("black","gray")) +
  # scale_fill_nejm() +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  xlab("MDS1") +
  ylab("MDS2")
  
# fit MDS model
env_fit <- envfit(mds_beta, res)

# correlation between dimension value and tissue ratio
cor_df <- data.frame(cor(res,mds_df[,c(1,2)]))
# adjust for plot scale
cor_df <- cor_df * ordiArrowMul(env_fit)
cor_df$distance <- sqrt(cor_df$V1 ** 2 + cor_df$V2 ** 2)

MDScorplot <- 
  ggplot(data = cor_df, aes(x = 0, y = 0, xend = V1, yend = V2, colour = distance)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm"))) +
  scale_color_gradient2(low = "gray",mid = "orange", high = "red") +
  ggrepel::geom_text_repel(data = cor_df, aes(x = V1, y = V2, label = rownames(cor_df)),
            color = "black", hjust = -0.2, vjust = -0.2,
            max.overlaps = Inf) +
  theme_bw() + 
  theme(aspect.ratio = 1)

dir.create("plots_R")
pdf("plots_R/001_betaDiversity_MDSplot.pdf", width = 12,height = 6)
MDSplot + MDScorplot
dev.off()



## *****************************************************************************
## CLR ----
## *****************************************************************************

library(compositions) # install.packages("compositions")
library(factoextra) # install.packages("factoextra")

res_clr <- clr(res)
pheatmap(res_clr)
pheatmap(t(res_clr))

d <- dist(res_clr)

# MDS
mds <- metaMDS(d)

mds_df2 <- as.data.frame(mds$points)
mds_df2 <- cbind(mds_df2,res)

MDSplot2 <- 
  ggplot(mds_df2, aes(x = MDS1, y = MDS2)) +
  geom_point() +
  theme_bw() + 
  theme(aspect.ratio = 1)


# fit MDS model
env_fit2 <- envfit(mds, res)

# correlation between dimension value and tissue ratio
cor_df2 <- data.frame(cor(res,mds_df2[,c(1,2)]))
# adjust for plot scale
cor_df2 <- cor_df2 * ordiArrowMul(env_fit2)
cor_df2$distance <- sqrt(cor_df2$MDS1 ** 2 + cor_df2$MDS2 ** 2)

MDScorplot2 <- 
  ggplot(data = cor_df2, aes(x = 0, y = 0, xend = MDS1, yend = MDS2, colour = distance)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm"))) +
  scale_color_gradient2(low = "gray",mid = "orange", high = "red") +
  ggrepel::geom_text_repel(data = cor_df2, aes(x = MDS1, y = MDS2, label = rownames(cor_df2)),
            color = "black", hjust = -0.2, vjust = -0.2,
            max.overlaps = Inf) +
  theme_bw() + 
  theme(aspect.ratio = 1)

MDSplot2 + MDScorplot2


## * CLRで組織間相関
res_clr2 <- as.matrix(res_clr)
res_clr_cor <- cor(res_clr2)
pheatmap(res_clr_cor)

## * CLRで患者間相関
res_clr2 <- as.matrix(res_clr)
res_clr_cor <- cor(t(res_clr2))
pheatmap(res_clr_cor, cutree_cols = 4, cutree_rows = 4,
         clustering_method = "ward.D2")
pheatmap(res_clr_cor, cutree_cols = 4, cutree_rows = 4,
         clustering_method = "average")


## * 比率で組織間相関
res_clr_cor <- cor(res)
diag(res_clr_cor) <- NA
pdf(file = "plots_R/005_Tissue_clr_pheatmap.pdf", width = 5, height = 5)
pheatmap(res_clr_cor)
dev.off()



## * 比率で患者間相関
res_clr_cor <- cor(t(res))
pheatmap(res_clr_cor, cutree_cols = 4, cutree_rows = 4,
         clustering_method = "ward.D2")


# ~~~~~~~~~~~~~~~~~~~ ----
## *****************************************************************************
## 接続パターン ----
## *****************************************************************************

Degree_df1 <- read.csv("~/東京大学 整形外科/ヒト滑膜/labels_QuPath6/Degree_df_241115.csv",
                       row.names = 1)
Degree_df2 <- read.csv("~/東京大学 整形外科/ヒト滑膜/labels_QuPath6/Degree_df2_241115.csv",
                       row.names = 1)

Degree_df1 <- t(Degree_df1)
Degree_df2 <- t(Degree_df2)
pheatmap(Degree_df1, annotation_col = anno)
pheatmap(Degree_df2, annotation_col = anno)

library(NMF)

estim_r1 <- nmf(Degree_df1, seq(2, 8), nrun = 10, seed = 123)
plot(estim_r1)
## --> 4で良さそう

estim_r2 <- nmf(Degree_df2, seq(2, 8), nrun = 10, seed = 123)
plot(estim_r2)
## --> 4で良さそう


nmf1 <- nmf(Degree_df1,rank = 4,seed=123, .options = "t")

# 結果の要約
summary(nmf1)
       #        rank sparseness.basis  sparseness.coef  silhouette.coef silhouette.basis        residuals 
       #   4.0000000        0.8371494        0.5041553        0.6221112        0.8846468        9.0539520 
       #       niter              cpu          cpu.all             nrun 
       # 570.0000000               NA               NA        1.0000000 

# 結果の可視化
plot(nmf1)

basismap(nmf1, scale="r1") # --> これがデフォルト 0-1に正規化される
basismap(nmf1, scale="row") # --> Z score
basismap(nmf1, scale="none")
coefmap(nmf1)
coefmap(nmf1, scale = "c1")
coefmap(nmf1, scale="none")
consensusmap(nmf1)

# 再構成された行列
W <- basis(nmf1)
H <- coef(nmf1)
V_reconstructed <- W %*% H

pheatmap(V_reconstructed)
# --> 入力dfが復元できている

range01 <- function(x){
  x <- (x-min(x))/(max(x)-min(x))
  return(x)
  }
W_scaled <- apply(W,1,range01)
pheatmap(t(W_scaled), cluster_cols = F)

H_scaled <- apply(H,1,range01)
pheatmap(t(H_scaled), cluster_rows = F)



nmf2 <- nmf(Degree_df2,rank = 4,seed=123, .options = "t")

# 結果の可視化
plot(nmf2)

basismap(nmf2, scale="r1") 
coefmap(nmf2)
coefmap(nmf2, scale = "c1")
consensusmap(nmf2)

# 再構成された行列
W2 <- basis(nmf2)
H2 <- coef(nmf2)
V_reconstructed2 <- W2 %*% H2

pheatmap(V_reconstructed2)
# --> 入力dfが復元できている

range01 <- function(x){
  x <- (x-min(x))/(max(x)-min(x))
  return(x)
  }
W2_scaled <- apply(W2,1,range01)
pheatmap(t(W2_scaled), cluster_cols = F)

H2_scaled <- apply(H2,1,range01)
pheatmap(t(H2_scaled), cluster_rows = F)


## k=6

nmf1 <- nmf(Degree_df1,rank = 6,seed=123, .options = "t")

basismap(nmf1)
coefmap(nmf1)
consensusmap(nmf1)

# 再構成された行列
W <- basis(nmf1)
H <- coef(nmf1)
V_reconstructed <- W %*% H

pheatmap(V_reconstructed)
# --> 入力dfが復元できている

range01 <- function(x){
  x <- (x-min(x))/(max(x)-min(x))
  return(x)
  }
W_scaled <- apply(W,1,range01)

pdf(file = "plots_R/006_NMF_basis.pdf", width = 6, height = 18)
pheatmap(t(W_scaled), cluster_cols = F)
dev.off()


H_scaled <- apply(H,1,range01)
pdf(file = "plots_R/007_NMF_coef.pdf", width = 18, height = 7)
pheatmap(t(H_scaled), cluster_rows = F, annotation_col = anno)
dev.off()




nmf2 <- nmf(Degree_df2,rank = 6,seed=123, .options = "t")

# 結果の可視化
plot(nmf2)

basismap(nmf2, scale="r1") 
coefmap(nmf2)
coefmap(nmf2, scale = "c1")
consensusmap(nmf2)

# 再構成された行列
W2 <- basis(nmf2)
H2 <- coef(nmf2)
V_reconstructed2 <- W2 %*% H2

pheatmap(V_reconstructed2)
# --> 入力dfが復元できている

range01 <- function(x){
  x <- (x-min(x))/(max(x)-min(x))
  return(x)
  }
W2_scaled <- apply(W2,1,range01)
pheatmap(t(W2_scaled), cluster_cols = F)

H2_scaled <- apply(H2,1,range01)
pheatmap(t(H2_scaled), cluster_rows = F)




pattern_cluster <- apply(W_scaled,2,which.max)
cluster_names <- sort(unique(gsub(pattern = "_.*", replacement = "", names(pattern_cluster))))
mat <- data.frame(matrix(nrow = length(cluster_names), ncol = length(cluster_names)))
dimnames(mat) <- list(cluster_names,cluster_names)

for( i in cluster_names){
  for( j in cluster_names){
    pattern <- paste(i,j,sep = "_")
    cluster <- pattern_cluster[pattern]
    if( !is.na(cluster) ){
      mat[i,j] <- mat[j,i] <- cluster
      }
  }
}

pheatmap(mat)


pattern_cluster2 <- apply(W2_scaled,2,which.max)
cluster_names2 <- sort(unique(gsub(pattern = "_.*", replacement = "", names(pattern_cluster2))))
mat2 <- data.frame(matrix(nrow = length(cluster_names2), ncol = length(cluster_names2)))
dimnames(mat2) <- list(cluster_names2,cluster_names2)

for( i in cluster_names2){
  for( j in cluster_names2){
    pattern <- paste(i,j,sep = "_")
    cluster <- pattern_cluster2[pattern]
    if( !is.na(cluster) ){
      mat2[i,j] <- mat2[j,i] <- cluster
      }
  }
}

pheatmap(mat2)



d <- dist(H_scaled)

# MDS
mds <- metaMDS(d)

colnames(H_scaled) <- paste0("Feature", seq(ncol(H_scaled)))
mds_df2 <- as.data.frame(mds$points)
mds_df2 <- cbind(mds_df2, H_scaled, anno)

ggplot(mds_df2, aes(x=MDS1,y=MDS2)) +
  geom_point(shape=21, aes(fill = Diagnosis,
                           size = Age, 
                           colour = Sex)) +
  scale_color_manual(values = c("black","gray")) +
  # scale_fill_nejm() +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  xlab("MDS1") +
  ylab("MDS2")
  
ggplot(mds_df2, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(colour = Feature1, size=Feature1)) +
  # scale_fill_nejm() +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  xlab("MDS1") +
  ylab("MDS2")

ggplot(mds_df2, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(colour = Feature2, size=Feature2)) +
  # scale_fill_nejm() +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  xlab("MDS1") +
  ylab("MDS2")

ggplot(mds_df2, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(colour = Feature3, size=Feature3)) +
  # scale_fill_nejm() +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  xlab("MDS1") +
  ylab("MDS2")

ggplot(mds_df2, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(colour = Feature4, size=Feature4)) +
  # scale_fill_nejm() +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  xlab("MDS1") +
  ylab("MDS2")

ggplot(mds_df2, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(colour = Feature5, size=Feature5)) +
  # scale_fill_nejm() +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  xlab("MDS1") +
  ylab("MDS2")

ggplot(mds_df2, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(colour = Feature6, size=Feature6)) +
  # scale_fill_nejm() +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  xlab("MDS1") +
  ylab("MDS2")


## *****************************************************************************
## * β diversity ----
## *****************************************************************************

library(vegan)
bc_dist <- vegdist( t(Degree_df1), method = "bray")

# MDS
mds_beta <- cmdscale(bc_dist, k = 2, eig = TRUE)
mdsData <- as.data.frame(mds_beta$points[,1:2])

mds_df <- cbind(mdsData, t(Degree_df1), anno)

MDSplot <- 
  ggplot(mds_df, aes(x=V1,y=V2)) +
  geom_point() +
  theme_bw() + 
  theme(aspect.ratio = 1)

MDSplot <- ggplot(mds_df, aes(x=V1,y=V2)) +
  geom_point(shape=21, aes(fill = Diagnosis,
                           size = Age, 
                           colour = Sex)) +
  scale_color_manual(values = c("black","gray")) +
  # scale_fill_nejm() +
  theme_bw() + 
  theme(aspect.ratio = 1) +
  xlab("MDS1") +
  ylab("MDS2")

MDSplot 


pheatmap(Degree_df2)
pheatmap(Degree_df2, scale = "row")
pheatmap(Degree_df2, scale = "column")

pheatmap(Degree_df2[rowSums(Degree_df2) > 0.1,])
pheatmap(Degree_df1[rowSums(Degree_df1) > 0.1,], annotation_col = anno)
pheatmap(Degree_df1[rowSums(Degree_df1) > 0.1,], annotation_col = anno,
         scale = "row", breaks = seq(-1.5,1.5,length.out=100),
         clustering_method = "ward.D2",
         cutree_cols = 3)


pdf(file = "plots_R/008_pattern_hierarchical.pdf", width = 10, height = 8)
pheatmap(Degree_df1[rowSums(Degree_df1) > 0.1,], annotation_col = anno,
         scale = "row", breaks = seq(-1.5,1.5,length.out=100),
         clustering_method = "ward.D2",
         cutree_cols = 3)
dev.off()

anno_row <- data.frame(
  A = gsub(pattern = "_.*", replacement = "", rownames(Degree_df1)),
  B = gsub(pattern = ".*_", replacement = "", rownames(Degree_df1)),
  row.names = rownames(Degree_df1)
  )
anno_row_color <- viridis::turbo(n = length(cluster_names))
names(anno_row_color) <- cluster_names

pdf(file = "plots_R/009_pattern_hierarchical.pdf", width = 10, height = 8)
pheatmap(Degree_df1[rowSums(Degree_df1) > 0.1,],
         annotation_col = anno,
         annotation_row = anno_row,
         annotation_colors = list(A = anno_row_color,
                                  B = anno_row_color),
         scale = "row", breaks = seq(-1.5,1.5,length.out=100),
         clustering_method = "ward.D2",
         cutree_cols = 3)
dev.off()

p <- pheatmap(Degree_df1[rowSums(Degree_df1) > 0.1,],
         annotation_col = anno,
         annotation_row = anno_row,
         annotation_colors = list(A = anno_row_color,
                                  B = anno_row_color),
         scale = "row", breaks = seq(-1.5,1.5,length.out=100),
         clustering_method = "ward.D2",
         cutree_cols = 3)


p$tree_col$labels[p$tree_col$order]
  #   [1] "D095"   "D011"   "D013"   "D005"   "D094"   "D021"   "D040"   "D036"   "D127"   "D129"  
  #  [11] "D046"   "D082"   "D071"   "D109"   "D117"   "D068"   "D092"   "D016"   "D054"   "D076"  
  #  [21] "D125"   "D018"   "D075"   "D055"   "D122"   "D093_L" "D020"   "D066"   "D118"   "D009"  
  #  [31] "D101"   "D025"   "D111"   "D010"   "D007"   "D060"   "D093_R" "D103"   "D026"   "D030"  
  #  [41] "D006"   "D024"   "D019"   "D096"   "D083"   "D084"   "D033"   "D073"   "D022"   "D085"  
  #  [51] "D014"   "D090"   "D034"   "D097"   "D112"   "D027"   "D087"   "D023"   "D002"   "D058"  
  #  [61] "D004"   "D038"   "D074_L" "D061"   "D050"   "D106"   "D031"   "D045"   "D044"   "D048"  
  #  [71] "D064"   "D028"   "D047"   "D063"   "D049"   "D108"   "D053"   "D119"   "D081"   "D080"  
  #  [81] "D113"   "D001"   "D035"   "D056"   "D098"   "D104"   "D123"   "D037"   "D079"   "D110"  
  #  [91] "D078"   "D042"   "D051"   "D077"   "D100"   "D065"   "D128"   "D067"   "D130"   "D029"  
  # [101] "D039"   "D099"   "D105"   "D041"   "D086"   "D102"   "D057"   "D089"   "D003"   "D072"  
  # [111] "D017"   "D059"   "D091"   "D015"   "D074_R" "D107"   "D069"   "D126" 



# ~~~~~~~~~~~~~~~~~~~ ----
## *****************************************************************************
## GNN Kmeans cluster ----
## *****************************************************************************


file_paths <- list.files(path = "../labels_QuPath6/Kmeans_ratio_241207/", full.names = TRUE)
samplenames <- sapply(file_paths, function(x) {
  x <- basename(x)
  gsub(pattern = "_HE.*",replacement = "",x)
  })
file_list <- lapply(file_paths, read.delim)
names(file_list) <- samplenames

file_list <- lapply(file_list, function(x){
  x[x$Object.type=="Tile",]
})

## * Kmeansクラスター数を集計 ----
tmp <- file_list[[1]]
tmp %>% 
  group_by(Classification) %>% 
  summarise(
    count = n()
  )



## * クラスター比率計算 ----

summary_list <- lapply(file_list, function(x){
  x <- x %>% 
    group_by(Classification) %>% 
    summarise(
      count = n()
    )
  return(x)
  })


summary_list <- lapply(summary_list, function(x){
  x$cluster_ratio <- x$count / sum(x$count)
  return(x)
  })


levels <- sort(unique(unlist(sapply(summary_list, function(x) unique(x$Classification)))))

res <- do.call(rbind, lapply(summary_list,function(x){
  value <- x$cluster_ratio
  names(value) <- x$Classification
  value <- value[levels]
  value[is.na(value)] <- 0
  return(value)
}))

pdf(file = "plots_R/010_KmeanCluster_ratio_pheatmap.pdf", width = 18, height = 4)
pheatmap(t(res))
dev.off()

pdf(file = "plots_R/011_KmeanCluster_ratio_pheatmap.pdf", width = 18, height = 4)
pheatmap(t(res), scale = "column")
dev.off()

pdf(file = "plots_R/012_KmeanCluster_ratio_pheatmap.pdf", width = 18, height = 4)
pheatmap(t(res), scale = "row")
dev.off()


# ~~~~~~~~~~~~~~~~~~~ ----
## *****************************************************************************
## Lining feature ----
## *****************************************************************************


file_paths <- list.files(path = "../labels_QuPath6/Kmeans_ratio_241207/", full.names = TRUE)
samplenames <- sapply(file_paths, function(x) {
  x <- basename(x)
  gsub(pattern = "_HE.*",replacement = "",x)
  })
file_list <- lapply(file_paths, data.table::fread)
names(file_list) <- samplenames

file_list <- lapply(file_list, function(x){
  x <- data.frame(x,row.names=1)
  x <- x[x$Object.type=="Tile",]
  x <- x[ x$Name == "lining", ]
  return(x)
})

x <- file_list[[1]]
x <- data.frame(x,row.names=1)
x <- x[x$Object.type=="Tile",]
x <- x[ x$Name == "lining", ] 
dput(colnames(x))

feature_columns <- c("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", 
"X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", 
"X20", "X21", "X22", "X23", "X24", "X25", "X26", "X27", "X28", 
"X29", "X30", "X31", "X32", "X33", "X34", "X35", "X36", "X37", 
"X38", "X39", "X40", "X41", "X42", "X43", "X44", "X45", "X46", 
"X47", "X48", "X49", "X50", "X51", "X52", "X53", "X54", "X55", 
"X56", "X57", "X58", "X59", "X60", "X61", "X62", "X63")

hist(x$X0)

df <- x
# ここでは各列の密度を計算
df_melted <- reshape2::melt(df, id.vars = NULL)

# ggplotを使った密度プロット
ggplot(df_melted, aes(x = value, color = variable)) +
  geom_density() +
  labs(title = "Density Plot of Variables", x = "Value", y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


feature_list <- lapply(file_list,function(x){
  x <- x[,feature_columns]
})

merge_feat <- do.call(rbind,feature_list)
dim(merge_feat)
  # [1] 355037     64 

df <- merge_feat
# ここでは各列の密度を計算
df_melted <- reshape2::melt(df, id.vars = NULL)

# ggplotを使った密度プロット
ggplot(df_melted, aes(x = value, color = variable)) +
  geom_density() +
  labs(title = "Density Plot of Variables", x = "Value", y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())

feat_mean <- colMeans(merge_feat)
feat_sd <- apply( merge_feat, 2, sd)

plot(x=feat_mean,y=feat_sd)

select_feature_names <- union(
  names(feat_mean[ abs(feat_mean) > 0.2]),
  names(feat_sd[ feat_sd > 0.25 ]))


feature_list2 <- lapply(file_list,function(x){
  x <- x[,select_feature_names]
})

feat_means <- do.call(rbind, lapply(feature_list2, colMeans))
feat_sds <- do.call(rbind, lapply(feature_list2, apply, 2, sd))

feat_means_melt <- reshape2::melt(feat_means)
feat_sds_melt <- reshape2::melt(feat_sds)

df <- data.frame( patient = feat_means_melt$Var1,
                  feature = feat_means_melt$Var2,
                  mean = feat_means_melt$value,
                  sd = feat_means_melt$value)

library(RColorBrewer)
ggplot(df, aes(x = feature, y = patient)) +
  geom_point(aes(fill=mean, size=sd), shape=22) +
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  # scale_size_continuous(range = c(10,1), limits = c(0,1)) +
  labs(fill = "Mean", size = "SD") +
  ylab("") +
  xlab("")


level1 <- unique(df$patient)
level2 <- unique(df$feature)
# 水準に基づき整数の値を割り振り。長方形の中心座標となる。
df$ycenter <- as.numeric(factor(df$patient,levels = level1))
df$xcenter <- as.numeric(factor(df$feature,levels = level2))

# 最大値で割って、0.1-1のサイズにclip
cellsize_vmax <- 
df$cellsize <- pmax(pmin(df$sd / cellsize_vmax, 1),0.1)

# 長方形の座標
df$xmin <- df$xcenter - df$cellsize/2
df$ymin <- df$ycenter - df$cellsize/2
df$xmax <- df$xmin + df$cellsize
df$ymax <- df$ymin + df$cellsize

ggplot(df, aes(x = feature, y = patient)) +
  geom_rect(aes(xmin=xmin,ymin=ymin,xmax = xmax,ymax = ymax,
                fill=mean),
            color="black" # 枠線
            ) +
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  labs(fill = "SD") +
  ylab("") +
  xlab("") +  
  theme(aspect.ratio = 1)


pheatmap(feat_means)



feat_means <- do.call(rbind, lapply(feature_list, colMeans))
feat_sds <- do.call(rbind, lapply(feature_list, apply, 2, sd))
pheatmap(feat_means)

feat_means2 <- feat_means[, apply(feat_means, 2, function(x){
  max(abs(x)) > 0.3
  })]

pheatmap(feat_means2)
