library(randomForestSRC)
library(pheatmap)
library(randomForest)
library(pROC)

df <- read.delim(file = "00_src/annotations_full.txt", sep = "\t", row.names = 1)
df_OA <- df[df$Diagnosis == "OA",]


# MLM
set.seed(123)

custom_labels <- c(
  "Fibrous_tissue__dense_irregular" = "Fibrous tissue (dense.irregular)",
  "Fibrous_tissue__dense_regular" = "Fibrous tissue (dense.regular)",
  "Fibrous_tissue__loose" = "Fibrous tissue (loose)",
  "Micro_vessel" = "Micro vessel",
  "Large_vessel" = "Large vessel",
  "Pre__KOOS_Pain" = "Pain (Pre Ope)",
  "Pre__KOOS_Symptom" = "Symptom (Pre Ope)",
  "Post_3Month__KOOS_Pain" = "Pain (Post-3M Ope)",
  "Post_3Month__KOOS_Symptom" = "Symptom (Post-3M Ope)",
  "Post_12Month__KOOS_Pain" = "Pain (Post-12M Ope)",
  "Post_12Month__KOOS_Symptom" = "Symptom (Post-12M Ope)"
)

dependent_vars <- c("Pre__KOOS_Pain", "Pre__KOOS_Symptom", "Post_3Month__KOOS_Pain",
                    "Post_3Month__KOOS_Symptom", "Post_12Month__KOOS_Pain", "Post_12Month__KOOS_Symptom")

data <- df_OA[,c(6:17,18,19,23,24,28,29)]
data <- na.omit(data)

obj <- rfsrc(formula = get.mv.formula(dependent_vars),
             data = data, ntree = 2000,
             importance = TRUE, 
             nsplit = 10, 
             splitrule = "mahalanobis", 
             seed = -37)

tmp <- get.mv.vimp(obj, standardize = TRUE)

rownames(tmp) <- ifelse(rownames(tmp) %in% names(custom_labels),
                        custom_labels[rownames(tmp)], rownames(tmp))
colnames(tmp) <- ifelse(colnames(tmp) %in% names(custom_labels),
                        custom_labels[colnames(tmp)], colnames(tmp))

range_value <- range(tmp, na.rm = TRUE)
max_abs <- max(abs(range_value))
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)

pdf("99_Fig/fig3/mlm_tissue_standardize.pdf" , width = 4.5, height = 3)
pheatmap(t(tmp),
         angle_col = 90,
         fontsize = 5,
         color = colors,
         breaks = breaks,
         cellwidth = 10,
         cellheight = 10,
         treeheight_row = 10,
         treeheight_col = 10,
)
dev.off()


# ROC
set.seed(123)

roc_data <- data[, c(1:12)]
predictors <- roc_data[,1:12]

roc_data$Pain_pre <- data$Pre__KOOS_Pain
roc_data$Symptom_pre <- data$Pre__KOOS_Symptom
roc_data$Pain_3M <- data$Post_3Month__KOOS_Pain
roc_data$Symptom_3M <- data$Post_3Month__KOOS_Symptom
roc_data$Pain_12M <- data$Post_12Month__KOOS_Pain
roc_data$Symptom_12M <- data$Post_12Month__KOOS_Symptom

## post12M symptom
median_val <- median(roc_data$Symptom_12M, na.rm = TRUE)

roc_data$binary <- ifelse(roc_data$Symptom_12M < median_val, 1, 0)

target <- as.factor(roc_data$binary)

rf_model <- randomForest(x = predictors, y = target)

oob_prob <- rf_model$votes[, "1"]

tapply(oob_prob, target, summary)

roc_oob <- roc(response = as.numeric(as.character(target)),
               predictor = oob_prob,
               levels = c(0, 1),
               direction = "<")

auc_val <- auc(roc_oob)
v <- var(roc_oob, method = "delong")
se <- sqrt(as.numeric(v))
z  <- (as.numeric(auc_val) - 0.5) / se
p_val <- 2 * (1 - pnorm(abs(z)))

### plot
lab_txt <- sprintf("AUC = %.2f\nP = %.3g", as.numeric(auc_val), p_val)

png("99_Fig/sup_fig3/roc_curve_post12_symptom.png", width = 6, height = 6, units = "cm", res = 600)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,
     legacy.axes = TRUE,                   
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()

pdf("99_Fig/sup_fig3/roc_curve_post12_symptom.pdf", width = 3, height = 3, useDingbats = FALSE)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,       
     legacy.axes = TRUE,            
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()


## post12M pain
median_val <- median(roc_data$Pain_12M, na.rm = TRUE)

roc_data$binary <- ifelse(roc_data$Pain_12M < median_val, 1, 0)

target <- as.factor(roc_data$binary)

rf_model <- randomForest(x = predictors, y = target)

oob_prob <- rf_model$votes[, "1"]

tapply(oob_prob, target, summary)

roc_oob <- roc(response = as.numeric(as.character(target)),
               predictor = oob_prob,
               levels = c(0, 1),
               direction = ">")

auc_val <- auc(roc_oob)
v <- var(roc_oob, method = "delong")
se <- sqrt(as.numeric(v))
z  <- (as.numeric(auc_val) - 0.5) / se
p_val <- 2 * (1 - pnorm(abs(z)))

### plot
lab_txt <- sprintf("AUC = %.2f\nP = %.3g", as.numeric(auc_val), p_val)

png("99_Fig/sup_fig3/roc_curve_post12_pain.png", width = 6, height = 6, units = "cm", res = 600)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,
     legacy.axes = TRUE,                   
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()

pdf("99_Fig/sup_fig3/roc_curve_post12_pain.pdf", width = 3, height = 3, useDingbats = FALSE)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,       
     legacy.axes = TRUE,            
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()


## post3M symptom
median_val <- median(roc_data$Symptom_3M, na.rm = TRUE)

roc_data$binary <- ifelse(roc_data$Symptom_3M < median_val, 1, 0)

target <- as.factor(roc_data$binary)

rf_model <- randomForest(x = predictors, y = target)

oob_prob <- rf_model$votes[, "1"]

tapply(oob_prob, target, summary)

roc_oob <- roc(response = as.numeric(as.character(target)),
               predictor = oob_prob,
               levels = c(0, 1),
               direction = ">")

auc_val <- auc(roc_oob)
v <- var(roc_oob, method = "delong")
se <- sqrt(as.numeric(v))
z  <- (as.numeric(auc_val) - 0.5) / se
p_val <- 2 * (1 - pnorm(abs(z)))

### plot
lab_txt <- sprintf("AUC = %.2f\nP = %.3g", as.numeric(auc_val), p_val)

png("99_Fig/sup_fig3/roc_curve_post3_symptom.png", width = 6, height = 6, units = "cm", res = 600)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,
     legacy.axes = TRUE,                   
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()

pdf("99_Fig/sup_fig3/roc_curve_post3_symptom.pdf", width = 3, height = 3, useDingbats = FALSE)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,       
     legacy.axes = TRUE,            
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()


## post3M pain
median_val <- median(roc_data$Pain_3M, na.rm = TRUE)

roc_data$binary <- ifelse(roc_data$Pain_3M < median_val, 1, 0)

target <- as.factor(roc_data$binary)

rf_model <- randomForest(x = predictors, y = target)

oob_prob <- rf_model$votes[, "1"]

tapply(oob_prob, target, summary)

roc_oob <- roc(response = as.numeric(as.character(target)),
               predictor = oob_prob,
               levels = c(0, 1),
               direction = ">")

auc_val <- auc(roc_oob)
v <- var(roc_oob, method = "delong")
se <- sqrt(as.numeric(v))
z  <- (as.numeric(auc_val) - 0.5) / se
p_val <- 2 * (1 - pnorm(abs(z)))

### plot
lab_txt <- sprintf("AUC = %.2f\nP = %.3g", as.numeric(auc_val), p_val)

png("99_Fig/sup_fig3/roc_curve_post3_pain.png", width = 6, height = 6, units = "cm", res = 600)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,
     legacy.axes = TRUE,                   
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()

pdf("99_Fig/sup_fig3/roc_curve_post3_pain.pdf", width = 3, height = 3, useDingbats = FALSE)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,       
     legacy.axes = TRUE,            
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()


## pre symptom
median_val <- median(roc_data$Symptom_pre, na.rm = TRUE)

roc_data$binary <- ifelse(roc_data$Symptom_pre < median_val, 1, 0)

target <- as.factor(roc_data$binary)

rf_model <- randomForest(x = predictors, y = target)

oob_prob <- rf_model$votes[, "1"]

tapply(oob_prob, target, summary)

roc_oob <- roc(response = as.numeric(as.character(target)),
               predictor = oob_prob,
               levels = c(0, 1),
               direction = ">")

auc_val <- auc(roc_oob)
v <- var(roc_oob, method = "delong")
se <- sqrt(as.numeric(v))
z  <- (as.numeric(auc_val) - 0.5) / se
p_val <- 2 * (1 - pnorm(abs(z)))

### plot
lab_txt <- sprintf("AUC = %.2f\nP = %.3g", as.numeric(auc_val), p_val)

png("99_Fig/sup_fig3/roc_curve_pre_symptom.png", width = 6, height = 6, units = "cm", res = 600)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,
     legacy.axes = TRUE,                   
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()

pdf("99_Fig/sup_fig3/roc_curve_pre_symptom.pdf", width = 3, height = 3, useDingbats = FALSE)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,       
     legacy.axes = TRUE,            
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()


## pre pain
median_val <- median(roc_data$Pain_pre, na.rm = TRUE)

roc_data$binary <- ifelse(roc_data$Pain_pre < median_val, 1, 0)

target <- as.factor(roc_data$binary)

rf_model <- randomForest(x = predictors, y = target)

oob_prob <- rf_model$votes[, "1"]

tapply(oob_prob, target, summary)

roc_oob <- roc(response = as.numeric(as.character(target)),
               predictor = oob_prob,
               levels = c(0, 1),
               direction = ">")

auc_val <- auc(roc_oob)
v <- var(roc_oob, method = "delong")
se <- sqrt(as.numeric(v))
z  <- (as.numeric(auc_val) - 0.5) / se
p_val <- 2 * (1 - pnorm(abs(z)))

### plot
lab_txt <- sprintf("AUC = %.2f\nP = %.3g", as.numeric(auc_val), p_val)

png("99_Fig/sup_fig3/roc_curve_pre_pain.png", width = 6, height = 6, units = "cm", res = 600)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,
     legacy.axes = TRUE,                   
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
dev.off()

pdf("99_Fig/sup_fig3/roc_curve_pre_pain.pdf", width = 3, height = 3, useDingbats = FALSE)
par(bg = "white")
plot(roc_oob,
     col = "red",             
     lwd = 1.5,       
     legacy.axes = TRUE,            
     main = "",       
     xlab = "1 - Specificity",    
     ylab = "Sensitivity",         
     cex.lab = 0.5,                
     cex.axis = 0.4,               
     font.lab = 2,                 
     las = 1)
auc_val <- auc(roc_oob)
text(0.65, 0.05, labels = lab_txt, cex = 0.5, col = "red", font = 1)
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
#  [1] pROC_1.18.5           randomForest_4.7-1.2  pheatmap_1.0.12       randomForestSRC_3.3.3
# 
# loaded via a namespace (and not attached):
#  [1] cli_3.6.4          data.tree_1.1.0    rlang_1.1.5        stringi_1.8.7      jsonlite_2.0.0     glue_1.8.0         colorspace_2.1-1   plyr_1.8.9         htmltools_0.5.8.1  scales_1.3.0       grid_4.3.3         munsell_0.5.1      visNetwork_2.1.2   fastmap_1.2.0     
# [15] lifecycle_1.0.4    DiagrammeR_1.0.11  compiler_4.3.3     RColorBrewer_1.1-3 Rcpp_1.0.14        htmlwidgets_1.6.4  rstudioapi_0.17.1  digest_0.6.37      R6_2.6.1           parallel_4.3.3     magrittr_2.0.3     tools_4.3.3        gtable_0.3.6      


