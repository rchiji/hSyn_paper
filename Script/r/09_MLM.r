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

