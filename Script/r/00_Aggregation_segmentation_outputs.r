library(dplyr)
library(ggsci)
library(ggplot2)
library(pheatmap)
library(patchwork)

file_paths <- list.files(path = "00_src/Annotation_ratio/", full.names = TRUE)
samplenames <- sapply(file_paths, function(x) {
  x <- basename(x)
  gsub(pattern = "_HE.*",replacement = "",x)
  })
file_list <- lapply(file_paths, read.delim)
names(file_list) <- samplenames

levels <- c("adipose", "Fibro(dense,irregular)", "Fibro(dense,regular)", 
"Fibro(loose)", "Immune cells", "lining", "muscle", "plasma", 
"RBC", "Stroma", "vessel", "vessel(large)")

# calculate_area_ratio
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

summary_list <- lapply(summary_list, function(x){
  x$Area_ratio <- x$Area_sum / sum(x$Area_sum)
  x$average_area <- x$Area_sum / x$count
  return(x)
  })

res <- do.call(rbind, lapply(summary_list,function(x){
  value <- x$Area_ratio
  names(value) <- x$Classification
  value <- value[levels]
  value[is.na(value)] <- 0
  return(value)
}))

write.table(res, "00_src/annotations_for_analysis.txt", sep = "\t", quote = FALSE, col.names = NA)
# -> Clinical information was manually added, and the file was saved as "annotations_full.txt".

