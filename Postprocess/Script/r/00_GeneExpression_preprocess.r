library(dplyr)
library(DESeq2)
library(edgeR)

data <- read.delim("00_src/gene_read_count.txt", sep = "\t", header = TRUE, row.names = 1)
# Gene ID duplication
# Gene IDs such as “PAR_Y_XX” appear to denote homologous genes located between the X and Y chromosomes (i.e., genes in the pseudoautosomal region).
# We identified 44 such duplicated entries; however, for all of them, the count values were identical between the homologous gene pairs. Therefore, we removed the PAR_Y entries in Excel.

data <- data[!is.na(data$Gene.type) & data$Gene.type == "protein_coding", ]
data <- data[,1:81]
data <- na.omit(data)

data$id <- rownames(data)
data <- data %>% 
  mutate(gene = paste(data$id, data$gene_name, sep = "_"))
row.names(data) <- data$gene
data <- data[,2:81]

write.table(data, "01_formatted/gene_read_count_preprocessed.txt", sep = "\t", row.names = T, col.names = NA)

