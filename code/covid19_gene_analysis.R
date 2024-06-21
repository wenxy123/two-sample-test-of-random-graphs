library(limma)
library(Glimma)
library(edgeR)
library(openxlsx)
data<- read.xlsx("C:/Users/86181/Desktop/real_data/sar_cov_2.xlsx", sheet = 1)
data <- data[,-2]
dim(data)
row.names(data) = data[,1]
data <- data[,-1]
samplenames <- colnames(data)

geneid <- rownames(data)
col.group <- group
design <- model.matrix(~0+group+lane)
