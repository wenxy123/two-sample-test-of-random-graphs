install.packages("Hmisc")
library(Hmisc)
library(WGCNA)
library("readxl")

data <- read_excel('C:/Users/86181/Desktop/sar_cov_2.xlsx',sheet = 1,col_names = TRUE)
data1 <- data[,-2]
data11 <- data1[,-1]
cor_mat = rcorr(t(as.matrix(data11)),type = "pearson")
cor_mat1 <- cor_mat$r
a <- data1[,1]
df <- data.frame(data11)
a<- list(data1[,1])
rownames(data11) = c(data1[,1])
colnames(cor_mat) = data1[,1]