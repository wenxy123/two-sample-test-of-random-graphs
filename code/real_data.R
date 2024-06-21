library(igraph)
library(ggplot2)
require(gridExtra)
library(sand)
library(igraphdata)
library(network)
library(fields)
library(lava)

data <- read.table('C:/Users/86181/Desktop/real_data.txt')
data1 <- read.csv('C:/Users/86181/Desktop/real_data1.csv')

edgelists = as.matrix(data[,1:2])
type = 1
ff = graph_from_edgelist(edgelists, directed=FALSE)
A <- as.matrix(get.adjacency(ff))

edgelists1 = as.matrix(data1[,1:2]+1)
type = 1
ff1 = graph_from_edgelist(edgelists1, directed=FALSE)
A1 <- as.matrix(get.adjacency(ff1))

r1 <- nrow(A)
r2 <- nrow(A1)

plot(ff,vertex.label=NA,vertex.size=5,vertex.color='lightblue')

plot(ff1,vertex.label=NA,vertex.size=5,vertex.color='lightblue')

nodes_list <- rep(50,100)

MT <- 100

repo_T1 <- matrix(NA,nrow=length(nodes_list),ncol=MT)
repo_T2 <- matrix(NA,nrow=length(nodes_list),ncol=MT)
repo_T3 <- matrix(NA,nrow=length(nodes_list),ncol=MT)

repo_T11 <- matrix(NA,nrow=length(nodes_list),ncol=MT)
repo_T22 <- matrix(NA,nrow=length(nodes_list),ncol=MT)
repo_T31 <- matrix(NA,nrow=length(nodes_list),ncol=MT)

for (i in 1:length(nodes_list)){
  sel_n <- nodes_list[i]
  st <- 100
  
  sample_node_index1 <- matrix(NA,nrow=st,ncol=sel_n)
  sample_node_index2 <- matrix(NA,nrow=st,ncol=sel_n)
  sample_node_index3 <- matrix(NA,nrow=st,ncol=sel_n)
  
  for(m in 1:st){
    sample_node_index1[m,] <- sample(1:500,sel_n)
    sample_node_index2[m,] <- sample(1:500,sel_n)
    sample_node_index3[m,] <- sample(1:500,sel_n)
  }
  
  for (k in 1:MT){
    idx1 <- sample(1:r1,500)
    idx2 <- sample(1:r1,500)
    idx3 <- sample(1:r2,500)
    
    sel_A <- sample_from_graph(A,length(idx1),idx1)
    sel_A1 <- sample_from_graph(A,length(idx2),idx2)
    sel_A2 <- sample_from_graph(A1,length(idx3),idx3)
    
    A11 <- get_mean_adjm(sel_A, sel_n, sample_node_index1)
    A22 <- get_mean_adjm(sel_A1, sel_n, sample_node_index2)
    A33 <- get_mean_adjm(sel_A2, sel_n, sample_node_index3)
    
    sigma1 <- get_sigma(A11, sel_n)
    sigma2 <- get_sigma(A22, sel_n)
    sigma3 <- get_sigma(A33, sel_n)
    
    #construct Z
    Z_diag <- NA
    for (n in 1:sel_n){
      Z_diag[n] <- sample(c(1/sqrt(sel_n),-1/sqrt(sel_n)), 1, prob=c(0.5,0.5))
    }
    
    Z1 <- construct_Z(A11, A22, sigma1, sigma2, st, sel_n, Z_diag)
    Z2 <- construct_Z(A11, A33, sigma1, sigma3, st, sel_n, Z_diag)
    
    #fill NaN value
    Z1[is.na(Z1)]=0
    Z2[is.na(Z2)]=0
    
    repo_T1[i,k] <- tr(Z1)
    repo_T2[i,k] <- tr(Z1%*%Z1)
    repo_T3[i,k] <- tr(Z1%*%Z1%*%Z1)
    
    repo_T11[i,k] <- tr(Z2)
    repo_T22[i,k] <- tr(Z2%*%Z2)
    repo_T31[i,k] <- tr(Z2%*%Z2%*%Z2)
  }
}

write.csv(repo_T1,"C:/Users/86181/Desktop/rd_repo_T1.csv")
write.csv(repo_T2,"C:/Users/86181/Desktop/rd_repo_T2.csv")
write.csv(repo_T3,"C:/Users/86181/Desktop/rd_repo_T3.csv")
write.csv(repo_T11,"C:/Users/86181/Desktop/rd_repo_T11.csv")
write.csv(repo_T22,"C:/Users/86181/Desktop/rd_repo_T22.csv")
write.csv(repo_T31,"C:/Users/86181/Desktop/rd_repo_T31.csv")




