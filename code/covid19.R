library(igraph)
library(ggplot2)
require(gridExtra)
library(sand)
library(igraphdata)
library(network)
library(fields)
library(lava)
library(readr)
library(data.table)

dir1 = "D:/桌面的数据/project1_network_hypothesis/real_data3/covid19"
dir2 = "D:/桌面的数据/project1_network_hypothesis/real_data3/covid19_1"
dir3 = "D:/桌面的数据/project1_network_hypothesis/real_data3/health"

file_list1 = list.files(path = dir1, recursive = TRUE,full.names = TRUE)
file_list2 = list.files(path = dir2, recursive = TRUE,full.names = TRUE)
file_list3 = list.files(path = dir3, recursive = TRUE,full.names = TRUE)

data1 <- read.csv(file_list1[5])
data2 <- read.csv(file_list2[5])
data3 <- read.csv(file_list3[5])

edgelists = as.matrix(data1[,2:3])
type = 1
ff = graph_from_edgelist(edgelists, directed=FALSE)
A <- as.matrix(get.adjacency(ff))

edgelists1 = as.matrix(data2[,2:3])
type = 1
ff1 = graph_from_edgelist(edgelists1, directed=FALSE)
A1 <- as.matrix(get.adjacency(ff1))

edgelists2 = as.matrix(data3[,2:3])
type = 1
ff2 = graph_from_edgelist(edgelists2, directed=FALSE)
A2 <- as.matrix(get.adjacency(ff2))

r1 <- nrow(A)
r2 <- nrow(A1)
r3 <- nrow(A2)

####mean
select_n = rep(10,500)
MT <- 100
repo_T2 <- matrix(NA,nrow=length(select_n),ncol=MT)
repo_T3 <- matrix(NA,nrow=length(select_n),ncol=MT)
repo_T21 <- matrix(NA,nrow=length(select_n),ncol=MT)
repo_T31 <- matrix(NA,nrow=length(select_n),ncol=MT)

for(l in 473:length(select_n)){
  sel_n <- select_n[l]
  for (i in 1:MT){
    st <- 100
    sample_node_index1 <- matrix(NA,nrow=st,ncol=sel_n)
    sample_node_index2 <- matrix(NA,nrow=st,ncol=sel_n)
    sample_node_index3 <- matrix(NA,nrow=st,ncol=sel_n)
    
    for(m in 1:st){
      sample_node_index1[m,] <- sample(1:r1,sel_n)
      sample_node_index2[m,] <- sample(1:r2,sel_n)
      sample_node_index3[m,] <- sample(1:r3,sel_n)
    }
    
    A11 <- get_mean_adjm(A, sel_n, sample_node_index1)
    A22 <- get_mean_adjm(A1, sel_n, sample_node_index2)
    A33 <- get_mean_adjm(A2, sel_n, sample_node_index3)
    
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
    Z1[is.infinite(Z1)]=100
    Z2[is.na(Z2)]=0
    Z2[is.infinite(Z2)]=100
    
    repo_T2[l,i] = tr(Z1%*%Z1)
    repo_T3[l,i] = tr(Z1%*%Z1%*%Z1)
    
    repo_T21[l,i] = tr(Z2%*%Z2)
    repo_T31[l,i] = tr(Z2%*%Z2%*%Z2)
  }
}


T2_var <- NA
T3_var <- NA
for (i in 1:15){
  T2_var[i] = var(repo_T2[i,])
  T3_var[i] = var(repo_T3[i,])
}
T2_var
T3_var

write.csv(repo_T2,"C:/Users/86181/Desktop/covid_mean_T2.csv")
write.csv(repo_T3,"C:/Users/86181/Desktop/covid_mean_T3.csv")
write.csv(repo_T21,"C:/Users/86181/Desktop/covid_mean_T21.csv")
write.csv(repo_T31,"C:/Users/86181/Desktop/covid_mean_T31.csv")

####embed
select_n = rep(20,500)
MT <- 100
repo_T2 <- matrix(NA,nrow=length(select_n),ncol=MT)
repo_T3 <- matrix(NA,nrow=length(select_n),ncol=MT)
repo_T21 <- matrix(NA,nrow=length(select_n),ncol=MT)
repo_T31 <- matrix(NA,nrow=length(select_n),ncol=MT)

for(l in 445:length(select_n)){
  sel_n <- select_n[l]
  for (i in 1:MT){
    st <- 100
    sample_node_index1 <- matrix(NA,nrow=st,ncol=sel_n)
    sample_node_index2 <- matrix(NA,nrow=st,ncol=sel_n)
    sample_node_index3 <- matrix(NA,nrow=st,ncol=sel_n)
    
    for(m in 1:st){
      sample_node_index1[m,] <- sample(1:r1,sel_n)
      sample_node_index2[m,] <- sample(1:r2,sel_n)
      sample_node_index3[m,] <- sample(1:r3,sel_n)
    }
    
    A11 <- get_embed_adjm(A, sel_n, sample_node_index1,3)
    A22 <- get_embed_adjm(A1, sel_n, sample_node_index2,3)
    A33 <- get_embed_adjm(A2, sel_n, sample_node_index3,3)
    
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
    Z1[is.infinite(Z1)]=100
    Z2[is.na(Z2)]=0
    Z2[is.infinite(Z2)]=100
    
    repo_T2[l,i] = tr(Z1%*%Z1)
    repo_T3[l,i] = tr(Z1%*%Z1%*%Z1)
    
    repo_T21[l,i] = tr(Z2%*%Z2)
    repo_T31[l,i] = tr(Z2%*%Z2%*%Z2)
  }
}

T2_var <- NA
T3_var <- NA
for (i in 1:15){
  T2_var[i] = var(repo_T2[i,])
  T3_var[i] = var(repo_T3[i,])
}
T2_var
T3_var

write.csv(repo_T2,"C:/Users/86181/Desktop/covid_embed_T2.csv")
write.csv(repo_T3,"C:/Users/86181/Desktop/covid_embed_T3.csv")
write.csv(repo_T21,"C:/Users/86181/Desktop/covid_embed_T21.csv")
write.csv(repo_T31,"C:/Users/86181/Desktop/covid_embed_T31.csv")

#############################################################
nodes_list <- rep(100,10)

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
    
    A11 <- get_embed_adjm(sel_A, sel_n, sample_node_index1,2)
    A22 <- get_embed_adjm(sel_A1, sel_n, sample_node_index2,2)
    A33 <- get_embed_adjm(sel_A2, sel_n, sample_node_index3,2)
    
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
    Z1[is.infinite(Z1)]=100
    Z2[is.na(Z2)]=0
    Z2[is.infinite(Z2)]=100
    
    repo_T1[i,k] <- tr(Z1)
    repo_T2[i,k] <- tr(Z1%*%Z1)
    repo_T3[i,k] <- tr(Z1%*%Z1%*%Z1)
    
    repo_T11[i,k] <- tr(Z2)
    repo_T22[i,k] <- tr(Z2%*%Z2)
    repo_T31[i,k] <- tr(Z2%*%Z2%*%Z2)
  }
}

