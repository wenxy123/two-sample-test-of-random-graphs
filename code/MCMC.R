library(igraph)
library(Matrix)
require(matrixcalc)
library(lava)
library(rBeta2009)
library(Rcpp)

Rcpp::sourceCpp('main_function.cpp')

#random dot product graph
generate_RDPG_pm <- function(gi,nodes){
  HW_curve <- matrix(NA,nrow=nodes,ncol=3)
  for (i in 1:nodes){
    HW_curve[i,1] <- gi[i]^2
    HW_curve[i,2] <- 2*gi[i]*(1-gi[i])
    HW_curve[i,3] <- (1-gi[i])^2
  }
  pm <- HW_curve%*%t(HW_curve)
  return(pm)
}

#random graph with kernel link
generate_kernel_pm <- function(x1,y1,nodes){
  pm <- matrix(NA,nrow=nodes,ncol=nodes)
  for(i in 1:nodes){
    for (j in 1:nodes){
      pm[i,j] <- exp(-(x1[i]-y1[j])^2)
    }
  }
  return(pm)
}

#stochastic block model
generate_sbm_pm <- function(p,q,nodes){
  pm <- matrix(NA,nrow=nodes,ncol=nodes)
  pi1 <- 0.6
  pi2 <- 0.4
  n1 <- pi1*nodes
  n2 <- pi2*nodes
  pm[1:n1,1:n1] = p^2
  pm[1:n1,(n1+1):nodes] = p*q
  pm[(n1+1):nodes,(n1+1):nodes] = q^2
  pm[(n1+1):nodes,1:n1] = p*q
  return(pm)
}

select_adjm <- function(pm,nodes){
  ibg <- matrix(NA,nrow=nodes,ncol=nodes)
  for(i in 1:nodes){
    for (j in i:nodes){
      ibg[i,j] = sample(c(1,0), 1, prob=c(pm[i,j],1-pm[i,j]))
      ibg[j,i] = ibg[i,j]
    }
  }
  return(ibg)
}

MT <- 300
nodes_list <- c(100,200,300,400,500,600,700,800,900,1000)

repo_T1 <- matrix(NA,nrow=length(nodes_list),ncol=MT)
repo_T2 <- matrix(NA,nrow=length(nodes_list),ncol=MT)
repo_T3 <- matrix(NA,nrow=length(nodes_list),ncol=MT)

repo_T11 <- matrix(NA,nrow=length(nodes_list),ncol=MT)
repo_T22 <- matrix(NA,nrow=length(nodes_list),ncol=MT)
repo_T31 <- matrix(NA,nrow=length(nodes_list),ncol=MT)


for (i in 1:length(nodes_list)){
  nodes <- nodes_list[i]
  sel_n <- nodes*0.1
  st <- sel_n
  
  #x1 <- rbeta(nodes, 1, 1)
  #y1 <- rbeta(nodes, 1, 1)
  #p1 <- 0.6
  #q1 <- 0.3
  #pm1 <- generate_kernel_pm(x1,y1,nodes)
  
  #x2 <- rbeta(nodes, 1, 2)
  #y2 <- rbeta(nodes, 1, 1)
  #p2 <- 0.7
  #q2 <- 0.4
  #pm2 <- generate_kernel_pm(x2,y2,nodes)
  gi1 <- rbeta(nodes, 1, 1)
  pm1 <- generate_RDPG_pm(gi1,nodes)
  
  gi2 <- rbeta(nodes, 1, 2)
  pm2 <- generate_RDPG_pm(gi2,nodes)
  
  sample_node_index <- matrix(NA,nrow=st,ncol=sel_n)
  for(m in 1:st){
    #sample_node_index[m,] <- sample(1:nodes,sel_n)
    sample_node_index[m,] <- sort(sample(1:nodes,sel_n))
  }
  
  for (k in 1:MT){
    ibg <- select_adjm(pm1,nodes)
    ibg1 <- select_adjm(pm1,nodes)
    ibg2 <- select_adjm(pm2,nodes)
    
    g=graph_from_adjacency_matrix(ibg,"undirected")
    g1=graph_from_adjacency_matrix(ibg1,"undirected")
    g2=graph_from_adjacency_matrix(ibg2,"undirected")
    
    A = as.matrix(get.adjacency(g))
    A_1 = as.matrix(get.adjacency(g1))
    A_2 = as.matrix(get.adjacency(g2))
    
    A11 <- get_mean_adjm(A, sel_n, sample_node_index)
    A22 <- get_mean_adjm(A_1, sel_n, sample_node_index)
    A33 <- get_mean_adjm(A_2, sel_n, sample_node_index)
    #A11 <- get_embed_adjm(A, sel_n, sample_node_index,3)
    #A22 <- get_embed_adjm(A_1, sel_n, sample_node_index,3)
    #A33 <- get_embed_adjm(A_2, sel_n, sample_node_index,3)
    
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
    print(k)
  }
  print(i)
}

write.csv(repo_T1,"/public3/home/scg5453/wenxy/lpvs_repo_T1.csv")
write.csv(repo_T2,"/public3/home/scg5453/wenxy/lpvs_repo_T2.csv")
write.csv(repo_T3,"/public3/home/scg5453/wenxy/lpvs_repo_T3.csv")
write.csv(repo_T11,"/public3/home/scg5453/wenxy/lpvs_repo_T11.csv")
write.csv(repo_T22,"/public3/home/scg5453/wenxy/lpvs_repo_T22.csv")
write.csv(repo_T31,"/public3/home/scg5453/wenxy/lpvs_repo_T31.csv")

#alpha <- 0.05
#qnorm(1-(alpha/2), mean = 0, sd = 4, lower.tail = TRUE, log.p = FALSE)




