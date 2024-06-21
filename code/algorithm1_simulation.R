library(igraph)
library(Matrix)
require(matrixcalc)
library(lava)
library(rBeta2009)


nodes <- 1000
sel_n <- 100
#generate beta distribution with a=1 and b=2.
gi <- rbeta(nodes, 1, 2)
HW_curve <- matrix(NA,nrow=nodes,ncol=3)
for (i in 1:nodes){
  HW_curve[i,1] <- gi[i]^2
  HW_curve[i,2] <- 2*gi[i]*(1-gi[i])
  HW_curve[i,3] <- (1-gi[i])^2
}
pm <- HW_curve%*%t(HW_curve)
#generate link matrix by gaussian kernel with sigma=1
nodes = 100
x1 <- rbeta(nodes, 1, 1)
y1 <- rbeta(nodes, 1, 1)
pm <- matrix(NA,nrow=nodes,ncol=nodes)
for(i in 1:nodes){
  for (j in 1:nodes){
    pm[i,j] <- exp(-(x1[i]-y1[j])^2)
  }
}
diag(pm) <- 0

x1 <- rbeta(nodes, 1, 2)
y1 <- rbeta(nodes, 1, 1)
pm1 <- matrix(NA,nrow=nodes,ncol=nodes)
for(i in 1:nodes){
  for (j in 1:nodes){
    pm1[i,j] <- exp(-(x1[i]-y1[j])^2)
  }
}
diag(pm1) <- 0

sample_time <- 100



T1 <- NA
rep_time <- 50
sample_time <- 100
nodes <- 1000
sel_n <- 100
sel_nodes_list <- sample(1:nodes,sel_n)

A1 <- sample_from_graph(adjm,sel_n,sel_nodes_list)
nodes <- 100

#let the sample time st=200
st <- 200
nodes <- 100
sel_n <- nodes

sample_node_index <- matrix(NA,nrow=st,ncol=sel_n)
for(i in 1:st){
  sample_node_index[i,] <- sample(1:nodes,sel_n)
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

rep_time <- 50
T1 <- NA
nodes <- 1000
st <- 100
sel_n <- 100

#generate beta distribution with a=1 and b=2.
gi <- rbeta(nodes, 1, 2)
HW_curve <- matrix(NA,nrow=nodes,ncol=3)
for (i in 1:nodes){
  HW_curve[i,1] <- gi[i]^2
  HW_curve[i,2] <- 2*gi[i]*(1-gi[i])
  HW_curve[i,3] <- (1-gi[i])^2
}
pm <- HW_curve%*%t(HW_curve)
diag(pm) <- 0

#stochastic block model
nodes <- 1000
p <- 0.6
q <- 0.3
pi1 <- 0.6
pi2 <- 0.4
pm <- matrix(NA,nrow=nodes,ncol=nodes)
n1 <- pi1*nodes
n2 <- pi2*nodes
pm[1:n1,1:n1] = p^2
pm[1:n1,(n1+1):nodes] = p*q
pm[(n1+1):nodes,(n1+1):nodes] = q^2
pm[(n1+1):nodes,1:n1] = p*q

sample_node_index <- matrix(NA,nrow=st,ncol=sel_n1)
for(i in 1:st){
  sample_node_index[i,] <- sample(1:nodes,sel_n1)
}

rep_time <- 50
T1 <- NA
nodes <- 1000
st <- 100
sel_n <- 100
sel_n1 <- 300
sel_n2 <- 500
sel_n3 <- 1000

T1_100 <- matrix(NA,nrow=5,ncol=rep_time)
T2_100 <- matrix(NA,nrow=5,ncol=rep_time)
T3_100 <- matrix(NA,nrow=5,ncol=rep_time)

T1_300 <- matrix(NA,nrow=5,ncol=rep_time)
T2_300 <- matrix(NA,nrow=5,ncol=rep_time)
T3_300 <- matrix(NA,nrow=5,ncol=rep_time)

T1_500 <- matrix(NA,nrow=5,ncol=rep_time)
T2_500 <- matrix(NA,nrow=5,ncol=rep_time)
T3_500 <- matrix(NA,nrow=5,ncol=rep_time)

T1_1000 <- matrix(NA,nrow=5,ncol=rep_time)
T2_1000 <- matrix(NA,nrow=5,ncol=rep_time)
T3_1000 <- matrix(NA,nrow=5,ncol=rep_time)

sel_n <- 100
sample_node_index <- matrix(NA,nrow=st,ncol=sel_n)
for(i in 1:st){
  sample_node_index[i,] <- sample(1:nodes,sel_n)
}

for (l in 1:5){
  sel_n <- sel_n
  for (k in 1:rep_time){
    ibg <- select_adjm(pm,nodes)
    ibg1 <- select_adjm(pm,nodes)
    
    g=graph_from_adjacency_matrix(ibg,"undirected")
    g1=graph_from_adjacency_matrix(ibg1,"undirected")
    
    A = as.matrix(get.adjacency(g))
    A_1 = as.matrix(get.adjacency(g1))
    
    A11 <- get_mean_adjm(A, sel_n, sample_node_index)
    A22 <- get_mean_adjm(A_1, sel_n, sample_node_index)
    
    sigma1 <- get_sigma(A11, sel_n)
    sigma2 <- get_sigma(A22, sel_n)
    
    #construct Z
    Z_diag <- NA
    for (i in 1:nodes){
      Z_diag[i] <- sample(c(1/sqrt(sel_n),-1/sqrt(sel_n)), 1, prob=c(0.5,0.5))
    }
    
    Z <- construct_Z(A11, A22, sigma1, sigma2, st, sel_n, Z_diag)
    
    T1_100[l,k] <- tr(Z)
    T2_100[l,k] <- tr(Z%*%Z)
    T3_100[l,k] <- tr(Z%*%Z%*%Z)
  }
}

mean(T3_100[1,])
var(T2_100[1,])
write.csv(T3_100,"C:/Users/86181/Desktop/sbm_T3_100.csv")

sample_node_index <- matrix(NA,nrow=st,ncol=sel_n1)
for(i in 1:st){
  sample_node_index[i,] <- sample(1:nodes,sel_n1)
}
for (l in 1:5){
  sel_n <- sel_n1
  for (k in 1:rep_time){
    ibg <- select_adjm(pm,nodes)
    ibg1 <- select_adjm(pm,nodes)
    
    g=graph_from_adjacency_matrix(ibg,"undirected")
    g1=graph_from_adjacency_matrix(ibg1,"undirected")
    
    A = as.matrix(get.adjacency(g))
    A_1 = as.matrix(get.adjacency(g1))
    
    A11 <- get_mean_adjm(A, sel_n, sample_node_index)
    A22 <- get_mean_adjm(A_1, sel_n, sample_node_index)
    
    sigma1 <- get_sigma(A11, sel_n)
    sigma2 <- get_sigma(A22, sel_n)
    
    #construct Z
    Z_diag <- NA
    for (i in 1:nodes){
      Z_diag[i] <- sample(c(1/sqrt(sel_n),-1/sqrt(sel_n)), 1, prob=c(0.5,0.5))
    }
    
    Z <- construct_Z(A11, A22, sigma1, sigma2, st, sel_n, Z_diag)
    
    T1_300[l,k] <- tr(Z)
    T2_300[l,k] <- tr(Z%*%Z)
    T3_300[l,k] <- tr(Z%*%Z%*%Z)
  }
}

mean(T3_300)
var(T2_300[5,])
write.csv(T3_300,"C:/Users/86181/Desktop/sbm_T3_300.csv")

sample_node_index <- matrix(NA,nrow=st,ncol=sel_n2)
for(i in 1:st){
  sample_node_index[i,] <- sample(1:nodes,sel_n2)
}

for (l in 1:5){
  sel_n <- sel_n2
  for (k in 1:rep_time){
    ibg <- select_adjm(pm,nodes)
    ibg1 <- select_adjm(pm,nodes)
    
    g=graph_from_adjacency_matrix(ibg,"undirected")
    g1=graph_from_adjacency_matrix(ibg1,"undirected")
    
    A = as.matrix(get.adjacency(g))
    A_1 = as.matrix(get.adjacency(g1))
    
    A11 <- get_mean_adjm(A, sel_n, sample_node_index)
    A22 <- get_mean_adjm(A_1, sel_n, sample_node_index)
    
    sigma1 <- get_sigma(A11, sel_n)
    sigma2 <- get_sigma(A22, sel_n)
    
    #construct Z
    Z_diag <- NA
    for (i in 1:nodes){
      Z_diag[i] <- sample(c(1/sqrt(sel_n),-1/sqrt(sel_n)), 1, prob=c(0.5,0.5))
    }
    
    Z <- construct_Z(A11, A22, sigma1, sigma2, st, sel_n, Z_diag)
    
    T1_500[l,k] <- tr(Z)
    T2_500[l,k] <- tr(Z%*%Z)
    T3_500[l,k] <- tr(Z%*%Z%*%Z)
  }
}

mean(T3_500)
var(T1_500[5,])
write.csv(T3_500,"C:/Users/86181/Desktop/sbm_T3_500.csv")


sample_node_index <- matrix(NA,nrow=st,ncol=sel_n3)
for(i in 1:st){
  sample_node_index[i,] <- sample(1:nodes,sel_n3)
}
for (l in 1:5){
  sel_n <- sel_n3
  for (k in 1:rep_time){
    ibg <- select_adjm(pm,nodes)
    ibg1 <- select_adjm(pm,nodes)
    
    g=graph_from_adjacency_matrix(ibg,"undirected")
    g1=graph_from_adjacency_matrix(ibg1,"undirected")
    
    A = as.matrix(get.adjacency(g))
    A_1 = as.matrix(get.adjacency(g1))
    
    A11 <- get_mean_adjm(A, sel_n, sample_node_index)
    A22 <- get_mean_adjm(A_1, sel_n, sample_node_index)
    
    sigma1 <- get_sigma(A11, sel_n)
    sigma2 <- get_sigma(A22, sel_n)
    
    #construct Z
    Z_diag <- NA
    for (i in 1:nodes){
      Z_diag[i] <- sample(c(1/sqrt(sel_n),-1/sqrt(sel_n)), 1, prob=c(0.5,0.5))
    }
    
    Z <- construct_Z(A11, A22, sigma1, sigma2, st, sel_n, Z_diag)
    
    T1_1000[l,k] <- tr(Z)
    T2_1000[l,k] <- tr(Z%*%Z)
    T3_1000[l,k] <- tr(Z%*%Z%*%Z)
  }
}
mean(T2_1000[2,])
var(T2_1000[5,])



write.csv(T3_1000,'C:/Users/86181/Desktop/sbm_T3_1000.csv')



