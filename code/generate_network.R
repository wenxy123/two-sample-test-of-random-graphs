library(igraph)
library(ggplot2)
require(gridExtra)
library(sand)
library(igraphdata)
library(network)
library(sna)
library(IRdisplay)
library(sbm)
require(matrixcalc)
library(asnipe)
library(fields)
library(clusterGeneration)
library(lava)
library(rBeta2009)
library(scatterplot3d)
library(Matrix)
library(RColorBrewer)

##method 1
g <- make_empty_graph(mode='undirect') + vertices(letters[1:10])
for (i in 1:5){
  g <- g + edges(i, 10, replace = TRUE, color = "gray")
}
plot(g)
g1 <- as.undirected(g)
plot(g1)

##method 2: generate stochastic block model
pm <- cbind( c(.42, .42), c(.42, .5) )
pm
g2 <- sample_sbm(100, pref.matrix=pm, block.sizes=c(60,40))
plot(g2)
A = as.matrix(get.adjacency(g2))
x<- NA

nodes <- 100
#should be multivariate distribution
lpvs <- matrix(rnorm(200), 100, 2)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
pm <- lpvs%*%t(lpvs)
lpvs <- sample_sphere_surface(dim=100, n=2)
lpvs <- sample_dirichlet(n=2, alpha=rep(1, 100))
lpvs1 <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2)))})
pm <- lpvs%*%t(lpvs)
diag(pm) <- 0
ibg=matrix(0,nrow=nodes, ncol=nodes)

for(i in 1:nodes){
  for (j in 1:nodes){
    pm <- lpvs[i,]%*%t(lpvs)[,j]
    pm1 <- abs(pm)/(sqrt(sum(lpvs[i,]^2))*sqrt(sum(lpvs[j,]^2)))
    ibg[i,j]=sample(c(1,0), 1, prob=c(pm1,1-pm1))
  }
}

adj=get_network(t(ibg), data_format="GBI", association_index="SRI")
g=graph_from_adjacency_matrix(adj, "undirected",weighted=T)
plot(g,vertex.label="", edge.width=E(g)$weight*5)

#######
lpvs <- matrix(rnorm(200), 100, 2)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
g <- sample_dot_product(lpvs)
plot(g,vertex.label="")
A = as.matrix(get.adjacency(g))
sum(A)

lpvs <- matrix(rnorm(200,mean=1,sd=2), 20, 10)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
g <- sample_dot_product(lpvs,directed = FALSE)
g <- sample_dot_product(lpvs1,directed = FALSE)
A = as.matrix(get.adjacency(g))
embed <- embed_adjacency_matrix(g, 5)
est <- embed$X
ev <- eigen(A)
V <- ev$vectors
L <- ev$values
length(L)
dim(V)
V1 <- V[,1:2]
L1 <- sqrt(L)
embed_X <- V1%*%diag(L1[1:2])



x<-NA
for (i in 1:1000){
  #lpvs <- matrix(rnorm(4,mean=2,sd=20), 2, 2)
  #pm <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
  #pm[2,1] <- pm[1,2]
  #g <- sample_dot_product(lpvs,directed = FALSE)
  #pm <- cbind( c(.3, .3), c(.3, .3) )
  pm <- cbind( c(.42, .42), c(.42, .5) )
  g2 <- sample_sbm(100, pref.matrix=pm, block.sizes=c(60,40))
  A = as.matrix(get.adjacency(g2))
  x[i] <- A[7,5]
}
var(x)
x
x1 <- c(0.63,-0.14)
x2 <- c(0.69,0.13)
#estimate latent position
pm <- cbind( c(.42, .42), c(.42, .5) )
g2 <- sample_sbm(100, pref.matrix=pm, block.sizes=c(60,40))
embed <- embed_adjacency_matrix(g2, 2)
est_X <- embed$X

est_X <- matrix(NA,nrow=10,ncol=2)
pm <- cbind( c(.02, .42), c(.42, .05) )
g2 <- sample_sbm(100, pref.matrix=pm, block.sizes=c(60,40))
A = as.matrix(get.adjacency(g2))
image.plot(matrix((data=A), ncol=100, nrow=100))
sum(A[1:60,60:100])/2400

ev <- eigen(pm)
V <- ev$vectors
L <- ev$values
length(L)
dim(V)
V1 <- V[,1:2]
L1 <- sqrt(L)
embed_X <- V1%*%diag(L1[1:2])


p <- 0.6
q <- 0.3
pi1 <- 0.6
pi2 <- 0.4
nodes <- 100
pm <- cbind( c(p^2, p*q), c(p*q, q^2) )
g2 <- sample_sbm(nodes, pref.matrix=pm, block.sizes=c(pi1*nodes,pi2*nodes))
plot(g2,vertex.label="", vertex.size=10,vertex.color='lightblue',edge.width=E(g2)$weight*5)
A = as.matrix(get.adjacency(g2))
sum(A)
cols <- brewer.pal(3, "Blues")
pal <- colorRampPalette(cols)
image.plot(matrix((data=A), ncol=100, nrow=100),col = pal(10))

#
sigma_p <- (pi1*p^4*(1-p^2)+pi2*p*q^3*(1-p*q))/(pi1*p^2+pi2*q^2)^2
sigma_q <- (pi1*p^3*q*(1-p*q)+pi2*q^4*(1-q^2))/(pi1*p^2+pi2*q^2)^2

est_p <- NA
est_q <- NA
p <- 0.6
q <- 0.3
pi1 <- 0.6
pi2 <- 0.4
nodes <- 100
pm <- cbind( c(p^2, p*q), c(p*q, q^2) )

for (i in 1:500){
  g2 <- sample_sbm(nodes, pref.matrix=pm, block.sizes=c(pi1*nodes,pi2*nodes))
  A = as.matrix(get.adjacency(g2))
  n1 <- pi1*nodes
  n2 <- pi2*nodes
  p1 <- sum(A[1:n1,1:n1])/(n1^2)
  p2 <- sum(A[1:n1,n1:n1+n2])/(n1*n2)
  p4 <- sum(A[n1:n1+n2,n1:n1+n2])/(n2^2)
  est_pm <- cbind(c(p1,p2),c(p2,p4))
  ev <- eigen(est_pm)
  V <- ev$vectors
  L <- ev$values
  L1 <- sqrt(L)
  embed_X <- V%*%diag(L1)
  est_p[i] <- embed_X[1,1]
  est_q[i] <- embed_X[2,1]
}

Z <- matrix(NA,nrow=nodes,ncol=nodes)

m1 <- 10
rep_time <- 100
nodes <- 100
sigma1 <- matrix(NA,nrow=nodes,ncol=nodes)
sigma2 <- matrix(NA,nrow=nodes,ncol=nodes)
n1 <- pi1*nodes
n2 <- pi2*nodes
sigma1[1:n1,1:n1] <- pm[1,1]*(1-pm[1,1])
sigma1[1:n1,(n1+1):nodes] <- pm[1,2]*(1-pm[1,2])
sigma1[(n1+1):nodes,1:n1] <- pm[1,2]*(1-pm[1,2])
sigma1[(n1+1):nodes,(n1+1):nodes] <- pm[2,2]*(1-pm[2,2])
sigma2 <- sigma1

T1 <- matrix(NA,nrow=10,ncol=rep_time)
construct_stat <- function(){
  T11 <- NA
  for (i in 1:rep_time){
    mean_A1 <- matrix(0,nrow=nodes,ncol=nodes)
    mean_A2 <- matrix(0,nrow=nodes,ncol=nodes)
    Z <- matrix(NA,nrow=nodes,ncol=nodes)
    
    for (j in 1:m1){
      g1 <- sample_sbm(nodes, pref.matrix=pm, block.sizes=c(pi1*nodes,pi2*nodes))
      A1 = as.matrix(get.adjacency(g1))
      mean_A1 <- mean_A1 + A1
      g2 <- sample_sbm(nodes, pref.matrix=pm, block.sizes=c(pi1*nodes,pi2*nodes))
      A2 = as.matrix(get.adjacency(g2))
      mean_A2 <- mean_A2 + A2
    }
    mean_A1 <- mean_A1/m1
    mean_A2 <- mean_A2/m1
    #construct Z
    for (k in 1:nodes){
      for (m in 1:nodes){
        if (k==m){
          Z[k,m] <- sample(c(1/sqrt(nodes),-1/sqrt(nodes)), 1, prob=c(0.5,0.5))
        }
        if (k!=m){
          Z[k,m] <- (mean_A1[k,m]-mean_A2[k,m])/(sqrt(nodes*(sigma1[k,m]/m1 + sigma2[k,m]/m1)))
        }
      }
    }
    #T1[i] <- (1/nodes)*tr(Z%*%Z%*%Z)
    T11[i] <- tr(Z%*%Z)
  }
  return(T11)
}
for(i in 1:10){
  T1[i,] <- construct_stat()
}
write.csv(T1,'C:/Users/86181/Desktop/nodes100_Z2.csv')
mean(T1)
var(T1[4,])

#latent position from beta distribution
gi <- rbeta(100, 1, 2)
HW_curve <- matrix(NA,nrow=100,ncol=3)
for (i in 1:100){
  HW_curve[i,1] <- gi[i]^2
  HW_curve[i,2] <- 2*gi[i]*(1-gi[i])
  HW_curve[i,3] <- (1-gi[i])^2
}
#write.csv(HW_curve, 'C:/Users/86181/Desktop/generate_lpvs.csv')
scatterplot3d(HW_curve[,1:3], pch = 6, color="steelblue")
pm <- HW_curve%*%t(HW_curve)
diag(pm) <- 0

rep_time <- 1000
A1 <- matrix(0,nrow=100,ncol=100)
for(r in 1:rep_time){
  ibg <- matrix(NA,nrow=100,ncol=100)
  for(i in 1:nodes){
    for (j in i:nodes){
      ibg[i,j] = sample(c(1,0), 1, prob=c(pm[i,j],1-pm[i,j]))
      ibg[j,i] = ibg[i,j]
    }
  }
  g=graph_from_adjacency_matrix(ibg,"undirected")
  A = as.matrix(get.adjacency(g))
  A1 <- A1 + A
}
A1 <- A1/rep_time


for(i in 1:nodes){
  for (j in i:nodes){
    ibg[i,j] = sample(c(1,0), 1, prob=c(pm[i,j],1-pm[i,j]))
    ibg[j,i] = ibg[i,j]
  }
}
sum(ibg)
g=graph_from_adjacency_matrix(ibg,"undirected")
plot(g,vertex.label="")
A = as.matrix(get.adjacency(g))
sum(A)
ev <- eigen(A1)
V <- ev$vectors
L <- ev$values
V1 <- V[,1:3]
L1 <- sqrt(L[1:3])
embed_X <- V1%*%diag(L1)
plot(embed_X)
plot(HW_curve[,1:2])
embed_X2 <- abs(embed_X)
ev <- eigen(A)
V <- ev$vectors
L <- ev$values
V1 <- V[,1:3]
L1 <- sqrt(L[1:3])
lpvs <- V1%*%diag(L1)
embed_X <- V1%*%diag(L1)
scatterplot3d(abs(embed_X[,1:3]), color="steelblue")
scatterplot3d(abs(lpvs[,1:3]), color="steelblue")

#calculate sigma
s <- svd(HW_curve[,1])
D0 <- NA
for (i in 1:100){
  x01 <- HW_curve[1,]
  item <- x01%*%t(HW_curve[i,])
  D0[i] <- item*(1-item)
}
D_n <- diag(D0)
delta <- (1/nodes)*t(HW_curve)%*%HW_curve
mid_term <- (1/nodes)*t(HW_curve)%*%D_n%*%HW_curve
sigma_x0 <- inv(delta)%*%mid_term%*%inv(delta)

sigma1 <- matrix(NA,nrow=100,ncol=100)
sigma2 <- matrix(NA,nrow=100,ncol=100)
for(i in 1:100){
  for(j in 1:100){
    sigma1[i,j] <- pm[i,j]*(1-pm[i,j])
  }
}
sigma2 <- sigma1

#treat A as bernoulli
select_adj <- function(pm){
  ibg <- matrix(NA,nrow=100,ncol=100)
  for(k in 1:nodes){
    for (m in i:nodes){
      ibg[k,m] = sample(c(1,0), 1, prob=c(pm[k,m],1-pm[k,m]))
      ibg[m,k] = ibg[k,m]
    }
  }
  return(ibg)
}
construct_z <- function(mean_A1,mean_A2,sigma1,sigma2){
  Z <- matrix(NA,nrow=nodes,ncol=nodes)
  for (k in 1:nodes){
    for (m in 1:nodes){
      if (k==m){
        Z[k,m] <- sample(c(1/sqrt(nodes),-1/sqrt(nodes)), 1, prob=c(0.5,0.5))
      }
      if (k!=m){
        Z[k,m] <- (mean_A1[k,m]-mean_A2[k,m])/(sqrt(nodes*(sigma1[k,m]/m1 + sigma2[k,m]/m1)))
      }
    }
  }
  return(Z)
}

T1 <- NA
rep_time <- 50
m1 <- 10
for (i in 1:rep_time){
  mean_A1 <- matrix(0,nrow=nodes,ncol=nodes)
  mean_A2 <- matrix(0,nrow=nodes,ncol=nodes)
  Z <- matrix(NA,nrow=nodes,ncol=nodes)
  ibg <- matrix(NA,nrow=nodes,ncol=nodes)
  
  for (j in 1:m1){
    ibg1 <- select_adj(pm)
    ibg2 <- select_adj(pm)
    g1 <- graph_from_adjacency_matrix(ibg1,"undirected")
    A1 = as.matrix(get.adjacency(g1))
    mean_A1 <- mean_A1 + A1
    g2 <- graph_from_adjacency_matrix(ibg2,"undirected")
    A2 = as.matrix(get.adjacency(g2))
    mean_A2 <- mean_A2 + A2
  }
  mean_A1 <- mean_A1/m1
  mean_A2 <- mean_A2/m1
  #construct Z
  Z <- construct_z(mean_A1,mean_A2,sigma1,sigma2)
  T1[i] <- tr(Z%*%Z%*%Z)
}
mean(T1)
var(T1)
#use estimated sigma


