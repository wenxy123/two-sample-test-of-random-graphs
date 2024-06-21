// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
mat sample_from_graph(mat A, int s_n, vec sel_nodes){
  mat sel_A(s_n,s_n);
  for (int i=0;i<s_n;i++){
    int row_index = sel_nodes.at(i)-1;
    for (int j=0;j<s_n;j++){
      int col_index = sel_nodes.at(j)-1;
      sel_A.at(i,j) = A.at(row_index,col_index);
    }
  }
  return sel_A;
}

// [[Rcpp::export]]
mat get_mean_adjm(mat A, int nodes, mat sel_nodes){
  mat A11(nodes,nodes,fill::zeros);
  int s_n = nodes;
  int st = sel_nodes.n_rows;
  int ncol = nodes;
  for(int r=0;r<st;r++){
    vec sel_nodes_list(ncol,fill::zeros);
    for(int k=0;k<ncol;k++){
      sel_nodes_list.at(k) = sel_nodes.at(r,k);
    }
    mat A1(nodes,nodes,fill::zeros);
    A1 = sample_from_graph(A, s_n, sel_nodes_list);
    A11 = A11 + A1;
  }
  A11 = A11/st;
  
  return A11;
}

// [[Rcpp::export]]
mat get_embed_adjm(mat A, int nodes, mat sel_nodes,int d){
  mat A11(nodes,nodes,fill::zeros);
  
  int s_n = nodes;
  int st = sel_nodes.n_rows;
  int ncol = nodes;
  for(int r=0;r<st;r++){
    vec sel_nodes_list(ncol,fill::zeros);
    for(int k=0;k<ncol;k++){
      sel_nodes_list.at(k) = sel_nodes.at(r,k);
    }
    mat A1(nodes,nodes,fill::zeros);
    A1 = sample_from_graph(A, s_n, sel_nodes_list);
    vec eigval;
    mat eigvec;
    mat p_hat;
    eig_sym(eigval, eigvec, A1);
    
    mat eigvec1 = eigvec(span::all, span(nodes-1-d,nodes-1));
    mat diag_eigval = diagmat(sqrt(eigval(span(nodes-1-d,nodes-1))));
    
    mat lpvs = eigvec1 * diag_eigval;
    
    p_hat = lpvs * lpvs.t();
    
    A11 = A11+p_hat;
  }
  A11 = A11/st;
  
  return A11;
}

// [[Rcpp::export]]
mat get_sigma(mat A, int sel_n){
  mat sigma(sel_n,sel_n,fill::zeros);
  for(int m=0;m<sel_n;m++){
    for(int n=0;n<sel_n;n++){
      sigma.at(m,n) = A.at(m,n)*(1-A.at(m,n));
    }
  }
  return sigma;
}

// [[Rcpp::export]]
mat construct_Z(mat A1, mat A2, mat sigma1, mat sigma2, int st, int sel_n, vec Z_diag){
  mat Z(sel_n, sel_n, fill::zeros);

  for (int a=0;a<sel_n;a++){
    for(int b=0;b<sel_n;b++){
      if(a==b){
        Z.at(a,b) = Z_diag.at(a);
      }
      else{
        Z.at(a,b) = (A1.at(a,b)-A2.at(a,b))/(sqrt(sel_n*((sigma1.at(a,b)/st) + (sigma2.at(a,b)/st))));
        // Z.at(a,b) = (A1.at(a,b)-A2.at(a,b))/(sqrt(((sigma1.at(a,b)/st) + (sigma2.at(a,b)/st))));
      }
    }
  }
  return Z;
}