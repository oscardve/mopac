/* MoPAC sgRNA sorting internal function
version: 1.0
date: May 2018
description: sorts sgRNAs of each gene by rank
author:  Oscar Villarreal
affiliation: University of Texas MD Anderson Cancer Center. Laboratory of Dr. Han Xu
contact: oscardvillarreal AT gmail.com
*/

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix get_sorted(NumericMatrix M) {
  for(int j=0;j<M.ncol();j++) {
    NumericVector sy = M(_,j);
    std::sort(sy.begin(),sy.end());
    M(_,j) = sy;
  }
  return(M);
}

// [[Rcpp::export]]
NumericMatrix get_sorted1(NumericVector V, int size) {
  int N = V.size() / size;
  for(int j=0;j<N;j++) {
    std::sort(V.begin()+j*size,V.begin()+(j+1)*size);
  }
  NumericMatrix M(size,N);
  for(int j=0;j<N;j++) {
    for(int i=0;i<size;i++) {
      M(i,j) = V(i+j*size);
    }
  }
  return(M);
}
