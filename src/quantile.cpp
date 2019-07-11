#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix get_quantile(NumericMatrix M) {
  int nrow = M.nrow(), ncol = M.ncol();
  // Sort and find ranks:
  NumericMatrix sorted(nrow,ncol);
  IntegerMatrix ranked(nrow,ncol);
  for(int j=0;j<ncol;j++){
    NumericVector sy = M(_,j);
    IntegerVector ry(nrow);
    std::sort(sy.begin(),sy.end());
    sorted(_,j) = sy;
    for(int i=0;i<nrow;i++)
      ry[i] = std::lower_bound(sy.begin(),sy.end(),M(i,j)) - sy.begin();
    ranked(_,j) = ry;
  }
  // Find row means of sorted values:
  NumericVector means(nrow);
  for(int i=0;i<nrow;i++){
    NumericVector x = sorted(i,_);
    means[i] = std::accumulate(x.begin(),x.end(),0.0) / ncol;
  }
  // Sort the row means according to ranks:
  NumericMatrix normed(nrow,ncol);
  for(int i=0;i<nrow;i++)
    for(int j=0;j<ncol;j++)
      normed(i,j) = means[ranked(i,j)];
  return(normed);
}

