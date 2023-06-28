#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sumC(NumericVector x){
  int n = x.size();
  int sum = 0;
  for(int i = 0; i < n; i ++){
    sum += x[i];
  }
  return(sum);
}
