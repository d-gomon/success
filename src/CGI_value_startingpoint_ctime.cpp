#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


/**
 * @param v - sorted vector instance
 * @param data - value to search
 * @return index of first value larger or equal than data, -1 otherwise
 */
int binary_search_larger_equal(NumericVector v, double data) {
  auto it = std::lower_bound(v.begin(), v.end(), data);
  if (it == v.end()) {
    return -1;
  }  else {
    std::size_t index = std::distance(v.begin(), it);
    return index;
  }
}



/**
 * @param v - sorted vector instance
 * @param data - value to search
 * @return index of first value larger than data, -1 otherwise
 */
int binary_search_larger(NumericVector v, double data) {
  auto it = std::upper_bound(v.begin(), v.end(), data);
  if (it == v.begin()) {
    return -1;
  } else {
    std::size_t index = std::distance(v.begin(), it)-1;
    return index;
  }
}


// [[Rcpp::export]]

NumericVector Active_Individuals_Lambdamat_rows(double helperstime, double ctime, NumericVector entrytimes){
  NumericVector out(2);
  out[0] = binary_search_larger_equal(entrytimes, helperstime);
  out[1] = binary_search_larger(entrytimes, ctime);
  return out;
}















