#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
double sum_rcpp(NumericVector log_h1, LogicalVector PE) {
  double total = 0.0;
  
  // Parallelized sum of exp(log_h1[PE])
#pragma omp parallel for reduction(+:total)
  for (int i = 0; i < log_h1.size(); ++i) {
    if (PE[i]) {
      total += std::exp(log_h1[i]);
    }
  }
  
  return total;
}
