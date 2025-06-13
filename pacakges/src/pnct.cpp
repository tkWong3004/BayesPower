// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp>

using namespace boost::math;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pnct(const NumericVector x, const double df, const NumericVector ncp = NumericVector::create(0.0), bool lower = true) {
  
  R_xlen_t nx = x.length();
  R_xlen_t nncp = ncp.length();
  
  NumericVector y(std::max(nx, nncp));
  
  if (nx == 1) {
    // Single x, multiple ncp
    double single_x = x[0];
    for (R_xlen_t i = 0; i < nncp; ++i) {
      non_central_t dist(df, ncp[i]);
      if (lower)
        y[i] = cdf(dist, single_x);
      else
        y[i] = cdf(complement(dist, single_x));
    }
  } else if (nncp == 1) {
    // Multiple x, single ncp
    double single_ncp = ncp[0];
    non_central_t dist(df, single_ncp);
    for (R_xlen_t i = 0; i < nx; ++i) {
      if (lower)
        y[i] = cdf(dist, x[i]);
      else
        y[i] = cdf(complement(dist, x[i]));
    }
  } else if (nx == nncp) {
    // Multiple x and ncp with same length
    for (R_xlen_t i = 0; i < nx; ++i) {
      non_central_t dist(df, ncp[i]);
      if (lower)
        y[i] = cdf(dist, x[i]);
      else
        y[i] = cdf(complement(dist, x[i]));
    }
  } else {
    stop("Lengths of x and ncp must either be the same or one of them must be 1.");
  }
  
  return y;
}
