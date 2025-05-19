// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp>
#include <cmath>

using namespace boost::math;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pnct(const NumericVector x, const double df,
                   const NumericVector ncp = NumericVector::create(0.0),
                   const bool lower = true) {
  R_xlen_t n = std::max(x.size(), ncp.size());
  NumericVector y(n);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi   = x.size()   == 1 ? x[0]   : x[i];
    double ncpi = ncp.size() == 1 ? ncp[0] : ncp[i];
    
    if (std::isinf(xi)) {
      y[i] = (xi > 0) ? 1.0 : 0.0;
      if (!lower) y[i] = 1.0 - y[i];
      continue;
    }
    
    try {
      non_central_t dist(df, ncpi);
      double val = cdf(dist, xi);
      y[i] = lower ? val : 1.0 - val;
    } catch (...) {
      y[i] = NA_REAL;
    }
  }
  
  return y;
}
