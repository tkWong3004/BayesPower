// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp>
#include <cmath> // For std::isnan

using namespace boost::math;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dnct(const NumericVector x, const double df,
                   const NumericVector ncp = NumericVector::create(0.0)) {
  R_xlen_t n = std::max(x.size(), ncp.size());
  NumericVector y(n);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi   = x.size()   == 1 ? x[0]   : x[i];
    double ncpi = ncp.size() == 1 ? ncp[0] : ncp[i];
    
    try {
      non_central_t dist(df, ncpi);
      double val = pdf(dist, xi);
      y[i] = std::isnan(val) ? 0.0 : val;
    } catch (...) {
      y[i] = 0.0;
    }
  }
  
  return y;
}
