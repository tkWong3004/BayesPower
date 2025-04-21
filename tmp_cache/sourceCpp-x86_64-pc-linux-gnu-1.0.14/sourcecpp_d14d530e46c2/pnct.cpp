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


#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pnct
NumericVector pnct(const NumericVector x, const double df, const NumericVector ncp, const bool lower);
RcppExport SEXP sourceCpp_3_pnct(SEXP xSEXP, SEXP dfSEXP, SEXP ncpSEXP, SEXP lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type ncp(ncpSEXP);
    Rcpp::traits::input_parameter< const bool >::type lower(lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(pnct(x, df, ncp, lower));
    return rcpp_result_gen;
END_RCPP
}
