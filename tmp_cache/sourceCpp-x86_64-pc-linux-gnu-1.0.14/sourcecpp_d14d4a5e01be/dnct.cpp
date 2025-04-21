// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp>
#include <cmath> // For std::isnan

using namespace boost::math;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dnct(const NumericVector x_raw, const double df,
                   const NumericVector ncp_raw = NumericVector::create(0.0)) {
  NumericVector x   = clone(x_raw);
  NumericVector ncp = clone(ncp_raw);
  
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



#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dnct
NumericVector dnct(const NumericVector x_raw, const double df, const NumericVector ncp_raw);
RcppExport SEXP sourceCpp_1_dnct(SEXP x_rawSEXP, SEXP dfSEXP, SEXP ncp_rawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x_raw(x_rawSEXP);
    Rcpp::traits::input_parameter< const double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type ncp_raw(ncp_rawSEXP);
    rcpp_result_gen = Rcpp::wrap(dnct(x_raw, df, ncp_raw));
    return rcpp_result_gen;
END_RCPP
}
