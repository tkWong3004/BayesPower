// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace boost::math;
using namespace boost::math::quadrature;
using namespace Rcpp;

// Prior density function
double prior_density(double delta, const std::string& model, double location, double scale, double dff) {
  if (model == "Normal") {
    double z = (delta - location) / scale;
    return R::dnorm(z, 0.0, 1.0, false) / scale;
  } else if (model == "Cauchy") {
    double z = (delta - location) / scale;
    return R::dcauchy(z, 0.0, 1.0, false) / scale;
  } else if (model == "t-distribution") {
    non_central_t dist(dff, 0.0);
    return pdf(dist, (delta - location) / scale) / scale;
  }
  return NA_REAL;
}

// Normalize prior over [a, b]
double prior_norm(const std::string& model, double location, double scale, double dff,
                  double a, double b) {
  auto f = [&](double delta) {
    return prior_density(delta, model, location, scale, dff);
  };
  gauss_kronrod<double, 61> integrator;
  return integrator.integrate(f, a, b, 15, 1e-8);
}

// [[Rcpp::export]]
NumericVector t1_BF10_cpp_vec(NumericVector t_vec, double df,
                              std::string model, double location, double scale, double dff,
                              std::string hypothesis) {
  R_xlen_t n = t_vec.size();
  NumericVector result(n);
  
  double a = (hypothesis == ">") ? 0.0 : -INFINITY;
  double b = (hypothesis == "<") ? 0.0 : INFINITY;
  
  double norm_const = (hypothesis == "!=") ? 1.0 : prior_norm(model, location, scale, dff, a, b);
  gauss_kronrod<double, 61> integrator;
  
  for (R_xlen_t i = 0; i < n; ++i) {
    double t = t_vec[i];
    auto integrand = [&](double delta) {
      non_central_t dist(df, delta * std::sqrt(df + 1.0));
      double lik = pdf(dist, t);
      double prior = prior_density(delta, model, location, scale, dff);
      return lik * prior / norm_const;
    };
    
    double marginal = integrator.integrate(integrand, a, b, 15, 1e-8);
    
    non_central_t null_dist(df, 0.0);
    double null = pdf(null_dist, t);
    result[i] = marginal / null;
  }
  
  return result;
}


#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// t1_BF10_cpp_vec
NumericVector t1_BF10_cpp_vec(NumericVector t_vec, double df, std::string model, double location, double scale, double dff, std::string hypothesis);
RcppExport SEXP sourceCpp_10_t1_BF10_cpp_vec(SEXP t_vecSEXP, SEXP dfSEXP, SEXP modelSEXP, SEXP locationSEXP, SEXP scaleSEXP, SEXP dffSEXP, SEXP hypothesisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t_vec(t_vecSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type dff(dffSEXP);
    Rcpp::traits::input_parameter< std::string >::type hypothesis(hypothesisSEXP);
    rcpp_result_gen = Rcpp::wrap(t1_BF10_cpp_vec(t_vec, df, model, location, scale, dff, hypothesis));
    return rcpp_result_gen;
END_RCPP
}
