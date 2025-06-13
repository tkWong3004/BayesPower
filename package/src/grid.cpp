#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
DataFrame BF_grid_rcpp(double D,
                       double a0, double b0,
                       double a1, double b1, int n1,
                       double a2, double b2, int n2,
                       std::string model1, double da1, double db1, double dp1,
                       std::string model2, double da2, double db2, double dp2) {
  
  int total = (n1 + 1) * (n2 + 1);
  IntegerVector k1(total), k2(total);
  NumericVector log_h1(total), log_h0(total), log_BF10(total);
  LogicalVector PE(total), NE(total);
  
  double logD = std::log(D);
  double logD_inv = -logD;
  
  // Precompute constants and beta function values
  double lbeta_a0b0 = R::lbeta(a0, b0);
  double lbeta_a1b1 = R::lbeta(a1, b1);
  double lbeta_a2b2 = R::lbeta(a2, b2);
  
  NumericVector lchoose_k1(n1 + 1), lbeta_k1(n1 + 1), d_k1(n1 + 1);
  for (int i = 0; i <= n1; ++i) {
    lchoose_k1[i] = R::lchoose(n1, i);
    lbeta_k1[i] = R::lbeta(a1 + i, b1 + n1 - i);
    d_k1[i] = lchoose_k1[i] + lbeta_k1[i] - lbeta_a1b1;
  }
  
  NumericVector lchoose_k2(n2 + 1), lbeta_k2(n2 + 1), d_k2(n2 + 1);
  for (int j = 0; j <= n2; ++j) {
    lchoose_k2[j] = R::lchoose(n2, j);
    lbeta_k2[j] = R::lbeta(a2 + j, b2 + n2 - j);
    d_k2[j] = lchoose_k2[j] + lbeta_k2[j] - lbeta_a2b2;
  }
  
  // Precompute lbeta terms for the null model
  int max_sum = n1 + n2;
  NumericVector lbeta_sum(max_sum + 1);
  for (int s = 0; s <= max_sum; ++s) {
    lbeta_sum[s] = R::lbeta(a0 + s, b0 + n1 + n2 - s);
  }
  
  // Parallelize this loop
#pragma omp parallel for collapse(2)
  for (int i = 0; i <= n1; ++i) {
    for (int j = 0; j <= n2; ++j) {
      int idx = i * (n2 + 1) + j;
      k1[idx] = i;
      k2[idx] = j;
      
      // Calculate log(h1) and log(h0)
      double logh1 = d_k1[i] + d_k2[j];
      double logh0 = lchoose_k1[i] + lchoose_k2[j] + lbeta_sum[i + j] - lbeta_a0b0;
      
      log_h1[idx] = logh1;
      log_h0[idx] = logh0;
      log_BF10[idx] = logh1 - logh0;
    }
  }
  
  // Vectorized logical comparisons for PE and NE
  PE = log_BF10 > logD;
  NE = log_BF10 < logD_inv;
  
  // Add dp based on models (adjusting d_k1 and d_k2 as needed)
  NumericVector d_k1_dp(n1 + 1), d_k2_dp(n2 + 1);
  
  // Handle model1 adjustments based on "same"
  if (model1 == "same") {
    d_k1_dp = d_k1; // Use d_k1 as is for "same" model
  } else if (model1 == "beta") {
    for (int i = 0; i <= n1; ++i) {
      d_k1_dp[i] = lchoose_k1[i] + R::lbeta(da1 + i, db1 + n1 - i) - R::lbeta(da1, db1);
    }
  } else if (model1 == "Point") {
    for (int i = 0; i <= n1; ++i) {
      d_k1_dp[i] = lchoose_k1[i] + i * std::log(dp1) + (n1 - i) * std::log(1 - dp1);
    }
  }
  
  // Handle model2 adjustments based on "same"
  if (model2 == "same") {
    d_k2_dp = d_k2; // Use d_k2 as is for "same" model
  } else if (model2 == "beta") {
    for (int j = 0; j <= n2; ++j) {
      d_k2_dp[j] = lchoose_k2[j] + R::lbeta(da2 + j, db2 + n2 - j) - R::lbeta(da2, db2);
    }
  } else if (model2 == "Point") {
    for (int j = 0; j <= n2; ++j) {
      d_k2_dp[j] = lchoose_k2[j] + j * std::log(dp2) + (n2 - j) * std::log(1 - dp2);
    }
  }
  
  // Update log(h1) for dp (if models are not "same")
  NumericVector log_h1_dp(total);
#pragma omp parallel for
  for (int i = 0; i <= n1; ++i) {
    for (int j = 0; j <= n2; ++j) {
      int idx = i * (n2 + 1) + j;
      double logh1_dp = d_k1_dp[i] + d_k2_dp[j];
      log_h1_dp[idx] = logh1_dp;
    }
  }
  
  // Create final result as DataFrame
  return DataFrame::create(
    Named("k1") = k1,
    Named("k2") = k2,
    Named("log_h1") = log_h1,
    Named("log_h0") = log_h0,
    Named("log_BF10") = log_BF10,
    Named("PE") = PE,
    Named("NE") = NE,
    Named("log_h1_dp") = log_h1_dp
  );
}
