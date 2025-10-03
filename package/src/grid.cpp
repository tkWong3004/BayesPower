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

  int rows = n1 + 1;
  int cols = n2 + 1;
  int total = rows * cols;

  IntegerVector k1(total), k2(total);
  NumericVector log_h1(total), log_h0(total), log_BF10(total);
  LogicalVector PE(total), NE(total);

  double logD = std::log(D);
  double logD_inv = -logD;

  // --- Precompute constants (single-threaded) ---
  double lbeta_a0b0 = R::lbeta(a0, b0);
  double lbeta_a1b1 = R::lbeta(a1, b1);
  double lbeta_a2b2 = R::lbeta(a2, b2);

  // precompute per-k values for group1
  std::vector<double> lchoose_k1(rows), lbeta_k1(rows), d_k1(rows);
  for (int i = 0; i <= n1; ++i) {
    lchoose_k1[i] = R::lchoose(n1, i);
    lbeta_k1[i] = R::lbeta(a1 + i, b1 + n1 - i);
    d_k1[i] = lchoose_k1[i] + lbeta_k1[i] - lbeta_a1b1;
  }

  // precompute per-k values for group2
  std::vector<double> lchoose_k2(cols), lbeta_k2(cols), d_k2(cols);
  for (int j = 0; j <= n2; ++j) {
    lchoose_k2[j] = R::lchoose(n2, j);
    lbeta_k2[j] = R::lbeta(a2 + j, b2 + n2 - j);
    d_k2[j] = lchoose_k2[j] + lbeta_k2[j] - lbeta_a2b2;
  }

  // lbeta sums for the null model
  int max_sum = n1 + n2;
  std::vector<double> lbeta_sum(max_sum + 1);
  for (int s = 0; s <= max_sum; ++s) {
    lbeta_sum[s] = R::lbeta(a0 + s, b0 + n1 + n2 - s);
  }

  // --- Prepare dp-adjusted d_k vectors (single-threaded) ---
  std::vector<double> d_k1_dp(rows), d_k2_dp(cols);

  // model1 adjustment
  if (model1 == "same") {
    for (int i = 0; i <= n1; ++i) d_k1_dp[i] = d_k1[i];
  } else if (model1 == "beta") {
    double lbeta_da1db1 = R::lbeta(da1, db1);
    for (int i = 0; i <= n1; ++i) {
      d_k1_dp[i] = lchoose_k1[i] + R::lbeta(da1 + i, db1 + n1 - i) - lbeta_da1db1;
    }
  } else if (model1 == "Point") {
    double log_dp1 = std::log(dp1);
    double log_1mdp1 = std::log(1.0 - dp1);
    for (int i = 0; i <= n1; ++i) {
      d_k1_dp[i] = lchoose_k1[i] + i * log_dp1 + (n1 - i) * log_1mdp1;
    }
  } else {
    // default fallback: use d_k1
    for (int i = 0; i <= n1; ++i) d_k1_dp[i] = d_k1[i];
  }

  // model2 adjustment
  if (model2 == "same") {
    for (int j = 0; j <= n2; ++j) d_k2_dp[j] = d_k2[j];
  } else if (model2 == "beta") {
    double lbeta_da2db2 = R::lbeta(da2, db2);
    for (int j = 0; j <= n2; ++j) {
      d_k2_dp[j] = lchoose_k2[j] + R::lbeta(da2 + j, db2 + n2 - j) - lbeta_da2db2;
    }
  } else if (model2 == "Point") {
    double log_dp2 = std::log(dp2);
    double log_1mdp2 = std::log(1.0 - dp2);
    for (int j = 0; j <= n2; ++j) {
      d_k2_dp[j] = lchoose_k2[j] + j * log_dp2 + (n2 - j) * log_1mdp2;
    }
  } else {
    for (int j = 0; j <= n2; ++j) d_k2_dp[j] = d_k2[j];
  }

  // --- Main parallel loop: single linear index for better cache behavior ---
  // Use schedule(static) so each thread gets contiguous ranges -> stable timings
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int idx = 0; idx < total; ++idx) {
    int i = idx / cols;       // row index (0..n1)
    int j = idx % cols;       // col index (0..n2)

    k1[idx] = i;
    k2[idx] = j;

    // compute logh1 (original model) and logh0 (null model)
    double logh1 = d_k1[i] + d_k2[j];
    double logh0 = lchoose_k1[i] + lchoose_k2[j] + lbeta_sum[i + j] - lbeta_a0b0;

    log_h1[idx] = logh1;
    log_h0[idx] = logh0;
    log_BF10[idx] = logh1 - logh0;
  }

  // Vectorized logical comparisons (single-threaded and cheap)
  for (int idx = 0; idx < total; ++idx) {
    PE[idx] = log_BF10[idx] > logD;
    NE[idx] = log_BF10[idx] < logD_inv;
  }

  // compute log_h1_dp (uses dp-adjusted d_k vectors)
  NumericVector log_h1_dp(total);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int idx = 0; idx < total; ++idx) {
    int i = idx / cols;
    int j = idx % cols;
    log_h1_dp[idx] = d_k1_dp[i] + d_k2_dp[j];
  }

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
