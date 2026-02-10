####
# This R script provides the codes for reproducing the results in Section 5 - Application
library(BayesPower)
#### Example 1 - Standardized mean difference

# given data
x1 <- 0.46; s1 <- 0.17; n1 <- 53
x2 <- 0.50; s2 <- 0.18; n2 <- 48

# pooled standard deviation
s_pooled <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))

# Cohen's d (Hedges' g would add small-sample correction)
d <- (x1 - x2) / s_pooled

# standard error of the difference under equal variances
se_diff_equal <- s_pooled * sqrt(1/n1 + 1/n2)

# t statistic
t_value <- (x1 - x2) / se_diff_equal

# degrees of freedom (equal variance case)
df <- n1 + n2 - 2

# two-sided p-value
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)

# Bayes factor
BF10.ttest.TwoSample(
  tval = -1.148,
  N1 = 53,
  N2 = 48,
  prior_analysis = "t-distribution",
  location = 0,
  scale = 0.707,
  dff = 1,
  alternative = "two.sided"
)

# interval Bayes Factor and equivalence test

BF10.ttest.TwoSample(
  tval = -1.148,
  N1 = 53,
  N2 = 48,
  prior_analysis = "t-distribution",
  location = 0,
  scale = 0.707,
  dff = 1,
  alternative = "two.sided",
  ROPE = c(-0.36,0.36)
)

# variance of d
var_d <- (n1 + n2) / (n1 * n2) + (d^2) / (2 * (n1 + n2))

# standard deviation of d
round(sqrt(var_d),3)

# power analysis for future study
BFpower.ttest.TwoSample(
  alternative = "two.sided",
  ROPE = c(-0.36, 0.36),
  threshold = 3,
  true_rate = 0.8,
  false_rate = 0.05,
  prior_analysis = "Normal",
  location = -0.23,
  scale = 0.2,
  dff = 1,
  type_rate = "negative",
  plot_power = TRUE,
  plot_rel = TRUE,
  r = 1
)

#### Example 2 - Correlation
#calculating the Bayes factor

BF10.cor(
  r = 0.393,
  n = 46,
  prior_analysis = "d_beta",
  k = 1,
  h0 = 0,
  alternative = "two.sided"
)

# if r is not rounded, the bayes factors between ours and Ly et al will be the same
BF10.cor(
  r = 0.3930924,
  n = 46,
  prior_analysis = "d_beta",
  k = 1,
  h0 = 0,
  alternative = "two.sided"
)

# power analysis
BFpower.cor(
  alternative = "greater",
  h0 = 0,
  threshold = 3,
  true_rate = 0.8,
  false_rate = 0.05,
  prior_analysis = "d_beta",
  k = 1,
  prior_design = "Point",
  location_d = 0.3,
  plot_power = TRUE,
  plot_rel = TRUE
)

#### Example 3 - ANOVA
BFpower.f.test(
  threshold = 3,
  true_rate = 0.8,
  false_rate = 0.05,
  p = 3,
  k = 4,
  prior_analysis = "effectsize",
  dff = 3,
  rscale = 0.18,
  f_m = 0.1,
  prior_design = "Point",
  f_m_d = 0.1,
  plot_power = TRUE,
  plot_rel = TRUE
)

#### Example 4 - one proportion
#calculating Bayes factor
BF10.bin.test(
  x = 42,
  n = 52,
  h0 = 0.5,
  prior_analysis = "beta",
  alternative = "greater",
  alpha = 1,
  beta = 1
)
# power analysis
BFpower.bin(
  alternative = "greater",
  threshold = 3,
  true_rate = 0.8,
  false_rate = 0.05,
  h0 = 0.5,
  prior_analysis = "beta",
  alpha = 1,
  beta = 1,
  plot_rel = TRUE,
  plot_power = TRUE
)

#### Example 5 - two proportions
#calculating Bayes factor
BF10.props(
  a0 = 1,
  b0 = 1,
  a1 = 1,
  b1 = 1,
  a2 = 1,
  b2 = 1,
  n1 = 493,
  n2 = 488,
  x1 = 155,
  x2 = 150
)
# power analysis
BFpower.props(
  threshold = 3,
  true_rate = 0.8,
  a0 = 1,
  b0 = 1,
  a1 = 156,
  b1 = 339,
  a2 = 151,
  b2 = 339,
  plot_power = TRUE,
  plot_rel = TRUE
)

