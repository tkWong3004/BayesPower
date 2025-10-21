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



# Bayes Factor

BF10.t.test.two_sample(
  tval = t_value,
  N1 = 53,
  N2 = 48,
  model = "t-distribution",
  location = 0,
  scale = 0.707,
  dff = 1,
  hypothesis = "!=",
  e = c(-0.36, 0.36)
)
1/ 0.1104477

# variance of d
var_d <- (n1 + n2) / (n1 * n2) + (d^2) / (2 * (n1 + n2))

# standard deviation of d
round(sqrt(var_d),3)

# power analysis for future study
BFpower.t.test_two_sample(
  hypothesis = "!=",
  e = c(-0.36, 0.36),
  interval = "2",
  D = 3,
  target = 0.8,
  alpha = 0.05,
  model = "Normal",
  location = -0.23,
  scale = 0.2,
  dff = 1,
  de_an_prior = 1,
  r = 1,
  mode_bf = 1,
  direct = "h0"
)

#### Example 2 - Correlation
#calculating the Bayes factor

BF10.cor(
  r = 0.393,
  n = 46,
  k = 1,
  alpha = 1,
  beta = 1,
  h0 = 0,
  hypothesis = "!=",
  location = 0,
  scale = 0.01,
  dff = 1,
  model = "d_beta"
)

# if r is not rounded, the bayes factors between ours and Ly et al will be the same

BF10.cor(
  r = 0.3930924,
  n = 46,
  k = 1,
  alpha = 1,
  beta = 1,
  h0 = 0,
  hypothesis = "!=",
  location = 0,
  scale = 0.01,
  dff = 1,
  model = "d_beta"
)

# power analysis
BFpower.cor(
  hypothesis = ">",
  h0 = 0,
  D = 3,
  target = 0.8,
  FP = 0.05,
  model = "d_beta",
  k = 1,
  model_d = "Point",
  location_d = 0.3,
  de_an_prior = 0,
  mode_bf = 1,
  direct = "h1"
)

#### Example 3 - ANOVA
BFpower.f(
  inter = "1",
  D = 3,
  target = 0.8,
  p = 3,
  k = 4,
  model = "effectsize",
  dff = 3,
  rscale = 0.18,
  f_m = 0.1,
  model_d = "Point",
  f_m_d = 0.1,
  de_an_prior = 0,
  mode_bf = 1,
  direct = "h1"
)

#### Example 4 - one proportion
#calculating Bayes factor
BF10.bin.test(
  x = 42,
  n = 52,
  alpha = 1,
  beta = 1,
  location = 0.5,
  scale = 1,
  model = "beta",
  hypothesis = ">"
)
# power analysis
BFpower.bin(
  hypothesis = ">",
  interval = "1",
  D = 3,
  target = 0.8,
  FP = 0.05,
  h0 = 0.5,
  location = 0.5,
  model = "beta",
  alpha = 1,
  beta = 1,
  de_an_prior = 1,
  mode_bf = 1,
  direct = "h1"
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
1/0.07604026
# power analysis
BFpower.props(
  D = 3,
  target = 0.8,
  a0 = 1,
  b0 = 1,
  model1 = "same",
  a1 = 156,
  b1 = 339,
  a2 = 151,
  b2 = 339,
  model2 = "same",
  mode_bf = 1,
  direct = "h1"
)
