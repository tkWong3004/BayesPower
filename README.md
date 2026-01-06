# BayesPower

**BayesPower** is a Shiny app for sample size determination and power calculation in Bayesian hypothesis testing for common testing problems.

## Available Statistical Tests

-   Standardized mean difference
-   Pearson’s correlation
-   Regression / ANOVA
-   One-proportion and two-proportion tests

## Usage

BayesPower can be downloaded and launched using the following R commands:

``` r
# CRAN version
install.packages("BayesPower")

# The current developmental version
install.packages("devtools")
devtools::install_github(repo = "tkWong3004/BayesPower", subdir = "package")

# launching the shiny app
BayesPower::BayesPower_BayesFactor()
```

## Citation

**Software paper**\
Wong, T. K., Pawel, S., & Tendeiro, J. (2025). BayesPower: A General Application of Power and Sample Size Calculation for the Bayes Factor. *PsyArXiv.* <https://doi.org/10.31234/osf.io/pgdac_v2>

**Software package**\
Wong, T. K., Pawel, S., & Tendeiro, J. (2025). BayesPower: Sample Size and Power Calculation for Bayesian Testing with Bayes Factor (Version 1.0.1) [R package]. <https://CRAN.R-project.org/package=BayesPower>

## Issues and Feedback

For bug reports, feature requests, or questions, please visit the GitHub Issues page above.
