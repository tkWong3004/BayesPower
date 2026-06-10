# BayesPower

**BayesPower** is a Shiny app for sample size determination and power calculation in Bayesian hypothesis testing for common testing problems.

## Available Statistical Tests

-   Standardized mean difference
-   Pearson’s correlation
-   Regression / ANOVA
-   One-proportion and two-proportion tests

## Installation instructions

To install and run **BayesPower**, please follow these steps.

### 1. Install R and RStudio

First, install the latest versions of **R** and **RStudio**.

Windows users also need to install **Rtools**, which is required for building R packages from source. Please install the version of Rtools that matches your version of R.

### 2. Install BayesPower

BayesPower can be installed either from CRAN or from GitHub.

#### Option A: Install the CRAN version

```r
install.packages("BayesPower")
```

#### Option B: Install the development version from GitHub

To install the latest development version from GitHub, first install the `pak` package and then install BayesPower from the GitHub repository:

```r
install.packages("pak")
pak::pak("tkWong3004/BayesPower/package")
```

### 3. Restart RStudio

After installation, close and reopen RStudio.

This step is required for the **RStudio Addin** and the package to appear correctly.

### 4. Launch BayesPower

BayesPower can be launched in two ways.

#### Option A: Launch from the R console

```r
BayesPower::BayesPower_BayesFactor()
```

#### Option B: Launch from the RStudio Addins menu

In RStudio, click **Addins** in the top main toolbar. You should see an entry for **BayesPower Shiny App**.

Clicking this entry launches the BayesPower Shiny app directly in your browser. 

## Citation

**Software paper**\
Wong, T. K., Pawel, S., & Tendeiro, J. (2025). BayesPower: A General Application of Power and Sample Size Calculation for the Bayes Factor. *PsyArXiv.* <https://doi.org/10.31234/osf.io/pgdac_v2>

**Software package**\
Wong, T. K., Pawel, S., & Tendeiro, J. (2025). BayesPower: Sample Size and Power Calculation for Bayesian Testing with Bayes Factor (Version 1.0.1) [R package]. <https://CRAN.R-project.org/package=BayesPower>

## Issues and Feedback

For bug reports, feature requests, or questions, please visit the GitHub Issues page above.
