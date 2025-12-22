#' Sample Size Determination for the One-Sample Bayesian t-Test
#'
#' Perform sample size determination or calculate the probability of obtaining 
#' compelling or misleading evidence for a one-sample Bayesian t-test. 
#' Can handle both point-null and interval-null hypothesis, and allows specifying
#' analysis and design priors.
#'
#' @param alternative Character. The alternative hypothesis being tested: \code{"two-sided"} (default), \code{"greater"}, or \code{"less"}.
#' @param e Optional numeric vector. Bounds for an interval null hypothesis:
#'   - For \code{alternative = "two.sided"}, must be a numeric vector of length 2 with distinct finite values.
#'   - For \code{alternative = "greater"}, must be a single numeric scalar > 0.
#'   - For \code{alternative = "less"}, must be a single numeric scalar < 0.
#' @param prior_analysis Character. Analysis prior under the alternative hypothesis:
#'   \code{"normal"} (default), \code{"moment"} (normal-moment prior), or \code{"t"} (t-Student distribution). 
#' @param location Numeric scalar. Location parameter for the analysis prior under the alternative hypothesis. Default is 0.
#' @param scale Numeric scalar. Scale parameter for the analysis prior under the alternative hypothesis (must be > 0). Default is 0.707.
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior under the alternative hypothesis (required if \code{prior_analysis = "t"}). Default is 1.
#' @param prior_design Optional character. Design prior under the alternative hypothesis:
#'   \code{"normal"}, \code{"moment"}, \code{"t"}, or \code{"point"}. If omitted, the design prior is the same as the analysis prior.
#' @param location_d Numeric scalar. Location parameter for the design prior under the alternative hypothesis. Default is 0. Used only if \code{prior_design} is chosen.
#' @param scale_d Numeric scalar. Scale parameter for the design prior under the alternative hypothesis. Default is 0.707. Used only if \code{prior_design} is chosen.
#' @param dff_d Numeric scalar. Degrees of freedom for the design prior under the alternative hypothesis (required if \code{prior_design = "t"}). Default is 1. Used only if \code{prior_design} is chosen. 
#' @param N Optional numeric scalar. Sample size. Only required if the goal is not sample size determination, but rather to calculate the probability of obtaining compelling or misleading evidence for a given sample size.
#' @param type_rate Character. Either \code{"positive"} (control true/false positive rates) or \code{"negative"} (control true/false negative rates). Default is \code{"positive"}.
#' @param true_rate Numeric scalar. Depending on \code{type_rate}, it is either the targeted true positive or true negative rate (between 0.6 and 0.999). Default is 0.8.
#' @param false_rate Numeric scalar. Depending on \code{type_rate}, it is either the targeted false positive or false negative rate (between 0.001 and 0.1). Default is 0.05.
#' @param D Numeric scalar. Threshold of compelling evidence (must be > 1). Default is 3.
#' @param plot_power Logical. If \code{TRUE}, plots power curve. Default is \code{FALSE}.
#' @param plot_rel Logical. If \code{TRUE}, plots relative likelihood curve. Default is \code{FALSE}.
#' 
#' @details
#' \strong{1. Sample size determination mode (when \code{N = NULL}):}
#'
#' If no sample size is provided, the function determines the minimum sample size required to meet the desired requirements. In this mode, the user must supply the following arguments:
#' \itemize{
#'   \item \code{type_rate} - either \code{"positive"} to control true/false positive rates,
#'         or \code{"negative"} to control true/false negative rates.
#'   \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#'   \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#'   \item \code{D} - the Bayes factor threshold for compelling evidence (must be > 1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability
#' of obtaining compelling evidence meets or exceeds \code{true_rate}, while the
#' probability of misleading evidence does not exceed \code{false_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{N} is supplied):}
#'
#' If a positive numeric sample size \code{N} is provided, the function computes
#' the probabilities of obtaining compelling or misleading evidence for that
#' fixed sample size. In this mode, the arguments \code{type_rate}, \code{true_rate},
#' and \code{false_rate} are ignored; only the Bayes factor threshold \code{D} is used.
#'
#' \strong{Analysis Priors:}
#'
#' The analysis prior specifies the prior distribution of the effect under the
#' alternative hypothesis. The user must provide:
#' \itemize{
#'   \item \code{prior_analysis} - the type of prior: \code{"normal"}, \code{"moment"} (normal-moment prior), or \code{"t"}.
#'   \item \code{location} - the mean or location of the prior.
#'   \item \code{scale} - the standard deviation or scale (must be positive).
#'   \item \code{dff} - degrees of freedom (required if \code{prior_analysis = "t"}).
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' A design prior can be supplied to reflect uncertainty about the effect size
#' during study planning. If provided, the following must be supplied:
#' \itemize{
#'   \item \code{prior_design} - the type of design prior: \code{"normal"}, \code{"moment"}, \code{"t"}, or \code{"point"}.
#'   \item \code{location_d} - the location of the design prior.
#'   \item \code{scale_d} - the scale parameter (positive for all models except \code{"point"}).
#'   \item \code{dff_d} - degrees of freedom for \code{"t"} design priors.
#' }
#'
#' \strong{Interval Null Hypothesis:}
#'
#' The argument \code{e} specifies the bounds of an interval null hypothesis.
#' If \code{e} is provided, the function evaluates the Bayes factor for an interval
#' null hypothesis. For a point-null hypothesis, \code{e} should be left as \code{NULL}.
#' 
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling
#' evidence as a function of sample size. 
#' If \code{plot_rel = TRUE}, a relative
#' likelihood curve is plotted to illustrate how evidence varies with effect size.
#' 
#' @return A list of class \code{BFpower_t} containing 10 elements:
#' \describe{
#'   \item{type}{Character string describing the test type, always "One-sample t-test".}
#'   \item{alternative}{Character, the tested alternative hypothesis.}
#'   \item{e}{Optional numeric vector for interval null bounds.}
#'   \item{analysis_h1}{List with the analysis prior parameters: \code{prior_analysis}, \code{location}, \code{scale}, and optionally \code{dff}.}
#'   \item{design_h1}{List with the design prior parameters: \code{prior_design}, \code{location}, \code{scale}, and optionally \code{dff} (or \code{NULL} if not provided).}
#'   \item{results}{Data frame with the probabilities of compelling/misleading evidence (\code{NaN} if calculation fails), and with the required sample size.}
#'   \item{D}{Numeric, threshold of compelling evidence.}
#'   \item{plot_power}{Logical, whether to plot the power curve.}
#'   \item{plot_rel}{Logical, whether the relationship between the BF and t-value is plotted.}
#' }
#' 
#' @examples
#' # Sample size determination with a point-null hypothesis:
#' BFpower.ttest.OneSample(
#'   alternative    = "two.sided",
#'   prior_analysis = "t",
#'   location       = 0,
#'   scale          = 0.707,
#'   dff            = 1,
#'   type_rate      = "positive",
#'   true_rate      = 0.8,
#'   false_rate     = 0.05,
#'   D              = 3
#' )
#' @export
BFpower.ttest.OneSample <- function(
    alternative = "two.sided", 
    e = NULL,
    prior_analysis = "normal", location   = 0, scale   = 0.707, dff   = 1,
    prior_design   = NULL,     location_d = 0, scale_d = 0.707, dff_d = 1,
    N = NULL,
    type_rate = "positive", true_rate = 0.8, false_rate = 0.05, 
    D = 3, 
    plot_power = FALSE, plot_rel = FALSE) 
{
  # Mode:
  if (is.null(N)) mode_bf <- 1 else mode_bf <- 0
  
  # Sample size:
  if (mode_bf == 0) {
    # Check that N is a positive numeric scalar #!#!# integer?
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 0) {
      stop("Argument [N] (sample size) must be a positive numeric scalar. Stop.")
    }
  } else {N = 2} #!#!# why?
  
  # Alternative hypothesis:
  if (alternative %in% c("less", "two.sided", "greater") == FALSE){
    stop("Argument [alternative] (the tested alternative hypothesis) must be either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test). Stop.")
  }
  
  # Equivalence test or not:
  interval <- if (is.null(e)) 1 else 0
  
  if (!is.null(e)) {
    
    if (alternative == "two.sided") {
      # e must be a numeric vector of length 2
      if (!is.numeric(e) || length(e) != 2 || any(!is.finite(e)) || e[1] == e[2]) {
        stop("For alternative hypothesis '!=', argument [e] (specifying the bounds for an interval null hypothesis) must be a numeric vector of length 2 with two distinct finite values. Stop.")
      }
    }
    
    if (alternative == "greater") {
      # e must be a numeric scalar > 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e <= 0) {
        stop("For alternative hypothesis '>', argument [e] (specifying the bounds for an interval null hypothesis) must be a numeric scalar > 0. Stop.")
      }
    }
    
    if (alternative == "less") {
      # e must be a numeric scalar < 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e >= 0) {
        stop("For alternative hypothesis '<', argument [e] (specifying the bounds for an interval null hypothesis) must be a numeric scalar < 0. Stop.")
      }
    }
    
  }
  
  # Analysis prior:
  # if (missing(prior_analysis)) {
  #   stop("Argument [prior_analysis] (the analysis prior under the alternative hypothesis) must be either `normal`, `moment` (normal-moment prior) or `t` (t-distribution). Stop.")
  # }
  if (prior_analysis %in% c("normal", "moment", "t") == FALSE){
    stop("Argument [prior_analysis] (the analysis prior under the alternative hypothesis) must be either `normal`,  `moment` (normal-moment prior) or `t` (t-distribution). Stop.")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("Argument [location] for the analysis prior must be a numeric scalar. Stop.")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale<=0) {
    stop("Argument [scale] for the analysis prior must be a positive numeric scalar. Stop.")
  }
  if (prior_analysis == "t") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("Argument [dff] (the degrees of freedom of the analysis prior under the alternative hypothesis) must be a positive numeric scalar when prior_analysis = 't'. Stop.")
    }
  } else {
    dff = 0
  }
  
  # Design prior:
  if (!is.null(prior_design)) {
    
    de_an_prior <- 0
    
    if (!(prior_design %in% c("normal", "moment", "t", "point"))) {
      stop("Argument [prior_design] (the design prior under the alternative hypothesis) must be either `normal`,  `moment` (normal-moment prior), `t` (t-distribution), or `point`. Stop.")
    }
    if (!is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d))
      stop("Argument [location_d] for the design prior must be a numeric scalar. Stop.")
    if (prior_design %in% c("normal", "moment", "t")) {
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0)
        stop("Argument [scale_d] for the design prior must be a positive numeric scalar. Stop.")
    }
    if (prior_design == "t") {
      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0)
        stop("Argument [dff_d] (the degrees of freedom of the design prior under the alternative hypothesis) must be a positive numeric scalar when prior_design = 't'. Stop.")
    } else {
      dff_d <- 0 #!#!# why?
    }
    
  } else {
    de_an_prior <- 1
  }
  
  # Desired power and strength of evidence:
  if (mode_bf == 1) {
    if (!(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate). Stop.")
    }
    direct <- switch (type_rate,            #!#!# why?
                      "positive" = "h1",
                      "negative" = "h0")
    
    if (!is.numeric(true_rate) || length(true_rate) != 1 || !is.finite(true_rate) || 
        true_rate <= 0.6 || true_rate >= 0.999){
      stop("Argument [true_rate] (targeted true positive or true negative rate) must be a numeric scalar strictly greater than 0.6 and smaller than 0.999. Stop.")
    }
    target <- true_rate
    
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("Argument [false_rate] (targeted false positive or false negative rate) must be a numeric scalar strictly greater than 0.001 and smaller than 0.1. Stop.")
    }
    alpha <- false_rate
    
    if (!is.numeric(D) || length(D) != 1 || !is.finite(D) || D <= 1) {
      stop("Argument [D] (threshold of compelling evidence) must be a numeric scalar greater than 1. Stop.")
    }
  } else {
    target <- 0
    alpha  <- 0
  }
  
  
  # Call appropriate table function with error handling:
  tryCatch(
    {
      if (interval == 1) {
        results = suppressWarnings(t1_Table(D, target, prior_analysis, location, scale, dff, alternative,
                                            prior_design, location_d, scale_d, dff_d, de_an_prior, N, mode_bf, alpha, direct))
      } else {
        results = suppressWarnings(t1e_table(D,target,prior_analysis,location,scale,dff, alternative,e ,
                                             prior_design,scale_d,dff_d, de_an_prior,N,mode_bf,location_d ,alpha,direct ))
      }
      
    },
    error = function(err) {
      message("Required sample size > 10,000")
      stop(NaN)
    }
  )
  
  analysis_h1 <- list(
    prior_analysis = prior_analysis,
    location       = location,
    scale          = scale
  )
  if (prior_analysis == "t") {
    analysis_h1$dff <- dff
  }
  
  if (!is.null(prior_design)) {
    
    # Base fields always included
    design_h1 <- list(
      prior_design = prior_design,
      location     = location_d,
      scale        = scale_d
    )
    
    # Only add dff if model is t-distribution
    if (prior_design == "t") {
      design_h1$dff <- dff_d
    }
    
  } else {
    # prior_design is NULL > fill all fields with NULL
    design_h1 <- list(
      prior_design = NULL,
      location     = NULL,
      scale        = NULL,
      dff          = NULL
    )
  }
  
  object <- list(
    type        = "One-sample t-test",
    alternative = alternative,
    e           = e,
    analysis_h1 = analysis_h1,
    design_h1   = design_h1,
    results     = results,
    D           = D,
    mode_bf     = mode_bf,
    plot_power  = plot_power,
    plot_rel    = plot_rel
  )
  class(object) <- "BFpower_t"
  plot(object)
  return(object)
}



#' Sample Size Determination for the Two-Sample Bayesian t-Test
#'
#' Perform sample size determination or calculate the probability of obtaining
#' compelling or misleading evidence for a two-sample Bayesian t-test.
#' Can handle both point-null and interval-null hypothesis, and allows specifying
#' analysis and design priors.
#'
#' @param alternative Character. The alternative hypothesis being tested: \code{"two-sided"} (default), \code{"greater"}, or \code{"less"}.
#' @param e Optional numeric vector. Bounds for an interval null hypothesis:
#'   - For \code{alternative = "two.sided"}, must be a numeric vector of length 2 with distinct finite values.
#'   - For \code{alternative = "greater"}, must be a single numeric scalar > 0.
#'   - For \code{alternative = "less"}, must be a single numeric scalar < 0.
#' @param prior_analysis Character. Analysis prior under the alternative hypothesis:
#'   \code{"normal"} (default), \code{"moment"} (normal-moment prior), or \code{"t"} (t-Student distribution). 
#' @param location Numeric scalar. Location parameter for the analysis prior under the alternative hypothesis. Default is 0.
#' @param scale Numeric scalar. Scale parameter for the analysis prior under the alternative hypothesis (must be > 0). Default is 0.707.
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior under the alternative hypothesis (required if \code{prior_analysis = "t"}). Default is 1.
#' @param prior_design Optional character. Design prior under the alternative hypothesis:
#'   \code{"normal"}, \code{"moment"}, \code{"t"}, or \code{"point"}. If omitted, the design prior is the same as the analysis prior.
#' @param location_d Numeric scalar. Location parameter for the design prior under the alternative hypothesis. Default is 0. Used only if \code{prior_design} is chosen.
#' @param scale_d Numeric scalar. Scale parameter for the design prior under the alternative hypothesis. Default is 0.707. Used only if \code{prior_design} is chosen.
#' @param dff_d Numeric scalar. Degrees of freedom for the design prior under the alternative hypothesis (required if \code{prior_design = "t"}). Default is 1. Used only if \code{prior_design} is chosen.
#' @param N1 Optional numeric scalar. Sample size for group 1 (used if \code{r = NULL}).
#' @param N2 Optional numeric scalar. Sample size for group 2 (used if \code{r = NULL}).
#' @param r Optional numeric scalar. Ratio of sample size \code{N2 / N1} (used if \code{N1} and \code{N2} are \code{NULL}).
#' @param type_rate Character. Either \code{"positive"} (control true/false positive rates) or \code{"negative"} (control true/false negative rates). Default is \code{"positive"}.
#' @param true_rate Numeric scalar. Depending on \code{type_rate}, it is either the targeted true positive or true negative rate (between 0.6 and 0.999). Default is 0.8.
#' @param false_rate Numeric scalar. Depending on \code{type_rate}, it is either the targeted false positive or false negative rate (between 0.001 and 0.1). Default is 0.05.
#' @param D Numeric scalar. Threshold of compelling evidence (must be > 1). Default is 3.
#' @param plot_power Logical. If \code{TRUE}, plots power curve. Default is \code{FALSE}.
#' @param plot_rel Logical. If \code{TRUE}, plots relative likelihood curve. Default is \code{FALSE}.
#'
#' @details
#' \strong{1. Sample size determination mode (when \code{N1 = NULL} and \code{N2 = NULL}, but \code{r} is provided):}
#'
#' If no sample sizes are provided, the function calculates the minimum sample size required for both groups to meet the desired requirements. In this mode, the user must supply the following arguments:
#' \itemize{
#'   \item \code{type_rate} - either \code{"positive"} to control true/false positive rates,
#'         or \code{"negative"} to control true/false negative rates.
#'   \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#'   \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#'   \item \code{D} - the Bayes factor threshold for compelling evidence (must be > 1).
#'   \item \code{r} - the allocation ratio of group 2 to group 1 sample sizes (\code{N2/N1}).
#' }
#'
#' The function iteratively finds the smallest sample sizes \code{N1} and \code{N2 = r * N1} for which the probability
#' of obtaining compelling evidence meets or exceeds \code{true_rate}, while the probability of misleading evidence
#' does not exceed \code{false_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{N1} and \code{N2} are supplied):}
#'
#' If positive numeric sample sizes \code{N1} and \code{N2} are provided, the function computes
#' the probabilities of obtaining compelling or misleading evidence for those fixed sample sizes. In this mode,
#' the arguments \code{type_rate}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{D} is used.
#'
#' \strong{Analysis Priors:}
#'
#' The analysis prior specifies the prior distribution of the effect under the alternative hypothesis. The user must provide:
#' \itemize{
#'   \item \code{prior_analysis} - the type of prior: \code{"normal"}, \code{"moment"} (normal-moment prior), or \code{"t"}.
#'   \item \code{location} - the mean or location of the prior.
#'   \item \code{scale} - the standard deviation or scale (must be positive).
#'   \item \code{dff} - degrees of freedom (required if \code{prior_analysis = "t"}).
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' A design prior can be supplied to reflect uncertainty about the effect size during study planning. If provided, the following must be supplied:
#' \itemize{
#'   \item \code{prior_design} - the type of design prior: \code{"normal"}, \code{"moment"}, \code{"t"}, or \code{"point"}.
#'   \item \code{location_d} - the location of the design prior.
#'   \item \code{scale_d} - the scale parameter (positive for all models except \code{"point"}).
#'   \item \code{dff_d} - degrees of freedom for \code{"t"} design priors.
#' }
#'
#' \strong{Interval Null Hypothesis:}
#'
#' The argument \code{e} specifies the bounds of an interval null hypothesis.
#' If \code{e} is provided, the function evaluates the Bayes factor for an interval
#' null hypothesis. For a point-null hypothesis, \code{e} should be left as \code{NULL}.
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling
#' evidence as a function of the sample sizes. If \code{plot_rel = TRUE}, a relative
#' likelihood curve is plotted to illustrate how evidence varies with effect size.
#'
#' @return An object of class \code{BFpower_t} containing 10 elements:
#' \describe{
#'   \item{type}{Character string describing the test type, always "Indepedent-samples t-test (equal variance)".}
##'   \item{alternative}{Character, the tested alternative hypothesis.}
#'   \item{e}{Optional numeric vector for interval null bounds.}
#'   \item{analysis_h1}{List with the analysis prior parameters: \code{prior_analysis}, \code{location}, \code{scale}, and optionally \code{dff}.}
#'   \item{design_h1}{List with the design prior parameters: \code{prior_design}, \code{location}, \code{scale}, and optionally \code{dff} (or \code{NULL} if not provided).}
#'   \item{results}{Data frame with the probabilities of compelling/misleading evidence (\code{NaN} if calculation fails), and with the required sample size.}
#'   \item{D}{Numeric, threshold of compelling evidence.}
#'   \item{plot_power}{Logical, whether to plot the power curve.}
#'   \item{plot_rel}{Logical, whether the relationship between the BF and t-value is plotted.}
#' }
#'
#' @examples
#' # Sample size determination with a point-null hypothesis:
#' BFpower.ttest.TwoSamples(
#'   alternative    = "two.sided",
#'   prior_analysis = "t",
#'   location       = 0,
#'   scale          = 0.707,
#'   dff            = 1,
#'   type_rate      = "positive",
#'   true_rate      = 0.8,
#'   false_rate     = 0.05,
#'   D              = 3, 
#'   r              = 1
#' )
#'
#' # Probability of obtaining compelling / misleading evidence using an interval null hypothesis:
#' BFpower.ttest.TwoSamples(
#'   alternative    = "greater",
#'   e              = 0.2,
#'   prior_analysis = "normal",
#'   location       = 0,
#'   scale          = 1,
#'   dff            = 0,
#'   N1             = 30,
#'   N2             = 25, 
#'   type_rate      = "positive",
#'   true_rate      = 0.8,
#'   false_rate     = 0.05,
#'   D              = 3
#' )
#'
#' @export
BFpower.ttest.TwoSamples <- function(
    alternative = "two.sided", 
    e = NULL,
    prior_analysis = "normal", location   = 0, scale   = 0.707, dff   = 1,
    prior_design   = NULL,     location_d = 0, scale_d = 0.707, dff_d = 1,
    N1 = NULL, N2 = NULL, r = NULL,
    type_rate = "positive", true_rate = 0.8, false_rate = 0.05, 
    D = 3, 
    plot_power = FALSE, plot_rel = FALSE) 
{
  # Checking N1, N2, r consistency:
  if (is.null(N1) && is.null(N2) && is.null(r)) {
    stop(
      "Argument [r] (ratio N2/N1) must be specified for sample size calculation.\n",
      "or\n", 
      "Arguments [N1], [N2] (groups sample sizes) must be specified for power calculation.\n", 
      "Stop."
    )
  }
  
  # If r provided -> N1 and N2 must be NULL:
  if (!is.null(r)) {
    mode_bf <- 1
    
    # r must be numeric scalar > 0
    if (!is.numeric(r) || length(r) != 1 || !is.finite(r) || r <= 0) {
      stop("Argument [r] (ratio N2/N1) must be a positive numeric scalar. Stop.")
    }
    
    if (!is.null(N1) || !is.null(N2)) {
      stop("If argument [r] is provided, both N1 and N2 must be NULL for sample size determination. Stop.")
    }
  }
  
  # If r is NULL -> N1 and N2 must both be valid numeric scalars:
  if (is.null(r)) {
    mode_bf <- 0
    if (is.null(N1) || is.null(N2)) {
      stop("If [r] (ratio N2/N1) is not provided, both N1 and N2 must be provided. Stop.")
    }
    
    if (!is.numeric(N1) || length(N1) != 1 || !is.finite(N1) || N1 <= 0) {
      stop("Argument [N1] (sample size for group 1) must be a positive numeric scalar. Stop")
    }
    if (!is.numeric(N2) || length(N2) != 1 || !is.finite(N2) || N2 <= 0) {
      stop("Argument [N2] (sample size for group 2) must be a positive numeric scalar. Stop")
    }
  }
  
  # Alternative hypothesis:
  if (alternative %in% c("less", "two.sided", "greater") == FALSE){
    stop("Argument [alternative] (the tested alternative hypothesis) must be either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test). Stop.")
  }
  
  # Equivalence test or not:
  interval <- if (is.null(e)) 1 else 0
  
  if (!is.null(e)) {
    
    if (alternative == "two.sided") {
      # e must be a numeric vector of length 2
      if (!is.numeric(e) || length(e) != 2 || any(!is.finite(e)) || e[1] == e[2]) {
        stop("For alternative hypothesis '!=', argument [e] (specifying the bounds for an interval null hypothesis) must be a numeric vector of length 2 with two distinct finite values. Stop.")
      }
    }
    
    if (alternative == "greater") {
      # e must be a numeric scalar > 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e <= 0) {
        stop("For alternative hypothesis '>', argument [e] (specifying the bounds for an interval null hypothesis) must be a numeric scalar > 0. Stop.")
      }
    }
    
    if (alternative == "less") {
      # e must be a numeric scalar < 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e >= 0) {
        stop("For alternative hypothesis '<', argument [e] (specifying the bounds for an interval null hypothesis) must be a numeric scalar < 0. Stop.")
      }
    }
    
  }
  
  # Analysis prior:
  # if (missing(prior_analysis)) {
  #   stop("Argument [prior_analysis] (the analysis prior under the alternative hypothesis) must be either `normal`, `moment` (normal-moment prior) or `t` (t-distribution). Stop.")
  # }
  if (prior_analysis %in% c("normal", "moment", "t") == FALSE){
    stop("Argument [prior_analysis] (the analysis prior under the alternative hypothesis) must be either `normal`,  `moment` (normal-moment prior) or `t` (t-distribution). Stop.")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("Argument [location] for the analysis prior must be a numeric scalar. Stop.")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale<=0) {
    stop("Argument [scale] for the analysis prior must be a positive numeric scalar. Stop.")
  }
  if (prior_analysis == "t") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("Argument [dff] (the degrees of freedom of the analysis prior under the alternative hypothesis) must be a positive numeric scalar when prior_analysis = 't'. Stop.")
    }
  } else {
    dff = 0
  }
  
  # Design prior:
  if (!is.null(prior_design)) {
    
    de_an_prior <- 0
    
    if (!(prior_design %in% c("normal", "moment", "t", "point"))) {
      stop("Argument [prior_design] (the design prior under the alternative hypothesis) must be either `normal`,  `moment` (normal-moment prior), `t` (t-distribution), or `point`. Stop.")
    }
    if (!is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d))
      stop("Argument [location_d] for the design prior must be a numeric scalar. Stop.")
    if (prior_design %in% c("normal", "moment", "t")) {
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0)
        stop("Argument [scale_d] for the design prior must be a positive numeric scalar. Stop.")
    }
    if (prior_design == "t") {
      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0)
        stop("Argument [dff_d] (the degrees of freedom of the design prior under the alternative hypothesis) must be a positive numeric scalar when prior_design = 't'. Stop.")
    } else {
      dff_d <- 0 #!#!# why?
    }
    
  } else {
    de_an_prior <- 1
  }
  
  # Desired power and strength of evidence:
  if (mode_bf == 1) {
    if (!(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate). Stop.")
    }
    direct <- switch (type_rate,            #!#!# why?
                      "positive" = "h1",
                      "negative" = "h0")
    
    if (!is.numeric(true_rate) || length(true_rate) != 1 || !is.finite(true_rate) || 
        true_rate <= 0.6 || true_rate >= 0.999){
      stop("Argument [true_rate] (targeted true positive or true negative rate) must be a numeric scalar strictly greater than 0.6 and smaller than 0.999. Stop.")
    }
    target <- true_rate
    
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("Argument [false_rate] (targeted false positive or false negative rate) must be a numeric scalar strictly greater than 0.001 and smaller than 0.1. Stop.")
    }
    alpha <- false_rate
    
    if (!is.numeric(D) || length(D) != 1 || !is.finite(D) || D <= 1) {
      stop("Argument [D] (threshold of compelling evidence) must be a numeric scalar greater than 1. Stop.")
    }
  } else {
    target <- 0
    alpha  <- 0
  }
  
  # Call appropriate table function with error handling:
  tryCatch(
    if (interval == 1) {
      results = suppressWarnings(t2_Table(D, r, target, prior_analysis, location, scale, dff, alternative,
                                          prior_design, location_d, scale_d, dff_d, de_an_prior, N1, N2, mode_bf, alpha, direct))
    } else {
      results = suppressWarnings(t2e_table(D, r, target, prior_analysis,location, scale, dff, alternative, e,
                                           prior_design,location_d, scale_d, dff_d, de_an_prior, N1, N2, mode_bf, alpha, direct))
    },
    error = function(err) {
      message("Required sample size > 10,000")
      stop(NaN)
    }
  )
  
  analysis_h1 <- list(
    prior_analysis = prior_analysis,
    location       = location,
    scale          = scale
  )
  if (prior_analysis == "t") {
    analysis_h1$dff <- dff
  }
  
  if (!is.null(prior_design)) {
    
    # Base fields always included
    design_h1 <- list(
      prior_design = prior_design,
      location     = location_d,
      scale        = scale_d
    )
    
    # Only add dff if model is t-distribution
    if (prior_design == "t") {
      design_h1$dff <- dff_d
    }
    
  } else {
    # prior_design is NULL > fill all fields with NULL
    design_h1 <- list(
      prior_design = NULL,
      location     = NULL,
      scale        = NULL,
      dff          = NULL
    )
  }
  
  object <- list(
    type        = "Indepedent-samples t-test (equal variance)",
    alternative = alternative,
    e           = e,
    analysis_h1 = analysis_h1,
    design_h1   = design_h1,
    results     = results,
    D           = D,
    mode_bf     = mode_bf,
    plot_power  = plot_power,
    plot_rel    = plot_rel
  )
  class(object) <- "BFpower_t"
  plot(object)
  return(object)
}



#' Sample size determination for Bayesian correlation test
#'
#' Perform sample size determination or calculate probabilities of compelling and misleading evidence
#' for a Bayesian correlation test.
#'
#' @param hypothesis The hypothesis being tested: two-sided (\code{"!="}), right-sided (\code{">"}), or left-sided (\code{"<"}).
#' @param h0 NaN value of the correlation; must be a numeric scalar between -0.8 and 0.8.
#' @param e Optional numeric vector or scalar specifying bounds for an interval null; used if interval Bayes factor is calculated.
#' @param D Threshold of compelling evidence (numeric scalar > 1).
#' @param true_rate Targeted true positive rate (if \code{positive = "positive"}) or true negative rate (if \code{positive = "negative"}).
#' @param false_rate Targeted false positive rate (if \code{positive = "positive"}) or false negative rate (if \code{positive = "negative"}).
#' @param model Analysis prior model under the alternative hypothesis: default beta (\code{"d_beta"}), beta (\code{"beta"}), or normal moment (\code{"moment"}).
#' @param k Parameter for the default beta prior (\code{"d_beta"}).
#' @param alpha Parameter for the beta prior (\code{"beta"}).
#' @param beta Parameter for the beta prior (\code{"beta"}).
#' @param scale Scale parameter for the normal moment prior (\code{"moment"}).
#' @param model_d Design prior model under the alternative hypothesis: default beta (\code{"d_beta"}), beta (\code{"beta"}), normal moment (\code{"moment"}), or point (\code{"point"}).
#' @param alpha_d Parameter for the design beta prior (\code{"beta"}).
#' @param beta_d Parameter for the design beta prior (\code{"beta"}).
#' @param location_d Location parameter for the design point prior (\code{"point"}).
#' @param k_d Parameter for the design default beta prior (\code{"d_beta"}).
#' @param scale_d Scale parameter for the design normal moment prior (\code{"moment"}).
#' @param N Sample size (numeric scalar).
#' @param positive Character indicating which rate to control: \code{"positive"} (true/false positive rates) or \code{"negative"} (true/false negative rates).
#' @param plot_power Logical; if TRUE, plots power curves.
#' @param plot_rel Logical; if TRUE, plots the relationship between the BF and data.
#'
#'@details
#' \strong{1. Sample size determination mode (when \code{N = NaN}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size. The user must provide:
#' \itemize{
#' \item \code{positive} - either \code{"positive"} to control true/false positive rates, or \code{"negative"} to control true/false negative rates.
#' \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#' \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#' \item \code{D} - the Bayes factor threshold for compelling evidence (must be > 1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence meets or exceeds \code{true_rate}, while the probability of misleading evidence does not exceed \code{false_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{N} is supplied):}
#'
#' If a positive numeric sample size \code{N} is provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, the arguments \code{positive}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{D} is used.
#'
#' \strong{Hypothesis specification:}
#'
#' The \code{hypothesis} argument defines the test type: \code{"!="} for two-sided, \code{">"} for right-sided, or \code{"<"} for left-sided tests. The optional \code{e} argument specifies bounds for an interval null hypothesis. If \code{e = NULL}, a point-NaN test is assumed.
#'
#' \strong{Analysis Priors:}
#'
#' The analysis prior specifies the prior distribution of the correlation under the alternative hypothesis. Depending on \code{model}, the user must supply:
#' \itemize{
#' \item \code{d_beta} (default beta): \code{k} > 0.
#' \item \code{beta} (stretched beta): \code{alpha} and \code{beta} > 0.
#' \item \code{Moment} (normal-moment prior): \code{scale} > 0.
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' A design prior can be specified for planning purposes. If provided, \code{model_d} must be one of \code{"d_beta"}, \code{"beta"}, \code{"moment"}, or \code{"point"}, and the corresponding parameters must be supplied:
#' \itemize{
#' \item \code{d_beta}: \code{k_d} > 0.
#' \item \code{beta}: \code{alpha_d} and \code{beta_d} > 0.
#' \item \code{Moment}: \code{scale_d} > 0.
#' \item \code{Point}: \code{location_d} numeric scalar.
#' }
#'
#' \strong{Interval null Hypothesis:}
#'
#' If \code{e} is provided, the function evaluates the Bayes factor for an interval null. Otherwise, a point-null hypothesis is assumed.
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling evidence as a function of sample size. If \code{plot_rel = TRUE}, the relationship between the BF and correlation is plotted.
#'
#' @return A list of class \code{BFpower_r} containing:
#' \itemize{
#'   \item \code{type}: Test type (always "Correlation").
#'   \item \code{hypothesis}: Hypothesis tested.
#'   \item \code{h0}: NaN correlation value.
#'   \item \code{e}: Bounds for interval null (if used).
#'   \item \code{analysis_h1}: Analysis prior parameters under the alternative hypothesis.
#'   \item \code{design_h1}: Design prior parameters under the alternative hypothesis.
#'   \item \code{results}: Data frame of probabilities of compelling and misleading evidence.
#'   \item \code{D}: Threshold of compelling evidence.
#'   \item \code{plot_power}: Logical, whether power curves are plotted.
#'   \item \code{plot_rel}: Logical, whether the relationship between the BF and correlation is plotted.
#' }
#'
#' @examples
#' BFpower.cor(
#'   hypothesis = "!=",
#'   h0 = 0,
#'   D = 3,
#'   true_rate = 0.8,
#'   false_rate = 0.05,
#'   model = "d_beta",
#'   k = 1,
#'   positive = "positive"
#' )
#'
#' @export
BFpower.cor<- function(hypothesis , h0, e = NULL,
                       D , true_rate, false_rate ,
                       model , k , alpha , beta , scale ,
                       model_d = NaN, alpha_d , beta_d , location_d ,
                       k_d , scale_d ,
                       N = NaN,  positive="positive",plot_power=FALSE,plot_rel=FALSE) {
  # mode
  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < -0.8 || h0 > 0.8) {
    stop("argument [h0] NaN value of rho must be a single numeric scalar between -0.8 and 0.8")
  }

  location <- h0
  dff <- dff_d <- 1
  if ( is.nan(N)) mode_bf=1 else mode_bf = 0


  # sample size
  if (mode_bf == 0) {
    # Check that N is a positive numeric scalar
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 0) {
      stop("argument [N] sample size must be a positive numeric scalar when mode_bf = 0")
    }
  }else {N=3}


  # hypothesis
  if(hypothesis %in% c("<", "!=", ">") == FALSE){
    stop("argument [hypothesis] should be set to either `<`  (left-sided test),  `!=` (two-sided test) or `>` (right-sided test)")
  }

  # Equivlance test or not
  interval <- if (is.null(e)) 1 else 0

  if (!is.null(e)) {

    if (hypothesis == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(e) || length(e) != 2 || any(!is.finite(e)) || e[1] == e[2]) {
        stop("For hypothesis '!=', argument [e] must be a numeric vector of length 2 with two distinct finite values")
      }
      # Additional bounds checks
      if (min(e) < -0.5 || max(e) > 0.5) {
        stop("For hypothesis '!=', e must satisfy min(e) >= -0.5 and max(e) <= 0.5")
      }
      if ((h0 + min(e)) <= -1 || (h0 + min(e)) >= 1) {
        stop("For hypothesis '!=', h0 + min(e) must be between -1 and 1")
      }

    } else if (hypothesis == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e <= 0) {
        stop("For hypothesis '>', argument [e] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (e > 0.5) stop("For hypothesis '>', e must be <= 0.5")
      if ((h0 + e) >= 1) stop("For hypothesis '>', h0 + e must be < 1")

    } else if (hypothesis == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e >= 0) {
        stop("For hypothesis '<', argument [e] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (e < -0.5) stop("For hypothesis '<', e must be >= -0.5")
      if ((h0 + e) <= -1) stop("For hypothesis '<', h0 + e must be > -1")
    }

  }


  # analysis prior model
  if (missing(model)) {
    stop("argument [model] for analysis prior must be one of `d_beta` (default stretched beta), `beta` (stretched beta), or `Moment` (normal-moment prior)")
  }

  # Analysis prior model validation
  if (!model %in% c("d_beta", "moment", "beta")) {
    stop("argument [model] for analysis prior must be one of `d_beta` (default stretched beta), `beta` (stretched beta), or `Moment` (normal-moment prior)")
  }

  # Model-specific checks
  if (model == "d_beta") {
    alpha=beta=scale=NaN
    # 'd_beta' requires k to be a single numeric scalar > 0
    if (!exists("k") || !is.numeric(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
      stop("For model 'd_beta', argument [k] must be a single numeric scalar > 0")
    }
  } else if (model == "beta") {
    k=scale=NaN
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (!exists("alpha") || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For model 'beta', argument [alpha] must be a single numeric scalar > 0")
    }
    if (!exists("beta") || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For model 'beta', argument [beta] must be a single numeric scalar > 0")
    }
  } else if (model == "moment") {
    k=alpha=beta=NaN
    # 'Moment' requires scale to be numeric scalar > 0
    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For model 'Moment', argument [scale] must be a numeric scalar > 0")
    }
  }


  ##  design prior

  if (!is.nan(model_d)) {

    de_an_prior <- 0

    # Validate model_d
    if (!model_d %in% c("d_beta", "moment", "beta", "point")) {
      stop("argument [model_d] for design prior must be one of `d_beta`, `beta`, `Moment`, or `Point`")
    }

    # Model-specific checks
    if (model_d == "d_beta") {
      alpha_d=beta_d=scale_d=location_d=NaN

      # 'd_beta' requires k_d to be a numeric scalar > 0
      if (!exists("k_d") || !is.numeric(k_d) || length(k_d) != 1 || !is.finite(k_d) || k_d <= 0) {
        stop("For design prior 'd_beta', argument [k_d] must be a single numeric scalar > 0")
      }
    } else if (model_d == "beta") {
      k_d=scale_d=location_d=NaN

      # 'beta' requires alpha_d and beta_d to be numeric scalars > 0
      if (!exists("alpha_d") || !is.numeric(alpha_d) || length(alpha_d) != 1 || !is.finite(alpha_d) || alpha_d <= 0) {
        stop("For design prior 'beta', argument [alpha_d] must be a single numeric scalar > 0")
      }
      if (!exists("beta_d") || !is.numeric(beta_d) || length(beta_d) != 1 || !is.finite(beta_d) || beta_d <= 0) {
        stop("For design prior 'beta', argument [beta_d] must be a single numeric scalar > 0")
      }
    } else if (model_d == "moment") {
      k_d=alpha_d=beta_d=NaN
      # 'Moment' requires scale_d numeric scalar > 0
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0) {
        stop("For design prior 'Moment', argument [scale_d] must be a numeric scalar > 0")
      }
    } else if (model_d == "point") {
      k_d=alpha_d=beta_d=scale_d=NaN

      # 'Point' requires location_d numeric scalar
      if (!is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d)) {
        stop("For design prior 'Point', argument [location_d] must be a numeric scalar")
      }
    }

  } else {
    de_an_prior <- 1
  }


  # desired power and strength of evidence
  if (mode_bf==1){
    if (!(positive %in% c("positive", "negative"))) {
      stop("argument [positive] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (positive,
                    "positive" = "h1",
                    "negative" = "h0"
    )
    if (!is.numeric(true_rate) || length(true_rate) != 1 || !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999) {
      stop("argument [true_rate] targeted true positive or negative rate must be a numeric scalar strictly greater than 0.6 and less than 0.999")
    }
    target = true_rate
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("argument [false_rate] must be a numeric scalar strictly greater than 0.001 and less than 0.1")
    }

    FP = false_rate
    if (!is.numeric(D) || length(D) != 1 || !is.finite(D) || D <= 1) {
      stop("argument [D] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=FP=0
  }


  tryCatch(
    suppressWarnings({
      if ( interval == 1) {
        results=r_table(D, target, model, k, alpha, beta, h0, location, scale, dff,
                        hypothesis, model_d, location_d, k_d, alpha_d, beta_d, scale_d,
                        dff_d, de_an_prior, N, mode_bf, FP, direct)
      } else {
        results=re_table(D, target, model, k, alpha, beta, h0, location, scale, dff,
                         hypothesis, model_d, location_d, k_d, alpha_d, beta_d, scale_d,
                         dff_d, de_an_prior, N, mode_bf, FP, e, direct)
      }
    }),
    error = function(e) {
      message("Required sample size > 5,000")
      stop(NaN)
    }
  )
  type = "correlation"
  analysis_h1 <- list(
    model = model,
    k = k,
    alpha=alpha,
    beta=beta,
    scale=scale
  )

  if (!is.nan(model_d)) {

    # Base fields always included
    design_h1 <-  list(
      model = model_d,
      location=location_d,
      k = k_d,
      alpha=alpha_d,
      beta=beta_d,
      scale=scale_d
    )


  } else {

    # model_d is NaN > fill all fields with NaN
    design_h1 <- list(
      model = NaN,
      location=NaN,
      k = NaN,
      alpha=NaN,
      beta=NaN,
      scale=NaN)

  }


  object <- list(
    type = "Correlation",
    hypothesis = hypothesis,
    h0=h0,
    e = e,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    D = D,
    mode_bf = mode_bf,
    plot_power=plot_power,
    plot_rel=plot_rel
  )
  class(object) <- "BFpower_r"
  plot(object)
  return(object)

}



#' Sample Size Determination for Bayesian F-test
#'
#' Computes required sample size or probabilities of compelling or misleading
#' evidence for a fixed sample size.
#'
#' @description
#' This function performs sample size determination (when \code{N = NaN}) or
#' calculates the probability of compelling/misleading evidence for a fixed sample
#' size.
#'
#' @param D Numeric scalar. Threshold for compelling evidence (must be > 1).
#'
#' @param true_rate Targeted true positive or true negative rate (used only when
#'   sample size determination is requested; \code{N = NaN}).
#'
#' @param false_rate Targeted false positive or false negative rate (used only when
#'   sample size determination is requested; \code{N = NaN}).
#'
#' @param p Number of predictors in the reduced model.
#'
#' @param k Number of predictors in the full model (must satisfy \code{k > p}).
#'
#' @param model Analysis prior model under the alternative hypothesis:
#'   \code{"effectsize"} or \code{"moment"}.
#'
#' @param dff Degrees of freedom for the analysis prior under the alternative
#'   hypothesis. Must be a positive scalar, and must be at least 3 if
#'   \code{model = "moment"}.
#'
#' @param rscale Scale parameter for the analysis effect-size prior (only used when
#'   \code{model = "effectsize"}).
#'
#' @param f_m Cohen's \eqn{f} effect-size parameter for the analysis prior (must be > 0).
#'
#' @param model_d Design prior model under the alternative hypothesis:
#'   \code{"effectsize"}, \code{"moment"}, or \code{"point"}.
#'
#' @param dff_d Degrees of freedom for the design prior. Must be a positive scalar,
#'   and at least 3 if \code{model_d = "moment"}.
#'
#' @param rscale_d Scale parameter for the design effect-size prior
#'   (only used when \code{model_d = "effectsize"}).
#'
#' @param f_m_d Cohen's \eqn{f} value for the design prior or the effect-size of the
#'   point design prior.
#'
#' @param N Sample size. If \code{NaN}, sample size determination is performed.
#'
#' @param positive Either `"positive"` (control true/false positive rates) or
#'   `"negative"` (control true/false negative rates).
#'
#' @param e Numeric bounds for the interval null (only used when interval
#'   Bayes factors are required).
#'
#' @param plot_power Logical. Whether to plot power curves when
#'   sample size determination is requested.
#'
#' @param plot_rel Logical. Whether to plot the relationship between the BF and data.
#'
#' @details
#'
#' \strong{1. Sample size determination mode (when \code{N = NaN}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size to achieve the desired configuration below. The user must provide:
#' \itemize{
#'   \item \code{positive} - either \code{"positive"} to control true/false positive rates, or \code{"negative"} to control true/false negative rates.
#'   \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#'   \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#'   \item \code{D} - the Bayes factor threshold for compelling evidence (must be > 1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence meets or exceeds \code{true_rate}, while the probability of misleading evidence does not exceed \code{false_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{N} is supplied):}
#'
#' If a positive numeric sample size \code{N} is provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, the arguments \code{positive}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{D} is used.
#'
#' \strong{Model specification:}
#'
#' The function requires the user to specify the full model (\code{k} predictors) and the reduced model (\code{p} predictors, \code{k > p}), and the analysis prior under the alternative hypothesis. Depending on the chosen \code{model}, different arguments are required:
#' \itemize{
#'   \item \code{model = "effectsize"}: requires \code{rscale} (scale parameter) and \code{f_m} (Cohen\'s f effect-size), and \code{dff} (degrees of freedom).
#'   \item \code{model = "moment"}: requires \code{f_m} (Cohen\'s f effect-size) and \code{dff} (degrees of freedom, must be >= 3); \code{rscale} is not used.
#' }
#' The design prior under the alternative hypothesis can optionally be specified using \code{model_d}, which can be:
#' \itemize{
#'   \item \code{"effectsize"}: requires \code{rscale_d}, \code{f_m_d}, and \code{dff_d}.
#'   \item \code{"moment"}: requires \code{f_m_d} and \code{dff_d} (>=3); \code{rscale_d} is not used.
#'   \item \code{"point"}: requires \code{f_m_d} only; \code{rscale_d} and \code{dff_d} are not used.
#' }
#'
#' \strong{Interval Null Hypothesis:}
#'
#' If \code{e} is provided, the function evaluates the Bayes factor for an interval null. Otherwise, a point-null hypothesis is assumed.
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling evidence as a function of sample size. If \code{plot_rel = TRUE}, the relationship between the Bayes factor and Cohen's \code{f} is plotted.
#'
#' @return A list of class \code{BFpower_f} containing:
#' \itemize{
#'   \item \code{type}: Test type ("Regression/ANOVA").
#'   \item \code{k}, \code{p}: Model sizes.
#'   \item \code{e}: Bounds for interval null (if used).
#'   \item \code{analysis_h1}: List describing the analysis prior.
#'   \item \code{design_h1}: List describing the design prior.
#'   \item \code{results}: Data frame of the probabilities of compelling/misleading evidence and the required or supplied sample size.
#'   \item \code{D}: Threshold of compelling evidence.
#'   \item \code{plot_power}: Logical, whether power curves are plotted.
#'   \item \code{plot_rel}: Logical, whether the relationship between the BF and data is plotted.
#' }
#'
#' If sample size determination fails, the function returns \code{NaN} and prints a message.
#'
#' @examples
#' BFpower.f.test(
#'   D = 3,
#'   true_rate = 0.8,
#'   false_rate = 0.05,
#'   p = 1,
#'   k = 2,
#'   model = "effectsize",
#'   dff = 1,
#'   rscale = 1,
#'   f_m = 0.316227766016838
#' )
#'
#' @export
BFpower.f.test <- function(D, true_rate, false_rate , p , k ,
                           model , dff , rscale , f_m ,
                           model_d = NaN, dff_d, rscale_d, f_m_d ,
                           N = NaN, positive="positive", e = NULL,plot_power=FALSE,plot_rel=FALSE) {

  ## mode
  if ( is.nan(N)) mode_bf=1 else mode_bf = 0

  ## model parameterchecks
  # Check p
  if (is.nan(p) || !is.numeric(p) || length(p) != 1 || is.na(p)) {
    stop("argument [p] (number of predictors in the reduced model) must be a positive numeric scalar")
  }

  # Check k
  if (is.nan(k) || !is.numeric(k) || length(k) != 1 || is.na(k)) {
    stop("argument [k] (number of predictors in the full model) must be a positive numeric scalar")
  }

  # Check relation
  if (k <= p) {
    stop("argument [k] (predictors in full model) must be greater than argument [p] (predictors in reduced model)")
  }

  # Equivlance test or not
  interval <- if (is.null(e)) 1 else 0

  # analysis prior model
  if (missing(model)) {
    stop("argument [model] for analysis prior should be set to either `effectsize`, or `Moment`")
  }
  if(model %in% c("effectsize","moment") == FALSE){
    stop("argument [model] for analysis prior should be set to either `effectsize`, or `Moment`")
  }

  if (model =="effectsize"){
    if (!is.numeric(rscale) || length(rscale) != 1 || !is.finite(rscale) || rscale <= 0) {
      stop("argument [rscale] scale parameter must be a positive numeric scalar")
    }
  }

  if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
    stop("argument [dff] degrees of freedom  for analysis prior must be a positive numeric scalar when model='t-distribution'")
  }

  if (!is.numeric(f_m) || length(f_m) != 1 || !is.finite(f_m) || f_m <= 0) {
    stop("argument [f_m] Cohen's f  for analysis prior must be a positive numeric scalar")
  }

  if (model == "moment"){
    rscale=NULL
    if (dff < 3) {
      stop("argument [dff] degrees of freedom for Moment analysis prior must be at least 3")
    }
  }

  # design prior

  if (!is.nan(model_d)) {

    de_an_prior <- 0

    # Validate model_d
    if (!(model_d %in% c("effectsize","moment","point"))) {
      stop("argument [model_d] for design prior must be either `effectsize`, or `Moment`")
    }



    if (model_d =="effectsize"){
      if (!is.numeric(rscale_d) || length(rscale_d) != 1 || !is.finite(rscale_d) || rscale_d <= 0) {
        stop("argument [rscale] scale parameter must be a positive numeric scalar")
      }

      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0) {
        stop("argument [dff] degrees of freedom  for design prior must be a positive numeric scalar when model='t-distribution'")
      }
      if (!is.numeric(f_m_d) || length(f_m_d) != 1 || !is.finite(f_m_d) || f_m_d <= 0) {
        stop("argument [f_m] Cohen's f  for design prior must be a positive numeric scalar")
      }


      }



    if (model_d == "moment"){
      rscale_d=NaN
      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0) {
        stop("argument [dff] degrees of freedom  for design prior must be a positive numeric scalar when model='t-distribution'")
      }
      if (!is.numeric(f_m_d) || length(f_m_d) != 1 || !is.finite(f_m_d) || f_m_d <= 0) {
        stop("argument [f_m] Cohen's f  for design prior must be a positive numeric scalar")
      }

      if (dff_d < 3) {
        stop("argument [dff_d] degrees of freedom for moment design prior must be at least 3")
      }
    }

    if (model_d == "point"){
      rscale_d=dff_d=NaN

      if (!is.numeric(f_m_d) || length(f_m_d) != 1 || !is.finite(f_m_d) || f_m_d <= 0) {
        stop("argument [f_m] Cohen's f  for design prior must be a positive numeric scalar")
      }

    }

  } else {
    de_an_prior <- 1
  }


  # desired power and strength of evidence
  if (mode_bf==1){
    if (!(positive %in% c("positive", "negative"))) {
      stop("argument [positive] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (positive,
                    "positive" = "h1",
                    "negative" = "h0"
    )
    if (!is.numeric(true_rate) || length(true_rate) != 1 ||
        !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999){
      stop("argument [true_rate] targeted true positive or negative rate must be a numeric scalar strictly greater than 0.6 and less than 0.999")
    }
    target = true_rate
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("argument [false_rate] must be a numeric scalar strictly greater than 0.001 and less than 0.1")
    }

    FP = false_rate
    if (!is.numeric(D) || length(D) != 1 || !is.finite(D) || D <= 1) {
      stop("argument [D] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=FP=0
  }









  results = tryCatch({
    suppressWarnings({
      if (!is.nan(interval) && interval == 1) {
        f_table(D, target, p, k, dff, rscale, f_m, model,
                dff_d, rscale_d, f_m_d, model_d, de_an_prior, N,
                mode_bf, FP, direct)
      } else {
        fe_table(D, target, p, k, dff, rscale, f_m, model,
                 dff_d, rscale_d, f_m_d, model_d, de_an_prior, N,
                 mode_bf, e, FP, direct)
      }
    })
  }, error = function(e) {

    if(dff<3|dff_d<3){
      stop(" Degrees of freedom[dff] for analysis prior or [dff_d] for design prior should be at least 3")

      }

    message("Required sample size > 10,000")

    return(invisible(NaN))
  })


  type = "Regression/ANOVA"
  analysis_h1 <- list(
    model = model,
    rscale = rscale,
    f_m = f_m,
    dff=dff
  )

  if (!is.nan(model_d)) {

    # Base fields always included
    design_h1 <- list(
      model = model_d,
      rscale = rscale_d,
      f_m = f_m_d,
      dff=dff_d)

  } else {

    # model_d is NaN > fill all fields with NaN
    design_h1 <- list(
      model  = NaN,
      rscale = NULL,
      f_m    = NaN,
      dff    = NaN)
  }


  object <- list(
    type = type,
    k=k,
    p=p,
    e = e,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    D = D,
    mode_bf = mode_bf,
    plot_power=plot_power,
    plot_rel=plot_rel
  )
  class(object) <- "BFpower_f"
  plot(object)
  return(object)

}


#' Sample size determination for Bayesian one-proportion test
#'
#' Perform sample size determination or the calculation of compelling and misleading evidence
#' for a Bayesian test of a single proportion.
#'
#' @param hypothesis The hypothesis being tested: two-sided (\code{"!="}), right-sided (\code{">"}), or left-sided (\code{"<"}).
#' @param D Numeric scalar. Threshold for compelling evidence (must be > 1).
#' @param h0 Null proportion value for the test (numeric scalar between 0.1 and 0.9).
#' @param true_rate Targeted true positive rate (if \code{direct = "h1"}) or true negative rate (if \code{direct = "h0"}). Used only when sample size determination is requested.
#' @param false_rate Targeted false positive rate (if \code{direct = "h1"}) or false negative rate (if \code{direct = "h0"}). Used only when sample size determination is requested.
#' @param model Analysis prior model under the alternative hypothesis: \code{"beta"} or \code{"moment"}.
#' @param alpha Parameter for the analysis beta prior (used when \code{model = "beta"}).
#' @param beta Parameter for the analysis beta prior (used when \code{model = "beta"}).
#' @param scale Scale parameter for the analysis moment prior (used when \code{model = "moment"}).
#' @param model_d Design prior model under the alternative hypothesis: \code{"beta"}, \code{"moment"}, or \code{"point"}.
#' @param alpha_d Parameter for the design beta prior (used when \code{model_d = "beta"}).
#' @param beta_d Parameter for the design beta prior (used when \code{model_d = "beta"}).
#' @param location_d Proportion value for the design point prior (\code{model_d = "point"}). Represents the true proportion under the alternative hypothesis.
#' @param scale_d Scale parameter for the design moment prior (used when \code{model_d = "moment"}).
#' @param N Sample size. If \code{NaN}, sample size determination is performed.
#' @param e Numeric bounds for the interval null (used when computing interval Bayes factors).
#' @param positive Either `"positive"` (controls true/false positive rates) or `"negative"` (controls true/false negative rates).
#' @param plot_power Logical. Whether to plot power curves when sample size determination is requested.
#' @param plot_rel Logical. Whether to plot probability of misleading evidence.
#'
#' @details
#'
#' \strong{1. Sample size determination mode (when \code{N = NaN}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size needed to achieve the desired configuration below. The user must provide:
#' \itemize{
#' \item \code{positive} - either \code{"positive"} to control true/false positive rates or \code{"negative"} to control true/false negative rates.
#' \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#' \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#' \item \code{D} - the Bayes factor threshold for compelling evidence (must be > 1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence meets or exceeds \code{true_rate}, while the probability of misleading evidence does not exceed \code{false_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{N} is supplied):}
#'
#' If a positive numeric sample size \code{N} is provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, \code{positive}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{D} is used.
#'
#' \strong{Model specification:}
#'
#' The user must specify the analysis prior under the alternative hypothesis using \code{model}:
#' \itemize{
#' \item \code{model = "beta"}: requires \code{alpha} and \code{beta} parameters (shape parameters of the beta distribution).
#' \item \code{model = "moment"}: requires \code{scale} parameter (scale of the moment prior).
#' }
#' The design prior under the alternative hypothesis can optionally be specified using \code{model_d}:
#' \itemize{
#' \item \code{"beta"}: requires \code{alpha_d} and \code{beta_d}.
#' \item \code{"moment"}: requires \code{scale_d}.
#' \item \code{"point"}: requires \code{location_d}, representing the true proportion under the alternative hypothesis.
#' }
#' If \code{model_d} is \code{NaN}, no design prior is used.
#'
#' \strong{interval null Hypothesis:}
#'
#' If \code{e} is provided, the function evaluates the Bayes factor for an interval null. Otherwise, a point-null hypothesis is assumed.
#'
#' \strong{Hypothesis:}
#'
#' The function supports one-sided (\code{">"} or \code{"<"}) and two-sided (\code{"!="}) tests. Design prior and interval null bounds must be consistent with the directionality of the hypothesis.
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling evidence as a function of sample size. If \code{plot_rel = TRUE}, the relationship between the Bayes factor and the number of successes (proportion) is plotted.
#'
#'
#' @return A list of class \code{"BFpower_bin"} containing:
#' \itemize{
#'   \item \code{type}: Test type ("One proportion").
#'   \item \code{hypothesis}: Tested hypothesis.
#'   \item \code{h0}: Null proportion.
#'   \item \code{analysis_h1}: List describing the analysis prior.
#'   \item \code{design_h1}: List describing the design prior.
#'   \item \code{results}: Data frame of probabilities of compelling/misleading evidence and the required or supplied sample size.
#'   \item \code{D}: Compelling-evidence threshold.
#'   \item \code{plot_power}: Logical, whether power curves are plotted.
#'   \item \code{plot_rel}: Logical, whether the relationship between the BF and data is plotted.
#' }
#'
#' If sample size determination fails, the function returns \code{NaN} and prints a message.
#'
#' @examples
#' BFpower.bin(
#'   hypothesis = "!=",
#'   D = 3,
#'   true_rate = 0.8,
#'   false_rate = 0.05,
#'   h0=.5,
#'   model = "beta",
#'   alpha = 1,
#'   beta = 1
#' )
#'
#'
#' @export
BFpower.bin <- function(hypothesis ,D , h0 ,
                        true_rate , false_rate ,
                        model , alpha , beta , scale ,
                        model_d = NaN, alpha_d , beta_d , location_d , scale_d ,
                        N = NaN, e = NULL, positive="positive",plot_power=FALSE,plot_rel=FALSE) {

  # mode
  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < .1 || h0 > 0.9) {
    stop("argument [h0] NaN value of proportion must be a single numeric scalar between .1 and 0.9")
  }

  location <- h0
  if ( is.nan(N)) mode_bf=1 else mode_bf = 0


  # sample size
  if (mode_bf == 0) {
    # Check that N is a positive numeric scalar
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 0) {
      stop("argument [N] sample size must be a positive numeric scalar when mode_bf = 0")
    }
  }else {N=3}


  # hypothesis
  if(hypothesis %in% c("<", "!=", ">") == FALSE){
    stop("argument [hypothesis] should be set to either `<`  (left-sided test),  `!=` (two-sided test) or `>` (right-sided test)")
  }

  # Equivlance test or not
  interval <- if (is.null(e)) 1 else 0

  if (!is.null(e)) {

    if (hypothesis == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(e) || length(e) != 2 || any(!is.finite(e)) || e[1] == e[2]) {
        stop("For hypothesis '!=', argument [e] must be a numeric vector of length 2 with two distinct finite values")
      }
      # Additional bounds checks
      if (min(e) < -0.5 || max(e) > 0.5) {
        stop("For hypothesis '!=', e must satisfy min(e) >= -0.5 and max(e) <= 0.5")
      }
      if ((h0 + min(e)) <= 0 || (h0 + min(e)) >= 1) {
        stop("For hypothesis '!=', h0 + min(e) must be between 0 and 1")
      }

    } else if (hypothesis == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e <= 0) {
        stop("For hypothesis '>', argument [e] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (e > 0.5) stop("For hypothesis '>', e must be <= 0.5")
      if ((h0 + e) >= 1) stop("For hypothesis '>', h0 + e must be < 1")

    } else if (hypothesis == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e >= 0) {
        stop("For hypothesis '<', argument [e] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (e < -0.5) stop("For hypothesis '<', e must be >= -0.5")
      if ((h0 + e) <= -1) stop("For hypothesis '<', h0 + e must be > 0")
    }

  }


  # analysis prior model
  if (missing(model)) {
    stop("argument [model] for analysis prior must be one of `beta`, or `Moment` (normal-moment prior)")
  }

  # Analysis prior model validation
  if (!model %in% c("moment", "beta")) {
    stop("argument [model] for analysis prior must be one of `beta` , or `Moment` (normal-moment prior)")
  }

  # Model-specific checks
  if (model == "beta") {
    scale=NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (!exists("alpha") || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For model 'beta', argument [alpha] must be a single numeric scalar > 0")
    }
    if (!exists("beta") || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For model 'beta', argument [beta] must be a single numeric scalar > 0")
    }
  } else if (model == "moment") {
    alpha=beta=NaN
    # 'Moment' requires scale to be numeric scalar > 0
    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For model 'Moment', argument [scale] must be a numeric scalar > 0")
    }
  }


  ##  design prior

  if (!is.nan(model_d)) {

    de_an_prior <- 0

    # Validate model_d
    if (!model_d %in% c("moment", "beta", "point")) {
      stop("argument [model_d] for design prior must be one of `beta`, `Moment`, or `Point`")
    }

    # Model-specific checks
    if (model_d == "beta") {
      scale_d=location_d=NaN

      # 'beta' requires alpha_d and beta_d to be numeric scalars > 0
      if (!exists("alpha_d") || !is.numeric(alpha_d) || length(alpha_d) != 1 || !is.finite(alpha_d) || alpha_d <= 0) {
        stop("For design prior 'beta', argument [alpha_d] must be a single numeric scalar > 0")
      }
      if (!exists("beta_d") || !is.numeric(beta_d) || length(beta_d) != 1 || !is.finite(beta_d) || beta_d <= 0) {
        stop("For design prior 'beta', argument [beta_d] must be a single numeric scalar > 0")
      }
    } else if (model_d == "moment") {
      alpha_d=beta_d=NaN
      # 'Moment' requires scale_d numeric scalar > 0
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0) {
        stop("For design prior 'Moment', argument [scale_d] must be a numeric scalar > 0")
      }
    } else if (model_d == "point") { alpha_d <- beta_d <- scale_d <- NaN  # Not needed for 'Point' prior

    # 'Point' prior requires location_d, which represents the true proportion
    # under the alternative hypothesis. It must be a numeric scalar.
    if (!is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d)) {
      stop("For design prior 'Point', argument [location_d] true proportion must be a numeric scalar")
    }

    # Validate location_d against hypothesis and h0
    # - "!=" : location_d must not equal h0
    # - ">"  : location_d must be greater than h0
    # - "<"  : location_d must be less than h0
    if (hypothesis == "!=" && location_d == h0) {
      stop("For hypothesis '!=', argument [location_d] true proportion must not equal h0")
    } else if (hypothesis == ">" && location_d <= h0) {
      stop("For hypothesis '>', argument [location_d] true proportion must be greater than h0")
    } else if (hypothesis == "<" && location_d >= h0) {
      stop("For hypothesis '<', argument [location_d] true proportion must be less than h0")
    }
    }

  } else {
    de_an_prior <- 1
  }


  # desired power and strength of evidence
  if (mode_bf==1){
    if (!(positive %in% c("positive", "negative"))) {
      stop("argument [positive] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (positive,
                    "positive" = "h1",
                    "negative" = "h0"
    )
    if (!is.numeric(true_rate) || length(true_rate) != 1 || !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999) {
      stop("argument [true_rate] targeted true positive or negative rate must be a numeric scalar strictly greater than 0.6 and less than 0.999")
    }
    target = true_rate
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("argument [false_rate] must be a numeric scalar strictly greater than 0.001 and less than 0.1")
    }

    FP = false_rate
    if (!is.numeric(D) || length(D) != 1 || !is.finite(D) || D <= 1) {
      stop("argument [D] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=FP=0
  }



  results <-tryCatch({
    suppressWarnings({
      if (is.null(e)) {
        bin_table(D, target, h0, alpha, beta, location, scale, model, hypothesis,
                  alpha_d, beta_d, location_d, scale_d, model_d, de_an_prior, N,
                  mode_bf, FP, direct)
      } else {
        bin_e_table(D, target, h0, alpha, beta, location, scale, model, hypothesis,
                    alpha_d, beta_d, location_d, scale_d, model_d, de_an_prior, N,
                    mode_bf, FP, e, direct)
      }
    })
  }, error = function(e) {
    message("Sample size cannot be determined")
    return(invisible(NaN))
  })

  type = "One proportion"
  analysis_h1 <- list(
    model = model,
    alpha=alpha,
    beta=beta,
    scale=scale
  )

  if (!is.nan(model_d)) {

    # Base fields always included
    design_h1 <-  list(
      model = model_d,
      alpha=alpha_d,
      beta=beta_d,
      scale=scale_d
    )


  } else {

    # model_d is NaN > fill all fields with NaN
    design_h1 <- list(
      model = NaN,
      alpha=NaN,
      beta=NaN,
      scale=NULL)

  }


  object <- list(
    type = type,
    hypothesis = hypothesis,
    h0=h0,
    e = e,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    D = D,
    mode_bf = mode_bf,
    plot_power=plot_power,
    plot_rel=plot_rel
  )
  class(object) <- "BFpower_bin"
  plot(object)
  return(object)

}

#' Sample size determination for Bayesian test of two proportions
#'
#' Perform sample size determination or calculate probabilities of compelling and misleading evidence
#' for a Bayesian comparison of two proportions.
#'
#' @param D Threshold of compelling evidence.
#' @param true_rate Targeted true positive rate (if \code{positive = "positive"}) or true negative rate (if \code{positive = "negative"}).
#' @param a0 Alpha parameter of the Beta prior under the null hypothesis.
#' @param b0 Beta parameter of the Beta prior under the null hypothesis.
#' @param a1 Alpha parameter of the Beta analysis prior for group 1 under the alternative hypothesis.
#' @param b1 Beta parameter of the Beta analysis prior for group 1 under the alternative hypothesis.
#' @param a2 Alpha parameter of the Beta analysis prior for group 2 under the alternative hypothesis.
#' @param b2 Beta parameter of the Beta analysis prior for group 2 under the alternative hypothesis.
#' @param model1 Model for the design prior of group 1: \code{"beta"}, \code{"point"}, or \code{"same"} (if \code{"same"}, the design prior is identical to the analysis prior).
#' @param a1d Alpha parameter of the design prior for group 1 (used if \code{model1 = "beta"}).
#' @param b1d Beta parameter of the design prior for group 1 (used if \code{model1 = "beta"}).
#' @param dp1 True proportion for group 1 in the design prior (used if \code{model1 = "point"}).
#' @param model2 Model for the design prior of group 2: \code{"beta"}, \code{"point"}, or \code{"same"} (if \code{"same"}, the design prior is identical to the analysis prior).
#' @param a2d Alpha parameter of the design prior for group 2 (used if \code{model2 = "beta"}).
#' @param b2d Beta parameter of the design prior for group 2 (used if \code{model2 = "beta"}).
#' @param dp2 True proportion for group 2 in the design prior (used if \code{model2 = "point"}).
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param positive Choose \code{"positive"} to control true/false positive rates or \code{"negative"} to control true/false negative rates.
#' @param plot_power Logical; if TRUE, plot the power curve.
#' @param plot_rel Logical; if TRUE, plot the grid for the values of BF across all possible combination of x1 and x2.
#'
#' @details
#'
#' \strong{1. Sample size determination mode (when \code{n1 = NaN} and \code{n2 = NaN}):}
#'
#' If no sample sizes are provided for the two groups, the function calculates the minimum sample sizes needed to achieve the desired configuration. The user must provide:
#' \itemize{
#' \item \code{positive} - either \code{"positive"} to control true/false positive rates or \code{"negative"} to control true/false negative rates.
#' \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#' \item \code{D} - the Bayes factor threshold for compelling evidence (must be > 1).
#' }
#'
#' The function iteratively finds the smallest sample sizes for which the probability of obtaining compelling evidence meets or exceeds \code{true_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{n1} and \code{n2} are supplied):}
#'
#' If positive numeric sample sizes \code{n1} and \code{n2} are provided, the function computes the probabilities of obtaining compelling or misleading evidence for these fixed sample sizes. In this mode, \code{positive} and \code{true_rate} are ignored; only the Bayes factor threshold \code{D} is used.
#'
#' \strong{Model specification:}
#'
#' The user must specify the analysis priors under the null and alternative hypotheses using Beta parameters:
#' \itemize{
#' \item \code{a0}, \code{b0} - Beta parameters for the null hypothesis prior.
#' \item \code{a1}, \code{b1} - Beta parameters for the analysis prior of group 1 under the alternative hypothesis.
#' \item \code{a2}, \code{b2} - Beta parameters for the analysis prior of group 2 under the alternative hypothesis.
#' }
#'
#' Design priors for the alternative hypothesis can optionally be specified:
#' \itemize{
#' \item \code{model1}, \code{a1d}, \code{b1d}, \code{dp1} - design prior for group 1 (\code{"same"} uses the analysis prior, \code{"beta"} requires Beta parameters, \code{"point"} uses a fixed proportion).
#' \item \code{model2}, \code{a2d}, \code{b2d}, \code{dp2} - design prior for group 2.
#' }
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, a power curve is plotted showing the probability of compelling evidence as a function of sample sizes. If \code{plot_rel = TRUE}, a grid of Bayes factors across possible outcomes is plotted.
#'
#' @return An object of class \code{BFpower_2p} containing:
#'   \itemize{
#'     \item Analysis priors under the null and alternative.
#'     \item Design priors for both groups.
#'     \item Computed probabilities of compelling and misleading evidence.
#'     \item Sample size results and other inputs.
#'   }
#'
#' @examples
#' BFpower.props(
#'   D = 3,
#'   true_rate = 0.8,
#'   a0 = 1,
#'   b0 = 1,
#'   a1 = 1,
#'   b1 = 1,
#'   a2 = 1,
#'   b2 = 1
#' )
#'
#' @export

BFpower.props <- function(D , true_rate , a0 , b0 , a1 , b1 ,
                          a2 , b2 , model1 = "same",
                          a1d , b1d , dp1 , model2 = "same",
                          a2d, b2d , dp2 ,
                          n1 = NaN, n2 = NaN,positive="positive",plot_power=FALSE,plot_rel=FALSE) {

  # Check NaN
  if (is.nan(n1) && is.nan(n2)) {
    mode_bf <- 1
  } else {
    mode_bf <- 0
  }

  # If both n1 and n2 are NaN > mode_bf = 1
  # If mode_bf = 1, check target range
  if (mode_bf==1){
    if (!(positive %in% c("positive", "negative"))) {
      stop("argument [positive] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (positive,
                    "positive" = "h1",
                    "negative" = "h0"
    )
    if (!is.numeric(true_rate) || length(true_rate) != 1 || !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999) {
      stop("argument [true_rate] targeted true positive or negative rate must be a numeric scalar strictly greater than 0.6 and less than 0.999")
    }
    target = true_rate

    if (!is.numeric(D) || length(D) != 1 || !is.finite(D) || D <= 1) {
      stop("argument [D] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=0
  }



  # If not NaN, check numeric, scalar, integer
  if (mode_bf == 0) {

    # Check n1
    if (!is.numeric(n1) || length(n1) != 1 || n1 %% 1 != 0 || n1 <= 0) {
      stop("arg [n1] sample size for group 1 must be a positive numeric scalar integer (> 0).")
    }

    # Check n2
    if (!is.numeric(n2) || length(n2) != 1 || n2 %% 1 != 0 || n2 <= 0) {
      stop("arg [n2] sample size for group 2 must be a positive numeric scalar integer (> 0).")
    }
  }

  # NaN hypothesis
  # Check a0 (alpha)
  if (!is.numeric(a0) || length(a0) != 1 || a0 <= 0) {
    stop("arg [a0] alpha for the Beta analysis prior under the null (\u03b8\u2080) must be a positive numeric scalar (> 0).")
  }

  # Check b0 (beta)
  if (!is.numeric(b0) || length(b0) != 1 || b0 <= 0) {
    stop("arg [b0] beta for the Beta analysis prior under the null (\u03b8\u2080) must be a positive numeric scalar (> 0).")
  }

  # alternative hypothesis theta1
  # Check a1 (alpha under the alternative)
  if (!is.numeric(a1) || length(a1) != 1 || a1 <= 0) {
    stop("arg [a1] alpha for the Beta analysis prior under the alternative (\u03b8\u2081) must be a positive numeric scalar (> 0).")
  }

  # Check b1 (beta under the alternative)
  if (!is.numeric(b1) || length(b1) != 1 || b1 <= 0) {
    stop("arg [b1] beta for the Beta analysis prior under the alternative (\u03b8\u2081) must be a positive numeric scalar (> 0).")
  }

  # alternative hypothesis theta2
  # Check a2 (alpha under the alternative)
  if (!is.numeric(a2) || length(a2) != 1 || a2 <= 0) {
    stop("arg [a2] alpha for the Beta analysis prior under the alternative (\u03b8\u2082) must be a positive numeric scalar (> 0).")
  }

  # Check b2 (beta under the alternative)
  if (!is.numeric(b2) || length(b2) != 1 || b2 <= 0) {
    stop("arg [b2] beta for the Beta analysis prior under the alternative (\u03b8\u2082) must be a positive numeric scalar (> 0).")
  }


  # --- Check model1 assumptions for design prior on theta1 ---

  if (model1 == "same") {

    # Automatically set all to NaN
    a1d <- 1
    b1d <- 1
    dp1 <- 0.5

  } else if (model1 == "beta") {

    # a1d and b1d must be valid Beta parameters
    if (!is.numeric(a1d) || length(a1d) != 1 || a1d <= 0) {
      stop("arg [a1d] alpha for the Beta design prior on \u03b8\u2081 must be a positive numeric scalar (> 0).")
    }
    if (!is.numeric(b1d) || length(b1d) != 1 || b1d <= 0) {
      stop("arg [b1d] beta for the Beta design prior on \u03b8\u2081 must be a positive numeric scalar (> 0).")
    }

    # dp1 irrelevant for beta prior > set to NaN automatically
    dp1 <- 0.5

  } else if (model1 == "point") {

    # Automatically set Beta parameters to NaN
    a1d <- 1
    b1d <- 1

    # dp1 must be numeric between 0 and 1
    if (!is.numeric(dp1) || length(dp1) != 1) {
      stop("arg [dp1] true \u03b8\u2081 must be a numeric scalar for model1 = 'Point'.")
    }
    if (dp1 <= 0 || dp1 >= 1) {
      stop("arg [dp1] must be > 0 and < 1 for model1 = 'Point'.")
    }

  } else {
    stop("arg [model1] must be one of: 'same', 'beta', 'Point'.")
  }

  # --- Check model2 assumptions for design prior on theta2 ---

  if (model2 == "same") {

    # Automatically set all to NaN
    a2d <- 1
    b2d <- 1
    dp2 <- .5

  } else if (model2 == "beta") {

    # a2d and b2d must be valid Beta parameters
    if (!is.numeric(a2d) || length(a2d) != 1 || a2d <= 0) {
      stop("arg [a2d] alpha for the Beta design prior on theta2 must be a positive numeric scalar (> 0).")
    }
    if (!is.numeric(b2d) || length(b2d) != 1 || b2d <= 0) {
      stop("arg [b2d] beta for the Beta design prior on theta2 must be a positive numeric scalar (> 0).")
    }

    # dp2 irrelevant for beta prior > set to NaN automatically
    dp2 <- .5

  } else if (model2 == "point") {

    # Automatically set Beta parameters to NaN
    a2d <- 1
    b2d <- 1

    # dp2 must be numeric between 0 and 1
    if (!is.numeric(dp2) || length(dp2) != 1) {
      stop("arg [dp2] must be a numeric scalar for model2 = 'Point'.")
    }
    if (dp2 <= 0 || dp2 >= 1) {
      stop("arg [dp2] must be > 0 and < 1 for model2 = 'Point'.")
    }

  } else {
    stop("arg [model2] must be one of: 'same', 'beta', 'Point'.")
  }



  r <- 1
  results=tryCatch({
    suppressWarnings({
      pro_table_p2(D, target, a0, b0, a1, b1, a2, b2, r, model1,
                   a1d, b1d, dp1, model2, a2d, b2d, dp2, mode_bf, n1, n2, direct)
    })
  }, error = function(e) {
    message("Required Sample size > 5000 per group")
    return(NaN)
  })



  type = "Two-proportions"
  analysis_h0 <- list(
    a = a0,
    b = b0
  )
  analysis_h1_theta_1 <- list(
    a = a1,
    b = b1
  )
  analysis_h1_theta_2 <- list(
    a = a2,
    b = b2
  )

  design_h1_theta_1 <- list(
    model=model1,
    a = a1d,
    b = b1d,
    p = dp1
  )
  design_h1_theta_2 <- list(
    model=model2,
    a = a2d,
    b = b2d,
    p = dp2
  )

  object <- list(
    type = type,
    analysis_h0=analysis_h0,
    analysis_h1_theta_1= analysis_h1_theta_1,
    analysis_h1_theta_2=analysis_h1_theta_2,
    design_h1_theta_1=design_h1_theta_1,
    design_h1_theta_2=design_h1_theta_2,
    results = results[[1]],
    grid=results[[2]],
    D = D,
    mode_bf = mode_bf,
    plot_power=plot_power,
    plot_rel=plot_rel
  )
  class(object) <- "BFpower_2p"
  plot(object)
  return(object)

}


#' Bayes Factor for a One-Sample Bayesian t-Test
#'
#' Computes the Bayes factor (BF10) for a one-sample t-test, comparing an observed t-value
#' against either a point NaN hypothesis or an interval null hypothesis.
#'
#' @param tval Numeric scalar. Observed t-value from the one-sample t-test.
#' @param df Numeric scalar. Degrees of freedom of the t-test (must be >= 1).
#' @param model Character string. Statistical model for the analysis prior under the alternative hypothesis.
#'   Choices are:
#'   \describe{
#'     \item{"normal"}{Normal distribution.}
#'     \item{"moment"}{Normal moment prior.}
#'     \item{"t"}{Scaled t-distribution.}
#'   }
#' @param location Numeric scalar. Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Numeric scalar. Scale parameter for the analysis prior under the alternative hypothesis (must be > 0).
#' @param dff Numeric scalar. Degrees of freedom for the t-distribution prior (only required if \code{model = "t"}; must be > 0). Ignored otherwise.
#' @param hypothesis Character string. Hypothesis being tested. One of:
#'   \describe{
#'     \item{"!="}{Two-sided (difference from 0).}
#'     \item{">"}{Right-sided (greater than 0).}
#'     \item{"<"}{Left-sided (less than 0).}
#'   }
#' @param e Optional numeric vector. Specifies bounds for an interval null hypothesis. For:
#'   \describe{
#'     \item{Two-sided (\code{"!="})}{Must be a numeric vector of length 2 with two distinct finite values.}
#'     \item{Right-sided (\code{">"})}{Must be a numeric scalar > 0.}
#'     \item{Left-sided (\code{"<"})}{Must be a numeric scalar < 0.}
#'   }
#'
#' @return A list of class \code{BFvalue_t} containing:
#'   \describe{
#'     \item{type}{Character, indicating "One-sample t-test".}
#'     \item{bf10}{Numeric, the Bayes factor (BF10).}
#'     \item{tval}{Observed t-value.}
#'     \item{df}{Degrees of freedom.}
#'     \item{analysis_h1}{List with the analysis prior parameters under H1: \code{model}, \code{location}, \code{scale}, and optionally \code{dff}.}
#'     \item{hypothesis}{Character, the tested hypothesis.}
#'     \item{e}{Optional numeric vector of interval null bounds.}
#'   }
#'
#' @examples
#' BF10.t.test.one.sample(
#'   tval = 2.31,
#'   df = 29,
#'   model = "t",
#'   location = 0,
#'   scale = 0.707,
#'   dff = 1,
#'   hypothesis = "!="
#' )
#'
#'
#' @export
BF10.t.test.one.sample <- function(tval, df, model, location, scale, dff, hypothesis, e = NULL) {

  ## -----------------------------
  ## Input validation
  ## -----------------------------

  # tval must be a numeric scalar
  if (!is.numeric(tval) || length(tval) != 1 || !is.finite(tval)) {
    stop("argument [tval] observed t-value must be a numeric scalar")
  }

  # df must be numeric >= 1
  if (!is.numeric(df) || length(df) != 1 || !is.finite(df) || df < 1) {
    stop("argument [df] degree of freedom must be a numeric scalar >= 1")
  }

  # Check e if provided
  if (!is.null(e)) {
    if (hypothesis == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(e) || length(e) != 2 || any(!is.finite(e)) || e[1] == e[2]) {
        stop("For hypothesis '!=', argument [e] must be a numeric vector of length 2 with two distinct finite values")
      }
    }
    if (hypothesis == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e <= 0) {
        stop("For hypothesis '>', argument [e] must be a numeric scalar > 0")
      }
    }
    if (hypothesis == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e >= 0) {
        stop("For hypothesis '<', argument [e] must be a numeric scalar < 0")
      }
    }
  }

  # Validate analysis prior
  if (missing(model)) {
    stop("argument [model] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }

  if (!model %in% c("normal","moment","t")) {
    stop("argument [model] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("argument [location] for analysis prior must be a numeric scalar")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
    stop("argument [scale] for analysis prior must be a positive numeric scalar (i.e., scale > 0)")
  }
  if (model == "t") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("argument [dff] degrees of freedom for analysis prior must be a positive numeric scalar when model='t-distribution'")
    }
  } else {
    dff <- 0
  }
  ## -----------------------------
  ## Call appropriate function
  ## -----------------------------
  suppressWarnings(
    if (is.null(e)) {
      bf10=t1_BF10(tval, df, model, location, scale, dff, hypothesis)
    } else {
      bf10=t1e_BF10(tval, df, model, location, scale, dff, hypothesis, e)
    }
  )

  type = "One-sample t-test"
  analysis_h1 <- list(
    model = model,
    location = location,
    scale = scale
  )
  if (model == "t") {
    analysis_h1$dff <- dff
  }
  object=list(type=type,bf10=bf10,tval=tval,df=df,analysis_h1=analysis_h1,hypothesis=hypothesis,e=e)

  class(object) <- "BFvalue_t"

  return(object)
}

#' Bayes Factor for Two-Sample Bayesian t-Test
#'
#' Compute the Bayes factor (BF10) for a two-sample independent-samples t-test
#' under a specified prior model. Supports both point-NaN and interval-NaN hypotheses.
#'
#' @param tval Numeric scalar. Observed t-value from the two-sample t-test.
#' @param N1 Numeric scalar. Sample size of group 1 (must be > 2, will be rounded to nearest integer).
#' @param N2 Numeric scalar. Sample size of group 2 (must be > 2, will be rounded to nearest integer).
#' @param model Character. Analysis prior under the alternative hypothesis:
#'   \code{"normal"}, \code{"moment"} (normal-moment prior), or \code{"t"}.
#' @param location Numeric scalar. Location parameter of the analysis prior.
#' @param scale Numeric scalar > 0. Scale parameter of the analysis prior.
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior (required if model = \code{"t"}; ignored otherwise).
#' @param hypothesis Character. Hypothesis to test: two-sided (\code{"!="}), right-sided (\code{">"}), or left-sided (\code{"<"}).
#' @param e Optional numeric. Bounds for an interval null:
#'   - For \code{hypothesis = "!="}, must be a numeric vector of length 2 with distinct finite values.
#'   - For \code{">"}, must be a single numeric scalar > 0.
#'   - For \code{"<"}, must be a single numeric scalar < 0.
#'
#' @return A list of class \code{BFvalue_t} containing:
#' \describe{
#'   \item{type}{Character string describing the test type.}
#'   \item{bf10}{Computed Bayes factor BF10.}
#'   \item{tval}{Observed t-value.}
#'   \item{df}{Degrees of freedom (currently NA / not computed).}
#'   \item{analysis_h1}{List of prior model parameters used for the alternative hypothesis.}
#'   \item{hypothesis}{Hypothesis tested (\code{"!="}, \code{">"}, or \code{"<"}).}
#'   \item{e}{Interval bounds used, if any.}
#'   \item{N1}{Sample size of group 1 (rounded integer).}
#'   \item{N2}{Sample size of group 2 (rounded integer).}
#' }
#'
#' @examples
#' BF10.t.test.two.sample(
#'   tval = 2.1,
#'   N1 = 30,
#'   N2 = 30,
#'   model = "t",
#'   location = 0,
#'   scale = 0.707,
#'   dff = 1,
#'   hypothesis = "!="
#' )
#'
#' # Using an interval null for a right-sided hypothesis
#' BF10.t.test.two.sample(
#'   tval = 1.8,
#'   N1 = 25,
#'   N2 = 25,
#'   model = "normal",
#'   location = 0,
#'   scale = 1,
#'   dff = 0,
#'   hypothesis = ">",
#'   e = 0.2
#' )
#'
#' @export
BF10.t.test.two.sample <- function(tval, N1, N2, model, location, scale, dff, hypothesis, e = NULL) {

  ## -----------------------------
  ## Input validation
  ## -----------------------------

  # tval must be a numeric scalar
  if (!is.numeric(tval) || length(tval) != 1 || !is.finite(tval)) {
    stop("argument [tval] observed t-value must be a numeric scalar")
  }

  # Example for N1
  if (!is.numeric(N1) || length(N1) != 1 || !is.finite(N1)) {
    stop("argument [N1] must be a numeric scalar")
  }
  if (N1 <= 2) {
    stop("argument [N1] must be greater than 2")
  }
  # Round to nearest integer
  N1 <- round(N1)

  # Similarly for N2
  if (!is.numeric(N2) || length(N2) != 1 || !is.finite(N2)) {
    stop("argument [N2] must be a numeric scalar")
  }
  if (N2 <= 2) {
    stop("argument [N2] must be greater than 2")
  }
  N2 <- round(N2)


  # Check e if provided
  if (!is.null(e)) {
    if (hypothesis == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(e) || length(e) != 2 || any(!is.finite(e)) || e[1] == e[2]) {
        stop("For hypothesis '!=', argument [e] must be a numeric vector of length 2 with two distinct finite values")
      }
    }
    if (hypothesis == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e <= 0) {
        stop("For hypothesis '>', argument [e] must be a numeric scalar > 0")
      }
    }
    if (hypothesis == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e >= 0) {
        stop("For hypothesis '<', argument [e] must be a numeric scalar < 0")
      }
    }
  }

  # Validate analysis prior
  if (missing(model)) {
    stop("argument [model] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }

  if (!model %in% c("normal","moment","t")) {
    stop("argument [model] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("argument [location] for analysis prior must be a numeric scalar")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
    stop("argument [scale] for analysis prior must be a positive numeric scalar (i.e., scale > 0)")
  }
  if (model == "t") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("argument [dff] degrees of freedom for analysis prior must be a positive numeric scalar when model='t-distribution'")
    }
  } else {
    dff <- 0
  }

  n1 <- N1
  n2 <- N2
  r <- n2 / n1
  df <-n1+n2-2

  suppressWarnings(
    if (is.null(e)) {
      bf10=t2_BF10(tval, n1, r, model, location, scale, dff, hypothesis)
    } else {
      bf10=t2e_BF10(tval, n1, r, model,location, scale, dff, hypothesis, e)
    }
  )

  type = "Indepedent-samples t-test (equal variance)"
  analysis_h1 <- list(
    model = model,
    location = location,
    scale = scale
  )
  if (model == "t") {
    analysis_h1$dff <- dff
  }
  object=list(type=type,bf10=bf10,tval=tval,df=df,analysis_h1=analysis_h1,hypothesis=hypothesis,e=e,N1=N1,N2=N2)

  class(object) <- "BFvalue_t"

  return(object)
}


#' Bayes factor for a Bayesian correlation test
#'
#' Calculate the Bayes factor (BF10) for a correlation coefficient, either against a point NaN
#' or an interval null hypothesis. Supports default beta (\code{"d_beta"}), stretched beta (\code{"beta"}),
#' and normal-moment (\code{"moment"}) priors for the alternative hypothesis.
#'
#' @param r Observed correlation coefficient. Must be a numeric scalar between -1 and 1.
#' @param n Sample size. Must be a numeric scalar greater than 3.
#' @param k Parameter for the analysis default beta prior (\code{"d_beta"}) under the alternative hypothesis.
#' @param alpha Parameter for the analysis beta prior (\code{"beta"}) under the alternative hypothesis.
#' @param beta Parameter for the analysis beta prior (\code{"beta"}) under the alternative hypothesis.
#' @param h0 NaN value of the correlation. Must be a numeric scalar between -0.8 and 0.8.
#' @param hypothesis The hypothesis being tested: two-sided (\code{"!="}), right-sided (\code{">"}), or left-sided (\code{"<"}).
#' @param scale Scale parameter for the analysis normal-moment prior (\code{"moment"}). Must be > 0.
#' @param model Statistical model for the analysis prior: default beta (\code{"d_beta"}), beta (\code{"beta"}), or normal-moment (\code{"moment"}).
#' @param e Optional numeric vector specifying bounds for an interval null hypothesis. For \code{"!="}, must be two distinct finite values between -0.5 and 0.5. For \code{">"} or \code{"<"}, must satisfy additional bounds relative to \code{h0}.
#'
#' @return A list with class \code{"BFvalue_r"} containing:
#' \itemize{
#'   \item \code{type}: "correlation"
#'   \item \code{bf10}: Calculated Bayes factor BF10
#'   \item \code{h0}: NaN value of the correlation
#'   \item \code{r}: Observed correlation coefficient
#'   \item \code{n}: Sample size
#'   \item \code{analysis_h1}: List of analysis prior parameters including \code{model}, \code{k}, \code{alpha}, \code{beta}, \code{scale}
#'   \item \code{hypothesis}: Hypothesis being tested
#'   \item \code{e}: Interval bounds if specified
#' }
#'
#' @examples
#' BF10.cor(
#'   r = 0.3,
#'   n = 50,
#'   k = 1,
#'   alpha = 0.05,
#'   beta = 0.2,
#'   h0 = 0,
#'   hypothesis = "!=",
#'   scale = 1,
#'   model = "d_beta"
#' )
#'
#' BF10.cor(
#'   r = 0.25,
#'   n = 40,
#'   alpha = 2,
#'   beta = 3,
#'   h0 = 0,
#'   hypothesis = ">",
#'   scale = 1,
#'   model = "beta"
#' )
#'
#' @export
BF10.cor <- function(r, n, k, alpha, beta, h0, hypothesis,  scale,  model, e = NULL) {

  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < -0.8 || h0 > 0.8) {
    stop("argument [h0] NaN value of rho must be a single numeric scalar between -0.8 and 0.8")
  }
  location = h0

  # hypothesis
  if(hypothesis %in% c("<", "!=", ">") == FALSE){
    stop("argument [hypothesis] should be set to either `<`  (left-sided test),  `!=` (two-sided test) or `>` (right-sided test)")
  }

  # Equivlance test or not
  interval <- if (is.null(e)) 1 else 0

  if (!is.null(e)) {

    if (hypothesis == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(e) || length(e) != 2 || any(!is.finite(e)) || e[1] == e[2]) {
        stop("For hypothesis '!=', argument [e] must be a numeric vector of length 2 with two distinct finite values")
      }
      # Additional bounds checks
      if (min(e) < -0.5 || max(e) > 0.5) {
        stop("For hypothesis '!=', e must satisfy min(e) >= -0.5 and max(e) <= 0.5")
      }
      if ((h0 + min(e)) <= -1 || (h0 + min(e)) >= 1) {
        stop("For hypothesis '!=', h0 + min(e) must be between -1 and 1")
      }

    } else if (hypothesis == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e <= 0) {
        stop("For hypothesis '>', argument [e] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (e > 0.5) stop("For hypothesis '>', e must be <= 0.5")
      if ((h0 + e) >= 1) stop("For hypothesis '>', h0 + e must be < 1")

    } else if (hypothesis == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e >= 0) {
        stop("For hypothesis '<', argument [e] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (e < -0.5) stop("For hypothesis '<', e must be >= -0.5")
      if ((h0 + e) <= -1) stop("For hypothesis '<', h0 + e must be > -1")
    }

  }

  # data
  if (!is.numeric(r) || length(r) != 1 || !is.finite(r) || r < -1 || r > 1) {
    stop("argument [r] observed correlation must be a single numeric scalar between -1 and 1")
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 3) {
    stop("argument [n] sample size must be a single numeric scalar greater than 3")
  }

  # analysis prior model
  if (missing(model)) {
    stop("argument [model] for analysis prior must be one of `d_beta` (default stretched beta), `beta` (stretched beta), or `Moment` (normal-moment prior)")
  }
  # Analysis prior model validation
  if (!model %in% c("d_beta", "moment", "beta")) {
    stop("argument [model] for analysis prior must be one of `d_beta` (default stretched beta), `beta` (stretched beta), or `Moment` (normal-moment prior)")
  }

  # Model-specific checks
  if (model == "d_beta") {
    alpha=beta=scale=NULL
    # 'd_beta' requires k to be a single numeric scalar > 0
    if (!exists("k") || !is.numeric(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
      stop("For model 'd_beta', argument [k] must be a single numeric scalar > 0")
    }
  } else if (model == "beta") {
    k=scale=NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (!exists("alpha") || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For model 'beta', argument [alpha] must be a single numeric scalar > 0")
    }
    if (!exists("beta") || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For model 'beta', argument [beta] must be a single numeric scalar > 0")
    }
  } else if (model == "moment") {
    k=alpha=beta=NaN
    # 'Moment' requires scale to be numeric scalar > 0
    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For model 'Moment', argument [scale] must be a numeric scalar > 0")
    }
  }

  suppressWarnings(
    if (is.null(e)) {
      bf10=r_BF10(r, n, k, alpha, beta, h0, hypothesis, location, scale, 1, model)
    } else {
      bf10=re_BF10(r, n, k, alpha, beta, h0, hypothesis, location, scale, 1, model, e)
    }
  )

  type = "correlation"
  analysis_h1 <- list(
    model = model,
    k = k,
    alpha=alpha,
    beta=beta,
    scale=scale
  )
  object=list(type=type,bf10=bf10,h0=h0,r=r,n=n,analysis_h1=analysis_h1,hypothesis=hypothesis,e=e)

  class(object) <- "BFvalue_r"
  return(object)
}



#' Bayes Factor for a Bayesian F-Test
#'
#' Computes the Bayes factor (BF10) for an F-test, comparing a full model to a
#' reduced model under either an effect-size prior or a Moment prior.
#' Optionally, an interval null hypothesis can be specified.
#'
#' @param fval Numeric scalar. Observed F statistic (must be ≥ 0).
#' @param df1 Numeric scalar. Numerator degrees of freedom (must be > 0).
#' @param df2 Numeric scalar. Denominator degrees of freedom (must be > 0).
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior under
#'   the alternative hypothesis. For the Moment prior, this must be \eqn{\ge 3}.
#' @param rscale Numeric scalar. Scale parameter for the effect-size prior
#'   (only used when \code{model = "effectsize"}).
#' @param f_m Numeric scalar. Cohen's f effect-size parameter for the
#'   analysis prior.
#' @param model Character string specifying the analysis prior under the
#'   alternative hypothesis. Must be either \code{"effectsize"} or
#'   \code{"moment"}.
#' @param e Optional numeric scalar specifying an upper bound for an interval
#'   NaN hypothesis. If provided, must be > 0.
#'
#' @return A list of class \code{"BFvalue_f"} containing:
#' \itemize{
#'   \item \code{fval} Input F-value
#'   \item \code{df1, df2} Degrees of freedom
#'   \item \code{e} Interval bound (if specified)
#'   \item \code{analysis_h1} Details of the analysis prior
#'   \item \code{bf10} The computed Bayes factor
#' }
#'
#' @examples
#' BF10.f.test(
#'   fval = 4.5,
#'   df1 = 2,
#'   df2 = 12,
#'   dff = 12,
#'   rscale = 0.707,
#'   f_m = 0.1,
#'   model = "effectsize"
#' )
#'
#' @export
BF10.f.test <- function(fval, df1, df2, dff, rscale, f_m, model, e = NULL) {


  ## Check fval
  if (!is.numeric(fval) || length(fval) != 1 || !is.finite(fval) || fval < 0) {
    stop("argument [fval] F-value must be a numeric scalar greater than or equal to 0")
  }

  ## Check df1
  if (!is.numeric(df1) || length(df1) != 1 || !is.finite(df1) || df1 <= 0) {
    stop("argument [df1] numerator degrees of freedom must be a positive numeric scalar")
  }

  ## Check df2
  if (!is.numeric(df2) || length(df2) != 1 || !is.finite(df2) || df2 <= 0) {
    stop("argument [df2] denominator degrees of freedom must be a positive numeric scalar")
  }


  # analysis prior model
  if (missing(model)) {
    stop("argument [model] for analysis prior should be set to either `effectsize`, or `Moment`")
  }
  if(model %in% c("effectsize","moment") == FALSE){
    stop("argument [model] for analysis prior should be set to either `effectsize`, or `Moment`")
  }

  if (model =="effectsize"){
    if (!is.numeric(rscale) || length(rscale) != 1 || !is.finite(rscale) || rscale <= 0) {
      stop("argument [rscale] scale parameter must be a positive numeric scalar")
    }
  }

  if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
    stop("argument [dff] degrees of freedom  for analysis prior must be a positive numeric scalar when model='t-distribution'")
  }

  if (!is.numeric(f_m) || length(f_m) != 1 || !is.finite(f_m) || f_m <= 0) {
    stop("argument [f_m] Cohen's f  for analysis prior must be a positive numeric scalar")
  }

  if (model == "moment"){
    rscale=NULL
    if (dff < 3) {
      stop("argument [dff] degrees of freedom for Moment analysis prior must be at least 3")
    }
  }

  if (!is.null(e)) {
    if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e <= 0) {
      stop("argument [e] interval bound must be a positive numeric scalar when specified")
    }
  }

  q <- df1
  m <- df1 + df2

  bf10= suppressWarnings(
    if (is.null(e)) {
      F_BF(fval, q, m, dff, rscale, f_m, model)
    } else {
      Fe_BF(fval, q, m, dff, rscale, f_m, model, e)
    }
  )

  type = "Regression/ANOVA"
  analysis_h1 <- list(
    model = model,
    rscale = rscale,
    f_m = f_m,
    dff=dff
  )
  object <- list(
    fval=fval,
    type = type,
    e = e,
    analysis_h1 = analysis_h1,
    df1=df1,
    df2=df2,
    bf10=bf10
  )
  class(object) <- "BFvalue_f"
  return(object)

}

#' Bayes Factor for a Bayesian One-Proportion Test
#'
#' Calculate the Bayes factor (BF10) for a single-proportion test, either against a point NaN
#' or an interval null hypothesis.
#'
#' @param x Observed number of successes (non-negative integer scalar, must be \eqn{\le n}).
#' @param n Sample size (positive integer scalar).
#' @param alpha Shape parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{model = "beta"}).
#' @param beta Shape parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{model = "beta"}).
#' @param h0 Null proportion value (numeric scalar between 0.1 and 0.9).
#' @param scale Scale parameter for the analysis prior (only used if \code{model = "moment"}).
#' @param model Statistical model for the analysis prior under the alternative hypothesis:
#'   \code{"beta"} (stretched beta) or \code{"moment"} (normal-moment prior).
#' @param hypothesis Hypothesis being tested: two-sided (\code{"!="}), right-sided (\code{">"}),
#'   or left-sided (\code{"<"}).
#' @param e Optional numeric vector specifying bounds for an interval null; used if interval BF is calculated.
#'
#' @return An object of class \code{"BFvalue_bin"} containing:
#'   \itemize{
#'     \item \code{bf10}: Bayes factor in favor of the alternative hypothesis.
#'     \item \code{type}: Test type ("one-proportion").
#'     \item \code{x}: Number of successes.
#'     \item \code{n}: Sample size.
#'     \item \code{h0}: Null proportion value.
#'     \item \code{analysis_h1}: List describing the analysis prior (\code{model}, \code{alpha}, \code{beta}, \code{scale}).
#'     \item \code{hypothesis}: Tested hypothesis.
#'     \item \code{e}: interval null bounds (if specified).
#'   }
#'
#' @examples
#' BF10.bin.test(
#'   x = 12,
#'   n = 50,
#'   alpha = 2,
#'   beta = 3,
#'   h0 = 0.5,
#'   scale = 1,
#'   model = "beta",
#'   hypothesis = "!="
#' )
#'
#' BF10.bin.test(
#'   x = 8,
#'   n = 20,
#'   h0 = 0.5,
#'   scale = 0.5,
#'   model = "moment",
#'   hypothesis = ">"
#' )
#'
#' @export

BF10.bin.test <- function(x, n, alpha, beta, h0, scale, model, hypothesis, e = NULL) {
  # Check n
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0 || n != floor(n)) {
    stop("argument [n] sample size must be a positive integer scalar")
  }

  # Check x
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x != floor(x)) {
    stop("argument [x] number of successes must be a non-negative integer scalar")
  }

  # Check relation
  if (x > n) {
    stop("argument [x] number of successes cannot exceed total sample size [n]")
  }
  # mode
  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < .1 || h0 > 0.9) {
    stop("argument [h0] NaN value of proportion must be a single numeric scalar between .1 and 0.9")
  }
  # hypothesis
  if(hypothesis %in% c("<", "!=", ">") == FALSE){
    stop("argument [hypothesis] should be set to either `<`  (left-sided test),  `!=` (two-sided test) or `>` (right-sided test)")
  }
  # Equivlance test or not
  interval <- if (is.null(e)) 1 else 0

  if (!is.null(e)) {

    if (hypothesis == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(e) || length(e) != 2 || any(!is.finite(e)) || e[1] == e[2]) {
        stop("For hypothesis '!=', argument [e] must be a numeric vector of length 2 with two distinct finite values")
      }
      # Additional bounds checks
      if (min(e) < -0.5 || max(e) > 0.5) {
        stop("For hypothesis '!=', e must satisfy min(e) >= -0.5 and max(e) <= 0.5")
      }
      if ((h0 + min(e)) <= 0 || (h0 + min(e)) >= 1) {
        stop("For hypothesis '!=', h0 + min(e) must be between 0 and 1")
      }

    } else if (hypothesis == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e <= 0) {
        stop("For hypothesis '>', argument [e] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (e > 0.5) stop("For hypothesis '>', e must be <= 0.5")
      if ((h0 + e) >= 1) stop("For hypothesis '>', h0 + e must be < 1")

    } else if (hypothesis == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(e) || length(e) != 1 || !is.finite(e) || e >= 0) {
        stop("For hypothesis '<', argument [e] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (e < -0.5) stop("For hypothesis '<', e must be >= -0.5")
      if ((h0 + e) <= -1) stop("For hypothesis '<', h0 + e must be > 0")
    }

  }
  # analysis prior model
  if (missing(model)) {
    stop("argument [model] for analysis prior must be one of `beta`, or `Moment` (normal-moment prior)")
  }

  # Analysis prior model validation
  if (!model %in% c("moment", "beta")) {
    stop("argument [model] for analysis prior must be one of `beta` , or `Moment` (normal-moment prior)")
  }

  # Model-specific checks
  if (model == "beta") {
    scale=NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (!exists("alpha") || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For model 'beta', argument [alpha] must be a single numeric scalar > 0")
    }
    if (!exists("beta") || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For model 'beta', argument [beta] must be a single numeric scalar > 0")
    }
  } else if (model == "moment") {
    alpha=beta=NaN
    # 'Moment' requires scale to be numeric scalar > 0
    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For model 'Moment', argument [scale] must be a numeric scalar > 0")
    }
  }


  bf10=suppressWarnings(
    if (is.null(e)) {
      bin_BF(x, n, alpha, beta, h0, scale, model, hypothesis)
    } else {
      bin_e_BF(x, n, alpha, beta, h0, scale, model, hypothesis, e)
    }
  )
  type = "one-proportion"
  analysis_h1 <- list(
    model = model,
    alpha=alpha,
    beta=beta,
    scale=scale
  )
  object=list(type=type,bf10=bf10,h0=h0,x=x,n=n,analysis_h1=analysis_h1,hypothesis=hypothesis,e=e)

  class(object) <- "BFvalue_bin"
  return(object)
}


#' Bayes factor for comparing two proportions
#'
#' Compute the Bayes factor (BF10) for a Bayesian test of two proportions.
#'
#' @param a0 Alpha parameter of the Beta prior under the null hypothesis.
#' @param b0 Beta parameter of the Beta prior under the null hypothesis.
#' @param a1 Alpha parameter of the Beta prior for group 1 under the alternative hypothesis.
#' @param b1 Beta parameter of the Beta prior for group 1 under the alternative hypothesis.
#' @param a2 Alpha parameter of the Beta prior for group 2 under the alternative hypothesis.
#' @param b2 Beta parameter of the Beta prior for group 2 under the alternative hypothesis.
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param x1 Number of successes observed in group 1.
#' @param x2 Number of successes observed in group 2.
#'
#' @return A list of class \code{BFvalue_2p} containing:
#' \itemize{
#'   \item \code{type}: the string "Two-proportions".
#'   \item \code{analysis_h0}: list with \code{a} and \code{b} for the null prior.
#'   \item \code{analysis_h1_theta_1}: list with \code{a} and \code{b} for group 1 prior under H1.
#'   \item \code{analysis_h1_theta_2}: list with \code{a} and \code{b} for group 2 prior under H1.
#'   \item \code{bf10}: the computed Bayes factor (BF10).
#'   \item \code{n1}, \code{x1}, \code{n2}, \code{x2}: the input sample sizes and observed successes.
#' }
#'
#' @examples
#' BF10.props(
#'   a0 = 2, b0 = 3,
#'   a1 = 2, b1 = 3,
#'   a2 = 2, b2 = 3,
#'   n1 = 50, n2 = 60,
#'   x1 = 25, x2 = 30
#' )
#'
#' @export
BF10.props <- function(a0, b0, a1, b1, a2, b2, n1, n2, x1, x2) {


  # null hypothesis
  # Check a0 (alpha)
  if (!is.numeric(a0) || length(a0) != 1 || a0 <= 0) {
    stop("arg [a0] alpha for the Beta analysis prior under the null (\u03b80) must be a positive numeric scalar (> 0).")
  }

  # Check b0 (beta)
  if (!is.numeric(b0) || length(b0) != 1 || b0 <= 0) {
    stop("arg [b0] beta for the Beta analysis prior under the null (\u03b80) must be a positive numeric scalar (> 0).")
  }

  # alternative hypothesis \u03b81
  # Check a1 (alpha under the alternative)
  if (!is.numeric(a1) || length(a1) != 1 || a1 <= 0) {
    stop("arg [a1] alpha for the Beta analysis prior under the alternative (\u03b81) must be a positive numeric scalar (> 0).")
  }

  # Check b1 (beta under the alternative)
  if (!is.numeric(b1) || length(b1) != 1 || b1 <= 0) {
    stop("arg [b1] beta for the Beta analysis prior under the alternative (\u03b81) must be a positive numeric scalar (> 0).")
  }

  # alternative hypothesis \u03b82
  # Check a2 (alpha under the alternative)
  if (!is.numeric(a2) || length(a2) != 1 || a2 <= 0) {
    stop("arg [a2] alpha for the Beta analysis prior under the alternative (\u03b82) must be a positive numeric scalar (> 0).")
  }

  # Check b2 (beta under the alternative)
  if (!is.numeric(b2) || length(b2) != 1 || b2 <= 0) {
    stop("arg [b2] beta for the Beta analysis prior under the alternative (\u03b82) must be a positive numeric scalar (> 0).")
  }

  # sample sizes
  if (!is.numeric(n1) || length(n1) != 1 || n1 %% 1 != 0 || n1 <= 0) {
    stop("arg [n1] sample size for group 1 must be a positive numeric scalar integer (> 0).")
  }
  if (!is.numeric(n2) || length(n2) != 1 || n2 %% 1 != 0 || n2 <= 0) {
    stop("arg [n2] sample size for group 2 must be a positive numeric scalar integer (> 0).")
  }

  # observed successes
  if (!is.numeric(x1) || length(x1) != 1 || x1 %% 1 != 0 || x1 < 0) {
    stop("arg [x1] for group 1 must be a non-negative numeric scalar integer (\u2265 0).")
  }
  if (!is.numeric(x2) || length(x2) != 1 || x2 %% 1 != 0 || x2 < 0) {
    stop("arg [x2] for group 2 must be a non-negative numeric scalar integer (\u2265 0).")
  }

  if (x1 > n1) {
    stop("arg [x1] number of successes in group 1 cannot exceed n1 (sample size).")
  }
  if (x2 > n2) {
    stop("arg [x2] number of successes in group 2 cannot exceed n2 (sample size).")
  }

  bf10=BF10_p2(a0, b0, a1, b1, a2, b2, n1, n2, x1, x2)
  type = "Two-proportions"
  analysis_h0 <- list(
    a = a0,
    b = b0
  )
  analysis_h1_theta_1 <- list(
    a = a1,
    b = b1
  )
  analysis_h1_theta_2 <- list(
    a = a2,
    b = b2
  )


  object <- list(
    type = type,
    analysis_h0=analysis_h0,
    analysis_h1_theta_1= analysis_h1_theta_1,
    analysis_h1_theta_2=analysis_h1_theta_2,
    bf10=bf10,
    n1=n1,
    x1=x1,
    n2=n2,
    x2=x2
  )
  class(object) <- "BFvalue_2p"
  return(object)

}

