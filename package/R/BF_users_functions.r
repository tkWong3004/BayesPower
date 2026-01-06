#' Sample Size Determination for One-Sample Bayesian t-Test
#'
#' Performs sample size determination or calculates the probabilities of obtaining compelling or misleading evidence
#' for a one-sample Bayesian t-test. Can handle both point-null and interval-null Bayes factors.
#'
#' @param alternative Character. The direction of the alternative hypothesis : two-sided (\code{"two.sided"}), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector. Bounds for an interval null hypothesis.
#' @param prior_analysis Character. The analysis prior under the alternative hypothesis:
#'   \code{"Normal"}, \code{"Moment"} (normal-moment prior), or \code{"t-distribution"}.
#' @param location Numeric scaler. Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Numeric scaler. Scale parameter for the analysis prior under the alternative hypothesis (must be > 0).
#' @param dff Numeric scaler. Degrees of freedom for the analysis prior under the alternative hypothesis (required if \code{prior_analysis = "t-distribution"}).
#' @param prior_design Optional character. The design prior under the alternative hypothesis:
#'   \code{"Normal"}, \code{"Moment"}, \code{"t-distribution"}, or \code{"Point"}.
#' @param location_d Numeric scaler. Location parameter for the design prior under the alternative hypothesis.
#' @param scale_d Numeric scaler. Scale parameter for the design prior under the alternative hypothesis.
#' @param dff_d Numeric scaler. Degrees of freedom for the design prior under the alternative hypothesis (required if \code{prior_design = "t-distribution"}).
#' @param N Numeric scaler. Sample size.
#' @param type_rate Character. Either \code{"positive"} (controls true/false positive rates) or \code{"negative"} (controls true/false negative rates).
#' @param true_rate Numeric scaler. Target true positive or negative rate (between 0.6 and 0.999).
#' @param false_rate Numeric scaler. Target false positive or false negative rate (between 0.001 and 0.1).
#' @param threshold Numeric scaler. Threshold of compelling evidence (must be > 1).
#' @param plot_power Logical. If \code{TRUE}, plots power curve.
#' @param plot_rel Logical. If \code{TRUE}, plots the relationship between the BF and data.
#'
#' @return An object of class \code{BFpower_t} (a list) containing:
#' \describe{
#'   \item{type}{Character, always "One-sample t-test".}
#'   \item{alternative}{Character, the direction of the alternative hypothesis.}
#'   \item{ROPE}{Optional numeric vector for interval null bounds.}
#'   \item{analysis_h1}{List with analysis prior parameters: \code{prior_analysis}, \code{location}, \code{scale}, and optionally \code{dff}.}
#'   \item{design_h1}{List with design prior parameters: \code{prior_analysis}, \code{location}, \code{scale}, and optionally \code{dff} (or \code{NULL} if not provided).}
#'   \item{results}{Data frame of probabilities: compelling/misleading evidence, or \code{NaN} if calculation fails.}
#'   \item{threshold}{Numeric, threshold of compelling evidence.}
#'   \item{plot_power}{Logical, whether to plot the power curve.}
#'   \item{plot_rel}{Logical, whether the relationship between the BF and t-value is plotted.}
#' }
#' @details
#' \strong{1. Sample size determination mode (when \code{N = NULL}):}
#'
#' If no sample size is provided, the function determines the minimum sample size. In this mode, the user
#' must supply the following arguments:
#' \itemize{
#'   \item \code{type_rate} - either \code{"positive"} to control true/false positive rates,
#'         or \code{"negative"} to control true/false negative rates.
#'   \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#'   \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#'   \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be > 1).
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
#' and \code{false_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Analysis Priors:}
#'
#' The analysis prior specifies the prior distribution of the effect under the
#' alternative hypothesis. The user must provide:
#' \itemize{
#'   \item \code{prior_analysis} - the type of prior: \code{"Normal"}, \code{"Moment"} (normal-moment prior), or \code{"t-distribution"}.
#'   \item \code{location} - the mean or location of the prior.
#'   \item \code{scale} - the standard deviation or scale (must be positive).
#'   \item \code{dff} - degrees of freedom (required if \code{prior_analysis = "t-distribution"}).
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' A design prior can be supplied to reflect uncertainty about the effect size
#' during study planning. If provided, the following must be supplied:
#' \itemize{
#'   \item \code{prior_design} - the type of design prior: \code{"Normal"}, \code{"Moment"}, \code{"t-distribution"}, or \code{"Point"}.
#'   \item \code{location_d} - the location of the design prior.
#'   \item \code{scale_d} - the scale parameter (positive for all models except \code{"Point"}).
#'   \item \code{dff_d} - degrees of freedom for \code{"t-distribution"} design priors.
#' }
#'
#'
#' \strong{interval null Hypothesis:}
#'
#' The argument \code{ROPE} specifies the bounds of an interval null hypothesis.
#' If \code{ROPE} is provided, the function evaluates the Bayes factor for an interval
#' null hypothesis. For a point-null hypothesis, \code{ROPE} should be left as \code{NULL}.
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling
#' evidence as a function of sample size. If \code{plot_rel = TRUE}, the relationship betwwen the BF and data is plotted.
#' @examples
#'BFpower.ttest.OneSample(
#'  alternative = "two.sided",
#'  threshold = 3,
#'  true_rate = 0.8,
#'  false_rate = 0.05,
#'  prior_analysis = "t-distribution",
#'  location = 0,
#'  scale = 0.707,
#'  dff = 1,
#'  plot_power = TRUE,
#'  plot_rel = TRUE
#')
#' @export
BFpower.ttest.OneSample <- function(
    alternative, ROPE=NULL,
    prior_analysis, location, scale, dff,
    prior_design=NULL, location_d, scale_d, dff_d,
    N=NULL,
    type_rate = "positive", true_rate, false_rate , threshold, plot_power = FALSE,plot_rel=FALSE
)  {
  # mode
  if ( is.null(N)) mode_bf=1 else mode_bf = 0

  # sample size
  if (mode_bf == 0) {
    # Check that N is a positive numeric scalar
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 0) {
      stop("Argument [N] sample size must be a positive numeric scalar when mode_bf = 0")
    }
  }else {N=2}

  # alternative
  if(alternative %in% c("two.sided", "less", "greater") == FALSE){
    stop("Argument [alternative] should be set to either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test)")
  }

  alternative <- switch(alternative,
                        "two.sided" = "!=",
                        "less"      = "<",
                        "greater"   = ">"
  )


  # Equivlance test or not
  interval <- if (is.null(ROPE)) 1 else 0

  if (!is.null(ROPE)) {

    if (alternative == "!=") {
      # e must be a numeric vector of length 2, both positive
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative '!=', Argument [ROPE] must be a numeric vector of length 2 with two distinct finite values")
      }

    }

    if (alternative == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative '>', Argument [ROPE] must be a numeric scalar > 0")
      }
    }

    if (alternative == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative '<', Argument [ROPE] must be a numeric scalar < 0")
      }
    }

  }


  # analysis prior prior_analysis
  if (missing(prior_analysis)) {
    stop("Argument [prior_analysis] for analysis prior should be set to either `Normal`, `Moment` (normal-moment prior) or `t-distribution`")
  }
  if(prior_analysis %in% c("Normal","Moment","t-distribution") == FALSE){
    stop("Argument [prior_analysis] for analysis prior should be set to either `Normal`,  `Moment` (normal-moment prior) or `t-distribution` ")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("Argument [location] for analysis prior  must be a numeric scalar")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale)||scale<=0) {
    stop("Argument [scale] for analysis prior must be a positive numeric scalar (i.e., scale > 0)")
  }
  if (prior_analysis == "t-distribution") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("Argument [dff] degrees of freedom  for analysis prior must be a positive numeric scalar when prior_analysis='t-distribution'")
    }
  }else{
    dff = 0
  }

  # design prior

  if (!is.null(prior_design)) {

    de_an_prior <- 0

    # Validate prior_design
    if (!(prior_design %in% c("Normal", "Moment", "t-distribution", "Point"))) {
      stop("Argument [prior_design] for design prior must be either `Normal`, `Moment`, `t-distribution`, or `Point`")
    }

    # Validate location_d in one line
    if (!is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d))
      stop("Argument [location_d] for design prior must be a numeric scalar")

    # Validate scale_d for prior_analysiss that require it
    if (prior_design %in% c("Normal", "Moment", "t-distribution")) {
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0)
        stop("Argument [scale_d] for design prior must be a positive numeric scalar (i.e., scale_d > 0)")
    }

    # Validate dff_d only when prior_design = t-distribution
    if (prior_design == "t-distribution") {
      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0)
        stop("Argument [dff_d] degrees of freedom for design prior must be a positive numeric scalar when prior_design='t-distribution'")
    } else {
      dff_d <- 0
    }

  } else {
    de_an_prior <- 1
  }



  # desired power and strength of evidence
  if (mode_bf==1){
    if (!(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (type_rate,
                    "positive" = "h1",
                    "negative" = "h0"
    )
    if (!is.numeric(true_rate) || length(true_rate) != 1 ||
        !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999){
      stop("Argument [true_rate] targeted true positive or negative rate must be a numeric scalar strictly greater than 0.6 and less than 0.999")
    }
    target = true_rate
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("Argument [false_rate] must be a numeric scalar strictly greater than 0.001 and less than 0.1")
    }

    alpha = false_rate
    if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold <= 1) {
      stop("Argument [threshold] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=alpha=0
  }
  ####


  # Call appropriate table function with error handling
  tryCatch(
    {
      if (interval == 1) {
        results = suppressWarnings(t1_Table(threshold, target, prior_analysis, location, scale, dff, alternative,
                                            prior_design, location_d, scale_d, dff_d, de_an_prior, N, mode_bf, alpha, direct))
      } else {
        results = suppressWarnings(t1e_table(threshold,target,prior_analysis,location,scale,dff, alternative,ROPE ,
                                             prior_design,scale_d,dff_d, de_an_prior,N,mode_bf,location_d ,alpha,direct ))
      }

    },
    error = function(err) {
      message("Required sample size > 10,000")
      stop(NaN)
    }
  )
  type = "One-sample t-test"
  analysis_h1 <- list(
    prior = prior_analysis,
    location = location,
    scale = scale
  )
  if (prior_analysis == "t-distribution") {
    analysis_h1$dff <- dff
  }

  if (!is.null(prior_design)) {

    # Base fields always included
    design_h1 <- list(
      prior    = prior_design,
      location = location_d,
      scale    = scale_d
    )

    # Only add dff if prior_analysis is t-distribution
    if (prior_design == "t-distribution") {
      design_h1$dff <- dff_d
    }

  } else {

    # prior_design is NULL > fill all fields with NULL
    design_h1 <- list(
      prior_analysis    = NULL,
      location = NULL,
      scale    = NULL,
      dff      = NULL
    )
  }


  object <- list(
    type = "One-sample t-test",
    alternative = alternative,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    threshold = threshold,
    mode_bf = mode_bf,
    plot_power=plot_power,
    plot_rel=plot_rel
  )
  class(object) <- "BFpower_t"
  plot(object)
  return(object)
}
#' Sample Size Determination for Two-Sample Bayesian t-Test
#'
#' Perform sample size determination or calculate the probabilities of obtaining
#' compelling or misleading evidence for a two-sample Bayesian t-test.
#' Supports point-null and interval-null hypotheses, and allows specifying
#' analysis and design priors.
#'
#' @param alternative Character. The direction of the alternative hypothesis: two-sided (\code{"two.sided"}),
#'   right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric. Bounds for an interval null:
#'   - For \code{hypothesis = "two.sided"}, must be a numeric vector of length 2 with distinct finite values.
#'   - For \code{"greater"}, must be a single numeric scalar > 0.
#'   - For \code{"less"}, must be a single numeric scalar < 0.
#' @param threshold Numeric scalar. Threshold for compelling evidence (must be > 1).
#' @param true_rate Numeric scalar. Target true positive or negative rate .
#' @param false_rate Numeric scalar. Target false positive or negative rate .
#' @param prior_analysis Character. Analysis prior under the alternative hypothesis:
#'   \code{"Normal"}, \code{"Moment"}, or \code{"t-distribution"}.
#' @param location Numeric scalar. Location parameter for the analysis prior.
#' @param scale Numeric scalar > 0. Scale parameter for the analysis prior.
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior (required if prior_analysis = \code{"t-distribution"}; ignored otherwise).
#' @param prior_design Optional character. Design prior under the alternative:
#'   \code{"Normal"}, \code{"Moment"}, \code{"t-distribution"}, or \code{"Point"}.
#' @param location_d Numeric scalar. Location parameter for the design prior.
#' @param scale_d Numeric scalar > 0. Scale parameter for the design prior.
#' @param dff_d Numeric scalar. Degrees of freedom for the design prior (required if \code{prior_design = "t-distribution"}; ignored otherwise).
#' @param N1 Sample size for group 1 (used if \code{r = NULL}).
#' @param N2 Sample size for group 2 (used if \code{r = NULL}).
#' @param r Optional numeric scalar. Ratio of sample size \code{N2 / N1} (used if \code{N1} and \code{N2} are NULL).
#' @param type_rate Character, either \code{"positive"} or \code{"negative"}; determines whether to control
#'   true/false positive or true/false negative rates .
#' @param plot_power Logical. If \code{TRUE}, a plot of the power or probability of compelling evidence is generated.
#' @param plot_rel Logical. Whether the relationship between the BF and data is plotted..
#'
#' @details
#' \strong{1. Sample size determination mode (when \code{N1 = NULL} and \code{N2 = NULL}, but \code{r} is provided):}
#'
#' If no sample sizes are provided, the function calculates the minimum required sample sizes for both groups. In this mode, the user
#' must supply:
#' \itemize{
#'   \item \code{type_rate} - either \code{"positive"} to control true/false positive rates,
#'         or \code{"negative"} to control true/false negative rates.
#'   \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#'   \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#'   \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be > 1).
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
#' the arguments \code{type_rate}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Analysis Priors:}
#'
#' The analysis prior specifies the prior distribution of the effect under the alternative hypothesis. The user must provide:
#' \itemize{
#'   \item \code{prior_analysis} - the type of prior: \code{"Normal"}, \code{"Moment"} (normal-moment prior), or \code{"t-distribution"}.
#'   \item \code{location} - the mean or location of the prior.
#'   \item \code{scale} - the standard deviation or scale (must be positive).
#'   \item \code{dff} - degrees of freedom (required if \code{prior_analysis = "t-distribution"}).
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' A design prior can be supplied to reflect uncertainty about the effect size during study planning. If provided, the following must be supplied:
#' \itemize{
#'   \item \code{prior_design} - the type of design prior: \code{"Normal"}, \code{"Moment"}, \code{"t-distribution"}, or \code{"Point"}.
#'   \item \code{location_d} - the location of the design prior.
#'   \item \code{scale_d} - the scale parameter (positive for all models except \code{"Point"}).
#'   \item \code{dff_d} - degrees of freedom for \code{"t-distribution"} design priors.
#' }
#'
#' \strong{interval null Hypothesis:}
#'
#' The argument \code{ROPE} specifies the bounds of an interval null hypothesis.
#' If \code{ROPE} is provided, the function evaluates the Bayes factor for an interval
#' null hypothesis. For a point-null hypothesis, \code{ROPE} should be left as \code{NULL}.
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling
#' evidence as a function of the sample sizes. If \code{plot_rel = TRUE}, the relationship between BF and data is plotted.
#'
#' @return An object of class \code{BFpower_t} containing:
#' \describe{
#'   \item{type}{Character string describing the test type.}
#'   \item{alternative}{alternative hypothesis (\code{"two.sided"}, \code{"greater"}, or \code{"less"}).}
#'   \item{ROPE}{Interval bounds under the null used, if any.}
#'   \item{analysis_h1}{List of analysis prior parameters used for the alternative hypothesis.}
#'   \item{design_h1}{List of design prior parameters used for the alternative hypothesis.}
#'   \item{results}{Data frame with probabilities of compelling/misleading evidence.}
#'   \item{threshold}{Threshold of compelling evidence.}
#'   \item{plot_power}{Logical flag for plotting power.}
#'   \item{plot_rel}{Logical flag for plotting the relationship between BF and t-value.}
#' }
#'
#' @examples
#'BFpower.ttest.TwoSample(
#'  alternative = "two.sided",
#'  ROPE = c(-0.36, 0.36),
#'  threshold = 3,
#'  true_rate = 0.8,
#'  false_rate = 0.05,
#'  prior_analysis = "Normal",
#'  location = -0.23,
#'  scale = 0.2,
#'  dff = 1,
#'  type_rate = "negative",
#'  plot_power = TRUE,
#'  plot_rel = TRUE,
#'  r = 1)
#' @export
BFpower.ttest.TwoSample <- function(alternative , ROPE = NULL,
                                    threshold , true_rate , false_rate ,
                                    prior_analysis , location , scale , dff ,
                                    prior_design = NULL, location_d , scale_d , dff_d ,
                                    N1 = NULL, N2 = NULL, r=NULL, type_rate="positive",plot_power=FALSE,plot_rel=FALSE) {
  ## ------------------------------
  ## CHECKING N1, N2, r CONSISTENCY
  ## ------------------------------

  if (is.null(N1) && is.null(N2) && is.null(r)) {
    stop(
      "argument [r] ratio of N2/N1 must be specified for sample size calculation\n",
      "or argument [N2] sample size for group 2 and [N1] sample size for group 1 for power calculation"
    )
  }

  # Case A: r provided > N1 and N2 must be NULL
  if (!is.null(r)) {
    mode_bf=1

    # r must be numeric scalar > 0
    if (!is.numeric(r) || length(r) != 1 || !is.finite(r) || r <= 0) {
      stop("argument [r] ratio of sample size for group 2 over 1 must be a positive numeric scalar")
    }

    if (!is.null(N1) || !is.null(N2)) {
      stop("If argument [r] is provided, both N1 and N2 must be NULL for sample size determination")
    }

  }

  # Case B: r is NULL > N1 and N2 must both be valid numeric scalars
  if (is.null(r)) {
    mode_bf=0
    if (is.null(N1) || is.null(N2)) {
      stop("If 'r' is NULL, both N1 and N2 must be provided")
    }

    if (!is.numeric(N1) || length(N1) != 1 || !is.finite(N1) || N1 <= 0) {
      stop("argument [N1] sample size for group 1 must be a positive numeric scalar")
    }
    if (!is.numeric(N2) || length(N2) != 1 || !is.finite(N2) || N2 <= 0) {
      stop("argument [N2] sample size for group 2 must be a positive numeric scalar")
    }
  }

  # alternative
  if(alternative %in% c("two.sided", "less", "greater") == FALSE){
    stop("Argument [alternative] should be set to either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test)")
  }

  alternative <- switch(alternative,
                        "two.sided" = "!=",
                        "less"      = "<",
                        "greater"   = ">"
  )

  # Equivlance test or not
  interval <- if (is.null(ROPE)) 1 else 0

  if (!is.null(ROPE)) {

    if (alternative == "!=") {
      # e must be a numeric vector of length 2, both positive
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative '!=', Argument [ROPE] must be a numeric vector of length 2 with two distinct finite values")
      }

    }

    if (alternative == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative '>', Argument [ROPE] must be a numeric scalar > 0")
      }
    }

    if (alternative == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative '<', Argument [ROPE] must be a numeric scalar < 0")
      }
    }

  }


  # analysis prior prior_analysis
  if (missing(prior_analysis)) {
    stop("argument [prior_analysis] for analysis prior should be set to either `Normal`, `Moment` (normal-moment prior) or `t-distribution`")
  }

  if(prior_analysis %in% c("Normal","Moment","t-distribution") == FALSE){
    stop("argument [prior_analysis] for analysis prior should be set to either `Normal`,  `Moment` (normal-moment prior) or `t-distribution` ")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("argument [location] for analysis prior  must be a numeric scalar")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale)||scale<=0) {
    stop("argument [scale] for analysis prior must be a positive numeric scalar (i.e., scale > 0)")
  }
  if (prior_analysis == "t-distribution") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("argument [dff] degrees of freedom  for analysis prior must be a positive numeric scalar when prior_analysis='t-distribution'")
    }
  }else{
    dff = 0
  }
  # design prior

  if (!is.null(prior_design)) {

    de_an_prior <- 0

    # Validate prior_design
    if (!(prior_design %in% c("Normal", "Moment", "t-distribution", "Point"))) {
      stop("argument [prior_design] for design prior must be either `Normal`, `Moment`, `t-distribution`, or `Point`")
    }

    # Validate location_d in one line
    if (!is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d))
      stop("argument [location_d] for design prior must be a numeric scalar")

    # Validate scale_d for prior_analysiss that require it
    if (prior_design %in% c("Normal", "Moment", "t-distribution")) {
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0)
        stop("argument [scale_d] for design prior must be a positive numeric scalar (i.e., scale_d > 0)")
    }

    # Validate dff_d only when prior_design = t-distribution
    if (prior_design == "t-distribution") {
      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0)
        stop("argument [dff_d] degrees of freedom for design prior must be a positive numeric scalar when prior_design='t-distribution'")
    } else {
      dff_d <- 0
    }

  } else {
    de_an_prior <- 1
  }


  # desired power and strength of evidence
  if (mode_bf==1){
    if (!(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (type_rate,
                    "positive" = "h1",
                    "negative" = "h0"
    )
    if (!is.numeric(true_rate) || length(true_rate) != 1 ||
        !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999){
      stop("Argument [true_rate] targeted true positive or negative rate must be a numeric scalar strictly greater than 0.6 and less than 0.999")
    }
    target = true_rate
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("Argument [false_rate] must be a numeric scalar strictly greater than 0.001 and less than 0.1")
    }

    alpha = false_rate
    if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold <= 1) {
      stop("Argument [threshold] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=alpha=0
  }

  tryCatch(
    suppressWarnings({
      if (interval == 1) {
        results=t2_Table(threshold, r, target, prior_analysis, location, scale, dff, alternative,
                         prior_design, location_d, scale_d, dff_d, de_an_prior, N1, N2, mode_bf, alpha, direct)
      } else {
        results=t2e_table(threshold, r, target, prior_analysis,location, scale, dff, alternative, ROPE,
                          prior_design,location_d, scale_d, dff_d, de_an_prior, mode_bf, N1, N2, alpha, direct)
      }
    }),
    error = function(err) {
      message("Required sample size > 10,000")
      stop(NaN)
    }
  )


  type = "Indepedent-samples t-test (equal variance)"
  analysis_h1 <- list(
    prior = prior_analysis,
    location = location,
    scale = scale
  )
  if (prior_analysis == "t-distribution") {
    analysis_h1$dff <- dff
  }

  if (!is.null(prior_design)) {

    # Base fields always included
    design_h1 <- list(
      prior    = prior_design,
      location = location_d,
      scale    = scale_d
    )

    # Only add dff if prior_analysis is t-distribution
    if (prior_design == "t-distribution") {
      design_h1$dff <- dff_d
    }

  } else {

    # prior_design is NULL > fill all fields with NULL
    design_h1 <- list(
      prior_analysis    = NULL,
      location = NULL,
      scale    = NULL,
      dff      = NULL
    )
  }


  object <- list(
    type = type,
    alternative = alternative,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    threshold = threshold,
    mode_bf = mode_bf,
    plot_power=plot_power,
    plot_rel=plot_rel
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
#' @param alternative The direction of the alternative hypothesis being tested: two-sided (\code{"two.sided"}), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param h0 Null value of the correlation; must be a numeric scalar between -0.8 and 0.8.
#' @param ROPE Optional numeric vector or scalar specifying bounds for an interval null; used if interval Bayes factor is calculated.
#' @param threshold Threshold of compelling evidence (numeric scalar > 1).
#' @param true_rate Targeted true positive rate (if \code{positive = "positive"}) or true negative rate (if \code{positive = "negative"}).
#' @param false_rate Targeted false positive rate (if \code{positive = "positive"}) or false negative rate (if \code{positive = "negative"}).
#' @param prior_analysis Analysis prior under the alternative hypothesis: default beta (\code{"d_beta"}), beta (\code{"beta"}), or normal moment (\code{"Moment"}).
#' @param k Parameter for the default beta prior (\code{"d_beta"}).
#' @param alpha Parameter for the beta prior (\code{"beta"}).
#' @param beta Parameter for the beta prior (\code{"beta"}).
#' @param scale Scale parameter for the normal moment prior (\code{"Moment"}).
#' @param prior_design Design prior  under the alternative hypothesis: default beta (\code{"d_beta"}), beta (\code{"beta"}), normal moment (\code{"Moment"}), or point (\code{"Point"}).
#' @param alpha_d Parameter for the design beta prior (\code{"beta"}).
#' @param beta_d Parameter for the design beta prior (\code{"beta"}).
#' @param location_d Location parameter for the design point prior (\code{"Point"}).
#' @param k_d Parameter for the design default beta prior (\code{"d_beta"}).
#' @param scale_d Scale parameter for the design normal moment prior (\code{"Moment"}).
#' @param N Sample size (numeric scalar).
#' @param type_rate Character indicating which rate to control: \code{"positive"} (true/false positive rates) or \code{"negative"} (true/false negative rates).
#' @param plot_power Logical; if TRUE, plots power curves.
#' @param plot_rel Logical; if TRUE, plots the relationship between the BF and data.
#'
#'@details
#' \strong{1. Sample size determination mode (when \code{N = NULL}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size. The user must provide:
#' \itemize{
#' \item \code{type_rate} - either \code{"positive"} to control true/false positive rates, or \code{"negative"} to control true/false negative rates.
#' \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#' \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#' \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be > 1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence meets or exceeds \code{true_rate}, while the probability of misleading evidence does not exceed \code{false_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{N} is supplied):}
#'
#' If a positive numeric sample size \code{N} is provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, the arguments \code{type_rate}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Hypothesis specification:}
#'
#' The \code{alternative} argument defines the direction of the alternative hypothesis : \code{"two.sided"} for two-sided, \code{"greater"} for right-sided, or \code{"less"} for left-sided tests. The optional \code{ROPE} argument specifies bounds for an interval null hypothesis. If \code{ROPE = NULL}, a point-null test is assumed.
#'
#' \strong{Analysis Priors:}
#'
#' The analysis prior specifies the prior distribution of the correlation under the alternative hypothesis. Depending on \code{prior_analysis}, the user must supply:
#' \itemize{
#' \item \code{d_beta} (default beta): \code{k} > 0.
#' \item \code{beta} (stretched beta): \code{alpha} and \code{beta} > 0.
#' \item \code{Moment} (normal-moment prior): \code{scale} > 0.
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' A design prior can be specified for planning purposes. If provided, \code{prior_design} must be one of \code{"d_beta"}, \code{"beta"}, \code{"Moment"}, or \code{"Point"}, and the corresponding parameters must be supplied:
#' \itemize{
#' \item \code{d_beta}: \code{k_d} > 0.
#' \item \code{beta}: \code{alpha_d} and \code{beta_d} > 0.
#' \item \code{Moment}: \code{scale_d} > 0.
#' \item \code{Point}: \code{location_d} numeric scalar.
#' }
#'
#' \strong{interval null Hypothesis:}
#'
#' If \code{ROPE} is provided, the function evaluates the Bayes factor for an interval null. Otherwise, a point-null hypothesis is assumed.
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling evidence as a function of sample size. If \code{plot_rel = TRUE}, the relationship between the BF and correlation is plotted.
#'
#' @return A list of class \code{BFpower_r} containing:
#' \itemize{
#'   \item \code{type}: Test type (always "Correlation").
#'   \item \code{alternative}: the direction of the alternative hypothesis.
#'   \item \code{h0}: the value of correlation under the null hypothesis.
#'   \item \code{ROPE}: Bounds for interval null (if used).
#'   \item \code{analysis_h1}: Analysis prior parameters under the alternative hypothesis.
#'   \item \code{design_h1}: Design prior parameters under the alternative hypothesis.
#'   \item \code{results}: Data frame of probabilities of compelling and misleading evidence.
#'   \item \code{threshold}: Threshold of compelling evidence.
#'   \item \code{plot_power}: Logical, whether power curves are plotted.
#'   \item \code{plot_rel}: Logical, whether the relationship between the BF and correlation is plotted.
#' }
#'
#' @examples
#' BFpower.cor(
#'  alternative = "greater",
#'  h0 = 0,
#'    threshold = 3,
#'    true_rate = 0.8,
#'    false_rate = 0.05,
#'    prior_analysis = "d_beta",
#'    k = 1,
#'    prior_design = "Point",
#'    location_d = 0.3,
#'    plot_power = TRUE,
#'    plot_rel = TRUE
#'  )
#'
#' @export
BFpower.cor<- function(alternative , h0, ROPE = NULL,
                       threshold , true_rate, false_rate ,
                       prior_analysis , k , alpha , beta , scale ,
                       prior_design = NULL, alpha_d , beta_d , location_d ,
                       k_d , scale_d ,
                       N = NULL,  type_rate="positive",plot_power=FALSE,plot_rel=FALSE) {
  # mode
  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < -0.8 || h0 > 0.8) {
    stop("argument [h0] NULL value of rho must be a single numeric scalar between -0.8 and 0.8")
  }

  location <- h0
  dff <- dff_d <- 1
  if ( is.null(N)) mode_bf=1 else mode_bf = 0


  # sample size
  if (mode_bf == 0) {
    # Check that N is a positive numeric scalar
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 0) {
      stop("argument [N] sample size must be a positive numeric scalar when mode_bf = 0")
    }
  }else {N=3}


  # alternative
  if(alternative %in% c("two.sided", "less", "greater") == FALSE){
    stop("Argument [alternative] should be set to either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test)")
  }

  alternative <- switch(alternative,
                        "two.sided" = "!=",
                        "less"      = "<",
                        "greater"   = ">"
  )


  # Equivlance test or not
  interval <- if (is.null(ROPE)) 1 else 0

  if (!is.null(ROPE)) {

    if (alternative == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative 'two.sided', argument [ROPE] must be a numeric vector of length 2 with two distinct finite values")
      }
      # Additional bounds checks
      if (min(ROPE) < -0.5 || max(ROPE) > 0.5) {
        stop("For alternative 'two.sided', ROPE must satisfy min(ROPE) >= -0.5 and max(ROPE) <= 0.5")
      }
      if ((h0 + min(ROPE)) <= -1 || (h0 + min(ROPE)) >= 1) {
        stop("For alternative '!=', h0 + min(ROPE) must be between -1 and 1")
      }

    } else if (alternative == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', argument [ROPE] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (ROPE > 0.5) stop("For alternative 'greater', ROPE must be <= 0.5")
      if ((h0 + ROPE) >= 1) stop("For alternative 'greater', h0 + ROPE must be < 1")

    } else if (alternative == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', argument [ROPE] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (ROPE < -0.5) stop("For alternative 'less', ROPE must be >= -0.5")
      if ((h0 + ROPE) <= -1) stop("For alternative 'less', h0 + ROPE must be > -1")
    }

  }


  # analysis prior prior_analysis
  if (missing(prior_analysis)) {
    stop("argument [prior_analysis] for analysis prior must be one of `d_beta` (default stretched beta), `beta` (stretched beta), or `Moment` (normal-moment prior)")
  }

  # Analysis prior prior_analysis validation
  if (!prior_analysis %in% c("d_beta", "Moment", "beta")) {
    stop("argument [prior_analysis] for analysis prior must be one of `d_beta` (default stretched beta), `beta` (stretched beta), or `Moment` (normal-moment prior)")
  }

  # prior_analysis-specific checks
  if (prior_analysis == "d_beta") {
    alpha=beta=scale=NULL
    # 'd_beta' requires k to be a single numeric scalar > 0
    if (!exists("k") || !is.numeric(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
      stop("For prior_analysis 'd_beta', argument [k] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "beta") {
    k=scale=NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (!exists("alpha") || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For prior_analysis 'beta', argument [alpha] must be a single numeric scalar > 0")
    }
    if (!exists("beta") || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For prior_analysis 'beta', argument [beta] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "Moment") {
    k=alpha=beta=NULL
    # 'Moment' requires scale to be numeric scalar > 0
    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For prior_analysis 'Moment', argument [scale] must be a numeric scalar > 0")
    }
  }


  ##  design prior

  if (!is.null(prior_design)) {

    de_an_prior <- 0

    # Validate prior_design
    if (!prior_design %in% c("d_beta", "Moment", "beta", "Point")) {
      stop("argument [prior_design] for design prior must be one of `d_beta`, `beta`, `Moment`, or `Point`")
    }

    # prior_analysis-specific checks
    if (prior_design == "d_beta") {
      alpha_d=beta_d=scale_d=location_d=NULL

      # 'd_beta' requires k_d to be a numeric scalar > 0
      if (!exists("k_d") || !is.numeric(k_d) || length(k_d) != 1 || !is.finite(k_d) || k_d <= 0) {
        stop("For design prior 'd_beta', argument [k_d] must be a single numeric scalar > 0")
      }
    } else if (prior_design == "beta") {
      k_d=scale_d=location_d=NULL

      # 'beta' requires alpha_d and beta_d to be numeric scalars > 0
      if (!exists("alpha_d") || !is.numeric(alpha_d) || length(alpha_d) != 1 || !is.finite(alpha_d) || alpha_d <= 0) {
        stop("For design prior 'beta', argument [alpha_d] must be a single numeric scalar > 0")
      }
      if (!exists("beta_d") || !is.numeric(beta_d) || length(beta_d) != 1 || !is.finite(beta_d) || beta_d <= 0) {
        stop("For design prior 'beta', argument [beta_d] must be a single numeric scalar > 0")
      }
    } else if (prior_design == "Moment") {
      k_d=alpha_d=beta_d=NULL
      # 'Moment' requires scale_d numeric scalar > 0
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0) {
        stop("For design prior 'Moment', argument [scale_d] must be a numeric scalar > 0")
      }
    } else if (prior_design == "Point") {
      k_d=alpha_d=beta_d=scale_d=NULL

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
    if (!(type_rate %in% c("positive", "negative"))) {
      stop("argument [positive] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (type_rate,
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
    if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold <= 1) {
      stop("argument [threshold] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=FP=0
  }


  tryCatch(
    suppressWarnings({
      if ( interval == 1) {
        results=r_table(threshold, target, prior_analysis, k, alpha, beta, h0, location, scale, dff,
                        alternative, prior_design, location_d, k_d, alpha_d, beta_d, scale_d,
                        dff_d, de_an_prior, N, mode_bf, FP, direct)
      } else {
        results=re_table(threshold, target, prior_analysis, k, alpha, beta, h0, location, scale, dff,
                         alternative, prior_design, location_d, k_d, alpha_d, beta_d, scale_d,
                         dff_d, de_an_prior, N, mode_bf, FP, ROPE, direct)
      }
    }),
    error = function(err) {
      message("Required sample size > 5,000")
      stop(NaN)
    }
  )
  type = "correlation"
  analysis_h1 <- list(
    prior = prior_analysis,
    k = k,
    alpha=alpha,
    beta=beta,
    scale=scale
  )

  if (!is.null(prior_design)) {

    # Base fields always included
    design_h1 <-  list(
      prior = prior_design,
      location=location_d,
      k = k_d,
      alpha=alpha_d,
      beta=beta_d,
      scale=scale_d
    )


  } else {

    # prior_design is NULL > fill all fields with NULL
    design_h1 <- list(
      prior = NULL,
      location=NULL,
      k = NULL,
      alpha=NULL,
      beta=NULL,
      scale=NULL)

  }


  object <- list(
    type = "Correlation",
    alternative = alternative,
    h0=h0,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    threshold = threshold,
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
#' This function performs sample size determination (when \code{N = NULL}) or
#' calculates the probability of compelling/misleading evidence for a fixed sample
#' size.
#'
#' @param threshold Numeric scalar. Threshold for compelling evidence (must be > 1).
#'
#' @param true_rate Targeted true positive or true negative rate (used only when
#'   sample size determination is requested; \code{N = NULL}).
#'
#' @param false_rate Targeted false positive or false negative rate (used only when
#'   sample size determination is requested; \code{N = NULL}).
#'
#' @param p Number of predictors in the reduced model.
#'
#' @param k Number of predictors in the full model (must satisfy \code{k > p}).
#'
#' @param prior_analysis Analysis prior model under the alternative hypothesis:
#'   \code{"effectsize"} or \code{"Moment"}.
#'
#' @param dff Degrees of freedom for the analysis prior under the alternative
#'   hypothesis. Must be a positive scalar, and must be at least 3 if
#'   \code{prior_analysis = "Moment"}.
#'
#' @param rscale Scale parameter for the analysis effect-size prior (only used when
#'   \code{prior_analysis = "effectsize"}).
#'
#' @param f_m Cohen's \eqn{f} effect-size parameter for the analysis prior (must be > 0).
#'
#' @param prior_design Design prior model under the alternative hypothesis:
#'   \code{"effectsize"}, \code{"Moment"}, or \code{"Point"}.
#'
#' @param dff_d Degrees of freedom for the design prior. Must be a positive scalar,
#'   and at least 3 if \code{prior_design = "Moment"}.
#'
#' @param rscale_d Scale parameter for the design effect-size prior
#'   (only used when \code{prior_design = "effectsize"}).
#'
#' @param f_m_d Cohen's \eqn{f} value for the design prior or the effect-size of the
#'   point design prior.
#'
#' @param N Sample size. If \code{NULL}, sample size determination is performed.
#'
#' @param type_rate Either `"positive"` (control true/false positive rates) or
#'   `"negative"` (control true/false negative rates).
#'
#' @param ROPE Numeric bounds for the interval null (only used when interval
#'   Bayes factors are required).
#'
#' @param plot_power Logical. Whether to plot power curves when
#'   sample size determination is requested.
#'
#' @param plot_rel Logical. Whether to plot the relationship between the BF and data.
#'
#' @details
#'
#' \strong{1. Sample size determination mode (when \code{N = NULL}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size to achieve the desired configuration below. The user must provide:
#' \itemize{
#'   \item \code{type_rate} - either \code{"positive"} to control true/false positive rates, or \code{"negative"} to control true/false negative rates.
#'   \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#'   \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#'   \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be > 1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence meets or exceeds \code{true_rate}, while the probability of misleading evidence does not exceed \code{false_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{N} is supplied):}
#'
#' If a positive numeric sample size \code{N} is provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, the arguments \code{type_rate}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Model specification:}
#'
#' The function requires the user to specify the full model (\code{k} predictors) and the reduced model (\code{p} predictors, \code{k > p}), and the analysis prior under the alternative hypothesis. Depending on the chosen \code{prior_analysis}, different arguments are required:
#' \itemize{
#'   \item \code{prior_analysis = "effectsize"}: requires \code{rscale} (scale parameter) and \code{f_m} (Cohen's f effect-size), and \code{dff} (degrees of freedom).
#'   \item \code{prior_analysis = "Moment"}: requires \code{f_m} (Cohen's f effect-size) and \code{dff} (degrees of freedom, must be >= 3); \code{rscale} is not used.
#' }
#' The design prior under the alternative hypothesis can optionally be specified using \code{prior_design}, which can be:
#' \itemize{
#'   \item \code{"effectsize"}: requires \code{rscale_d}, \code{f_m_d}, and \code{dff_d}.
#'   \item \code{"Moment"}: requires \code{f_m_d} and \code{dff_d} (>=3); \code{rscale_d} is not used.
#'   \item \code{"Point"}: requires \code{f_m_d} only; \code{rscale_d} and \code{dff_d} are not used.
#' }
#'
#' \strong{interval null Hypothesis:}
#'
#' If \code{ROPE} is provided, the function evaluates the Bayes factor for an interval null. Otherwise, a point-null hypothesis is assumed.
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling evidence as a function of sample size. If \code{plot_rel = TRUE}, the relationship between the Bayes factor and Cohen's \code{f} is plotted.
#'
#' @return A list of class \code{BFpower_f} containing:
#' \itemize{
#'   \item \code{type}: Test type ("Regression/ANOVA").
#'   \item \code{k}, \code{p}: Model sizes.
#'   \item \code{ROPE}: Bounds for interval null (if used).
#'   \item \code{analysis_h1}: List describing the analysis prior.
#'   \item \code{design_h1}: List describing the design prior.
#'   \item \code{results}: Data frame of the probabilities of compelling/misleading evidence and the required or supplied sample size.
#'   \item \code{threshold}: Threshold of compelling evidence.
#'   \item \code{plot_power}: Logical, whether power curves are plotted.
#'   \item \code{plot_rel}: Logical, whether the relationship between the BF and data is plotted.
#' }
#'
#' If sample size determination fails, the function returns \code{NaN} and prints a message.
#'
#' @examples
#'BFpower.f.test(
#'  threshold = 3,
#'  true_rate = 0.8,
#'  false_rate = 0.05,
#'  p = 3,
#'  k = 4,
#'  prior_analysis = "effectsize",
#'  dff = 3,
#'  rscale = 0.18,
#'  f_m = 0.1,
#'  prior_design = "Point",
#'  f_m_d = 0.1,
#'  plot_power = TRUE,
#'  plot_rel = TRUE)
#'
#' @export
BFpower.f.test <- function(threshold, true_rate, false_rate , p , k ,
                           prior_analysis , dff , rscale , f_m ,
                           prior_design = NULL, dff_d, rscale_d, f_m_d ,
                           N = NULL, type_rate="positive", ROPE = NULL,plot_power=FALSE,plot_rel=FALSE) {

  ## mode
  if ( is.null(N)) mode_bf=1 else mode_bf = 0

  ## prior_analysis parameterchecks
  # Check p
  if (is.null(p) || !is.numeric(p) || length(p) != 1 || is.na(p)) {
    stop("argument [p] (number of predictors in the reduced prior_analysis) must be a positive numeric scalar")
  }

  # Check k
  if (is.null(k) || !is.numeric(k) || length(k) != 1 || is.na(k)) {
    stop("argument [k] (number of predictors in the full prior_analysis) must be a positive numeric scalar")
  }

  # Check relation
  if (k <= p) {
    stop("argument [k] (predictors in full prior_analysis) must be greater than argument [p] (predictors in reduced prior_analysis)")
  }

  # Equivlance test or not
  interval <- if (is.null(ROPE)) 1 else 0

  # analysis prior prior_analysis
  if (missing(prior_analysis)) {
    stop("argument [prior_analysis] for analysis prior should be set to either `effectsize`, or `Moment`")
  }
  if(prior_analysis %in% c("effectsize","Moment") == FALSE){
    stop("argument [prior_analysis] for analysis prior should be set to either `effectsize`, or `Moment`")
  }

  if (prior_analysis =="effectsize"){
    if (!is.numeric(rscale) || length(rscale) != 1 || !is.finite(rscale) || rscale <= 0) {
      stop("argument [rscale] scale parameter must be a positive numeric scalar")
    }
  }

  if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
    stop("argument [dff] degrees of freedom  for analysis prior must be a positive numeric scalar when prior_analysis='t-distribution'")
  }

  if (!is.numeric(f_m) || length(f_m) != 1 || !is.finite(f_m) || f_m <= 0) {
    stop("argument [f_m] Cohen's f  for analysis prior must be a positive numeric scalar")
  }

  if (prior_analysis == "Moment"){
    rscale=NULL
    if (dff < 3) {
      stop("argument [dff] degrees of freedom for Moment analysis prior must be at least 3")
    }
  }

  # design prior

  if (!is.null(prior_design)) {

    de_an_prior <- 0

    # Validate prior_design
    if (!(prior_design %in% c("effectsize","Moment","Point"))) {
      stop("argument [prior_design] for design prior must be either `effectsize`, or `Moment`")
    }



    if (prior_design =="effectsize"){
      if (!is.numeric(rscale_d) || length(rscale_d) != 1 || !is.finite(rscale_d) || rscale_d <= 0) {
        stop("argument [rscale] scale parameter must be a positive numeric scalar")
      }

      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0) {
        stop("argument [dff] degrees of freedom  for design prior must be a positive numeric scalar when prior_analysis='t-distribution'")
      }
      if (!is.numeric(f_m_d) || length(f_m_d) != 1 || !is.finite(f_m_d) || f_m_d <= 0) {
        stop("argument [f_m] Cohen's f  for design prior must be a positive numeric scalar")
      }


    }



    if (prior_design == "Moment"){
      rscale_d=NULL
      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0) {
        stop("argument [dff] degrees of freedom  for design prior must be a positive numeric scalar when prior_analysis='t-distribution'")
      }
      if (!is.numeric(f_m_d) || length(f_m_d) != 1 || !is.finite(f_m_d) || f_m_d <= 0) {
        stop("argument [f_m] Cohen's f  for design prior must be a positive numeric scalar")
      }

      if (dff_d < 3) {
        stop("argument [dff_d] degrees of freedom for moment design prior must be at least 3")
      }
    }

    if (prior_design == "Point"){
      rscale_d=dff_d=NULL

      if (!is.numeric(f_m_d) || length(f_m_d) != 1 || !is.finite(f_m_d) || f_m_d <= 0) {
        stop("argument [f_m] Cohen's f  for design prior must be a positive numeric scalar")
      }

    }

  } else {
    de_an_prior <- 1
  }


  # desired power and strength of evidence
  if (mode_bf==1){
    if (!(type_rate %in% c("positive", "negative"))) {
      stop("argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (type_rate,
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
    if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold <= 1) {
      stop("argument [threshold] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=FP=0
  }









  results = tryCatch({
    suppressWarnings({
      if (!is.null(interval) && interval == 1) {
        f_table(threshold, target, p, k, dff, rscale, f_m, prior_analysis,
                dff_d, rscale_d, f_m_d, prior_design, de_an_prior, N,
                mode_bf, FP, direct)
      } else {
        fe_table(threshold, target, p, k, dff, rscale, f_m, prior_analysis,
                 dff_d, rscale_d, f_m_d, prior_design, de_an_prior, N,
                 mode_bf, ROPE, FP, direct)
      }
    })
  }, error = function(err) {

    if(dff<3|dff_d<3){
      stop(" Degrees of freedom[dff] for analysis prior or [dff_d] for design prior should be at least 3")

    }

    message("Required sample size > 10,000")

    return(NaN)
  })


  type = "Regression/ANOVA"
  analysis_h1 <- list(
    prior = prior_analysis,
    rscale = rscale,
    f_m = f_m,
    dff=dff
  )

  if (!is.null(prior_design)) {

    # Base fields always included
    design_h1 <- list(
      prior = prior_design,
      rscale = rscale_d,
      f_m = f_m_d,
      dff=dff_d)

  } else {

    # prior_design is NULL > fill all fields with NULL
    design_h1 <- list(
      prior_analysis  = NULL,
      rscale = NULL,
      f_m    = NULL,
      dff    = NULL)
  }


  object <- list(
    type = type,
    k=k,
    p=p,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    threshold = threshold,
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
#' @param alternative The direction of the alternative hypothesis : two-sided (\code{"two.sided"}), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param threshold Numeric scalar. Threshold for compelling evidence (must be > 1).
#' @param h0 Null proportion value for the test (numeric scalar between 0.1 and 0.9).
#' @param true_rate Targeted true positive rate  or true negative rate .
#' @param false_rate Targeted false positive rate  or false negative rate .
#' @param prior_analysis Analysis prior under the alternative hypothesis: \code{"beta"} or \code{"Moment"}.
#' @param alpha Parameter for the analysis beta prior (used when \code{prior_analysis = "beta"}).
#' @param beta Parameter for the analysis beta prior (used when \code{prior_analysis = "beta"}).
#' @param scale Scale parameter for the analysis moment prior (used when \code{prior_analysis = "Moment"}).
#' @param prior_design Design prior under the alternative hypothesis: \code{"beta"}, \code{"Moment"}, or \code{"Point"}.
#' @param alpha_d Parameter for the design beta prior (used when \code{prior_design = "beta"}).
#' @param beta_d Parameter for the design beta prior (used when \code{prior_design = "beta"}).
#' @param location_d Proportion value for the design point prior (\code{prior_design = "Point"}). Represents the true proportion under the alternative hypothesis.
#' @param scale_d Scale parameter for the design moment prior (used when \code{prior_design = "Moment"}).
#' @param N Sample size. If \code{NULL}, sample size determination is performed.
#' @param ROPE Numeric bounds for the interval null (used when computing interval Bayes factors).
#' @param type_rate Either `"positive"` (controls true/false positive rates) or `"negative"` (controls true/false negative rates).
#' @param plot_power Logical. Whether to plot power curves when sample size determination is requested.
#' @param plot_rel Logical. Whether to plot probability of misleading evidence.
#'
#' @details
#'
#' \strong{1. Sample size determination mode (when \code{N = NULL}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size needed to achieve the desired configuration below. The user must provide:
#' \itemize{
#' \item \code{type_rate} - either \code{"positive"} to control true/false positive rates or \code{"negative"} to control true/false negative rates.
#' \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#' \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#' \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be > 1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence meets or exceeds \code{true_rate}, while the probability of misleading evidence does not exceed \code{false_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{N} is supplied):}
#'
#' If a positive numeric sample size \code{N} is provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, \code{type_rate}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Model specification:}
#'
#' The user must specify the analysis prior under the alternative hypothesis using \code{prior_analysis}:
#' \itemize{
#' \item \code{prior_analysis = "beta"}: requires \code{alpha} and \code{beta} parameters (shape parameters of the beta distribution).
#' \item \code{prior_analysis = "Moment"}: requires \code{scale} parameter (scale of the moment prior).
#' }
#' The design prior under the alternative hypothesis can optionally be specified using \code{prior_design}:
#' \itemize{
#' \item \code{"beta"}: requires \code{alpha_d} and \code{beta_d}.
#' \item \code{"Moment"}: requires \code{scale_d}.
#' \item \code{"Point"}: requires \code{location_d}, representing the true proportion under the alternative hypothesis.
#' }
#' If \code{prior_design} is \code{NULL}, no design prior is used.
#'
#' \strong{interval null Hypothesis:}
#'
#' If \code{ROPE} is provided, the function evaluates the Bayes factor for an interval null. Otherwise, a point-null hypothesis is assumed.
#'
#' \strong{Hypothesis:}
#'
#' The function supports one-sided (\code{"greater"} or \code{"less"}) and two-sided (\code{"two.sided"}) tests. Design prior and interval null bounds must be consistent with the directionality of the hypothesis.
#'
#' \strong{Plotting:}
#'
#' If \code{plot_power = TRUE}, the function plots the probability of compelling evidence as a function of sample size. If \code{plot_rel = TRUE}, the relationship between the Bayes factor and the number of successes (proportion) is plotted.
#'
#'
#' @return A list of class \code{"BFpower_bin"} containing:
#' \itemize{
#'   \item \code{type}: Test type ("One proportion").
#'   \item \code{alternative}: alternative hypothesis.
#'   \item \code{h0}: The proportion under the null hypothesis.
#'   \item \code{analysis_h1}: List describing the analysis prior.
#'   \item \code{design_h1}: List describing the design prior.
#'   \item \code{results}: Data frame of probabilities of compelling/misleading evidence and the required or supplied sample size.
#'   \item \code{threshold}: Compelling-evidence threshold.
#'   \item \code{plot_power}: Logical, whether power curves are plotted.
#'   \item \code{plot_rel}: Logical, whether the relationship between the BF and data is plotted.
#' }
#'
#' If sample size determination fails, the function returns \code{NaN} and prints a message.
#'
#' @examples
#' BFpower.bin(
#'   alternative = "greater",
#'   threshold = 3,
#'   true_rate = 0.8,
#'   false_rate = 0.05,
#'   h0 = 0.5,
#'   prior_analysis = "beta",
#'   alpha = 1,
#'   beta = 1,
#'   plot_rel = TRUE,
#'   plot_power = TRUE)
#'
#'
#' @export
BFpower.bin <- function(alternative ,threshold , h0 ,
                        true_rate , false_rate ,
                        prior_analysis , alpha , beta , scale ,
                        prior_design = NULL, alpha_d , beta_d , location_d , scale_d ,
                        N = NULL, ROPE = NULL, type_rate="positive",plot_power=FALSE,plot_rel=FALSE) {

  # mode
  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < .1 || h0 > 0.9) {
    stop("argument [h0] NULL value of proportion must be a single numeric scalar between .1 and 0.9")
  }

  location <- h0
  if ( is.null(N)) mode_bf=1 else mode_bf = 0


  # sample size
  if (mode_bf == 0) {
    # Check that N is a positive numeric scalar
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 0) {
      stop("argument [N] sample size must be a positive numeric scalar when mode_bf = 0")
    }
  }else {N=3}


  # alternative
  if(alternative %in% c("two.sided", "less", "greater") == FALSE){
    stop("Argument [alternative] should be set to either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test)")
  }

  alternative <- switch(alternative,
                        "two.sided" = "!=",
                        "less"      = "<",
                        "greater"   = ">"
  )


  # Equivlance test or not
  interval <- if (is.null(ROPE)) 1 else 0

  if (!is.null(ROPE)) {

    if (alternative == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative '!=', argument [ROPE] must be a numeric vector of length 2 with two distinct finite values")
      }
      # Additional bounds checks
      if (min(ROPE) < -0.5 || max(ROPE) > 0.5) {
        stop("For alternative '!=', ROPE must satisfy min(ROPE) >= -0.5 and max(ROPE) <= 0.5")
      }
      if ((h0 + min(ROPE)) <= 0 || (h0 + min(ROPE)) >= 1) {
        stop("For alternative '!=', h0 + min(ROPE) must be between 0 and 1")
      }

    } else if (alternative == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative '>', argument [ROPE] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (ROPE > 0.5) stop("For alternative '>', ROPE must be <= 0.5")
      if ((h0 + ROPE) >= 1) stop("For alternative '>', h0 + ROPE must be < 1")

    } else if (alternative == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative '<', argument [ROPE] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (ROPE < -0.5) stop("For alternative '<', ROPE must be >= -0.5")
      if ((h0 + ROPE) <= -1) stop("For alternative '<', h0 + ROPE must be > 0")
    }

  }


  # analysis prior prior_analysis
  if (missing(prior_analysis)) {
    stop("argument [prior_analysis] for analysis prior must be one of `beta`, or `Moment` (normal-moment prior)")
  }

  # Analysis prior prior_analysis validation
  if (!prior_analysis %in% c("Moment", "beta")) {
    stop("argument [prior_analysis] for analysis prior must be one of `beta` , or `Moment` (normal-moment prior)")
  }

  # prior_analysis-specific checks
  if (prior_analysis == "beta") {
    scale=NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (!exists("alpha") || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For prior_analysis 'beta', argument [alpha] must be a single numeric scalar > 0")
    }
    if (!exists("beta") || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For prior_analysis 'beta', argument [beta] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "Moment") {
    alpha=beta=NULL
    # 'Moment' requires scale to be numeric scalar > 0
    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For prior_analysis 'Moment', argument [scale] must be a numeric scalar > 0")
    }
  }


  ##  design prior

  if (!is.null(prior_design)) {

    de_an_prior <- 0

    # Validate prior_design
    if (!prior_design %in% c("Moment", "beta", "Point")) {
      stop("argument [prior_design] for design prior must be one of `beta`, `Moment`, or `Point`")
    }

    # prior_analysis-specific checks
    if (prior_design == "beta") {
      scale_d=location_d=NULL

      # 'beta' requires alpha_d and beta_d to be numeric scalars > 0
      if (!exists("alpha_d") || !is.numeric(alpha_d) || length(alpha_d) != 1 || !is.finite(alpha_d) || alpha_d <= 0) {
        stop("For design prior 'beta', argument [alpha_d] must be a single numeric scalar > 0")
      }
      if (!exists("beta_d") || !is.numeric(beta_d) || length(beta_d) != 1 || !is.finite(beta_d) || beta_d <= 0) {
        stop("For design prior 'beta', argument [beta_d] must be a single numeric scalar > 0")
      }
    } else if (prior_design == "Moment") {
      alpha_d=beta_d=NULL
      # 'Moment' requires scale_d numeric scalar > 0
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0) {
        stop("For design prior 'Moment', argument [scale_d] must be a numeric scalar > 0")
      }
    } else if (prior_design == "Point") { alpha_d <- beta_d <- scale_d <- NULL  # Not needed for 'Point' prior

    # 'Point' prior requires location_d, which represents the true proportion
    # under the alternative alternative. It must be a numeric scalar.
    if (!is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d)) {
      stop("For design prior 'Point', argument [location_d] true proportion must be a numeric scalar")
    }

    # Validate location_d against alternative and h0
    # - "!=" : location_d must not equal h0
    # - ">"  : location_d must be greater than h0
    # - "<"  : location_d must be less than h0
    if (alternative == "!=" && location_d == h0) {
      stop("For alternative '!=', argument [location_d] true proportion must not equal h0")
    } else if (alternative == ">" && location_d <= h0) {
      stop("For alternative '>', argument [location_d] true proportion must be greater than h0")
    } else if (alternative == "<" && location_d >= h0) {
      stop("For alternative '<', argument [location_d] true proportion must be less than h0")
    }
    }

  } else {
    de_an_prior <- 1
  }


  # desired power and strength of evidence
  if (mode_bf==1){
    if (!(type_rate %in% c("positive", "negative"))) {
      stop("argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (type_rate,
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
    if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold <= 1) {
      stop("argument [threshold] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=FP=0
  }



  results <-tryCatch({
    suppressWarnings({
      if (is.null(ROPE)) {
        bin_table(threshold, target, h0, alpha, beta, location, scale, prior_analysis, alternative,
                  alpha_d, beta_d, location_d, scale_d, prior_design, de_an_prior, N,
                  mode_bf, FP, direct)
      } else {
        bin_e_table(threshold, target, h0, alpha, beta, location, scale, prior_analysis, alternative,
                    alpha_d, beta_d, location_d, scale_d, prior_design, de_an_prior, N,
                    mode_bf, FP, ROPE, direct)
      }
    })
  }, error = function(err) {
    message("Sample size cannot be determined")
    return(NaN)
  })

  type = "One proportion"
  analysis_h1 <- list(
    prior = prior_analysis,
    alpha=alpha,
    beta=beta,
    scale=scale
  )

  if (!is.null(prior_design)) {

    # Base fields always included
    design_h1 <-  list(
      prior = prior_design,
      alpha=alpha_d,
      beta=beta_d,
      scale=scale_d
    )


  } else {

    # prior_design is NULL > fill all fields with NULL
    design_h1 <- list(
      prior = NULL,
      alpha=NULL,
      beta=NULL,
      scale=NULL)

  }


  object <- list(
    type = type,
    alternative = alternative,
    h0=h0,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    threshold = threshold,
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
#' @param threshold Threshold of compelling evidence.
#' @param true_rate Targeted true positive rate (if \code{positive = "positive"}) or true negative rate (if \code{positive = "negative"}).
#' @param a0 Alpha parameter of the Beta prior under the null hypothesis.
#' @param b0 Beta parameter of the Beta prior under the null hypothesis.
#' @param a1 Alpha parameter of the Beta analysis prior for group 1 under the alternative hypothesis.
#' @param b1 Beta parameter of the Beta analysis prior for group 1 under the alternative hypothesis.
#' @param a2 Alpha parameter of the Beta analysis prior for group 2 under the alternative hypothesis.
#' @param b2 Beta parameter of the Beta analysis prior for group 2 under the alternative hypothesis.
#' @param prior_design_1 the design prior of group 1: \code{"beta"}, \code{"Point"}, or \code{"same"} (if \code{"same"}, the design prior is identical to the analysis prior).
#' @param a1d Alpha parameter of the design prior for group 1 (used if \code{model1 = "beta"}).
#' @param b1d Beta parameter of the design prior for group 1 (used if \code{model1 = "beta"}).
#' @param dp1 True proportion for group 1 in the design prior (used if \code{model1 = "Point"}).
#' @param prior_design_2 the design prior of group 2: \code{"beta"}, \code{"Point"}, or \code{"same"} (if \code{"same"}, the design prior is identical to the analysis prior).
#' @param a2d Alpha parameter of the design prior for group 2 (used if \code{model2 = "beta"}).
#' @param b2d Beta parameter of the design prior for group 2 (used if \code{model2 = "beta"}).
#' @param dp2 True proportion for group 2 in the design prior (used if \code{model2 = "Point"}).
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param type_rate Choose \code{"positive"} to control true/false positive rates or \code{"negative"} to control true/false negative rates.
#' @param plot_power Logical; if TRUE, plot the power curve.
#' @param plot_rel Logical; if TRUE, plot the grid for the values of BF across all possible combination of x1 and x2.
#'
#' @details
#'
#' \strong{1. Sample size determination mode (when \code{n1 = NULL} and \code{n2 = NULL}):}
#'
#' If no sample sizes are provided for the two groups, the function calculates the minimum sample sizes needed to achieve the desired configuration. The user must provide:
#' \itemize{
#' \item \code{type_rate} - either \code{"positive"} to control true/false positive rates or \code{"negative"} to control true/false negative rates.
#' \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#' \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be > 1).
#' }
#'
#' The function iteratively finds the smallest sample sizes for which the probability of obtaining compelling evidence meets or exceeds \code{true_rate}.
#'
#' \strong{2. Fixed-sample analysis mode (when \code{n1} and \code{n2} are supplied):}
#'
#' If positive numeric sample sizes \code{n1} and \code{n2} are provided, the function computes the probabilities of obtaining compelling or misleading evidence for these fixed sample sizes. In this mode, \code{type_rate} and \code{true_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
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
#' \item \code{prior_design_1}, \code{a1d}, \code{b1d}, \code{dp1} - design prior for group 1 (\code{"same"} uses the analysis prior, \code{"beta"} requires Beta parameters, \code{"Point"} uses a fixed proportion).
#' \item \code{prior_design_2}, \code{a2d}, \code{b2d}, \code{dp2} - design prior for group 2.
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
#' threshold = 3,
#' true_rate = 0.8,
#' a0 = 1,
#' b0 = 1,
#' a1 = 156,
#' b1 = 339,
#' a2 = 151,
#' b2 = 339,
#' plot_power = TRUE,
#' plot_rel = TRUE)
#'
#' @export

BFpower.props <- function(threshold , true_rate , a0 , b0 , a1 , b1 ,
                          a2 , b2 , prior_design_1 = "same",
                          a1d , b1d , dp1 , prior_design_2 = "same",
                          a2d, b2d , dp2 ,
                          n1 = NULL, n2 = NULL,type_rate="positive",plot_power=FALSE,plot_rel=FALSE) {

  # Check NULL
  if (is.null(n1) && is.null(n2)) {
    mode_bf <- 1
  } else {
    mode_bf <- 0
  }

  # If both n1 and n2 are NULL > mode_bf = 1
  # If mode_bf = 1, check target range
  if (mode_bf==1){
    if (!(type_rate %in% c("positive", "negative"))) {
      stop("argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rate)")
    }
    direct= switch (type_rate,
                    "positive" = "h1",
                    "negative" = "h0"
    )
    if (!is.numeric(true_rate) || length(true_rate) != 1 || !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999) {
      stop("argument [true_rate] targeted true positive or negative rate must be a numeric scalar strictly greater than 0.6 and less than 0.999")
    }
    target = true_rate

    if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold <= 1) {
      stop("argument [threshold] threshold of compelling evidence must be a numeric scalar greater than 1")
    }
  } else{
    target=0
  }



  # If not NULL, check numeric, scalar, integer
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

  # NULL hypothesis
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


  # --- Check prior_design_1 assumptions for design prior on theta1 ---

  if (prior_design_1 == "same") {

    # Automatically set all to NULL
    a1d <- 1
    b1d <- 1
    dp1 <- 0.5

  } else if (prior_design_1 == "beta") {

    # a1d and b1d must be valid Beta parameters
    if (!is.numeric(a1d) || length(a1d) != 1 || a1d <= 0) {
      stop("arg [a1d] alpha for the Beta design prior on \u03b8\u2081 must be a positive numeric scalar (> 0).")
    }
    if (!is.numeric(b1d) || length(b1d) != 1 || b1d <= 0) {
      stop("arg [b1d] beta for the Beta design prior on \u03b8\u2081 must be a positive numeric scalar (> 0).")
    }

    # dp1 irrelevant for beta prior > set to NULL automatically
    dp1 <- 0.5

  } else if (prior_design_1 == "Point") {

    # Automatically set Beta parameters to NULL
    a1d <- 1
    b1d <- 1

    # dp1 must be numeric between 0 and 1
    if (!is.numeric(dp1) || length(dp1) != 1) {
      stop("arg [dp1] true \u03b8\u2081 must be a numeric scalar for prior_design_1 = 'Point'.")
    }
    if (dp1 <= 0 || dp1 >= 1) {
      stop("arg [dp1] must be > 0 and < 1 for prior_design_1 = 'Point'.")
    }

  } else {
    stop("arg [prior_design_1] must be one of: 'same', 'beta', 'Point'.")
  }

  # --- Check prior_design_2 assumptions for design prior on theta2 ---

  if (prior_design_2 == "same") {

    # Automatically set all to NULL
    a2d <- 1
    b2d <- 1
    dp2 <- .5

  } else if (prior_design_2 == "beta") {

    # a2d and b2d must be valid Beta parameters
    if (!is.numeric(a2d) || length(a2d) != 1 || a2d <= 0) {
      stop("arg [a2d] alpha for the Beta design prior on theta2 must be a positive numeric scalar (> 0).")
    }
    if (!is.numeric(b2d) || length(b2d) != 1 || b2d <= 0) {
      stop("arg [b2d] beta for the Beta design prior on theta2 must be a positive numeric scalar (> 0).")
    }

    # dp2 irrelevant for beta prior > set to NULL automatically
    dp2 <- .5

  } else if (prior_design_2 == "Point") {

    # Automatically set Beta parameters to NULL
    a2d <- 1
    b2d <- 1

    # dp2 must be numeric between 0 and 1
    if (!is.numeric(dp2) || length(dp2) != 1) {
      stop("arg [dp2] must be a numeric scalar for prior_design_2 = 'Point'.")
    }
    if (dp2 <= 0 || dp2 >= 1) {
      stop("arg [dp2] must be > 0 and < 1 for prior_design_2 = 'Point'.")
    }

  } else {
    stop("arg [prior_design_2] must be one of: 'same', 'beta', 'Point'.")
  }



  r <- 1
  results=tryCatch({
    suppressWarnings({
      pro_table_p2(threshold, target, a0, b0, a1, b1, a2, b2, r, prior_design_1,
                   a1d, b1d, dp1, prior_design_2, a2d, b2d, dp2, mode_bf, n1, n2, direct)
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
    prior=prior_design_1,
    a = a1d,
    b = b1d,
    p = dp1
  )
  design_h1_theta_2 <- list(
    prior=prior_design_2,
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
    threshold = threshold,
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
#' against either a point null hypothesis or an interval null hypothesis.
#'
#' @param tval Numeric scalar. Observed t-value from the one-sample t-test.
#' @param df Numeric scalar. Degrees of freedom of the t-test (must be >= 1).
#' @param prior_analysis Character string. Statistical model for the analysis prior under the alternative hypothesis.
#'   Choices are:
#'   \describe{
#'     \item{"Normal"}{Normal distribution.}
#'     \item{"Moment"}{Normal moment prior.}
#'     \item{"t-distribution"}{Scaled t-distribution.}
#'   }
#' @param location Numeric scalar. Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Numeric scalar. Scale parameter for the analysis prior under the alternative hypothesis (must be > 0).
#' @param dff Numeric scalar. Degrees of freedom for the t-distribution prior (only required if \code{prior_analysis = "t-distribution"}; must be > 0). Ignored otherwise.
#' @param alternative Character string. The direction of the alternative hypothesis. One of:
#'   \describe{
#'     \item{"!="}{Two-sided (difference from 0).}
#'     \item{">"}{Right-sided (greater than 0).}
#'     \item{"<"}{Left-sided (less than 0).}
#'   }
#' @param ROPE Optional numeric vector. Specifies bounds for an interval null hypothesis. For:
#'   \describe{
#'     \item{Two-sided (\code{"two.sided"})}{Must be a numeric vector of length 2 with two distinct finite values.}
#'     \item{Right-sided (\code{"greater"})}{Must be a numeric scalar > 0.}
#'     \item{Left-sided (\code{"less"})}{Must be a numeric scalar < 0.}
#'   }
#'
#' @return A list of class \code{BFvalue_t} containing:
#'   \describe{
#'     \item{type}{Character, indicating "One-sample t-test".}
#'     \item{bf10}{Numeric, the Bayes factor (BF10).}
#'     \item{tval}{Observed t-value.}
#'     \item{df}{Degrees of freedom.}
#'     \item{analysis_h1}{List with the analysis prior parameters under H1: \code{prior_analysis}, \code{location}, \code{scale}, and optionally \code{dff}.}
#'     \item{alternative}{Character, the direction of the alternative hypothesis.}
#'     \item{ROPE}{Optional numeric vector of interval null bounds.}
#'   }
#'
#' @examples
#' BF10.ttest.OneSample(
#' tval = 2,
#' df = 50,
#' prior_analysis = "t-distribution",
#' location = 0,
#' scale = 0.707,
#' dff = 1,
#' alternative = "two.sided")
#'
#'
#' @export
BF10.ttest.OneSample <- function(tval, df, prior_analysis, location, scale, dff, alternative, ROPE = NULL) {

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
  # alternative
  if(alternative %in% c("two.sided", "less", "greater") == FALSE){
    stop("Argument [alternative] should be set to either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test)")
  }

  alternative <- switch(alternative,
                        "two.sided" = "!=",
                        "less"      = "<",
                        "greater"   = ">"
  )
  # Check e if provided
  if (!is.null(ROPE)) {
    if (alternative == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative 'two.sided', argument [e] must be a numeric vector of length 2 with two distinct finite values")
      }
    }
    if (alternative == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'less', argument [e] must be a numeric scalar > 0")
      }
    }
    if (alternative == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'greater', argument [e] must be a numeric scalar < 0")
      }
    }
  }

  # Validate analysis prior
  if (missing(prior_analysis)) {
    stop("argument [prior_analysis] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }

  if (!prior_analysis %in% c("Normal","Moment","t-distribution")) {
    stop("argument [prior_analysis] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("argument [location] for analysis prior must be a numeric scalar")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
    stop("argument [scale] for analysis prior must be a positive numeric scalar (i.e., scale > 0)")
  }
  if (prior_analysis == "t-distribution") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("argument [dff] degrees of freedom for analysis prior must be a positive numeric scalar when prior_analysis='t-distribution'")
    }
  } else {
    dff <- 0
  }
  ## -----------------------------
  ## Call appropriate function
  ## -----------------------------
  suppressWarnings(
    if (is.null(ROPE)) {
      bf10=t1_BF10(tval, df, prior_analysis, location, scale, dff, alternative)
    } else {
      bf10=t1e_BF10(tval, df, prior_analysis, location, scale, dff, alternative, ROPE)
    }
  )

  type = "One-sample t-test"
  analysis_h1 <- list(
    prior = prior_analysis,
    location = location,
    scale = scale
  )
  if (prior_analysis == "t-distribution") {
    analysis_h1$dff <- dff
  }
  object=list(type=type,bf10=bf10,tval=tval,df=df,analysis_h1=analysis_h1,alternative=alternative,ROPE=ROPE)

  class(object) <- "BFvalue_t"

  return(object)
}



#' Bayes Factor for Two-Sample Bayesian t-Test
#'
#' Compute the Bayes factor (BF10) for a two-sample independent-samples t-test. Supports both point-null and interval-null hypotheses.
#'
#' @param tval Numeric scalar. Observed t-value from the two-sample t-test.
#' @param N1 Numeric scalar. Sample size of group 1 (must be > 2, will be rounded to nearest integer).
#' @param N2 Numeric scalar. Sample size of group 2 (must be > 2, will be rounded to nearest integer).
#' @param prior_analysis Character. Analysis prior under the alternative hypothesis:
#'   \code{"Normal"}, \code{"Moment"} (normal-moment prior), or \code{"t-distribution"}.
#' @param location Numeric scalar. Location parameter of the analysis prior.
#' @param scale Numeric scalar > 0. Scale parameter of the analysis prior.
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior (required if prior_analysis = \code{"t-distribution"}; ignored otherwise).
#' @param alternative Character. The direction of the alternative hypothesis two-sided (\code{"two.sided"}), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric. Bounds for an interval null:
#'   - For \code{"two.sided"}, must be a numeric vector of length 2 with distinct finite values.
#'   - For \code{"greater"}, must be a single numeric scalar > 0.
#'   - For \code{"less"}, must be a single numeric scalar < 0.
#'
#' @return A list of class \code{BFvalue_t} containing:
#' \describe{
#'   \item{type}{Character string describing the test type.}
#'   \item{bf10}{Computed Bayes factor BF10.}
#'   \item{tval}{Observed t-value.}
#'   \item{df}{Degrees of freedom (currently NA / not computed).}
#'   \item{analysis_h1}{List of prior model parameters used for the alternative hypothesis.}
#'   \item{alternative}{Hypothesis tested (\code{"two.sided"}, \code{"greater"}, or \code{"less"}).}
#'   \item{ROPE}{Interval bounds used, if any.}
#'   \item{N1}{Sample size of group 1 (rounded integer).}
#'   \item{N2}{Sample size of group 2 (rounded integer).}
#' }
#'
#' @examples
#'BF10.ttest.TwoSample(
#'  tval = -1.148,
#'  N1 = 53,
#'  N2 = 48,
#'  prior_analysis = "t-distribution",
#'  location = 0,
#'  scale = 0.707,
#'  dff = 1,
#'  alternative = "two.sided",
#'  ROPE = c(-0.36,0.36))
#'
#' @export
BF10.ttest.TwoSample <- function(tval, N1, N2, prior_analysis, location, scale, dff,alternative, ROPE = NULL) {

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

  # alternative
  if(alternative %in% c("two.sided", "less", "greater") == FALSE){
    stop("Argument [alternative] should be set to either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test)")
  }

  alternative <- switch(alternative,
                        "two.sided" = "!=",
                        "less"      = "<",
                        "greater"   = ">"
  )
  # Check e if provided
  if (!is.null(ROPE)) {
    if (alternative == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative 'two.sided', argument [e] must be a numeric vector of length 2 with two distinct finite values")
      }
    }
    if (alternative == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'less', argument [e] must be a numeric scalar > 0")
      }
    }
    if (alternative == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'greater', argument [e] must be a numeric scalar < 0")
      }
    }
  }

  # Validate analysis prior
  if (missing(prior_analysis)) {
    stop("argument [prior_analysis] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }

  if (!prior_analysis %in% c("Normal","Moment","t-distribution")) {
    stop("argument [prior_analysis] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("argument [location] for analysis prior must be a numeric scalar")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
    stop("argument [scale] for analysis prior must be a positive numeric scalar (i.e., scale > 0)")
  }
  if (prior_analysis == "t-distribution") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("argument [dff] degrees of freedom for analysis prior must be a positive numeric scalar when prior_analysis='t-distribution'")
    }
  } else {
    dff <- 0
  }

  n1 <- N1
  n2 <- N2
  r <- n2 / n1
  df <-n1+n2-2

  suppressWarnings(
    if (is.null(ROPE)) {
      bf10=t2_BF10(tval, n1, r, prior_analysis, location, scale, dff,alternative)
    } else {
      bf10=t2e_BF10(tval, n1, r, prior_analysis,location, scale, dff,alternative, ROPE)
    }
  )

  type = "Indepedent-samples t-test (equal variance)"
  analysis_h1 <- list(
    prior = prior_analysis,
    location = location,
    scale = scale
  )
  if (prior_analysis == "t-distribution") {
    analysis_h1$dff <- dff
  }
  object=list(type=type,bf10=bf10,tval=tval,df=df,analysis_h1=analysis_h1,alternative=alternative,ROPE=ROPE,N1=N1,N2=N2)

  class(object) <- "BFvalue_t"

  return(object)
}



#' Bayes factor for a Bayesian correlation test
#'
#' Calculate the Bayes factor (BF10) for a correlation coefficient, either against a point null
#' or an interval null hypothesis. Supports default beta (\code{"d_beta"}), stretched beta (\code{"beta"}),
#' and normal-moment (\code{"Moment"}) priors for the alternative hypothesis.
#'
#' @param r Observed correlation coefficient. Must be a numeric scalar between -1 and 1.
#' @param n Sample size. Must be a numeric scalar greater than 3.
#' @param k Parameter for the analysis default beta prior (\code{"d_beta"}) under the alternative hypothesis.
#' @param alpha Parameter for the analysis beta prior (\code{"beta"}) under the alternative hypothesis.
#' @param beta Parameter for the analysis beta prior (\code{"beta"}) under the alternative hypothesis.
#' @param h0 Null value of the correlation. Must be a numeric scalar between -0.8 and 0.8.
#' @param alternative The direction of the alternative hypothesis being tested: two-sided (\code{"two.sided"}), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param scale Scale parameter for the analysis normal-moment prior (\code{"Moment"}). Must be > 0.
#' @param prior_analysis Analysis prior: default beta (\code{"d_beta"}), beta (\code{"beta"}), or normal-moment (\code{"Moment"}).
#' @param ROPE Optional numeric vector specifying bounds for an interval null hypothesis. For \code{"two.sided"}, must be two distinct finite values between -0.5 and 0.5. For \code{"greater"} or \code{"less"}, must satisfy additional bounds relative to \code{h0}.
#'
#' @return A list with class \code{"BFvalue_r"} containing:
#' \itemize{
#'   \item \code{type}: "correlation"
#'   \item \code{bf10}: Calculated Bayes factor BF10
#'   \item \code{h0}: Null value of the correlation
#'   \item \code{r}: Observed correlation coefficient
#'   \item \code{n}: Sample size
#'   \item \code{analysis_h1}: List of analysis prior parameters including \code{prior_analysis}, \code{k}, \code{alpha}, \code{beta}, \code{scale}
#'   \item \code{alternative}: the direction of the alternative hypothesis
#'   \item \code{ROPE}: Interval bounds if specified
#' }
#'
#' @examples
#' BF10.cor(
#'   r = 0.3930924,
#'   n = 46,
#'   prior_analysis = "d_beta",
#'   k = 1,
#'   h0 = 0,
#'   alternative = "two.sided")
#' @export

BF10.cor <- function(r, n, k, alpha, beta, h0, alternative,  scale,  prior_analysis, ROPE = NULL) {

  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < -0.8 || h0 > 0.8) {
    stop("argument [h0] null value of rho must be a single numeric scalar between -0.8 and 0.8")
  }
  location = h0


  # alternative
  if(alternative %in% c("two.sided", "less", "greater") == FALSE){
    stop("Argument [alternative] should be set to either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test)")
  }

  alternative <- switch(alternative,
                        "two.sided" = "!=",
                        "less"      = "<",
                        "greater"   = ">"
  )


  # Equivlance test or not
  interval <- if (is.null(ROPE)) 1 else 0

  if (!is.null(ROPE)) {

    if (alternative == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative 'two.sided', argument [ROPE] must be a numeric vector of length 2 with two distinct finite values")
      }
      # Additional bounds checks
      if (min(ROPE) < -0.5 || max(ROPE) > 0.5) {
        stop("For alternative 'two.sided', ROPE must satisfy min(ROPE) >= -0.5 and max(ROPE) <= 0.5")
      }
      if ((h0 + min(ROPE)) <= -1 || (h0 + min(ROPE)) >= 1) {
        stop("For alternative '!=', h0 + min(ROPE) must be between -1 and 1")
      }

    } else if (alternative == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', argument [ROPE] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (ROPE > 0.5) stop("For alternative 'greater', ROPE must be <= 0.5")
      if ((h0 + ROPE) >= 1) stop("For alternative 'greater', h0 + ROPE must be < 1")

    } else if (alternative == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', argument [ROPE] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (ROPE < -0.5) stop("For alternative 'less', ROPE must be >= -0.5")
      if ((h0 + ROPE) <= -1) stop("For alternative 'less', h0 + ROPE must be > -1")
    }

  }


  # analysis prior prior_analysis
  if (missing(prior_analysis)) {
    stop("argument [prior_analysis] for analysis prior must be one of `d_beta` (default stretched beta), `beta` (stretched beta), or `Moment` (normal-moment prior)")
  }

  # Analysis prior prior_analysis validation
  if (!prior_analysis %in% c("d_beta", "Moment", "beta")) {
    stop("argument [prior_analysis] for analysis prior must be one of `d_beta` (default stretched beta), `beta` (stretched beta), or `Moment` (normal-moment prior)")
  }

  # prior_analysis-specific checks
  if (prior_analysis == "d_beta") {
    alpha=beta=scale=NULL
    # 'd_beta' requires k to be a single numeric scalar > 0
    if (!exists("k") || !is.numeric(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
      stop("For prior_analysis 'd_beta', argument [k] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "beta") {
    k=scale=NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (!exists("alpha") || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For prior_analysis 'beta', argument [alpha] must be a single numeric scalar > 0")
    }
    if (!exists("beta") || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For prior_analysis 'beta', argument [beta] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "Moment") {
    k=alpha=beta=NULL
    # 'Moment' requires scale to be numeric scalar > 0
    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For prior_analysis 'Moment', argument [scale] must be a numeric scalar > 0")
    }
  }


  suppressWarnings(
    if (is.null(ROPE)) {
      bf10=r_BF10(r, n, k, alpha, beta, h0, alternative, location, scale, 1, prior_analysis)
    } else {
      bf10=re_BF10(r, n, k, alpha, beta, h0, alternative, location, scale, 1, prior_analysis, ROPE)
    }
  )

  type = "correlation"
  analysis_h1 <- list(
    prior = prior_analysis,
    k = k,
    alpha=alpha,
    beta=beta,
    scale=scale
  )
  object=list(type=type,bf10=bf10,h0=h0,r=r,n=n,analysis_h1=analysis_h1,alternative=alternative,ROPE=ROPE)

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
#'   (only used when \code{prior_analysis = "effectsize"}).
#' @param f_m Numeric scalar. Cohen's f effect-size parameter for the
#'   analysis prior.
#' @param prior_analysis Analysis prior under the
#'   alternative hypothesis. Must be either \code{"effectsize"} or
#'   \code{"Moment"}.
#' @param ROPE Optional numeric scalar specifying an upper bound for an interval
#'   null hypothesis. If provided, must be > 0.
#'
#' @return A list of class \code{"BFvalue_f"} containing:
#' \itemize{
#'   \item \code{fval} Input F-value
#'   \item \code{df1, df2} Degrees of freedom
#'   \item \code{ROPE} Interval bound (if specified)
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
#'   prior_analysis = "effectsize"
#' )
#'
#' @export
BF10.f.test <- function(fval, df1, df2, dff, rscale, f_m, prior_analysis, ROPE = NULL) {


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


  # analysis prior prior_analysis
  if (missing(prior_analysis)) {
    stop("argument [prior_analysis] for analysis prior should be set to either `effectsize`, or `Moment`")
  }
  if(prior_analysis %in% c("effectsize","Moment") == FALSE){
    stop("argument [prior_analysis] for analysis prior should be set to either `effectsize`, or `Moment`")
  }

  if (prior_analysis =="effectsize"){
    if (!is.numeric(rscale) || length(rscale) != 1 || !is.finite(rscale) || rscale <= 0) {
      stop("argument [rscale] scale parameter must be a positive numeric scalar")
    }
  }

  if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
    stop("argument [dff] degrees of freedom  for analysis prior must be a positive numeric scalar when prior_analysis='t-distribution'")
  }

  if (!is.numeric(f_m) || length(f_m) != 1 || !is.finite(f_m) || f_m <= 0) {
    stop("argument [f_m] Cohen's f  for analysis prior must be a positive numeric scalar")
  }

  if (prior_analysis == "Moment"){
    rscale=NULL
    if (dff < 3) {
      stop("argument [dff] degrees of freedom for Moment analysis prior must be at least 3")
    }
  }

  if (!is.null(ROPE)) {
    if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
      stop("argument [ROPE] interval bound must be a positive numeric scalar when specified")
    }
  }

  q <- df1
  m <- df1 + df2

  bf10= suppressWarnings(
    if (is.null(ROPE)) {
      F_BF(fval, q, m, dff, rscale, f_m, prior_analysis)
    } else {
      Fe_BF(fval, q, m, dff, rscale, f_m, prior_analysis, ROPE)
    }
  )

  type = "Regression/ANOVA"
  analysis_h1 <- list(
    prior = prior_analysis,
    rscale = rscale,
    f_m = f_m,
    dff=dff
  )
  object <- list(
    fval=fval,
    type = type,
    ROPE = ROPE,
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
#' Calculate the Bayes factor (BF10) for a single-proportion test, either against a point null
#' or an interval null hypothesis.
#'
#' @param x Observed number of successes (non-negative integer scalar, must be \eqn{\le n}).
#' @param n Sample size (positive integer scalar).
#' @param alpha Shape parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param beta Shape parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param h0 Null proportion value (numeric scalar between 0.1 and 0.9).
#' @param scale Scale parameter for the analysis prior (only used if \code{prior_analysis = "Moment"}).
#' @param prior_analysis the analysis prior under the alternative hypothesis:
#'   \code{"beta"} (stretched beta) or \code{"Moment"} (normal-moment prior).
#' @param alternative Hypothesis being tested: two-sided (\code{"two.sided"}), right-sided (\code{"greater"}),
#'   or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector specifying bounds for an interval null; used if interval BF is calculated.
#'
#' @return An object of class \code{"BFvalue_bin"} containing:
#'   \itemize{
#'     \item \code{bf10}: Bayes factor in favor of the alternative hypothesis.
#'     \item \code{type}: Test type ("one-proportion").
#'     \item \code{x}: Number of successes.
#'     \item \code{n}: Sample size.
#'     \item \code{h0}: Null proportion value.
#'     \item \code{analysis_h1}: List describing the analysis prior (\code{prior_analysis}, \code{alpha}, \code{beta}, \code{scale}).
#'     \item \code{alternative}: the direction of the alternative hypothesis.
#'     \item \code{ROPE}: interval null bounds (if specified).
#'   }
#'
#' @examples
#'BF10.bin.test(
#'  x = 42,
#'  n = 52,
#'  h0 = 0.5,
#'  prior_analysis = "beta",
#'  alternative = "greater",
#'  alpha = 1,
#'  beta = 1)
#'
#' @export
BF10.bin.test <- function(x, n, alpha, beta, h0, scale, prior_analysis, alternative, ROPE = NULL) {
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
    stop("argument [h0] NULL value of proportion must be a single numeric scalar between .1 and 0.9")
  }
  # alternative
  if(alternative %in% c("two.sided", "less", "greater") == FALSE){
    stop("Argument [alternative] should be set to either `less`  (left-sided test),  `two.sided` (two-sided test) or `greater` (right-sided test)")
  }

  alternative <- switch(alternative,
                        "two.sided" = "!=",
                        "less"      = "<",
                        "greater"   = ">"
  )


  # Equivlance test or not
  interval <- if (is.null(ROPE)) 1 else 0

  if (!is.null(ROPE)) {

    if (alternative == "!=") {
      # e must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative '!=', argument [ROPE] must be a numeric vector of length 2 with two distinct finite values")
      }
      # Additional bounds checks
      if (min(ROPE) < -0.5 || max(ROPE) > 0.5) {
        stop("For alternative '!=', ROPE must satisfy min(ROPE) >= -0.5 and max(ROPE) <= 0.5")
      }
      if ((h0 + min(ROPE)) <= 0 || (h0 + min(ROPE)) >= 1) {
        stop("For alternative '!=', h0 + min(ROPE) must be between 0 and 1")
      }

    } else if (alternative == ">") {
      # e must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative '>', argument [ROPE] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (ROPE > 0.5) stop("For alternative '>', ROPE must be <= 0.5")
      if ((h0 + ROPE) >= 1) stop("For alternative '>', h0 + ROPE must be < 1")

    } else if (alternative == "<") {
      # e must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative '<', argument [ROPE] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (ROPE < -0.5) stop("For alternative '<', ROPE must be >= -0.5")
      if ((h0 + ROPE) <= -1) stop("For alternative '<', h0 + ROPE must be > 0")
    }

  }


  # analysis prior prior_analysis
  if (missing(prior_analysis)) {
    stop("argument [prior_analysis] for analysis prior must be one of `beta`, or `Moment` (normal-moment prior)")
  }

  # Analysis prior prior_analysis validation
  if (!prior_analysis %in% c("Moment", "beta")) {
    stop("argument [prior_analysis] for analysis prior must be one of `beta` , or `Moment` (normal-moment prior)")
  }

  # prior_analysis-specific checks
  if (prior_analysis == "beta") {
    scale=NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (!exists("alpha") || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For prior_analysis 'beta', argument [alpha] must be a single numeric scalar > 0")
    }
    if (!exists("beta") || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For prior_analysis 'beta', argument [beta] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "Moment") {
    alpha=beta=NULL
    # 'Moment' requires scale to be numeric scalar > 0
    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For prior_analysis 'Moment', argument [scale] must be a numeric scalar > 0")
    }
  }


  bf10=suppressWarnings(
    if (is.null(ROPE)) {
      bin_BF(x, n, alpha, beta, h0, scale, prior_analysis, alternative)
    } else {
      bin_e_BF(x, n, alpha, beta, h0, scale, prior_analysis, alternative, ROPE)
    }
  )
  type = "one-proportion"
  analysis_h1 <- list(
    prior = prior_analysis,
    alpha=alpha,
    beta=beta,
    scale=scale
  )
  object=list(type=type,bf10=bf10,h0=h0,x=x,n=n,analysis_h1=analysis_h1,alternative=alternative,ROPE=ROPE)

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
#' a0 = 1,
#' b0 = 1,
#' a1 = 1,
#' b1 = 1,
#' a2 = 1,
#' b2 = 1,
#' n1 = 493,
#' n2 = 488,
#' x1 = 155,
#' x2 = 150)
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

