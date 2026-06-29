#' Sample Size Determination for the One-Sample Bayesian t-Test
#'
#' Perform sample size determination or power calculation of compelling and misleading evidence for a one-sample Bayesian t-test.
#' Can handle both point-null and interval-null hypothesis, and allows specifying
#' analysis and design priors.
#'
#'
#' @param threshold Numeric scalar. Threshold of compelling evidence (must be \eqn{\ge 1}).
#' @param type_rate Character. Either \code{"positive"} (controls true/false positive rates) or \code{"negative"} (controls true/false negative rates).
#' @param true_rate Numeric scalar. Target true positive or negative rate (between 0.6 and 0.999) for sample size determination.
#' @param false_rate Numeric scalar. Target false positive or false negative rate (between 0.001 and 0.1) for sample size determination.
#' @param N Numeric integer. Sample size for power calculation. If \code{NULL}, sample size determination is performed.
#' @param alternative Character. The direction of the alternative hypothesis : two-sided (\code{"two.sided"} ), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector or scalar. Specifies the region of practical
#'   equivalence relative to the point null value of zero.
#'   Thus, \code{ROPE} defines the interval of values considered
#'   practically equivalent to the null value.
#'
#'   For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two distinct finite values such that the first element
#'   is negative and the second element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). Example: If \code{alternative = "two.sided"} and \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[-0.2, 0.2]}.
#'
#'   For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0. Example: If \code{alternative = "greater"} and \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0, 0.2]}.
#'
#'   For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0.Example: If \code{alternative = "less"} and \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[-0.2, 0]}.
#'
#' @param prior_analysis Character. The analysis prior under the alternative hypothesis:
#'   \code{"Normal"}, \code{"Moment"} (normal-moment prior), or \code{"t-distribution"}.
#' @param location Numeric scalar. Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Numeric scalar. Scale parameter for the analysis prior under the alternative hypothesis (must be > 0).
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior under the alternative hypothesis (required if \code{prior_analysis = "t-distribution"}).
#' @param prior_design Optional Character. The design prior under the alternative hypothesis:
#'   \code{"Normal"}, \code{"Moment"} (normal-moment prior), \code{"t-distribution"}, or \code{"Point"}.
#' @param location_d Numeric scalar. Location parameter for the design prior under the alternative hypothesis.
#' @param scale_d Numeric scalar. Scale parameter for the design prior under the
#'   alternative hypothesis. Required if \code{prior_design} is
#'   \code{"Normal"}, \code{"Moment"}, or \code{"t-distribution"}; must be > 0.
#'   Not used when \code{prior_design = "Point"}.
#' @param dff_d Numeric scalar. Degrees of freedom for the design prior under the alternative hypothesis (required if \code{prior_design = "t-distribution"}).
#' @details
#' \strong{Sample Size Determination Mode (when \code{N = NULL}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size needed to achieve the desired configuration below. The user must provide:
#' \itemize{
#'   \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be \eqn{\ge 1}).
#'   \item \code{type_rate} - either \code{"positive"} to control true/false positive rates,
#'         or \code{"negative"} to control true/false negative rates.
#'   \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#'   \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence (i.e., true positive/negative rate) meets or exceeds \code{true_rate}, while the probability of misleading evidence (i.e., false positive/negative rate) does not exceed \code{false_rate}.
#'
#' \strong{Fixed-sample Analysis Mode (when \code{N} is supplied):}
#'
#' If a positive integer sample size \code{N} is provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, the arguments \code{type_rate}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Direction of the Alternative Hypothesis:}
#'
#' The argument \code{alternative} specifies the direction of the test and can be set to \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#'
#'
#' \strong{Interval Null Hypothesis:}
#'
#' The interval null hypothesis can be specified using the argument \code{ROPE},
#' which defines a region of practical equivalence around the null value of 0.
#'
#' The required form of \code{ROPE} depends on the direction of \code{alternative}:
#' \itemize{
#'
#'   \item For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two distinct finite values such that the first element
#'   is negative and the second element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). Example: If \code{alternative = "two.sided"} and \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[-0.2, 0.2]}.
#'
#'   \item For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0. Example: If \code{alternative = "greater"} and \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0, 0.2]}.
#'
#'   \item For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0.Example: If \code{alternative = "less"} and \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[-0.2, 0]}.
#'
#' }
#'
#' If \code{ROPE = NULL}, a point-null hypothesis is assumed.
#'
#' \strong{Analysis Priors:}
#'
#' The user must specify the analysis prior under the alternative hypothesis using \code{prior_analysis}:
#' \itemize{
#' \item \code{Normal} (normal prior): \code{location} with \code{scale} > 0.
#' \item \code{Moment} (normal-moment prior): \code{location} with \code{scale} > 0.
#' \item \code{t-distribution} (scaled t prior): \code{location}, \code{scale} > 0, and \code{dff} > 0.
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' The design prior under the alternative hypothesis can optionally be specified using \code{prior_design}:
#' \itemize{
#' \item \code{Normal} (normal prior): \code{location_d} with \code{scale_d} > 0.
#' \item \code{Moment} (normal-moment prior): \code{location_d} with \code{scale_d} > 0.
#' \item \code{t-distribution} (scaled t prior): \code{location_d} with \code{scale_d} > 0, and \code{dff_d} > 0.
#' \item \code{Point} (point prior): \code{location_d}.
#' }
#'
#' If \code{prior_design} is \code{NULL}, the analysis prior is used as the design prior.
#'
#' @return An object of class \code{BFpower} containing:
#'   \itemize{
#'     \item \code{type}: Character. Test type (always "One-sample t-test").
#'     \item \code{threshold}: Numeric scalar. threshold of compelling evidence.
#'     \item \code{alternative}: Character. The direction of the alternative hypothesis (\code{"two.sided"}, \code{"greater"}, or \code{"less"}).
#'     \item \code{ROPE}: Optional numeric vector or scalar for interval null bounds.
#'     \item \code{analysis_h1}: List with the analysis prior parameters:
#'       \code{prior}, \code{location}, \code{scale}, and optionally \code{dff}.
#'     \item \code{design_h1}: List with the design prior parameters:
#'       \code{prior}, \code{location}, \code{scale}, and optionally \code{dff}.
#'     \item \code{results}: Data frame of probabilities: compelling/misleading evidence.
#'     \item \code{setting}: List containing \code{mode_bf}, indicating whether
#'       sample size determination (\code{1}) or power calculation (\code{0}) is
#'       performed, and \code{same.priors}, indicating whether the design and
#'       analysis priors are the same (\code{1}) or not the same (\code{0}).
#'          }
#' @examples
#'BFpower.ttest.OneSample(
#'  threshold = 3,
#'  true_rate = 0.8,
#'  false_rate = 0.05,
#'  alternative = "two.sided",
#'  prior_analysis = "t-distribution",
#'  location = 0,
#'  scale = 0.707,
#'  dff = 1
#')
#' @export
BFpower.ttest.OneSample <- function(
    threshold,type_rate = "positive", true_rate, false_rate , N=NULL,
    alternative, ROPE=NULL,
    prior_analysis, location, scale, dff,
    prior_design=NULL, location_d, scale_d, dff_d)  {
  if (is.null(N)) {
    # mode: sample size determination
    mode_bf <- 1
  } else {
    # mode: power calculation for a fixed sample size
    mode_bf <- 0
  }

  if (mode_bf == 0) {
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 0 || N != floor(N)) {
      stop("Argument [N] sample size must be a positive integer")
    }
  }

  if (missing(alternative) || !is.character(alternative) || length(alternative) != 1 ||
      !(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Argument [alternative] should be set to either `less` (left-sided test), `two.sided` (two-sided test), or `greater` (right-sided test)")
  }


  if (!is.null(ROPE)) {

    if (alternative ==  "two.sided") {
      # ROPE must be numeric length 2, finite, distinct, with negative lower and positive upper bound
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) ||
          ROPE[1] == ROPE[2] || ROPE[1] >= 0 || ROPE[2] <= 0) {
        stop("For alternative 'two.sided', Argument [ROPE] must be a numeric vector of length 2 with ROPE[1] < 0 and ROPE[2] > 0")
      }
    }
    if (alternative == "greater") {
      # ROPE must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', Argument [ROPE] must be a numeric scalar > 0")
      }
    }

    if (alternative == "less") {
      # ROPE must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', Argument [ROPE] must be a numeric scalar < 0")
      }
    }

  }


  # analysis prior prior_analysis
  if (missing(prior_analysis) || !is.character(prior_analysis) || length(prior_analysis) != 1 ||
      !(prior_analysis %in% c("Normal", "Moment", "t-distribution"))) {
    stop("Argument [prior_analysis] for analysis prior should be set to either `Normal`, `Moment` (normal-moment prior), or `t-distribution`")
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
    dff <- NULL
  }

  # design prior

  if (!is.null(prior_design)) {
    # de_an_prior == 0 mean design prior differs from analysis prior
    # de_an_prior == 1 mean design and analysis priors are the same
    de_an_prior <- 0

    # Validate prior_design
    if (!is.character(prior_design) || length(prior_design) != 1 ||
        !(prior_design %in% c("Normal", "Moment", "t-distribution", "Point"))) {
      stop("Argument [prior_design] for design prior must be either `Normal`, `Moment`, `t-distribution`, or `Point`")
    }

    # Validate location_d in one line
    if (!is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d))
      stop("Argument [location_d] for design prior must be a numeric scalar")

    # Validate scale_d for prior_analysis that require it
    if (prior_design %in% c("Normal", "Moment", "t-distribution")) {
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0)
        stop("Argument [scale_d] for design prior must be a positive numeric scalar (i.e., scale_d > 0)")
    }

    # Validate dff_d only when prior_design = t-distribution
    if (prior_design == "t-distribution") {
      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0)
        stop("Argument [dff_d] degrees of freedom for design prior must be a positive numeric scalar when prior_design='t-distribution'")
    } else {
      dff_d <- NULL
    }
    if (prior_design == "Point") {
      scale_d <- NULL
    }
  } else {
    # de_an_prior == 1 mean design and analysis priors are the same
    de_an_prior <- 1
  }

 # desired strength of evidence
  if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold < 1) {
    stop("Argument [threshold] threshold of compelling evidence must be a numeric scalar being at least 1")
  }
  # desired power and strength of evidence
  if (mode_bf==1){
    if (!is.character(type_rate) || length(type_rate) != 1 ||
        !(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rates)")
    }
    if (!is.numeric(true_rate) || length(true_rate) != 1 ||
        !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999){
      stop("Argument [true_rate] (targeted true positive or true negative rate) must be a numeric scalar strictly greater than 0.6 and smaller than 0.999.")
    }
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("Argument [false_rate] (targeted false positive or false negative rate) must be a numeric scalar strictly greater than 0.001 and smaller than 0.1")
    }

  } else{
    # true_rate and false_rate are not needed when calculating power for a fixed sample size
    true_rate<-false_rate<-NULL
  }
  ####

  # Call appropriate table function with error handling
  results<-tryCatch(
    {
      if (is.null(ROPE)) {
        suppressWarnings(t1_Table(threshold, true_rate, prior_analysis, location, scale, dff, alternative,
                                            prior_design, location_d, scale_d, dff_d, de_an_prior, N, mode_bf, false_rate, type_rate))
      } else {
         suppressWarnings(t1e_table(threshold,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                                             prior_design,scale_d,dff_d, de_an_prior,N,mode_bf,location_d ,false_rate,type_rate ))
      }

    },
    error = function(err) {
      message("Error: Required sample size > 10,000")
      return(NaN)
    }
  )
  if (is.numeric(results) && length(results) == 1 && is.nan(results)) {
    return(NaN)
  }
  type = "One-sample t-test"
  setting  <- list(
    same.priors = de_an_prior,
    mode_bf     = mode_bf
  )
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

    # prior_design is NULL prior_design = prior_analysis
    design_h1 <- analysis_h1
  }


  object <- list(
    type = type,
    threshold = threshold,
    alternative = alternative,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    setting = setting
  )
  class(object) <- "BFpower"

  return(object)
}
#' Sample Size Determination for the Two-Sample Bayesian t-Test
#'
#' Perform sample size determination or power calculation of compelling and misleading evidence for a two-sample Bayesian t-test with equal variances.
#' Can handle both point-null and interval-null hypothesis, and allows specifying
#' analysis and design priors.
#'
#' @param threshold Numeric scalar. Threshold for compelling evidence (must be \eqn{\ge 1}).
#' @param type_rate Character. either \code{"positive"} or \code{"negative"}; determines whether to control
#'   true/false positive or true/false negative rates .
#' @param true_rate Numeric scalar. Target true positive or negative rate .
#' @param false_rate Numeric scalar. Target false positive or false negative rate (between 0.001 and 0.1) for sample size determination.
#' @param N1 Positive numeric integer. Sample size for group 1 for power calculation, used if \code{r = NULL} (must be \eqn{\ge 2}).
#' @param N2 Positive numeric integer. Sample size for group 2 for power calculation, used if \code{r = NULL} (must be \eqn{\ge 2}).
#' @param r Optional numeric scalar. Ratio of sample size \code{N2 / N1} for sample size determination (used if \code{N1} and \code{N2} are NULL).
#' @param alternative Character. The direction of the alternative hypothesis : two-sided (\code{"two.sided"} ), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector or scalar. Specifies the region of practical
#'   equivalence relative to the point null value of zero.
#'   Thus, \code{ROPE} defines the interval of values considered
#'   practically equivalent to the null value.
#'
#'   For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two distinct finite values such that the first element
#'   is negative and the second element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). Example: If \code{alternative = "two.sided"} and \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[-0.2, 0.2]}.
#'
#'   For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0. Example: If \code{alternative = "greater"} and \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0, 0.2]}.
#'
#'   For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0.Example: If \code{alternative = "less"} and \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[-0.2, 0]}.
#'
#' @param prior_analysis Character. The analysis prior under the alternative hypothesis:
#'   \code{"Normal"}, \code{"Moment"} (normal-moment prior), or \code{"t-distribution"}.
#' @param location Numeric scalar. Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Numeric scalar > 0. Scale parameter for the analysis prior under the alternative hypothesis.
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior (required if prior_analysis = \code{"t-distribution"}; ignored otherwise).
#' @param prior_design Optional Character. Design prior under the alternative:
#'   \code{"Normal"}, \code{"Moment"}(normal-moment prior), \code{"t-distribution"}, or \code{"Point"}.
#' @param location_d Numeric scalar. Location parameter for the design prior under the alternative hypothesis.
#' @param scale_d Numeric scalar. Scale parameter for the design prior under the
#'   alternative hypothesis. Required only if \code{prior_design} is
#'   \code{"Normal"}, \code{"Moment"}, or \code{"t-distribution"}; must be > 0.
#'   Not used when \code{prior_design = "Point"}.
#' @param dff_d Numeric scalar. Degrees of freedom for the design prior under the alternative hypothesis (required if \code{prior_design = "t-distribution"}; ignored otherwise).
#' @details
#' \strong{Sample size determination mode (when \code{N1 = NULL} and \code{N2 = NULL}, but \code{r} is provided):}
#'
#' If no sample size is provided, the function calculates the minimum sample size needed to achieve the desired configuration below. The user must provide:
#' \itemize{
#'   \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be at least 1).
#'   \item \code{type_rate} - either \code{"positive"} to control true/false positive rates,
#'         or \code{"negative"} to control true/false negative rates.
#'   \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#'   \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#'   \item \code{r} - the allocation ratio of group 2 to group 1 sample sizes (\code{N2/N1}).
#' }
#'
#' The function iteratively finds the smallest sample size \code{N1} and \code{N2 = r * N1} for which the probability of obtaining compelling evidence (i.e., true positive/negative rate) meets or exceeds \code{true_rate}, while the probability of misleading evidence (i.e., false positive/negative rate) does not exceed \code{false_rate}.
#'
#' \strong{Fixed-sample analysis mode (when \code{N1} and \code{N2} are supplied):}
#'
#'
#' If positive integer sample sizes \code{N1} and \code{N2} are provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, the arguments \code{type_rate}, \code{r}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Direction of the Alternative Hypothesis:}
#'
#' The argument \code{alternative} specifies the direction of the test and can be set to \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#'

#' \strong{Interval Null Hypothesis:}
#'
#' The interval null hypothesis can be specified using the argument \code{ROPE},
#' which defines a region of practical equivalence around the null value of 0.
#'
#' The required form of \code{ROPE} depends on the direction of \code{alternative}:
#' \itemize{
#'
#'   \item For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two distinct finite values such that the first element
#'   is negative and the second element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). Example: If \code{alternative = "two.sided"} and \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[-0.2, 0.2]}.
#'
#'   \item For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0. Example: If \code{alternative = "greater"} and \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0, 0.2]}.
#'
#'   \item For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0.Example: If \code{alternative = "less"} and \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[-0.2, 0]}.
#'
#' }
#'
#' If \code{ROPE = NULL}, a point-null hypothesis is assumed.
#'
#' \strong{Analysis Priors:}
#'
#' The user must specify the analysis prior under the alternative hypothesis using \code{prior_analysis}:
#' \itemize{
#' \item \code{Normal} (normal prior): \code{location} with \code{scale} > 0.
#' \item \code{Moment} (normal-moment prior): \code{scale} > 0.
#' \item \code{t-distribution} (scaled t prior): \code{location} with \code{scale} > 0, and \code{dff} > 0.
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' The design prior under the alternative hypothesis can optionally be specified using \code{prior_design}:
#' \itemize{
#' \item \code{Normal} (normal prior): \code{location_d} with \code{scale_d} > 0.
#' \item \code{Moment} (normal-moment prior): \code{location_d} with \code{scale_d} > 0.
#' \item \code{t-distribution} (scaled t prior): \code{location_d} with \code{scale_d} > 0, and \code{dff_d} > 0.
#' \item \code{Point} (point prior): \code{location_d}.
#' }
#'
#' If \code{prior_design} is \code{NULL}, the analysis prior is used as the design prior.
#'
#' @return An object of class \code{BFpower} containing:
#'   \itemize{
#'     \item \code{type}: Character. Test type (always "Independent-samples t-test (equal variance)").
#'     \item \code{threshold}: Numeric scalar. Threshold of compelling evidence.
#'     \item \code{alternative}: Character. The direction of the alternative hypothesis (\code{"two.sided"}, \code{"greater"}, or \code{"less"}).
#'     \item \code{ROPE}: Optional numeric vector or scalar. Interval bounds under the null, if any.
#'     \item \code{analysis_h1}: List with the analysis prior parameters:
#'       \code{prior}, \code{location}, \code{scale}, and optionally \code{dff}.
#'     \item \code{design_h1}:  List with the design prior parameters:
#'       \code{prior}, \code{location}, \code{scale}, and optionally \code{dff}.
#'     \item \code{results}: Data frame with probabilities of compelling/misleading evidence.
#'     \item \code{setting}: List containing \code{mode_bf}, indicating whether
#'       sample size determination (\code{1}) or power calculation (\code{0}) is
#'       performed, and \code{same.priors}, indicating whether the design and
#'       analysis priors are the same (\code{1}) or not the same (\code{0}).
#'
#'   }
#'
#' @examples
#'BFpower.ttest.TwoSample(
#'  threshold = 3,
#'  type_rate = "negative",
#'  true_rate = 0.8,
#'  false_rate = 0.05,
#'  r = 1,
#'  alternative = "two.sided",
#'  ROPE = c(-0.36, 0.36),
#'  prior_analysis = "Normal",
#'  location = -0.23,
#'  scale = 0.2,
#'  dff = 1)
#' @export
BFpower.ttest.TwoSample <- function(threshold ,type_rate = "positive", true_rate , false_rate ,
                                    N1 = NULL, N2 = NULL, r=NULL,
                                    alternative , ROPE = NULL,
                                    prior_analysis , location , scale , dff ,
                                    prior_design = NULL, location_d , scale_d , dff_d ) {
  ## ------------------------------
  ## CHECKING N1, N2, r CONSISTENCY
  ## ------------------------------

  if (is.null(N1) && is.null(N2) && is.null(r)) {
    stop(
      "Argument [r] ratio of N2/N1 must be specified for sample size calculation\n",
      "or Argument [N2] sample size for group 2 and [N1] sample size for group 1 for power calculation"
    )
  }

  # Case A: r provided > N1 and N2 must be NULL
  if (!is.null(r)) {
    # mode_bf == 1: sample size determination
    mode_bf = 1

    # r must be numeric scalar > 0
    if (!is.numeric(r) || length(r) != 1 || !is.finite(r) || r <= 0) {
      stop("Argument [r] ratio of sample size for group 2 over 1 must be a positive numeric scalar")
    }

    if (!is.null(N1) || !is.null(N2)) {
      stop("If Argument [r] is provided, both N1 and N2 must be NULL for sample size determination")
    }

  }

  # Case B: r is NULL > N1 and N2 must both be valid numeric scalars
  if (is.null(r)) {
    # mode_bf == 0: power calculation

    mode_bf = 0
    if (is.null(N1) || is.null(N2)) {
      stop("If 'r' is NULL, both N1 and N2 must be provided")
    }

    if (!is.numeric(N1) || length(N1) != 1 || !is.finite(N1) || N1 < 2 || N1 != floor(N1)) {
      stop("Argument [N1] sample size for group 1 must be a positive integer being at least 2")
    }
    if (!is.numeric(N2) || length(N2) != 1 || !is.finite(N2) || N2 < 2 || N2 != floor(N2)) {
      stop("Argument [N2] sample size for group 2 must be a positive integer being at least 2")
    }
  }

  # alternative
  if (missing(alternative) || !is.character(alternative) || length(alternative) != 1 ||
      !(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Argument [alternative] should be set to either `less` (left-sided test), `two.sided` (two-sided test), or `greater` (right-sided test)")
  }


  if (!is.null(ROPE)) {

    if (alternative ==  "two.sided") {
      # ROPE must be numeric length 2, finite, distinct, with negative lower and positive upper bound
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) ||
          ROPE[1] == ROPE[2] || ROPE[1] >= 0 || ROPE[2] <= 0) {
        stop("For alternative 'two.sided', Argument [ROPE] must be a numeric vector of length 2 with ROPE[1] < 0 and ROPE[2] > 0")
      }
    }
    if (alternative == "greater") {
      # ROPE must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', Argument [ROPE] must be a numeric scalar > 0")
      }
    }

    if (alternative == "less") {
      # ROPE must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', Argument [ROPE] must be a numeric scalar < 0")
      }
    }

  }


  # analysis prior prior_analysis
  if (missing(prior_analysis) || !is.character(prior_analysis) || length(prior_analysis) != 1 ||
      !(prior_analysis %in% c("Normal", "Moment", "t-distribution"))) {
    stop("Argument [prior_analysis] for analysis prior should be set to either `Normal`, `Moment` (normal-moment prior), or `t-distribution`")
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
    dff <- NULL
  }
  # design prior

  if (!is.null(prior_design)) {
    # design prior and analysis prior are not the same
    de_an_prior <- 0

    # Validate prior_design
    if (!is.character(prior_design) || length(prior_design) != 1 ||
        !(prior_design %in% c("Normal", "Moment", "t-distribution", "Point"))) {
      stop("Argument [prior_design] for design prior must be either `Normal`, `Moment`, `t-distribution`, or `Point`")
    }

    # Validate location_d in one line
    if (!is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d))
      stop("Argument [location_d] for design prior must be a numeric scalar")

    # Validate scale_d for prior_analysis that require it
    if (prior_design %in% c("Normal", "Moment", "t-distribution")) {
      if (!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0)
        stop("Argument [scale_d] for design prior must be a positive numeric scalar (i.e., scale_d > 0)")
    }

    # Validate dff_d only when prior_design = t-distribution
    if (prior_design == "t-distribution") {
      if (!is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0)
        stop("Argument [dff_d] degrees of freedom for design prior must be a positive numeric scalar when prior_design='t-distribution'")
    } else {
      dff_d <- NULL
    }


    # fix scale_d to 0 for point design prior
    if ( prior_design == "Point"){
      scale_d <- NULL

    }


  } else {
    # design prior and analysis prior are the same
    de_an_prior <- 1
  }

  # desired strength of evidence
  if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold < 1) {
    stop("Argument [threshold] threshold of compelling evidence must be a numeric scalar being at least 1")
  }
  # desired true/false positive rates or true/false negative rates
  if (mode_bf==1){
    if (!is.character(type_rate) || length(type_rate) != 1 ||
        !(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rates)")
    }

    if (!is.numeric(true_rate) || length(true_rate) != 1 ||
        !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999){
      stop("Argument [true_rate] (targeted true positive or true negative rate) must be a numeric scalar strictly greater than 0.6 and smaller than 0.999.")
    }
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("Argument [false_rate] (targeted false positive or false negative rate) must be a numeric scalar strictly greater than 0.001 and smaller than 0.1")
    }

  } else{
    true_rate<-false_rate<-NULL
  }

  results <- tryCatch(
    {
      if (is.null(ROPE)) {
        suppressWarnings(
          t2_Table(
            threshold, r, true_rate, prior_analysis, location, scale, dff, alternative,
            prior_design, location_d, scale_d, dff_d, de_an_prior,
            N1, N2, mode_bf, false_rate, type_rate
          )
        )
      } else {
        suppressWarnings(
          t2e_table(
            threshold, r, true_rate, prior_analysis, location, scale, dff,
            alternative, ROPE, prior_design, location_d, scale_d, dff_d,
            de_an_prior, mode_bf, N1, N2, false_rate, type_rate
          )
        )
      }
    },
    error = function(err) {
      message("Error: Required sample size > 10,000")
      return(NaN)
    }
  )
  if (is.numeric(results) && length(results) == 1 && is.nan(results)) {
    return(NaN)
  }


  type = "Independent-samples t-test (equal variance)"
  setting  <- list(
    same.priors = de_an_prior,
    mode_bf     = mode_bf
  )
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

    # prior_design is NULL, design prior = analyis prior
    design_h1 <- analysis_h1
  }


  object <- list(
    type = type,
    threshold = threshold,
    alternative = alternative,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    setting = setting
  )
  class(object) <- "BFpower"

  return(object)

}


#' Sample Size Determination for the Bayesian Correlation Test
#'
#' Perform sample size determination or power calculation of compelling and misleading evidence for a Bayesian correlation test.
#' Can handle both point-null and interval-null hypothesis, and allows specifying
#' analysis and design priors.
#'
#' @param threshold Numeric scalar. Threshold for compelling evidence (must be \eqn{\ge 1}).
#' @param type_rate Character. Either `"positive"` (controls true/false positive rates) or `"negative"` (controls true/false negative rates).
#' @param true_rate Numeric scalar. Target true positive or negative rate (between 0.6 and 0.999) for sample size determination.
#' @param false_rate Numeric scalar. Target false positive or false negative rate (between 0.001 and 0.1) for sample size determination.
#' @param N Numeric integer. Sample size for power calculation. If \code{NULL}, sample size determination is performed.
#' @param h0 Numeric scalar. Null value of the correlation. Must be a numeric scalar between -0.8 and 0.8.
#' @param alternative Character. The direction of the alternative hypothesis : two-sided (\code{"two.sided"} ), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector or scalar. Specifies the region of practical
#'   equivalence relative to the point null value of \code{h0}.
#'   Thus, \code{ROPE} defines the interval of values considered
#'   practically equivalent to the null value by \code{h0 + ROPE}.
#'
#'   For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two
#'   distinct finite values such that the first element is negative and the second
#'   element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). The resulting region of practical equivalence
#'   is \code{[h0 + ROPE[1], h0 + ROPE[2]]}. Example: If \code{h0 = 0.1}, \code{alternative = "two.sided"}, \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[-0.1, 0.3]}.
#'
#'   For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0,
#'   the interval null extends from \code{h0} to \code{h0 + ROPE}.
#'   Example: If \code{h0 = 0.1}, \code{alternative = "greater"}, \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0.1, 0.3]}.
#'
#'   For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0,
#'   the interval null extends from \code{h0 + ROPE} to \code{h0}.
#'   Example: If \code{h0 = 0.1}, \code{alternative = "less"}, \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[-0.1, 0.1]}.
#'
#' @param prior_analysis Character. The analysis prior under the alternative hypothesis:
#'        default beta (\code{"d_beta"}), beta (\code{"beta"}), or normal-moment prior (\code{"Moment"}).
#' @param k Numeric scalar. Shape parameter of the analysis default beta prior under the alternative hypothesis given  \eqn{\alpha = \beta = \frac{1}{\kappa}}{alpha = beta = 1/kappa} (required if \code{prior_analysis = "d_beta"}).
#' @param alpha Numeric scalar.  Alpha parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param beta Numeric scalar.  Beta parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param scale Numeric scalar.  Scale parameter of the analysis prior under the alternative hypothesis (required if \code{prior_analysis = "Moment"}).
#' @param prior_design Character. Design prior  under the alternative hypothesis: default beta (\code{"d_beta"}), beta (\code{"beta"}), normal-moment prior (\code{"Moment"}), or point (\code{"Point"}).
#' @param k_d Numeric scalar. Shape parameter of the design default beta prior under the alternative hypothesis given  \eqn{\alpha = \beta = \frac{1}{\kappa}}{alpha = beta = 1/kappa} (required if \code{prior_design = "d_beta"}).
#' @param alpha_d Numeric scalar. Alpha parameter of the design beta prior under the alternative hypothesis(\code{"beta"}).
#' @param beta_d Numeric scalar. Beta Parameter of the design beta prior under the alternative hypothesis(\code{"beta"}).
#' @param location_d Numeric scalar. Location parameter of the design prior under the alternative hypothesis.
#'   Required for \code{prior_design = "Moment"} and \code{prior_design = "Point"}.
#'   For \code{"Moment"}, it must satisfy \code{-1 < location_d < 1}.
#'   For \code{"Point"}, it represents the true correlation and must satisfy
#'   direction-specific constraints: for \code{alternative = "greater"},
#'   \code{h0 < location_d < 1}; for \code{alternative = "less"},
#'   \code{-1 < location_d < h0}; and for \code{alternative = "two.sided"},
#'   \code{-1 < location_d < 1} and \code{location_d != h0}.
#' @param scale_d Numeric scalar. Scale parameter of the design normal-moment prior (\code{"Moment"}) under the alternative hypothesis.
#'
#'
#'
#'@details
#' \strong{Sample Size Determination Mode (when \code{N = NULL}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size needed to achieve the desired configuration below. The user must provide:
#' \itemize{
#' \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be \eqn{\ge 1}).
#' \item \code{type_rate} - either \code{"positive"} to control true/false positive rates, or \code{"negative"} to control true/false negative rates.
#' \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#' \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence (i.e., true positive/negative rate) meets or exceeds \code{true_rate}, while the probability of misleading evidence (i.e., false positive/negative rate) does not exceed \code{false_rate}.
#'
#' \strong{Fixed-sample Analysis Mode (when \code{N} is supplied):}
#'
#' If a positive integer sample size \code{N} is provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, the arguments \code{type_rate}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Direction of the Alternative Hypothesis:}
#'
#' The argument \code{alternative} specifies the direction of the test and can be set to \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#'
#' \strong{Interval Null Hypothesis:}
#'
#' The interval null hypothesis can be specified using the argument \code{ROPE},
#' which defines a region of practical equivalence around the null value of \code{h0}.
#' Thus, \code{ROPE} defines the interval of values considered
#' practically equivalent to the null value by \code{h0 + ROPE}.
#'
#'
#' The required form of \code{ROPE} depends on the direction of \code{alternative}:
#' \itemize{
#'
#'   \item For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two
#'   distinct finite values such that the first element is negative and the second
#'   element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). The resulting region of practical equivalence
#'   is \code{[h0 + ROPE[1], h0 + ROPE[2]]}. Example: If \code{h0 = 0.1}, \code{alternative = "two.sided"}, \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[-0.1, 0.3]}.
#'
#'   \item For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0,
#'   the interval null extends from \code{h0} to \code{h0 + ROPE}.
#'   Example: If \code{h0 = 0.1}, \code{alternative = "greater"}, \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0.1, 0.3]}.
#'
#'   \item For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0,
#'   the interval null extends from \code{h0 + ROPE} to \code{h0}.
#'   Example: If \code{h0 = 0.1}, \code{alternative = "less"}, \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[-0.1, 0.1]}.
#'
#' }
#'
#'
#' If \code{ROPE = NULL}, a point-null hypothesis is assumed.
#'
#' \strong{Analysis Priors:}
#'
#' The user must specify the analysis prior under the alternative hypothesis using \code{prior_analysis}:
#' \itemize{
#' \item \code{d_beta} (default beta): \code{k} > 0.
#' \item \code{beta} (stretched beta): \code{alpha} and \code{beta} > 0.
#' \item \code{Moment} (normal-moment prior): \code{scale} > 0.
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' The design prior under the alternative hypothesis can optionally be specified
#' using \code{prior_design}:
#' \itemize{
#' \item \code{d_beta} (default beta): requires \code{k_d > 0}.
#' \item \code{beta} (stretched beta): requires \code{alpha_d > 0} and \code{beta_d > 0}.
#' \item \code{Moment} (normal-moment prior): requires \code{scale_d > 0} and
#'   \code{-1 < location_d < 1}.
#' \item \code{Point}: requires direction-specific constraints on \code{location_d}:
#'   for \code{"greater"}, \code{h0 < location_d < 1}; for \code{"less"},
#'   \code{-1 < location_d < h0}; and for \code{"two.sided"},
#'   \code{-1 < location_d < 1} and \code{location_d != h0}.
#' }
#'
#'
#' If \code{prior_design} is \code{NULL}, the analysis prior is used as the design prior.
#'
#' @return A list of class \code{BFpower} containing:
#' \itemize{
#'   \item \code{type}: Character. Test type (always "Correlation").
#'   \item \code{threshold}: Numeric scalar. Threshold of compelling evidence.
#'   \item \code{h0}: Numeric scalar. the value of correlation under the null hypothesis.
#'   \item \code{alternative}: Character. The direction of the alternative hypothesis (\code{"two.sided"}, \code{"greater"}, or \code{"less"}).
#'   \item \code{ROPE}: Optional numeric vector or scalar. Interval bounds under the null, if any.
#'   \item \code{analysis_h1}: List with the analysis prior parameters:
#'   \code{prior}, \code{k}, \code{alpha}, \code{beta},\code{location} (being the same as \code{h0} for the moment-normal prior, otherwise it is \code{NULL}), and \code{scale}.
#'   \item \code{design_h1}: List with the design prior parameters:
#'   \code{prior}, \code{k}, \code{alpha}, \code{beta}, \code{location}, and \code{scale}.
#'   \item \code{results}: Data frame with the probabilities of compelling/misleading evidence, and with the required sample size.
#'     \item \code{setting}: List containing \code{mode_bf}, indicating whether
#'       sample size determination (\code{1}) or power calculation (\code{0}) is
#'       performed, and \code{same.priors}, indicating whether the design and
#'       analysis priors are the same (\code{1}) or not the same (\code{0}).
#'         }
#'
#' @examples
#' BFpower.cor(
#'    threshold = 3,
#'    true_rate = 0.8,
#'    false_rate = 0.05,
#'    h0 = 0,
#'    alternative = "greater",
#'    prior_analysis = "d_beta",
#'    k = 1,
#'    prior_design = "Point",
#'    location_d = 0.3)
#'
#' @export
BFpower.cor<- function(threshold , type_rate="positive",true_rate, false_rate ,N = NULL,
                       h0, alternative ,  ROPE = NULL,
                       prior_analysis , k , alpha , beta , scale ,
                       prior_design = NULL,k_d , alpha_d , beta_d , location_d , scale_d) {
  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < -0.8 || h0 > 0.8) {
    stop("Argument [h0] null value of rho must be a single numeric scalar between -0.8 and 0.8")
  }
  # the correlation under h0
  location <- h0

  # sample size determination or power calculation
  if (is.null(N)) mode_bf=1 else mode_bf = 0


  # sample size
  if (mode_bf == 0) {
    # Check that N is a positive numeric integer
    # sample size
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 3 || N != floor(N)) {
      stop("Argument [N] sample size must be an integer greater than 3")
    }
  }

  # alternative
  if (missing(alternative) || !is.character(alternative) || length(alternative) != 1 ||
      !(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Argument [alternative] should be set to either `less` (left-sided test), `two.sided` (two-sided test), or `greater` (right-sided test)")
  }



  if (!is.null(ROPE)) {

    if (alternative ==  "two.sided") {

      # Basic structure + finiteness
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE))) {
        stop("For alternative 'two.sided', Argument [ROPE] must be a numeric vector of length 2 with finite values")
      }

      # Enforce sign structure: ROPE[1] < 0 < ROPE[2]
      if (ROPE[1] >= 0 || ROPE[2] <= 0) {
        stop("For alternative 'two.sided', ROPE must satisfy ROPE[1] < 0 and ROPE[2] > 0")
      }

      # Additional bounds checks
      if (ROPE[1] < -0.5 || ROPE[2] > 0.5) {
        stop("For alternative 'two.sided', ROPE must satisfy ROPE[1] >= -0.5 and ROPE[2] <= 0.5")
      }

      if ((h0 + ROPE[1]) <= -1 || (h0 + ROPE[2]) >= 1) {
        stop("For alternative 'two.sided', h0 + ROPE must lie strictly within (-1, 1)")
      }

    } else if (alternative == "greater") {
      # ROPE must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', Argument [ROPE] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (ROPE > 0.5) stop("For alternative 'greater', ROPE must be <= 0.5")
      if ((h0 + ROPE) >= 1) stop("For alternative 'greater', h0 + ROPE must be < 1")

    } else if (alternative == "less") {
      # ROPE must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', Argument [ROPE] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (ROPE < -0.5) stop("For alternative 'less', ROPE must be >= -0.5")
      if ((h0 + ROPE) <= -1) stop("For alternative 'less', h0 + ROPE must be > -1")
    }

  }

 # analysis prior
  if (missing(prior_analysis) || !is.character(prior_analysis) || length(prior_analysis) != 1 ||
      !(prior_analysis %in% c("d_beta", "Moment", "beta"))) {
    stop("Argument [prior_analysis] for analysis prior must be one of `d_beta`, `beta`, or `Moment`")
  }

  # prior_analysis-specific checks
  if (prior_analysis == "d_beta") {
    alpha<-beta<-scale<-NULL
    # 'd_beta' requires k to be a single numeric scalar > 0
    if (missing(k) || !is.numeric(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
      stop("For prior_analysis 'd_beta', Argument [k] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "beta") {
    k<-scale<-NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (missing(alpha) || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For prior_analysis 'beta', Argument [alpha] must be a single numeric scalar > 0")
    }
    if (missing(beta) || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For prior_analysis 'beta', Argument [beta] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "Moment") {
    k<-alpha<-beta<-NULL
    # 'Moment' requires scale to be numeric scalar > 0
    if (missing(scale)||!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For prior_analysis 'Moment', Argument [scale] must be a numeric scalar > 0")
    }
  }


  ##  design prior

  if (!is.null(prior_design)) {
    # design prior and analysis prior are not the same
    de_an_prior <- 0

    # Validate prior_design
    if (!is.character(prior_design) || length(prior_design) != 1 ||
        !(prior_design %in% c("d_beta", "Moment", "beta", "Point"))) {
      stop("Argument [prior_design] for design prior must be one of `d_beta`, `beta`, `Moment`, or `Point`")
    }

    # prior_analysis-specific checks
    if (prior_design == "d_beta") {
      alpha_d<-beta_d<-scale_d<-location_d<-NULL

      # 'd_beta' requires k_d to be a numeric scalar > 0
      if (missing(k_d) || !is.numeric(k_d) || length(k_d) != 1 || !is.finite(k_d) || k_d <= 0) {
        stop("For design prior 'd_beta', Argument [k_d] must be a single numeric scalar > 0")
      }
    } else if (prior_design == "beta") {
      k_d<-scale_d<-location_d<-NULL

      # 'beta' requires alpha_d and beta_d to be numeric scalars > 0
      if (missing(alpha_d) || !is.numeric(alpha_d) || length(alpha_d) != 1 || !is.finite(alpha_d) || alpha_d <= 0) {
        stop("For design prior 'beta', Argument [alpha_d] must be a single numeric scalar > 0")
      }
      if (missing(beta_d) || !is.numeric(beta_d) || length(beta_d) != 1 || !is.finite(beta_d) || beta_d <= 0) {
        stop("For design prior 'beta', Argument [beta_d] must be a single numeric scalar > 0")
      }
    } else if (prior_design == "Moment") {

      k_d <- alpha_d <- beta_d <- NULL

      if (missing(location_d) ||
          !is.numeric(location_d) ||
          length(location_d) != 1 ||
          !is.finite(location_d) ||
          location_d <= -1 ||
          location_d >= 1) {
        stop("For design prior 'Moment', Argument [location_d] must be a numeric scalar strictly between -1 and 1")
      }

      if (missing(scale_d) ||
          !is.numeric(scale_d) ||
          length(scale_d) != 1 ||
          !is.finite(scale_d) ||
          scale_d <= 0) {
        stop("For design prior 'Moment', Argument [scale_d] must be a numeric scalar > 0")
      }




    } else if (prior_design == "Point") {


      k_d <- alpha_d <- beta_d <- scale_d <- NULL

      if (missing(location_d) ||
          !is.numeric(location_d) ||
          length(location_d) != 1 ||
          !is.finite(location_d)) {
        stop("For design prior 'Point', Argument [location_d] true correlation must be a numeric scalar")
      }

      if (alternative == "two.sided") {
        if (location_d <= -1 || location_d >= 1 || location_d == h0) {
          stop("For alternative 'two.sided', Argument [location_d] must satisfy -1 < location_d < 1 and location_d != h0")
        }

      } else if (alternative == "greater") {
        if (location_d <= h0 || location_d >= 1) {
          stop("For alternative 'greater', Argument [location_d] must satisfy h0 < location_d < 1")
        }

      } else if (alternative == "less") {
        if (location_d <= -1 || location_d >= h0) {
          stop("For alternative 'less', Argument [location_d] must satisfy -1 < location_d < h0")
        }
      }



    }

  } else {
    # design and analaysis priors are the same
    de_an_prior <- 1
  }

  # desired strength of evidence
  if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold < 1) {
    stop("Argument [threshold] threshold of compelling evidence must be a numeric scalar being at least 1")
  }

  # desired true/false positive rates or negative rates
  if (mode_bf == 1) {
    if (!is.character(type_rate) || length(type_rate) != 1 ||
        !(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rates)")
    }

    if (!is.numeric(true_rate) || length(true_rate) != 1 || !is.finite(true_rate) ||
        true_rate <= 0.6 || true_rate >= 0.999) {
      stop("Argument [true_rate] (targeted true positive or true negative rate) must be a numeric scalar strictly greater than 0.6 and smaller than 0.999.")
    }

    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("Argument [false_rate] (targeted false positive or false negative rate) must be a numeric scalar strictly greater than 0.001 and smaller than 0.1")
    }

  } else {
    true_rate <- false_rate <- NULL
  }
  # degree of freedom for t-distribution, it is not used for this version of package
  dff <- dff_d <- NULL

  results <- tryCatch(
    {
      if (is.null(ROPE)) {
        suppressWarnings(
          r_table(
            threshold, true_rate, prior_analysis, k, alpha, beta, h0, location, scale, dff,
            alternative, prior_design, location_d, k_d, alpha_d, beta_d, scale_d,
            dff_d, de_an_prior, N, mode_bf, false_rate, type_rate
          )
        )
      } else {
        suppressWarnings(
          re_table(
            threshold, true_rate, prior_analysis, k, alpha, beta, h0, location, scale, dff,
            alternative, prior_design, location_d, k_d, alpha_d, beta_d, scale_d,
            dff_d, de_an_prior, N, mode_bf, false_rate, ROPE, type_rate
          )
        )
      }
    },
    error = function(err) {
      message("Error: Required sample size > 5,000")
      return(NaN)
    }
  )
  if (is.numeric(results) && length(results) == 1 && is.nan(results)) {
    return(NaN)
  }
  type = "Correlation"
  setting  <- list(
    same.priors = de_an_prior,
    mode_bf     = mode_bf
  )
  analysis_h1 <- list(
    prior = prior_analysis,
    k = k,
    alpha=alpha,
    beta=beta,
    location = if(prior_analysis=="Moment") h0 else NULL,
    scale=scale
  )

  if (!is.null(prior_design)) {

    # Base fields always included
    design_h1 <-  list(
      prior = prior_design,
      k = k_d,
      alpha=alpha_d,
      beta=beta_d,
      location=location_d,
      scale=scale_d
    )


  } else {

    # prior_design is NULL, design and analysis priors are the same
    design_h1 <- analysis_h1

  }
  object <- list(
    type =  type,
    threshold = threshold,
    h0=h0,
    alternative = alternative,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    setting = setting
  )
  class(object) <- "BFpower"

  return(object)

}

#' Sample Size Determination for the Bayesian F-Test
#'
#'
#' @description
#' Perform sample size determination or power calculation of compelling and misleading evidence for a Bayesian F-test comparing a full model to a nested reduced model.
#' Can handle both point-null and interval-null hypothesis, and allows specifying
#' analysis and design priors.
#'
#' @param threshold Numeric scalar. Threshold for compelling evidence (must be \eqn{\ge 1}).
#' @param type_rate Character. Either `"positive"` (control true/false positive rates) or
#'   `"negative"` (control true/false negative rates).
#' @param true_rate Numeric scalar. Target true positive or negative rate (between 0.6 and 0.999) for sample size determination.
#' @param false_rate Numeric scalar. Target false positive or false negative rate (between 0.001 and 0.1) for sample size determination.
#' @param N Optional integer. Sample size for power calculation. If \code{NULL}, sample size determination is performed.
#'   If \code{N} of at least \code{k + 1} is supplied, power calculation for a fixed sample size is performed.
#' @param p Numeric integer. Number of predictors in the reduced model.
#' @param k Numeric integer. Number of predictors in the full model (must satisfy \code{k > p}).
#' @param ROPE Optional numeric scalar. Specifies the upper bound of the region
#'   of practical equivalence, whose lower bound is fixed at zero. Thus,
#'   \code{ROPE} defines the interval of values considered practically
#'   equivalent to the null value. If provided, it must be positive.
#'
#'   For example, if \code{ROPE = 0.2}, then the region of practical equivalence
#'   is \code{[0, 0.2]}.
#' @param prior_analysis Character. The analysis prior model under the alternative hypothesis:
#'   \code{"effectsize"} or \code{"Moment"}.
#' @param rscale Numeric scalar. Scale parameter for the effect-size analysis prior under the alternative hypothesis (required if \code{prior_analysis = "effectsize"}).
#' @param f_m Numeric scalar. Cohen's f location parameter for the analysis prior under the alternative hypothesis.
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior under
#'   the alternative hypothesis. For the Moment prior, this must be \eqn{\ge 3}.
#' @param prior_design Character. Design prior model under the alternative hypothesis:
#'   \code{"effectsize"}, \code{"Moment"}, or \code{"Point"}.
#' @param rscale_d Numeric scalar. Scale parameter for the effect-size design prior under the alternative hypothesis (required if \code{prior_design = "effectsize"}).
#' @param f_m_d Numeric scalar.  Cohen's f location parameter for the design prior under the alternative hypothesis.
#' @param dff_d Numeric scalar. Degrees of freedom for the design prior under
#'   the alternative hypothesis. For the Moment prior, this must be \eqn{\ge 3}.
#'
#' @details
#'
#' \strong{Sample Size Determination Mode (when \code{N = NULL}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size needed to achieve the desired configuration below. The user must provide:
#' \itemize{
#'   \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be \eqn{\ge 1}).
#'   \item \code{type_rate} - either \code{"positive"} to control true/false positive rates, or \code{"negative"} to control true/false negative rates.
#'   \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#'   \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence (i.e., true positive/negative rate) meets or exceeds \code{true_rate}, while the probability of misleading evidence (i.e., false positive/negative rate) does not exceed \code{false_rate}.
#'
#' \strong{Fixed-sample Analysis Mode (when \code{N} is supplied):}
#'
#' If a positive integer sample size \code{N} is provided, the function computes the probabilities
#' of obtaining compelling or misleading evidence for that fixed sample size. In this mode,
#' \code{type_rate}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor
#' threshold \code{threshold} is used. The supplied \code{N} must satisfy
#' \code{N >= k + 1}.
#'
#' \strong{Interval Null Hypothesis:}
#'
#' The interval null hypothesis can be specified using the argument \code{ROPE},
#' which defines  the upper bound of the region
#' of practical equivalence, whose lower bound is fixed at zero. Thus,
#' \code{ROPE} defines the interval of values considered practically
#' equivalent to the null value. If provided, it must be positive.
#'
#' Example: If \code{ROPE =  0.2}, then the effective null interval is \code{[0, 0.2]}.
#'
#' If \code{ROPE = NULL}, a point-null hypothesis is assumed.

#'
#' \strong{Analysis Priors:}
#'
#' The user must specify the analysis prior under the alternative hypothesis using \code{prior_analysis}:
#' \itemize{
#' \item \code{effectsize} (effect size prior): \code{rscale > 0}, \code{f_m > 0}, and \code{dff > 0}.
#' \item \code{Moment} (normal-moment prior): \code{f_m > 0} and \code{dff >= 3}.
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' The design prior under the alternative hypothesis can optionally be specified using \code{prior_design}:
#' \itemize{
#' \item \code{effectsize} (effect size prior): \code{rscale_d > 0}, \code{f_m_d > 0}, and \code{dff_d > 0}.
#' \item \code{Moment} (normal-moment prior): \code{f_m_d > 0} and \code{dff_d >= 3}.
#' \item \code{Point} (point prior): \code{f_m_d > 0}.
#' }
#'
#' If \code{prior_design} is \code{NULL}, the analysis prior is used as the design prior.
#'
#' @return A list of class \code{BFpower} containing:
#'   \itemize{
#'     \item \code{type}: Character. Test type (always "Regression/ANOVA").
#'     \item \code{threshold}: Numeric scalar. Threshold of compelling evidence.
#'     \item \code{p}: Numeric integer. Number of predictors in the reduced model.
#'     \item \code{k}: Numeric integer. Number of predictors in the full model (must satisfy \code{k > p}).
#'     \item \code{ROPE}: Optional numeric scalar. Interval bounds under the null, if any.
#'     \item \code{analysis_h1}: List containing the analysis prior specification, including
#'       the prior distribution \code{prior}, the scale \code{rscale}, \code{f_m}, and degrees of freedom \code{dff}.
#'     \item \code{design_h1}: List containing the design prior specification, including
#'       the prior distribution \code{prior}, the scale \code{rscale}, \code{f_m}, and degrees of freedom \code{dff}.
#'     \item \code{results}: Data frame of probabilities of compelling/misleading evidence and
#'       the required or supplied sample size.
#'     \item \code{setting}: List containing \code{mode_bf}, indicating whether
#'       sample size determination (\code{1}) or power calculation (\code{0}) is
#'       performed, and \code{same.priors}, indicating whether the design and
#'       analysis priors are the same (\code{1}) or not the same (\code{0}).
#'          }
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
#'  rscale = 0.18,
#'  f_m = 0.1,
#'  dff = 3,
#'  prior_design = "Point",
#'  f_m_d = 0.1)
#'
#' @export
BFpower.f.test <- function(threshold, type_rate="positive",true_rate, false_rate ,N = NULL,
                           p , k ,ROPE = NULL,
                           prior_analysis , rscale , f_m , dff ,
                           prior_design = NULL,  rscale_d, f_m_d,dff_d) {

  ## mode
  # mode_bf == 1 : sample size determination
  # mode_bf == 0 : power calculation for fixed sample size
  if ( is.null(N)) mode_bf=1 else mode_bf = 0


  # Check p
  if (missing(p) || !is.numeric(p) || length(p) != 1 || !is.finite(p) || p < 0 || p != as.integer(p)) {
    stop("Argument [p] number of predictors in the reduced model must be a non-negative integer scalar")
  }

  # Check k
  if (missing(k) || !is.numeric(k) || length(k) != 1 || !is.finite(k) || k <= 0 || k != as.integer(k)) {
    stop("Argument [k] number of predictors in the full model must be a positive integer scalar")
  }

  # Check relation
  if (k <= p) {
    stop("Argument [k] number of predictors in the full model must be greater than argument [p] number of predictors in the reduced model")
  }


  ## Check N for fixed-sample mode
  if (mode_bf == 0) {
    #smallest possible sample size
    min_N <- k + 1

    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) ||
        N != as.integer(N) || N < min_N) {
      stop(sprintf(
        "Argument [N] sample size must be an integer scalar at least k + 1; with p = %s and k = %s, N must be at least %s",
        p, k, min_N
      ))
    }
  }


  # checking ROPE
  if (!is.null(ROPE)) {
    if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
      stop("Argument [ROPE] interval bound must be a positive numeric scalar when specified")
    }
  }


  # analysis prior prior_analysis
  if (missing(prior_analysis) || !is.character(prior_analysis) || length(prior_analysis) != 1 ||
      !(prior_analysis %in% c("effectsize", "Moment"))) {
    stop("Argument [prior_analysis] for analysis prior should be set to either `effectsize`, or `Moment`")
  }

  if (prior_analysis == "effectsize") {
    if (missing(rscale) || !is.numeric(rscale) || length(rscale) != 1 || !is.finite(rscale) || rscale <= 0) {
      stop("Argument [rscale] scale parameter must be a positive numeric scalar")
    }
  }

  if (missing(dff) || !is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
    stop("Argument [dff] degrees of freedom for analysis prior must be a positive numeric scalar")
  }

  if (missing(f_m) || !is.numeric(f_m) || length(f_m) != 1 || !is.finite(f_m) || f_m <= 0) {
    stop("Argument [f_m] Cohen's f for analysis prior must be a positive numeric scalar")
  }

  if (prior_analysis == "Moment"){
    rscale=NULL
    if (dff < 3) {
      stop("Argument [dff] degrees of freedom for Moment analysis prior must be at least 3")
    }
  }

  # design prior

  if (!is.null(prior_design)) {
    # design prior and analysis prior are not the same
    de_an_prior <- 0

    # Validate prior_design
    if (!is.character(prior_design) || length(prior_design) != 1 ||
        !(prior_design %in% c("effectsize", "Moment", "Point"))) {
      stop("Argument [prior_design] for design prior must be one of `effectsize`, `Moment`, or `Point`")
    }

    if (prior_design == "effectsize") {

      if (missing(rscale_d) || !is.numeric(rscale_d) || length(rscale_d) != 1 ||
          !is.finite(rscale_d) || rscale_d <= 0) {
        stop("Argument [rscale_d] scale parameter for design prior must be a positive numeric scalar")
      }

      if (missing(dff_d) || !is.numeric(dff_d) || length(dff_d) != 1 ||
          !is.finite(dff_d) || dff_d <= 0) {
        stop("Argument [dff_d] degrees of freedom for design prior must be a positive numeric scalar")
      }

      if (missing(f_m_d) || !is.numeric(f_m_d) || length(f_m_d) != 1 ||
          !is.finite(f_m_d) || f_m_d <= 0) {
        stop("Argument [f_m_d] Cohen's f for design prior must be a positive numeric scalar")
      }
    }

    if (prior_design == "Moment"){
      rscale_d <- NULL
      if (missing(dff_d) || !is.numeric(dff_d) || length(dff_d) != 1 || !is.finite(dff_d) || dff_d <= 0) {
        stop("Argument [dff_d] degrees of freedom for design prior must be a positive numeric scalar")
        }
      if (missing(f_m_d) ||!is.numeric(f_m_d) || length(f_m_d) != 1 || !is.finite(f_m_d) || f_m_d <= 0) {
        stop("Argument [f_m_d] Cohen's f  for design prior must be a positive numeric scalar")
      }

      if (dff_d < 3) {
        stop("Argument [dff_d] degrees of freedom for moment design prior must be at least 3")
      }
    }

    if (prior_design == "Point"){
      rscale_d<-dff_d<-NULL

      if (missing(f_m_d) || !is.numeric(f_m_d) || length(f_m_d) != 1 || !is.finite(f_m_d) || f_m_d <= 0) {
        stop("Argument [f_m_d] Cohen's f for design prior must be a positive numeric scalar")
      }

    }

  } else {
    # design prior and analysis prior are the same
    de_an_prior <- 1
  }
  # desired strength of evidence
  if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold < 1) {
    stop("Argument [threshold] threshold of compelling evidence must be a numeric scalar at least 1")
  }
  # desired true/false positive rates or negative rates
  if (mode_bf==1){
    if (!is.character(type_rate) || length(type_rate) != 1 ||
        !(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rates)")
    }

    if (!is.numeric(true_rate) || length(true_rate) != 1 ||
        !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999){
      stop("Argument [true_rate] (targeted true positive or true negative rate) must be a numeric scalar strictly greater than 0.6 and smaller than 0.999.")
    }
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("Argument [false_rate] (targeted false positive or false negative rate) must be a numeric scalar strictly greater than 0.001 and smaller than 0.1")
    }


  } else{
    true_rate <- false_rate <- NULL
  }


  results <- tryCatch({
    suppressWarnings({
      if (is.null(ROPE)) {
        f_table(threshold, true_rate, p, k, dff, rscale, f_m, prior_analysis,
                dff_d, rscale_d, f_m_d, prior_design, de_an_prior, N,
                mode_bf, false_rate, type_rate)
      } else {
        fe_table(threshold, true_rate, p, k, dff, rscale, f_m, prior_analysis,
                 dff_d, rscale_d, f_m_d, prior_design, de_an_prior, N,
                 mode_bf, ROPE, false_rate, type_rate)
      }
    })
  }, error = function(err) {
    message("Error: Required sample size > 10,000")
    return(NaN)
  })
  if (is.numeric(results) && length(results) == 1 && is.nan(results)) {
    return(NaN)
  }

  type = "Regression/ANOVA"
  setting  <- list(
    same.priors = de_an_prior,
    mode_bf     = mode_bf
  )
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

    # prior_design is NULL, design prior is the same as analysis prior
    design_h1 <- analysis_h1
  }

  object <- list(
    type = type,
    threshold = threshold,
    k=k,
    p=p,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    setting = setting
  )
  class(object) <- "BFpower"

  return(object)

}


#' Sample Size Determination for the Bayesian One-Proportion Test
#'
#' Perform sample size determination or power calculation of compelling and misleading evidence for a Bayesian test of a single proportion.
#' Can handle both point-null and interval-null hypothesis, and allows specifying
#' analysis and design priors.
#'
#' @param threshold Numeric scalar. Threshold for compelling evidence (must be \eqn{\ge 1}).
#' @param type_rate Character. Either `"positive"` (controls true/false positive rates) or `"negative"` (controls true/false negative rates).
#' @param true_rate Numeric scalar. Target true positive or negative rate (between 0.6 and 0.999) for sample size determination.
#' @param false_rate Numeric scalar. Target false positive or false negative rate (between 0.001 and 0.1) for sample size determination.
#' @param N Numeric integer. Sample size for power calculation. If \code{NULL}, sample size determination is performed.
#' @param h0 Numeric scalar. Null proportion value (numeric scalar between 0.1 and 0.9).
#' @param alternative Character. The direction of the alternative hypothesis : two-sided (\code{"two.sided"} ), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector or scalar. Specifies the region of practical
#'   equivalence relative to the point null value of \code{h0}.
#'   Thus, \code{ROPE} defines the interval of values considered
#'   practically equivalent to the null value by \code{h0 + ROPE}.
#'
#'   For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two
#'   distinct finite values such that the first element is negative and the second
#'   element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). The resulting region of practical equivalence
#'   is \code{[h0 + ROPE[1], h0 + ROPE[2]]}. Example: If \code{h0 = 0.5}, \code{alternative = "two.sided"}, \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[0.3, 0.7]}.
#'
#'   For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0,
#'   the interval null extends from \code{h0} to \code{h0 + ROPE}.
#'   Example: If \code{h0 = 0.5}, \code{alternative = "greater"}, \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0.5, 0.7]}.
#'
#'   For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0,
#'   the interval null extends from \code{h0 + ROPE} to \code{h0}.
#'   Example: If \code{h0 = 0.5}, \code{alternative = "less"}, \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[0.3, 0.5]}.
#' @param prior_analysis Character. The analysis prior under the alternative hypothesis: \code{"beta"} or \code{"Moment"} (normal-moment prior).
#' @param alpha Numeric scalar.  Alpha parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param beta Numeric scalar.  Beta parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param scale Numeric scalar.  Scale parameter of the analysis prior under the alternative hypothesis (required if \code{prior_analysis = "Moment"}).
#' @param prior_design Character. Design prior under the alternative hypothesis: \code{"beta"}, \code{"Moment"}(normal-moment prior), or \code{"Point"}.
#' @param alpha_d Numeric scalar. Alpha parameter of the design beta prior under the alternative hypothesis (required if \code{prior_design = "beta"}).
#' @param beta_d Numeric scalar. Beta Parameter of the design beta prior under the alternative hypothesis (required if \code{prior_design = "beta"}).
#' @param location_d Numeric scalar. Location parameter for the design prior under the alternative hypothesis.
#'   Required for \code{prior_design = "Moment"} and \code{prior_design = "Point"}.
#'   For \code{"Moment"}, it must satisfy \code{0 < location_d < 1}.
#'   For \code{"Point"}, it represents the true proportion and must satisfy direction-specific
#'   constraints: for \code{alternative = "greater"}, \code{h0 < location_d < 1};
#'   for \code{alternative = "less"}, \code{0 < location_d < h0}; and for
#'   \code{alternative = "two.sided"}, \code{0 < location_d < 1} and
#'   \code{location_d != h0}.
#' @param scale_d Numeric scalar.  Scale parameter of the design prior under the alternative hypothesis (required if \code{prior_design = "Moment"}).
#' @details
#'
#' \strong{Sample Size Determination Mode (when \code{N = NULL}):}
#'
#' If no sample size is provided, the function calculates the minimum sample size needed to achieve the desired configuration below. The user must provide:
#' \itemize{
#' \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be \eqn{\ge 1}).
#' \item \code{type_rate} - either \code{"positive"} to control true/false positive rates or \code{"negative"} to control true/false negative rates.
#' \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#' \item \code{false_rate} - the acceptable false positive or false negative rate (between 0.001 and 0.1).
#' }
#'
#' The function iteratively finds the smallest sample size for which the probability of obtaining compelling evidence (i.e., true positive/negative rate) meets or exceeds \code{true_rate}, while the probability of misleading evidence (i.e., false positive/negative rate) does not exceed \code{false_rate}.
#'
#' \strong{Fixed-sample Analysis Mode (when \code{N} is supplied):}
#'
#' If a positive integer sample size \code{N} is provided, the function computes the probabilities of obtaining compelling or misleading evidence for that fixed sample size. In this mode, \code{type_rate}, \code{true_rate}, and \code{false_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Direction of the Alternative Hypothesis:}
#'
#' The argument \code{alternative} specifies the direction of the test and can be set to \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#'
#' \strong{Interval Null Hypothesis:}
#'
#' The interval null hypothesis can be specified using the argument \code{ROPE},
#' which defines a region of practical equivalence around the null value of \code{h0}.
#' Thus, \code{ROPE} defines the interval of values considered
#' practically equivalent to the null value by \code{h0 + ROPE}.
#'
#' The required form of \code{ROPE} depends on the direction of \code{alternative}:
#' \itemize{
#'   \item For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two
#'   distinct finite values such that the first element is negative and the second
#'   element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). The resulting region of practical equivalence
#'   is \code{[h0 + ROPE[1], h0 + ROPE[2]]}. Example: If \code{h0 = 0.5}, \code{alternative = "two.sided"}, \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[0.3, 0.7]}.
#'
#'   \item For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0,
#'   the interval null extends from \code{h0} to \code{h0 + ROPE}.
#'   Example: If \code{h0 = 0.5}, \code{alternative = "greater"}, \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0.5, 0.7]}.
#'
#'   \item For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0,
#'   the interval null extends from \code{h0 + ROPE} to \code{h0}.
#'   Example: If \code{h0 = 0.5}, \code{alternative = "less"}, \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[0.3, 0.5]}.
#'
#' }
#'
#' If \code{ROPE = NULL}, a point-null hypothesis is assumed.
#'
#' \strong{Analysis Priors:}
#'
#' The user must specify the analysis prior under the alternative hypothesis using \code{prior_analysis}:
#' \itemize{
#' \item \code{beta} (beta prior): \code{alpha} and \code{beta} > 0.
#' \item \code{Moment} (normal-moment prior) : \code{scale} > 0.
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' The design prior under the alternative hypothesis can optionally be specified using \code{prior_design}:
#' \itemize{
#' \item \code{beta}: requires \code{alpha_d > 0} and \code{beta_d > 0}.
#' \item \code{Moment}: requires \code{scale_d > 0} and \code{0 < location_d < 1}.
#' \item \code{Point}: requires direction-specific constraints on \code{location_d}: for \code{"greater"}, \code{h0 < location_d < 1}; for \code{"less"}, \code{0 < location_d < h0}; and for \code{"two.sided"}, \code{0 < location_d < 1} and \code{location_d != h0}.
#' }
#' If \code{prior_design} is \code{NULL}, the analysis prior is used as the design prior.
#'
#' @return A list of class \code{BFpower} containing:
#' \itemize{
#'   \item \code{type}: Character. Test type (always "One-proportion").
#'   \item \code{threshold}: Numeric scalar. Compelling-evidence threshold.
#'   \item \code{h0}: Numeric scalar. Null proportion value.
#'   \item \code{alternative}: Character. The direction of the alternative hypothesis (\code{"two.sided"}, \code{"greater"}, or \code{"less"}).
#'   \item \code{ROPE}: Optional numeric vector or scalar. Interval bounds under the null, if any.
#'   \item \code{analysis_h1}: List describing the analysis prior, containing
#'     \code{prior} (prior distribution) \code{alpha} (alpha parameter),
#'     \code{beta} (beta parameter),  \code{location} (location parameter being the same as \code{h0} for the moment-normal prior),
#'      and \code{scale} (scale parameter).
#'   \item \code{design_h1}: List describing the design prior, containing, the list contains
#'     \code{prior} (prior distribution), \code{alpha} (alpha parameter), \code{beta} (beta parameter),
#'     \code{location} (location parameter), and  \code{scale} (scale parameter).
#'   \item \code{results}: Data frame of probabilities of compelling/misleading evidence and the required or supplied sample size.
#'     \item \code{setting}: List containing \code{mode_bf}, indicating whether
#'       sample size determination (\code{1}) or power calculation (\code{0}) is
#'       performed, and \code{same.priors}, indicating whether the design and
#'       analysis priors are the same (\code{1}) or not the same (\code{0}).
#'          }
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
#'   beta = 1)
#'
#' @export
BFpower.bin <- function(threshold, type_rate="positive"  ,true_rate , false_rate , N = NULL,
                        h0 ,alternative ,ROPE = NULL,
                        prior_analysis , alpha , beta , scale ,
                        prior_design = NULL, alpha_d , beta_d , location_d , scale_d ) {

  # mode
  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < .1 || h0 > 0.9) {
    stop("Argument [h0] null value of proportion must be a single numeric scalar between .1 and 0.9")
  }
  # the proportion under the null hypothesis
  # location here is the parameter of the moment analaysis prior, if used.
  location <- h0

  # mode_bf == 1 : sample size determination
  # mode_bf == 0 : power calculation for a fixed sample size
  if ( is.null(N)) mode_bf=1 else mode_bf = 0


  # sample size
  if (mode_bf == 0) {
    # Check that N is a positive numeric scalar
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 0 || N != floor(N)) {
      stop("Argument [N] sample size must be a positive integer")
    }
  }else {
    N <- NULL
    }
  # alternative
  if (missing(alternative) || !is.character(alternative) || length(alternative) != 1 ||
      !(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Argument [alternative] should be set to either `less` (left-sided test), `two.sided` (two-sided test), or `greater` (right-sided test)")
  }


  if (!is.null(ROPE)) {

    if (alternative ==  "two.sided") {
      # ROPE must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative 'two.sided', Argument [ROPE] must be a numeric vector of length 2 with two distinct finite values")
      }
      if (ROPE[1] >= 0 || ROPE[2] <= 0) {
        stop("For alternative 'two.sided', ROPE must satisfy ROPE[1] < 0 and ROPE[2] > 0")
      }
      # Additional bounds checks
      if (min(ROPE) < -0.5 || max(ROPE) > 0.5) {
        stop("For alternative 'two.sided', ROPE must satisfy min(ROPE) >= -0.5 and max(ROPE) <= 0.5")
      }
      if ((h0 + ROPE[1]) <= 0 || (h0 + ROPE[2]) >= 1) {
        stop("For alternative 'two.sided', h0 + ROPE must be between 0 and 1")
      }

    } else if (alternative == "greater") {
      # ROPE must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', Argument [ROPE] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (ROPE > 0.5) stop("For alternative 'greater', ROPE must be <= 0.5")
      if ((h0 + ROPE) >= 1) stop("For alternative 'greater', h0 + ROPE must be < 1")

    } else if (alternative == "less") {
      # ROPE must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', Argument [ROPE] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (ROPE < -0.5) stop("For alternative 'less', ROPE must be >= -0.5")
      if ((h0 + ROPE) <= 0) stop("For alternative 'less', h0 + ROPE must be > 0")
    }

  }


  # analysis prior prior_analysis
  if (missing(prior_analysis) || !is.character(prior_analysis) || length(prior_analysis) != 1 ||
      !(prior_analysis %in% c("Moment", "beta"))) {
    stop("Argument [prior_analysis] for analysis prior must be either `beta` or `Moment`")
  }
  # prior_analysis-specific checks
  if (prior_analysis == "beta") {
    scale<-NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (missing(alpha)  || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For prior_analysis 'beta', Argument [alpha] must be a single numeric scalar > 0")
    }
    if (missing(beta)  || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For prior_analysis 'beta', Argument [beta] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "Moment") {
    alpha<-beta<-NULL
    # 'Moment' requires scale to be numeric scalar > 0
    if (missing(scale) ||!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For prior_analysis 'Moment', Argument [scale] must be a numeric scalar > 0")
    }
  }
  ##  design prior

  if (!is.null(prior_design)) {
    # design prior and analysis prior are not the same
    de_an_prior <- 0

    # Validate prior_design
    if (!is.character(prior_design) || length(prior_design) != 1 ||
        !(prior_design %in% c("Moment", "beta", "Point"))) {
      stop("Argument [prior_design] for design prior must be one of `beta`, `Moment`, or `Point`")
    }
    # prior_analysis-specific checks
    if (prior_design == "beta") {
      scale_d<-location_d<-NULL

      # 'beta' requires alpha_d and beta_d to be numeric scalars > 0
      if (missing(alpha_d)  || !is.numeric(alpha_d) || length(alpha_d) != 1 || !is.finite(alpha_d) || alpha_d <= 0) {
        stop("For design prior 'beta', Argument [alpha_d] must be a single numeric scalar > 0")
      }
      if (missing(beta_d)  || !is.numeric(beta_d) || length(beta_d) != 1 || !is.finite(beta_d) || beta_d <= 0) {
        stop("For design prior 'beta', Argument [beta_d] must be a single numeric scalar > 0")
      }
    } else if (prior_design == "Moment") {
      alpha_d<-beta_d<-NULL
      if (missing(location_d) || !is.numeric(location_d) || length(location_d) != 1 ||
          !is.finite(location_d) || location_d <= 0 || location_d >= 1) {
        stop("Argument [location_d] for design prior must be a numeric scalar strictly between 0 and 1")
      }
      # 'Moment' requires scale_d numeric scalar > 0
      if (missing(scale_d) ||!is.numeric(scale_d) || length(scale_d) != 1 || !is.finite(scale_d) || scale_d <= 0) {
        stop("For design prior 'Moment', Argument [scale_d] must be a numeric scalar > 0")
      }
    } else if (prior_design == "Point") {
      alpha_d <- beta_d <- scale_d <- NULL  # Not needed for 'Point' prior

      # 'Point' prior requires location_d, which represents the true proportion
      # under the alternative hypothesis.
      if (missing(location_d) || !is.numeric(location_d) || length(location_d) != 1 || !is.finite(location_d)) {
        stop("For design prior 'Point', Argument [location_d] true proportion must be a numeric scalar")
      }

      if (alternative == "two.sided") {
        if (location_d <= 0 || location_d >= 1 || location_d == h0) {
          stop("For alternative 'two.sided', Argument [location_d] must satisfy 0 < location_d < 1 and location_d != h0")
        }
      } else if (alternative == "greater") {
        if (location_d <= h0 || location_d >= 1) {
          stop("For alternative 'greater', Argument [location_d] must satisfy h0 < location_d < 1")
        }
      } else if (alternative == "less") {
        if (location_d <= 0 || location_d >= h0) {
          stop("For alternative 'less', Argument [location_d] must satisfy 0 < location_d < h0")
        }
      }

    }

  } else {
    # design and analysis priors are the same
    de_an_prior <- 1
  }
  # desired strength of evidence
  if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold < 1) {
    stop("Argument [threshold] threshold of compelling evidence must be a numeric scalar being at least 1")
  }
  # desired power
  if (mode_bf==1){
    if (!is.character(type_rate) || length(type_rate) != 1 ||
        !(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true/false positive rates) or `negative` (controlling true/false negative rates)")
    }

    if (!is.numeric(true_rate) || length(true_rate) != 1 || !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999) {
      stop("Argument [true_rate] (targeted true positive or true negative rate) must be a numeric scalar strictly greater than 0.6 and smaller than 0.999.")
    }
    if (!is.numeric(false_rate) || length(false_rate) != 1 || !is.finite(false_rate) ||
        false_rate <= 0.001 || false_rate >= 0.1) {
      stop("Argument [false_rate] (targeted false positive or false negative rate) must be a numeric scalar strictly greater than 0.001 and smaller than 0.1")
    }

  } else{
    true_rate<-false_rate<-NULL
  }

  results <-tryCatch({
    suppressWarnings({
      if (is.null(ROPE)) {
        bin_table(threshold, true_rate, h0, alpha, beta, location, scale, prior_analysis, alternative,
                  alpha_d, beta_d, location_d, scale_d, prior_design, de_an_prior, N,
                  mode_bf, false_rate, type_rate)
      } else {
        bin_e_table(threshold, true_rate, h0, alpha, beta, location, scale, prior_analysis, alternative,
                    alpha_d, beta_d, location_d, scale_d, prior_design, de_an_prior, N,
                    mode_bf, false_rate, ROPE, type_rate)
      }
    })
  }, error = function(err) {
    message("Error: Required sample size > 10,000")
    return(NaN)
  })
  if (is.numeric(results) && length(results) == 1 && is.nan(results)) {
    return(NaN)
  }
  type = "One-proportion"
  setting  <- list(
    same.priors = de_an_prior,
    mode_bf     = mode_bf
  )
  analysis_h1 <- list(
    prior = prior_analysis,
    alpha=alpha,
    beta=beta,
    location = if(prior_analysis == "Moment") location else NULL,
    scale=scale
  )

  if (!is.null(prior_design)) {

    # Base fields always included
    design_h1 <-  list(
      prior = prior_design,
      alpha=alpha_d,
      beta=beta_d,
      location = location_d,
      scale=scale_d
    )


  } else {

    # prior_design is NULL , design prior = analysis prior
    design_h1 <- analysis_h1

  }


  object <- list(
    type = type,
    threshold = threshold,
    h0=h0,
    alternative = alternative,
    ROPE = ROPE,
    analysis_h1 = analysis_h1,
    design_h1 = design_h1,
    results = results,
    setting = setting
  )
  class(object) <- "BFpower"

  return(object)

}


#' Sample Size Determination for the Bayesian Test of Two Proportions
#'
#'
#' Perform sample size determination or power calculation of compelling and misleading evidence for a Bayesian test of two proportions.
#' Under the null hypothesis, \eqn{\theta_1 = \theta_2} and it is
#' assigned a shared analysis beta prior. Under the alternative hypothesis, \eqn{\theta_1} and
#' \eqn{\theta_2} are treated as distinct parameters and are assigned independent beta analysis priors.
#' The function supports the specification of beta and point design priors.
#'
#' @param threshold Numeric scalar. Threshold of compelling evidence (must be \eqn{\ge 1}).
#' @param type_rate Character. Choose \code{"positive"} to control true positive rate or \code{"negative"} to control true negative rate.
#' @param true_rate Numeric scalar. Target true positive rate (between 0.6 and 0.999) for sample size determination.
#' @param N1 Optional positive integer. Sample size for group 1 for power calculation. Must be supplied together with \code{N2}; if both are \code{NULL}, sample size determination is performed.
#' @param N2 Optional positive integer. Sample size for group 2 for power calculation. Must be supplied together with \code{N1}; if both are \code{NULL}, sample size determination is performed.
#' @param a0 Positive numeric scalar. Alpha parameter of the Beta analysis prior under the null hypothesis.
#' @param b0 Positive numeric scalar. Beta parameter of the Beta analysis prior under the null hypothesis.
#' @param a1 Positive numeric scalar. Alpha parameter of the Beta analysis prior for group 1 under the alternative hypothesis.
#' @param b1 Positive numeric scalar. Beta parameter of the Beta analysis prior for group 1 under the alternative hypothesis.
#' @param a2 Positive numeric scalar. Alpha parameter of the Beta analysis prior for group 2 under the alternative hypothesis.
#' @param b2 Positive numeric scalar. Beta parameter of the Beta analysis prior for group 2 under the alternative hypothesis.
#' @param prior_design_1 Character. The design prior of group 1: \code{"beta"}, \code{"Point"}, or \code{"same"} (if \code{"same"}, the design prior is identical to the analysis prior).
#' @param a1d Positive numeric scalar. Alpha parameter of the design prior for group 1 (required if \code{prior_design_1 = "beta"}).
#' @param b1d Positive numeric scalar. Beta parameter of the design prior for group 1 (required if \code{prior_design_1 = "beta"}).
#' @param dp1 Numeric scalar. True proportion for group 1 in the design prior (required if \code{prior_design_1 = "Point"}).
#' @param prior_design_2 Character. The design prior of group 2: \code{"beta"}, \code{"Point"}, or \code{"same"} (if \code{"same"}, the design prior is identical to the analysis prior).
#' @param a2d Positive numeric scalar. Alpha parameter of the design prior for group 2 (required if \code{prior_design_2 = "beta"}).
#' @param b2d Positive numeric scalar. Beta parameter of the design prior for group 2 (required if \code{prior_design_2 = "beta"}).
#' @param dp2 Numeric scalar. True proportion for group 2 in the design prior (required if \code{prior_design_2 = "Point"}).
#' @details
#'
#' \strong{Sample Size Determination Mode (when \code{N1 = NULL} and \code{N2 = NULL}):}
#'
#' If no sample sizes are provided for the two groups, the function calculates the minimum sample sizes needed to achieve the desired configuration. The user must provide:
#' \itemize{
#' \item \code{threshold} - the Bayes factor threshold for compelling evidence (must be \eqn{\ge 1}).
#' \item \code{type_rate} - either \code{"positive"} to control true positive rates or \code{"negative"} to control true negative rates.
#' \item \code{true_rate} - the targeted true positive or true negative rate (between 0.6 and 0.999).
#' }
#'
#' The function iteratively finds the smallest sample sizes for which the probability of obtaining compelling evidence (i.e., true positive/negative rate) meets or exceeds \code{true_rate}.
#'
#' \strong{Fixed-sample Analysis Mode (when \code{N1} and \code{N2} are supplied):}
#'
#' If positive integer sample sizes \code{N1} and \code{N2} are provided, the function computes the probabilities of obtaining compelling or misleading evidence for these fixed sample sizes. In this mode, \code{type_rate} and \code{true_rate} are ignored; only the Bayes factor threshold \code{threshold} is used.
#'
#' \strong{Analysis Priors:}
#'
#' The user must specify the analysis priors under the null and alternative hypotheses:
#' \itemize{
#' \item Null hypothesis: Beta prior with parameters \code{a0} and \code{b0}.
#' \item Alternative hypothesis:
#'   \itemize{
#'   \item Group 1: Beta prior with parameters \code{a1} and \code{b1}.
#'   \item Group 2: Beta prior with parameters \code{a2} and \code{b2}.
#'   }
#' }
#'
#' \strong{Design Priors (optional):}
#'
#' Design priors for the alternative hypothesis can optionally be specified:
#' \itemize{
#' \item Group 1 design prior (\code{prior_design_1}):
#'   \itemize{
#'   \item \code{"same"}: uses the corresponding analysis prior (\code{a1}, \code{b1}).
#'   \item \code{"beta"} (beta prior): requires parameters \code{a1d} and \code{b1d}.
#'   \item \code{"Point"} (point prior): requires fixed proportion \code{dp1}.
#'   }
#' \item Group 2 design prior (\code{prior_design_2}):
#'   \itemize{
#'   \item \code{"same"}: uses the corresponding analysis prior (\code{a2}, \code{b2}).
#'   \item \code{"beta"} (beta prior): requires parameters \code{a2d} and \code{b2d}.
#'   \item \code{"Point"} (point prior): requires fixed proportion \code{dp2}.
#'   }
#' }
#'
#' @return A list of class \code{BFpower} containing:
#'   \itemize{
#'     \item \code{type}: Character. Test type (always "Two-proportions").
#'     \item \code{threshold}: Numeric scalar. Threshold of compelling evidence.
#'     \item \code{analysis_h0}: List of analysis prior parameters under the null, containing \code{a} and \code{b}.
#'     \item \code{analysis_h1_theta_1}: List of analysis prior parameters for group 1 under the alternative, containing \code{a} and \code{b}.
#'     \item \code{analysis_h1_theta_2}: List of analysis prior parameters for group 2 under the alternative, containing \code{a} and \code{b}.
#'     \item \code{design_h1_theta_1}: List describing the design prior for
#'       group 1 under the alternative hypothesis. The list contains \code{prior} (prior distribution),
#'       \code{a} (alpha parameter), \code{b} (beta parameter), and
#'       \code{p} (point-prior proportion).
#'     \item \code{design_h1_theta_2}: List describing the design prior for
#'       group 2 under the alternative hypothesis. The list contains \code{prior} (prior distribution),
#'       \code{a} (alpha parameter), \code{b} (beta parameter), and
#'       \code{p} (point-prior proportion).
#'     \item \code{results}: Data frame of probabilities of compelling and misleading evidence.
#'     \item \code{grid}: Grid used internally for the computation of the results (i.e., true/false positive and negative rates) and the plot method.
#'   \item \code{mode_bf}: Numeric scalar. Indicates whether sample size determination (\code{1}) or power calculation (\code{0}) is performed. This output is only used internally in the print method.
#'   }
#' @examples
#' BFpower.props(
#' threshold = 3,
#' true_rate = 0.8,
#' a0 = 1,
#' b0 = 1,
#' a1 = 156,
#' b1 = 339,
#' a2 = 151,
#' b2 = 339)
#'
#' @export
BFpower.props <- function(threshold ,type_rate="positive", true_rate ,
                          N1 = NULL, N2 = NULL,
                          a0 , b0 , a1 , b1 ,a2 , b2 ,
                          prior_design_1 = "same",
                          a1d , b1d , dp1 ,
                          prior_design_2 = "same",
                          a2d, b2d , dp2 ) {

  if (is.null(N1) && is.null(N2)) {
    # mode : sample size determination
    mode_bf <- 1
    # ratio of sample size per group, it is always 1
    r <- 1
  } else {
    if (xor(is.null(N1), is.null(N2))) {
      stop("Arguments [N1] and [N2] must either both be NULL or both be supplied.")
    }
    # mode: power calculation for a fixed sample size
    mode_bf <- 0
  }

  if (mode_bf == 0) {
    if (!is.numeric(N1) || length(N1) != 1 || !is.finite(N1) || N1 %% 1 != 0 || N1 <= 0) {
      stop("arg [N1] sample size for group 1 must be a positive numeric scalar integer (> 0).")
    }

    if (!is.numeric(N2) || length(N2) != 1 || !is.finite(N2) || N2 %% 1 != 0 || N2 <= 0) {
      stop("arg [N2] sample size for group 2 must be a positive numeric scalar integer (> 0).")
    }

    r <- N2 / N1
  }

  if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold < 1) {
    stop("Argument [threshold] threshold of compelling evidence must be a numeric scalar being at least 1")
  }

  # If both N1 and N2 are NULL > mode_bf = 1
  # If mode_bf = 1, check target range
  if (mode_bf==1){
    if (!is.character(type_rate) || length(type_rate) != 1 ||
        !(type_rate %in% c("positive", "negative"))) {
      stop("Argument [type_rate] must be `positive` (controlling true positive rates) or `negative` (controlling true negative rates)")
    }

    if (!is.numeric(true_rate) || length(true_rate) != 1 || !is.finite(true_rate) || true_rate <= 0.6 || true_rate >= 0.999) {
      stop("Argument [true_rate] (targeted true positive or true negative rate) must be a numeric scalar strictly greater than 0.6 and smaller than 0.999.")
    }

  } else{
    true_rate <- 0
  }




  # Null hypothesis
  # Check a0 (alpha)
  if (!is.numeric(a0) || length(a0) != 1 || !is.finite(a0) || a0 <= 0) {
    stop("arg [a0] alpha for the Beta analysis prior under the null (\u03b8\u2080) must be a positive numeric scalar (> 0).")
  }

  # Check b0 (beta)
  if (!is.numeric(b0) || length(b0) != 1 || !is.finite(b0) || b0 <= 0) {
    stop("arg [b0] beta for the Beta analysis prior under the null (\u03b8\u2080) must be a positive numeric scalar (> 0).")
  }


  # alternative hypothesis theta1
  # Check a1 (alpha under the alternative)
  if (!is.numeric(a1) || length(a1) != 1 || !is.finite(a1) || a1 <= 0) {
    stop("arg [a1] alpha for the Beta analysis prior under the alternative (\u03b8\u2081) must be a positive numeric scalar (> 0).")
  }

  # Check b1 (beta under the alternative)
  if (!is.numeric(b1) || length(b1) != 1 || !is.finite(b1) || b1 <= 0) {
    stop("arg [b1] beta for the Beta analysis prior under the alternative (\u03b8\u2081) must be a positive numeric scalar (> 0).")
  }

  # alternative hypothesis theta2
  # Check a2 (alpha under the alternative)
  if (!is.numeric(a2) || length(a2) != 1 || !is.finite(a2) || a2 <= 0) {
    stop("arg [a2] alpha for the Beta analysis prior under the alternative (\u03b8\u2082) must be a positive numeric scalar (> 0).")
  }

  # Check b2 (beta under the alternative)
  if (!is.numeric(b2) || length(b2) != 1 || !is.finite(b2) || b2 <= 0) {
    stop("arg [b2] beta for the Beta analysis prior under the alternative (\u03b8\u2082) must be a positive numeric scalar (> 0).")
  }


  # Checking design priors
  if (!is.character(prior_design_1) || length(prior_design_1) != 1 ||
      !(prior_design_1 %in% c("same", "beta", "Point"))) {
    stop("arg [prior_design_1] must be one of: 'same', 'beta', 'Point'.")
  }

  if (!is.character(prior_design_2) || length(prior_design_2) != 1 ||
      !(prior_design_2 %in% c("same", "beta", "Point"))) {
    stop("arg [prior_design_2] must be one of: 'same', 'beta', 'Point'.")
  }


  # --- Check prior_design_1 assumptions for design prior on theta1 ---

  if (prior_design_1 == "same") {

    a1d <- a1
    b1d <- b1
    dp1 <- 0.5  # Dummy value to prevent error when running grid.cpp

  } else if (prior_design_1 == "beta") {

    if (!is.numeric(a1d) || length(a1d) != 1 || !is.finite(a1d) || a1d <= 0) {
      stop("arg [a1d] alpha for the Beta design prior on \u03b8\u2081 must be a positive numeric scalar (> 0).")
    }

    if (!is.numeric(b1d) || length(b1d) != 1 || !is.finite(b1d) || b1d <= 0) {
      stop("arg [b1d] beta for the Beta design prior on \u03b8\u2081 must be a positive numeric scalar (> 0).")
    }

    dp1 <- 0.5 # Dummy value to prevent error when running grid.cpp

  } else if (prior_design_1 == "Point") {

    if (!is.numeric(dp1) || length(dp1) != 1 || !is.finite(dp1)) {
      stop("arg [dp1] true \u03b8\u2081 must be a numeric scalar for prior_design_1 = 'Point'.")
    }

    if (dp1 <= 0 || dp1 >= 1) {
      stop("arg [dp1] must be > 0 and < 1 for prior_design_1 = 'Point'.")
    }

    a1d <- 1 # Dummy value to prevent error when running grid.cpp
    b1d <- 1 # Dummy value to prevent error when running grid.cpp

  } else {
    stop("arg [prior_design_1] must be one of: 'same', 'beta', 'Point'.")
  }


  # --- Check prior_design_2 assumptions for design prior on theta2 ---

  if (prior_design_2 == "same") {

    a2d <- a2
    b2d <- b2
    dp2 <- 0.5 # Dummy value to prevent error when running grid.cpp

  } else if (prior_design_2 == "beta") {

    if (!is.numeric(a2d) || length(a2d) != 1 || !is.finite(a2d) || a2d <= 0) {
      stop("arg [a2d] alpha for the Beta design prior on theta2 must be a positive numeric scalar (> 0).")
    }

    if (!is.numeric(b2d) || length(b2d) != 1 || !is.finite(b2d) || b2d <= 0) {
      stop("arg [b2d] beta for the Beta design prior on theta2 must be a positive numeric scalar (> 0).")
    }

    dp2 <- 0.5 # Dummy value to prevent error when running grid.cpp

  } else if (prior_design_2 == "Point") {

    if (!is.numeric(dp2) || length(dp2) != 1 || !is.finite(dp2)) {
      stop("arg [dp2] must be a numeric scalar for prior_design_2 = 'Point'.")
    }

    if (dp2 <= 0 || dp2 >= 1) {
      stop("arg [dp2] must be > 0 and < 1 for prior_design_2 = 'Point'.")
    }

    a2d <- 1 # Dummy value to prevent error when running grid.cpp
    b2d <- 1 # Dummy value to prevent error when running grid.cpp

  } else {
    stop("arg [prior_design_2] must be one of: 'same', 'beta', 'Point'.")
  }




  results=tryCatch({
    suppressWarnings({
      pro_table_p2(threshold, true_rate, a0, b0, a1, b1, a2, b2, r, prior_design_1,
                   a1d, b1d, dp1, prior_design_2, a2d, b2d, dp2, mode_bf, N1, N2, type_rate)
    })
  }, error = function(e) {
    message("Error: Required Sample size > 5000 per group")
    return(NaN)
  })

  if (is.numeric(results) && length(results) == 1 && is.nan(results)) {
    results_out <- NaN
    grid_out <- NULL
    return(NaN)
  } else {
    results_out <- results[[1]]
    grid_out <- results[[2]]
  }

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

  if (prior_design_1 == "same") {
    design_h1_theta_1 <- list(
      prior = "same",  # "same" is used instead of "beta" for speeding up the calculation in the back-end
      a = a1,
      b = b1,
      p = NULL
    )
  } else if (prior_design_1 == "beta") {
    design_h1_theta_1 <- list(
      prior = "beta",
      a = a1d,
      b = b1d,
      p = NULL
    )
  } else if (prior_design_1 == "Point") {
    design_h1_theta_1 <- list(
      prior = "Point",
      a = NULL,
      b = NULL,
      p = dp1
    )
  }

  if (prior_design_2 == "same") {
    design_h1_theta_2 <- list(
      prior = "same",  # "same" is used instead of "beta" for speeding up the calculation in the back-end
      a = a2,
      b = b2,
      p = NULL
    )
  } else if (prior_design_2 == "beta") {
    design_h1_theta_2 <- list(
      prior = "beta",
      a = a2d,
      b = b2d,
      p = NULL
    )
  } else if (prior_design_2 == "Point") {
    design_h1_theta_2 <- list(
      prior = "Point",
      a = NULL,
      b = NULL,
      p = dp2
    )
  }

  object <- list(
    type = type,
    threshold = threshold,
    analysis_h0=analysis_h0,
    analysis_h1_theta_1= analysis_h1_theta_1,
    analysis_h1_theta_2=analysis_h1_theta_2,
    design_h1_theta_1=design_h1_theta_1,
    design_h1_theta_2=design_h1_theta_2,
    results = results_out,
    grid= grid_out,
    mode_bf = mode_bf
  )
  class(object) <- "BFpower"

  return(object)

}

#' Bayes Factor for a One-Sample Bayesian t-Test
#'
#' Calculate the Bayes factor (BF10) for a one-sample t-test, comparing an observed t-value
#' against either a point null hypothesis or an interval null hypothesis.
#'
#' @param tval Numeric scalar. Observed t-value from the one-sample t-test.
#' @param df Numeric scalar. Degrees of freedom of the t-test (must be \eqn{\ge 1}).
#' @param alternative Character. The direction of the alternative hypothesis : two-sided (\code{"two.sided"} ), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector or scalar. Specifies the region of practical
#'   equivalence relative to the point null value of zero.
#'   Thus, \code{ROPE} defines the interval of values considered
#'   practically equivalent to the null value.
#'
#'   For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two distinct finite values such that the first element
#'   is negative and the second element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). Example: If \code{alternative = "two.sided"} and \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[-0.2, 0.2]}.
#'
#'   For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0. Example: If \code{alternative = "greater"} and \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0, 0.2]}.
#'
#'   For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0.Example: If \code{alternative = "less"} and \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[-0.2, 0]}.
#'
#' @param prior_analysis Character. The analysis prior under the alternative hypothesis: \code{"Normal"} (normal distribution),
#'   \code{"Moment"} (normal-moment prior), or \code{"t-distribution"} (t-distribution).
#' @param location Numeric scalar. Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Numeric scalar. Scale parameter for the analysis prior under the alternative hypothesis (must be > 0).
#' @param dff Numeric scalar. Degrees of freedom for the t-distribution prior under the alternative hypothesis (required if \code{prior_analysis = "t-distribution"}; must be > 0).
#'
#' @return A list of class \code{BFvalue} containing:
#'   \itemize{
#'     \item \code{type}: Character. Test type (always "One-sample t-test").
#'     \item \code{bf10}: Numeric scalar. The computed Bayes factor in favor of the alternative hypothesis relative to the null hypothesis.
#'     \item \code{tval}: Numeric scalar. Observed t-value.
#'     \item \code{df}: Numeric scalar. Degrees of freedom.
#'     \item \code{analysis_h1}: List with the analysis prior parameters:
#'       \code{prior} (prior distribution), \code{location}, \code{scale}, and
#'       optionally \code{dff}.
#'     \item \code{alternative}: Character. the direction of the alternative hypothesis.
#'     \item \code{ROPE}: Optional numeric vector or scalar. Interval bounds under the null, if any.
#'     \item \code{d}: Numeric scalar. Observed Cohen's d.
#'     \item \code{p.value}: Numeric scalar. p-value.
#'   }
#' @examples
#' BF10.ttest.OneSample(
#' tval = 2,
#' df = 50,
#' alternative = "two.sided",
#' prior_analysis = "t-distribution",
#' location = 0,
#' scale = 0.707,
#' dff = 1)
#'
#'
#' @export
BF10.ttest.OneSample <- function(tval, df, alternative, ROPE = NULL, prior_analysis, location, scale, dff) {

  ## -----------------------------
  ## Input validation
  ## -----------------------------


  # tval must be a numeric scalar
  if (!is.numeric(tval) || length(tval) != 1 || !is.finite(tval)) {
    stop("Argument [tval] observed t-value must be a numeric scalar")
  }

  # df must be numeric >= 1
  if (!is.numeric(df) || length(df) != 1 || !is.finite(df) || df < 1) {
    stop("Argument [df] degree of freedom must be a numeric scalar >= 1")
  }
  # alternative
  if (missing(alternative) || !is.character(alternative) || length(alternative) != 1 ||
      !(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Argument [alternative] should be set to either `less` (left-sided test), `two.sided` (two-sided test), or `greater` (right-sided test)")
  }

  # Check ROPE if provided
  if (!is.null(ROPE)) {
    if (alternative ==  "two.sided") {
      # ROPE must be numeric length 2, finite, distinct, with negative lower and positive upper bound
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) ||
          ROPE[1] == ROPE[2] || ROPE[1] >= 0 || ROPE[2] <= 0) {
        stop("For alternative 'two.sided', Argument [ROPE] must be a numeric vector of length 2 with ROPE[1] < 0 and ROPE[2] > 0")
      }
    }
    if (alternative == "greater") {
      # ROPE must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', Argument [ROPE] must be a numeric scalar > 0")
      }
    }
    if (alternative == "less") {
      # ROPE must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', Argument [ROPE] must be a numeric scalar < 0")
      }
    }
  }

  # Validate analysis prior
  if (missing(prior_analysis)) {
    stop("Argument [prior_analysis] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }

  if (!prior_analysis %in% c("Normal","Moment","t-distribution")) {
    stop("Argument [prior_analysis] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("Argument [location] for analysis prior must be a numeric scalar")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
    stop("Argument [scale] for analysis prior must be a positive numeric scalar (i.e., scale > 0)")
  }
  if (prior_analysis == "t-distribution") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("Argument [dff] degrees of freedom for analysis prior must be a positive numeric scalar when prior_analysis='t-distribution'")
    }
  } else {
    dff <- NULL
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

  p.value <- t.pval(tval=tval, n1=df+1, n2 = NULL, alternative, ROPE = ROPE, type = "One-sample t-test")
  analysis_h1 <- list(
    prior = prior_analysis,
    location = location,
    scale = scale
  )
  if (prior_analysis == "t-distribution") {
    analysis_h1$dff <- dff
  }
  object=list(type=type,bf10=bf10,tval=tval,df=df,analysis_h1=analysis_h1,alternative=alternative,ROPE=ROPE,d=tval/sqrt(df+1),p.value=p.value)

  class(object) <- "BFvalue"

  return(object)
}



#' Bayes Factor for a Two-Sample Bayesian t-Test
#'
#' Calculate the Bayes factor (BF10) for a two-sample independent t-test with equal variances. Supports both point-null and interval-null hypotheses.
#'
#' @param tval Numeric scalar. Observed t-value from the two-sample t-test.
#' @param N1 Numeric integer. Sample size of group 1 (must be > 2).
#' @param N2 Numeric integer. Sample size of group 2 (must be > 2).
#' @param alternative Character. The direction of the alternative hypothesis : two-sided (\code{"two.sided"} ), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector or scalar. Specifies the region of practical
#'   equivalence relative to the point null value of zero.
#'   Thus, \code{ROPE} defines the interval of values considered
#'   practically equivalent to the null value.
#'
#'   For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two distinct finite values such that the first element
#'   is negative and the second element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). Example: If \code{alternative = "two.sided"} and \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[-0.2, 0.2]}.
#'
#'   For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0. Example: If \code{alternative = "greater"} and \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0, 0.2]}.
#'
#'   For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0.Example: If \code{alternative = "less"} and \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[-0.2, 0]}.
#'
#'
#' @param prior_analysis Character. Analysis prior under the alternative hypothesis:
#'   \code{"Normal"}, \code{"Moment"} (normal-moment prior), or \code{"t-distribution"}.
#' @param location Numeric scalar. Location parameter of the analysis prior under the alternative hypothesis.
#' @param scale Numeric scalar > 0. Scale parameter of the analysis prior under the alternative hypothesis.
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior under the alternative hypothesis (required if prior_analysis = \code{"t-distribution"}; ignored otherwise).
#'
#' @return A list of class \code{BFvalue} containing:
#'   \itemize{
#'     \item \code{type}: Character. Test type (always  "Independent-samples t-test (equal variance)").
#'     \item \code{bf10}: Numeric scalar. The computed Bayes factor in favor of the alternative hypothesis relative to the null hypothesis.
#'     \item \code{tval}: Numeric scalar. Observed t-value.
#'     \item \code{df}: Numeric scalar. Degrees of freedom.
#'     \item \code{analysis_h1}: List with the analysis prior parameters:
#'       \code{prior} (prior distribution), \code{location}, \code{scale}, and
#'       optionally \code{dff}.
#'     \item \code{alternative}: Character. The direction of the alternative hypothesis (\code{"two.sided"}, \code{"greater"}, or \code{"less"}).
#'     \item \code{ROPE}: Optional numeric vector or scalar. Interval bounds under the null, if any.
#'     \item \code{N1}: Positive integer scalar. Sample size of group 1.
#'     \item \code{N2}: Positive integer scalar. Sample size of group 2.
#'     \item \code{d}: Numeric scalar. Observed Cohen's d.
#'     \item \code{p.value}: Numeric scalar. p-value.
#'   }
#' @examples
#'BF10.ttest.TwoSample(
#'  tval = -1.148,
#'  N1 = 53,
#'  N2 = 48,
#'  alternative = "two.sided",
#'  ROPE = c(-0.36,0.36),
#'  prior_analysis = "t-distribution",
#'  location = 0,
#'  scale = 0.707,
#'  dff = 1)
#'
#' @export
BF10.ttest.TwoSample <- function(tval, N1, N2,alternative, ROPE = NULL, prior_analysis, location, scale, dff) {

  ## -----------------------------
  ## Input validation
  ## -----------------------------

  # tval must be a numeric scalar
  if (!is.numeric(tval) || length(tval) != 1 || !is.finite(tval)) {
    stop("Argument [tval] observed t-value must be a numeric scalar")
  }

  if (!is.numeric(N1) || length(N1) != 1 || !is.finite(N1) || N1 <= 2 || N1 != floor(N1)) {
    stop("Argument [N1] sample size for group 1 must be an integer greater than 2")
  }
  if (!is.numeric(N2) || length(N2) != 1 || !is.finite(N2) || N2 <= 2 || N2 != floor(N2)) {
    stop("Argument [N2] sample size for group 2 must be an integer greater than 2")
  }


  # alternative
  if (missing(alternative) || !is.character(alternative) || length(alternative) != 1 ||
      !(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Argument [alternative] should be set to either `less` (left-sided test), `two.sided` (two-sided test), or `greater` (right-sided test)")
  }

  # Check ROPE if provided
  if (!is.null(ROPE)) {
    if (alternative ==  "two.sided") {
      # ROPE must be numeric length 2, finite, distinct, with negative lower and positive upper bound
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) ||
          ROPE[1] == ROPE[2] || ROPE[1] >= 0 || ROPE[2] <= 0) {
        stop("For alternative 'two.sided', Argument [ROPE] must be a numeric vector of length 2 with ROPE[1] < 0 and ROPE[2] > 0")
      }
    }

    if (alternative == "greater") {
      # ROPE must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', Argument [ROPE] must be a numeric scalar > 0")
      }
    }
    if (alternative == "less") {
      # ROPE must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', Argument [ROPE] must be a numeric scalar < 0")
      }
    }
  }

  # Validate analysis prior
  if (missing(prior_analysis)) {
    stop("Argument [prior_analysis] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }

  if (!prior_analysis %in% c("Normal","Moment","t-distribution")) {
    stop("Argument [prior_analysis] for analysis prior should be either `Normal`,  `Moment` (normal-moment prior) or `t-distribution`")
  }
  if (!is.numeric(location) || length(location) != 1 || !is.finite(location)) {
    stop("Argument [location] for analysis prior must be a numeric scalar")
  }
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
    stop("Argument [scale] for analysis prior must be a positive numeric scalar (i.e., scale > 0)")
  }
  if (prior_analysis == "t-distribution") {
    if (!is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
      stop("Argument [dff] degrees of freedom for analysis prior must be a positive numeric scalar when prior_analysis='t-distribution'")
    }
  } else {
    dff <- NULL
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

  type = "Independent-samples t-test (equal variance)"
  p.value <- t.pval(tval=tval, n1=n1, n2 = n2, alternative, ROPE = ROPE, type = "Independent-samples t-test (equal variance)")

  analysis_h1 <- list(
    prior = prior_analysis,
    location = location,
    scale = scale
  )
  if (prior_analysis == "t-distribution") {
    analysis_h1$dff <- dff
  }
  object=list(type=type,bf10=bf10,tval=tval,df=df,analysis_h1=analysis_h1,alternative=alternative,ROPE=ROPE,N1=N1,N2=N2,d=tval/sqrt((n1*n2)/(n1+n2)),p.value=p.value)

  class(object) <- "BFvalue"

  return(object)
}



#' Bayes Factor for a Bayesian Correlation Test
#'
#' Calculate the Bayes factor (BF10) for a correlation coefficient, either against a point null
#' or an interval null hypothesis. Supports default beta (\code{"d_beta"}), stretched beta (\code{"beta"}),
#' and normal-moment (\code{"Moment"}) priors for the alternative hypothesis.
#'
#' @param r Numeric scalar. Observed correlation coefficient. Must be a numeric scalar between -1 and 1.
#' @param n Numeric integer. Sample size. Must be a numeric scalar greater than 3.
#' @param h0 Numeric scalar. Null value of the correlation. Must be a numeric scalar between -0.8 and 0.8.
#' @param alternative Character. The direction of the alternative hypothesis : two-sided (\code{"two.sided"} ), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector or scalar. Specifies the region of practical
#'   equivalence relative to the point null value of \code{h0}.
#'   Thus, \code{ROPE} defines the interval of values considered
#'   practically equivalent to the null value by \code{h0 + ROPE}.
#'
#'   For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two
#'   distinct finite values such that the first element is negative and the second
#'   element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). The resulting region of practical equivalence
#'   is \code{[h0 + ROPE[1], h0 + ROPE[2]]}. Example: If \code{h0 = 0.1}, \code{alternative = "two.sided"}, \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[-0.1, 0.3]}.
#'
#'   For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0,
#'   the interval null extends from \code{h0} to \code{h0 + ROPE}.
#'   Example: If \code{h0 = 0.1}, \code{alternative = "greater"}, \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0.1, 0.3]}.
#'
#'   For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0,
#'   the interval null extends from \code{h0 + ROPE} to \code{h0}.
#'   Example: If \code{h0 = 0.1}, \code{alternative = "less"}, \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[-0.1, 0.1]}.
#'
#' @param prior_analysis Character. The analysis prior under the alternative hypothesis: default beta (\code{"d_beta"}), beta (\code{"beta"}), or normal-moment (\code{"Moment"}).
#' @param k Numeric scalar. Shape parameter for the analysis default beta prior under the alternative hypothesis given  \eqn{\alpha = \beta = \frac{1}{\kappa}}{alpha = beta = 1/kappa} (required if \code{prior_analysis = "d_beta"}).
#' @param alpha Numeric scalar.  Alpha parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param beta Numeric scalar.  Beta parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param scale Numeric scalar.  Scale parameter for the analysis prior (required if \code{prior_analysis = "Moment"}).
#'
#' @return A list with class \code{BFvalue} containing:
#' \itemize{
#'   \item \code{type}: Character. Test type (always "Correlation").
#'   \item \code{bf10}: Numeric scalar. The computed Bayes factor in favor of the alternative hypothesis relative to the null hypothesis.
#'   \item \code{h0}: Numeric scalar. Null value of the correlation.
#'   \item \code{r}: Numeric scalar. Observed correlation coefficient.
#'   \item \code{n}: Positive integer scalar. Sample size.
#'   \item \code{analysis_h1}: List with the analysis prior parameters:
#'   \code{prior}, \code{k}, \code{alpha}, \code{beta},\code{location} (being the same as \code{h0} for the moment-normal prior, otherwise it is \code{NULL}), and \code{scale}.
#'   \item \code{alternative}: Character. The direction of the alternative hypothesis (\code{"two.sided"}, \code{"greater"}, or \code{"less"}).
#'   \item \code{ROPE}: Optional numeric vector or scalar. Interval bounds under the null, if any.
#'   \item \code{p.value}: Numeric scalar. p-value.
#' }
#'
#' @examples
#' BF10.cor(
#'   r = 0.393,
#'   n = 46,
#'   h0 = 0,
#'   alternative = "two.sided",
#'   prior_analysis = "d_beta",
#'   k = 1)
#' @export

BF10.cor <- function(r, n, h0, alternative, ROPE = NULL,  prior_analysis, k, alpha, beta,  scale) {

  # checking the obsered correlation
  if (!is.numeric(r) || length(r) != 1 || is.na(r) || !is.finite(r) || r <= -1 || r >= 1) {
    stop("Argument [r] observed correlation must be a single numeric scalar strictly between -1 and 1")
  }

  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < -0.8 || h0 > 0.8) {
    stop("Argument [h0] null value of rho must be a single numeric scalar between -0.8 and 0.8")
  }
  location = h0
  # sample size
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 3 || n != floor(n)) {
    stop("Argument [n] sample size must be an integer greater than 3")
  }
  # alternative
  if (missing(alternative) || !is.character(alternative) || length(alternative) != 1 ||
      !(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Argument [alternative] should be set to either `less` (left-sided test), `two.sided` (two-sided test), or `greater` (right-sided test)")
  }



  if (!is.null(ROPE)) {

    if (alternative ==  "two.sided") {
      # Basic structure + finiteness
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE))) {
        stop("For alternative 'two.sided', Argument [ROPE] must be a numeric vector of length 2 with finite values")
      }


      # Enforce sign structure: ROPE[1] < 0 < ROPE[2]
      if (ROPE[1] >= 0 || ROPE[2] <= 0) {
        stop("For alternative 'two.sided', ROPE must satisfy ROPE[1] < 0 and ROPE[2] > 0")
      }

      # Additional bounds checks
      if (ROPE[1] < -0.5 || ROPE[2] > 0.5) {
        stop("For alternative 'two.sided', ROPE must satisfy ROPE[1] >= -0.5 and ROPE[2] <= 0.5")
      }

      if ((h0 + ROPE[1]) <= -1 || (h0 + ROPE[2]) >= 1) {
        stop("For alternative 'two.sided', h0 + ROPE must lie strictly within (-1, 1)")
      }

    } else if (alternative == "greater") {
      # ROPE must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', Argument [ROPE] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (ROPE > 0.5) stop("For alternative 'greater', ROPE must be <= 0.5")
      if ((h0 + ROPE) >= 1) stop("For alternative 'greater', h0 + ROPE must be < 1")

    } else if (alternative == "less") {
      # e must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', Argument [ROPE] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (ROPE < -0.5) stop("For alternative 'less', ROPE must be >= -0.5")
      if ((h0 + ROPE) <= -1) stop("For alternative 'less', h0 + ROPE must be > -1")
    }

  }


  # analysis prior
  if (missing(prior_analysis) || !is.character(prior_analysis) || length(prior_analysis) != 1) {
    stop("Argument [prior_analysis] for analysis prior must be one of `d_beta`, `beta`, or `Moment`")
  }

  if (!prior_analysis %in% c("d_beta", "Moment", "beta")) {
    stop("Argument [prior_analysis] for analysis prior must be one of `d_beta`, `beta`, or `Moment`")
  }

  # prior_analysis-specific checks
  if (prior_analysis == "d_beta") {
    alpha<-beta<-scale<-NULL
    # 'd_beta' requires k to be a single numeric scalar > 0
    if (missing(k) || !is.numeric(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
      stop("For prior_analysis 'd_beta', Argument [k] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "beta") {
    k<-scale<-NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (missing(alpha) || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For prior_analysis 'beta', Argument [alpha] must be a single numeric scalar > 0")
    }
    if (missing(beta) || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For prior_analysis 'beta', Argument [beta] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "Moment") {
    k<-alpha<-beta<-NULL
    # 'Moment' requires scale to be numeric scalar > 0
    if (missing(scale)||!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For prior_analysis 'Moment', Argument [scale] must be a numeric scalar > 0")
    }
  }
  suppressWarnings(
    if (is.null(ROPE)) {
      bf10=r_BF10(r, n, k, alpha, beta, h0, alternative, location, scale, 1, prior_analysis)
    } else {
      bf10=re_BF10(r, n, k, alpha, beta, h0, alternative, location, scale, 1, prior_analysis, ROPE)
    }
  )

  type = "Correlation"
  p.value <- r.pval(r, n,h0, alternative , ROPE)

  analysis_h1 <- list(
    prior = prior_analysis,
    location = location,
    k = k,
    alpha=alpha,
    beta=beta,
    scale=scale
  )
  object=list(type=type,bf10=bf10,h0=h0,r=r,n=n,analysis_h1=analysis_h1,alternative=alternative,ROPE=ROPE,p.value=p.value)

  class(object) <- "BFvalue"
  return(object)
}




#' Bayes Factor for a Bayesian F-Test
#'
#' Calculate the Bayes factor (BF10) for an F-test, comparing a full model to a
#' reduced model under either an effect-size prior or a moment prior.
#' Optionally, an interval null hypothesis can be specified.
#'
#' @param fval Numeric scalar. Observed F statistic (must be \eqn{\ge 0}).
#' @param df1 Numeric scalar. Numerator degrees of freedom (must be > 0).
#' @param df2 Numeric scalar. Denominator degrees of freedom (must be > 0).
#' @param ROPE Optional numeric scalar. Specifies the upper bound of the region
#'   of practical equivalence, whose lower bound is fixed at zero. Thus,
#'   \code{ROPE} defines the interval of values considered practically
#'   equivalent to the null value. If provided, it must be positive.
#'
#'   For example, if \code{ROPE = 0.2}, then the region of practical equivalence
#'   is \code{[0, 0.2]}.
#' @param prior_analysis Character. The analysis prior under the alternative hypothesis: \code{"effectsize"} or \code{"Moment"}.
#' @param rscale Numeric scalar. Scale parameter for the effect-size analysis prior under the alternative hypothesis (required if \code{prior_analysis = "effectsize"}).
#' @param f_m Numeric scalar. Cohen's f location parameter for the analysis prior under the alternative hypothesis.
#' @param dff Numeric scalar. Degrees of freedom for the analysis prior under
#'   the alternative hypothesis. For the Moment prior, this must be \eqn{\ge 3}.
#'
#' @return A list of class \code{BFvalue} containing:
#'   \itemize{
#'     \item \code{type}: Character. Test type (always "Regression/ANOVA").
#'     \item \code{bf10}: Numeric scalar. The computed Bayes factor in favor of the alternative hypothesis relative to the null hypothesis.
#'     \item \code{fval}: Numeric scalar. Input F-value.
#'     \item \code{df1}, \code{df2}: Numeric scalar. Degrees of freedom.
#'     \item \code{analysis_h1}: List containing the analysis prior specification, including
#'       the prior distribution, the scale \code{rscale}, \code{f_m}, and degrees of freedom \code{dff}.
#'     \item \code{ROPE}: Optional numeric scalar. Interval bounds under the null, if any.
#'     \item \code{p.value}: Numeric scalar. p-value.
#'   }
#' @examples
#' BF10.f.test(
#'   fval = 4.5,
#'   df1 = 2,
#'   df2 = 12,
#'   dff = 12,
#'   prior_analysis = "effectsize",
#'   rscale = 0.707,
#'   f_m = 0.1)
#'
#' @export
BF10.f.test <- function(fval, df1, df2, ROPE = NULL, prior_analysis, rscale, f_m, dff) {


  ## Check fval
  if (!is.numeric(fval) || length(fval) != 1 || !is.finite(fval) || fval < 0) {
    stop("Argument [fval] F-value must be a numeric scalar greater than or equal to 0")
  }

  ## Check df1
  if (!is.numeric(df1) || length(df1) != 1 || !is.finite(df1) || df1 <= 0) {
    stop("Argument [df1] numerator degrees of freedom must be a positive numeric scalar")
  }

  ## Check df2
  if (!is.numeric(df2) || length(df2) != 1 || !is.finite(df2) || df2 <= 0) {
    stop("Argument [df2] denominator degrees of freedom must be a positive numeric scalar")
  }


  # analysis prior prior_analysis
  if (missing(prior_analysis) || !is.character(prior_analysis) || length(prior_analysis) != 1 ||
      !(prior_analysis %in% c("effectsize", "Moment"))) {
    stop("Argument [prior_analysis] for analysis prior should be set to either `effectsize`, or `Moment`")
  }

  if (prior_analysis == "effectsize") {
    if (missing(rscale) || !is.numeric(rscale) || length(rscale) != 1 || !is.finite(rscale) || rscale <= 0) {
      stop("Argument [rscale] scale parameter must be a positive numeric scalar")
    }
  }

  if (missing(dff) || !is.numeric(dff) || length(dff) != 1 || !is.finite(dff) || dff <= 0) {
    stop("Argument [dff] degrees of freedom for analysis prior must be a positive numeric scalar")
  }

  if (missing(f_m) || !is.numeric(f_m) || length(f_m) != 1 || !is.finite(f_m) || f_m <= 0) {
    stop("Argument [f_m] Cohen's f for analysis prior must be a positive numeric scalar")
  }

  if (prior_analysis == "Moment"){
    rscale <- NULL
    if (dff < 3) {
      stop("Argument [dff] degrees of freedom for Moment analysis prior must be at least 3")
    }
  }

  if (!is.null(ROPE)) {
    if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
      stop("Argument [ROPE] interval bound must be a positive numeric scalar when specified")
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
  p.value <- f.pval(fval, df1,df2,ROPE)
  analysis_h1 <- list(
    prior = prior_analysis,
    rscale = rscale,
    f_m = f_m,
    dff=dff
  )
  object <- list(
    type = type,
    bf10=bf10,
    fval=fval,
    df1=df1,
    df2=df2,
    analysis_h1 = analysis_h1,
    ROPE = ROPE,
    p.value=p.value
  )
  class(object) <- "BFvalue"
  return(object)

}


#' Bayes Factor for a Bayesian One-Proportion Test
#'
#' Calculate the Bayes factor (BF10) for a one-proportion test, either against a point null
#' or an interval null hypothesis.
#'
#' @param n Numeric integer. Sample size (positive integer scalar).
#' @param x Numeric integer. Observed number of successes (non-negative integer scalar, must be \eqn{\ge 0} and \eqn{\le n} ).
#' @param h0 Numeric scalar.  Null proportion value (numeric scalar between 0.1 and 0.9).
#' @param alternative Character. The direction of the alternative hypothesis : two-sided (\code{"two.sided"} ), right-sided (\code{"greater"}), or left-sided (\code{"less"}).
#' @param ROPE Optional numeric vector or scalar. Specifies the region of practical
#'   equivalence relative to the point null value of \code{h0}.
#'   Thus, \code{ROPE} defines the interval of values considered
#'   practically equivalent to the null value by \code{h0 + ROPE}.
#'
#'   For \code{alternative = "two.sided"}, argument \code{ROPE} must be a numeric vector of length 2 with two
#'   distinct finite values such that the first element is negative and the second
#'   element is positive (i.e., \code{ROPE[1] < 0 < ROPE[2]}). The resulting region of practical equivalence
#'   is \code{[h0 + ROPE[1], h0 + ROPE[2]]}. Example: If \code{h0 = 0.5}, \code{alternative = "two.sided"}, \code{ROPE = c(-0.2, 0.2)},
#'   then the region of practical equivalence is \code{[0.3, 0.7]}.
#'
#'   For \code{alternative = "greater"}, argument \code{ROPE} must be a numeric scalar > 0,
#'   the interval null extends from \code{h0} to \code{h0 + ROPE}.
#'   Example: If \code{h0 = 0.5}, \code{alternative = "greater"}, \code{ROPE = 0.2},
#'   then the region of practical equivalence is \code{[0.5, 0.7]}.
#'
#'   For \code{alternative = "less"}, argument \code{ROPE} must be a numeric scalar < 0,
#'   the interval null extends from \code{h0 + ROPE} to \code{h0}.
#'   Example: If \code{h0 = 0.5}, \code{alternative = "less"}, \code{ROPE = -0.2},
#'   then the region of practical equivalence is \code{[0.3, 0.5]}.
#'
#' @param prior_analysis Character. The analysis prior under the alternative hypothesis:
#'   \code{"beta"} or \code{"Moment"} (normal-moment prior).
#' @param alpha Numeric scalar.  Alpha parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param beta Numeric scalar.  Beta parameter of the analysis beta prior under the alternative hypothesis
#'   (required if \code{prior_analysis = "beta"}).
#' @param scale Numeric scalar.  Scale parameter for the analysis prior (required if \code{prior_analysis = "Moment"}).
#'
#' @return An object of class \code{BFvalue} containing:
#'   \itemize{
#'     \item \code{type}: Character. Test type (always "One-proportion").
#'     \item \code{bf10}: Numeric scalar. The computed Bayes factor in favor of the alternative hypothesis relative to the null hypothesis.
#'     \item \code{h0}: Numeric scalar. Null proportion value.
#'     \item \code{x}: Non-negative integer scalar.  Number of successes.
#'     \item \code{n}: Positive integer scalar. Sample size.
#'    \item \code{analysis_h1}: List describing the analysis prior, containing
#'     \code{prior} (prior distribution) \code{alpha} (alpha parameter),
#'     \code{beta} (beta parameter),  \code{location} (location parameter being the same as \code{h0} for the moment-normal prior),
#'     \item \code{alternative}: Character. The direction of the alternative hypothesis (\code{"two.sided"}, \code{"greater"}, or \code{"less"}).
#'     \item \code{ROPE}: Optional numeric vector or scalar. Interval bounds under the null, if any.
#'     \item \code{p.value}: Numeric scalar. p-value.
#'   }
#'
#' @examples
#'BF10.bin.test(
#'  n = 52,
#'  x = 42,
#'  h0 = 0.5,
#'  alternative = "greater",
#'  prior_analysis = "beta",
#'  alpha = 1,
#'  beta = 1)
#'
#' @export
BF10.bin.test <- function( n,x, h0, alternative, ROPE = NULL, prior_analysis, alpha, beta,  scale) {
  # Check n
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0 || n != floor(n)) {
    stop("Argument [n] sample size must be a positive integer")
  }

  # Check x
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x != floor(x)) {
    stop("Argument [x] number of successes must be a non-negative integer")
  }

  # Check relation
  if (x > n) {
    stop("Argument [x] number of successes cannot exceed total sample size [n]")
  }
  # mode
  # Check h0
  if (!is.numeric(h0) || length(h0) != 1 || !is.finite(h0) || h0 < .1 || h0 > 0.9) {
    stop("Argument [h0] null value of proportion must be a single numeric scalar between .1 and 0.9")
  }
  # alternative
  if (missing(alternative) || !is.character(alternative) || length(alternative) != 1 ||
      !(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Argument [alternative] should be set to either `less` (left-sided test), `two.sided` (two-sided test), or `greater` (right-sided test)")
  }


  if (!is.null(ROPE)) {

    if (alternative ==  "two.sided") {

      # ROPE must be a numeric vector of length 2, both finite and distinct
      if (!is.numeric(ROPE) || length(ROPE) != 2 || any(!is.finite(ROPE)) || ROPE[1] == ROPE[2]) {
        stop("For alternative 'two.sided', Argument [ROPE] must be a numeric vector of length 2 with two distinct finite values")
      }
      if (ROPE[1] >= 0 || ROPE[2] <= 0) {
        stop("For alternative 'two.sided', ROPE must satisfy ROPE[1] < 0 and ROPE[2] > 0")
      }

      # Additional bounds checks
      if (min(ROPE) < -0.5 || max(ROPE) > 0.5) {
        stop("For alternative 'two.sided', ROPE must satisfy min(ROPE) >= -0.5 and max(ROPE) <= 0.5")
      }
      if ((h0 + ROPE[1]) <= 0 || (h0 + ROPE[2]) >= 1) {
        stop("For alternative 'two.sided', h0 + ROPE must be between 0 and 1")
      }

    } else if (alternative == "greater") {
      # ROPE must be a numeric scalar > 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE <= 0) {
        stop("For alternative 'greater', Argument [ROPE] must be a numeric scalar > 0")
      }
      # Additional bounds checks
      if (ROPE > 0.5) stop("For alternative 'greater', ROPE must be <= 0.5")
      if ((h0 + ROPE) >= 1) stop("For alternative 'greater', h0 + ROPE must be < 1")

    } else if (alternative == "less") {
      # ROPE must be a numeric scalar < 0
      if (!is.numeric(ROPE) || length(ROPE) != 1 || !is.finite(ROPE) || ROPE >= 0) {
        stop("For alternative 'less', Argument [ROPE] must be a numeric scalar < 0")
      }
      # Additional bounds checks
      if (ROPE < -0.5) stop("For alternative 'less', ROPE must be >= -0.5")
      if ((h0 + ROPE) <= 0) stop("For alternative 'less', h0 + ROPE must be > 0")
    }

  }
  # analysis prior prior_analysis
  if (missing(prior_analysis) || !is.character(prior_analysis) || length(prior_analysis) != 1 ||
      !(prior_analysis %in% c("Moment", "beta"))) {
    stop("Argument [prior_analysis] for analysis prior must be either `beta` or `Moment`")
  }

  # prior_analysis-specific checks
  if (prior_analysis == "beta") {
    scale <- NULL
    # 'beta' requires alpha and beta to be numeric scalars > 0
    if (missing(alpha) || !is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("For prior_analysis 'beta', Argument [alpha] must be a single numeric scalar > 0")
    }
    if (missing(beta) || !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
      stop("For prior_analysis 'beta', Argument [beta] must be a single numeric scalar > 0")
    }
  } else if (prior_analysis == "Moment") {
    alpha <- beta <- NULL
    # 'Moment' requires scale to be numeric scalar > 0
    if (missing(scale)||!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("For prior_analysis 'Moment', Argument [scale] must be a numeric scalar > 0")
    }
  }


  bf10=suppressWarnings(
    if (is.null(ROPE)) {
      bin_BF(x, n, alpha, beta, h0, scale, prior_analysis, alternative)
    } else {
      bin_e_BF(x, n, alpha, beta, h0, scale, prior_analysis, alternative, ROPE)
    }
  )

  type = "One-proportion"
  p.value <- bin.pval(x,n,h0,alternative,ROPE)

  analysis_h1 <- list(
    prior = prior_analysis,
    location = h0,
    alpha=alpha,
    beta=beta,
    scale=scale
  )
  object=list(type=type,
              bf10=bf10,
              h0=h0,
              x=x,
              n=n,
              analysis_h1=analysis_h1,
              alternative=alternative,
              ROPE=ROPE,
              p.value=p.value)

  class(object) <- "BFvalue"
  return(object)
}


#' Bayes Factor for Comparing Two Proportions
#'
#' Compute the Bayes factor (BF10) for a Bayesian test of two proportions.
#'
#' @param N1 Positive numeric integer. Sample size for group 1 (must be > 0).
#' @param x1 Non-negative numeric integer. Number of successes observed in group 1 (must be \eqn{\ge 0} and \eqn{\le N_1} ).
#' @param N2 Positive numeric integer. Sample size for group 2 (must be > 0).
#' @param x2 Non-negative numeric integer. Number of successes observed in group 2 (must be \eqn{\ge 0} and \eqn{\le N_2}).
#' @param a0 Positive numeric scalar. Alpha parameter of the Beta analysis prior under the null hypothesis.
#' @param b0 Positive numeric scalar. Beta parameter of the Beta analysis prior under the null hypothesis.
#' @param a1 Positive numeric scalar. Alpha parameter of the Beta analysis prior for group 1 under the alternative hypothesis.
#' @param b1 Positive numeric scalar. Beta parameter of the Beta analysis prior for group 1 under the alternative hypothesis.
#' @param a2 Positive numeric scalar. Alpha parameter of the Beta analysis prior for group 2 under the alternative hypothesis.
#' @param b2 Positive numeric scalar. Beta parameter of the Beta analysis prior for group 2 under the alternative hypothesis.
#'
#' @return A list of class \code{BFvalue} containing:
#' \itemize{
#'   \item \code{type}: Character. Test type (always "Two-proportions").
#'   \item \code{bf10}: Numeric scalar. The computed Bayes factor in favor of the alternative hypothesis relative to the null hypothesis.
#'   \item \code{N1}: Positive integer scalar. Sample size for group 1.
#'   \item \code{x1}: Non-negative integer scalar. Observed successes for group 1.
#'   \item \code{N2}: Positive integer scalar. Sample size for group 2.
#'   \item \code{x2}: Non-negative integer scalar. Observed successes for group 2.
#'   \item \code{analysis_h0}: list with \code{a} (alpha parameter) and \code{b} (beta parameter) for the null prior.
#'   \item \code{analysis_h1_theta_1}: list with \code{a} (alpha parameter) and \code{b} (beta parameter) for group 1 prior under H1.
#'   \item \code{analysis_h1_theta_2}: list with \code{a} (alpha parameter) and \code{b} (beta parameter) for group 2 prior under H1.
#'    \item \code{OddsRatio}: Numeric scalar. Observed odds ratio.
#'    \item \code{p.value}: Numeric scalar. p-value.
#' }
#' @examples
#' BF10.props(
#' N1 = 493,
#' x1 = 155,
#' N2 = 488,
#' x2 = 150,
#' a0 = 1,
#' b0 = 1,
#' a1 = 1,
#' b1 = 1,
#' a2 = 1,
#' b2 = 1)
#' @export
BF10.props <- function(N1, x1, N2, x2, a0, b0, a1, b1, a2, b2) {


  # NULL hypothesis
  # Check a0 (alpha)
  if (!is.numeric(a0) || length(a0) != 1 || !is.finite(a0) || a0 <= 0) {
    stop("arg [a0] alpha for the Beta analysis prior under the null (\u03b8\u2080) must be a positive numeric scalar (> 0).")
  }

  # Check b0 (beta)
  if (!is.numeric(b0) || length(b0) != 1 || !is.finite(b0) || b0 <= 0) {
    stop("arg [b0] beta for the Beta analysis prior under the null (\u03b8\u2080) must be a positive numeric scalar (> 0).")
  }


  # alternative hypothesis theta1
  # Check a1 (alpha under the alternative)
  if (!is.numeric(a1) || length(a1) != 1 || !is.finite(a1) || a1 <= 0) {
    stop("arg [a1] alpha for the Beta analysis prior under the alternative (\u03b8\u2081) must be a positive numeric scalar (> 0).")
  }

  # Check b1 (beta under the alternative)
  if (!is.numeric(b1) || length(b1) != 1 || !is.finite(b1) || b1 <= 0) {
    stop("arg [b1] beta for the Beta analysis prior under the alternative (\u03b8\u2081) must be a positive numeric scalar (> 0).")
  }

  # alternative hypothesis theta2
  # Check a2 (alpha under the alternative)
  if (!is.numeric(a2) || length(a2) != 1 || !is.finite(a2) || a2 <= 0) {
    stop("arg [a2] alpha for the Beta analysis prior under the alternative (\u03b8\u2082) must be a positive numeric scalar (> 0).")
  }

  # Check b2 (beta under the alternative)
  if (!is.numeric(b2) || length(b2) != 1 || !is.finite(b2) || b2 <= 0) {
    stop("arg [b2] beta for the Beta analysis prior under the alternative (\u03b8\u2082) must be a positive numeric scalar (> 0).")
  }


  # sample sizes
  if (!is.numeric(N1) || length(N1) != 1 || N1 %% 1 != 0 ||!is.finite(N1) || N1 <= 0) {
    stop("arg [N1] sample size for group 1 must be a positive numeric scalar integer (> 0).")
  }
  if (!is.numeric(N2) || length(N2) != 1 || N2 %% 1 != 0 ||!is.finite(N2) || N2 <= 0) {
    stop("arg [N2] sample size for group 2 must be a positive numeric scalar integer (> 0).")
  }

  # observed successes
  if (!is.numeric(x1) || length(x1) != 1 || x1 %% 1 != 0 ||!is.finite(x1) || x1 < 0) {
    stop("arg [x1] for group 1 must be a non-negative numeric scalar integer (\u2265 0).")
  }
  if (!is.numeric(x2) || length(x2) != 1 || x2 %% 1 != 0 ||!is.finite(x2) || x2 < 0) {
    stop("arg [x2] for group 2 must be a non-negative numeric scalar integer (\u2265 0).")
  }

  if (x1 > N1) {
    stop("arg [x1] number of successes in group 1 cannot exceed N1 (sample size).")
  }
  if (x2 > N2) {
    stop("arg [x2] number of successes in group 2 cannot exceed N2 (sample size).")
  }

  bf10=BF10_p2(a0, b0, a1, b1, a2, b2, N1, N2, x1, x2)
  tab <- matrix(
    c(x1, N1 - x1,
      x2, N2 - x2),
    nrow = 2,
    byrow = TRUE
  )
  results <- stats::fisher.test(tab)

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
    bf10=bf10,
    N1=N1,
    x1=x1,
    N2=N2,
    x2=x2,
    analysis_h0=analysis_h0,
    analysis_h1_theta_1= analysis_h1_theta_1,
    analysis_h1_theta_2=analysis_h1_theta_2,
    OddsRatio = unname(results$estimate),
    p.value=results$p.value
  )
  class(object) <- "BFvalue"
  return(object)

}

