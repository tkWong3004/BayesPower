#' Sample size determination for one-sample Bayesian t-test
#'
#' Perform sample size determination or the calculation of compelling and misleading evidence.
#'
#' @param hypothesis The hypothesis being tested (e.g., two-sided \code{"!="}, right-sided \code{">"}, left-sided \code{"<"}).
#' @param e The bounds for the interval Bayes factor (used when \code{interval = 0}).
#' @param interval Integer (1 or 0). If \code{1}, Bayes factor with a point null against a composite alternative hypothesis;
#'   otherwise Bayes factor with interval null and alternative hypotheses.
#' @param D The bound of compelling evidence.
#' @param target The targeted true positive rate (if \code{direct = "h1"}) or true negative rate (if \code{direct = "h0"}).
#' @param alpha The targeted false positive rate (if \code{direct = "h1"}) or false negative rate (if \code{direct = "h0"}).
#' @param model Statistical model of the analysis prior under the alternative hypothesis: Normal distribution (\code{"Normal"}), Normal moment (\code{"NLP"}), or scaled t (\code{"t-distribution"}).
#' @param location Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Scale parameter for the analysis prior under the alternative hypothesis.
#' @param dff Degrees of freedom for the analysis prior under the alternative hypothesis (if applicable).
#' @param model_d Statistical model of the design prior under the alternative hypothesis: Normal distribution (\code{"Normal"}), Normal moment (\code{"NLP"}), or scaled t (\code{"t-distribution"}).
#' @param location_d Location parameter for the design prior under the alternative hypothesis.
#' @param scale_d Scale parameter for the design prior under the alternative hypothesis.
#' @param dff_d Degrees of freedom parameter for the design prior under the alternative hypothesis.
#' @param de_an_prior Integer (0 or 1). If 1, analysis and design priors under the alternative are the same; if 0, they are not.
#' @param N Sample size.
#' @param mode_bf Integer (0 or 1). If \code{1}, sample size determination; if \code{2}, \code{N} is needed for the calculation of probabilities of compelling and misleading evidence.
#' @param direct If \code{"h1"}, controlling true/false positive rates; if \code{"h0"}, controlling true/false negative rates.
#'
#' @examples
#' \dontrun{
#' BFpower.t.test_one_sample(
#'   hypothesis = "!=",
#'   interval = "1",
#'   D = 3,
#'   target = 0.8,
#'   alpha = 0.05,
#'   model = "t-distribution",
#'   location = 0,
#'   scale = 0.707,
#'   dff = 1,
#'   de_an_prior = 1,
#'   N = NULL,
#'   mode_bf = 1,
#'   direct = "h1"
#' )
#' }
#'
#' @export
BFpower.t.test_one_sample <- function(hypothesis = NULL, e = NULL, interval = NULL,
                                 D = NULL, target = NULL, alpha = NULL,
                                 model = NULL, location = NULL, scale = NULL, dff = NULL,
                                 model_d = NULL, location_d = NULL, scale_d = NULL, dff_d = NULL,
                                 de_an_prior = NULL,
                                 N = NULL, mode_bf = NULL, direct = NULL) {

  # Automatically assign 1 to any NULL or "NULL" argument
  args <- list(D = D, target = target, model = model, location = location, scale = scale,
               dff = dff, hypothesis = hypothesis, model_d = model_d, location_d = location_d,
               scale_d = scale_d, dff_d = dff_d, de_an_prior = de_an_prior, N = N,
               mode_bf = mode_bf, alpha = alpha, direct = direct, e = e)

  args <- lapply(args, function(x) {
    if (is.null(x) || (is.character(x) && toupper(x) == "NULL")) 1 else x
  })

  # Extract updated arguments
  D <- args$D
  target <- args$target
  model <- args$model
  location <- args$location
  scale <- args$scale
  dff <- args$dff
  hypothesis <- args$hypothesis
  model_d <- args$model_d
  location_d <- args$location_d
  scale_d <- args$scale_d
  dff_d <- args$dff_d
  de_an_prior <- args$de_an_prior
  N <- args$N
  mode_bf <- args$mode_bf
  alpha <- args$alpha
  direct <- args$direct
  e <- args$e

  # Call appropriate table function with error handling
  tryCatch(
    suppressWarnings({
      if (interval == 1) {
        t1_Table(D, target, model, location, scale, dff, hypothesis,
                 model_d, location_d, scale_d, dff_d, de_an_prior, N, mode_bf, alpha, direct)
      } else {
        t1e_table(D, target, model, scale, dff, hypothesis, e,
                  model_d, scale_d, dff_d, de_an_prior, N, mode_bf, location_d, alpha, direct)
      }
    }),
    error = function(e) {
      message("Sample size cannot be determined")
      return(NULL)
    }
  )
}

#' Sample size determination for two-sample Bayesian t-test
#'
#' Perform sample size determination or the calculation of compelling and misleading evidence
#' for a two-sample Bayesian t-test.
#'
#' @param hypothesis The hypothesis being tested (e.g., two-sided \code{"!="}, right-sided \code{">"}, left-sided \code{"<"}).
#' @param e The bounds for the interval Bayes factor (used when \code{interval = 0}).
#' @param interval Character or integer (0 or 1). If \code{"1"}, Bayes factor with a point null against a composite alternative hypothesis;
#'   otherwise Bayes factor with interval null and alternative hypotheses.
#' @param D The bound of compelling evidence.
#' @param target The targeted true positive rate (if \code{direct = "h1"}) or true negative rate (if \code{direct = "h0"}).
#' @param alpha The targeted false positive rate (if \code{direct = "h1"}) or false negative rate (if \code{direct = "h0"}).
#' @param model Statistical model of the analysis prior under the alternative hypothesis: Normal distribution (\code{"Normal"}), Normal moment (\code{"NLP"}), or scaled t (\code{"t-distribution"}).
#' @param location Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Scale parameter for the analysis prior under the alternative hypothesis.
#' @param dff Degrees of freedom for the analysis prior under the alternative hypothesis (if applicable).
#' @param model_d Statistical model of the design prior under the alternative hypothesis: Normal distribution (\code{"Normal"}), Normal moment (\code{"NLP"}), or scaled t (\code{"t-distribution"}).
#' @param location_d Location parameter for the design prior under the alternative hypothesis.
#' @param scale_d Scale parameter for the design prior under the alternative hypothesis.
#' @param dff_d Degrees of freedom parameter for the design prior under the alternative hypothesis.
#' @param de_an_prior Integer (0 or 1). If 1, analysis and design priors under the alternative are the same; if 0, they are not.
#' @param N1 Sample size of group 1.
#' @param N2 Sample size of group 2.
#' @param r Ratio of the sample size of group 2 over group 1 (N2 / N1).
#' @param mode_bf Integer (0 or 1). If \code{1}, sample size determination; if \code{2}, \code{N1} and \code{N2} are needed for the calculation of probabilities of compelling and misleading evidence.
#' @param direct If \code{"h1"}, BF10; if \code{"h0"}, BF01.
#'
#' @examples
#' \dontrun{
#' BFpower.t.test_two_sample(
#'   hypothesis = "!=",
#'   e = NULL,
#'   interval = "1",
#'   D = 3,
#'   target = 0.8,
#'   alpha = 0.05,
#'   model = "t-distribution",
#'   location = 0,
#'   scale = 0.707,
#'   dff = 1,
#'   de_an_prior = 1,
#'   r = 1,
#'   mode_bf = 1,
#'   direct = "h1"
#' )
#' }
#'
#' @export
BFpower.t.test_two_sample <- function(hypothesis = NULL, e = NULL, interval = NULL,
                                 D = NULL, target = NULL, alpha = NULL,
                                 model = NULL, location = NULL, scale = NULL, dff = NULL,
                                 model_d = NULL, location_d = NULL, scale_d = NULL, dff_d = NULL,
                                 de_an_prior = NULL,
                                 N1 = NULL, N2 = NULL, r = NULL, mode_bf = NULL, direct = NULL) {

  # Automatically assign 1 to any NULL or "NULL" argument
  args <- list(D = D, target = target, model = model, location = location, scale = scale,
               dff = dff, hypothesis = hypothesis, model_d = model_d, location_d = location_d,
               scale_d = scale_d, dff_d = dff_d, de_an_prior = de_an_prior, N1 = N1, N2 = N2,
               r = r, mode_bf = mode_bf, alpha = alpha, direct = direct, e = e, interval = interval)

  args <- lapply(args, function(x) {
    if (is.null(x) || (is.character(x) && toupper(x) == "NULL")) 1 else x
  })

  # Extract updated arguments
  D <- args$D
  target <- args$target
  model <- args$model
  location <- args$location
  scale <- args$scale
  dff <- args$dff
  hypothesis <- args$hypothesis
  model_d <- args$model_d
  location_d <- args$location_d
  scale_d <- args$scale_d
  dff_d <- args$dff_d
  de_an_prior <- args$de_an_prior
  N1 <- args$N1
  N2 <- args$N2
  r <- args$r
  mode_bf <- args$mode_bf
  alpha <- args$alpha
  direct <- args$direct
  e <- args$e
  interval <- args$interval

  tryCatch(
    suppressWarnings({
      if (!is.null(interval) && interval == 1) {
        t2_Table(D, r, target, model, location, scale, dff, hypothesis,
                 model_d, location_d, scale_d, dff_d, de_an_prior, N1, N2, mode_bf, alpha, direct)
      } else {
        t2e_table(D, r, target, model, scale, dff, hypothesis, e,
                  model_d, scale_d, dff_d, de_an_prior, mode_bf, location, N1, N2, alpha, direct)
      }
    }),
    error = function(e) {
      message("Sample size cannot be determined")
      return(NULL)
    }
  )
}


#' Sample size determination for Bayesian correlation test
#'
#' Perform sample size determination or the calculation of compelling and misleading evidence
#' for a Bayesian correlation test.
#'
#' @param hypothesis The hypothesis being tested (e.g., two-sided \code{"!="}, right-sided \code{">"}, left-sided \code{"<"}).
#' @param h0 Null value of the correlation.
#' @param e The bounds for the interval Bayes factor (used when \code{interval = 0}).
#' @param interval Character or integer (0 or 1). If \code{"1"}, Bayes factor with a point null against a composite alternative hypothesis;
#'   otherwise Bayes factor with interval null and alternative hypotheses.
#' @param D The bound of compelling evidence.
#' @param target The targeted true positive rate (if \code{direct = "h1"}) or true negative rate (if \code{direct = "h0"}).
#' @param FP The targeted false positive rate (if \code{direct = "h1"}) or false negative rate (if \code{direct = "h0"}).
#' @param model Statistical model of the analysis prior under the alternative hypothesis: default beta (\code{"d_beta"}), beta (\code{"beta"}), or normal moment (\code{"NLP"})
#' @param k Parameter for the analysis default beta prior under the alternative hypothesis.
#' @param alpha Parameter for the analysis beta prior under the alternative hypothesis.
#' @param beta Parameter for the analysis beta prior under the alternative hypothesis.
#' @param scale Scale parameter for the analysis normal moment prior under the alternative hypothesis.
#' @param model_d Statistical model of the design prior under the alternative hypothesis:default beta (\code{"d_beta"}), beta (\code{"beta"}), normal moment (\code{"NLP"} , or point \code{"Point"})
#' @param alpha_d Parameter for the design beta prior under the alternative hypothesis.
#' @param beta_d Parameter for the design beta prior under the alternative hypothesis.
#' @param location_d Location parameter for the design point prior under the alternative hypothesis.
#' @param k_d Parameter for the design default beta prior under the alternative hypothesis.
#' @param scale_d Scale parameter for the design normal moment prior under the alternative hypothesis.
#' @param de_an_prior Integer (0 or 1). If 1, analysis and design priors under the alternative are the same; if 0, they are not.
#' @param N Sample size.
#' @param mode_bf Integer (0 or 1). If \code{1}, sample size determination; if \code{2}, \code{N} is needed for the calculation of probabilities of compelling and misleading evidence.
#' @param direct If \code{"h1"}, BF10; if \code{"h0"}, BF01.
#'
#' @examples
#' \dontrun{
#' BFpower.cor(
#'   hypothesis = "!=",
#'   h0 = 0,
#'   e = NULL,
#'   interval = "1",
#'   D = 3,
#'   target = 0.8,
#'   FP = 0.05,
#'   model = "d_beta",
#'   k = 1,
#'   de_an_prior = 1,
#'   mode_bf = 1,
#'   direct = "h1"
#' )
#' }
#'
#' @export
BFpower.cor <- function(hypothesis = NULL, h0 = NULL, e = NULL, interval = NULL,
                   D = NULL, target = NULL, FP = NULL,
                   model = NULL, k = NULL, alpha = NULL, beta = NULL, scale = NULL,
                   model_d = NULL, alpha_d = NULL, beta_d = NULL, location_d = NULL,
                   k_d = NULL, scale_d = NULL,
                   de_an_prior = NULL,
                   N = NULL, mode_bf = NULL, direct = NULL) {

  location <- h0
  dff <- dff_d <- 1

  # Handle NULL and "NULL"
  args <- list(
    hypothesis = hypothesis, h0 = h0, e = e, interval = interval, D = D, target = target,
    FP = FP, model = model, k = k, alpha = alpha, beta = beta, scale = scale,
    model_d = model_d, alpha_d = alpha_d, beta_d = beta_d, location_d = location_d,
    k_d = k_d, scale_d = scale_d, de_an_prior = de_an_prior, N = N,
    mode_bf = mode_bf, direct = direct
  )

  args <- lapply(args, function(x) {
    if (is.null(x) || (is.character(x) && toupper(x) == "NULL")) 1 else x
  })

  # Extract updated arguments
  hypothesis <- args$hypothesis
  h0 <- args$h0
  e <- args$e
  interval <- args$interval
  D <- args$D
  target <- args$target
  FP <- args$FP
  model <- args$model
  k <- args$k
  alpha <- args$alpha
  beta <- args$beta
  scale <- args$scale
  model_d <- args$model_d
  alpha_d <- args$alpha_d
  beta_d <- args$beta_d
  location_d <- args$location_d
  k_d <- args$k_d
  scale_d <- args$scale_d
  de_an_prior <- args$de_an_prior
  N <- args$N
  mode_bf <- args$mode_bf
  direct <- args$direct

  tryCatch(
    suppressWarnings({
      if (!is.null(interval) && interval == 1) {
        r_table(D, target, model, k, alpha, beta, h0, location, scale, dff,
                hypothesis, model_d, location_d, k_d, alpha_d, beta_d, scale_d,
                dff_d, de_an_prior, N, mode_bf, FP, direct)
      } else {
        re_table(D, target, model, k, alpha, beta, h0, location, scale, dff,
                 hypothesis, model_d, location_d, k_d, alpha_d, beta_d, scale_d,
                 dff_d, de_an_prior, N, mode_bf, FP, e, direct)
      }
    }),
    error = function(e) {
      message("Sample size cannot be determined")
      return(NULL)
    }
  )
}

#' Sample size determination for Bayesian F-test
#'
#' Perform sample size determination or the calculation of compelling and misleading evidence
#' for a Bayesian F-test.
#'
#' @param interval Character or integer (0 or 1). If \code{"1"}, Bayes factor with a point null against a composite alternative hypothesis;
#'   otherwise Bayes factor with interval null and alternative hypotheses.
#' @param e The bounds for the interval Bayes factor (used when \code{interval = 0}).
#' @param D The bound of compelling evidence.
#' @param target The targeted true positive rate (if \code{direct = "h1"}) or true negative rate (if \code{direct = "h0"}).
#' @param FP The targeted false positive rate (if \code{direct = "h1"}) or false negative rate (if \code{direct = "h0"}).
#' @param p Number of predictors in the reduced model.
#' @param k Number of predictors in the full model.
#' @param model Statistical model of the analysis prior under the alternative hypothesis: effect size prior  (\code{"effectsize"}) or Moment prior (\code{"Moment"})
#' @param dff Degrees of freedom for the analysis prior under the alternative hypothesis.(must be >3 if moment prior is used)
#' @param rscale Scaling parameter for the analysis effect size prior.
#' @param f_m Cohen's f effect size parameter for the analysis prior.
#' @param model_d Statistical model of the design prior under the alternative hypothesis:: effect size prior  (\code{"effectsize"}), Moment prior (\code{"Moment"}), or Point prior (\code{"Point"})
#' @param dff_d Degrees of freedom for the design prior under the alternative hypothesis. (must be >3 if moment prior is used)
#' @param rscale_d Scaling parameter for the design effect size prior.
#' @param f_m_d Cohen's f effect size parameter for the design prior or the point design prior.
#' @param de_an_prior Integer (0 or 1). If 1, analysis and design priors under the alternative are the same; if 0, they are not.
#' @param N Sample size.
#' @param mode_bf Integer (0 or 1). If \code{1}, sample size determination; if \code{2}, \code{N} is needed for the calculation of probabilities of compelling and misleading evidence.
#' @param direct If \code{"h1"}, BF10; if \code{"h0"}, BF01.
#'
#' @examples
#' \dontrun{
#' BFpower.f(
#'   interval = "1",
#'   e = NULL,
#'   D = 3,
#'   target = 0.8,
#'   FP = 0.05,
#'   p = 2,
#'   k = 3,
#'   model = "f-distribution",
#'   dff = 1,
#'   rscale = 0.5,
#'   f_m = 0.15,
#'   de_an_prior = 1,
#'   mode_bf = 1,
#'   direct = "h1"
#' )
#' }
#'
#' @export
BFpower.f <- function(interval = NULL,
                 D = NULL, target = NULL, FP = NULL, p = NULL, k = NULL,
                 model = NULL, dff = NULL, rscale = NULL, f_m = NULL,
                 model_d = NULL, dff_d = NULL, rscale_d = NULL, f_m_d = NULL,
                 de_an_prior = NULL,
                 N = NULL, mode_bf = NULL, direct = NULL, e = NULL) {

  # Handle NULL and "NULL"
  args <- list(
    interval = interval, D = D, target = target, FP = FP, p = p, k = k,
    model = model, dff = dff, rscale = rscale, f_m = f_m,
    model_d = model_d, dff_d = dff_d, rscale_d = rscale_d, f_m_d = f_m_d,
    de_an_prior = de_an_prior, N = N, mode_bf = mode_bf,
    direct = direct, e = e
  )

  args <- lapply(args, function(x) {
    if (is.null(x) || (is.character(x) && toupper(x) == "NULL")) {
      1
    } else if (is.character(x) && x == "1") {
      1
    } else {
      x
    }
  })
  # Extract updated arguments
  interval   <- args$interval
  D          <- args$D
  target     <- args$target
  FP         <- args$FP
  p          <- args$p
  k          <- args$k
  model      <- args$model
  dff        <- args$dff
  rscale     <- args$rscale
  f_m        <- args$f_m
  model_d    <- args$model_d
  dff_d      <- args$dff_d
  rscale_d   <- args$rscale_d
  f_m_d      <- args$f_m_d
  de_an_prior<- args$de_an_prior
  N          <- args$N
  mode_bf    <- args$mode_bf
  direct     <- args$direct
  e          <- args$e

  tryCatch({
    suppressWarnings({
      if (!is.null(interval) && interval == 1) {
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
    message("Sample size cannot be determined")
    return(invisible(NULL))
  })
}


#' Sample size determination for Bayesian one-proportion test
#'
#' Perform sample size determination or the calculation of compelling and misleading evidence
#' for a Bayesian test of a single proportion.
#'
#' @param hypothesis The hypothesis being tested (e.g., two-sided \code{"!="}, right-sided \code{">"}, left-sided \code{"<"}).
#' @param interval Character or integer (0 or 1). If \code{"1"}, Bayes factor with a point null against a composite alternative hypothesis;
#'   otherwise Bayes factor with interval null and alternative hypotheses.
#' @param D The bound of compelling evidence.
#' @param target The targeted true positive rate (if \code{direct = "h1"}) or true negative rate (if \code{direct = "h0"}).
#' @param FP The targeted false positive rate (if \code{direct = "h1"}) or false negative rate (if \code{direct = "h0"}).
#' @param location Null proportion value.
#' @param model Statistical model of the analysis prior under the alternative hypothesis: beta prior  (\code{"beta"}) or Moment prior (\code{"Moment"})
#' @param alpha Parameter for the analysis prior under the alternative hypothesis.
#' @param beta Parameter for the analysis prior under the alternative hypothesis.
#' @param scale Scale parameter for the analysis prior under the alternative hypothesis.
#' @param model_d Statistical model of the design prior under the alternative hypothesis:beta prior  (\code{"beta"}) , Moment prior (\code{"Moment"}), or Point prior (\code{"Point"})
#' @param alpha_d Parameter for the design prior under the alternative hypothesis.
#' @param beta_d Parameter for the design prior under the alternative hypothesis.
#' @param location_d The proportion value for the design point prior.
#' @param scale_d Scale parameter for the design prior under the alternative hypothesis.
#' @param de_an_prior Integer (0 or 1). If 1, analysis and design priors under the alternative are the same; if 0, they are not.
#' @param N Sample size.
#' @param mode_bf Integer (0 or 1). If \code{1}, sample size determination; if \code{2}, \code{N} is needed for the calculation of probabilities of compelling and misleading evidence.
#' @param e The bounds for the interval Bayes factor (used when \code{interval = 0}).
#' @param direct If \code{"h1"}, BF10; if \code{"h0"}, BF01.
#'
#' @examples
#' \dontrun{
#' BFpower.bin(
#'   hypothesis = "!=",
#'   interval = "1",
#'   D = 3,
#'   target = 0.8,
#'   FP = 0.05,
#'   location = 0.5,
#'   model = "beta",
#'   alpha = 1,
#'   beta = 1,
#'   de_an_prior = 1,
#'   mode_bf = 1,
#'   direct = "h1"
#' )
#' }
#'
#' @export
BFpower.bin <- function(hypothesis = NULL, interval = NULL,
                   D = NULL, target = NULL, FP = NULL, location = NULL,
                   model = NULL, alpha = NULL, beta = NULL, scale = NULL,
                   model_d = NULL, alpha_d = NULL, beta_d = NULL, location_d = NULL, scale_d = NULL,
                   de_an_prior = NULL,
                   N = NULL, mode_bf = NULL, e = NULL, direct = NULL, h0 = NULL) {

  # Handle NULL and "NULL"
  args <- list(
    hypothesis = hypothesis, interval = interval, D = D, target = target, FP = FP,
    location = location, model = model, alpha = alpha, beta = beta, scale = scale,
    model_d = model_d, alpha_d = alpha_d, beta_d = beta_d,
    location_d = location_d, scale_d = scale_d,
    de_an_prior = de_an_prior, N = N, mode_bf = mode_bf, e = e, direct = direct,
    h0 = h0
  )

  args <- lapply(args, function(x) {
    if (is.null(x) || (is.character(x) && toupper(x) == "NULL")) {
      1
    } else if (is.character(x) && x == "1") {
      1
    } else {
      x
    }
  })

  # Extract updated arguments
  hypothesis  <- args$hypothesis
  interval    <- args$interval
  D           <- args$D
  target      <- args$target
  FP          <- args$FP
  location    <- args$location
  model       <- args$model
  alpha       <- args$alpha
  beta        <- args$beta
  scale       <- args$scale
  model_d     <- args$model_d
  alpha_d     <- args$alpha_d
  beta_d      <- args$beta_d
  location_d  <- args$location_d
  scale_d     <- args$scale_d
  de_an_prior <- args$de_an_prior
  N           <- args$N
  mode_bf     <- args$mode_bf
  e           <- args$e
  direct      <- args$direct
  h0          <- args$h0

  tryCatch({
    suppressWarnings({
      if (!is.null(interval) && interval == 1) {
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
    return(invisible(NULL))
  })
}

#' Sample size determination for Bayesian test of two proportions
#'
#' Perform sample size determination or the calculation of compelling and misleading evidence
#' for a Bayesian comparison of two proportions.
#'
#' @param D The bound of compelling evidence.
#' @param target The targeted true positive rate (if \code{direct = "h1"}) or true negative rate (if \code{direct = "h0"}).
#' @param a0 Alpha parameter of the beta distribution under the null .
#' @param b0 Beta parameter of the beta distribution  under the null.
#' @param model1 Statistical model of the design prior for group 1: beta (\code{"beta"}), Point prior (\code{"Point"}, or same as analysis prior \code{"same"})
#' @param a1 Alpha parameter of the analysis beta prior  distribution for group 1 under the alternative hypothesis.
#' @param b1 Beta parameter of the analysis beta prior  distribution for group 1 under the alternative hypothesis.
#' @param a2 Alpha parameter of the analysis beta prior distribution for group 2 under the alternative hypothesis.
#' @param b2 Beta parameter  of the analysis beta prior  distribution for group 2 under the alternative hypothesis.
#' @param model2 Statistical model of the design prior for group 1: beta (\code{"beta"}), or Point prior (\code{"Point"}, or same as analysis prior \code{"same"})
#' @param a1d Alpha parameter for the design prior of group 1.
#' @param b1d Beta parameter for the design prior of group 1.
#' @param dp1 True proportion for group 1 in the design prior.
#' @param a2d Alpha parameter for the design prior of group 2.
#' @param b2d Beta parameter for the design prior of group 2.
#' @param dp2 True proportion for group 2 in the design prior.
#' @param mode_bf Integer (0 or 1). If \code{1}, sample size determination; if \code{2}, \code{n1} and \code{n2} are used for the calculation of probabilities of compelling and misleading evidence.
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param direct If \code{"h1"}, BF10; if \code{"h0"}, BF01.
#'
#' @examples
#' \dontrun{
#' BFpower.props(
#'   D = 3,
#'   target = 0.8,
#'   a0 = 1,
#'   b0 = 1,
#'   model1 = "same",
#'   a1 = 1,
#'   b1 = 1,
#'   a2 = 1,
#'   b2 = 1,
#'   model2 = "same",
#'   mode_bf = 1,
#'   direct = "h1"
#' )
#' }
#'
#' @export
BFpower.props <- function(D = NULL, target = NULL, a0 = NULL, b0 = NULL, a1 = NULL, b1 = NULL,
                     a2 = NULL, b2 = NULL, model1 = NULL,
                     a1d = NULL, b1d = NULL, dp1 = NULL, model2 = NULL,
                     a2d = NULL, b2d = NULL, dp2 = NULL,
                     mode_bf = NULL, n1 = NULL, n2 = NULL, direct = NULL) {

  r <- 1

  # Handle NULL and "NULL"
  args <- list(
    D = D, target = target, a0 = a0, b0 = b0, a1 = a1, b1 = b1,
    a2 = a2, b2 = b2, model1 = model1, a1d = a1d, b1d = b1d, dp1 = dp1,
    model2 = model2, a2d = a2d, b2d = b2d, dp2 = dp2,
    mode_bf = mode_bf, n1 = n1, n2 = n2, direct = direct
  )

  args <- lapply(args, function(x) {
    if (is.null(x) || (is.character(x) && toupper(x) == "NULL")) {
      1
    } else if (is.character(x) && x == "1") {
      1
    } else {
      x
    }
  })

  # Extract updated arguments
  D        <- args$D
  target   <- args$target
  a0       <- args$a0
  b0       <- args$b0
  a1       <- args$a1
  b1       <- args$b1
  a2       <- args$a2
  b2       <- args$b2
  model1   <- args$model1
  a1d      <- args$a1d
  b1d      <- args$b1d
  dp1      <- args$dp1
  model2   <- args$model2
  a2d      <- args$a2d
  b2d      <- args$b2d
  dp2      <- args$dp2
  mode_bf  <- args$mode_bf
  n1       <- args$n1
  n2       <- args$n2
  direct   <- args$direct

  tryCatch({
    suppressWarnings({
      pro_table_p2(D, target, a0, b0, a1, b1, a2, b2, r, model1,
                   a1d, b1d, dp1, model2, a2d, b2d, dp2, mode_bf, n1, n2, direct)
    })
  }, error = function(e) {
    message("Sample size cannot be determined")
    return(invisible(NULL))
  })
}




#' Bayes factor for one-sample Bayesian t-test
#'
#' Calculate the Bayes factor (BF10) for a one-sample Bayesian t-test, either against a point null or an interval null hypothesis.
#'
#' @param tval Observed t-value from the one-sample t-test.
#' @param df Degrees of freedom for the t-test.
#' @param model Statistical model of the analysis prior under the alternative hypothesis: Normal distribution (\code{"Normal"}), Normal moment (\code{"NLP"}), or scaled t (\code{"t-distribution"}).
#' @param location Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Scale parameter for the analysis prior under the alternative hypothesis.
#' @param dff Degrees of freedom for the analysis prior under the alternative hypothesis (if applicable).
#' @param hypothesis The hypothesis being tested: two-sided (\code{"!="}), right-sided (\code{">"}), or left-sided (\code{"<"}).
#' @param e Optional numeric vector specifying bounds for an interval null; used if interval BF is calculated.
#'
#' @return The Bayes factor (BF10) for the one-sample t-test.
#'
#' @examples
#' \dontrun{
#' BF10.t.test.one_sample(
#'   tval = 2.31,
#'   df = 29,
#'   model = "t-distribution",
#'   location = 0,
#'   scale = 0.707,
#'   dff = 1,
#'   hypothesis = "!="
#' )
#' }
#' @export
BF10.t.test.one_sample <- function(tval, df, model, location, scale, dff, hypothesis, e = NULL) {
  suppressWarnings(
    if (is.null(e)) {
      t1_BF10(tval, df, model, location, scale, dff, hypothesis)
    } else {
      t1e_BF10(tval, df,model,scale,dff , hypothesis,e )
    }
  )
}


#' Bayes factor for two-sample Bayesian t-test
#'
#' Calculate the Bayes factor (BF10) for a two-sample Bayesian t-test, either against a point null or an interval null hypothesis.
#'
#' @param tval Observed t-value from the two-sample t-test.
#' @param N1 Sample size of group 1.
#' @param N2 Sample size of group 2.
#' @param model Statistical model of the analysis prior under the alternative hypothesis: Normal distribution (\code{"Normal"}), Normal moment (\code{"NLP"}), or scaled t (\code{"t-distribution"}).
#' @param location Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Scale parameter for the analysis prior under the alternative hypothesis.
#' @param dff Degrees of freedom for the analysis prior under the alternative hypothesis (if applicable).
#' @param hypothesis The hypothesis being tested: two-sided (\code{"!="}), right-sided (\code{">"}), or left-sided (\code{"<"}).
#' @param e Optional numeric vector specifying bounds for an interval null; used if interval BF is calculated.
#'
#' @return The Bayes factor (BF10) for the two-sample t-test.
#'
#' @examples
#' \dontrun{
#' BF10.t.test.two_sample(
#'   tval = 2.1,
#'   N1 = 30,
#'   N2 = 30,
#'   model = "t-distribution",
#'   location = 0,
#'   scale = 0.707,
#'   dff = 1,
#'   hypothesis = "!="
#' )
#' }
#' @export
BF10.t.test.two_sample <- function(tval, N1, N2, model, location, scale, dff, hypothesis, e = NULL) {
  n1 <- N1
  n2 <- N2
  r <- n2 / n1
  suppressWarnings(
    if (is.null(e)) {
      t2_BF10(tval, n1, r, model, location, scale, dff, hypothesis)
    } else {
      t2e_BF10(tval, n1, r, model, scale, dff, hypothesis, e)
    }
  )
}

#' Bayes factor for a Bayesian correlation test
#'
#' Calculate the Bayes factor (BF10) for a correlation, either against a point null or an interval null hypothesis.
#'
#' @param r Observed correlation coefficient.
#' @param n Sample size.
#' @param k Parameter for the analysis default beta prior under the alternative hypothesis.
#' @param alpha Parameter for the analysis beta prior under the alternative hypothesis.
#' @param beta Parameter for the analysis beta prior under the alternative hypothesis.
#' @param h0 Null value of the correlation.
#' @param hypothesis The hypothesis being tested: two-sided (\code{"!="}), right-sided (\code{">"}), or left-sided (\code{"<"}).
#' @param location Location parameter for the analysis prior under the alternative hypothesis.
#' @param scale Scale parameter for the analysis normal moment prior under the alternative hypothesis.
#' @param dff Degrees of freedom for the analysis prior under the alternative hypothesis (if applicable).
#' @param model Statistical model of the analysis prior under the alternative hypothesis: default beta (\code{"d_beta"}), beta (\code{"beta"}), or normal moment (\code{"NLP"}).
#' @param e Optional numeric vector specifying bounds for an interval null; used if interval BF is calculated.
#'
#' @return The Bayes factor (BF10) for the correlation test.
#'
#' @examples
#' \dontrun{
#' BF10.cor(
#'   r = 0.3,
#'   n = 50,
#'   k = 1,
#'   alpha = 0.05,
#'   beta = 0.2,
#'   h0 = 0,
#'   hypothesis = "!=",
#'   location = 0,
#'   scale = 1,
#'   dff = 49,
#'   model = "d_beta"
#' )
#' }
#' @export
BF10.cor <- function(r, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model, e = NULL) {
  suppressWarnings(
    if (is.null(e)) {
      r_BF10(r, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model)
    } else {
      re_BF10(r, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model, e)
    }
  )
}


#' Bayes factor for a Bayesian F-test
#'
#' Calculate the Bayes factor (BF10) for an F-test, either against a point null or an interval null hypothesis.
#'
#' @param fval Observed F-value from the F-test.
#' @param df1 Degrees of freedom for the numerator of the F-test.
#' @param df2 Degrees of freedom for the denominator of the F-test.
#' @param dff Degrees of freedom for the analysis prior under the alternative hypothesis (if applicable).
#' @param rscale Scaling parameter for the analysis effect size prior.
#' @param f_m Cohen's f effect size parameter for the analysis prior.
#' @param model Statistical model of the analysis prior under the alternative hypothesis: effect size prior (\code{"effectsize"}) or Moment prior (\code{"Moment"}).
#' @param e Optional numeric vector specifying bounds for an interval null; used if interval BF is calculated.
#'
#' @return The Bayes factor (BF10) for the F-test.
#'
#' @examples
#' \dontrun{
#' BF10.f.test(
#'   fval = 4.5,
#'   df1 = 2,
#'   df2 = 12,
#'   dff = 12,
#'   rscale = 0.707,
#'   f_m = "medium",
#'   model = "effectsize"
#' )
#' }
#' @export
BF10.f.test <- function(fval, df1, df2, dff, rscale, f_m, model, e = NULL) {

  q <- df1
  m <- df1 + df2

  suppressWarnings(
    if (is.null(e)) {
      F_BF(fval, q, m, dff, rscale, f_m, model)
    } else {
      Fe_BF(fval, q, m, dff, rscale, f_m, model, e)
    }
  )
}


#' Bayes factor for a Bayesian one-proportion test
#'
#' Calculate the Bayes factor (BF10) for a test of a single proportion, either against a point null or an interval null hypothesis.
#'
#' @param x Observed number of successes.
#' @param n Sample size.
#' @param alpha Parameter for the analysis beta prior under the alternative hypothesis.
#' @param beta Parameter for the analysis beta prior under the alternative hypothesis.
#' @param location Null proportion value.
#' @param scale Scale parameter for the analysis prior (if applicable, e.g., for Moment prior).
#' @param model Statistical model of the analysis prior under the alternative hypothesis: beta prior (\code{"beta"}) or Moment prior (\code{"Moment"}).
#' @param hypothesis The hypothesis being tested: two-sided (\code{"!="}), right-sided (\code{">"}), or left-sided (\code{"<"}).
#' @param e Optional numeric vector specifying bounds for an interval null; used if interval BF is calculated.
#'
#' @return The Bayes factor (BF10) for the one-proportion test.
#'
#' @examples
#' \dontrun{
#' BF10.bin.test(
#'   x = 12,
#'   n = 50,
#'   alpha = 2,
#'   beta = 3,
#'   location = 0.5,
#'   scale = 1,
#'   model = "beta",
#'   hypothesis = "!="
#' )
#' }
#' @export
BF10.bin.test <- function(x, n, alpha, beta, location, scale, model, hypothesis, e = NULL) {
  suppressWarnings(
    if (is.null(e)) {
      bin_BF(x, n, alpha, beta, location, scale, model, hypothesis)
    } else {
      bin_e_BF(x, n, alpha, beta, location, scale, model, hypothesis, e)
    }
  )
}


#' Bayes factor for a Bayesian test of two proportions
#'
#' Calculate the Bayes factor (BF10) for comparing two proportions using a Bayesian framework.
#'
#' @param a0 Alpha parameter of the beta distribution under the null hypothesis.
#' @param b0 Beta parameter of the beta distribution under the null hypothesis.
#' @param a1 Alpha parameter of the analysis beta prior distribution for group 1 under the alternative hypothesis.
#' @param b1 Beta parameter of the analysis beta prior distribution for group 1 under the alternative hypothesis.
#' @param a2 Alpha parameter of the analysis beta prior distribution for group 2 under the alternative hypothesis.
#' @param b2 Beta parameter of the analysis beta prior distribution for group 2 under the alternative hypothesis.
#' @param n1 Sample size for group 1.
#' @param n2 Sample size for group 2.
#' @param y1 Observed number of successes for group 1.
#' @param y2 Observed number of successes for group 2.
#'
#' @return The Bayes factor (BF10) for comparing two proportions.
#'
#' @examples
#' \dontrun{
#' BF10.props(
#'   a0 = 2, b0 = 3,
#'   a1 = 2, b1 = 3,
#'   a2 = 2, b2 = 3,
#'   n1 = 50, n2 = 60,
#'   y1 = 25, y2 = 30
#' )
#' }
#' @export
BF10.props <- function(a0, b0, a1, b1, a2, b2, n1, n2, x1, x2) {
  BF10_p2(a0, b0, a1, b1, a2, b2, n1, n2, x1, x2)
}

