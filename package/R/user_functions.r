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
#' bp_t.test_one_sample(
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
bp_t.test_one_sample <- function(hypothesis = NULL, e = NULL, interval = NULL,
                                 D = NULL, target = NULL, alpha = NULL,
                                 model = NULL, location = NULL, scale = NULL, dff = NULL,
                                 model_d = NULL, location_d = NULL, scale_d = NULL, dff_d = NULL,
                                 de_an_prior = NULL,
                                 N = NULL, mode_bf = NULL, direct = NULL) {

  # Automatically assign 1 to any NULL argument
  args <- list(D = D, target = target, model = model, location = location, scale = scale,
               dff = dff, hypothesis = hypothesis, model_d = model_d, location_d = location_d,
               scale_d = scale_d, dff_d = dff_d, de_an_prior = de_an_prior, N = N,
               mode_bf = mode_bf, alpha = alpha, direct = direct, e = e)

  args <- lapply(args, function(x) if (is.null(x)) 1 else x)

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

  # Call appropriate table function
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
#' bp_t.test_two_sample(
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
bp_t.test_two_sample <- function(hypothesis = NULL, e = NULL, interval = NULL,
                                 D = NULL, target = NULL, alpha = NULL,
                                 model = NULL, location = NULL, scale = NULL, dff = NULL,
                                 model_d = NULL, location_d = NULL, scale_d = NULL, dff_d = NULL,
                                 de_an_prior = NULL,
                                 N1 = NULL, N2 = NULL, r = NULL, mode_bf = NULL, direct = NULL) {

  tryCatch(
    suppressWarnings({
      if (!is.null(interval) && interval == "1") {
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
#' bp_cor(
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
bp_cor <- function(hypothesis = NULL, h0 = NULL, e = NULL, interval = NULL,
                   D = NULL, target = NULL, FP = NULL,
                   model = NULL, k = NULL, alpha = NULL, beta = NULL, scale = NULL,
                   model_d = NULL, alpha_d = NULL, beta_d = NULL, location_d = NULL,
                   k_d = NULL, scale_d = NULL,
                   de_an_prior = NULL,
                   N = NULL, mode_bf = NULL, direct = NULL) {

  location <- h0
  dff <- dff_d <- 1

  tryCatch(
    suppressWarnings({
      if (!is.null(interval) && interval == "1") {
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
#' bp_f(
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
bp_f <- function(interval = NULL,
                 D = NULL, target = NULL, FP = NULL, p = NULL, k = NULL,
                 model = NULL, dff = NULL, rscale = NULL, f_m = NULL,
                 model_d = NULL, dff_d = NULL, rscale_d = NULL, f_m_d = NULL,
                 de_an_prior = NULL,
                 N = NULL, mode_bf = NULL, direct = NULL, e = NULL) {

  tryCatch({
    suppressWarnings({
      if (interval == "1") {
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
#' bp_bin(
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
bp_bin <- function(hypothesis = NULL, interval = NULL,
                   D = NULL, target = NULL, FP = NULL, location = NULL,
                   model = NULL, alpha = NULL, beta = NULL, scale = NULL,
                   model_d = NULL, alpha_d = NULL, beta_d = NULL, location_d = NULL, scale_d = NULL,
                   de_an_prior = NULL,
                   N = NULL, mode_bf = NULL, e = NULL, direct = NULL) {

  tryCatch({
    suppressWarnings({
      if (!is.null(interval) && interval == "1") {
        bin_table(D, target, alpha, beta, location, scale, model, hypothesis,
                  alpha_d, beta_d, location_d, scale_d, model_d, de_an_prior, N,
                  mode_bf, FP, direct)
      } else {
        bin_e_table(D, target, alpha, beta, location, scale, model, hypothesis,
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
#' bp_props(
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
bp_props <- function(D = NULL, target = NULL, a0 = NULL, b0 = NULL, a1 = NULL, b1 = NULL,
                     a2 = NULL, b2 = NULL, model1 = NULL,
                     a1d = NULL, b1d = NULL, dp1 = NULL, model2 = NULL,
                     a2d = NULL, b2d = NULL, dp2 = NULL,
                     mode_bf = NULL, n1 = NULL, n2 = NULL, direct = NULL) {

  r <- 1
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
