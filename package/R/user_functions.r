#' @export
bp_t.test_one_sample <- function(interval,
                                         D, target, model, location, scale, dff, hypothesis,
                                         model_d, location_d, scale_d, dff_d, de_an_prior,
                                         N, mode_bf, alpha, direct, e) {

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
  if (interval == 1) {
    t1_Table(D, target, model, location, scale, dff, hypothesis,
             model_d, location_d, scale_d, dff_d, de_an_prior, N, mode_bf, alpha, direct)
  } else {
    t1e_table(D, target, model, scale, dff, hypothesis, e,
              model_d, scale_d, dff_d, de_an_prior, N, mode_bf, location_d, alpha, direct)
  }
}


#' @export
bp_t.test_two_sample <- function(D, r, target, model, location, scale, dff, hypothesis,
                                         model_d, location_d, scale_d, dff_d, de_an_prior,
                                         N1, N2, mode_bf, alpha, direct, e , interval) {

  if (interval == "1") {
    t2_Table(D, r, target, model, location, scale, dff, hypothesis,
             model_d, location_d, scale_d, dff_d, de_an_prior, N1, N2, mode_bf, alpha, direct)
  } else {
    t2e_table(D, r, target, model, scale, dff, hypothesis, e,
              model_d, scale_d, dff_d, de_an_prior, mode_bf, location, N1, N2, alpha, direct)
  }
}
#' @export
bp_cor <- function(interval,
                           D, target, model, k, alpha, beta, h0,  scale, dff,
                           hypothesis, model_d, location_d, k_d, alpha_d, beta_d, scale_d, dff_d,
                           de_an_prior, N, mode_bf, FP, e, direct) {
  location = h0
  if (interval == "1") {
    r_table(D, target, model, k, alpha, beta, h0, location, scale, dff,
            hypothesis, model_d, location_d, k_d, alpha_d, beta_d, scale_d,
            dff_d, de_an_prior, N, mode_bf, FP, direct)

  } else {
    re_table(D, target, model, k, alpha, beta, h0, location, scale, dff,
             hypothesis, model_d, location_d, k_d, alpha_d, beta_d, scale_d,
             dff_d, de_an_prior, N, mode_bf, FP, e, direct)
  }
}
#' @export
bp_f <- function(interval,
                         D, target, p, k, dff, rscale, f_m, model,
                         dff_d, rscale_d, f_m_d, model_d, de_an_prior, n,
                         mode_bf, FP, e, direct) {

  if (interval == "1") {
    f_table(D, target, p, k, dff, rscale, f_m, model,
            dff_d, rscale_d, f_m_d, model_d, de_an_prior, n,
            mode_bf, FP, direct)

  } else {
    fe_table(D, target, p, k, dff, rscale, f_m, model,
             dff_d, rscale_d, f_m_d, model_d, de_an_prior, n,
             mode_bf, e, FP, direct)
  }
}

#' @export
bp_bin <- function(interval,
                           D, target, alpha, beta, location, scale, model, hypothesis,
                           alpha_d, beta_d, location_d, scale_d, model_d, de_an_prior, N,
                           mode_bf, FP, e = NULL, direct) {

  if (interval == "1") {
    bin_table(D, target, alpha, beta, location, scale, model, hypothesis,
              alpha_d, beta_d, location_d, scale_d, model_d, de_an_prior, N,
              mode_bf, FP, direct)

  } else {
    bin_e_table(D, target, alpha, beta, location, scale, model, hypothesis,
                alpha_d, beta_d, location_d, scale_d, model_d, de_an_prior, N,
                mode_bf, FP, e, direct)
  }
}

#' @export
bp_props <- function(D,target, a0, b0, a1, b1, a2, b2,model1,
                             da1,db1,dp1,model2,da2,db2,dp2,mode_bf,n1,n2,direct) {
r=1
  pro_table_p2(D,target, a0, b0, a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2,mode_bf,n1,n2,direct)

  }

