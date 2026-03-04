### Bivariate Correlation
library("BayesPower")
library("crch")
library("ExtDist")
library("mvtnorm")
#### Setting up the simulation


# number of iteration
iter <- 500

# Specification of analysis and design prior
alternative <- "two.sided"; h0 <- 0; ROPE <- NULL
prior_analysis <- "beta";k <- 1; alpha <- 1; beta <- 1; scale <- .5
prior_design <- "Moment"; alpha_d <- 2; beta_d <- 2;
location_d <- 0 ;  k_d <- 2; scale_d <- .3
threshold <- 3

# input validation
results <- BFpower.cor(
  alternative=alternative, h0=h0, ROPE=ROPE,
  prior_analysis=prior_analysis,k=k,alpha=alpha,beta=beta,scale=scale,
  prior_design=prior_design ,alpha_d=alpha_d,beta_d=beta_d,
  location_d=location_d, k_d=k_d, scale_d =scale_d,
  N=50,  threshold=threshold)


# -------------------------------
# Helpers
# -------------------------------

cdf_fun <- function(x, prior, k,alpha,beta,location,scale) {
  switch(prior,
         "d_beta"   = {pBeta_ab(x,shape1 = 1/k,shape2 = 1/k,
                                a = -1,b = 1)},
         "Moment"   = {  BayesPower:::pmom(x-location,
                                           tau = scale^2)},
         "beta" = {pBeta_ab(x,shape1 = alpha,shape2 = beta,
                            a = -1,b = 1)}
  )
}

q_fun <- function(p,  prior, k,alpha,beta,location,scale) {
  switch(prior,
         "d_beta"   = {  qBeta_ab(p,shape1 = 1/k,shape2 = 1/k,a = -1,b = 1)},
         "Moment"   = {  qmom(p, tau = scale^2)+h0},
         "beta"     = {  qBeta_ab(p,shape1 = alpha,shape2 = beta,a =-1,b = 1)}
  )
}

draw_trunc <- function(n, lower, upper, prior, k,alpha,beta,location,scale) {
  if (prior == "Point") return(rep(location, iter))

  p_lower <- cdf_fun(lower, prior,k,alpha,beta,location, scale)
  p_upper <- cdf_fun(upper, prior, k,alpha,beta,location, scale)

  u <- runif(n, p_lower, p_upper)
  q_fun(u, prior, k,alpha,beta,location,scale)
}

# -------------------------------
# CASE 1: Point Null
# -------------------------------

if (is.null(ROPE)) {

  rho_h0 <- rep(h0, iter)

  bounds <- switch(alternative,
                   "greater"   = c(h0, 1),
                   "less"      = c(-1, h0),
                   "two.sided" = c(-1, 1)
  )

  rho_h1 <- draw_trunc(iter,
                       lower = bounds[1],
                       upper = bounds[2],
                       prior = prior_design,
                       k  = k_d,
                       alpha = alpha_d,
                       beta = beta_d,
                       location = location_d,
                       scale = scale_d)

}

# -------------------------------
# CASE 2: Interval Null
# -------------------------------

if (!is.null(ROPE)) {

  # -------- H0 --------
  # Note: for interval Bayes factor, the null has the
  #       same analysis prior as the H1. However,
  #       the delta under H1 made use of the design prior.

  if (length(ROPE) == 2) {


    rho_h0   <- draw_trunc(iter,
                           lower = min(ROPE),
                           upper = max(ROPE),
                           prior = prior_analysis,
                           k  = k,
                           alpha = alpha,
                           beta = beta,
                           location = h0,
                           scale = scale)

  } else {

    if (alternative == "greater") {

      rho_h0   <- draw_trunc(iter,
                             lower = h0,
                             upper = h0+ROPE,
                             prior = prior_analysis,
                             k  = k,
                             alpha = alpha,
                             beta = beta,
                             location = h0,
                             scale = scale_d)


    } else {
      rho_h0   <- draw_trunc(iter,
                             lower = h0+ROPE,
                             upper = h0,
                             prior = prior_analysis,
                             k  = k,
                             alpha = alpha,
                             beta = beta,
                             location = location,
                             scale = scale)
    }
  }

  # -------- H1 --------


  if (length(ROPE) == 2) {

    # proportional tail sampling
    lower <- h0+min(ROPE)
    upper <- h0+max(ROPE)

    rho_h0   <- draw_trunc(iter,
                           lower = lower,
                           upper = upper,
                           prior = prior_design,
                           k  = k_d,
                           alpha = alpha_d,
                           beta = beta_d,
                           location = location_d,
                           scale = scale_d)



    p_lower <- cdf_fun(lower, prior_design, k_d,
                       alpha_d,beta_d,location = location_d,scale_d)
    p_upper <- cdf_fun(upper, prior_design, k_d,
                       alpha_d,beta_d,location = location_d,scale_d)

    mass_lower <- p_lower
    mass_upper <- 1 - p_upper

    n_upper <- round(iter * mass_upper / (mass_lower + mass_upper))
    n_lower <- iter - n_upper


    rho_h1 <- c(draw_trunc(n_upper,
                           lower = upper,
                           upper = 1,
                           prior = prior_design,
                           k  = k_d,
                           alpha = alpha_d,
                           beta = beta_d,
                           location = location_d,
                           scale = scale_d),
                draw_trunc(n_lower,
                           lower = -1,
                           upper =lower,
                           prior = prior_design,
                           k  = k_d,
                           alpha = alpha_d,
                           beta = beta_d,
                           location = location_d,
                           scale = scale_d)
    )


  } else {

    if (alternative == "greater") {
      rho_h1 <- draw_trunc(iter,
                           lower = h0+ROPE,
                           upper = 1,
                           prior = prior_design,
                           k  = k_d,
                           alpha = alpha_d,
                           beta = beta_d,
                           location = location_d,
                           scale = scale_d)

    } else {
      rho_h1 <- draw_trunc(iter,
                           lower = -1,
                           upper = h0+ROPE,
                           prior = prior_design,
                           k  = k_d,
                           alpha = alpha_d,
                           beta = beta_d,
                           location = location_d,
                           scale = scale_d)
    }
  }
}



sample_size  <- seq(10,1000,50)
simulated    <- array(NA, dim=c(length(sample_size),4))
numeric_ones <- array(NA, dim=c(length(sample_size),4))

# the location parameter of the analysis prior is always equal to h0
location <- h0

for (i in 1:length(sample_size)){
  N <- sample_size[i]
  # --- simulate r  ---
  rval_h0 <- sapply(1:iter, function(i) {
    r <- rho_h0[i]
    Sigma <- matrix(c(1, r,
                      r, 1), nrow = 2, byrow = TRUE)
    # mean vector
    mu <- rep(0, 2)
    # generate multivariate normal
    X <- rmvnorm(N, mean = mu, sigma = Sigma)
    # correlation between first two columns
    cor(X)[1, 2]
  })
  rval_h1 <- sapply(1:iter, function(i) {
    r <- rho_h1[i]
    Sigma <- matrix(c(1, r,
                      r, 1), nrow = 2, byrow = TRUE)
    # mean vector
    mu <- rep(0, 2)
    # generate multivariate normal
    X <- rmvnorm(N, mean = mu, sigma = Sigma)
    # correlation between first two columns
    cor(X)[1, 2]
  })


  # --- critical correlation ---
  ## We made use of the r.cri to find how many BF10>k and BF01>k
  ## We did not calculate BF for every r to speed thing up


  r.cri.h1 <- if (is.null(ROPE)){ BayesPower:::r_BF_bound_10(threshold,N,k,alpha,beta,h0,alternative ,location,scale ,1 , prior_analysis ) }else{
    BayesPower:::re_BF_bound_10(threshold,N,k,alpha,beta,h0,alternative ,location,scale ,dff ,prior_analysis  ,ROPE)}

  r.cri.h0 <- if (is.null(ROPE)){ BayesPower:::r_BF_bound_01(threshold,N,k,alpha,beta,h0,alternative ,location,scale ,1 , prior_analysis ) } else {
    BayesPower:::re_BF_bound_01(threshold,N,k,alpha,beta,h0,alternative ,location,scale ,dff ,prior_analysis  ,ROPE)}


  # --- operational characteristics ---
  TPR <- switch(alternative,
                "two.sided" = mean(rval_h1 > max(r.cri.h1) | rval_h1 <
                                     min(r.cri.h1)),
                "greater"  = mean(rval_h1 > max(r.cri.h1)),
                "less"     = mean(rval_h1 < max(r.cri.h1)))
  FPR <- switch(alternative,
                "two.sided" = mean(rval_h0 > max(r.cri.h1) | rval_h0 <
                                     min(r.cri.h1)),
                "greater"  = mean(rval_h0 > max(r.cri.h1)),
                "less"     = mean(rval_h0 < max(r.cri.h1)))
  if (any(r.cri.h0 == "bound cannot be found")){
    TNR <- FNR <- 0
  } else {
  TNR <- switch(alternative,
                "two.sided" = mean(rval_h0 < max(r.cri.h0) & rval_h0 >
                                     min(r.cri.h0)),
                "greater"  = mean(rval_h0 < max(r.cri.h0)),
                "less"  = mean(rval_h0 > max(r.cri.h0)))
  FNR <- switch(alternative,
                "two.sided" = mean(rval_h1 < max(r.cri.h0) & rval_h1 >
                                     min(r.cri.h0)),
                "greater"  = mean(rval_h1 < max(r.cri.h0)),
                "less"  = mean(rval_h1 > max(r.cri.h0)))
  }

  simulated[i,] <- c(TPR,FPR,TNR,FNR)

  numeric_results<- BFpower.cor(
    alternative=alternative, h0=h0, ROPE=ROPE,
    prior_analysis=prior_analysis,k=k,alpha=alpha,beta=beta,scale=scale,
    prior_design=prior_design ,alpha_d=alpha_d,beta_d=beta_d,
    location_d=location_d, k_d=k_d, scale_d =scale_d,
    N=N,  threshold=threshold)

  numeric_ones[i,] <-  c(numeric_results$results[[1]],
                         numeric_results$results[[4]],
                         numeric_results$results[[3]],
                         numeric_results$results[[2]])
}

par(mfrow = c(1, 2))
plot(sample_size, numeric_ones[,1], type = "l",
     xlab = "Total sample size",
     ylab = "Probability",
     ylim = c(0,1),
     frame.plot = FALSE,
     lwd = 3,
     main = bquote(bold("Power curve for BF"[10]~">"~.(threshold))))
lines(sample_size, simulated[,1], col = "gray", lwd = 2, lty = 3)
lines(sample_size, numeric_ones[,2], col = "black", lwd = 3)
lines(sample_size, simulated[,2], col = "grey", lwd = 2, lty = 3)


plot(sample_size, numeric_ones[,3], type = "l",
     xlab = "Total sample size",
     ylab = "Probability",
     ylim = c(0,1),
     frame.plot = FALSE,
     lwd = 3,
     main = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold))))
lines(sample_size, simulated[,3], col = "gray", lwd = 2, lty = 3)
lines(sample_size, numeric_ones[,4], col = "black", lwd = 3)
lines(sample_size, simulated[,4], col = "grey", lwd = 2, lty = 3)
