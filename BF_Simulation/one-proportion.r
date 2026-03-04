### One proportion
library("BayesPower")
library("crch")
library("ExtDist")
library("mvtnorm")
#### Setting up the simulation

# number of iteration
iter <- 10000

# Specification of analysis and design prior
alternative <- "two.sided"; threshold <- 3; h0 <- .5
prior_analysis <- "beta"
alpha <- 1; beta <- 1; scale <- .25
prior_design <- "beta"
alpha_d <- 1; beta_d <- 2; location_d <- .4; scale_d <- .1
ROPE <- NULL

# input validation
results <- BFpower.bin(
  alternative=alternative,threshold=threshold,h0=h0,
  prior_analysis=prior_analysis,
  alpha=alpha,beta=beta,scale=scale,
  prior_design=prior_design,
  alpha_d=alpha_d,beta_d=beta_d,location_d=location_d,scale_d=scale_d,
  N =50,ROPE=ROPE)


# -------------------------------
# Helpers
# -------------------------------

cdf_fun <- function(x, prior, location, scale, alpha, beta) {
  switch(prior,
         "Moment" = BayesPower:::pmom(x - location, tau = scale^2),
         "beta"   = pBeta_ab(x, shape1 = alpha, shape2 = beta, a = 0, b = 1)
  )
}

q_fun <- function(p, prior, location, scale, alpha, beta) {
  switch(prior,
         "Moment" = qmom(p, tau = scale^2) + location,
         "beta"   = qBeta_ab(p, shape1 = alpha, shape2 = beta, a = 0, b = 1)
  )
}

draw_trunc <- function(n, lower, upper, prior, location, scale, alpha, beta) {

  if (prior == "Point")
    return(rep(location, n))

  p_lower <- cdf_fun(lower, prior, location, scale, alpha, beta)
  p_upper <- cdf_fun(upper, prior, location, scale, alpha, beta)

  u <- runif(n, p_lower, p_upper)

  q_fun(u, prior, location, scale, alpha, beta)
}
# -------------------------------
# CASE 1: Point Null
# -------------------------------

if (is.null(ROPE)) {

  theta_h0 <- rep(h0, iter)

  bounds <- switch(alternative,
                   "greater"   = c(h0, 1),
                   "less"      = c(0, h0),
                   "two.sided" = c(0, 1)
  )

  theta_h1 <- draw_trunc(iter, min(bounds), max(bounds),
                         prior_design, location_d, scale_d, alpha_d, beta_d)
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

    theta_h0 <- draw_trunc(iter, h0+min(ROPE), h0+max(ROPE),
                           prior_analysis, h0, scale, alpha, beta)

  } else {

    if (alternative == "greater") {
      theta_h0 <- draw_trunc(iter, h0, h0+ROPE,
                             prior_analysis, h0, scale, alpha, beta)

    }
    if (alternative == "less"){
      theta_h0 <- draw_trunc(iter, h0+ROPE, h0,
                             prior_analysis, h0, scale, alpha, beta)
    }
  }

  # -------- H1 --------


  if (length(ROPE) == 2) {

    # proportional tail sampling
    lower <- h0+min(ROPE)
    upper <- h0+max(ROPE)

    p_lower <- cdf_fun(lower, prior_design,
                       location_d, scale_d, alpha_d, beta_d)
    p_upper <- cdf_fun(upper, prior_design,
                       location_d, scale_d, alpha_d, beta_d)

    mass_lower <- p_lower
    mass_upper <- 1 - p_upper

    n_upper <- round(iter * mass_upper / (mass_lower + mass_upper))
    n_lower <- iter - n_upper

    theta_h1 <- c(
      draw_trunc(n_upper, upper, 1,
                 prior_design,
                 location_d, scale_d, alpha_d, beta_d),
      draw_trunc(n_lower, 0, lower,
                 prior_design,
                 location_d, scale_d, alpha_d, beta_d)
    )

  } else {

    if (alternative == "greater") {
      theta_h1 <- draw_trunc(iter,h0+ROPE, 1,
                             prior_design,
                             location_d, scale_d, alpha_d, beta_d)

    } else {
      theta_h1 <- draw_trunc(iter,0, h0+ROPE,
                             prior_design,
                             location_d, scale_d, alpha_d, beta_d)
    }
  }
}

sample_size  <- seq(10,1000,50)
simulated    <- array(NA, dim=c(length(sample_size),4))
numeric_ones <- array(NA, dim=c(length(sample_size),4))
for (i in 1:length(sample_size)){
  n = sample_size[i]

  # --- simulate number of success ---
  suc_h0 <-  sapply(1:iter, function(ii) {
    rbinom(1,n,theta_h0[ii])
  })
  suc_h1 <-  sapply(1:iter, function(ii) {
    rbinom(1,n,theta_h1[ii])
  })

  # --- operational characteristics ---
  BF10_h0 <-  sapply(1:iter, function(ii) {
    #BF10.bin.test(suc_h0[ii],n, alpha,beta,h0,scale,
    #                      prior_analysis,alternative,ROPE)$bf10

    # we use the BF function from the backend instead of
    # the exported one to speed thing up
    if(is.null(ROPE)) {
      BayesPower:::bin_BF(suc_h0[ii],n,alpha,beta,h0,
                          scale,prior_analysis,alternative)
    } else{

      BayesPower:::bin_e_BF(suc_h0[ii],n,alpha,beta,h0,
                            scale,prior_analysis,alternative,ROPE)
    }
  })
  BF10_h1 <-  sapply(1:iter, function(ii) {
    #BF10.bin.test(suc_h1[ii],n, alpha,beta,h0,scale,
    #                     prior_analysis,alternative,ROPE)$bf10
    if(is.null(ROPE)) {
      BayesPower:::bin_BF(suc_h1[ii],n,alpha,beta,h0,
                          scale,prior_analysis,alternative)
    } else{

      BayesPower:::bin_e_BF(suc_h1[ii],n,alpha,beta,h0,
                            scale,prior_analysis,alternative,ROPE)
    }
  })



  # --- operational characteristics ---
  TPR <- sum(BF10_h1>=threshold)/iter
  FPR <- sum(BF10_h0>=threshold)/iter

  TNR <- sum((1/BF10_h0)>=threshold)/iter
  FNR <- sum((1/BF10_h1)>=threshold)/iter


  simulated[i,] <- c(TPR,FPR,TNR,FNR)

  numeric_results<- BFpower.bin(
    alternative=alternative,threshold=threshold,h0=h0,
    prior_analysis=prior_analysis,
    alpha=alpha,beta=beta,scale=scale,
    prior_design=prior_design,
    alpha_d=alpha_d,beta_d=beta_d,location_d=location_d,scale_d=scale_d,
    N =n,ROPE=ROPE)

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



