### Independent samples t-tests (equal variance)
library("BayesPower")
library("crch")
library("ExtDist")
library("mvtnorm")
# number of iteration
iter <- 5000

# Specification of analysis and design prior
alternative <- "greater"; ROPE <- .2
prior_analysis <- "t-distribution"; location <- 0; scale <- .707; dff <- 3
prior_design <- "Moment"; location_d <- 0; scale_d <- .707; dff_d <- .3
threshold <- 3

# input validation
results <-BFpower.ttest.TwoSample(
  alternative = alternative,ROPE = ROPE,
  prior_analysis = prior_analysis,location = location, scale = scale,dff= dff,
  prior_design = prior_design,location_d = location_d,scale_d = scale_d,
  dff_d=dff_d,threshold = threshold,N1=50,N2=50
)
# -------------------------------
# Helpers
# -------------------------------

cdf_fun <- function(x, prior, location, scale, dff) {
  switch(prior,
         "Normal" = pnorm(x, location, scale),
         "t-distribution" = pt((x - location)/scale, df = dff),
         "Moment" = BayesPower:::pmom(x - location, tau = scale^2)
  )
}

q_fun <- function(p, prior, location, scale, dff) {
  switch(prior,
         "Normal" = qnorm(p, location, scale),
         "t-distribution" = location + scale * qt(p, df = dff),
         "Moment" = qmom(p, tau = scale^2) + location
  )
}

draw_trunc <- function(n, lower, upper, prior, location, scale, dff) {
  if (prior == "Point") return(rep(location, iter))
  p_lower <- if (is.infinite(lower)) 0 else cdf_fun(lower, prior, location, scale, dff)
  p_upper <- if (is.infinite(upper)) 1 else cdf_fun(upper, prior, location, scale, dff)

  u <- runif(n, p_lower, p_upper)
  q_fun(u, prior, location, scale, dff)
}

# -------------------------------
# CASE 1: Point Null
# -------------------------------

if (is.null(ROPE)) {

  delta_h0 <- rep(0, iter)

  bounds <- switch(alternative,
                   "greater"   = c(0, Inf),
                   "less"      = c(-Inf, 0),
                   "two.sided" = c(-Inf, Inf)
  )

  delta_h1 <- draw_trunc(iter,
                         lower = bounds[1],
                         upper = bounds[2],
                         prior = prior_design,
                         location = location_d,
                         scale = scale_d,
                         dff = dff_d)
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

    delta_h0 <- draw_trunc(iter,
                           lower = min(ROPE),
                           upper = max(ROPE),
                           prior = prior_analysis,
                           location = location,
                           scale = scale,
                           dff = dff)

  } else {

    if (alternative == "greater") {
      delta_h0 <- draw_trunc(iter,
                             lower = 0,
                             upper = ROPE,
                             prior = prior_analysis,
                             location = location,
                             scale = scale,
                             dff = dff)

    } else {
      delta_h0 <- draw_trunc(iter,
                             lower = ROPE,
                             upper = 0,
                             prior = prior_analysis,
                             location = location,
                             scale = scale,
                             dff = dff)
    }
  }

  # -------- H1 --------


  if (length(ROPE) == 2) {

    # proportional tail sampling
    lower <- min(ROPE)
    upper <- max(ROPE)

    p_lower <- cdf_fun(lower, prior_design, location_d, scale_d, dff_d)
    p_upper <- cdf_fun(upper, prior_design, location_d, scale_d, dff_d)

    mass_lower <- p_lower
    mass_upper <- 1 - p_upper

    n_upper <- round(iter * mass_upper / (mass_lower + mass_upper))
    n_lower <- iter - n_upper

    delta_h1 <- c(
      draw_trunc(n_upper, upper, Inf, prior_design, location_d, scale_d, dff_d),
      draw_trunc(n_lower, -Inf, lower, prior_design, location_d, scale_d, dff_d)
    )

  } else {

    if (alternative == "greater") {
      delta_h1 <- draw_trunc(iter,
                             lower = ROPE,
                             upper = Inf,
                             prior = prior_design,
                             location = location_d,
                             scale = scale_d,
                             dff = dff_d)

    } else {
      delta_h1 <- draw_trunc(iter,
                             lower = -Inf,
                             upper = ROPE,
                             prior = prior_design,
                             location = location_d,
                             scale = scale_d,
                             dff = dff_d)
    }
  }
}


sample_size  <- seq(10,1000,50)
r            <- 1    # ratio of sample size n2/n1
simulated    <- array(NA, dim=c(length(sample_size),4))
numeric_ones <- array(NA, dim=c(length(sample_size),4))
total   <- array(NA, dim = length(sample_size))
for (i in 1:length(sample_size)){
  n1 = sample_size[i]
  n2 = n1*r
  total[i] = n1 + n2

  df  <- n1 + n2 - 2
  ess <- (n1 * n2)/(n1 + n2)

  # --- simulate t statistics ---
  tval_h0 <- rt(iter, df, ncp = sqrt(ess) * delta_h0)
  tval_h1 <- rt(iter, df, ncp = sqrt(ess) * delta_h1)

  # --- critical t-values ---
  ## We made use of the t.cri to find how many BF10>k and BF01>k
  ## We did not calculate BF for every t-value to speed thing up
  suppressWarnings(
    if(is.null(ROPE)){
      t.cri.h1 <- BayesPower:::t2_BF10_bound(threshold, n1, r,  prior_analysis,
                                             location, scale, dff, alternative)
      t.cri.h0 <- BayesPower:::t2_BF01_bound(threshold, n1,r, prior_analysis,
                                             location, scale, dff, alternative)
    } else{

      t.cri.h1 <- BayesPower:::t2e_BF10_bound(threshold, n1,r, prior_analysis,
                                              location, scale, dff,
                                              alternative,ROPE)
      t.cri.h0 <- BayesPower:::t2e_BF01_bound(threshold, n1,r, prior_analysis,
                                              location, scale, dff,
                                              alternative,ROPE)

    }
  )
  # --- operational characteristics ---
  TPR <- switch(alternative,
                "two.sided" = mean(tval_h1 > max(t.cri.h1) | tval_h1 <
                                     min(t.cri.h1)),
                "greater"  = mean(tval_h1 > max(t.cri.h1)),
                "less"     = mean(tval_h1 < max(t.cri.h1)))
  FPR <- switch(alternative,
                "two.sided" = mean(tval_h0 > max(t.cri.h1) | tval_h0 <
                                     min(t.cri.h1)),
                "greater"  = mean(tval_h0 > max(t.cri.h1)),
                "less"     = mean(tval_h0 < max(t.cri.h1)))

  if (any(t.cri.h0 == "bound cannot be found")) {

    TNR <- 0
    FNR <- 0

  } else {

    TNR <- switch(alternative,
                  "two.sided" = mean(tval_h0 < max(t.cri.h0) &
                                       tval_h0 > min(t.cri.h0)),
                  "greater"   = mean(tval_h0 < max(t.cri.h0)),
                  "less"      = mean(tval_h0 > min(t.cri.h0))
    )

    FNR <- switch(alternative,
                  "two.sided" = mean(tval_h1 < max(t.cri.h0) &
                                       tval_h1 > min(t.cri.h0)),
                  "greater"   = mean(tval_h1 < max(t.cri.h0)),
                  "less"      = mean(tval_h1 > min(t.cri.h0))
    )
  }


  simulated[i,] <- c(TPR,FPR,TNR,FNR)

  numeric_results<- BFpower.ttest.TwoSample(
    alternative = alternative,ROPE = ROPE,
    prior_analysis = prior_analysis,location = location,
    scale = scale,dff= dff,
    prior_design = prior_design,location_d = location_d,
    scale_d = scale_d,dff_d=dff_d,
    threshold = threshold, N1 =n1,N2=n2)

  numeric_ones[i,] <-  c(numeric_results$results[[1]],
                         numeric_results$results[[4]],
                         numeric_results$results[[3]],
                         numeric_results$results[[2]])
}

par(mfrow = c(1, 2))
plot(total, numeric_ones[,1], type = "l",
     xlab = "Total sample size",
     ylab = "Probability",
     ylim = c(0,1),
     frame.plot = FALSE,
     lwd = 3,
     main = bquote(bold("Power curve for BF"[10]~">"~.(threshold))))
lines(total, simulated[,1], col = "gray", lwd = 2, lty = 3)
lines(total, numeric_ones[,2], col = "black", lwd = 3)
lines(total, simulated[,2], col = "grey", lwd = 2, lty = 3)

plot(total, numeric_ones[,3], type = "l",
     xlab = "Total sample size",
     ylab = "Probability",
     ylim = c(0,1),
     frame.plot = FALSE,
     lwd = 3,
     main = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold))))
lines(total, simulated[,3], col = "gray", lwd = 2, lty = 3)
lines(total, numeric_ones[,4], col = "black", lwd = 3)
lines(total, simulated[,4], col = "grey", lwd = 2, lty = 3)

