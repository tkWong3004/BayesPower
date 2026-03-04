### Regression/ANOVA
library("BayesPower")
library("crch")
library("ExtDist")
library("mvtnorm")
#### Setting up the simulation
# number of iteration
iter <- 1000

# Specification of analysis and design prior
threshold <- 3; p <- 3; k <- 20
prior_analysis <- "effectsize"; dff <- 3 ; rscale <- .2; f_m <- sqrt(.1)
prior_design <- "Moment"; dff_d <- 5; rscale_d = .15; f_m_d <- sqrt(.01)
ROPE = NULL

# the number of additional predictors
q  <- k - p

# input validation
results <- BFpower.f.test(
  threshold=threshold,p=p,k=k,
  prior_analysis=prior_analysis,dff=dff,rscale=rscale, f_m=f_m,
  prior_design=prior_design ,dff_d=dff_d,rscale_d=rscale_d,f_m_d=f_m_d,
  N=100 , ROPE=ROPE
)

# density function of effect size prior
d_fes <- function(fsq, q, dff, rscale, f_m) {
  gamma((q + dff) / 2) / (gamma(dff / 2) * gamma(q / 2)) *
    (dff * rscale^2)^(dff / 2) * fsq^(q / 2 - 1) *
    (dff * rscale^2 + f_m^2 + fsq)^(-dff / 2 - q / 2) *
    hypergeo::genhypergeo(
      c((dff + q) / 4, (2 + dff + q) / 4),
      q / 2,
      4 * f_m^2 * fsq / (dff * rscale^2 + f_m^2 + fsq)^2
    )
}
# cumulative distribution function of effect size prior
p_fes <- function(fsq, q, dff, rscale, f_m) {
  sapply(fsq, function(x) {
    if (x == 0) {
      return(0)  # probability is 0 if fsq = 0
    } else if (is.infinite(x)) {
      return(1)  # probability is 1 if fsq = Inf
    } else {
      integrate(
        d_fes,
        lower = 0,
        upper = x,
        q = q,
        dff = dff,
        rscale = rscale,
        f_m = f_m,
        rel.tol = 1e-10,
        stop.on.error = FALSE
      )$value
    }
  })
}
# inverse of cumulative distribution function of effect size prior

q_fes <- function(p, q, dff, rscale, f_m) {
  sapply(p, function(prob) {
    if (prob == 0) return(0)
    if (prob == 1) return(Inf)

    objfun <- function(fsq) {
      p_fes(fsq, q, dff, rscale, f_m)- prob

    }

    opt <- uniroot(objfun,lower=0,upper=10000)$root

    return(opt)
  })
}
# density function of moment prior

d_fmoment<-function(fsq,q,dff,f_m){
  temp <- f_m^2 * (dff + q - 2)/2

  gamma((q + dff) / 2) / gamma(dff / 2) / gamma(q / 2) *
    2 * (dff - 2) / q / (dff-2 + q) / f_m^2 *
    fsq^(q/2) * temp^(dff/2) * (temp + fsq)^(-(dff+q)/2)
}
# cumulative distribution function of moment prior

p_fmoment <- function(fsq, q, dff, f_m) {
  sapply(fsq, function(x) {
    integrate(
      d_fmoment,
      lower = 0,
      upper = x,
      q = q,
      dff = dff,
      f_m = f_m,
      rel.tol = 1e-10,
      stop.on.error = F
    )$value
  })
}
# inverse of cumulative distribution function of moment prior

q_fmoment <- function(p, q, dff, f_m) {
  sapply(p, function(prob) {
    if (prob == 0) return(0)
    if (round(prob,10) == 1) return(Inf)

    objfun <- function(fsq) {
      p_fmoment(fsq, q, dff, f_m)- prob
    }

    # give a wider search space
    opt <- uniroot(objfun,lower=0,upper=10000)$root
    return(opt)
  })
}


# -------------------------------
# Helpers
# -------------------------------

cdf_fun <- function(x, prior, q, dff, rscale, f_m) {
  switch(prior,
         "effectsize"   = {p_fes(x, q, dff, rscale, f_m)},

         "Moment"   = {  p_fmoment(x, q, dff,  f_m)}
  )
}

q_fun <- function(p,  prior, q, dff, rscale, f_m) {
  switch(prior,
         "effectsize"   = {q_fes(p, q, dff, rscale, f_m)},
         "Moment"   = {  q_fmoment(p, q, dff,  f_m)}
  )
}

draw_trunc <- function(n, lower, upper, prior, q, dff, rscale, f_m) {

  if(prior == "Point"){return(rep(f_m^2,iter))}
  p_lower <- if(lower == 0) 0 else {cdf_fun(lower, prior, q, dff, rscale, f_m)}
  p_upper <- if(upper == Inf) 1 else { cdf_fun(upper, prior, q, dff, rscale, f_m)}

  u <- runif(n, p_lower, p_upper)
  q_fun(u, prior, q, dff, rscale, f_m)
}

# -------------------------------
# CASE 1: Point Null
# -------------------------------

if (is.null(ROPE)) {

  fsq_h0 <- rep(0, iter)
  fsq_h1 <- draw_trunc(iter, 0, Inf, prior_design, q,
                       dff_d, rscale_d, f_m_d)

}

# -------------------------------
# CASE 2: Interval Null
# -------------------------------

if (!is.null(ROPE)) {

  fsq_h0 <- draw_trunc(iter, 0, ROPE, prior_analysis, q,
                       dff, rscale, f_m)
  fsq_h1 <- draw_trunc(iter, ROPE, Inf, prior_design, q,
                       dff_d, rscale_d, f_m_d)
}


sample_size  <- seq(2*k-p+2,1000,50)
simulated    <- array(NA, dim=c(length(sample_size),4))
numeric_ones <- array(NA, dim=c(length(sample_size),4))

for (i in 1:length(sample_size)){
  n <- sample_size[i]

  q  =  k - p
  m  =  n - p

  # --- simulate f-values  ---
  fval_h0  <- sapply(1:iter, function(ii) {rf(1, q, n - k, m * fsq_h0[ii])})
  fval_h1  <- sapply(1:iter, function(ii) {rf(1, q, n - k, m * fsq_h1[ii])})



  # --- critical f-values ---
  ## We made use of the f.cri to find how many BF10>k and BF01>k
  ## We did not calculate BF for every f to speed thing up


  f.cri.h1 <- if(is.null(ROPE)) {BayesPower:::F_BF_bound_10(threshold, q, m, dff, rscale, f_m, prior_analysis )} else {
    BayesPower:::Fe_BF_bound_10(threshold, q, m, dff, rscale, f_m, prior_analysis ,ROPE)
    }

  f.cri.h0 <- if(is.null(ROPE)) {BayesPower:::F_BF_bound_01(threshold, q, m, dff, rscale, f_m, prior_analysis )} else {
    BayesPower:::Fe_BF_bound_01(threshold, q, m, dff, rscale, f_m, prior_analysis ,ROPE)}


  # --- operational characteristics ---
  TPR <- sum(fval_h1>=f.cri.h1)/iter
  FPR <- sum(fval_h0>=f.cri.h1)/iter
  # Check if f.cri.h0 is valid
  if (is.character(f.cri.h0) && f.cri.h0 == "no bound is found") {
    TNR <- 0
    FNR <- 0
  } else {
    # regular calculation
    TNR <- sum(fval_h0 <= f.cri.h0) / iter
    FNR <- sum(fval_h1 <= f.cri.h0) / iter
  }


  simulated[i,] <- c(TPR,FPR,TNR,FNR)

  numeric_results<-  BFpower.f.test(
    threshold=threshold,p=p,k=k,
    prior_analysis=prior_analysis,dff=dff,rscale=rscale, f_m=f_m,
    prior_design=prior_design ,dff_d=dff_d,rscale_d=rscale_d,f_m_d=f_m_d,
    N=n , ROPE=ROPE
  )

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


