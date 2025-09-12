library(BayesPower)
library(mombf)
library("crch")

# functions to simulate delta
sim_delta <-function(iter, location, scale, dff, model,hypothesis){
  if(model=="Point"){
    return(rep(location, iter))
  }
  # set bounds
  p_bound <- switch(
    hypothesis,
    ">"  = c(a = "a", b = 1),
    "<"  = c(a = 0, b = "a"),
    "!=" = c(a = 0, b = 1)
  )

  # function to replace 0 with probability
  replace_zero <- function(bound, model, location , scale , dff = NULL) {
    if (bound == "a") {
      if (model == "Normal") {
        return(pnorm(0, mean = location, sd = scale))
      } else if (model == "t-distribution") {
        return(pt((0-location)/scale, df = dff, ncp = 0, lower.tail = TRUE))
      } else if (model == "NLP") {
        return(BayesPower:::pmom(0 - location, tau = scale^2))
      }
    } else {
      return(as.numeric(bound))
    }
  }

  # apply to both bounds
  p_bound <- sapply(p_bound, replace_zero, model = model, location = location,
                    scale = scale, dff = dff)

  sim_p=runif(iter,min(p_bound),max(p_bound))


  switch(model,
         "Normal"         = qnorm (sim_p, location, scale),
         "NLP"            = qmom(sim_p,  tau = scale^2)+location,
         #"t-distribution" = (qt(sim_p,df=dff)+location)/(scale)
         "t-distribution" = qtt(sim_p, location = location, scale = scale, dff, left = -Inf, right = Inf,
                                lower.tail = TRUE, log.p = FALSE)
  )

}


sim_delta_e_out <- function(iter, location, scale, dff, model, hypothesis, e = NULL) {
  if(model=="Point"){
    return(rep(location, iter))
  }
  # --- Helper: convert bounds ---
  replace_zero <- function(bound, e, model, location, scale, dff, lower = TRUE) {
    if (bound == "f") {
      if (model == "Normal") {
        return(pnorm(e, mean = location, sd = scale, lower.tail = lower))
      } else if (model == "t-distribution") {
        return(pt((e-location)/scale, df = dff, ncp = 0, lower.tail = lower))
      } else if (model == "NLP") {
        return(if (lower) BayesPower:::pmom(e - location, tau = scale^2)
               else 1 - BayesPower:::pmom(e - location, tau = scale^2))
      }
    }
    return(as.numeric(bound))
  }

  # --- Bounds setup ---
  if (length(e) == 2) {
    # two-sided interval
    upper_bounds <- c("f", 1)
    lower_bounds <- c(0, "f")

    p_upper <- sapply(upper_bounds, replace_zero, e = max(e),
                      model = model, location = location, scale = scale, dff = dff, lower = TRUE)
    p_lower <- sapply(lower_bounds, replace_zero, e = min(e),
                      model = model, location = location, scale = scale, dff = dff, lower = TRUE)

    ratio <-(1- p_upper[1]) / p_lower[2]
    n1 <- round(iter * ratio / (1 + ratio))
    n2 <- iter - n1

    sim_p <- c(
      runif(n1, min(p_upper), max(p_upper)),
      runif(n2, min(p_lower), max(p_lower))
    )

  } else {
    # one-sided interval
    bounds <- switch(hypothesis,
                     ">" = c("f", 1),
                     "<" = c(0, "f"))

    p_bound <- sapply(bounds, replace_zero,
                      e = e, # placeholder if e not given
                      model = model, location = location, scale = scale, dff = dff)

    sim_p <- runif(iter, min(p_bound), max(p_bound))
  }

  # --- Transform back ---
  switch(model,
         "Normal"         = qnorm(sim_p, mean = location, sd = scale),
         "NLP"            = qmom(sim_p , tau = scale^2)+location,
         "t-distribution" = location+ scale * qt(sim_p, df = dff)
  )
}

sim_delta_e_in <- function(iter, location, scale, dff, model, hypothesis, e = NULL) {
  if(model=="Point"){
    return(rep(location, iter))
  }
  # --- Helper: convert bounds ---
  replace_zero <- function(bound, e, model, location, scale, dff, lower = TRUE) {
    if (bound == "f") {
      if (model == "Normal") {
        return(pnorm(e, mean = location, sd = scale, lower.tail = lower))
      } else if (model == "t-distribution") {
        return(pt((e-location)/(scale), df = dff, ncp = 0, lower.tail = lower))
      } else if (model == "NLP") {
        return(if (lower) BayesPower:::pmom(e - location, tau = scale^2)
               else 1 - BayesPower:::pmom(e - location, tau = scale^2))
      }
    }
    return(as.numeric(bound))
  }

  # --- Bounds setup ---
  bounds <- switch(hypothesis,
                   "!=" = c("f", "f"),
                   ">"  = c(0.5, "f"),
                   "<"  = c("f", 0.5))

  p_bound <- sapply(bounds, replace_zero,e=e,
                    model = model, location = location, scale = scale, dff = dff, lower = TRUE)

  # then simulate
  sim_p <- runif(iter, min(p_bound), max(p_bound))

  # --- Transform back ---
  switch(model,
         "Normal"         = qnorm(sim_p, mean = location, sd = scale),
         "NLP"            = qmom(sim_p, tau = scale^2)+location,
         "t-distribution" = location+ scale * qt(sim_p, df = dff)
  )
}

# functions to simulate data and  the pro(BF>k|h_i)
sim_t <- function(iter, D, n, r = NULL,
                  model, location, scale, dff, hypothesis,
                  model_d, location_d, scale_d, dff_d,
                  de_an_prior, t1_T2, e = NULL,h0 =F) {

  # --- choose prior ---

  if(h0 && !is.null(e)){
    delta <- sim_delta_e_in(iter, location, scale, dff, model, hypothesis,e)
  } else{

    if (de_an_prior == 1) {
      delta <- if (is.null(e)) sim_delta(iter, location, scale, dff, model, hypothesis)
      else sim_delta_e_out(iter, location, scale, dff, model, hypothesis, e)
    } else {
      delta <- if (is.null(e)) sim_delta(iter, location_d, scale_d, dff_d, model_d, hypothesis)
      else sim_delta_e_out(iter, location_d, scale_d, dff_d, model_d, hypothesis, e)
    }

  }

  # --- df and scaling ---
  if (t1_T2) {
    df <- n - 1
    constant <- n
  } else {
    n2 <- n * r
    df <- n + n2 - 2
    constant <- (n * n2)/(n + n2)
  }
  # --- simulate t statistics ---
  sim_vals <- rt(iter, df, ncp = sqrt(constant) * delta)
  # --- bounds ---
  BF10_B <- if (t1_T2){
    if (is.null(e)) {BayesPower:::t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)}
    else {BayesPower:::t1e_BF10_bound(D, df, model, scale, dff, hypothesis, e)}
  } else{
    if (is.null(e)) {BayesPower:::t2_BF10_bound(D, n, r, model, location, scale, dff, hypothesis)}
    else {BayesPower:::t2e_BF10_bound(D, n, r, model, scale, dff, hypothesis, e)}
  }
  BF01_B <- if (t1_T2){
    if (is.null(e)) BayesPower:::t1_BF01_bound(D, df, model, location, scale, dff, hypothesis)
    else BayesPower:::t1e_BF01_bound(D, df, model, scale, dff, hypothesis, e)} else{
      if (is.null(e)) BayesPower:::t2_BF01_bound(D, n, r, model, location, scale, dff, hypothesis)
      else BayesPower:::t2e_BF01_bound(D, n, r, model, scale, dff, hypothesis, e)
    }

  # --- P(BF>k) ---
  PE_sim <- switch(hypothesis,
                   "!=" = mean(sim_vals > max(BF10_B) | sim_vals < min(BF10_B)),
                   ">"  = mean(sim_vals > max(BF10_B)),
                   "<"  = mean(sim_vals < max(BF10_B)))

  NE_sim <- switch(hypothesis,
                   "!=" = mean(sim_vals < max(BF01_B) & sim_vals > min(BF01_B)),
                   ">"  = mean(sim_vals < max(BF01_B)),
                   "<"  = mean(sim_vals > max(BF01_B)))

  return(c(PE_sim, NE_sim))
}

# comparing the results with the numeric method and simulation
sim_method_t<-function(iter,D,  model, location, scale, dff, hypothesis,
                       model_d, location_d, scale_d, dff_d, de_an_prior, n,r,t1_T2,e=NULL){

  results <-
    if(t1_T2){

      if (is.null(e)) {
        BayesPower:::t1_Table(
          D, target = 0.8, model, location, scale, dff, hypothesis,
          model_d, location_d, scale_d, dff_d, de_an_prior,
          N = n, mode_bf = 0, alpha = 0.05, direct = "h1"
        )
      } else {
        BayesPower:::t1e_table(
          D, target = 0.8, model, scale, dff, hypothesis, e,
          model_d, scale_d, dff_d, de_an_prior, N = n,
          mode_bf = 0, location_d = location_d,
          alpha = 0.05, direct = "h1"
        )
      }
    }else{


      if (is.null(e)) {
        N1=n
        N2 = N1 * r
        BayesPower:::t2_Table(D, r, target=.8, model, location, scale, dff, hypothesis,
                              model_d, location_d, scale_d, dff_d, de_an_prior, N1, N2,
                              mode_bf=0, alpha=.8, direct="h1")
      } else {
        N1=n
        N2 = N1 * r
        BayesPower:::t2e_table(D, r, target=.8, model, scale, dff, hypothesis, e, model_d,
                               scale_d, dff_d, de_an_prior, mode_bf=0, location_d, N1, N2, alpha=.05,
                               direct="h1")
      }

    }


  pro_h1<-sim_t(iter, D, n, r, model, location, scale, dff, hypothesis,
                model_d, location_d, scale_d, dff_d, de_an_prior, t1_T2,e,h0=F)
  if(is.null(e)) model_d ="Point"
  pro_h0<-sim_t(iter, D, n, r, model, location, scale, dff, hypothesis,
                model_d=model_d, location_d=0, scale_d, dff_d, de_an_prior=0, t1_T2,e,h0=T)
  pro_h0<-c(pro_h0[2],pro_h0[1])
  results = as.vector(unlist(results))
  out_matrix = c(results[1:4], pro_h1, pro_h0)
  return(out_matrix)
}

# making the plot and doing it for seqence of N
t_ver<-function(iter,D,  model, location, scale, dff, hypothesis,
                model_d, location_d, scale_d, dff_d, de_an_prior, r,t1_T2,e){
  if(t1_T2){
    n = seq(100,1000,100)
    total = n
  }else{
    n = seq(100,1000,100)
    n2= n*r
    total = n+n2
  }

  probs <- matrix(NA, nrow = length(n), ncol = 8)  # simpler than array
  for(i in 1:length(n)){
    probs[i,] =  as.vector(unlist(sim_method_t(iter,D,  model, location, scale, dff, hypothesis,
                                               model_d, location_d, scale_d, dff_d, de_an_prior, n[i],r,t1_T2,e)))

  }

  par(mfrow = c(1, 2))
  plot(total, probs[,1], type = "l",
       xlab = "Total sample size",
       ylab = "Probability",
       ylim = c(0,1),
       frame.plot = FALSE,
       lwd = 3,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(total, probs[,1+4], col = "gray", lwd = 2, lty = 3)
  lines(total, probs[,2], col = "black", lwd = 3)
  lines(total, probs[,2+4], col = "grey", lwd = 2, lty = 3)


  plot(total, probs[,3], type = "l",
       xlab = "Total sample size",
       ylab = "Probability",
       ylim = c(0,1),
       frame.plot = FALSE,
       lwd = 3,
       main = bquote(bold("Power curve for BF"[0][1]~">"~.(D))))
  lines(total, probs[,3+4], col = "gray", lwd = 2, lty = 3)
  lines(total, probs[,4], col = "black", lwd = 3)
  lines(total, probs[,4+4], col = "grey", lwd = 2, lty = 3)

}


# input
iter =10000 # number of iteration

D =3 # bound of evidence
model = "t-distribution"   # "Normal" or "t-distribution" or "NLP"
location = 0
scale = .707
dff = 1
hypothesis = "!="  # "!="   ">"   "<"
e= c(-.2,.2)
# e for equivalence testing
# length have to 2 when hypothesis = "!="
# one should be >0 another one should be <0
# e should <0 when hypothesis = "<"
# e should >0 when hypothesis = ">"
de_an_prior=1     # whether design and analysis priors are the same
model_d="Point" # "Normal" or "t-distribution" or "NLP" or "Point"
location_d = 2
scale_d = 1
dff_d = 1

r =1         # ratio of sample size
t1_T2 =T    # TRUE: one sample / paired t-test FALSE : two sample t-test


t_ver(iter,D,  model, location, scale, dff, hypothesis,
      model_d, location_d, scale_d, dff_d, de_an_prior, r,t1_T2,e)
# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


