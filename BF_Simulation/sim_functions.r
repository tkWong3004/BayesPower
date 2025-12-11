library(BayesPower)
library(mombf)
library(crch)
library(ExtDist)
library(mvtnorm)

#### t-test ####
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

# functions to simulate delta for equivalence test

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
              else {BayesPower:::t1e_BF10_bound(D, df, model,location, scale, dff, hypothesis, e)}
    } else{
    if (is.null(e)) {BayesPower:::t2_BF10_bound(D, n, r, model, location, scale, dff, hypothesis)}
    else {BayesPower:::t2e_BF10_bound(D, n, r, model, location,scale, dff, hypothesis, e)}
  }
  BF01_B <- if (t1_T2){
    if (is.null(e)) BayesPower:::t1_BF01_bound(D, df, model, location, scale, dff, hypothesis)
    else BayesPower:::t1e_BF01_bound(D, df, model,location, scale, dff, hypothesis, e)} else{
    if (is.null(e)) BayesPower:::t2_BF01_bound(D, n, r, model, location, scale, dff, hypothesis)
    else BayesPower:::t2e_BF01_bound(D, n, r, model, location,scale, dff, hypothesis, e)
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
      D, target = 0.8, model, location,scale, dff, hypothesis, e,
      model_d, location_d = location_d,scale_d, dff_d, de_an_prior, N = n,
      mode_bf = 0,
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
        BayesPower:::t2e_table(D, r, target=.8, model,location, scale, dff, hypothesis, e, model_d,location_d,
                               scale_d, dff_d, de_an_prior, mode_bf=0,  N1, N2, alpha=.05,
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
#### correlation ####



# functions to simulate rho
sim_rho <-function(iter, h0,k,alpha,beta,scale, model,hypothesis){
  if(model=="Point"){
    return(rep(h0, iter))
  }
  # set bounds
  rho_bound <- switch(
    hypothesis,
    ">"  = c(a = h0, b = 1),
    "<"  = c(a = -1, b = h0),
    "!=" = c(a = -1, b = 1)
  )
  p_bound<-switch(model,
                  "d_beta"   = {pBeta_ab(rho_bound,shape1 = 1/k,shape2 = 1/k,a = -1,b = 1)},
                  "NLP"   = {  BayesPower:::pmom(rho_bound-h0, tau = scale^2)},
                  "beta" = {pBeta_ab(rho_bound,shape1 = alpha,shape2 = beta,a = -1,b = 1)})


  sim_p=runif(iter,min(p_bound),max(p_bound))
  switch(model,
         "d_beta"   = {qBeta_ab(sim_p,shape1 = 1/k,shape2 = 1/k,a = -1,b = 1)},
         "NLP"   = {  qmom(sim_p, tau = scale^2)+h0},
         "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = beta,a =-1,b = 1)})

}
# functions to simulate rho for equivalence test under h1
sim_rho_out <-function(iter, h0,k,alpha,beta,scale, model,hypothesis, e = NULL){
  if(model=="Point"){
    return(rep(h0, iter))
  }

  rho_bound <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2]))

  p_bound<-switch(model,
                  "d_beta"   = {pBeta_ab(rho_bound,shape1 = 1/k,shape2 = 1/k,a = -1,b = 1)},
                  "NLP"   = {  BayesPower:::pmom(rho_bound-h0, tau = scale^2)},
                  "beta" = {pBeta_ab(rho_bound,shape1 = alpha,shape2 = beta,a =-1,b = 1)})



  if (model =="NLP"){

    sim_p <- switch(hypothesis,
                    ">" = runif(iter, max(p_bound), BayesPower:::pmom(1-h0, tau = scale^2)),
                    "<" = runif(iter, BayesPower:::pmom(-1-h0, tau = scale^2), min(p_bound)),
                    "!=" = {
                      ratio <- 1/((1 - p_bound[2]) / p_bound[1])
                      n1 <- round(iter * ratio / (1 + ratio))
                      n2 <- iter - n1
                      c(
                        runif(n1, BayesPower:::pmom(-1-h0, tau = scale^2), min(p_bound)),
                        runif(n2, max(p_bound), BayesPower:::pmom(1-h0, tau = scale^2))
                      )
                    })
  } else {
    sim_p <- switch(hypothesis,
                    ">" = runif(iter, max(p_bound), 1),
                    "<" = runif(iter, 0, min(p_bound)),
                    "!=" = {
                      ratio <- 1/((1 - p_bound[2]) / p_bound[1])
                      n1 <- round(iter * ratio / (1 + ratio))
                      n2 <- iter - n1
                      c(
                        runif(n1, 0, min(p_bound)),
                        runif(n2, max(p_bound), 1)
                      )
                    })
  }

  switch(model,
         "d_beta"   = {qBeta_ab(sim_p,shape1 = 1/k,shape2 = 1/k,a = -1,b = 1)},
         "NLP"   = {  qmom(sim_p, tau = scale^2)+h0},
         "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = beta,a = -1,b = 1)})

}

# functions to simulate rho for equivalence test under h0

sim_rho_in <-function(iter, h0,k,alpha,beta,scale, model,hypothesis, e = NULL){
  if(model=="Point"){
    return(rep(h0, iter))
  }

  rho_bound <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2]))

  p_bound<-switch(model,
                  "d_beta"   = {pBeta_ab(rho_bound,shape1 = 1/k,shape2 = 1/k,a = -1,b = 1)},
                  "NLP"   = {  BayesPower:::pmom(rho_bound-h0, tau = scale^2)},
                  "beta" = {pBeta_ab(rho_bound,shape1 = alpha,shape2 = beta,a =-1,b = 1)})

  sim_p <- runif(iter, min(p_bound), max(p_bound))

  switch(model,
         "d_beta"   = {qBeta_ab(sim_p,shape1 = 1/k,shape2 = 1/k,a = -1,b = 1)},
         "NLP"   = {  qmom(sim_p, tau = scale^2)+h0},
         "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = beta,a = -1,b = 1)})
}

# functions to simulate dataset (to get bivraite correlation of the sample) and  the pro(BF>k|h_i)
sim_cor <- function(iter, D,model,k, alpha, beta,h0,scale,dff, hypothesis ,model_d,
                    location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N,e=NULL, IN =FALSE) {
  location =h0
  # --- choose prior ---

  if(IN && !is.null(e)){
    rho <- sim_rho_in(iter, h0,k,alpha,beta,scale, model,hypothesis,e)
  } else{

    if (de_an_prior == 1) {
      rho <- if (is.null(e))sim_rho(iter, h0,k,alpha,beta,scale, model,hypothesis)
      else sim_rho_out(iter, h0,k,alpha,beta,scale, model,hypothesis,e)
    } else {
      rho <- if (is.null(e))sim_rho(iter, location_d,k_d,alpha_d,beta_d,scale_d, model_d,hypothesis)
      else sim_rho_out(iter, location_d,k_d,alpha_d,beta_d,scale_d, model_d,hypothesis,e)
    }

  }

  # --- simulate cor statistics ---

  sim_vals <- sapply(1:iter, function(i) {
    r <- rho[i]
    Sigma <- matrix(c(1, r,
                      r, 1), nrow = 2, byrow = TRUE)
    # mean vector
    mu <- rep(0, 2)
    # generate multivariate normal
    X <- rmvnorm(N, mean = mu, sigma = Sigma)
    # correlation between first two columns
    cor(X)[1, 2]
  })


  # --- bounds ---
  BF10_B <- if (is.null(e)) BayesPower:::r_BF_bound_10(D,N,k,alpha,beta,h0,hypothesis ,location,scale ,dff ,model )else
    BayesPower:::re_BF_bound_10(D,N,k,alpha,beta,h0,hypothesis ,location,scale ,dff ,model ,e)


  BF01_B <- if (is.null(e)) BayesPower:::r_BF_bound_01(D,N,k,alpha,beta,h0,hypothesis ,location,scale ,dff ,model )else
    BayesPower:::re_BF_bound_01(D,N,k,alpha,beta,h0,hypothesis ,location,scale ,dff ,model ,e)


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
sim_method_R<-function(iter, D,model,k, alpha, beta,h0,scale, hypothesis ,model_d,
                       location_d,k_d, alpha_d, beta_d,scale_d,de_an_prior,N,e=NULL){
  location = h0
  dff =dff_d=1
  results <-if (is.null(e)) {
    BayesPower:::r_table(D, target=.8, model, k, alpha, beta, h0, location=h0, scale,
                         dff, hypothesis, model_d, location_d, k_d, alpha_d, beta_d,
                         scale_d, dff_d, de_an_prior, N, mode_bf=0, FP=,05, direct)

  } else {
    BayesPower:::re_table(D, target=.8, model, k, alpha, beta, h0, location=h0, scale,
                          dff, hypothesis, model_d, location_d, k_d, alpha_d, beta_d,
                          scale_d, dff_d, de_an_prior, N, mode_bf=0, FP=,05,e, direct)
  }


  pro_h1<-sim_cor(iter, D,model,k, alpha, beta,h0,scale,dff, hypothesis ,model_d,
                  location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N,e=e, IN =FALSE)
  if(is.null(e)) model_d ="Point"
  pro_h0<-sim_cor(iter, D,model,k, alpha, beta,h0,scale,dff, hypothesis ,model_d,
                  location_d=h0,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior=0,N,e, IN =T)

  pro_h0<-c(pro_h0[2],pro_h0[1])
  results = as.vector(unlist(results))
  out_matrix = c(results[1:4], pro_h1, pro_h0)
  return(out_matrix)
}

# making the plot and doing it for seqence of N
r_ver<-function(iter, D,model,k, alpha, beta,h0,scale, hypothesis ,model_d,
                location_d,k_d, alpha_d, beta_d,scale_d,de_an_prior,e){
  n = seq(50,1000,by=(1000-50)/10)
  total = n

  probs <- matrix(NA, nrow = length(n), ncol = 8)  # simpler than array
  for(i in 1:length(n)){

    probs[i,] =  as.vector(unlist(sim_method_R(iter, D,model,k, alpha, beta,h0,scale, hypothesis ,model_d,
                                               location_d,k_d, alpha_d, beta_d,scale_d,de_an_prior,N=n[i],e)))

  }

  par(mfrow = c(1, 2))
  plot(total, probs[,1], type = "l",
       xlab = "Total sample size",
       ylab = "Probability",
       ylim = c(0,1),
       frame.plot = FALSE,
       lwd = 3,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(total, probs[,5], col = "gray", lwd = 2, lty = 3)
  lines(total, probs[,2], col = "black", lwd = 3)
  lines(total, probs[,6], col = "grey", lwd = 2, lty = 3)


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
####  f-test ####

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

    opt <- uniroot(objfun,lower=0,upper=Inf)$root

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


# functions to simulate f effect size under h0
sim_fsq_h0 <-function(iter,q, dff, rscale, f_m,model,e ){

  if (is.null(e)) return(rep(0,iter))
  if (is.null(e)){e=1}
  f_bound <- c(0,e)
  p_bound<-switch(model,
                  "effectsize"   = {p_fes(f_bound, q, dff, rscale, f_m)},
                  "Moment"   = {  p_fmoment(f_bound, q, dff,  f_m)}
  )


  sim_p=runif(iter,min(p_bound),max(p_bound))


  switch(model,
         "effectsize"   = {q_fes(sim_p, q, dff, rscale, f_m)},
         "Moment"   = {  q_fmoment(sim_p, q, dff,  f_m)}
  )

}
# functions to simulate f effect size under h1

sim_fsq_h1 <-function(iter,q, dff, rscale, f_m,model,e){

  if (model=="Point") return(rep(f_m,iter))

  if (is.null(e)){e =0}
  f_bound <- c(e,Inf)
  p_bound<-switch(model,
                  "effectsize"   = {p_fes(f_bound, q, dff, rscale, f_m)},
                  "Moment"   = {  p_fmoment(f_bound, q, dff,  f_m)}
  )

  sim_p=runif(iter,min(p_bound),max(p_bound))

  switch(model,
         "effectsize"   = {q_fes(sim_p, q, dff, rscale, f_m)},
         "Moment"   = {  q_fmoment(sim_p, q, dff,  f_m)}
  )
}

# functions to simulate F value (the data)

sim_f_value<-function(iter,D,p,k,dff,rscale,f_m,model,
                      dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, h0=F,e){
  q= k-p
  m=n-k

  if(h0) lambda=sim_fsq_h0(iter,q, dff, rscale, f_m,model,e ) else{

    if (de_an_prior ==1 ) lambda=sim_fsq_h1(iter,q, dff, rscale, f_m,model,e )else
      lambda=sim_fsq_h1(iter,q, dff_d, rscale_d, f_m_d,model_d,e )
  }
  h0=T
  if (length(unique(lambda)) == 1) {
    f_value <- sapply(1:iter, function(i) {rf(1, q, n - k, (n - k) * lambda[i]^2)})
  } else {
    f_value <- sapply(1:iter, function(i) {rf(1, q, n - k, (n - k) * lambda[i])})
  }




  BF10 = if(is.null(e)) BayesPower:::F_BF_bound_10(D, q, m, dff, rscale, f_m, model) else
    BayesPower:::Fe_BF_bound_10(D, q, m, dff, rscale, f_m, model,e)
  BF01 = if(is.null(e)) BayesPower:::F_BF_bound_01(D, q, m, dff, rscale, f_m, model) else
    BayesPower:::Fe_BF_bound_01(D, q, m, dff, rscale, f_m, model,e)

  # --- Pro(BF>k) ---
  PE_sim <-  if(BF10=="no bound is found"){0} else sum(f_value>=BF10)/iter
  NE_sim <- if(BF01=="no bound is found"){0} else sum(f_value<=BF01)/iter

  return(c(PE_sim, NE_sim))

}

# comparing the results with the numeric method and simulation
sim_method_F<-function(iter,D,p,k,dff,rscale,f_m,model,
                       dff_d,rscale_d,f_m_d,model_d,de_an_prior,n,e){

  results <-if(is.null(e)) BayesPower:::f_table(D, target=.8, p, k, dff, rscale, f_m, model, dff_d, rscale_d,
                                                f_m_d, model_d, de_an_prior, n, mode_bf=0, FP=.05, direct="h1") else
                                                  BayesPower:::fe_table(D, target=.8, p, k, dff, rscale, f_m, model, dff_d, rscale_d,
                                                                        f_m_d, model_d, de_an_prior, n, mode_bf=0, e, FP=.05, direct="h1")

  pro_h1 <- sim_f_value(iter,D,p,k,dff,rscale,f_m,model,
                        dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, h0=F,e)
  if(is.null(e)) model_d ="Point"
  pro_h0 <- sim_f_value(iter,D,p,k,dff,rscale,f_m,model,
                        dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, h0=T,e)

  pro_h0<-c(pro_h0[2],pro_h0[1])
  results = as.vector(unlist(results))
  out_matrix = c(results[1:4], pro_h1, pro_h0)
  return(out_matrix)
}

# making the plot and doing it for sequence of N
f_ver<-function(iter,D,p,k,dff,rscale,f_m,model,
                dff_d,rscale_d,f_m_d,model_d,de_an_prior,e){


  n = round(seq(50,1000,by=(1000-50)/10))
  total = n

  probs <- matrix(NA, nrow = length(n), ncol = 8)  # simpler than array
  for(i in 1:length(n)){

    probs[i,] =  as.vector(unlist(sim_method_F(iter,D,p,k,dff,rscale,f_m,model,
                                               dff_d,rscale_d,f_m_d,model_d,de_an_prior,n[i],e)))

  }


  par(mfrow = c(1, 2))
  plot(total, probs[,1], type = "l",
       xlab = "Total sample size",
       ylab = "Probability",
       ylim = c(0,1),
       frame.plot = FALSE,
       lwd = 3,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(total, probs[,5], col = "gray", lwd = 2, lty = 3)
  lines(total, probs[,2], col = "black", lwd = 3)
  lines(total, probs[,6], col = "grey", lwd = 2, lty = 3)


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


#### one proportion ####
# functions to simulate theta true proportion
sim_pro <-function(iter, h0,location, alpha,beta,scale, model,hypothesis){
  if(model=="Point"){
    return(rep(location, iter))
  }
  # set bounds
  pro_bound <- switch(
    hypothesis,
    ">"  = c(a = h0, b = 1),
    "<"  = c(a = 0, b = h0),
    "!=" = c(a = 0, b = 1)
  )

  p_bound<-switch(model,
                  "Moment"   = {  pmom(pro_bound-location, tau = scale^2)},
                  "beta" = {pBeta_ab(pro_bound,shape1 = alpha,shape2 = beta,a = 0,b = 1)})

  sim_p=runif(iter,min(p_bound),max(p_bound))

  switch(model,
         "Moment"   = {
           qmom(sim_p, tau = scale^2)+location
         },
         "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = beta,a =0,b = 1)})

}
# functions to simulate theta true proportion under h0 for equivalence test
sim_pro_out <-function(iter, h0,location, alpha,beta,scale, model,hypothesis, e = NULL){
  if(model=="Point"){
    return(rep(location, iter))
  }

  pro_bound <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2]))

  p_bound<-switch(model,
                  "Moment"   = {  BayesPower:::pmom(pro_bound-location, tau = scale^2)},
                  "beta" = {pBeta_ab(pro_bound,shape1 = alpha,shape2 = beta,a =0,b = 1)})

  if (model =="Moment"){

    sim_p <- switch(hypothesis,
                    ">" = runif(iter, max(p_bound), BayesPower:::pmom(1-location, tau = scale^2)),
                    "<" = runif(iter, BayesPower:::pmom(0-location, tau = scale^2), min(p_bound)),
                    "!=" = {
                      ratio <- 1/((1 - p_bound[2]) / p_bound[1])
                      n1 <- round(iter * ratio / (1 + ratio))
                      n2 <- iter - n1
                      c(
                        runif(n1, BayesPower:::pmom(0-location, tau = scale^2), min(p_bound)),
                        runif(n2, max(p_bound), BayesPower:::pmom(1-location, tau = scale^2))
                      )
                    })
  } else {
    sim_p <- switch(hypothesis,
                    ">" = runif(iter, max(p_bound), 1),
                    "<" = runif(iter, 0, min(p_bound)),
                    "!=" = {
                      ratio <- 1/((1 - p_bound[2]) / p_bound[1])
                      n1 <- round(iter * ratio / (1 + ratio))
                      n2 <- iter - n1
                      c(
                        runif(n1, 0, min(p_bound)),
                        runif(n2, max(p_bound), 1)
                      )
                    })
  }

  switch(model,
         "Moment"   = {  qmom(sim_p, tau = scale^2)+location},
         "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = beta,a = 0,b = 1)})

}
# functions to simulate theta true proportion under h1 for equivalence test

sim_pro_in <-function(iter,h0, location, alpha,beta,scale, model,hypothesis, e = NULL){
  if(model=="Point"){
    return(rep(location, iter))
  }

  pro_bound <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2]))

  p_bound<-switch(model,
                  "Moment"   = {  BayesPower:::pmom(pro_bound-location, tau = scale^2)},
                  "beta" = {pBeta_ab(pro_bound,shape1 = alpha,shape2 = beta,a =0,b = 1)})

  sim_p <- runif(iter, min(p_bound), max(p_bound))

  switch(model,
         "Moment"   = {  qmom(sim_p, tau = scale^2)+location},
         "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = beta,a = 0,b = 1)})
}

# functions to simulate data (number of success) and  the pro(BF>k|h_i)
sim_success <- function(iter,h0, D,alpha,beta,location,scale,model,hypothesis,
                        alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,n,e=NULL,IN=F) {
  # --- choose prior ---
  if(IN && !is.null(e)){
    pro <- sim_pro_in(iter,h0, location, alpha,beta,scale, model,hypothesis,e)
  } else{

    if (de_an_prior == 1) {
      pro <- if (is.null(e))sim_pro(iter, h0,location, alpha,beta,scale, model,hypothesis)
      else sim_pro_out(iter, h0,location, alpha,beta,scale, model,hypothesis,e)
    } else {
      pro <- if (is.null(e)) sim_pro(iter,  h0,location_d, alpha_d,beta_d,scale_d, model_d,hypothesis)
      else sim_pro_out(iter, h0,location_d,alpha_d,beta_d,scale_d, model_d,hypothesis,e)
    }

  }
  # --- simulate cor statistics ---
  sim_vals <- sapply(1:iter, function(i) {
    rbinom(1,n,pro[i])
  })

  BF10 = if(is.null(e)) BayesPower:::bin_BF(sim_vals,n,alpha,beta,location,scale,model,hypothesis) else
    BayesPower:::bin_e_BF(sim_vals,n,alpha,beta,location,scale,model,hypothesis,e)

  # --- Pro(BF>k) ---
  PE_sim <- sum(BF10>=D)/iter
  NE_sim <- sum((1/BF10)>=D)/iter
  return(c(PE_sim, NE_sim))
}

# comparing the results with the numeric method and simulation
sim_method_pro<-function(iter,h0, D,alpha,beta,location,scale,model,hypothesis,
                         alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,n,e){

  results <-if (is.null(e)) {
    BayesPower:::bin_table(D,target=.8,h0,alpha,beta,location,scale,model,hypothesis,
                           alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,n, mode_bf=0,FP=0.05,direct="h1")

  } else {
    BayesPower:::bin_e_table(D, target=.8, h0,alpha, beta, location, scale, model, hypothesis,
                             alpha_d, beta_d, location_d, scale_d, model_d, de_an_prior,
                             n, mode_bf=0, FP=.05, e, direct="h1")
  }


  pro_h1<- sim_success(iter,h0, D,alpha,beta,location,scale,model,hypothesis,
                       alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,n,e,IN=F)
  if (is.null(e)){
    model_d="Point"
    location_d=h0
  }

  pro_h0<- sim_success(iter, h0,D,alpha,beta,location,scale,model,hypothesis,
                       alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior=0,n,e,IN=T)

  pro_h0<-c(pro_h0[2],pro_h0[1])
  results = as.vector(unlist(results))
  out_matrix = c(results[1:4], pro_h1, pro_h0)
  return(out_matrix)
}

# making the plot and doing it for seqence of N
bin_ver<-function(iter,h0, D,model,  alpha, beta,location,scale, hypothesis ,model_d,
                  location_d,k_d, alpha_d, beta_d,scale_d,de_an_prior,e){
  n = round(seq(50,1000,by=(1000-50)/10))
  total = n

  probs <- matrix(NA, nrow = length(n), ncol = 8)  # simpler than array
  for(i in 1:length(n)){

    probs[i,] =  as.vector(unlist(sim_method_pro(iter,h0, D,alpha,beta,location,scale,model,hypothesis,
                                                 alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,n[i],e)))

  }


  par(mfrow = c(1, 2))
  plot(total, probs[,1], type = "l",
       xlab = "Total sample size",
       ylab = "Probability",
       ylim = c(0,1),
       frame.plot = FALSE,
       lwd = 3,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(total, probs[,5], col = "gray", lwd = 2, lty = 3)
  lines(total, probs[,2], col = "black", lwd = 3)
  lines(total, probs[,6], col = "grey", lwd = 2, lty = 3)


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


#### two proportion ####
# functions to simulate theta under h0
sim_pro_h0 <-function(iter, a0,b0){
  pro1 <- rbeta(iter, a0, b0)
  pro2 <- pro1

  # Make a 2-column matrix
  cbind(pro1, pro2)

}
# functions to simulate theta under h1
sim_pro_h1 <-function(iter, a1,b1,p1,a2,b2,p2,model1,model2){

  pro1 <- switch(model1,
                 "beta"={rbeta(iter, a1, b1)},
                 "Point"= rep(p1,iter))
  pro2 <- switch(model2,
                 "beta"={rbeta(iter, a2, b2)},
                 "Point"= rep(p2,iter))
  cbind(pro1, pro2)

}

# function to simulate data (number of success)

sim_success_ps<-function(iter, D,a0, b0, n1,n2,a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2,h0=F){
  if(h0) pro=sim_pro_h0(iter, a0,b0) else{

    if (model1==model2&model1=="same"){
      pro=sim_pro_h1(iter, a1,b1,p1,a2,b2,p2,model1,model2)
    } else {
      pro=sim_pro_h1(iter, da1,db1,dp1,da2,db2,dp2,model1,model2)
    }
  }


  x1 <- sapply(1:iter, function(i) {
    rbinom(1,n1,pro[i,1])
  })
  x2<-sapply(1:iter, function(i) {
    rbinom(1,n2,pro[i,2])
  })
  BF10 = BayesPower:::BF10_p2(a0, b0, a1, b1, a2, b2,n1,n2,x1,x2)

  # --- Pro(BF>k) ---
  PE_sim <- sum(BF10>=D)/iter
  NE_sim <- sum((1/BF10)>=D)/iter

  return(c(PE_sim, NE_sim))

}


# comparing the results with the numeric method and simulation
sim_method_pros<-function(iter,D ,a0, b0, a1, b1, a2, b2, model1,da1,db1,dp1,model2,da2,db2,dp2,n){
  n1=n2=n
  results <-
    BayesPower:::pro_table_p2(D,target=.8, a0, b0, a1, b1, a2, b2, r=1,
                              model1,da1,db1,dp1,model2,da2,db2,dp2,mode_bf=0,n1,n2,direct="h1")[[1]]


  pro_h1<- sim_success_ps(iter, D,a0, b0, n1,n2,a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2)
  pro_h0<- sim_success_ps(iter, D,a0, b0, n1,n2,a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2,h0=T)

  pro_h0<-c(pro_h0[2],pro_h0[1])
  results = as.vector(unlist(results))
  out_matrix = c(results[1:4], pro_h1, pro_h0)
  return(out_matrix)
}

# making the plot and doing it for seqence of N
p2_ver<-function(iter,D, a0, b0, a1, b1, a2, b2, model1,da1,db1,dp1,model2,da2,db2,dp2){
  if(model1==model2&model1=="same"){

    model1=model2="beta"
    da1=a1
    db1=b1
    da2=a2
    db2=b2
  }

  n = round(seq(50,1000,by=(1000-50)/10))
  total = n

  probs <- matrix(NA, nrow = length(n), ncol = 8)  # simpler than array
  for(i in 1:length(n)){

    probs[i,] =  as.vector(unlist(sim_method_pros(iter,D, a0, b0, a1, b1, a2, b2, model1,da1,db1,dp1,model2,da2,db2,dp2,n[i])))

  }


  par(mfrow = c(1, 2))
  plot(total, probs[,1], type = "l",
       xlab = "Total sample size",
       ylab = "Probability",
       ylim = c(0,1),
       frame.plot = FALSE,
       lwd = 3,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(total, probs[,5], col = "gray", lwd = 2, lty = 3)
  lines(total, probs[,2], col = "black", lwd = 3)
  lines(total, probs[,6], col = "grey", lwd = 2, lty = 3)


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

