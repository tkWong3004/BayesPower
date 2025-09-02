library(BayesPower)
library(mombf)
library(ExtDist)
library(mvtnorm)


# functions to simulate delta
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
         "beta" = {pBeta_ab(rho_bound,shape1 = alpha,shape2 = alpha,a = -1,b = 1)})


  sim_p=runif(iter,min(p_bound),max(p_bound))
  switch(model,
                  "d_beta"   = {qBeta_ab(sim_p,shape1 = 1/k,shape2 = 1/k,a = -1,b = 1)},
                  "NLP"   = {  qmom(sim_p, tau = scale^2)+h0},
                  "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = alpha,a =-1,b = 1)})

  }

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
                  "beta" = {pBeta_ab(rho_bound,shape1 = alpha,shape2 = alpha,a =-1,b = 1)})



 if (model =="NLP"){

   sim_p <- switch(hypothesis,
                   ">" = runif(iter, max(p_bound), BayesPower:::pmom(1-h0, tau = scale^2)),
                   "<" = runif(iter, BayesPower:::pmom(-1-h0, tau = scale^2), min(p_bound)),
                   "!=" = {
                     ratio <- (1 - p_bound[2]) / p_bound[1]
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
                     ratio <- (1 - p_bound[2]) / p_bound[1]
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
         "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = alpha,a = -1,b = 1)})

  }

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
                  "beta" = {pBeta_ab(rho_bound,shape1 = alpha,shape2 = alpha,a =-1,b = 1)})

  sim_p <- runif(iter, min(p_bound), max(p_bound))

  switch(model,
         "d_beta"   = {qBeta_ab(sim_p,shape1 = 1/k,shape2 = 1/k,a = -1,b = 1)},
         "NLP"   = {  qmom(sim_p, tau = scale^2)+h0},
         "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = alpha,a = -1,b = 1)})
}

# functions to simulate data and  the pro(BF>k|h_i)
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


  # --- posterior & null evidence ---
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
# input
iter = 1000 # number of iteration

D =3 # bound of evidence
model = "NLP"   # "d_beta" or "beta" or "NLP"
k =2
alpha=1
beta=1
h0=0.1
scale=2
hypothesis = "!="  # "!="   ">"   "<"
e= c(-.2,.2)
# Let it NULL when testing a point -null
# e for equivalence testing
# length have to 2 when hypothesis = "!="
# one should be >0 another one should be <0
# e should <0 when hypothesis = "<"
# e should >0 when hypothesis = ">"
de_an_prior=1     # whether design and analysis priors are the same
model_d = "Point" ##  # "d_beta" or "beta" or "NLP" "Point"
location_d=.2
k_d=1
alpha_d=1
beta_d=1
scale_d=1



r_ver(iter=1000, D,model,k, alpha, beta,h0,scale, hypothesis ,model_d,
                location_d,k_d, alpha_d, beta_d,scale_d,de_an_prior,e)
# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


