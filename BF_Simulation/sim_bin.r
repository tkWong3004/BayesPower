library(BayesPower)
library(mombf)
library(ExtDist)
library(mvtnorm)

# functions to simulate delta
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

sim_pro_in <-function(iter,h0, location, alpha,beta,scale, model,hypothesis, e = NULL){
  if(model=="Point"){
    return(rep(location, iter))
  }

  pro_bound <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = location),
                      "!=" = c(a = h0+e[1], b = h0+e[2]))

  p_bound<-switch(model,
                  "Moment"   = {  BayesPower:::pmom(pro_bound-location, tau = scale^2)},
                  "beta" = {pBeta_ab(pro_bound,shape1 = alpha,shape2 = beta,a =0,b = 1)})

  sim_p <- runif(iter, min(p_bound), max(p_bound))

  switch(model,
         "Moment"   = {  qmom(sim_p, tau = scale^2)+location},
         "beta" = {qBeta_ab(sim_p,shape1 = alpha,shape2 = beta,a = 0,b = 1)})
}

# functions to simulate data and  the pro(BF>k|h_i)
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
# input
iter = 1000 # number of iteration

D =3 # bound of evidence
model = "beta"   #  "beta" or "Moment"

alpha=3
beta=3
h0=location=0.5
scale=2
hypothesis = "!="  # "!="   ">"   "<"
e= NULL
# Let it NULL when testing a point -null
# e for equivalence testing
# length have to be 2 when hypothesis = "!="
# one should be >0 another one should be <0
# e should <0 when hypothesis = "<"
# e should >0 when hypothesis = ">"
de_an_prior=0     # whether design and analysis priors are the same
model_d = "Point" ##  # "beta" or "Moment" or "Point"
location_d=.38
alpha_d=2
beta_d=2
scale_d=1


bin_ver(iter,h0, D,model,  alpha, beta,location,scale, hypothesis ,model_d,
                  location_d,k_d, alpha_d, beta_d,scale_d,de_an_prior,e)
# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


