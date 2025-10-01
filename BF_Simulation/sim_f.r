library(BayesPower)


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

d_fmoment<-function(fsq,q,dff,f_m){
  temp <- f_m^2 * (dff + q - 2)/2

gamma((q + dff) / 2) / gamma(dff / 2) / gamma(q / 2) *
  2 * (dff - 2) / q / (dff-2 + q) / f_m^2 *
  fsq^(q/2) * temp^(dff/2) * (temp + fsq)^(-(dff+q)/2)
}

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

# functions to simulate delta
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

# making the plot and doing it for seqence of N
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
### !!! computationally demanding simulation !!! ###
# input
iter = 500 # number of iteration

D=3
p=3       # number of predictor from the reduced model
k=4       # number of predictor from the full model

# the analysis prior
model="effectsize"    # "Moment" "effectsize"
dff=3
rscale=1
f_m=sqrt(.01)
e = NULL


#the design prior
de_an_prior=  1   # design and analysis prior being the same :1 ,not the same :2
model_d="Point"      # "Moment" "effectsize" "Point"
dff_d=7                   # must be > 3 if "Moment"
rscale_d=1.8
f_m_d= .10



f_ver(iter,D,p,k,dff,rscale,f_m,model,
                dff_d,rscale_d,f_m_d,model_d,de_an_prior,e)
# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


