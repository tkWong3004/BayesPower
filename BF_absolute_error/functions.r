#### t-test ####
t1_TPE_abs_error <- function(t, df, prior_analysis, location, scale, dff){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (prior_analysis == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(BayesPower:::pnct(min(t), df, ncp) + (1 - BayesPower:::pnct(max(t), df, ncp)))
    # Length 1:
    return(if (t >= 0) 1 - BayesPower:::pnct(t, df, ncp) else BayesPower:::pnct(t, df, ncp))
  }

  bound  <- switch(alternative,
                   "greater"  = c(a = 0,    b = Inf),
                   "less"  = c(a = -Inf, b = 0),
                   "two.sided" = c(a = -Inf, b = Inf))

  normalization <- if (alternative == "two.sided") 1 else
    switch(prior_analysis,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) BayesPower:::pmom(bound[2]-location, tau=scale^2) else 1-BayesPower:::pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

  int <- if (length(t) == 2) { # two-sided test
    function(delta) {
      pro1 <- 1 - BayesPower:::pnct(max(t), df, delta * sqrt(df + 1))
      pro2 <-     BayesPower:::pnct(min(t), df, delta * sqrt(df + 1))
      (pro1 + pro2) * BayesPower:::t1_prior(delta, location, scale, dff, prior_analysis) / normalization
    }
  } else if (t >= 0) { # one-sided test with delta > 0
    function(delta) (1 - BayesPower:::pnct(t, df, delta * sqrt(df + 1))) * BayesPower:::t1_prior(delta, location, scale, dff, prior_analysis) / normalization
  } else {             # one-sided test with delta < 0
    function(delta) BayesPower:::pnct(t, df, delta * sqrt(df + 1)) * BayesPower:::t1_prior(delta, location, scale, dff, prior_analysis) / normalization
  }

  error = 1e-4
  stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$abs.error

}

t2_TPE_abs_error <-function(t,n1,r,prior_analysis ,location ,scale,dff , alternative ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (prior_analysis == "Point"){
    pro = switch(alternative,
                 "two.sided"= BayesPower:::pnct(min(t),df,ncp = location*constant,lower  = T)+BayesPower:::pnct(max(t),df,ncp = location*constant,lower  = F),
                 "greater" = BayesPower:::pnct(t,df,ncp = location*constant,lower  = F),
                 "less" = BayesPower:::pnct(t,df,ncp = location*constant,lower  = T))
    return(pro)
  }

  bound  <- switch(alternative,
                   "greater" = c(a = 0, b = Inf),
                   "less" = c(a = -Inf, b = 0),
                   "two.sided" = c(a = -Inf, b = Inf)
  )



  x = NULL

  normalization <- if (alternative == "two.sided") 1 else
    switch(prior_analysis,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) BayesPower:::pmom(bound[2]-location, tau=scale^2) else 1-BayesPower:::pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(alternative,
                  "two.sided" = {
                    pro1 <- BayesPower:::pnct(max(t), df, ncp = ncp, lower = FALSE)
                    pro2 <- BayesPower:::pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 + pro2
                  },
                  "greater" = BayesPower:::pnct(t, df, ncp = ncp, lower = FALSE),
                  "less" = BayesPower:::pnct(t, df, ncp = ncp, lower = TRUE)
    )

    pro * BayesPower:::t1_prior(delta, location, scale, dff, prior_analysis) / normalization
  }

  error = 1e-4
  if (prior_analysis == "Moment" & scale <.3 ){
    error = 1e-14
  }
  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = error,stop.on.error=FALSE)$abs.error

  return(x)

}

#### correlation ####
r_TPE_abs_error <-function(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)

  if (prior_analysis =="Point"){
    x = switch(alternative,
               "two.sided" = {BayesPower:::p_cor(max(r),location,n,lower.tail = F)+ BayesPower:::p_cor(min(r),location,n,lower.tail = T)},
               "greater"  = {BayesPower:::p_cor(r,location,n,lower.tail =F)},
               "less"  = {BayesPower:::p_cor(r,location,n,lower.tail =T)}
    )
    return(x)
  }

  bound  <- switch(alternative,
                   "greater" = c(a = h0, b = 1),
                   "less" = c(a = -1, b = h0),
                   "two.sided" = c(a = -1, b = 1)
  )
  normalization <-   normalization <- if (alternative == "two.sided") {
    switch(prior_analysis,
           "d_beta"   = 1,
           "beta" = 1,
           "Moment"   = { BayesPower:::pmom(bound[2]-location, tau=scale^2)-BayesPower:::pmom(bound[1]-location, tau=scale^2)})

  }else{
    switch(prior_analysis,
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,-1,1)-p_beta(bound[1], 1/k,1/k,-1,1) ,
           "beta" = p_beta(bound[2], alpha, beta,-1,1)-p_beta(bound[1], alpha, beta,-1,1),
           "Moment"   = {BayesPower:::pmom(bound[2]-location, tau=scale^2)-BayesPower:::pmom(bound[1]-location, tau=scale^2)})
  }
  int <- function(rho) {
    prob <- switch(alternative,
                   "two.sided" = BayesPower:::p_cor(max(r), rho, n, lower.tail = FALSE) +
                     BayesPower:::p_cor(min(r), rho, n, lower.tail = TRUE),
                   "greater"  = BayesPower:::p_cor(r, rho, n, lower.tail = FALSE),
                   "less"  = BayesPower:::p_cor(r, rho, n, lower.tail = TRUE)
    )

    prob * BayesPower:::r_prior(rho, k, location, scale, dff, prior_analysis, alpha, beta,min(bound),max(bound)) / normalization
  }
  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-4)$abs.error
  return(x)

}
#### f-test ####
F_TPE_abs_error<-function(f,q,m,dff,rscale,f_m,prior_analysis){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  if (prior_analysis == "Point"){
    x = stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = F)
    return(x)
  }
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*BayesPower:::F_prior(fsq,q,dff,rscale,f_m,prior_analysis)
  }
  x = stats::integrate(int,lower = 0,upper = Inf)$abs.error
  return(x)
}
#### one-proportion ####

bin_TPE_abs_error<-function(x,n,h0,alpha,beta,location,scale,prior_analysis,alternative){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  if (prior_analysis =="Point"){
    TPE = switch(alternative,
                 "two.sided" = {

                   switch(length(x)==2,
                          "1" ={stats::pbinom(min(x),n,location,lower.tail = T)+ stats::pbinom(max(x)-1,n,location,lower.tail = F)},
                          "0"=  {
                            switch(x/n>location,
                                   "1" = stats::pbinom(x-1,n,location,lower.tail = F),
                                   "0" = stats::pbinom(x,n,location,lower.tail = T))

                          })
                 },
                 "greater"  = {stats::pbinom(x-1,n,location,lower.tail = F)},
                 "less"  = {stats::pbinom(x,n,location,lower.tail = T)}
    )
    return(TPE)
  }

  bound  <- switch(alternative,
                   "greater" = c(a = h0, b = 1),
                   "less" = c(a = 0, b = h0),
                   "two.sided" = c(a = 0, b = 1)
  )
  normalization <- if (alternative == "two.sided") {
    switch(prior_analysis,
           "Moment"   = BayesPower:::pmom(bound[2]-location, tau=scale^2)-BayesPower:::pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(prior_analysis,
           "Moment"   = BayesPower:::pmom(bound[2]-location, tau=scale^2)-BayesPower:::pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
  }
  int <- function(prop) {
    pro <- switch(alternative,
                  "two.sided" = {
                    if (length(x) == 2) {
                      stats::pbinom(min(x), n, prop, lower.tail = TRUE) +
                        stats::pbinom(max(x) - 1, n, prop, lower.tail = FALSE)
                    } else {
                      mapply(function(x_i, n_i, p_i) {
                        if (x_i / n_i > location) {
                          stats::pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                        } else {
                          stats::pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                        }
                      }, x, n, prop)
                    }
                  },
                  "greater" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "less" = stats::pbinom(x, n, prop, lower.tail = TRUE)
    )

    pro * BayesPower:::bin_prior(prop, alpha, beta, location, scale, prior_analysis) / normalization
  }

  TPE = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-5)$abs.error

  return(TPE)

}
