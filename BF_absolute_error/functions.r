#### t-test ####
t1_TPE_abs_error <- function(t, df, model, location, scale, dff){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  hypothesis <- if (length(t) == 2) "!=" else if (t >= 0) ">" else "<"

  if (model == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(min(t), df, ncp) + (1 - pnct(max(t), df, ncp)))
    # Length 1:
    return(if (t >= 0) 1 - pnct(t, df, ncp) else pnct(t, df, ncp))
  }

  bound  <- switch(hypothesis,
                   ">"  = c(a = 0,    b = Inf),
                   "<"  = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf))

  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "NLP"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

  int <- if (length(t) == 2) { # two-sided test
    function(delta) {
      pro1 <- 1 - pnct(max(t), df, delta * sqrt(df + 1))
      pro2 <-     pnct(min(t), df, delta * sqrt(df + 1))
      (pro1 + pro2) * t1_prior(delta, location, scale, dff, model) / normalization
    }
  } else if (t >= 0) { # one-sided test with delta > 0
    function(delta) (1 - pnct(t, df, delta * sqrt(df + 1))) * t1_prior(delta, location, scale, dff, model) / normalization
  } else {             # one-sided test with delta < 0
    function(delta) pnct(t, df, delta * sqrt(df + 1)) * t1_prior(delta, location, scale, dff, model) / normalization
  }

  # setting error value such that error are prevented:
  #error <- if (model == "NLP" && scale < 0.3) 1e-14 else if (scale > 0.3) .Machine$double.eps^0.25 else 1e-8
  error = 1e-4
  stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$abs.error

}

t2_TPE_abs_error <-function(t,n1,r,model ,location ,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (model == "Point"){
    pro = switch(hypothesis,
                 "!="= pnct(min(t),df,ncp = location*constant,lower  = T)+pnct(max(t),df,ncp = location*constant,lower  = F),
                 ">" = pnct(t,df,ncp = location*constant,lower  = F),
                 "<" = pnct(t,df,ncp = location*constant,lower  = T))
    return(pro)
  }

  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )



  x = NULL

  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "NLP"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(hypothesis,
                  "!=" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = FALSE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 + pro2
                  },
                  ">" = pnct(t, df, ncp = ncp, lower = FALSE),
                  "<" = pnct(t, df, ncp = ncp, lower = TRUE)
    )

    pro * t1_prior(delta, location, scale, dff, model) / normalization
  }

  error = 1e-4
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = error,stop.on.error=FALSE)$abs.error

  return(x)

}

#### correlation ####
r_TPE_abs_error <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)

  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {p_cor(max(r),location,n,lower.tail = F)+ p_cor(min(r),location,n,lower.tail = T)},
               ">"  = {p_cor(r,location,n,lower.tail =F)},
               "<"  = {p_cor(r,location,n,lower.tail =T)}
    )
    return(x)
  }

  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  normalization <-   normalization <- if (hypothesis == "!=") {
    switch(model,
           "d_beta"   = 1,
           "beta" = 1,
           "NLP"   = { pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})

  }else{
    switch(model,
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,-1,1)-p_beta(bound[1], 1/k,1/k,-1,1) ,
           "beta" = p_beta(bound[2], alpha, beta,-1,1)-p_beta(bound[1], alpha, beta,-1,1),
           "NLP"   = {pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})
  }
  int <- function(rho) {
    prob <- switch(hypothesis,
                   "!=" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                     p_cor(min(r), rho, n, lower.tail = TRUE),
                   ">"  = p_cor(r, rho, n, lower.tail = FALSE),
                   "<"  = p_cor(r, rho, n, lower.tail = TRUE)
    )

    prob * r_prior(rho, k, location, scale, dff, model, alpha, beta,min(bound),max(bound)) / normalization
  }
  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-4)$abs.error
  return(x)

}
#### f-test ####
F_TPE_abs_error<-function(f,q,m,dff,rscale,f_m,model){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  if (model == "Point"){
    x = stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = F)
    return(x)
  }
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,model)
  }
  x = stats::integrate(int,lower = 0,upper = Inf)$abs.error
  return(x)
}
#### one-proportion ####

bin_TPE_abs_error<-function(x,n,h0,alpha,beta,location,scale,model,hypothesis){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  if (model =="Point"){
    TPE = switch(hypothesis,
                 "!=" = {

                   switch(length(x)==2,
                          "1" ={stats::pbinom(min(x),n,location,lower.tail = T)+ stats::pbinom(max(x)-1,n,location,lower.tail = F)},
                          "0"=  {
                            switch(x/n>location,
                                   "1" = stats::pbinom(x-1,n,location,lower.tail = F),
                                   "0" = stats::pbinom(x,n,location,lower.tail = T))

                          })
                 },
                 ">"  = {stats::pbinom(x-1,n,location,lower.tail = F)},
                 "<"  = {stats::pbinom(x,n,location,lower.tail = T)}
    )
    return(TPE)
  }

  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = 0, b = h0),
                   "!=" = c(a = 0, b = 1)
  )
  normalization <- if (hypothesis == "!=") {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
  }
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
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
                  ">" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "<" = stats::pbinom(x, n, prop, lower.tail = TRUE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalization
  }

  TPE = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-5)$abs.error

  return(TPE)

}
