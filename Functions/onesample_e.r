robust_uniroot <- function(f, lower, upper_start = 200, max_attempts = 20, step = 500, ...) {
  upper <- upper_start
  attempt <- 1
  
  repeat {
    result <- tryCatch(
      {
        uniroot(f, lower = lower, upper = upper, ...)$root
      },
      error = function(e) {
        NULL
      }
    )
    
    if (!is.null(result)) {
      return(result)
    }
    
    if (attempt >= max_attempts) {
      stop("Failed to find root after increasing upper bound ", max_attempts, " times.")
    }
    
    upper <- upper + step
    attempt <- attempt + 1
  }
}

te_prior<- function(delta,scale,dff,model){
  
  switch(model,
        "Cauchy"        = tstude(delta,0, scale,1),
        "Normal"        = dnorm(delta,0,scale),
        "NLP"            = dnlp(delta,0,scale),
        "t-distribution" = tstude(delta,0,scale,dff))
  
  
}

t1e_BF10i <-function(t,df,model ,scale,dff , hypothesis,e ){
  bound_h1  <- switch(hypothesis,
                      "!=" = c(a = e[1], b = e[2]),
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = e)
                      )
  
  bound_h0  <- switch(hypothesis,
                      "!=" = c(a = e[1], b = e[2]),
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0)
                      )

  normalizationh1 <- switch(hypothesis,
                            "!=" = 1 - switch(model,
                                              "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                              "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                              "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                              "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff)),
                            "<"  = switch(model,
                                          "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                          "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                          "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                          "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff)),
                            ">"  = switch(model,
                                          "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                          "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                          "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                          "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff))
  )
  
  
  # H0 Normalization
  normalizationh0 <- switch(model,
                            "Cauchy"         = pcauchy(bound_h0[2], 0, scale) - pcauchy(bound_h0[1], 0, scale),
                            "Normal"         = pnorm  (bound_h0[2], 0, scale) - pnorm  (bound_h0[1], 0, scale),
                            "NLP"            = pmom   (bound_h0[2], tau = scale^2) - pmom   (bound_h0[1], tau = scale^2),
                            "t-distribution" = pt     (bound_h0[2] / scale, df = dff) - pt  (bound_h0[1] / scale, df = dff)
  )
  
  int  <- function(delta){
    dt(t,df,ncp = delta *sqrt(df+1))* te_prior(delta,scale,dff,model)/normalizationh1}
  
   error = 1e-4

  if (hypothesis == "!="){
  lh1 = integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value 
  }else{
    lh1 = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value
    
  }

  
  int  <- function(delta){
    dt(t,df,ncp = delta *sqrt(df+1))* te_prior(delta,scale,dff,model)/normalizationh0}
  
  lh0 = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(lh1/lh0)
}

t1e_BF10 <-function(t,df,model,scale,dff , hypothesis,e ){

  x <- sapply(t, function(ti) t1e_BF10i(ti, df, model, scale, dff, hypothesis, e))
  return(x)
}
#
t1e_BF10_bound <-function(D, df,model,scale,dff , hypothesis,e){
  y <- numeric(0)
  Bound_finding <-function(t){
    t1e_BF10(t,df,model,scale,dff , hypothesis,e )- D
  }
  
  switch(hypothesis,
         "!=" ={
           x <- tryCatch(uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
           y <- tryCatch(uniroot(Bound_finding, lower =  0, upper = 20)$root, error = function(e) NA)
         },
         ">"={
           x <- tryCatch(uniroot(Bound_finding, lower = 0, upper = 20)$root, error = function(e) NA)
         },
         "<" = {
           x <- tryCatch(uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
         })


  results <- c(x, y)
  
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")
  
  BF.vals  <- t1e_BF10(results,df,model,scale,dff , hypothesis,e )
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
  
}


t1e_BF01_bound <-function(D, df,model,scale,dff , hypothesis,e){
  t1e_BF10_bound(1/D, df,model,scale,dff , hypothesis,e)
}



t1e_TPE <-function(t,df,model ,scale,dff , hypothesis ,e,location){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (model == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(min(t), df, ncp) + (1 - pnct(max(t), df, ncp)))
    # Length 1:
    return(if (t >= 0) 1 - pnct(t, df, ncp) else pnct(t, df, ncp))
  }

  
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = e),
                      "!=" = c(a = e[1], b = e[2])
  )
  
  normalizationh1 <- switch(hypothesis,
                            "!=" = 1 - switch(model,
                                              "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                              "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                              "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                              "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff)),
                            "<"  = switch(model,
                                          "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                          "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                          "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                          "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff)),
                            ">"  = switch(model,
                                          "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                          "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                          "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                          "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff))
  )
  
  
  x = NULL

  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)
    
    pro <- switch(hypothesis,
                  "!=" = pnct(max(t), df, ncp = ncp, lower  = FALSE) +
                    pnct(min(t), df, ncp = ncp, lower  = TRUE),
                  ">"  = pnct(t, df, ncp = ncp, lower  = FALSE),
                  "<"  = pnct(t, df, ncp = ncp, lower  = TRUE)
    )
    
    pro * te_prior(delta, scale, dff, model) / normalizationh1
  }
  
   error = 1e-4

  x <- switch(hypothesis,
              "!=" = integrate(int, -Inf, bound_h1[1], rel.tol = error)$value +
                integrate(int, bound_h1[2], Inf, rel.tol = error)$value,
              "<"  = ,
              ">"  = integrate(int, bound_h1[1], bound_h1[2], rel.tol = error)$value
             
  )
  return(x) 
  
} 

t1e_FNE <-function(t,df,model ,scale,dff , hypothesis ,e,location){
  
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  if (model == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(max(t), df, ncp) - pnct(min(t), df, ncp))
    # Length 1:
    return(if (t >= 0) pnct(t, df, ncp) else 1 - pnct(t, df, ncp))
  }
  
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = e),
                      "!=" = c(a = e[1], b = e[2])
  )

  normalizationh1 <- switch(hypothesis,
                            "!=" = 1 - switch(model,
                                              "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                              "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                              "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                              "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff)),
                            "<"  = switch(model,
                                          "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                          "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                          "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                          "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff)),
                            ">"  = switch(model,
                                          "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                          "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                          "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                          "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff))) 
  x = NULL
  
  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)
    
    pro <- switch(hypothesis,
                   "!=" = pnct(max(t), df, ncp = ncp, lower  = TRUE) -
                     pnct(min(t), df, ncp = ncp, lower  = TRUE),
                   ">"  = pnct(t, df, ncp = ncp, lower  = TRUE),
                   "<"  = pnct(t,df, ncp = ncp, lower  = FALSE)
    )
    
    pro * te_prior(delta, scale, dff, model) / normalizationh1
  }
  
   error = 1e-4
  
  x <- switch(hypothesis,
              "!=" = integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value ,
              "<"  = ,
              ">"  = integrate(int, bound_h1[1], bound_h1[2], rel.tol = error)$value)
  return(x) 
  
} 

t1e_TNE <-function(t,df,model ,scale,dff , hypothesis ,e){
  
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )
  normalizationh0 <- switch(model,
                            "Cauchy"         = pcauchy(bound_h0[2], 0, scale) - pcauchy(bound_h0[1], 0, scale),
                            "Normal"         = pnorm  (bound_h0[2], 0, scale) - pnorm  (bound_h0[1], 0, scale),
                            "NLP"            = pmom   (bound_h0[2], tau = scale^2) - pmom   (bound_h0[1], tau = scale^2),
                            "t-distribution" = pt     (bound_h0[2] / scale, df = dff) - pt  (bound_h0[1] / scale, df = dff)
  )
  
  x = NULL
  
  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)
    pro <- switch(hypothesis,
                   "!=" = pnct(max(t), df, ncp, lower  = TRUE) - 
                     pnct(min(t), df, ncp, lower  = TRUE),
                   ">"  = pnct(t, df, ncp, lower  = TRUE),
                   "<"  = pnct(t, df, ncp, lower  = FALSE),
                   stop("Unsupported hypothesis")
    )
    
    pro * te_prior(delta, scale, dff, model) / normalizationh0
  }
  
   error = 1e-4
  
  x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error)$value

  return(x) 
  
} 

t1e_FPE <-function(t,df,model ,scale,dff , hypothesis ,e){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )
  
  normalizationh0 <- switch(model,
                            "Cauchy"         = pcauchy(bound_h0[2], 0, scale) - pcauchy(bound_h0[1], 0, scale),
                            "Normal"         = pnorm  (bound_h0[2], 0, scale) - pnorm  (bound_h0[1], 0, scale),
                            "NLP"            = pmom   (bound_h0[2], tau = scale^2) - pmom   (bound_h0[1], tau = scale^2),
                            "t-distribution" = pt     (bound_h0[2] / scale, df = dff) - pt  (bound_h0[1] / scale, df = dff)
  )
  x = NULL
  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)
    
    pro <- switch(hypothesis,
                   "!=" = pnct(max(t), df, ncp, lower  = FALSE) + pnct(min(t), df, ncp, lower  = TRUE),
                   ">"  = pnct(t, df, ncp, lower  = FALSE),
                   "<"  = pnct(t, df, ncp, lower  = TRUE),
                   stop("Unsupported hypothesis")
    )
    
    pro * te_prior(delta, scale, dff, model) / normalizationh0
  }

   error = 1e-4
  
  x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  
  return(x) 
  
} 

t1e_N_finder<-function(D,target,model,scale,dff, hypothesis,e ,
                   model_d,scale_d,dff_d, de_an_prior,location_d  ,alpha){
  
  lower <- 2
  t2 <-t1e_BF10_bound(D, lower,model,scale,dff , hypothesis,e)
  p2 <- if (de_an_prior == 1)
    t1e_TPE(t2,lower,model ,scale,dff , hypothesis ,e,location) else
    t1e_TPE(t2,lower,model ,scale_d,dff_d , hypothesis ,e,location_d)
  if (p2 > target) return(lower)

  Power_root <- function(df) {

    t <- t1e_BF10_bound(D, df, model, scale, dff, hypothesis, e)
    
    pro <- if (de_an_prior == 1) {
      t1e_TPE(t, df, model, scale, dff, hypothesis, e)
    } else {
      t1e_TPE(t, df, model_d, scale_d, dff_d, hypothesis, e, location_d)
    }
    
    target - pro
  }
  
  df.power <- robust_uniroot(Power_root, lower = 2)
  t <- t1e_BF10_bound(D,df.power,model,scale,dff,hypothesis ,e )
  FPE <-t1e_FPE(t,df.power,model ,scale,dff , hypothesis ,e)
  if (FPE <= alpha) return(df.power + 1)
  
  alpha.root <- function(df) {
      t <- t1e_BF10_bound(D,df,model,scale,dff,hypothesis ,e )
      pro <- t1e_FPE(t , df , model=model , scale=scale,dff=dff, hypothesis,e)
      return(pro - alpha)
    }
  df.alpha <- uniroot(alpha.root, lower = df.power, upper = upper)$root
  return(df.alpha+1)
  
}


  
t1e_table<-function(D,target,model,scale,dff, hypothesis,e ,
                    model_d,scale_d,dff_d, de_an_prior,df,mode_bf,location_d ,alpha ){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  if (mode_bf == "0") df <- N-1 else {
    N  <- ceiling(t1e_N_finder(D,target,model,scale,dff, hypothesis,e ,
                               model_d,scale_d,dff_d, de_an_prior ,location_d,alpha ))
    df <- N -1
  }
  # t bounds:
  t10 <- t1e_BF10_bound(D, df,model,scale,dff , hypothesis,e)
  t01 <- t1e_BF01_bound(D, df,model,scale,dff , hypothesis,e)
  
  # max BF10 possible:
  max_BF <- 1 / t1e_BF10i(0,df,model ,scale,dff , hypothesis,e )
  BF_D   <- t10
  
  # FPE and TPE:
  FPE       <- t1e_FPE(t10,df,model ,scale,dff , hypothesis ,e)
  if (de_an_prior == 1) {
    TPE       <- t1e_TPE(t10,df,model ,scale,dff , hypothesis ,e,location_d)
    TPR_model <- model
    TPR_scale <- scale
    TPR_dff   <- dff
  } else {
    TPE       <- t1e_TPE(t10,df,model_d ,scale_d,dff_d , hypothesis ,e,location_d)
    TPR_model <- model_d
    TPR_location   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }
  
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t1e_FNE(t01,df,TPR_model ,TPR_scale,TPR_dff , hypothesis ,e,TPR_location)
    TNE <-  t1e_TNE(t01,df,model ,scale,dff , hypothesis ,e)
  }
  # table:
  tab.names <- c(
    sprintf("p(BF10 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H0)", D),
    sprintf("p(BF10 > %0.f | H0)", D),
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, N, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
  
}

compute.prior.density.te.h1 <- function(tt, model, scale, dff, hypothesis,e,location) {
  if (model == "Point") return(rep(NA, length(tt)))
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = -e),
                      "!=" = c(a = e[1], b = e[2])
  )
  
  normalizationh1 <- switch(hypothesis,
                            "!=" = 1 - switch(model,
                                              "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                              "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                              "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                              "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff)),
                            "<"  = switch(model,
                                          "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                          "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                          "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                          "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff)),
                            ">"  = switch(model,
                                          "Cauchy"         = pcauchy(bound_h1[2], 0, scale) - pcauchy(bound_h1[1], 0, scale),
                                          "Normal"         = pnorm (bound_h1[2], 0, scale) - pnorm (bound_h1[1], 0, scale),
                                          "NLP"            = pmom(bound_h1[2], tau = scale^2) - pmom(bound_h1[1], tau = scale^2),
                                          "t-distribution" = pt(bound_h1[2] / scale, df = dff) - pt(bound_h1[1] / scale, df = dff))
  )
  
  
  #prior_h1<-te_prior(tt,scale,dff,model) / normalizationh1
  prior_h1<-te_prior(tt,scale,dff,model)
  switch(hypothesis,
         "!=" = { prior_h1[tt>min(bound_h1)&tt<max(bound_h1)]=0 },
         ">" = { prior_h1[tt<bound_h1[1]]=0 },
         "<" = { prior_h1[tt>bound_h1[2]]=0 }
         )
  prior_h1
}
compute.prior.density.te.h0 <- function(tt, model, scale, dff, hypothesis,e,location) {
  if (model == "Point") return(rep(NA, length(tt)))
  bound_h0  <- switch(hypothesis,
                      "!=" = c(a = e[1], b = e[2]),
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0)
  )
  # H0 Normalization
  normalizationh0 <- switch(model,
                            "Cauchy"         = pcauchy(bound_h0[2], 0, scale) - pcauchy(bound_h0[1], 0, scale),
                            "Normal"         = pnorm  (bound_h0[2], 0, scale) - pnorm  (bound_h0[1], 0, scale),
                            "NLP"            = pmom   (bound_h0[2], tau = scale^2) - pmom   (bound_h0[1], tau = scale^2),
                            "t-distribution" = pt     (bound_h0[2] / scale, df = dff) - pt  (bound_h0[1] / scale, df = dff)
  )
  
  
  #prior_h0 <- te_prior(tt,scale,dff,model) / normalizationh0
  prior_h0 <- te_prior(tt,scale,dff,model)
  switch(hypothesis,
         "!=" = { prior_h0[!(tt>min(bound_h0)&tt<max(bound_h0))]=0},
         ">" = { prior_h0[tt>bound_h0[2]]=0 },
         "<" = { prior_h0[tt<bound_h0[1]]=0 }
         
  )
  prior_h0
}


t1e_prior_plot <- function(model, scale, dff, hypothesis, e,
                           de_an_prior, model_d, scale_d, dff_d, location_d) {
  par(mfrow = c(1, 1))
  
  plot.bounds <- switch(hypothesis,
                        ">"  = c(0, 5),
                        "<"  = c(-5, 0),
                        "!=" = c(-5, 5))
  tt <- seq(plot.bounds[1], plot.bounds[2], 0.01)
  
  prior.analysis.h1 <- compute.prior.density.te.h1(tt, model, scale, dff, hypothesis, e, location_d)
  prior.analysis.h0 <- compute.prior.density.te.h0(tt, model, scale, dff, hypothesis, e)
  prior.design <- if (de_an_prior == 0 && model_d != "Point") {
    compute.prior.density.te.h1(tt, model_d, scale_d, dff_d, hypothesis, e, location_d)
  } else {
    rep(NA, length(tt))
  }
  
  ylim.max <- max(prior.analysis.h1, prior.analysis.h0, prior.design, na.rm = TRUE)
  
  # Base plot
  plot(tt, prior.analysis.h1, type = "l", lwd = 2,
       xlab = expression(bold(delta)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~delta~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))
  
  lines(tt, prior.analysis.h0, lty = 2, col = "black", lwd = 2)
  
  # Optional: design prior
  legend.labels <- c("H1 - Analysis Prior", "H0 - Analysis Prior")
  legend.cols   <- c("black", "black")
  legend.lty    <- c(1, 2)
  legend.lwd    <- c(2, 2)
  
  if (de_an_prior == 0) {
    if (model_d == "Point") {
      arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim.max,
             length = 0.1, col = "gray", lty = 2)
    } else {
      lines(tt, prior.design, lty = 1, col = "gray", lwd = 3)
    }
    
    # Add design prior to legend
    legend.labels <- c(legend.labels, "Design prior")
    legend.cols   <- c(legend.cols, "gray")
    legend.lty    <- c(legend.lty, 1)
    legend.lwd    <- c(legend.lwd, 2)
  }
  
  legend("topleft",
         legend = legend.labels,
         col = legend.cols,
         lty = legend.lty,
         lwd = legend.lwd,
         bty = "n")
}

te1_BF <-function(D,df,model ,scale,dff , hypothesis ,e){
  tt <- seq(-5, 5, 0.2)
  par(mfrow = c(1, 2))
  # Compute BF10 and t-bounds:
  BF10   <- t1e_BF10(tt,df,model,scale,dff , hypothesis,e )
  t.BF10 <- t1e_BF10_bound(D, df,model,scale,dff , hypothesis,e)
  
  # Left plot - BF10:
  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(D)~"when t = "~.(format(t.BF10, digits = 4))))
  } else {
    bquote(bold("BF"[10]~"="~.(D)~"when t = "~.(format(t.BF10[1], digits = 4))~"or"~.(format(t.BF10[2], digits = 4))))
  }
  plot(tt, BF10, type = "l", log = "y", xlab = "t-value", ylab = expression("BF"[10]),
       main = main.bf10, frame.plot = FALSE, xaxt = "n",xlim = c(-5,5))
  abline(v = t.BF10)
  axis(1, c(-5, 5))
  if (length(t.BF10)) axis(1, round(t.BF10, 2))
  
  # right plot - BF01:
  BF01   <- 1 / BF10
  t.BF01 <- t1e_BF01_bound(D, df,model,scale,dff , hypothesis,e)
  # Check if BF01 = D is possible:
  max.BF01   <- 1 / t1_BF10(0, df, model, location, scale, dff, hypothesis = "!=")
  impossible <- (max.BF01 < D || identical(t.BF01, "bound cannot be found"))
  
  plot(tt, BF01, type = "l", log = "y", xlab = "t-value", ylab = bquote("BF"['01']),
       main = "", frame.plot = FALSE, xaxt = "n")
  axis(1, c(-5, 5))
  if (impossible) {
    title(main = bquote(bold("It is impossible to have BF"[01]~"="~.(D))))
  } else {
    abline(v = t.BF01)
    axis(1, round(t.BF01, 2))
    main.bf01 <- if (length(t.BF01) == 1) {
      bquote(bold("BF"['01']~"="~.(D)~"when t = "~.(format(t.BF01, digits = 4))))
    } else {
      bquote(bold("BF"['01']~"="~.(D)~"when t = "~.(format(t.BF01[1], digits = 4))~"or"~.(format(t.BF01[2], digits = 4))))
    }
    title(main = main.bf01)
  }
  
}
Power_t1e<-function(D,model,location,scale,dff, hypothesis,
                   model_d,location_d,scale_d,dff_d, de_an_prior,N,e){
  par(mfrow = c(1, 1))
  # df range to evaluate power:
  df.min     <- 2
  df.max     <- ceiling(N * 1.2)
  dfs        <- seq(df.min, df.max, length.out = 30)
  power.vals <- numeric(length(dfs))
  
  for (i in seq_along(dfs)) {
    t <- t1e_BF10_bound(D, dfs[i],model,scale,dff , hypothesis,e)
    # Choose correct design prior:
    power.vals[i] <- if (de_an_prior == 1)
      t1e_TPE(t,dfs[i],model ,scale,dff , hypothesis ,e,location_d) else
        t1e_TPE(t,dfs[i],model_d ,scale_d,dff_d , hypothesis ,e,location_d)
  }
  
  plot(dfs + 1, power.vals, type = "l",
       xlab = "Sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  
}

