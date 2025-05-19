
t2e_BF10i <-function(t,n1,r,model ,scale,dff , hypothesis,e ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = -e),
                      "!=" = c(a = e[1], b = e[2])
  )
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = -e, b = 0),
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
  
  
  
  normalizationh0 <- switch(model,
                            "Cauchy"         = pcauchy(bound_h0[2], 0, scale) - pcauchy(bound_h0[1], 0, scale),
                            "Normal"         = pnorm  (bound_h0[2], 0, scale) - pnorm  (bound_h0[1], 0, scale),
                            "NLP"            = pmom   (bound_h0[2], tau = scale^2) - pmom   (bound_h0[1], tau = scale^2),
                            "t-distribution" = pt     (bound_h0[2] / scale, df = dff) - pt  (bound_h0[1] / scale, df = dff))
  
  

  int  <- function(delta){
    dt(t,df,ncp = delta *constant)* te_prior(delta,scale,dff,model)/normalizationh1}

  error = 1e-4

  if (hypothesis == "!="){
  lh1 = integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value 
  }else{
    lh1 = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value}
    
  
  int  <- function(delta){
    dt(t,df,ncp = delta *constant)* te_prior(delta,scale,dff,model)/normalizationh0}
  
  lh0 = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(lh1/lh0)
}

t2e_BF10 <-function(t,n1,r,model,scale,dff , hypothesis,e ){
  sapply(t, function(ti) t2e_BF10i(ti,n1,r,model ,scale,dff , hypothesis,e ))
}

t2e_BF10_bound <-function(D, n1,r,model,scale,dff , hypothesis,e){
  
  y <- numeric(0)
  Bound_finding <-function(t){
    t2e_BF10(t,n1,r,model,scale,dff , hypothesis,e )- D
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
  BF.vals  <- t2e_BF10(results,n1,r,model,scale,dff , hypothesis,e )
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

# finding the t that correspond to BF10=D
t2e_BF01_bound <-function(D, n1,r,model,scale,dff , hypothesis,e){
  t2e_BF10_bound(1/D, n1,r,model,scale,dff , hypothesis,e)
  
}

t2e_TPE <-function(t,n1,r,model ,scale,dff , hypothesis ,e,location){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {pnct(min(t),df,ncp= location*constant,lower = T)+ pnct(max(t),df,ncp=location*constant,lower = F)},
               "<"  = {pnct(t,df,ncp = location *constant,lower  = T)},
               ">"  = {pnct(t,df,ncp = location *constant,lower  = F)}
               )
    return(x)
  }
  
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
    
    pro * te_prior(delta, scale, dff, model) / normalizationh1
  }
  
  
  error = 1e-4
  
  if (hypothesis == "!="){
    x = integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value 
  }else{
    x = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value
    
  }
  return(x) 
  
} 

t2e_FNE <-function(t,n1,r,model ,scale,dff , hypothesis ,e,location){
  
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {pnct(min(t),df,ncp= location*constant,lower = T)- pnct(max(t),df,ncp=location*constant,lower = T)},
               "<"  = {pnct(t,df,ncp = location *constant,lower  = F)},
               ">"  = {pnct(t,df,ncp = location *constant,lower  = T)}
    )
    return(x)
  }
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
  
  int <- function(delta) {
    ncp <- delta * constant
    
    pro <- switch(hypothesis,
                  "!=" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = TRUE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 - pro2
                  },
                  ">" = pnct(t, df, ncp = ncp, lower = TRUE),
                  "<" = pnct(t, df, ncp = ncp, lower = FALSE)
    )
    
    pro * te_prior(delta, scale, dff, model) / normalizationh1
  }
  
  error = 1e-4
  if (hypothesis == "!="){
    x = integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value 
  }else{
    x = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value
    
  }
  return(x) 
  
} 


t2e_TNE <-function(t,n1,r,model ,scale,dff , hypothesis ,e){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = -e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )
  
  normalizationh0 <- switch(model,
                            "Cauchy"         = pcauchy(bound_h0[2], 0, scale) - pcauchy(bound_h0[1], 0, scale),
                            "Normal"         = pnorm  (bound_h0[2], 0, scale) - pnorm  (bound_h0[1], 0, scale),
                            "NLP"            = pmom   (bound_h0[2], tau = scale^2) - pmom   (bound_h0[1], tau = scale^2),
                            "t-distribution" = pt     (bound_h0[2] / scale, df = dff) - pt  (bound_h0[1] / scale, df = dff))
  
  
  int <- function(delta) {
    ncp <- delta * constant
    
    pro <- switch(hypothesis,
                  "!=" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = TRUE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 - pro2
                  },
                  ">" = pnct(t, df, ncp = ncp, lower = TRUE),
                  "<" = pnct(t, df, ncp = ncp, lower = FALSE)
    )
    
    pro * te_prior(delta, scale, dff, model) / normalizationh0
  }
  error = 1e-4
    x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(x) 
  
} 

t2e_FPE <-function(t,n1,r,model ,scale,dff , hypothesis ,e){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = -e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )
  
  normalizationh0 <- switch(model,
                            "Cauchy"         = pcauchy(bound_h0[2], 0, scale) - pcauchy(bound_h0[1], 0, scale),
                            "Normal"         = pnorm  (bound_h0[2], 0, scale) - pnorm  (bound_h0[1], 0, scale),
                            "NLP"            = pmom   (bound_h0[2], tau = scale^2) - pmom   (bound_h0[1], tau = scale^2),
                            "t-distribution" = pt     (bound_h0[2] / scale, df = dff) - pt  (bound_h0[1] / scale, df = dff))
  
  
  
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
    
    pro * te_prior(delta, scale, dff, model) / normalizationh0
  }
  error = 1e-4
  x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  
  return(x) 
  
} 



t2e_N_finder<-function(D,r,target,model,scale,dff, hypothesis,e ,
                   model_d,scale_d,dff_d, de_an_prior,location,alpha ){
  
  lower <- 2
  t2 <-t2e_BF10_bound(D, lower,r,model,scale,dff , hypothesis,e)
  p2 <- if (de_an_prior == 1)
    t2e_TPE (t2,lower,r,model ,scale,dff , hypothesis,e ,location) else
      t2e_TPE (t2,lower,r,model_d ,scale_d,dff_d , hypothesis,e ,location)
  if (p2 > target) return(lower)
  
  Power_root <- function(n1) {
    
    t <- t2e_BF10_bound(D, n1,r,model,scale,dff , hypothesis,e)
    
    pro <- if (de_an_prior == 1) {
      t2e_TPE (t,n1,r,model ,scale,dff , hypothesis,e )
    } else {
      t2e_TPE (t,n1,r,model_d ,scale_d,dff_d , hypothesis,e ,location)
    }
    
    target - pro
  }
  N1.power <- robust_uniroot(Power_root, lower = 2)
  t <- t2e_BF10_bound(D, N1.power,r,model,scale,dff , hypothesis,e)
  FPE <-t2e_FPE(t,N1.power,r,model ,scale,dff , hypothesis ,e)

  if (FPE <= alpha) return(N1.power + 1)
  
  alpha.root <- function(n1) {
    t <- t2e_BF10_bound(D, n1,r,model,scale,dff , hypothesis,e)
    pro <- t2e_FPE(t,n1,r,model ,scale,dff , hypothesis ,e)
    return(pro - alpha)
  }
  N1.alpha <- robust_uniroot(alpha.root , lower = N1.power)
  return(N1.alpha)
  }

  
t2e_table<-function(D,r,target,model,scale,dff, hypothesis,e ,
                    model_d,scale_d,dff_d, de_an_prior,mode_bf,location,N1,N2,alpha ){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  
  if (mode_bf == 1){
    
    n1 = ceiling(t2e_N_finder(D,r,target,model,scale,dff, hypothesis,e ,
                      model_d,scale_d,dff_d, de_an_prior,location,alpha ))
    n2 = n1*r
  } else {
    n1 = N1
    n2 = N2
    r= n2/n1
  }
  # t bounds:
  t10 <- t2e_BF10_bound(D, n1,r,model,scale,dff , hypothesis,e)
  t01 <- t2e_BF01_bound(D, n1,r,model,scale,dff , hypothesis,e)
  
  # max BF10 possible:
  max_BF <- 1 / t2e_BF10i(0,n1,r,model ,scale,dff , hypothesis,e )
  BF_D   <- t10
  
  # FPE and TPE:
  FPE       <- t2e_FPE(t10,n1,r,model ,scale,dff , hypothesis ,e)
  if (de_an_prior == 1) {
    TPE       <- t2e_TPE(t10,n1,r,model ,scale,dff , hypothesis,e ,location)
    TPR_model <- model
    TPR_loc   <- location
    TPR_scale <- scale
    TPR_dff   <- dff
  } else {
    TPE       <- t2e_TPE(t10,n1,r,model_d ,scale_d,dff_d , hypothesis,e ,location)
    TPR_model <- model_d
    TPR_loc   <- location
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }
  
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t2e_FNE(t01,n1,r,TPR_model ,TPR_scale,TPR_dff , hypothesis ,e,TPR_location)
    TNE <- t2e_TNE(t01,n1,r,model ,scale,dff , hypothesis ,e)
  }
  
  # table:
  tab.names <- c(
    sprintf("p(BF10 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H0)", D),
    sprintf("p(BF10 > %0.f | H0)", D),
    "Required N1",
    "Required N2"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n1,n2, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}


t2e_BF <-function(D,n1,r,model ,scale,dff , hypothesis ,e){
  tt <- seq(-5, 5, 0.2)
  par(mfrow = c(1, 2))
  # Compute BF10 and t-bounds:
  BF10   <- t2e_BF10(tt, n1,r,model,scale,dff , hypothesis,e)
  t.BF10 <- t2e_BF10_bound(D, n1,r,model,scale,dff , hypothesis,e)
  
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
  t.BF01 <- t2e_BF01_bound(D, N1.power,r,model,scale,dff , hypothesis,e)
  # Check if BF01 = D is possible:
  max.BF01   <- 1 / t2e_BF10i(0,n1,r,model ,scale,dff , hypothesis,e )
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


Power_t2e<-function(D,model,location,scale,dff, hypothesis,
                   model_d,location_d,scale_d,dff_d, de_an_prior,n1,r){
  par(mfrow = c(1, 1))
  Total_ = n1 + n1*r
  smin = 4
  smax = Total_*1.2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  sn1 = sdf/(1+r)
  power =  array(NA, dim = c(length(sdf)))
  
  for ( i in 1:length(sdf)){
    t = t2e_BF10_bound(D, sn1[i],r,model,scale,dff , hypothesis,e)
    power[i] = switch(de_an_prior,
                      "1" = t2e_TPE(t,sn1[i],r,model ,scale,dff , hypothesis ,e,location),
                      "0" = t2e_TPE(t,sn1[i],r,model_d ,scale_d,dff_d , hypothesis ,e,location_d))
    
    
  }
  plot(sdf,power, type = "l",
       xlab = "Total sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  
}
