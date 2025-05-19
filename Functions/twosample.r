
# the Bayes Factor

t2_BF10 <-function(t,n1,r,model ,location,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )

  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Cauchy"         = pcauchy(bound[2], location, scale)     - pcauchy(bound[1], location, scale),
           "Normal"         = pnorm (bound[2], location, scale)      - pnorm (bound[1], location, scale),
           "NLP"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = pt((bound[2] - location) / scale, dff, 0) - pt((bound[1] - location) / scale, dff, 0))
  
  error = 1e-5
  x <- sapply(t, function(ti) {
    int <- function(delta) {
      dt(ti, df, ncp = delta * constant) * t1_prior(delta, location, scale, dff, model) / normalization
    }
    
    integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value / 
      dt(ti, df, ncp = 0)
  })
  
  return(x)
}



# finding the t that correspond to BF10=D
t2_BF10_bound <-function(D, n1,r,model ,location ,scale,dff , hypothesis){
  y <- numeric(0)
  Bound_finding <-function(t){
    t2_BF10(t,n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )- D
  }
  x <- tryCatch(uniroot(Bound_finding, lower = -7, upper = 0)$root, error = function(e) NA)
  y <- tryCatch(uniroot(Bound_finding, lower =  0, upper = 7)$root, error = function(e) NA)
  results <- c(x, y)
  
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")
  
  BF.vals  <- t2_BF10(results,n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  
  return(results[BF.close])
}


# finding the t that correspond to BF01=D
t2_BF01_bound <-function(D , n1,r,model ,location ,scale,dff , hypothesis){
  t2_BF10_bound(1/D, n1,r,model ,location ,scale,dff , hypothesis)
}


# p(BF01>D|H0)
t2_TNE <- function(t , n1,r,hypothesis){
  n2 = n1*r
  df = n1+n2-2

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  pro <- switch(hypothesis,
                "!=" = pt(max(t), df) - pt(min(t), df),
                ">"  = pt(t, df),
                "<"  = 1 - pt(t, df)
  )
  
  return(pro)
  
}

# p(BF10>D|H1)
t2_TPE <-function(t,n1,r,model ,location ,scale,dff , hypothesis ){
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
           "Cauchy"         = pcauchy(bound[2], location, scale)     - pcauchy(bound[1], location, scale),
           "Normal"         = pnorm (bound[2], location, scale)      - pnorm (bound[1], location, scale),
           "NLP"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = pt((bound[2] - location) / scale, dff, 0) - pt((bound[1] - location) / scale, dff, 0))
  
  
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
  x = integrate(int,lower = bound[1],upper = bound[2], rel.tol = error,stop.on.error=FALSE)$value
  
  return(x) 
  
} 


# p(BF01>D|H1)
t2_FNE<-function(t,n1,r,model ,location ,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  if (model == "Point"){
    pro = switch(hypothesis,
                 "!="=  pnct(max(t),df,ncp = location*constant,lower  = T) - pnct(min(t),df,ncp = location*constant,lower  = T),
                 ">" = pnct(t,df,ncp = location*constant,lower  = T),
                 "<" = pnct(t,df,ncp = location*constant,lower  = F))
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
           "Cauchy"         = pcauchy(bound[2], location, scale)     - pcauchy(bound[1], location, scale),
           "Normal"         = pnorm (bound[2], location, scale)      - pnorm (bound[1], location, scale),
           "NLP"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = pt((bound[2] - location) / scale, dff, 0) - pt((bound[1] - location) / scale, dff, 0))
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
    
    pro * t1_prior(delta, location, scale, dff, model) / normalization
  }
  
  
  x = integrate(int,lower = bound[1],upper = bound[2],stop.on.error = FALSE)$value

  return(x) 
} 


# p(BF10>D|H0)
t2_FPE <- function(t,n1,r, hypothesis){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  pro <- switch(hypothesis,
                "!=" = pt(max(t), df = df, lower.tail = FALSE) +
                  pt(min(t), df = df, lower.tail = TRUE),
                ">"  = pt(t, df = df, lower.tail = FALSE),
                "<"  = pt(t, df = df, lower.tail = TRUE)
  )
  return(pro)

}

# Finding the degree of freedom that ensure p(BF10>D|H1) > targeted probability

t2_N_finder<-function(D,r,target,model,location,scale,dff, hypothesis ,
                   model_d,location_d,scale_d,dff_d,de_an_prior ,alpha){
 
  lower <- 2
  upper <- 10000
  t2 <- t2_BF10_bound(D, lower,r,model ,location ,scale,dff , hypothesis)
  p2 <- if (de_an_prior == 1)
    t2_TPE(t2 , n1=lower,r , model , location ,scale,dff , hypothesis) else
      t2_TPE(t2 , n1=lower,r , model_d , location_d,scale_d,dff_d, hypothesis)
  if (p2 > target) return(lower)
  Power_root <- function(n1) {
    i<<-i+1
    t <- t2_BF10_bound(D, n1,r,model ,location ,scale,dff , hypothesis)
    if (de_an_prior == 1)
      t2_TPE(t , n1,r , model , location ,scale,dff , hypothesis) - target else
        t2_TPE(t , n1,r , model_d , location_d,scale_d,dff_d, hypothesis) - target
  }
  N1.power <-  uniroot(Power_root,lower = lower,upper =  upper)$root
  t  <-  t2_BF10_bound(D,  N1.power,r,model ,location ,scale,dff , hypothesis)
  FPE <- t2_FPE(t,N1.power,r, hypothesis)
  if (FPE <= alpha) return(N1.power)
  alpha.root <- function(n1) {
    t <- t2_BF10_bound(D,  n1,r,model ,location ,scale,dff , hypothesis)
    t2_FPE(t,n1,r, hypothesis) - alpha
  }
  
  N1.alpha <- uniroot(alpha.root, lower = N1.power, upper = upper)$root
  return(N1.alpha  )
}
  
  


# probability table
t2_Table <- function(D,r,target,model,location,scale,dff, hypothesis,
                  model_d,location_d,scale_d,dff_d, de_an_prior,N1,N2, mode_bf ,alpha ){

  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  
  if (mode_bf == 1) {
    n1 <- ceiling(t2_N_finder(D, r, target, model, location, scale, dff, 
                              hypothesis, model_d, location_d, scale_d, dff_d, 
                              de_an_prior, alpha))
    n2 <- n1 * r
  } else {
    n1 <- N1
    n2 <- N2
    r  <- n2 / n1
  }
  
  # t bounds:
  t10 <- t2_BF10_bound(D, n1,r,model,location,scale,dff , hypothesis)
  t01 <- t2_BF01_bound(D, n1,r,model,location,scale,dff , hypothesis)
  
  # max BF10 possible:
  max_BF <- 1 / t2_BF10(0,n1,r,model ,location,scale,dff , hypothesis )
  BF_D   <- t10
  
  # FPE and TPE:
  FPE       <- t2_FPE(t10,n1,r, hypothesis)
  if (de_an_prior == 1) {
    TPE       <- t2_TPE(t10,n1,r,model ,location ,scale,dff , hypothesis )
    TPR_model <- model
    TPR_loc   <- location
    TPR_scale <- scale
    TPR_dff   <- dff
  } else {
    TPE       <- t2_TPE(t10,n1,r,model_d ,location_d ,scale_d,dff_d , hypothesis )
    TPR_model <- model_d
    TPR_loc   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }
  
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t2_FNE(t01, n1,r,TPR_model, TPR_loc, TPR_scale, TPR_dff, hypothesis )
    TNE <- t2_TNE(t01 , n1,r,hypothesis)
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


# plots for showing the relationship between BF and t-values 

t2_bf10 <-function(D ,n1,r, target,model ,location ,scale,dff  , hypothesis ){
  par(mfrow = c(1, 2))
  tt <- seq(-5, 5, 0.2)
  
  # Compute BF10 and t-bounds:
  BF10   <- t2_BF10(tt,n1,r,model ,location,scale,dff,hypothesis)
  t.BF10 <- t2_BF10_bound(D, n1,r,model ,location ,scale,dff , hypothesis)
  # Left plot - BF10:
  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(D)~"when t = "~.(format(t.BF10, digits = 4))))
  } else {
    bquote(bold("BF"[10]~"="~.(D)~"when t = "~.(format(t.BF10[1], digits = 4))~"or"~.(format(t.BF10[2], digits = 4))))
  }
  plot(tt, BF10, type = "l", log = "y", xlab = "t-value", ylab = expression("BF"[10]),
       main = main.bf10, frame.plot = FALSE, xaxt = "n")
  abline(v = t.BF10)
  axis(1, c(-5, 5))
  if (length(t.BF10)) axis(1, round(t.BF10, 2))
  
  # Left plot - BF01:
  BF01   <- 1 / BF10
  t.BF01 <- t2_BF01_bound(D, n1,r,model ,location ,scale,dff , hypothesis)
  
  # Check if BF01 = D is possible:
  max.BF01   <- 1 / t2_BF10 (0,n1,r,model ,location,scale,dff ,"!=")
  impossible <- (hypothesis == "!=") && (max.BF01 < D || identical(t.BF01, "bound cannot be found"))
  
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

Power_t2<-function(D,model,location,scale,dff, hypothesis,
                   model_d,location_d,scale_d,dff_d, de_an_prior,n1,r){
  par(mfrow = c(1, 1))
  Total_ = n1 + n1*r
  smin = 4
  smax = Total_*1.2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  sn1 = sdf/(1+r)
  power.vals =  array(NA, dim = c(length(sdf)))
  
  for ( i in 1:length(sdf)){
    t = t2_BF10_bound(D , sn1[i],r,model ,location ,scale,dff , hypothesis)
    power.vals[i] = switch(de_an_prior,
                      "1" = t2_TPE(t,sn1[i],r,model ,location ,scale,dff , hypothesis ),
                      "0" = t2_TPE(t,sn1[i],r,model_d ,location_d ,scale_d,dff_d , hypothesis ))
    
    
  }
  plot(sdf, power.vals, type = "l",
       xlab = "Total sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  
}

