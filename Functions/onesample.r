# The app is riddled with dt() warnings of the type:
#    In dt(t[i], df, ncp = delta * sqrt(df + 1)) :
#      full precision may not have been achieved in 'pnt{final}'
# This is why I used dnct() in the past. 
# I am bringing it back, slightly optimized.
# We can later talk of possible alternatives.
# Since sourcing these functions takes time, what we can do it create a mini-package
# just with these functions. Loading this package would make these function immediately
# available, no more loading time required.
sourceCpp("Functions/dnct.cpp")
sourceCpp("Functions/pnct.cpp")

############# prior density function #############
# Probability density function of non-local prior:
dnlp <-function(delta, mu, ta){
  #((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
  dmom(delta-mu, ta^2) # 'mombf' package (includes pmom(), etc)
}
# Probability density function of informed t prior:
tstude <- function(t, location = 0, scale = sqrt(2)/2, df = 1) {
  # gamma((df+1)/2) * ((df+((t-location)/scale)^2)/df)^(-((df+1)/2)) / (scale*sqrt(df*pi)*gamma(df/2))
  dnct((t-location)/scale, df, ncp = 0)/scale 
}

t1_prior<- function(delta, location, scale, dff, model){
  switch(model,
         "Cauchy"         = tstude(delta, location, scale, 1),
         "Normal"         = dnorm (delta, location, scale),
         "NLP"            = dnlp  (delta, location, scale),
         "t-distribution" = tstude(delta, location, scale, dff))
}

############# the Bayes Factor #############

t1_BF10 <-function(t, df, model, location, scale, dff, hypothesis){
  bound  <- switch(hypothesis,
                   ">"  = c(a = 0, b = Inf),
                   "<"  = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x <- numeric(length(t))
  
  # Normalize the prior outside the for-loop:
  # normalization  <- integrate(function(delta) t1_prior(delta, location, scale, dff, model),lower = bound[1], upper = bound[2])$value
  # For all priors, the prior integrates to 1 when a = -Inf, b = Inf.
  # For all priors, we use their CDFs when either a = 0 or b = 0 and minimize manual integrations.
  # Note: pmom() errors at -Inf and Inf, so we avoid it below.
  normalization <- if (hypothesis == "!=") 1 else 
    switch(model,
           "Cauchy"         = pcauchy(bound[2], location, scale)     - pcauchy(bound[1], location, scale),
           "Normal"         = pnorm (bound[2], location, scale)      - pnorm (bound[1], location, scale),
           "NLP"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = pnct((bound[2] - location) / scale, dff, 0) - pnct((bound[1] - location) / scale, dff, 0))
  
  for(i in 1:length(t)){
    int  <- function(delta)
      dnct(t[i], df, ncp = delta * sqrt(df+1)) * t1_prior(delta, location, scale, dff, model)/normalization

    # Removed stop.on.error = FALSE as it is bad form:
    x[i] <- integrate(int, lower = bound[1], upper = bound[2], rel.tol = 1e-8)$value / dt(t[i], df, ncp = 0)
    
  }
  return(x)
}



############# bound function  #############
# for finding the t value such that BF10 = D (code optimized):
t1_BF10_bound <- function(D, df, model, location, scale, dff, hypothesis) {
  Bound_finding <- function(t) t1_BF10(t, df, model, location, scale, dff, hypothesis) - D
  
  if (hypothesis == "!=") {
    x <- tryCatch(uniroot(Bound_finding, lower = -8, upper = 0)$root, error = function(e) NA)
    y <- tryCatch(uniroot(Bound_finding, lower =  0, upper = 8)$root, error = function(e) NA)
    results <- c(x, y)
  } else {
    x <- tryCatch(uniroot.all(Bound_finding, lower = -8, upper = 8), error = function(e) NA)
    results <- x
  }
  
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")
  
  BF.vals  <- t1_BF10(results, df, model, location, scale, dff, hypothesis)
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  
  return(results[BF.close])
}

# finding the t that correspond to BF01 = D is the same as 
# finding the t that corresponds to BF10 = 1/D:
t1_BF01_bound <- function(D, df, model, location, scale, dff, hypothesis) {
  t1_BF10_bound(1 / D, df, model, location, scale, dff, hypothesis)
}

# p(BF01>D|H0)
# t is the t-value lead to BF = b based on the bound functions (optimized):
t1_TNE <- function(t, df) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  if (length(t) == 2) return(pt(max(t), df) - pt(min(t), df))
  
  # length(t) = 1:
  return(if (t > 0) pt(t, df) else 1 - pt(t, df))
}


# p(BF10>D|H1)
# Argument 'hypothesis' is fully determined by the length and sign of the t values.
# I removed it as a function argument and compute it inside t1_TPE() instead.
t1_TPE <- function(t, df, model, location, scale, dff) {
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
  
  normalization <- integrate(function(delta) t1_prior(delta, location, scale, dff, model),
                             lower = bound[1], upper = bound[2])$value
  
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
  error <- if (model == "NLP" && scale < 0.3) 1e-14 else if (scale > 0.3) .Machine$double.eps^0.25 else 1e-8
  
  integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value
}


# p(BF01>D|H1)
# Similar as above:
t1_FNE <- function(t, df, model, location, scale, dff){
  
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  hypothesis <- if (length(t) == 2) "!=" else if (t >= 0) ">" else "<"
  
  if (model == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(max(t), df, ncp) - pnct(min(t), df, ncp))
    # Length 1:
    return(if (t >= 0) pnct(t, df, ncp) else 1 - pnct(t, df, ncp))
  }
  
  bound  <- switch(hypothesis,
                   ">"  = c(a = 0,    b = Inf),
                   "<"  = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf))
  
  # in some weird situation uniroot.all gives 4 t-values leading to BF =b 
  # but after some checking, only the middle two are the t-values leading to BF =b.
  # Jorge: This would be better dealt with inside t1_BF10_bound() and t1_BF01_bound()
  # I'll remove it and deal with the problem when I see it:
  # t <- as.numeric(t)
  # if (length(t) == 4) {  # Corrected condition check
  #   t = t[2:3]
  # }
  
  normalization <- integrate(function(delta) t1_prior(delta, location, scale, dff, model),
                             lower = bound[1], upper = bound[2])$value
  
  int <- if (length(t) == 2) { # two-sided test
    function(delta) {
      pro1 <- pnct(max(t), df, delta * sqrt(df + 1))
      pro2 <- pnct(min(t), df, delta * sqrt(df + 1))
      (pro1 - pro2) * t1_prior(delta, location, scale, dff, model) / normalization
    }
  } else if (t >= 0) { # one-sided test with delta > 0
    function(delta) pnct(t, df, delta * sqrt(df + 1)) * t1_prior(delta, location, scale, dff, model) / normalization
  } else {             # one-sided test with delta < 0
    function(delta) (1 - pnct(t, df, delta * sqrt(df + 1))) * t1_prior(delta, location, scale, dff, model) / normalization
  }
  
  # setting error value such that error are prevented:
  error <- if (model == "NLP" && scale < 0.3) 1e-14 else if (scale > 0.3) .Machine$double.eps^0.25 else 1e-8
  
  integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value
} 


# p(BF10>D|H0)
t1_FPE <- function(t, df) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  # if (length(t) == 4) t <- t[2:3]
  
  if (length(t) == 2) return(pt(min(t), df) + (1 - pt(max(t), df)))
  
  # length(t) = 1:
  return(if (t > 0) 1 - pt(t, df) else pt(t, df))
}



# Finding the degree of freedom that ensure p(BF10>D|H1) > targeted probability:
t1_N_finder <- function(D, target, model, location, scale, dff, hypothesis, 
                        model_d, location_d, scale_d, dff_d, de_an_prior, alpha) {
  #de_an_prior: 1 = design prior and analysis priors are the same, 0 otherwise.
  
  # error prevention 
  # sometimes, power can go higher than .8 with N= 2 already.
  # So, N should be returned now, otherwise, error will occur later.
  # Jorge: Below, I added design prior for N=2 too.
  lower <- 2
  upper <- 100000
  
  t2 <- t1_BF10_bound(D, df = lower, model, location, scale, dff, hypothesis)
  p2 <- if (de_an_prior == 1) 
    t1_TPE(t2, df = lower, model, location, scale, dff) else
      t1_TPE(t2, df = lower, model_d, location_d, scale_d, dff_d)
  if (p2 > target) return(lower)
  
  Power_root <- function(df) {
    t <- t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)
    if (de_an_prior == 1) 
      t1_TPE(t, df, model, location, scale, dff) - target else
        t1_TPE(t, df, model_d, location_d, scale_d, dff_d) - target
  }
  
  ## finding the required df, i will do the plus one to get the N in the later function.
  # Jorge: I'm not sure whether you meant later inside t1_N_finder(), or elsewhere.
  #        Anyway, it's a pity if we don't fix it here already, super easy to do right now.
  #        So I did it, see below.
  
  # Jorge: 'df.power' makes for a more accurate name than 'N'.
  df.power <- uniroot(Power_root, lower = lower, upper = upper)$root
  
  ## checking if the N lead to an acceptable alpha level
  t   <- t1_BF10_bound(D, df.power, model, location, scale, dff, hypothesis)
  FPE <- t1_FPE(t, df.power, model, location, scale, dff)
  if (FPE <= alpha) return(df.power + 1)
  
  # if the FPE > alpha, then we search for another df 
  # Jorge: 'alpha.root' is better than 'alpha_bound'.
  alpha.root <- function(df) {
    t <- t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)
    t1_FPE(t, df, model, location, scale, dff) - alpha
  }
  
  # Jorge: 'df.alpha' is better than 'NN'.
  df.alpha <- uniroot(alpha.root, lower = df.power, upper = upper)$root
  return(df.alpha + 1)

  #N = auto_uniroot_fixed_lower_going_up(Power_root,fixed_lower = 2, upper = 2000, step = 500, max_attempts = 16)
}



############ probability table
t1_Table <- function(D,target,model,location,scale,dff, hypothesis,
                  model_d,location_d,scale_d,dff_d, de_an_prior,N, mode_bf ,alpha ){
  # mode_bf == "0" means that the design analysis is done for a fixed N
  # Otherwise, it searches N where power > targeted power with FPE < FP
  if (mode_bf == "0"){
  df = N-1
 }else {
  
  df =  t1_N_finder(D,target,model,location,scale,dff, hypothesis,
                model_d,location_d,scale_d,dff_d, de_an_prior ,alpha )
  # here is the exact df that we need, required N is df +1 , which is done at the very end
  df = ceiling(df+1) -1
 }
  
  
  
  t10 = t1_BF10_bound(D ,df ,model ,location ,scale ,dff ,hypothesis )
  BF_D= t01 = t1_BF01_bound(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)
  
  
  if (de_an_prior == 1 ){
  TPE = t1_TPE(t10,df ,model ,location ,scale ,dff ,hypothesis )
  FPE = t1_FPE(t10,df,model ,location ,scale,dff, hypothesis)
  
  max_BF = 1/t1_BF10(0,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )

  if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found"  )) {
    FNE = 0
    TNR = 0
  }else{
    FNE = t1_FNE(t01,df,model ,location ,scale,dff, hypothesis)
    TNR = t1_TNE(t01 , df,model ,location ,scale,dff , hypothesis)
  }
  table <- data.frame(
    H11 = TPE,
    H10 = FNE,
    H00 = TNR,
    H01 = FPE,
    N =  ceiling(df+1))
  colnames(table) <- c(sprintf("p(BF10> %0.f|H1)",D), sprintf("p(BF01> %0.f|H1)",D), sprintf("p(BF01> %0.f|H0)",D), sprintf("p(BF10> %0.f|H0)",D), " Reqiured N")
  
  } else {
    TPE = t1_TPE(t10,df,model_d ,location_d ,scale_d ,dff_d ,hypothesis )
    FPE = t1_FPE(t10,df,model ,location ,scale,dff, hypothesis)
    
    max_BF = 1/t1_BF10(0,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )
    BF_D = t1_BF10_bound(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)
    
    if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found"  )) {
      FNE = 0
      TNR = 0
    }else{
      
      FNE = t1_FNE(t01,df,model_d,location_d ,scale_d,dff_d, hypothesis)
      TNR = t1_TNE(t01 , df,model ,location ,scale,dff , hypothesis)
    }
    table <- data.frame(
      H11 = TPE,
      H10 = FNE,
      H00 = TNR,
      H01 = FPE,
      N =  ceiling(df+1))
    colnames(table) <- c(sprintf("p(BF10> %0.f|H1)",D), sprintf("p(BF01> %0.f|H1)",D), sprintf("p(BF01> %0.f|H0)",D), sprintf("p(BF10> %0.f|H0)",D), " Reqiured N")
    
    }
  return(table)
}




# plot for the selected prior 
t1_prior_plot<-function(D ,target,model ,location ,scale,dff , hypothesis,model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior){
  par(mfrow = c(1, 1))
  # the bound for the plot
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = 5),
                   "<" = c(a = -5, b = 0),
                   "!=" = c(a = -5, b = 5)
  )
  tt= seq(bound[1],bound[2],.01)
  
  # the bound for the intergrals
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  
  normalization  <- integrate(function(delta)t1_prior(delta, location,scale,dff,model),lower = bound[1],upper = bound[2])$value
  # anlaysis prior
  prior_DELTA = NA
  prior_DELTA = t1_prior(tt, location,scale,dff,model)/normalization
  
  plot(tt,prior_DELTA,xlab= bquote(bold(delta)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~delta~" under the alternative hypothesis")),frame.plot = FALSE)
  if (de_an_prior ==0){
    
    if (model_d != "Point"){
      # design prior
      normalization_d  <- integrate(function(delta)t1_prior(delta,location_d,scale_d,dff_d,model_d),lower = bound[1],upper = bound[2])$value
      
      prior_DELTA_D = t1_prior(tt,location_d,scale_d,dff_d,model_d)/normalization_d
      
      plot(tt,prior_DELTA,xlab= bquote(bold(delta)),ylim = c(0,max(max(prior_DELTA_D),max(prior_DELTA))),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~delta~" under the alternative hypothesis")),frame.plot = FALSE)
      
      lines(tt,prior_DELTA_D ,lty = 2)
      
    } else{
      # point prior
      plot(tt,prior_DELTA,xlab= bquote(bold(delta)),ylim = c(0,max(prior_DELTA)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~delta~" under the alternative hypothesis")),frame.plot = FALSE)
      arrows(x0=location_d, y0=0, x1=location_d, y1=max(prior_DELTA), length=0.2, code=2, col="black", lwd=1,,lty = 2)
    }
    
    legend("topright", 
           legend = c("Analysis prior", "Design prior"), 
           lty = c(1, 2), 
           col = c("black", "black"),
           bty = "n") 
    
  }
}

# plots for showing the relationship between BF and t-values 

bf10_t1 <-function(D =3,df, target,model = "NA",location =0,scale=.707,dff = 1, hypothesis ){
  
  tt= seq(from = -5,to = 5,.2)
  # finding the t leading to BF10 = b
  BF_D = t1_BF10_bound(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)
  BF10 = t1_BF10(tt,df,model ,location,scale,dff,hypothesis)
  
  # setting the title of the plot
  if (length(BF_D) == 1){
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D, digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f ",D,BF_D)
  } else {
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f or %.3f ",D,BF_D[1],BF_D[2])
  }
  par(mfrow = c(1, 2))
  plot(tt,log10(BF10),xlab= "t-value",type="l", ylab = expression("logarithm of BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
  abline(v = BF_D)
  axis(1, c(-5,5))
  if (length(BF_D) != 0 ){
  axis(1, round(BF_D,2))}
  
  max_BF = 1/t1_BF10(0,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis ="!=" )
  BF_D = t1_BF01_bound(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)

  
  
  plot(tt,log10(1/BF10),xlab= "t-value",type="l",main = "",frame.plot = FALSE,ylab = bquote("logarithm of BF"[0][1]),xaxt = "n")
  axis(1, c(-5,5))
  if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found" ) ) {
    main = bquote(bold("It is impossible to have BF"[0][1]~"="~.(D)))
    title(main = main)
    #sprintf("It is impossible to have BF01 = %.3f ",D)
  } else      {
    abline(v = BF_D)
    axis(1, round(BF_D,2))
    if (length(BF_D) == 1){
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when t = "~.(format(BF_D, digits = 4))))
      title(main = main)
    } else {
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when t = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
      title(main = main)
    }}


  

  
}

# Power curve functions
Power_t1<-function(D,model,location,scale,dff, hypothesis,
                   model_d,location_d,scale_d,dff_d, de_an_prior,N){
  
  smin = 2
  smax = N*1.2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  power =  array(NA, dim = c(length(sdf)))
  
  for ( i in 1:length(sdf)){
    t = t1_BF10_bound(D ,sdf[i],model ,location ,scale  ,dff ,hypothesis )
    power[i] = switch(de_an_prior,
                      "1" = t1_TPE(t,sdf[i],model  ,location  ,scale,dff, hypothesis),
                      "0" = t1_TPE(t,sdf[i],model_d  ,location_d  ,scale_d,dff_d, hypothesis))
    
    
  }
  plot(sdf+1,power,type="l",main = "",frame.plot = FALSE,xlab = "Sample size", ylab = "Probability of True positive evidence",xlim = c(1,max(sdf)), 
       ylim = c(0,1) )
  
}

