############# prior density function #############
# likelihood of non-local prior
dnlp <-function(delta,mu,ta){
  ((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
}
# likelihood of informed t prior
tstude <- function(t, location = 0, scale = sqrt(2)/2, df = 1) {
  gamma((df+1)/2) * ((df+((t-location)/scale)^2)/df)^(-((df+1)/2)) / (scale*sqrt(df*pi)*gamma(df/2))
  
  #dt((t-location)/scale,df,ncp = 0)/scale
}

t1_prior<- function(delta, location,scale,dff,model){
  
  switch(model,
         "Cauchy" = tstude(delta,location, scale,1),
         "Normal"   = dnorm(delta,location,scale),
         "NLP"   = dnlp(delta,location,scale),
         "t-distribution" = tstude(delta,location,scale,dff))
  
  
}

############# the Bayes Factor #############

t1_BF10 <-function(t,df,model ,location,scale,dff , hypothesis ){
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x <- numeric(length(t))
  
  for(i in 1:length(t)){
    
    normalization  <- integrate(function(delta)t1_prior(delta, location,scale,dff,model),lower = bound[1],upper = bound[2])$value
    
    int  <- function(delta){
      dt(t[i],df,ncp = delta *sqrt(df+1))* t1_prior(delta, location,scale,dff,model)/normalization}

    error = 1e-8
    x[i]= integrate(int,lower = bound[1],upper = bound[2], rel.tol=error,stop.on.error = F)$value/dt(t[i],df,ncp = 0)
    
  }
  return(x)
}



############# bound function  #############
# for finding the t value being equal to a value of BF

t1_BF10_bound <-function(D, df,model ,location ,scale,dff , hypothesis){

  y <- numeric(0)
  Bound_finding <-function(t){
    t1_BF10(t,df,model,location,scale,dff, hypothesis  )- D
  }
  
  if (hypothesis=="!="){
    x <- tryCatch( uniroot(Bound_finding, lower = -8, upper = 0)$root, error=function(e){})
    
    y <- tryCatch( uniroot(Bound_finding, lower = 0, upper = 8)$root, error=function(e){})
    
    
  } else {
    x <- tryCatch( uniroot.all(Bound_finding, lower = -8, upper = 8), error=function(e){})

  }
  if (length(x) ==0){
    x = "bound cannot be found"
    return(x)
  } 
  if (length(y) == 1){
    x = cbind(x,y)
  }
  ##preventing possible wrong t-value
  BF = t1_BF10(x,df,model,location,scale,dff,hypothesis )
  x = x[round(BF,2)==round(D,2)]
  
  
  if (length(x) ==0){
    x = "bound cannot be found"
    return(x)
  } 
  return(x)
}



# finding the t that correspond to BF01=D
t1_BF01_bound <-function(D , df,model ,location ,scale,dff , hypothesis){
  
  Bound_finding <-function(t){
    1/t1_BF10(t,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )-D
  }
  
   x = uniroot.all(Bound_finding, lower = -8, upper = 8)
  
  if (length(x) == 0 ){
    x = "bound cannot be found"
    return(x)
  }
  
  BF = 1/t1_BF10(x ,df,model,location,scale,dff, hypothesis )
  x = x[round(BF,1)== round(D,1)]
  
  return(x)
}



# p(BF01>D|H0)
# t is the t-value lead to BF = b based on the bound functions
t1_TNE <- function(t , df,model,location ,scale,dff , hypothesis){
  
  # where there is no t-value, then the probability is of course is 0 since there is no bound
  
  if (any(t == "bound cannot be found")){
    return(t)
  }
  
  if (length(t)>1){
    pro = pt(max(t),df) - pt(min(t),df)
  }else{
    # when t>0, the alternative hypothesis is delta >0
    # so , the true negative rate is based on left tail.
    if (t >0){
      pro = pt(t,df)
    } else {
      pro = 1-pt(t,df)
    }
  }
  
   
  
  return(pro)
}

# p(BF10>D|H1)
t1_TPE <-function(t,df,model ,location ,scale,dff , hypothesis ){
  
  if (any(t =="bound cannot be found" | length(t)==0)){
    t=0
    return(t)
  }
  
  if (model == "Point"){
    pro = switch(hypothesis,
                 "!="= pnct(t[t<0],df,ncp = location *sqrt(df+1),lower  = T)+pnct(t[t>0],df,ncp = location *sqrt(df+1),lower  = F),
                 ">" = pnct(t,df,ncp = location *sqrt(df+1),lower  = F),
                 "<" = pnct(t,df,ncp = location *sqrt(df+1),lower  = T))
    return(pro)
  }
  
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x = NULL

  normalization  <- integrate(function(delta)t1_prior(delta, location,scale,dff,model),lower = bound[1],upper = bound[2])$value
  
  if (length(t) > 1){
    # two-sided test
    int  <- function(delta){
      pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = F)
      pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1+pro2)* t1_prior(delta, location,scale,dff,model)/normalization }
    
    
    
  }  else if (t >= 0) { 
    # one-sided test with delta >0
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
      (pro1)* t1_prior(delta, location,scale,dff,model)/normalization }
    
  } else if (t < 0){
    # one-sided test with delta <0
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1)* t1_prior(delta, location,scale,dff,model)/normalization }
    
  }
  # setting error value such that error are prevented
  if (scale >.3){
    error = .Machine$double.eps^0.25
  }else {
    error = 1e-8
  }
  
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
  x = integrate(int,lower = bound[1],upper = bound[2], rel.tol = error,stop.on.error=FALSE)$value
  
  return(x) 
  
} 

# p(BF01>D|H1)
t1_FNE<-function(t,df,model ,location ,scale,dff , hypothesis ){
  
  if (any(t =="bound cannot be found" | length(t)==0)){
    t=0
    return(t)
  }
  
  if (model == "Point"){
    pro = switch(hypothesis,
                 "!="=  pnct(t[t>0],df,ncp = location *sqrt(df+1),lower  = T) - pnct(t[t<0],df,ncp = location *sqrt(df+1),lower  = T),
                 ">" = pnct(t,df,ncp = location *sqrt(df+1),lower  = T),
                 "<" = pnct(t,df,ncp = location *sqrt(df+1),lower  = F))
    return(pro)
  }
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x = NULL
  # in some weird situation uniroot.all gives 4 t-values leading to BF =b 
  # but after some checking, only the middle two are the t-values leading to BF =b.
  t <- as.numeric(t)
  if (length(t) == 4) {  # Corrected condition check
    t = t[2:3]
  }
  
  normalization  <- integrate(function(delta)t1_prior(delta, location,scale,dff,model),lower = bound[1],upper = bound[2])$value
  if (length(t)>1){
    
    int  <- function(delta){
      pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = T)
      pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1-pro2)* t1_prior(delta, location,scale,dff,model)/normalization }
    
  } else if (t>0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1)* t1_prior(delta, location,scale,dff,model)/normalization }
    
  } else if (t<0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
      (pro1)* t1_prior(delta, location,scale,dff,model)/normalization }
    
  }
  x = integrate(int,lower = bound[1],upper = bound[2],stop.on.error = FALSE)$value
  
  return(x) 
} 

# p(BF10>D|H0)
t1_FPE <- function(t,df,model ,location ,scale,dff, hypothesis){
  if (any(t =="bound cannot be found" | length(t)==0)){
    return(0)
  }
  

  if (hypothesis=="!="){
    
    if (any(t[1] == "bound cannot be found"|t[2] == "bound cannot be found")){
    return(0)
  }
  } else if (t == "bound cannot be found"){
    return(0)
  }
  
  if (length(t) == 4) {  # Corrected condition check
    t = t[2:3]
  }
  if (length(t)>1){
    pro = pt(max(t),df=df,lower.tail = F) +pt(min(t),df=df,lower.tail = T)
  }else if (t>0){
    pro = pt(t[t>0],df=df,lower.tail = F)
  }else if (t<0){
    pro = pt(t[t<0],df=df,lower.tail = T)
  } 
    return(pro)
    
  
  
}



# Finding the degree of freedom that ensure p(BF10>D|H1) > targeted probability

t1_N_finder<-function(D,target,model,location,scale,dff, hypothesis ,
                   model_d,location_d,scale_d,dff_d,de_an_prior,alpha){
  # error prevention 
  #sometimes, power can go higher than .8 with N= 2 already.
  # So, N should be returned now, otherwise , error will occur later
  
  lo = 2
  t = t1_BF10_bound(D=D,df=lo,model=model,location=location ,scale=scale ,dff = dff, hypothesis =hypothesis)
  
  pro = t1_TPE(t , df=lo , model=model , location=location ,scale=scale,dff=dff, hypothesis=hypothesis)
  
  if (pro>target){
    N = 2
    return(N)
  }
  
  
  Power_root <- function(df){
    
    #de_an_prior: 1 means design prior and analysis priors are the same
    # otherwise, they are different, so the calculate of power differs as well here.
    
    t = t1_BF10_bound(D,df,model,location ,scale ,dff, hypothesis )
    if (de_an_prior == 1){
    pro = t1_TPE(t , df , model , location ,scale,dff, hypothesis)
    }else{
      pro = t1_TPE(t , df=df , model=model_d , location=location_d ,scale=scale_d,dff=dff_d, hypothesis )
    }
    
    return(pro-target)}
  
  
  up= 100000
  ## finding the required df , i will do the plus one to get the N in the later function.
  N = uniroot(Power_root,lower = lo,upper =  up)$root
  
  
  ## checking if the N lead to an acceptable alpha level
  t = t1_BF10_bound(D=D,df=N,model=model,location=location ,scale=scale ,dff=dff, hypothesis =hypothesis)
  FPE =t1_FPE(t , df=N , model=model , location=location ,scale=scale,dff=dff, hypothesis=hypothesis )
  # if the PFE > alpha, then we search for another df 
  
  if (FPE<alpha){
    
    return(N)
  }else{
    alpha_bound <- function(df) {
      t <- t1_BF10_bound(D=D, df=df, model=model, location=location, scale=scale, dff=dff, hypothesis=hypothesis)
      pro <- t1_FPE(t, df=df, model=model, location=location, scale=scale, dff=dff, hypothesis=hypothesis)
      return(pro - alpha)
      
    }
    NN = uniroot(alpha_bound,lower = N,upper =  up)$root
    return(NN)
    }

  
  
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

