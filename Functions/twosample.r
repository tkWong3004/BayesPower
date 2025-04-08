
# likelihood of non-local prior
dnlp <-function(delta,mu,ta){
  ((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
}
# likelihood of informed t prior
tstude <- function(t, location = 0, scale = sqrt(2)/2, df = 1) {
  gamma((df+1)/2) * ((df+((t-location)/scale)^2)/df)^(-((df+1)/2)) / (scale*sqrt(df*pi)*gamma(df/2))
  
  #dt((t-location)/scale,df,ncp = 0)/scale
}

t2_prior<- function(delta, location,scale,dff,model){
  
  switch(model,
         "Cauchy" = tstude(delta,location, scale,1),
         "Normal"   = dnorm(delta,location,scale),
         "NLP"   = dnlp(delta,location,scale),
         "t-distribution" = tstude(delta,location,scale,dff))
  
  
}



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
  x <- numeric(length(t))
  
  for(i in 1:length(t)){
    
    
    normalization  <- integrate(function(delta)t2_prior(delta, location,scale,dff,model),lower = bound[1],upper = bound[2])$value
    
    int  <- function(delta){
      dt(t[i],df,ncp = delta *constant)* t2_prior(delta, location,scale,dff,model)/normalization}

    error = 1e-8
    if (model == "NLP" & scale <.3 ){
      error = 1e-14
    }
    x[i]= integrate(int,lower = bound[1],upper = bound[2], rel.tol=error,stop.on.error = F)$value/dt(t[i],df,ncp = 0)
    
  }
  return(x)
}



# finding the t that correspond to BF10=D
t2_BF10_bound <-function(D, n1,r,model ,location ,scale,dff , hypothesis){
  y <- numeric(0)
  Bound_finding <-function(t){
    t2_BF10(t,n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )- D
  }
  
  if (hypothesis=="!="){
    x <- tryCatch({
      uniroot(Bound_finding, lower = -8, upper = 0)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    y <- tryCatch({
      uniroot(Bound_finding, lower = 0, upper = 8)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
  } else {
    x <- tryCatch({
      uniroot.all(Bound_finding, lower = -8, upper = 8)
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
    
  }
  
  if (length(x) ==0){
    x = "bound cannot be found"
    return(x)
  } 
  if (length(y) == 1){
    x = cbind(x,y)
  }
  BF = t2_BF10(x,n1,r,model,location,scale,dff,hypothesis )
  x = x[round(BF,2)==round(D,2)]
  if (length(x) ==0){
    x = "bound cannot be found"
    return(x)
  } 
  return(x)
}


# finding the t that correspond to BF01=D

t2_BF01_bound <-function(D , n1,r,model ,location ,scale,dff , hypothesis){
  y <-numeric(0)
  Bound_finding <-function(t){
    1/t2_BF10(t,n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )-D
  }
  
  if (hypothesis=="!="){
    x <- tryCatch({
      uniroot(Bound_finding, lower = -8, upper = 0)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    y <- tryCatch({
      uniroot(Bound_finding, lower = 0, upper = 8)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
  } else {
    x <- tryCatch({
      uniroot.all(Bound_finding, lower = -8, upper = 8)
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
    
  }
  
  if (length(x) == 0 ){
    x = "bound cannot be found"
    return(x)
  }
  if (length(y) == 1){
    x = cbind(x,y)
  }
  BF = 1/t2_BF10(x,n1,r,model,location,scale,dff,hypothesis )
  x = x[round(BF,1)== round(D,1)]
  
  return(x)
}

# p(BF01>D|H0)
t2_TNE <- function(t , n1,r,model,location ,scale,dff , hypothesis){
  n2 = n1*r
  df = n1+n2-2

  if (any(t == "bound cannot be found")){
    return(t)
  }
  if (length(t)>1){
    pro = pt(max(t),df) - pt(min(t),df)
  }else{
    if (t >0){
      pro = pt(t,df)
    } else {
      pro = 1-pt(t,df)
    }
  }
  return(pro)
}

# p(BF10>D|H1)
t2_TPE <-function(t,n1,r,model ,location ,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t =="bound cannot be found" | length(t)==0)){
    t=0
    return(t)
  }
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

  normalization  <- integrate(function(delta)t2_prior(delta, location,scale,dff,model),lower = bound[1],upper = bound[2])$value
  
  if (length(t) > 1){
    
    int  <- function(delta){
      pro1 = pnct(t[t>0],df,ncp = delta *constant,lower  = F)
      pro2 = pnct(t[t<0],df,ncp = delta *constant,lower  = T)
      (pro1+pro2)* t2_prior(delta, location,scale,dff,model)/normalization }
    
    
    
  }  else if (t >= 0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = F)
      (pro1)* t2_prior(delta, location,scale,dff,model)/normalization }
    
  } else if (t < 0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = T)
      (pro1)* t2_prior(delta, location,scale,dff,model)/normalization }
    
  }
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
t2_FNE<-function(t,n1,r,model ,location ,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t =="bound cannot be found" | length(t)==0)){
    t=0
    return(t)
  }
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
  t <- as.numeric(t)
  if (length(t) == 4) {  # Corrected condition check
    t = t[2:3]
  }
  
  normalization  <- integrate(function(delta)t2_prior(delta, location,scale,dff,model),lower = bound[1],upper = bound[2])$value
  if (length(t)>1){
    
    int  <- function(delta){
      pro1 = pnct(max(t),df,ncp = delta *constant,lower  = T)
      pro2 = pnct(min(t),df,ncp = delta *constant,lower  = T)
      (pro1-pro2)* t2_prior(delta, location,scale,dff,model)/normalization }
    
  } else if (t>0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = T)
      (pro1)* t2_prior(delta, location,scale,dff,model)/normalization }
    
  } else if (t<0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = F)
      (pro1)* t2_prior(delta, location,scale,dff,model)/normalization }
    
  }
  x = integrate(int,lower = bound[1],upper = bound[2],stop.on.error = FALSE)$value
  
  if (length(t) == 1 ){
    if (t>0 & hypothesis == "<"){
      x = 1-x
    }
    if (t<0 & hypothesis == ">") {
      x = 1-x
    }}
  return(x) 
} 

# p(BF10>D|H0)
t2_FPE <- function(t,n1,r,model ,location ,scale,dff, hypothesis){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t =="bound cannot be found" | length(t)==0)){
    t=0
    return(t)
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

t2_N_finder<-function(D,r,target,model,location,scale,dff, hypothesis ,
                   model_d,location_d,scale_d,dff_d,de_an_prior ,alpha){
 
  lo = 2
  t = t2_BF10_bound(D, lo,r,model ,location ,scale,dff , hypothesis)
  pro = t2_TPE(t , n1=lo,r , model=model , location=location ,scale=scale,dff = dff, hypothesis=hypothesis)
  
  if (pro>target){
    N = 2
    return(N)
  }
  
  up= 100000
  if (de_an_prior == 1){
  Power_root <- function(n1){
    
    t = t2_BF10_bound(D, n1,r,model ,location ,scale,dff , hypothesis)
    pro = t2_TPE(t , n1,r , model=model , location=location ,scale=scale,dff = dff, hypothesis=hypothesis)
    return(target- pro)}
  }else {

    Power_root <- function(n1){
      
      t = t2_BF10_bound(D, n1,r,model ,location ,scale,dff , hypothesis)
      pro = t2_TPE(t , n1,r , model_d , location_d ,scale_d,dff_d, hypothesis=hypothesis)
      return(target- pro)}
    
    
  }
  N = uniroot(Power_root,lower = 2,upper =  up)$root
  t = t2_BF10_bound(D, N,r,model ,location ,scale,dff , hypothesis)
  FPE =t2_FPE(t,N,r,model ,location ,scale,dff, hypothesis)
 
  if (FPE<alpha){
    
    return(N)
  }else{
    alpha_bound <- function(df) {
      t <- t2_BF10_bound(D, df,r,model ,location ,scale,dff , hypothesis)
      pro <- t2_FPE(t,df,r,model ,location ,scale,dff, hypothesis)
      return(pro - alpha)
    }
    NN = uniroot(alpha_bound,lower = N,upper =  up)$root
    return(NN)
  }
  
  }


# probability table
t2_Table <- function(D,r,target,model,location,scale,dff, hypothesis,
                  model_d,location_d,scale_d,dff_d, de_an_prior,N1,N2, mode_bf ,alpha ){

  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  
  if (mode_bf == 1){
    
    n1 = ceiling(t2_N_finder(D,r,target,model,location,scale,dff, hypothesis ,
                     model_d,location_d,scale_d,dff_d, de_an_prior ,alpha))
    n2 = n1*r
  } else {
    n1 = N1
    n2 = N2
    r=n2/n1
  }
  bound01 = t2_BF01_bound(D, n1,r,model,location,scale,dff , hypothesis)
  
  
  if (length(bound01)==0|any(bound01 == "bound cannot be found")){
    FNE = 0
    TNE = 0
  }else{
    
    if (de_an_prior ==1){
      FNE = t2_FNE(bound01,n1,r,model ,location,scale,dff , hypothesis )
      TNE = t2_TNE(bound01,n1,r,model ,location,scale,dff , hypothesis )
    }
    if (de_an_prior ==0){
      FNE = t2_FNE(bound01,n1,r,model_d,location_d,scale_d,dff_d , hypothesis )
      TNE = t2_TNE(bound01,n1,r,model,location ,scale,dff , hypothesis )
    }
    
    
    
  }
  bound10 = t2_BF10_bound(D, n1,r,model,location,scale,dff , hypothesis)
  
  if (de_an_prior ==1){
    TPE = t2_TPE(bound10,n1,r,model ,location,scale,dff , hypothesis )
    FPE = t2_FPE(bound10,n1,r,model ,location,scale,dff , hypothesis )
  }else{
    TPE = t2_TPE(bound10,n1,r,model_d ,location_d,scale_d,dff_d , hypothesis)
    FPE = t2_FPE(bound10,n1,r,model ,location,scale,dff , hypothesis )
  }
  
  
  table <- data.frame(
    TPE = TPE,
    FNE = FNE,
    TNE = TNE,
    FPE = FPE,
    N1 =  n1,
    N2 = n2 )
  return(table)
}


# plots for showing the relationship between BF and t-values 

t2_bf10 <-function(D ,n1,r, target,model ,location ,scale,dff  , hypothesis ){
  
  tt= seq(from = -5,to = 5,.2)
  BF_D = t2_BF10_bound(D, n1,r,model ,location ,scale,dff , hypothesis)
  BF10 = t2_BF10(tt,n1,r,model ,location,scale,dff,hypothesis)
  
  
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
  
  max_BF = 1/t2_BF10(0,n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis ="!=" )
  BF_D = t2_BF01_bound(D , n1,r,model ,location ,scale,dff , hypothesis)

  
  
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

Power_t2<-function(D,model,location,scale,dff, hypothesis,
                   model_d,location_d,scale_d,dff_d, de_an_prior,n1,r){
  Total_ = n1 + n1*r
  smin = 4
  smax = Total_*1.2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  sn1 = sdf/(1+r)
  power =  array(NA, dim = c(length(sdf)))
  
  for ( i in 1:length(sdf)){
    t = t2_BF10_bound(D , sn1[i],r,model ,location ,scale,dff , hypothesis)
    power[i] = switch(de_an_prior,
                      "1" = t2_TPE(t,sn1[i],r,model ,location ,scale,dff , hypothesis ),
                      "0" = t2_TPE(t,sn1[i],r,model_d ,location_d ,scale_d,dff_d , hypothesis ))
    
    
  }
  plot(sdf+1,power,type="l",main = "",frame.plot = FALSE,xlab = "Total sample size", ylab = "Probability of True positive evidence",xlim = c(1,max(sdf)), 
       ylim = c(0,1) )
  
}

