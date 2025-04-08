

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

  normalizationh1<- switch(hypothesis,
                           "!=" = 1-integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h1[1],upper = bound_h1[2])$value,
                           "<"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h1[1],upper = bound_h1[2])$value,
                           ">"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h1[1],upper = bound_h1[2])$value)

  int  <- function(delta){
    dt(t,df,ncp = delta *constant)* te_prior(delta,scale,dff,model)/normalizationh1}
  
  error = 1e-8
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
  if (hypothesis == "!="){
  lh1 = integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value 
  }else{
    lh1 = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value
    
  }
  normalizationh0<- switch(hypothesis,
                           "!=" = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           "<"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           ">"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value)
  
  
  
  int  <- function(delta){
    dt(t,df,ncp = delta *constant)* te_prior(delta,scale,dff,model)/normalizationh0}
  
  lh0 = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(lh1/lh0)
}

t2e_BF10 <-function(t,n1,r,model,scale,dff , hypothesis,e ){
  

  x <- numeric(length(t))
  
  for(i in 1:length(t)){
    
    x[i]= t2e_BF10i(t[i],n1,r,model ,scale,dff , hypothesis,e )
    
  }
  return(x)
}

# finding the t that correspond to BF10=D
t2e_BF01_bound <-function(D, n1,r,model,scale,dff , hypothesis,e){
  y <- numeric(0)
  Bound_finding <-function(t){
    1/t2e_BF10(t,n1,r,model,scale,dff , hypothesis,e )- D
  }
  
  if (hypothesis=="!="){
    x <- tryCatch({
      uniroot(Bound_finding, lower = -10, upper = 0)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    y <- tryCatch({
      uniroot(Bound_finding, lower = 0, upper = 10)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
  } else {
    x <- tryCatch({
      uniroot.all(Bound_finding, lower = -10, upper = 10)
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
    
  }
  
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  if (length(y) == 1){
    x = cbind(x,y)
  }
  BF = 1/t2e_BF10(x,n1,r,model,scale,dff , hypothesis,e )
  x = x[round(BF,2)==round(D,2)]
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  return(x)
}


t2e_BF10_bound <-function(D, n1,r,model,scale,dff , hypothesis,e){
  
  y <- numeric(0)
  Bound_finding <-function(t){
    t2e_BF10(t,n1,r,model,scale,dff , hypothesis,e )- D
  }
  
  if (hypothesis=="!="){
    x <- tryCatch({
      uniroot(Bound_finding, lower = -10, upper = 0)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    y <- tryCatch({
      uniroot(Bound_finding, lower = 0, upper = 10)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
  } else {
    x <- tryCatch({
      uniroot.all(Bound_finding, lower = -10, upper = 10)
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
    
  }
  
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  if (length(y) == 1){
    x = cbind(x,y)
  }
  BF = t2e_BF10(x,n1,r,model,scale,dff , hypothesis,e )
  x = x[round(BF,2)==round(D,2)]
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  return(x)
}



t2e_TPE <-function(t,n1,r,model ,scale,dff , hypothesis ,e,location){
   n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t =="no bound is found" | length(t)==0)){
    t=0
    return(t)
  }
  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {pnct(min(t),df,ncp= location*constant,lower = T)+ pnct(max(t),df,ncp=scale,lower = F)},
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
  
  normalizationh1<- switch(hypothesis,
                           "!=" = 1-integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h1[1],upper = bound_h1[2])$value,
                           "<"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h1[1],upper = bound_h1[2])$value,
                           ">"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h1[1],upper = bound_h1[2])$value)
  x = NULL

  if (length(t) > 1){
    
    int  <- function(delta){
      pro1 = pnct(max(t),df,ncp = delta *constant,lower  = F)
      pro2 = pnct(min(t),df,ncp = delta *constant,lower  = T)
      (pro1+pro2)* te_prior(delta,scale,dff,model)/normalizationh1 }
    
    
    
  }  else if (t >= 0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = F)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh1}
    
  } else if (t < 0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = T)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh1 }
    
  }
  if (scale >.3){
    error = .Machine$double.eps^0.25
  }else {
    error = 1e-8
  }
  
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
  
  if (hypothesis == "!="){
    x = integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value 
  }else{
    x = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value
    
  }
  return(x) 
  
} 

t2e_FNE <-function(t,n1,r,model ,scale,dff , hypothesis ,e,location){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  
  if (any(t =="no bound is found" | length(t)==0)){
    t=0
    return(t)
  }
  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {pnct(min(t),df,ncp= location*constant,lower = T)- pnct(max(t),df,ncp=scale,lower = T)},
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

  normalizationh1<- switch(hypothesis,
                           "!=" = 1-integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h1[1],upper = bound_h1[2])$value,
                           "<"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h1[1],upper = bound_h1[2])$value,
                           ">"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h1[1],upper = bound_h1[2])$value)
  x = NULL
  
  if (length(t) > 1){
    
    int  <- function(delta){
      pro1 = pnct(max(t),df,ncp = delta *constant,lower  = T)
      pro2 = pnct(min(t),df,ncp = delta *constant,lower  = T)
      (pro1-pro2)* te_prior(delta,scale,dff,model)/normalizationh1 }
    
    
    
  }  else if (t >= 0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = T)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh1}
    
  } else if (t < 0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = F)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh1 }
    
  }
  if (scale >.3){
    error = .Machine$double.eps^0.25
  }else {
    error = 1e-8
  }
  
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
  
  if (hypothesis == "!="){
    x = integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value 
  }else{
    x = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value
    
  }
  return(x) 
  
} 

t2e_TNE <-function(t,n1,r,model ,scale,dff , hypothesis ,e){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  
  if (any(t =="no bound is found" | length(t)==0)){
    t=0
    return(t)
  }
  
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = -e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )
  
  normalizationh0<- switch(hypothesis,
                           "!=" = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           "<"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           ">"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value)
  
  
  x = NULL
  
  if (length(t) > 1){
    
    int  <- function(delta){
      pro1 = pnct(max(t),df,ncp = delta *constant,lower  = T)
      pro2 = pnct(min(t),df,ncp = delta *constant,lower  = T)
      (pro1-pro2)* te_prior(delta,scale,dff,model)/normalizationh0 }
    
    
    
  }  else if (t >= 0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = T)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh0}
    
  } else if (t < 0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = F)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh0 }
    
  }
  if (scale >.3){
    error = .Machine$double.eps^0.25
  }else {
    error = 1e-8
  }
  
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
  

    x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  1-x
  return(x) 
  
} 

t2e_FPE <-function(t,n1,r,model ,scale,dff , hypothesis ,e){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t =="no bound is found" | length(t)==0)){
    t=0
    return(t)
  }

  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = -e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )
  
  normalizationh0<- switch(hypothesis,
                           "!=" = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           "<"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           ">"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value)

  x = NULL
  
  if (length(t) > 1){
    
    int  <- function(delta){
      pro1 = pnct(max(t),df,ncp = delta *constant,lower  = F)
      pro2 = pnct(min(t),df,ncp = delta *constant,lower  = T)
      (pro1+pro2)* te_prior(delta,scale,dff,model)/normalizationh0 }
    
    
    
  }  else if (t >= 0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = F)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh0}
    
  } else if (t < 0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *constant,lower  = T)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh0 }
    
  }
  if (scale >.3){
    error = .Machine$double.eps^0.25
  }else {
    error = 1e-8
  }
  
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
  
  
  x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  
  return(x) 
  
} 



t2e_N_finder<-function(D,r,target,model,scale,dff, hypothesis,e ,
                   model_d,scale_d,dff_d, de_an_prior,location,alpha ){
  
  lo = 2
  t = t2e_BF10_bound(D, lo,r,model,scale,dff , hypothesis,e)
  pro =t2e_TPE (t,lo,r,model ,scale,dff , hypothesis,e ,location)

  if (pro>target){
    N = 2
    return(N)
  }
  


  if (de_an_prior == 1){
    Power_root <- function(n1){
      
      t = t2e_BF10_bound(D, n1,r,model,scale,dff , hypothesis,e)
      pro = t2e_TPE (t,n1,r,model ,scale,dff , hypothesis,e )
      return(target- pro)}
  }else {
    
    Power_root <- function(n1){
      
      t = t2e_BF10_bound(D, n1,r,model,scale,dff , hypothesis,e)
      pro = t2e_TPE (t,n1,r,model_d ,scale_d,dff_d , hypothesis,e ,location)
      return(target- pro)}
    
    
  }

  N = auto_uniroot_fixed_lower(Power_root,lower = 2)
  
  t = t2e_BF10_bound(D, N,r,model,scale,dff , hypothesis,e)
  FPE =t2e_FPE(t,N,r,model ,scale,dff , hypothesis ,e)
  
  if (FPE<alpha){
    
    return(N)
  }else{
    alpha_bound <- function(df) {
      t <- t2e_BF10_bound(D, df,r,model,scale,dff , hypothesis,e)
      pro <- t2e_FPE(t,df,r,model ,scale,dff , hypothesis ,e)
      return(pro - alpha)
    }
    
    
    NN = uniroot(alpha_bound,lower = N,upper =  up)$root
    return(NN)
  }
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
  
  max_BF = 1/t2e_BF10i(0,n1,r,model ,scale,dff , hypothesis,e )
  if (max_BF<D){
    TPE = 0
    FPE = 0
  }
  
  bound01 = t2e_BF01_bound(D, n1,r,model,scale,dff , hypothesis,e)
  
  
  if (length(bound01)==0){
    FNE = 0
    TNE = 0
  }else{
    
    if (de_an_prior ==1){
      FNE = t2e_FNE(bound01,n1,r,model ,scale,dff , hypothesis ,e)
      TNE = t2e_TNE(bound01,n1,r,model ,scale,dff , hypothesis ,e)
    }
    if (de_an_prior ==0){
      FNE = t2e_FNE(bound01,n1,r,model_d,scale_d,dff_d , hypothesis ,e,location)
      TNE = t2e_TNE(bound01,n1,r,model ,scale,dff , hypothesis ,e)
    }
    
    
    
  }
  
  bound10 = t2e_BF10_bound(D, n1,r,model,scale,dff , hypothesis,e)
  
  if (de_an_prior ==1){
    TPE = t2e_TPE(bound10,n1,r,model ,scale,dff , hypothesis ,e)
    FPE = t2e_FPE(bound10,n1,r,model ,scale,dff , hypothesis ,e)
  }else{
    TPE = t2e_TPE(bound10,n1,r,model_d ,scale_d,dff_d , hypothesis ,e,location)
    FPE = t2e_FPE(bound10,n1,r,model ,scale,dff , hypothesis ,e)
  }
  
  
  table <- data.frame(
    TPE = TPE,
    FNE = FNE,
    TNE = TNE,
    FPE = FPE,
    N1 =  n1,
    N2 = n1*r )
  return(table)
  
}


t2e_BF <-function(D,n1,r,model ,scale,dff , hypothesis ,e){
  
  tt= seq(from = -5,to = 5,.1)
  BF_D = t2e_BF10_bound(D, n1,r,model,scale,dff , hypothesis,e)
  BF10 = t2e_BF10(tt, n1,r,model,scale,dff , hypothesis,e)
  
  
  if (length(BF_D) == 1){
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D, digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f ",D,BF_D)
  } else {
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f or %.3f ",D,BF_D[1],BF_D[2])
  }
  par(mfrow = c(1, 2))
  plot(tt,log10(BF10),xlim=c(-5,5),xlab= "t-value",type="l", ylab = expression("logarithm of BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
  abline(v = BF_D,)
  axis(1, c(-5,5))
  if (length(BF_D) != 0 ){
    axis(1, round(BF_D,2))}
  
  max_BF = 1/ t2e_BF10(0, n1,r,model,scale,dff , hypothesis,e)
  BF_D = t2e_BF01_bound(D, n1,r,model,scale,dff , hypothesis,e)
  
  
  
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


Power_t2e<-function(D,model,location,scale,dff, hypothesis,
                   model_d,location_d,scale_d,dff_d, de_an_prior,n1,r){
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
  plot(sdf+1,power,type="l",main = "",frame.plot = FALSE,xlab = "Total sample size", ylab = "Probability of True positive evidence",xlim = c(1,max(sdf)), 
       ylim = c(0,1) )
  
}
