#### error prevention functions
# sometimes, when the bound is too large,the uniroot function cannot find the root.
# so the upper bound needed to be extended incrementally 
auto_uniroot_fixed_lower <- function(f,lower, upper = 1000, step = 1000, max_attempts = 10, ...) {
  attempts <- 0
  while (attempts < max_attempts) {
    attempts <- attempts + 1
    # Try to find the root with the current bounds
    result <- tryCatch({
      uniroot(f, lower = lower, upper = upper)$root
    }, error = function(e) {
      # If there's an error, return NA to indicate no root found
      return(NA)
    })
    
    # If a root is found (not NA), return the result
    if (!is.na(result)) {
      return(result)
    }
    
    # If no root is found, expand the search range and try again
    upper <- upper +step
  }
  
}

auto_uniroot_fixed_upper <- function(f,upper, lower = 1000, step = 1000, max_attempts = 10, ...) {
  attempts <- 0
  while (attempts < max_attempts) {
    attempts <- attempts + 1
    # Try to find the root with the current bounds
    result <- tryCatch({
      uniroot(f, lower = lower, upper = upper)$root
    }, error = function(e) {
      # If there's an error, return NA to indicate no root found
      return(NA)
    })
    
    # If a root is found (not NA), return the result
    if (!is.na(result)) {
      return(result)
    }
    
    # If no root is found, expand the search range and try again
    lower <- lower - step
  }
  
}






te_prior<- function(delta,scale,dff,model){
  
  switch(model,
         "Cauchy" = tstude(delta,0, scale,1),
         "Normal"   = dnorm(delta,0,scale),
         "NLP"   = dnlp(delta,0,scale),
         "t-distribution" = tstude(delta,0,scale,dff))
  
  
}

t1e_BF10i <-function(t,df,model ,scale,dff , hypothesis,e ){
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
  normalizationh0<- switch(hypothesis,
                           "!=" = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           "<"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           ">"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value)

  
  int  <- function(delta){
    dt(t,df,ncp = delta *sqrt(df+1))* te_prior(delta,scale,dff,model)/normalizationh1}
  
  error = 1e-8
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
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
  

  x <- numeric(length(t))
  
  for(i in 1:length(t)){
    
    x[i]= t1e_BF10i(t[i],df,model ,scale,dff , hypothesis,e )
    
  }
  return(x)
}

# finding the t that correspond to BF10=D
t1e_BF01_bound <-function(D, df,model,scale,dff , hypothesis,e){
  y <- numeric(0)
  Bound_finding <-function(t){
    1/t1e_BF10(t,df,model,scale,dff , hypothesis,e )- D
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
  BF = 1/t1e_BF10(x,df,model,scale,dff , hypothesis,e )
  x = x[round(BF,2)==round(D,2)]
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  return(x)
}


t1e_BF10_bound <-function(D, df,model,scale,dff , hypothesis,e){
  y <- numeric(0)
  Bound_finding <-function(t){
    t1e_BF10(t,df,model,scale,dff , hypothesis,e )- D
  }
  
  if (hypothesis=="!="){
    x <- tryCatch({
      uniroot(Bound_finding, lower = -15, upper = 0)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    y <- tryCatch({
      uniroot(Bound_finding, lower = 0, upper = 15)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
  } else {
    x <- tryCatch({
      uniroot.all(Bound_finding, lower = -15, upper = 15)
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
  BF = t1e_BF10(x,df,model,scale,dff , hypothesis,e )
  x = x[round(BF,2)==round(D,2)]
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  return(x)
}



t1e_TPE <-function(t,df,model ,scale,dff , hypothesis ,e,location){
  
  if (any(t =="no bound is found" | length(t)==0)){
    t=0
    return(t)
  }
  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {pnct(min(t),df,ncp= location*sqrt(df+1),lower = T)+ pnct(max(t),df,ncp=location*sqrt(df+1),lower = F)},
               "<"  = {pnct(t,df,ncp = location *sqrt(df+1),lower  = T)},
               ">"  = {pnct(t,df,ncp = location *sqrt(df+1),lower  = F)}
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
      pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = F)
      pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1+pro2)* te_prior(delta,scale,dff,model)/normalizationh1 }
    
    
    
  }  else if (hypothesis  == ">") { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh1}
    
  } else if (hypothesis  == "<"){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
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

t1e_FNE <-function(t,df,model ,scale,dff , hypothesis ,e,location){
  
  if (any(t =="no bound is found" | length(t)==0)){
    t=0
    return(t)
  }
  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {pnct(max(t),df,ncp= location*sqrt(df+1),lower = T)- pnct(min(t),df,ncp=location*sqrt(df+1),lower = T)},
               "<"  = {pnct(t,df,ncp = location *sqrt(df+1),lower  = F)},
               ">"  = {pnct(t,df,ncp = location *sqrt(df+1),lower  = T)}
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
      pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = T)
      pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1-pro2)* te_prior(delta,scale,dff,model)/normalizationh1 }
    
    
    
  }  else if (t >= 0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh1}
    
  } else if (t < 0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
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

t1e_TNE <-function(t,df,model ,scale,dff , hypothesis ,e){
  
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
      pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = T)
      pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1-pro2)* te_prior(delta,scale,dff,model)/normalizationh0 }
    
    
    
  }  else if (t >= 0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh0}
    
  } else if (t < 0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh0 }
    
  }
  error = 1e-8
  
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
  

    x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  1-x
  return(x) 
  
} 

t1e_FPE <-function(t,df,model ,scale,dff , hypothesis ,e){
  
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
      pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = F)
      pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
      (pro1+pro2)* te_prior(delta,scale,dff,model)/normalizationh0 }
    
    
    
  }  else if (t >= 0) { 
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
      (pro1)* te_prior(delta,scale,dff,model)/normalizationh0}
    
  } else if (t < 0){
    int  <- function(delta){
      pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
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

t1e_N_finder<-function(D,target,model,scale,dff, hypothesis,e ,
                   model_d,scale_d,dff_d, de_an_prior,location_d  ,alpha){
  
  lo = 2
  t = t1e_BF10_bound(D, lo,model,scale,dff , hypothesis,e)
  pro =t1e_TPE (t,lo,model ,scale,dff , hypothesis,e ,location_d)

  if (pro>target){
    N = 2
    return(N)
  }
  


  if (de_an_prior == 1){
    Power_root <- function(df){
      
      t = t1e_BF10_bound(D, df,model,scale,dff , hypothesis,e)
      pro = t1e_TPE (t,df,model ,scale,dff , hypothesis,e )
      return(target- pro)}
  }else {
    
    Power_root <- function(df){
      
      t = t1e_BF10_bound(D, df,model,scale,dff , hypothesis,e)
      pro = t1e_TPE (t,df,model_d ,scale_d,dff_d , hypothesis,e ,location_d)
      return(target- pro)}
    
    
  }
  
  N = auto_uniroot_fixed_lower(Power_root,lower = 2)
    
  t = t1e_BF10_bound(D,N,model,scale,dff,hypothesis ,e )
  
  FPE =t1e_FPE(t,N,model ,scale,dff , hypothesis ,e)
  
  
  
  if (FPE<alpha){
    
    return(N)
  }else{
    alpha_bound <- function(df) {
      t <- t1e_BF10_bound(D,df,model,scale,dff,hypothesis ,e )
      pro <- t1e_FPE(t , df , model=model , scale=scale,dff=dff, hypothesis,e)
      return(pro - alpha)
    }
    NN = auto_uniroot_fixed_lower(alpha_bound,lower = N)
    return(NN)
  }
        }
    
    
  

  
  
  

  
t1e_table<-function(D,target,model,scale,dff, hypothesis,e ,
                    model_d,scale_d,dff_d, de_an_prior,df,mode_bf,location_d ,alpha ){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  
  if (mode_bf == 1){
    
    df = t1e_N_finder(D,target,model,scale,dff, hypothesis,e ,
                      model_d,scale_d,dff_d, de_an_prior ,location_d,alpha )
    df = ceiling(df+1) -1
  } else {
    df = df
  }
  
  max_BF = 1/t1e_BF10i(0,df,model ,scale,dff , hypothesis,e )
  if (max_BF<D){
    TPE = 0
    FPE = 0
  }
  
  bound01 = t1e_BF01_bound(D, df,model,scale,dff , hypothesis,e)
  
  
  if (length(bound01)==0){
    FNE = 0
    TNE = 0
  }else{
    
    if (de_an_prior ==1){
      FNE = t1e_FNE(bound01,df,model ,scale,dff , hypothesis ,e)
      TNE = t1e_TNE(bound01,df,model ,scale,dff , hypothesis ,e)
    }
    if (de_an_prior ==0){
      FNE = t1e_FNE(bound01,df,model_d,scale_d,dff_d , hypothesis ,e,location_d)
      TNE = t1e_TNE(bound01,df,model ,scale,dff , hypothesis ,e)
    }
    
    
    
  }
  
  bound10 = t1e_BF10_bound(D, df,model,scale,dff , hypothesis,e)
  
  if (de_an_prior ==1){
    TPE = t1e_TPE(bound10,df,model ,scale,dff , hypothesis ,e)
    FPE = t1e_FPE(bound10,df,model ,scale,dff , hypothesis ,e)
  }else{
    TPE = t1e_TPE(bound10,df,model_d ,scale_d,dff_d , hypothesis ,e,location_d)
    FPE = t1e_FPE(bound10,df,model ,scale,dff , hypothesis ,e)
  }
  
  
  table <- data.frame(
    TPE = TPE,
    FNE = FNE,
    TNE = TNE,
    FPE = FPE,
    N =  df)
  return(table)
  
}

t1e_prior_plot <-function(model,scale,dff , hypothesis,e,de_an_prior,model_d,scale_d,dff_d,location_d ){
  par(mfrow = c(1, 1))
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
  normalizationh0<- switch(hypothesis,
                           "!=" = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           "<"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value,
                           ">"  = integrate(function(delta) te_prior(delta, scale,dff,model),lower = bound_h0[1],upper = bound_h0[2])$value)
  
  
  dd = switch(hypothesis,
              "!=" = seq(-3,3,.01),
              ">" = seq(0,3,.01),
              "<" = seq(-3,0,.01)
  )
  prior_h1 = te_prior(dd, scale,dff,model)/normalizationh1
  prior_h0 = te_prior(dd, scale,dff,model)/normalizationh0

  switch(hypothesis,
         "!=" = { prior_h1[dd>min(e)&dd<max(e)]=-1 },
         ">" = { prior_h1[dd<e]=-1 },
         "<" = { prior_h1[dd>e]=-1 }
         
         )
  
  switch(hypothesis,
         "!=" = { prior_h0[!(dd>min(e)&dd<max(e))]=-1},
         ">" = { prior_h0[dd>e]=-1 },
         "<" = { prior_h0[dd<e]=-1 }
         
  )
  
  
  
  plot(dd, prior_h0,
       type = "l",
       lty = 3,
       lwd = 4,
       xlab = bquote(bold(delta)),
       ylim = c(0,max(prior_h0)*1.1),
       ylab = "density",
       main = bquote(bold("prior distribution on " ~ delta ~ " under the two hypotheses")),
       frame.plot = FALSE)
  lines(dd,prior_h1, type="l", lty=1, lwd=4)
  legend("topright", 
         legend = c("H1", "H0"),  # Labels for the lines
         col = c("black", "black"),  # Line colors
         lty = c(1, 3),  # Line types
         lwd = c(4, 4),  # Line widths
         bty = "n",
         title = "Analysis Prior")  # No box around the legend
  
 
  if (de_an_prior ==0){
    if (model_d == "Point"){

      arrows(x0=location_d, y0=0, x1=location_d, y1=max(prior_h0), length=0.2, code=2, col="gray", lwd=4,,lty = 1)
    }else{
    normalizationd<- switch(hypothesis,
                             "!=" = 1-integrate(function(delta) te_prior(delta, scale_d,dff_d,model_d),lower = bound_h1[1],upper = bound_h1[2])$value,
                             "<"  = integrate(function(delta) te_prior(delta, scale_d,dff_d,model_d),lower = bound_h1[1],upper = bound_h1[2])$value,
                             ">"  = integrate(function(delta) te_prior(delta, scale_d,dff_d,model_d),lower = bound_h1[1],upper = bound_h1[2])$value)
    
    prior_d = te_prior(dd, scale_d,dff_d,model_d)/normalizationd
    
    switch(hypothesis,
           "!=" = { prior_d[dd>min(e)&dd<max(e)]=-1 },
           ">" = { prior_d[dd<e]=-1 },
           "<" = { prior_d[dd>e]=-1 }
           
    )
    lines(dd,prior_d, type="l", lty=1, lwd=4,col="gray")
}
    legend("topleft", 
           legend = c("Analysis prior", "Design prior"), 
           lty = c(1, 1), 
           lwd = c(4, 4),
           col = c("black", "gray"),
           bty = "n") 
    
  }
}
te1_BF <-function(D,df,model ,scale,dff , hypothesis ,e){
  tt= seq(from = -10,to = 10,.5)
  BF_D = t1e_BF10_bound(D, df,model,scale,dff , hypothesis,e)
  BF10 = t1e_BF10(tt, df,model,scale,dff , hypothesis,e)
  
  
  if (length(BF_D) == 1){
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D, digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f ",D,BF_D)
  } else {
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f or %.3f ",D,BF_D[1],BF_D[2])
  }
  par(mfrow = c(1, 2))
  plot(tt,log10(BF10),xlim=c(min(tt),max(tt)),xlab= "t-value",type="l", ylab = expression("logarithm of BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
  abline(v = BF_D,)
  axis(1, c(min(tt),max(tt)))
  if (length(BF_D) != 0 ){
    axis(1, round(BF_D,2))}
  
  max_BF = 1/ t1e_BF10(0, df,model,scale,dff , hypothesis,e)
  BF_D = t1e_BF01_bound(D, df,model,scale,dff , hypothesis,e)
  
  
  
  plot(tt,log10(1/BF10),xlab= "t-value",xlim = c(min(tt),max(tt)),type="l",main = "",frame.plot = FALSE,ylab = bquote("logarithm of BF"[0][1]),xaxt = "n")
  axis(1, c(min(tt),max(tt)))
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
Power_t1e<-function(D,model,location,scale,dff, hypothesis,
                   model_d,location_d,scale_d,dff_d, de_an_prior,N,e){
  
  smin = 2
  smax = N*1.2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  power =  array(NA, dim = c(length(sdf)))
  
  for ( i in 1:length(sdf)){
    t = t1e_BF10_bound(D, sdf[i],model,scale,dff , hypothesis,e)
    
    power[i] = switch(de_an_prior,
                      "1" = t1e_TPE(t,sdf[i],model ,scale,dff , hypothesis ,e,location),
                      "0" = t1e_TPE(t,sdf[i],model_d,scale_d,dff_d , hypothesis ,e,location))
    
    
  }
  plot(sdf+1,power,type="l",main = "",frame.plot = FALSE,xlab = "Sample size", ylab = "Probability of True positive evidence",xlim = c(1,max(sdf)), 
       ylim = c(0,1) )
  
}

