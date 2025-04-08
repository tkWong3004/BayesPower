library("gsl")
#########################


re_BF10i<-function(r,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
    x = NA
    bound_h1  <- switch(hypothesis,
                        ">" = c(a = h0+e, b = 1),
                        "<" = c(a = -1, b = h0+e),
                        "!=" = c(a = h0+e[1], b = h0+e[2])
    )
    bound_h0  <- switch(hypothesis,
                        ">" = c(a = h0, b = h0+e),
                        "<" = c(a = h0+e, b = h0),
                        "!=" = c(a = h0+e[1], b = h0+e[2])
    )
    
    
    normalizationh1  <- switch(hypothesis,
                               "!=" = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = -1,upper = bound_h1[1],rel.tol = 1e-10)$value+integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[2],upper = 1,rel.tol = 1e-10)$value,
                               ">"  = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                               "<"  = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
    
  
    int  <- function(rho){d_cor(r,rho,n)*r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh1
    }
    
    if (hypothesis == "!="){
      lh1 = integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-10,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-10,stop.on.error = F)$value 
    }else{
      lh1 = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-10,stop.on.error = F)$value
      
    }
    normalizationh0 <- integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h0[1],upper = bound_h0[2])$value
    
    lh0 = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-10,stop.on.error = F)$value
  
    x = lh1/lh0
    
    return(x)
  }

re_BF10<-function(r,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  x = NA
  for (i in 1:length(r)){
    x[i] = re_BF10i(r[i],n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  }
  return(x)
}



re_BF_bound_10 <-function(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  y =x<- numeric(0)
  Bound_finding <-function(r){
    re_BF10(r,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)- D
  }
  
  if (hypothesis=="!="){
    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)  
    y <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
    
    
  } 
  if (hypothesis == ">"){
    
    x <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
  }
  
  if (hypothesis == "<"){
    
    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)  
  }
  
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  if (length(y) == 1){
    x = cbind(x,y)
  }
 # BF = r_BF10(x,n,k,hypothesis,location,scale,dff,model)
 # x = x[round(BF,1)==round(D,1)]
  #if (length(x) ==0){
   # x = "no bound is found"
    #return(x)
#  } 
  return(x)
}

re_BF_bound_01 <-function(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  y =x<- numeric(0)
  Bound_finding <-function(r){
    1/re_BF10(r,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)- D
  }
  
  if (hypothesis=="!="){
    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)  
    y <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
    
    
  } 
  if (hypothesis == ">"){
    
    x <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
  }
  
  if (hypothesis == "<"){
    
    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)  
  }
  
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  if (length(y) == 1){
    x = cbind(x,y)
  }
  # BF = r_BF10(x,n,k,hypothesis,location,scale,dff,model)
  # x = x[round(BF,1)==round(D,1)]
  #if (length(x) ==0){
  # x = "no bound is found"
  #return(x)
  #  } 
  return(x)
}


re_TPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){

  if (any(r =="no bound is found" | length(r)==0)){
    r=0
    return(r)
  }
  
  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {p_cor(max(r),location,n,lower.tail = F)+ p_cor(min(r),location,n,lower.tail = T)},
               ">"  = {p_cor(r,location,n,lower.tail =F)},
               "<"  = {p_cor(r,location,n,lower.tail =T)}
    )
    return(x)
  }
  
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = -1, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh1  <- switch(hypothesis,
                             "!=" = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = -1,upper = bound_h1[1],rel.tol = 1e-10)$value+integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[2],upper = 1,rel.tol = 1e-10)$value,
                             ">"  = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                             "<"  = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
  
  
  if (hypothesis == "!=" ){
    
    int  <- function(rho){
      pro1 = p_cor(max(r),rho,n,lower.tail = F)
      pro2 = p_cor(min(r),rho,n,lower.tail = T)

      (pro1+pro2)* r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh1
      }
    
  } else if (r >= 0) { 
    int  <- function(rho){(p_cor(r,rho,n,lower.tail =F))* r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh1 } 
  } else if (r <= 0) {
    
    int  <- function(rho){(p_cor(r,rho,n,lower.tail =T))* r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh1 }
    
    
  }
  
  if (hypothesis == "!="){
    x = integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-10,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-10,stop.on.error = F)$value 
  }else{
    x = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-10,stop.on.error = F)$value
    
  }
  return(x)
  
}


re_FNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  
  if (any(r =="no bound is found" | length(r)==0)){
    r=0
    return(r)
  }
  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {p_cor(max(r),location,n,lower.tail = T)- p_cor(min(r),location,n,lower.tail = T)},
               ">"  = {p_cor(r,location,n,lower.tail =T)},
               "<"  = {p_cor(r,location,n,lower.tail =F)}
    )
    return(x)
  }
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = -1, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  
  normalizationh1  <- switch(hypothesis,
                             "!=" = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = -1,upper = bound_h1[1],rel.tol = 1e-10)$value+integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[2],upper = 1,rel.tol = 1e-10)$value,
                             ">"  = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                             "<"  = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
  
  
  if (hypothesis == "!=" ){
    
    int  <- function(rho){
      pro1 = p_cor(max(r),rho,n,lower.tail = T)
      pro2 = p_cor(min(r),rho,n,lower.tail = T)
      (pro1-pro2)* r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh1}
    
  } else if (r >= 0) { 
    int  <- function(rho){(p_cor(r,rho,n,lower.tail =T))* r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh1 } 
  } else if (r <= 0) {
    
    int  <- function(rho){(p_cor(r,rho,n,lower.tail =F))* r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh1 }
    
    
  }
  
  if (hypothesis == "!="){
    x = integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-10,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-10,stop.on.error = F)$value 
  }else{
    x = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-10,stop.on.error = F)$value
    
  }
  
  
  return(x)
  
}


re_FPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  
  if (any(r =="no bound is found" | length(r)==0)){
    r=0
    return(r)
  }
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh0 <- integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h0[1],upper = bound_h0[2])$value
  
  
  
  if (hypothesis == "!=" ){
    

    int <-function(rho){
    pro1 = p_cor(max(r),rho,n,lower.tail = F)
    pro2 = p_cor(min(r),rho,n,lower.tail = T)
    (pro1+pro2)* r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh0
    }
    
  } else if (r >= 0) { 
    int <-function(rho){
   pro = p_cor(r,rho,n,lower.tail =F)
   pro*r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh0
    }
   
  } else if (r <= 0) {
    int <-function(rho){
    pro = p_cor(r,rho,n,lower.tail =T)
    pro*r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh0
    }
    
  }
  
  x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-10,stop.on.error = F)$value
  
  
  
  return(x)
  
}

re_TNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  
  if (any(r =="no bound is found" | length(r)==0)){
    r=0
    return(r)
  }
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  
  normalizationh0 <- integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h0[1],upper = bound_h0[2])$value
  
  
  if (hypothesis == "!=" ){
    
    
    int <-function(rho){
      pro1 = p_cor(max(r),rho,n,lower.tail = T)
      pro2 = p_cor(min(r),rho,n,lower.tail = T)
      (pro1-pro2)* r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh0
    }
    
  } else if (r >= 0) { 
    int <-function(rho){
      pro = p_cor(r,rho,n,lower.tail =T)
      pro*r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh0
    }
    
  } else if (r <= 0) {
    int <-function(rho){
      pro = p_cor(r,rho,n,lower.tail =F)
      pro*r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh0
    }
    
  }
  x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-10,stop.on.error = F)$value
  
  
  return(x)
  
}



re_N_finder<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                      location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP,e){

  lo = 20
  r = re_BF_bound_10(D,lo,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  
  if (de_an_prior == 0 ){
    pro = re_TPE(r,lo,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
  } else {
    pro = re_TPE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  }
  
  
  if (pro>target){
    b10 = re_BF_bound_10(D,lo,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    
    FPE = re_FPE(b10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    
    
    if ( FPE<FP){
      return(lo)
    } else{
      alpha_bound <- function(n) {
        b10 <-  re_BF_bound_10(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
        pro <- re_FPE(b10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
        return(pro - FP)
        
      }
      NN = uniroot(alpha_bound,lower = lo,upper = 100000)$root
      #auto_uniroot_fixed_lower (alpha_bound,N, upper = 1000, step = 1000, max_attempts = 100)
      
      
      #uniroot(alpha_bound,lower = N,upper = 10000)$root
    }
    
    
  }
  
  up= 5000
  Power_root <- function(N){

    r = re_BF_bound_10(D,N,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)

    if (de_an_prior == 0 ){
      pro = re_TPE(r,N,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    }else {
      pro = re_TPE(r,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
        
      }
return(pro-target)
    }

  
  N = auto_uniroot_fixed_lower(Power_root ,lower=lo,upper = up,step = 50,max_attempts = 100)
  b10 = re_BF_bound_10(D,N,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  FPE = re_FPE(b10,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  
  if ( FPE<FP){
    return(N)
  } else{
    alpha_bound <- function(n) {
      b10 <-  re_BF_bound_10(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
      pro <- re_FPE(b10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
      return(pro - FP)
      
    }
    NN = uniroot(alpha_bound,lower = N,upper = 100000)$root
    return(NN)
    #auto_uniroot_fixed_lower (alpha_bound,N, upper = 1000, step = 1000, max_attempts = 100)
    
    
    #uniroot(alpha_bound,lower = N,upper = 10000)$root
  }
}


re_table<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                   location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,FP ,e){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  
  if (mode_bf == 1){
    
    n = ceiling (re_N_finder(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                    location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP,e))
    
  } else {
   n=N
 }
  
  max_BF = 1/re_BF10(h0,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  if (max_BF<D){
    TPE = 0
    FPE = 0
  }
  
  bound01 = re_BF_bound_01(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)

  
  if (length(bound01)==0){
    FNE = 0
    TNE = 0
  }else{
    
    if (de_an_prior ==1){
      FNE = re_FNE(bound01,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
      TNE = re_TNE(bound01,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    }
    if (de_an_prior ==0){
      FNE = re_FNE(bound01,n,k, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
      TNE = re_TNE(bound01,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    }
    
  
  
  }
  
  bound10 = re_BF_bound_10(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  
  if (de_an_prior ==1){
  TPE = re_TPE(bound10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  FPE = re_FPE(bound10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  }else{
    TPE = re_TPE(bound01,n,k, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    FPE = re_FPE(bound10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  }
  
  
  table <- data.frame(
    TPE = TPE,
    FNE = FNE,
    TNE = TNE,
    FPE = FPE,
    N =  ceiling(n))
  return(table)
  
}



re_prior_plot <-function(k, alpha, beta,h0,location,scale,dff,model,de_an_prior,
                         k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d,hypothesis,e){
  par(mfrow = c(1, 1))
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = -1, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )

  normalizationh1  <- switch(hypothesis,
                             "!=" = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = -1,upper = bound_h1[1],rel.tol = 1e-10)$value+integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[2],upper = 1,rel.tol = 1e-10)$value,
                             ">"  = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                             "<"  = integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
  
  normalizationh0 <- integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound_h0[1],upper = bound_h0[2])$value
  
  rho = switch(hypothesis,
              "!=" = seq(-1,1,.002),
              ">" = seq(h0,1,.002),
              "<" = seq(-1,h0,.002)
  )
  
  
  
prior_h1 = r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh1
 prior_h0 = r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalizationh0
  switch(hypothesis,
         "!=" = { prior_h1[rho>=min(bound_h1)&rho<=max(bound_h1)]=-1 },
         ">" = { prior_h1[rho<min(bound_h1)]=-1 },
         "<" = { prior_h1[rho>max(bound_h1)]=-1 }
         
  )
  
  switch(hypothesis,
         "!=" = { prior_h0[!(rho>min(bound_h0)&rho<max(bound_h0))]=-1},
         ">" = { prior_h0[rho>max(bound_h0)]=-1 },
         "<" = { prior_h0[rho<min(bound_h0)]=-1 }
         
  )
  
  plot(rho,prior_h0,
       type = "l",
       lty = 3,
       lwd = 4,
       xlab = bquote(bold(rho)),
       ylim = c(0,max(prior_h0)*1.1),
       ylab = "density",
       main = bquote(bold("prior distribution on " ~ rho ~ " under the two hypotheses")),
       frame.plot = FALSE)
  lines(rho,prior_h1, type="l", lty=1, lwd=4)
  legend("topright", 
         legend = c("H1", "H0"),  # Labels for the lines
         col = c("black", "black"),  # Line colors
         lty = c(1, 3),  # Line types
         lwd = c(4, 4),  # Line widths
         bty = "n",
         title = "Analysis Prior")  # No box around the legend
  
  if (de_an_prior ==0){
    if (model_d == "Point"){

      arrows(x0=h0, y0=0, x1=h0, y1=max(prior_h1,prior_h0), length=0.2, code=2, col="gray", lwd=4,lty = 2)
    }else{
      normalizationd  <- switch(hypothesis,
                                 "!=" = 1-integrate(function(rho)r_prior(rho,k_d,d_location,d_scale,d_dff,model_d ,alpha_d, beta_d),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                                 ">"  = integrate(function(rho)r_prior(rho,k_d,d_location,d_scale,d_dff,model_d ,alpha_d, beta_d),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                                 "<"  = integrate(function(rho)r_prior(rho,k_d,d_location,d_scale,d_dff,model_d ,alpha_d, beta_d),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
    
    prior_rho_D = r_prior(rho,k_d,d_location,d_scale,d_dff,model_d ,alpha_d, beta_d)/normalizationd
    switch(hypothesis,
           "!=" = { prior_rho_D[rho>=min(bound_h1)&rho<=max(bound_h1)]=-1 },
           ">" = { prior_rho_D[rho<min(bound_h1)]=-1 },
           "<" = { prior_rho_D[rho>max(bound_h1)]=-1 }
           
    )

    lines(rho,prior_rho_D, type="l", lty=1, lwd=4,col="gray")
    
    
    }
    legend("topleft", 
           legend = c("Analysis prior", "Design prior"), 
           lty = c(1, 1), 
           lwd = c(4, 4),
           col = c("black", "gray"),
           bty = "n") 
  
  }
}




re_bf10_p <-function(D,n,k,h0,hypothesis,location,scale,dff,model,e){
  
  rr= seq(from = -.99,to = .99,.01)
  BF_D = re_BF_bound_10(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  BF10 = re_BF10(rr,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  
  
  if (length(BF_D) == 1){
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when r = "~.(format(BF_D, digits = 4))))
    #sprintf("BF10 = %.0f when r = %.3f ",D,BF_D)
  } else {
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when r = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
    #sprintf("BF10 = %.0f when r = %.3f or %.3f ",D,BF_D[1],BF_D[2])
  }
  par(mfrow = c(1, 2))
  plot(rr,log10(BF10),xlim=c(-1,1),xlab= "Correlation",type="l", ylab = expression("logarithm of BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
  abline(v = BF_D,)
  axis(1, c(-1,1))
  if (length(BF_D) != 0 ){
    axis(1, round(BF_D,2))}
  if (hypothesis =="!="){
  max_BF = 1/ re_BF10(h0,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  }else{
    max_BF =Inf
    }
  BF_D = re_BF_bound_01(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  
  
  
  plot(rr,log10(1/BF10),xlab= "Correlation",type="l",main = "",frame.plot = FALSE,ylab = bquote("logarithm of BF"[0][1]),xaxt = "n")
  axis(1, c(-1,1))
  if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found" ) ) {
    main = bquote(bold("It is impossible to have BF"[0][1]~"="~.(D)))
    title(main = main)
    #sprintf("It is impossible to have BF01 = %.3f ",D)
  } else      {
    abline(v = BF_D)
    axis(1, round(BF_D,2))
    if (length(BF_D) == 1){
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when r = "~.(format(BF_D, digits = 4))))
      title(main = main)
    } else {
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when r = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
      title(main = main)
    }}

}


Power_re<-function(D,k, alpha, beta,h0,hypothesis,location,scale,dff,model, 
                  k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d, de_an_prior,N,e){
  par(mfrow = c(1, 1))
  smin = 4
  smax = N*1.2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  power =  array(NA, dim = c(length(sdf)))
  
  for ( i in 1:length(sdf)){
    r = re_BF_bound_10(D,sdf[i],k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    
    
    if (de_an_prior == 1){
      power[i] = re_TPE(r,sdf[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
      
    } else {
      power[i] =re_TPE(r,sdf[i],k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    }

    
  }
  plot(sdf+1,power,type="l",main = "",frame.plot = FALSE,xlab = "Sample size", ylab = "Probability of True positive evidence",xlim = c(1,max(sdf)), 
       ylim = c(0,1) )
  
}




