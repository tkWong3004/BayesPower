library("gsl")
###################3

r_auto_uniroot_fixed_lower <- function(f,lower, upper = 1, step = .05, max_attempts = 25) {
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
    upper <- upper - step
  }
  
}

r_auto_uniroot_fixed_upper <- function(f,upper, lower = -1, step = .05, max_attempts = 25, ...) {
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
    lower <- lower + step
  }
  
}


#Fisher
r_mean <-function(r){
  (1/2)*log((1+r)/(1-r))
}

r_sd <-function(N){
  1/sqrt(N-3)
}
#prior
d_strechted_beta <-function(rho,k){
  2^((k-2)/k)*(1-rho^2)^((1-k)/k)/beta(1/k,1/k)
  
}


beta_pdf <- function(rho, alpha, beta) {

  a = -1
  c = 1
  # Beta function
  B_ab <- beta(alpha, beta)
  
  # Compute the PDF
  pdf_value <- ((rho - a)^(alpha - 1) * (c - rho)^(beta - 1)) / ((c - a)^(alpha + beta - 1) * B_ab)
  
  return(pdf_value)
}


# likelihood of non-local prior
dnlp <-function(delta,mu,ta){
  ((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
}

r_prior<- function(rho,k,location,scale,dff,model, alpha, beta){
  
  switch(model,
         "Normal" = dnorm(rho,location,scale),
         "d_beta"   = d_strechted_beta(rho,k),
          "NLP"   = dnlp(rho,location,scale),
          "t_dis" = tstude(rho,location,scale,dff),
         "beta" = beta_pdf (rho, alpha, beta))
}

#########################
d_cor <- function(r, rho, n) {
  n=n-1
  
  # Calculate the logarithmic terms
  log_gamma_n <- lgamma(n)
  log_gamma_n_plus_half <- lgamma(n + 0.5)
  
  # Calculate the logarithmic difference
  log_difference <- log_gamma_n - log_gamma_n_plus_half
  
  # Exponentiate to get the ratio
  ratio <- exp(log_difference)
  
  # Logarithmic version of the rest of the terms
  log_likelihood_value <- log(n - 1) - 0.5 * log(2 * pi) + log(ratio) +
    0.5 * n * log(1 - rho^2) +
    0.5 * (n - 3) * log(1 - r^2) +
    (-n + 0.5) * log(1 - rho * r)  # This term might go to infinity
  
  # Exponentiate the result to return it in original scale
  likelihood_value <- exp(log_likelihood_value) *
    hyperg_2F1(0.5, 0.5, n + 0.5, 0.5 * (r * rho + 1))
  
  return(likelihood_value)
}


r_BF10<-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  x = NA
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  
  normalization  <- integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound[1],upper = bound[2],rel.tol = 1e-10)$value
  
  int  <- function(rho){d_cor(r[i],rho,n)*r_prior(rho,k,location,scale,dff,model, alpha, beta)
  }
  
  for (i in 1:length(r)){
    lh1 =integrate(int,lower = bound[1],upper = bound[2],stop.on.error = F,rel.tol = 1e-10)$value/normalization
    lh0 = d_cor(r[i],h0,n)
    x[i] = lh1/lh0
  }
  return(x)
}




r_BF_bound_10 <-function(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  y <- numeric(0)
  Bound_finding <-function(r){
    r_BF10(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)- D
  }
  
  if (hypothesis=="!="){

      
    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)  
    y <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)

    
  } 
  if (hypothesis == ">"){
    
    x <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
  }
  if (hypothesis == "<"){
    
    x <-r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)  
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

r_BF_bound_01 <-function(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  x <- numeric(0)
  y <- numeric(0)
  Bound_finding <-function(r){
    1/r_BF10(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)- D
  }
  

  if (hypothesis=="!="){
    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25) 
    y <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
    
    
    
  } 
  if (hypothesis == "<"){
    
    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)
  }
  if (hypothesis == ">"){
    
    x <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
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

p_cor<-function(limit,rho,n,lower.tail){
  
  pnorm(r_mean(limit),r_mean(rho),sd = r_sd(n),lower.tail =  lower.tail)
  
  
}

r_TPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){

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
  
  
  
  
  
  
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  normalization  <- integrate(function(rho)r_prior(rho ,k ,location ,scale ,dff ,model , alpha, beta),lower = bound[1],upper = bound[2])$value
  
  if (hypothesis == "!=" ){
    
    int  <- function(rho){
      pro1 = p_cor(max(r),rho,n,lower.tail = F)
      pro2 = p_cor(min(r),rho,n,lower.tail = T)

      (pro1+pro2)* r_prior(rho ,k ,location ,scale ,dff ,model, alpha, beta )/normalization
      }
    
  } else if (hypothesis == ">") { 
    int  <- function(rho){(p_cor(r,rho,n,lower.tail =F))* r_prior(rho ,k ,location ,scale ,dff ,model , alpha, beta)/normalization } 
  } else if (hypothesis == "<") {
    
    int  <- function(rho){(p_cor(r,rho,n,lower.tail =T))* r_prior(rho ,k ,location ,scale ,dff ,model, alpha, beta )/normalization }
    
    
  }
  
  x = integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-8)$value
  return(x)
  
}



r_FNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  
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
  
  
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  
  normalization  <- integrate(function(rho)r_prior(rho ,k ,location ,scale ,dff ,model , alpha, beta),lower = bound[1],upper = bound[2])$value
  
  if (hypothesis == "!=" ){
    
    int  <- function(rho){
      pro1 = p_cor(max(r),rho,n,lower.tail = T)
      pro2 = p_cor(min(r),rho,n,lower.tail = T)
      (pro1-pro2)* r_prior(rho ,k ,location ,scale ,dff ,model , alpha, beta)/normalization}
    
  } else if (hypothesis == ">") { 
    int  <- function(rho){(p_cor(r,rho,n,lower.tail =T))* r_prior(rho ,k ,location ,scale ,dff ,model, alpha, beta )/normalization } 
  } else if (hypothesis == "<") {
    
    int  <- function(rho){(p_cor(r,rho,n,lower.tail =F))* r_prior(rho ,k ,location ,scale ,dff ,model, alpha, beta )/normalization }
    
    
  }
  
  x = integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-8, subdivisions=10000000)$value
  return(x)
  
}



r_FPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  
  if (any(r =="no bound is found" | length(r)==0)){
    r=0
    return(r)
  }
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  
  
  
  if (hypothesis == "!=" ){
    

      pro1 = p_cor(max(r),h0,n,lower.tail = F)
      pro2 = p_cor(min(r),h0,n,lower.tail = T)
      x = pro1+pro2
    
  } else if (hypothesis == ">") { 
   x = p_cor(r,h0,n,lower.tail =F)
  } else if (hypothesis == "<") {
    
    x = p_cor(r,h0,n,lower.tail =T)
    
    
  }
  
  
  return(x)
  
}



r_TNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  
  if (any(r =="no bound is found" | length(r)==0)){
    r=0
    return(r)
  }
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  
  
  
  if (hypothesis == "!=" ){
    
    
    pro1 = p_cor(max(r),h0,n,lower.tail = T)
    pro2 = p_cor(min(r),h0,n,lower.tail = T)
    x = pro1-pro2
    
  } else if (hypothesis == ">") { 
    x = p_cor(r,h0,n,lower.tail =T)
  } else if (hypothesis == "<") {
    
    x = p_cor(r,h0,n,lower.tail =F)
    
    
  }
  
  
  return(x)
  
}



r_N_finder<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                       location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP){

  lo = 15
  r = r_BF_bound_10(D,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  if (de_an_prior == 0 ){
    pro = r_TPE(r,lo,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
  } else {
    pro = r_TPE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  }
  
  if (pro>target){
    b10 = r_BF_bound_10(D,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    
    FPE = r_FPE(b10,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    
    
    if ( FPE<FP){
      return(lo)
    } else{
      alpha_bound <- function(n) {
        b10 <- r_BF_bound_10(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
        pro <- r_FPE(b10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
        return(pro - FP)
        
      }
      NN = uniroot(alpha_bound,lower = lo,upper = 100000)$root
      #auto_uniroot_fixed_lower (alpha_bound,N, upper = 1000, step = 1000, max_attempts = 100)
      
      
      #uniroot(alpha_bound,lower = N,upper = 10000)$root
    }
    

  }
  
  up= 5000
  Power_root <- function(N){
    r = numeric(0)
    r = r_BF_bound_10(D,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    
    if (de_an_prior == 0 ){
      pro = r_TPE(r,N,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
    } else {
      pro = r_TPE(r,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    }
    #pro = switch(de_an_prior,
    #             "1" = r_TPE(r,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model),
     #            "0" = r_TPE(r,N,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
     #            )
    
                 return(pro-target)
    }

  N = auto_uniroot_fixed_lower(Power_root ,lower=lo,upper = up,step = 50,max_attempts = 100)
  
  
  b10 = r_BF_bound_10(D,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  
  FPE = r_FPE(b10,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  
  if ( FPE<FP){
    return(N)
  } else{
    alpha_bound <- function(n) {
      b10 <- r_BF_bound_10(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
      pro <- r_FPE(b10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
      return(pro - alpha)
      
    }
    NN = uniroot(alpha_bound,lower = N,upper = 100000)$root
      #auto_uniroot_fixed_lower (alpha_bound,N, upper = 1000, step = 1000, max_attempts = 100)
      
      return(NN)
      #uniroot(alpha_bound,lower = N,upper = 10000)$root
  }

  }


r_table<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                    location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,FP ){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  
  if (mode_bf == 1){
    
    n = ceiling(r_N_finder(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                           location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP))
  } else {
   n = N
 }
  
  max_BF = 1/r_BF10(h0,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  if (max_BF<D){
    TPE = 0
    FPE = 0
  }
  
  bound01 = r_BF_bound_01(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)

  
  if (length(bound01)==0){
    FNE = 0
    TNE = 0
  }else{
    
    if (de_an_prior ==1){
      FNE = r_FNE(bound01,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
      TNE = r_TNE(bound01,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    }
    if (de_an_prior ==0){
      FNE = r_FNE(bound01,n,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
      TNE = r_TNE(bound01,n,, alpha, betak,h0,hypothesis,location,scale,dff,model)
    }
    
  
  
  }
  
  bound10 = r_BF_bound_10(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  
  if (de_an_prior ==1){
  TPE = r_TPE(bound10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  FPE = r_FPE(bound10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  }else{
    TPE = r_TPE(bound10,n,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
    FPE = r_FPE(bound10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    
 
  }
  
  
  table <- data.frame(
    TPE = TPE,
    FNE = FNE,
    TNE = TNE,
    FPE = FPE,
    N =  n)
  return(table)
  
}




r_prior_plot <-function(k, alpha, beta,h0,location,scale,dff,model,de_an_prior,
                        k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d,hypothesis){
  par(mfrow = c(1, 1))
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  rho= seq(bound[1],bound[2],.01)
  

  normalization  <- integrate(function(rho)r_prior(rho,k,location,scale,dff,model, alpha, beta),lower = bound[1],upper = bound[2])$value
  
  
  prior_rho = NA
  prior_rho = r_prior(rho,k,location,scale,dff,model, alpha, beta)/normalization
  

    
  plot(rho,prior_rho,xlab= bquote(bold(rho)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~rho~" under the alternative hypothesis")),frame.plot = FALSE)
  if (de_an_prior ==0){
    
    if (model_d != "Point"){
      
      normalization_d  <- integrate(function(rho)r_prior(rho,k_d,location_d,scale_d,dff_d,model_d, alpha, beta),lower = bound[1],upper = bound[2])$value
      
      prior_rho_D = r_prior(rho,k_d,location_d,scale_d,dff_d,model_d, alpha_d, beta_d)/normalization_d
      clean_vec <- prior_rho_D[is.finite(prior_rho_D) & !is.na(prior_rho_D)]
      plot(rho,prior_rho,ylim =c(0,max(clean_vec)),xlab= bquote(bold(rho)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~rho~" under the alternative hypothesis")),frame.plot = FALSE)
      
      lines(rho,prior_rho_D ,lty = 2)
      
    } else {
      plot(rho,prior_rho,xlab= bquote(bold(rho)),ylim = c(0,max(prior_rho)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~rho~" under the alternative hypothesis")),frame.plot = FALSE)
      arrows(x0=location_d, y0=0, x1=location_d, y1=max(prior_rho), length=0.2, code=2, col="black", lwd=1,,lty = 2)
      
      
    }

    
    legend("topright", 
           legend = c("Analysis prior", "Design prior"), 
           lty = c(1, 2), 
           col = c("black", "black"),
           bty = "n") 
    
  }
  
  
}


r_bf10_p <-function(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model){
  
  rr= seq(from = -.99,to = .99,.01)
  BF_D = r_BF_bound_10(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model)
  BF10 = r_BF10(rr,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model)
  
  
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
  
  max_BF = 1/ r_BF10(h0,n,k,alpha, beta,h0,hypothesis="!=",location,scale,dff,model)
  BF_D = r_BF_bound_01(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model)
  
  
  
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

Power_r<-function(D,k, alpha, beta,h0,hypothesis,location,scale,dff,model, 
                  k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d, de_an_prior,N){
  
  smin = 4
  smax = N*1.2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  power =  array(NA, dim = c(length(sdf)))
  
  for ( i in 1:length(sdf)){
    r = r_BF_bound_10(D,sdf[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    power[i] = switch(de_an_prior,
                      "1" = r_TPE(r,sdf[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model),
                      "0" = r_TPE(r,sdf[i],k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d))
    
   
    
    
  }
  plot(sdf+1,power,type="l",main = "",frame.plot = FALSE,xlab = "Sample size", ylab = "Probability of True positive evidence",xlim = c(1,max(sdf)), 
       ylim = c(0,1) )
  
}

