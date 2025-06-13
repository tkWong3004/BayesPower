r_auto_uniroot_fixed_lower <- function(f,lower, upper = 1, step = .05, max_attempts = 25) {
  attempts <- 0
  while (attempts < max_attempts) {
    attempts <- attempts + 1
    # Try to find the root with the current bounds
    result <- tryCatch({
      uniroot(f, lower = lower, upper = upper , tol = 1e-10)$root
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
      uniroot(f, lower = lower, upper = upper, tol = 1e-10)$root
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
    
    normalizationh1 <- switch(hypothesis,
                              "!=" = switch(model,
                                            "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                            "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                            "NLP"          = {
                                              (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                                (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                            }),
                              
                              "<"  = switch(model,
                                            "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                            "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                            "NLP"          = {
                                              (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                            }),
                              ">"  = switch(model,
                                            "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                            "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                            "NLP"          = {
                                              (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                            })
    )
    normalizationh0 <- switch(model,
                              "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                              "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                              "NLP"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                              }
    )
  
    int  <- function(rho){d_cor(r,rho,n)*r_prior(rho,k,location,scale,dff,model, alpha, beta,-1,1)/normalizationh1
    }
    
    if (hypothesis == "!="){
      lh1 = integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-5,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-5,stop.on.error = F)$value 
    }else{
      lh1 = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-5,stop.on.error = F)$value
      
    }
    lh0 = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value
  
    x = (lh1/normalizationh1)/(lh0/normalizationh0)
    return(x)
  }
re_BF10<-function(r,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
   sapply(r, re_BF10i, n = n, k = k, alpha = alpha, beta = beta,
              h0 = h0, hypothesis = hypothesis, location = location,
              scale = scale, dff = dff, model = model, e = e)
}



re_BF_bound_10 <-function(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  y <- numeric(0)
  Bound_finding <-function(r){
    re_BF10(r,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)- D
  }
  opt_result <- optimize(Bound_finding, interval = c(-1, 1))$minimum

  if (hypothesis=="!="){
    x <- r_auto_uniroot_fixed_upper (Bound_finding,opt_result, lower = -1, step = .05, max_attempts = 25)  
    y <- r_auto_uniroot_fixed_lower(Bound_finding,opt_result, upper = 1, step = .05, max_attempts = 25)
  } 
  if (hypothesis == ">"){
    
    x <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
  }
  
  if (hypothesis == "<"){
    
    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)  
  }
  results <- c(x, y)
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")
  
  BF.vals  <- re_BF10(results,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

re_BF_bound_01 <-function(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  re_BF_bound_10 (1/D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
}

re_TPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){

  if (any(r =="bound cannot be found" | length(r)==0)){
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
  normalizationh1 <- switch(hypothesis,
                            "!=" = switch(model,
                                          "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "NLP"          = {
                                            (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                              (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                          }),
                            
                            "<"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "NLP"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          }),
                            ">"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "NLP"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          })
  )
  
  int <- function(rho) {
    pro <- switch(hypothesis,
                  "!=" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  ">"  = p_cor(r, rho, n, lower.tail = FALSE),
                  "<"  = p_cor(r, rho, n, lower.tail = TRUE)
    )
    
    pro * r_prior(rho, k, location, scale, dff, model, alpha, beta,-1,1) / normalizationh1
  }
  

  
  if (hypothesis == "!="){
    x = integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-5,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-5,stop.on.error = F)$value 
  }else{
    x = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-5,stop.on.error = F)$value
    
  }
  return(x)
  
}

re_FNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  
  if (any(r =="bound cannot be found" | length(r)==0)){
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
  
  normalizationh1 <- switch(hypothesis,
                            "!=" = switch(model,
                                          "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "NLP"          = {
                                            (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                              (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                          }),
                            
                            "<"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "NLP"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          }),
                            ">"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "NLP"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          })
  )
  
  int <- function(rho) {
    pro <- switch(hypothesis,
                  "!=" = p_cor(max(r), rho, n, lower.tail = T) -
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  ">"  = p_cor(r, rho, n, lower.tail = T),
                  "<"  = p_cor(r, rho, n, lower.tail = F)
    )
    
    pro * r_prior(rho, k, location, scale, dff, model, alpha, beta,-1,1) / normalizationh1
  }
  
  
  if (hypothesis == "!="){
    x = integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-10,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-10,stop.on.error = F)$value 
  }else{
    x = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-10,stop.on.error = F)$value
    
  }

  return(x)
  
}

re_FPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  
  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh0 <- switch(model,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "NLP"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )
  
  
  int <- function(rho) {
    pro <- switch(hypothesis,
                  "!=" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                         p_cor(min(r), rho, n, lower.tail = TRUE),
                  ">"  = p_cor(r, rho, n, lower.tail = FALSE),
                  "<"  = p_cor(r, rho, n, lower.tail = TRUE)
    )
    
    pro * r_prior(rho, k, location, scale, dff, model, alpha, beta,-1,1) / normalizationh0
  }
  
  x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value
  
  
  return(x)
  
}

re_TNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  
  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  
  normalizationh0 <- switch(model,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "NLP"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )
  
  int <- function(rho) {
    pro <- switch(hypothesis,
                  "!=" = p_cor(max(r), rho, n, lower.tail = TRUE) -
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  ">"  = p_cor(r, rho, n, lower.tail = TRUE),
                  "<"  = p_cor(r, rho, n, lower.tail = FALSE)
    )
    
    pro * r_prior(rho, k, location, scale, dff, model, alpha, beta,-1,1) / normalizationh0
  }
  
  x = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value
  
  
  return(x)
  
}


re_N_finder<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                      location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP,e){
  lo = 10
  upper = 5000
  
  r = re_BF_bound_10(D,lo,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  TPE_lo <- if (de_an_prior == 1)
    re_TPE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e) else
      re_TPE(r,lo,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
  FPE_lo <-  re_FPE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  
  if (TPE_lo > target && FPE_lo < FP) {
    return(lo)
  } else if (TPE_lo > target) {
    alpha.root <- function(n) {
      r <- re_BF_bound_10(D, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model, e)
      re_FPE(r, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model, e) - FP
    }
    return(uniroot(alpha.root, lower = lo, upper = upper)$root)
  }
  
  
  
  Power_root <- function(N){

    r = re_BF_bound_10(D,N,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)

    if (de_an_prior == 0 ){
      pro = re_TPE(r,N,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    }else {
      pro = re_TPE(r,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
        
      }
return(pro-target)
    }

  N.power <- robust_uniroot(Power_root, lower = lo)
  r = re_BF_bound_10(D, N.power,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  FPE = re_FPE(r, N.power,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  if (FPE <= FP) return(N.power)
  
  alpha.root <- function(n) {
    r <- re_BF_bound_10(D, n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    re_FPE(r, n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)-FP
  }
  N.alpha = uniroot(alpha.root,lower = N.power,upper = upper)$root
  return(N.alpha)
}


re_table<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                   location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,FP ,e){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  
  if (mode_bf == 1) n = ceiling (re_N_finder(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                                             location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP,e)) else  n = N
  
  
  # r bounds:
  r10 = re_BF_bound_10(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  r01 = re_BF_bound_01(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  
  # max BF10 possible:
  max_BF <- 1 / re_BF10(h0,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  BF_D   <- r10
  
  # FPE and TPE:
  FPE       <- re_FPE(r10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)

  if (de_an_prior == 1) {
      TPE           <- re_TPE(r10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
      TPR_k         <- k
      TPR_alpha     <- alpha
      TPR_beta      <- beta
      TPR_location  <- location
      TPR_scale     <- scale
      TPR_dff       <- dff
      TPR_model     <- model
        
  } else {
    TPE           <- re_TPE(r10,n,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    TPR_k         <- k_d
    TPR_alpha     <- alpha_d
    TPR_beta      <- beta_d
    TPR_location  <- location_d
    TPR_scale     <- scale_d
    TPR_dff       <- dff_d
    TPR_model     <- model_d
  } 
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- re_FNE(r01,n,TPR_k, TPR_alpha, TPR_beta,h0,hypothesis,TPR_location,TPR_scale,TPR_dff,TPR_model,e)
    TNE <- re_TNE(r01,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  }
  
  # table:
  tab.names <- c(
    sprintf("p(BF10 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H0)", D),
    sprintf("p(BF10 > %0.f | H0)", D),
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}


compute.prior.density.re.h1 <- function(rho,h0, k,location,scale,dff,model, alpha, beta,hypothesis,e) {
  if (model == "Point") return(rep(NA, length(rho)))
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = -1, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh1 <- switch(hypothesis,
                            "!=" = switch(model,
                                          "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "NLP"          = {
                                            (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                              (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                          }),
                            
                            "<"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "NLP"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          }),
                            ">"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "NLP"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          })
  )
  
  #r_prior(rho,k,location,scale,dff,model, alpha, beta,min(bound),max(bound)) / normalization

  prior_h1<-r_prior(rho,k,location,scale,dff,model, alpha, beta,-1,1)
  switch(hypothesis,
         "!=" = { prior_h1[rho>min(bound_h1)&rho<max(bound_h1)]=0 },
         ">" = { prior_h1[rho<bound_h1[1]]=0 },
         "<" = { prior_h1[rho>bound_h1[2]]=0 }
  )
  prior_h1
  
  }
compute.prior.density.re.h0 <- function(rho,h0, k,location,scale,dff,model, alpha, beta,hypothesis,e) {
  if (model == "Point") return(rep(NA, length(rho)))
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh0 <- switch(model,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "NLP"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )
  
  #r_prior(rho,k,location,scale,dff,model, alpha, beta,min(bound),max(bound)) / normalization
  
  prior_h0<-r_prior(rho,k,location,scale,dff,model, alpha, beta,-1,1)
  switch(hypothesis,
         "!=" = { prior_h0[rho<min(bound_h0)|rho>max(bound_h0)]=0 },
         ">" = { prior_h0[rho>bound_h0[2]]=0 },
         "<" = { prior_h0[rho<bound_h0[1]]=0 }
  )
  prior_h0
  
}



re_prior_plot <-function(k, alpha, beta,h0,location,scale,dff,model,de_an_prior,
                         k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d,hypothesis,e){
  par(mfrow = c(1, 1))

  
  plot.bounds <- switch(hypothesis,
                        ">"  = c(h0, 1),
                        "<"  = c(-1, h0),
                        "!=" = c(-1, 1))
  rr <- seq(plot.bounds[1], plot.bounds[2], 0.0025)
  
  prior.analysis.h1 <- compute.prior.density.re.h1(rr,h0, k,location,scale,dff,model, alpha, beta,hypothesis,e) 
  prior.analysis.h0<- compute.prior.density.re.h0(rr,h0, k,location,scale,dff,model, alpha, beta,hypothesis,e) 
  prior.design <- if (de_an_prior == 0 && model_d != "Point") {
    compute.prior.density.re.h1(rr, k_d,location_d,scale_d,dff_d,model_d, alpha_d, beta_d,hypothesis,e) 
  } else {
    rep(NA, length(rr))
  }
  # Combine all values into one vector
  all_vals <- c(prior.analysis.h1, prior.analysis.h0, prior.design)
  
  # Filter out NA and infinite values
  finite_vals <- all_vals[is.finite(all_vals)]
  
  # Get the max from finite values only
  ylim.max <- max(finite_vals)
  
  
  # Base plot
  plot(rr, prior.analysis.h1, type = "l", lwd = 2,
       xlab = expression(bold(rho)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~rho~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))
  
  lines(rr, prior.analysis.h0, lty = 2, col = "black", lwd = 2)
  
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
      lines(rr, prior.design, lty = 1, col = "gray", lwd = 3)
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



re_bf10_p <-function(D,n,k,h0,hypothesis,location,scale,dff,model,e){
  
  rr  <- seq(from = -.99,to = .99,.01)
  
  # Compute BF10 and t-bounds:
  r.BF10 <- re_BF_bound_10(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  BF10 <- re_BF10(rr,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  
  par(mfrow = c(1, 2))
  # Left plot - BF10:
  main.bf10 <- if (length(r.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(D)~"when r = "~.(format(r.BF10, digits = 4))))
  } else {
    bquote(bold("BF"[10]~"="~.(D)~"when r = "~.(format(r.BF10[1], digits = 4))~"or"~.(format(r.BF10[2], digits = 4))))
  }
  plot(rr, BF10, type = "l", log = "y", xlab = "Correlation", ylab = expression("BF"[10]),
       main = main.bf10, frame.plot = FALSE, xaxt = "n")
  abline(v = r.BF10)
  axis(1, c(-1, 1))
  if (length(r.BF10)) axis(1, round(r.BF10, 2))
  
  # Left plot - BF01:
  BF01   <- 1 / BF10
  r.BF01 <- re_BF_bound_01(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  
  # Check if BF01 = D is possible:
  max.BF01   <- 1 / re_BF10(h0,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  impossible <- (hypothesis == "!=") && (max.BF01 < D || identical(r.BF01, "bound cannot be found"))
  
  plot(rr, BF01, type = "l", log = "y", xlab = "Correlation", ylab = bquote("BF"['01']),
       main = "", frame.plot = FALSE, xaxt = "n")
  axis(1, c(-1, 1))
  if (impossible) {
    title(main = bquote(bold("It is impossible to have BF"[01]~"="~.(D))))
  } else {
    abline(v = r.BF01)
    axis(1, round(r.BF01, 2))
    main.bf01 <- if (length(r.BF01) == 1) {
      bquote(bold("BF"['01']~"="~.(D)~"when r = "~.(format(r.BF01, digits = 4))))
    } else {
      bquote(bold("BF"['01']~"="~.(D)~"when r = "~.(format(r.BF01[1], digits = 4))~"or"~.(format(r.BF01[2], digits = 4))))
    }
    title(main = main.bf01)
  }
}


Power_re<-function(D,k, alpha, beta,h0,hypothesis,location,scale,dff,model, 
                  k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d, de_an_prior,N,e){
  # N range to evaluate power:
  N.min     <- 4
  N.max     <- ceiling(N * 1.2)
  Ns        <- seq(N.min, N.max, length.out = 31)
  TPE <- numeric(length(Ns))
  FPE <- numeric(length(Ns))
  TNE <- numeric(length(Ns))
  FNE <- numeric(length(Ns))
  
  for (i in seq_along(Ns)) {
    r10 <-  re_BF_bound_10(D,Ns[i],k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    r01 <-  re_BF_bound_01(D,Ns[i],k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    
    # Choose correct design prior:
    TPE[i] <- if (de_an_prior == 1)
      re_TPE(r10,Ns[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model,e) else
        re_TPE(r10,Ns[i],k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    FNE[i] <- if (de_an_prior == 1)
      re_FNE(r01,Ns[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model,e) else
        re_FNE(r01,Ns[i],k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    FPE[i]       <- re_FPE(r10,Ns[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    TNE[i]       <- re_TNE(r01,Ns[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    
    
      }
  
  par(mfrow = c(1, 2))
  plot(Ns, TPE, type = "l",
       xlab = "Total sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(Ns,FPE,col = "grey")
  legend(x = N.max*.4,y=.5,              # position of the legend
         legend = c("True positive", "False positive"),  # labels
         col = c("black", "grey"),   # colors matching lines
         lty = 1,                   # line type (solid)
         bty = "n")  
  
  plot(Ns, TNE, type = "l",
       xlab = "Total sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[0][1]~">"~.(D))))
  lines(Ns,FNE,col = "grey")
  legend(x = N.max*.4,y=.5,              # position of the legend
         legend = c("True negative", "False negative"),  # labels
         col = c("black", "grey"),   # colors matching lines
         lty = 1,                   # line type (solid)
         bty = "n")   

}




