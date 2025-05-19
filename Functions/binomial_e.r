adjust_root_10_e <- function(root, n, alpha, beta, location, scale, model, hypothesis, D,e) {
  # Evaluate BF at the root
  BF_val <- bin_e_BF(root,n,alpha,beta,location,scale,model,hypothesis,e)
  
  if (BF_val <= D) {
    # Try root - 1
    BF_prev <- bin_e_BF(root-1,n,alpha,beta,location,scale,model,hypothesis,e)
    if (BF_prev > D) return(root - 1)
    
    # Try root + 1
    BF_next <- bin_e_BF(root+1,n,alpha,beta,location,scale,model,hypothesis,e)
    if (BF_next > D) return(root + 1)
  }
  
  # Return original if already valid or no better nearby found
  return(root)
}

adjust_root_01_e <- function(root, n, alpha, beta, location, scale, model, hypothesis, D,e) {
  # Evaluate BF at the root
  BF_val <- 1/bin_e_BF(root,n,alpha,beta,location,scale,model,hypothesis,e)
  
  if (BF_val <= D) {
    # Try root - 1
    BF_prev <- 1/bin_e_BF(root-1,n,alpha,beta,location,scale,model,hypothesis,e)
    if (BF_prev > D) return(root - 1)
    
    # Try root + 1
    BF_next <- 1/bin_e_BF(root+1,n,alpha,beta,location,scale,model,hypothesis,e)
    if (BF_next > D) return(root + 1)
  }
  
  # Return original if already valid or no better nearby found
  return(root)
}

bin_e_BF<-function(x,n,alpha,beta,location,scale,model,hypothesis,e){
  BF = NA
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = location+e, b = 1),
                      "<" = c(a = 0, b = location+e),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = location, b = location+e),
                      "<" = c(a = location+e, b = location),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  
  normalizationh1 <- switch(hypothesis,
                            "!=" = {
                              if (model == "beta") {
                                1 - (pbeta(bound_h1[2], alpha, beta) - pbeta(bound_h1[1], alpha, beta))
                              } else if (model == "Moment") {
                                (pmom(1 - location, tau = scale^2) - pmom(bound_h1[2] - location, tau = scale^2)) +
                                  (pmom(bound_h1[1] - location, tau = scale^2) - pmom(-1 - location, tau = scale^2))
                              }
                            },
                            "<" = ,
                            ">" = {
                              if (model == "beta") {
                                pbeta(bound_h1[2], alpha, beta) - pbeta(bound_h1[1], alpha, beta)
                              } else if (model == "Moment") {
                                pmom(bound_h1[2] - location, tau = scale^2) - pmom(bound_h1[1] - location, tau = scale^2)
                              }
                            }
  )
  
  normalizationh0 <- switch(model,
                            "beta"      =   pbeta(bound_h0[2], alpha, beta) - pbeta(bound_h0[1], alpha, beta),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )
  
    for (i in 1:length(x)){
  int  <- function(prop){dbinom(x[i], size=n, prob=prop) *bin_prior(prop,alpha,beta,location,scale,model)
  }
  
  if (hypothesis == "!="){
    lh1 = integrate(int,lower = 0,upper = bound_h1[1], rel.tol=1e-5,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-5,stop.on.error = F)$value 
  }else{
    lh1 = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-5,stop.on.error = F)$value
    
  }
  lh0 = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value

    
    BF[i] = (lh1/normalizationh1)/(lh0/normalizationh0)
  }
  return(BF)
  
}


bin_e_BF_bound_10 <-function(D,n,alpha,beta,location,scale,model,hypothesis,e){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)- D
  }
  
  x <- tryCatch(uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)
  
  results <- round(results[!is.na(results)])
  if (length(results) == 0) return("bound cannot be found")
  
  
  results <- sapply(results, function(root) {
    adjust_root_10_e(root, n, alpha, beta, location, scale, model, hypothesis, D,e)
  })
  
  BF.vals  <- bin_e_BF(results,n,alpha,beta,location,scale,model,hypothesis,e)
  BF.close <- which(BF.vals > D)
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}


bin_e_BF_bound_01 <-function(D,n,alpha,beta,location,scale,model,hypothesis,e){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    1/bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)- D
  }
  
  x <- tryCatch(uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)
  
  results <- round(results[!is.na(results)])
  if (length(results) == 0) return("bound cannot be found")
  
  
  results <- sapply(results, function(root) {
    adjust_root_01_e(root, n, alpha, beta, location, scale, model, hypothesis, D,e)
  })
  
  BF.vals  <- 1/bin_e_BF(results,n,alpha,beta,location,scale,model,hypothesis,e)
  BF.close <- which(BF.vals > D)
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}


bin_e_TPE<-function(x,n,alpha,beta,location,scale,model,hypothesis,e){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(x)
  

  if (model =="Point"){
    TPE = switch(hypothesis,
               "!=" = {
                 
                 switch(length(x)==2,
                        "1" ={pbinom(min(x),n,location,lower.tail = T)+ pbinom(max(x)-1,n,location,lower.tail = F)},
                        "0"=  {    
                          switch(x/n>location,
                                 "1" = pbinom(x,n,location,lower.tail = F),
                                 "0" = pbinom(x-1,n,location,lower.tail = T))
                          
                        })
                 },
               ">"  = {pbinom(x,n,location,lower.tail = F)},
               "<"  = {pbinom(x-1,n,location,lower.tail = T)}
    )
    return(TPE)
  }

  
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = location+e, b = 1),
                      "<" = c(a = 0, b = location+e),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  normalizationh1 <- switch(hypothesis,
                            "!=" = {
                              if (model == "beta") {
                                1 - (pbeta(bound_h1[2], alpha, beta) - pbeta(bound_h1[1], alpha, beta))
                              } else if (model == "Moment") {
                                (pmom(1 - location, tau = scale^2) - pmom(bound_h1[2] - location, tau = scale^2)) +
                                  (pmom(bound_h1[1] - location, tau = scale^2) - pmom(-1 - location, tau = scale^2))
                              }
                            },
                            "<" = ,
                            ">" = {
                              if (model == "beta") {
                                pbeta(bound_h1[2], alpha, beta) - pbeta(bound_h1[1], alpha, beta)
                              } else if (model == "Moment") {
                                pmom(bound_h1[2] - location, tau = scale^2) - pmom(bound_h1[1] - location, tau = scale^2)
                              }
                            }
  )  
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      pbinom(min(x), n, prop, lower.tail = TRUE) +
                        pbinom(max(x) - 1, n, prop, lower.tail = FALSE)
                    } else {
                      mapply(function(x_i, n_i, p_i) {
                        if (x_i / n_i > p_i) {
                          pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                        } else {
                          pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                        }
                      }, x, n, prop)
                    }
                  },
                  ">" = pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "<" = pbinom(x, n, prop, lower.tail = TRUE)
    )
    
    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalizationh1
  }
  
  if(hypothesis == "!="){
    TPE = integrate(int,lower = 0,upper = bound_h1[1], rel.tol = 1e-5)$value + integrate(int,lower = bound_h1[2],upper = 1, rel.tol = 1e-5)$value
  }else{
    TPE = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol = 1e-5)$value
  }
  return(TPE)
  
}


bin_e_FNE<-function(x,n,alpha,beta,location,scale,model,hypothesis,e){
  
  if (length(x) == 0 || any(x == "bound cannot be found")) return(x)
  
  
  if (model =="Point"){
    FNE = switch(hypothesis,
                 "!=" = {
                   
                   switch(length(x)==2,
                          "1" ={pbinom(max(x),n,location,lower.tail = T)- pbinom(min(x)-1,n,location,lower.tail = T)},
                          "0"=  {    
                            switch(x/n>location,
                                   "1" = pbinom(x,n,location,lower.tail = T),
                                   "0" = pbinom(x-1,n,location,lower.tail = F))
                            
                          })},
                 ">"  = {pbinom(x,n,location,lower.tail = T)},
                 "<"  = {pbinom(x-1,n,location,lower.tail = F)}
    )
    return(FNE)
  }

  
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = location+e, b = 1),
                      "<" = c(a = 0, b = location+e),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  normalizationh1 <- switch(hypothesis,
                            "!=" = {
                              if (model == "beta") {
                                1 - (pbeta(bound_h1[2], alpha, beta) - pbeta(bound_h1[1], alpha, beta))
                              } else if (model == "Moment") {
                                (pmom(1 - location, tau = scale^2) - pmom(bound_h1[2] - location, tau = scale^2)) +
                                  (pmom(bound_h1[1] - location, tau = scale^2) - pmom(-1 - location, tau = scale^2))
                              }
                            },
                            "<" = ,
                            ">" = {
                              if (model == "beta") {
                                pbeta(bound_h1[2], alpha, beta) - pbeta(bound_h1[1], alpha, beta)
                              } else if (model == "Moment") {
                                pmom(bound_h1[2] - location, tau = scale^2) - pmom(bound_h1[1] - location, tau = scale^2)
                              }
                            }
  )
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      pbinom(max(x), n, prop, lower.tail = TRUE) - pbinom(min(x) - 1, n, prop, lower.tail = TRUE)
                    } else {
                      if ((x / n) > prop) {
                        pbinom(x, n, prop, lower.tail = TRUE)
                      } else {
                        pbinom(x - 1, n, prop, lower.tail = FALSE)
                      }
                    }
                  },
                  ">" = pbinom(x , n, prop, lower.tail = TRUE),
                  "<" = pbinom(x - 1, n, prop, lower.tail = FALSE)
    )
    
    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalizationh1
  }
  if(hypothesis == "!="){
    FNE = integrate(int,lower = 0,upper = bound_h1[1], rel.tol = 1e-5)$value + integrate(int,lower = bound_h1[2],upper = 1, rel.tol = 1e-5)$value
  }else{
    FNE = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol = 1e-5)$value
  }
  
  
  return(FNE)
  
}

bin_e_FPE<-function(x,n,alpha,beta,location,scale,model,hypothesis,e){
  
  if (length(x) == 0 || any(x == "bound cannot be found")) return(x)
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = location, b = location+e),
                      "<" = c(a = location+e, b = location),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  normalizationh0 <- switch(model,
                            "beta"      =   pbeta(bound_h0[2], alpha, beta) - pbeta(bound_h0[1], alpha, beta),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )
  
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      pbinom(min(x), n, prop, lower.tail = TRUE) +
                        pbinom(max(x) - 1, n, prop, lower.tail = FALSE)
                    } else {
                      mapply(function(x_i, n_i, p_i) {
                        if (x_i / n_i > p_i) {
                          pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                        } else {
                          pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                        }
                      }, x, n, prop)
                    }
                  },
                  ">" = pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "<" = pbinom(x, n, prop, lower.tail = TRUE)
    )
    
    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalizationh0
  }
  
  
    FPE = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol = 1e-5)$value
    return(FPE)

}

bin_e_TNE<-function(x,n,alpha,beta,location,scale,model,hypothesis,e){
  
  
  if (length(x) == 0 || any(x == "bound cannot be found")) return(x)
  
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = location, b = location+e),
                      "<" = c(a = location+e, b = location),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  
  normalizationh0 <- switch(model,
                            "beta"      =   pbeta(bound_h0[2], alpha, beta) - pbeta(bound_h0[1], alpha, beta),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      pbinom(max(x), n, prop, lower.tail = TRUE) - pbinom(min(x) - 1, n, prop, lower.tail = TRUE)
                    } else {
                      if ((x / n) > prop) {
                        pbinom(x, n, prop, lower.tail = TRUE)
                      } else {
                        pbinom(x - 1, n, prop, lower.tail = FALSE)
                      }
                    }
                  },
                  ">" = pbinom(x , n, prop, lower.tail = TRUE),
                  "<" = pbinom(x - 1, n, prop, lower.tail = FALSE)
    )
    
    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalizationh0
  }

    TNE = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol = 1e-5)$value
  
  
  
  return(TNE)
}

bin_e_N_finder <-function(D,target,alpha,beta,location,scale,model,hypothesis,
                        alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP,e){
  lower = 10
  upper = 10000
  
  b10 =  bin_e_BF_bound_10(D,lower,alpha,beta,location,scale,model,hypothesis,e)
  
  TPE_lo <- if (de_an_prior == 1)
    bin_e_TPE(b10,lower,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e) else
      bin_e_TPE(b10,lower,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
  if (TPE_lo > target) return(lower)
  FPE_lo <-  bin_e_FPE(b10,lower,alpha,beta,location,scale,model,hypothesis,e)
  if (TPE_lo > target&FPE_lo<FP) return(lower)
  
  Power_root <- function(N){
    N =round(N)
    x = bin_e_BF_bound_10(D,N,alpha,beta,location,scale,model,hypothesis,e)

    if(de_an_prior == 1){
      pro = bin_e_TPE(x,N,alpha,beta,location,scale,model,hypothesis,e)
    } else{
      pro = bin_e_TPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
    }
    return(pro-target)
  }

  N.power = round(uniroot(Power_root,lower = lower,upper = upper)$root)+1

    while(TRUE) {
    b10 <- bin_e_BF_bound_10(D,N.power,alpha,beta,location,scale,model,hypothesis,e)
    pro <- if (de_an_prior == 0) {
      bin_e_TPE(b10, N.power, alpha_d, beta_d, location_d, scale_d, model_d, hypothesis,e)
    } else {
      bin_e_TPE(b10,N.power,alpha,beta,location,scale,model,hypothesis,e)
    }
    
    if (pro > target) break
    N.power <- N.power + 1
  }
  b10 = bin_e_BF_bound_10(D,N.power,alpha,beta,location,scale,model,hypothesis,e)
  FPE =  bin_e_FPE(b10,N.power,alpha,beta,location,scale,model,hypothesis,e)
  if (FPE <= FP) return(N.power)
  
  alpha.root <- function(n) {
    n=round(n)
    b10 <- bin_e_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis,e)
    bin_e_FPE(b10,n,alpha,beta,location,scale,model,hypothesis,e)-FP
  }
  N.alpha = round(uniroot(alpha.root,lower = N.power,upper = upper)$root)
  return(N.alpha)
    
  }


bin_e_table<-function(D,target,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,N, mode_bf,FP,e){
  if (mode_bf == "0") n = N else n =  bin_e_N_finder(D,target,alpha,beta,location,scale,model,hypothesis,
                                                     alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP,e)
  
  # b bounds:
  b10 <- bin_e_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis,e)
  b01 <-  bin_e_BF_bound_01(D,n,alpha,beta,location,scale,model,hypothesis,e)
  
  max_BF <- 1 /bin_e_BF(ceiling(n/2),n,alpha,beta,location,scale,model,hypothesis,e)
  
  # FPE and TPE:
  FPE       <- bin_e_FPE(b10,n,alpha,beta,location,scale,model,hypothesis,e)
  if (de_an_prior == 1) {
    TPE          <- bin_e_TPE(b10,n,alpha,beta,location,scale,model,hypothesis,e)
    TPR_alpha    <- alpha
    TPR_beta     <- beta
    TPR_location <- location
    TPR_scale    <- scale
    TPR_model    <- model
    
  } else {
    TPE          <- bin_e_TPE(b10,n,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
    TPR_alpha    <- alpha_d
    TPR_beta     <- beta_d
    TPR_location <- location_d
    TPR_scale    <- scale_d
    TPR_model    <- model_d
  }
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | b01 == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- bin_e_FNE(b01,n,TPR_alpha,TPR_beta,TPR_location,TPR_scale,TPR_model,hypothesis,e)
    TNE <- bin_e_TNE(b01,n,alpha,beta,location,scale,model,hypothesis,e)
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

bin_e_bf10 <-function(DD,n,alpha,beta,location,scale,model,hypothesis,e){
  x= seq(from = 0,to =n,by= 3)
  
  # Compute BF10 and x-bounds:
  b.BF10     <- bin_e_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis,e)
  BF10_at_b  <- round(bin_e_BF(b.BF10,n,alpha,beta,location,scale,model,hypothesis,e),2)
  BF10       <- bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)

  if (length(b.BF10)== 2){  part1 = bquote(bold("BF"[10] ~ "=" ~ .(BF10_at_b[1]) / .(BF10_at_b[2])))}else{part1 = bquote(bold("BF"[10] ~ "=" ~ .(BF10_at_b[1])))}
  
  if (length(b.BF10)== 2){  part2 = bquote("when x = " ~ .(b.BF10[1]) / .(b.BF10[2]))}else{  part2 <- bquote("when x = " ~ .(b.BF10[1]))}
  
  main <- bquote(bold(.(part1) ~ .(part2)))
  
  par(mfrow = c(1, 2))
  plot(x, BF10, type = "l", log = "y", xlab = "Number of success", ylab = expression("BF"[10]),
       main = main, frame.plot = FALSE, xaxt = "n")
  abline(v = b.BF10)
  axis(1, c(0, n))
  if (length(b.BF10)) axis(1, round(b.BF10, 2))
  
  # right plot - BF01:
  BF01   <- 1 / BF10
  b.BF01 <- bin_e_BF_bound_01(D,n,alpha,beta,location,scale,model,hypothesis,e)
  BF01_at_b <- round(1/bin_e_BF(b.BF01,n,alpha,beta,location,scale,model,hypothesis,e),2)
  
  # Check if BF01 = D is possible:
  max.BF01 <- 1 / bin_e_BF(round(n/2),n,alpha,beta,location,scale,model,hypothesis,e)
  impossible <- (hypothesis == "!=") && (max.BF01 < D || identical(b.BF01, "bound cannot be found"))
  
  plot(x, BF01, type = "l", log = "y", xlab = "Number of success", ylab = bquote("BF"['01']),
       main = "", frame.plot = FALSE, xaxt = "n")
  axis(1, c(0, n))
  
  if (impossible) {
    title(main = bquote(bold("It is impossible to have BF"[01]~"="~.(D))))
  } else {
    abline(v = b.BF01)
    axis(1, round(b.BF01, 2))
    
    if (length(b.BF10) == 2) {
      part1 <- bquote("BF"[10] == bold(.(BF10_at_b[1])) / bold(.(BF10_at_b[2])))
      part2 <- bquote(bold("when x = " ~ bold(.(b.BF10[1])) / bold(.(b.BF10[2]))))
    } else {
      part1 <- bquote("BF"[10] == bold(.(BF10_at_b[1])))
      part2 <- bquote(bold("when x = " ~ bold(.(b.BF10[1]))))
    }
    main.bf01 = bquote(bold(.(part1) ~ .(part2)))
    title(main = main.bf01)
  }

}


Power_e_bin<-function(D,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d, de_an_prior,N,e){

  smin = 10
  smax = N*1.2
  sN = ceiling(seq(smin,smax , by = (smax-smin)/50))
  
  power =  array(NA, dim = c(length(sN)))
  
  for ( i in 1:length(sN)){
    x = bin_e_BF_bound_10 (D,sN[i],alpha,beta,location,scale,model,hypothesis,e)

 
    if(de_an_prior ==1){
      power[i] = bin_e_TPE(x,sN[i],alpha,beta,location,scale,model,hypothesis,e)
    
    }else{
      power[i] = bin_e_TPE(x,sN[i],alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
    }

  }
  par(mfrow = c(1, 1))
  
  plot(sN,power,type="l",main = "",frame.plot = FALSE,xlab = "Total sample size", ylab = "Probability of True positive evidence",xlim = c(10,max(sN)), 
       ylim = c(0,1) )
  
}

compute.prior.density.be.h1 <- function(prop,alpha,beta,location,scale,model,hypothesis,e) {
  if (model == "Point") return(rep(NA, length(prop)))
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = location+e, b = 1),
                      "<" = c(a = 0, b = location+e),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  
  prior_h1<- bin_prior(prop,alpha,beta,location,scale,model)
  switch(hypothesis,
         "!=" = { prior_h1[prop>min(bound_h1)&prop<max(bound_h1)]=0 },
         ">" = { prior_h1[prop<bound_h1[1]]=0 },
         "<" = { prior_h1[prop>bound_h1[2]]=0 }
  )
  prior_h1
}


compute.prior.density.be.h0 <- function(prop,alpha,beta,location,scale,model,hypothesis,e) {
  if (model == "Point") return(rep(NA, length(prop)))
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = location, b = location+e),
                      "<" = c(a = location+e, b = location),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  
  prior_h0<- bin_prior(prop,alpha,beta,location,scale,model)
  switch(hypothesis,
         "!=" = { prior_h0[prop<min(bound_h0)|prop>max(bound_h0)]=0 },
         ">" = { prior_h0[prop>bound_h0[2]]=0 },
         "<" = { prior_h0[prop<bound_h0[1]]=0 }
  )
  prior_h0
}





bin_e_prior_plot <-function(alpha,beta,location,scale,model,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,de_an_prior,e){
  par(mfrow = c(1, 1))
  bounds <- switch(hypothesis,
                        ">"  = c(location, 1),
                        "<"  = c(0, location),
                        "!=" = c(0, 1))
  prop           <- seq(bounds[1],bounds[2],.002)
  
  
  prior.analysis.h1 <- compute.prior.density.be.h1(prop,alpha,beta,location,scale,model,hypothesis,e) 
  prior.analysis.h0<-  compute.prior.density.be.h0(prop,alpha,beta,location,scale,model,hypothesis,e) 
  prior.design <- if (de_an_prior == 0 && model_d != "Point") {
    compute.prior.density.be.h1(prop,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e) 
  } else {
    rep(NA, length(prop))
  }
  # Combine all values into one vector
  all_vals <- c(prior.analysis.h1, prior.analysis.h0, prior.design)
  
  # Filter out NA and infinite values
  finite_vals <- all_vals[is.finite(all_vals)]
  
  # Get the max from finite values only
  ylim.max <- max(finite_vals)
  
  
  # Base plot
  plot(prop, prior.analysis.h1, type = "l", lwd = 2,
       xlab =bquote(bold(rho)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~rho~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))
  
  lines(prop, prior.analysis.h0, lty = 2, col = "black", lwd = 2)
  
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
      lines(prop, prior.design, lty = 1, col = "gray", lwd = 3)
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


