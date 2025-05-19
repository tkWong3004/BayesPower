adjust_root_10 <- function(root, n, alpha, beta, location, scale, model, hypothesis, D) {
  # If root is less than 0, return NA
  if (root < 0) return(NA)
  
  # Evaluate BF at the root
  BF_val <- bin_BF(root, n, alpha, beta, location, scale, model, hypothesis)
  
  if (BF_val <= D) {
    # Try root - 1 only if root > 0
    if (root > 0) {
      BF_prev <- bin_BF(root - 1, n, alpha, beta, location, scale, model, hypothesis)
      if (BF_prev > D) return(root - 1)
    }
    
    # Try root + 1
    BF_next <- bin_BF(root + 1, n, alpha, beta, location, scale, model, hypothesis)
    if (BF_next > D) return(root + 1)
  }
  
  # Return original if already valid or no better nearby found
  return(root)
}


adjust_root_01 <- function(root, n, alpha, beta, location, scale, model, hypothesis, D) {
  # Evaluate BF at the root
  BF_val <- 1/bin_BF(root, n, alpha, beta, location, scale, model, hypothesis)
  
  if (BF_val <= D) {
    # Try root - 1
    BF_prev <- 1/bin_BF(root - 1, n, alpha, beta, location, scale, model, hypothesis)
    if (BF_prev > D) return(root - 1)
    
    # Try root + 1
    BF_next <- 1/bin_BF(root + 1, n, alpha, beta, location, scale, model, hypothesis)
    if (BF_next > D) return(root + 1)
  }
  
  # Return original if already valid or no better nearby found
  return(root)
}


bin_prior <-function(prop,alpha,beta,location,scale,model){
  
  switch(model,
         "beta" = dbeta(prop, alpha,beta),
         "Moment" = dnlp(prop,location,scale))
}
bin_BF<-function(x,n,alpha,beta,location,scale,model,hypothesis){
  BF = NA
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = 0, b = location),
                   "!=" = c(a = 0, b = 1)
  )
  

  normalization <- if (hypothesis == "!=") {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)
    
  } else {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = pbeta(bound[2],alpha,beta)-pbeta(bound[1],alpha,beta))
    }
  for( i in 1:length(x)){
    int  <- function(prop){dbinom(x[i], size=n, prob=prop) *bin_prior(prop,alpha,beta,location,scale,model)}
    lh1 <- integrate(int, lower = bound[1], upper = bound[2], rel.tol = 1e-5)$value / normalization
    lh0 <- dbinom(x[i], size = n, prob = location)
    BF[i] = lh1 / lh0
  }

  
  return(BF)
  
}

bin_BF_bound_10 <-function(D,n,alpha,beta,location,scale,model,hypothesis){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)- D
  }
  
  x <- tryCatch(uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)
  #results <- tryCatch(uniroot.all(Bound_finding, lower = 0 ,upper = n), error = function(e) NA)
  results <- round(results[!is.na(results) & is.finite(results)])
  
  if (length(results) == 0) return("bound cannot be found")
  

  results <- sapply(results, function(root) {
    adjust_root_10(root, n, alpha, beta, location, scale, model, hypothesis, D)
  })
  
  
  BF.vals  <- bin_BF(results,n,alpha,beta,location,scale,model,hypothesis)
  
  BF.close <- which(BF.vals > D)
  if (length(BF.close) == 0 || all(!is.finite(BF.close))) return("bound cannot be found")
  return(results[BF.close])
}

bin_BF_bound_01 <-function(D,n,alpha,beta,location,scale,model,hypothesis){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    1/bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)- D
  }
  
  x <- tryCatch(uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)
  
  results <- round(results[!is.na(results)])
  if (length(results) == 0) return("bound cannot be found")
  
  
  results <- sapply(results, function(root) {
    adjust_root_01(root, n, alpha, beta, location, scale, model, hypothesis, D)
  })
  
  
  BF.vals  <- 1/bin_BF(results,n,alpha,beta,location,scale,model,hypothesis)
  
  BF.close <- which(BF.vals > D)
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}
 
bin_TPE<-function(x,n,alpha,beta,location,scale,model,hypothesis){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  if (model =="Point"){
    TPE = switch(hypothesis,
               "!=" = {
                 
                 switch(length(x)==2,
                        "1" ={pbinom(min(x),n,location,lower.tail = T)+ pbinom(max(x)-1,n,location,lower.tail = F)},
                        "0"=  {    
                          switch(x/n>location,
                                 "1" = pbinom(x-1,n,location,lower.tail = F),
                                 "0" = pbinom(x,n,location,lower.tail = T))
                          
                        })
                 },
               ">"  = {pbinom(x-1,n,location,lower.tail = F)},
               "<"  = {pbinom(x,n,location,lower.tail = T)}
    )
    return(TPE)
  }
  
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = 0, b = location),
                   "!=" = c(a = 0, b = 1)
  )
  normalization <- if (hypothesis == "!=") {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)
    
  } else {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = pbeta(bound[2],alpha,beta)-pbeta(bound[1],alpha,beta))
  }
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
    
    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalization
  }
  
  TPE = integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-5)$value
  
  return(TPE)
  
}


bin_FNE<-function(x,n,alpha,beta,location,scale,model,hypothesis){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)
  
  if (model == "Point") {
    FNE <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      pbinom(max(x), n, location, lower.tail = TRUE) - pbinom(min(x) - 1, n, location, lower.tail = TRUE)
                    } else {
                      if ((x / n) > location) {
                        pbinom(x, n, location, lower.tail = TRUE)
                      } else {
                        pbinom(x - 1, n, location, lower.tail = FALSE)
                      }
                    }
                  },
                  ">" = pbinom(x, n, location, lower.tail = TRUE),
                  "<" = pbinom(x - 1, n, location, lower.tail = FALSE)
    )
    return(FNE)
  }
  

  
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = 0, b = location),
                   "!=" = c(a = 0, b = 1)
  )
  
  normalization <- if (hypothesis == "!=") {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)
    
  } else {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = pbeta(bound[2],alpha,beta)-pbeta(bound[1],alpha,beta))
  }
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
    
    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalization
  }
  FNE = integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-5)$value
  return(FNE)
  
}


bin_FPE<-function(x,n,location,hypothesis){
  
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)
  
  FPE <- switch(hypothesis,
                "!=" = {
                  if (length(x) == 2) {
                    pbinom(min(x), n, location, lower.tail = TRUE) +
                      pbinom(max(x) - 1, n, location, lower.tail = FALSE)
                  } else {
                    mapply(function(x_i, n_i, p_i) {
                      if (x_i / n_i > p_i) {
                        pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                      } else {
                        pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                      }
                    }, x, n, location)
                  }
                },
                ">" = pbinom(x - 1, n, location, lower.tail = FALSE),
                "<" = pbinom(x, n, location, lower.tail = TRUE)
  )
  
    return(FPE)

}

bin_TNE<-function(x,n,location,hypothesis){
  
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)
  
  
  TNE <- switch(hypothesis,
                "!=" = {
                  if (length(x) == 2) {
                    pbinom(max(x), n, location, lower.tail = TRUE) - pbinom(min(x) - 1, n, location, lower.tail = TRUE)
                  } else {
                    if ((x / n) > location) {
                      pbinom(x, n, location, lower.tail = TRUE)
                    } else {
                      pbinom(x - 1, n, location, lower.tail = FALSE)
                    }
                  }
                },
                ">" = pbinom(x, n, location, lower.tail = TRUE),
                "<" = pbinom(x - 1, n, location, lower.tail = FALSE)
  )
  
  return(TNE)
  
}

bin_N_finder <-function(D,target,alpha,beta,location,scale,model,hypothesis,
                        alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP){
  lower = 10
  upper = 10000
  
  b10 = bin_BF_bound_10(D,lower,alpha,beta,location,scale,model,hypothesis)
  TPE_lo <- if (de_an_prior == 1)
    bin_TPE(b10,lower,alpha,beta,location,scale,model,hypothesis) else
      bin_TPE(b10,lower,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
  if (TPE_lo > target) return(lower)
  FPE_lo <-  bin_FPE(b10,lower,location,hypothesis)
  if (TPE_lo > target&FPE_lo<FP) return(lower)

  
  Power_root <- function(N){
    N =round(N)
    b10 = bin_BF_bound_10 (D,N,alpha,beta,location,scale,model,hypothesis)
    pro <- if (de_an_prior==0){
      bin_TPE(b10,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
      }else bin_TPE(b10,N,alpha,beta,location,scale,model,hypothesis)
    
    pro-target
  }

  N.power = round(uniroot(Power_root,lower = lower,upper = upper)$root)+1
  while(TRUE) {
    b10 <- bin_BF_bound_10(D, N.power, alpha, beta, location, scale, model, hypothesis)
    pro <- if (de_an_prior == 0) {
      bin_TPE(b10, N.power, alpha_d, beta_d, location_d, scale_d, model_d, hypothesis)
    } else {
      bin_TPE(b10, N.power, alpha, beta, location, scale, model, hypothesis)
    }
    
    if (pro > target) break
    N.power <- N.power + 1
  }
  
  
  b10 = bin_BF_bound_10(D,N.power,alpha,beta,location,scale,model,hypothesis)
  FPE = bin_FPE(b10,N.power,location,hypothesis)
  if (FPE <= FP) return(N.power)
  
  
  alpha.root <- function(n) {
    n=round(n)
    b10 <- bin_BF_bound_10 (D,n,alpha,beta,location,scale,model,hypothesis)
    bin_FPE(b10,n,location,hypothesis)-FP
  }
  N.alpha = round(uniroot(alpha.root,lower = N.power,upper = upper)$root)
  return(N.alpha)
}


bin_table<-function(D,target,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,N, mode_bf,FP){
  if (mode_bf == "0") n = N else n =  bin_N_finder(D,target,alpha,beta,location,scale,model,hypothesis,
                                 alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP)
  
  
  # b bounds:
  b10 <- bin_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis)
  b01 <-  bin_BF_bound_01(D,n,alpha,beta,location,scale,model,hypothesis)
  
  
  # max BF10 possible:
  max_BF <- 1 / bin_BF(round(location*n),n,alpha,beta,location,scale,model,hypothesis)
  BF_D   <- b10
    
  # FPE and TPE:
  FPE       <- bin_FPE(b10,n,location,hypothesis)
  if (de_an_prior == 1) {
    TPE          <- bin_TPE(b10,n,alpha,beta,location,scale,model,hypothesis)
    TPR_alpha    <- alpha
    TPR_beta     <- beta
    TPR_location <- location
    TPR_scale    <- scale
    TPR_model    <- model
    
  } else {
    TPE          <- bin_TPE(b10,n,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
    TPR_alpha    <- alpha_d
    TPR_beta     <- beta_d
    TPR_location <- location_d
    TPR_scale    <- scale_d
    TPR_model    <- model_d
  }
  
  
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- bin_FNE(b01,n,TPR_alpha,TPR_beta,TPR_location,TPR_scale,TPR_model,hypothesis)
    TNE <- bin_TNE(b10,n,location,hypothesis)
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

bin_bf10 <-function(D,n,alpha,beta,location,scale,model,hypothesis){
  x= seq(from = 0,to =n,by= 3)
  
  # Compute BF10 and x-bounds:
  b.BF10 <- bin_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis)
  BF10_at_b <- round(bin_BF(b.BF10,n,alpha,beta,location,scale,model,hypothesis),2)
  BF10 <- bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)

  par(mfrow = c(1, 2))
  
  if (length(b.BF10)== 2){  part1 = bquote(bold("BF"[10] ~ "=" ~ .(BF10_at_b[1]) / .(BF10_at_b[2])))}else{part1 = bquote(bold("BF"[10] ~ "=" ~ .(BF10_at_b[1])))}
  
  if (length(b.BF10)== 2){  part2 = bquote("when x = " ~ .(b.BF10[1]) / .(b.BF10[2]))}else{  part2 <- bquote("when x = " ~ .(b.BF10[1]))}
  
  main <- bquote(bold(.(part1) ~ .(part2)))
  
  plot(x, BF10, type = "l", log = "y", xlab = "Number of success", ylab = expression("BF"[10]),
       main = main, frame.plot = FALSE, xaxt = "n")
  abline(v = b.BF10)
  axis(1, c(0, n))
  if (length(b.BF10)) axis(1, round(b.BF10, 2))
  
  
  # right plot - BF01:
  BF01   <- 1 / BF10
  b.BF01 <- bin_BF_bound_01(D,n,alpha,beta,location,scale,model,hypothesis)
  BF01_at_b <- round(1/bin_BF(b.BF01,n,alpha,beta,location,scale,model,hypothesis),2)
  
  # Check if BF01 = D is possible:
  max.BF01 <- 1 / bin_BF(round(location*n),n,alpha,beta,location,scale,model,hypothesis)
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

Power_bin<-function(D,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d, de_an_prior,N){
  par(mfrow = c(1, 1))
  smin = 10
  smax = N*1.2
  Ns = ceiling(seq(smin,smax , by = (smax-smin)/50))
  
  power =  array(NA, dim = c(length(Ns)))
  
  for ( i in 1:length(Ns)){
    x = bin_BF_bound_10(D,Ns[i],alpha,beta,location,scale,model,hypothesis)
 
    if(de_an_prior ==1){
      power[i] = bin_TPE(x,Ns[i],alpha,beta,location,scale,model,hypothesis)
    }else{
      power[i] =bin_TPE(x,Ns[i],alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
    }

  }
  plot(Ns, power, type = "l",
       xlab = "Sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  
  
}


compute.prior.density.b <- function(prop,alpha,beta,location,scale,model,hypothesis) {
  if (model == "Point") return(rep(NA, length(prop)))
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = 0, b = location),
                   "!=" = c(a = 0, b = 1)
  )
  normalization <- if (hypothesis == "!=") {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)
    
  } else {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = pbeta(bound[2],alpha,beta)-pbeta(bound[1],alpha,beta))
  }
  bin_prior(prop,alpha,beta,location,scale,model)/ normalization
}


bin_prior_plot <-function(alpha,beta,location,scale,model,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,de_an_prior){
  par(mfrow = c(1, 1))
  bound          <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = 0, b = location),
                   "!=" = c(a = 0, b = 1)
  )
  prop           <- seq(bound[1],bound[2],.01)
  normalization  <- integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound[1],upper = bound[2],rel.tol = 1e-10)$value
  prior.analysis <- compute.prior.density.b(prop,alpha,beta,location,scale,model,hypothesis)
  prior.design   <- if (de_an_prior == 0 && model_d != "Point")
    compute.prior.density.b(prop,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis) else
      rep(NA, length(prop))

  # Combine all values into one vector
  all_vals <- c(prior.analysis, prior.design)
  
  # Filter out NA and infinite values
  finite_vals <- all_vals[is.finite(all_vals)]
  
  # Get the max from finite values only
  ylim.max <- max(finite_vals)
  # Base plot:
  plot(prop, prior.analysis, type = "l", lwd = 2,
       xlab = expression(bold(delta)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~rho~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))
  
  # If design prior != analysis prior:
  if (de_an_prior == 0) {
    if (model_d == "Point")
      arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim.max, length = 0.1, col = "black", lty = 2) else
        lines(prop, prior.design, lty = 2)
    
    # Add legend:
    legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }
}




