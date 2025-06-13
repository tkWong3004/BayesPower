
#Fisher
r_mean <-function(r){
  as.numeric(r)
  (1/2)*log((1+r)/(1-r))
}

r_sd <-function(N){
  1/sqrt(N-3)
}
#prior
d_strechted_beta <-function(rho,k,a,b){
  alpha = beta=1/k
  d_beta(rho, alpha, beta,a,b)
  #2^((k-2)/k)*(1-rho^2)^((1-k)/k)/beta(1/k,1/k)
  
}

p_beta <-function(rho, alpha, beta,a,b){
  pBeta_ab(
    rho,
    shape1 = alpha,
    shape2 = beta,
    a = a,
    b = b
  )
}


d_beta <- function(rho, alpha, beta,a,b) {

  # Beta function
  B_ab <- beta(alpha, beta)
  
  # Compute the PDF
  pdf_value <- ((rho - a)^(alpha - 1) * (b - rho)^(beta - 1)) / ((b - a)^(alpha + beta - 1) * B_ab)
  
  return(pdf_value)
}


# likelihood of non-local prior
dnlp <-function(delta,mu,ta){
  ((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
}

r_prior<- function(rho,k,location,scale,dff,model, alpha, beta,a,b){
  
  switch(model,
         "Normal" = dnorm(rho,location,scale),
         "d_beta"   = d_strechted_beta(rho,k,a,b),
          "NLP"   = dnlp(rho,location,scale),
          "t_dis" = tstude(rho,location,scale,dff),
         "beta" = d_beta(rho, alpha, beta,a,b))
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
  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,-1,1)-p_beta(bound[1], 1/k,1/k,-1,1) , 
           "beta" = p_beta(bound[2], alpha, beta,-1,1)-p_beta(bound[1], alpha, beta,-1,1),
           "NLP"   = {pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})
  

  # Define the integrand function for marginal likelihood under H1
  int <- function(rho, ri) {
    d_cor(ri, rho, n) * r_prior(rho, k, location, scale, dff, model, alpha, beta, min(bound), max(bound))
  }
  
  # Compute Bayes factors for each observed correlation ri
  x <- sapply(r, function(ri) {
    # Marginal likelihood under H1 (integrated over rho)
    lh1 <- integrate(int, ri = ri, lower = bound[1], upper = bound[2],
                     stop.on.error = FALSE, rel.tol = 1e-4)$value / normalization
    # Likelihood under H0 (fixed rho = h0)
    lh0 <- d_cor(ri, h0, n)
    # Bayes factor
    lh1 / lh0
  })
  
  return(x)
}

r_BF_bound_10 <-function(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  y <- numeric(0)
  Bound_finding <-function(r)r_BF10(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)- D
  
  x <- tryCatch(uniroot(Bound_finding, lower = -.99, upper = h0,tol = 1e-5)$root, error = function(e) NA)
  y <- tryCatch(uniroot(Bound_finding, lower =  h0, upper = .99,tol = 1e-5)$root, error = function(e) NA)
  results <- c(x, y)
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")
  
  BF.vals  <- r_BF10(results,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

r_BF_bound_01 <-function(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  r_BF_bound_10(1/D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
}

p_cor<-function(limit,rho,n,lower.tail){
  
  pnorm(r_mean(limit),r_mean(rho),sd = r_sd(n),lower.tail =  lower.tail)
  
  
}

r_TPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)
  
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
  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Normal" = pnorm(bound[2],location,scale)-pnorm(bound[1],location,scale),
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,min(bound),max(bound))-p_beta(bound[1], 1/k,1/k,min(bound),max(bound)) , 
           "NLP"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "t_dis" = pt((bound[2] - location) / scale, dff, 0) - pt((bound[1] - location) / scale, dff, 0),
           "beta" = p_beta(bound[2], alpha, beta,min(bound),max(bound))-p_beta(bound[1], alpha, beta,min(bound),max(bound)))
  
  int <- function(rho) {
    prob <- switch(hypothesis,
                   "!=" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                     p_cor(min(r), rho, n, lower.tail = TRUE),
                   ">"  = p_cor(r, rho, n, lower.tail = FALSE),
                   "<"  = p_cor(r, rho, n, lower.tail = TRUE)
    )
    
    prob * r_prior(rho, k, location, scale, dff, model, alpha, beta,min(bound),max(bound)) / normalization
  }
  x = integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-4)$value
  return(x)
  
}

r_FNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  
  if (any(r == "bound cannot be found") || length(r) == 0) return(0)
  
  
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
  
  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Normal" = pnorm(bound[2],location,scale)-pnorm(bound[1],location,scale),
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,min(bound),max(bound))-p_beta(bound[1], 1/k,1/k,min(bound),max(bound)) , 
           "NLP"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "t_dis" = pt((bound[2] - location) / scale, dff, 0) - pt((bound[1] - location) / scale, dff, 0),
           "beta" = p_beta(bound[2], alpha, beta,min(bound),max(bound))-p_beta(bound[1], alpha, beta,min(bound),max(bound)))
  
  int <- function(rho) {
    prob <- switch(hypothesis,
                   "!=" = p_cor(max(r), rho, n, lower.tail = TRUE) -
                     p_cor(min(r), rho, n, lower.tail = TRUE),
                   ">"  = p_cor(r, rho, n, lower.tail = TRUE),
                   "<"  = p_cor(r, rho, n, lower.tail = FALSE)
    )
    
    prob * r_prior(rho, k, location, scale, dff, model, alpha, beta,min(bound),max(bound)) / normalization
  }
  
  
  x = integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-8, subdivisions=10000000)$value
  return(x)
  
}

r_FPE <-function(r,n,h0,hypothesis){
  
  if (any(r == "bound cannot be found") || length(r) == 0) return(0)
  
  x <- switch(hypothesis,
              "!=" = p_cor(max(r), h0, n, lower.tail = FALSE) +
                p_cor(min(r), h0, n, lower.tail = TRUE),
              ">"  = p_cor(r, h0, n, lower.tail = FALSE),
              "<"  = p_cor(r, h0, n, lower.tail = TRUE)
  )
  return(x)
  
}


r_TNE <-function(r,n,h0,hypothesis){
  
  if (any(r == "bound cannot be found") || length(r) == 0) return(0)
  
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  
  x <- switch(hypothesis,
              "!=" = p_cor(max(r), h0, n, lower.tail = TRUE) -
                p_cor(min(r), h0, n, lower.tail = TRUE),
              ">"  = p_cor(r, h0, n, lower.tail = TRUE),
              "<"  = p_cor(r, h0, n, lower.tail = FALSE)
  )
  
  return(x)
  
}



r_N_finder<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                       location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP){

  lo = 10
  upper = 5000

  r = r_BF_bound_10(D,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  TPE_lo <- if (de_an_prior == 1)
    r_TPE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model) else
      r_TPE(r,lo,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
  FPE_lo <-  r_FPE(r,lo,h0,hypothesis )
   
  if (TPE_lo > target&FPE_lo<FP) return(lo)

  Power_root <- function(N) {
    r <- r_BF_bound_10(D, N, k, alpha, beta, h0, hypothesis, location, scale, dff, model)
    pro <- if (de_an_prior==0){ r_TPE(r, N, k_d, alpha_d, beta_d, h0, hypothesis, location_d, scale_d, dff_d, model_d) }else r_TPE(r, N, k, alpha, beta, h0, hypothesis, location, scale, dff, model)

    pro - target
  }

  N.power = uniroot(Power_root,lower = lo,upper = upper)$root
  
  ## checking if the N lead to an acceptable alpha level
  r = r_BF_bound_10(D,N.power,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  
  FPE = r_FPE(r,N.power,h0,hypothesis)
  if (FPE <= FP) return(N.power)
  
  alpha.root <- function(n) {
    r <- r_BF_bound_10(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
   r_FPE(r,n,h0,hypothesis)-FP
  }
  N.alpha = uniroot(alpha.root,lower = N.power,upper = upper)$root
  return(N.alpha)
  }


r_table<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                    location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,FP ){
  
  if (mode_bf == 1) n = ceiling(r_N_finder(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                           location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP)) else  n = N

  # r bounds:
  r10 <- r_BF_bound_10(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  r01 <-  r_BF_bound_01(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  
  # max BF10 possible:
  max_BF <- 1 / r_BF10(h0,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  BF_D   <- r10
  
  # FPE and TPE:
  FPE       <- r_FPE(r10,n,h0,hypothesis)
  if (de_an_prior == 1) {
    TPE         <- r_TPE(r10, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model)
    TPR_model   <- model
    TPR_k       <- k
    TPR_alpha   <- alpha
    TPR_beta    <- beta
    TPR_location<- location
    TPR_scale   <- scale
    
  } else {
    TPE       <- r_TPE(r10, n, k_d, alpha_d, beta_d, h0, hypothesis, location_d, scale_d, dff_d, model_d)
    TPR_model   <- model_d
    TPR_k       <- k_d
    TPR_alpha   <- alpha_d
    TPR_beta    <- beta_d
    TPR_location<- location_d
    TPR_scale   <- scale_d
  }
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- r_FNE(r01,n,TPR_k, TPR_alpha, TPR_beta,h0,hypothesis,TPR_location,TPR_scale,TPR_dff,TPR_model)
    TNE <- r_TNE(r01,n,h0,hypothesis)
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

compute.prior.density.r <- function(rho, k,location,scale,dff,model, alpha, beta,hypothesis) {
  if (model == "Point") return(rep(NA, length(rho)))
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = -1, b = location),
                   "!=" = c(a = -1, b = 1)
  )
  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Normal" = pnorm(bound[2],location,scale)-pnorm(bound[1],location,scale),
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,min(bound),max(bound))-p_beta(bound[1], 1/k,1/k,min(bound),max(bound)) , 
           "NLP"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "t_dis" = pt((bound[2] - location) / scale, dff, 0) - pt((bound[1] - location) / scale, dff, 0),
           "beta" = p_beta(bound[2], alpha, beta,min(bound),max(bound))-p_beta(bound[1], alpha, beta,min(bound),max(bound)))
  
  
  r_prior(rho,k,location,scale,dff,model, alpha, beta,min(bound),max(bound)) / normalization
}


r_prior_plot <-function(k, alpha, beta,h0,location,scale,dff,model,de_an_prior,
                        k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d,hypothesis){
  par(mfrow = c(1, 1))
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  rho <- seq(bound[1],bound[2],.01)
  prior.analysis <-compute.prior.density.r(rho, k,location,scale,dff,model, alpha, beta,hypothesis)
  prior.design   <- if (de_an_prior == 0 && model_d != "Point")
    compute.prior.density.r(rho, k_d,location_d,scale_d,dff_d,model_d, alpha_d, beta_d,hypothesis) else
      rep(NA, length(rho))
  # Combine all values into one vector
  all_vals <- c(prior.analysis,prior.design)
  
  # Filter out NA and infinite values
  finite_vals <- all_vals[is.finite(all_vals)]
  
  # Get the max from finite values only
  ylim.max <- max(finite_vals)
  # Base plot:
  plot(rho, prior.analysis, type = "l", lwd = 2,
       xlab = expression(bold(delta)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~rho~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))
  
  # If design prior != analysis prior:
  if (de_an_prior == 0) {
    if (model_d == "Point")
      arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim.max, length = 0.1, col = "black", lty = 2) else
        lines(rho, prior.design, lty = 2)
    
    # Add legend:
    legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }

}


r_bf10_p <-function(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model){
  
  rr  <- seq(from = -.99,to = .99,.01)
  
  # Compute BF10 and t-bounds:
  r.BF10 <- r_BF_bound_10(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model)
  BF10 <- r_BF10(rr,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model)
  
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
  r.BF01 <- r_BF_bound_01(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model)
  
  # Check if BF01 = D is possible:
  max.BF01   <- 1 / r_BF10(h0,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
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

Power_r<-function(D,k, alpha, beta,h0,hypothesis,location,scale,dff,model, 
                  k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d, de_an_prior,N){
  
  # N range to evaluate power:
  N.min     <- 4
  N.max     <- ceiling(N * 1.2)
  Ns        <- seq(N.min, N.max, length.out = 31)
  TPE <- numeric(length(Ns))
  FPE <- numeric(length(Ns))
  TNE <- numeric(length(Ns))
  FNE <- numeric(length(Ns))
  

  
  for (i in seq_along(Ns)) {
    r10 <- r_BF_bound_10(D,Ns[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    r01 <- r_BF_bound_01(D,Ns[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    
    # Choose correct design prior:
    TPE[i] <- if (de_an_prior == 1)
      r_TPE(r10,Ns[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model) else
        r_TPE(r10,Ns[i],k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
    FPE[i] <-  r_FPE(r10,Ns[i],h0,hypothesis)
    FNE[i] <- if (de_an_prior == 1)
      r_FNE(r01,Ns[i],k, alpha, beta,h0,hypothesis,location,scale,dff,model) else
        r_FNE(r01,Ns[i],k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
    TNE[i] <- r_TNE(r01,Ns[i],h0,hypothesis)
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

