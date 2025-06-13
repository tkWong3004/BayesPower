
# k = number of predictor in the full model
# p = number  of predictor in the reduced model
# m = N-p 
# q = k-p

F_prior<- function(fsq,q,dff,rscale,f,model) {
  
  
  switch(model,
         "tdis" = {gamma((q + dff) / 2) / gamma(dff / 2) /gamma(q / 2) *
    (dff * rscale^2)^(dff / 2) * fsq^(q / 2 - 1) *
    (dff * rscale^2 + f^2 + fsq)^(-dff / 2 - q / 2) *
    hypergeo::genhypergeo(c((dff + q) / 4, (2 + dff + q) / 4),
                          q / 2, 4 * f^2 * fsq / (dff * rscale^2 + f^2 + fsq)^2)},
          "Moment" = { temp <- f^2 * (dff + q - 2)/2
          
          gamma((q + dff) / 2) / gamma(dff / 2) / gamma(q / 2) *
            2 * (dff - 2) / q / (dff-2 + q) / f^2 *
            fsq^(q/2) * temp^(dff/2) * (temp + fsq)^(-(dff+q)/2)})
  
}


F_BF <- function(f, q, m, dff, rscale, f_m, model) {
  sapply(f, function(fi) {
    int <- function(fsq) {
      df(fi, q, m - q, ncp = m * fsq) * F_prior(fsq, q, dff, rscale, f_m, model)
    }
    lh1 <- integrate(int, lower = 0, upper = Inf, stop.on.error = FALSE, rel.tol = 1e-4)$value
    lh0 <- df(fi, q, m - q)
    lh1 / lh0
  })
}


F_BF_bound_10 <-function(D,q,m,dff,rscale,f_m,model){
  x = numeric(0)
  Bound_finding <-function(f){
    F_BF(f,q,m,dff,rscale,f_m,model)-D
  }

  x = tryCatch( uniroot(Bound_finding,lower=0.01,upper = 40 )$root, error=function(e){})
  if (length(x) == 0) return("no bound is found")

  return(x)
}

F_BF_bound_01 <-function(D,q,m,dff,rscale,f_m,model){
  F_BF_bound_10(1/D,q,m,dff,rscale,f_m,model)
}


F_TPE<-function(f,q,m,dff,rscale,f_m,model){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)
  
  if (model == "Point"){
    x = pf(f,q,m-q,ncp =m*f_m,lower.tail = F)
    return(x)
  }
  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,model)
  }
  x = integrate(int,lower = 0,upper = Inf)$value
    return(x)
}

F_FNE<-function(f,q,m,dff,rscale,f_m,model){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)
  
  if (model == "Point"){
    
    x = pf(f,q,m-q,ncp =m*f_m,lower.tail = T)
    return(x)
  }
  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,model)
  }
  x = integrate(int,lower = 0,upper = Inf)$value
  return(x)
}

F_TNE<-function(f,q,m){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)
    x = pf(f,q,m-q,ncp =0,lower.tail = T)
  return(x)
}

F_FPE<-function(f,q,m){
  
  if (length(f) == 0 || any(f == "no bound is found")) return(0)
  
  x = pf(f,q,m-q,ncp =0,lower.tail = F)
  return(x)
}

f_N_finder<-function(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,FP){
  q= k-p
  lower = 2*k-p+1
  m= lower-p
  f <- F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  p2 <- if (de_an_prior == 1)
    F_TPE(f,q,m,dff,rscale,f_m,model) else
      F_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d)
  if (p2 > target) return(lower)

  Power_root <- function(n){
    m= n-p
    f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
    
    pro <- if (de_an_prior == 1)
      F_TPE(f,q,m,dff,rscale,f_m,model) else
        F_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d)
    
    return(pro-target)
  }
  
  N.power = uniroot(Power_root,lower = lower,upper =  10000)$root
  m= N.power-p
  f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  FPE = F_FPE(f,q,m)
  
  if (FPE <= FP) return(N.power + 1)
  alpha_root <- function(n){
    m= n-p
    f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
    pro = F_FPE(f,q,m)
    
    
    return(pro-FP)
  }
  N.alpha = uniroot(alpha_root,lower = N.power,upper =  10000)$root
  return(N.alpha)
}

f_table<-function(D,target,p,k,dff,rscale,f_m,model,
                  dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, mode_bf,FP ){

  
  if (mode_bf == 1){
    
    n = ceiling(f_N_finder(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,FP ))
  } else {
    n=n
  }
  q= k-p
  m= n-p
  
  # f bounds:
  f10 <- F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  f01 <- F_BF_bound_01(D,q,m,dff,rscale,f_m,model)
  
  # max BF10 possible:
  max_BF <- 1/F_BF(0.00001,q,m,dff,rscale,f_m,model)
  BF_D   <- f10
  
  # FPE and TPE:
  FPE       <- F_FPE(f10,q,m)
  if (de_an_prior == 1) {
    TPE       <- F_TPE(f10,q,m,dff,rscale,f_m,model)
    TPR_dff   <- dff
    TPR_rscale<- rscale
    TPR_f_m   <- f_m
    TPR_model <- model
  } else {
    TPE       <- F_TPE(f10,q,m,dff_d,rscale_d,f_m_d,model_d)
    TPR_dff   <- dff_d
    TPR_rscale<- rscale_d
    TPR_f_m   <- f_m_d
    TPR_model <- model_d
  }
  
  # FNE and TNE:
  if ( max_BF < D | BF_D == "no bound is found") {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- F_FNE(f01,q,m,TPR_dff,TPR_rscale,TPR_f_m,TPR_model)
    TNE <- F_TNE(f01,q,m)
     
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

prior_plot_f <-function(q,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior ){
  par(mfrow = c(1, 1))
  fsq = seq(0.001,15,.2)
  
  prior.analysis = F_prior(fsq,q,dff,rscale,f_m,model)
  prior.design   <- if (de_an_prior == 0 && model_d != "Point")
    F_prior(fsq,q,dff_d,rscale_d,f_m_d,model_d) else
      rep(NA, length(fsq))
  ylim.max <- max(prior.analysis, prior.design, na.rm = TRUE)
  
  
  # Base plot:
  plot(fsq, prior.analysis, type = "l", lwd = 2,
       xlab = expression(bold(lambda^2)),
       ylab = "density",
       main = bquote(bold("prior distribution on "~lambda^2~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))
  
  # If design prior != analysis prior:
  if (de_an_prior == 0) {
    if (model_d == "Point")
      arrows(x0 = f_m, y0 = 0, x1 = f_m, y1 = ylim.max, length = 0.1, col = "black", lty = 2) else
        lines(fsq, prior.design, lty = 2)
    
    # Add legend:
    legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }
  
}

bf10_f <-function(D,n,k,p,dff,rscale,f_m,model){
  q       <- k-p
  m       <- n-p
  ff      <- seq(from = .01,to = 10,.05)
  f.BF10  <- F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  BF10    <- F_BF(ff,q,m,dff,rscale,f_m,model)
  
  
  if (length(f.BF10) == 1){
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when f = "~.(format(f.BF10, digits = 4))))
  } else {
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when f = "~.(format(f.BF10[1], digits = 4))~"or"~.(format(f.BF10[2], digits = 4))))
  }
  
  par(mfrow = c(1, 2))
  plot(ff,BF10,log = "y",xlim=c(0,10),xlab= "f-value",type="l", ylab = expression("logarithm of BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
  abline(v = f.BF10)
  axis(1, c(0,10))
  
  if (length(f.BF10) != 0 ){
    axis(1, round(f.BF10,2))}
  
  max_BF = 1/  F_BF(.001,q,m,dff,rscale,f_m,model)
  f.BF01 =  F_BF_bound_01(D,q,m,dff,rscale,f_m,model)
  
  
  
  plot(ff,1/BF10,log = "y",xlab= "f-value",type="l",main = "",frame.plot = FALSE,ylab = bquote("logarithm of BF"[0][1]),xaxt = "n")
  axis(1, c(0,10))
  if (any(max_BF<D |f.BF01 == "bound cannot be found" ) ) {
    main = bquote(bold("It is impossible to have BF"[0][1]~"="~.(D)))
    title(main = main)
    #sprintf("It is impossible to have BF01 = %.3f ",D)
  } else      {
    abline(v = f.BF01)
    axis(1, round(f.BF01,2))
    if (length(f.BF01) == 1){
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when f = "~.(format(f.BF01, digits = 4))))
      title(main = main)
    } else {
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when f = "~.(format(f.BF01[1], digits = 4))~"or"~.(format(f.BF01[2], digits = 4))))
      title(main = main)
    }}
  
}
Power_f<-function(D,k,p,dff,rscale,f_m,model,k_d,p_d,dff_d,rscale_d,f_m_d,model_d,de_an_prior,N){
  smin = (2*k-p+1)
  smax = N*1.2
  n = seq(smin,smax , by = (smax-smin)/30)
  q = k-p
  m = n-p
  TPE =  array(NA, dim = c(length(n)))
  FPE =  array(NA, dim = c(length(n)))
  TNE =  array(NA, dim = c(length(n)))
  FNE =  array(NA, dim = c(length(n)))
  
  for ( i in 1:length(n)){
    f10 = F_BF_bound_10(D,q,m[i],dff,rscale,f_m,model)
    f01 = F_BF_bound_01(D,q,m[i],dff,rscale,f_m,model)
    
   
    if(de_an_prior == 1){
      TPE[i] = F_TPE(f10,q,m[i],dff,rscale,f_m,model)
    }else{
      TPE[i] = F_TPE(f10,q,m[i],dff_d,rscale_d,f_m_d,model_d)
    }

    FPE[i] = F_FPE(f10,q,m[i])
    TNE[i] = F_TNE(f01,q,m[i])
    if(de_an_prior == 1){
      FNE[i] = F_FNE(f01,q,m[i],dff,rscale,f_m,model)
    }else{
      FNE[i] = F_FNE(f01,q,m[i],dff_d,rscale_d,f_m_d,model_d)
    }
    
  }

  par(mfrow = c(1, 2))
  plot(n, TPE, type = "l",
       xlab = "Total sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(n,FPE,col = "grey")
  legend(x = smax*.4,y=.5,              # position of the legend
         legend = c("True positive", "False positive"),  # labels
         col = c("black", "grey"),   # colors matching lines
         lty = 1,                   # line type (solid)
         bty = "n")  
  
  plot(n, TNE, type = "l",
       xlab = "Total sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[0][1]~">"~.(D))))
  lines(n,FNE,col = "grey")
  legend(x = smax*.4,y=.5,              # position of the legend
         legend = c("True negative", "False negative"),  # labels
         col = c("black", "grey"),   # colors matching lines
         lty = 1,                   # line type (solid)
         bty = "n")   
  
}
