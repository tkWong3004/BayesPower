Fe_BF <- function(f, q, m, dff, rscale, f_m, model, e) {
  # Compute normalizations once
  normalizationh1  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = e,upper = Inf,rel.tol = 1e-10)$value
  normalizationh0 <- 1 - normalizationh1
  
  # Define likelihood ratio function
  sapply(f, function(fi) {
    int1 <- function(fsq) {
      df(fi, q, m - q, ncp = m * fsq) * F_prior(fsq,q,dff,rscale,f_m,model) / normalizationh1
    }
    int0 <- function(fsq) {
      df(fi, q, m - q, ncp = m * fsq) * F_prior(fsq,q,dff,rscale,f_m,model) / normalizationh0
    }
    lh1 <- integrate(int1, lower = e, upper = Inf, stop.on.error = FALSE)$value
    lh0 <- integrate(int0, lower = 0, upper = e,   stop.on.error = FALSE)$value
    lh1 / lh0
  })
}

Fe_BF_bound_10 <-function(D,q,m,dff,rscale,f_m,model,e){
  x = numeric(0)
  Bound_finding <-function(f){
    Fe_BF(f,q,m,dff,rscale,f_m,model,e)-D
  }
  x = tryCatch( uniroot(Bound_finding,lower=0.01,upper = 500 )$root, error=function(e){})
  if (length(x) == 0) return("no bound is found")
  
  return(x)
}

Fe_BF_bound_01 <-function(D,q,m,dff,rscale,f_m,model,e){
  Fe_BF_bound_10(1/D,q,m,dff,rscale,f_m,model,e)
}


Fe_TPE<-function(f,q,m,dff,rscale,f_m,model,e){
  
  
  if (length(f) == 0 || any(f == "no bound is found")) return(0)
  
  
  if (model == "Point") return(pf(f,q,m-q,ncp =m*f_m,lower.tail = F))
  
  normalizationh1  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = e,upper = Inf,rel.tol = 1e-5)$value
  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh1
  }
  x = integrate(int,lower = e,upper = Inf)$value
    return(x)
}

Fe_FNE<-function(f,q,m,dff,rscale,f_m,model,e){
  
  if (length(f) == 0 || any(f == "no bound is found")) return(0)
  
  
  if (model == "Point") pf(f,q,m-q,ncp =m*f_m,lower.tail = T)
  
  normalizationh1  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = e,upper = Inf,rel.tol = 1e-10)$value
  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh1
  }
  x = integrate(int,lower = e,upper = Inf)$value
  return(x)
}

Fe_TNE<-function(f,q,m,dff,rscale,f_m,model,e){
  
  if (length(f) == 0 || any(f == "no bound is found")) return(0)
  
  normalizationh0  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = 0,upper = e,rel.tol = 1e-5)$value

  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh0
  }
  x = integrate(int,lower = 0,upper = e)$value
  return(x)
}

Fe_FPE<-function(f,q,m,dff,rscale,f_m,model,e){
  
  if (length(f) == 0 || any(f == "no bound is found")) return(0)
  
  normalizationh0  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = 0,upper = e,rel.tol = 1e-10)$value
  

  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh0
  }
  x = integrate(int,lower = 0,upper = e)$value
  return(x)
}

fe_N_finder<-function(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,de_an_prior,e,FP ){
  q     <- k-p
  lower <- 2*k-p+1
  m     <- lower-p
  f     <- Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
  
  p2    <- if (de_an_prior == 1)
    Fe_TPE(f,q,m,dff,rscale,f_m,model,e) else
      Fe_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d,e)
  if (p2 > target) return(lower)
  
  Power_root <- function(n){
    m   <- n-p
    f   <- Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
    pro <- if (de_an_prior == 1)
      Fe_TPE(f,q,m,dff,rscale,f_m,model,e) else
        Fe_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d,e)
    return(pro-target)
  }
  
  #N.power <- uniroot(Power_root,lower = lower,upper =  5000)$root
  N.power <-robust_uniroot(Power_root,lower=lower)
  m       <- N.power-p
  f       <- Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
  FPE     <- Fe_FPE(f,q,m,dff,rscale,f_m,model,e)
  
  if (FPE <= FP) return(N.power)
    
    alpha_root <- function(n){
      m= n-p
      f = Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
      pro = Fe_FPE(f,q,m,dff,rscale,f_m,model,e)
      return(pro-FP)
    }
  N.alpha <- robust_uniroot(alpha_root,lower=lower)
    #uniroot(alpha_root,lower = N.power,upper =  5000)$root
    return(N.alpha)
}

fe_table<-function(D,target,p,k,dff,rscale,f_m,model,
                  dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, mode_bf,e ,FP){
  
  if (mode_bf == 1){
    
    n = ceiling(fe_N_finder(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,de_an_prior,e ,FP))
  } else {
    n=n
  }
  q= k-p
  m= n-p
  
  # f bounds:
  f10 <- Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
  f01 <- Fe_BF_bound_01(D,q,m,dff,rscale,f_m,model,e)
  
  # max BF10 possible:
  max_BF <- 1/Fe_BF(0.00001,q,m,dff,rscale,f_m,model,e)
  BF_D   <- f10
  
  # FPE and TPE:
  FPE       <- Fe_FPE(f10,q,m,dff,rscale,f_m,model,e)
  if (de_an_prior == 1) {
    TPE       <- Fe_TPE(f10,q,m,dff,rscale,f_m,model,e)
    TPR_dff   <- dff
    TPR_rscale<- rscale
    TPR_f_m   <- f_m
    TPR_model <- model
  } else {
    TPE       <- Fe_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d,e)
    TPR_dff   <- dff_d
    TPR_rscale<- rscale_d
    TPR_f_m   <- f_m_d
    TPR_model <- model_d
  }
  
  # FNE and TNE:
  if (max_BF < D | BF_D == "no bound is found") {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- Fe_FNE(f01,q,m,TPR_dff,TPR_rscale,TPR_f_m,TPR_model,e)
    TNE <-  Fe_TNE(f01,q,m,dff,rscale,f_m,model,e)
    
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




prior_plot_fe <-function(q,dff,rscale,f,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior ,e){
  par(mfrow = c(1, 1))
  
  normalization  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f,model),lower = e,upper = Inf,rel.tol = 1e-10)$value

  fsq = seq(0.001,20,.05)
  #prior.analysis.h1 = F_prior(fsq,q,dff,rscale,f,model)/normalization
  prior.analysis.h1 = F_prior(fsq,q,dff,rscale,f,model)
  prior.analysis.h1[fsq<e]=0
  #prior.analysis.h0 = F_prior(fsq,q,dff,rscale,f,model)/(1-normalization)
  prior.analysis.h0 = F_prior(fsq,q,dff,rscale,f,model)
  prior.analysis.h0[fsq>e]=0
  prior.design <- if (de_an_prior == 0 && model_d != "Point"){
    #normalization_d  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f,model),lower = e,upper = Inf,rel.tol = 1e-10)$value
    
    F_prior(fsq,q,dff_d,rscale_d,f_m_d,model_d)/normalization_d
    F_prior(fsq,q,dff_d,rscale_d,f_m_d,model_d)}else
      rep(NA, length(fsq))
  ylim.max <- max(prior.analysis.h1, prior.analysis.h0, prior.design, na.rm = TRUE)
  
  plot(fsq, prior.analysis.h1, type = "l", lwd = 2,
       xlab = bquote(bold(lambda^2)),
       ylab = "density",
       main  = bquote(bold("prior distribution on "~lambda^2~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))
  
  lines(fsq, prior.analysis.h0, lty = 2, col = "black", lwd = 2)
  
  # Optional: design prior
  legend.labels <- c("H1 - Analysis Prior", "H0 - Analysis Prior")
  legend.cols   <- c("black", "black")
  legend.lty    <- c(1, 2)
  legend.lwd    <- c(2, 2)
  
  if (de_an_prior == 0) {
    if (model_d == "Point") {
      arrows(x0 = f_m_d, y0 = 0, x1 = f_m_d, y1 = ylim.max,
             length = 0.1, col = "gray", lty = 2,lwd = 3)
    } else {
      lines(fsq, prior.design, lty = 1, col = "gray", lwd = 3)
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

bf10_fe <-function(D,n,k,p,dff,rscale,f_m,model,e){
  q       <- k-p
  m       <- n-p
  ff      <- seq(from = .01,to = 10,.05)
  f.BF10  <- Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
  BF10    <- Fe_BF(ff,q,m,dff,rscale,f_m,model,e) 
  
  if (length(f.BF10) == 1){
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when f = "~.(format(f.BF10, digits = 4))))
  } else {
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when f = "~.(format(f.BF10[1], digits = 4))~"or"~.(format(f.BF10[2], digits = 4))))
  }
  par(mfrow = c(1, 2))
  plot(ff,BF10,log = "y",xlim=c(0,10),xlab= "f-value",type="l", ylab = expression("BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
  abline(v = f.BF10)
  axis(1, c(0,10))
  
  if (length(f.BF10) != 0 ){
    axis(1, round(f.BF10,2))}
  
  max_BF <- 1/  Fe_BF(.001,q,m,dff,rscale,f_m,model,e)
  f.BF01   <-  Fe_BF_bound_01(D,q,m,dff,rscale,f_m,model,e)
  
  
  
  plot(ff,1/BF10,log = "y",xlab= "f-value",type="l",main = "",frame.plot = FALSE,ylab = bquote("BF"[0][1]),xaxt = "n")
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

Power_fe<-function(D,k,p,dff,rscale,f_m,model,k_d,p_d,dff_d,rscale_d,f_m_d,model_d,de_an_prior,N,e){
  par(mfrow = c(1, 1))
  smin = (2*k-p+1)
  smax = N*2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  
  q = k-p
  m = sdf-p
  TPE =  array(NA, dim = c(length(sdf)))
  FPE =  array(NA, dim = c(length(sdf)))
  TNE =  array(NA, dim = c(length(sdf)))
  FNE =  array(NA, dim = c(length(sdf)))
  for ( i in 1:length(sdf)){
    f10 = Fe_BF_bound_10(D,q,m[i],dff,rscale,f_m,model,e)
    f01 = Fe_BF_bound_10(D,q,m[i],dff,rscale,f_m,model,e)
    
    
    if(de_an_prior == 1){
      TPE[i] = Fe_TPE(f10,q,m[i],dff,rscale,f_m,model,e)
      
    }else{
      TPE[i] = Fe_TPE(f10,q,m[i],dff_d,rscale_d,f_m_d,model_d,e)
    }
    
    if(de_an_prior == 1){
      FNE[i] = Fe_FNE(f01,q,m[i],dff,rscale,f_m,model,e)
      
    }else{
      FNE[i] = Fe_FNE(f01,q,m[i],dff_d,rscale_d,f_m_d,model_d,e)
    }
    
    FPE[i]       <- Fe_FPE(f10,q,m[i],dff,rscale,f_m,model,e)
    TNE[i]       <- Fe_TNE(f01,q,m[i],dff,rscale,f_m,model,e)
    
    
  }
  par(mfrow = c(1, 2))
  plot(sdf, TPE, type = "l",
       xlab = "Total sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(sdf,FPE,col = "grey")
  legend(x = smax*.4,y=.5,              # position of the legend
         legend = c("True positive", "False positive"),  # labels
         col = c("black", "grey"),   # colors matching lines
         lty = 1,                   # line type (solid)
         bty = "n")  
  
  plot(sdf, TNE, type = "l",
       xlab = "Total sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[0][1]~">"~.(D))))
  lines(sdf,FNE,col = "grey")
  legend(x = smax*.4,y=.5,              # position of the legend
         legend = c("True negative", "False negative"),  # labels
         col = c("black", "grey"),   # colors matching lines
         lty = 1,                   # line type (solid)
         bty = "n")
  
}

