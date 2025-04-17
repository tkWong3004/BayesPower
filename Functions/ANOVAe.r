library(rootSolve)

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

Fe_BF <-function(f,q,m,dff,rscale,f_m,model,e){

  
  normalizationh1  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = e,upper = Inf,rel.tol = 1e-10)$value
  normalizationh0 = 1-normalizationh1
  x = NA
  for(i in 1:length(f)){
  
  int1  <- function(fsq){df(f[i],q,m-q,ncp =m*fsq)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh1
   
  }
  int0  <- function(fsq){df(f[i],q,m-q,ncp =m*fsq)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh0
    
  }
  lh1 = integrate(int1,lower = e,upper = Inf,stop.on.error = F)$value
  lh0 = integrate(int0,lower = 0,upper = e,stop.on.error = F)$value
  x[i] = lh1/lh0
  }
  return(x)
}


Fe_BF_bound_10 <-function(D,q,m,dff,rscale,f_m,model,e){
  x = numeric(0)
  Bound_finding <-function(f){
    Fe_BF(f,q,m,dff,rscale,f_m,model,e)-D
  }
 
  x = tryCatch( uniroot(Bound_finding,lower=0.001,upper = 200 )$root, error=function(e){})
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  

  return(x)
}

Fe_BF_bound_01 <-function(D,q,m,dff,rscale,f_m,model,e){
  Bound_finding <-function(f){
    Fe_BF(f,q,m,dff,rscale,f_m,model,e)-1/D
  }
  
  x = tryCatch( uniroot(Bound_finding,lower=0.001,upper = 200 )$root, error=function(e){})
  
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  
  return(x)
}


Fe_TPE<-function(f,q,m,dff,rscale,f_m,model,e){
  
  
  if (any(f =="no bound is found" | length(f)==0)){
    t=0
    return(t)
  }
  
  if (model == "Point"){
    
    x = pf(f,q,m-q,ncp =m*f_m,lower.tail = F)
    return(x)
  }
  
  normalizationh1  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = e,upper = Inf,rel.tol = 1e-10)$value
  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh1
  }
  x = integrate(int,lower = e,upper = Inf)$value
    return(x)
}

Fe_FNE<-function(f,q,m,dff,rscale,f_m,model,e){
  
  
  if (any(f =="no bound is found" | length(f)==0)){
    t=0
    return(t)
  }
  
  if (model == "Point"){
    
    x = pf(f,q,m-q,ncp =m*f_m,lower.tail = T)
    return(x)
  }
  
  normalizationh1  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = e,upper = Inf,rel.tol = 1e-10)$value
  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh1
  }
  x = integrate(int,lower = e,upper = Inf)$value
  return(x)
}

Fe_TNE<-function(f,q,m,dff,rscale,f_m,model,e){
  normalizationh0  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = 0,upper = e,rel.tol = 1e-10)$value
  
  if (any(f =="no bound is found" | length(f)==0)){
    t=0
    return(t)
  }
  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh0
  }
  x = integrate(int,lower = 0,upper = e)$value
  return(x)
}

Fe_FPE<-function(f,q,m,dff,rscale,f_m,model,e){
  normalizationh0  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = 0,upper = e,rel.tol = 1e-10)$value
  
  if (any(f =="no bound is found" | length(f)==0)){
    t=0
    return(t)
  }
  int  <- function(fsq){
    
    pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh0
  }
  x = integrate(int,lower = 0,upper = e)$value
  return(x)
}

fe_N_finder<-function(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,de_an_prior,e,FP ){
  q= k-p
  lo = 2*k-p+1
  m= lo-p
  f = Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
  if ( de_an_prior ==1 ){
    power = Fe_TPE(f,q,m,dff,rscale,f_m,model,e)
  }else{
    power = Fe_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d,e)
      
  }  
  if(any(power >target)){
    return(lo)
  }
  
  
  
  Power_root <- function(n){
    m= n-p
    f = Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
    if ( de_an_prior ==1 ){
      pro = Fe_TPE(f,q,m,dff,rscale,f_m,model,e)
    }else{
      pro = Fe_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d,e)
      
    }  
    return(pro-target)
  }
  
  N = uniroot(Power_root,lower = lo,upper =  5000)$root
  m= N-p
  f = Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
  FPE = Fe_FPE(f,q,m,dff,rscale,f_m,model,e)
  
  
  if (FP>FPE){
    return(N)
  } else{
    
    alpha_root <- function(n){
      m= n-p
      f = Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
      pro = Fe_FPE(f,q,m,dff,rscale,f_m,model,e)
      
      
      return(pro-FP)
    }
    NN = uniroot(alpha_root,lower = N,upper =  10000)$root
    return(NN)
  }


}

fe_table<-function(D,target,p,k,dff,rscale,f_m,model,
                  dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, mode_bf,e ,FP){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  
  if (mode_bf == 1){
    
    n = fe_N_finder(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,de_an_prior,e ,FP)
  } else {
    n=n
  }
  q= k-p
  m= n-p
  
  max_BF = 1/Fe_BF(0.00001,q,m,dff,rscale,f_m,model,e)
  if (max_BF<D){
    TPE = 0
    FPE = 0
  }
  
  bound01 = Fe_BF_bound_01(D,q,m,dff,rscale,f_m,model,e)
  
  
  if (length(bound01)==0){
    FNE = 0
    TNE = 0
  }else{
    
    if (de_an_prior ==1){
      FNE = Fe_FNE(bound01,q,m,dff,rscale,f_m,model,e)
      TNE = Fe_TNE(bound01,q,m,dff,rscale,f_m,model,e)
    }
    if (de_an_prior ==0){
      FNE = Fe_FNE(bound01,q,m,dff_d,rscale_d,f_m_d,model_d,e)
      TNE = Fe_TNE(bound01,q,m,dff,rscale,f_m,model,e)
    }
    
    
    
  }
  
  bound10 = Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
  
  if (de_an_prior ==1){
    TPE = Fe_TPE(bound10,q,m,dff,rscale,f_m,model,e)
    FPE = Fe_FPE(bound10,q,m,dff,rscale,f_m,model,e)
  }else{
    TPE = Fe_TPE(bound10,q,m,dff_d,rscale_d,f_m_d,model_d,e)
    FPE = Fe_FPE(bound10,q,m,dff,rscale,f_m,model,e)
  }
  
  
  table <- data.frame(
    TPE = TPE,
    FNE = FNE,
    TNE = TNE,
    FPE = FPE,
    N =  ceiling(n))
  return(table)
  
}




prior_plot_fe <-function(q,dff,rscale,f,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior ,e){
  par(mfrow = c(1, 1))
  normalization  <- integrate(function(fsq)F_prior(fsq,q,dff,rscale,f,model),lower = e,upper = Inf,rel.tol = 1e-10)$value
  

  fsq = seq(0.001,20,.05)
  priorh1 = F_prior(fsq,q,dff,rscale,f,model)/normalization
  priorh1[fsq<e]=-1
  priorh0 = F_prior(fsq,q,dff,rscale,f,model)/(1-normalization)
  priorh0[fsq>e]=-1
  

  plot(fsq,priorh0,
       lty = 3,
       lwd = 4,
       xlab= bquote(bold(lambda^2)),
       ylab= "density",
       type = "l",
       ylim = c(0,max(max(priorh1),max(priorh0))),
       main  = bquote(bold("prior distribution on "~lambda^2~" under the alternative hypothesis")),frame.plot = FALSE)
  lines(fsq,priorh1, type="l", lty=1, lwd=4)

  
  if (de_an_prior ==0){
    if (model_d != "Point"){
      
      normalization  <- integrate(function(fsq)F_prior(fsq,q,dff_d,rscale_d,f_m_d,model_d),lower = e,upper = Inf,rel.tol = 1e-10)$value
      
      
      
    prior_d = F_prior(fsq,q,dff_d,rscale_d,f_m_d,model)/normalization
    prior_d[fsq<e]=-1
    
    
    plot(fsq,priorh0,
         lty = 3,
         lwd = 4,
         xlab= bquote(bold(lambda^2)),
         ylab= "density",
         type = "l",
         ylim = c(0,max(max(priorh1),max(prior_d),max(priorh0))),
         main  = bquote(bold("prior distribution on "~lambda^2~" under the alternative hypothesis")),frame.plot = FALSE)
    lines(fsq,priorh1, type="l", lty=1, lwd=4)
    
    
    
  
    lines(fsq,prior_d, type="l", lty=1, lwd=2,col="gray")}else{
      
      arrows(x0=f_m_d, y0=0, x1=f_m_d, y1=max(priorh1,priorh0), length=0.2, code=2, col="gray", lwd=4,lty = 2)
    }

    legend("topleft", 
           legend = c("Analysis prior", "Design prior"), 
           lty = c(1, 1), 
           lwd = c(4, 2),
           col = c("black", "gray"),
           bty = "n") 
    
  }
  
  
  
  legend("topright", 
         legend = c("H1", "H0"),  # Labels for the lines
         col = c("black", "black"),  # Line colors
         lty = c(1, 3),  # Line types
         lwd = c(4, 4),  # Line widths
         bty = "n",
         title = "Analysis Prior")  # No box around the legend
  
  
}

bf10_fe <-function(D,n,k,p,dff,rscale,f_m,model,e){
  q= k-p
  m= n-p
  ff = seq(from = .01,to = 10,.05)
  BF_D = Fe_BF_bound_10(D,q,m,dff,rscale,f_m,model,e)
  BF10 = Fe_BF(ff,q,m,dff,rscale,f_m,model,e)
  
  
  if (length(BF_D) == 1){
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when f = "~.(format(BF_D, digits = 4))))
    #sprintf("BF10 = %.0f when f = %.3f ",D,BF_D)
  } else {
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when f = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
    #sprintf("BF10 = %.0f when f = %.3f or %.3f ",D,BF_D[1],BF_D[2])
  }
  par(mfrow = c(1, 2))
  plot(ff,log10(BF10),xlim=c(0,10),xlab= "f-value",type="l", ylab = expression("logarithm of BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
  abline(v = BF_D)
  axis(1, c(0,10))
  if (length(BF_D) != 0 ){
    axis(1, round(BF_D,2))}
  
  max_BF = 1/  Fe_BF(.001,q,m,dff,rscale,f_m,model,e)
  BF_D =  Fe_BF_bound_01(D,q,m,dff,rscale,f_m,model,e)
  
  
  
  plot(ff,log10(1/BF10),xlab= "f-value",type="l",main = "",frame.plot = FALSE,ylab = bquote("logarithm of BF"[0][1]),xaxt = "n")
  axis(1, c(0,10))
  if (any(max_BF<D |BF_D == "bound cannot be found" ) ) {
    main = bquote(bold("It is impossible to have BF"[0][1]~"="~.(D)))
    title(main = main)
    #sprintf("It is impossible to have BF01 = %.3f ",D)
  } else      {
    abline(v = BF_D)
    axis(1, round(BF_D,2))
    if (length(BF_D) == 1){
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when f = "~.(format(BF_D, digits = 4))))
      title(main = main)
    } else {
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when f = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
      title(main = main)
    }}
  
}

Power_fe<-function(D,k,p,dff,rscale,f_m,model,k_d,p_d,dff_d,rscale_d,f_m_d,model_d,de_an_prior,N,e){
  par(mfrow = c(1, 1))
  smin = (2*k-p+1)
  smax = N*2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  power =  array(NA, dim = c(length(sdf)))
  q = k-p
  m = sdf-p
  for ( i in 1:length(sdf)){
    f = Fe_BF_bound_10(D,q,m[i],dff,rscale,f_m,model,e)
    
    
    if(de_an_prior == 1){
      power[i] = Fe_TPE(f,q,m[i],dff,rscale,f_m,model,e)
      
    }else{
      power[i] = Fe_TPE(f,q,m[i],dff_d,rscale_d,f_m_d,model_d,e)
    }
    
    
    
  }
  plot(sdf+1,power,type="l",main = "",frame.plot = FALSE,xlab = "Sample size", ylab = "Probability of True positive evidence",xlim = c(smin,max(sdf)), 
       ylim = c(0,1) )
  
}

