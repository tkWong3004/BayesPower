library(rootSolve)
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



F_BF <-function(f,q,m,dff,rscale,f_m,model){
  x = NA
  for(i in 1:length(f)){
  
  int  <- function(fsq){df(f[i],q,m-q,ncp =m*fsq)*F_prior(fsq,q,dff,rscale,f_m,model)
 
     }
  lh1 = integrate(int,lower = 0,upper = Inf,stop.on.error = F)$value
  lh0 = df(f[i],q,m-q)
  x[i] = lh1/lh0
  }
  return(x)
}


F_BF_bound_10 <-function(D,q,m,dff,rscale,f_m,model){
  x = numeric(0)
  Bound_finding <-function(f){
    F_BF(f,q,m,dff,rscale,f_m,model)-D
  }
 
  x = tryCatch( uniroot(Bound_finding,lower=0.01,upper = 40 )$root, error=function(e){})
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  

  return(x)
}

F_BF_bound_01 <-function(D,q,m,dff,rscale,f_m,model){
  Bound_finding <-function(f){
    1/F_BF(f,q,m,dff,rscale,f_m,model)-D
  }
  
  x = tryCatch( uniroot(Bound_finding,lower=0.01,upper = 40 )$root, error=function(e){})
  
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  
  return(x)
}


F_TPE<-function(f,q,m,dff,rscale,f_m,model){
  if (any(f =="no bound is found" | length(f)==0)){
    t=0
    return(t)
  }
  
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
  if (any(f =="no bound is found" | length(f)==0)){
    t=0
    return(t)
  }
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

F_TNE<-function(f,q,m,dff,rscale,f_m,model){
  if (any(f =="no bound is found" | length(f)==0)){
    t=0
    return(t)
  }
  
    
    x = pf(f,q,m-q,ncp =0,lower.tail = T)
  return(x)
}

F_FPE<-function(f,q,m,dff,rscale,f_m,model){
  
  if (any(f =="no bound is found" | length(f)==0)){
    t=0
    return(t)
  }
  
  x = pf(f,q,m-q,ncp =0,lower.tail = F)
  return(x)
}

f_N_finder<-function(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,FP){
  q= k-p
  lo = 2*k-p+1
  m= lo-p
  f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
 if ( de_an_prior ==1 ){
   power = F_TPE(f,q,m,dff,rscale,f_m,model)
 }else{
   power = F_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d)
 }
     
  if(any(power >target)){
    return(lo)
  }
  
  Power_root <- function(n){
    m= n-p
    f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
    if(de_an_prior == 1){
      pro = F_TPE(f,q,m,dff,rscale,f_m,model)
      
    }else{
    
    pro = F_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d)}
    
    return(pro-target)
  }
  
  N = uniroot(Power_root,lower = lo,upper =  10000)$root
  m= N-p
  f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  FPE = F_FPE(f,q,m,dff,rscale,f_m,model)
  
  if (FP>FPE){
    return(N)
  } else{
    
    alpha_root <- function(n){
      m= n-p
      f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
      pro = F_FPE(f,q,m,dff,rscale,f_m,model)

      
      return(pro-FP)
    }
    NN = uniroot(alpha_root,lower = N,upper =  10000)$root
    
  }
  
}

f_table<-function(D,target,p,k,dff,rscale,f_m,model,
                  dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, mode_bf,FP ){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  
  if (mode_bf == 1){
    
    n = ceiling(f_N_finder(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,FP ))
  } else {
    n=n
  }
  q= k-p
  m= n-p
  
  max_BF = 1/F_BF(0.00001,q,m,dff,rscale,f_m,model)
  if (max_BF<D){
    TPE = 0
    FPE = 0
  }
  
  bound01 = F_BF_bound_01(D,q,m,dff,rscale,f_m,model)
  
  
  if (length(bound01)==0){
    FNE = 0
    TNE = 0
  }else{
    
    if (de_an_prior ==1){
      FNE = F_FNE(bound01,q,m,dff,rscale,f_m,model)
      TNE = F_TNE(bound01,q,m,dff,rscale,f_m,model)
    }
    if (de_an_prior ==0){
      FNE = F_FNE(bound01,q,m,dff_d,rscale_d,f_m_d,model_d)
      TNE = F_TNE(bound01,q,m,dff,rscale,f_m,model)
    }
    
    
    
  }
  
  bound10 = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  
  if (de_an_prior ==1){
    TPE = F_TPE(bound10,q,m,dff,rscale,f_m,model)
    FPE = F_FPE(bound10,q,m,dff,rscale,f_m,model)
  }else{
    TPE = F_TPE(bound10,q,m,dff_d,rscale_d,f_m_d,model_d)
    FPE = F_FPE(bound10,q,m,dff,rscale,f_m,model)
  }
  
  
  table <- data.frame(
    TPE = TPE,
    FNE = FNE,
    TNE = TNE,
    FPE = FPE,
    N =  ceiling(n))
  return(table)
  
}

prior_plot_f <-function(q,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior ){
  par(mfrow = c(1, 1))

  fsq = seq(0.001,15,.2)
  prior = F_prior(fsq,q,dff,rscale,f_m,model)
  
  
  plot(fsq,prior,xlab= bquote(bold(lambda^2)),
       ylab= "density",type = "l",
       main  = bquote(bold("prior distribution on "~lambda^2~" under the alternative hypothesis")),
       frame.plot = FALSE)
  if (de_an_prior == 0){
    
    if (model_d != "Point"){
    prior_d = F_prior(fsq,q,dff_d,rscale_d,f_m_d,model)
    plot(fsq,prior,ylim=c(0,max(prior_d,prior)),xlab= bquote(bold(lambda^2)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~lambda^2~" under the alternative hypothesis")),frame.plot = FALSE)
    lines(fsq,prior_d ,lty = 2)
    }else{
      
      plot(fsq,prior,xlab= bquote(bold(lambda^2)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~lambda^2~" under the alternative hypothesis")),frame.plot = FALSE)
      arrows(x0=f_m_d, y0=0, x1=f_m_d, y1=max(prior), length=0.2, code=2, col="black", lwd=1,,lty = 2)
    }
    legend("topright", 
           legend = c("Analysis prior", "Design prior"), 
           lty = c(1, 2), 
           col = c("black", "black"),
           bty = "n") 
    
  }
  
  
}

bf10_f <-function(D,n,k,p,dff,rscale,f_m,model){
  q= k-p
  m= n-p
  ff = seq(from = .01,to = 10,.05)
  BF_D = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  BF10 = F_BF(ff,q,m,dff,rscale,f_m,model)
  
  
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
  
  max_BF = 1/  F_BF(.001,q,m,dff,rscale,f_m,model)
  BF_D =  F_BF_bound_01(D,q,m,dff,rscale,f_m,model)
  
  
  
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
Power_f<-function(D,k,p,dff,rscale,f_m,model,k_d,p_d,dff_d,rscale_d,f_m_d,model_d,de_an_prior,N){
  par(mfrow = c(1, 1))
  smin = (2*k-p+1)
  smax = N*1.2
  sdf = seq(smin,smax , by = (smax-smin)/30)
  power =  array(NA, dim = c(length(sdf)))
  q = k-p
  m = sdf-p
  for ( i in 1:length(sdf)){
    f = F_BF_bound_10(D,q,m[i],dff,rscale,f_m,model)
      
   
    if(de_an_prior == 1){
      power[i] = F_TPE(f,q,m[i],dff,rscale,f_m,model)
    }else{
      power[i] = F_TPE(f,q,m[i],dff_d,rscale_d,f_m_d,model_d)
    }

    
    
  }
  plot(sdf+1,power,type="l",main = "",frame.plot = FALSE,xlab = "Sample size", ylab = "Probability of True positive evidence",xlim = c(1,max(sdf)), 
       ylim = c(0,1) )
  
}
