# likelihood of non-local prior
dnlp <-function(delta,mu,ta){
  ((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
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

  normalization  <- integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound[1],upper = bound[2],rel.tol = 1e-10)$value

  
  int  <- function(prop){dbinom(x[i], size=n, prob=prop) *bin_prior(prop,alpha,beta,location,scale,model)}

  for (i in 1:length(x)){
    lh1 =integrate(int,lower = bound[1],upper = bound[2],stop.on.error = F,rel.tol = 1e-10)$value/normalization
    lh0 = dbinom(x[i], size=n, prob=location)
    BF[i] = lh1/lh0
  }
  return(BF)
  
}


bin_BF_bound_10 <-function(D,n,alpha,beta,location,scale,model,hypothesis){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)- D
  }
  switch(hypothesis,
         "!="= {
           
      x <- tryCatch( round(uniroot(Bound_finding,lower = 0 ,upper = round(location*n))$root), error=function(e){})
      y <-   tryCatch(  round(uniroot(Bound_finding,lower = round(location*n) ,upper = n)$root), error=function(e){})
      },
      ">" = {
        x <- tryCatch(  round(uniroot(Bound_finding,lower = round(location*n) ,upper = n)$root), error=function(e){})
      },
      "<" = {
        x <- tryCatch( round(uniroot(Bound_finding,lower = 0 ,upper = round(location*n))$root), error=function(e){})
      } ) 

  
  
  switch(hypothesis,
         "!=" ={
           if (length(x) ==0&length(y)==0){
             x = "no bound is found"
             return(x)
           } 
         },
           "<" = {if (length(x) ==0){
             x = "no bound is found"
             return(x)
           }},
         ">" = {if (length(x) ==0){
           x = "no bound is found"
           return(x)
         }})
  
  
  if (length(y) == 1){
    x = cbind(x,y)
  }
  
  BF = bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)
  
  while (any(D>BF)){
    
    if (length(x) ==2){
      
      x =c(x[1] - 1, x[2] + 1)
    }else{
      
      if (x>=location*n){
        x=x+1
      }else{
        x=x-1
      }
      
    }
    if (all(x == 0)) {  
      break
    }
    BF = bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)
  }  
 
    
    
    return(x)
}


bin_BF_bound_01 <-function(D,n,alpha,beta,location,scale,model,hypothesis){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    1/bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)- D
  }
  
  switch(hypothesis,
         "!="= {
           
           x <- tryCatch( round(uniroot(Bound_finding,lower = 0 ,upper = round(location*n))$root), error=function(e){})
           y <-   tryCatch(  round(uniroot(Bound_finding,lower = round(location*n) ,upper = n)$root), error=function(e){})
         },
         ">" = {
           x <- tryCatch(  round(uniroot(Bound_finding,lower = round(location*n) ,upper = n)$root), error=function(e){})
         },
         "<" = {
           x <- tryCatch( round(uniroot(Bound_finding,lower = 0 ,upper = round(location*n))$root), error=function(e){})
         } ) 
  
  switch(hypothesis,
         "!=" ={
           if (length(x) ==0&length(y)==0){
             x = "no bound is found"
             return(x)
           } 
         },
         "<" = {if (length(x) ==0){
           x = "no bound is found"
           return(x)
         }},
         ">" = {if (length(x) ==0){
           x = "no bound is found"
           return(x)
         }})
  
  
  if (length(y) == 1){
    x = cbind(x,y)
  }
  
  BF = 1/bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)
  
  while (any(D>BF)){
    
    if (length(x) ==2){
      
      x =c(x[1] + 1, x[2] - 1)
    }else{
      if(x>location*n){
        x=x-1
      }else{
        x=x+1
      }
      
    }
    BF = 1/bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)
  }  
  
  return(x)
}


bin_TPE<-function(x,n,alpha,beta,location,scale,model,hypothesis){
  if (any(x =="no bound is found" | length(x)==0)){
    x=0
    return(x)
  }
 
  
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
  normalization  <- integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound[1],upper = bound[2],rel.tol = 1e-10)$value
  
  switch(hypothesis,
         
         "!=" =
           
           int  <- function(prop){
             if (length(x)==2){
               pro = pbinom(min(x),n,prop,lower.tail = T)+ pbinom(max(x)-1,n,prop,lower.tail = F)
             }else{
               
               pro = mapply(function(x_i, n_i, p_i) {
                 if (x_i / n_i > p_i) {
                   pbinom(x_i-1, n_i, p_i, lower.tail = FALSE)
                 } else {
                   pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                 }
               }, x, n, prop)
               #if(any(x/n>prop)){
              # pro = pbinom(x,n,prop,lower.tail = F)
              # }else{
               #  pro = pbinom(x,n,prop,lower.tail = T)
                # }
             }
             
             
             
         #  pro = switch(length(x)==2,
        #                "1" ={pbinom(min(x),n,prop,lower.tail = T)+ pbinom(max(x),n,prop,lower.tail = F)},
        #                "0"=  {    
        #                  switch(x/n>prop,
        #                         "1" = pbinom(x,n,prop,lower.tail = F),
        #                         "0" = pbinom(x,n,prop,lower.tail = T)
        #                         )
                          
        #                })
           
           pro*bin_prior(prop,alpha,beta,location,scale,model)/normalization
           },
         ">" = 
           int  <- function(prop){
             pro = pbinom(x-1,n,prop,lower.tail = F)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalization
           },
         "<" = 
           int  <- function(prop){
             pro = pbinom(x,n,prop,lower.tail = T)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalization
           })
  
  TPE = integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-10)$value
  return(TPE)
  
}


bin_FNE<-function(x,n,alpha,beta,location,scale,model,hypothesis){
  if (any(x =="no bound is found" | length(x)==0)){
    x=0
    return(x)
  }
  
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

  
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = 0, b = location),
                   "!=" = c(a = 0, b = 1)
  )
  normalization  <- integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound[1],upper = bound[2],rel.tol = 1e-10)$value
  switch(hypothesis,
         
         "!=" =
           
           int  <- function(prop){
             pro = switch(length(x)==2,
                          "1" ={pbinom(max(x),n,prop,lower.tail = T)- pbinom(min(x)-1,n,prop,lower.tail = T)},
                          "0"=  {    
                            switch(x/n>prop,
                                   "1" = pbinom(x,n,prop,lower.tail = T),
                                   "0" = pbinom(x-1,n,prop,lower.tail = F))
                            
                          })
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalization
           },
         ">" = 
           int  <- function(prop){
             pro = pbinom(x-1,n,prop,lower.tail = T)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalization
           },
         "<" = 
           int  <- function(prop){
             pro = pbinom(x-1,n,prop,lower.tail = F)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalization
           })
  FNE = integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-8)$value
  return(FNE)
  
}


bin_FPE<-function(x,n,alpha,beta,location,scale,model,hypothesis){
  
  if (any(x =="no bound is found" | length(x)==0)){
    x=0
    return(x)
  }
  

    FPE = switch(hypothesis,
                 "!=" = {
              
                   if (length(x)==2){
                     pbinom(min(x),n,location,lower.tail = T)+ pbinom(max(x)-1,n,location,lower.tail = F)
                   }else{  
                     #switch(x/n>location,
                    #              "1" = pbinom(x,n,location,lower.tail = F),
                    #              "0" = pbinom(x,n,location,lower.tail = T))}
                     mapply(function(x_i, n_i, p_i) {
                       if (x_i / n_i > p_i) {
                         pbinom(x_i-1, n_i, p_i, lower.tail = FALSE)
                       } else {
                         pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                       }
                     }, x, n, location)    
                            
                            
                          
                            
                          }},
                 ">"  = {pbinom(x-1,n,location,lower.tail = F)},
                 "<"  = {pbinom(x,n,location,lower.tail = T)}
    )
    return(FPE)

}


bin_TNE<-function(x,n,alpha,beta,location,scale,model,hypothesis){
  
  if (any(x =="no bound is found" | length(x)==0)){
    x=0
    return(x)
  }
  
  TNE = switch(hypothesis,
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
  return(TNE)
}

bin_N_finder <-function(D,target,alpha,beta,location,scale,model,hypothesis,
                        alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP){
  lo = 10
  x = bin_BF_bound_10(D,lo,alpha,beta,location,scale,model,hypothesis)
  pro = bin_TPE(x,lo,alpha,beta,location,scale,model,hypothesis)
  
  if ( pro>target){

    FPE = bin_FPE(x,lo,alpha,beta,location,scale,model,hypothesis)
    attempt = 0 
    if(FPE<FP){
      return(lo)
    }else{
      while(FP<pro){
        attempt = attempt +1
        lo = lo +1
        x = bin_BF_bound_10 (D,lo,alpha,beta,location,scale,model,hypothesis)
        pro = bin_FPE(x,lo,alpha,beta,location,scale,model,hypothesis)
        
        if(attempt>100){
          return(lo)
          break
        }
      }

       return(lo)

      }
      
      
    
    
    
    
    
  }
  
  Power_root <- function(N){
  N =round(N)
    x = bin_BF_bound_10 (D,N,alpha,beta,location,scale,model,hypothesis)

    
    if(de_an_prior == 1){
      pro = bin_TPE(x,N,alpha,beta,location,scale,model,hypothesis)
    } else{
      pro = bin_TPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
    }
    #pro = switch(de_an_prior,
    #             "1" = bin_TPE(x,N,alpha,beta,location,scale,model,hypothesis),
    #            "0" = bin_TPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
    #            )
    
    return(pro-target)
  }
  
  up= 50000
  N = round(uniroot(Power_root,lower = lo,upper = up)$root)
  x10 = bin_BF_bound_10(D,N,alpha,beta,location,scale,model,hypothesis)
  FPE = bin_FPE(x10,N,alpha,beta,location,scale,model,hypothesis)
  
  if(FPE<FP){
    return(N)
  }else{
    alpha_root <- function(N){
      N =round(N)
      x = bin_BF_bound_10(D,N,alpha,beta,location,scale,model,hypothesis)
      pro = bin_FPE(x,N,alpha,beta,location,scale,model,hypothesis)

      #pro = switch(de_an_prior,
      #             "1" = bin_TPE(x,N,alpha,beta,location,scale,model,hypothesis),
      #            "0" = bin_TPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
      #            )
      
      return(pro-FP)
    }
    
    N = round(uniroot(alpha_root,lower = N,upper = up)$root)
    
    
    return(N)
    
    
    
    
    
    
    #attempt = 0
    #while(FP<pro){
     # attempt = attempt +1
     # N = N +5
     # x = bin_BF_bound_10 (D,N,alpha,beta,location,scale,model,hypothesis)
      #pro = switch(de_an_prior,
      #             "1" = bin_FPE(x,N,alpha,beta,location,scale,model,hypothesis),
      #             "0" = bin_FPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
      #)
      #if(attempt>100){
      #  return(N)
      #  break
      #}
    }
    
    #return(N)
    
  #}
    
    
  
}

bin_table<-function(D,target,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,N, mode_bf,FP){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  if (mode_bf == "0"){
    n = N
  }else {
    
    n =  bin_N_finder(D,target,alpha,beta,location,scale,model,hypothesis,
                                 alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP)
  
  }
  bound01 = bin_BF_bound_01(D,n,alpha,beta,location,scale,model,hypothesis)
  
  if (length(bound01)==0|any(bound01 == "bound cannot be found")){
    FNE = 0
    TNE = 0
  }else{
    
    if (de_an_prior ==1){
      FNE = bin_FNE(bound01,n,alpha,beta,location,scale,model,hypothesis)
      TNE = bin_TNE(bound01,n,alpha,beta,location,scale,model,hypothesis)
    }
    if (de_an_prior ==0){
      FNE = bin_FNE(bound01,n,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
      TNE = bin_TNE(bound01,n,alpha,beta,location,scale,model,hypothesis)
    }
}
  bound10 = bin_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis)
  

  if (de_an_prior ==1){
    TPE = bin_TPE(bound10,n,alpha,beta,location,scale,model,hypothesis)
    FPE = bin_FPE(bound10,n,alpha,beta,location,scale,model,hypothesis)
  }else{
    TPE = bin_TPE(bound10,n,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
    FPE = bin_TNE(bound10,n,alpha,beta,location,scale,model,hypothesis)
  }
  
  
  table <- data.frame(
    TPE = TPE,
    FNE = FNE,
    TNE = TNE,
    FPE = FPE,
    N =  n)
  return(table)
}

bin_bf10 <-function(DD,n,alpha,beta,location,scale,model,hypothesis){
  x= seq(from = 0,to =n,by= 5)
  BF_D = bin_BF_bound_10(DD,n,alpha,beta,location,scale,model,hypothesis)
  D = round(bin_BF(BF_D,n,alpha,beta,location,scale,model,hypothesis),2)
  BF10 = bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)
  
  
  if (length(D)== 2){part1 = bquote(bold("BF"[10] ~ "=" ~ .(D[1]) / .(D[2])))}else{part1 = bquote(bold("BF"[10] ~ "=" ~ .(D[1])))}
  
  if (length(D)== 2){  part2 = bquote("when x = " ~ .(BF_D[1]) / .(BF_D[2]))}else{  part2 <- bquote("when x = " ~ .(BF_D[1]))}
  
  

  
  main = bquote(.(part1) ~ .(part2))
 
  
  par(mfrow = c(1, 2))
  plot(x,log10(BF10),xlab= "Number of success",type="l", ylab = expression("logarithm of BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
  abline(v = BF_D)
  axis(1, c(0,n))
  if (length(BF_D) != 0 ){
    for(i in 1:length(BF_D))
      axis(1, BF_D[i])
      }
  
  max_BF = 1/bin_BF(ceiling(location*n),n,alpha,beta,location,scale,model,hypothesis)
  
  BF_D = bin_BF_bound_01(DD,n,alpha,beta,location,scale,model,hypothesis)

  
  plot(x,log10(1/BF10),xlab= "Number of success",type="l",main = "",frame.plot = FALSE,ylab = bquote("logarithm of BF"[0][1]),xaxt = "n")
  axis(1, c(0,n))
  
  if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found" |BF_D == "no bound is found" ) ) {
    main = bquote(bold("It is impossible to have BF"[0][1]~"="~.(D)))
    title(main = main)
    #sprintf("It is impossible to have BF01 = %.3f ",D)
  } else      {
    D = round(1/bin_BF(BF_D,n,alpha,beta,location,scale,model,hypothesis),2)
    
    abline(v = BF_D)
    
    if (length(BF_D) != 0 ){
      for(i in 1:length(BF_D))
        axis(1, BF_D[i])
    }
    
    if (length(D)== 2){part1 = bquote(bold("BF"[10] ~ "=" ~ .(D[1]) / .(D[2])))}else{part1 = bquote(bold("BF"[10] ~ "=" ~ .(D[1])))}
    
    if (length(D)== 2){  part2 = bquote("when x = " ~ .(BF_D[1]) / .(BF_D[2]))}else{  part2 <- bquote("when x = " ~ .(BF_D[1]))}
    
    main = bquote(.(part1) ~ .(part2))
    title(main = main)
    }
}


Power_bin<-function(D,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d, de_an_prior,N){

  smin = 10
  smax = N*1.2
  sdf = ceiling(seq(smin,smax , by = (smax-smin)/50))
  
  power =  array(NA, dim = c(length(sdf)))
  
  for ( i in 1:length(sdf)){
    x = bin_BF_bound_10(D,sdf[i],alpha,beta,location,scale,model,hypothesis)
 
    if(de_an_prior ==1){
      power[i] = bin_TPE(x,sdf[i],alpha,beta,location,scale,model,hypothesis)
    }else{
      power[i] =bin_TPE(x,sdf[i],alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
    }
    #power[i] = switch(de_an_prior,
    #                  "1" =  bin_TPE(x,sdf[i],alpha,beta,location,scale,model,hypothesis), 
    #                  "0" = bin_TPE(x,sdf[i],alpha_d,beta_d,location_d,scale_d,model_d,hypothesis))
    #
    
  }
  plot(sdf,power,type="l",main = "",frame.plot = FALSE,xlab = "Total sample size", ylab = "Probability of True positive evidence",xlim = c(10,max(sdf)), 
       ylim = c(0,1) )
  
}



bin_prior_plot <-function(alpha,beta,location,scale,model,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,de_an_prior){
  par(mfrow = c(1, 1))
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = 0, b = location),
                   "!=" = c(a = 0, b = 1)
  )
  propp = seq(bound[1],bound[2],.01)
  
  normalization  <- integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound[1],upper = bound[2],rel.tol = 1e-10)$value

  
  
  prior_prop = NA
  prior_prop = bin_prior(propp,alpha,beta,location,scale,model)/normalization
  
  
  
  plot(propp,prior_prop,xlab= bquote(bold(p)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~p~" under the alternative hypothesis")),frame.plot = FALSE)
  
  
  if (de_an_prior ==0){
    
    if (model_d != "Point"){
      
      normalization_d  <- integrate(function(prop)bin_prior(prop,alpha_d,beta_d,location_d,scale_d,model_d),lower = bound[1],upper = bound[2],rel.tol = 1e-10)$value
      
      prior_pro_D = bin_prior(propp,alpha_d,beta_d,location_d,scale_d,model_d)/normalization_d
      plot(propp,prior_prop,ylim =c(0,max(prior_pro_D,prior_prop)),xlab= bquote(bold(rho)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~rho~" under the alternative hypothesis")),frame.plot = FALSE)
      
      lines(propp,prior_pro_D ,lty = 2)
      
    } else {
      plot(propp,prior_prop,
           xlab= bquote(bold(p)),
           ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~p~" under the alternative hypothesis")),frame.plot = FALSE)
      
      arrows(x0=location_d, y0=0, x1=location_d, y1=max(prior_prop), length=0.2, code=2, col="black", lwd=1,,lty = 2)
      
      
    }
    
    
    legend("topright", 
           legend = c("Analysis prior", "Design prior"), 
           lty = c(1, 2), 
           col = c("black", "black"),
           bty = "n") 
    
  }
  
  
}




