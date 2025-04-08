
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
  
  normalizationh1  <- switch(hypothesis,
                             "!=" = {integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = 0,upper = bound_h1[1],rel.tol = 1e-10)$value+integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[2],upper = 1,rel.tol = 1e-10)$value},
                             ">"  = integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                             "<"  = integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
  
  normalizationh0 <- integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h0[1],upper = bound_h0[2])$value
  for (i in 1:length(x)){
  int  <- function(prop){dbinom(x[i], size=n, prob=prop) *bin_prior(prop,alpha,beta,location,scale,model)}
  
  if (hypothesis == "!="){
    lh1 = integrate(int,lower = 0,upper = bound_h1[1], rel.tol=1e-10,stop.on.error = F)$value+integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-10,stop.on.error = F)$value 
  }else{
    lh1 = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-10,stop.on.error = F)$value
    
  }
  lh0 = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-10,stop.on.error = F)$value

    
    BF[i] = lh1/lh0
  }
  return(BF)
  
}



bin_e_BF_bound_10 <-function(D,n,alpha,beta,location,scale,model,hypothesis,e){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)- D
  }
  switch(hypothesis,
         "!="= {

           
           
      x <-  tryCatch( round(uniroot(Bound_finding,lower = 0 ,upper = ceiling(n/2))$root), error=function(e){})
      y <-tryCatch( round(uniroot(Bound_finding,lower = ceiling(n/2) ,upper =n)$root), error=function(e){})

      },
      ">" = {
        x <- tryCatch( round(uniroot(Bound_finding,lower = ceiling(n/2) ,upper =n)$root), error=function(e){})
      },
      "<" = {
        x <- tryCatch( round(uniroot(Bound_finding,lower = 0 ,upper = ceiling(n/2))$root), error=function(e){})

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
  
  
  BF = bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)
  

  
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
    BF = bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)
  }  
    
    
    
    return(x)
}


bin_e_BF_bound_01 <-function(D,n,alpha,beta,location,scale,model,hypothesis,e){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    1/bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)- D
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
  
  BF = 1/bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)
  
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
    BF = 1/bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)
  }  
 
  
  return(x)
}

bin_e_TPE<-function(x,n,alpha,beta,location,scale,model,hypothesis,e){
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
  normalizationh1  <- switch(hypothesis,
                             "!=" = {integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = 0,upper = bound_h1[1],rel.tol = 1e-10)$value+integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[2],upper = 1,rel.tol = 1e-10)$value},
                             ">"  = integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                             "<"  = integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
  
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
           
           pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh1
           },
         ">" = 
           int  <- function(prop){
             pro = pbinom(x-1,n,prop,lower.tail = F)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh1
           },
         "<" = 
           int  <- function(prop){
             pro = pbinom(x,n,prop,lower.tail = T)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh1
           })
  
  
  
  if(hypothesis == "!="){
    TPE = integrate(int,lower = 0,upper = bound_h1[1], rel.tol = 1e-10)$value + integrate(int,lower = bound_h1[2],upper = 1, rel.tol = 1e-10)$value
  }else{
    TPE = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol = 1e-10)$value
  }
  
  
  return(TPE)
  
}


bin_e_FNE<-function(x,n,alpha,beta,location,scale,model,hypothesis,e){
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

  
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = location+e, b = 1),
                      "<" = c(a = 0, b = location+e),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  normalizationh1  <- switch(hypothesis,
                             "!=" = {integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = 0,upper = bound_h1[1],rel.tol = 1e-10)$value+integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[2],upper = 1,rel.tol = 1e-10)$value},
                             ">"  = integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                             "<"  = integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
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
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh1
           },
         ">" = 
           int  <- function(prop){
             pro = pbinom(x,n,prop,lower.tail = T)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh1
           },
         "<" = 
           int  <- function(prop){
             pro = pbinom(x-1,n,prop,lower.tail = F)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh1
           })
  
  if(hypothesis == "!="){
    FNE = integrate(int,lower = 0,upper = bound_h1[1], rel.tol = 1e-10)$value + integrate(int,lower = bound_h1[2],upper = 1, rel.tol = 1e-10)$value
  }else{
    FNE = integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol = 1e-10)$value
  }
  
  
  return(FNE)
  
}


bin_e_FPE<-function(x,n,alpha,beta,location,scale,model,hypothesis,e){
  
  if (any(x =="no bound is found" | length(x)==0)){
    x=0
    return(x)
  }
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = location, b = location+e),
                      "<" = c(a = location+e, b = location),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  normalizationh0 <- integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h0[1],upper = bound_h0[2])$value
  
  
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
             
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh0
           },
         ">" = 
           int  <- function(prop){
             pro = pbinom(x-1,n,prop,lower.tail = F)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh0
           },
         "<" = 
           int  <- function(prop){
             pro = pbinom(x,n,prop,lower.tail = T)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh0
           })
  
  
  

    FPE = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol = 1e-10)$value
    return(FPE)

}


bin_e_TNE<-function(x,n,alpha,beta,location,scale,model,hypothesis,e){
  
  if (any(x =="no bound is found" | length(x)==0)){
    x=0
    return(x)
  }
  
  
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = location, b = location+e),
                      "<" = c(a = location+e, b = location),
                      "!=" = c(a = location+e[1], b = location+e[2])
  )
  
  normalizationh0 <- integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h0[1],upper = bound_h0[2])$value
  
  
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
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh0
           },
         ">" = 
           int  <- function(prop){
             pro = pbinom(x,n,prop,lower.tail = T)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh0
           },
         "<" = 
           int  <- function(prop){
             pro = pbinom(x-1,n,prop,lower.tail = F)
             pro*bin_prior(prop,alpha,beta,location,scale,model)/normalizationh0
           })
  

    TNE = integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol = 1e-10)$value
  
  
  
  return(TNE)
}

bin_e_N_finder <-function(D,target,alpha,beta,location,scale,model,hypothesis,
                        alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP,e){
  lo = 10
  x = bin_e_BF_bound_10(D,lo,alpha,beta,location,scale,model,hypothesis,e)
  pro = bin_e_TPE(x,lo,alpha,beta,location,scale,model,hypothesis,e)
  
  if ( pro>target){

    FPE = bin_e_FPE(x,lo,alpha,beta,location,scale,model,hypothesis,e)
    #attempt = 0 
    if(FPE<FP){
      return(lo)
    }else{
      while(FP<pro){
        alpha_root <- function(N){
          N =round(N)
          x = bin_e_BF_bound_10(D,N,alpha,beta,location,scale,model,hypothesis,e)
          
          
          if(de_an_prior == 1){
            pro = bin_e_FPE(x,N,alpha,beta,location,scale,model,hypothesis,e)
          } else{
            pro = bin_e_FPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
          }
          return(pro-FP)
        }
        
        N = round(uniroot(alpha_root,lower = lo,upper = up)$root)
        
        
        return(N)}}}
        
        
        
        
        
        
        
        #attempt = attempt +1
        #lo = lo +1
        #x = bin_e_BF_bound_10(D,lo,alpha,beta,location,scale,model,hypothesis,e)
        #pro = switch(de_an_prior,
        #             "1" = bin_e_FPE(x,lo,alpha,beta,location,scale,model,hypothesis,e),
        #             "0" = bin_e_FPE(x,lo,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
        #)
        #if(attempt>100){
        #  return(lo)
        #  break
        #}
      

       #return(lo)
 

  
  Power_root <- function(N){
    N =round(N)
    x = bin_e_BF_bound_10(D,N,alpha,beta,location,scale,model,hypothesis,e)

    
    if(de_an_prior == 1){
      pro = bin_e_TPE(x,N,alpha,beta,location,scale,model,hypothesis,e)
    } else{
      pro = bin_e_TPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
    }
    #pro = switch(de_an_prior,
    #             "1" = bin_TPE(x,N,alpha,beta,location,scale,model,hypothesis),
    #            "0" = bin_TPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
    #            )
    
    return(pro-target)
  }
  
  up= 50000
  N = round(uniroot(Power_root,lower = lo,upper = up)$root)
  x10 = bin_e_BF_bound_10(D,N,alpha,beta,location,scale,model,hypothesis,e)
  FPE = bin_e_FPE(x10,N,alpha,beta,location,scale,model,hypothesis,e)
  
  if(FPE<FP){
    return(N)
  }else{
    
    
    alpha_root <- function(N){
      N =round(N)
      x = bin_e_BF_bound_10(D,N,alpha,beta,location,scale,model,hypothesis,e)
      
      
      if(de_an_prior == 1){
        pro = bin_e_FPE(x,N,alpha,beta,location,scale,model,hypothesis,e)
      } else{
        pro = bin_e_FPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
      }
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
    #  attempt = attempt +1
    #  N = N +5
    #  x = bin_e_BF_bound_10(D,N,alpha,beta,location,scale,model,hypothesis,e)
    #  pro = switch(de_an_prior,
    #               "1" = bin_e_FPE(x,N,alpha,beta,location,scale,model,hypothesis,e),
    #               "0" = bin_e_FPE(x,N,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
    #  )
    #  if(attempt>100){
    #    
    #    return(N)
    #    break
    # }
    #}
    
    
    
  }
    
    
  
}

bin_e_table<-function(D,target,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,N, mode_bf,FP,e){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)
  if (mode_bf == "0"){
    n = N
  }else {
    
    n =  bin_e_N_finder(D,target,alpha,beta,location,scale,model,hypothesis,
                                 alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP,e)
  
  }
  bound01 = bin_e_BF_bound_01(D,n,alpha,beta,location,scale,model,hypothesis,e)
  
  if (length(bound01)==0|any(bound01 == "bound cannot be found")){
    FNE = 0
    TNE = 0
  }else{
    
    if (de_an_prior ==1){
      FNE = bin_e_FNE(bound01,n,alpha,beta,location,scale,model,hypothesis,e)
      TNE = bin_e_TNE(bound01,n,alpha,beta,location,scale,model,hypothesis,e)
    }
    if (de_an_prior ==0){
      FNE = bin_e_FNE(bound01,n,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
      TNE = bin_e_TNE(bound01,n,alpha,beta,location,scale,model,hypothesis,e)
    }
}
  bound10 = bin_e_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis,e)
  

  if (de_an_prior ==1){
    TPE = bin_e_TPE(bound10,n,alpha,beta,location,scale,model,hypothesis,e)
    FPE = bin_e_FPE(bound10,n,alpha,beta,location,scale,model,hypothesis,e)
  }else{
    TPE = bin_e_TPE(bound10,n,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
    FPE = bin_e_FPE(bound10,n,alpha,beta,location,scale,model,hypothesis,e)
  }
  
  
  table <- data.frame(
    TPE = TPE,
    FNE = FNE,
    TNE = TNE,
    FPE = FPE,
    N =  n)
  return(table)
}

bin_e_bf10 <-function(DD,n,alpha,beta,location,scale,model,hypothesis,e){
  x= seq(from = 0,to =n,by= 5)
  BF_D = bin_e_BF_bound_10(DD,n,alpha,beta,location,scale,model,hypothesis,e)
  D = round(bin_e_BF(BF_D,n,alpha,beta,location,scale,model,hypothesis,e),2)
  BF10 = bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)
  
  
  
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
  
  max_BF = 1/bin_e_BF(ceiling(location*n),n,alpha,beta,location,scale,model,hypothesis,e)
    
  
  BF_D = bin_e_BF_bound_01(DD,n,alpha,beta,location,scale,model,hypothesis,e)

  
  plot(x,log10(1/BF10),xlab= "Number of success",type="l",main = "",frame.plot = FALSE,ylab = bquote("logarithm of BF"[0][1]),xaxt = "n")
  axis(1, c(0,n))
  
  if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found" |BF_D == "no bound is found" ) ) {
    main = bquote(bold("It is impossible to have BF"[0][1]~"="~.(D)))
    title(main = main)
    #sprintf("It is impossible to have BF01 = %.3f ",D)
  } else      {
    D = round(1/bin_e_BF(BF_D,n,alpha,beta,location,scale,model,hypothesis,e),2)
    
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


Power_e_bin<-function(D,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d, de_an_prior,N,e){

  smin = 10
  smax = N*1.2
  sdf = ceiling(seq(smin,smax , by = (smax-smin)/50))
  
  power =  array(NA, dim = c(length(sdf)))
  
  for ( i in 1:length(sdf)){
    x = bin_e_BF_bound_10 (D,sdf[i],alpha,beta,location,scale,model,hypothesis,e)

 
    if(de_an_prior ==1){
      power[i] = bin_e_TPE(x,sdf[i],alpha,beta,location,scale,model,hypothesis,e)
        
        
        
    }else{
      power[i] = bin_e_TPE(x,sdf[i],alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
    }
    #power[i] = switch(de_an_prior,
    #                  "1" =  bin_TPE(x,sdf[i],alpha,beta,location,scale,model,hypothesis), 
    #                  "0" = bin_TPE(x,sdf[i],alpha_d,beta_d,location_d,scale_d,model_d,hypothesis))
    #
    
  }
  plot(sdf,power,type="l",main = "",frame.plot = FALSE,xlab = "Total sample size", ylab = "Probability of True positive evidence",xlim = c(10,max(sdf)), 
       ylim = c(0,1) )
  
}



bin_e_prior_plot <-function(alpha,beta,location,scale,model,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,de_an_prior,e){
  par(mfrow = c(1, 1))

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
  
  normalizationh1  <- switch(hypothesis,
                             "!=" = {integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = 0,upper = bound_h1[1],rel.tol = 1e-10)$value+integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[2],upper = 1,rel.tol = 1e-10)$value},
                             ">"  = integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                             "<"  = integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
  normalizationh0 <- integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound_h0[1],upper = bound_h0[2])$value
  
  
  
  
  
  propp = switch(hypothesis,
               "!=" = seq(0,1,.001),
               ">" = seq(location,1,.001),
               "<" = seq(0,location,.001)
  )
  
  


  

  prior_h1 = bin_prior(propp,alpha,beta,location,scale,model)/normalizationh1
  prior_h0 = bin_prior(propp,alpha,beta,location,scale,model)/normalizationh0
  
  
  switch(hypothesis,
         "!=" = { prior_h1[propp>=min(bound_h1)&propp<=max(bound_h1)]=-1 },
         ">" = { prior_h1[propp<min(bound_h1)]=-1 },
         "<" = { prior_h1[propp>max(bound_h1)]=-1 }
         
  )
  
  switch(hypothesis,
         "!=" = { prior_h0[!(propp>min(bound_h0)&propp<max(bound_h0))]=-1},
         ">" = { prior_h0[propp>max(bound_h0)]=-1 },
         "<" = { prior_h0[propp<min(bound_h0)]=-1 }
         
  )
  
  
  plot(propp,prior_h0,
       type = "l",
       lty = 3,
       lwd = 4,
       xlab = bquote(bold(p)),
       ylim = c(0,max(prior_h0)*1.1),
       ylab = "density",
       main = bquote(bold("prior distribution on " ~ rho ~ " under the two hypotheses")),
       frame.plot = FALSE)
  
  lines(propp,prior_h1, type="l", lty=1, lwd=4)
  legend("topright", 
         legend = c("H1", "H0"),  # Labels for the lines
         col = c("black", "black"),  # Line colors
         lty = c(1, 3),  # Line types
         lwd = c(4, 4),  # Line widths
         bty = "n",
         title = "Analysis Prior")  # No box around the legend
  
  
  
  if (de_an_prior ==0){
    if (model_d == "Point"){
      
      arrows(x0=location, y0=0, x1=location, y1=max(prior_h1,prior_h0), length=0.2, code=2, col="gray", lwd=4,lty = 2)
    }else{
      
      normalizationd  <- switch(hypothesis,
                                 "!=" = 1-integrate(function(prop)bin_prior(prop,alpha_d,beta_d,location_d,scale_d,model_d),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                                 ">"  = integrate(function(prop)bin_prior(prop,alpha_d,beta_d,location_d,scale_d,model_d),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value,
                                 "<"  = integrate(function(prop)bin_prior(prop,alpha_d,beta_d,location_d,scale_d,model_d),lower = bound_h1[1],upper = bound_h1[2],rel.tol = 1e-10)$value)
      
      prior_d = bin_prior(propp,alpha_d,beta_d,location_d,scale_d,model_d)/normalizationd
      
      
      switch(hypothesis,
             "!=" = { prior_d[propp>=min(bound_h1)&propp<=max(bound_h1)]=-1 },
             ">" = { prior_d[propp<min(bound_h1)]=-1 },
             "<" = { prior_d[propp>max(bound_h1)]=-1 }
             
      )
      
      lines(propp,prior_d, type="l", lty=1, lwd=4,col="gray")
      
      
    }
    legend("topleft", 
           legend = c("Analysis prior", "Design prior"), 
           lty = c(1, 1), 
           lwd = c(4, 4),
           col = c("black", "gray"),
           bty = "n") 
    
  }
  
  
}




