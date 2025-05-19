BF10_p2<-function(a0, b0, a1, b1, a2, b2,n1,n2,k1,k2){
  
  logBF = lbeta(k1 + k2 + a0, n1 + n2 - k1 - k2 + b0) -
    lbeta(k1 + a1, n1 - k1 + b1) -
    lbeta(k2 + a2, n2 - k2 + b2) +
    lbeta(a1, b1) +
    lbeta(a2, b2) -
    lbeta(a0, b0)
  1/exp(logBF)
}

ps_N_finder<-function(D,target, a0, b0, a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2) {
  
  lo_n1 <- 10
  n2 <- round(lo_n1)*r
  grid <- BF_grid_rcpp(D, a0, b0, a1, b1, lo_n1, a2, b2, n2,model1,da1,db1,dp1,model2,da2,db2,dp2) 
  pro <- sum_rcpp(grid$log_h1_dp,grid$PE)
  
  if ( pro>target){
    return(list(grid,lo_n1))
  }
  power<-function(n1){
    n1 = round(n1)
    n2 = n1*r
    grid <<- BF_grid_rcpp(D, a0, b0, a1, b1, n1, a2, b2, n2,model1,da1,db1,dp1,model2,da2,db2,dp2) 
    pro <- sum_rcpp(grid$log_h1_dp,grid$PE)
    return(pro-target)
  }
  n1 <- suppressWarnings(round(uniroot(power, lower = lo_n1, upper = 5000,maxiter = 10)$root))
  
  return(list(grid,n1))
  
}

pro_table_p2<-function(D,target, a0, b0, a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2,mode_bf,n1,n2) {
  
  if (mode_bf==1){
    x = ps_N_finder(D,target, a0, b0, a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2)
    grid = x[[1]]
    n1  = x[[2]]
    n2 =  x[[2]]*r
  }else{
    grid = BF_grid_rcpp(D, a0, b0, a1, b1, n1, a2, b2, n2,model1,da1,db1,dp1,model2,da2,db2,dp2) 
  }
  table <- data.frame(
    TPE = sum_rcpp(grid$log_h1_dp,grid$PE),
    FNE = sum_rcpp(grid$log_h1_dp,grid$NE),
    TNE = sum_rcpp(grid$log_h0,grid$NE),
    FPE = sum_rcpp(grid$log_h0,grid$PE),
    n1  = n1,
    n2 =  n2)
  colnames(table) <- c(sprintf("p(BF10> %0.f|H1)",D), sprintf("p(BF01> %0.f|H1)",D), sprintf("p(BF01> %0.f|H0)",D), sprintf("p(BF10> %0.f|H0)",D), " N1", "N2")
  table
}


p2_prior_plot<-function(a,b,ad,bd,dp,model,nu){

  par(mfrow = c(1, 1))
  
  prop    <- seq( 0,1,.001)
  
  prior.analysis <- dbeta(prop,a,b)
  prior.design   <- switch(model,
                           "same" = dbeta(prop,a,b),
                           "beta" = dbeta(prop,ad,bd),
                           "Point" = rep(NA, length(prop)))
  
  ylim.max <- max(prior.analysis, prior.design, na.rm = TRUE)
  
  plot(prop, prior.analysis, type = "l", lwd = 2,
       xlab = substitute(bold(rho[nu_val]), list(nu_val = nu)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~rho[.(nu)]~"")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))
  
  if (model != "same") {
    if (model == "Point")
      arrows(x0 = dp, y0 = 0, x1 = dp, y1 = ylim.max, length = 0.1, col = "black", lty = 2) else
        lines(prop, prior.design, lty = 2)
    
    # Add legend:
    legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }
  
}

