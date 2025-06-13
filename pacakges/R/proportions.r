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
    return(pro - target - 1e-3)
  }
  n1 <- suppressWarnings(round(uniroot(power, lower = lo_n1, upper = 5000,maxiter = 10)$root))
  grid_power <- grid
  
 
  
  return(list(grid_power,n1))
  
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
  list(table,grid)
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
       xlab = substitute(bold(p[nu_val]), list(nu_val = nu)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~p[.(nu)]~"")),
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



Power_p2<-function(D,n1, a0, b0, a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2){
  smax = n1*1.2
  Ns = ceiling(seq(10,smax,by = (smax-10)/20))
  n2 = Ns*r
  Nt= Ns+n2
  TPE =  array(NA, dim = c(length(Ns)))
  FPE =  array(NA, dim = c(length(Ns)))
  TNE =  array(NA, dim = c(length(Ns)))
  FNE =  array(NA, dim = c(length(Ns)))
    
    for(i in 1:length(Ns)){
      
      n2 = Ns[i]*r
      grid <- BF_grid_rcpp(D, a0, b0, a1, b1, Ns[i], a2, b2, n2,model1,da1,db1,dp1,model2,da2,db2,dp2) 
      TPE[i] <- sum_rcpp(grid$log_h1_dp,grid$PE)
      FPE[i] <- sum_rcpp(grid$log_h0,grid$PE)
      TNE[i] <- sum_rcpp(grid$log_h0,grid$NE)
      FNE[i] <- sum_rcpp(grid$log_h1_dp,grid$NE)
  
    }

  par(mfrow = c(1, 2))
  plot(Nt, TPE, type = "l",
       xlab = "Total sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(Nt,FPE,col = "grey")
  legend(x = max(Nt)*.4,y=.5,              # position of the legend
         legend = c("True positive", "False positive"),  # labels
         col = c("black", "grey"),   # colors matching lines
         lty = 1,                   # line type (solid)
         bty = "n") 
  plot(Nt, TNE, type = "l",
       xlab = "Total ample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[0][1]~">"~.(D))))
  lines(Nt,FNE,col = "grey")
  legend(x = max(Nt)*.4,y=.5,              # position of the legend
         legend = c("True negative", "False negative"),  # labels
         col = c("black", "grey"),   # colors matching lines
         lty = 1,                   # line type (solid)
         bty = "n") 
  
  
  
  
}



heatmap_p2_0<-function(x,D){
  # Prepare data
  df <- data.frame(
    k1 = x$k1,
    k2 = x$k2,
    PE = x$PE,
    NE = x$NE,
    BF = x$log_BF10
  )
  
  # Derive effect type
  df$effect <- with(df, ifelse(PE == 1, "PE",
                               ifelse(NE == 1, "NE", "None")))
  
  # Use bquote to fix the math expressions in legend labels
  labels <- c(
    "PE"   = bquote(BF[10] > .(D)),
    "NE"   = bquote(BF[0][1] > .(D)),
    "None" = bquote(1 / .(D) < BF[10] ~ "<" ~ .(D))
  )
  
  # Plot
  ggplot(df, aes(x = k1, y = k2, fill = effect)) +
    geom_tile() +
    scale_fill_manual(name ="classification",
                      values = c("PE" = "#440154", "NE" = "#21908C", "None" = "#FDE725"),
      
      labels = labels
    ) +
    labs(title = "BF and number of success",
         x = "x1", y = "x2") +
    coord_fixed() +
    theme_minimal()
  
  
  
  
  
  ggplot(df, aes(x = k1, y = k2, fill = BF)) +
    geom_tile() +
    scale_fill_viridis_c(name = "BF") +   # continuous Viridis color scale
    labs(
      title = "Heatmap of BF",
      x = "k1",
      y = "k2"
    ) +
    coord_fixed() +  # equal aspect ratio
    theme_minimal()
  
  
  
  
}



heatmap_p2<-function(x,D){
  # Prepare data
  df <- data.frame(
    k1 = x$k1,
    k2 = x$k2,
    PE = x$PE,
    NE = x$NE,
    BF = x$log_BF10
  )
  
  # Derive effect type
  df$effect <- with(df, ifelse(PE == 1, "PE",
                               ifelse(NE == 1, "NE", "None")))
  
  # Use bquote to fix the math expressions in legend labels
  labels <- c(
    "PE"   = bquote(BF[10] > .(D)),
    "NE"   = bquote(BF[0][1] > .(D)),
    "None" = bquote(1 / .(D) < BF[10] ~ "<" ~ .(D))
  )
  
  # First plot: categorical heatmap with custom colors and math labels
  p1 <- ggplot(df, aes(x = k1, y = k2, fill = effect)) +
    geom_tile() +
    scale_fill_manual(
      name = "Classification",
      values = c("PE" = "#FDE725", "NE" = "#440154", "None" = "#21908C"),
      labels = labels
    ) +
    labs(title = "BF and number of success", x = "x1", y = "x2") +
    coord_fixed() +
    theme_minimal()
  
  # Second plot: continuous heatmap of BF
  p2 <- ggplot(df, aes(x = k1, y = k2, fill = BF)) +
    geom_tile() +
    scale_fill_viridis_c(name = "log BF10") +
    labs(title = "Heatmap of BF", x = "k1", y = "k2") +
    coord_fixed() +
    theme_minimal()
  
  # Combine side by side
  combined_plot <- p1 + p2 + plot_layout(ncol = 2)
  
  # Display
  print(combined_plot)
  
  
  
}

