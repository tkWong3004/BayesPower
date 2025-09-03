library(BayesPower)
library(mombf)
library(ExtDist)
library(mvtnorm)
# functions to simulate delta
sim_pro_h0 <-function(iter, a0,b0){
  pro1 <- rbeta(iter, a0, b0)
  pro2 <- pro1

  # Make a 2-column matrix
  cbind(pro1, pro2)

}

sim_pro_h1 <-function(iter, a1,b1,p1,a2,b2,p2,model1,model2){

  pro1 <- switch(model1,
                 "beta"={rbeta(iter, a1, b1)},
                 "Point"= rep(p1,iter))
  pro2 <- switch(model2,
                 "beta"={rbeta(iter, a2, b2)},
                 "Point"= rep(p2,iter))
  cbind(pro1, pro2)

}


sim_success_ps<-function(iter, D,a0, b0, n1,n2,a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2,h0=F){
  if(h0) pro=sim_pro_h0(iter, a0,b0) else{

  if (model1==model2&model1=="same"){
    pro=sim_pro_h1(iter, a1,b1,p1,a2,b2,p2,model1,model2)
  } else {
    pro=sim_pro_h1(iter, da1,db1,dp1,da2,db2,dp2,model1,model2)
  }
  }


  x1 <- sapply(1:iter, function(i) {
  rbinom(1,n1,pro[i,1])
})
  x2<-sapply(1:iter, function(i) {
    rbinom(1,n2,pro[i,2])
  })
  BF10 = BayesPower:::BF10_p2(a0, b0, a1, b1, a2, b2,n1,n2,x1,x2)

  # --- Pro(BF>k) ---
  PE_sim <- sum(BF10>=D)/iter
  NE_sim <- sum((1/BF10)>=D)/iter

  return(c(PE_sim, NE_sim))

}


# comparing the results with the numeric method and simulation
sim_method_pros<-function(iter,D ,a0, b0, a1, b1, a2, b2, model1,da1,db1,dp1,model2,da2,db2,dp2,n){
  n1=n2=n
  results <-
    BayesPower:::pro_table_p2(D,target=.8, a0, b0, a1, b1, a2, b2, r=1,
                            model1,da1,db1,dp1,model2,da2,db2,dp2,mode_bf=0,n1,n2,direct="h1")[[1]]


  pro_h1<- sim_success_ps(iter, D,a0, b0, n1,n2,a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2)
  pro_h0<- sim_success_ps(iter, D,a0, b0, n1,n2,a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2,h0=T)

  pro_h0<-c(pro_h0[2],pro_h0[1])
  results = as.vector(unlist(results))
  out_matrix = c(results[1:4], pro_h1, pro_h0)
  return(out_matrix)
}

# making the plot and doing it for seqence of N
p2_ver<-function(iter,D, a0, b0, a1, b1, a2, b2, model1,da1,db1,dp1,model2,da2,db2,dp2){
  if(model1==model2&model1=="same"){

    model1=model2="beta"
    da1=a1
    db1=b1
    da2=a2
    db2=b2
  }

  n = round(seq(50,1000,by=(1000-50)/10))
  total = n

  probs <- matrix(NA, nrow = length(n), ncol = 8)  # simpler than array
  for(i in 1:length(n)){

    probs[i,] =  as.vector(unlist(sim_method_pros(iter,D, a0, b0, a1, b1, a2, b2, model1,da1,db1,dp1,model2,da2,db2,dp2,n[i])))

  }


  par(mfrow = c(1, 2))
  plot(total, probs[,1], type = "l",
       xlab = "Total sample size",
       ylab = "Probability",
       ylim = c(0,1),
       frame.plot = FALSE,
       lwd = 3,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
  lines(total, probs[,5], col = "gray", lwd = 2, lty = 3)
  lines(total, probs[,2], col = "black", lwd = 3)
  lines(total, probs[,6], col = "grey", lwd = 2, lty = 3)


  plot(total, probs[,3], type = "l",
       xlab = "Total sample size",
       ylab = "Probability",
       ylim = c(0,1),
       frame.plot = FALSE,
       lwd = 3,
       main = bquote(bold("Power curve for BF"[0][1]~">"~.(D))))
  lines(total, probs[,3+4], col = "gray", lwd = 2, lty = 3)
  lines(total, probs[,4], col = "black", lwd = 3)
  lines(total, probs[,4+4], col = "grey", lwd = 2, lty = 3)

}
# input
iter = 5000 # number of iteration
# beta prior under h0
a0=b0=5
### beta prior under h1
#p1
a1=44
b1=44
#p2
a2=99
b2=99

# if analysis prior and design prior are the same under h1
model1=model2="same"
# else specify the things below


###design prior under h1

# p1
model1="beta"    #"Point" or "beta"
da1=db1=2
dp1=.2
#p2
model2="beta"    #"Point" or "beta"
da2=db2=2
dp2=.2


p2_ver(iter,D, a0, b0, a1, b1, a2, b2, model1,da1,db1,dp1,model2,da2,db2,dp2)
# black solid lines are the estimated ones
# dotted gray lines are the simulated ones
# the top lines are the true positive(left)/negative(right)
# the lowers lines are the false positive(left)/negative(right)


