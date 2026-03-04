### Two proportions
library("BayesPower")
library("crch")
library("ExtDist")
library("mvtnorm")
#### Setting up the simulation

# number of iteration
iter <- 5000

# Specification of analysis and design prior
threshold <- 3
a0 <- 1; b0 <- 1
a1 <- 1; b1 <- 1
a2 <- 1; b2 <- 1
prior_design_1 <- "beta"
a1d <- 2; b1d <- 2; dp1 <- .7
prior_design_2 <- "beta"
a2d <- 2; b2d <- 2; dp2 <- .4

# input validation
results <- BFpower.props(
  threshold= threshold,
  a0 = a0, b0 = b0,
  a1 = a1, b1 = b1,
  a2 = a2, b2 = b2,
  prior_design_1 = prior_design_1,
  a1d = a1d, b1d = b1d, dp1 = dp1,
  prior_design_2 = prior_design_2,
  a2d = a2d, b2d = b2d, dp2 = dp2,
  n1 = 50, n2 = 50
)


#### Under the null
theta1_h0 <- rbeta(iter, a0, b0)
theta2_h0 <- theta1_h0

#### Under the alternative
# "same" means analysis and design priors being the same
theta1_h1 <- switch(prior_design_1,
                    "beta" = {rbeta(iter, a1d, b1d)},
                    "Point"= rep(dp1,iter),
                    "same" = rbeta(iter, a1, b1))
theta2_h1 <- switch(prior_design_2,
                    "beta"={rbeta(iter, a2d, b2d)},
                    "Point"= rep(dp2,iter),
                    "same" = rbeta(iter, a2, b2))


sample_size  <- seq(10,1000,50)
r            <- 1
simulated    <- array(NA, dim=c(length(sample_size),4))
numeric_ones <- array(NA, dim=c(length(sample_size),4))
total   <- array(NA, dim = length(sample_size))

for (i in 1:length(sample_size)){
  n1 = sample_size[i]
  n2 = n1*r
  total[i] = n1 + n2

  ### simulating number of success per group and calculate BF
  # under the null
  x1_h0 <- sapply(1:iter, function(ii) {
    rbinom(1,n1,theta1_h0[ii])
  })
  x2_h0<-sapply(1:iter, function(ii) {
    rbinom(1,n2,theta2_h0[ii])
  })

  BF10_h0 <-sapply(1:iter, function(ii) {
    # BF10.props below is the exported function from the package
    # BF10.props(
    #  a0, b0, a1, b1, a2, b2,
    #  n1, n2 ,
    #  x1 = x1_h0[i],x2 = x2_h0[ii])$bf10
    # however using the underlying function for calculating BF is faster
    # since it would not calculate other irrelevent things e.g., p-value
    BayesPower:::BF10_p2(a0, b0, a1, b1, a2, b2,n1,n2,x1_h0[ii],x2_h0[ii])

  })

  # under the alternative
  x1_h1 <- sapply(1:iter, function(ii) {
    rbinom(1,n1,theta1_h1[ii])
  })
  x2_h1<-sapply(1:iter, function(ii) {
    rbinom(1,n2,theta2_h1[ii])
  })

  BF10_h1 <-sapply(1:iter, function(ii) {
    #BF10.props(
    #  a0, b0, a1, b1, a2, b2,
    #  n1, n2 ,
    #  x1 = x1_h1[ii],x2 = x2_h1[ii])$bf10

    BayesPower:::BF10_p2(a0, b0, a1, b1, a2, b2,n1,n2,x1_h1[ii],x2_h1[ii])

  })

  # --- operational characteristics ---
  TPR <- sum(BF10_h1>=threshold)/iter
  FPR <- sum(BF10_h0>=threshold)/iter

  TNR <- sum((1/BF10_h0)>=threshold)/iter
  FNR <- sum((1/BF10_h1)>=threshold)/iter


  simulated[i,] <- c(TPR,FPR,TNR,FNR)

  numeric_results<- BFpower.props(
    threshold= threshold,
    a0 = a0, b0 = b0,
    a1 = a1, b1 = b1,
    a2 = a2, b2 = b2,
    prior_design_1 = prior_design_1,
    a1d = a1d, b1d = b1d, dp1 = dp1,
    prior_design_2 = prior_design_2,
    a2d = a2d, b2d = b2d, dp2 = dp2,
    n1 = n1, n2 = n2
  )

  numeric_ones[i,] <-  c(numeric_results$results[[1]],
                         numeric_results$results[[4]],
                         numeric_results$results[[3]],
                         numeric_results$results[[2]])
}


par(mfrow = c(1, 2))
plot(total, numeric_ones[,1], type = "l",
     xlab = "Total sample size",
     ylab = "Probability",
     ylim = c(0,1),
     frame.plot = FALSE,
     lwd = 3,
     main = bquote(bold("Power curve for BF"[10]~">"~.(threshold))))
lines(total, simulated[,1], col = "gray", lwd = 2, lty = 3)
lines(total, numeric_ones[,2], col = "black", lwd = 3)
lines(total, simulated[,2], col = "grey", lwd = 2, lty = 3)

plot(total, numeric_ones[,3], type = "l",
     xlab = "Total sample size",
     ylab = "Probability",
     ylim = c(0,1),
     frame.plot = FALSE,
     lwd = 3,
     main = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold))))
lines(total, simulated[,3], col = "gray", lwd = 2, lty = 3)
lines(total, numeric_ones[,4], col = "black", lwd = 3)
lines(total, simulated[,4], col = "grey", lwd = 2, lty = 3)

