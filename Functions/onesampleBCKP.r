# The app is riddled with dt() warnings of the type:
#    In dt(t[i], df, ncp = delta * sqrt(df + 1)) :
#      full precision may not have been achieved in 'pnt{final}'
# This is why I used dnct() in the past. 
# I am bringing it back, slightly optimized. See dnct and pnct in app.r.
# We can later talk of possible alternatives.
# Since sourcing these functions takes time, what we can do it create a mini-package
# just with these functions. Loading this package would make these function immediately
# available, no more loading time required.

############# prior density function #############
# Probability density function of non-local prior:
dnlp <-function(delta, mu, ta){
  #((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
  dmom(delta-mu, ta^2) # 'mombf' package (includes pmom(), etc)
}
# Probability density function of informed t prior:
tstude <- function(t, location = 0, scale = sqrt(2)/2, df = 1) {
  # gamma((df+1)/2) * ((df+((t-location)/scale)^2)/df)^(-((df+1)/2)) / (scale*sqrt(df*pi)*gamma(df/2))
  dnct((t-location)/scale, df, ncp = 0)/scale 
}

t1_prior<- function(delta, location, scale, dff, model){
  switch(model,
         "Cauchy"         = tstude(delta, location, scale, 1),
         "Normal"         = dnorm (delta, location, scale),
         "NLP"            = dnlp  (delta, location, scale),
         "t-distribution" = tstude(delta, location, scale, dff))
}

############# the Bayes Factor #############

t1_BF10 <-function(t, df, model, location, scale, dff, hypothesis){
  bound  <- switch(hypothesis,
                   ">"  = c(a = 0, b = Inf),
                   "<"  = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x <- numeric(length(t))
  
  # Normalize the prior outside the for-loop:
  # normalization  <- integrate(function(delta) t1_prior(delta, location, scale, dff, model),lower = bound[1], upper = bound[2])$value
  # For all priors, the prior integrates to 1 when a = -Inf, b = Inf.
  # For all priors, we use their CDFs when either a = 0 or b = 0 and minimize manual integrations.
  # Note: pmom() errors at -Inf and Inf, so we avoid it below.
  normalization <- if (hypothesis == "!=") 1 else 
    switch(model,
           "Cauchy"         = pcauchy(bound[2], location, scale)     - pcauchy(bound[1], location, scale),
           "Normal"         = pnorm (bound[2], location, scale)      - pnorm (bound[1], location, scale),
           "NLP"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = pnct((bound[2] - location) / scale, dff, 0) - pnct((bound[1] - location) / scale, dff, 0))
  
  for(i in 1:length(t)){
    int  <- dnct(t[i], df, ncp = delta * sqrt(df+1)) * t1_prior(delta, location, scale, dff, model)/normalization

    # Removed stop.on.error = FALSE as it is bad form:
    x[i] <- integrate(int, lower = bound[1], upper = bound[2], rel.tol = 1e-8)$value / dt(t[i], df, ncp = 0)
  }
  return(x)
}



############# bound function  #############
# for finding the t value such that BF10 = D (code optimized):
t1_BF10_bound <- function(D, df, model, location, scale, dff, hypothesis) {
  Bound_finding <- function(t) t1_BF10(t, df, model, location, scale, dff, hypothesis) - D
  
  if (hypothesis == "!=") {
    x <- tryCatch(uniroot(Bound_finding, lower = -8, upper = 0)$root, error = function(e) NA)
    y <- tryCatch(uniroot(Bound_finding, lower =  0, upper = 8)$root, error = function(e) NA)
    results <- c(x, y)
  } else {
    x <- tryCatch(uniroot.all(Bound_finding, lower = -8, upper = 8), error = function(e) NA)
    results <- x
  }
  
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")
  
  BF.vals  <- t1_BF10(results, df, model, location, scale, dff, hypothesis)
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  
  return(results[BF.close])
}

# finding the t that correspond to BF01 = D is the same as 
# finding the t that corresponds to BF10 = 1/D:
t1_BF01_bound <- function(D, df, model, location, scale, dff, hypothesis) {
  t1_BF10_bound(1 / D, df, model, location, scale, dff, hypothesis)
}

# p(BF01>D|H0)
# t is the t-value lead to BF = b based on the bound functions (optimized):
t1_TNE <- function(t, df) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  if (length(t) == 2) return(pt(max(t), df) - pt(min(t), df))
  
  # length(t) = 1:
  return(if (t > 0) pt(t, df) else 1 - pt(t, df))
}


# p(BF10>D|H1)
# Argument 'hypothesis' is fully determined by the length and sign of the t values.
# I removed it as a function argument and compute it inside t1_TPE() instead.
t1_TPE <- function(t, df, model, location, scale, dff) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  hypothesis <- if (length(t) == 2) "!=" else if (t >= 0) ">" else "<"
  
  if (model == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(min(t), df, ncp) + (1 - pnct(max(t), df, ncp)))
    # Length 1:
    return(if (t >= 0) 1 - pnct(t, df, ncp) else pnct(t, df, ncp))
  }
  
  bound  <- switch(hypothesis,
                   ">"  = c(a = 0,    b = Inf),
                   "<"  = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf))
  
  normalization <- integrate(function(delta) t1_prior(delta, location, scale, dff, model),
                             lower = bound[1], upper = bound[2])$value
  
  int <- if (length(t) == 2) { # two-sided test
    function(delta) {
      pro1 <- 1 - pnct(max(t), df, delta * sqrt(df + 1))
      pro2 <-     pnct(min(t), df, delta * sqrt(df + 1))
      (pro1 + pro2) * t1_prior(delta, location, scale, dff, model) / normalization
    }
  } else if (t >= 0) { # one-sided test with delta > 0
    function(delta) (1 - pnct(t, df, delta * sqrt(df + 1))) * t1_prior(delta, location, scale, dff, model) / normalization
  } else {             # one-sided test with delta < 0
    function(delta) pnct(t, df, delta * sqrt(df + 1)) * t1_prior(delta, location, scale, dff, model) / normalization
  }
  
  # setting error value such that error are prevented:
  error <- if (model == "NLP" && scale < 0.3) 1e-14 else if (scale > 0.3) .Machine$double.eps^0.25 else 1e-8
  
  integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value
}


# p(BF01>D|H1)
# Similar as above:
t1_FNE <- function(t, df, model, location, scale, dff){
  
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  hypothesis <- if (length(t) == 2) "!=" else if (t >= 0) ">" else "<"
  
  if (model == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(max(t), df, ncp) - pnct(min(t), df, ncp))
    # Length 1:
    return(if (t >= 0) pnct(t, df, ncp) else 1 - pnct(t, df, ncp))
  }
  
  bound  <- switch(hypothesis,
                   ">"  = c(a = 0,    b = Inf),
                   "<"  = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf))
  
  # in some weird situation uniroot.all gives 4 t-values leading to BF =b 
  # but after some checking, only the middle two are the t-values leading to BF =b.
  # Jorge: This would be better dealt with inside t1_BF10_bound() and t1_BF01_bound()
  # I'll remove it and deal with the problem when I see it:
  # t <- as.numeric(t)
  # if (length(t) == 4) {  # Corrected condition check
  #   t = t[2:3]
  # }
  
  normalization <- integrate(function(delta) t1_prior(delta, location, scale, dff, model),
                             lower = bound[1], upper = bound[2])$value
  
  int <- if (length(t) == 2) { # two-sided test
    function(delta) {
      pro1 <- pnct(max(t), df, delta * sqrt(df + 1))
      pro2 <- pnct(min(t), df, delta * sqrt(df + 1))
      (pro1 - pro2) * t1_prior(delta, location, scale, dff, model) / normalization
    }
  } else if (t >= 0) { # one-sided test with delta > 0
    function(delta) pnct(t, df, delta * sqrt(df + 1)) * t1_prior(delta, location, scale, dff, model) / normalization
  } else {             # one-sided test with delta < 0
    function(delta) (1 - pnct(t, df, delta * sqrt(df + 1))) * t1_prior(delta, location, scale, dff, model) / normalization
  }
  
  # setting error value such that error are prevented:
  error <- if (model == "NLP" && scale < 0.3) 1e-14 else if (scale > 0.3) .Machine$double.eps^0.25 else 1e-8
  
  integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value
} 


# p(BF10>D|H0)
t1_FPE <- function(t, df) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  
  # if (length(t) == 4) t <- t[2:3]
  
  if (length(t) == 2) return(pt(min(t), df) + (1 - pt(max(t), df)))
  
  # length(t) = 1:
  return(if (t > 0) 1 - pt(t, df) else pt(t, df))
}



# Finding the degree of freedom that ensure p(BF10>D|H1) > targeted probability:
t1_N_finder <- function(D, target, model, location, scale, dff, hypothesis,
                        model_d, location_d, scale_d, dff_d, de_an_prior, alpha) {
  #de_an_prior: 1 = design prior and analysis priors are the same, otherwise different.

  # error prevention
  # sometimes, power can go higher than .8 with N= 2 already.
  # So, N should be returned now, otherwise, error will occur later.
  # Jorge: Below, I added design prior for N=2 too.
  lower <- 2
  upper <- 100000

  t2 <- t1_BF10_bound(D, df = lower, model, location, scale, dff, hypothesis)
  p2 <- if (de_an_prior == 1)
    t1_TPE(t2, df = lower, model, location, scale, dff) else
      t1_TPE(t2, df = lower, model_d, location_d, scale_d, dff_d)
  if (p2 > target) return(lower)

  Power_root <- function(df) {
    t <- t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)
    if (de_an_prior == 1)
      t1_TPE(t, df, model, location, scale, dff) - target else
        t1_TPE(t, df, model_d, location_d, scale_d, dff_d) - target
  }

  ## finding the required df, i will do the plus one to get the N in the later function.
  # Jorge: It's a pity if we don't fix it here already, super easy to do right now.
  #        So I did it, see below. I'll adapt latter functions if needed be.

  # Jorge: 'df.power' makes for a more accurate name than 'N'.message("Power at lower = ", Power_root(lower))
  message("lower = ", lower)
  message("upper = ", upper)
  message("Power at lower = ", Power_root(lower))
  message("Power at upper = ", Power_root(upper))
  df.power <- uniroot(Power_root, lower = lower, upper = upper)$root

  ## checking if the N lead to an acceptable alpha level
  t   <- t1_BF10_bound(D, df.power, model, location, scale, dff, hypothesis)
  FPE <- t1_FPE(t, df.power)
  if (FPE <= alpha) return(df.power + 1)

  # if the FPE > alpha, then we search for another df
  # Jorge: 'alpha.root' is better than 'alpha_bound'.
  alpha.root <- function(df) {
    t <- t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)
    t1_FPE(t, df) - alpha
  }

  # Jorge: 'df.alpha' is better than 'NN'.
  df.alpha <- uniroot(alpha.root, lower = df.power, upper = upper)$root
  return(df.alpha + 1)

  #N = auto_uniroot_fixed_lower_going_up(Power_root,fixed_lower = 2, upper = 2000, step = 500, max_attempts = 16)
}



############ probability table
# Jorge: I edited so that it used N returned by t1_N_finder().
t1_Table <- function(D, target, model, location, scale, dff, hypothesis,
                     model_d, location_d, scale_d, dff_d, de_an_prior, N, mode_bf, alpha) {
  # mode_bf == "0" means that the design analysis is done for a fixed N
  # Otherwise, it searches N where power > targeted power with FPE < FP
  if (mode_bf == "0") df <- N-1 else {
    N  <- ceiling(t1_N_finder(D, target, model, location, scale, dff, hypothesis,
                      model_d, location_d, scale_d, dff_d, de_an_prior, alpha))
    df <- N -1
  }

  # t bounds:
  t10 <- t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)
  t01 <- t1_BF01_bound(D, df, model, location, scale, dff, hypothesis)

  # max BF10 possible:
  max_BF <- 1 / t1_BF10(0, df, model, location, scale, dff, hypothesis)
  BF_D   <- t10

  # FPE and TPE:
  FPE       <- t1_FPE(t10, df)
  if (de_an_prior == 1) {
    TPE       <- t1_TPE(t10, df, model, location, scale, dff)
    TPR_model <- model
    TPR_loc   <- location
    TPR_scale <- scale
    TPR_dff   <- dff
  } else {
    TPE       <- t1_TPE(t10, df, model_d, location_d, scale_d, dff_d)
    TPR_model <- model_d
    TPR_loc   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }

  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNR <- 0
  } else {
    FNE <- t1_FNE(t01, df, TPR_model, TPR_loc, TPR_scale, TPR_dff)
    TNR <- t1_TNE(t01, df)
  }

  # table:
  tab.names <- c(
    sprintf("p(BF10 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H0)", D),
    sprintf("p(BF10 > %0.f | H0)", D),
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNR, FPE, N, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}

# For plotting, compute normalized prior density over tt:
compute.prior.density <- function(tt, model, location, scale, dff, hypothesis) {
  if (model == "Point") return(rep(NA, length(tt)))
  bounds <- switch(hypothesis,
                   ">"  = c(0, Inf),
                   "<"  = c(-Inf, 0),
                   "!=" = c(-Inf, Inf))
  norm <- integrate(function(delta) t1_prior(delta, location, scale, dff, model),
                    lower = bounds[1], upper = bounds[2])$value
  t1_prior(tt, location, scale, dff, model) / norm
}

# plot for the selected prior 
t1_prior_plot <- function(D, target, model, location, scale, dff, hypothesis,
                          model_d, location_d, scale_d, dff_d, de_an_prior) {
  par(mfrow = c(1, 1))
  
  plot.bounds    <- switch(hypothesis,
                           ">"  = c(0, 5),
                           "<"  = c(-5, 0),
                           "!=" = c(-5, 5))
  tt             <- seq(plot.bounds[1], plot.bounds[2], 0.01)
  prior.analysis <- compute.prior.density(tt, model, location, scale, dff, hypothesis)
  prior.design   <- if (de_an_prior == 0 && model_d != "Point") 
    compute.prior.density(tt, model_d, location_d, scale_d, dff_d, hypothesis) else 
      rep(NA, length(tt))
  ylim.max <- max(prior.analysis, prior.design, na.rm = TRUE)
  
  # Base plot:
  plot(tt, prior.analysis, type = "l", lwd = 2,
       xlab = expression(bold(delta)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~delta~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))
  
  # If design prior != analysis prior:
  if (de_an_prior == 0) {
    if (model_d == "Point") 
      arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim_max, length = 0.1, col = "black", lty = 2) else 
        lines(tt, prior.design, lty = 2)
    
    # Add legend:
    legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }
}


# Plot BF10 and BF01 vs. t-values:
bf10_t1 <-function(D = 3, df, target, model = "NA", location = 0, scale = 0.707, dff = 1, hypothesis) {
  tt <- seq(-5, 5, 0.2)
  
  # Compute BF10 and t-bounds:
  BF10   <- t1_BF10(tt, df, model, location, scale, dff, hypothesis)
  t.BF10 <- t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)
  
  par(mfrow = c(1, 2))
  # Left plot - BF10:
  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(D)~"when t = "~.(format(t.BF10, digits = 4))))
  } else {
    bquote(bold("BF"[10]~"="~.(D)~"when t = "~.(format(t.BF10[1], digits = 4))~"or"~.(format(t.BF10[2], digits = 4))))
  }
  plot(tt, BF10, type = "l", log = "y", xlab = "t-value", ylab = expression("BF"[10]),
       main = main.bf10, frame.plot = FALSE, xaxt = "n")
  abline(v = t.BF10)
  axis(1, c(-5, 5))
  if (length(t.BF10)) axis(1, round(t.BF10, 2))
  
  # Left plot - BF01:
  BF01   <- 1 / BF10
  t.BF01 <- t1_BF01_bound(D, df, model, location, scale, dff, hypothesis)
  
  # Check if BF01 = D is possible:
  max.BF01   <- 1 / t1_BF10(0, df, model, location, scale, dff, hypothesis = "!=")
  impossible <- (hypothesis == "!=") && (max.BF01 < D || identical(t.BF01, "bound cannot be found"))
  
  plot(tt, BF01, type = "l", log = "y", xlab = "t-value", ylab = bquote("BF"[01]),
       main = "", frame.plot = FALSE, xaxt = "n")
  axis(1, c(-5, 5))
  if (impossible) {
    title(main = bquote(bold("It is impossible to have BF"[01]~"="~.(D))))
  } else {
    abline(v = t.BF01)
    axis(1, round(t.BF01, 2))
    main.bf01 <- if (length(t.BF01) == 1) {
      bquote(bold("BF"[01]~"="~.(D)~"when t = "~.(format(t.BF01, digits = 4))))
    } else {
      bquote(bold("BF"[01]~"="~.(D)~"when t = "~.(format(t.BF01[1], digits = 4))~"or"~.(format(t.BF01[2], digits = 4))))
    }
    title(main = main.bf01)
  }
}

# Power curve function for BF10 > D under H1:
Power_t1 <- function(D, model, location, scale, dff, hypothesis,
                     model_d, location_d, scale_d, dff_d,
                     de_an_prior, N) {
  
  # df range to evaluate power:
  df.min     <- 2
  df.max     <- ceiling(N * 1.2)
  dfs        <- seq(df.min, df.max, length.out = 31)
  power.vals <- numeric(length(dfs))
  
  for (i in seq_along(dfs)) {
    t <- t1_BF10_bound(D, dfs[i], model, location, scale, dff, hypothesis)
    # Choose correct design prior:
    power.vals[i] <- if (de_an_prior == 1) 
      t1_TPE(t, dfs[i], model, location, scale, dff) else 
        t1_TPE(t, dfs[i], model_d, location_d, scale_d, dff_d)
  }

  plot(dfs + 1, power.vals, type = "l",
       xlab = "Sample size",
       ylab = bquote("p(BF"[10]~">"~.(D)~"| H1)"),
       ylim = c(0, 1), frame.plot = FALSE,
       main = bquote(bold("Power curve for BF"[10]~">"~.(D))))
}
