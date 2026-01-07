#' @useDynLib BayesPower, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom rlang .data
#' @importFrom utils globalVariables
NULL

# Avoid R CMD check notes
utils::globalVariables(c(
  "SampleSize", "Probability", "Type",
  "True Positive", "False Positive",
  "True Negative", "False Negative",
  "x", "BF","r","f"
))

fmt_val <- function(x) {
  if (is.numeric(x) && length(x) == 1) return(as.character(x))
  if (is.numeric(x) && length(x) > 1) return(paste(x, collapse = ", "))
  if (is.character(x)) return(shQuote(x))
  return(as.character(x))
}
# ---- ANOVA.r ----

# k = number of predictor in the full model
# p = number  of predictor in the reduced model
# m = N-p
# q = k-p

F_prior<- function(fsq,q,dff,rscale,f,model) {


  switch(model,
         "effectsize" = {gamma((q + dff) / 2) / gamma(dff / 2) /gamma(q / 2) *
    (dff * rscale^2)^(dff / 2) * fsq^(q / 2 - 1) *
    (dff * rscale^2 + f^2 + fsq)^(-dff / 2 - q / 2) *
    hypergeo::genhypergeo(c((dff + q) / 4, (2 + dff + q) / 4),
                          q / 2, 4 * f^2 * fsq / (dff * rscale^2 + f^2 + fsq)^2)},
          "Moment" = { temp <- f^2 * (dff + q - 2)/2

          gamma((q + dff) / 2) / gamma(dff / 2) / gamma(q / 2) *
            2 * (dff - 2) / q / (dff-2 + q) / f^2 *
            fsq^(q/2) * temp^(dff/2) * (temp + fsq)^(-(dff+q)/2)})

}


F_BF <- function(f, q, m, dff, rscale, f_m, model) {
  sapply(f, function(fi) {
    int <- function(fsq) {
      stats::df(fi, q, m - q, ncp = m * fsq) * F_prior(fsq, q, dff, rscale, f_m, model)
    }
    lh1 <- stats::integrate(int, lower = 0, upper = Inf, stop.on.error = FALSE, rel.tol = 1e-4)$value
    lh0 <- stats::df(fi, q, m - q)
    lh1 / lh0
  })
}


F_BF_bound_10 <-function(D,q,m,dff,rscale,f_m,model){
  x = numeric(0)
  Bound_finding <-function(f){
    F_BF(f,q,m,dff,rscale,f_m,model)-D
  }

  x = tryCatch( stats::uniroot(Bound_finding,lower=0.01,upper = 40 )$root, error=function(e){})
  if (length(x) == 0) return("no bound is found")

  return(x)
}

F_BF_bound_01 <-function(D,q,m,dff,rscale,f_m,model){
  F_BF_bound_10(1/D,q,m,dff,rscale,f_m,model)
}


F_TPE<-function(f,q,m,dff,rscale,f_m,model){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  if (model == "Point"){
    x = stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = F)
    return(x)
  }
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,model)
  }
  x = stats::integrate(int,lower = 0,upper = Inf)$value
    return(x)
}

F_FNE<-function(f,q,m,dff,rscale,f_m,model){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  if (model == "Point"){

    x = stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = T)
    return(x)
  }
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,model)
  }
  x = stats::integrate(int,lower = 0,upper = Inf)$value
  return(x)
}

F_TNE<-function(f,q,m){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)
    x = stats::pf(f,q,m-q,ncp =0,lower.tail = T)
  return(x)
}

F_FPE<-function(f,q,m){

  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  x = stats::pf(f,q,m-q,ncp =0,lower.tail = F)
  return(x)
}

f_N_finder<-function(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,FP){
  q= k-p
  lower = 2*k-p+1
  m= lower-p
  f <- F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  p2 <- if (de_an_prior == 1)
    F_TPE(f,q,m,dff,rscale,f_m,model) else
      F_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d)
  if (p2 > target) return(lower)

  Power_root <- function(n){
    m= n-p
    f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)

    pro <- if (de_an_prior == 1)
      F_TPE(f,q,m,dff,rscale,f_m,model) else
        F_TPE(f,q,m,dff_d,rscale_d,f_m_d,model_d)

    return(pro-target)
  }

  N.power = stats::uniroot(Power_root,lower = lower,upper =  10000)$root
  m= N.power-p
  f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  FPE = F_FPE(f,q,m)

  if (FPE <= FP) return(N.power + 1)
  alpha_root <- function(n){
    m= n-p
    f = F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
    pro = F_FPE(f,q,m)


    return(pro-FP)
  }
  N.alpha = stats::uniroot(alpha_root,lower = N.power,upper =  10000)$root
  return(N.alpha)
}

f_N_01_finder<-function(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,FP){
  q= k-p
  lower = 2*k-p+1
  m= lower-p
  upper =  10000
  f <- F_BF_bound_01(D,q,m,dff,rscale,f_m,model)
  TNE_lo <- F_TNE(f,q,m)
  FNE_lo <- if (de_an_prior == 1)
    F_FNE(f,q,m,dff,rscale,f_m,model) else
      F_FNE(f,q,m,dff_d,rscale_d,f_m_d,model_d)

  if (TNE_lo > target && FNE_lo < FP) {
    return(lower)
  } else if (TNE_lo > target) {
    FN.root <- function(n) {
      q= k-p
      m= n-p
      f <- F_BF_bound_01(D,q,m,dff,rscale,f_m,model)
      FNE <- if (de_an_prior == 1)
        F_FNE(f,q,m,dff,rscale,f_m,model) else
          F_FNE(f,q,m,dff_d,rscale_d,f_m_d,model_d)
      FNE - FP
    }
    return(stats::uniroot(FN.root, lower = lower, upper = upper)$root)
  }

  TN_root <- function(n){
    m= n-p
    f = F_BF_bound_01(D,q,m,dff,rscale,f_m,model)
    TNE <- F_TNE(f,q,m)

    return(TNE-target)
  }

  N.TN = stats::uniroot(TN_root,lower = lower,upper =  upper)$root
  m= N.TN-p
  f = F_BF_bound_01(D,q,m,dff,rscale,f_m,model)
  FNE = if (de_an_prior == 1)
    F_FNE(f,q,m,dff,rscale,f_m,model) else
      F_FNE(f,q,m,dff_d,rscale_d,f_m_d,model_d)

  if (FNE <= FP) return(N.TN + 1)

  FN.root <- function(n) {
    q= k-p
    m= n-p
    f <- F_BF_bound_01(D,q,m,dff,rscale,f_m,model)
    FNE <- if (de_an_prior == 1)
      F_FNE(f,q,m,dff,rscale,f_m,model) else
        F_FNE(f,q,m,dff_d,rscale_d,f_m_d,model_d)
    FNE - FP
  }

  N.FN = stats::uniroot(FN.root, lower = N.TN, upper = upper)$root
  return(N.FN)
}

f_table<-function(D,target,p,k,dff,rscale,f_m,model,
                  dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, mode_bf,FP,direct ){


  if (mode_bf == 1){

    n = switch(direct,
               "h1"= ceiling(f_N_finder(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,FP )),
               "h0" = ceiling(f_N_01_finder(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,FP )))
  } else {
    n=n
  }
  q= k-p
  m= n-p

  # f bounds:
  f10 <- F_BF_bound_10(D,q,m,dff,rscale,f_m,model)
  f01 <- F_BF_bound_01(D,q,m,dff,rscale,f_m,model)

  # max BF10 possible:
  max_BF <- 1/F_BF(0.00001,q,m,dff,rscale,f_m,model)
  BF_D   <- f10

  # FPE and TPE:
  FPE       <- F_FPE(f10,q,m)
  if (de_an_prior == 1) {
    TPE       <- F_TPE(f10,q,m,dff,rscale,f_m,model)
    TPR_dff   <- dff
    TPR_rscale<- rscale
    TPR_f_m   <- f_m
    TPR_model <- model
  } else {
    TPE       <- F_TPE(f10,q,m,dff_d,rscale_d,f_m_d,model_d)
    TPR_dff   <- dff_d
    TPR_rscale<- rscale_d
    TPR_f_m   <- f_m_d
    TPR_model <- model_d
  }

  # FNE and TNE:
  if ( max_BF < D | BF_D == "no bound is found") {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- F_FNE(f01,q,m,TPR_dff,TPR_rscale,TPR_f_m,TPR_model)
    TNE <- F_TNE(f01,q,m)

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

prior_plot_f <-function(q,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior ){
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 1))
  fsq = seq(0.001,15,.2)

  prior.analysis = F_prior(fsq,q,dff,rscale,f_m,model)
  prior.design   <- if (de_an_prior == 0 && model_d != "Point")
    F_prior(fsq,q,dff_d,rscale_d,f_m_d,model_d) else
      rep(NA, length(fsq))
  ylim.max <- max(prior.analysis, prior.design, na.rm = TRUE)


  # Base plot:
  plot(fsq, prior.analysis, type = "l", lwd = 2,
       xlab = expression(bold(lambda^2)),
       ylab = "density",
       main = bquote(bold("prior distribution on "~lambda^2~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))

  # If design prior != analysis prior:
  if (de_an_prior == 0) {
    if (model_d == "Point")
      graphics::arrows(x0 = f_m, y0 = 0, x1 = f_m, y1 = ylim.max, length = 0.1, col = "black", lty = 2) else
        graphics::lines(fsq, prior.design, lty = 2)

    # Add legend:
    graphics::legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }

}

bf10_f <- function(D, n, k, p, dff, rscale, f_m, model) {

  q <- k - p
  m <- n - p
  ff <- seq(0.01, 10, 0.05)

  # BF10 values and bounds
  BF10 <- F_BF(ff, q, m, dff, rscale, f_m, model)
  f.BF10 <- F_BF_bound_10(D, q, m, dff, rscale, f_m, model)

  df10 <- data.frame(f = ff, BF = BF10)

  main.bf10 <- if (length(f.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(D)~" when f = "~.(round(f.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(D)~" when f = "~.(round(f.BF10[1], 2))~" or "~.(round(f.BF10[2], 2))))
  }

  # x-axis breaks for BF10
  x_breaks_10 <- sort(unique(c(0, 10, round(f.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(f, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = f.BF10, linetype = "dashed") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(0, 10), breaks = x_breaks_10) +
    ggplot2::labs(
      x = "f-value",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # BF01 values and bounds
  BF01 <- 1 / BF10
  f.BF01 <- F_BF_bound_01(D, q, m, dff, rscale, f_m, model)
  max_BF01 <- 1 / F_BF(0.001, q, m, dff, rscale, f_m, model)
  impossible <- any(max_BF01 < D | f.BF01 == "bound cannot be found")

  df01 <- data.frame(f = ff, BF = BF01)

  if (impossible) {
    p2 <- ggplot2::ggplot(df01, ggplot2::aes(f, BF)) +
      ggplot2::geom_line(size = 1.2, color = "black") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(0, 10), breaks = c(0, 10)) +
      ggplot2::labs(
        x = "f-value",
        y = expression("BF"[0][1] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[0][1]~"="~.(D)))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  } else {
    main.bf01 <- if (length(f.BF01) == 1) {
      bquote(bold("BF"[0][1]~"="~.(D)~" when f = "~.(round(f.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(D)~" when f = "~.(round(f.BF01[1], 2))~" or "~.(round(f.BF01[2], 2))))
    }

    x_breaks_01 <- sort(unique(c(0, 10, round(f.BF01, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(f, BF)) +
      ggplot2::geom_line(size = 1.2, color = "black") +
      ggplot2::geom_vline(xintercept = f.BF01, linetype = "dashed") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(0, 10), breaks = x_breaks_01) +
      ggplot2::labs(
        x = "f-value",
        y = expression("BF"[0][1] * " (log scale)"),
        title = main.bf01
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }

  # Combine plots side by side
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

Power_f <- function(D, k, p, dff, rscale, f_m, model,
                    k_d, p_d, dff_d, rscale_d, f_m_d, model_d,
                    de_an_prior, N) {

  # Sample size range
  smin <- (2 * k - p + 1)
  smax <- N * 1.2
  n <- seq(smin, smax, length.out = 31)
  q <- k - p
  m <- n - p

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(n))

  for (i in seq_along(n)) {
    f10 <- F_BF_bound_10(D, q, m[i], dff, rscale, f_m, model)
    f01 <- F_BF_bound_01(D, q, m[i], dff, rscale, f_m, model)

    # True Positive
    TPE[i] <- if (de_an_prior == 1) {
      F_TPE(f10, q, m[i], dff, rscale, f_m, model)
    } else {
      F_TPE(f10, q, m[i], dff_d, rscale_d, f_m_d, model_d)
    }

    # False Positive / True Negative
    FPE[i] <- F_FPE(f10, q, m[i])
    TNE[i] <- F_TNE(f01, q, m[i])

    # False Negative
    FNE[i] <- if (de_an_prior == 1) {
      F_FNE(f01, q, m[i], dff, rscale, f_m, model)
    } else {
      F_FNE(f01, q, m[i], dff_d, rscale_d, f_m_d, model_d)
    }
  }

  # Prepare data for ggplot
  df_bf10 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = n,
      `True Positive` = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df_bf10$Type <- factor(df_bf10$Type, levels = c("True Positive", "False Positive"))

  df_bf01 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = n,
      `True Negative` = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df_bf01$Type <- factor(df_bf01$Type, levels = c("True Negative", "False Negative"))

  # Colors for lines
  type_colors <- c(
    "True Positive" = "black",
    "False Positive" = "grey50",
    "True Negative" = "black",
    "False Negative" = "grey50"
  )

  # Clean theme
  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Legend theme
  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  # BF10 plot
  p1 <- ggplot2::ggplot(df_bf10, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # BF01 plot
  p2 <- ggplot2::ggplot(df_bf01, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots side by side
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

# ---- ANOVAe.r ----
Fe_BF <- function(f, q, m, dff, rscale, f_m, model, e) {
  # Compute normalizations once
  normalizationh1  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = e,upper = Inf,rel.tol = 1e-10)$value
  normalizationh0 <- 1 - normalizationh1

  # Define likelihood ratio function
  sapply(f, function(fi) {
    int1 <- function(fsq) {
      stats::df(fi, q, m - q, ncp = m * fsq) * F_prior(fsq,q,dff,rscale,f_m,model) / normalizationh1
    }
    int0 <- function(fsq) {
      stats::df(fi, q, m - q, ncp = m * fsq) * F_prior(fsq,q,dff,rscale,f_m,model) / normalizationh0
    }
    lh1 <- stats::integrate(int1, lower = e, upper = Inf, stop.on.error = FALSE)$value
    lh0 <- stats::integrate(int0, lower = 0, upper = e,   stop.on.error = FALSE)$value
    lh1 / lh0
  })
}

Fe_BF_bound_10 <-function(D,q,m,dff,rscale,f_m,model,e){
  x = numeric(0)
  Bound_finding <-function(f){
    Fe_BF(f,q,m,dff,rscale,f_m,model,e)-D
  }
  #x = tryCatch( stats::uniroot(Bound_finding,lower=0.01,upper = 100 )$root, error=function(e){})
  x = tryCatch( rootSolve::uniroot.all(Bound_finding,lower=0.01,upper = 500 ), error=function(e){})

  if (length(x) == 0) return("no bound is found")

  return(x)
}

Fe_BF_bound_01 <-function(D,q,m,dff,rscale,f_m,model,e){
  Fe_BF_bound_10(1/D,q,m,dff,rscale,f_m,model,e)
}

Fe_TPE<-function(f,q,m,dff,rscale,f_m,model,e){


  if (length(f) == 0 || any(f == "no bound is found")) return(0)


  if (model == "Point") return(stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = F))

  normalizationh1  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = e,upper = Inf,rel.tol = 1e-5)$value
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh1
  }
  x = stats::integrate(int,lower = e,upper = Inf)$value
    return(x)
}

Fe_FNE<-function(f,q,m,dff,rscale,f_m,model,e){

  if (length(f) == 0 || any(f == "no bound is found")) return(0)


  if (model == "Point") return(stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = T))

  normalizationh1  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = e,upper = Inf,rel.tol = 1e-10)$value
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh1
  }
  x = stats::integrate(int,lower = e,upper = Inf)$value
  return(x)
}

Fe_TNE<-function(f,q,m,dff,rscale,f_m,model,e){

  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  normalizationh0  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = 0,upper = e,rel.tol = 1e-5)$value

  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh0
  }
  x = stats::integrate(int,lower = 0,upper = e)$value
  return(x)
}

Fe_FPE<-function(f,q,m,dff,rscale,f_m,model,e){

  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  normalizationh0  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,model),lower = 0,upper = e,rel.tol = 1e-10)$value


  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,model)/normalizationh0
  }
  x = stats::integrate(int,lower = 0,upper = e)$value
  return(x)
}

fe_N_finder<-function(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,e,FP ){
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

  #N.power <- stats::uniroot(Power_root,lower = lower,upper =  5000)$root
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
  N.alpha <- robust_uniroot(alpha_root,lower=N.power)
    #stats::uniroot(alpha_root,lower = N.power,upper =  5000)$root
    return(N.alpha)
}

fe_N_01_finder<-function(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,e,FP ){
  q     <- k-p
  lower <- 2*k-p+1
  m     <- lower-p
  upper <-  5000
  f     <- Fe_BF_bound_01(D,q,m,dff,rscale,f_m,model,e)
  TNE_lo <- Fe_TNE(f,q,m,dff,rscale,f_m,model,e)
  FNE_lo <- if (de_an_prior == 1)
    Fe_FNE(f,q,m,dff,rscale,f_m,model,e) else
      Fe_FNE(f,q,m,dff_d,rscale_d,f_m_d,model_d,e)

  if (TNE_lo > target && FNE_lo < FP) {
    return(lower)
  } else if (TNE_lo > target) {
    FN.root <- function(n) {
      q     <- k-p
      m     <- n-p
      f <- Fe_BF_bound_01(D,q,m,dff,rscale,f_m,model,e)
      FNE <- if (de_an_prior == 1)
        Fe_FNE(f,q,m,dff,rscale,f_m,model,e) else
          Fe_FNE(f,q,m,dff_d,rscale_d,f_m_d,model_d,e)
      FNE - FP
    }
    return(stats::uniroot(FN.root, lower = lower, upper = upper)$root)
  }

  TN_root <- function(n){
    m   <- n-p
    f   <- Fe_BF_bound_01(D,q,m,dff,rscale,f_m,model,e)
    pro <- Fe_TNE(f,q,m,dff,rscale,f_m,model,e)
    return(pro-target)
  }

  #N.TN <- stats::uniroot(Power_root,lower = lower,upper =  5000)$root
  N.TN <-robust_uniroot(TN_root,lower=lower)
  m       <- N.TN-p
  f       <- Fe_BF_bound_01(D,q,m,dff,rscale,f_m,model,e)
  FNE     <- if (de_an_prior == 1)
    Fe_FNE(f,q,m,dff,rscale,f_m,model,e) else
      Fe_FNE(f,q,m,dff_d,rscale_d,f_m_d,model_d,e)

  if (FNE <= FP) return(N.TN)

  FN.root <- function(n) {
    q     <- k-p
    m     <- n-p
    f <- Fe_BF_bound_01(D,q,m,dff,rscale,f_m,model,e)
    FNE <- if (de_an_prior == 1)
      Fe_FNE(f,q,m,dff,rscale,f_m,model,e) else
        Fe_FNE(f,q,m,dff_d,rscale_d,f_m_d,model_d,e)
    FNE - FP
  }
  N.FN <- robust_uniroot(FN.root,lower=N.TN)
  #stats::uniroot(alpha_root,lower = N.TN,upper =  5000)$root
  return(N.FN)
}


fe_table<-function(D,target,p,k,dff,rscale,f_m,model,
                  dff_d,rscale_d,f_m_d,model_d,de_an_prior,n, mode_bf,e ,FP,direct){

  if (mode_bf == 1){

    n = switch(direct,
               "h1" = ceiling(fe_N_finder(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,e ,FP)),
               "h0" = ceiling(fe_N_01_finder(D,target,p,k,dff,rscale,f_m,model,dff_d,rscale_d,f_m_d,model_d,de_an_prior,e ,FP)))
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
    TPE       <- Fe_TPE(f10,q,m,dff_d,rscale_d,f_m_d,model_d,e)
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
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 1))

  normalization  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f,model),lower = e,upper = Inf,rel.tol = 1e-10)$value

  fsq = seq(0.001,20,.05)
  #prior.analysis.h1 = F_prior(fsq,q,dff,rscale,f,model)/normalization
  prior.analysis.h1 = F_prior(fsq,q,dff,rscale,f,model)
  prior.analysis.h1[fsq<e]=0
  #prior.analysis.h0 = F_prior(fsq,q,dff,rscale,f,model)/(1-normalization)
  prior.analysis.h0 = F_prior(fsq,q,dff,rscale,f,model)
  prior.analysis.h0[fsq>e]=0
  prior.design <- if (de_an_prior == 0 && model_d != "Point"){
    #normalization_d  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f,model),lower = e,upper = Inf,rel.tol = 1e-10)$value

    #F_prior(fsq,q,dff_d,rscale_d,f_m_d,model_d)/normalization_d
    F_prior(fsq,q,dff_d,rscale_d,f_m_d,model_d)}else
      rep(NA, length(fsq))
  ylim.max <- max(prior.analysis.h1, prior.analysis.h0, prior.design, na.rm = TRUE)

  plot(fsq, prior.analysis.h1, type = "l", lwd = 2,
       xlab = bquote(bold(lambda^2)),
       ylab = "density",
       main  = bquote(bold("prior distribution on "~lambda^2~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))

  graphics::lines(fsq, prior.analysis.h0, lty = 2, col = "black", lwd = 2)

  # Optional: design prior
  legend.labels <- c("H1 - Analysis Prior", "H0 - Analysis Prior")
  legend.cols   <- c("black", "black")
  legend.lty    <- c(1, 2)
  legend.lwd    <- c(2, 2)

  if (de_an_prior == 0) {
    if (model_d == "Point") {
      graphics::arrows(x0 = f_m_d, y0 = 0, x1 = f_m_d, y1 = ylim.max,
             length = 0.1, col = "gray", lty = 2,lwd = 3)
    } else {
      graphics::lines(fsq, prior.design, lty = 1, col = "gray", lwd = 3)
    }

    # Add design prior to legend
    legend.labels <- c(legend.labels, "Design prior")
    legend.cols   <- c(legend.cols, "gray")
    legend.lty    <- c(legend.lty, 1)
    legend.lwd    <- c(legend.lwd, 2)
  }
  graphics::legend("topleft",
         legend = legend.labels,
         col = legend.cols,
         lty = legend.lty,
         lwd = legend.lwd,
         bty = "n")


}

bf10_fe <- function(D, n, k, p, dff, rscale, f_m, model, e) {

  q <- k - p
  m <- n - p
  ff <- seq(0.01, 10, 0.05)

  # BF10 values and bounds
  BF10 <- Fe_BF(ff, q, m, dff, rscale, f_m, model, e)
  f.BF10 <- Fe_BF_bound_10(D, q, m, dff, rscale, f_m, model, e)

  df10 <- data.frame(f = ff, BF = BF10)

  main.bf10 <- if (length(f.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(D)~" when f = "~.(round(f.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(D)~" when f = "~.(round(f.BF10[1], 2))~" or "~.(round(f.BF10[2], 2))))
  }

  # x-axis breaks for BF10
  x_breaks_10 <- sort(unique(c(0, 10, round(f.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(f, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = f.BF10, linetype = "dashed") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(0, 10), breaks = x_breaks_10) +
    ggplot2::labs(
      x = "f-value",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # BF01 values and bounds
  BF01 <- 1 / BF10
  f.BF01 <- Fe_BF_bound_01(D, q, m, dff, rscale, f_m, model, e)
  max_BF01 <- 1 / Fe_BF(0.001, q, m, dff, rscale, f_m, model, e)
  impossible <- any(max_BF01 < D | f.BF01 == "bound cannot be found")

  df01 <- data.frame(f = ff, BF = BF01)

  if (impossible) {
    p2 <- ggplot2::ggplot(df01, ggplot2::aes(f, BF)) +
      ggplot2::geom_line(size = 1.2, color = "black") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(0, 10), breaks = c(0, 10)) +
      ggplot2::labs(
        x = "f-value",
        y = expression("BF"[0][1] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[0][1]~"="~.(D)))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  } else {
    main.bf01 <- if (length(f.BF01) == 1) {
      bquote(bold("BF"[0][1]~"="~.(D)~" when f = "~.(round(f.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(D)~" when f = "~.(round(f.BF01[1], 2))~" or "~.(round(f.BF01[2], 2))))
    }

    x_breaks_01 <- sort(unique(c(0, 10, round(f.BF01, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(f, BF)) +
      ggplot2::geom_line(size = 1.2, color = "black") +
      ggplot2::geom_vline(xintercept = f.BF01, linetype = "dashed") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(0, 10), breaks = x_breaks_01) +
      ggplot2::labs(
        x = "f-value",
        y = expression("BF"[0][1] * " (log scale)"),
        title = main.bf01
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }

  # Combine plots side by side
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}


Power_fe <- function(D, k, p, dff, rscale, f_m, model,
                     k_d, p_d, dff_d, rscale_d, f_m_d, model_d,
                     de_an_prior, N, e) {

  # Sample size range
  smin <- (2 * k - p + 1)
  smax <- N * 2
  sdf <- seq(smin, smax, length.out = 31)
  q <- k - p
  m <- sdf - p

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(sdf))

  for (i in seq_along(sdf)) {
    f10 <- Fe_BF_bound_10(D, q, m[i], dff, rscale, f_m, model, e)
    f01 <- Fe_BF_bound_01(D, q, m[i], dff, rscale, f_m, model, e)

    # True Positive
    TPE[i] <- if (de_an_prior == 1) {
      Fe_TPE(f10, q, m[i], dff, rscale, f_m, model, e)
    } else {
      Fe_TPE(f10, q, m[i], dff_d, rscale_d, f_m_d, model_d, e)
    }

    # False Negative
    FNE[i] <- if (de_an_prior == 1) {
      Fe_FNE(f01, q, m[i], dff, rscale, f_m, model, e)
    } else {
      Fe_FNE(f01, q, m[i], dff_d, rscale_d, f_m_d, model_d, e)
    }

    # False Positive / True Negative
    FPE[i] <- Fe_FPE(f10, q, m[i], dff, rscale, f_m, model, e)
    TNE[i] <- Fe_TNE(f01, q, m[i], dff, rscale, f_m, model, e)
  }

  # Prepare data for ggplot
  df_bf10 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = sdf,
      `True Positive` = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df_bf10$Type <- factor(df_bf10$Type, levels = c("True Positive", "False Positive"))

  df_bf01 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = sdf,
      `True Negative` = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df_bf01$Type <- factor(df_bf01$Type, levels = c("True Negative", "False Negative"))

  # Colors for lines
  type_colors <- c(
    "True Positive" = "black",
    "False Positive" = "grey50",
    "True Negative" = "black",
    "False Negative" = "grey50"
  )

  # Clean theme
  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Legend theme
  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  # BF10 plot
  p1 <- ggplot2::ggplot(df_bf10, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # BF01 plot
  p2 <- ggplot2::ggplot(df_bf01, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots side by side
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

# ---- binomial.r ----
adjust_root_10 <- function(root, n, alpha, beta, location, scale, model, hypothesis, D) {
  # If root is less than 0, return NA
  if (root < 0) return(NA)

  # Evaluate BF at the root
  BF_val <- bin_BF(root, n, alpha, beta, location, scale, model, hypothesis)

  if (BF_val <= D) {
    # Try root - 1 only if root > 0
    if (root > 0) {
      BF_prev <- bin_BF(root - 1, n, alpha, beta, location, scale, model, hypothesis)
      if (BF_prev > D) return(root - 1)
    }

    # Try root + 1
    BF_next <- bin_BF(root + 1, n, alpha, beta, location, scale, model, hypothesis)
    if (BF_next > D) return(root + 1)
  }

  # Return original if already valid or no better nearby found
  return(root)
}


adjust_root_01 <- function(root, n, alpha, beta, location, scale, model, hypothesis, D) {
  # Evaluate BF at the root
  BF_val <- 1/bin_BF(root, n, alpha, beta, location, scale, model, hypothesis)

  if (BF_val <= D) {
    # Try root - 1
    BF_prev <- 1/bin_BF(root - 1, n, alpha, beta, location, scale, model, hypothesis)
    if (BF_prev > D) return(root - 1)

    # Try root + 1
    BF_next <- 1/bin_BF(root + 1, n, alpha, beta, location, scale, model, hypothesis)
    if (BF_next > D) return(root + 1)
  }

  # Return original if already valid or no better nearby found
  return(root)
}


bin_prior <-function(prop,alpha,beta,location,scale,model){

  switch(model,
         "beta" = stats::dbeta(prop, alpha,beta),
         "Moment" = dMoment(prop,location,scale))
}
bin_BF<-function(x,n,alpha,beta,location,scale,model,hypothesis){
  BF = NA
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = 0, b = location),
                   "!=" = c(a = 0, b = 1)
  )


  normalization <- if (hypothesis == "!=") {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
    }
  for( i in 1:length(x)){
    int  <- function(prop){stats::dbinom(x[i], size=n, prob=prop) *bin_prior(prop,alpha,beta,location,scale,model)}
    lh1 <- stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = 1e-5)$value / normalization
    lh0 <- stats::dbinom(x[i], size = n, prob = location)
    BF[i] = lh1 / lh0
  }


  return(BF)

}

bin_BF_bound_10 <-function(D,n,alpha,beta,location,scale,model,hypothesis){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)- D
  }

  x <- tryCatch(stats::uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)
  #results <- tryCatch(uniroot.all(Bound_finding, lower = 0 ,upper = n), error = function(e) NA)
  results <- round(results[!is.na(results) & is.finite(results)])

  if (length(results) == 0) return("bound cannot be found")


  results <- sapply(results, function(root) {
    adjust_root_10(root, n, alpha, beta, location, scale, model, hypothesis, D)
  })


  BF.vals  <- bin_BF(results,n,alpha,beta,location,scale,model,hypothesis)

  BF.close <- which(BF.vals > D)
  if (length(BF.close) == 0 || all(!is.finite(BF.close))) return("bound cannot be found")
  return(results[BF.close])
}

bin_BF_bound_01 <-function(D,n,alpha,beta,location,scale,model,hypothesis){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    1/bin_BF(x,n,alpha,beta,location,scale,model,hypothesis)- D
  }

  x <- tryCatch(stats::uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)

  results <- round(results[!is.na(results)])
  if (length(results) == 0) return("bound cannot be found")


  results <- sapply(results, function(root) {
    adjust_root_01(root, n, alpha, beta, location, scale, model, hypothesis, D)
  })


  BF.vals  <- 1/bin_BF(results,n,alpha,beta,location,scale,model,hypothesis)

  BF.close <- which(BF.vals > D)
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}


bin_TPE<-function(x,n,h0,alpha,beta,location,scale,model,hypothesis){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  if (model =="Point"){
    TPE = switch(hypothesis,
               "!=" = {

                 switch(length(x)==2,
                        "1" ={stats::pbinom(min(x),n,location,lower.tail = T)+ stats::pbinom(max(x)-1,n,location,lower.tail = F)},
                        "0"=  {
                          switch(x/n>location,
                                 "1" = stats::pbinom(x-1,n,location,lower.tail = F),
                                 "0" = stats::pbinom(x,n,location,lower.tail = T))

                        })
                 },
               ">"  = {stats::pbinom(x-1,n,location,lower.tail = F)},
               "<"  = {stats::pbinom(x,n,location,lower.tail = T)}
    )
    return(TPE)
  }

  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = 0, b = h0),
                   "!=" = c(a = 0, b = 1)
  )
  normalization <- if (hypothesis == "!=") {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
  }
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      stats::pbinom(min(x), n, prop, lower.tail = TRUE) +
                        stats::pbinom(max(x) - 1, n, prop, lower.tail = FALSE)
                    } else {
                      mapply(function(x_i, n_i, p_i) {
                        if (x_i / n_i > location) {
                          stats::pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                        } else {
                          stats::pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                        }
                      }, x, n, prop)
                    }
                  },
                  ">" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "<" = stats::pbinom(x, n, prop, lower.tail = TRUE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalization
  }

  TPE = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-5)$value

  return(TPE)

}

bin_FNE<-function(x,n,h0,alpha,beta,location,scale,model,hypothesis){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  if (model == "Point") {
    FNE <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      stats::pbinom(max(x), n, location, lower.tail = TRUE) - stats::pbinom(min(x) - 1, n, location, lower.tail = TRUE)
                    } else {
                      if ((x / n) > location) {
                        stats::pbinom(x, n, location, lower.tail = TRUE)
                      } else {
                        stats::pbinom(x - 1, n, location, lower.tail = FALSE)
                      }
                    }
                  },
                  ">" = stats::pbinom(x, n, location, lower.tail = TRUE),
                  "<" = stats::pbinom(x - 1, n, location, lower.tail = FALSE)
    )
    return(FNE)
  }



  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = 0, b = h0),
                   "!=" = c(a = 0, b = 1)
  )

  normalization <- if (hypothesis == "!=") {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
  }
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      stats::pbinom(max(x), n, prop, lower.tail = TRUE) - stats::pbinom(min(x) - 1, n, prop, lower.tail = TRUE)
                    } else {
                      if ((x / n) > location) {
                        stats::pbinom(x, n, prop, lower.tail = TRUE)
                      } else {
                        stats::pbinom(x - 1, n, prop, lower.tail = FALSE)
                      }
                    }
                  },
                  ">" = stats::pbinom(x , n, prop, lower.tail = TRUE),
                  "<" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalization
  }
  FNE = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-5)$value
  return(FNE)

}

bin_FPE<-function(x,n,location,hypothesis){

  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  FPE <- switch(hypothesis,
                "!=" = {
                  if (length(x) == 2) {
                    stats::pbinom(min(x), n, location, lower.tail = TRUE) +
                      stats::pbinom(max(x) - 1, n, location, lower.tail = FALSE)
                  } else {
                    mapply(function(x_i, n_i, p_i) {
                      if (x_i / n_i > location) {
                        stats::pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                      } else {
                        stats::pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                      }
                    }, x, n, location)
                  }
                },
                ">" = stats::pbinom(x - 1, n, location, lower.tail = FALSE),
                "<" = stats::pbinom(x, n, location, lower.tail = TRUE)
  )

    return(FPE)

}

bin_TNE<-function(x,n,location,hypothesis){

  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)


  TNE <- switch(hypothesis,
                "!=" = {
                  if (length(x) == 2) {
                    stats::pbinom(max(x), n, location, lower.tail = TRUE) - stats::pbinom(min(x) - 1, n, location, lower.tail = TRUE)
                  } else {
                    if ((x / n) > location) {
                      stats::pbinom(x, n, location, lower.tail = TRUE)
                    } else {
                      stats::pbinom(x - 1, n, location, lower.tail = FALSE)
                    }
                  }
                },
                ">" = stats::pbinom(x, n, location, lower.tail = TRUE),
                "<" = stats::pbinom(x - 1, n, location, lower.tail = FALSE)
  )

  return(TNE)

}

bin_N_finder <-function(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                        alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP){
  lower = 10
  upper = 10000

  b10 = bin_BF_bound_10(D,lower,alpha,beta,location,scale,model,hypothesis)
  TPE_lo <- if (de_an_prior == 1)
    bin_TPE(b10,lower,h0,alpha,beta,location,scale,model,hypothesis) else
      bin_TPE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
  if (TPE_lo > target) return(lower)
  FPE_lo <-  bin_FPE(b10,lower,location,hypothesis)
  if (TPE_lo > target&FPE_lo<FP) return(lower)


  Power_root <- function(N){
    N =round(N)
    b10 = bin_BF_bound_10 (D,N,alpha,beta,location,scale,model,hypothesis)
    pro <- if (de_an_prior==0){
      bin_TPE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
      }else bin_TPE(b10,N,h0,alpha,beta,location,scale,model,hypothesis)

    pro-target
  }

  N.power = round(stats::uniroot(Power_root,lower = lower,upper = upper)$root)+1

  N.extended = seq(N.power,N.power+20,2)
  Power.extended = unlist(lapply(N.extended, Power_root))

  if (any(Power.extended<0)){
    lower = which(Power.extended < 0)[1]
    N.power = round(stats::uniroot(Power_root,lower = lower,upper = upper)$root)+1

  }



  while(TRUE) {
    b10 <- bin_BF_bound_10(D, N.power, alpha, beta, location, scale, model, hypothesis)
    pro <- if (de_an_prior == 0) {
      bin_TPE(b10, N.power,h0, alpha_d, beta_d, location_d, scale_d, model_d, hypothesis)
    } else {
      bin_TPE(b10, N.power,h0, alpha, beta, location, scale, model, hypothesis)
    }

    if (pro > target) break
    N.power <- N.power + 1
  }


  b10 = bin_BF_bound_10(D,N.power,alpha,beta,location,scale,model,hypothesis)
  FPE = bin_FPE(b10,N.power,location,hypothesis)
  if (FPE <= FP) return(N.power)


  alpha.root <- function(n) {
    n=round(n)
    b10 <- bin_BF_bound_10 (D,n,alpha,beta,location,scale,model,hypothesis)
    bin_FPE(b10,n,location,hypothesis)-FP
  }
  N.alpha = round(stats::uniroot(alpha.root,lower = N.power,upper = upper)$root)
  return(N.alpha)
}

bin_N_01_finder <-function(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                           alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP){
  lower = 10
  upper = 10000

  b10 = bin_BF_bound_01(D,lower,alpha,beta,location,scale,model,hypothesis)
  TNE_lo = bin_TNE(b10,lower,location,hypothesis)
  FNE_lo <- if (de_an_prior == 1)
    bin_FNE(b10,lower,h0,alpha,beta,location,scale,model,hypothesis) else
      bin_FNE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)

  if (TNE_lo > target && FNE_lo < FP) {
    return(lower)
  } else if (TNE_lo > target) {
    FN_root <- function(N){
      N =round(N)
      b10 = bin_BF_bound_01 (D,N,alpha,beta,location,scale,model,hypothesis)
      pro <- if (de_an_prior == 1)
        bin_FNE(b10,N,h0,alpha,beta,location,scale,model,hypothesis) else
          bin_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)

      pro-FP
    }
    return(round(stats::uniroot(FN_root, lower = lower, upper = upper)$root))
  }
  TN_root <- function(N){
    N =round(N)
    b10 = bin_BF_bound_01 (D,N,alpha,beta,location,scale,model,hypothesis)
    pro <-  bin_TNE(b10,N,location,hypothesis)

    pro-target
  }

  N.TN = round(stats::uniroot(TN_root,lower = lower,upper = upper)$root)+1
  N.extended = seq(N.TN,N.TN+20,2)
  Power.extended = unlist(lapply(N.extended, TN_root))

  if (any(Power.extended<0)){
    lower = which(Power.extended < 0)[1]
    N.TN = round(stats::uniroot(TN_root,lower = lower,upper = upper)$root)+1

  }



  while(TRUE) {
    b10 <- bin_BF_bound_01(D, N.TN, alpha, beta, location, scale, model, hypothesis)
    pro <- bin_TNE(b10,N.TN,location,hypothesis)

    if (pro > target) break
    N.TN <- N.TN + 1
  }


  b10 = bin_BF_bound_01(D,N.TN,alpha,beta,location,scale,model,hypothesis)
  FNE = if (de_an_prior == 1)
    bin_FNE(b10,N.TN,h0,alpha,beta,location,scale,model,hypothesis) else
      bin_FNE(b10,N.TN,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
  if (FNE <= FP) return(N.TN)
  FN_root <- function(N){
    N =round(N)
    b10 = bin_BF_bound_01 (D,N,alpha,beta,location,scale,model,hypothesis)
    pro <- if (de_an_prior == 1)
      bin_FNE(b10,N,h0,alpha,beta,location,scale,model,hypothesis) else
        bin_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)

    pro-FP
  }
  N.FN = round(stats::uniroot(FN_root,lower = N.TN,upper = upper)$root)
  return(N.FN)
}

bin_table<-function(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,N, mode_bf,FP,direct){
  if (mode_bf == "0") n = N else n = switch(
    direct,
    "h1" = bin_N_finder(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                                 alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP),
    "h0" = bin_N_01_finder(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                        alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP))


  # b bounds:
  b10 <- bin_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis)
  b01 <-  bin_BF_bound_01(D,n,alpha,beta,location,scale,model,hypothesis)


  # max BF10 possible:
  max_BF <- 1 / bin_BF(round(location*n),n,alpha,beta,location,scale,model,hypothesis)
  BF_D   <- b10

  # FPE and TPE:
  FPE       <- bin_FPE(b10,n,location,hypothesis)
  if (de_an_prior == 1) {
    TPE          <- bin_TPE(b10,n,h0,alpha,beta,location,scale,model,hypothesis)
    TPR_alpha    <- alpha
    TPR_beta     <- beta
    TPR_location <- location
    TPR_scale    <- scale
    TPR_model    <- model

  } else {
    TPE          <- bin_TPE(b10,n,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis)
    TPR_alpha    <- alpha_d
    TPR_beta     <- beta_d
    TPR_location <- location_d
    TPR_scale    <- scale_d
    TPR_model    <- model_d
  }


  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- bin_FNE(b01,n,h0,TPR_alpha,TPR_beta,TPR_location,TPR_scale,TPR_model,hypothesis)
    TNE <- bin_TNE(b01,n,location,hypothesis)
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

bin_bf10<- function(D, n, alpha, beta, location, scale, model, hypothesis) {

  # Sequence of successes
  x <- seq(0, n, by = 3)

  # Compute BF10 and bounds
  BF10 <- bin_BF(x, n, alpha, beta, location, scale, model, hypothesis)
  b.BF10 <- bin_BF_bound_10(D, n, alpha, beta, location, scale, model, hypothesis)
  BF10_at_b <- bin_BF(b.BF10, n, alpha, beta, location, scale, model, hypothesis)

  BF01 <- 1 / BF10
  b.BF01 <- bin_BF_bound_01(D, n, alpha, beta, location, scale, model, hypothesis)
  BF01_at_b <- 1 / bin_BF(b.BF01, n, alpha, beta, location, scale, model, hypothesis)

  # Check if BF01 = D is impossible
  max.BF01 <- 1 / bin_BF(round(location * n), n, alpha, beta, location, scale, model, hypothesis)
  impossible <- (hypothesis == "!=") && (max.BF01 < D || identical(b.BF01, "bound cannot be found"))
  # Titles for BF10
  main.bf10 <- if (length(b.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(round(BF10_at_b, 2)) ~ " when x = " ~ .(round(b.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(round(BF10_at_b[1], 2)) ~ "/" ~ .(round(BF10_at_b[2], 2)) ~
                  " when x = " ~ .(round(b.BF10[1], 2)) ~ " or " ~ .(round(b.BF10[2], 2))))
  }

  # Titles for BF01
  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(D)))
  } else if (length(b.BF01) == 1) {
    bquote(bold("BF"[0][1] ~ "=" ~ .(round(BF01_at_b, 2)) ~ " when x = " ~ .(round(b.BF01, 2))))
  } else {
    bquote(bold("BF"[0][1] ~ "=" ~ .(round(BF01_at_b[1], 2)) ~ "/" ~ .(round(BF01_at_b[2], 2)) ~
                  " when x = " ~ .(round(b.BF01[1], 2)) ~ " or " ~ .(round(b.BF01[2], 2))))
  }


  # Data frames for ggplot
  df_bf10 <- data.frame(x = x, BF = BF10)
  df_bf01 <- data.frame(x = x, BF = BF01)

  # Clean theme
  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text  = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  ## ---------- BF10 ----------
  x_breaks_10 <- sort(unique(c(0, n, round(b.BF10, 2))))

  p1 <- ggplot2::ggplot(df_bf10, ggplot2::aes(x = x, y = BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = b.BF10, linetype = "dashed") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(0, n), breaks = x_breaks_10) +
    ggplot2::labs(
      x = "Number of successes",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    clean_theme

  ## ---------- BF01 ----------
  x_breaks_01 <- if (impossible) c(0, n)
  else sort(unique(c(0, n, round(b.BF01, 2))))

  p2 <- ggplot2::ggplot(df_bf01, ggplot2::aes(x = x, y = BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(
      xintercept = if (!impossible) b.BF01 else NA,
      linetype = "dashed"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(0, n), breaks = x_breaks_01) +
    ggplot2::labs(
      x = "Number of successes",
      y = expression("BF"[0][1] * " (log scale)"),
      title = main.bf01
    ) +
    clean_theme

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

Power_bin <- function(D, h0, alpha, beta, location, scale, model, hypothesis,
                         alpha_d, beta_d, location_d, scale_d, model_d,
                         de_an_prior, N) {

  # Sample size range
  Ns <- ceiling(seq(10, N*1.2, length.out = 31))

  # Initialize probability vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(Ns))

  # Compute bounds and probabilities
  for (i in seq_along(Ns)) {
    x10 <- bin_BF_bound_10(D, Ns[i], alpha, beta, location, scale, model, hypothesis)
    x01 <- bin_BF_bound_01(D, Ns[i], alpha, beta, location, scale, model, hypothesis)

    TPE[i] <- if (de_an_prior == 1) {
      bin_TPE(x10, Ns[i], h0, alpha, beta, location, scale, model, hypothesis)
    } else {
      bin_TPE(x10, Ns[i], h0, alpha_d, beta_d, location_d, scale_d, model_d, hypothesis)
    }

    FNE[i] <- if (de_an_prior == 1) {
      bin_FNE(x01, Ns[i], h0, alpha, beta, location, scale, model, hypothesis)
    } else {
      bin_FNE(x01, Ns[i], h0, alpha_d, beta_d, location_d, scale_d, model_d, hypothesis)
    }

    FPE[i] <- bin_FPE(x10, Ns[i], location, hypothesis)
    TNE[i] <- bin_TNE(x01, Ns[i], location, hypothesis)
  }

  # Prepare data for ggplot
  df_BF10 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = Ns,
      `True Positive` = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df_BF10$Type <- factor(df_BF10$Type, levels = c("True Positive", "False Positive"))

  df_BF01 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = Ns,
      `True Negative` = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df_BF01$Type <- factor(df_BF01$Type, levels = c("True Negative", "False Negative"))

  # Colors for lines
  type_colors <- c(
    "True Positive" = "black",
    "False Positive" = "grey50",
    "True Negative" = "black",
    "False Negative" = "grey50"
  )

  # ---------- Theme for axes, text, and grid ----------
  axis_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),      # remove background grid
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # ---------- Theme for legend ----------
  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),           # inside top-left corner
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  # ---------- BF10 Plot ----------
  p1 <- ggplot2::ggplot(df_BF10, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    axis_theme +
    legend_theme

  # ---------- BF01 Plot ----------
  p2 <- ggplot2::ggplot(df_BF01, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    axis_theme +
    legend_theme

  # Combine side-by-side
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}


compute.prior.density.b <- function(prop,alpha,beta,location,scale,model,hypothesis) {
  if (model == "Point") return(rep(NA, length(prop)))
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = 0, b = location),
                   "!=" = c(a = 0, b = 1)
  )
  normalization <- if (hypothesis == "!=") {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(model,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
  }
  bin_prior(prop,alpha,beta,location,scale,model)/ normalization
}


bin_prior_plot <-function(h0,alpha,beta,location,scale,model,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,de_an_prior){
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 1))
  bound          <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = 0, b = h0),
                   "!=" = c(a = 0, b = 1)
  )
  prop           <- seq(bound[1],bound[2],.01)
  normalization  <- stats::integrate(function(prop)bin_prior(prop,alpha,beta,location,scale,model),lower = bound[1],upper = bound[2],rel.tol = 1e-10)$value
  prior.analysis <- compute.prior.density.b(prop,alpha,beta,location,scale,model,hypothesis)
  prior.design   <- if (de_an_prior == 0 && model_d != "Point")
    compute.prior.density.b(prop,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis) else
      rep(NA, length(prop))

  # Combine all values into one vector
  all_vals <- c(prior.analysis, prior.design)

  # Filter out NA and infinite values
  finite_vals <- all_vals[is.finite(all_vals)]

  # Get the max from finite values only
  ylim.max <- max(finite_vals)
  # Base plot:
  plot(prop, prior.analysis, type = "l", lwd = 2,
       xlab = bquote(atop(italic(theta))),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~italic(theta)~"under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))

  # If design prior != analysis prior:
  if (de_an_prior == 0) {
    if (model_d == "Point")
      graphics::arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim.max, length = 0.1, col = "black", lty = 2) else
        graphics::lines(prop, prior.design, lty = 2)

    # Add legend:
    graphics::legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }
}





# ---- binomial_e.r ----
adjust_root_10_e <- function(root, n, alpha, beta, location, scale, model, hypothesis, D,e) {
  # Evaluate BF at the root
  BF_val <- bin_e_BF(root,n,alpha,beta,location,scale,model,hypothesis,e)

  if (BF_val <= D) {
    # Try root - 1
    BF_prev <- bin_e_BF(root-1,n,alpha,beta,location,scale,model,hypothesis,e)
    if (BF_prev > D) return(root - 1)

    # Try root + 1
    BF_next <- bin_e_BF(root+1,n,alpha,beta,location,scale,model,hypothesis,e)
    if (BF_next > D) return(root + 1)
  }

  # Return original if already valid or no better nearby found
  return(root)
}

adjust_root_01_e <- function(root, n, alpha, beta, location, scale, model, hypothesis, D,e) {
  # Evaluate BF at the root
  BF_val <- 1/bin_e_BF(root,n,alpha,beta,location,scale,model,hypothesis,e)

  if (BF_val <= D) {
    # Try root - 1
    BF_prev <- 1/bin_e_BF(root-1,n,alpha,beta,location,scale,model,hypothesis,e)
    if (!is.nan(BF_prev) && !is.na(BF_prev) && BF_prev > D) return(root - 1)

    # Try root + 1
    BF_next <- 1/bin_e_BF(root+1,n,alpha,beta,location,scale,model,hypothesis,e)
    if (!is.nan(BF_next) && !is.na(BF_next) &&BF_next > D) return(root + 1)
  }

  # Return original if already valid or no better nearby found
  return(root)
}

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

  normalizationh1 <- switch(hypothesis,
                            "!=" = {
                              if (model == "beta") {
                                1 - (stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta))
                              } else if (model == "Moment") {
                                (pmom(1 - location, tau = scale^2) - pmom(bound_h1[2] - location, tau = scale^2)) +
                                  (pmom(bound_h1[1] - location, tau = scale^2) - pmom(0 - location, tau = scale^2))
                              }
                            },
                            "<" = ,
                            ">" = {
                              if (model == "beta") {
                                stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta)
                              } else if (model == "Moment") {
                                pmom(bound_h1[2] - location, tau = scale^2) - pmom(bound_h1[1] - location, tau = scale^2)
                              }
                            }
  )

  normalizationh0 <- switch(model,
                            "beta"      =   stats::pbeta(bound_h0[2], alpha, beta) - stats::pbeta(bound_h0[1], alpha, beta),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )

    for (i in 1:length(x)){
  int  <- function(prop){stats::dbinom(x[i], size=n, prob=prop) *bin_prior(prop,alpha,beta,location,scale,model)
  }

  if (hypothesis == "!="){
    lh1 = stats::integrate(int,lower = 0,upper = bound_h1[1], rel.tol=1e-5,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-5,stop.on.error = F)$value
  }else{
    lh1 = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-5,stop.on.error = F)$value

  }
  lh0 = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value


    BF[i] = (lh1/normalizationh1)/(lh0/normalizationh0)
  }
  return(BF)

}


bin_e_BF_bound_10 <-function(D,n,alpha,beta,location,scale,model,hypothesis,e){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)- D
  }
  x <- tryCatch(stats::uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)

  results <- round(results[!is.na(results)])
  if (length(results) == 0) return("bound cannot be found")


  results <- sapply(results, function(root) {
    adjust_root_10_e(root, n, alpha, beta, location, scale, model, hypothesis, D,e)
  })

  BF.vals  <- bin_e_BF(results,n,alpha,beta,location,scale,model,hypothesis,e)
  BF.close <- which(BF.vals > D)
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

bin_e_BF_bound_01 <-function(D,n,alpha,beta,location,scale,model,hypothesis,e){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    1/bin_e_BF(x,n,alpha,beta,location,scale,model,hypothesis,e)- D
  }

  x <- tryCatch(stats::uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)

  results <- round(results[!is.na(results)])
  if (length(results) == 0) return("bound cannot be found")


  results <- sapply(results, function(root) {
    adjust_root_01_e(root, n, alpha, beta, location, scale, model, hypothesis, D,e)
  })

  BF.vals  <- 1/bin_e_BF(results,n,alpha,beta,location,scale,model,hypothesis,e)
  BF.close <- which(BF.vals > D)
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}


bin_e_TPE<-function(x,n,h0,alpha,beta,location,scale,model,hypothesis,e){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)


  if (model =="Point"){
    TPE = switch(hypothesis,
               "!=" = {

                 switch(length(x)==2,
                        "1" ={stats::pbinom(min(x),n,location,lower.tail = T)+ stats::pbinom(max(x)-1,n,location,lower.tail = F)},
                        "0"=  {
                          switch(x/n>location,
                                 "1" = stats::pbinom(x-1,n,location,lower.tail = F),
                                 "0" = stats::pbinom(x,n,location,lower.tail = T))

                        })
                 },
               ">"  = {stats::pbinom(x-1,n,location,lower.tail = F)},
               "<"  = {stats::pbinom(x,n,location,lower.tail = T)}
    )
    return(TPE)
  }

  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = 0, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh1 <- switch(hypothesis,
                            "!=" = {
                              if (model == "beta") {
                                1 - (stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta))
                              } else if (model == "Moment") {
                                (pmom(1 - location, tau = scale^2) - pmom(bound_h1[2] - location, tau = scale^2)) +
                                  (pmom(bound_h1[1] - location, tau = scale^2) - pmom(0 - location, tau = scale^2))
                              }
                            },
                            "<" = ,
                            ">" = {
                              if (model == "beta") {
                                stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta)
                              } else if (model == "Moment") {
                                pmom(bound_h1[2] - location, tau = scale^2) - pmom(bound_h1[1] - location, tau = scale^2)
                              }
                            }
  )
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      stats::pbinom(min(x), n, prop, lower.tail = TRUE) +
                        stats::pbinom(max(x) - 1, n, prop, lower.tail = FALSE)
                    } else {
                      mapply(function(x_i, n_i, p_i) {
                        if (x_i / n_i > location) {
                          stats::pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                        } else {
                          stats::pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                        }
                      }, x, n, prop)
                    }
                  },
                  ">" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "<" = stats::pbinom(x, n, prop, lower.tail = TRUE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalizationh1
  }

  if(hypothesis == "!="){
    TPE = stats::integrate(int,lower = 0,upper = bound_h1[1], rel.tol = 1e-5)$value + stats::integrate(int,lower = bound_h1[2],upper = 1, rel.tol = 1e-5)$value
  }else{
    TPE = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol = 1e-5)$value
  }
  return(TPE)

}


bin_e_FNE<-function(x,n,h0,alpha,beta,location,scale,model,hypothesis,e){

  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)


  if (model =="Point"){
    FNE = switch(hypothesis,
                 "!=" = {

                   switch(length(x)==2,
                          "1" ={stats::pbinom(max(x),n,location,lower.tail = T)- stats::pbinom(min(x)-1,n,location,lower.tail = T)},
                          "0"=  {
                            switch(x/n>location,
                                   "1" = stats::pbinom(x,n,location,lower.tail = T),
                                   "0" = stats::pbinom(x-1,n,location,lower.tail = F))

                          })},
                 ">"  = {stats::pbinom(x,n,location,lower.tail = T)},
                 "<"  = {stats::pbinom(x-1,n,location,lower.tail = F)}
    )
    return(FNE)
  }


  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = 0, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh1 <- switch(hypothesis,
                            "!=" = {
                              if (model == "beta") {
                                1 - (stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta))
                              } else if (model == "Moment") {
                                (pmom(1 - location, tau = scale^2) - pmom(bound_h1[2] - location, tau = scale^2)) +
                                  (pmom(bound_h1[1] - location, tau = scale^2) - pmom(0 - location, tau = scale^2))
                              }
                            },
                            "<" = ,
                            ">" = {
                              if (model == "beta") {
                                stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta)
                              } else if (model == "Moment") {
                                pmom(bound_h1[2] - location, tau = scale^2) - pmom(bound_h1[1] - location, tau = scale^2)
                              }
                            }
  )
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      stats::pbinom(max(x), n, prop, lower.tail = TRUE) - stats::pbinom(min(x) - 1, n, prop, lower.tail = TRUE)
                    } else {
                      mapply(function(x_i, n_i, p_i) {
                      if ((x_i / n_i) > location) {
                        stats::pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                      } else {
                        stats::pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                      }
                    }, x, n, prop)
                    }
                  },
                  ">" = stats::pbinom(x , n, prop, lower.tail = TRUE),
                  "<" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalizationh1
  }
  if(hypothesis == "!="){
    FNE = stats::integrate(int,lower = 0,upper = bound_h1[1], rel.tol = 1e-5)$value + stats::integrate(int,lower = bound_h1[2],upper = 1, rel.tol = 1e-5)$value
  }else{
    FNE = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol = 1e-5)$value
  }


  return(FNE)

}

bin_e_FPE<-function(x,n,h0,alpha,beta,location,scale,model,hypothesis,e){

  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)


  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh0 <- switch(model,
                            "beta"      =   stats::pbeta(bound_h0[2], alpha, beta) - stats::pbeta(bound_h0[1], alpha, beta),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )

  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      stats::pbinom(min(x), n, prop, lower.tail = TRUE) +
                        stats::pbinom(max(x) - 1, n, prop, lower.tail = FALSE)
                    } else {
                      mapply(function(x_i, n_i, p_i) {
                        if (x_i / n_i > location) {
                          stats::pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                        } else {
                          stats::pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                        }
                      }, x, n, prop)
                    }
                  },
                  ">" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "<" = stats::pbinom(x, n, prop, lower.tail = TRUE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalizationh0
  }


    FPE = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol = 1e-5)$value
    return(FPE)

}

bin_e_TNE<-function(x,n,h0,alpha,beta,location,scale,model,hypothesis,e){


  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )

  normalizationh0 <- switch(model,
                            "beta"      =   stats::pbeta(bound_h0[2], alpha, beta) - stats::pbeta(bound_h0[1], alpha, beta),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )
  int <- function(prop) {
    pro <- switch(hypothesis,
                  "!=" = {
                    if (length(x) == 2) {
                      stats::pbinom(max(x), n, prop, lower.tail = TRUE) - stats::pbinom(min(x) - 1, n, prop, lower.tail = TRUE)
                    } else {


                      mapply(function(x_i, n_i, p_i) {
                        if ((x_i / n_i) > location) {
                          stats::pbinom(x_i, n_i, p_i, lower.tail = TRUE)
                        } else {
                          stats::pbinom(x_i - 1, n_i, p_i, lower.tail = FALSE)
                        }
                      }, x, n, prop)



                    }
                  },
                  ">" = stats::pbinom(x , n, prop, lower.tail = TRUE),
                  "<" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, model) / normalizationh0
  }

    TNE = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol = 1e-5)$value



  return(TNE)
}

bin_e_N_finder <-function(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                        alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP,e){
  lower = 10
  upper = 10000

  b10 =  bin_e_BF_bound_10(D,lower,alpha,beta,location,scale,model,hypothesis,e)

  TPE_lo <- if (de_an_prior == 1)
    bin_e_TPE(b10,lower,h0,alpha,beta,location,scale,model,hypothesis,e) else
      bin_e_TPE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
  FPE_lo <-  bin_e_FPE(b10,lower,h0,alpha,beta,location,scale,model,hypothesis,e)
  if (TPE_lo > target&FPE_lo<FP) return(lower)

  Power_root <- function(N){
    N =round(N)
    x = bin_e_BF_bound_10(D,N,alpha,beta,location,scale,model,hypothesis,e)

    if(de_an_prior == 1){
      pro = bin_e_TPE(x,N,h0,alpha,beta,location,scale,model,hypothesis,e)
    } else{
      pro = bin_e_TPE(x,N,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
    }
    return(pro-target)
  }

  N.power = round(stats::uniroot(Power_root,lower = lower,upper = upper)$root)+1

  N.extended = seq(N.power,N.power+20,2)
  Power.extended = unlist(lapply(N.extended, Power_root))

  if (any(Power.extended<0)){
    lower = which(Power.extended < 0)[1]
    N.power = round(stats::uniroot(Power_root,lower = lower,upper = upper)$root)+1

  }




    while(TRUE) {
    b10 <- bin_e_BF_bound_10(D,N.power,alpha,beta,location,scale,model,hypothesis,e)
    pro <- if (de_an_prior == 0) {
      bin_e_TPE(b10, N.power, h0,alpha_d, beta_d, location_d, scale_d, model_d, hypothesis,e)
    } else {
      bin_e_TPE(b10,N.power,h0,alpha,beta,location,scale,model,hypothesis,e)
    }

    if (pro > target) break
    N.power <- N.power + 1
  }
  b10 = bin_e_BF_bound_10(D,N.power,alpha,beta,location,scale,model,hypothesis,e)
  FPE =  bin_e_FPE(b10,N.power,h0,alpha,beta,location,scale,model,hypothesis,e)
  if (FPE <= FP) return(N.power)

  alpha.root <- function(n) {
    n=round(n)
    b10 <- bin_e_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis,e)
    bin_e_FPE(b10,n,h0,alpha,beta,location,scale,model,hypothesis,e)-FP
  }
  N.alpha = round(stats::uniroot(alpha.root,lower = N.power,upper = upper)$root)
  return(N.alpha)

  }
bin_e_N_01_finder <-function(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                             alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP,e){
  lower = 10
  upper = 10000

  b10 =  bin_e_BF_bound_01(D,lower,alpha,beta,location,scale,model,hypothesis,e)
  TNE_lo =  bin_e_TPE(b10,lower,h0,alpha,beta,location,scale,model,hypothesis,e)
  FNE_lo <-  if (de_an_prior == 1)
    bin_e_FNE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e) else
      bin_e_FNE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)

  if (TNE_lo > target && FNE_lo < FP) {
    return(lower)
  } else if (TNE_lo > target) {
    FN_root <- function(N){
      N =round(N)
      b10 =  bin_e_BF_bound_01(D,N,alpha,beta,location,scale,model,hypothesis,e)
      pro <- if (de_an_prior == 1)
        bin_e_FNE(b10,N,h0,alpha,beta,location,scale,model,hypothesis,e) else
          bin_e_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)

      pro-FP
    }
    return(round(stats::uniroot(FN_root, lower = lower, upper = upper)$root))
  }

  TN_root <- function(N){
    N =round(N)
    x = bin_e_BF_bound_01(D,N,alpha,beta,location,scale,model,hypothesis,e)

    pro = bin_e_TNE(x,N,h0,alpha,beta,location,scale,model,hypothesis,e)
    return(pro-target)
  }

  N.TN = round(stats::uniroot(TN_root,lower = lower,upper = upper)$root)+1

  N.extended = seq(N.TN,N.TN+20,2)
  Power.extended = unlist(lapply(N.extended, TN_root))

  if (any(Power.extended<0)){
    lower = which(Power.extended < 0)[1]
    N.TN = round(stats::uniroot(TN_root,lower = lower,upper = upper)$root)+1

  }




  while(TRUE) {
    b10 <- bin_e_BF_bound_01(D,N.TN,alpha,beta,location,scale,model,hypothesis,e)
    pro <- bin_e_TNE(b10,N.TN,h0,alpha,beta,location,scale,model,hypothesis,e)

    if (pro > target) break
    N.TN <- N.TN + 1
  }
  b10 = bin_e_BF_bound_01(D,N.TN,alpha,beta,location,scale,model,hypothesis,e)
  FNE =  if (de_an_prior == 1)
    bin_e_FNE(b10,N.TN,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e) else
      bin_e_FNE(b10,N.TN,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
  if (FNE <= FP) return(N.TN)

  FN_root <- function(N){
    N =round(N)
    b10 =  bin_e_BF_bound_01(D,N,alpha,beta,location,scale,model,hypothesis,e)
    pro <- if (de_an_prior == 1)
      bin_e_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e) else
        bin_e_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)

    pro-FP
  }
  N.FN = round(stats::uniroot(FN_root,lower = N.TN,upper = upper)$root)
  return(N.FN)

}

bin_e_table<-function(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                    alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,N, mode_bf,FP,e,direct){
  if (mode_bf == "0") n = N else n = switch(
    direct,
    "h1" = bin_e_N_finder(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                                                     alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP,e),
    "h0" = bin_e_N_01_finder(D,target,h0,alpha,beta,location,scale,model,hypothesis,
                          alpha_d,beta_d,location_d,scale_d,model_d,de_an_prior,FP,e))

  # b bounds:
  b10 <- bin_e_BF_bound_10(D,n,alpha,beta,location,scale,model,hypothesis,e)
  b01 <-  bin_e_BF_bound_01(D,n,alpha,beta,location,scale,model,hypothesis,e)

  max_BF <- 1 /bin_e_BF(round(location*n),n,alpha,beta,location,scale,model,hypothesis,e)

  # FPE and TPE:
  FPE       <- bin_e_FPE(b10,n,h0,alpha,beta,location,scale,model,hypothesis,e)
  if (de_an_prior == 1) {
    TPE          <- bin_e_TPE(b10,n,h0,alpha,beta,location,scale,model,hypothesis,e)
    TPR_alpha    <- alpha
    TPR_beta     <- beta
    TPR_location <- location
    TPR_scale    <- scale
    TPR_model    <- model

  } else {
    TPE          <- bin_e_TPE(b10,n,h0,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
    TPR_alpha    <- alpha_d
    TPR_beta     <- beta_d
    TPR_location <- location_d
    TPR_scale    <- scale_d
    TPR_model    <- model_d
  }
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | b01 == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- bin_e_FNE(b01,n,h0,TPR_alpha,TPR_beta,TPR_location,TPR_scale,TPR_model,hypothesis,e)
    TNE <- bin_e_TNE(b01,n,h0,alpha,beta,location,scale,model,hypothesis,e)
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

bin_e_bf10 <- function(D, n, alpha, beta, location, scale, model, hypothesis, e) {

  # Sequence of successes
  x <- seq(0, n, by = 3)

  # Compute BF10 and bounds
  BF10 <- bin_e_BF(x, n, alpha, beta, location, scale, model, hypothesis, e)
  b.BF10 <- bin_e_BF_bound_10(D, n, alpha, beta, location, scale, model, hypothesis, e)
  BF10_at_b <- bin_e_BF(b.BF10, n, alpha, beta, location, scale, model, hypothesis, e)

  # Compute BF01 and bounds
  BF01 <- 1 / BF10
  b.BF01 <- bin_e_BF_bound_01(D, n, alpha, beta, location, scale, model, hypothesis, e)
  BF01_at_b <- 1 / bin_e_BF(b.BF01, n, alpha, beta, location, scale, model, hypothesis, e)

  # Check if BF01 = D is impossible
  max.BF01 <- 1 / bin_e_BF(round(n / 2), n, alpha, beta, location, scale, model, hypothesis, e)
  impossible <- (hypothesis == "!=") && (max.BF01 < D || identical(b.BF01, "bound cannot be found"))

  # Titles for BF10
  main.bf10 <- if (length(b.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(round(BF10_at_b, 2)) ~ " when x = " ~ .(round(b.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(round(BF10_at_b[1], 2)) ~ "/" ~ .(round(BF10_at_b[2], 2)) ~
                  " when x = " ~ .(round(b.BF10[1], 2)) ~ " or " ~ .(round(b.BF10[2], 2))))
  }

  # Titles for BF01
  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(D)))
  } else if (length(b.BF01) == 1) {
    bquote(bold("BF"[01] ~ "=" ~ .(round(BF01_at_b, 2)) ~ " when x = " ~ .(round(b.BF01, 2))))
  } else {
    bquote(bold("BF"[01] ~ "=" ~ .(round(BF01_at_b[1], 2)) ~ "/" ~ .(round(BF01_at_b[2], 2)) ~
                  " when x = " ~ .(round(b.BF01[1], 2)) ~ " or " ~ .(round(b.BF01[2], 2))))
  }


  # Data frames for ggplot
  df_bf10 <- data.frame(x = x, BF = BF10)
  df_bf01 <- data.frame(x = x, BF = BF01)

  # Clean theme
  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text  = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  ## ---------- BF10 ----------
  x_breaks_10 <- sort(unique(c(0, n, round(b.BF10, 2))))

  p1 <- ggplot2::ggplot(df_bf10, ggplot2::aes(x = x, y = BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = b.BF10, linetype = "dashed") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(0, n), breaks = x_breaks_10) +
    ggplot2::labs(
      x = "Number of successes",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    clean_theme

  ## ---------- BF01 ----------
  x_breaks_01 <- if (impossible) c(0, n)
  else sort(unique(c(0, n, round(b.BF01, 2))))

  p2 <- ggplot2::ggplot(df_bf01, ggplot2::aes(x = x, y = BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(
      xintercept = if (!impossible) b.BF01 else NA,
      linetype = "dashed"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(0, n), breaks = x_breaks_01) +
    ggplot2::labs(
      x = "Number of successes",
      y = expression("BF"[0][1] * " (log scale)"),
      title = main.bf01
    ) +
    clean_theme

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

Power_e_bin <- function(D, h0, alpha, beta, location, scale, model, hypothesis,
                        alpha_d, beta_d, location_d, scale_d, model_d, de_an_prior, N, e) {

  # Sample size range
  smin <- 10
  smax <- N * 1.2
  sN <- ceiling(seq(smin, smax, length.out = 51))  # 51 points for smooth curves

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(sN))

  for (i in seq_along(sN)) {

    x10 <- bin_e_BF_bound_10(D, sN[i], alpha, beta, location, scale, model, hypothesis, e)
    x01 <- bin_e_BF_bound_01(D, sN[i], alpha, beta, location, scale, model, hypothesis, e)

    # True Positive
    TPE[i] <- if (de_an_prior == 1) {
      bin_e_TPE(x10, sN[i], h0, alpha, beta, location, scale, model, hypothesis, e)
    } else {
      bin_e_TPE(x10, sN[i], h0, alpha_d, beta_d, location_d, scale_d, model_d, hypothesis, e)
    }

    # False Negative
    FNE[i] <- if (de_an_prior == 1) {
      bin_e_FNE(x01, sN[i], h0, alpha, beta, location, scale, model, hypothesis, e)
    } else {
      bin_e_FNE(x01, sN[i], h0, alpha_d, beta_d, location_d, scale_d, model_d, hypothesis, e)
    }

    # False Positive & True Negative
    FPE[i] <- bin_e_FPE(x10, sN[i], h0, alpha, beta, location, scale, model, hypothesis, e)
    TNE[i] <- bin_e_TNE(x01, sN[i], h0, alpha, beta, location, scale, model, hypothesis, e)
  }

  # Prepare data for ggplot
  df_bf10 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = sN,
      `True Positive` = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df_bf10$Type <- factor(df_bf10$Type, levels = c("True Positive", "False Positive"))

  df_bf01 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = sN,
      `True Negative` = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df_bf01$Type <- factor(df_bf01$Type, levels = c("True Negative", "False Negative"))

  # Colors
  type_colors <- c(
    "True Positive" = "black",
    "False Positive" = "grey50",
    "True Negative" = "black",
    "False Negative" = "grey50"
  )

  # Clean theme
  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Legend theme
  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  # BF10 plot
  p1 <- ggplot2::ggplot(df_bf10, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # BF01 plot
  p2 <- ggplot2::ggplot(df_bf01, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

compute.prior.density.be.h1 <- function(h0,prop,alpha,beta,location,scale,model,hypothesis,e) {
  if (model == "Point") return(rep(NA, length(prop)))
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = 0, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )

  prior_h1<- bin_prior(prop,alpha,beta,location,scale,model)
  switch(hypothesis,
         "!=" = { prior_h1[prop>min(bound_h1)&prop<max(bound_h1)]=0 },
         ">" = { prior_h1[prop<bound_h1[1]]=0 },
         "<" = { prior_h1[prop>bound_h1[2]]=0 }
  )
  prior_h1
}


compute.prior.density.be.h0 <- function(h0,prop,alpha,beta,location,scale,model,hypothesis,e) {
  if (model == "Point") return(rep(NA, length(prop)))
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )

  prior_h0<- bin_prior(prop,alpha,beta,location,scale,model)
  switch(hypothesis,
         "!=" = { prior_h0[prop<min(bound_h0)|prop>max(bound_h0)]=0 },
         ">" = { prior_h0[prop>bound_h0[2]]=0 },
         "<" = { prior_h0[prop<bound_h0[1]]=0 }
  )
  prior_h0
}





bin_e_prior_plot <-function(h0,alpha,beta,location,scale,model,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,de_an_prior,e){
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 1))
  bounds <- switch(hypothesis,
                        ">"  = c(h0, 1),
                        "<"  = c(0, h0),
                        "!=" = c(0, 1))
  prop           <- seq(bounds[1],bounds[2],.002)


  prior.analysis.h1 <- compute.prior.density.be.h1(h0,prop,alpha,beta,location,scale,model,hypothesis,e)
  prior.analysis.h0<-  compute.prior.density.be.h0(h0,prop,alpha,beta,location,scale,model,hypothesis,e)
  prior.design <- if (de_an_prior == 0 && model_d != "Point") {
    compute.prior.density.be.h1(h0,prop,alpha_d,beta_d,location_d,scale_d,model_d,hypothesis,e)
  } else {
    rep(NA, length(prop))
  }
  # Combine all values into one vector
  all_vals <- c(prior.analysis.h1, prior.analysis.h0, prior.design)

  # Filter out NA and infinite values
  finite_vals <- all_vals[is.finite(all_vals)]

  # Get the max from finite values only
  ylim.max <- max(finite_vals)


  # Base plot
  plot(prop, prior.analysis.h1, type = "l", lwd = 2,
       xlab =bquote(atop(italic(theta))),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~italic(theta)~"under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))

  graphics::lines(prop, prior.analysis.h0, lty = 2, col = "black", lwd = 2)

  # Optional: design prior
  legend.labels <- c("H1 - Analysis Prior", "H0 - Analysis Prior")
  legend.cols   <- c("black", "black")
  legend.lty    <- c(1, 2)
  legend.lwd    <- c(2, 2)

  if (de_an_prior == 0) {
    if (model_d == "Point") {
      graphics::arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim.max,
             length = 0.1, col = "gray", lty = 2)
    } else {
      graphics::lines(prop, prior.design, lty = 1, col = "gray", lwd = 3)
    }

    # Add design prior to legend
    legend.labels <- c(legend.labels, "Design prior")
    legend.cols   <- c(legend.cols, "gray")
    legend.lty    <- c(legend.lty, 1)
    legend.lwd    <- c(legend.lwd, 2)
  }

  graphics::legend("topleft",
         legend = legend.labels,
         col = legend.cols,
         lty = legend.lty,
         lwd = legend.lwd,
         bty = "n")

}



# ---- Correlation.r ----

#Fisher
r_mean <-function(r){
  as.numeric(r)
  (1/2)*log((1+r)/(1-r))
}

r_sd <-function(N){
  1/sqrt(N-3)
}
#prior
d_strechted_beta <-function(rho,k,a,b){
  alpha = beta=1/k
  d_beta(rho, alpha, beta,-1,1)
  #2^((k-2)/k)*(1-rho^2)^((1-k)/k)/beta(1/k,1/k)

}

p_beta <-function(rho, alpha, beta,a,b){
  ExtDist::pBeta_ab(
    rho,
    shape1 = alpha,
    shape2 = beta,
    a = -1,
    b = 1
  )
}


d_beta <- function(rho, alpha, beta,a,b) {

  # Beta function
  B_ab <- beta(alpha, beta)
  #
  a=-1
  b=1
  # Compute the PDF
  pdf_value <- ((rho - a)^(alpha - 1) * (b - rho)^(beta - 1)) / ((b - a)^(alpha + beta - 1) * B_ab)

  return(pdf_value)
}


# likelihood of non-local prior
dMoment <-function(delta,mu,ta){
  ((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
}

r_prior<- function(rho,k,location,scale,dff,model, alpha, beta,a,b){

  switch(model,
         "Normal" = stats::dnorm(rho,location,scale),
         "d_beta"   = d_strechted_beta(rho,k,a,b),
          "Moment"   = dMoment(rho,location,scale),
          "t_dis" = tstude(rho,location,scale,dff),
         "beta" = d_beta(rho, alpha, beta,a,b))
}


d_cor <- function(r, rho, n) {
  n=n-1

  # Calculate the logarithmic terms
  log_gamma_n <- lgamma(n)
  log_gamma_n_plus_half <- lgamma(n + 0.5)

  # Calculate the logarithmic difference
  log_difference <- log_gamma_n - log_gamma_n_plus_half

  # Exponentiate to get the ratio
  ratio <- exp(log_difference)

  # Logarithmic version of the rest of the terms
  log_likelihood_value <- log(n - 1) - 0.5 * log(2 * pi) + log(ratio) +
    0.5 * n * log(1 - rho^2) +
    0.5 * (n - 3) * log(1 - r^2) +
    (-n + 0.5) * log(1 - rho * r)  # This term might go to infinity

  # Exponentiate the result to return it in original scale
  likelihood_value <- exp(log_likelihood_value) *
    gsl::hyperg_2F1(0.5, 0.5, n + 0.5, 0.5 * (r * rho + 1))

  return(likelihood_value)
}


r_BF10<-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  x = NA
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  normalization <- if (hypothesis == "!=") {
    switch(model,
           "d_beta"   = 1,
           "beta" = 1,
           "Moment"   = { pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})

  }else{
    switch(model,
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,-1,1)-p_beta(bound[1], 1/k,1/k,-1,1) ,
           "beta" = p_beta(bound[2], alpha, beta,-1,1)-p_beta(bound[1], alpha, beta,-1,1),
           "Moment"   = {pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})
}

  # Define the integrand function for marginal likelihood under H1
  int <- function(rho, ri) {
    d_cor(ri, rho, n) * r_prior(rho, k, location, scale, dff, model, alpha, beta, min(bound), max(bound))
  }

  # Compute Bayes factors for each observed correlation ri
  x <- sapply(r, function(ri) {
    # Marginal likelihood under H1 (integrated over rho)
    lh1 <- stats::integrate(int, ri = ri, lower = bound[1], upper = bound[2],
                     stop.on.error = FALSE, rel.tol = 1e-4)$value / normalization
    # Likelihood under H0 (fixed rho = h0)
    lh0 <- d_cor(ri, h0, n)
    # Bayes factor
    lh1 / lh0
  })

  return(x)
}

r_BF_bound_10 <-function(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  y <- numeric(0)
  Bound_finding <-function(r)r_BF10(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)- D

  x <- tryCatch(stats::uniroot(Bound_finding, lower = -.99, upper = h0,tol = 1e-5)$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower =  h0, upper = .99,tol = 1e-5)$root, error = function(e) NA)
  results <- c(x, y)
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")

  BF.vals  <- r_BF10(results,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

r_BF_bound_01 <-function(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){
  r_BF_bound_10(1/D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
}

p_cor<-function(limit,rho,n,lower.tail){

  stats::pnorm(r_mean(limit),r_mean(rho),sd = r_sd(n),lower.tail =  lower.tail)


}

r_TPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)

  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {p_cor(max(r),location,n,lower.tail = F)+ p_cor(min(r),location,n,lower.tail = T)},
               ">"  = {p_cor(r,location,n,lower.tail =F)},
               "<"  = {p_cor(r,location,n,lower.tail =T)}
    )
    return(x)
  }

  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  normalization <-   normalization <- if (hypothesis == "!=") {
    switch(model,
           "d_beta"   = 1,
           "beta" = 1,
           "Moment"   = { pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})

  }else{
    switch(model,
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,-1,1)-p_beta(bound[1], 1/k,1/k,-1,1) ,
           "beta" = p_beta(bound[2], alpha, beta,-1,1)-p_beta(bound[1], alpha, beta,-1,1),
           "Moment"   = {pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})
  }
  int <- function(rho) {
    prob <- switch(hypothesis,
                   "!=" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                     p_cor(min(r), rho, n, lower.tail = TRUE),
                   ">"  = p_cor(r, rho, n, lower.tail = FALSE),
                   "<"  = p_cor(r, rho, n, lower.tail = TRUE)
    )

    prob * r_prior(rho, k, location, scale, dff, model, alpha, beta,min(bound),max(bound)) / normalization
  }
  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-4)$value
  return(x)

}

r_FNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)


  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {p_cor(max(r),location,n,lower.tail = T)- p_cor(min(r),location,n,lower.tail = T)},
               ">"  = {p_cor(r,location,n,lower.tail =T)},
               "<"  = {p_cor(r,location,n,lower.tail =F)}
    )
    return(x)
  }


  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )

  normalization <-  normalization <- if (hypothesis == "!=") {
    switch(model,
           "d_beta"   = 1,
           "beta" = 1,
           "Moment"   = { pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})

  }else{
    switch(model,
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,-1,1)-p_beta(bound[1], 1/k,1/k,-1,1) ,
           "beta" = p_beta(bound[2], alpha, beta,-1,1)-p_beta(bound[1], alpha, beta,-1,1),
           "Moment"   = {pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})
  }
  int <- function(rho) {
    prob <- switch(hypothesis,
                   "!=" = p_cor(max(r), rho, n, lower.tail = TRUE) -
                     p_cor(min(r), rho, n, lower.tail = TRUE),
                   ">"  = p_cor(r, rho, n, lower.tail = TRUE),
                   "<"  = p_cor(r, rho, n, lower.tail = FALSE)
    )

    prob * r_prior(rho, k, location, scale, dff, model, alpha, beta,min(bound),max(bound)) / normalization
  }


  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-8, subdivisions=10000000)$value
  return(x)

}

r_FPE <-function(r,n,h0,hypothesis){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)

  x <- switch(hypothesis,
              "!=" = p_cor(max(r), h0, n, lower.tail = FALSE) +
                p_cor(min(r), h0, n, lower.tail = TRUE),
              ">"  = p_cor(r, h0, n, lower.tail = FALSE),
              "<"  = p_cor(r, h0, n, lower.tail = TRUE)
  )
  return(x)

}


r_TNE <-function(r,n,h0,hypothesis){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)

  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )

  x <- switch(hypothesis,
              "!=" = p_cor(max(r), h0, n, lower.tail = TRUE) -
                p_cor(min(r), h0, n, lower.tail = TRUE),
              ">"  = p_cor(r, h0, n, lower.tail = TRUE),
              "<"  = p_cor(r, h0, n, lower.tail = FALSE)
  )

  return(x)

}



r_N_finder<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                       location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP){

  lo = 10
  upper = 5000

  r = r_BF_bound_10(D,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  TPE_lo <- if (de_an_prior == 1)
    r_TPE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model) else
      r_TPE(r,lo,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
  FPE_lo <-  r_FPE(r,lo,h0,hypothesis )

  if (TPE_lo > target&FPE_lo<FP) return(lo)

  Power_root <- function(N) {
    r <- r_BF_bound_10(D, N, k, alpha, beta, h0, hypothesis, location, scale, dff, model)
    pro <- if (de_an_prior==0){ r_TPE(r, N, k_d, alpha_d, beta_d, h0, hypothesis, location_d, scale_d, dff_d, model_d) }else r_TPE(r, N, k, alpha, beta, h0, hypothesis, location, scale, dff, model)

    pro - target
  }

  N.power = stats::uniroot(Power_root,lower = lo,upper = upper)$root

  ## checking if the N lead to an acceptable alpha level
  r = r_BF_bound_10(D,N.power,k, alpha, beta,h0,hypothesis,location,scale,dff,model)

  FPE = r_FPE(r,N.power,h0,hypothesis)
  if (FPE <= FP) return(N.power)

  alpha.root <- function(n) {
    r <- r_BF_bound_10(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
   r_FPE(r,n,h0,hypothesis)-FP
  }
  N.alpha = stats::uniroot(alpha.root,lower = N.power,upper = upper)$root
  return(N.alpha)
  }
r_N_01_finder<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                        location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP){

  lo = 10
  upper = 5000

  r = r_BF_bound_01(D,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  TNE_lo <- r_TNE(r,lo,h0,hypothesis )
  FNE_lo <-  if (de_an_prior == 1)
    r_FNE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model) else
      r_FNE(r,lo,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)

  if (TNE_lo > target && TNE_lo < FP) {
    return(lo)
  } else if (TNE_lo > target) {
    FN.root <- function(n) {
     r <- r_BF_bound_01(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    FNE<-  if (de_an_prior == 1)
        r_FNE(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model) else
          r_FNE(r,n,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d)
    FNE- FP
    }
    return(stats::uniroot(FN.root, lower = lo, upper = upper)$root)
  }

  TN_root <- function(N) {
    r <- r_BF_bound_01(D, N, k, alpha, beta, h0, hypothesis, location, scale, dff, model)
    pro <- r_TNE(r,N,h0,hypothesis )
    pro - target
  }

  N.TN <- tryCatch(
    stats::uniroot(TN_root, lower = lo, upper = upper)$root,
    error = function(e) 20
  )

  ## checking if the N lead to an acceptable alpha level
  r = r_BF_bound_01(D,N.TN,k, alpha, beta,h0,hypothesis,location,scale,dff,model)

  FNE =   if (de_an_prior==0){ r_FNE(r, N.TN, k_d, alpha_d, beta_d, h0, hypothesis, location_d, scale_d, dff_d, model_d) }else r_FNE(r, N.TN, k, alpha, beta, h0, hypothesis, location, scale, dff, model)

  if (FNE <= FP) return(N.TN)

  FN.root <- function(n) {
    r   <- r_BF_bound_10(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
    FNE <- if (de_an_prior==0){ r_FNE(r, n, k_d, alpha_d, beta_d, h0, hypothesis, location_d, scale_d, dff_d, model_d)
      }else r_FNE(r, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model)
    FNE-FP
  }
  N.FN = stats::uniroot(FN.root,lower = N.TN,upper = upper)$root
  return(N.FN)
}


r_table<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                    location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,FP,direct ){

  n <- if (mode_bf == 1) {
    switch(direct,
           "h1" = ceiling(r_N_finder(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                           location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP)),
           "h0" = ceiling(r_N_01_finder(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                                     location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP)))} else  n = N

  # r bounds:
  r10 <- r_BF_bound_10(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  r01 <-  r_BF_bound_01(D,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)

  # max BF10 possible:
  max_BF <- 1 / r_BF10(h0,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model)
  BF_D   <- r10

  # FPE and TPE:
  FPE       <- r_FPE(r10,n,h0,hypothesis)
  if (de_an_prior == 1) {
    TPE         <- r_TPE(r10, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model)
    TPR_model   <- model
    TPR_k       <- k
    TPR_alpha   <- alpha
    TPR_beta    <- beta
    TPR_location<- location
    TPR_scale   <- scale
    TPR_dff     <- dff

  } else {
    TPE       <- r_TPE(r10, n, k_d, alpha_d, beta_d, h0, hypothesis, location_d, scale_d, dff_d, model_d)
    TPR_model   <- model_d
    TPR_k       <- k_d
    TPR_alpha   <- alpha_d
    TPR_beta    <- beta_d
    TPR_location<- location_d
    TPR_scale   <- scale_d
    TPR_dff     <- dff_d
  }
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- r_FNE(r01,n,TPR_k, TPR_alpha, TPR_beta,h0,hypothesis,TPR_location,TPR_scale,TPR_dff,TPR_model)
    TNE <- r_TNE(r01,n,h0,hypothesis)
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

compute.prior.density.r <- function(rho, k,location,scale,dff,model, alpha, beta,hypothesis) {
  if (model == "Point") return(rep(NA, length(rho)))
  bound  <- switch(hypothesis,
                   ">" = c(a = location, b = 1),
                   "<" = c(a = -1, b = location),
                   "!=" = c(a = -1, b = 1)
  )
  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Normal" = stats::pnorm(bound[2],location,scale)-stats::pnorm(bound[1],location,scale),
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,min(bound),max(bound))-p_beta(bound[1], 1/k,1/k,min(bound),max(bound)) ,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "t_dis" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0),
           "beta" = p_beta(bound[2], alpha, beta,min(bound),max(bound))-p_beta(bound[1], alpha, beta,min(bound),max(bound)))


  r_prior(rho,k,location,scale,dff,model, alpha, beta,min(bound),max(bound)) / normalization
}


r_prior_plot <-function(k, alpha, beta,h0,location,scale,dff,model,de_an_prior,
                        k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d,hypothesis){
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 1))
  bound  <- switch(hypothesis,
                   ">" = c(a = h0, b = 1),
                   "<" = c(a = -1, b = h0),
                   "!=" = c(a = -1, b = 1)
  )
  rho <- seq(bound[1],bound[2],.01)
  prior.analysis <-compute.prior.density.r(rho, k,location,scale,dff,model, alpha, beta,hypothesis)
  prior.design   <- if (de_an_prior == 0 && model_d != "Point")
    compute.prior.density.r(rho, k_d,location_d,scale_d,dff_d,model_d, alpha_d, beta_d,hypothesis) else
      rep(NA, length(rho))
  # Combine all values into one vector
  all_vals <- c(prior.analysis,prior.design)

  # Filter out NA and infinite values
  finite_vals <- all_vals[is.finite(all_vals)]

  # Get the max from finite values only
  ylim.max <- max(finite_vals)
  # Base plot:
  plot(rho, prior.analysis, type = "l", lwd = 2,
       xlab = expression(bold(rho)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~rho~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))

  # If design prior != analysis prior:
  if (de_an_prior == 0) {
    if (model_d == "Point")
      graphics::arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim.max, length = 0.1, col = "black", lty = 2) else
        graphics::lines(rho, prior.design, lty = 2)

    # Add legend:
    graphics::legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }

}


r_bf10_p <- function(D, n, k, alpha, beta, h0, hypothesis,
                     location, scale, dff, model) {

  rr <- seq(-0.99, 0.99, 0.01)

  BF10   <- r_BF10(rr, n, k, alpha, beta, h0, hypothesis,
                   location, scale, dff, model)
  r.BF10 <- r_BF_bound_10(D, n, k, alpha, beta, h0, hypothesis,
                          location, scale, dff, model)

  BF01   <- 1 / BF10
  r.BF01 <- r_BF_bound_01(D, n, k, alpha, beta, h0, hypothesis,
                          location, scale, dff, model)

  max.BF01 <- 1 / r_BF10(h0, n, k, alpha, beta, h0, hypothesis,
                         location, scale, dff, model)

  impossible <- (hypothesis == "!=") &&
    (max.BF01 < D || identical(r.BF01, "bound cannot be found"))

  ## ---------- Titles ----------
  main.bf10 <- if (length(r.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(D) ~ " when r = " ~ .(round(r.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(D) ~ " when r = " ~ .(round(r.BF10[1], 2)) ~
                  " or " ~ .(round(r.BF10[2], 2))))
  }

  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(D)))
  } else if (length(r.BF01) == 1) {
    bquote(bold("BF"[0][1] ~ "=" ~ .(D) ~ " when r = " ~ .(round(r.BF01, 2))))
  } else {
    bquote(bold("BF"[0][1] ~ "=" ~ .(D) ~ " when r = " ~ .(round(r.BF01[1], 2)) ~
                  " or " ~ .(round(r.BF01[2], 2))))
  }

  df10 <- data.frame(r = rr, BF = BF10)
  df01 <- data.frame(r = rr, BF = BF01)

  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text  = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  ## ---------- BF10 ----------
  x_breaks_10 <- sort(unique(c(-1, 1, round(r.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(r, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = r.BF10, linetype = "dashed") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = x_breaks_10) +
    ggplot2::labs(
      x = "Correlation",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    clean_theme

  ## ---------- BF01 ----------
  x_breaks_01 <- if (impossible) c(-1, 1)
  else sort(unique(c(-1, 1, round(r.BF01, 2))))

  p2 <- ggplot2::ggplot(df01, ggplot2::aes(r, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +   # ← all black now
    ggplot2::geom_vline(
      xintercept = if (!impossible) r.BF01 else NA,
      linetype = "dashed"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = x_breaks_01) +
    ggplot2::labs(
      x = "Correlation",
      y = expression("BF"[0][1] * " (log scale)"),
      title = main.bf01
    ) +
    clean_theme

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}


Power_r <- function(D, k, alpha, beta, h0, hypothesis,
                    location, scale, dff, model,
                    k_d, alpha_d, beta_d,
                    location_d, scale_d, dff_d, model_d,
                    de_an_prior, N) {

  Ns <- seq(4, ceiling(N * 1.2), length.out = 31)

  TPE <- FPE <- TNE <- FNE <- numeric(length(Ns))

  for (i in seq_along(Ns)) {

    r10 <- r_BF_bound_10(D, Ns[i], k, alpha, beta, h0, hypothesis,
                         location, scale, dff, model)
    r01 <- r_BF_bound_01(D, Ns[i], k, alpha, beta, h0, hypothesis,
                         location, scale, dff, model)

    TPE[i] <- if (de_an_prior == 1)
      r_TPE(r10, Ns[i], k, alpha, beta, h0, hypothesis,
            location, scale, dff, model)
    else
      r_TPE(r10, Ns[i], k_d, alpha_d, beta_d, h0, hypothesis,
            location_d, scale_d, dff_d, model_d)

    FPE[i] <- r_FPE(r10, Ns[i], h0, hypothesis)

    TNE[i] <- r_TNE(r01, Ns[i], h0, hypothesis)

    FNE[i] <- if (de_an_prior == 1)
      r_FNE(r01, Ns[i], k, alpha, beta, h0, hypothesis,
            location, scale, dff, model)
    else
      r_FNE(r01, Ns[i], k_d, alpha_d, beta_d, h0, hypothesis,
            location_d, scale_d, dff_d, model_d)
  }

  ## ---------- Data ----------
  df1 <- tidyr::pivot_longer(
    data.frame(
      SampleSize = Ns,
      `True Positive`  = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df1$Type <- factor(df1$Type, levels = c("True Positive", "False Positive"))

  df2 <- tidyr::pivot_longer(
    data.frame(
      SampleSize = Ns,
      `True Negative`  = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df2$Type <- factor(df2$Type, levels = c("True Negative", "False Negative"))

  ## ---------- Style ----------
  type_colors <- c(
    "True Positive"  = "black",
    "False Positive" = "grey50",
    "True Negative"  = "black",
    "False Negative" = "grey50"
  )

  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x  = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 12),
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  ## ---------- Plots ----------
  p1 <- ggplot2::ggplot(df1,
                        ggplot2::aes(SampleSize, Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  p2 <- ggplot2::ggplot(df2,
                        ggplot2::aes(SampleSize, Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

# ---- Correlation_e.r ----
r_auto_uniroot_fixed_lower <- function(f,lower, upper = 1, step = .05, max_attempts = 25) {
  attempts <- 0
  while (attempts < max_attempts) {
    attempts <- attempts + 1
    # Try to find the root with the current bounds
    result <- tryCatch({
      stats::uniroot(f, lower = lower, upper = upper , tol = 1e-10)$root
    }, error = function(e) {
      # If there's an error, return NA to indicate no root found
      return(NA)
    })

    # If a root is found (not NA), return the result
    if (!is.na(result)) {
      return(result)
    }

    # If no root is found, expand the search range and try again
    upper <- upper - step
  }

}

r_auto_uniroot_fixed_upper <- function(f,upper, lower = -1, step = .05, max_attempts = 25, ...) {
  attempts <- 0
  while (attempts < max_attempts) {
    attempts <- attempts + 1
    # Try to find the root with the current bounds
    result <- tryCatch({
      stats::uniroot(f, lower = lower, upper = upper, tol = 1e-10)$root
    }, error = function(e) {
      # If there's an error, return NA to indicate no root found
      return(NA)
    })

    # If a root is found (not NA), return the result
    if (!is.na(result)) {
      return(result)
    }

    # If no root is found, expand the search range and try again
    lower <- lower + step
  }

}


re_BF10i<-function(r,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
    x = NA

    bound_h1  <- switch(hypothesis,
                        ">" = c(a = h0+e, b = 1),
                        "<" = c(a = -1, b = h0+e),
                        "!=" = c(a = h0+e[1], b = h0+e[2])
    )
    bound_h0  <- switch(hypothesis,
                        ">" = c(a = h0, b = h0+e),
                        "<" = c(a = h0+e, b = h0),
                        "!=" = c(a = h0+e[1], b = h0+e[2])
    )

    normalizationh1 <- switch(hypothesis,
                              "!=" = switch(model,
                                            "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                            "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                            "Moment"          = {
                                              (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                                (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                            }),

                              "<"  = switch(model,
                                            "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                            "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                            "Moment"          = {
                                              (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                            }),
                              ">"  = switch(model,
                                            "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                            "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                            "Moment"          = {
                                              (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                            })
    )
    normalizationh0 <- switch(model,
                              "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                              "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                              "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                              }
    )

    int  <- function(rho){d_cor(r,rho,n)*r_prior(rho,k,location,scale,dff,model, alpha, beta,-1,1)/normalizationh1
    }

    if (hypothesis == "!="){
      lh1 = stats::integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-5,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-5,stop.on.error = F)$value
    }else{
      lh1 = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-5,stop.on.error = F)$value

    }
    lh0 = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value

    x = (lh1/normalizationh1)/(lh0/normalizationh0)
    return(x)
  }
re_BF10<-function(r,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
   sapply(r, re_BF10i, n = n, k = k, alpha = alpha, beta = beta,
              h0 = h0, hypothesis = hypothesis, location = location,
              scale = scale, dff = dff, model = model, e = e)
}



re_BF_bound_10 <-function(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  y <- numeric(0)
  Bound_finding <-function(r){
    re_BF10(r,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)- D
  }
  opt_result <- stats::optimize(Bound_finding, interval = c(-1, 1))$minimum

  if (hypothesis=="!="){
    x <- r_auto_uniroot_fixed_upper (Bound_finding,opt_result, lower = -1, step = .05, max_attempts = 25)
    y <- r_auto_uniroot_fixed_lower(Bound_finding,opt_result, upper = 1, step = .05, max_attempts = 25)
  }
  if (hypothesis == ">"){

    x <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
  }

  if (hypothesis == "<"){

    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)
  }
  results <- c(x, y)
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")

  BF.vals  <- re_BF10(results,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

re_BF_bound_01 <-function(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e){
  re_BF_bound_10 (1/D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
}

re_TPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){

  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }

  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {p_cor(max(r),location,n,lower.tail = F)+ p_cor(min(r),location,n,lower.tail = T)},
               ">"  = {p_cor(r,location,n,lower.tail =F)},
               "<"  = {p_cor(r,location,n,lower.tail =T)}
    )
    return(x)
  }

  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = -1, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh1 <- switch(hypothesis,
                            "!=" = switch(model,
                                          "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "Moment"          = {
                                            (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                              (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                          }),

                            "<"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "Moment"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          }),
                            ">"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "Moment"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          })
  )

  int <- function(rho) {
    pro <- switch(hypothesis,
                  "!=" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  ">"  = p_cor(r, rho, n, lower.tail = FALSE),
                  "<"  = p_cor(r, rho, n, lower.tail = TRUE)
    )

    pro * r_prior(rho, k, location, scale, dff, model, alpha, beta,-1,1) / normalizationh1
  }



  if (hypothesis == "!="){
    x = stats::integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-5,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-5,stop.on.error = F)$value
  }else{
    x = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-5,stop.on.error = F)$value

  }
  return(x)

}

re_FNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){

  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }
  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {p_cor(max(r),location,n,lower.tail = T)- p_cor(min(r),location,n,lower.tail = T)},
               ">"  = {p_cor(r,location,n,lower.tail =T)},
               "<"  = {p_cor(r,location,n,lower.tail =F)}
    )
    return(x)
  }
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = -1, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )

  normalizationh1 <- switch(hypothesis,
                            "!=" = switch(model,
                                          "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "Moment"          = {
                                            (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                              (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                          }),

                            "<"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "Moment"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          }),
                            ">"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "Moment"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          })
  )

  int <- function(rho) {
    pro <- switch(hypothesis,
                  "!=" = p_cor(max(r), rho, n, lower.tail = T) -
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  ">"  = p_cor(r, rho, n, lower.tail = T),
                  "<"  = p_cor(r, rho, n, lower.tail = F)
    )

    pro * r_prior(rho, k, location, scale, dff, model, alpha, beta,-1,1) / normalizationh1
  }


  if (hypothesis == "!="){
    x = stats::integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-10,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-10,stop.on.error = F)$value
  }else{
    x = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-10,stop.on.error = F)$value

  }

  return(x)

}

re_FPE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){

  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh0 <- switch(model,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )


  int <- function(rho) {
    pro <- switch(hypothesis,
                  "!=" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                         p_cor(min(r), rho, n, lower.tail = TRUE),
                  ">"  = p_cor(r, rho, n, lower.tail = FALSE),
                  "<"  = p_cor(r, rho, n, lower.tail = TRUE)
    )

    pro * r_prior(rho, k, location, scale, dff, model, alpha, beta,-1,1) / normalizationh0
  }

  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value


  return(x)

}

re_TNE <-function(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e){

  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )

  normalizationh0 <- switch(model,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )

  int <- function(rho) {
    pro <- switch(hypothesis,
                  "!=" = p_cor(max(r), rho, n, lower.tail = TRUE) -
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  ">"  = p_cor(r, rho, n, lower.tail = TRUE),
                  "<"  = p_cor(r, rho, n, lower.tail = FALSE)
    )

    pro * r_prior(rho, k, location, scale, dff, model, alpha, beta,-1,1) / normalizationh0
  }

  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value


  return(x)

}


re_N_finder<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                      location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP,e){
  lo = 10
  upper = 5000

  r = re_BF_bound_10(D,lo,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  TPE_lo <- if (de_an_prior == 1)
    re_TPE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e) else
      re_TPE(r,lo,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
  FPE_lo <-  re_FPE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)

  if (TPE_lo > target && FPE_lo < FP) {
    return(lo)
  } else if (TPE_lo > target) {
    alpha.root <- function(n) {
      r <- re_BF_bound_10(D, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model, e)
      re_FPE(r, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model, e) - FP
    }
    return(stats::uniroot(alpha.root, lower = lo, upper = upper)$root)
  }



  Power_root <- function(N){

    r = re_BF_bound_10(D,N,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)

    if (de_an_prior == 0 ){
      pro = re_TPE(r,N,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    }else {
      pro = re_TPE(r,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)

      }
return(pro-target)
    }

  N.power <- robust_uniroot(Power_root, lower = lo)
  r = re_BF_bound_10(D, N.power,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  FPE = re_FPE(r, N.power,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  if (FPE <= FP) return(N.power)

  alpha.root <- function(n) {
    r <- re_BF_bound_10(D, n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    re_FPE(r, n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)-FP
  }
  N.alpha = stats::uniroot(alpha.root,lower = N.power,upper = upper)$root
  return(N.alpha)
}

re_N_01_finder<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                         location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP,e){
  lo = 10
  upper = 5000

  r = re_BF_bound_01(D,lo,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  TNE_lo <- re_TNE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  FNE_lo <-  if (de_an_prior == 1)
    re_FNE(r,lo,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e) else
      re_FNE(r,lo,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)


  if (TNE_lo > target && FNE_lo < FP) {
    return(lo)
  } else if (TNE_lo > target) {
    FN.root <- function(n) {
      r  <- re_BF_bound_01(D, n, k, alpha, beta, h0, hypothesis, location, scale, dff, model, e)
      FNE<-  if (de_an_prior == 0 ){
        pro = re_FNE(r,n,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
      }else {
        pro = re_FNE(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)

      }
      FNE - FP
    }
    return(stats::uniroot(FN.root, lower = lo, upper = upper)$root)
  }


  TN_root <- function(N){

    r = re_BF_bound_10(D,N,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    pro = re_TNE(r,N,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    return(pro-target)
  }

  N.TN <- robust_uniroot(TN_root, lower = lo)
  r    <- re_BF_bound_01(D, N.TN,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  FNE  <- re_FNE(r, N.TN,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  if (FNE <= FP) return(N.TN)

  FN.root <- function(n) {
    r <- re_BF_bound_01(D, n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
    FNE<-  if (de_an_prior == 0 ){
      pro = re_FNE(r,n,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    }else {
      pro = re_FNE(r,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)

    }
    FNE - FP
  }
  N.FN = stats::uniroot(FN.root,lower = N.TN,upper = upper)$root
  return(N.FN)
}

re_table<-function(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                   location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,FP ,e,direct){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)

  n <- if (mode_bf == 1) {
    switch(direct,
           "h1" = ceiling(re_N_finder(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                                      location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP,e)),
           "h0" = ceiling(re_N_01_finder(D,target,model,k, alpha, beta,h0,location,scale,dff, hypothesis ,model_d,
                                      location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,FP,e)))} else  n = N

  # r bounds:
  r10 = re_BF_bound_10(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  r01 = re_BF_bound_01(D,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)

  # max BF10 possible:
  max_BF <- 1 / re_BF10(h0,n,k,alpha, beta,h0,hypothesis,location,scale,dff,model,e)
  BF_D   <- r10

  # FPE and TPE:
  FPE       <- re_FPE(r10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)

  if (de_an_prior == 1) {
      TPE           <- re_TPE(r10,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
      TPR_k         <- k
      TPR_alpha     <- alpha
      TPR_beta      <- beta
      TPR_location  <- location
      TPR_scale     <- scale
      TPR_dff       <- dff
      TPR_model     <- model

  } else {
    TPE           <- re_TPE(r10,n,k_d, alpha_d, beta_d,h0,hypothesis,location_d,scale_d,dff_d,model_d,e)
    TPR_k         <- k_d
    TPR_alpha     <- alpha_d
    TPR_beta      <- beta_d
    TPR_location  <- location_d
    TPR_scale     <- scale_d
    TPR_dff       <- dff_d
    TPR_model     <- model_d
  }
  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- re_FNE(r01,n,TPR_k, TPR_alpha, TPR_beta,h0,hypothesis,TPR_location,TPR_scale,TPR_dff,TPR_model,e)
    TNE <- re_TNE(r01,n,k, alpha, beta,h0,hypothesis,location,scale,dff,model,e)
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


compute.prior.density.re.h1 <- function(rho,h0, k,location,scale,dff,model, alpha, beta,hypothesis,e) {
  if (model == "Point") return(rep(NA, length(rho)))
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = h0+e, b = 1),
                      "<" = c(a = -1, b = h0+e),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh1 <- switch(hypothesis,
                            "!=" = switch(model,
                                          "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "Moment"          = {
                                            (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                              (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                          }),

                            "<"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "Moment"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          }),
                            ">"  = switch(model,
                                          "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                          "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                          "Moment"          = {
                                            (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                          })
  )

  #r_prior(rho,k,location,scale,dff,model, alpha, beta,min(bound),max(bound)) / normalization

  prior_h1<-r_prior(rho,k,location,scale,dff,model, alpha, beta,-1,1)
  switch(hypothesis,
         "!=" = { prior_h1[rho>min(bound_h1)&rho<max(bound_h1)]=0 },
         ">" = { prior_h1[rho<bound_h1[1]]=0 },
         "<" = { prior_h1[rho>bound_h1[2]]=0 }
  )
  prior_h1

  }
compute.prior.density.re.h0 <- function(rho,h0, k,location,scale,dff,model, alpha, beta,hypothesis,e) {
  if (model == "Point") return(rep(NA, length(rho)))
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = h0, b = h0+e),
                      "<" = c(a = h0+e, b = h0),
                      "!=" = c(a = h0+e[1], b = h0+e[2])
  )
  normalizationh0 <- switch(model,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )

  #r_prior(rho,k,location,scale,dff,model, alpha, beta,min(bound),max(bound)) / normalization

  prior_h0<-r_prior(rho,k,location,scale,dff,model, alpha, beta,-1,1)
  switch(hypothesis,
         "!=" = { prior_h0[rho<min(bound_h0)|rho>max(bound_h0)]=0 },
         ">" = { prior_h0[rho>bound_h0[2]]=0 },
         "<" = { prior_h0[rho<bound_h0[1]]=0 }
  )
  prior_h0

}



re_prior_plot <-function(k, alpha, beta,h0,location,scale,dff,model,de_an_prior,
                         k_d, alpha_d, beta_d,location_d,scale_d,dff_d,model_d,hypothesis,e){
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 1))


  plot.bounds <- switch(hypothesis,
                        ">"  = c(h0, 1),
                        "<"  = c(-1, h0),
                        "!=" = c(-1, 1))
  rr <- seq(plot.bounds[1], plot.bounds[2], 0.0025)

  prior.analysis.h1 <- compute.prior.density.re.h1(rr,h0, k,location,scale,dff,model, alpha, beta,hypothesis,e)
  prior.analysis.h0<- compute.prior.density.re.h0(rr,h0, k,location,scale,dff,model, alpha, beta,hypothesis,e)
  prior.design <- if (de_an_prior == 0 && model_d != "Point") {
    compute.prior.density.re.h1(rr, h0,k_d,location_d,scale_d,dff_d,model_d, alpha_d, beta_d,hypothesis,e)
  } else {
    rep(NA, length(rr))
  }
  # Combine all values into one vector
  all_vals <- c(prior.analysis.h1, prior.analysis.h0, prior.design)

  # Filter out NA and infinite values
  finite_vals <- all_vals[is.finite(all_vals)]

  # Get the max from finite values only
  ylim.max <- max(finite_vals)


  # Base plot
  plot(rr, prior.analysis.h1, type = "l", lwd = 2,
       xlab = expression(bold(rho)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~rho~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))

  graphics::lines(rr, prior.analysis.h0, lty = 2, col = "black", lwd = 2)

  # Optional: design prior
  legend.labels <- c("H1 - Analysis Prior", "H0 - Analysis Prior")
  legend.cols   <- c("black", "black")
  legend.lty    <- c(1, 2)
  legend.lwd    <- c(2, 2)

  if (de_an_prior == 0) {
    if (model_d == "Point") {
      graphics::arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim.max,
             length = 0.1, col = "gray", lty = 2)
    } else {
      graphics::lines(rr, prior.design, lty = 1, col = "gray", lwd = 3)
    }

    # Add design prior to legend
    legend.labels <- c(legend.labels, "Design prior")
    legend.cols   <- c(legend.cols, "gray")
    legend.lty    <- c(legend.lty, 1)
    legend.lwd    <- c(legend.lwd, 2)
  }

  graphics::legend("topleft",
         legend = legend.labels,
         col = legend.cols,
         lty = legend.lty,
         lwd = legend.lwd,
         bty = "n")
}



re_bf10_p <- function(D, n, k, alpha, beta, h0, hypothesis,
                      location, scale, dff, model, e) {

  rr <- seq(-0.99, 0.99, 0.01)

  # Compute BF10 and bounds
  BF10   <- re_BF10(rr, n, k, alpha, beta, h0, hypothesis,
                    location, scale, dff, model, e)
  r.BF10 <- re_BF_bound_10(D, n, k, alpha, beta, h0, hypothesis,
                           location, scale, dff, model, e)

  BF01   <- 1 / BF10
  r.BF01 <- re_BF_bound_01(D, n, k, alpha, beta, h0, hypothesis,
                           location, scale, dff, model, e)

  max.BF01 <- 1 / re_BF10(h0, n, k, alpha, beta, h0, hypothesis,
                          location, scale, dff, model, e)
  impossible <- (hypothesis == "!=") &&
    (max.BF01 < D || identical(r.BF01, "bound cannot be found"))

  ## ---------- Titles ----------
  main.bf10 <- if (length(r.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(D) ~ " when r = " ~ .(round(r.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(D) ~ " when r = " ~ .(round(r.BF10[1], 2)) ~
                  " or " ~ .(round(r.BF10[2], 2))))
  }

  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(D)))
  } else if (length(r.BF01) == 1) {
    bquote(bold("BF"[0][1] ~ "=" ~ .(D) ~ " when r = " ~ .(round(r.BF01, 2))))
  } else {
    bquote(bold("BF"[0][1] ~ "=" ~ .(D) ~ " when r = " ~ .(round(r.BF01[1], 2)) ~
                  " or " ~ .(round(r.BF01[2], 2))))
  }



  df10 <- data.frame(r = rr, BF = BF10)
  df01 <- data.frame(r = rr, BF = BF01)

  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text  = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # BF10 plot
  p1 <- ggplot2::ggplot(df10, ggplot2::aes(r, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = r.BF10, linetype = "dashed", color = "black") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = sort(unique(c(-1, 1, round(r.BF10, 2))))) +
    ggplot2::labs(
      x = "Correlation",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    clean_theme

  # BF01 plot
  x_breaks_01 <- if (impossible) c(-1, 1) else sort(unique(c(-1, 1, round(r.BF01, 2))))

  p2 <- ggplot2::ggplot(df01, ggplot2::aes(r, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +   # all black line
    ggplot2::geom_vline(
      xintercept = if (!impossible) r.BF01 else NA,
      linetype = "dashed", color = "black"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = x_breaks_01) +
    ggplot2::labs(
      x = "Correlation",
      y = expression("BF"[0][1] * " (log scale)"),
      title = main.bf01
    ) +
    clean_theme

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

Power_re <- function(D, k, alpha, beta, h0, hypothesis,
                     location, scale, dff, model,
                     k_d, alpha_d, beta_d, location_d, scale_d, dff_d, model_d,
                     de_an_prior, N, e) {

  # N range
  N.min <- 4
  N.max <- ceiling(N * 1.2)
  Ns <- seq(N.min, N.max, length.out = 31)

  TPE <- FPE <- TNE <- FNE <- numeric(length(Ns))

  for (i in seq_along(Ns)) {
    r10 <- re_BF_bound_10(D, Ns[i], k, alpha, beta, h0, hypothesis, location, scale, dff, model, e)
    r01 <- re_BF_bound_01(D, Ns[i], k, alpha, beta, h0, hypothesis, location, scale, dff, model, e)

    TPE[i] <- if (de_an_prior == 1)
      re_TPE(r10, Ns[i], k, alpha, beta, h0, hypothesis, location, scale, dff, model, e) else
        re_TPE(r10, Ns[i], k_d, alpha_d, beta_d, h0, hypothesis, location_d, scale_d, dff_d, model_d, e)

    FPE[i] <- re_FPE(r10, Ns[i], k, alpha, beta, h0, hypothesis, location, scale, dff, model, e)
    TNE[i] <- re_TNE(r01, Ns[i], k, alpha, beta, h0, hypothesis, location, scale, dff, model, e)

    FNE[i] <- if (de_an_prior == 1)
      re_FNE(r01, Ns[i], k, alpha, beta, h0, hypothesis, location, scale, dff, model, e) else
        re_FNE(r01, Ns[i], k_d, alpha_d, beta_d, h0, hypothesis, location_d, scale_d, dff_d, model_d, e)
  }

  # Prepare data frames for ggplot
  df1 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = Ns,
      `True Positive` = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df1$Type <- factor(df1$Type, levels = c("True Positive", "False Positive"))

  df2 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = Ns,
      `True Negative` = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df2$Type <- factor(df2$Type, levels = c("True Negative", "False Negative"))

  type_colors <- c(
    "True Positive" = "black",
    "False Positive" = "grey50",
    "True Negative" = "black",
    "False Negative" = "grey50"
  )

  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  p1 <- ggplot2::ggplot(df1, ggplot2::aes(SampleSize, Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme + legend_theme

  p2 <- ggplot2::ggplot(df2, ggplot2::aes(SampleSize, Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme + legend_theme

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}



# ---- onesample.r ----
pmom <- function(q,V1=1,tau=1) {

  z <- .5-(stats::pnorm(abs(q)/sqrt(V1*tau)) - abs(q)/sqrt(2*pi*V1*tau) * exp(-.5*q^2/(tau*V1)) - .5)
  return(z*(q<=0)+(1-z)*(q>0))
}

# Probability density function of non-local prior:
dMoment <-function(delta,mu,ta){
  ((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
}
# Probability density function of informed t prior:
tstude <- function(t, location = 0, scale = sqrt(2)/2, df = 1) {
  gamma((df+1)/2) * ((df+((t-location)/scale)^2)/df)^(-((df+1)/2)) / (scale*sqrt(df*pi)*gamma(df/2))
  #stats::dt((t-location)/scale, df, ncp = 0)/scale
}

t1_prior<- function(delta, location, scale, dff, model){
  switch(model,
         "Cauchy"         = tstude(delta, location, scale, 1),
         "Normal"         = stats::dnorm (delta, location, scale),
         "Moment"            = dMoment  (delta, location, scale),
         "t-distribution" = tstude(delta, location, scale, dff))
}

t1_BF10 <-function(t, df, model, location, scale, dff, hypothesis){
  bound  <- switch(hypothesis,
                   ">"  = c(a = 0, b = Inf),
                   "<"  = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x <- numeric(length(t))

  # Normalize the prior outside the for-loop:
  # normalization  <- stats::integrate(function(delta) t1_prior(delta, location, scale, dff, model),lower = bound[1], upper = bound[2])$value
  # For all priors, the prior integrates to 1 when a = -Inf, b = Inf.
  # For all priors, we use their CDFs when either a = 0 or b = 0 and minimize manual integrations.
  # Note: pmom() errors at -Inf and Inf, so we avoid it below.
  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

  for(i in 1:length(t)){
    # int  <- function(delta) stats::dt(t[i], df, ncp = delta * sqrt(df+1)) * t1_prior(delta, location, scale, dff, model)/normalization
    int <- function(delta) stats::dt(t[i], df, ncp = delta * sqrt(df + 1)) * t1_prior(delta, location, scale, dff, model) / normalization

    # Removed stop.on.error = FALSE as it is bad form.
    # Below, I increased rel.tol. It gives good enough precision, and the app becomes quite faster:
    x[i] <- stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = 1e-5)$value / stats::dt(t[i], df, ncp = 0)
  }
  return(x)
}



# for finding the t value such that BF10 = D (code stats::optimized):
t1_BF10_bound <- function(D, df, model, location, scale, dff, hypothesis) {
  Bound_finding <- function(t) t1_BF10(t, df, model, location, scale, dff, hypothesis) - D

    x <- tryCatch(stats::uniroot(Bound_finding, lower = -7, upper = 0)$root, error = function(e) NA)
    y <- tryCatch(stats::uniroot(Bound_finding, lower =  0, upper = 7)$root, error = function(e) NA)
    results <- c(x, y)

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
# t is the t-value lead to BF = b based on the bound functions (stats::optimized):
t1_TNE <- function(t, df,hypothesis) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (length(t) == 2) return(stats::pt(max(t), df) - stats::pt(min(t), df))

  # length(t) = 1:
  return(if (hypothesis==">") stats::pt(t, df,0) else 1 - stats::pt(t, df))
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

  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

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
  #error <- if (model == "Moment" && scale < 0.3) 1e-14 else if (scale > 0.3) .Machine$double.eps^0.25 else 1e-8
  error = 1e-4
  stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value
}

# p(BF01>D|H1)
# Similar as above:
t1_FNE <- function(t, df, model, location, scale, dff,hypothesis){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (model == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(max(t), df, ncp) - pnct(min(t), df, ncp))
    # Length 1:
    return(if (hypothesis==">") pnct(t, df, ncp) else 1 - pnct(t, df, ncp))
  }

  bound  <- switch(hypothesis,
                   ">"  = c(a = 0,    b = Inf),
                   "<"  = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf))


  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

  int <- if (hypothesis =="!=") { # two-sided test
    function(delta) {
      pro1 <- pnct(max(t), df, delta * sqrt(df + 1))
      pro2 <- pnct(min(t), df, delta * sqrt(df + 1))
      (pro1 - pro2) * t1_prior(delta, location, scale, dff, model) / normalization
    }
  } else if (hypothesis ==">") { # one-sided test with delta > 0
    function(delta) pnct(t, df, delta * sqrt(df + 1)) * t1_prior(delta, location, scale, dff, model) / normalization
  } else {             # one-sided test with delta < 0
    function(delta) (1 - pnct(t, df, delta * sqrt(df + 1))) * t1_prior(delta, location, scale, dff, model) / normalization
  }

  # setting error value such that error are prevented:
  #error <- if (model == "Moment" && scale < 0.3) 1e-14 else if (scale > 0.3) .Machine$double.eps^0.25 else 1e-8
  error = 1e-4
  stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value
}

# p(BF10>D|H0)
t1_FPE <- function(t, df,hypothesis) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  # if (length(t) == 4) t <- t[2:3]

  if (length(t) == 2) return(stats::pt(min(t), df) + (1 - stats::pt(max(t), df)))

  # length(t) = 1:
  return(if (hypothesis==">") 1 - stats::pt(t, df) else stats::pt(t, df))
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
  upper <- 10000

  t2 <- t1_BF10_bound(D, df = lower, model, location, scale, dff,hypothesis)
  p2 <- if (de_an_prior == 1)
    t1_TPE(t2, df = lower, model, location, scale, dff ) else
      t1_TPE(t2, df = lower, model_d, location_d, scale_d, dff_d)
  if (p2 > target) return(lower)
  Power_root <- function(df) {
    t <- t1_BF10_bound(D, df, model, location, scale, dff,hypothesis)
    pro <- if (de_an_prior == 1)
      t1_TPE(t, df, model, location, scale, dff)  else
        t1_TPE(t, df, model_d, location_d, scale_d, dff_d)
    pro-target
  }

  ## finding the required df, i will do the plus one to get the N in the later function.
  # Jorge: It's a pity if we don't fix it here already, super easy to do right now.
  #        So I did it, see below. I'll adapt latter functions if needed be.

  # Jorge: 'df.power' makes for a more accurate name than 'N'.message("Power at lower = ", Power_root(lower))
  df.power <- stats::uniroot(Power_root, lower = lower, upper = upper)$root

  ## checking if the N lead to an acceptable alpha level
  t   <- t1_BF10_bound(D, df.power, model, location, scale, dff,hypothesis)
  FPE <- t1_FPE(t, df.power,hypothesis)
  if (FPE <= alpha) return(df.power + 1)

  # if the FPE > alpha, then we search for another df
  # Jorge: 'alpha.root' is better than 'alpha_bound'.
  alpha.root <- function(df) {
    t <- t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)
    t1_FPE(t, df,hypothesis) - alpha
  }

  # Jorge: 'df.alpha' is better than 'NN'.
  df.alpha <- stats::uniroot(alpha.root, lower = df.power, upper = upper)$root
  return(df.alpha + 1)
  }

t1_N_01_finder <- function(D, target, model, location, scale, dff, hypothesis,
                           model_d, location_d, scale_d, dff_d, de_an_prior, alpha) {
  lower <- 2
  upper <- 10000

  t2 <- t1_BF01_bound(D, df = lower, model, location, scale, dff,hypothesis)
  TNE_lo <- t1_TNE(t2, df = lower,hypothesis)
  if (TNE_lo > target) return(lower)

  FNE_lo <-  if (de_an_prior == 1)
    t1_FNE(t2, df = lower, model, location, scale, dff,hypothesis ) else
      t1_FNE(t2, df = lower, model_d, location_d, scale_d, dff_d,hypothesis)
  if (TNE_lo > target&FNE_lo<alpha) return(lower)


  TN_root <- function(df) {
    t <- t1_BF01_bound(D, df, model, location, scale, dff,hypothesis)
    t1_TNE(t, df = df,hypothesis) - target
  }

  df.TN <- stats::uniroot(TN_root, lower = lower, upper = upper)$root

  t   <- t1_BF01_bound(D, df.TN, model, location, scale, dff,hypothesis)
  FNE <- if (de_an_prior == 1)
    t1_FNE(t, df = df.TN, model, location, scale, dff,hypothesis ) else
      t1_FNE(t, df = df.TN, model_d, location_d, scale_d, dff_d,hypothesis)

  if (FNE <= alpha) return(df.TN + 1)

  FN.root <- function(df) {
    t <- t1_BF01_bound(D, df, model, location, scale, dff, hypothesis)
    pro = if (de_an_prior == 1)
      t1_FNE(t, df = df, model, location, scale, dff ,hypothesis) - alpha else
        t1_FNE(t, df = df, model_d, location_d, scale_d, dff_d,hypothesis)
    pro- alpha
  }
  df.FN <- stats::uniroot(FN.root, lower = df.TN, upper = upper)$root
  return(df.FN + 1)
}




############ probability table
# Jorge: I edited so that it used N returned by t1_N_finder().
t1_Table <- function(D, target, model, location, scale, dff, hypothesis,
                     model_d, location_d, scale_d, dff_d, de_an_prior, N, mode_bf, alpha,direct) {

    df <- if (mode_bf == "0") {
    N - 1
  } else {
    fun <- if (direct == "h1") t1_N_finder else t1_N_01_finder
    ceiling(fun(D, target, model, location, scale, dff, hypothesis,
                model_d, location_d, scale_d, dff_d, de_an_prior, alpha)) - 1
  }

  # t bounds:
  t10 <- t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)
  t01 <- t1_BF01_bound(D, df, model, location, scale, dff, hypothesis)

  # max BF10 possible:
  max_BF <- 1 / t1_BF10(0, df, model, location, scale, dff, hypothesis)
  BF_D   <- t10

  # FPE and TPE:
  FPE       <- t1_FPE(t10, df,hypothesis)
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
    TNE <- 0
  } else {
    FNE <- t1_FNE(t01, df, TPR_model, TPR_loc, TPR_scale, TPR_dff,hypothesis)
    TNE <- t1_TNE(t01, df,hypothesis)
  }

  # table:
  tab.names <- c(
    sprintf("p(BF10 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H0)", D),
    sprintf("p(BF10 > %0.f | H0)", D),
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, df+1, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}

# For plotting, compute normalized prior density over tt:
compute.prior.density.t <- function(tt, model, location, scale, dff, hypothesis) {
  if (model == "Point") return(rep(NA, length(tt)))
  bounds <- switch(hypothesis,
                   ">"  = c(0, Inf),
                   "<"  = c(-Inf, 0),
                   "!=" = c(-Inf, Inf))
  norm <- stats::integrate(function(delta) t1_prior(delta, location, scale, dff, model),
                    lower = bounds[1], upper = bounds[2])$value
  t1_prior(tt, location, scale, dff, model) / norm
}

# plot for the selected prior
t1_prior_plot <- function(D, target, model, location, scale, dff, hypothesis,
                          model_d, location_d, scale_d, dff_d, de_an_prior) {
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 1))

  plot.bounds    <- switch(hypothesis,
                           ">"  = c(0, 5),
                           "<"  = c(-5, 0),
                           "!=" = c(-5, 5))
  tt             <- seq(plot.bounds[1], plot.bounds[2], 0.01)
  prior.analysis <- compute.prior.density.t(tt, model, location, scale, dff, hypothesis)
  prior.design   <- if (de_an_prior == 0 && model_d != "Point")
    compute.prior.density.t(tt, model_d, location_d, scale_d, dff_d, hypothesis) else
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
      graphics::arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim.max, length = 0.1, col = "black", lty = 2) else
        graphics::lines(tt, prior.design, lty = 2)

    # Add legend:
    graphics::legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }
}

# plots for showing the relationship between BF and t-values

bf10_t1 <- function(D = 3, df, target = 0.8,
                    model = "NA", location = 0, scale = 0.707,
                    dff = 1, hypothesis) {

  tt <- seq(-5, 5, 0.2)

  # Compute BF10 and bounds
  BF10   <- t1_BF10(tt, df, model, location, scale, dff, hypothesis)
  t.BF10 <- t1_BF10_bound(D, df, model, location, scale, dff, hypothesis)

  BF01   <- 1 / BF10
  t.BF01 <- t1_BF01_bound(D, df, model, location, scale, dff, hypothesis)

  # BF10 title
  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(D) ~ " when t = " ~ .(round(t.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(D) ~ " when t = " ~ .(round(t.BF10[1], 2)) ~
                  " or " ~ .(round(t.BF10[2], 2))))
  }

  # Check if BF01 is impossible
  impossible <- (hypothesis == "!=") &&
    (max(BF01) < D || identical(t.BF01, "bound cannot be found"))

  # BF01 title
  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(D)))
  } else if (length(t.BF01) == 1) {
    bquote(bold("BF"[0][1] ~ "=" ~ .(D) ~ " when t = " ~ .(round(t.BF01, 2))))
  } else {
    bquote(bold("BF"[0][1] ~ "=" ~ .(D) ~ " when t = " ~ .(round(t.BF01[1], 2)) ~
                  " or " ~ .(round(t.BF01[2], 2))))
  }



  df_bf10 <- data.frame(tt = tt, BF = BF10)
  df_bf01 <- data.frame(tt = tt, BF = BF01)

  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text  = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  ## ---------- BF10 ----------
  x_breaks_10 <- sort(unique(c(-5, 5, round(t.BF10, 2))))

  p1 <- ggplot2::ggplot(df_bf10, ggplot2::aes(tt, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = t.BF10, linetype = "dashed") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = x_breaks_10) +
    ggplot2::labs(
      x = "t-value",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    clean_theme

  ## ---------- BF01 ----------
  x_breaks_01 <- if (impossible) c(-5, 5)
  else sort(unique(c(-5, 5, round(t.BF01, 2))))

  p2 <- ggplot2::ggplot(df_bf01, ggplot2::aes(tt, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +   # ← FIXED LINE COLOR
    ggplot2::geom_vline(
      xintercept = if (!impossible) t.BF01 else NA,
      linetype = "dashed"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = x_breaks_01) +
    ggplot2::labs(
      x = "t-value",
      y = expression("BF"[0][1] * " (log scale)"),
      title = main.bf01
    ) +
    clean_theme

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

# Power curve function for BF10 > D under H1:
Power_t1 <- function(D, model, location, scale, dff, hypothesis,
                        model_d, location_d, scale_d, dff_d,
                        de_an_prior, N) {

  # df range
  df.min <- 2
  df.max <- ceiling(N * 1.2)
  dfs <- seq(df.min, df.max, length.out = 31)

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(dfs))

  for (i in seq_along(dfs)) {
    t10 <- t1_BF10_bound(D, dfs[i], model, location, scale, dff, hypothesis)
    t01 <- t1_BF01_bound(D, dfs[i], model, location, scale, dff, hypothesis)

    TPE[i] <- if (de_an_prior == 1) {
      t1_TPE(t10, dfs[i], model, location, scale, dff)
    } else {
      t1_TPE(t10, dfs[i], model_d, location_d, scale_d, dff_d)
    }

    FPE[i] <- t1_FPE(t10, dfs[i], hypothesis)
    TNE[i] <- t1_TNE(t01, dfs[i], hypothesis)
    FNE[i] <- if (de_an_prior == 1) {
      t1_FNE(t01, dfs[i], model, location, scale, dff, hypothesis)
    } else {
      t1_FNE(t01, dfs[i], model_d, location_d, scale_d, dff_d, hypothesis)
    }
  }

  # Prepare data for ggplot (check.names = FALSE avoids column name issues)
  df1 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = dfs + 1,
      `True Positive` = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df1$Type <- factor(df1$Type, levels = c("True Positive", "False Positive"))

  df2 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = dfs + 1,
      `True Negative` = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df2$Type <- factor(df2$Type, levels = c("True Negative", "False Negative"))

  # Colors for lines
  type_colors <- c(
    "True Positive" = "black",
    "False Positive" = "grey50",
    "True Negative" = "black",
    "False Negative" = "grey50"
  )

  # Clean theme with larger axis labels and no background grid
  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"), # increased font size for x-axis
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Legend theme: inside top-left, no border, black first then gray
  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  # Plot BF10
  p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # Plot BF01
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots side by side
  combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 2)

  print(combined_plot)
}










# ---- onesample_e.r ----
robust_uniroot <- function(f, lower, upper_start = 500, max_attempts = 20, step = 500, ...) {
  upper <- upper_start
  attempt <- 1

  repeat {
    result <- tryCatch(
      {
        stats::uniroot(f, lower = lower, upper = upper, ...)$root
      },
      error = function(e) {
        NULL
      }
    )

    if (!is.null(result)) {
      return(result)
    }

    if (attempt >= max_attempts) {
      stop("Failed to find root after increasing upper bound ", max_attempts, " times.")
    }

    upper <- upper + step
    attempt <- attempt + 1
  }
}

te_prior<- function(delta,location,scale,dff,model){

  switch(model,
        "Cauchy"        = tstude(delta,location, scale,1),
        "Normal"        = stats::dnorm(delta,location,scale),
        "Moment"            = dMoment(delta,location,scale),
        "t-distribution" = tstude(delta,location,scale,dff))


}
norm_h1 <- function(hypothesis, model, bound_h1, location, scale, dff = NULL) {
  normalizationh1 <- switch(
    hypothesis,

    "!=" = 1 - switch(
      model,
      "Cauchy"         = stats::pcauchy(bound_h1[2], location, scale) -
        stats::pcauchy(bound_h1[1], location, scale),
      "Normal"         = stats::pnorm(bound_h1[2], location, scale) -
        stats::pnorm(bound_h1[1], location, scale),
      "Moment"            = pmom(bound_h1[2] - location, tau = scale^2) -
        pmom(bound_h1[1] - location, tau = scale^2),
      "t-distribution" = stats::pt((bound_h1[2] - location)/scale, df = dff) -
        stats::pt((bound_h1[1] - location)/scale, df = dff)
    ),

    "<" = switch(
      model,
      "Cauchy"         = stats::pcauchy(bound_h1[2], location, scale) -
        stats::pcauchy(bound_h1[1], location, scale),
      "Normal"         = stats::pnorm(bound_h1[2], location, scale) -
        stats::pnorm(bound_h1[1], location, scale),
      "Moment"            = pmom(bound_h1[2] - location, tau = scale^2),
      "t-distribution" = stats::pt((bound_h1[2] - location)/scale, df = dff) -
        stats::pt((bound_h1[1] - location)/scale, df = dff)
    ),

    ">" = switch(
      model,
      "Cauchy"         = stats::pcauchy(bound_h1[2], location, scale) -
        stats::pcauchy(bound_h1[1], location, scale),
      "Normal"         = stats::pnorm(bound_h1[2], location, scale) -
        stats::pnorm(bound_h1[1], location, scale),
      "Moment"            = 1 - pmom(bound_h1[1] - location, tau = scale^2),
      "t-distribution" = stats::pt((bound_h1[2] - location)/scale, df = dff) -
        stats::pt((bound_h1[1] - location)/scale, df = dff)
    )
  )

  return(normalizationh1)
}
norm_h0 <- function(model, bound_h0, location, scale, dff = NULL) {
  normalizationh0 <- switch(
    model,

    "Cauchy"         = stats::pcauchy(bound_h0[2], location, scale) -
      stats::pcauchy(bound_h0[1], location, scale),

    "Normal"         = stats::pnorm(bound_h0[2], location, scale) -
      stats::pnorm(bound_h0[1], location, scale),

    "Moment"            = pmom(bound_h0[2] - location, tau = scale^2) -
      pmom(bound_h0[1] - location, tau = scale^2),

    "t-distribution" = stats::pt((bound_h0[2] - location)/scale, df = dff) -
      stats::pt((bound_h0[1] - location)/scale, df = dff)
  )

  return(normalizationh0)
}


t1e_BF10i <-function(t,df,model ,location, scale,dff , hypothesis,e ){
  bound_h1  <- switch(hypothesis,
                      "!=" = c(a = e[1], b = e[2]),
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = e)
                      )

  bound_h0  <- switch(hypothesis,
                      "!=" = c(a = e[1], b = e[2]),
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0)
                      )

  normalizationh1 <- norm_h1(hypothesis, model, bound_h1, location, scale, dff)



  # H0 Normalization
  normalizationh0 <- norm_h0(model, bound_h0, location, scale, dff)

  int  <- function(delta){
    stats::dt(t,df,ncp = delta *sqrt(df+1))* te_prior(delta,location,scale,dff,model)/normalizationh1
    }

   error = 1e-4

  if (hypothesis == "!="){
  lh1 = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value
  }else{
    lh1 = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value

  }


  int  <- function(delta){
    stats::dt(t,df,ncp = delta *sqrt(df+1))* te_prior(delta,location,scale,dff,model)/normalizationh0}

  lh0 = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(lh1/lh0)
}

t1e_BF10 <-function(t,df,model,location, scale,dff , hypothesis,e ){

  x <- sapply(t, function(ti) t1e_BF10i(ti, df, model,location, scale, dff, hypothesis, e))
  return(x)
}
#
t1e_BF10_bound <-function(D, df,model,location,scale,dff , hypothesis,e){
  y <- numeric(0)
  Bound_finding <-function(t){
    t1e_BF10(t,df,model,location,scale,dff , hypothesis,e )- D
  }

  switch(hypothesis,
         "!=" ={
           x <- tryCatch(stats::uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
           y <- tryCatch(stats::uniroot(Bound_finding, lower =  0, upper = 20)$root, error = function(e) NA)
         },
         ">"={
           x <- tryCatch(stats::uniroot(Bound_finding, lower = 0, upper = 20)$root, error = function(e) NA)
         },
         "<" = {
           x <- tryCatch(stats::uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
         })


  results <- c(x, y)

  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")

  BF.vals  <- t1e_BF10(results,df,model,location,scale,dff , hypothesis,e )
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])

}


t1e_BF01_bound <-function(D, df,model,location,scale,dff , hypothesis,e){
  t1e_BF10_bound(1/D, df,model,location,scale,dff , hypothesis,e)
}



t1e_TPE <-function(t,df,model ,location,scale,dff , hypothesis ,e){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (model == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(min(t), df, ncp) + (1 - pnct(max(t), df, ncp)))
    # Length 1:
    return(if (t >= 0) 1 - pnct(t, df, ncp) else pnct(t, df, ncp))
  }


  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = e),
                      "!=" = c(a = e[1], b = e[2])
  )

  normalizationh1 <- norm_h1(hypothesis, model, bound_h1, location, scale, dff)


  x = NULL

  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)

    pro <- switch(hypothesis,
                  "!=" = pnct(max(t), df, ncp = ncp, lower  = FALSE) +
                    pnct(min(t), df, ncp = ncp, lower  = TRUE),
                  ">"  = pnct(t, df, ncp = ncp, lower  = FALSE),
                  "<"  = pnct(t, df, ncp = ncp, lower  = TRUE)
    )

    pro * te_prior(delta, location,scale, dff, model) / normalizationh1
  }

   error = 1e-4

  x <- switch(hypothesis,
              "!=" = stats::integrate(int, -Inf, bound_h1[1], rel.tol = error)$value +
                stats::integrate(int, bound_h1[2], Inf, rel.tol = error)$value,
              "<"  = ,
              ">"  = stats::integrate(int, bound_h1[1], bound_h1[2], rel.tol = error)$value

  )
  return(x)

}

t1e_FNE <-function(t,df,model ,location,scale,dff , hypothesis ,e){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (model == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(max(t), df, ncp) - pnct(min(t), df, ncp))
    # Length 1:
    return(if (t >= 0) pnct(t, df, ncp) else 1 - pnct(t, df, ncp))
  }

  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = e),
                      "!=" = c(a = e[1], b = e[2])
  )

  normalizationh1 <- norm_h1(hypothesis, model, bound_h1, location, scale, dff)

  x = NULL

  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)

    pro <- switch(hypothesis,
                   "!=" = pnct(max(t), df, ncp = ncp, lower  = TRUE) -
                     pnct(min(t), df, ncp = ncp, lower  = TRUE),
                   ">"  = pnct(t, df, ncp = ncp, lower  = TRUE),
                   "<"  = pnct(t,df, ncp = ncp, lower  = FALSE)
    )

    pro * te_prior(delta,location, scale, dff, model) / normalizationh1
  }

   error = 1e-4

  x <- switch(hypothesis,
              "!=" = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value ,
              "<"  = ,
              ">"  = stats::integrate(int, bound_h1[1], bound_h1[2], rel.tol = error)$value)
  return(x)

}

t1e_TNE <-function(t,df,model ,location,scale,dff , hypothesis ,e){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )
  normalizationh0 <- norm_h0(model, bound_h0, location, scale, dff)

  x = NULL

  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)
    pro <- switch(hypothesis,
                   "!=" = pnct(max(t), df, ncp, lower  = TRUE) -
                     pnct(min(t), df, ncp, lower  = TRUE),
                   ">"  = pnct(t, df, ncp, lower  = TRUE),
                   "<"  = pnct(t, df, ncp, lower  = FALSE),
                   stop("Unsupported hypothesis")
    )

    pro * te_prior(delta,location, scale, dff, model) / normalizationh0
  }

   error = 1e-4

  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error)$value

  return(x)

}

t1e_FPE <-function(t,df,model ,location,scale,dff , hypothesis ,e){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )

  normalizationh0 <- norm_h0(model, bound_h0, location, scale, dff)
  x = NULL
  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)

    pro <- switch(hypothesis,
                   "!=" = pnct(max(t), df, ncp, lower  = FALSE) + pnct(min(t), df, ncp, lower  = TRUE),
                   ">"  = pnct(t, df, ncp, lower  = FALSE),
                   "<"  = pnct(t, df, ncp, lower  = TRUE),
                   stop("Unsupported hypothesis")
    )

    pro * te_prior(delta,location, scale, dff, model) / normalizationh0
  }

   error = 1e-4

  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value

  return(x)

}

t1e_N_finder<-function(D,target,model,location,scale,dff, hypothesis,e ,
                   model_d,scale_d,dff_d, de_an_prior,location_d  ,alpha){

  lower <- 2
  upper <- 10000
  t2 <-t1e_BF10_bound(D, lower,model,location,scale,dff , hypothesis,e)
  p2 <- if (de_an_prior == 1)
    t1e_TPE(t2,lower,model ,location,scale,dff , hypothesis ,e) else
    t1e_TPE(t2,lower,model,location_d ,scale_d,dff_d , hypothesis ,e)
  if (p2 > target) return(lower)

  Power_root <- function(df) {

    t <- t1e_BF10_bound(D, df, model,location, scale, dff, hypothesis, e)

    pro <- if (de_an_prior == 1) {
      t1e_TPE(t, df, model,location, scale, dff, hypothesis, e)
    } else {
      t1e_TPE(t, df, model_d,location_d, scale_d, dff_d, hypothesis, e)
    }

    target - pro
  }

  df.power <- robust_uniroot(Power_root, lower = 2)
  t <- t1e_BF10_bound(D,df.power,model,location,scale,dff,hypothesis ,e )
  FPE <-t1e_FPE(t,df.power,model,location ,scale,dff , hypothesis ,e)
  if (FPE <= alpha) return(df.power + 1)

  alpha.root <- function(df) {
      t <- t1e_BF10_bound(D,df,model,location,scale,dff,hypothesis ,e )
      pro <- t1e_FPE(t , df , model,location , scale,dff, hypothesis,e)
      return(pro - alpha)
    }
  df.alpha <- stats::uniroot(alpha.root, lower = df.power, upper = upper)$root
  return(df.alpha+1)

}

t1e_N_01_finder<-function(D,target,model,location,scale,dff, hypothesis,e ,
                          model_d,scale_d,dff_d, de_an_prior,location_d  ,alpha){

  lower <- 10
  upper <- 10000
  t2 <-t1e_BF01_bound(D, lower,model,location,scale,dff , hypothesis,e)
  TNE_lo <-t1e_TNE(t2,lower,model,location ,scale,dff , hypothesis ,e)
  if (TNE_lo > target) return(lower)

  FNE_lo <- if (de_an_prior == 1)
    t1e_FPE(t2,lower,model,location ,scale,dff , hypothesis ,e) else
      t1e_FPE(t2,lower,model,location_d ,scale_d,dff_d , hypothesis ,e)
  if (TNE_lo > target&FNE_lo<alpha) return(lower)

  TN_root <- function(df) {
    t <- t1e_BF01_bound(D, df, model,location, scale, dff, hypothesis, e)

    pro <-t1e_TNE(t,df,model,location ,scale,dff , hypothesis ,e)
    target-pro
  }

  df.TN <- robust_uniroot(TN_root, lower = lower)
  t <- t1e_BF01_bound(D,df.TN,model,location,scale,dff,hypothesis ,e )
  FNE <-t1e_FNE(t,df.TN,model,location ,scale,dff , hypothesis ,e)
  if (FNE <= alpha) return(df.TN + 1)

  FN.root <- function(df) {
    t <- t1e_BF01_bound(D,df,model,location,scale,dff,hypothesis ,e )
    pro <-   if (de_an_prior == 1) {
      t1e_FNE(t, df, model,location, scale, dff, hypothesis, e)
    } else {
      t1e_FNE(t, df, model_d,location_d, scale_d, dff_d, hypothesis, e)
    }
    return(pro - alpha)
  }
  df.FN <- robust_uniroot(FN.root, lower = df.TN )
  return(df.FN+1)

}
t1e_table<-function(D,target,model,location,scale,dff, hypothesis,e ,
                    model_d,scale_d,dff_d, de_an_prior,N,mode_bf,location_d ,alpha,direct ){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)

  df <- if (mode_bf == 0) {
    N - 1
  } else {
    fun <- if (direct == "h1") t1e_N_finder else t1e_N_01_finder
    ceiling(fun(D,target,model,location,scale,dff, hypothesis,e ,
                model_d,scale_d,dff_d, de_an_prior ,location_d,alpha )) - 1
  }



  # t bounds:
  t10 <- t1e_BF10_bound(D, df,model,location,scale,dff , hypothesis,e)
  t01 <- t1e_BF01_bound(D, df,model,location,scale,dff , hypothesis,e)

  # max BF10 possible:
  max_BF <- 1 / t1e_BF10(0,df,model,location ,scale,dff , hypothesis,e )
  BF_D   <- t10

  # FPE and TPE:
  FPE       <- t1e_FPE(t10,df,model ,location,scale,dff , hypothesis ,e)
  if (de_an_prior == 1) {
    TPE         <- t1e_TPE(t10,df,model ,location,scale,dff , hypothesis ,e)
    TPR_model   <- model
    TPR_scale   <- scale
    TPR_dff     <- dff
    TPR_location<- location
  } else {
    TPE       <- t1e_TPE(t10,df,model_d ,location_d,scale_d,dff_d , hypothesis ,e)
    TPR_model <- model_d
    TPR_location   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }

  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t1e_FNE(t01,df,TPR_model,TPR_location ,TPR_scale,TPR_dff , hypothesis ,e)
    TNE <-  t1e_TNE(t01,df,model,location ,scale,dff , hypothesis ,e)
  }
  # table:
  tab.names <- c(
    sprintf("p(BF10 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H0)", D),
    sprintf("p(BF10 > %0.f | H0)", D),
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, df+1, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table

}

compute.prior.density.te.h1 <- function(tt, model,location, scale, dff, hypothesis,e) {
  if (model == "Point") return(rep(NA, length(tt)))
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = -e),
                      "!=" = c(a = e[1], b = e[2])
  )

  normalizationh1 <- norm_h1(hypothesis, model, bound_h1, location, scale, dff)


  #prior_h1<-te_prior(tt,scale,dff,model) / normalizationh1
  prior_h1<-te_prior(tt,location,scale,dff,model)
  switch(hypothesis,
         "!=" = { prior_h1[tt>min(bound_h1)&tt<max(bound_h1)]=0 },
         ">" = { prior_h1[tt<bound_h1[1]]=0 },
         "<" = { prior_h1[tt>bound_h1[2]]=0 }
         )
  prior_h1
}

compute.prior.density.te.h0 <- function(tt, model,location, scale, dff, hypothesis,e) {
  if (model == "Point") return(rep(NA, length(tt)))
  bound_h0  <- switch(hypothesis,
                      "!=" = c(a = e[1], b = e[2]),
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0)
  )
  # H0 Normalization
  normalizationh0 <- norm_h0(model, bound_h0, location, scale, dff)


  #prior_h0 <- te_prior(tt,scale,dff,model) / normalizationh0
  prior_h0 <- te_prior(tt,location,scale,dff,model)
  switch(hypothesis,
         "!=" = { prior_h0[!(tt>min(bound_h0)&tt<max(bound_h0))]=0},
         ">" = { prior_h0[tt>bound_h0[2]]=0 },
         "<" = { prior_h0[tt<bound_h0[1]]=0 }

  )
  prior_h0
}


t1e_prior_plot <- function(model,location, scale, dff, hypothesis, e,
                           de_an_prior, model_d, scale_d, dff_d, location_d) {
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 1))

  plot.bounds <- switch(hypothesis,
                        ">"  = c(0, 5),
                        "<"  = c(-5, 0),
                        "!=" = c(-5, 5))
  tt <- seq(plot.bounds[1], plot.bounds[2], 0.01)

  prior.analysis.h1 <- compute.prior.density.te.h1(tt, model,location, scale, dff, hypothesis, e)
  prior.analysis.h0 <- compute.prior.density.te.h0(tt, model,location, scale, dff, hypothesis, e)
  prior.design <- if (de_an_prior == 0 && model_d != "Point") {
    compute.prior.density.te.h1(tt, model_d,location_d, scale_d, dff_d, hypothesis, e)
  } else {
    rep(NA, length(tt))
  }

  ylim.max <- max(prior.analysis.h1, prior.analysis.h0, prior.design, na.rm = TRUE)

  # Base plot
  plot(tt, prior.analysis.h1, type = "l", lwd = 2,
       xlab = expression(bold(delta)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~delta~" under the alternative hypothesis")),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))

  graphics::lines(tt, prior.analysis.h0, lty = 2, col = "black", lwd = 2)

  # Optional: design prior
  legend.labels <- c("H1 - Analysis Prior", "H0 - Analysis Prior")
  legend.cols   <- c("black", "black")
  legend.lty    <- c(1, 2)
  legend.lwd    <- c(2, 2)

  if (de_an_prior == 0) {
    if (model_d == "Point") {
      graphics::arrows(x0 = location_d, y0 = 0, x1 = location_d, y1 = ylim.max,
             length = 0.1, col = "gray", lty = 2)
    } else {
      graphics::lines(tt, prior.design, lty = 1, col = "gray", lwd = 3)
    }

    # Add design prior to legend
    legend.labels <- c(legend.labels, "Design prior")
    legend.cols   <- c(legend.cols, "gray")
    legend.lty    <- c(legend.lty, 1)
    legend.lwd    <- c(legend.lwd, 2)
  }

  graphics::legend("topleft",
         legend = legend.labels,
         col = legend.cols,
         lty = legend.lty,
         lwd = legend.lwd,
         bty = "n")
}

te1_BF<- function(D, df,
                          model, location, scale, dff, hypothesis, e) {

  tt <- seq(-5, 5, 0.2)

  ## ---------- BF10 ----------
  BF10   <- t1e_BF10(tt, df, model, location, scale, dff, hypothesis, e)
  t.BF10 <- t1e_BF10_bound(D, df, model, location, scale, dff, hypothesis, e)

  df10 <- data.frame(t = tt, BF = BF10)

  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(D)~" when t = "~.(round(t.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(D)~" when t = "~.(round(t.BF10[1], 2))~
                  " or "~.(round(t.BF10[2], 2))))
  }

  x_breaks_10 <- sort(unique(c(-5, 5, round(t.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(t, BF)) +
    ggplot2::geom_line(size = 1.1) +
    ggplot2::scale_y_log10() +
    ggplot2::geom_vline(xintercept = t.BF10) +
    ggplot2::scale_x_continuous(
      limits = c(-5, 5),
      breaks = x_breaks_10
    ) +
    ggplot2::labs(
      x = "t-value",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 13, face = "bold"),
      axis.text  = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  ## ---------- BF01 ----------
  BF01   <- 1 / BF10
  t.BF01 <- t1e_BF01_bound(D, df, model, location, scale, dff, hypothesis, e)

  df01 <- data.frame(t = tt, BF = BF01)

  max.BF01   <- 1 / t1e_BF10(0, df, model, location, scale, dff, hypothesis, e)
  impossible <- (max.BF01 < D || identical(t.BF01, "bound cannot be found"))

  if (impossible) {

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(size = 1.1) +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(
        limits = c(-5, 5),
        breaks = c(-5, 5)
      ) +
      ggplot2::labs(
        x = "t-value",
        y = bquote("BF"['01'] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[01]~"="~.(D)))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.text  = ggplot2::element_text(size = 11),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

  } else {

    main.bf01 <- if (length(t.BF01) == 1) {
      bquote(bold("BF"[0][1]~"="~.(D)~" when t = "~.(round(t.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(D)~" when t = "~.(round(t.BF01[1], 2))~
                    " or "~.(round(t.BF01[2], 2))))
    }

    x_breaks_01 <- sort(unique(c(-5, 5, round(t.BF01, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(size = 1.1) +
      ggplot2::scale_y_log10() +
      ggplot2::geom_vline(xintercept = t.BF01) +
      ggplot2::scale_x_continuous(
        limits = c(-5, 5),
        breaks = x_breaks_01
      ) +
      ggplot2::labs(
        x = "t-value",
        y = bquote("BF"[0][1] * " (log scale)"),
        title = main.bf01
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.text  = ggplot2::element_text(size = 11),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

Power_t1e<- function(D, model, location, scale, dff, hypothesis,
                         model_d, location_d, scale_d, dff_d,
                         de_an_prior, N, e) {

  # df range
  df.min <- 2
  df.max <- ceiling(N * 1.2)
  dfs <- seq(df.min, df.max, length.out = 30)

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(dfs))

  for (i in seq_along(dfs)) {
    t10 <- t1e_BF10_bound(D, dfs[i], model, location, scale, dff, hypothesis, e)
    t01 <- t1e_BF01_bound(D, dfs[i], model, location, scale, dff, hypothesis, e)

    # Choose correct design prior
    TPE[i] <- if (de_an_prior == 1) {
      t1e_TPE(t10, dfs[i], model, location, scale, dff, hypothesis, e)
    } else {
      t1e_TPE(t10, dfs[i], model_d, location_d, scale_d, dff_d, hypothesis, e)
    }

    FPE[i] <- t1e_FPE(t10, dfs[i], model, location, scale, dff, hypothesis, e)
    TNE[i] <- t1e_TNE(t01, dfs[i], model, location, scale, dff, hypothesis, e)
    FNE[i] <- if (de_an_prior == 1) {
      t1e_FNE(t01, dfs[i], model, location, scale, dff, hypothesis, e)
    } else {
      t1e_FNE(t01, dfs[i], model_d, location_d, scale_d, dff_d, hypothesis, e)
    }
  }

  # Prepare data frames for ggplot
  df1 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = dfs + 1,
      `True Positive` = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df1$Type <- factor(df1$Type, levels = c("True Positive", "False Positive"))

  df2 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = dfs + 1,
      `True Negative` = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df2$Type <- factor(df2$Type, levels = c("True Negative", "False Negative"))

  # Colors
  type_colors <- c(
    "True Positive" = "black",
    "False Positive" = "grey50",
    "True Negative" = "black",
    "False Negative" = "grey50"
  )

  # Clean theme
  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"), # increased font size for x-axis
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Legend theme
  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  # BF10 plot
  p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # BF01 plot
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots side by side
  patchwork::wrap_plots(p1, p2, ncol = 2)
}


# ---- proportions.r ----
BF10_p2<-function(a0, b0, a1, b1, a2, b2,n1,n2,k1,k2){

  logBF = lbeta(k1 + k2 + a0, n1 + n2 - k1 - k2 + b0) -
    lbeta(k1 + a1, n1 - k1 + b1) -
    lbeta(k2 + a2, n2 - k2 + b2) +
    lbeta(a1, b1) +
    lbeta(a2, b2) -
    lbeta(a0, b0)
  1/exp(logBF)
}

ps_N_finder <- function(D, target, a0, b0, a1, b1, a2, b2, r,
                        model1, da1, db1, dp1, model2, da2, db2, dp2) {

  lo_n1 <- 10
  n2 <- round(lo_n1 * r)
  grid <- BF_grid_rcpp(D, a0, b0, a1, b1, lo_n1, a2, b2, n2,
                       model1, da1, db1, dp1, model2, da2, db2, dp2)

  pro <- sum_rcpp(grid$log_h1_dp, grid$PE)
  if (pro > target) return(list(grid, lo_n1))

  # Function for uniroot
  power_fun <- function(n1){
    n1 <- round(n1)
    n2 <- n1 * r
    g <- BF_grid_rcpp(D, a0, b0, a1, b1, n1, a2, b2, n2,
                      model1, da1, db1, dp1, model2, da2, db2, dp2)
    sum_rcpp(g$log_h1_dp, g$PE) - target - 0.005
  }

  n1 <- suppressWarnings(round(stats::uniroot(power_fun, lower = lo_n1, upper = 5000, maxiter = 10)$root))
  n2 <- round(n1 * r)
  grid <- BF_grid_rcpp(D, a0, b0, a1, b1, n1, a2, b2, n2,
                       model1, da1, db1, dp1, model2, da2, db2, dp2)

  list(grid, n1)
}


ps_N_01_finder<-function(D,target, a0, b0, a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2) {

  lo_n1 <- 10
  n2 <- round(lo_n1)*r
  grid <- BF_grid_rcpp(D, a0, b0, a1, b1, lo_n1, a2, b2, n2,model1,da1,db1,dp1,model2,da2,db2,dp2)
  pro <- sum_rcpp(grid$log_h0,grid$NE)

  if ( pro>target){
    return(list(grid,lo_n1))
  }
  TN<-function(n1){
    n1 = round(n1)
    n2 = n1*r
    grid <<- BF_grid_rcpp(D, a0, b0, a1, b1, n1, a2, b2, n2,model1,da1,db1,dp1,model2,da2,db2,dp2)
    pro <- sum_rcpp(grid$log_h0,grid$NE)
    return(pro - target - .01)
  }
  n1 <- suppressWarnings(round(stats::uniroot(TN, lower = lo_n1, upper = 5000,maxiter = 10)$root))
  grid_power <- grid



  return(list(grid_power,n1))

}


pro_table_p2<-function(D,target, a0, b0, a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2,mode_bf,n1,n2,direct) {

  if (mode_bf==1){
    x = switch(direct,
               "h1" = ps_N_finder(D,target, a0, b0, a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2),
               "h0" = ps_N_01_finder(D,target, a0, b0, a1, b1, a2, b2, r,model1,da1,db1,dp1,model2,da2,db2,dp2))
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
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1, 1))

  prop    <- seq( 0,1,.001)

  prior.analysis <- stats::dbeta(prop,a,b)
  prior.design   <- switch(model,
                           "same" = stats::dbeta(prop,a,b),
                           "beta" = stats::dbeta(prop,ad,bd),
                           "Point" = rep(NA, length(prop)))

  ylim.max <- max(prior.analysis, prior.design, na.rm = TRUE)

  plot(prop, prior.analysis, type = "l", lwd = 2,
       xlab = substitute(bold(theta[nu_val]), list(nu_val = nu)),
       ylab = "density",
       main = bquote(bold("Prior distribution on "~theta[.(nu)])),
       frame.plot = FALSE,
       ylim = c(0, ylim.max))


  if (model != "same") {
    if (model == "Point")
      graphics::arrows(x0 = dp, y0 = 0, x1 = dp, y1 = ylim.max, length = 0.1, col = "black", lty = 2) else
        graphics::lines(prop, prior.design, lty = 2)

    # Add legend:
    graphics::legend("topright",
           legend = c("Analysis prior", "Design prior"),
           lty = c(1, 2),
           col = c("black", "black"),
           bty = "n")
  }

}


Power_p2 <- function(D, n1, a0, b0, a1, b1, a2, b2, r,
                     model1, da1, db1, dp1,
                     model2, da2, db2, dp2) {

  # Define sample size range
  smax <- n1 * 1.2
  Ns <- ceiling(seq(10, smax, by = (smax - 10) / 20))
  n2 <- Ns * r
  Nt <- Ns + n2

  # Initialize probability vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(Ns))

  # Compute probabilities for each Ns
  for (i in seq_along(Ns)) {
    n2 <- Ns[i] * r
    grid <- BF_grid_rcpp(D, a0, b0, a1, b1, Ns[i],
                         a2, b2, n2,
                         model1, da1, db1, dp1,
                         model2, da2, db2, dp2)
    TPE[i] <- sum_rcpp(grid$log_h1_dp, grid$PE)
    FPE[i] <- sum_rcpp(grid$log_h0, grid$PE)
    TNE[i] <- sum_rcpp(grid$log_h0, grid$NE)
    FNE[i] <- sum_rcpp(grid$log_h1_dp, grid$NE)
  }

  # Prepare data frames for ggplot
  df1 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = Nt,
      `True Positive` = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df1$Type <- factor(df1$Type, levels = c("True Positive", "False Positive"))

  df2 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = Nt,
      `True Negative` = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df2$Type <- factor(df2$Type, levels = c("True Negative", "False Negative"))

  # Colors
  type_colors <- c(
    "True Positive" = "black",
    "False Positive" = "grey50",
    "True Negative" = "black",
    "False Negative" = "grey50"
  )

  # Clean theme
  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Legend theme
  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  # Plot BF10
  p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # Plot BF01
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots
  combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 2)

  print(combined_plot)
}


heatmap_p2 <- function(x, D) {
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

  # Legend labels (math expressions)
  labels <- c(
    "PE"   = bquote(BF[10] > .(D)),
    "NE"   = bquote(BF[0][1] > .(D)),
    "None" = bquote(1 / .(D) < BF[10] ~ "<" ~ .(D))
  )

  # First plot: categorical heatmap
  p1 <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$k1, y = .data$k2, fill = .data$effect)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(
      name = "Classification",
      values = c("PE" = "#FDE725", "NE" = "#440154", "None" = "#21908C"),
      labels = labels
    ) +
    ggplot2::labs(
      title = "BF and number of success",
      x = "x1", y = "x2"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal()

  # Second plot: continuous heatmap of BF
  p2 <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$k1, y = .data$k2, fill = .data$BF)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "log BF10") +
    ggplot2::labs(
      title = "Heatmap of BF",
      x = "x1", y = "x2"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal()

  # Combine plots side by side
  combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 2)

  # Return combined plot (print when called interactively)
  combined_plot
}



# ---- Server_bin.r ----

server_bin<- function(input, output, session) {
input_bin <- shiny::reactive({
  mode_bf <- switch(input$Modebin,
                    "1" = 1,
                    "2" = 0,
                    "3" = 0)# mode

 direct <- switch(input$bin_direct,
                  "1" = "h1",
                  "0" = "h0")

  interval <- input$h0bin # point null or interval

  hypothesis <- switch(interval,
                       "1" =   switch(input$h1bin,        # direction of the test
                                      "1" = "!=",
                                      "2" =  ">",
                                      "3" =  "<"),
                       "2" = switch(input$h1bine,        # direction of the test
                                    "1" = "!=",
                                    "2" =  ">",
                                    "3" =  "<"))
  h0       <-  input$h0prop
  location <- input$h0prop
  lbbin <- input$lbbine
  ubbin <- input$ubbine

  if ((location+lbbin)<(0)){
    lbbin = lbbin+-1-(location+lbbin)

  }

  if ((location+ubbin)>(+1)){
    ubbin = ubbin+1-(location+ubbin)

  }


  e <- switch(input$h1bine,        # bound for interval test
              "1" = c(lbbin, ubbin),
              "2" = ubbin,
              "3" = lbbin)

  inter <- switch(interval,
                  "1" = input$h1bin,
                  "2" = input$h1bine)


  model <- switch(input$modelbin,
                  "1" = "beta",
                  "2" = "Moment")
  alpha <- input$alphabin
  beta <- input$betabin
  scale <- input$sbin
  de_an_prior <- switch(input$priorbin,
                        "1" = 1,
                        "2" = 0)
  alpha_d <- input$alphabind
  beta_d <- input$betabind
  scale_d <- input$sbind
  model_d <- switch(input$modelbind,
                  "1" = "beta",
                  "2" = "Moment",
                  "3" = "Point")
  location_d <- input$h0bind
  target <- input$powerbin
  FP <- input$FP_bin
  D <- input$bbin
  N <- input$nbin
  Suc <- input$xbin
  pc   <- "1" %in% input$o_plot_bin
  rela <- "2" %in% input$o_plot_bin

  if (mode_bf!=1){
    pc=rela=F
  }
 ############


  # Add all variables to the final list
  list(
    mode_bf = mode_bf,
    direct = direct,
    interval = interval,
    hypothesis =hypothesis,
    h0=h0,
    location = location,
    e = e,
    lbbin = lbbin,
    ubbin = ubbin,
    inter = inter,
    model = model,
    alpha = alpha,
    beta = beta,
    scale = scale,
    de_an_prior = de_an_prior,
    alpha_d = alpha_d,
    beta_d = beta_d,
    scale_d = scale_d,
    location_d = location_d,
    model_d = model_d,
    target = target,
    FP = FP,
    D = D,
    N = N,
    Suc =Suc,
    pc=pc,
    rela=rela

  )
})



output$bin_lower<-shiny::renderUI({
  bin = input_bin()


  table_html <-  paste0('
                        \\theta_0 - \\epsilon = ', bin$location+bin$lbbin,'')

  shiny::tagList(
    # Render the table using MathJax
    shiny::withMathJax(
      shiny::em('$$', table_html, '$$')
    )
  )

})


output$bin_upper<-shiny::renderUI({
 bin = input_bin()

  table_html <-  paste0('
                        \\theta_0 - \\epsilon = ', bin$location+bin$ubbin,'')

  shiny::tagList(
    # Render the table using MathJax
    shiny::withMathJax(
      shiny::em('$$', table_html, '$$')
    )
  )


})




shiny::observeEvent(input$runbin, {
  bin = input_bin()

  dat <- tryCatch({
    suppressWarnings(switch(bin$interval,
                "1" = {bin_table(bin$D,bin$target,bin$h0,bin$alpha,bin$beta,bin$location,
                  bin$scale,bin$model,bin$hypothesis,
                  bin$alpha_d,bin$beta_d,bin$location_d,bin$scale_d,
                  bin$model_d,bin$de_an_prior,bin$N, bin$mode_bf,bin$FP,bin$direct)},
                "2" = {
                  bin_e_table(bin$D,bin$target,bin$h0,bin$alpha,bin$beta,bin$location,
                                 bin$scale,bin$model,bin$hypothesis,
                                 bin$alpha_d,bin$beta_d,bin$location_d,bin$scale_d,
                                 bin$model_d,bin$de_an_prior,bin$N, bin$mode_bf,bin$FP,bin$e,bin$direct)


                  }))}, error = function(e) {
                    "Error"
                  })

  output$result_bin <-  shiny::renderText({
    paste("# Function to be used in R", show_bin_code(bin), sep = "\n")
  })
  output$prior_bin <- shiny::renderPlot({
    switch(bin$interval,
           "1" = {bin_prior_plot(bin$h0,bin$alpha,bin$beta,bin$location,bin$scale,bin$model,
                                 bin$alpha_d,bin$beta_d,bin$location_d,
                                 bin$scale_d,bin$model_d,bin$hypothesis,
                                 bin$de_an_prior)},
           "2" = bin_e_prior_plot (bin$h0,bin$alpha,bin$beta,bin$location,bin$scale,
                                  bin$model,bin$alpha_d,bin$beta_d,bin$location_d,
                                  bin$scale_d,bin$model_d,
                                  bin$hypothesis,bin$de_an_prior,bin$e)
           )



  })

  output$resultbin <- shiny::renderUI({
    if (identical(dat, "Error")){
      table_html <- shiny::span("\\(\\text{Note: Error when the required N > 10,000}\\)", style = "color: red;")
    }else{
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('$$', '
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    p\\text{(BF}_{10} > ', bin$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{01} > ', bin$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[3], 3), nsmall = 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', bin$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[2], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{10} > ', bin$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[4], 3), nsmall = 3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N} & ', dat[5], ' \\\\
    \\end{array}
  ', '$$')
    }
    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em(table_html)
      )
    )
  })
  # Define reactive containers OUTSIDE the if blocks
  pc_bin   <- shiny::reactiveVal(NULL)
  rela_bin <- shiny::reactiveVal(NULL)


  # ===================================================
  #               POWER CURVE (bin)
  # ===================================================

  if (isTRUE(bin$pc)) {

    output$plot_power_bin_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Power Curve}$$"))
      )
    })

    output$plot_power_bin <- shiny::renderPlot({

      suppressWarnings(
        switch(
          bin$interval,

          "1" = Power_bin(
            bin$D, bin$h0, bin$alpha, bin$beta, bin$location, bin$scale,
            bin$model, bin$hypothesis,
            bin$alpha_d, bin$beta_d, bin$location_d,
            bin$scale_d, bin$model_d, bin$de_an_prior,
            dat[1,5]
          ),

          "2" = Power_e_bin(
            bin$D, bin$h0, bin$alpha, bin$beta, bin$location, bin$scale,
            bin$model, bin$hypothesis,
            bin$alpha_d, bin$beta_d, bin$location_d,
            bin$scale_d, bin$model_d, bin$de_an_prior,
            dat[1,5], bin$e
          )
        )
      )

      pc_bin(grDevices::recordPlot())
    })

  } else {

    pc_bin(NULL)
    output$plot_power_bin_text <- shiny::renderUI(NULL)
    output$plot_power_bin      <- shiny::renderPlot(NULL)
  }


  # ===================================================
  #           RELATIONSHIP BETWEEN BF & DATA (bin)
  # ===================================================

  if (isTRUE(bin$rela)) {

    output$plot_rel_bin_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Relationship between BF and data}$$"))
      )
    })

    output$plot_rel_bin <- shiny::renderPlot({

      suppressWarnings(
        switch(
          bin$interval,

          "1" = bin_bf10(
            bin$D, dat[1,5], bin$alpha, bin$beta, bin$location, bin$scale,
            bin$model, bin$hypothesis
          ),

          "2" = bin_e_bf10(
            bin$D, dat[1,5], bin$alpha, bin$beta, bin$location, bin$scale,
            bin$model, bin$hypothesis, bin$e
          )
        )
      )

      rela_bin(grDevices::recordPlot())
    })

  } else {

    rela_bin(NULL)
    output$plot_rel_bin_text <- shiny::renderUI(NULL)
    output$plot_rel_bin      <- shiny::renderPlot(NULL)
  }

  output$export_bin <- shiny::downloadHandler(
    filename = function() {
      "BayesPower-report.html"
    },
    content = function(file) {
      template_path <- system.file("report_templates", "report_bin.Rmd", package = "BayesPower")

      tempReport <- file.path(tempdir(), "report_bin.Rmd")
      file.copy(template_path, tempReport, overwrite = TRUE)

      rmarkdown::render(
        input = tempReport,output_format ="html_document",
        output_file = file,
        params = list(bin = bin, dat = dat,pc_bin=pc_bin(),rela_bin=rela_bin()),  # ✅ pass to `params`
        envir = new.env(parent = globalenv())  # environment still required
      )
    }
  )




})


shiny::observeEvent(input$calbin, {
  bin = input_bin()
  BF10 <- switch(bin$interval,
                 "1" = bin_BF(bin$Suc,bin$N,bin$alpha,bin$beta,bin$location,bin$scale,bin$model,bin$hypothesis),
                 "2" = bin_e_BF(bin$Suc,bin$N,bin$alpha,bin$beta,bin$location,bin$scale,bin$model,bin$hypothesis,bin$e))

  output$BFbin <- shiny::renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    N = ', bin$N, ', x = ', bin$Suc, '; \\textit{BF}_{10} = ', round(BF10, 4), ', \\textit{BF}_{01} = ',round(1/BF10, 4),'
')

    output$result_bin <- shiny::renderText({
      args <- list(
        x = bin$Suc,
        n = bin$N,
        h0 = bin$location,
        model = bin$model,
        hypothesis = bin$hypothesis
      )

      # Add model-specific parameters
      if (bin$model == "beta") {
        args$alpha <- bin$alpha
        args$beta  <- bin$beta
      } else if (bin$model == "Moment") {
        args$scale <- bin$scale
      }

      # Include e only if interval != 1
      if (!is.null(bin$interval) && bin$interval != 1) {
        args$e <- bin$e
      }

      fmt_val <- function(x) {
        if (is.numeric(x) && length(x) == 1) return(as.character(x))
        if (is.numeric(x) && length(x) > 1) return(paste(x, collapse = ", "))
        if (is.character(x)) return(shQuote(x))
        return(as.character(x))
      }

      # Build string with each argument on a new line
      arg_strings <- sapply(names(args), function(arg) {

        val <- args[[arg]]
        arg_print <- arg  # default printed name

        ## model → prior_analysis
        if (arg == "model") arg_print <- "prior_analysis"

        ## e → ROPE
        if (arg == "e") arg_print <- "ROPE"

        ## hypothesis → alternative
        if (arg == "hypothesis") {

          arg_print <- "alternative"

          val <- switch(val,
                        "<"  = "less",
                        "!=" = "two.sided",
                        ">"  = "greater",
                        stop("Invalid hypothesis")
          )
        }

        # Format e as c(...) if vector
        if (arg == "e" && length(val) > 1) {
          sprintf("  %s = c(%s)", arg_print, paste(fmt_val(val), collapse = ", "))
        } else {
          sprintf("  %s = %s", arg_print, fmt_val(val))
        }
      })

      call_string <- paste0(
        "# Function to be used in R\n",
        "BF10.bin.test(\n",
        paste(arg_strings, collapse = ",\n"),
        "\n)"
      )

      call_string

    })




    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em('$$', table_html, '$$')
      )
    )
  })

})
}




# ---- Server_f.r ----

server_f<- function(input, output, session) {
input_f <- shiny::reactive({
  mode_bf <- switch(input$Modef,
                    "1" = 1,
                    "2" = 0,
                    "3" = 0)
  direct<- switch(input$f_direct,
                  "1" = "h1",
                  "0" = "h0")
 anovareg <- input$ANOREG
 reduced_model <-input$redf
 f1 <-input$f1
 f2 <-input$f2
  if (input$ANOREG == 2){
    p <- input$pf
    k <- input$kf
  }else{
    p <-switch(input$redf,
               "1" = 1,
               "2" = input$f1-1+1,
               "3" = input$f1-1 +input$f2-1 +1
    )
    full_model <-switch(input$redf,
                        "1"=input$full1,
                        "2"=input$full2,
                        "3"=input$full3)
    k <-switch(full_model,
               "2" = input$f1-1+1,
               "3" = input$f1-1 +input$f2-1 +1,
               "4" = input$f1-1 +input$f2-1 +1 + (input$f1-1)*(input$f2-1)
    )
  }



  inter <- input$h0f

  e <- switch(inter,
              "1" = input$epsilinff,
              "2" = input$epsilinff)


  model <- switch(input$modelf,
                  "1" = "effectsize",
                  "2" = "Moment")

  rscale <- input$rf
  f_m <- sqrt(input$fsdf)
  dff <- input$dff
  de_an_prior <- switch(input$priorf,
                        "1" = 1,
                        "2" = 0)

  model_d <- switch(input$modelfd,
                    "1" = "effectsize",
                    "2" = "Moment",
                    "3" = "Point")

  rscale_d <- input$rfd
  f_m_d <- sqrt(input$fsdfd)

  if ( input$modelfd == "3"){
    f_m_d <-sqrt(input$lfd)
  }

  dff_d <- input$dffd
  target <- input$powerf
  alpha <- input$alphaf
  N <- input$nf
  D <- input$bff
  fval <- input$fval
  df1 <- input$df1f
  df2 <- input$df2f
  q = k -p
  pc   <- "1" %in% input$o_plot_f
  rela <- "2" %in% input$o_plot_f
  if (mode_bf!=1){
    pc=rela=F
  }
  # Add all variables to the final list


  if (input$ANOREG == 1){
    list(
    mode_bf = mode_bf,
    direct = direct,
    p = p,
    k = k,
    q = q,
    inter=inter,
    e=e,
    model=model,
    rscale=rscale,
    f_m=f_m,
    dff=dff,
    de_an_prior=de_an_prior,
    model_d=model_d,
    rscale_d=rscale_d,
    f_m_d=f_m_d,
    dff_d=dff_d,
    target=target,
    alpha=alpha,
    N=N,
    D=D,
    fval = fval,
    df1=df1,
    df2=df2,
    pc = pc,
    rela = rela,
    anovareg=anovareg,
    full_model = full_model,
    reduced_model=reduced_model,
    f1 =f1,f2=f2
  )}else{
    list(
    mode_bf = mode_bf,
    direct=direct,
    p = p,
    k = k,
    q = q,
    inter=inter,
    e=e,
    model=model,
    rscale=rscale,
    f_m=f_m,
    dff=dff,
    de_an_prior=de_an_prior,
    model_d=model_d,
    rscale_d=rscale_d,
    f_m_d=f_m_d,
    dff_d=dff_d,
    target=target,
    alpha=alpha,
    N=N,
    D=D,
    fval = fval,
    df1=df1,
    df2=df2,
    pc = pc,
    rela = rela,
    anovareg=anovareg
  )}

})


output$prior_suggest <- shiny::renderUI({
  ff = input_f()
  if (ff$model == "effectsize"){

  table_html <- paste0('
\\textit{df} = ', 3, ', \\textit{r} = \\sqrt{\\frac{df - 2}{dfq}} \\times f = ',round(sqrt((3 - 2) / 3*ff$q) * sqrt(ff$f_m),2),'
')
  }else{

    table_html <- paste0('
\\textit{df = 5+(q-1)} = ', 5+ff$q-1,'
')


}
  # Render the table using MathJax
  shiny::tagList(
    # Render the table using MathJax
    shiny::withMathJax(
      shiny::em('$$', table_html, '$$')
    )
  )
})


shiny::observeEvent(input$runf, {
  ff = input_f()

  output$result_f <-  shiny::renderText({
    paste("# Function to be used in R", show_f_code(ff), sep = "\n")
  })


  dat = tryCatch({ switch(ff$inter,
               "1" = f_table(ff$D,ff$target,ff$p,ff$k,ff$dff,ff$rscale,ff$f_m,ff$model,
                ff$dff_d,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior,ff$N, ff$mode_bf,ff$alpha ,ff$direct),
               "2" = fe_table(ff$D,ff$target,ff$p,ff$k,ff$dff,ff$rscale,ff$f_m,ff$model,
                              ff$dff_d,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior,ff$N, ff$mode_bf,ff$e ,ff$alpha,ff$direct))
  }, error = function(e) {
    "Error"
  })

  output$priorff <- shiny::renderPlot({

    switch(ff$inter,
           "1" =prior_plot_f(ff$q,ff$dff,ff$rscale,ff$f_m,ff$model,ff$dff_d
                 ,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior),
           "2" = prior_plot_fe(ff$q,ff$dff,ff$rscale,ff$f_m,ff$model,ff$dff_d
                               ,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior,ff$e))

  })

  output$resultf <- shiny::renderUI({
    if (identical(dat, "Error")){
      table_html <- shiny::span("\\(\\text{Error: the required N > 5,000} \\)", style = "color: red;")

    }else{
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('$$','
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    p\\text{(BF}_{10} > ', ff$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{01} > ', ff$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[3], 3), nsmall = 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', ff$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[2], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{10} > ', ff$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[4], 3), nsmall = 3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N} & ', dat[5], ' \\\\
    \\end{array}
  ', '$$')
    }
    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em(table_html)
      )
    )
  })
  # Define reactive containers OUTSIDE the if blocks
  pc_f   <- shiny::reactiveVal(NULL)
  rela_f <- shiny::reactiveVal(NULL)


  # ===================================================
  #               POWER CURVE (F)
  # ===================================================

  if (isTRUE(ff$pc)) {

    output$plot_power_f_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$${\\text{Power Curve}}$$"))
      )
    })

    output$plot_power_f <- shiny::renderPlot({

      suppressWarnings(
        switch(
          ff$inter,

          "1" = Power_f(
            ff$D, ff$k, ff$p, ff$dff, ff$rscale,
            ff$f_m, ff$model,
            ff$k_d, ff$p_d, ff$dff_d, ff$rscale_d, ff$f_m_d, ff$model_d,
            ff$de_an_prior,
            dat[1,5]
          ),

          "2" = Power_fe(
            ff$D, ff$k, ff$p, ff$dff, ff$rscale,
            ff$f_m, ff$model,
            ff$k_d, ff$p_d, ff$dff_d, ff$rscale_d, ff$f_m_d, ff$model_d,
            ff$de_an_prior,
            dat[1,5],
            ff$e
          )
        )
      )

      pc_f(grDevices::recordPlot())
    })

  } else {

    pc_f(NULL)
    output$plot_power_f_text <- shiny::renderUI(NULL)
    output$plot_power_f      <- shiny::renderPlot(NULL)
  }



  # ===================================================
  #           RELATIONSHIP BETWEEN BF & DATA (F)
  # ===================================================

  if (isTRUE(ff$rela)) {

    output$plot_rel_f_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$${\\text{Relationship between BF and data}}$$"))
      )
    })

    output$plot_rel_f <- shiny::renderPlot({

      suppressWarnings(
        switch(
          ff$inter,

          "1" = bf10_f(
            ff$D, dat[1,5], ff$k, ff$p, ff$dff, ff$rscale, ff$f_m, ff$model
          ),

          "2" = bf10_fe(
            ff$D, dat[1,5], ff$k, ff$p, ff$dff, ff$rscale, ff$f_m, ff$model, ff$e
          )
        )
      )

      rela_f(grDevices::recordPlot())
    })

  } else {

    rela_f(NULL)
    output$plot_rel_f_text <- shiny::renderUI(NULL)
    output$plot_rel_f      <- shiny::renderPlot(NULL)
  }


  output$export_f <- shiny::downloadHandler(
    filename = function() {
      "BayesPower-report.html"
    },
    content = function(file) {
      template_path <- system.file("report_templates", "report_f.Rmd", package = "BayesPower")

      tempReport <- file.path(tempdir(), "report_f.Rmd")
      file.copy(template_path, tempReport, overwrite = TRUE)

      rmarkdown::render(
        input = tempReport,output_format ="html_document",
        output_file = file,
        params = list(ff = ff, dat = dat,pc_f=pc_f(),rela_f=rela_f()),  # ✅ pass to `params`
        envir = new.env(parent = globalenv())  # environment still required
      )
    }
  )




})






shiny::observeEvent(input$calf, {
  ff = input_f()
  m = ff$df2+ff$df1
  output$priorff <- shiny::renderPlot({

    switch(ff$inter,
           "1" =prior_plot_f(ff$q,ff$dff,ff$rscale,ff$f_m,ff$model,ff$dff_d
                             ,ff$rscale_d,ff$f_m_d,ff$model_d,1),
           "2" = prior_plot_fe(ff$q,ff$dff,ff$rscale,ff$f_m,ff$model,ff$dff_d
                               ,ff$rscale_d,ff$f_m_d,ff$model_d,1,ff$e))

  })
  BF10 <- if (ff$inter==1) F_BF(ff$fval,ff$df1,m,ff$dff,ff$rscale,ff$f_m,ff$model) else
    Fe_BF(ff$fval,ff$df1,m,ff$dff,ff$rscale,ff$f_m,ff$model,ff$e)

  output$result_f <- shiny::renderText({

    args <- list(
      fval = ff$fval,
      df1 = ff$df1,
      df2 = ff$df2,
      dff = ff$dff,
      rscale = ff$rscale,
      f_m = ff$f_m,
      model = ff$model
    )

    if (ff$inter!=1) {
      args$e <- ff$e
    }

    # Build string with each argument on a new line
    arg_strings <- sapply(names(args), function(arg) {

      val <- args[[arg]]
      arg_print <- arg

      ## model > prior_analysis
      if (arg == "model") {
        arg_print <- "prior_analysis"
      }

      ## e > ROPE
      if (arg == "e") {
        arg_print <- "ROPE"
      }

      sprintf("  %s = %s", arg_print, fmt_val(val))
    })


    call_string <- paste0(
      "# Function to be used in R\n",
      "BF10.f.test(\n",
      paste(arg_strings, collapse = ",\n"),
      "\n)"
    )

    call_string
  })


  output$BFcalf <- shiny::renderUI({
    # Create the LaTeX formatted strings for the table
    output$BFcalf <- shiny::renderUI({
      # Create the LaTeX formatted string with proper escaping
      table_latex <- paste0(
        "$$ \\textit{F}(", ff$df1, ",", ff$df2,
        ") = ", round(ff$fval, 3),
        ",\\ \\textit{BF}_{10} = ", round(BF10, 4),", \\textit{BF}_{01} =" ,round(1/BF10, 4)," $$"
      )

      shiny::tagList(
        shiny::withMathJax(
          shiny::em(table_latex)
        )
      )
    })

  })

})
}

# ---- Server_p2.r ----

server_p2<- function(input, output, session) {
input_p2 <- shiny::reactive({


  mode_bf <- switch(input$Modep2,
                    "1" = 1,
                    "2" = 0,
                    "3" = 0)# mode
 direct <- switch(input$p2_direct,
                  "1" = "h1",
                  "0" = "h0")
  a0 <- input$alpha0
  b0 <- input$beta0

  a1 <- input$alpha1
  b1 <- input$beta1

  a2 <- input$alpha2
  b2 <- input$beta2

  a1d <- input$alpha1d
  b1d <- input$beta1d

  a2d <- input$alpha2d
  b2d <- input$beta2d

  dp1 <- input$location1d
  dp2 <- input$location2d


  model_p1 <-switch(input$model_p1,
                    "1" = "Point",
                    "2" = "beta")
  model_p2 <-switch(input$model_p2,
                    "1" = "Point",
                    "2" = "beta")

  if (input$priorp2 == 1){
    model_p1 = model_p2 = "same"
  }
  de_an_prior<-input$priorp2
  D <- input$bp2

  n1 <- input$n1p2
  n2 <- input$n2p2

  x1 <- input$x1p2
  x2 <- input$x2p2

  target <- input$powerp2
  pc   <- "1" %in% input$o_plot_p2
  rela <- "2" %in% input$o_plot_p2

  if (mode_bf!=1){
    pc=rela=F
  }
 ############

  list(
    mode_bf = mode_bf,
    direct = direct,
    a0 = a0,
    b0 = b0,
    a1 = a1,
    b1 = b1,
    a2 = a2,
    b2 = b2,
    a1d =  a1d,
    b1d = b1d,
    a2d = a2d,
    b2d = b2d,
    model1 = model_p1,
    model2 = model_p2,
    dp1 = dp1,
    dp2 = dp2,
    D = D,
    n1 = round(n1),
    n2 = round(n2),
    k1 = round(x1),
    k2 = round(x2),
    target = target,
    r=1,
    pc=pc,
    rela=rela,
    de_an_prior=de_an_prior

  )
})



shiny::observeEvent(input$runp2, {

  p2 <- input_p2()

  # Compute data safely

  dat <- tryCatch({
    pro_table_p2(
      p2$D, p2$target, p2$a0, p2$b0,
      p2$a1, p2$b1, p2$a2, p2$b2, p2$r,
      p2$model1, p2$a1d, p2$b1d, p2$dp1,
      p2$model2, p2$a2d, p2$b2d, p2$dp2,
      p2$mode_bf, p2$n1, p2$n2, p2$direct
    )
  }, error = function(e) "Error")


  # If dat is NULL, skip processing
  shiny::req(!is.null(dat))
  table <- dat[[1]]

  # Render function code
  output$result_p2 <- shiny::renderText({
    paste("# Function to be used in R", show_props_code(p2), sep = "\n")
  })

  # Render priors
  output$prior_p0 <- shiny::renderPlot({
    p2_prior_plot(p2$a0, p2$b0, 1, 1, 0, "same", 0)
  })
  output$prior_p1 <- shiny::renderPlot({
    p2_prior_plot(p2$a1, p2$b1, p2$a1d, p2$b1d, p2$dp1, p2$model1, 1)
  })
  output$prior_p2 <- shiny::renderPlot({
    p2_prior_plot(p2$a2, p2$b2, p2$a2d, p2$b2d, p2$dp2, p2$model2, 2)
  })

  # Render results table
  output$resultp2 <- shiny::renderUI({
    if (identical(dat, "Error")) {
      shiny::withMathJax(
        shiny::span("\\(\\text{Note: Error when required } N > 5,000\\)", style = "color: red;")
      )
    } else {
      table <- dat[[1]]
      table_html <- paste0('$$',
                           '\\begin{array}{l c}
      \\textbf{Probability of Compelling Evidence} & \\\\
      \\hline
      p\\text{(BF}_{10} > ', p2$D, '\\, | \\, \\mathcal{H}_1) & ', format(round(table[1,1],3), nsmall = 3), ' \\\\
      p\\text{(BF}_{01} > ', p2$D, '\\, | \\, \\mathcal{H}_0) & ', format(round(table[1,3],3), nsmall = 3), ' \\\\
      \\textbf{Probability of Misleading Evidence} & \\\\
      \\hline
      p\\text{(BF}_{01} > ', p2$D, '\\, | \\, \\mathcal{H}_1) & ', format(round(table[1,2],3), nsmall = 3), ' \\\\
      p\\text{(BF}_{10} > ', p2$D, '\\, | \\, \\mathcal{H}_0) & ', format(round(table[1,4],3), nsmall = 3), ' \\\\
      \\textbf{Required Sample Size} & \\\\
      \\hline
      \\text{N}_1 & ', table[1,5], ' \\\\
      \\text{N}_2 & ', table[1,6], ' \\\\
      \\end{array}
      $$'
      )
      shiny::withMathJax(shiny::em(table_html))
    }
  })

  # Define reactive containers OUTSIDE the if blocks
  pc_p2   <- shiny::reactiveVal(NULL)
  rela_p2 <- shiny::reactiveVal(NULL)

  # ===================================================
  #               POWER CURVE
  # ===================================================

  if (isTRUE(p2$pc)) {

    output$plot_power_p2_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Power Curve}$$"))
      )
    })

    output$plot_power_p2 <- shiny::renderPlot({

      suppressWarnings(
        Power_p2(
          p2$D, table[1,5], p2$a0, p2$b0, p2$a1, p2$b1, p2$a2,
          p2$b2, table[1,6] / table[1,5], p2$model1, p2$a1d, p2$b1d, p2$dp1,
          p2$model2, p2$a2d, p2$b2d, p2$dp2
        )
      )

      pc_p2(grDevices::recordPlot())
    })

  } else {

    pc_p2(NULL)
    output$plot_power_p2_text <- shiny::renderUI(NULL)
    output$plot_power_p2      <- shiny::renderPlot(NULL)
  }

  # ===================================================
  #           RELATIONSHIP BETWEEN BF & DATA
  # ===================================================

  if (isTRUE(p2$rela)) {

    output$plot_rel_p2_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Relationship between BF and data}$$"))
      )
    })

    output$plot_rel_p2 <- shiny::renderPlot({
      # Explicitly assign and print the ggplot
      plt <- heatmap_p2(dat[[2]], p2$D)
      print(plt)

      # Save the plot for later
      rela_p2(grDevices::recordPlot())
    })

  } else {

    rela_p2(NULL)
    output$plot_rel_p2_text <- shiny::renderUI(NULL)
    output$plot_rel_p2      <- shiny::renderPlot(NULL)
  }

  # Download handler
  output$export_p2 <- shiny::downloadHandler(
    filename = function() "BayesPower-report.html",
    content = function(file) {
      template_path <- system.file("report_templates", "report_2p.Rmd", package = "BayesPower")
      tempReport <- file.path(tempdir(), "report_2p.Rmd")
      file.copy(template_path, tempReport, overwrite = TRUE)
      rmarkdown::render(
        input = tempReport,
        output_format = "html_document",
        output_file = file,
        params = list(p2 = p2, dat = dat, pc_p2 = pc_p2(), rela_p2 = rela_p2()),
        envir = new.env(parent = globalenv())
      )
    }
  )

})



shiny::observeEvent(input$calp2, {
  p2 = input_p2()
  BF10 <- BF10_p2(p2$a0, p2$b0, p2$a1, p2$b1, p2$a2, p2$b2,p2$n1,p2$n2,p2$k1,p2$k2)

  output$BFp2 <- shiny::renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0(
      'n_1 = ', p2$n1, ', ',
      'n_2 = ', p2$n2, ', ',
      'x_1 = ', p2$k1, ', ',
      'x_2 = ', p2$k2, ' \\\\ ',
      '\\textit{BF}_{10} = ', round(BF10, 4),
      ', \\textit{BF}_{01} = ', round(1/BF10, 4)
    )


    output$result_p2 <- shiny::renderText({


      args <- list(
        a0 = p2$a0,
        b0 = p2$b0,
        a1 = p2$a1,
        b1 = p2$b1,
        a2 = p2$a2,
        b2 = p2$b2,
        n1 = p2$n1,
        n2 = p2$n2,
        x1 = p2$k1,
        x2 = p2$k2
      )
      fmt_val <- function(x) {
        if (is.numeric(x) && length(x) == 1) return(as.character(x))
        if (is.numeric(x) && length(x) > 1) return(paste(x, collapse = ", "))
        if (is.character(x)) return(shQuote(x))
        return(as.character(x))
      }
      # Build string with each argument on a new line
      arg_strings <- sapply(names(args), function(nm) {
        sprintf("  %s = %s", nm, fmt_val(args[[nm]]))
      })

      call_string <- paste0(
        "# Function to be used in R\n",
        "BF10.props(\n",
        paste(arg_strings, collapse = ",\n"),
        "\n)"
      )

      call_string
    })




    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em('$$', table_html, '$$')
      )
    )
  })

})
}

# ---- Server_r.r ----

server_r<- function(input, output, session) {
input_r <- shiny::reactive({
  mode_bf <- switch(input$Moder,
                    "1" = 1,
                    "2" = 0,
                    "3" = 0)# mode
  direct <- switch(input$r_direct,
                   "1" = "h1",
                   "0" = "h0")
  interval <- input$h0r # point null or interval
  h0 <- input$h0pho
  lbre <- input$lbre
  ubre <- input$ubre

  if ((h0+lbre)<(-1)){
    lbre = lbre+-1-(h0+lbre)

  }

  if ((h0+ubre)>(+1)){
    ubre = ubre+1-(h0+ubre)

  }



  e <- switch(input$h1re,        # bound for interval test
              "1" = c(lbre, ubre),
              "2" = ubre,
              "3" = lbre)
  inter <- switch(interval,
                  "1" = input$h1r,
                  "2" = input$h1re)


  hypothesis <- switch(interval,
                       "1" =   switch(input$h1r,        # direction of the test
                                      "1" = "!=",
                                      "2" =  ">",
                                      "3" =  "<"),
                       "2" = switch(input$h1re,        # direction of the test
                                    "1" = "!=",
                                    "2" =  ">",
                                    "3" =  "<"))





  model <- switch(input$modelr,
                  "1" = "d_beta",
                  "2" = "beta",
                  "3" = "Moment")
  k <- input$kr
  scale <- input$sr
  alpha <- input$ralpha
  beta <- input$rbeta
  de_an_prior <- switch(input$priorr,
                        "1" = 1,
                        "2" = 0)
  model_d <- switch(input$modelrd,
                    "1" = "d_beta",
                    "2" = "beta",
                    "3" = "Moment",
                    "4" = "Point")
  location_d <- input$h0phod
  k_d <- input$rkd
  scale_d <- input$rsd
  alpha_d <- input$ralphad
  beta_d<- input$rbetad
  target <- input$powerr
  FP <- input$alphapr
  D <- input$br
  N <-  switch(input$Moder,
               "1" = 2,
               "2" = input$nr,
               "3" = input$rdf)
  rval <- input$rval
  pc   <- 1 %in% input$o_plot_r
  rela <- 2 %in% input$o_plot_r
  if (mode_bf!=1){
    pc=rela=F
  }
  ###########
  location <- h0
  dff <- 1

  dff_d <- 1



 ############


  # Add all variables to the final list
  list(
    mode_bf = mode_bf,
    direct = direct,
    interval = interval,
    e = e,
    lbre = lbre,
    ubre = ubre,
    inter = inter,
    hypothesis = hypothesis,
    h0 = h0,
    model = model,
    k =k,
    scale = scale,
    alpha = alpha,
    beta = beta,
    de_an_prior = de_an_prior,
    model_d =model_d,
    location_d = location_d,
    k_d = k_d,
    scale_d = scale_d,
    alpha_d = alpha_d,
    beta_d = beta_d,
    target = target,
    FP = FP,
    D = D,
    N = N,
    rval = rval,
    location = location,
    dff = dff ,
    dff_d = dff_d,
    pc = pc,
    rela = rela

  )
})

output$r_lower<-shiny::renderUI({
 rr = input_r()


  table_html <-  paste0('
                        \\rho_0 - \\epsilon = ', rr$h0+rr$lbre,'')

  shiny::tagList(
    # Render the table using MathJax
    shiny::withMathJax(
      shiny::em('$$', table_html, '$$')
    )
  )

})


output$r_upper<-shiny::renderUI({
  rr = input_r()

  table_html <-  paste0('
                        \\rho_0 - \\epsilon = ', rr$h0+rr$ubre,'')

  shiny::tagList(
    # Render the table using MathJax
    shiny::withMathJax(
      shiny::em('$$', table_html, '$$')
    )
  )


})

shiny::observeEvent(input$runr, {
  rr = input_r()

  dat <- tryCatch({
    switch(rr$interval,

    "1" = r_table(rr$D,rr$target,rr$model,rr$k,
            rr$alpha, rr$beta,rr$h0,rr$location,
            rr$scale,rr$dff, rr$hypothesis ,rr$model_d,
            rr$location_d,rr$k_d, rr$alpha_d, rr$beta_d,
            rr$scale_d,rr$dff_d,rr$de_an_prior,rr$N,
            rr$mode_bf,rr$FP,rr$direct ),
    "2" = re_table(rr$D,rr$target,rr$model,rr$k,
                  rr$alpha, rr$beta,rr$h0,rr$location,
                  rr$scale,rr$dff, rr$hypothesis ,rr$model_d,
                  rr$location_d,rr$k_d, rr$alpha_d, rr$beta_d,
                  rr$scale_d,rr$dff_d,rr$de_an_prior,rr$N,
                  rr$mode_bf,rr$FP,rr$e,rr$direct ))
  }, error = function(e) {
    "Error"
  })

  output$result_r <-  shiny::renderText({
    paste("# Function to be used in R", show_cor_code(rr), sep = "\n")
  })
  output$prior_r <- shiny::renderPlot({

    switch(rr$interval,
           "1" = r_prior_plot(rr$k, rr$alpha, rr$beta,
                              rr$h0,rr$location,rr$scale,
                              rr$dff,rr$model,rr$de_an_prior,
                              rr$k_d, rr$alpha_d, rr$beta_d,
                              rr$location_d,rr$scale_d,rr$dff_d,
                              rr$model_d,rr$hypothesis),
           "2" = re_prior_plot(rr$k, rr$alpha, rr$beta,
                              rr$h0,rr$location,rr$scale,
                              rr$dff,rr$model,rr$de_an_prior,
                              rr$k_d, rr$alpha_d, rr$beta_d,
                              rr$location_d,rr$scale_d,rr$dff_d,
                              rr$model_d,rr$hypothesis,rr$e))



  })

  output$resultr <- shiny::renderUI({
    if (identical(dat, "Error")){
      table_html <- shiny::span("\\(\\text{Note: Potential Error when the required N > 5,000} \\)", style = "color: red;")

    }else{
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('$$','
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    p\\text{(BF}_{10} > ', rr$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{01} > ', rr$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[3], 3), nsmall = 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', rr$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[2], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{10} > ', rr$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[4], 3), nsmall = 3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N} & ', dat[5], ' \\\\
    \\end{array}
  ','$$')
    }
    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em(table_html)
      )
    )
  })


  # Define reactive containers OUTSIDE the if blocks
  pc_r   <- shiny::reactiveVal(NULL)
  rela_r <- shiny::reactiveVal(NULL)


  # ===================================================
  #                POWER CURVE (r)
  # ===================================================

  if (isTRUE(rr$pc)) {

    output$plot_power_r_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Power Curve}$$"))
      )
    })

    output$plot_power_r <- shiny::renderPlot({

      suppressWarnings(
        switch(
          rr$interval,

          "1" = Power_r(
            rr$D, rr$k, rr$alpha, rr$beta, rr$h0, rr$hypothesis,
            rr$location, rr$scale, rr$dff, rr$model,
            rr$k_d, rr$alpha_d, rr$beta_d, rr$location_d, rr$scale_d,
            rr$dff_d, rr$model_d, rr$de_an_prior, dat[1,5]
          ),

          "2" = Power_re(
            rr$D, rr$k, rr$alpha, rr$beta, rr$h0, rr$hypothesis,
            rr$location, rr$scale, rr$dff, rr$model,
            rr$k_d, rr$alpha_d, rr$beta_d, rr$location_d, rr$scale_d,
            rr$dff_d, rr$model_d, rr$de_an_prior, dat[1,5], rr$e
          )
        )
      )

      pc_r(grDevices::recordPlot())
    })

  } else {

    pc_r(NULL)
    output$plot_power_r_text <- shiny::renderUI(NULL)
    output$plot_power_r      <- shiny::renderPlot(NULL)
  }



  # ===================================================
  #                RELATIONSHIP (r)
  # ===================================================

  if (isTRUE(rr$rela)) {

    output$plot_rel_r_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Relationship between BF and data}$$"))
      )
    })

    output$plot_rel_r <- shiny::renderPlot({

      suppressWarnings(
        switch(
          rr$interval,

          "1" = r_bf10_p(
            rr$D, dat[1,5], rr$k, rr$alpha, rr$beta, rr$h0,
            rr$hypothesis, rr$location, rr$scale, rr$dff, rr$model
          ),

          "2" = re_bf10_p(
            rr$D, dat[1,5], rr$k, rr$alpha, rr$beta, rr$h0,
            rr$hypothesis, rr$location, rr$scale, rr$dff, rr$model, rr$e
          )
        )
      )

      rela_r(grDevices::recordPlot())
    })

  } else {

    rela_r(NULL)
    output$plot_rel_r_text <- shiny::renderUI(NULL)
    output$plot_rel_r      <- shiny::renderPlot(NULL)
  }

  output$export_r <- shiny::downloadHandler(
    filename = function() {
      "BayesPower-report.html"
    },
    content = function(file) {
      template_path <- system.file("report_templates", "report_r.Rmd", package = "BayesPower")

      tempReport <- file.path(tempdir(), "report_r.Rmd")
      file.copy(template_path, tempReport, overwrite = TRUE)

      rmarkdown::render(
        input = tempReport,output_format ="html_document",
        output_file = file,
        params = list(rr = rr, dat = dat,pc_r=pc_r(),rela_r=rela_r()),  # ✅ pass to `params`
        envir = new.env(parent = globalenv())  # environment still required
      )
    }
  )



})

shiny::observeEvent(input$calr, {
  rr = input_r()
  output$prior_r <- shiny::renderPlot({

    switch(rr$interval,
           "1" = r_prior_plot(rr$k, rr$alpha, rr$beta,
                              rr$h0,rr$location,rr$scale,
                              rr$dff,rr$model,1,
                              rr$k_d, rr$alpha_d, rr$beta_d,
                              rr$location_d,rr$scale_d,rr$dff_d,
                              rr$model_d,rr$hypothesis),
           "2" = re_prior_plot(rr$k, rr$alpha, rr$beta,
                               rr$h0,rr$location,rr$scale,
                               rr$dff,rr$model,1,
                               rr$k_d, rr$alpha_d, rr$beta_d,
                               rr$location_d,rr$scale_d,rr$dff_d,
                               rr$model_d,rr$hypothesis,rr$e))



  })
  BF10 <- switch(rr$interval ,
                 "1" = r_BF10(rr$rval,rr$N,rr$k, rr$alpha, rr$beta,rr$h0,rr$hypothesis,rr$location,rr$scale,rr$dff,rr$model),
                 "2" = re_BF10(rr$rval,rr$N,rr$k, rr$alpha, rr$beta,rr$h0,rr$hypothesis,rr$location,rr$scale,rr$dff,rr$model,rr$e))

  output$result_r <- shiny::renderText({

    build_BF10_call <- function(rr) {

      fmt_val <- function(x) {
        if (is.numeric(x) && length(x) == 1) return(as.character(x))
        if (is.numeric(x) && length(x) > 1) return(paste(x, collapse = ", "))
        if (is.character(x)) return(shQuote(x))
        return(as.character(x))
      }

      args <- list(
        r = rr$rval,
        n = rr$N,
        model = rr$model
      )

      # Add model-specific arguments
      if (!is.null(rr$model)) {
        if (rr$model == "d_beta") {
          args$k <- rr$k
        } else if (rr$model == "beta") {
          args$alpha <- rr$alpha
          args$beta  <- rr$beta
        } else if (rr$model == "Moment") {
          args$scale <- rr$scale
        }
      }

      # Common arguments
      args$h0 <- rr$h0
      args$hypothesis <- rr$hypothesis

      # Optional ROPE (only if interval != 1)
      if (!is.null(rr$interval) && rr$interval != 1) {
        args$e <- rr$e
      }

      # Build string with renaming & mapping rules
      arg_strings <- sapply(names(args), function(arg) {

        val <- args[[arg]]
        arg_print <- arg

        ## model > prior_analysis
        if (arg == "model") {
          arg_print <- "prior_analysis"
        }

        ## e > ROPE
        if (arg == "e") {
          arg_print <- "ROPE"
          if (length(val) > 1) {
            return(sprintf(
              "  %s = c(%s)",
              arg_print,
              paste(fmt_val(val), collapse = ", ")
            ))
          }
        }

        ## hypothesis > alternative
        if (arg == "hypothesis") {

          arg_print <- "alternative"

          val <- switch(val,
                        "<"  = "less",
                        "!=" = "two.sided",
                        ">"  = "greater",
                        stop("Invalid hypothesis")
          )

          val <- shQuote(val)
        } else {
          val <- fmt_val(val)
        }

        sprintf("  %s = %s", arg_print, val)
      })

      paste0(
        "# Function to be used in R\n",
        "BF10.cor(\n",
        paste(arg_strings, collapse = ",\n"),
        "\n)"
      )
    }



    build_BF10_call(rr)
  })





  output$BFrv <- shiny::renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    \\textit{r}(n = ', rr$N , ') = ',rr$rval,', \\textit{BF}_{10} = ', round(BF10, 4),", \\textit{BF}_{01} = ",round(1/BF10, 4), '
')


    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em('$$', table_html, '$$')
      )
    )
  })

})
}




# ---- Server_t1.r ----

server_t1<- function(input, output, session) {
input_t1 <- shiny::reactive({
  mode_bf <- switch(input$Modet1,
                    "1" = 1,
                    "2" = 0,
                    "3" = 0)# mode
  N <-  switch(input$Modet1,
               "1" = 2,
               "2" = input$nt1,
               "3" = input$t1df)
  direct <- switch(input$t1_direct,
                   "1" = "h1",
                   "0" = "h0")

  interval <- input$h0t1 # point null or interval

  e <- switch(input$h1t1e,        # bound for interval test
              "1" = c(input$lbt1e, input$ubt1e),
              "2" = input$ubt1e,
              "3" = input$lbt1e)
  inter <- switch(interval,
                  "1" = input$h1t1,
                  "2" = input$h1t1e)

  hypothesis <- switch(interval,
                       "1" =   switch(input$h1t1,        # direction of the test
                                      "1" = "!=",
                                      "2" =  ">",
                                      "3" =  "<"),
                       "2" = switch(input$h1t1e,        # direction of the test
                                    "1" = "!=",
                                    "2" =  ">",
                                    "3" =  "<"))

  model <- switch(input$modelt1,
                  "1" = "t-distribution",
                  "2" = "Normal",
                  "3" = "Moment")

  location <- input$lt1
  scale <- input$st1
  dff <- input$dft1
  de_an_prior <- switch(input$prior,
                        "1" = 1,
                        "2" = 0)
  model_d <- switch(input$modelt1d,
                    "1" = "t-distribution",
                    "2" = "Normal",
                    "3" = "Moment",
                    "4" = "Point")
  location_d <- input$lt1d

  scale_d <- input$st1d
  dff_d <- input$dft1d
  D <- input$bt1
  type <- input$typet1
  target <- input$powert1
  alpha <- input$alphat1

  tval <- input$t1tval
  pc   <- "1" %in% input$o_plot_t1
  rela <- "2" %in% input$o_plot_t1

  if (mode_bf!=1){
    pc=rela=F
  }
  # Add all variables to the final list
  list(
    mode_bf = mode_bf,
    direct = direct,
    interval = interval,
    hypothesis = hypothesis ,
    e = e,
    model = model,
    location = location,
    scale = scale,
    dff = dff,
    de_an_prior = de_an_prior,
    model_d = model_d,
    location_d = location_d,
    scale_d = scale_d,
    dff_d = dff_d,
    type = type,
    D = D,
    target = target,
    alpha = alpha ,
    N = N,
    tval = tval,
    pc = pc,
    rela = rela
  )
})

shiny::observeEvent(input$runt1, {
  x = input_t1()

  dat = tryCatch({suppressWarnings(switch(x$interval, "1" =  t1_Table(x$D,x$target,x$model,x$location,x$scale,x$dff, x$hypothesis,
                                                        x$model_d,x$location_d,x$scale_d,x$dff_d, x$de_an_prior,x$N, x$mode_bf ,
                                                        x$alpha,x$direct),
                                          "2" = t1e_table(x$D,x$target,x$model,x$location,x$scale,x$dff, x$hypothesis,x$e ,
                                                                                 x$model_d,x$scale_d,x$dff_d, x$de_an_prior,x$N,x$mode_bf,x$location_d ,x$alpha,x$direct)))
  }, error = function(e) {
    "Error"
  })

  output$result_t1 <- shiny::renderText({
    paste("# Function to be used in R", show_t1_code(x), sep = "\n")
  })
  output$priort1 <- shiny::renderPlot({
    suppressWarnings(switch(x$interval,
           "1"= t1_prior_plot(
      D = x$D,                  # Access 'D' explicitly
      target = x$target,        # Access 'target' explicitly
      model = x$model,          # Access 'model' explicitly
      location = x$location,    # Access 'location' explicitly
      scale = x$scale,          # Access 'scale' explicitly
      dff = x$dff,              # Access 'dff' explicitly
      hypothesis = x$hypothesis,  # Access 'hypothesis' explicitly
      model_d = x$model_d,        # Access 'model_d' explicitly
      location_d = x$location_d,  # Access 'location_d' explicitly
      scale_d = x$scale_d,        # Access 'scale_d' explicitly
      dff_d = x$dff_d,            # Access 'dff_d' explicitly
      de_an_prior = x$de_an_prior   # Access 'de_an_prior' explicitly
    ),
    "2" = t1e_prior_plot(x$model,
                  x$location,
                   x$scale,
                   x$dff ,
                   x$hypothesis,
                   x$e,
                   x$de_an_prior,
                   x$model_d,
                   x$scale_d,
                   x$dff_d,
                   x$location_d )

  ))

  })



  output$resultt1 <- shiny::renderUI({
    if (identical(dat, "Error")){
      table_html <- shiny::span("\\(\\text{Error: the required N > 10,000}\\)", style = "color: red;")
    }else{
    # Create the LaTeX formatted strings for the table
    table_html <- paste0( '$$','
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    p\\text{(BF}_{10} > ', x$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{01} > ', x$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[3], 3), nsmall = 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', x$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[2], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{10} > ', x$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[4], 3), nsmall = 3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N} & ', dat[5], ' \\\\
    \\end{array}
  ', '$$')
    }
    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em( table_html)
      )
    )
  })
  # --- reactive containers should be defined OUTSIDE the if blocks ----
  pc_t1   <- shiny::reactiveVal(NULL)
  rela_t1 <- shiny::reactiveVal(NULL)


  # ===============================
  #   POWER CURVE SECTION (pc)
  # ===============================

  if (isTRUE(x$pc)) {

    output$plot_power_t1_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Power Curve}$$"))
      )
    })

    output$plot_power_t1 <- shiny::renderPlot({

      plt <-suppressWarnings(
        switch(
          x$interval,
          "1" = Power_t1(
            x$D, x$model, x$location, x$scale, x$dff, x$hypothesis,
            x$model_d, x$location_d, x$scale_d, x$dff_d,
            x$de_an_prior, dat[1,5]
          ),
          "2" = Power_t1e(
            x$D, x$model, x$location, x$scale, x$dff, x$hypothesis,
            x$model_d, x$location_d, x$scale_d, x$dff_d,
            x$de_an_prior, dat[1,5], x$e
          )
        )
      )
      print(plt)
      pc_t1(grDevices::recordPlot())
    })

  } else {

    pc_t1(NULL)
    output$plot_power_t1_text <- shiny::renderUI(NULL)
    output$plot_power_t1      <- shiny::renderPlot(NULL)
  }



  # ===============================
  #   RELATIONSHIP SECTION (rela)
  # ===============================

  if (isTRUE(x$rela)) {

    output$plot_rel_t1_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Relationship between BF and data}$$"))
      )
    })

    output$plot_rel_t1 <- shiny::renderPlot({

      plt<-suppressWarnings(
        switch(
          x$interval,
          "1" =
            bf10_t1(
              D         = x$D,
              df        = dat[1,5],
              target    = x$target,
              model     = x$model,
              location  = x$location,
              scale     = x$scale,
              dff       = x$dff,
              hypothesis = x$hypothesis
            ),
          "2" =
            te1_BF(
              x$D, dat[1,5], x$model, x$location, x$scale, x$dff,
              x$hypothesis, x$e
            )
        )
      )

      rela_t1(grDevices::recordPlot())
    })

  } else {

    rela_t1(NULL)
    output$plot_rel_t1_text <- shiny::renderUI(NULL)
    output$plot_rel_t1      <- shiny::renderPlot(NULL)
  }


  output$export_t1 <- shiny::downloadHandler(
    filename = function() {
      "BayesPower-report.html"
    },
    content = function(file) {

      template_path <- system.file("report_templates", "report_t1.Rmd", package = "BayesPower")

      tempReport <- file.path(tempdir(), "report_t1.Rmd")
      file.copy( template_path, tempReport, overwrite = TRUE)

      rmarkdown::render(
        input = tempReport,output_format ="html_document",
        output_file = file,
        params = list(x = x, dat = dat,pc_t1=pc_t1(),rela_t1=rela_t1()),  # ✅ pass to `params`
        envir = new.env(parent = globalenv())  # environment still required
      )
    }
  )


})

shiny::observeEvent(input$cal1, {
  x = input_t1()

  output$result_t1 <- shiny::renderText({args <- list(
    tval = x$tval,
    df = x$N,
    model = x$model,
    location = x$location,
    scale = x$scale,
    dff = x$dff,
    hypothesis = x$hypothesis
  )

  if (!is.null(x$e) && x$interval != 1) {
    args$e <- x$e
  }

  # Build string with each argument on a new line
  arg_strings <- sapply(names(args), function(arg) {

    # fmt_val defined inside sapply, same style as before
    fmt_val <- function(x) {
      if (is.character(x)) {
        # Wrap character strings in quotes
        return(shQuote(x))
      } else if (is.numeric(x) && length(x) > 1) {
        # Wrap numeric vectors in c(...)
        return(paste0("c(", paste(x, collapse = ", "), ")"))
      } else if (is.numeric(x)) {
        return(as.character(x))
      } else if (is.logical(x)) {
        return(ifelse(x, "TRUE", "FALSE"))
      } else {
        stop("Unsupported type in fmt_val")
      }
    }

    val <- args[[arg]]
    arg_print <- arg

    ## model > prior_analysis
    if (arg == "model") {
      arg_print <- "prior_analysis"
    }

    ## e > ROPE
    if (arg == "e") {
      arg_print <- "ROPE"
    }

    ## hypothesis > alternative
    if (arg == "hypothesis") {
      arg_print <- "alternative"

      val <- switch(val,
                    "<"  = "less",
                    "!=" = "two.sided",
                    ">"  = "greater",
                    stop("Invalid hypothesis")
      )

      val <- shQuote(val)

    } else {
      val <- fmt_val(val)
    }

    sprintf("  %s = %s", arg_print, val)
  })

  call_string <- paste0(
    "# Function to be used in R\n",
    "BF10.ttest.OneSample(\n",
    paste(arg_strings, collapse = ",\n"),
    "\n)"
  )

  call_string
  })





  output$priort1 <- shiny::renderPlot({
    suppressWarnings(switch(x$interval,
                            "1"= t1_prior_plot(
                              D = x$D,                  # Access 'D' explicitly
                              target = x$target,        # Access 'target' explicitly
                              model = x$model,          # Access 'model' explicitly
                              location = x$location,    # Access 'location' explicitly
                              scale = x$scale,          # Access 'scale' explicitly
                              dff = x$dff,              # Access 'dff' explicitly
                              hypothesis = x$hypothesis,  # Access 'hypothesis' explicitly
                              model_d = x$model_d,        # Access 'model_d' explicitly
                              location_d = x$location_d,  # Access 'location_d' explicitly
                              scale_d = x$scale_d,        # Access 'scale_d' explicitly
                              dff_d = x$dff_d,            # Access 'dff_d' explicitly
                              de_an_prior = 1   # Access 'de_an_prior' explicitly
                            ),
                            "2" = t1e_prior_plot(x$model,
                                                 x$location,
                                                 x$scale,
                                                 x$dff ,
                                                 x$hypothesis,
                                                 x$e,
                                                 1,
                                                 x$model_d,
                                                 x$scale_d,
                                                 x$dff_d,
                                                 x$location_d )

    ))

  })
  BF10 <- suppressWarnings(switch(x$interval,
                 "1" = t1_BF10(x$tval,x$N,x$model ,x$location,x$scale,x$dff , x$hypothesis ),
                 "2" = t1e_BF10(x$tval,x$N,x$model,x$location,x$scale,x$dff , x$hypothesis,x$e )))

  output$BFt1 <- shiny::renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    \\textit{t}(', x$N , ') = ',x$tval,', \\textit{BF}_{10} = ', round(BF10, 4),", \\textit{BF}_{01} = ",round(1/BF10, 4), '
')


    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em('$$', table_html, '$$')
      )
    )
  })

})
# Reactive expression to calculate t-value
t_value <- shiny::reactive({
  # Extract inputs
  x_bar <- input$t1_s_mean
  mu <- input$t1_mean
  s <- input$t1_sd
  n <- input$t1_sample_size

  # Avoid division by zero
  if (s <= 0 || n <= 0) return(NA)

  # Compute t-value
  t <- (x_bar - mu) / (s / sqrt(n))
  t
})

# Render LaTeX output
output$cal_t1 <- shiny::renderUI({
  t <- t_value()
  n <- input$t1_sample_size
  df <- n - 1  # Degrees of freedom

  if (is.na(t)) return(shiny::HTML("Invalid input"))

  # LaTeX formula with df
  shiny::withMathJax(
    shiny::HTML(
      paste0(
        "\\( t = \\frac{\\bar{x} - \\mu}{s / \\sqrt{n}} = ",
        round(t, 4),
        ", \\quad df = ", df,
        "\\)"
      )
    )
  )
})

}


# ---- Server_t2.r ----

server_t2<- function(input, output, session) {
input_t2 <- shiny::reactive({
  mode_bf <- switch(input$Modet2,
                    "1" = 1,
                    "2" = 0,
                    "3" = 0)# mode
  direct <- switch(input$t2_direct,
                   "1" = "h1",
                   "0" = "h0")
  interval <- input$h0t2 # point null or interval

  e <- switch(input$h1t2e,        # bound for interval test
              "1" = c(input$lbt2e, input$ubt2e),
              "2" = input$ubt2e,
              "3" = input$lbt2e)
  inter <- switch(interval,
                  "1" = input$h1t2,
                  "2" = input$h1t2e)

  hypothesis <- switch(interval,
                       "1" =   switch(input$h1t2,        # direction of the test
                                      "1" = "!=",
                                      "2" =  ">",
                                      "3" =  "<"),
                       "2" = switch(input$h1t2e,        # direction of the test
                                    "1" = "!=",
                                    "2" =  ">",
                                    "3" =  "<"))



  model <- switch(input$modelt2,
                  "1" = "t-distribution",
                  "2" = "Normal",
                  "3" = "Moment")

  location <- input$lt2
  scale <- input$st2
  dff <- input$dft2
  de_an_prior <- switch(input$priort2,
                        "1" = 1,
                        "2" = 0)
  model_d <- switch(input$modelt2d,
                    "1" = "t-distribution",
                    "2" = "Normal",
                    "3" = "Moment",
                    "4" = "Point")
  location_d <- input$lt2d
  scale_d <- input$st2d
  dff_d <- input$dft2d
  D <- input$bt2
  type <- input$typet2
  target <- input$powert2
  alpha <- input$alphat2
  tval <- input$t2tval
  r <- switch(input$Modet2,
              "1" = input$rt2,
              "2" = input$n2t2/input$n1t2,
              "3" = input$rt2)
  N1 = input$n1t2
  N2 = input$n2t2
  pc   <- "1" %in% input$o_plot_t2
  rela <- "2" %in% input$o_plot_t2
  if (mode_bf!=1){
    pc=rela=F
  }
  # Add all variables to the final list
  list(
    mode_bf = mode_bf,
    direct = direct,
    interval = interval,
    hypothesis = hypothesis,
    e = e,
    model = model,
    location = location,
    scale = scale,
    dff = dff,
    de_an_prior = de_an_prior,
    model_d = model_d,
    location_d = location_d,
    scale_d = scale_d,
    dff_d = dff_d,
    type = type,
    D = D,
    target = target,
    alpha = alpha,
    tval = tval,
    r = r,
    N1=N1,
    N2=N2,
    df = df,
    pc = pc,
    rela = rela
  )
})

shiny::observeEvent(input$runt2, {
  t2 = input_t2()


  dat <- tryCatch({
    suppressWarnings(switch(t2$interval,
                            "1" = t2_Table(t2$D, t2$r, t2$target, t2$model, t2$location, t2$scale, t2$dff, t2$hypothesis,
                                           t2$model_d, t2$location_d, t2$scale_d, t2$dff_d, t2$de_an_prior, t2$N1, t2$N2, t2$mode_bf, t2$alpha,t2$direct),
                            "2" = t2e_table(t2$D, t2$r, t2$target, t2$model,t2$location, t2$scale, t2$dff, t2$hypothesis, t2$e,
                                            t2$model_d,t2$location_d, t2$scale_d, t2$dff_d, t2$de_an_prior, t2$mode_bf,  t2$N1, t2$N2, t2$alpha,t2$direct)
    ))
  }, error = function(e) {
    "Error"
  })
  output$result_t2 <- shiny::renderText({ paste("# Function to be used in R", show_t2_code(t2), sep = "\n") })

  output$priort2 <- shiny::renderPlot({
    suppressWarnings(switch(t2$interval,
           "1"=
             t1_prior_plot(
               D = t2$D,                  # Access 'D' explicitly
               target = t2$target,        # Access 'target' explicitly
               model = t2$model,          # Access 'model' explicitly
               location = t2$location,    # Access 'location' explicitly
               scale = t2$scale,          # Access 'scale' explicitly
               dff = t2$dff,              # Access 'dff' explicitly
               hypothesis = t2$hypothesis,  # Access 'hypothesis' explicitly
               model_d = t2$model_d,        # Access 'model_d' explicitly
               location_d = t2$location_d,  # Access 'location_d' explicitly
               scale_d = t2$scale_d,        # Access 'scale_d' explicitly
               dff_d = t2$dff_d,            # Access 'dff_d' explicitly
               de_an_prior = t2$de_an_prior   # Access 'de_an_prior' explicitly
             ), "2" =
             t1e_prior_plot(t2$model,
                            t2$location,
                            t2$scale,
                            t2$dff ,
                            t2$hypothesis,
                            t2$e,
                            t2$de_an_prior,
                            t2$model_d,
                            t2$scale_d,
                            t2$dff_d,
                            t2$location )

    ))

  })



  output$resultt2 <- shiny::renderUI({
    # Create the LaTeX formatted strings for the table
    if (identical(dat, "Error")){
      table_html <- shiny::em(shiny::span("\\(\\text{Error: the required } N > 10,000\\)", style = "color: red;"))

    }else{
    table_html <- paste0("$$",'
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    p\\text{(BF}_{10} > ', t2$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1,1], 3),nsmall=3), ' \\\\
    p\\text{(BF}_{01} > ', t2$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[1,3], 3),nsmall=3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', t2$D, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1,2], 3),nsmall=3), ' \\\\
    p\\text{(BF}_{10} > ', t2$D, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[1,4], 3),nsmall=3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N}_1 & ', dat[1,5], ' \\\\
    \\text{N}_2 & ', dat[1,6], ' \\\\
    \\end{array}
  ',"$$")}

    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em(table_html)
      )
    )
  })
  # Define reactive containers OUTSIDE the if blocks
  pc_t2   <- shiny::reactiveVal(NULL)
  rela_t2 <- shiny::reactiveVal(NULL)


  # ===================================================
  #                POWER CURVE (t2)
  # ===================================================

  if (isTRUE(t2$pc)) {

    output$plot_power_t2_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Power Curve}$$"))
      )
    })

    output$plot_power_t2 <- shiny::renderPlot({

      suppressWarnings(
        switch(
          t2$interval,

          "1" = Power_t2(
            t2$D, t2$model, t2$location, t2$scale, t2$dff, t2$hypothesis,
            t2$model_d, t2$location_d, t2$scale_d, t2$dff_d,
            t2$de_an_prior,
            dat[1,5],
            dat[1,6] / dat[1,5]
          ),

          "2" = Power_t2e(
            t2$D, t2$model, t2$location, t2$scale, t2$dff, t2$hypothesis,
            t2$model_d, t2$location_d, t2$scale_d, t2$dff_d,
            t2$de_an_prior,
            dat[1,5],
            dat[1,6] / dat[1,5],
            t2$e
          )
        )
      )

      pc_t2(grDevices::recordPlot())
    })

  } else {

    pc_t2(NULL)
    output$plot_power_t2_text <- shiny::renderUI(NULL)
    output$plot_power_t2      <- shiny::renderPlot(NULL)
  }


  # ===================================================
  #                RELATIONSHIP (t2)
  # ===================================================

  if (isTRUE(t2$rela)) {

    output$plot_rel_t2_text <- shiny::renderUI({
      shiny::tagList(
        shiny::withMathJax(shiny::em("$$\\text{Relationship between BF and data}$$"))
      )
    })

    output$plot_rel_t2 <- shiny::renderPlot({

      suppressWarnings(
        switch(
          t2$interval,

          "1" = t2_BF(
            t2$D, dat[1,5], t2$r, t2$target,
            t2$model, t2$location, t2$scale, t2$dff,
            t2$hypothesis
          ),

          "2" = t2e_BF(
            t2$D, dat[1,5], t2$r,
            t2$model, t2$location, t2$scale, t2$dff,
            t2$hypothesis, t2$e
          )
        )
      )

      rela_t2(grDevices::recordPlot())
    })

  } else {

    rela_t2(NULL)
    output$plot_rel_t2_text <- shiny::renderUI(NULL)
    output$plot_rel_t2      <- shiny::renderPlot(NULL)
  }


  output$export_t2 <- shiny::downloadHandler(
    filename = function() {
      "BayesPower-report.html"
    },
    content = function(file) {
      template_path <- system.file("report_templates", "report_t2.Rmd", package = "BayesPower")

      tempReport <- file.path(tempdir(), "report_t2.Rmd")
      file.copy(template_path, tempReport, overwrite = TRUE)

      rmarkdown::render(
        input = tempReport,output_format ="html_document",
        output_file = file,
        params = list(t2 = t2, dat = dat,pc_t2=pc_t2(),rela_t2=rela_t2()),  # ✅ pass to `params`
        envir = new.env(parent = globalenv())  # environment still required
      )
    }
  )







})

shiny::observeEvent(input$cal2, {
  t2 = input_t2()




  output$priort2 <- shiny::renderPlot({
    suppressWarnings(switch(t2$interval,
                            "1"=
                              t1_prior_plot(
                                D = t2$D,                  # Access 'D' explicitly
                                target = t2$target,        # Access 'target' explicitly
                                model = t2$model,          # Access 'model' explicitly
                                location = t2$location,    # Access 'location' explicitly
                                scale = t2$scale,          # Access 'scale' explicitly
                                dff = t2$dff,              # Access 'dff' explicitly
                                hypothesis = t2$hypothesis,  # Access 'hypothesis' explicitly
                                model_d = t2$model_d,        # Access 'model_d' explicitly
                                location_d = t2$location_d,  # Access 'location_d' explicitly
                                scale_d = t2$scale_d,        # Access 'scale_d' explicitly
                                dff_d = t2$dff_d,            # Access 'dff_d' explicitly
                                de_an_prior = 1   # Access 'de_an_prior' explicitly
                              ), "2" =
                              t1e_prior_plot(t2$model,
                                             t2$location,
                                             t2$scale,
                                             t2$dff ,
                                             t2$hypothesis,
                                             t2$e,
                                             1,
                                             t2$model_d,
                                             t2$scale_d,
                                             t2$dff_d,
                                             t2$location_d )

    ))

  })
  r = t2$N2/t2$N1
  N1 = t2$N1
  ddff = t2$N1+t2$N2-2

  BF10 <- suppressWarnings(switch(t2$interval,
                 "1" = t2_BF10(t2$tval,N1,r,t2$model ,t2$location,t2$scale,t2$dff , t2$hypothesis ),
                 "2" = t2e_BF10(t2$tval,N1,r,t2$model,t2$location,t2$scale,t2$dff , t2$hypothesis,t2$e )))


  output$result_t2 <- shiny::renderText({


    fmt_val <- function(val) {
      if (is.character(val)) {
        sprintf('"%s"', val)
      } else if (is.numeric(val) && length(val) > 1) {
        paste0("c(", paste(val, collapse = ","), ")")
      } else {
        as.character(val)
      }
    }



    args <- list(
      tval = t2$tval,
      N1 = t2$N1,
      N2 = t2$N2,
      model = t2$model,
      location = t2$location,
      scale = t2$scale,
      dff = t2$dff,
      hypothesis = t2$hypothesis
    )

    if (t2$interval != 1) {
      args$e <- t2$e
    }

    # Build string with each argument on a new line
    arg_strings <- sapply(names(args), function(arg) {

      val <- args[[arg]]
      arg_print <- arg

      ## model > prior_analysis
      if (arg == "model") {
        arg_print <- "prior_analysis"
      }

      ## e > ROPE
      if (arg == "e") {
        arg_print <- "ROPE"
      }

      ## hypothesis > alternative
      if (arg == "hypothesis") {

        arg_print <- "alternative"

        val <- switch(val,
                      "<"  = "less",
                      "!=" = "two.sided",
                      ">"  = "greater",
                      stop("Invalid hypothesis")
        )

        val <- shQuote(val)

      } else {
        val <- fmt_val(val)
      }

      sprintf("  %s = %s", arg_print, val)
    })

    call_string <- paste0(
      "# Function to be used in R\n",
      "BF10.ttest.TwoSample(\n",
      paste(arg_strings, collapse = ",\n"),
      "\n)"
    )

    call_string
  })



  output$BFt2 <- shiny::renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    \\textit{t}(', ddff, ') = ',t2$tval,', \\textit{BF}_{10} = ', round(BF10, 4),", \\textit{BF}_{01} = " ,round(1/BF10, 4),'
')


    # Render the table using MathJax
    shiny::tagList(
      # Render the table using MathJax
      shiny::withMathJax(
        shiny::em('$$', table_html, '$$')
      )
    )
  })

})

# Reactive calculation for independent t-test (equal variance)
t2_value <- shiny::reactive({
  x1 <- input$t2_mean1
  x2 <- input$t2_mean2
  s1 <- input$t2_sd1
  s2 <- input$t2_sd2
  n1 <- input$t2_n1
  n2 <- input$t2_n2

  # Pooled standard deviation
  s_p <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2))

  # t-value
  t <- (x1 - x2) / (s_p * sqrt(1/n1 + 1/n2))
  t
})

# Degrees of freedom
df <- shiny::reactive({
  input$t2_n1 + input$t2_n2 - 2
})

# Render LaTeX output
output$cal_t2 <- shiny::renderUI({
  t <- t2_value()
  d <- df()

  shiny::withMathJax(
    shiny::HTML(
      paste0(
        "\\( t = \\frac{\\bar{x}_1 - \\bar{x}_2}{s_p \\sqrt{1/n_1 + 1/n_2}} = ",
        round(t, 4),
        ", \\quad df = ", d,
        "\\)"
      )
    )
  )
})
}













# ---- twosample.r ----

# the Bayes Factor

t2_BF10 <-function(t,n1,r,model ,location,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )

  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

  error = 1e-10
  x <- sapply(t, function(ti) {
    int <- function(delta) {
      stats::dt(ti, df, ncp = delta * constant) * t1_prior(delta, location, scale, dff, model) / normalization
    }

    stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value /
      stats::dt(ti, df, ncp = 0)
  })

  return(x)
}



# finding the t that correspond to BF10=D
t2_BF10_bound <-function(D, n1,r,model ,location ,scale,dff , hypothesis){
  y <- numeric(0)
  Bound_finding <-function(t){
    t2_BF10(t,n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )- D
  }
  x <- tryCatch(stats::uniroot(Bound_finding, lower = -6, upper = 0)$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower =  0, upper = 6)$root, error = function(e) NA)
  results <- c(x, y)

  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")

  BF.vals  <- t2_BF10(results,n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")

  return(results[BF.close])
}


# finding the t that correspond to BF01=D
t2_BF01_bound <-function(D , n1,r,model ,location ,scale,dff , hypothesis){
  t2_BF10_bound(1/D, n1,r,model ,location ,scale,dff , hypothesis)
}


# p(BF01>D|H0)
t2_TNE <- function(t , n1,r,hypothesis){
  n2 = n1*r
  df = n1+n2-2

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  pro <- switch(hypothesis,
                "!=" = stats::pt(max(t), df) - stats::pt(min(t), df),
                ">"  = stats::pt(t, df),
                "<"  = 1 - stats::pt(t, df)
  )

  return(pro)

}

# p(BF10>D|H1)
t2_TPE <-function(t,n1,r,model ,location ,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (model == "Point"){
    pro = switch(hypothesis,
                 "!="= pnct(min(t),df,ncp = location*constant,lower  = T)+pnct(max(t),df,ncp = location*constant,lower  = F),
                 ">" = pnct(t,df,ncp = location*constant,lower  = F),
                 "<" = pnct(t,df,ncp = location*constant,lower  = T))
    return(pro)
  }

  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )



  x = NULL

  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(hypothesis,
                  "!=" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = FALSE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 + pro2
                  },
                  ">" = pnct(t, df, ncp = ncp, lower = FALSE),
                  "<" = pnct(t, df, ncp = ncp, lower = TRUE)
    )

    pro * t1_prior(delta, location, scale, dff, model) / normalization
  }

  error = 1e-4
  if (model == "Moment" & scale <.3 ){
    error = 1e-14
  }
  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = error,stop.on.error=FALSE)$value

  return(x)

}


# p(BF01>D|H1)
t2_FNE<-function(t,n1,r,model ,location ,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (model == "Point"){
    pro = switch(hypothesis,
                 "!="=  pnct(max(t),df,ncp = location*constant,lower  = T) - pnct(min(t),df,ncp = location*constant,lower  = T),
                 ">" = pnct(t,df,ncp = location*constant,lower  = T),
                 "<" = pnct(t,df,ncp = location*constant,lower  = F))
    return(pro)
  }
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x = NULL


  normalization <- if (hypothesis == "!=") 1 else
    switch(model,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))
  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(hypothesis,
                  "!=" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = TRUE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 - pro2
                  },
                  ">" = pnct(t, df, ncp = ncp, lower = TRUE),
                  "<" = pnct(t, df, ncp = ncp, lower = FALSE)
    )

    pro * t1_prior(delta, location, scale, dff, model) / normalization
  }


  x = stats::integrate(int,lower = bound[1],upper = bound[2],stop.on.error = FALSE)$value

  return(x)
}


# p(BF10>D|H0)
t2_FPE <- function(t,n1,r, hypothesis){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  pro <- switch(hypothesis,
                "!=" = stats::pt(max(t), df = df, lower.tail = FALSE) +
                  stats::pt(min(t), df = df, lower.tail = TRUE),
                ">"  = stats::pt(t, df = df, lower.tail = FALSE),
                "<"  = stats::pt(t, df = df, lower.tail = TRUE)
  )
  return(pro)

}

# Finding the degree of freedom that ensure p(BF10>D|H1) > targeted probability

t2_N_finder<-function(D,r,target,model,location,scale,dff, hypothesis ,
                   model_d,location_d,scale_d,dff_d,de_an_prior ,alpha){

  lower <- 2
  upper <- 10000
  t2 <- t2_BF10_bound(D, lower,r,model ,location ,scale,dff , hypothesis)
  p2 <- if (de_an_prior == 1)
    t2_TPE(t2 , n1=lower,r , model , location ,scale,dff , hypothesis) else
      t2_TPE(t2 , n1=lower,r , model_d , location_d,scale_d,dff_d, hypothesis)
  if (p2 > target) return(lower)
  Power_root <- function(n1) {
    t <- t2_BF10_bound(D, n1,r,model ,location ,scale,dff , hypothesis)
    if (de_an_prior == 1)
      t2_TPE(t , n1,r , model , location ,scale,dff , hypothesis) - target else
        t2_TPE(t , n1,r , model_d , location_d,scale_d,dff_d, hypothesis) - target
  }
  N1.power <-  stats::uniroot(Power_root,lower = lower,upper =  upper)$root
  #N1.power <- robust_uniroot(Power_root, lower = 2)
  t  <-  t2_BF10_bound(D,  N1.power,r,model ,location ,scale,dff , hypothesis)
  FPE <- t2_FPE(t,N1.power,r, hypothesis)
  if (FPE <= alpha) return(N1.power)
  alpha.root <- function(n1) {
    t <- t2_BF10_bound(D,  n1,r,model ,location ,scale,dff , hypothesis)
    t2_FPE(t,n1,r, hypothesis) - alpha
  }

  N1.alpha <- stats::uniroot(alpha.root, lower = N1.power, upper = upper)$root
  return(N1.alpha  )
}

t2_N_01_finder<-function(D,r,target,model,location,scale,dff, hypothesis ,
                         model_d,location_d,scale_d,dff_d,de_an_prior ,alpha){

  lower <- 2
  upper <- 10000
  t2 <- t2_BF01_bound(D, lower,r,model ,location ,scale,dff , hypothesis)
  TNE_lo <- t2_TNE(t2,lower,r,hypothesis)
  if (TNE_lo > target) return(lower)
  FNE_lo <-  if (de_an_prior == 1)
    t2_FNE(t2,lower,r,model ,location ,scale,dff , hypothesis ) else
      t2_FNE(t2,lower,r,model_d ,location_d ,scale_d,dff_d , hypothesis )
  if (TNE_lo > target&FNE_lo<alpha) return(lower)

  TN_root <- function(n1) {
    t <- t2_BF01_bound(D, n1,r,model ,location ,scale,dff , hypothesis)
    t2_TNE(t,lower,r,hypothesis)-target
  }
  N1.TN <-  stats::uniroot(TN_root,lower = lower,upper =  upper)$root
  t  <-  t2_BF01_bound(D,  N1.TN,r,model ,location ,scale,dff , hypothesis)
  FNE <- if (de_an_prior == 1)
    t2_FNE(t,N1.TN,r,model ,location ,scale,dff , hypothesis ) else
      t2_FNE(t,N1.TN,r,model_d ,location_d ,scale_d,dff_d , hypothesis )
  if (FNE <= alpha) return(N1.TN)

  FN.root <- function(n1) {
    t <- t2_BF01_bound(D,  n1,r,model ,location ,scale,dff , hypothesis)
    FNE <- if (de_an_prior == 1)
      t2_FNE(t,n1,r,model ,location ,scale,dff , hypothesis ) else
        t2_FNE(t,n1,r,model_d ,location_d ,scale_d,dff_d , hypothesis )
    FNE -alpha
  }
  N1.FN <- stats::uniroot(FN.root, lower = N1.TN, upper = upper)$root
  return( N1.FN )
}

# probability table
t2_Table <- function(D,r,target,model,location,scale,dff, hypothesis,
                  model_d,location_d,scale_d,dff_d, de_an_prior,N1,N2, mode_bf ,alpha ,direct){

  bound01 = as.numeric(0)
  bound10 = as.numeric(0)

  if (mode_bf == 1) {
    n1 <- switch(direct,
                 "h1" = {ceiling(t2_N_finder(D, r, target, model, location, scale, dff,
                              hypothesis, model_d, location_d, scale_d, dff_d,
                              de_an_prior, alpha))},
                 "h0" = {ceiling(t2_N_01_finder(D, r, target, model, location, scale, dff,
                                             hypothesis, model_d, location_d, scale_d, dff_d,
                                             de_an_prior, alpha))}     )
    n2 <- n1 * r
  } else {
    n1 <- N1
    n2 <- N2
    r  <- n2 / n1
  }

  # t bounds:
  t10 <- t2_BF10_bound(D, n1,r,model,location,scale,dff , hypothesis)
  t01 <- t2_BF01_bound(D, n1,r,model,location,scale,dff , hypothesis)

  # max BF10 possible:
  max_BF <- 1 / t2_BF10(0,n1,r,model ,location,scale,dff , hypothesis )
  BF_D   <- t10

  # FPE and TPE:
  FPE       <- t2_FPE(t10,n1,r, hypothesis)
  if (de_an_prior == 1) {
    TPE       <- t2_TPE(t10,n1,r,model ,location ,scale,dff , hypothesis )
    TPR_model <- model
    TPR_loc   <- location
    TPR_scale <- scale
    TPR_dff   <- dff
  } else {
    TPE       <- t2_TPE(t10,n1,r,model_d ,location_d ,scale_d,dff_d , hypothesis )
    TPR_model <- model_d
    TPR_loc   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }

  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t2_FNE(t01, n1,r,TPR_model, TPR_loc, TPR_scale, TPR_dff, hypothesis )
    TNE <- t2_TNE(t01 , n1,r,hypothesis)
  }

  # table:
  tab.names <- c(
    sprintf("p(BF10 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H0)", D),
    sprintf("p(BF10 > %0.f | H0)", D),
    "Required N1",
    "Required N2"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n1,n2, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}


# plots for showing the relationship between BF and t-values

t2_BF <- function(D, n1, r, target,
                  model, location, scale, dff, hypothesis) {

  tt <- seq(-5, 5, 0.2)

  ## ---------- BF10 ----------
  BF10   <- t2_BF10(tt, n1, r, model, location, scale, dff, hypothesis)
  t.BF10 <- t2_BF10_bound(D, n1, r, model, location, scale, dff, hypothesis)

  df10 <- data.frame(t = tt, BF = BF10)

  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(D)~" when t = "~.(round(t.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(D)~" when t = "~.(round(t.BF10[1], 2))~
                  " or "~.(round(t.BF10[2], 2))))
  }

  x_breaks_10 <- sort(unique(c(-5, 5, round(t.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(t, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = t.BF10, linetype = "solid") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = x_breaks_10) +
    ggplot2::labs(
      x = "t-value",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 13, face = "bold"),
      axis.text  = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  ## ---------- BF01 ----------
  BF01   <- 1 / BF10
  t.BF01 <- t2_BF01_bound(D, n1, r, model, location, scale, dff, hypothesis)

  df01 <- data.frame(t = tt, BF = BF01)

  max.BF01   <- 1 / t2_BF10(0, n1, r, model, location, scale, dff, "!=")
  impossible <- (hypothesis == "!=") &&
    (max.BF01 < D || identical(t.BF01, "bound cannot be found"))

  if (impossible) {

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(size = 1.2, color = "black") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = c(-5, 5)) +
      ggplot2::labs(
        x = "t-value",
        y = bquote("BF"['01'] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[01]~"="~.(D)))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.text  = ggplot2::element_text(size = 11),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

  } else {

    main.bf01 <- if (length(t.BF01) == 1) {
      bquote(bold("BF"[0][1]~"="~.(D)~" when t = "~.(round(t.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(D)~" when t = "~.(round(t.BF01[1], 2))~
                    " or "~.(round(t.BF01[2], 2))))
    }

    x_breaks_01 <- sort(unique(c(-5, 5, round(t.BF01, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(size = 1.2, color = "black") +
      ggplot2::geom_vline(xintercept = t.BF01, linetype = "solid") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = x_breaks_01) +
      ggplot2::labs(
        x = "t-value",
        y = bquote("BF"['01'] * " (log scale)"),
        title = main.bf01
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.text  = ggplot2::element_text(size = 11),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}



Power_t2 <- function(D, model, location, scale, dff, hypothesis,
                            model_d, location_d, scale_d, dff_d,
                            de_an_prior, n1, r) {

  Total_ <- n1 + n1 * r
  smin <- 4
  smax <- Total_ * 1.2
  sdf  <- seq(smin, smax, length.out = 31)
  sn1  <- sdf / (1 + r)

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(sdf))

  for (i in seq_along(sdf)) {

    t10 <- t2_BF10_bound(D, sn1[i], r, model, location, scale, dff, hypothesis)
    t01 <- t2_BF01_bound(D, sn1[i], r, model, location, scale, dff, hypothesis)

    TPE[i] <- if (de_an_prior == 1) {
      t2_TPE(t10, sn1[i], r, model, location, scale, dff, hypothesis)
    } else {
      t2_TPE(t10, sn1[i], r, model_d, location_d, scale_d, dff_d, hypothesis)
    }

    FPE[i] <- t2_FPE(t10, sn1[i], r, hypothesis)

    FNE[i] <- if (de_an_prior == 1) {
      t2_FNE(t01, sn1[i], r, model, location, scale, dff, hypothesis)
    } else {
      t2_FNE(t01, sn1[i], r, model_d, location_d, scale_d, dff_d, hypothesis)
    }

    TNE[i] <- t2_TNE(t01, sn1[i], r, hypothesis)
  }

  ## ---------- Data for ggplot ----------

  df1 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = sdf,
      `True Positive`  = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df1$Type <- factor(df1$Type, levels = c("True Positive", "False Positive"))

  df2 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = sdf,
      `True Negative`  = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df2$Type <- factor(df2$Type, levels = c("True Negative", "False Negative"))

  ## ---------- Styling ----------

  type_colors <- c(
    "True Positive"  = "black",
    "False Positive" = "grey50",
    "True Negative"  = "black",
    "False Negative" = "grey50"
  )

  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x  = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 12),
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  ## ---------- Plots ----------

  p1 <- ggplot2::ggplot(df1,
                        ggplot2::aes(x = SampleSize,
                                     y = Probability,
                                     color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  p2 <- ggplot2::ggplot(df2,
                        ggplot2::aes(x = SampleSize,
                                     y = Probability,
                                     color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  ## ---------- Combine ----------

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}


# ---- twosample_e.r ----

t2e_BF10i <-function(t,n1,r,model ,location,scale,dff , hypothesis,e ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = e),
                      "!=" = c(a = e[1], b = e[2])
  )
  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )
  normalizationh1 <- norm_h1(hypothesis, model, bound_h1, location, scale, dff)


  normalizationh0 <- norm_h0(model, bound_h0, location, scale, dff)


  int  <- function(delta){
    stats::dt(t,df,ncp = delta *constant)* te_prior(delta,location,scale,dff,model)/normalizationh1}

  error = 1e-4

  if (hypothesis == "!="){
  lh1 = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value
  }else{
    lh1 = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value}


  int  <- function(delta){
    stats::dt(t,df,ncp = delta *constant)* te_prior(delta,location,scale,dff,model)/normalizationh0}

  lh0 = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(lh1/lh0)
}

t2e_BF10 <-function(t,n1,r,model,location,scale,dff , hypothesis,e ){
  sapply(t, function(ti) t2e_BF10i(ti,n1,r,model ,location,scale,dff , hypothesis,e ))
}

t2e_BF10_bound <-function(D, n1,r,model,location,scale,dff , hypothesis,e){

  y <- numeric(0)
  Bound_finding <-function(t){
    t2e_BF10(t,n1,r,model,location,scale,dff , hypothesis,e )- D
  }

  switch(hypothesis,
         "!=" ={
           x <- tryCatch(stats::uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
           y <- tryCatch(stats::uniroot(Bound_finding, lower =  0, upper = 20)$root, error = function(e) NA)
         },
         ">"={
           x <- tryCatch(stats::uniroot(Bound_finding, lower = 0, upper = 20)$root, error = function(e) NA)
         },
         "<" = {
           x <- tryCatch(stats::uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
         })

  results <- c(x, y)
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")
  BF.vals  <- t2e_BF10(results,n1,r,model,location,scale,dff , hypothesis,e )
  BF.close <- which(round(BF.vals, 2) == round(D, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

# finding the t that correspond to BF10=D
t2e_BF01_bound <-function(D, n1,r,model,location,scale,dff , hypothesis,e){
  t2e_BF10_bound(1/D, n1,r,model,location,scale,dff , hypothesis,e)

}

t2e_TPE <-function(t,n1,r,model ,location,scale,dff , hypothesis ,e){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {pnct(min(t),df,ncp= location*constant,lower = T)+ pnct(max(t),df,ncp=location*constant,lower = F)},
               "<"  = {pnct(t,df,ncp = location *constant,lower  = T)},
               ">"  = {pnct(t,df,ncp = location *constant,lower  = F)}
               )
    return(x)
  }

  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = e),
                      "!=" = c(a = e[1], b = e[2])
  )

  normalizationh1 <- norm_h1(hypothesis, model, bound_h1, location, scale, dff)



  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(hypothesis,
                  "!=" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = FALSE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 + pro2
                  },
                  ">" = pnct(t, df, ncp = ncp, lower = FALSE),
                  "<" = pnct(t, df, ncp = ncp, lower = TRUE)
    )

    pro * te_prior(delta,location, scale, dff, model) / normalizationh1
  }


  error = 1e-4

  if (hypothesis == "!="){
    x = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value
  }else{
    x = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value

  }
  return(x)

}

t2e_FNE <-function(t,n1,r,model ,location,scale,dff , hypothesis ,e){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  if (model =="Point"){
    x = switch(hypothesis,
               "!=" = {pnct(max(t),df,ncp= location*constant,lower = T)- pnct(min(t),df,ncp=location*constant,lower = T)},
               "<"  = {pnct(t,df,ncp = location *constant,lower  = F)},
               ">"  = {pnct(t,df,ncp = location *constant,lower  = T)}
    )
    return(x)
  }
  bound_h1  <- switch(hypothesis,
                      ">" = c(a = e, b = Inf),
                      "<" = c(a = -Inf, b = e),
                      "!=" = c(a = e[1], b = e[2])
  )

  normalizationh1 <- norm_h1(hypothesis, model, bound_h1, location, scale, dff)


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(hypothesis,
                  "!=" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = TRUE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 - pro2
                  },
                  ">" = pnct(t, df, ncp = ncp, lower = TRUE),
                  "<" = pnct(t, df, ncp = ncp, lower = FALSE)
    )

    pro * te_prior(delta, location,scale, dff, model) / normalizationh1
  }

  error = 1e-4
  if (hypothesis == "!="){
    x = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value
  }else{
    x = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value

  }
  return(x)

}


t2e_TNE <-function(t,n1,r,model ,location,scale,dff , hypothesis ,e){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))


  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )

  normalizationh0 <- norm_h0(model, bound_h0, location, scale, dff)


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(hypothesis,
                  "!=" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = TRUE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 - pro2
                  },
                  ">" = pnct(t, df, ncp = ncp, lower = TRUE),
                  "<" = pnct(t, df, ncp = ncp, lower = FALSE)
    )

    pro * te_prior(delta,location, scale, dff, model) / normalizationh0
  }
  error = 1e-4
    x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(x)

}

t2e_FPE <-function(t,n1,r,model ,location,scale,dff , hypothesis ,e){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  bound_h0  <- switch(hypothesis,
                      ">" = c(a = 0, b = e),
                      "<" = c(a = e, b = 0),
                      "!=" = c(a = e[1], b = e[2])
  )

  normalizationh0 <- norm_h0(model, bound_h0, location, scale, dff)


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(hypothesis,
                  "!=" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = FALSE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 + pro2
                  },
                  ">" = pnct(t, df, ncp = ncp, lower = FALSE),
                  "<" = pnct(t, df, ncp = ncp, lower = TRUE)
    )

    pro * te_prior(delta,location, scale, dff, model) / normalizationh0
  }
  error = 1e-4
  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value

  return(x)

}



t2e_N_finder<-function(D,r,target,model,location,scale,dff, hypothesis,e ,
                   model_d,location_d,scale_d,dff_d, de_an_prior,alpha ){

  lower <- 2
  t2 <-t2e_BF10_bound(D, lower,r,model,location,scale,dff , hypothesis,e)
  p2 <- if (de_an_prior == 1)
    t2e_TPE (t2,lower,r,model ,location,scale,dff , hypothesis,e) else
      t2e_TPE (t2,lower,r,model_d ,location_d,scale_d,dff_d , hypothesis,e)
  if (p2 > target) return(lower)

  Power_root <- function(n1) {

    t <- t2e_BF10_bound(D, n1,r,model,location,scale,dff , hypothesis,e)

    pro <- if (de_an_prior == 1) {
      t2e_TPE (t,n1,r,model ,location,scale,dff , hypothesis,e )
    } else {
      t2e_TPE (t,n1,r,model_d ,location_d,scale_d,dff_d , hypothesis,e)
    }

    target - pro
  }
  N1.power <- robust_uniroot(Power_root, lower = 2)
  t <- t2e_BF10_bound(D, N1.power,r,model,location,scale,dff , hypothesis,e)
  FPE <-t2e_FPE(t,N1.power,r,model,location,scale,dff , hypothesis ,e)

  if (FPE <= alpha) return(N1.power + 1)

  alpha.root <- function(n1) {
    t <- t2e_BF10_bound(D, n1,r,model,location,scale,dff , hypothesis,e)
    pro <- t2e_FPE(t,n1,r,model,location ,scale,dff , hypothesis ,e)
    return(pro - alpha)
  }
  N1.alpha <- robust_uniroot(alpha.root , lower = N1.power)
  return(N1.alpha)
  }
t2e_N_01_finder<-function(D,r,target,model,location,scale,dff, hypothesis,e ,
                          model_d,location_d,scale_d,dff_d, de_an_prior,alpha ){

  lower <- 10
  t2 <-t2e_BF01_bound(D, lower,r,mode,location,scale,dff , hypothesis,e)
  TNE_lo <- t2e_TNE(t2,lower,r,model ,location,scale,dff , hypothesis,e)
  if (TNE_lo > target) return(lower)
  FNE_lo <- if (de_an_prior == 1)
    t2e_FNE (t2,lower,r,model,location ,scale,dff , hypothesis,e ) else
      t2e_FNE (t2,lower,r,model_d ,location_d,scale_d,dff_d , hypothesis,e )
  if (TNE_lo > target&FNE_lo<alpha) return(lower)

  TN_root <- function(n1) {

    t <- t2e_BF01_bound(D, n1,r,model,location,scale,dff , hypothesis,e)

    pro <- t2e_TNE(t,n1,r,model ,location,scale,dff , hypothesis,e)

    target - pro
  }
  N1.TN <- robust_uniroot(TN_root, lower = 2)
  t <- t2e_BF01_bound(D, N1.TN,r,model,location,scale,dff , hypothesis,e)
  FNE <-t2e_FNE(t,N1.TN,r,model ,location,scale,dff , hypothesis ,e)

  if (FNE <= alpha) return(N1.TN + 1)

  FN.root <- function(n1) {
    t <- t2e_BF01_bound(D, n1,r,model,location,scale,dff , hypothesis,e)
    pro <- if (de_an_prior == 1) {
      t2e_FNE (t,n1,r,model ,location,scale,dff , hypothesis,e )
    } else {
      t2e_FNE (t,n1,r,model_d ,location_d,scale_d,dff_d , hypothesis,e)
    }
    return(pro - alpha)
  }
  N1.FN <- robust_uniroot(FN.root , lower = N1.TN)
  return(N1.FN)
}

t2e_table<-function(D,r,target,model,location,scale,dff, hypothesis,e ,
                    model_d,location_d,scale_d,dff_d, de_an_prior,mode_bf,N1,N2,alpha ,direct){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)

  if (mode_bf == 1){

    n1 = switch(direct,
                "h1" = ceiling(t2e_N_finder(D,r,target,model,location,scale,dff, hypothesis,e ,
                                                      model_d,location_d,scale_d,dff_d, de_an_prior,alpha )),
                "h0" = ceiling(t2e_N_01_finder(D,r,target,model,location,scale,dff, hypothesis,e ,
                                                         model_d,location_d,scale_d,dff_d, de_an_prior,alpha ) ))
    n2 = n1*r
  } else {
    n1 = N1
    n2 = N2
    r= n2/n1
  }
  # t bounds:
  t10 <- t2e_BF10_bound(D, n1,r,model,location,scale,dff , hypothesis,e)
  t01 <- t2e_BF01_bound(D, n1,r,model,location,scale,dff , hypothesis,e)

  # max BF10 possible:
  max_BF <- 1 / t2e_BF10i(0,n1,r,model ,location,scale,dff , hypothesis,e )
  BF_D   <- t10

  # FPE and TPE:
  FPE       <- t2e_FPE(t10,n1,r,model,location ,scale,dff , hypothesis ,e)
  if (de_an_prior == 1) {
    TPE       <- t2e_TPE(t10,n1,r,model ,location,scale,dff , hypothesis,e )
    TPR_model <- model
    TPR_location   <- location
    TPR_scale <- scale
    TPR_dff   <- dff
  } else {
    TPE       <- t2e_TPE(t10,n1,r,model_d ,location_d,scale_d,dff_d , hypothesis,e )
    TPR_model <- model_d
    TPR_location   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }

  # FNE and TNE:
  if (any(hypothesis == "!=" & max_BF < D | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t2e_FNE(t01,n1,r,TPR_model ,TPR_location,TPR_scale,TPR_dff , hypothesis ,e)
    TNE <- t2e_TNE(t01,n1,r,model,location ,scale,dff , hypothesis ,e)
  }

  # table:
  tab.names <- c(
    sprintf("p(BF10 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H1)", D),
    sprintf("p(BF01 > %0.f | H0)", D),
    sprintf("p(BF10 > %0.f | H0)", D),
    "Required N1",
    "Required N2"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n1,n2, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}

t2e_BF <- function(D, n1, r,
                   model, location, scale, dff, hypothesis, e) {

  tt <- seq(-5, 5, 0.2)

  ## ---------- BF10 ----------
  BF10   <- t2e_BF10(tt, n1, r, model, location, scale, dff, hypothesis, e)
  t.BF10 <- t2e_BF10_bound(D, n1, r, model, location, scale, dff, hypothesis, e)

  df10 <- data.frame(t = tt, BF = BF10)

  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(D)~" when t = "~.(round(t.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(D)~" when t = "~.(round(t.BF10[1], 2))~
                  " or "~.(round(t.BF10[2], 2))))
  }

  x_breaks_10 <- sort(unique(c(-5, 5, round(t.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(t, BF)) +
    ggplot2::geom_line(size = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = t.BF10, linetype = "dashed") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = x_breaks_10) +
    ggplot2::labs(
      x = "t-value",
      y = expression("BF"[10] * " (log scale)"),
      title = main.bf10
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 13, face = "bold"),
      axis.text  = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  ## ---------- BF01 ----------
  BF01   <- 1 / BF10
  t.BF01 <- t2e_BF01_bound(D, n1, r, model, location, scale, dff, hypothesis, e)

  df01 <- data.frame(t = tt, BF = BF01)

  max.BF01   <- 1 / t2e_BF10i(0, n1, r, model, location, scale, dff, hypothesis, e)
  impossible <- (max.BF01 < D || identical(t.BF01, "bound cannot be found"))

  if (impossible) {

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(size = 1.2, color = "black") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = c(-5, 5)) +
      ggplot2::labs(
        x = "t-value",
        y = bquote("BF"['01'] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[01]~"="~.(D)))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.text  = ggplot2::element_text(size = 11),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

  } else {

    main.bf01 <- if (length(t.BF01) == 1) {
      bquote(bold("BF"[0][1]~"="~.(D)~" when t = "~.(round(t.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(D)~" when t = "~.(round(t.BF01[1], 2))~
                    " or "~.(round(t.BF01[2], 2))))
    }

    x_breaks_01 <- sort(unique(c(-5, 5, round(t.BF01, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(size = 1.2, color = "black") +
      ggplot2::geom_vline(xintercept = t.BF01, linetype = "dashed") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = x_breaks_01) +
      ggplot2::labs(
        x = "t-value",
        y = bquote("BF"['01'] * " (log scale)"),
        title = main.bf01
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.text  = ggplot2::element_text(size = 11),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}




Power_t2e <- function(D, model, location, scale, dff, hypothesis,
                      model_d, location_d, scale_d, dff_d,
                      de_an_prior, n1, r, e) {

  Total_ <- n1 + n1 * r
  smin   <- 4
  smax   <- Total_ * 1.2
  sdf    <- seq(smin, smax, length.out = 31)
  sn1    <- sdf / (1 + r)

  TPE <- FPE <- TNE <- FNE <- numeric(length(sdf))

  for (i in seq_along(sdf)) {

    t10 <- t2e_BF10_bound(D, sn1[i], r, model, location, scale, dff, hypothesis, e)
    t01 <- t2e_BF01_bound(D, sn1[i], r, model, location, scale, dff, hypothesis, e)

    TPE[i] <- if (de_an_prior == 1) {
      t2e_TPE(t10, sn1[i], r, model, location, scale, dff, hypothesis, e)
    } else {
      t2e_TPE(t10, sn1[i], r, model_d, location_d, scale_d, dff_d, hypothesis, e)
    }

    FPE[i] <- t2e_FPE(t10, sn1[i], r, model, location, scale, dff, hypothesis, e)
    TNE[i] <- t2e_TNE(t01, sn1[i], r, model, location, scale, dff, hypothesis, e)

    FNE[i] <- if (de_an_prior == 1) {
      t2e_FNE(t01, sn1[i], r, model, location, scale, dff, hypothesis, e)
    } else {
      t2e_FNE(t01, sn1[i], r, model_d, location_d, scale_d, dff_d, hypothesis, e)
    }
  }

  ## ---------- Data frames ----------
  df1 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = sdf,
      `True Positive`  = TPE,
      `False Positive` = FPE,
      check.names = FALSE
    ),
    cols = c(`True Positive`, `False Positive`),
    names_to = "Type",
    values_to = "Probability"
  )
  df1$Type <- factor(df1$Type, levels = c("True Positive", "False Positive"))

  df2 <- tidyr::pivot_longer(
    data = data.frame(
      SampleSize = sdf,
      `True Negative`  = TNE,
      `False Negative` = FNE,
      check.names = FALSE
    ),
    cols = c(`True Negative`, `False Negative`),
    names_to = "Type",
    values_to = "Probability"
  )
  df2$Type <- factor(df2$Type, levels = c("True Negative", "False Negative"))

  ## ---------- Colors ----------
  type_colors <- c(
    "True Positive"  = "black",
    "False Positive" = "grey50",
    "True Negative"  = "black",
    "False Negative" = "grey50"
  )

  ## ---------- Themes ----------
  clean_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.x  = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 12),
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  legend_theme <- ggplot2::theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 12)
  )

  ## ---------- BF10 plot ----------
  p1 <- ggplot2::ggplot(df1,
                        ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  ## ---------- BF01 plot ----------
  p2 <- ggplot2::ggplot(df2,
                        ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(D)))
    ) +
    clean_theme +
    legend_theme

  ## ---------- Combine ----------
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

