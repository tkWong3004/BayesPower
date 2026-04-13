#' @useDynLib BayesPower, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom rlang .data
#' @importFrom utils globalVariables
NULL

# Avoid R CMD check notes regarding the used variables in ggplot functions.
utils::globalVariables(c(
  "SampleSize", "Probability", "Type",
  "True Positive", "False Positive",
  "True Negative", "False Negative",
  "x", "BF","r","f","Density", "Prior", "lambda2", "theta", "fsq"
))

fmt_val <- function(x) {
  if (is.numeric(x) && length(x) == 1) return(as.character(x))
  if (is.numeric(x) && length(x) > 1) return(paste(x, collapse = ", "))
  if (is.character(x)) return(shQuote(x))
  return(as.character(x))
}
# ---- ANOVA.r ----

# k = number of predictor in the full prior_analysis
# p = number  of predictor in the reduced prior_analysis
# m = N-p
# q = k-p

F_prior<- function(fsq,q,dff,rscale,f,prior_analysis) {


  switch(prior_analysis,
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


F_BF <- function(f, q, m, dff, rscale, f_m, prior_analysis) {
  sapply(f, function(fi) {
    int <- function(fsq) {
      stats::df(fi, q, m - q, ncp = m * fsq) * F_prior(fsq, q, dff, rscale, f_m, prior_analysis)
    }
    lh1 <- stats::integrate(int, lower = 0, upper = Inf, stop.on.error = FALSE, rel.tol = 1e-4)$value
    lh0 <- stats::df(fi, q, m - q)
    lh1 / lh0
  })
}


F_BF_bound_10 <-function(threshold,q,m,dff,rscale,f_m,prior_analysis){
  x = numeric(0)
  Bound_finding <-function(f){
    F_BF(f,q,m,dff,rscale,f_m,prior_analysis)-threshold
  }

  x = tryCatch( stats::uniroot(Bound_finding,lower=0.01,upper = 40 )$root, error=function(e){})
  if (length(x) == 0) return("no bound is found")

  return(x)
}

F_BF_bound_01 <-function(threshold,q,m,dff,rscale,f_m,prior_analysis){
  F_BF_bound_10(1/threshold,q,m,dff,rscale,f_m,prior_analysis)
}


F_TPE<-function(f,q,m,dff,rscale,f_m,prior_analysis){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  if (prior_analysis == "Point"){
    x = stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = F)
    return(x)
  }
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,prior_analysis)
  }
  x = stats::integrate(int,lower = 0,upper = Inf)$value
  return(x)
}

F_FNE<-function(f,q,m,dff,rscale,f_m,prior_analysis){
  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  if (prior_analysis == "Point"){

    x = stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = T)
    return(x)
  }
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,prior_analysis)
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

f_N_finder<-function(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,dff_d,rscale_d,f_m_d,prior_design,de_an_prior,false_rate){
  q= k-p
  lower = 2*k-p+1
  m= lower-p
  f <- F_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis)
  p2 <- if (de_an_prior == 1)
    F_TPE(f,q,m,dff,rscale,f_m,prior_analysis) else
      F_TPE(f,q,m,dff_d,rscale_d,f_m_d,prior_design)
  if (p2 > true_rate) return(lower)

  Power_root <- function(n){
    m= n-p
    f = F_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis)

    pro <- if (de_an_prior == 1)
      F_TPE(f,q,m,dff,rscale,f_m,prior_analysis) else
        F_TPE(f,q,m,dff_d,rscale_d,f_m_d,prior_design)

    return(pro-true_rate)
  }

  N.power = stats::uniroot(Power_root,lower = lower,upper =  10000)$root
  m= N.power-p
  f = F_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis)
  FPE = F_FPE(f,q,m)

  if (FPE <= false_rate) return(N.power + 1)
  alpha_root <- function(n){
    m= n-p
    f = F_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis)
    pro = F_FPE(f,q,m)


    return(pro-false_rate)
  }
  N.alpha = stats::uniroot(alpha_root,lower = N.power,upper =  10000)$root
  return(N.alpha)
}

f_N_01_finder<-function(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,dff_d,rscale_d,f_m_d,prior_design,de_an_prior,false_rate){
  q= k-p
  lower = 2*k-p+1
  m= lower-p
  upper =  10000
  f <- F_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis)
  TNE_lo <- F_TNE(f,q,m)
  FNE_lo <- if (de_an_prior == 1)
    F_FNE(f,q,m,dff,rscale,f_m,prior_analysis) else
      F_FNE(f,q,m,dff_d,rscale_d,f_m_d,prior_design)

  if (TNE_lo > true_rate && FNE_lo < false_rate) {
    return(lower)
  } else if (TNE_lo > true_rate) {
    FN.root <- function(n) {
      q= k-p
      m= n-p
      f <- F_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis)
      FNE <- if (de_an_prior == 1)
        F_FNE(f,q,m,dff,rscale,f_m,prior_analysis) else
          F_FNE(f,q,m,dff_d,rscale_d,f_m_d,prior_design)
      FNE - false_rate
    }
    return(stats::uniroot(FN.root, lower = lower, upper = upper)$root)
  }

  TN_root <- function(n){
    m= n-p
    f = F_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis)
    TNE <- F_TNE(f,q,m)

    return(TNE-true_rate)
  }

  N.TN = stats::uniroot(TN_root,lower = lower,upper =  upper)$root
  m= N.TN-p
  f = F_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis)
  FNE = if (de_an_prior == 1)
    F_FNE(f,q,m,dff,rscale,f_m,prior_analysis) else
      F_FNE(f,q,m,dff_d,rscale_d,f_m_d,prior_design)

  if (FNE <= false_rate) return(N.TN + 1)

  FN.root <- function(n) {
    q= k-p
    m= n-p
    f <- F_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis)
    FNE <- if (de_an_prior == 1)
      F_FNE(f,q,m,dff,rscale,f_m,prior_analysis) else
        F_FNE(f,q,m,dff_d,rscale_d,f_m_d,prior_design)
    FNE - false_rate
  }

  N.FN = stats::uniroot(FN.root, lower = N.TN, upper = upper)$root
  return(N.FN)
}

f_table<-function(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,
                  dff_d,rscale_d,f_m_d,prior_design,de_an_prior,n, mode_bf,false_rate,type_rate ){


  if (mode_bf == 1){

    n = switch(type_rate,
               "positive"= ceiling(f_N_finder(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,dff_d,rscale_d,f_m_d,prior_design,de_an_prior,false_rate )),
               "negative" = ceiling(f_N_01_finder(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,dff_d,rscale_d,f_m_d,prior_design,de_an_prior,false_rate )))
  } else {
    n=n
  }
  q= k-p
  m= n-p

  # f bounds:
  f10 <- F_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis)
  f01 <- F_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis)

  # max BF10 possible:
  max_BF <- 1/F_BF(0.00001,q,m,dff,rscale,f_m,prior_analysis)
  BF_D   <- f10

  # FPE and TPE:
  FPE       <- F_FPE(f10,q,m)
  if (de_an_prior == 1) {
    TPE       <- F_TPE(f10,q,m,dff,rscale,f_m,prior_analysis)
    TPR_dff   <- dff
    TPR_rscale<- rscale
    TPR_f_m   <- f_m
    TPR_prior <- prior_analysis
  } else {
    TPE       <- F_TPE(f10,q,m,dff_d,rscale_d,f_m_d,prior_design)
    TPR_dff   <- dff_d
    TPR_rscale<- rscale_d
    TPR_f_m   <- f_m_d
    TPR_prior <- prior_design
  }

  # FNE and TNE:
  if ( max_BF < threshold | BF_D == "no bound is found") {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- F_FNE(f01,q,m,TPR_dff,TPR_rscale,TPR_f_m,TPR_prior)
    TNE <- F_TNE(f01,q,m)

  }

  # table:
  tab.names <- c(
    "TruePositve",
    "FalseNegative",
    "TrueNegative",
    "FalsePositive",
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table

}

prior_plot_f <- function(q, dff, rscale, f_m, prior_analysis,
                         dff_d, rscale_d, f_m_d, prior_design,
                         de_an_prior) {

  # ---- Sequence ----
  fsq <- seq(0.01, 3, 0.025)
  plot.bounds <- c(0.01, 3)

  # ---- Compute Analysis Prior ----
  prior_analysis_dens <- F_prior(fsq, q, dff, rscale, f_m, prior_analysis)

  # ---- Base data frame ----
  df <- data.frame(
    lambda2 = fsq,
    Density = prior_analysis_dens,
    Prior = "H1 - Analysis Prior"
  )

  # ---- Conditionally add Design Prior ----
  if (de_an_prior == 0) {

    if (prior_design == "Point") {
      # Dummy line for legend only
      df_design <- data.frame(
        lambda2 = c(NA, NA),
        Density = c(NA, NA),
        Prior = "H1 - Design Prior"
      )
      df <- rbind(df, df_design)

    } else {
      prior_design_dens <- F_prior(fsq, q, dff_d, rscale_d, f_m_d, prior_design)

      df_design <- data.frame(
        lambda2 = fsq,
        Density = prior_design_dens,
        Prior = "H1 - Design Prior"
      )

      df <- rbind(df, df_design)
    }
  }

  # ---- Legend position ----
  legend_pos <- c(0.65, 0.95)

  # ---- Base ggplot ----
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = lambda2,
                                    y = Density,
                                    color = Prior,
                                    linetype = Prior)) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(
      "H1 - Analysis Prior" = "black",
      "H1 - Design Prior"   = "gray"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "H1 - Analysis Prior" = "solid",
      "H1 - Design Prior"   = "dashed"
    )) +
    ggplot2::labs(
      x = expression(bold(lambda^2)),
      y = "density",
      title = bquote(bold("Prior distribution on "~lambda^2~
                            " under the alternative hypothesis"))
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.justification = c(0, 1),
      legend.background =
        ggplot2::element_rect(fill = scales::alpha("white", 0.8),
                              color = NA),
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 1.5)
      )
    )

  # ---- Add vertical arrow for Point design prior ----
  if (de_an_prior == 0 && prior_design == "Point") {

    ylim_max <- max(prior_analysis_dens, na.rm = TRUE)

    p <- p +
      ggplot2::annotate("segment",
                        x = f_m_d, xend = f_m_d,
                        y = 0, yend = ylim_max,
                        color = "gray",
                        linetype = "dashed",
                        arrow = ggplot2::arrow(
                          length = grid::unit(0.1, "inches")
                        ))
  }

  # ---- Set limits ----
  ylim_max <- max(df$Density, na.rm = TRUE)

  p <- p +
    ggplot2::coord_cartesian(
      xlim = plot.bounds,
      ylim = c(0, ylim_max)
    )

  return(p)
}
bf10_f <- function(threshold, n, k, p, dff, rscale, f_m, prior_analysis) {

  q <- k - p
  m <- n - p
  ff <- seq(0.01, 10, 0.05)

  # BF10 values and bounds
  BF10 <- F_BF(ff, q, m, dff, rscale, f_m, prior_analysis)
  f.BF10 <- F_BF_bound_10(threshold, q, m, dff, rscale, f_m, prior_analysis)

  df10 <- data.frame(f = ff, BF = BF10)

  main.bf10 <- if (length(f.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(threshold)~" when f = "~.(round(f.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(threshold)~" when f = "~.(round(f.BF10[1], 2))~" or "~.(round(f.BF10[2], 2))))
  }

  # x-axis breaks for BF10
  x_breaks_10 <- sort(unique(c(0, 10, round(f.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(f, BF)) +
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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
  f.BF01 <- F_BF_bound_01(threshold, q, m, dff, rscale, f_m, prior_analysis)
  max_BF01 <- 1 / F_BF(0.001, q, m, dff, rscale, f_m, prior_analysis)
  impossible <- any(max_BF01 < threshold | f.BF01 == "bound cannot be found")

  df01 <- data.frame(f = ff, BF = BF01)

  if (impossible) {
    p2 <- ggplot2::ggplot(df01, ggplot2::aes(f, BF)) +
      ggplot2::geom_line(linewidth = 1.2, color = "black") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(0, 10), breaks = c(0, 10)) +
      ggplot2::labs(
        x = "f-value",
        y = expression("BF"[0][1] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[0][1]~"="~.(threshold)))
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
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when f = "~.(round(f.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when f = "~.(round(f.BF01[1], 2))~" or "~.(round(f.BF01[2], 2))))
    }

    x_breaks_01 <- sort(unique(c(0, 10, round(f.BF01, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(f, BF)) +
      ggplot2::geom_line(linewidth = 1.2, color = "black") +
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

Power_f <- function(threshold, k, p, dff, rscale, f_m, prior_analysis,
                    k_d, p_d, dff_d, rscale_d, f_m_d, prior_design,
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
    f10 <- F_BF_bound_10(threshold, q, m[i], dff, rscale, f_m, prior_analysis)
    f01 <- F_BF_bound_01(threshold, q, m[i], dff, rscale, f_m, prior_analysis)

    # True Positive
    TPE[i] <- if (de_an_prior == 1) {
      F_TPE(f10, q, m[i], dff, rscale, f_m, prior_analysis)
    } else {
      F_TPE(f10, q, m[i], dff_d, rscale_d, f_m_d, prior_design)
    }

    # False Positive / True Negative
    FPE[i] <- F_FPE(f10, q, m[i])
    TNE[i] <- F_TNE(f01, q, m[i])

    # False Negative
    FNE[i] <- if (de_an_prior == 1) {
      F_FNE(f01, q, m[i], dff, rscale, f_m, prior_analysis)
    } else {
      F_FNE(f01, q, m[i], dff_d, rscale_d, f_m_d, prior_design)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # BF01 plot
  p2 <- ggplot2::ggplot(df_bf01, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots side by side
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

# ---- ANOVAe.r ----
Fe_BF <- function(f, q, m, dff, rscale, f_m, prior_analysis, ROPE) {
  # Compute normalizations once
  normalizationh1  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,prior_analysis),lower = ROPE,upper = Inf,rel.tol = 1e-10)$value
  normalizationh0 <- 1 - normalizationh1

  # Define likelihood ratio function
  sapply(f, function(fi) {
    int1 <- function(fsq) {
      stats::df(fi, q, m - q, ncp = m * fsq) * F_prior(fsq,q,dff,rscale,f_m,prior_analysis) / normalizationh1
    }
    int0 <- function(fsq) {
      stats::df(fi, q, m - q, ncp = m * fsq) * F_prior(fsq,q,dff,rscale,f_m,prior_analysis) / normalizationh0
    }
    lh1 <- stats::integrate(int1, lower = ROPE, upper = Inf, stop.on.error = FALSE)$value
    lh0 <- stats::integrate(int0, lower = 0, upper = ROPE,   stop.on.error = FALSE)$value
    lh1 / lh0
  })
}

Fe_BF_bound_10 <-function(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE){
  x = numeric(0)
  Bound_finding <-function(f){
    Fe_BF(f,q,m,dff,rscale,f_m,prior_analysis,ROPE)-threshold
  }
  #x = tryCatch( stats::uniroot(Bound_finding,lower=0.01,upper = 100 )$root, error=function(e){})
  x = tryCatch( rootSolve::uniroot.all(Bound_finding,lower=0.01,upper = 500 ), error=function(e){})

  if (length(x) == 0) return("no bound is found")

  return(x)
}

Fe_BF_bound_01 <-function(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE){
  Fe_BF_bound_10(1/threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
}

Fe_TPE<-function(f,q,m,dff,rscale,f_m,prior_analysis,ROPE){


  if (length(f) == 0 || any(f == "no bound is found")) return(0)


  if (prior_analysis == "Point") return(stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = F))

  normalizationh1  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,prior_analysis),lower = ROPE,upper = Inf,rel.tol = 1e-5)$value
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,prior_analysis)/normalizationh1
  }
  x = stats::integrate(int,lower = ROPE,upper = Inf)$value
  return(x)
}

Fe_FNE<-function(f,q,m,dff,rscale,f_m,prior_analysis,ROPE){

  if (length(f) == 0 || any(f == "no bound is found")) return(0)


  if (prior_analysis == "Point") return(stats::pf(f,q,m-q,ncp =m*f_m^2,lower.tail = T))

  normalizationh1  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,prior_analysis),lower = ROPE,upper = Inf,rel.tol = 1e-10)$value
  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,prior_analysis)/normalizationh1
  }
  x = stats::integrate(int,lower = ROPE,upper = Inf)$value
  return(x)
}

Fe_TNE<-function(f,q,m,dff,rscale,f_m,prior_analysis,ROPE){

  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  normalizationh0  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,prior_analysis),lower = 0,upper = ROPE,rel.tol = 1e-5)$value

  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = T)*F_prior(fsq,q,dff,rscale,f_m,prior_analysis)/normalizationh0
  }
  x = stats::integrate(int,lower = 0,upper = ROPE)$value
  return(x)
}

Fe_FPE<-function(f,q,m,dff,rscale,f_m,prior_analysis,ROPE){

  if (length(f) == 0 || any(f == "no bound is found")) return(0)

  normalizationh0  <- stats::integrate(function(fsq)F_prior(fsq,q,dff,rscale,f_m,prior_analysis),lower = 0,upper = ROPE,rel.tol = 1e-10)$value


  int  <- function(fsq){

    stats::pf(f,q,m-q,ncp =m*fsq,lower.tail = F)*F_prior(fsq,q,dff,rscale,f_m,prior_analysis)/normalizationh0
  }
  x = stats::integrate(int,lower = 0,upper = ROPE)$value
  return(x)
}

fe_N_finder<-function(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,dff_d,rscale_d,f_m_d,prior_design,de_an_prior,ROPE,false_rate ){
  q     <- k-p
  lower <- 2*k-p+1
  m     <- lower-p
  f     <- Fe_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)

  p2    <- if (de_an_prior == 1)
    Fe_TPE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE) else
      Fe_TPE(f,q,m,dff_d,rscale_d,f_m_d,prior_design,ROPE)
  if (p2 > true_rate) return(lower)

  Power_root <- function(n){
    m   <- n-p
    f   <- Fe_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
    pro <- if (de_an_prior == 1)
      Fe_TPE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE) else
        Fe_TPE(f,q,m,dff_d,rscale_d,f_m_d,prior_design,ROPE)
    return(pro-true_rate)
  }

  #N.power <- stats::uniroot(Power_root,lower = lower,upper =  5000)$root
  N.power <-robust_uniroot(Power_root,lower=lower)
  m       <- N.power-p
  f       <- Fe_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
  FPE     <- Fe_FPE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE)

  if (FPE <= false_rate) return(N.power)

  alpha_root <- function(n){
    m= n-p
    f = Fe_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
    pro = Fe_FPE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE)
    return(pro-false_rate)
  }
  N.alpha <- robust_uniroot(alpha_root,lower=N.power)
  #stats::uniroot(alpha_root,lower = N.power,upper =  5000)$root
  return(N.alpha)
}

fe_N_01_finder<-function(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,dff_d,rscale_d,f_m_d,prior_design,de_an_prior,ROPE,false_rate ){
  q     <- k-p
  lower <- 2*k-p+1
  m     <- lower-p
  upper <-  5000
  f     <- Fe_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
  TNE_lo <- Fe_TNE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE)
  FNE_lo <- if (de_an_prior == 1)
    Fe_FNE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE) else
      Fe_FNE(f,q,m,dff_d,rscale_d,f_m_d,prior_design,ROPE)

  if (TNE_lo > true_rate && FNE_lo < false_rate) {
    return(lower)
  } else if (TNE_lo > true_rate) {
    FN.root <- function(n) {
      q     <- k-p
      m     <- n-p
      f <- Fe_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
      FNE <- if (de_an_prior == 1)
        Fe_FNE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE) else
          Fe_FNE(f,q,m,dff_d,rscale_d,f_m_d,prior_design,ROPE)
      FNE - false_rate
    }
    return(stats::uniroot(FN.root, lower = lower, upper = upper)$root)
  }

  TN_root <- function(n){
    m   <- n-p
    f   <- Fe_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
    pro <- Fe_TNE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE)
    return(pro-true_rate)
  }

  #N.TN <- stats::uniroot(Power_root,lower = lower,upper =  5000)$root
  N.TN <-robust_uniroot(TN_root,lower=lower)
  m       <- N.TN-p
  f       <- Fe_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
  FNE     <- if (de_an_prior == 1)
    Fe_FNE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE) else
      Fe_FNE(f,q,m,dff_d,rscale_d,f_m_d,prior_design,ROPE)

  if (FNE <= false_rate) return(N.TN)

  FN.root <- function(n) {
    q     <- k-p
    m     <- n-p
    f <- Fe_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
    FNE <- if (de_an_prior == 1)
      Fe_FNE(f,q,m,dff,rscale,f_m,prior_analysis,ROPE) else
        Fe_FNE(f,q,m,dff_d,rscale_d,f_m_d,prior_design,ROPE)
    FNE - false_rate
  }
  N.FN <- robust_uniroot(FN.root,lower=N.TN)
  #stats::uniroot(alpha_root,lower = N.TN,upper =  5000)$root
  return(N.FN)
}


fe_table<-function(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,
                   dff_d,rscale_d,f_m_d,prior_design,de_an_prior,n, mode_bf,ROPE ,false_rate,type_rate){

  if (mode_bf == 1){

    n = switch(type_rate,
               "positive" = ceiling(fe_N_finder(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,dff_d,rscale_d,f_m_d,prior_design,de_an_prior,ROPE ,false_rate)),
               "negative" = ceiling(fe_N_01_finder(threshold,true_rate,p,k,dff,rscale,f_m,prior_analysis,dff_d,rscale_d,f_m_d,prior_design,de_an_prior,ROPE ,false_rate)))
  } else {
    n=n
  }
  q= k-p
  m= n-p

  # f bounds:
  f10 <- Fe_BF_bound_10(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)
  f01 <- Fe_BF_bound_01(threshold,q,m,dff,rscale,f_m,prior_analysis,ROPE)

  # max BF10 possible:
  max_BF <- 1/Fe_BF(0.00001,q,m,dff,rscale,f_m,prior_analysis,ROPE)
  BF_D   <- f10

  # FPE and TPE:
  FPE       <- Fe_FPE(f10,q,m,dff,rscale,f_m,prior_analysis,ROPE)
  if (de_an_prior == 1) {
    TPE       <- Fe_TPE(f10,q,m,dff,rscale,f_m,prior_analysis,ROPE)
    TPR_dff   <- dff
    TPR_rscale<- rscale
    TPR_f_m   <- f_m
    TPR_prior <- prior_analysis
  } else {
    TPE       <- Fe_TPE(f10,q,m,dff_d,rscale_d,f_m_d,prior_design,ROPE)
    TPR_dff   <- dff_d
    TPR_rscale<- rscale_d
    TPR_f_m   <- f_m_d
    TPR_prior <- prior_design
  }

  # FNE and TNE:
  if (max_BF < threshold | BF_D == "no bound is found") {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- Fe_FNE(f01,q,m,TPR_dff,TPR_rscale,TPR_f_m,TPR_prior,ROPE)
    TNE <-  Fe_TNE(f01,q,m,dff,rscale,f_m,prior_analysis,ROPE)

  }
  # table:
  tab.names <- c(
    "TruePositve",
    "FalseNegative",
    "TrueNegative",
    "FalsePositive",
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table

}


prior_plot_fe <- function(q, dff, rscale, f_m, prior_analysis,
                          dff_d, rscale_d, f_m_d, prior_design,
                          de_an_prior, ROPE) {

  fsq <- seq(0.01, 3, 0.01)

  prior_h1 <- F_prior(fsq, q, dff, rscale, f_m, prior_analysis)
  prior_h1[fsq < ROPE] <- 0

  prior_h0 <- F_prior(fsq, q, dff, rscale, f_m, prior_analysis)
  prior_h0[fsq > ROPE] <- 0

  df_lines <- data.frame(
    fsq = rep(fsq, 2),
    Density = c(prior_h1, prior_h0),
    Prior = rep(c("H1 - Analysis Prior", "H0 - Analysis Prior"),
                each = length(fsq))
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = df_lines,
                       ggplot2::aes(x = fsq, y = Density,
                                    color = Prior,
                                    linetype = Prior,
                                    linewidth = Prior)) +
    ggplot2::scale_color_manual(values = c(
      "H1 - Analysis Prior" = "black",
      "H0 - Analysis Prior" = "black",
      "H1 - Design Prior"   = "gray"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "H1 - Analysis Prior" = "solid",
      "H0 - Analysis Prior" = "dashed",
      "H1 - Design Prior"   = "solid"
    )) +
    ggplot2::scale_linewidth_manual(values = c(
      "H1 - Analysis Prior" = 1.2,
      "H0 - Analysis Prior" = 1.2,
      "H1 - Design Prior"   = 2
    )) +
    ggplot2::labs(
      x = expression(bold(lambda^2)),
      y = "Density",
      title = bquote(bold("Prior distribution on "~lambda^2~
                            " under the alternative hypothesis"))
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = c(0.80, 0.95),
      legend.justification = c("right", "top"),
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  # Non-point design prior
  if (de_an_prior == 0 && prior_design != "Point") {
    prior_design_vals <- F_prior(fsq, q, dff_d, rscale_d, f_m_d, prior_design)
    prior_design_vals[fsq < ROPE] <- 0

    df_design <- data.frame(
      fsq = fsq,
      Density = prior_design_vals,
      Prior = "H1 - Design Prior"
    )

    p <- p + ggplot2::geom_line(
      data = df_design,
      ggplot2::aes(x = fsq, y = Density,
                   color = Prior,
                   linetype = Prior,
                   linewidth = Prior)
    )
  }

  # Point design prior
  if (de_an_prior == 0 && prior_design == "Point") {
    ylim_max <- max(prior_h1, prior_h0, na.rm = TRUE)

    df_dummy <- data.frame(fsq = NA, Density = NA,
                           Prior = "H1 - Design Prior")

    p <- p +
      ggplot2::geom_line(data = df_dummy,
                         ggplot2::aes(x = fsq, y = Density,
                                      color = Prior,
                                      linetype = Prior,
                                      linewidth = Prior),
                         na.rm = TRUE) +
      ggplot2::geom_segment(
        ggplot2::aes(x = f_m_d, xend = f_m_d,
                     y = 0, yend = ylim_max),
        color = "gray",
        linetype = "dashed",
        arrow = ggplot2::arrow(length = grid::unit(0.1, "inches"))
      )
  }

  return(p)
}
bf10_fe <- function(threshold, n, k, p, dff, rscale, f_m, prior_analysis, ROPE) {

  q <- k - p
  m <- n - p
  ff <- seq(0.01, 10, 0.05)

  # BF10 values and bounds
  BF10 <- Fe_BF(ff, q, m, dff, rscale, f_m, prior_analysis, ROPE)
  f.BF10 <- Fe_BF_bound_10(threshold, q, m, dff, rscale, f_m, prior_analysis, ROPE)

  df10 <- data.frame(f = ff, BF = BF10)

  main.bf10 <- if (length(f.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(threshold)~" when f = "~.(round(f.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(threshold)~" when f = "~.(round(f.BF10[1], 2))~" or "~.(round(f.BF10[2], 2))))
  }

  # x-axis breaks for BF10
  x_breaks_10 <- sort(unique(c(0, 10, round(f.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(f, BF)) +
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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
  f.BF01 <- Fe_BF_bound_01(threshold, q, m, dff, rscale, f_m, prior_analysis, ROPE)
  max_BF01 <- 1 / Fe_BF(0.001, q, m, dff, rscale, f_m, prior_analysis, ROPE)
  impossible <- any(max_BF01 < threshold | f.BF01 == "bound cannot be found")

  df01 <- data.frame(f = ff, BF = BF01)

  if (impossible) {
    p2 <- ggplot2::ggplot(df01, ggplot2::aes(f, BF)) +
      ggplot2::geom_line(linewidth = 1.2, color = "black") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(0, 10), breaks = c(0, 10)) +
      ggplot2::labs(
        x = "f-value",
        y = expression("BF"[0][1] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[0][1]~"="~.(threshold)))
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
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when f = "~.(round(f.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when f = "~.(round(f.BF01[1], 2))~" or "~.(round(f.BF01[2], 2))))
    }

    x_breaks_01 <- sort(unique(c(0, 10, round(f.BF01, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(f, BF)) +
      ggplot2::geom_line(linewidth = 1.2, color = "black") +
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


Power_fe <- function(threshold, k, p, dff, rscale, f_m, prior_analysis,
                     k_d, p_d, dff_d, rscale_d, f_m_d, prior_design,
                     de_an_prior, N, ROPE) {

  # Sample size range
  smin <- (2 * k - p + 1)
  smax <- N * 2
  sdf <- seq(smin, smax, length.out = 31)
  q <- k - p
  m <- sdf - p

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(sdf))

  for (i in seq_along(sdf)) {
    f10 <- Fe_BF_bound_10(threshold, q, m[i], dff, rscale, f_m, prior_analysis, ROPE)
    f01 <- Fe_BF_bound_01(threshold, q, m[i], dff, rscale, f_m, prior_analysis, ROPE)

    # True Positive
    TPE[i] <- if (de_an_prior == 1) {
      Fe_TPE(f10, q, m[i], dff, rscale, f_m, prior_analysis, ROPE)
    } else {
      Fe_TPE(f10, q, m[i], dff_d, rscale_d, f_m_d, prior_design, ROPE)
    }

    # False Negative
    FNE[i] <- if (de_an_prior == 1) {
      Fe_FNE(f01, q, m[i], dff, rscale, f_m, prior_analysis, ROPE)
    } else {
      Fe_FNE(f01, q, m[i], dff_d, rscale_d, f_m_d, prior_design, ROPE)
    }

    # False Positive / True Negative
    FPE[i] <- Fe_FPE(f10, q, m[i], dff, rscale, f_m, prior_analysis, ROPE)
    TNE[i] <- Fe_TNE(f01, q, m[i], dff, rscale, f_m, prior_analysis, ROPE)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # BF01 plot
  p2 <- ggplot2::ggplot(df_bf01, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots side by side
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

# ---- binomial.r ----
adjust_root_10 <- function(root, n, alpha, beta, location, scale, prior_analysis, alternative, threshold) {
  # If root is less than 0, return NA
  if (root < 0) return(NA)

  # Evaluate BF at the root
  BF_val <- bin_BF(root, n, alpha, beta, location, scale, prior_analysis, alternative)

  if (BF_val <= threshold) {
    # Try root - 1 only if root > 0
    if (root > 0) {
      BF_prev <- bin_BF(root - 1, n, alpha, beta, location, scale, prior_analysis, alternative)
      if (BF_prev > threshold) return(root - 1)
    }

    # Try root + 1
    BF_next <- bin_BF(root + 1, n, alpha, beta, location, scale, prior_analysis, alternative)
    if (BF_next > threshold) return(root + 1)
  }

  # Return original if already valid or no better nearby found
  return(root)
}


adjust_root_01 <- function(root, n, alpha, beta, location, scale, prior_analysis, alternative, threshold) {
  # Evaluate BF at the root
  BF_val <- 1/bin_BF(root, n, alpha, beta, location, scale, prior_analysis, alternative)

  if (BF_val <= threshold) {
    # Try root - 1
    BF_prev <- 1/bin_BF(root - 1, n, alpha, beta, location, scale, prior_analysis, alternative)
    if (BF_prev > threshold) return(root - 1)

    # Try root + 1
    BF_next <- 1/bin_BF(root + 1, n, alpha, beta, location, scale, prior_analysis, alternative)
    if (BF_next > threshold) return(root + 1)
  }

  # Return original if already valid or no better nearby found
  return(root)
}


bin_prior <-function(prop,alpha,beta,location,scale,prior_analysis){

  switch(prior_analysis,
         "beta" = stats::dbeta(prop, alpha,beta),
         "Moment" = dMoment(prop,location,scale))
}
bin_BF<-function(x,n,alpha,beta,location,scale,prior_analysis,alternative){
  BF = NA
  bound  <- switch(alternative,
                   "greater" = c(a = location, b = 1),
                   "less" = c(a = 0, b = location),
                   "two.sided" = c(a = 0, b = 1)
  )


  normalization <- if (alternative == "two.sided") {
    switch(prior_analysis,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(prior_analysis,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
  }
  for( i in 1:length(x)){
    int  <- function(prop){stats::dbinom(x[i], size=n, prob=prop) *bin_prior(prop,alpha,beta,location,scale,prior_analysis)}
    lh1 <- stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = 1e-5)$value / normalization
    lh0 <- stats::dbinom(x[i], size = n, prob = location)
    BF[i] = lh1 / lh0
  }


  return(BF)

}

bin_BF_bound_10 <-function(threshold,n,alpha,beta,location,scale,prior_analysis,alternative){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    bin_BF(x,n,alpha,beta,location,scale,prior_analysis,alternative)- threshold
  }

  x <- tryCatch(stats::uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)
  #results <- tryCatch(uniroot.all(Bound_finding, lower = 0 ,upper = n), error = function(e) NA)
  results <- round(results[!is.na(results) & is.finite(results)])

  if (length(results) == 0) return("bound cannot be found")


  results <- sapply(results, function(root) {
    adjust_root_10(root, n, alpha, beta, location, scale, prior_analysis, alternative, threshold)
  })


  BF.vals  <- bin_BF(results,n,alpha,beta,location,scale,prior_analysis,alternative)

  BF.close <- which(BF.vals > threshold)
  if (length(BF.close) == 0 || all(!is.finite(BF.close))) return("bound cannot be found")
  return(results[BF.close])
}

bin_BF_bound_01 <-function(threshold,n,alpha,beta,location,scale,prior_analysis,alternative){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    1/bin_BF(x,n,alpha,beta,location,scale,prior_analysis,alternative)- threshold
  }

  x <- tryCatch(stats::uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)

  results <- round(results[!is.na(results)])
  if (length(results) == 0) return("bound cannot be found")


  results <- sapply(results, function(root) {
    adjust_root_01(root, n, alpha, beta, location, scale, prior_analysis, alternative, threshold)
  })


  BF.vals  <- 1/bin_BF(results,n,alpha,beta,location,scale,prior_analysis,alternative)

  BF.close <- which(BF.vals > threshold)
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}


bin_TPE<-function(x,n,h0,alpha,beta,location,scale,prior_analysis,alternative){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  if (prior_analysis =="Point"){
    TPE = switch(alternative,
                 "two.sided" = {

                   switch(length(x)==2,
                          "1" ={stats::pbinom(min(x),n,location,lower.tail = T)+ stats::pbinom(max(x)-1,n,location,lower.tail = F)},
                          "0"=  {
                            switch(x/n>location,
                                   "1" = stats::pbinom(x-1,n,location,lower.tail = F),
                                   "0" = stats::pbinom(x,n,location,lower.tail = T))

                          })
                 },
                 "greater"  = {stats::pbinom(x-1,n,location,lower.tail = F)},
                 "less"  = {stats::pbinom(x,n,location,lower.tail = T)}
    )
    return(TPE)
  }

  bound  <- switch(alternative,
                   "greater" = c(a = h0, b = 1),
                   "less" = c(a = 0, b = h0),
                   "two.sided" = c(a = 0, b = 1)
  )
  normalization <- if (alternative == "two.sided") {
    switch(prior_analysis,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(prior_analysis,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
  }
  int <- function(prop) {
    pro <- switch(alternative,
                  "two.sided" = {
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
                  "greater" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "less" = stats::pbinom(x, n, prop, lower.tail = TRUE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, prior_analysis) / normalization
  }

  TPE = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-5)$value

  return(TPE)

}

bin_FNE<-function(x,n,h0,alpha,beta,location,scale,prior_analysis,alternative){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  if (prior_analysis == "Point") {
    FNE <- switch(alternative,
                  "two.sided" = {
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
                  "greater" = stats::pbinom(x, n, location, lower.tail = TRUE),
                  "less" = stats::pbinom(x - 1, n, location, lower.tail = FALSE)
    )
    return(FNE)
  }



  bound  <- switch(alternative,
                   "greater" = c(a = h0, b = 1),
                   "less" = c(a = 0, b = h0),
                   "two.sided" = c(a = 0, b = 1)
  )

  normalization <- if (alternative == "two.sided") {
    switch(prior_analysis,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(prior_analysis,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
  }
  int <- function(prop) {
    pro <- switch(alternative,
                  "two.sided" = {
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
                  "greater" = stats::pbinom(x , n, prop, lower.tail = TRUE),
                  "less" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, prior_analysis) / normalization
  }
  FNE = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-5)$value
  return(FNE)

}

bin_FPE<-function(x,n,location,alternative){

  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  FPE <- switch(alternative,
                "two.sided" = {
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
                "greater" = stats::pbinom(x - 1, n, location, lower.tail = FALSE),
                "less" = stats::pbinom(x, n, location, lower.tail = TRUE)
  )

  return(FPE)

}

bin_TNE<-function(x,n,location,alternative){

  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)


  TNE <- switch(alternative,
                "two.sided" = {
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
                "greater" = stats::pbinom(x, n, location, lower.tail = TRUE),
                "less" = stats::pbinom(x - 1, n, location, lower.tail = FALSE)
  )

  return(TNE)

}

bin_N_finder <-function(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                        alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,false_rate){
  lower = 10
  upper = 10000

  b10 = bin_BF_bound_10(threshold,lower,alpha,beta,location,scale,prior_analysis,alternative)
  TPE_lo <- if (de_an_prior == 1)
    bin_TPE(b10,lower,h0,alpha,beta,location,scale,prior_analysis,alternative) else
      bin_TPE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative)
  if (TPE_lo > true_rate) return(lower)
  FPE_lo <-  bin_FPE(b10,lower,location,alternative)
  if (TPE_lo > true_rate&FPE_lo<false_rate) return(lower)


  Power_root <- function(N){
    N =round(N)
    b10 = bin_BF_bound_10 (threshold,N,alpha,beta,location,scale,prior_analysis,alternative)
    pro <- if (de_an_prior==0){
      bin_TPE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative)
    }else bin_TPE(b10,N,h0,alpha,beta,location,scale,prior_analysis,alternative)

    pro-true_rate
  }

  N.power = round(stats::uniroot(Power_root,lower = lower,upper = upper)$root)+1

  N.extended = seq(N.power,N.power+20,2)
  Power.extended = unlist(lapply(N.extended, Power_root))

  if (any(Power.extended<0)){
    lower = which(Power.extended < 0)[1]
    N.power = round(stats::uniroot(Power_root,lower = lower,upper = upper)$root)+1

  }



  while(TRUE) {
    b10 <- bin_BF_bound_10(threshold, N.power, alpha, beta, location, scale, prior_analysis, alternative)
    pro <- if (de_an_prior == 0) {
      bin_TPE(b10, N.power,h0, alpha_d, beta_d, location_d, scale_d, prior_design, alternative)
    } else {
      bin_TPE(b10, N.power,h0, alpha, beta, location, scale, prior_analysis, alternative)
    }

    if (pro > true_rate) break
    N.power <- N.power + 1
  }


  b10 = bin_BF_bound_10(threshold,N.power,alpha,beta,location,scale,prior_analysis,alternative)
  FPE = bin_FPE(b10,N.power,location,alternative)
  if (FPE <= false_rate) return(N.power)


  alpha.root <- function(n) {
    n=round(n)
    b10 <- bin_BF_bound_10 (threshold,n,alpha,beta,location,scale,prior_analysis,alternative)
    bin_FPE(b10,n,location,alternative)-false_rate
  }
  N.alpha = round(stats::uniroot(alpha.root,lower = N.power,upper = upper)$root)
  return(N.alpha)
}

bin_N_01_finder <-function(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                           alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,false_rate){
  lower = 10
  upper = 10000

  b10 = bin_BF_bound_01(threshold,lower,alpha,beta,location,scale,prior_analysis,alternative)
  TNE_lo = bin_TNE(b10,lower,location,alternative)
  FNE_lo <- if (de_an_prior == 1)
    bin_FNE(b10,lower,h0,alpha,beta,location,scale,prior_analysis,alternative) else
      bin_FNE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative)

  if (TNE_lo > true_rate && FNE_lo < false_rate) {
    return(lower)
  } else if (TNE_lo > true_rate) {
    FN_root <- function(N){
      N =round(N)
      b10 = bin_BF_bound_01 (threshold,N,alpha,beta,location,scale,prior_analysis,alternative)
      pro <- if (de_an_prior == 1)
        bin_FNE(b10,N,h0,alpha,beta,location,scale,prior_analysis,alternative) else
          bin_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative)

      pro-false_rate
    }
    return(round(stats::uniroot(FN_root, lower = lower, upper = upper)$root))
  }
  TN_root <- function(N){
    N =round(N)
    b10 = bin_BF_bound_01 (threshold,N,alpha,beta,location,scale,prior_analysis,alternative)
    pro <-  bin_TNE(b10,N,location,alternative)

    pro-true_rate
  }

  N.TN = round(stats::uniroot(TN_root,lower = lower,upper = upper)$root)+1
  N.extended = seq(N.TN,N.TN+20,2)
  Power.extended = unlist(lapply(N.extended, TN_root))

  if (any(Power.extended<0)){
    lower = which(Power.extended < 0)[1]
    N.TN = round(stats::uniroot(TN_root,lower = lower,upper = upper)$root)+1

  }



  while(TRUE) {
    b10 <- bin_BF_bound_01(threshold, N.TN, alpha, beta, location, scale, prior_analysis, alternative)
    pro <- bin_TNE(b10,N.TN,location,alternative)

    if (pro > true_rate) break
    N.TN <- N.TN + 1
  }


  b10 = bin_BF_bound_01(threshold,N.TN,alpha,beta,location,scale,prior_analysis,alternative)
  FNE = if (de_an_prior == 1)
    bin_FNE(b10,N.TN,h0,alpha,beta,location,scale,prior_analysis,alternative) else
      bin_FNE(b10,N.TN,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative)
  if (FNE <= false_rate) return(N.TN)
  FN_root <- function(N){
    N =round(N)
    b10 = bin_BF_bound_01 (threshold,N,alpha,beta,location,scale,prior_analysis,alternative)
    pro <- if (de_an_prior == 1)
      bin_FNE(b10,N,h0,alpha,beta,location,scale,prior_analysis,alternative) else
        bin_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative)

    pro-false_rate
  }
  N.FN = round(stats::uniroot(FN_root,lower = N.TN,upper = upper)$root)
  return(N.FN)
}

bin_table<-function(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                    alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,N, mode_bf,false_rate,type_rate){
  if (mode_bf == "0") n = N else n = switch(
    type_rate,
    "positive" = bin_N_finder(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                              alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,false_rate),
    "negative" = bin_N_01_finder(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                                 alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,false_rate))


  # b bounds:
  b10 <- bin_BF_bound_10(threshold,n,alpha,beta,location,scale,prior_analysis,alternative)
  b01 <-  bin_BF_bound_01(threshold,n,alpha,beta,location,scale,prior_analysis,alternative)


  # max BF10 possible:
  max_BF <- 1 / bin_BF(round(location*n),n,alpha,beta,location,scale,prior_analysis,alternative)
  BF_D   <- b10

  # FPE and TPE:
  FPE       <- bin_FPE(b10,n,location,alternative)
  if (de_an_prior == 1) {
    TPE          <- bin_TPE(b10,n,h0,alpha,beta,location,scale,prior_analysis,alternative)
    TPR_alpha    <- alpha
    TPR_beta     <- beta
    TPR_location <- location
    TPR_scale    <- scale
    TPR_prior    <- prior_analysis

  } else {
    TPE          <- bin_TPE(b10,n,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative)
    TPR_alpha    <- alpha_d
    TPR_beta     <- beta_d
    TPR_location <- location_d
    TPR_scale    <- scale_d
    TPR_prior    <- prior_design
  }


  # FNE and TNE:
  if (any(alternative == "two.sided" & max_BF < threshold | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- bin_FNE(b01,n,h0,TPR_alpha,TPR_beta,TPR_location,TPR_scale,TPR_prior,alternative)
    TNE <- bin_TNE(b01,n,location,alternative)
  }

  # table:
  tab.names <- c(
    "TruePositve",
    "FalseNegative",
    "TrueNegative",
    "FalsePositive",
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}



bin_bf10<- function(threshold, n, alpha, beta, location, scale, prior_analysis, alternative) {

  # Sequence of successes
  x <- seq(0, n, by = 3)

  # Compute BF10 and bounds
  BF10 <- bin_BF(x, n, alpha, beta, location, scale, prior_analysis, alternative)
  b.BF10 <- bin_BF_bound_10(threshold, n, alpha, beta, location, scale, prior_analysis, alternative)
  BF10_at_b <- bin_BF(b.BF10, n, alpha, beta, location, scale, prior_analysis, alternative)

  BF01 <- 1 / BF10
  b.BF01 <- bin_BF_bound_01(threshold, n, alpha, beta, location, scale, prior_analysis, alternative)
  BF01_at_b <- 1 / bin_BF(b.BF01, n, alpha, beta, location, scale, prior_analysis, alternative)

  # Check if BF01 = D is impossible
  max.BF01 <- 1 / bin_BF(round(location * n), n, alpha, beta, location, scale, prior_analysis, alternative)
  impossible <- (alternative == "two.sided") && (max.BF01 < threshold || identical(b.BF01, "bound cannot be found"))
  # Titles for BF10
  main.bf10 <- if (length(b.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(round(BF10_at_b, 2)) ~ " when x = " ~ .(round(b.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(round(BF10_at_b[1], 2)) ~ "/" ~ .(round(BF10_at_b[2], 2)) ~
                  " when x = " ~ .(round(b.BF10[1], 2)) ~ " or " ~ .(round(b.BF10[2], 2))))
  }

  # Titles for BF01
  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(threshold)))
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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

Power_bin <- function(threshold, h0, alpha, beta, location, scale, prior_analysis, alternative,
                      alpha_d, beta_d, location_d, scale_d, prior_design,
                      de_an_prior, N) {

  # Sample size range
  Ns <- ceiling(seq(10, N*1.2, length.out = 31))

  # Initialize probability vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(Ns))

  # Compute bounds and probabilities
  for (i in seq_along(Ns)) {
    x10 <- bin_BF_bound_10(threshold, Ns[i], alpha, beta, location, scale, prior_analysis, alternative)
    x01 <- bin_BF_bound_01(threshold, Ns[i], alpha, beta, location, scale, prior_analysis, alternative)

    TPE[i] <- if (de_an_prior == 1) {
      bin_TPE(x10, Ns[i], h0, alpha, beta, location, scale, prior_analysis, alternative)
    } else {
      bin_TPE(x10, Ns[i], h0, alpha_d, beta_d, location_d, scale_d, prior_design, alternative)
    }

    FNE[i] <- if (de_an_prior == 1) {
      bin_FNE(x01, Ns[i], h0, alpha, beta, location, scale, prior_analysis, alternative)
    } else {
      bin_FNE(x01, Ns[i], h0, alpha_d, beta_d, location_d, scale_d, prior_design, alternative)
    }

    FPE[i] <- bin_FPE(x10, Ns[i], location, alternative)
    TNE[i] <- bin_TNE(x01, Ns[i], location, alternative)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    axis_theme +
    legend_theme

  # ---------- BF01 Plot ----------
  p2 <- ggplot2::ggplot(df_BF01, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
    ) +
    axis_theme +
    legend_theme

  # Combine side-by-side
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}


compute.prior.density.b <- function(prop,alpha,beta,location,scale,prior_analysis,alternative) {
  if (prior_analysis == "Point") return(rep(NA, length(prop)))
  bound  <- switch(alternative,
                   "greater" = c(a = location, b = 1),
                   "less" = c(a = 0, b = location),
                   "two.sided" = c(a = 0, b = 1)
  )
  normalization <- if (alternative == "two.sided") {
    switch(prior_analysis,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = 1)

  } else {
    switch(prior_analysis,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "beta"     = stats::pbeta(bound[2],alpha,beta)-stats::pbeta(bound[1],alpha,beta))
  }
  bin_prior(prop,alpha,beta,location,scale,prior_analysis)/ normalization
}


bin_prior_plot <- function(h0, alpha, beta,
                           location, scale, prior_analysis,
                           alpha_d, beta_d,
                           location_d, scale_d,
                           prior_design,
                           alternative, de_an_prior) {

  # ---- Determine bounds ----
  plot.bounds <- switch(alternative,
                        "greater"   = c(h0, 1),
                        "less"      = c(0, h0),
                        "two.sided" = c(0, 1))

  prop <- seq(plot.bounds[1], plot.bounds[2], 0.01)

  # ---- Compute analysis prior ----
  prior.analysis <- compute.prior.density.b(
    prop, alpha, beta, location, scale,
    prior_analysis, alternative
  )

  # ---- Base data frame ----
  df <- data.frame(
    theta = prop,
    Density = prior.analysis,
    Prior = "H1 - Analysis Prior"
  )

  # ---- Add design prior if needed ----
  if (de_an_prior == 0) {

    if (prior_design == "Point") {

      df_design <- data.frame(
        theta = c(NA, NA),
        Density = c(NA, NA),
        Prior = "H1 - Design Prior"
      )

      df <- rbind(df, df_design)

    } else {

      prior.design <- compute.prior.density.b(
        prop, alpha_d, beta_d,
        location_d, scale_d,
        prior_design, alternative
      )

      df_design <- data.frame(
        theta = prop,
        Density = prior.design,
        Prior = "H1 - Design Prior"
      )

      df <- rbind(df, df_design)
    }
  }

  # ---- Y limits ----
  ylim_max <- max(df$Density[is.finite(df$Density)], na.rm = TRUE)

  # ---- Legend position (match t1_prior_plot) ----
  legend_pos <- switch(alternative,
                       "greater"   = c(0.65, 0.95),
                       "two.sided" = c(0.65, 0.95),
                       "less"      = c(0.05, 0.95))

  # ---- Build ggplot ----
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = theta,
                                    y = Density,
                                    color = Prior,
                                    linetype = Prior)) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(
      "H1 - Analysis Prior" = "black",
      "H1 - Design Prior"   = "gray"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "H1 - Analysis Prior" = "solid",
      "H1 - Design Prior"   = "dashed"
    )) +
    ggplot2::labs(
      x = expression(bold(theta)),
      y = "density",
      title = bquote(bold("Prior distribution on "~theta~" under the alternative"))
    ) +
    ggplot2::coord_cartesian(
      xlim = plot.bounds,
      ylim = c(0, ylim_max)
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.justification = c(0, 1),
      legend.background =
        ggplot2::element_rect(fill = scales::alpha("white", 0.8), color = NA),
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 1.5))
    )

  # ---- Handle Point prior ----
  if (de_an_prior == 0 && prior_design == "Point") {

    p <- p +
      ggplot2::annotate("segment",
                        x = location_d,
                        xend = location_d,
                        y = 0,
                        yend = ylim_max,
                        color = "gray",
                        linetype = "dashed",
                        arrow = ggplot2::arrow(
                          length = grid::unit(0.1, "inches")
                        ))
  }

  return(p)
}

# ---- binomial_e.r ----
adjust_root_10_e <- function(root, n, alpha, beta, location, scale, prior_analysis, alternative, threshold,ROPE) {
  # Evaluate BF at the root
  BF_val <- bin_e_BF(root,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)

  if (BF_val <=  threshold) {
    # Try root - 1
    BF_prev <- bin_e_BF(root-1,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    if (BF_prev > threshold) return(root - 1)

    # Try root + 1
    BF_next <- bin_e_BF(root+1,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    if (BF_next > threshold) return(root + 1)
  }

  # Return original if already valid or no better nearby found
  return(root)
}

adjust_root_01_e <- function(root, n, alpha, beta, location, scale, prior_analysis, alternative, threshold,ROPE) {
  # Evaluate BF at the root
  BF_val <- 1/bin_e_BF(root,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)

  if (BF_val <=  threshold) {
    # Try root - 1
    BF_prev <- 1/bin_e_BF(root-1,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    if (!is.nan(BF_prev) && !is.na(BF_prev) && BF_prev > threshold) return(root - 1)

    # Try root + 1
    BF_next <- 1/bin_e_BF(root+1,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    if (!is.nan(BF_next) && !is.na(BF_next) &&BF_next > threshold) return(root + 1)
  }

  # Return original if already valid or no better nearby found
  return(root)
}

bin_e_BF<-function(x,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE){
  BF = NA
  bound_h1  <- switch(alternative,
                      "greater" = c(a = location+ROPE, b = 1),
                      "less" = c(a = 0, b = location+ROPE),
                      "two.sided" = c(a = location+ROPE[1], b = location+ROPE[2])
  )
  bound_h0  <- switch(alternative,
                      "greater" = c(a = location, b = location+ROPE),
                      "less" = c(a = location+ROPE, b = location),
                      "two.sided" = c(a = location+ROPE[1], b = location+ROPE[2])
  )

  normalizationh1 <- switch(alternative,
                            "two.sided" = {
                              if (prior_analysis == "beta") {
                                1 - (stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta))
                              } else if (prior_analysis == "Moment") {
                                (pmom(1 - location, tau = scale^2) - pmom(bound_h1[2] - location, tau = scale^2)) +
                                  (pmom(bound_h1[1] - location, tau = scale^2) - pmom(0 - location, tau = scale^2))
                              }
                            },
                            "less" = ,
                            "greater" = {
                              if (prior_analysis == "beta") {
                                stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta)
                              } else if (prior_analysis == "Moment") {
                                pmom(bound_h1[2] - location, tau = scale^2) - pmom(bound_h1[1] - location, tau = scale^2)
                              }
                            }
  )

  normalizationh0 <- switch(prior_analysis,
                            "beta"      =   stats::pbeta(bound_h0[2], alpha, beta) - stats::pbeta(bound_h0[1], alpha, beta),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )

  for (i in 1:length(x)){
    int  <- function(prop){stats::dbinom(x[i], size=n, prob=prop) *bin_prior(prop,alpha,beta,location,scale,prior_analysis)
    }

    if (alternative == "two.sided"){
      lh1 = stats::integrate(int,lower = 0,upper = bound_h1[1], rel.tol=1e-5,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-5,stop.on.error = F)$value
    }else{
      lh1 = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-5,stop.on.error = F)$value

    }
    lh0 = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value


    BF[i] = (lh1/normalizationh1)/(lh0/normalizationh0)
  }
  return(BF)

}


bin_e_BF_bound_10 <-function(threshold,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    bin_e_BF(x,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)- threshold
  }
  x <- tryCatch(stats::uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)

  results <- round(results[!is.na(results)])
  if (length(results) == 0) return("bound cannot be found")


  results <- sapply(results, function(root) {
    adjust_root_10_e(root, n, alpha, beta, location, scale, prior_analysis, alternative, threshold,ROPE)
  })

  BF.vals  <- bin_e_BF(results,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  BF.close <- which(BF.vals > threshold)
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

bin_e_BF_bound_01 <-function(threshold,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE){
  y =x= numeric(0)
  Bound_finding <-function(x){
    x = round(x)
    1/bin_e_BF(x,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)- threshold
  }

  x <- tryCatch(stats::uniroot(Bound_finding, lower = 0 ,upper = round(location*n))$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower = round(location*n) ,upper = n)$root, error = function(e) NA)
  results <- c(x, y)

  results <- round(results[!is.na(results)])
  if (length(results) == 0) return("bound cannot be found")


  results <- sapply(results, function(root) {
    adjust_root_01_e(root, n, alpha, beta, location, scale, prior_analysis, alternative, threshold,ROPE)
  })

  BF.vals  <- 1/bin_e_BF(results,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  BF.close <- which(BF.vals > threshold)
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}


bin_e_TPE<-function(x,n,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE){
  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)


  if (prior_analysis =="Point"){
    TPE = switch(alternative,
                 "two.sided" = {

                   switch(length(x)==2,
                          "1" ={stats::pbinom(min(x),n,location,lower.tail = T)+ stats::pbinom(max(x)-1,n,location,lower.tail = F)},
                          "0"=  {
                            switch(x/n>location,
                                   "1" = stats::pbinom(x-1,n,location,lower.tail = F),
                                   "0" = stats::pbinom(x,n,location,lower.tail = T))

                          })
                 },
                 "greater"  = {stats::pbinom(x-1,n,location,lower.tail = F)},
                 "less"  = {stats::pbinom(x,n,location,lower.tail = T)}
    )
    return(TPE)
  }

  bound_h1  <- switch(alternative,
                      "greater" = c(a = h0+ROPE, b = 1),
                      "less" = c(a = 0, b = h0+ROPE),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )
  normalizationh1 <- switch(alternative,
                            "two.sided" = {
                              if (prior_analysis == "beta") {
                                1 - (stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta))
                              } else if (prior_analysis == "Moment") {
                                (pmom(1 - location, tau = scale^2) - pmom(bound_h1[2] - location, tau = scale^2)) +
                                  (pmom(bound_h1[1] - location, tau = scale^2) - pmom(0 - location, tau = scale^2))
                              }
                            },
                            "less" = ,
                            "greater" = {
                              if (prior_analysis == "beta") {
                                stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta)
                              } else if (prior_analysis == "Moment") {
                                pmom(bound_h1[2] - location, tau = scale^2) - pmom(bound_h1[1] - location, tau = scale^2)
                              }
                            }
  )
  int <- function(prop) {
    pro <- switch(alternative,
                  "two.sided" = {
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
                  "greater" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "less" = stats::pbinom(x, n, prop, lower.tail = TRUE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, prior_analysis) / normalizationh1
  }

  if(alternative == "two.sided"){
    TPE = stats::integrate(int,lower = 0,upper = bound_h1[1], rel.tol = 1e-5)$value + stats::integrate(int,lower = bound_h1[2],upper = 1, rel.tol = 1e-5)$value
  }else{
    TPE = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol = 1e-5)$value
  }
  return(TPE)

}


bin_e_FNE<-function(x,n,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE){

  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)


  if (prior_analysis =="Point"){
    FNE = switch(alternative,
                 "two.sided" = {

                   switch(length(x)==2,
                          "1" ={stats::pbinom(max(x),n,location,lower.tail = T)- stats::pbinom(min(x)-1,n,location,lower.tail = T)},
                          "0"=  {
                            switch(x/n>location,
                                   "1" = stats::pbinom(x,n,location,lower.tail = T),
                                   "0" = stats::pbinom(x-1,n,location,lower.tail = F))

                          })},
                 "greater"  = {stats::pbinom(x,n,location,lower.tail = T)},
                 "less"  = {stats::pbinom(x-1,n,location,lower.tail = F)}
    )
    return(FNE)
  }


  bound_h1  <- switch(alternative,
                      "greater" = c(a = h0+ROPE, b = 1),
                      "less" = c(a = 0, b = h0+ROPE),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )
  normalizationh1 <- switch(alternative,
                            "two.sided" = {
                              if (prior_analysis == "beta") {
                                1 - (stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta))
                              } else if (prior_analysis == "Moment") {
                                (pmom(1 - location, tau = scale^2) - pmom(bound_h1[2] - location, tau = scale^2)) +
                                  (pmom(bound_h1[1] - location, tau = scale^2) - pmom(0 - location, tau = scale^2))
                              }
                            },
                            "less" = ,
                            "greater" = {
                              if (prior_analysis == "beta") {
                                stats::pbeta(bound_h1[2], alpha, beta) - stats::pbeta(bound_h1[1], alpha, beta)
                              } else if (prior_analysis == "Moment") {
                                pmom(bound_h1[2] - location, tau = scale^2) - pmom(bound_h1[1] - location, tau = scale^2)
                              }
                            }
  )
  int <- function(prop) {
    pro <- switch(alternative,
                  "two.sided" = {
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
                  "greater" = stats::pbinom(x , n, prop, lower.tail = TRUE),
                  "less" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, prior_analysis) / normalizationh1
  }
  if(alternative == "two.sided"){
    FNE = stats::integrate(int,lower = 0,upper = bound_h1[1], rel.tol = 1e-5)$value + stats::integrate(int,lower = bound_h1[2],upper = 1, rel.tol = 1e-5)$value
  }else{
    FNE = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol = 1e-5)$value
  }


  return(FNE)

}

bin_e_FPE<-function(x,n,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE){

  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)


  bound_h0  <- switch(alternative,
                      "greater" = c(a = h0, b = h0+ROPE),
                      "less" = c(a = h0+ROPE, b = h0),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )
  normalizationh0 <- switch(prior_analysis,
                            "beta"      =   stats::pbeta(bound_h0[2], alpha, beta) - stats::pbeta(bound_h0[1], alpha, beta),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )

  int <- function(prop) {
    pro <- switch(alternative,
                  "two.sided" = {
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
                  "greater" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE),
                  "less" = stats::pbinom(x, n, prop, lower.tail = TRUE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, prior_analysis) / normalizationh0
  }


  FPE = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol = 1e-5)$value
  return(FPE)

}

bin_e_TNE<-function(x,n,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE){


  if (length(x) == 0 || any(x == "bound cannot be found")) return(0)

  bound_h0  <- switch(alternative,
                      "greater" = c(a = h0, b = h0+ROPE),
                      "less" = c(a = h0+ROPE, b = h0),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )

  normalizationh0 <- switch(prior_analysis,
                            "beta"      =   stats::pbeta(bound_h0[2], alpha, beta) - stats::pbeta(bound_h0[1], alpha, beta),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )
  int <- function(prop) {
    pro <- switch(alternative,
                  "two.sided" = {
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
                  "greater" = stats::pbinom(x , n, prop, lower.tail = TRUE),
                  "less" = stats::pbinom(x - 1, n, prop, lower.tail = FALSE)
    )

    pro * bin_prior(prop, alpha, beta, location, scale, prior_analysis) / normalizationh0
  }

  TNE = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol = 1e-5)$value



  return(TNE)
}

bin_e_N_finder <-function(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                          alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,false_rate,ROPE){
  lower = 10
  upper = 10000

  b10 =  bin_e_BF_bound_10(threshold,lower,alpha,beta,location,scale,prior_analysis,alternative,ROPE)

  TPE_lo <- if (de_an_prior == 1)
    bin_e_TPE(b10,lower,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE) else
      bin_e_TPE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE)
  FPE_lo <-  bin_e_FPE(b10,lower,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  if (TPE_lo > true_rate&FPE_lo<false_rate) return(lower)

  Power_root <- function(N){
    N =round(N)
    x = bin_e_BF_bound_10(threshold,N,alpha,beta,location,scale,prior_analysis,alternative,ROPE)

    if(de_an_prior == 1){
      pro = bin_e_TPE(x,N,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    } else{
      pro = bin_e_TPE(x,N,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE)
    }
    return(pro-true_rate)
  }

  N.power = round(stats::uniroot(Power_root,lower = lower,upper = upper)$root)+1

  N.extended = seq(N.power,N.power+20,2)
  Power.extended = unlist(lapply(N.extended, Power_root))

  if (any(Power.extended<0)){
    lower = which(Power.extended < 0)[1]
    N.power = round(stats::uniroot(Power_root,lower = lower,upper = upper)$root)+1

  }




  while(TRUE) {
    b10 <- bin_e_BF_bound_10(threshold,N.power,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    pro <- if (de_an_prior == 0) {
      bin_e_TPE(b10, N.power, h0,alpha_d, beta_d, location_d, scale_d, prior_design, alternative,ROPE)
    } else {
      bin_e_TPE(b10,N.power,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    }

    if (pro > true_rate) break
    N.power <- N.power + 1
  }
  b10 = bin_e_BF_bound_10(threshold,N.power,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  FPE =  bin_e_FPE(b10,N.power,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  if (FPE <= false_rate) return(N.power)

  alpha.root <- function(n) {
    n=round(n)
    b10 <- bin_e_BF_bound_10(threshold,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    bin_e_FPE(b10,n,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)-false_rate
  }
  N.alpha = round(stats::uniroot(alpha.root,lower = N.power,upper = upper)$root)
  return(N.alpha)

}



bin_e_N_01_finder <-function(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                             alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,false_rate,ROPE){
  lower = 10
  upper = 10000

  b10 =  bin_e_BF_bound_01(threshold,lower,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  TNE_lo =  bin_e_TPE(b10,lower,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  FNE_lo <-  if (de_an_prior == 1)
    bin_e_FNE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE) else
      bin_e_FNE(b10,lower,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE)

  if (TNE_lo > true_rate && FNE_lo < false_rate) {
    return(lower)
  } else if (TNE_lo > true_rate) {
    FN_root <- function(N){
      N =round(N)
      b10 =  bin_e_BF_bound_01(threshold,N,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
      pro <- if (de_an_prior == 1)
        bin_e_FNE(b10,N,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE) else
          bin_e_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE)

      pro-false_rate
    }
    return(round(stats::uniroot(FN_root, lower = lower, upper = upper)$root))
  }

  TN_root <- function(N){
    N =round(N)
    x = bin_e_BF_bound_01(threshold,N,alpha,beta,location,scale,prior_analysis,alternative,ROPE)

    pro = bin_e_TNE(x,N,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    return(pro-true_rate)
  }

  N.TN = round(stats::uniroot(TN_root,lower = lower,upper = upper)$root)+1

  N.extended = seq(N.TN,N.TN+20,2)
  Power.extended = unlist(lapply(N.extended, TN_root))

  if (any(Power.extended<0)){
    lower = which(Power.extended < 0)[1]
    N.TN = round(stats::uniroot(TN_root,lower = lower,upper = upper)$root)+1

  }




  while(TRUE) {
    b10 <- bin_e_BF_bound_01(threshold,N.TN,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    pro <- bin_e_TNE(b10,N.TN,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)

    if (pro > true_rate) break
    N.TN <- N.TN + 1
  }
  b10 = bin_e_BF_bound_01(threshold,N.TN,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  FNE =  if (de_an_prior == 1)
    bin_e_FNE(b10,N.TN,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE) else
      bin_e_FNE(b10,N.TN,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE)
  if (FNE <= false_rate) return(N.TN)

  FN_root <- function(N){
    N =round(N)
    b10 =  bin_e_BF_bound_01(threshold,N,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    pro <- if (de_an_prior == 1)
      bin_e_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE) else
        bin_e_FNE(b10,N,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE)

    pro- false_rate
  }
  N.FN = round(stats::uniroot(FN_root,lower = N.TN,upper = upper)$root)
  return(N.FN)

}

bin_e_table<-function(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                      alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,N, mode_bf,false_rate,ROPE,type_rate){
  if (mode_bf == "0") n = N else n = switch(
    type_rate,
    "positive" = bin_e_N_finder(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                                alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,false_rate,ROPE),
    "negative" = bin_e_N_01_finder(threshold,true_rate,h0,alpha,beta,location,scale,prior_analysis,alternative,
                                   alpha_d,beta_d,location_d,scale_d,prior_design,de_an_prior,false_rate,ROPE))

  # b bounds:
  b10 <- bin_e_BF_bound_10(threshold,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  b01 <-  bin_e_BF_bound_01(threshold,n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)

  max_BF <- 1 /bin_e_BF(round(location*n),n,alpha,beta,location,scale,prior_analysis,alternative,ROPE)

  # FPE and TPE:
  FPE       <- bin_e_FPE(b10,n,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  if (de_an_prior == 1) {
    TPE          <- bin_e_TPE(b10,n,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
    TPR_alpha    <- alpha
    TPR_beta     <- beta
    TPR_location <- location
    TPR_scale    <- scale
    TPR_prior    <- prior_analysis

  } else {
    TPE          <- bin_e_TPE(b10,n,h0,alpha_d,beta_d,location_d,scale_d,prior_design,alternative,ROPE)
    TPR_alpha    <- alpha_d
    TPR_beta     <- beta_d
    TPR_location <- location_d
    TPR_scale    <- scale_d
    TPR_prior    <- prior_design
  }
  # FNE and TNE:
  if (any(alternative == "two.sided" & max_BF < threshold | b01 == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- bin_e_FNE(b01,n,h0,TPR_alpha,TPR_beta,TPR_location,TPR_scale,TPR_prior,alternative,ROPE)
    TNE <- bin_e_TNE(b01,n,h0,alpha,beta,location,scale,prior_analysis,alternative,ROPE)
  }
  # table:
  tab.names <- c(
    "TruePositve",
    "FalseNegative",
    "TrueNegative",
    "FalsePositive",
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}

bin_e_bf10 <- function(threshold, n, alpha, beta, location, scale, prior_analysis, alternative, ROPE) {

  # Sequence of successes
  x <- seq(0, n, by = 3)

  # Compute BF10 and bounds
  BF10 <- bin_e_BF(x, n, alpha, beta, location, scale, prior_analysis, alternative, ROPE)
  b.BF10 <- bin_e_BF_bound_10(threshold, n, alpha, beta, location, scale, prior_analysis, alternative, ROPE)
  BF10_at_b <- bin_e_BF(b.BF10, n, alpha, beta, location, scale, prior_analysis, alternative, ROPE)

  # Compute BF01 and bounds
  BF01 <- 1 / BF10
  b.BF01 <- bin_e_BF_bound_01(threshold, n, alpha, beta, location, scale, prior_analysis, alternative, ROPE)
  BF01_at_b <- 1 / bin_e_BF(b.BF01, n, alpha, beta, location, scale, prior_analysis, alternative, ROPE)

  # Check if BF01 = D is impossible
  max.BF01 <- 1 / bin_e_BF(round(n / 2), n, alpha, beta, location, scale, prior_analysis, alternative, ROPE)
  impossible <- (alternative == "two.sided") && (max.BF01 < threshold || identical(b.BF01, "bound cannot be found"))

  # Titles for BF10
  main.bf10 <- if (length(b.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(round(BF10_at_b, 2)) ~ " when x = " ~ .(round(b.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(round(BF10_at_b[1], 2)) ~ "/" ~ .(round(BF10_at_b[2], 2)) ~
                  " when x = " ~ .(round(b.BF10[1], 2)) ~ " or " ~ .(round(b.BF10[2], 2))))
  }

  # Titles for BF01
  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(threshold)))
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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

Power_e_bin <- function(threshold, h0, alpha, beta, location, scale, prior_analysis, alternative,
                        alpha_d, beta_d, location_d, scale_d, prior_design, de_an_prior, N, ROPE) {

  # Sample size range
  smin <- 10
  smax <- N * 1.2
  sN <- ceiling(seq(smin, smax, length.out = 51))  # 51 points for smooth curves

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(sN))

  for (i in seq_along(sN)) {

    x10 <- bin_e_BF_bound_10(threshold, sN[i], alpha, beta, location, scale, prior_analysis, alternative, ROPE)
    x01 <- bin_e_BF_bound_01(threshold, sN[i], alpha, beta, location, scale, prior_analysis, alternative, ROPE)

    # True Positive
    TPE[i] <- if (de_an_prior == 1) {
      bin_e_TPE(x10, sN[i], h0, alpha, beta, location, scale, prior_analysis, alternative, ROPE)
    } else {
      bin_e_TPE(x10, sN[i], h0, alpha_d, beta_d, location_d, scale_d, prior_design, alternative, ROPE)
    }

    # False Negative
    FNE[i] <- if (de_an_prior == 1) {
      bin_e_FNE(x01, sN[i], h0, alpha, beta, location, scale, prior_analysis, alternative, ROPE)
    } else {
      bin_e_FNE(x01, sN[i], h0, alpha_d, beta_d, location_d, scale_d, prior_design, alternative, ROPE)
    }

    # False Positive & True Negative
    FPE[i] <- bin_e_FPE(x10, sN[i], h0, alpha, beta, location, scale, prior_analysis, alternative, ROPE)
    TNE[i] <- bin_e_TNE(x01, sN[i], h0, alpha, beta, location, scale, prior_analysis, alternative, ROPE)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # BF01 plot
  p2 <- ggplot2::ggplot(df_bf01, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

compute.prior.density.be.h1 <- function(h0,prop,alpha,beta,location,scale,prior_analysis,alternative,ROPE) {
  if (prior_analysis == "Point") return(rep(NA, length(prop)))
  bound_h1  <- switch(alternative,
                      "greater" = c(a = h0+ROPE, b = 1),
                      "less" = c(a = 0, b = h0+ROPE),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )

  prior_h1<- bin_prior(prop,alpha,beta,location,scale,prior_analysis)
  switch(alternative,
         "two.sided" = { prior_h1[prop>min(bound_h1)&prop<max(bound_h1)]=0 },
         "greater" = { prior_h1[prop<bound_h1[1]]=0 },
         "less" = { prior_h1[prop>bound_h1[2]]=0 }
  )
  prior_h1
}


compute.prior.density.be.h0 <- function(h0,prop,alpha,beta,location,scale,prior_analysis,alternative,ROPE) {
  if (prior_analysis == "Point") return(rep(NA, length(prop)))
  bound_h0  <- switch(alternative,
                      "greater" = c(a = h0, b = h0+ROPE),
                      "less" = c(a = h0+ROPE, b = h0),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )

  prior_h0<- bin_prior(prop,alpha,beta,location,scale,prior_analysis)
  switch(alternative,
         "two.sided" = { prior_h0[prop<min(bound_h0)|prop>max(bound_h0)]=0 },
         "greater" = { prior_h0[prop>bound_h0[2]]=0 },
         "less" = { prior_h0[prop<bound_h0[1]]=0 }
  )
  prior_h0
}





bin_e_prior_plot <- function(h0,
                             alpha, beta, location, scale, prior_analysis,
                             alpha_d, beta_d, location_d, scale_d, prior_design,
                             alternative, de_an_prior, ROPE) {

  # ---- Plot bounds ----
  plot.bounds <- switch(alternative,
                        "greater" = c(h0, 1),
                        "less" = c(0, h0),
                        "two.sided" = c(0, 1))

  theta <- seq(plot.bounds[1], plot.bounds[2], 0.002)

  # ---- Compute H1 and H0 priors ----
  prior_h1 <- compute.prior.density.be.h1(
    h0, theta, alpha, beta, location, scale,
    prior_analysis, alternative, ROPE
  )

  prior_h0 <- compute.prior.density.be.h0(
    h0, theta, alpha, beta, location, scale,
    prior_analysis, alternative, ROPE
  )

  # ---- Long format data (H1/H0) ----
  df_lines <- data.frame(
    theta = rep(theta, 2),
    Density = c(prior_h1, prior_h0),
    Prior = rep(c("H1 - Analysis Prior",
                  "H0 - Analysis Prior"),
                each = length(theta))
  )

  # ---- Legend position ----
  legend_pos <- switch(alternative,
                       "greater"   = c(0.75, 0.95),
                       "two.sided" = c(0.75, 0.95),
                       "less"      = c(0.2, 0.95))

  # ---- Base ggplot ----
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = df_lines,
      ggplot2::aes(x = theta,
                   y = Density,
                   color = Prior,
                   linetype = Prior,
                   linewidth = Prior)
    ) +
    ggplot2::scale_color_manual(values = c(
      "H1 - Analysis Prior" = "black",
      "H0 - Analysis Prior" = "black",
      "H1 - Design Prior"   = "gray"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "H1 - Analysis Prior" = "solid",
      "H0 - Analysis Prior" = "dashed",
      "H1 - Design Prior"   = "solid"
    )) +
    ggplot2::scale_linewidth_manual(values = c(
      "H1 - Analysis Prior" = 1.2,
      "H0 - Analysis Prior" = 1.2,
      "H1 - Design Prior"   = 2
    )) +
    ggplot2::labs(
      x = expression(bold(theta)),
      y = "Density",
      title = bquote(bold("Prior distribution on "~theta~
                            " under the alternative"))
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  # ---- Add design prior line (non-point) ----
  if (de_an_prior == 0 && prior_design != "Point") {

    prior_design_vals <- compute.prior.density.be.h1(
      h0, theta,
      alpha_d, beta_d,
      location_d, scale_d,
      prior_design,
      alternative, ROPE
    )

    df_design <- data.frame(
      theta = theta,
      Density = prior_design_vals,
      Prior = "H1 - Design Prior"
    )

    df_design <- df_design[!is.na(df_design$Density), ]

    p <- p +
      ggplot2::geom_line(
        data = df_design,
        ggplot2::aes(x = theta,
                     y = Density,
                     color = Prior,
                     linetype = Prior,
                     linewidth = Prior)
      )
  }

  # ---- Add vertical arrow for Point design prior ----
  if (de_an_prior == 0 && prior_design == "Point") {

    ylim_max <- max(prior_h1, prior_h0, na.rm = TRUE)

    # Invisible dummy line for legend
    df_dummy <- data.frame(
      theta = c(NA, NA),
      Density = c(NA, NA),
      Prior = "H1 - Design Prior"
    )

    p <- p +
      ggplot2::geom_line(
        data = df_dummy,
        ggplot2::aes(x = theta,
                     y = Density,
                     color = Prior,
                     linetype = Prior,
                     linewidth = Prior),
        na.rm = TRUE,
        show.legend = TRUE
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = location_d,
                     xend = location_d,
                     y = 0,
                     yend = ylim_max),
        color = "gray",
        linetype = "dashed",
        arrow = ggplot2::arrow(
          length = grid::unit(0.1, "inches")
        )
      )
  }

  return(p)
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

r_prior<- function(rho,k,location,scale,dff,prior_analysis, alpha, beta,a,b){

  switch(prior_analysis,
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


r_BF10<-function(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis){
  x = NA
  bound  <- switch(alternative,
                   "greater" = c(a = h0, b = 1),
                   "less" = c(a = -1, b = h0),
                   "two.sided" = c(a = -1, b = 1)
  )
  normalization <- if (alternative == "two.sided") {
    switch(prior_analysis,
           "d_beta"   = 1,
           "beta" = 1,
           "Moment"   = { pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})

  }else{
    switch(prior_analysis,
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,-1,1)-p_beta(bound[1], 1/k,1/k,-1,1) ,
           "beta" = p_beta(bound[2], alpha, beta,-1,1)-p_beta(bound[1], alpha, beta,-1,1),
           "Moment"   = {pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})
  }

  # Define the integrand function for marginal likelihood under H1
  int <- function(rho, ri) {
    d_cor(ri, rho, n) * r_prior(rho, k, location, scale, dff, prior_analysis, alpha, beta, min(bound), max(bound))
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

r_BF_bound_10 <-function(threshold,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis){
  y <- numeric(0)
  Bound_finding <-function(r)r_BF10(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)- threshold

  x <- tryCatch(stats::uniroot(Bound_finding, lower = -.99, upper = h0,tol = 1e-5)$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower =  h0, upper = .99,tol = 1e-5)$root, error = function(e) NA)
  results <- c(x, y)
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")

  BF.vals  <- r_BF10(results,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)
  BF.close <- which(round(BF.vals, 2) == round(threshold, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

r_BF_bound_01 <-function(threshold,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis){
  r_BF_bound_10(1/threshold,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)
}

p_cor<-function(limit,rho,n,lower.tail){

  stats::pnorm(r_mean(limit),r_mean(rho),sd = r_sd(n),lower.tail =  lower.tail)


}

r_TPE <-function(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)

  if (prior_analysis =="Point"){
    x = switch(alternative,
               "two.sided" = {p_cor(max(r),location,n,lower.tail = F)+ p_cor(min(r),location,n,lower.tail = T)},
               "greater"  = {p_cor(r,location,n,lower.tail =F)},
               "less"  = {p_cor(r,location,n,lower.tail =T)}
    )
    return(x)
  }

  bound  <- switch(alternative,
                   "greater" = c(a = h0, b = 1),
                   "less" = c(a = -1, b = h0),
                   "two.sided" = c(a = -1, b = 1)
  )
  normalization <-   normalization <- if (alternative == "two.sided") {
    switch(prior_analysis,
           "d_beta"   = 1,
           "beta" = 1,
           "Moment"   = { pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})

  }else{
    switch(prior_analysis,
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,-1,1)-p_beta(bound[1], 1/k,1/k,-1,1) ,
           "beta" = p_beta(bound[2], alpha, beta,-1,1)-p_beta(bound[1], alpha, beta,-1,1),
           "Moment"   = {pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})
  }
  int <- function(rho) {
    prob <- switch(alternative,
                   "two.sided" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                     p_cor(min(r), rho, n, lower.tail = TRUE),
                   "greater"  = p_cor(r, rho, n, lower.tail = FALSE),
                   "less"  = p_cor(r, rho, n, lower.tail = TRUE)
    )

    prob * r_prior(rho, k, location, scale, dff, prior_analysis, alpha, beta,min(bound),max(bound)) / normalization
  }
  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-4)$value
  return(x)

}

r_FNE <-function(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)


  if (prior_analysis =="Point"){
    x = switch(alternative,
               "two.sided" = {p_cor(max(r),location,n,lower.tail = T)- p_cor(min(r),location,n,lower.tail = T)},
               "greater"  = {p_cor(r,location,n,lower.tail =T)},
               "less"  = {p_cor(r,location,n,lower.tail =F)}
    )
    return(x)
  }


  bound  <- switch(alternative,
                   "greater" = c(a = h0, b = 1),
                   "less" = c(a = -1, b = h0),
                   "two.sided" = c(a = -1, b = 1)
  )

  normalization <-  normalization <- if (alternative == "two.sided") {
    switch(prior_analysis,
           "d_beta"   = 1,
           "beta" = 1,
           "Moment"   = { pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})

  }else{
    switch(prior_analysis,
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,-1,1)-p_beta(bound[1], 1/k,1/k,-1,1) ,
           "beta" = p_beta(bound[2], alpha, beta,-1,1)-p_beta(bound[1], alpha, beta,-1,1),
           "Moment"   = {pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2)})
  }
  int <- function(rho) {
    prob <- switch(alternative,
                   "two.sided" = p_cor(max(r), rho, n, lower.tail = TRUE) -
                     p_cor(min(r), rho, n, lower.tail = TRUE),
                   "greater"  = p_cor(r, rho, n, lower.tail = TRUE),
                   "less"  = p_cor(r, rho, n, lower.tail = FALSE)
    )

    prob * r_prior(rho, k, location, scale, dff, prior_analysis, alpha, beta,min(bound),max(bound)) / normalization
  }


  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = 1e-8, subdivisions=10000000)$value
  return(x)

}

r_FPE <-function(r,n,h0,alternative){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)

  x <- switch(alternative,
              "two.sided" = p_cor(max(r), h0, n, lower.tail = FALSE) +
                p_cor(min(r), h0, n, lower.tail = TRUE),
              "greater"  = p_cor(r, h0, n, lower.tail = FALSE),
              "less"  = p_cor(r, h0, n, lower.tail = TRUE)
  )
  return(x)

}


r_TNE <-function(r,n,h0,alternative){

  if (any(r == "bound cannot be found") || length(r) == 0) return(0)

  bound  <- switch(alternative,
                   "greater" = c(a = h0, b = 1),
                   "less" = c(a = -1, b = h0),
                   "two.sided" = c(a = -1, b = 1)
  )

  x <- switch(alternative,
              "two.sided" = p_cor(max(r), h0, n, lower.tail = TRUE) -
                p_cor(min(r), h0, n, lower.tail = TRUE),
              "greater"  = p_cor(r, h0, n, lower.tail = TRUE),
              "less"  = p_cor(r, h0, n, lower.tail = FALSE)
  )

  return(x)

}



r_N_finder<-function(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                     location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,false_rate){

  lo = 10
  upper = 5000

  r = r_BF_bound_10(threshold,lo,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)
  TPE_lo <- if (de_an_prior == 1)
    r_TPE(r,lo,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis) else
      r_TPE(r,lo,k_d, alpha_d, beta_d,h0,alternative,location_d,scale_d,dff_d,prior_design)
  FPE_lo <-  r_FPE(r,lo,h0,alternative )

  if (TPE_lo > true_rate&FPE_lo<false_rate) return(lo)

  Power_root <- function(N) {
    r <- r_BF_bound_10(threshold, N, k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis)
    pro <- if (de_an_prior==0){ r_TPE(r, N, k_d, alpha_d, beta_d, h0, alternative, location_d, scale_d, dff_d, prior_design) }else r_TPE(r, N, k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis)

    pro - true_rate
  }

  N.power = stats::uniroot(Power_root,lower = lo,upper = upper)$root

  ## checking if the N lead to an acceptable alpha level
  r = r_BF_bound_10(threshold,N.power,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)

  FPE = r_FPE(r,N.power,h0,alternative)
  if (FPE <= false_rate) return(N.power)

  alpha.root <- function(n) {
    r <- r_BF_bound_10(threshold,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)
    r_FPE(r,n,h0,alternative)-false_rate
  }
  N.alpha = stats::uniroot(alpha.root,lower = N.power,upper = upper)$root
  return(N.alpha)
}
r_N_01_finder<-function(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                        location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,false_rate){

  lo = 10
  upper = 5000

  r = r_BF_bound_01(threshold,lo,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)
  TNE_lo <- r_TNE(r,lo,h0,alternative )
  FNE_lo <-  if (de_an_prior == 1)
    r_FNE(r,lo,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis) else
      r_FNE(r,lo,k_d, alpha_d, beta_d,h0,alternative,location_d,scale_d,dff_d,prior_design)

  if (TNE_lo > true_rate && TNE_lo < false_rate) {
    return(lo)
  } else if (TNE_lo > true_rate) {
    FN.root <- function(n) {
      r <- r_BF_bound_01(threshold,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)
      FNE<-  if (de_an_prior == 1)
        r_FNE(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis) else
          r_FNE(r,n,k_d, alpha_d, beta_d,h0,alternative,location_d,scale_d,dff_d,prior_design)
      FNE- false_rate
    }
    return(stats::uniroot(FN.root, lower = lo, upper = upper)$root)
  }

  TN_root <- function(N) {
    r <- r_BF_bound_01(threshold, N, k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis)
    pro <- r_TNE(r,N,h0,alternative )
    pro - true_rate
  }

  N.TN <- tryCatch(
    stats::uniroot(TN_root, lower = lo, upper = upper)$root,
    error = function(e) 20
  )

  ## checking if the N lead to an acceptable alpha level
  r = r_BF_bound_01(threshold,N.TN,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)

  FNE =   if (de_an_prior==0){ r_FNE(r, N.TN, k_d, alpha_d, beta_d, h0, alternative, location_d, scale_d, dff_d, prior_design) }else r_FNE(r, N.TN, k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis)

  if (FNE <= false_rate) return(N.TN)

  FN.root <- function(n) {
    r   <- r_BF_bound_01(threshold,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)
    FNE <- if (de_an_prior==0){
      r_FNE(r, n, k_d, alpha_d, beta_d, h0, alternative, location_d, scale_d, dff_d, prior_design)
    } else {
      r_FNE(r, n, k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis)
      }
    FNE-false_rate
  }
  N.FN = stats::uniroot(FN.root,lower = N.TN,upper = upper)$root
  return(N.FN)
}

r_table<-function(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                  location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,false_rate,type_rate ){

  n <- if (mode_bf == 1) {
    switch(type_rate,
           "positive" = ceiling(r_N_finder(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                                           location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,false_rate)),
           "negative" = ceiling(r_N_01_finder(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                                              location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,false_rate)))} else  n = N

         # r bounds:
         r10 <- r_BF_bound_10(threshold,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)
         r01 <-  r_BF_bound_01(threshold,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)

         # max BF10 possible:
          max_BF <- 1 / r_BF10(h0,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis)
          BF_D   <- r10

         # FPE and TPE:
          FPE       <- r_FPE(r10,n,h0,alternative)
           if (de_an_prior == 1) {
           TPE         <- r_TPE(r10, n, k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis)
           TPR_prior   <- prior_analysis
           TPR_k       <- k
           TPR_alpha   <- alpha
           TPR_beta    <- beta
           TPR_location<- location
           TPR_scale   <- scale
           TPR_dff     <- dff
           } else {
           TPE       <- r_TPE(r10, n, k_d, alpha_d, beta_d, h0, alternative, location_d, scale_d, dff_d, prior_design)
           TPR_prior   <- prior_design
           TPR_k       <- k_d
           TPR_alpha   <- alpha_d
           TPR_beta    <- beta_d
           TPR_location<- location_d
           TPR_scale   <- scale_d
           TPR_dff     <- dff_d
           }
           # FNE and TNE:
            if (any(alternative == "two.sided" & max_BF < threshold | BF_D == "bound cannot be found")) {
                 FNE <- 0
                 TNE <- 0
             } else {
             FNE <- r_FNE(r01,n,TPR_k, TPR_alpha, TPR_beta,h0,alternative,TPR_location,TPR_scale,TPR_dff,TPR_prior)
             TNE <- r_TNE(r01,n,h0,alternative)
             }


            # table:
             tab.names <- c("TruePositve",
                            "FalseNegative",
                            "TrueNegative",
                            "FalsePositive",
                            "Required N"
                            )
            table <- data.frame(TPE, FNE, TNE, FPE, n, check.names = FALSE, row.names = NULL)
            colnames(table) <- tab.names
            table
}


compute.prior.density.r <- function(rho, k,location,scale,dff,prior_analysis, alpha, beta,alternative) {
  if (prior_analysis == "Point") return(rep(NA, length(rho)))
  bound  <- switch(alternative,
                   "greater" = c(a = location, b = 1),
                   "less" = c(a = -1, b = location),
                   "two.sided" = c(a = -1, b = 1)
  )
  normalization <- if (alternative == "two.sided") 1 else
    switch(prior_analysis,
           "Normal" = stats::pnorm(bound[2],location,scale)-stats::pnorm(bound[1],location,scale),
           "d_beta"   = p_beta(bound[2], 1/k, 1/k,min(bound),max(bound))-p_beta(bound[1], 1/k,1/k,min(bound),max(bound)) ,
           "Moment"   = pmom(bound[2]-location, tau=scale^2)-pmom(bound[1]-location, tau=scale^2),
           "t_dis" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0),
           "beta" = p_beta(bound[2], alpha, beta,min(bound),max(bound))-p_beta(bound[1], alpha, beta,min(bound),max(bound)))


  r_prior(rho,k,location,scale,dff,prior_analysis, alpha, beta,min(bound),max(bound)) / normalization
}


r_prior_plot <- function(k, alpha, beta, h0,
                         location, scale, dff, prior_analysis, de_an_prior,
                         k_d, alpha_d, beta_d,
                         location_d, scale_d, dff_d,
                         prior_design, alternative) {

  # ---- Determine bounds ----
  bound <- switch(alternative,
                  "greater"   = c(h0, 1),
                  "less"      = c(-1, h0),
                  "two.sided" = c(-1, 1))

  rho <- seq(bound[1], bound[2], .01)

  # ---- Compute priors ----
  prior.analysis <- compute.prior.density.r(
    rho, k, location, scale, dff,
    prior_analysis, alpha, beta, alternative
  )

  # Base data frame
  df <- data.frame(
    rho = rho,
    Density = prior.analysis,
    Prior = "H1 - Analysis Prior"
  )

  # ---- Add design prior if needed ----
  if (de_an_prior == 0) {

    if (prior_design == "Point") {
      # Dummy row for legend only
      df_design <- data.frame(
        rho = c(NA, NA),
        Density = c(NA, NA),
        Prior = "H1 - Design Prior"
      )
      df <- rbind(df, df_design)

    } else {
      prior.design <- compute.prior.density.r(
        rho, k_d, location_d, scale_d, dff_d,
        prior_design, alpha_d, beta_d, alternative
      )

      df_design <- data.frame(
        rho = rho,
        Density = prior.design,
        Prior = "H1 - Design Prior"
      )

      df <- rbind(df, df_design)
    }
  }

  # ---- Y limits ----
  ylim_max <- max(df$Density[is.finite(df$Density)], na.rm = TRUE)

  # ---- Legend position (match t1_prior_plot logic) ----
  legend_pos <- switch(alternative,
                       "greater"   = c(0.65, 0.95),
                       "two.sided" = c(0.65, 0.95),
                       "less"      = c(0.05, 0.95))

  # ---- Build ggplot ----
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = rho,
                                    y = Density,
                                    color = Prior,
                                    linetype = Prior)) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(
      "H1 - Analysis Prior" = "black",
      "H1 - Design Prior"   = "gray"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "H1 - Analysis Prior" = "solid",
      "H1 - Design Prior"   = "dashed"
    )) +
    ggplot2::labs(
      x = expression(bold(rho)),
      y = "density",
      title = bquote(bold("Prior distribution on "~rho~" under the alternative"))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, ylim_max),
                             xlim = bound) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.justification = c(0, 1),
      legend.background =
        ggplot2::element_rect(fill = scales::alpha("white", 0.8), color = NA),
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 1.5))
    )

  # ---- Add vertical arrow for point prior ----
  if (de_an_prior == 0 && prior_design == "Point") {

    p <- p +
      ggplot2::annotate("segment",
                        x = location_d,
                        xend = location_d,
                        y = 0,
                        yend = ylim_max,
                        color = "gray",
                        linetype = "dashed",
                        arrow = ggplot2::arrow(
                          length = grid::unit(0.1, "inches")
                        ))
  }

  return(p)
}
r_bf10_p <- function(threshold, n, k, alpha, beta, h0, alternative,
                     location, scale, dff, prior_analysis) {

  rr <- seq(-0.99, 0.99, 0.01)

  BF10   <- r_BF10(rr, n, k, alpha, beta, h0, alternative,
                   location, scale, dff, prior_analysis)
  r.BF10 <- r_BF_bound_10(threshold, n, k, alpha, beta, h0, alternative,
                          location, scale, dff, prior_analysis)

  BF01   <- 1 / BF10
  r.BF01 <- r_BF_bound_01(threshold, n, k, alpha, beta, h0, alternative,
                          location, scale, dff, prior_analysis)

  max.BF01 <- 1 / r_BF10(h0, n, k, alpha, beta, h0, alternative,
                         location, scale, dff, prior_analysis)

  impossible <- (alternative == "two.sided") &&
    (max.BF01 < threshold || identical(r.BF01, "bound cannot be found"))

  ## ---------- Titles ----------
  main.bf10 <- if (length(r.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(threshold) ~ " when r = " ~ .(round(r.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(threshold) ~ " when r = " ~ .(round(r.BF10[1], 2)) ~
                  " or " ~ .(round(r.BF10[2], 2))))
  }

  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(threshold)))
  } else if (length(r.BF01) == 1) {
    bquote(bold("BF"[0][1] ~ "=" ~ .(threshold) ~ " when r = " ~ .(round(r.BF01, 2))))
  } else {
    bquote(bold("BF"[0][1] ~ "=" ~ .(threshold) ~ " when r = " ~ .(round(r.BF01[1], 2)) ~
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +   # ← all black now
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


Power_r <- function(threshold, k, alpha, beta, h0, alternative,
                    location, scale, dff, prior_analysis,
                    k_d, alpha_d, beta_d,
                    location_d, scale_d, dff_d, prior_design,
                    de_an_prior, N) {

  Ns <- seq(4, ceiling(N * 1.2), length.out = 31)

  TPE <- FPE <- TNE <- FNE <- numeric(length(Ns))

  for (i in seq_along(Ns)) {

    r10 <- r_BF_bound_10(threshold, Ns[i], k, alpha, beta, h0, alternative,
                         location, scale, dff, prior_analysis)
    r01 <- r_BF_bound_01(threshold, Ns[i], k, alpha, beta, h0, alternative,
                         location, scale, dff, prior_analysis)

    TPE[i] <- if (de_an_prior == 1)
      r_TPE(r10, Ns[i], k, alpha, beta, h0, alternative,
            location, scale, dff, prior_analysis)
    else
      r_TPE(r10, Ns[i], k_d, alpha_d, beta_d, h0, alternative,
            location_d, scale_d, dff_d, prior_design)

    FPE[i] <- r_FPE(r10, Ns[i], h0, alternative)

    TNE[i] <- r_TNE(r01, Ns[i], h0, alternative)

    FNE[i] <- if (de_an_prior == 1)
      r_FNE(r01, Ns[i], k, alpha, beta, h0, alternative,
            location, scale, dff, prior_analysis)
    else
      r_FNE(r01, Ns[i], k_d, alpha_d, beta_d, h0, alternative,
            location_d, scale_d, dff_d, prior_design)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  p2 <- ggplot2::ggplot(df2,
                        ggplot2::aes(SampleSize, Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
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

    # If no root is found, ROPExpand the search range and try again
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

    # If no root is found, ROPExpand the search range and try again
    lower <- lower + step
  }

}


re_BF10i<-function(r,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE){
  x = NA

  bound_h1  <- switch(alternative,
                      "greater" = c(a = h0+ROPE, b = 1),
                      "less" = c(a = -1, b = h0+ROPE),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )
  bound_h0  <- switch(alternative,
                      "greater" = c(a = h0, b = h0+ROPE),
                      "less" = c(a = h0+ROPE, b = h0),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )

  normalizationh1 <- switch(alternative,
                            "two.sided" = switch(prior_analysis,
                                                 "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                                 "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                                 "Moment"          = {
                                                   (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                                     (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                                 }),

                            "less"  = switch(prior_analysis,
                                             "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                             "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                             "Moment"          = {
                                               (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                             }),
                            "greater"  = switch(prior_analysis,
                                                "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                                "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                                "Moment"          = {
                                                  (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                                })
  )
  normalizationh0 <- switch(prior_analysis,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )

  int  <- function(rho){d_cor(r,rho,n)*r_prior(rho,k,location,scale,dff,prior_analysis, alpha, beta,-1,1)/normalizationh1
  }

  if (alternative == "two.sided"){
    lh1 = stats::integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-5,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-5,stop.on.error = F)$value
  }else{
    lh1 = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-5,stop.on.error = F)$value

  }
  lh0 = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value

  x = (lh1/normalizationh1)/(lh0/normalizationh0)
  return(x)
}
re_BF10<-function(r,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE){
  sapply(r, re_BF10i, n = n, k = k, alpha = alpha, beta = beta,
         h0 = h0, alternative = alternative, location = location,
         scale = scale, dff = dff, prior_analysis = prior_analysis, ROPE = ROPE)
}



re_BF_bound_10 <-function(threshold,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE){
  y <- numeric(0)
  Bound_finding <-function(r){
    re_BF10(r,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)- threshold
  }
  opt_result <- stats::optimize(Bound_finding, interval = c(-.999, .999))$minimum

  if (alternative=="two.sided"){
    x <- r_auto_uniroot_fixed_upper (Bound_finding,opt_result, lower = -1, step = .05, max_attempts = 25)
    y <- r_auto_uniroot_fixed_lower(Bound_finding,opt_result, upper = 1, step = .05, max_attempts = 25)
  }
  if (alternative == "greater"){

    x <- r_auto_uniroot_fixed_lower(Bound_finding,h0, upper = 1, step = .05, max_attempts = 25)
  }

  if (alternative == "less"){

    x <- r_auto_uniroot_fixed_upper (Bound_finding,h0, lower = -1, step = .05, max_attempts = 25)
  }
  results <- c(x, y)
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")

  BF.vals  <- re_BF10(results,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  BF.close <- which(round(BF.vals, 2) == round(threshold, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

re_BF_bound_01 <-function(threshold,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE){
  re_BF_bound_10 (1/threshold,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
}

re_TPE <-function(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE){

  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }

  if (prior_analysis =="Point"){
    x = switch(alternative,
               "two.sided" = {p_cor(max(r),location,n,lower.tail = F)+ p_cor(min(r),location,n,lower.tail = T)},
               "greater"  = {p_cor(r,location,n,lower.tail =F)},
               "less"  = {p_cor(r,location,n,lower.tail =T)}
    )
    return(x)
  }

  bound_h1  <- switch(alternative,
                      "greater" = c(a = h0+ROPE, b = 1),
                      "less" = c(a = -1, b = h0+ROPE),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )
  normalizationh1 <- switch(alternative,
                            "two.sided" = switch(prior_analysis,
                                                 "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                                 "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                                 "Moment"          = {
                                                   (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                                     (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                                 }),

                            "less"  = switch(prior_analysis,
                                             "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                             "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                             "Moment"          = {
                                               (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                             }),
                            "greater"  = switch(prior_analysis,
                                                "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                                "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                                "Moment"          = {
                                                  (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                                })
  )

  int <- function(rho) {
    pro <- switch(alternative,
                  "two.sided" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  "greater"  = p_cor(r, rho, n, lower.tail = FALSE),
                  "less"  = p_cor(r, rho, n, lower.tail = TRUE)
    )

    pro * r_prior(rho, k, location, scale, dff, prior_analysis, alpha, beta,-1,1) / normalizationh1
  }



  if (alternative == "two.sided"){
    x = stats::integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-5,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-5,stop.on.error = F)$value
  }else{
    x = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-5,stop.on.error = F)$value

  }
  return(x)

}

re_FNE <-function(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE){

  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }
  if (prior_analysis =="Point"){
    x = switch(alternative,
               "two.sided" = {p_cor(max(r),location,n,lower.tail = T)- p_cor(min(r),location,n,lower.tail = T)},
               "greater"  = {p_cor(r,location,n,lower.tail =T)},
               "less"  = {p_cor(r,location,n,lower.tail =F)}
    )
    return(x)
  }
  bound_h1  <- switch(alternative,
                      "greater" = c(a = h0+ROPE, b = 1),
                      "less" = c(a = -1, b = h0+ROPE),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )

  normalizationh1 <- switch(alternative,
                            "two.sided" = switch(prior_analysis,
                                                 "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                                 "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                                 "Moment"          = {
                                                   (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                                     (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                                 }),

                            "less"  = switch(prior_analysis,
                                             "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                             "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                             "Moment"          = {
                                               (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                             }),
                            "greater"  = switch(prior_analysis,
                                                "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                                "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                                "Moment"          = {
                                                  (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                                })
  )

  int <- function(rho) {
    pro <- switch(alternative,
                  "two.sided" = p_cor(max(r), rho, n, lower.tail = T) -
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  "greater"  = p_cor(r, rho, n, lower.tail = T),
                  "less"  = p_cor(r, rho, n, lower.tail = F)
    )

    pro * r_prior(rho, k, location, scale, dff, prior_analysis, alpha, beta,-1,1) / normalizationh1
  }


  if (alternative == "two.sided"){
    x = stats::integrate(int,lower = -1,upper = bound_h1[1], rel.tol=1e-10,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = 1, rel.tol=1e-10,stop.on.error = F)$value
  }else{
    x = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=1e-10,stop.on.error = F)$value

  }

  return(x)

}

re_FPE <-function(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE){

  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }
  bound_h0  <- switch(alternative,
                      "greater" = c(a = h0, b = h0+ROPE),
                      "less" = c(a = h0+ROPE, b = h0),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )
  normalizationh0 <- switch(prior_analysis,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )


  int <- function(rho) {
    pro <- switch(alternative,
                  "two.sided" = p_cor(max(r), rho, n, lower.tail = FALSE) +
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  "greater"  = p_cor(r, rho, n, lower.tail = FALSE),
                  "less"  = p_cor(r, rho, n, lower.tail = TRUE)
    )

    pro * r_prior(rho, k, location, scale, dff, prior_analysis, alpha, beta,-1,1) / normalizationh0
  }

  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value


  return(x)

}

re_TNE <-function(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE){

  if (any(r =="bound cannot be found" | length(r)==0)){
    r=0
    return(r)
  }
  bound_h0  <- switch(alternative,
                      "greater" = c(a = h0, b = h0+ROPE),
                      "less" = c(a = h0+ROPE, b = h0),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )

  normalizationh0 <- switch(prior_analysis,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )

  int <- function(rho) {
    pro <- switch(alternative,
                  "two.sided" = p_cor(max(r), rho, n, lower.tail = TRUE) -
                    p_cor(min(r), rho, n, lower.tail = TRUE),
                  "greater"  = p_cor(r, rho, n, lower.tail = TRUE),
                  "less"  = p_cor(r, rho, n, lower.tail = FALSE)
    )

    pro * r_prior(rho, k, location, scale, dff, prior_analysis, alpha, beta,-1,1) / normalizationh0
  }

  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=1e-5,stop.on.error = F)$value


  return(x)

}




re_N_finder<-function(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                      location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,false_rate,ROPE){
  lo = 10
  upper = 5000

  r = re_BF_bound_10(threshold,lo,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  TPE_lo <- if (de_an_prior == 1)
    re_TPE(r,lo,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE) else
      re_TPE(r,lo,k_d, alpha_d, beta_d,h0,alternative,location_d,scale_d,dff_d,prior_design,ROPE)
  FPE_lo <-  re_FPE(r,lo,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)

  if (TPE_lo > true_rate && FPE_lo < false_rate) {
    return(lo)
  } else if (TPE_lo > true_rate) {
    alpha.root <- function(n) {
      r <- re_BF_bound_10(threshold, n, k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis, ROPE)
      re_FPE(r, n, k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis, ROPE) - false_rate
    }
    return(stats::uniroot(alpha.root, lower = lo, upper = upper)$root)
  }



  Power_root <- function(N){

    r = re_BF_bound_10(threshold,N,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)

    if (de_an_prior == 0 ){
      pro = re_TPE(r,N,k_d, alpha_d, beta_d,h0,alternative,location_d,scale_d,dff_d,prior_design,ROPE)
    }else {
      pro = re_TPE(r,N,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)

    }
    return(pro-true_rate)
  }

  N.power <- robust_uniroot(Power_root, lower = lo)
  r = re_BF_bound_10(threshold, N.power,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  FPE = re_FPE(r, N.power,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  if (FPE <= false_rate) return(N.power)

  alpha.root <- function(n) {
    r <- re_BF_bound_10(threshold, n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
    re_FPE(r, n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)-false_rate
  }
  N.alpha = stats::uniroot(alpha.root,lower = N.power,upper = upper)$root
  return(N.alpha)
}

re_N_01_finder<-function(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                         location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,false_rate,ROPE){
  lo = 10
  upper = 5000

  r = re_BF_bound_01(threshold,lo,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  TNE_lo <- re_TNE(r,lo,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  FNE_lo <-  if (de_an_prior == 1)
    re_FNE(r,lo,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE) else
      re_FNE(r,lo,k_d, alpha_d, beta_d,h0,alternative,location_d,scale_d,dff_d,prior_design,ROPE)


  if (TNE_lo > true_rate && FNE_lo < false_rate) {
    return(lo)
  } else if (TNE_lo > true_rate) {
    FN.root <- function(n) {
      r  <- re_BF_bound_01(threshold, n, k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis, ROPE)
      FNE<-  if (de_an_prior == 0 ){
        pro = re_FNE(r,n,k_d, alpha_d, beta_d,h0,alternative,location_d,scale_d,dff_d,prior_design,ROPE)
      }else {
        pro = re_FNE(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)

      }
      FNE - false_rate
    }
    return(stats::uniroot(FN.root, lower = lo, upper = upper)$root)
  }


  TN_root <- function(N){

    r = re_BF_bound_01(threshold,N,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
    pro = re_TNE(r,N,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
    return(pro-true_rate)
  }

  N.TN <- robust_uniroot(TN_root, lower = lo)
  r    <- re_BF_bound_01(threshold, N.TN,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  FNE  <- re_FNE(r, N.TN,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  if (FNE <= false_rate) return(N.TN)

  FN.root <- function(n) {
    r <- re_BF_bound_01(threshold, n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
    FNE<-  if (de_an_prior == 0 ){
      pro = re_FNE(r,n,k_d, alpha_d, beta_d,h0,alternative,location_d,scale_d,dff_d,prior_design,ROPE)
    }else {
      pro = re_FNE(r,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)

    }
    FNE - false_rate
  }
  N.FN = stats::uniroot(FN.root,lower = N.TN,upper = upper)$root
  return(N.FN)
}

re_table<-function(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                   location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior,N, mode_bf,false_rate,ROPE,type_rate){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)

  n <- if (mode_bf == 1) {
    switch(type_rate,
           "positive" = ceiling(re_N_finder(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                                            location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,false_rate,ROPE)),
           "negative" = ceiling(re_N_01_finder(threshold,true_rate,prior_analysis,k, alpha, beta,h0,location,scale,dff, alternative ,prior_design,
                                               location_d,k_d, alpha_d, beta_d,scale_d,dff_d,de_an_prior ,false_rate,ROPE)))} else  n = N

  # r bounds:
  r10 = re_BF_bound_10(threshold,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  r01 = re_BF_bound_01(threshold,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)

  # max BF10 possible:
  max_BF <- 1 / re_BF10(h0,n,k,alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  BF_D   <- r10

  # FPE and TPE:
  FPE       <- re_FPE(r10,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)

  if (de_an_prior == 1) {
    TPE           <- re_TPE(r10,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
    TPR_k         <- k
    TPR_alpha     <- alpha
    TPR_beta      <- beta
    TPR_location  <- location
    TPR_scale     <- scale
    TPR_dff       <- dff
    TPR_prior     <- prior_analysis

  } else {
    TPE           <- re_TPE(r10,n,k_d, alpha_d, beta_d,h0,alternative,location_d,scale_d,dff_d,prior_design,ROPE)
    TPR_k         <- k_d
    TPR_alpha     <- alpha_d
    TPR_beta      <- beta_d
    TPR_location  <- location_d
    TPR_scale     <- scale_d
    TPR_dff       <- dff_d
    TPR_prior     <- prior_design
  }
  # FNE and TNE:
  if (any(alternative == "two.sided" & max_BF < threshold | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- re_FNE(r01,n,TPR_k, TPR_alpha, TPR_beta,h0,alternative,TPR_location,TPR_scale,TPR_dff,TPR_prior,ROPE)
    TNE <- re_TNE(r01,n,k, alpha, beta,h0,alternative,location,scale,dff,prior_analysis,ROPE)
  }

  # table:
  tab.names <- c(
    "TruePositve",
    "FalseNegative",
    "TrueNegative",
    "FalsePositive",
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}


compute.prior.density.re.h1 <- function(rho,h0, k,location,scale,dff,prior_analysis, alpha, beta,alternative,ROPE) {
  if (prior_analysis == "Point") return(rep(NA, length(rho)))
  bound_h1  <- switch(alternative,
                      "greater" = c(a = h0+ROPE, b = 1),
                      "less" = c(a = -1, b = h0+ROPE),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )
  normalizationh1 <- switch(alternative,
                            "two.sided" = switch(prior_analysis,
                                                 "d_beta"       = 1-(p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                                 "beta"         = 1-(p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                                 "Moment"          = {
                                                   (pmom(1-location, tau = scale^2)-pmom(bound_h1[2]-location, tau = scale^2))+
                                                     (pmom(bound_h1[1]-location, tau = scale^2)-pmom(-1-location, tau = scale^2))
                                                 }),

                            "less"  = switch(prior_analysis,
                                             "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                             "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                             "Moment"          = {
                                               (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                             }),
                            "greater"  = switch(prior_analysis,
                                                "d_beta"       = (p_beta(bound_h1[2], 1/k, 1/k,-1,1) - p_beta(bound_h1[1], 1/k, 1/k,-1,1)),
                                                "beta"         = (p_beta(bound_h1[2], alpha, beta,-1,1) - p_beta(bound_h1[1], alpha, beta,-1,1)),
                                                "Moment"          = {
                                                  (pmom(bound_h1[2]-location, tau = scale^2)-pmom(bound_h1[1]-location, tau = scale^2))
                                                })
  )

  #r_prior(rho,k,location,scale,dff,prior_analysis, alpha, beta,min(bound),max(bound)) / normalization

  prior_h1<-r_prior(rho,k,location,scale,dff,prior_analysis, alpha, beta,-1,1)
  switch(alternative,
         "two.sided" = { prior_h1[rho>min(bound_h1)&rho<max(bound_h1)]=0 },
         "greater" = { prior_h1[rho<bound_h1[1]]=0 },
         "less" = { prior_h1[rho>bound_h1[2]]=0 }
  )
  prior_h1

}
compute.prior.density.re.h0 <- function(rho,h0, k,location,scale,dff,prior_analysis, alpha, beta,alternative,ROPE) {
  if (prior_analysis == "Point") return(rep(NA, length(rho)))
  bound_h0  <- switch(alternative,
                      "greater" = c(a = h0, b = h0+ROPE),
                      "less" = c(a = h0+ROPE, b = h0),
                      "two.sided" = c(a = h0+ROPE[1], b = h0+ROPE[2])
  )
  normalizationh0 <- switch(prior_analysis,
                            "d_beta" = p_beta(bound_h0[2], 1/k, 1/k, -1, 1) - p_beta(bound_h0[1], 1/k, 1/k, -1, 1),
                            "beta"   = p_beta(bound_h0[2], alpha, beta, -1, 1) - p_beta(bound_h0[1], alpha, beta, -1, 1),
                            "Moment"    = {pmom(bound_h0[2]-location, tau = scale^2) - pmom(bound_h0[1]-location, tau = scale^2)
                            }
  )

  #r_prior(rho,k,location,scale,dff,prior_analysis, alpha, beta,min(bound),max(bound)) / normalization

  prior_h0<-r_prior(rho,k,location,scale,dff,prior_analysis, alpha, beta,-1,1)
  switch(alternative,
         "two.sided" = { prior_h0[rho<min(bound_h0)|rho>max(bound_h0)]=0 },
         "greater" = { prior_h0[rho>bound_h0[2]]=0 },
         "less" = { prior_h0[rho<bound_h0[1]]=0 }
  )
  prior_h0

}



re_prior_plot <- function(k, alpha, beta, h0,
                          location, scale, dff, prior_analysis, de_an_prior,
                          k_d, alpha_d, beta_d,
                          location_d, scale_d, dff_d,
                          prior_design, alternative, ROPE) {

  # ---- Plot bounds ----
  plot.bounds <- switch(alternative,
                        "greater"   = c(h0, 1),
                        "less"      = c(-1, h0),
                        "two.sided" = c(-1, 1))

  rr <- seq(plot.bounds[1], plot.bounds[2], 0.0025)

  # ---- Compute priors ----
  prior.analysis.h1 <- compute.prior.density.re.h1(
    rr, h0, k, location, scale, dff,
    prior_analysis, alpha, beta, alternative, ROPE
  )

  prior.analysis.h0 <- compute.prior.density.re.h0(
    rr, h0, k, location, scale, dff,
    prior_analysis, alpha, beta, alternative, ROPE
  )

  # ---- Base long-format data ----
  df_lines <- data.frame(
    rr = rep(rr, 2),
    Density = c(prior.analysis.h1, prior.analysis.h0),
    Prior = rep(c("H1 - Analysis Prior",
                  "H0 - Analysis Prior"),
                each = length(rr))
  )

  # ---- Compute design prior if needed ----
  if (de_an_prior == 0 && prior_design != "Point") {

    prior.design <- compute.prior.density.re.h1(
      rr, h0, k_d, location_d, scale_d, dff_d,
      prior_design, alpha_d, beta_d, alternative, ROPE
    )

    df_design <- data.frame(
      rr = rr,
      Density = prior.design,
      Prior = "H1 - Design Prior"
    )

    df_lines <- rbind(df_lines, df_design)
  }

  # ---- Y limits ----
  ylim_max <- max(df_lines$Density[is.finite(df_lines$Density)], na.rm = TRUE)

  # ---- Legend position (match your other functions) ----
  legend_pos <- switch(alternative,
                       "greater"   = c(0.65, 0.95),
                       "two.sided" = c(0.65, 0.95),
                       "less"      = c(0.05, 0.95))

  # ---- Build ggplot ----
  p <- ggplot2::ggplot(df_lines,
                       ggplot2::aes(x = rr,
                                    y = Density,
                                    color = Prior,
                                    linetype = Prior,
                                    linewidth = Prior)) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(
      "H1 - Analysis Prior" = "black",
      "H0 - Analysis Prior" = "black",
      "H1 - Design Prior"   = "gray"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "H1 - Analysis Prior" = "solid",
      "H0 - Analysis Prior" = "dashed",
      "H1 - Design Prior"   = "solid"
    )) +
    ggplot2::scale_linewidth_manual(values = c(
      "H1 - Analysis Prior" = 1.2,
      "H0 - Analysis Prior" = 1.2,
      "H1 - Design Prior"   = 2
    )) +
    ggplot2::labs(
      x = expression(bold(rho)),
      y = "Density",
      title = bquote(bold("Prior distribution on "~rho~" under the alternative"))
    ) +
    ggplot2::coord_cartesian(
      xlim = plot.bounds,
      ylim = c(0, ylim_max)
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.justification = c(0, 1),
      legend.background =
        ggplot2::element_rect(fill = scales::alpha("white", 0.8), color = NA),
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  # ---- Handle Point design prior ----
  if (de_an_prior == 0 && prior_design == "Point") {

    # dummy row for legend
    df_dummy <- data.frame(
      rr = c(NA, NA),
      Density = c(NA, NA),
      Prior = "H1 - Design Prior"
    )

    p <- p +
      ggplot2::geom_line(data = df_dummy,
                         ggplot2::aes(x = rr, y = Density,
                                      color = Prior,
                                      linetype = Prior,
                                      linewidth = Prior),
                         na.rm = TRUE,
                         show.legend = TRUE) +
      ggplot2::annotate("segment",
                        x = location_d,
                        xend = location_d,
                        y = 0,
                        yend = ylim_max,
                        color = "gray",
                        linetype = "dashed",
                        arrow = ggplot2::arrow(
                          length = grid::unit(0.1, "inches")
                        ))
  }

  return(p)
}

re_bf10_p <- function(threshold, n, k, alpha, beta, h0, alternative,
                      location, scale, dff, prior_analysis, ROPE) {

  rr <- seq(-0.99, 0.99, 0.01)

  # Compute BF10 and bounds
  BF10   <- re_BF10(rr, n, k, alpha, beta, h0, alternative,
                    location, scale, dff, prior_analysis, ROPE)
  r.BF10 <- re_BF_bound_10(threshold, n, k, alpha, beta, h0, alternative,
                           location, scale, dff, prior_analysis, ROPE)

  BF01   <- 1 / BF10
  r.BF01 <- re_BF_bound_01(threshold, n, k, alpha, beta, h0, alternative,
                           location, scale, dff, prior_analysis, ROPE)

  max.BF01 <- 1 / re_BF10(h0, n, k, alpha, beta, h0, alternative,
                          location, scale, dff, prior_analysis, ROPE)
  impossible <- (alternative == "two.sided") &&
    (max.BF01 < threshold || identical(r.BF01, "bound cannot be found"))

  ## ---------- Titles ----------
  main.bf10 <- if (length(r.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(threshold) ~ " when r = " ~ .(round(r.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(threshold) ~ " when r = " ~ .(round(r.BF10[1], 2)) ~
                  " or " ~ .(round(r.BF10[2], 2))))
  }

  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(threshold)))
  } else if (length(r.BF01) == 1) {
    bquote(bold("BF"[0][1] ~ "=" ~ .(threshold) ~ " when r = " ~ .(round(r.BF01, 2))))
  } else {
    bquote(bold("BF"[0][1] ~ "=" ~ .(threshold) ~ " when r = " ~ .(round(r.BF01[1], 2)) ~
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +   # all black line
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

Power_re <- function(threshold, k, alpha, beta, h0, alternative,
                     location, scale, dff, prior_analysis,
                     k_d, alpha_d, beta_d, location_d, scale_d, dff_d, prior_design,
                     de_an_prior, N, ROPE) {

  # N range
  N.min <- 4
  N.max <- ceiling(N * 1.2)
  Ns <- seq(N.min, N.max, length.out = 31)

  TPE <- FPE <- TNE <- FNE <- numeric(length(Ns))

  for (i in seq_along(Ns)) {
    r10 <- re_BF_bound_10(threshold, Ns[i], k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis, ROPE)
    r01 <- re_BF_bound_01(threshold, Ns[i], k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis, ROPE)

    TPE[i] <- if (de_an_prior == 1)
      re_TPE(r10, Ns[i], k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis, ROPE) else
        re_TPE(r10, Ns[i], k_d, alpha_d, beta_d, h0, alternative, location_d, scale_d, dff_d, prior_design, ROPE)

    FPE[i] <- re_FPE(r10, Ns[i], k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis, ROPE)
    TNE[i] <- re_TNE(r01, Ns[i], k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis, ROPE)

    FNE[i] <- if (de_an_prior == 1)
      re_FNE(r01, Ns[i], k, alpha, beta, h0, alternative, location, scale, dff, prior_analysis, ROPE) else
        re_FNE(r01, Ns[i], k_d, alpha_d, beta_d, h0, alternative, location_d, scale_d, dff_d, prior_design, ROPE)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme + legend_theme

  p2 <- ggplot2::ggplot(df2, ggplot2::aes(SampleSize, Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
    ) +
    clean_theme + legend_theme

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}



# ---- onesample.r ----

pmom <- function(q,V1=1,tau=1) {

  z <- .5-(stats::pnorm(abs(q)/sqrt(V1*tau)) - abs(q)/sqrt(2*pi*V1*tau) * exp(-.5*q^2/(tau*V1)) - .5)
  return(z*(q<=0)+(1-z)*(q>0))
}
qmom <- function(p, tau) {
  sapply(p, function(prob) {
    if (prob == 0) return(0)
    if (prob == 1) return(Inf)

    objfun <- function(x) {
      pmom(x, tau = tau)- prob
    }
    opt <- stats::uniroot(objfun,lower=-5,upper=5)$root

    return(opt)
  })
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

t1_prior<- function(delta, location, scale, dff, prior_analysis){
  switch(prior_analysis,
         "Cauchy"         = tstude(delta, location, scale, 1),
         "Normal"         = stats::dnorm (delta, location, scale),
         "Moment"            = dMoment  (delta, location, scale),
         "t-distribution" = tstude(delta, location, scale, dff))
}

t1_BF10 <-function(t, df, prior_analysis, location, scale, dff, alternative){
  bound  <- switch(alternative,
                   "greater"  = c(a = 0, b = Inf),
                   "less"  = c(a = -Inf, b = 0),
                   "two.sided" = c(a = -Inf, b = Inf)
  )
  x <- numeric(length(t))

  # Normalize the prior outside the for-loop:
  # normalization  <- stats::integrate(function(delta) t1_prior(delta, location, scale, dff, prior_analysis),lower = bound[1], upper = bound[2])$value
  # For all priors, the prior integrates to 1 when a = -Inf, b = Inf.
  # For all priors, we use their CDFs when either a = 0 or b = 0 and minimize manual integrations.
  # Note: pmom() errors at -Inf and Inf, so we avoid it below.
  normalization <- if (alternative == "two.sided") 1 else
    switch(prior_analysis,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

  for(i in 1:length(t)){
    # int  <- function(delta) stats::dt(t[i], df, ncp = delta * sqrt(df+1)) * t1_prior(delta, location, scale, dff, prior_analysis)/normalization
    int <- function(delta) stats::dt(t[i], df, ncp = delta * sqrt(df + 1)) * t1_prior(delta, location, scale, dff, prior_analysis) / normalization

    # Removed stop.on.error = FALSE as it is bad form.
    # Below, I increased rel.tol. It gives good enough precision, and the app becomes quite faster:
    x[i] <- stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = 1e-5)$value / stats::dt(t[i], df, ncp = 0)
  }
  return(x)
}



# for finding the t value such that BF10 = D (code stats::optimized):
t1_BF10_bound <- function(threshold, df, prior_analysis, location, scale, dff, alternative) {
  Bound_finding <- function(t) t1_BF10(t, df, prior_analysis, location, scale, dff, alternative) - threshold

  x <- tryCatch(stats::uniroot(Bound_finding, lower = -6, upper = 0)$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower =  0, upper = 6)$root, error = function(e) NA)
  results <- c(x, y)

  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")

  BF.vals  <- t1_BF10(results, df, prior_analysis, location, scale, dff, alternative)
  BF.close <- which(round(BF.vals, 2) == round(threshold, 2))
  if (length(BF.close) == 0) return("bound cannot be found")


  return(results[BF.close])
}

# finding the t that correspond to BF01 = D is the same as
# finding the t that corresponds to BF10 = 1/threshold:
t1_BF01_bound <- function(threshold, df, prior_analysis, location, scale, dff, alternative) {
  t1_BF10_bound(1 / threshold, df, prior_analysis, location, scale, dff, alternative)
}

# p(BF01>D|H0)
# t is the t-value lead to BF = b based on the bound functions (stats::optimized):
t1_TNE <- function(t, df,alternative) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (length(t) == 2) return(stats::pt(max(t), df) - stats::pt(min(t), df))

  # length(t) = 1:
  return(if (alternative=="greater") stats::pt(t, df,0) else 1 - stats::pt(t, df))
}

# p(BF10>D|H1)
# Argument 'alternative' is fully determined by the length and sign of the t values.
# I removed it as a function argument and compute it inside t1_TPE() instead.
t1_TPE <- function(t, df, prior_analysis, location, scale, dff) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  alternative <- if (length(t) == 2) "two.sided" else if (t >= 0) "greater" else "less"

  if (prior_analysis == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(min(t), df, ncp) + (1 - pnct(max(t), df, ncp)))
    # Length 1:
    return(if (t >= 0) 1 - pnct(t, df, ncp) else pnct(t, df, ncp))
  }

  bound  <- switch(alternative,
                   "greater"  = c(a = 0,    b = Inf),
                   "less"  = c(a = -Inf, b = 0),
                   "two.sided" = c(a = -Inf, b = Inf))

  normalization <- if (alternative == "two.sided") 1 else
    switch(prior_analysis,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

  int <- if (length(t) == 2) { # two-sided test
    function(delta) {
      pro1 <- 1 - pnct(max(t), df, delta * sqrt(df + 1))
      pro2 <-     pnct(min(t), df, delta * sqrt(df + 1))
      (pro1 + pro2) * t1_prior(delta, location, scale, dff, prior_analysis) / normalization
    }
  } else if (t >= 0) { # one-sided test with delta > 0
    function(delta) (1 - pnct(t, df, delta * sqrt(df + 1))) * t1_prior(delta, location, scale, dff, prior_analysis) / normalization
  } else {             # one-sided test with delta < 0
    function(delta) pnct(t, df, delta * sqrt(df + 1)) * t1_prior(delta, location, scale, dff, prior_analysis) / normalization
  }

  # setting error value such that error are prevented:
  #error <- if (prior_analysis == "Moment" && scale < 0.3) 1e-14 else if (scale > 0.3) .Machine$double.eps^0.25 else 1e-8
  error = 1e-4
  stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value
}

# p(BF01>D|H1)
# Similar as above:
t1_FNE <- function(t, df, prior_analysis, location, scale, dff,alternative){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (prior_analysis == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(max(t), df, ncp) - pnct(min(t), df, ncp))
    # Length 1:
    return(if (alternative=="greater") pnct(t, df, ncp) else 1 - pnct(t, df, ncp))
  }

  bound  <- switch(alternative,
                   "greater"  = c(a = 0,    b = Inf),
                   "less"  = c(a = -Inf, b = 0),
                   "two.sided" = c(a = -Inf, b = Inf))


  normalization <- if (alternative == "two.sided") 1 else
    switch(prior_analysis,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

  int <- if (length(t) == 2) { # two-sided test
    function(delta) {
      pro1 <- pnct(max(t), df, delta * sqrt(df + 1))
      pro2 <- pnct(min(t), df, delta * sqrt(df + 1))
      (pro1 - pro2) * t1_prior(delta, location, scale, dff, prior_analysis) / normalization
    }
  } else if (alternative =="greater") { # one-sided test with delta > 0
    function(delta) pnct(t, df, delta * sqrt(df + 1)) * t1_prior(delta, location, scale, dff, prior_analysis) / normalization
  } else if (alternative =="less"){             # one-sided test with delta < 0
    function(delta) (1 - pnct(t, df, delta * sqrt(df + 1))) * t1_prior(delta, location, scale, dff, prior_analysis) / normalization
  }

  # setting error value such that error are prevented:
  #error <- if (prior_analysis == "Moment" && scale < 0.3) 1e-14 else if (scale > 0.3) .Machine$double.eps^0.25 else 1e-8
  error = 1e-4
  stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value
}

# p(BF10>D|H0)
t1_FPE <- function(t, df,alternative) {
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  # if (length(t) == 4) t <- t[2:3]

  if (length(t) == 2) return(stats::pt(min(t), df) + (1 - stats::pt(max(t), df)))

  # length(t) = 1:
  return(if (alternative=="greater") 1 - stats::pt(t, df) else stats::pt(t, df))
}




# Finding the degree of freedom that ensure p(BF10>D|H1) > targeted probability:
t1_N_finder <- function(threshold, true_rate, prior_analysis, location, scale, dff, alternative,
                        prior_design, location_d, scale_d, dff_d, de_an_prior, false_rate) {
  #de_an_prior: 1 = design prior and analysis priors are the same, otherwise different.

  # error prevention
  # sometimes, power can go higher than .8 with N= 2 already.
  # So, N should be returned now, otherwise, error will occur later.
  # Jorge: Below, I added design prior for N=2 too.
  lower <- 2
  upper <- 10000

  t2 <- t1_BF10_bound(threshold, df = lower, prior_analysis, location, scale, dff,alternative)
  p2 <- if (de_an_prior == 1)
    t1_TPE(t2, df = lower, prior_analysis, location, scale, dff ) else
      t1_TPE(t2, df = lower, prior_design, location_d, scale_d, dff_d)
  if (p2 > true_rate) return(lower)
  Power_root <- function(df) {
    t <- t1_BF10_bound(threshold, df, prior_analysis, location, scale, dff,alternative)
    pro <- if (de_an_prior == 1)
      t1_TPE(t, df, prior_analysis, location, scale, dff)  else
        t1_TPE(t, df, prior_design, location_d, scale_d, dff_d)
    pro-true_rate
  }

  ## finding the required df, i will do the plus one to get the N in the later function.
  # Jorge: It's a pity if we don't fix it here already, super easy to do right now.
  #        So I did it, see below. I'll adapt latter functions if needed be.

  # Jorge: 'df.power' makes for a more accurate name than 'N'.message("Power at lower = ", Power_root(lower))
  df.power <- stats::uniroot(Power_root, lower = lower, upper = upper)$root

  ## checking if the N lead to an acceptable alpha level
  t   <- t1_BF10_bound(threshold, df.power, prior_analysis, location, scale, dff,alternative)
  FPE <- t1_FPE(t, df.power,alternative)
  if (FPE <= false_rate) return(df.power + 1)

  # if the FPE > false_rate, then we search for another df
  # Jorge: 'alpha.root' is better than 'alpha_bound'.
  alpha.root <- function(df) {
    t <- t1_BF10_bound(threshold, df, prior_analysis, location, scale, dff, alternative)
    t1_FPE(t, df,alternative) - false_rate
  }

  # Jorge: 'df.alpha' is better than 'NN'.
  df.alpha <- stats::uniroot(alpha.root, lower = df.power, upper = upper)$root
  return(df.alpha + 1)
}

t1_N_01_finder <- function(threshold, true_rate, prior_analysis, location, scale, dff, alternative,
                           prior_design, location_d, scale_d, dff_d, de_an_prior, false_rate) {
  lower <- 2
  upper <- 10000

  t2 <- t1_BF01_bound(threshold, df = lower, prior_analysis, location, scale, dff,alternative)
  TNE_lo <- t1_TNE(t2, df = lower,alternative)
  if (TNE_lo > true_rate) return(lower)

  FNE_lo <-  if (de_an_prior == 1)
    t1_FNE(t2, df = lower, prior_analysis, location, scale, dff,alternative ) else
      t1_FNE(t2, df = lower, prior_design, location_d, scale_d, dff_d,alternative)
  if (TNE_lo > true_rate&FNE_lo<false_rate) return(lower)


  TN_root <- function(df) {
    t <- t1_BF01_bound(threshold, df, prior_analysis, location, scale, dff,alternative)
    t1_TNE(t, df = df,alternative) - true_rate
  }

  df.TN <- stats::uniroot(TN_root, lower = lower, upper = upper)$root

  t   <- t1_BF01_bound(threshold, df.TN, prior_analysis, location, scale, dff,alternative)
  FNE <- if (de_an_prior == 1)
    t1_FNE(t, df = df.TN, prior_analysis, location, scale, dff,alternative ) else
      t1_FNE(t, df = df.TN, prior_design, location_d, scale_d, dff_d,alternative)

  if (FNE <= false_rate) return(df.TN + 1)

  FN.root <- function(df) {
    t <- t1_BF01_bound(threshold, df, prior_analysis, location, scale, dff, alternative)
    pro = if (de_an_prior == 1)
      t1_FNE(t, df = df, prior_analysis, location, scale, dff ,alternative) else
        t1_FNE(t, df = df, prior_design, location_d, scale_d, dff_d,alternative)
    pro- false_rate
  }
  df.FN <- stats::uniroot(FN.root, lower = df.TN, upper = upper)$root
  return(df.FN + 1)
}




############ probability table
# Jorge: I edited so that it used N returned by t1_N_finder().
t1_Table <- function(threshold, true_rate, prior_analysis, location, scale, dff, alternative,
                     prior_design, location_d, scale_d, dff_d, de_an_prior, N, mode_bf, false_rate,type_rate) {

  df <- if (mode_bf == "0") {
    N - 1
  } else {
    fun <- if (type_rate == "positive") t1_N_finder else t1_N_01_finder
    ceiling(fun(threshold, true_rate, prior_analysis, location, scale, dff, alternative,
                prior_design, location_d, scale_d, dff_d, de_an_prior, false_rate)) - 1
  }

  # t bounds:
  t10 <- t1_BF10_bound(threshold, df, prior_analysis, location, scale, dff, alternative)
  t01 <- t1_BF01_bound(threshold, df, prior_analysis, location, scale, dff, alternative)

  # max BF10 possible:
  max_BF <- 1 / t1_BF10(0, df, prior_analysis, location, scale, dff, alternative)
  BF_D   <- t10

  # FPE and TPE:
  FPE       <- t1_FPE(t10, df,alternative)
  if (de_an_prior == 1) {
    TPE       <- t1_TPE(t10, df, prior_analysis, location, scale, dff)
    TPR_prior <- prior_analysis
    TPR_loc   <- location
    TPR_scale <- scale
    TPR_dff   <- dff
  } else {
    TPE       <- t1_TPE(t10, df, prior_design, location_d, scale_d, dff_d)
    TPR_prior <- prior_design
    TPR_loc   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }

  # FNE and TNE:
  if (any(alternative == "two.sided" & max_BF < threshold | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t1_FNE(t01, df, TPR_prior, TPR_loc, TPR_scale, TPR_dff,alternative)
    TNE <- t1_TNE(t01, df,alternative)
  }

  # table:
  tab.names <- c(
    "TruePositve",
    "FalseNegative",
    "TrueNegative",
    "FalsePositive",
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, df+1, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}



# For plotting, compute normalized prior density over tt:
compute.prior.density.t <- function(tt, prior_analysis, location, scale, dff, alternative) {
  if (prior_analysis == "Point") return(rep(NA, length(tt)))
  bounds <- switch(alternative,
                   "greater"  = c(0, Inf),
                   "less"  = c(-Inf, 0),
                   "two.sided" = c(-Inf, Inf))
  norm <- stats::integrate(function(delta) t1_prior(delta, location, scale, dff, prior_analysis),
                           lower = bounds[1], upper = bounds[2])$value
  t1_prior(tt, location, scale, dff, prior_analysis) / norm
}

# plot for the selected prior
t1_prior_plot <- function( prior_analysis, location, scale, dff, alternative,
                          prior_design, location_d, scale_d, dff_d, de_an_prior) {
  # ---- Determine plotting bounds ----
  plot.bounds <- switch(alternative,
                        "greater" = c(0, 5),
                        "less" = c(-5, 0),
                        "two.sided" = c(-5, 5))
  tt <- seq(plot.bounds[1], plot.bounds[2], 0.01)

  # ---- Compute Analysis prior ----
  prior_analysis_dens <- compute.prior.density.t(tt, prior_analysis, location, scale, dff, alternative)

  # ---- Base data frame ----
  df <- data.frame(
    t = tt,
    Density = prior_analysis_dens,
    Prior = "H1 - Analysis Prior"
  )

  # ---- Conditionally add Design prior ----
  if (de_an_prior == 0) {
    if (prior_design == "Point") {
      # Dummy line for legend only
      df_design <- data.frame(
        t = c(NA, NA),
        Density = c(NA, NA),
        Prior = "H1 - Design Prior"
      )
      df <- rbind(df, df_design)
    } else {
      prior_design_dens <- compute.prior.density.t(tt, prior_design, location_d, scale_d, dff_d, alternative)
      df_design <- data.frame(
        t = tt,
        Density = prior_design_dens,
        Prior = "H1 - Design Prior"
      )
      df <- rbind(df, df_design)
    }
  }

  # ---- Legend position (match t1e_prior_plot) ----
  legend_pos <- switch(alternative,
                       "greater" = c(0.65, 0.95),
                       "two.sided" = c(0.65, 0.95),
                       "less" = c(0.05, 0.95))

  # ---- Base ggplot ----
  p <- ggplot2::ggplot(df, ggplot2::aes(x = t, y = Density, color = Prior, linetype = Prior)) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("H1 - Analysis Prior" = "black",
                                           "H1 - Design Prior" = "gray")) +
    ggplot2::scale_linetype_manual(values = c("H1 - Analysis Prior" = "solid",
                                              "H1 - Design Prior" = "dashed")) +
    ggplot2::labs(
      x = expression(bold(delta)),
      y = "density",
      title = bquote(bold("Prior distribution on "~delta~" under the alternative"))
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.8), color = NA),
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linewidth = 1.5))
    )

  # ---- Add vertical segment for Point prior (like t1e_prior_plot) ----
  if (de_an_prior == 0 && prior_design == "Point") {
    ylim_max <- max(df$Density, na.rm = TRUE)
    # Dummy line for legend
    p <- p +
      ggplot2::geom_line(data = data.frame(t = c(NA, NA), Density = c(NA, NA), Prior = "H1 - Design Prior"),
                         ggplot2::aes(x = t, y = Density, color = Prior, linetype = Prior),
                         na.rm = TRUE, show.legend = TRUE) +
      # Vertical gray line using annotate
      ggplot2::annotate("segment",
                        x = location_d, xend = location_d,
                        y = 0, yend = ylim_max,
                        color = "gray",
                        linetype = "dashed",
                        arrow = ggplot2::arrow(length = grid::unit(0.1, "inches")))
  }

  # ---- Set x/y limits ----
  ylim_max <- max(df$Density, na.rm = TRUE)
  p <- p + ggplot2::coord_cartesian(ylim = c(0, ylim_max), xlim = plot.bounds)

  return(p)
}



# plots for showing the relationship between BF and t-values

bf10_t1 <- function(threshold = 3, df,
                    prior_analysis = "NA", location = 0, scale = 0.707,
                    dff = 1, alternative) {

  tt <- seq(-5, 5, 0.2)

  # Compute BF10 and bounds
  BF10   <- t1_BF10(tt, df, prior_analysis, location, scale, dff, alternative)
  t.BF10 <- t1_BF10_bound(threshold, df, prior_analysis, location, scale, dff, alternative)

  BF01   <- 1 / BF10
  t.BF01 <- t1_BF01_bound(threshold, df, prior_analysis, location, scale, dff, alternative)

  # BF10 title
  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10] ~ "=" ~ .(threshold) ~ " when t = " ~ .(round(t.BF10, 2))))
  } else {
    bquote(bold("BF"[10] ~ "=" ~ .(threshold) ~ " when t = " ~ .(round(t.BF10[1], 2)) ~
                  " or " ~ .(round(t.BF10[2], 2))))
  }

  # Check if BF01 is impossible
  impossible <- (alternative == "two.sided") &&
    (max(BF01) < threshold || identical(t.BF01, "bound cannot be found"))

  # BF01 title
  main.bf01 <- if (impossible) {
    bquote(bold("It is impossible to have BF"[01] ~ "=" ~ .(threshold)))
  } else if (length(t.BF01) == 1) {
    bquote(bold("BF"[0][1] ~ "=" ~ .(threshold) ~ " when t = " ~ .(round(t.BF01, 2))))
  } else {
    bquote(bold("BF"[0][1] ~ "=" ~ .(threshold) ~ " when t = " ~ .(round(t.BF01[1], 2)) ~
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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
    ggplot2::geom_line(linewidth = 1.2, color = "black") +   # ← FIXED LINE COLOR
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


# Power curve function for BF10 > threshold under H1:
Power_t1 <- function(threshold, prior_analysis, location, scale, dff, alternative,
                     prior_design, location_d, scale_d, dff_d,
                     de_an_prior, N) {

  # df range
  df.min <- 2
  df.max <- ceiling(N * 1.2)
  dfs <- seq(df.min, df.max, length.out = 31)

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(dfs))

  for (i in seq_along(dfs)) {
    t10 <- t1_BF10_bound(threshold, dfs[i], prior_analysis, location, scale, dff, alternative)
    t01 <- t1_BF01_bound(threshold, dfs[i], prior_analysis, location, scale, dff, alternative)

    TPE[i] <- if (de_an_prior == 1) {
      t1_TPE(t10, dfs[i], prior_analysis, location, scale, dff)
    } else {
      t1_TPE(t10, dfs[i], prior_design, location_d, scale_d, dff_d)
    }

    FPE[i] <- t1_FPE(t10, dfs[i], alternative)
    TNE[i] <- t1_TNE(t01, dfs[i], alternative)
    FNE[i] <- if (de_an_prior == 1) {
      t1_FNE(t01, dfs[i], prior_analysis, location, scale, dff, alternative)
    } else {
      t1_FNE(t01, dfs[i], prior_design, location_d, scale_d, dff_d, alternative)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # Plot BF01
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
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

te_prior<- function(delta,location,scale,dff,prior_analysis){

  switch(prior_analysis,
         "Cauchy"        = tstude(delta,location, scale,1),
         "Normal"        = stats::dnorm(delta,location,scale),
         "Moment"            = dMoment(delta,location,scale),
         "t-distribution" = tstude(delta,location,scale,dff))


}
norm_h1 <- function(alternative, prior_analysis, bound_h1, location, scale, dff = NULL) {
  normalizationh1 <- switch(
    alternative,

    "two.sided" = 1 - switch(
      prior_analysis,
      "Cauchy"         = stats::pcauchy(bound_h1[2], location, scale) -
        stats::pcauchy(bound_h1[1], location, scale),
      "Normal"         = stats::pnorm(bound_h1[2], location, scale) -
        stats::pnorm(bound_h1[1], location, scale),
      "Moment"            = pmom(bound_h1[2] - location, tau = scale^2) -
        pmom(bound_h1[1] - location, tau = scale^2),
      "t-distribution" = stats::pt((bound_h1[2] - location)/scale, df = dff) -
        stats::pt((bound_h1[1] - location)/scale, df = dff)
    ),

    "less" = switch(
      prior_analysis,
      "Cauchy"         = stats::pcauchy(bound_h1[2], location, scale) -
        stats::pcauchy(bound_h1[1], location, scale),
      "Normal"         = stats::pnorm(bound_h1[2], location, scale) -
        stats::pnorm(bound_h1[1], location, scale),
      "Moment"            = pmom(bound_h1[2] - location, tau = scale^2),
      "t-distribution" = stats::pt((bound_h1[2] - location)/scale, df = dff) -
        stats::pt((bound_h1[1] - location)/scale, df = dff)
    ),

    "greater" = switch(
      prior_analysis,
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
norm_h0 <- function(prior_analysis, bound_h0, location, scale, dff = NULL) {
  normalizationh0 <- switch(
    prior_analysis,

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


t1e_BF10i <-function(t,df,prior_analysis ,location, scale,dff , alternative,ROPE ){
  bound_h1  <- switch(alternative,
                      "two.sided" = c(a = ROPE[1], b = ROPE[2]),
                      "greater" = c(a = ROPE, b = Inf),
                      "less" = c(a = -Inf, b = ROPE)
  )

  bound_h0  <- switch(alternative,
                      "two.sided" = c(a = ROPE[1], b = ROPE[2]),
                      "greater" = c(a = 0, b = ROPE),
                      "less" = c(a = ROPE, b = 0)
  )

  normalizationh1 <- norm_h1(alternative, prior_analysis, bound_h1, location, scale, dff)



  # H0 Normalization
  normalizationh0 <- norm_h0(prior_analysis, bound_h0, location, scale, dff)

  int  <- function(delta){
    stats::dt(t,df,ncp = delta *sqrt(df+1))* te_prior(delta,location,scale,dff,prior_analysis)/normalizationh1
  }

  error = 1e-4

  if (alternative == "two.sided"){
    lh1 = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value
  }else{
    lh1 = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value

  }


  int  <- function(delta){
    stats::dt(t,df,ncp = delta *sqrt(df+1))* te_prior(delta,location,scale,dff,prior_analysis)/normalizationh0}

  lh0 = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(lh1/lh0)
}

t1e_BF10 <-function(t,df,prior_analysis,location, scale,dff , alternative,ROPE ){

  x <- sapply(t, function(ti) t1e_BF10i(ti, df, prior_analysis,location, scale, dff, alternative, ROPE))
  return(x)
}
#
t1e_BF10_bound <-function(threshold, df,prior_analysis,location,scale,dff , alternative,ROPE){
  y <- numeric(0)
  Bound_finding <-function(t){
    t1e_BF10(t,df,prior_analysis,location,scale,dff , alternative,ROPE )- threshold
  }

  switch(alternative,
         "two.sided" ={
           x <- tryCatch(stats::uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
           y <- tryCatch(stats::uniroot(Bound_finding, lower =  0, upper = 20)$root, error = function(e) NA)
         },
         "greater"={
           x <- tryCatch(stats::uniroot(Bound_finding, lower = 0, upper = 20)$root, error = function(e) NA)
         },
         "less" = {
           x <- tryCatch(stats::uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
         })


  results <- c(x, y)

  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")

  BF.vals  <- t1e_BF10(results,df,prior_analysis,location,scale,dff , alternative,ROPE )
  BF.close <- which(round(BF.vals, 2) == round(threshold, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])

}


t1e_BF01_bound <-function(threshold, df,prior_analysis,location,scale,dff , alternative,ROPE){
  t1e_BF10_bound(1/threshold, df,prior_analysis,location,scale,dff , alternative,ROPE)
}



t1e_TPE <-function(t,df,prior_analysis ,location,scale,dff , alternative ,ROPE){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (prior_analysis == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(min(t), df, ncp) + (1 - pnct(max(t), df, ncp)))
    # Length 1:
    return(if (t >= 0) 1 - pnct(t, df, ncp) else pnct(t, df, ncp))
  }


  bound_h1  <- switch(alternative,
                      "greater" = c(a = ROPE, b = Inf),
                      "less" = c(a = -Inf, b = ROPE),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )

  normalizationh1 <- norm_h1(alternative, prior_analysis, bound_h1, location, scale, dff)


  x = NULL

  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)

    pro <- switch(alternative,
                  "two.sided" = pnct(max(t), df, ncp = ncp, lower  = FALSE) +
                    pnct(min(t), df, ncp = ncp, lower  = TRUE),
                  "greater"  = pnct(t, df, ncp = ncp, lower  = FALSE),
                  "less"  = pnct(t, df, ncp = ncp, lower  = TRUE)
    )

    pro * te_prior(delta, location,scale, dff, prior_analysis) / normalizationh1
  }

  error = 1e-4

  x <- switch(alternative,
              "two.sided" = stats::integrate(int, -Inf, bound_h1[1], rel.tol = error)$value +
                stats::integrate(int, bound_h1[2], Inf, rel.tol = error)$value,
              "less"  = ,
              "greater"  = stats::integrate(int, bound_h1[1], bound_h1[2], rel.tol = error)$value

  )
  return(x)

}

t1e_FNE <-function(t,df,prior_analysis ,location,scale,dff , alternative ,ROPE){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (prior_analysis == "Point") {
    ncp <- location * sqrt(df + 1)
    if (length(t) == 2) return(pnct(max(t), df, ncp) - pnct(min(t), df, ncp))
    # Length 1:
    return(if (t >= 0) pnct(t, df, ncp) else 1 - pnct(t, df, ncp))
  }

  bound_h1  <- switch(alternative,
                      "greater" = c(a = ROPE, b = Inf),
                      "less" = c(a = -Inf, b = ROPE),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )

  normalizationh1 <- norm_h1(alternative, prior_analysis, bound_h1, location, scale, dff)

  x = NULL

  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)

    pro <- switch(alternative,
                  "two.sided" = pnct(max(t), df, ncp = ncp, lower  = TRUE) -
                    pnct(min(t), df, ncp = ncp, lower  = TRUE),
                  "greater"  = pnct(t, df, ncp = ncp, lower  = TRUE),
                  "less"  = pnct(t,df, ncp = ncp, lower  = FALSE)
    )

    pro * te_prior(delta,location, scale, dff, prior_analysis) / normalizationh1
  }

  error = 1e-4

  x <- switch(alternative,
              "two.sided" = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value ,
              "less"  = ,
              "greater"  = stats::integrate(int, bound_h1[1], bound_h1[2], rel.tol = error)$value)
  return(x)

}

t1e_TNE <-function(t,df,prior_analysis ,location,scale,dff , alternative ,ROPE){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  bound_h0  <- switch(alternative,
                      "greater" = c(a = 0, b = ROPE),
                      "less" = c(a = ROPE, b = 0),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )
  normalizationh0 <- norm_h0(prior_analysis, bound_h0, location, scale, dff)

  x = NULL

  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)
    pro <- switch(alternative,
                  "two.sided" = pnct(max(t), df, ncp, lower  = TRUE) -
                    pnct(min(t), df, ncp, lower  = TRUE),
                  "greater"  = pnct(t, df, ncp, lower  = TRUE),
                  "less"  = pnct(t, df, ncp, lower  = FALSE),
                  stop("Unsupported alternative")
    )

    pro * te_prior(delta,location, scale, dff, prior_analysis) / normalizationh0
  }

  error = 1e-4

  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error)$value

  return(x)

}

t1e_FPE <-function(t,df,prior_analysis ,location,scale,dff , alternative ,ROPE){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)
  bound_h0  <- switch(alternative,
                      "greater" = c(a = 0, b = ROPE),
                      "less" = c(a = ROPE, b = 0),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )

  normalizationh0 <- norm_h0(prior_analysis, bound_h0, location, scale, dff)
  x = NULL
  int <- function(delta) {
    ncp <- delta * sqrt(df + 1)

    pro <- switch(alternative,
                  "two.sided" = pnct(max(t), df, ncp, lower  = FALSE) + pnct(min(t), df, ncp, lower  = TRUE),
                  "greater"  = pnct(t, df, ncp, lower  = FALSE),
                  "less"  = pnct(t, df, ncp, lower  = TRUE),
                  stop("Unsupported alternative")
    )

    pro * te_prior(delta,location, scale, dff, prior_analysis) / normalizationh0
  }

  error = 1e-4

  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value

  return(x)

}

t1e_N_finder<-function(threshold,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                       prior_design,scale_d,dff_d, de_an_prior,location_d  ,false_rate){

  lower <- 2
  upper <- 10000
  t2 <-t1e_BF10_bound(threshold, lower,prior_analysis,location,scale,dff , alternative,ROPE)
  p2 <- if (de_an_prior == 1)
    t1e_TPE(t2,lower,prior_analysis ,location,scale,dff , alternative ,ROPE) else
      t1e_TPE(t2,lower,prior_analysis,location_d ,scale_d,dff_d , alternative ,ROPE)
  if (p2 > true_rate) return(lower)

  Power_root <- function(df) {

    t <- t1e_BF10_bound(threshold, df, prior_analysis,location, scale, dff, alternative, ROPE)

    pro <- if (de_an_prior == 1) {
      t1e_TPE(t, df, prior_analysis,location, scale, dff, alternative, ROPE)
    } else {
      t1e_TPE(t, df, prior_design,location_d, scale_d, dff_d, alternative, ROPE)
    }

    true_rate - pro
  }

  df.power <- robust_uniroot(Power_root, lower = 2)
  t <- t1e_BF10_bound(threshold,df.power,prior_analysis,location,scale,dff,alternative ,ROPE )
  FPE <-t1e_FPE(t,df.power,prior_analysis,location ,scale,dff , alternative ,ROPE)
  if (FPE <= false_rate) return(df.power + 1)

  alpha.root <- function(df) {
    t <- t1e_BF10_bound(threshold,df,prior_analysis,location,scale,dff,alternative ,ROPE )
    pro <- t1e_FPE(t , df , prior_analysis,location , scale,dff, alternative,ROPE)
    return(pro - false_rate)
  }
  df.alpha <- stats::uniroot(alpha.root, lower = df.power, upper = upper)$root
  return(df.alpha+1)

}

t1e_N_01_finder<-function(threshold,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                          prior_design,scale_d,dff_d, de_an_prior,location_d  ,false_rate){

  lower <- 10
  upper <- 10000
  t2 <-t1e_BF01_bound(threshold, lower,prior_analysis,location,scale,dff , alternative,ROPE)
  TNE_lo <-t1e_TNE(t2,lower,prior_analysis,location ,scale,dff , alternative ,ROPE)
  if (TNE_lo > true_rate) return(lower)

  FNE_lo <- if (de_an_prior == 1)
    t1e_FPE(t2,lower,prior_analysis,location ,scale,dff , alternative ,ROPE) else
      t1e_FPE(t2,lower,prior_analysis,location_d ,scale_d,dff_d , alternative ,ROPE)
  if (TNE_lo > true_rate&FNE_lo<false_rate) return(lower)

  TN_root <- function(df) {
    t <- t1e_BF01_bound(threshold, df, prior_analysis,location, scale, dff, alternative, ROPE)

    pro <-t1e_TNE(t,df,prior_analysis,location ,scale,dff , alternative ,ROPE)
    true_rate-pro
  }

  df.TN <- robust_uniroot(TN_root, lower = lower)
  t <- t1e_BF01_bound(threshold,df.TN,prior_analysis,location,scale,dff,alternative ,ROPE )
  FNE <-t1e_FNE(t,df.TN,prior_analysis,location ,scale,dff , alternative ,ROPE)
  if (FNE <= false_rate) return(df.TN + 1)

  FN.root <- function(df) {
    t <- t1e_BF01_bound(threshold,df,prior_analysis,location,scale,dff,alternative ,ROPE )
    pro <-   if (de_an_prior == 1) {
      t1e_FNE(t, df, prior_analysis,location, scale, dff, alternative, ROPE)
    } else {
      t1e_FNE(t, df, prior_design,location_d, scale_d, dff_d, alternative, ROPE)
    }
    return(pro - false_rate)
  }
  df.FN <- robust_uniroot(FN.root, lower = df.TN )
  return(df.FN+1)

}
t1e_table<-function(threshold,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                    prior_design,scale_d,dff_d, de_an_prior,N,mode_bf,location_d ,false_rate,type_rate ){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)

  df <- if (mode_bf == 0) {
    N - 1
  } else {
    fun <- if (type_rate == "positive") t1e_N_finder else t1e_N_01_finder
    ceiling(fun(threshold,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                prior_design,scale_d,dff_d, de_an_prior ,location_d,false_rate )) - 1
  }



  # t bounds:
  t10 <- t1e_BF10_bound(threshold, df,prior_analysis,location,scale,dff , alternative,ROPE)
  t01 <- t1e_BF01_bound(threshold, df,prior_analysis,location,scale,dff , alternative,ROPE)

  # max BF10 possible:
  max_BF <- 1 / t1e_BF10(0,df,prior_analysis,location ,scale,dff , alternative,ROPE )
  BF_D   <- t10

  # FPE and TPE:
  FPE       <- t1e_FPE(t10,df,prior_analysis ,location,scale,dff , alternative ,ROPE)
  if (de_an_prior == 1) {
    TPE         <- t1e_TPE(t10,df,prior_analysis ,location,scale,dff , alternative ,ROPE)
    TPR_prior   <- prior_analysis
    TPR_scale   <- scale
    TPR_dff     <- dff
    TPR_location<- location
  } else {
    TPE       <- t1e_TPE(t10,df,prior_design ,location_d,scale_d,dff_d , alternative ,ROPE)
    TPR_prior <- prior_design
    TPR_location   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }

  # FNE and TNE:
  if (any(alternative == "two.sided" & max_BF < threshold | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t1e_FNE(t01,df,TPR_prior,TPR_location ,TPR_scale,TPR_dff , alternative ,ROPE)
    TNE <-  t1e_TNE(t01,df,prior_analysis,location ,scale,dff , alternative ,ROPE)
  }
  # table:
  tab.names <- c(
    "TruePositve",
    "FalseNegative",
    "TrueNegative",
    "FalsePositive",
    "Required N"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, df+1, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table

}

compute.prior.density.te.h1 <- function(tt, prior_analysis,location, scale, dff, alternative,ROPE) {
  if (prior_analysis == "Point") return(rep(NA, length(tt)))
  bound_h1  <- switch(alternative,
                      "greater" = c(a = ROPE, b = Inf),
                      "less" = c(a = -Inf, b = -ROPE),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )

  normalizationh1 <- norm_h1(alternative, prior_analysis, bound_h1, location, scale, dff)


  #prior_h1<-te_prior(tt,scale,dff,prior_analysis) / normalizationh1
  prior_h1<-te_prior(tt,location,scale,dff,prior_analysis)
  switch(alternative,
         "two.sided" = { prior_h1[tt>min(bound_h1)&tt<max(bound_h1)]=0 },
         "greater" = { prior_h1[tt<bound_h1[1]]=0 },
         "less" = { prior_h1[tt>bound_h1[2]]=0 }
  )
  prior_h1
}

compute.prior.density.te.h0 <- function(tt, prior_analysis,location, scale, dff, alternative,ROPE) {
  if (prior_analysis == "Point") return(rep(NA, length(tt)))
  bound_h0  <- switch(alternative,
                      "two.sided" = c(a = ROPE[1], b = ROPE[2]),
                      "greater" = c(a = 0, b = ROPE),
                      "less" = c(a = ROPE, b = 0)
  )
  # H0 Normalization
  normalizationh0 <- norm_h0(prior_analysis, bound_h0, location, scale, dff)


  #prior_h0 <- te_prior(tt,scale,dff,prior_analysis) / normalizationh0
  prior_h0 <- te_prior(tt,location,scale,dff,prior_analysis)
  switch(alternative,
         "two.sided" = { prior_h0[!(tt>min(bound_h0)&tt<max(bound_h0))]=0},
         "greater" = { prior_h0[tt>bound_h0[2]]=0 },
         "less" = { prior_h0[tt<bound_h0[1]]=0 }

  )
  prior_h0
}



t1e_prior_plot <- function(prior_analysis, location, scale, dff, alternative, ROPE,
                           de_an_prior, prior_design, scale_d, dff_d, location_d) {

  # Plot bounds
  plot.bounds <- switch(alternative,
                        "greater" = c(0, 5),
                        "less" = c(-5, 0),
                        "two.sided" = c(-5, 5))
  tt <- seq(plot.bounds[1], plot.bounds[2], 0.01)

  # Compute H1 and H0 priors
  prior_h1 <- compute.prior.density.te.h1(tt, prior_analysis, location, scale, dff, alternative, ROPE)
  prior_h0 <- compute.prior.density.te.h0(tt, prior_analysis, location, scale, dff, alternative, ROPE)

  # Base long-format data for H1/H0
  df_lines <- data.frame(
    tt = rep(tt, 2),
    Density = c(prior_h1, prior_h0),
    Prior = rep(c("H1 - Analysis Prior", "H0 - Analysis Prior"), each = length(tt))
  )

  # Determine legend position
  legend_pos <- switch(alternative,
                       "greater" = c(0.75, 0.95),
                       "two.sided" = c(0.75, 0.95),
                       "less" = c(0.05, 0.95))

  # Start ggplot with H1/H0 lines
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = df_lines,
      ggplot2::aes(x = tt, y = Density, color = Prior, linetype = Prior, linewidth = Prior)
    ) +
    ggplot2::scale_color_manual(values = c(
      "H1 - Analysis Prior" = "black",
      "H0 - Analysis Prior" = "black",
      "H1 - Design Prior" = "gray"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "H1 - Analysis Prior" = "solid",
      "H0 - Analysis Prior" = "dashed",
      "H1 - Design Prior" = "solid"
    )) +
    ggplot2::scale_linewidth_manual(values = c(
      "H1 - Analysis Prior" = 1.2,
      "H0 - Analysis Prior" = 1.2,
      "H1 - Design Prior" = 2
    ))+
    ggplot2::labs(
      x = expression(bold(delta)),
      y = "Density",
      title = bquote(bold("Prior distribution on "~delta~" under the alternative"))
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank() # remove all grid lines
    )

  # Add design prior line if non-point
  if (de_an_prior == 0 && prior_design != "Point") {
    prior_design_vals <- compute.prior.density.te.h1(tt, prior_design, location_d, scale_d, dff_d, alternative, ROPE)
    df_design <- data.frame(
      tt = tt,
      Density = prior_design_vals,
      Prior = "H1 - Design Prior"
    )
    df_design <- df_design[!is.na(df_design$Density), ]
    p <- p + ggplot2::geom_line(data = df_design,
                                ggplot2::aes(x = tt, y = Density, color = Prior, linetype = Prior, linewidth = Prior))
  }

  # Add vertical arrow for point prior without affecting H1/H0 lines
  if (de_an_prior == 0 && prior_design == "Point") {
    ylim_max <- max(prior_h1, prior_h0, na.rm = TRUE)
    # Add invisible line for legend only
    df_dummy <- data.frame(tt = c(NA, NA), Density = c(NA, NA), Prior = "H1 - Design Prior")
    p <- p +
      ggplot2::geom_line(data = df_dummy,
                         ggplot2::aes(x = tt, y = Density, color = Prior, linetype = Prior, linewidth = Prior),
                         na.rm = TRUE, show.legend = TRUE) +
      # vertical arrow
      ggplot2::geom_segment(ggplot2::aes(x = location_d, xend = location_d, y = 0, yend = ylim_max),
                            color = "gray", linetype = "dashed",
                            arrow = ggplot2::arrow(length = grid::unit(0.1, "inches")))
  }

  return(p)
}


te1_BF<- function(threshold, df,
                  prior_analysis, location, scale, dff, alternative, ROPE) {

  tt <- seq(-5, 5, 0.2)

  ## ---------- BF10 ----------
  BF10   <- t1e_BF10(tt, df, prior_analysis, location, scale, dff, alternative, ROPE)
  t.BF10 <- t1e_BF10_bound(threshold, df, prior_analysis, location, scale, dff, alternative, ROPE)

  df10 <- data.frame(t = tt, BF = BF10)

  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(threshold)~" when t = "~.(round(t.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(threshold)~" when t = "~.(round(t.BF10[1], 2))~
                  " or "~.(round(t.BF10[2], 2))))
  }

  x_breaks_10 <- sort(unique(c(-5, 5, round(t.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(t, BF)) +
    ggplot2::geom_line(linewidth = 1.1) +
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
  t.BF01 <- t1e_BF01_bound(threshold, df, prior_analysis, location, scale, dff, alternative, ROPE)

  df01 <- data.frame(t = tt, BF = BF01)

  max.BF01   <- 1 / t1e_BF10(0, df, prior_analysis, location, scale, dff, alternative, ROPE)
  impossible <- (max.BF01 < threshold || identical(t.BF01, "bound cannot be found"))

  if (impossible) {

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(
        limits = c(-5, 5),
        breaks = c(-5, 5)
      ) +
      ggplot2::labs(
        x = "t-value",
        y = bquote("BF"['01'] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[01]~"="~.(threshold)))
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
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when t = "~.(round(t.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when t = "~.(round(t.BF01[1], 2))~
                    " or "~.(round(t.BF01[2], 2))))
    }

    x_breaks_01 <- sort(unique(c(-5, 5, round(t.BF01, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(linewidth = 1.1) +
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

Power_t1e<- function(threshold, prior_analysis, location, scale, dff, alternative,
                     prior_design, location_d, scale_d, dff_d,
                     de_an_prior, N, ROPE) {

  # df range
  df.min <- 2
  df.max <- ceiling(N * 1.2)
  dfs <- seq(df.min, df.max, length.out = 30)

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(dfs))

  for (i in seq_along(dfs)) {
    t10 <- t1e_BF10_bound(threshold, dfs[i], prior_analysis, location, scale, dff, alternative, ROPE)
    t01 <- t1e_BF01_bound(threshold, dfs[i], prior_analysis, location, scale, dff, alternative, ROPE)

    # Choose correct design prior
    TPE[i] <- if (de_an_prior == 1) {
      t1e_TPE(t10, dfs[i], prior_analysis, location, scale, dff, alternative, ROPE)
    } else {
      t1e_TPE(t10, dfs[i], prior_design, location_d, scale_d, dff_d, alternative, ROPE)
    }

    FPE[i] <- t1e_FPE(t10, dfs[i], prior_analysis, location, scale, dff, alternative, ROPE)
    TNE[i] <- t1e_TNE(t01, dfs[i], prior_analysis, location, scale, dff, alternative, ROPE)
    FNE[i] <- if (de_an_prior == 1) {
      t1e_FNE(t01, dfs[i], prior_analysis, location, scale, dff, alternative, ROPE)
    } else {
      t1e_FNE(t01, dfs[i], prior_design, location_d, scale_d, dff_d, alternative, ROPE)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # BF01 plot
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
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

ps_N_finder <- function(threshold, true_rate, a0, b0, a1, b1, a2, b2, r,
                        prior_design_1, da1, db1, dp1, prior_design_2, da2, db2, dp2) {

  lo_n1 <- 10
  n2 <- round(lo_n1 * r)
  grid <- BF_grid_rcpp(threshold, a0, b0, a1, b1, lo_n1, a2, b2, n2,
                       prior_design_1, da1, db1, dp1, prior_design_2, da2, db2, dp2)

  pro <- sum_rcpp(grid$log_h1_dp, grid$PE)
  if (pro > true_rate) return(list(grid, lo_n1))

  # Function for uniroot
  power_fun <- function(n1){
    n1 <- round(n1)
    n2 <- n1 * r
    g <- BF_grid_rcpp(threshold, a0, b0, a1, b1, n1, a2, b2, n2,
                      prior_design_1, da1, db1, dp1, prior_design_2, da2, db2, dp2)
    sum_rcpp(g$log_h1_dp, g$PE) - true_rate - 0.005
  }

  n1 <- suppressWarnings(round(stats::uniroot(power_fun, lower = lo_n1, upper = 5000, maxiter = 10)$root))
  n2 <- round(n1 * r)
  grid <- BF_grid_rcpp(threshold, a0, b0, a1, b1, n1, a2, b2, n2,
                       prior_design_1, da1, db1, dp1, prior_design_2, da2, db2, dp2)

  list(grid, n1)
}


ps_N_01_finder<-function(threshold,true_rate, a0, b0, a1, b1, a2, b2, r,prior_design_1,da1,db1,dp1,prior_design_2,da2,db2,dp2) {

  lo_n1 <- 10
  n2 <- round(lo_n1)*r
  grid <- BF_grid_rcpp(threshold, a0, b0, a1, b1, lo_n1, a2, b2, n2,prior_design_1,da1,db1,dp1,prior_design_2,da2,db2,dp2)
  pro <- sum_rcpp(grid$log_h0,grid$NE)

  if ( pro>true_rate){
    return(list(grid,lo_n1))
  }
  TN<-function(n1){
    n1 = round(n1)
    n2 = n1*r
    grid <<- BF_grid_rcpp(threshold, a0, b0, a1, b1, n1, a2, b2, n2,prior_design_1,da1,db1,dp1,prior_design_2,da2,db2,dp2)
    pro <- sum_rcpp(grid$log_h0,grid$NE)
    return(pro - true_rate - .01)
  }
  n1 <- suppressWarnings(round(stats::uniroot(TN, lower = lo_n1, upper = 5000,maxiter = 10)$root))
  grid_power <- grid



  return(list(grid_power,n1))

}


pro_table_p2<-function(threshold,true_rate, a0, b0, a1, b1, a2, b2, r,prior_design_1,da1,db1,dp1,prior_design_2,da2,db2,dp2,mode_bf,n1,n2,type_rate) {

  if (mode_bf==1){
    x = switch(type_rate,
               "positive" = ps_N_finder(threshold,true_rate, a0, b0, a1, b1, a2, b2, r,prior_design_1,da1,db1,dp1,prior_design_2,da2,db2,dp2),
               "negative" = ps_N_01_finder(threshold,true_rate, a0, b0, a1, b1, a2, b2, r,prior_design_1,da1,db1,dp1,prior_design_2,da2,db2,dp2))
    grid = x[[1]]
    n1  = x[[2]]
    n2 =  x[[2]]*r
  }else{
    grid = BF_grid_rcpp(threshold, a0, b0, a1, b1, n1, a2, b2, n2,prior_design_1,da1,db1,dp1,prior_design_2,da2,db2,dp2)
  }
  table <- data.frame(
    TPE = sum_rcpp(grid$log_h1_dp,grid$PE),
    FNE = sum_rcpp(grid$log_h1_dp,grid$NE),
    TNE = sum_rcpp(grid$log_h0,grid$NE),
    FPE = sum_rcpp(grid$log_h0,grid$PE),
    n1  = n1,
    n2 =  n2)
  colnames(table) <- c(sprintf("p(BF10> %0.f|H1)",threshold),
                       sprintf("p(BF01> %0.f|H1)",threshold),
                       sprintf("p(BF01> %0.f|H0)",threshold),
                       sprintf("p(BF10> %0.f|H0)",threshold), " N1", "N2")
  list(table,grid)
}




p2_prior_plot <- function(a, b, ad, bd, dp, prior_analysis, nu) {

  # ---- Sequence of probabilities ----
  prop <- seq(0, 1, 0.001)

  # ---- Compute analysis prior ----
  prior_analysis_dens <- stats::dbeta(prop, a, b)

  # ---- Compute design prior ----
  prior_design_dens <- switch(prior_analysis,
                              "same"  = stats::dbeta(prop, a, b),
                              "beta"  = stats::dbeta(prop, ad, bd),
                              "Point" = rep(NA, length(prop)))

  # ---- Combine into a long-format data frame ----
  df <- data.frame(
    prop = prop,
    Density = prior_analysis_dens,
    Prior = "H1 - Analysis Prior"
  )

  if (prior_analysis != "same" && prior_analysis != "Point") {
    df_design <- data.frame(
      prop = prop,
      Density = prior_design_dens,
      Prior = "H1 - Design Prior"
    )
    df <- rbind(df, df_design)
  }

  # ---- Base ggplot ----
  p <- ggplot2::ggplot(df, ggplot2::aes(x = prop, y = Density, color = Prior, linetype = Prior)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = c("H1 - Analysis Prior" = "black",
                                           "H1 - Design Prior" = "gray")) +
    ggplot2::scale_linetype_manual(values = c("H1 - Analysis Prior" = "solid",
                                              "H1 - Design Prior" = "dashed")) +
    ggplot2::labs(
      x = substitute(bold(theta[nu_val]), list(nu_val = nu)),
      y = "density",
      title = bquote(bold("Prior distribution on "~theta[.(nu)]))
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = c(0.65, 0.95),
      legend.justification = c(0, 1),
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1.5)))

  # ---- Add vertical line for Point prior ----
  if (prior_analysis == "Point") {
    ylim_max <- max(prior_analysis_dens, na.rm = TRUE)

    # Dummy line for legend
    df_dummy <- data.frame(prop = c(NA, NA), Density = c(NA, NA),
                           Prior = "H1 - Design Prior")

    p <- p +
      ggplot2::geom_line(data = df_dummy,
                         ggplot2::aes(x = prop, y = Density,
                                      color = Prior, linetype = Prior),
                         na.rm = TRUE, show.legend = TRUE) +
      ggplot2::annotate("segment",
                        x = dp, xend = dp,
                        y = 0, yend = ylim_max,
                        colour = "gray",
                        linetype = "dashed",
                        arrow = ggplot2::arrow(length = grid::unit(0.1, "inches")))
  }

  return(p)
}
Power_p2 <- function(threshold, n1, a0, b0, a1, b1, a2, b2, r,
                     prior_design_1, da1, db1, dp1,
                     prior_design_2, da2, db2, dp2) {

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
    grid <- BF_grid_rcpp(threshold, a0, b0, a1, b1, Ns[i],
                         a2, b2, n2,
                         prior_design_1, da1, db1, dp1,
                         prior_design_2, da2, db2, dp2)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # Plot BF01
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  # Combine plots
  combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 2)

  print(combined_plot)
}


heatmap_p2 <- function(x, threshold) {
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
    "PE"   = bquote(BF[10] > .(threshold)),
    "NE"   = bquote(BF[0][1] > .(threshold)),
    "None" = bquote(1 / .(threshold) < BF[10] ~ "<" ~ .(threshold))
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

    type_rate <- switch(input$bin_type_rate,
                        "1" = "positive",
                        "0" = "negative")

    interval <- input$h0bin # point null or interval

    alternative <- switch(interval,
                          "1" =   switch(input$h1bin,
                                         "1" = "two.sided",
                                         "2" =  "greater",
                                         "3" =  "less"),
                          "2" = switch(input$h1bine,
                                       "1" = "two.sided",
                                       "2" =  "greater",
                                       "3" =  "less"))
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


    ROPE <- switch(input$h1bine,        # bound for interval test
                   "1" = c(lbbin, ubbin),
                   "2" = ubbin,
                   "3" = lbbin)

    inter <- switch(interval,
                    "1" = input$h1bin,
                    "2" = input$h1bine)


    prior_analysis <- switch(input$modelbin,
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
    prior_design <- switch(input$modelbind,
                           "1" = "beta",
                           "2" = "Moment",
                           "3" = "Point")
    location_d <- input$h0bind
    true_rate <- input$true_rate_bin
    false_rate <- input$false_rate_bin
    threshold <- input$threshold_bin
    N <- input$nbin
    Suc <- input$xbin
    pc   <- "1" %in% input$o_plot_bin
    rela <- "2" %in% input$o_plot_bin

    ############


    # Add all variables to the final list
    list(
      mode_bf = mode_bf,
      type_rate = type_rate,
      interval = interval,
      alternative =alternative,
      h0=h0,
      location = location,
      ROPE = ROPE,
      lbbin = lbbin,
      ubbin = ubbin,
      inter = inter,
      prior_analysis = prior_analysis,
      alpha = alpha,
      beta = beta,
      scale = scale,
      de_an_prior = de_an_prior,
      alpha_d = alpha_d,
      beta_d = beta_d,
      scale_d = scale_d,
      location_d = location_d,
      prior_design = prior_design,
      true_rate = true_rate,
      false_rate = false_rate,
      threshold = threshold,
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
                        \\theta_0 + \\epsilon = ', bin$location+bin$ubbin,'')

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
                              "1" = {bin_table(bin$threshold,bin$true_rate,bin$h0,bin$alpha,bin$beta,bin$location,
                                               bin$scale,bin$prior_analysis,bin$alternative,
                                               bin$alpha_d,bin$beta_d,bin$location_d,bin$scale_d,
                                               bin$prior_design,bin$de_an_prior,bin$N, bin$mode_bf,bin$false_rate,bin$type_rate)},
                              "2" = {
                                bin_e_table(bin$threshold,bin$true_rate,bin$h0,bin$alpha,bin$beta,bin$location,
                                            bin$scale,bin$prior_analysis,bin$alternative,
                                            bin$alpha_d,bin$beta_d,bin$location_d,bin$scale_d,
                                            bin$prior_design,bin$de_an_prior,bin$N, bin$mode_bf,bin$false_rate, bin$ROPE,bin$type_rate)


                              }))}, error = function(e) {
                                "Error"
                              })

    output$result_bin <-  shiny::renderText({
      paste("# Function to be used in R", show_bin_code(bin), sep = "\n")
    })
    output$prior_bin <- shiny::renderPlot({
      switch(bin$interval,
             "1" = {bin_prior_plot(bin$h0,bin$alpha,bin$beta,bin$location,bin$scale,bin$prior_analysis,
                                   bin$alpha_d,bin$beta_d,bin$location_d,
                                   bin$scale_d,bin$prior_design,bin$alternative,
                                   bin$de_an_prior)},
             "2" = bin_e_prior_plot (bin$h0,bin$alpha,bin$beta,bin$location,bin$scale,
                                     bin$prior_analysis,bin$alpha_d,bin$beta_d,bin$location_d,
                                     bin$scale_d,bin$prior_design,
                                     bin$alternative,bin$de_an_prior, bin$ROPE)
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
    p\\text{(BF}_{10} > ', bin$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{01} > ', bin$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[3], 3), nsmall = 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', bin$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[2], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{10} > ', bin$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[4], 3), nsmall = 3), ' \\\\
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

      if (identical(dat, "Error")) {

        output$plot_power_bin_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Power curve is not shown due to an error}$$")
          )
        })

        output$plot_power_bin <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_power_bin_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Power Curve}$$")
            )
          )
        })

        output$plot_power_bin <- shiny::renderPlot({

          suppressWarnings(
            switch(
              bin$interval,

              "1" = Power_bin(
                bin$threshold, bin$h0, bin$alpha, bin$beta, bin$location, bin$scale,
                bin$prior_analysis, bin$alternative,
                bin$alpha_d, bin$beta_d, bin$location_d,
                bin$scale_d, bin$prior_design, bin$de_an_prior,
                dat[1, 5]
              ),

              "2" = Power_e_bin(
                bin$threshold, bin$h0, bin$alpha, bin$beta, bin$location, bin$scale,
                bin$prior_analysis, bin$alternative,
                bin$alpha_d, bin$beta_d, bin$location_d,
                bin$scale_d, bin$prior_design, bin$de_an_prior,
                dat[1, 5],  bin$ROPE
              )
            )
          )

          pc_bin(grDevices::recordPlot())
        })
      }

    } else {

      pc_bin(NULL)
      output$plot_power_bin_text <- shiny::renderUI(NULL)
      output$plot_power_bin      <- shiny::renderPlot(NULL)
    }


    # ===================================================
    #           RELATIONSHIP BETWEEN BF & DATA (bin)
    # ===================================================

    if (isTRUE(bin$rela)) {

      if (identical(dat, "Error")) {

        output$plot_rel_bin_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Relationship plot is not shown due to an error}$$")
          )
        })

        output$plot_rel_bin <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_rel_bin_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Relationship between BF and data}$$")
            )
          )
        })

        output$plot_rel_bin <- shiny::renderPlot({

          n <- dat[1, 5]

          suppressWarnings(
            switch(
              bin$interval,

              "1" = bin_bf10(
                bin$threshold, n, bin$alpha, bin$beta, bin$location, bin$scale,
                bin$prior_analysis, bin$alternative
              ),

              "2" = bin_e_bf10(
                bin$threshold, n, bin$alpha, bin$beta, bin$location, bin$scale,
                bin$prior_analysis, bin$alternative,  bin$ROPE
              )
            )
          )

          rela_bin(grDevices::recordPlot())
        })
      }

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
                   "1" = bin_BF(bin$Suc,bin$N,bin$alpha,bin$beta,bin$location,bin$scale,bin$prior_analysis,bin$alternative),
                   "2" = bin_e_BF(bin$Suc,bin$N,bin$alpha,bin$beta,bin$location,bin$scale,bin$prior_analysis,bin$alternative, bin$ROPE))

    output$BFbin <- shiny::renderUI({
      # Create the LaTeX formatted strings for the table
      ROPE    <- switch(bin$interval,"1"=NULL,"2"= bin$ROPE)
      p.value <- bin.pval(bin$Suc,bin$N,bin$h0,bin$alternative,ROPE)
      table_html <- paste0('
    N = ', bin$N, ', x = ', bin$Suc,', \\textit{p} = ',round(p.value,4), ',\\\\ \\textit{BF}_{10} = ', round(BF10, 4), ', \\textit{BF}_{01} = ',round(1/BF10, 4),'
')

      output$result_bin <- shiny::renderText({
        args <- list(
          x = bin$Suc,
          n = bin$N,
          h0 = bin$location,
          prior_analysis = bin$prior_analysis,
          alternative = bin$alternative
        )

        # Add prior_analysis-specific parameters
        if (bin$prior_analysis == "beta") {
          args$alpha <- bin$alpha
          args$beta  <- bin$beta
        } else if (bin$prior_analysis == "Moment") {
          args$scale <- bin$scale
        }

        # Include e only if interval != 1
        if (!is.null(bin$interval) && bin$interval != 1) {
          args$ROPE <-  bin$ROPE
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

          ## prior_analysis → prior_analysis
          if (arg == "prior_analysis") arg_print <- "prior_analysis"

          ## e → ROPE
          if (arg == "e") arg_print <- "ROPE"

          ## alternative → alternative
          if (arg == "alternative") {

            arg_print <- "alternative"

            val <- switch(val,
                          "less"  = "less",
                          "two.sided" = "two.sided",
                          "greater"  = "greater",
                          stop("Invalid alternative")
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
  shiny::observeEvent(input$modelfd, {
    if (input$modelfd == 2) {  # Moment prior
      # If current df < 3, set it to 3
      if (input$dffd < 3) {
        shiny::updateSliderInput(session, "dffd", value = 3, min = 3)
      } else {
        shiny::updateSliderInput(session, "dffd", min = 3)
      }
    } else if (input$modelfd == 1) {  # Effect size prior
      shiny::updateSliderInput(session, "dffd", min = 1)
    }
  })
  shiny::observeEvent(input$modelf, {
    if (input$modelf == 2) {  # Moment prior
      # Update the slider to enforce df >= 3
      if (input$dff < 3) {
        shiny::updateSliderInput(session, "dff", value = 3, min = 3)
      } else {
        shiny::updateSliderInput(session, "dff", min = 3)
      }
    } else {  # Effect size prior
      shiny::updateSliderInput(session, "dff", min = 1)
    }
  })



  input_f <- shiny::reactive({
    mode_bf <- switch(input$Modef,
                      "1" = 1,
                      "2" = 0,
                      "3" = 0)
    type_rate<- switch(input$f_type_rate,
                       "1" = "positive",
                       "0" = "negative")
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

    ROPE <- switch(inter,
                   "1" = input$epsilinff,
                   "2" = input$epsilinff)


    prior_analysis <- switch(input$modelf,
                             "1" = "effectsize",
                             "2" = "Moment")

    rscale <- input$rf
    f_m <- sqrt(input$fsdf)
    dff <- input$dff
    de_an_prior <- switch(input$priorf,
                          "1" = 1,
                          "2" = 0)

    prior_design <- switch(input$modelfd,
                           "1" = "effectsize",
                           "2" = "Moment",
                           "3" = "Point")

    rscale_d <- input$rfd
    f_m_d <- sqrt(input$fsdfd)

    if ( input$modelfd == "3"){
      f_m_d <-sqrt(input$lfd)
    }

    dff_d <- input$dffd
    true_rate <- input$true_rate_f
    false_rate <- input$false_rate_f
    N <- input$nf
    threshold <- input$threshold_f
    fval <- input$fval
    df1 <- input$df1f
    df2 <- input$df2f
    q = k -p
    pc   <- "1" %in% input$o_plot_f
    rela <- "2" %in% input$o_plot_f

    # Add all variables to the final list


    if (input$ANOREG == 1){
      list(
        mode_bf = mode_bf,
        type_rate = type_rate,
        p = p,
        k = k,
        q = q,
        inter=inter,
        ROPE=ROPE,
        prior_analysis=prior_analysis,
        rscale=rscale,
        f_m=f_m,
        dff=dff,
        de_an_prior=de_an_prior,
        prior_design=prior_design,
        rscale_d=rscale_d,
        f_m_d=f_m_d,
        dff_d=dff_d,
        true_rate=true_rate,
        false_rate=false_rate,
        N=N,
        threshold = threshold,
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
          type_rate=type_rate,
          p = p,
          k = k,
          q = q,
          inter=inter,
          ROPE=ROPE,
          prior_analysis=prior_analysis,
          rscale=rscale,
          f_m=f_m,
          dff=dff,
          de_an_prior=de_an_prior,
          prior_design=prior_design,
          rscale_d=rscale_d,
          f_m_d=f_m_d,
          dff_d=dff_d,
          true_rate=true_rate,
          false_rate=false_rate,
          N=N,
          threshold = threshold,
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
    if (ff$prior_analysis == "effectsize"){

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
                            "1" = f_table(ff$threshold,ff$true_rate,ff$p,ff$k,ff$dff,ff$rscale,ff$f_m,ff$prior_analysis,
                                          ff$dff_d,ff$rscale_d,ff$f_m_d,ff$prior_design,ff$de_an_prior,ff$N, ff$mode_bf,ff$false_rate ,ff$type_rate),
                            "2" = fe_table(ff$threshold,ff$true_rate,ff$p,ff$k,ff$dff,ff$rscale,ff$f_m,ff$prior_analysis,
                                           ff$dff_d,ff$rscale_d,ff$f_m_d,ff$prior_design,ff$de_an_prior,ff$N, ff$mode_bf,ff$ROPE ,ff$false_rate,ff$type_rate))
    }, error = function(e) {
      "Error"
    })

    output$priorff <- shiny::renderPlot({

      switch(ff$inter,
             "1" =prior_plot_f(ff$q,ff$dff,ff$rscale,ff$f_m,ff$prior_analysis,ff$dff_d
                               ,ff$rscale_d,ff$f_m_d,ff$prior_design,ff$de_an_prior),
             "2" = prior_plot_fe(ff$q,ff$dff,ff$rscale,ff$f_m,ff$prior_analysis,ff$dff_d
                                 ,ff$rscale_d,ff$f_m_d,ff$prior_design,ff$de_an_prior,ff$ROPE))

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
    p\\text{(BF}_{10} > ', ff$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{01} > ', ff$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[3], 3), nsmall = 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', ff$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[2], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{10} > ', ff$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[4], 3), nsmall = 3), ' \\\\
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

      if (identical(dat, "Error")) {

        output$plot_power_f_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Power curve is not shown due to an error}$$")
          )
        })

        output$plot_power_f <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_power_f_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Power Curve}$$")
            )
          )
        })

        output$plot_power_f <- shiny::renderPlot({

          suppressWarnings(
            switch(
              ff$inter,

              "1" = Power_f(
                ff$threshold, ff$k, ff$p, ff$dff, ff$rscale,
                ff$f_m, ff$prior_analysis,
                ff$k_d, ff$p_d, ff$dff_d, ff$rscale_d, ff$f_m_d, ff$prior_design,
                ff$de_an_prior,
                dat[1, 5]
              ),

              "2" = Power_fe(
                ff$threshold, ff$k, ff$p, ff$dff, ff$rscale,
                ff$f_m, ff$prior_analysis,
                ff$k_d, ff$p_d, ff$dff_d, ff$rscale_d, ff$f_m_d, ff$prior_design,
                ff$de_an_prior,
                dat[1, 5],
                ff$ROPE
              )
            )
          )

          pc_f(grDevices::recordPlot())
        })
      }

    } else {

      pc_f(NULL)
      output$plot_power_f_text <- shiny::renderUI(NULL)
      output$plot_power_f      <- shiny::renderPlot(NULL)
    }


    # ===================================================
    #           RELATIONSHIP BETWEEN BF & DATA (F)
    # ===================================================
    if (isTRUE(ff$rela)) {

      if (identical(dat, "Error")) {

        output$plot_rel_f_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Relationship plot is not shown due to an error}$$")
          )
        })

        output$plot_rel_f <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_rel_f_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Relationship between BF and data}$$")
            )
          )
        })

        output$plot_rel_f <- shiny::renderPlot({

          n <- dat[1, 5]

          suppressWarnings(
            switch(
              ff$inter,

              "1" = bf10_f(
                ff$threshold, n, ff$k, ff$p, ff$dff, ff$rscale, ff$f_m, ff$prior_analysis
              ),

              "2" = bf10_fe(
                ff$threshold, n, ff$k, ff$p, ff$dff, ff$rscale, ff$f_m, ff$prior_analysis, ff$ROPE
              )
            )
          )

          rela_f(grDevices::recordPlot())
        })
      }

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
             "1" =prior_plot_f(ff$q,ff$dff,ff$rscale,ff$f_m,ff$prior_analysis,ff$dff_d
                               ,ff$rscale_d,ff$f_m_d,ff$prior_design,1),
             "2" = prior_plot_fe(ff$q,ff$dff,ff$rscale,ff$f_m,ff$prior_analysis,ff$dff_d
                                 ,ff$rscale_d,ff$f_m_d,ff$prior_design,1,ff$ROPE))

    })
    BF10 <- if (ff$inter==1) F_BF(ff$fval,ff$df1,m,ff$dff,ff$rscale,ff$f_m,ff$prior_analysis) else
      Fe_BF(ff$fval,ff$df1,m,ff$dff,ff$rscale,ff$f_m,ff$prior_analysis,ff$ROPE)

    output$result_f <- shiny::renderText({

      args <- list(
        fval = ff$fval,
        df1 = ff$df1,
        df2 = ff$df2,
        dff = ff$dff,
        rscale = ff$rscale,
        f_m = ff$f_m,
        prior_analysis = ff$prior_analysis
      )

      if (ff$inter!=1) {
        args$e <- ff$ROPE
      }

      # Build string with each argument on a new line
      arg_strings <- sapply(names(args), function(arg) {

        val <- args[[arg]]
        arg_print <- arg

        ## prior_analysis > prior_analysis
        if (arg == "prior_analysis") {
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
      ROPE <- switch(ff$inter,"1" = NULL,"2" = ff$ROPE)
      p.value <- f.pval(ff$fval, ff$df1,ff$df2,ROPE=ROPE)

      # Create the LaTeX formatted strings for the table
      output$BFcalf <- shiny::renderUI({
        # Create the LaTeX formatted string with proper escaping
        table_latex <- paste0(
          "$$ \\textit{F}(", ff$df1, ",", ff$df2,
          ") = ", round(ff$fval, 3),", \\textit{p} = ",round(p.value,4) ,
          ",\\\\ \\textit{BF}_{10} = ", round(BF10, 4),", \\textit{BF}_{01} =" ,round(1/BF10, 4)," $$"
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
    type_rate <- switch(input$p2_type_rate,
                        "1" = "positive",
                        "0" = "negative")
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


    prior_design_1 <-switch(input$model_p1,
                            "1" = "Point",
                            "2" = "beta")
    prior_design_2 <-switch(input$model_p2,
                            "1" = "Point",
                            "2" = "beta")

    if (input$priorp2 == 1){
      prior_design_1 = prior_design_2 = "same"
    }
    de_an_prior<-input$priorp2
    threshold <- input$threshold_p2

    N1 <- input$n1p2
    N2 <- input$n2p2

    x1 <- input$x1p2
    x2 <- input$x2p2

    true_rate <- input$true_rate_p2
    pc   <- "1" %in% input$o_plot_p2
    rela <- "2" %in% input$o_plot_p2

    ############

    list(
      mode_bf = mode_bf,
      type_rate = type_rate,
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
      prior_design_1 = prior_design_1,
      prior_design_2 = prior_design_2,
      dp1 = dp1,
      dp2 = dp2,
      threshold = threshold,
      N1 = round(N1),
      N2 = round(N2),
      k1 = round(x1),
      k2 = round(x2),
      true_rate = true_rate,
      r=1,
      pc=pc,
      rela=rela,
      de_an_prior=de_an_prior

    )
  })



  shiny::observeEvent(input$runp2, {

    p2 <- input_p2()

    # Compute data safely
    dat <-  tryCatch({pro_table_p2(
      p2$threshold, p2$true_rate, p2$a0, p2$b0,
      p2$a1, p2$b1, p2$a2, p2$b2, p2$r,
      p2$prior_design_1, p2$a1d, p2$b1d, p2$dp1,
      p2$prior_design_2, p2$a2d, p2$b2d, p2$dp2,
      p2$mode_bf, p2$N1, p2$N2, p2$type_rate
    )}, error = function(e) {
      "Error"
    })


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
      p2_prior_plot(p2$a1, p2$b1, p2$a1d, p2$b1d, p2$dp1, p2$prior_design_1, 1)
    })
    output$prior_p2 <- shiny::renderPlot({
      p2_prior_plot(p2$a2, p2$b2, p2$a2d, p2$b2d, p2$dp2, p2$prior_design_2, 2)
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
      p\\text{(BF}_{10} > ', p2$threshold, '\\, | \\, \\mathcal{H}_1) & ', format(round(table[1,1],3), nsmall = 3), ' \\\\
      p\\text{(BF}_{01} > ', p2$threshold, '\\, | \\, \\mathcal{H}_0) & ', format(round(table[1,3],3), nsmall = 3), ' \\\\
      \\textbf{Probability of Misleading Evidence} & \\\\
      \\hline
      p\\text{(BF}_{01} > ', p2$threshold, '\\, | \\, \\mathcal{H}_1) & ', format(round(table[1,2],3), nsmall = 3), ' \\\\
      p\\text{(BF}_{10} > ', p2$threshold, '\\, | \\, \\mathcal{H}_0) & ', format(round(table[1,4],3), nsmall = 3), ' \\\\
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

      if (identical(table, "Error")) {

        output$plot_power_p2_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Power curve is not shown due to an error}$$")
          )
        })

        output$plot_power_p2 <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_power_p2_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Power Curve}$$")
            )
          )
        })

        output$plot_power_p2 <- shiny::renderPlot({

          suppressWarnings(
            Power_p2(
              p2$threshold, table[1, 5], p2$a0, p2$b0, p2$a1, p2$b1, p2$a2,
              p2$b2, table[1, 6] / table[1, 5], p2$prior_design_1, p2$a1d, p2$b1d, p2$dp1,
              p2$prior_design_2, p2$a2d, p2$b2d, p2$dp2
            )
          )

          pc_p2(grDevices::recordPlot())
        })
      }

    } else {

      pc_p2(NULL)
      output$plot_power_p2_text <- shiny::renderUI(NULL)
      output$plot_power_p2      <- shiny::renderPlot(NULL)
    }

    # ===================================================
    #           RELATIONSHIP BETWEEN BF & DATA
    # ===================================================

    if (isTRUE(p2$rela)) {

      if (identical(dat, "Error")) {

        output$plot_rel_p2_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Relationship plot is not shown due to an error}$$")
          )
        })

        output$plot_rel_p2 <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_rel_p2_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Relationship between BF and data}$$")
            )
          )
        })

        output$plot_rel_p2 <- shiny::renderPlot({

          # Explicitly assign and print the ggplot
          plt <- heatmap_p2(dat[[2]], p2$threshold)
          print(plt)

          # Save the plot for later
          rela_p2(grDevices::recordPlot())
        })
      }

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
    BF10 <- BF10_p2(p2$a0, p2$b0, p2$a1, p2$b1, p2$a2, p2$b2,p2$N1,p2$N2,p2$k1,p2$k2)
    tab <- matrix(
      c(p2$k1, p2$N1 - p2$k1,
        p2$k2, p2$N2 - p2$k2),
      nrow = 2,
      byrow = TRUE
    )


    # Render priors
    output$prior_p0 <- shiny::renderPlot({
      p2_prior_plot(p2$a0, p2$b0, 1, 1, 0, "same", 0)
    })
    output$prior_p1 <- shiny::renderPlot({
      p2_prior_plot(p2$a1, p2$b1, p2$a1d, p2$b1d, p2$dp1, p2$prior_design_1, 1)
    })
    output$prior_p2 <- shiny::renderPlot({
      p2_prior_plot(p2$a2, p2$b2, p2$a2d, p2$b2d, p2$dp2, p2$prior_design_2, 2)
    })


    results <-stats::fisher.test(tab)
    output$BFp2 <- shiny::renderUI({
      # Create the LaTeX formatted strings for the table
      table_html <- paste0(
        'n_1 = ', p2$N1, ', ',
        'n_2 = ', p2$N2, ', ',
        'x_1 = ', p2$k1, ', ',
        'x_2 = ', p2$k2, ' \\\\ ',
        '\\textit{Odd Ratio = }',round(results$estimate,4),', \\textit{p} = ',round(results$p.value,4),' \\\\ ',
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
          N1 = p2$N1,
          N2 = p2$N2,
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
    type_rate <- switch(input$r_type_rate,
                        "1" = "positive",
                        "0" = "negative")
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



    ROPE <- switch(input$h1re,        # bound for interval test
                   "1" = c(lbre, ubre),
                   "2" = ubre,
                   "3" = lbre)
    inter <- switch(interval,
                    "1" = input$h1r,
                    "2" = input$h1re)


    alternative <- switch(interval,
                          "1" =   switch(input$h1r,
                                         "1" = "two.sided",
                                         "2" =  "greater",
                                         "3" =  "less"),
                          "2" = switch(input$h1re,
                                       "1" = "two.sided",
                                       "2" =  "greater",
                                       "3" =  "less"))





    prior_analysis <- switch(input$modelr,
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
    prior_design <- switch(input$modelrd,
                           "1" = "d_beta",
                           "2" = "beta",
                           "3" = "Moment",
                           "4" = "Point")
    location_d <- input$h0phod
    k_d <- input$rkd
    scale_d <- input$rsd
    alpha_d <- input$ralphad
    beta_d<- input$rbetad
    true_rate <- input$true_rate_r
    false_rate <- input$false_rate_r
    threshold <- input$threshold_r
    N <-  switch(input$Moder,
                 "1" = 2,
                 "2" = input$nr,
                 "3" = input$rdf)
    rval <- input$rval
    pc   <- 1 %in% input$o_plot_r
    rela <- 2 %in% input$o_plot_r

    ###########
    location <- h0
    dff <- 1

    dff_d <- 1



    ############


    # Add all variables to the final list
    list(
      mode_bf = mode_bf,
      type_rate = type_rate,
      interval = interval,
      ROPE = ROPE,
      lbre = lbre,
      ubre = ubre,
      inter = inter,
      alternative = alternative,
      h0 = h0,
      prior_analysis = prior_analysis,
      k =k,
      scale = scale,
      alpha = alpha,
      beta = beta,
      de_an_prior = de_an_prior,
      prior_design =prior_design,
      location_d = location_d,
      k_d = k_d,
      scale_d = scale_d,
      alpha_d = alpha_d,
      beta_d = beta_d,
      true_rate = true_rate,
      false_rate = false_rate,
      threshold = threshold,
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

             "1" = r_table(rr$threshold,rr$true_rate,rr$prior_analysis,rr$k,
                           rr$alpha, rr$beta,rr$h0,rr$location,
                           rr$scale,rr$dff, rr$alternative ,rr$prior_design,
                           rr$location_d,rr$k_d, rr$alpha_d, rr$beta_d,
                           rr$scale_d,rr$dff_d,rr$de_an_prior,rr$N,
                           rr$mode_bf,rr$false_rate,rr$type_rate ),
             "2" = re_table(rr$threshold,rr$true_rate,rr$prior_analysis,rr$k,
                            rr$alpha, rr$beta,rr$h0,rr$location,
                            rr$scale,rr$dff, rr$alternative ,rr$prior_design,
                            rr$location_d,rr$k_d, rr$alpha_d, rr$beta_d,
                            rr$scale_d,rr$dff_d,rr$de_an_prior,rr$N,
                            rr$mode_bf,rr$false_rate,rr$ROPE,rr$type_rate ))
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
                                rr$dff,rr$prior_analysis,rr$de_an_prior,
                                rr$k_d, rr$alpha_d, rr$beta_d,
                                rr$location_d,rr$scale_d,rr$dff_d,
                                rr$prior_design,rr$alternative),
             "2" = re_prior_plot(rr$k, rr$alpha, rr$beta,
                                 rr$h0,rr$location,rr$scale,
                                 rr$dff,rr$prior_analysis,rr$de_an_prior,
                                 rr$k_d, rr$alpha_d, rr$beta_d,
                                 rr$location_d,rr$scale_d,rr$dff_d,
                                 rr$prior_design,rr$alternative,rr$ROPE))



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
    p\\text{(BF}_{10} > ', rr$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{01} > ', rr$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[3], 3), nsmall = 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', rr$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[2], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{10} > ', rr$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[4], 3), nsmall = 3), ' \\\\
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

      if (identical(dat, "Error")) {

        output$plot_power_r_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Power curve is not shown due to an error}$$")
          )
        })

        output$plot_power_r <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_power_r_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Power Curve}$$")
            )
          )
        })

        output$plot_power_r <- shiny::renderPlot({

          suppressWarnings(
            switch(
              rr$interval,

              "1" = Power_r(
                rr$threshold, rr$k, rr$alpha, rr$beta, rr$h0, rr$alternative,
                rr$location, rr$scale, rr$dff, rr$prior_analysis,
                rr$k_d, rr$alpha_d, rr$beta_d, rr$location_d, rr$scale_d,
                rr$dff_d, rr$prior_design, rr$de_an_prior, dat[1, 5]
              ),

              "2" = Power_re(
                rr$threshold, rr$k, rr$alpha, rr$beta, rr$h0, rr$alternative,
                rr$location, rr$scale, rr$dff, rr$prior_analysis,
                rr$k_d, rr$alpha_d, rr$beta_d, rr$location_d, rr$scale_d,
                rr$dff_d, rr$prior_design, rr$de_an_prior, dat[1, 5], rr$ROPE
              )
            )
          )

          pc_r(grDevices::recordPlot())
        })
      }

    } else {

      pc_r(NULL)
      output$plot_power_r_text <- shiny::renderUI(NULL)
      output$plot_power_r      <- shiny::renderPlot(NULL)
    }




    # ===================================================
    #                RELATIONSHIP (r)
    # ===================================================

    if (isTRUE(rr$rela)) {

      if (identical(dat, "Error")) {

        output$plot_rel_r_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Relationship plot is not shown due to an error}$$")
          )
        })

        output$plot_rel_r <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_rel_r_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Relationship between BF and data}$$")
            )
          )
        })

        output$plot_rel_r <- shiny::renderPlot({

          n <- dat[1, 5]

          suppressWarnings(
            switch(
              rr$interval,

              "1" = r_bf10_p(
                rr$threshold, n, rr$k, rr$alpha, rr$beta, rr$h0,
                rr$alternative, rr$location, rr$scale, rr$dff, rr$prior_analysis
              ),

              "2" = re_bf10_p(
                rr$threshold, n, rr$k, rr$alpha, rr$beta, rr$h0,
                rr$alternative, rr$location, rr$scale, rr$dff, rr$prior_analysis, rr$ROPE
              )
            )
          )

          rela_r(grDevices::recordPlot())
        })
      }

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
                                rr$dff,rr$prior_analysis,1,
                                rr$k_d, rr$alpha_d, rr$beta_d,
                                rr$location_d,rr$scale_d,rr$dff_d,
                                rr$prior_design,rr$alternative),
             "2" = re_prior_plot(rr$k, rr$alpha, rr$beta,
                                 rr$h0,rr$location,rr$scale,
                                 rr$dff,rr$prior_analysis,1,
                                 rr$k_d, rr$alpha_d, rr$beta_d,
                                 rr$location_d,rr$scale_d,rr$dff_d,
                                 rr$prior_design,rr$alternative,rr$ROPE))



    })
    BF10 <- switch(rr$interval ,
                   "1" = r_BF10(rr$rval,rr$N,rr$k, rr$alpha, rr$beta,rr$h0,rr$alternative,rr$location,rr$scale,rr$dff,rr$prior_analysis),
                   "2" = re_BF10(rr$rval,rr$N,rr$k, rr$alpha, rr$beta,rr$h0,rr$alternative,rr$location,rr$scale,rr$dff,rr$prior_analysis,rr$ROPE))

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
          prior_analysis = rr$prior_analysis
        )

        # Add prior_analysis-specific arguments
        if (!is.null(rr$prior_analysis)) {
          if (rr$prior_analysis == "d_beta") {
            args$k <- rr$k
          } else if (rr$prior_analysis == "beta") {
            args$alpha <- rr$alpha
            args$beta  <- rr$beta
          } else if (rr$prior_analysis == "Moment") {
            args$scale <- rr$scale
          }
        }

        # Common arguments
        args$h0 <- rr$h0
        args$alternative <- rr$alternative

        # Optional ROPE (only if interval != 1)
        if (!is.null(rr$interval) && rr$interval != 1) {
          args$e <- rr$ROPE
        }

        # Build string with renaming & mapping rules
        arg_strings <- sapply(names(args), function(arg) {

          val <- args[[arg]]
          arg_print <- arg

          ## prior_analysis > prior_analysis
          if (arg == "prior_analysis") {
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

          if (arg == "alternative") {
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

      ROPE <- switch(rr$interval,"1" = NULL,"2" = rr$ROPE)
      p.value <- r.pval(rr$rval, rr$N,rr$h0, rr$alternative , ROPE = ROPE)

      # Create the LaTeX formatted strings for the table
      table_html <- paste0('
    \\textit{r}(n = ', rr$N , ') = ',rr$rval,', \\textit{p} = ',round(p.value,4),', \\\\ \\textit{BF}_{10} = ', round(BF10, 4),", \\textit{BF}_{01} = ",round(1/BF10, 4), '
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
    type_rate <- switch(input$t1_type_rate,
                        "1" = "positive",
                        "0" = "negative")

    interval <- input$h0t1 # point null or interval

    ROPE <- switch(input$h1t1e,        # bound for interval test
                   "1" = c(input$lbt1e, input$ubt1e),
                   "2" = input$ubt1e,
                   "3" = input$lbt1e)
    inter <- switch(interval,
                    "1" = input$h1t1,
                    "2" = input$h1t1e)

    alternative <- switch(interval,
                          "1" =   switch(input$h1t1,
                                         "1" = "two.sided",
                                         "2" =  "greater",
                                         "3" =  "less"),
                          "2" = switch(input$h1t1e,
                                       "1" = "two.sided",
                                       "2" =  "greater",
                                       "3" =  "less"))

    prior_analysis <- switch(input$modelt1,
                             "1" = "t-distribution",
                             "2" = "Normal",
                             "3" = "Moment")

    location <- input$lt1
    scale <- input$st1
    dff <- input$dft1
    de_an_prior <- switch(input$prior,
                          "1" = 1,
                          "2" = 0)
    prior_design <- switch(input$modelt1d,
                           "1" = "t-distribution",
                           "2" = "Normal",
                           "3" = "Moment",
                           "4" = "Point")
    location_d <- input$lt1d

    scale_d <- input$st1d
    dff_d <- input$dft1d
    threshold <- input$threshold_t1
    type <- input$typet1
    true_rate <- input$true_rate_t1
    false_rate <- input$false_rate_t1

    tval <- input$t1tval
    pc   <- "1" %in% input$o_plot_t1
    rela <- "2" %in% input$o_plot_t1


    # Add all variables to the final list
    list(
      mode_bf = mode_bf,
      type_rate = type_rate,
      interval = interval,
      alternative = alternative ,
      ROPE = ROPE,
      prior_analysis = prior_analysis,
      location = location,
      scale = scale,
      dff = dff,
      de_an_prior = de_an_prior,
      prior_design = prior_design,
      location_d = location_d,
      scale_d = scale_d,
      dff_d = dff_d,
      type = type,
      threshold = threshold,
      true_rate = true_rate,
      false_rate = false_rate ,
      N = N,
      tval = tval,
      pc = pc,
      rela = rela
    )
  })

  shiny::observeEvent(input$runt1, {
    x = input_t1()

    dat = tryCatch({suppressWarnings(switch(x$interval, "1" =  t1_Table(x$threshold,x$true_rate,x$prior_analysis,x$location,x$scale,x$dff, x$alternative,
                                                                        x$prior_design,x$location_d,x$scale_d,x$dff_d, x$de_an_prior,x$N, x$mode_bf ,
                                                                        x$false_rate,x$type_rate),
                                            "2" = t1e_table(x$threshold,x$true_rate,x$prior_analysis,x$location,x$scale,x$dff, x$alternative,x$ROPE ,
                                                            x$prior_design,x$scale_d,x$dff_d, x$de_an_prior,x$N,x$mode_bf,x$location_d ,x$false_rate,x$type_rate)))
    }, error = function(e) {
      "Error"
    })

    output$result_t1 <- shiny::renderText({
      paste("# Function to be used in R", show_t1_code(x), sep = "\n")
    })
    output$priort1 <- shiny::renderPlot({
      suppressWarnings(switch(x$interval,
                              "1"= t1_prior_plot(        # Access 'target' explicitly
                                prior_analysis = x$prior_analysis,          # Access 'prior_analysis' explicitly
                                location = x$location,    # Access 'location' explicitly
                                scale = x$scale,          # Access 'scale' explicitly
                                dff = x$dff,              # Access 'dff' explicitly
                                alternative = x$alternative,  # Access 'alternative' explicitly
                                prior_design = x$prior_design,        # Access 'prior_design' explicitly
                                location_d = x$location_d,  # Access 'location_d' explicitly
                                scale_d = x$scale_d,        # Access 'scale_d' explicitly
                                dff_d = x$dff_d,            # Access 'dff_d' explicitly
                                de_an_prior = x$de_an_prior   # Access 'de_an_prior' explicitly
                              ),
                              "2" = t1e_prior_plot(x$prior_analysis,
                                                   x$location,
                                                   x$scale,
                                                   x$dff ,
                                                   x$alternative,
                                                   x$ROPE,
                                                   x$de_an_prior,
                                                   x$prior_design,
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
    p\\text{(BF}_{10} > ', x$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{01} > ', x$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[3], 3), nsmall = 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', x$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[2], 3), nsmall = 3), ' \\\\
    p\\text{(BF}_{10} > ', x$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[4], 3), nsmall = 3), ' \\\\
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
      if (identical(dat, "Error")) {

        output$plot_power_t1_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Power curve is not shown due to an error}$$")
          )
        })

        output$plot_power_t1 <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_power_t1_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Power Curve}$$")
            )
          )
        })

        output$plot_power_t1 <- shiny::renderPlot({

          plt <- suppressWarnings(
            switch(
              x$interval,
              "1" = Power_t1(
                x$threshold, x$prior_analysis, x$location, x$scale, x$dff, x$alternative,
                x$prior_design, x$location_d, x$scale_d, x$dff_d,
                x$de_an_prior, dat[1,5]
              ),
              "2" = Power_t1e(
                x$threshold, x$prior_analysis, x$location, x$scale, x$dff, x$alternative,
                x$prior_design, x$location_d, x$scale_d, x$dff_d,
                x$de_an_prior, dat[1,5], x$ROPE
              )
            )
          )

          print(plt)
          pc_t1(grDevices::recordPlot())
        })
      }

    } else {

      pc_t1(NULL)
      output$plot_power_t1_text <- shiny::renderUI(NULL)
      output$plot_power_t1      <- shiny::renderPlot(NULL)
    }



    # ===============================
    #   RELATIONSHIP SECTION (rela)
    # ===============================

    if (isTRUE(x$rela)) {

      if (identical(dat, "Error")) {

        output$plot_rel_t1_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Relationship plot is not shown due to an error}$$")
          )
        })

        output$plot_rel_t1 <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_rel_t1_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Relationship between BF and data}$$")
            )
          )
        })

        output$plot_rel_t1 <- shiny::renderPlot({

          plt <- suppressWarnings(
            switch(
              x$interval,
              "1" =
                bf10_t1(
                  threshold          = x$threshold,
                  df         = dat[1,5],
                  prior_analysis      = x$prior_analysis,
                  location   = x$location,
                  scale      = x$scale,
                  dff        = x$dff,
                  alternative = x$alternative
                ),
              "2" =
                te1_BF(
                  x$threshold, dat[1,5], x$prior_analysis, x$location, x$scale, x$dff,
                  x$alternative, x$ROPE
                )
            )
          )

          rela_t1(grDevices::recordPlot())
        })
      }

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
      prior_analysis = x$prior_analysis,
      location = x$location,
      scale = x$scale,
      dff = x$dff,
      alternative = x$alternative
    )

    if (!is.null(x$ROPE) && x$interval != 1) {
      args$ROPE <- x$ROPE
    }

    # Build string with each argument on a new line
    arg_strings <- sapply(names(args), function(arg) {

      val <- args[[arg]]
      arg_print <- arg

      # ---- Omission rules ----
      if ((x$prior_analysis %in% c("Normal", "Moment")) && arg == "dff") {
        return(NULL)
      }

      # Handle alternative
      if (arg == "alternative") {
        val <- shQuote(val)
      } else if (arg == "ROPE") {
        # Wrap vector ROPE in c(...)
        if (length(val) > 1) {
          val <- paste0("c(", paste(val, collapse = ", "), ")")
        } else {
          val <- as.character(val)
        }
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
                              "1"= t1_prior_plot(       # Access 'target' explicitly
                                prior_analysis = x$prior_analysis,          # Access 'prior_analysis' explicitly
                                location = x$location,    # Access 'location' explicitly
                                scale = x$scale,          # Access 'scale' explicitly
                                dff = x$dff,              # Access 'dff' explicitly
                                alternative = x$alternative,  # Access 'alternative' explicitly
                                prior_design = x$prior_design,        # Access 'prior_design' explicitly
                                location_d = x$location_d,  # Access 'location_d' explicitly
                                scale_d = x$scale_d,        # Access 'scale_d' explicitly
                                dff_d = x$dff_d,            # Access 'dff_d' explicitly
                                de_an_prior = 1   # Access 'de_an_prior' explicitly
                              ),
                              "2" = t1e_prior_plot(x$prior_analysis,
                                                   x$location,
                                                   x$scale,
                                                   x$dff ,
                                                   x$alternative,
                                                   x$ROPE,
                                                   1,
                                                   x$prior_design,
                                                   x$scale_d,
                                                   x$dff_d,
                                                   x$location_d )

      ))

    })
    BF10 <- suppressWarnings(switch(x$interval,
                                    "1" = t1_BF10(x$tval,x$N,x$prior_analysis ,x$location,x$scale,x$dff , x$alternative ),
                                    "2" = t1e_BF10(x$tval,x$N,x$prior_analysis,x$location,x$scale,x$dff , x$alternative,x$ROPE )))
    d.obs <- x$tval/sqrt(x$N)
    ROPE <- switch(x$interval,"1" = NULL,"2" = x$ROPE)
    p.value <- t.pval(x$tval, x$N+1, n2 = NULL, x$alternative, ROPE = ROPE, type = "One-sample t-test")
    output$BFt1 <- shiny::renderUI({
      # Create the LaTeX formatted strings for the table
      table_html <- paste0('
    \\textit{t}(', x$N , ') = ',x$tval,', \\textit{p} = ',round(p.value,4),', \\textit{d} = ',round(d.obs,4),',\\\\ \\textit{BF}_{10} = ', round(BF10, 4),", \\textit{BF}_{01} = ",round(1/BF10, 4), '
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
    type_rate <- switch(input$t2_type_rate,
                        "1" = "positive",
                        "0" = "negative")
    interval <- input$h0t2 # point null or interval

    ROPE <- switch(input$h1t2e,        # bound for interval test
                   "1" = c(input$lbt2e, input$ubt2e),
                   "2" = input$ubt2e,
                   "3" = input$lbt2e)
    inter <- switch(interval,
                    "1" = input$h1t2,
                    "2" = input$h1t2e)

    alternative <- switch(interval,
                          "1" =   switch(input$h1t2,
                                         "1" = "two.sided",
                                         "2" =  "greater",
                                         "3" =  "less"),
                          "2" = switch(input$h1t2e,
                                       "1" = "two.sided",
                                       "2" =  "greater",
                                       "3" =  "less"))



    prior_analysis <- switch(input$modelt2,
                             "1" = "t-distribution",
                             "2" = "Normal",
                             "3" = "Moment")

    location <- input$lt2
    scale <- input$st2
    dff <- input$dft2
    de_an_prior <- switch(input$priort2,
                          "1" = 1,
                          "2" = 0)
    prior_design <- switch(input$modelt2d,
                           "1" = "t-distribution",
                           "2" = "Normal",
                           "3" = "Moment",
                           "4" = "Point")
    location_d <- input$lt2d
    scale_d <- input$st2d
    dff_d <- input$dft2d
    threshold <- input$threshold_t2
    type <- input$typet2
    true_rate <- input$true_rate_t2
    false_rate <- input$false_rate_t2
    tval <- input$t2tval
    r <- switch(input$Modet2,
                "1" = input$rt2,
                "2" = input$n2t2/input$n1t2,
                "3" = input$rt2)
    N1 = input$n1t2
    N2 = input$n2t2
    pc   <- "1" %in% input$o_plot_t2
    rela <- "2" %in% input$o_plot_t2

    # Add all variables to the final list
    list(
      mode_bf = mode_bf,
      type_rate = type_rate,
      interval = interval,
      alternative = alternative,
      ROPE = ROPE,
      prior_analysis = prior_analysis,
      location = location,
      scale = scale,
      dff = dff,
      de_an_prior = de_an_prior,
      prior_design = prior_design,
      location_d = location_d,
      scale_d = scale_d,
      dff_d = dff_d,
      type = type,
      threshold = threshold,
      true_rate = true_rate,
      false_rate = false_rate,
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
                              "1" = t2_Table(t2$threshold, t2$r, t2$true_rate, t2$prior_analysis, t2$location, t2$scale, t2$dff, t2$alternative,
                                             t2$prior_design, t2$location_d, t2$scale_d, t2$dff_d, t2$de_an_prior, t2$N1, t2$N2, t2$mode_bf, t2$false_rate,t2$type_rate),
                              "2" = t2e_table(t2$threshold, t2$r, t2$true_rate, t2$prior_analysis,t2$location, t2$scale, t2$dff, t2$alternative, t2$ROPE,
                                              t2$prior_design,t2$location_d, t2$scale_d, t2$dff_d, t2$de_an_prior, t2$mode_bf,  t2$N1, t2$N2, t2$false_rate,t2$type_rate)
      ))
    }, error = function(e) {
      "Error"
    })
    output$result_t2 <- shiny::renderText({ paste("# Function to be used in R", show_t2_code(t2), sep = "\n") })

    output$priort2 <- shiny::renderPlot({
      suppressWarnings(switch(t2$interval,
                              "1"=
                                t1_prior_plot(       # Access 'true_rate' explicitly
                                  prior_analysis = t2$prior_analysis,          # Access 'prior_analysis' explicitly
                                  location = t2$location,    # Access 'location' explicitly
                                  scale = t2$scale,          # Access 'scale' explicitly
                                  dff = t2$dff,              # Access 'dff' explicitly
                                  alternative = t2$alternative,  # Access 'alternative' explicitly
                                  prior_design = t2$prior_design,        # Access 'prior_design' explicitly
                                  location_d = t2$location_d,  # Access 'location_d' explicitly
                                  scale_d = t2$scale_d,        # Access 'scale_d' explicitly
                                  dff_d = t2$dff_d,            # Access 'dff_d' explicitly
                                  de_an_prior = t2$de_an_prior   # Access 'de_an_prior' explicitly
                                ), "2" =
                                t1e_prior_plot(t2$prior_analysis,
                                               t2$location,
                                               t2$scale,
                                               t2$dff ,
                                               t2$alternative,
                                               t2$ROPE,
                                               t2$de_an_prior,
                                               t2$prior_design,
                                               t2$scale_d,
                                               t2$dff_d,
                                               t2$location_d )

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
    p\\text{(BF}_{10} > ', t2$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1,1], 3),nsmall=3), ' \\\\
    p\\text{(BF}_{01} > ', t2$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[1,3], 3),nsmall=3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    p\\text{(BF}_{01} > ', t2$threshold, '\\, | \\, \\mathcal{H}_1)\\ & ', format(round(dat[1,2], 3),nsmall=3), ' \\\\
    p\\text{(BF}_{10} > ', t2$threshold, '\\, | \\, \\mathcal{H}_0)\\ & ', format(round(dat[1,4], 3),nsmall=3), ' \\\\
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

      if (identical(dat, "Error")) {

        output$plot_power_t2_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Power curve is not shown due to an error}$$")
          )
        })

        output$plot_power_t2 <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_power_t2_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Power Curve}$$")
            )
          )
        })

        output$plot_power_t2 <- shiny::renderPlot({

          plt<-suppressWarnings(
            switch(
              t2$interval,

              "1" = Power_t2(
                t2$threshold, t2$prior_analysis, t2$location, t2$scale, t2$dff, t2$alternative,
                t2$prior_design, t2$location_d, t2$scale_d, t2$dff_d,
                t2$de_an_prior,
                unlist(dat[1,5]),
                unlist(dat[1,6]) / unlist(dat[1,5])
              ),

              "2" = Power_t2e(
                t2$threshold, t2$prior_analysis, t2$location, t2$scale, t2$dff, t2$alternative,
                t2$prior_design, t2$location_d, t2$scale_d, t2$dff_d,
                t2$de_an_prior,
                dat[1,5],
                dat[1,6] / dat[1,5],
                t2$ROPE
              )
            )
          )
          print(plt)
          pc_t2(grDevices::recordPlot())
        })
      }

    } else {

      pc_t2(NULL)
      output$plot_power_t2_text <- shiny::renderUI(NULL)
      output$plot_power_t2      <- shiny::renderPlot(NULL)
    }


    # ===================================================
    #                RELATIONSHIP (t2)
    # ===================================================

    if (isTRUE(t2$rela)) {

      if (identical(dat, "Error")) {

        output$plot_rel_t2_text <- shiny::renderUI({
          shiny::withMathJax(
            shiny::em("$$\\text{Relationship plot is not shown due to an error}$$")
          )
        })

        output$plot_rel_t2 <- shiny::renderPlot({
          NULL
        })

      } else {

        output$plot_rel_t2_text <- shiny::renderUI({
          shiny::tagList(
            shiny::withMathJax(
              shiny::em("$$\\text{Relationship between BF and data}$$")
            )
          )
        })

        output$plot_rel_t2 <- shiny::renderPlot({

          plt<-suppressWarnings(
            switch(
              t2$interval,

              "1" = t2_BF(
                t2$threshold, dat[1, 5], t2$r, t2$true_rate,
                t2$prior_analysis, t2$location, t2$scale, t2$dff,
                t2$alternative
              ),

              "2" = t2e_BF(
                t2$threshold, dat[1, 5], t2$r,
                t2$prior_analysis, t2$location, t2$scale, t2$dff,
                t2$alternative, t2$ROPE
              )
            )
          )
          print(plt)
          rela_t2(grDevices::recordPlot())
        })
      }

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
                                t1_prior_plot(        # Access 'true_rate' explicitly
                                  prior_analysis = t2$prior_analysis,          # Access 'prior_analysis' explicitly
                                  location = t2$location,    # Access 'location' explicitly
                                  scale = t2$scale,          # Access 'scale' explicitly
                                  dff = t2$dff,              # Access 'dff' explicitly
                                  alternative = t2$alternative,  # Access 'alternative' explicitly
                                  prior_design = t2$prior_design,        # Access 'prior_design' explicitly
                                  location_d = t2$location_d,  # Access 'location_d' explicitly
                                  scale_d = t2$scale_d,        # Access 'scale_d' explicitly
                                  dff_d = t2$dff_d,            # Access 'dff_d' explicitly
                                  de_an_prior = 1   # Access 'de_an_prior' explicitly
                                ), "2" =
                                t1e_prior_plot(t2$prior_analysis,
                                               t2$location,
                                               t2$scale,
                                               t2$dff ,
                                               t2$alternative,
                                               t2$ROPE,
                                               1,
                                               t2$prior_design,
                                               t2$scale_d,
                                               t2$dff_d,
                                               t2$location_d )

      ))

    })
    r = t2$N2/t2$N1
    N1 = t2$N1
    ddff = t2$N1+t2$N2-2

    BF10 <- suppressWarnings(switch(t2$interval,
                                    "1" = t2_BF10(t2$tval,N1,r,t2$prior_analysis ,t2$location,t2$scale,t2$dff , t2$alternative ),
                                    "2" = t2e_BF10(t2$tval,N1,r,t2$prior_analysis,t2$location,t2$scale,t2$dff , t2$alternative,t2$ROPE )))


    output$result_t2 <- shiny::renderText({

      fmt_val <- function(x) {
        if (is.numeric(x) && length(x) == 1) return(as.character(x))
        if (is.numeric(x) && length(x) > 1) return(paste(x, collapse = ", "))
        if (is.character(x)) return(shQuote(x))
        return(as.character(x))
      }

      args <- list(
        tval = t2$tval,
        N1 = t2$N1,
        N2 = t2$N2,
        prior_analysis = t2$prior_analysis,
        location = t2$location,
        scale = t2$scale,
        dff = t2$dff,
        alternative = t2$alternative
      )

      if (t2$interval != 1) {
        args$ROPE <- t2$ROPE
      }

      # Build string with each argument on a new line, applying omission rules
      arg_strings <- sapply(names(args), function(arg) {

        val <- args[[arg]]
        arg_print <- arg

        # ---- Omission rules ----
        if ((t2$prior_analysis %in% c("Normal", "Moment")) && arg == "dff") {
          return(NULL)
        }

        # Handle alternative
        if (arg == "alternative") {
          val <- shQuote(val)
        } else if (arg == "ROPE") {
          # Wrap vector ROPE in c(...)
          if (length(val) > 1) {
            val <- paste0("c(", paste(val, collapse = ", "), ")")
          } else {
            val <- as.character(val)
          }
        } else {
          val <- fmt_val(val)
        }

        sprintf("  %s = %s", arg_print, val)
      })

      # Remove NULLs from omitted arguments
      arg_strings <- arg_strings[!sapply(arg_strings, is.null)]

      call_string <- paste0(
        "# Function to be used in R\n",
        "BF10.ttest.TwoSample(\n",
        paste(arg_strings, collapse = ",\n"),
        "\n)"
      )

      call_string
    })
    d.obs <-  t2$tval / sqrt((t2$N1 * t2$N2) / (t2$N1 + t2$N2))
    ROPE <- switch(t2$interval,"1" = NULL,"2" = t2$ROPE)
    p.value <- t.pval(t2$tval, t2$N1, t2$N2, t2$alternative, ROPE = ROPE, type = "two")
    output$BFt2 <- shiny::renderUI({
      # Create the LaTeX formatted strings for the table
      table_html <- paste0(
        '\\textit{t}(', ddff, ') = ', t2$tval,', \\textit{p} = ',round(p.value,4),
        ', \\textit{d} = ', round(d.obs, 4), ', \\\\ ',
        '\\textit{BF}_{10} = ', round(BF10, 4),
        ', \\textit{BF}_{01} = ', round(1/BF10, 4)
      )



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

t2_BF10 <-function(t,n1,r,prior_analysis ,location,scale,dff , alternative ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  bound  <- switch(alternative,
                   "greater" = c(a = 0, b = Inf),
                   "less" = c(a = -Inf, b = 0),
                   "two.sided" = c(a = -Inf, b = Inf)
  )

  normalization <- if (alternative == "two.sided") 1 else
    switch(prior_analysis,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))

  error = 1e-10
  x <- sapply(t, function(ti) {
    int <- function(delta) {
      stats::dt(ti, df, ncp = delta * constant) * t1_prior(delta, location, scale, dff, prior_analysis) / normalization
    }

    stats::integrate(int, lower = bound[1], upper = bound[2], rel.tol = error, stop.on.error = FALSE)$value /
      stats::dt(ti, df, ncp = 0)
  })

  return(x)
}



# finding the t that correspond to BF10=D
t2_BF10_bound <-function(threshold, n1,r,prior_analysis ,location ,scale,dff , alternative){
  y <- numeric(0)
  Bound_finding <-function(t){
    t2_BF10(t,n1,r,prior_analysis=prior_analysis,location=location,scale=scale,dff=dff, alternative =alternative )- threshold
  }
  x <- tryCatch(stats::uniroot(Bound_finding, lower = -6, upper = 0)$root, error = function(e) NA)
  y <- tryCatch(stats::uniroot(Bound_finding, lower =  0, upper = 6)$root, error = function(e) NA)
  results <- c(x, y)

  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")

  BF.vals  <- t2_BF10(results,n1,r,prior_analysis=prior_analysis,location=location,scale=scale,dff=dff, alternative =alternative )
  BF.close <- which(round(BF.vals, 2) == round(threshold, 2))
  if (length(BF.close) == 0) return("bound cannot be found")

  return(results[BF.close])
}


# finding the t that correspond to BF01=D
t2_BF01_bound <-function(threshold , n1,r,prior_analysis ,location ,scale,dff , alternative){
  t2_BF10_bound(1/threshold, n1,r,prior_analysis ,location ,scale,dff , alternative)
}


# p(BF01>D|H0)
t2_TNE <- function(t , n1,r,alternative){
  n2 = n1*r
  df = n1+n2-2

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  pro <- switch(alternative,
                "two.sided" = stats::pt(max(t), df) - stats::pt(min(t), df),
                "greater"  = stats::pt(t, df),
                "less"  = 1 - stats::pt(t, df)
  )

  return(pro)

}

# p(BF10>D|H1)
t2_TPE <-function(t,n1,r,prior_analysis ,location ,scale,dff , alternative ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (prior_analysis == "Point"){
    pro = switch(alternative,
                 "two.sided"= pnct(min(t),df,ncp = location*constant,lower  = T)+pnct(max(t),df,ncp = location*constant,lower  = F),
                 "greater" = pnct(t,df,ncp = location*constant,lower  = F),
                 "less" = pnct(t,df,ncp = location*constant,lower  = T))
    return(pro)
  }

  bound  <- switch(alternative,
                   "greater" = c(a = 0, b = Inf),
                   "less" = c(a = -Inf, b = 0),
                   "two.sided" = c(a = -Inf, b = Inf)
  )



  x = NULL

  normalization <- if (alternative == "two.sided") 1 else
    switch(prior_analysis,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(alternative,
                  "two.sided" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = FALSE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 + pro2
                  },
                  "greater" = pnct(t, df, ncp = ncp, lower = FALSE),
                  "less" = pnct(t, df, ncp = ncp, lower = TRUE)
    )

    pro * t1_prior(delta, location, scale, dff, prior_analysis) / normalization
  }

  error = 1e-4
  if (prior_analysis == "Moment" & scale <.3 ){
    error = 1e-14
  }
  x = stats::integrate(int,lower = bound[1],upper = bound[2], rel.tol = error,stop.on.error=FALSE)$value

  return(x)

}


# p(BF01>D|H1)
t2_FNE<-function(t,n1,r,prior_analysis ,location ,scale,dff , alternative ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  if (prior_analysis == "Point"){
    pro = switch(alternative,
                 "two.sided"=  pnct(max(t),df,ncp = location*constant,lower  = T) - pnct(min(t),df,ncp = location*constant,lower  = T),
                 "greater" = pnct(t,df,ncp = location*constant,lower  = T),
                 "less" = pnct(t,df,ncp = location*constant,lower  = F))
    return(pro)
  }
  bound  <- switch(alternative,
                   "greater" = c(a = 0, b = Inf),
                   "less" = c(a = -Inf, b = 0),
                   "two.sided" = c(a = -Inf, b = Inf)
  )
  x = NULL


  normalization <- if (alternative == "two.sided") 1 else
    switch(prior_analysis,
           "Cauchy"         = stats::pcauchy(bound[2], location, scale)     - stats::pcauchy(bound[1], location, scale),
           "Normal"         = stats::pnorm (bound[2], location, scale)      - stats::pnorm (bound[1], location, scale),
           "Moment"            = if (bound[2] == 0) pmom(bound[2]-location, tau=scale^2) else 1-pmom(bound[1]-location, tau=scale^2),
           "t-distribution" = stats::pt((bound[2] - location) / scale, dff, 0) - stats::pt((bound[1] - location) / scale, dff, 0))
  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(alternative,
                  "two.sided" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = TRUE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 - pro2
                  },
                  "greater" = pnct(t, df, ncp = ncp, lower = TRUE),
                  "less" = pnct(t, df, ncp = ncp, lower = FALSE)
    )

    pro * t1_prior(delta, location, scale, dff, prior_analysis) / normalization
  }


  x = stats::integrate(int,lower = bound[1],upper = bound[2],stop.on.error = FALSE)$value

  return(x)
}


# p(BF10>D|H0)
t2_FPE <- function(t,n1,r, alternative){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  pro <- switch(alternative,
                "two.sided" = stats::pt(max(t), df = df, lower.tail = FALSE) +
                  stats::pt(min(t), df = df, lower.tail = TRUE),
                "greater"  = stats::pt(t, df = df, lower.tail = FALSE),
                "less"  = stats::pt(t, df = df, lower.tail = TRUE)
  )
  return(pro)

}

# Finding the degree of freedom that ensure p(BF10>D|H1) > targeted probability


t2_N_finder<-function(threshold,r,true_rate,prior_analysis,location,scale,dff, alternative ,
                      prior_design,location_d,scale_d,dff_d,de_an_prior ,false_rate){

  lower <- 2
  upper <- 10000
  t2 <- t2_BF10_bound(threshold, lower,r,prior_analysis ,location ,scale,dff , alternative)
  p2 <- if (de_an_prior == 1)
    t2_TPE(t2 , n1=lower,r , prior_analysis , location ,scale,dff , alternative) else
      t2_TPE(t2 , n1=lower,r , prior_design , location_d,scale_d,dff_d, alternative)
  if (p2 > true_rate) return(lower)
  Power_root <- function(n1) {
    t <- t2_BF10_bound(threshold, n1,r,prior_analysis ,location ,scale,dff , alternative)
    if (de_an_prior == 1)
      t2_TPE(t , n1,r , prior_analysis , location ,scale,dff , alternative) - true_rate else
        t2_TPE(t , n1,r , prior_design , location_d,scale_d,dff_d, alternative) - true_rate
  }
  N1.power <-  stats::uniroot(Power_root,lower = lower,upper =  upper)$root
  #N1.power <- robust_uniroot(Power_root, lower = 2)
  t  <-  t2_BF10_bound(threshold,  N1.power,r,prior_analysis ,location ,scale,dff , alternative)
  FPE <- t2_FPE(t,N1.power,r, alternative)
  if (FPE <= false_rate) return(N1.power)
  alpha.root <- function(n1) {
    t <- t2_BF10_bound(threshold,  n1,r,prior_analysis ,location ,scale,dff , alternative)
    t2_FPE(t,n1,r, alternative) - false_rate
  }

  N1.alpha <- stats::uniroot(alpha.root, lower = N1.power, upper = upper)$root
  return(N1.alpha  )
}

t2_N_01_finder<-function(threshold,r,true_rate,prior_analysis,location,scale,dff, alternative ,
                         prior_design,location_d,scale_d,dff_d,de_an_prior ,false_rate){

  lower <- 2
  upper <- 10000
  t2 <- t2_BF01_bound(threshold, lower,r,prior_analysis ,location ,scale,dff , alternative)
  TNE_lo <- t2_TNE(t2,lower,r,alternative)
  if (TNE_lo > true_rate) return(lower)
  FNE_lo <-  if (de_an_prior == 1)
    t2_FNE(t2,lower,r,prior_analysis ,location ,scale,dff , alternative ) else
      t2_FNE(t2,lower,r,prior_design ,location_d ,scale_d,dff_d , alternative )
  if (TNE_lo > true_rate&FNE_lo<false_rate) return(lower)

  TN_root <- function(n1) {
    t <- t2_BF01_bound(threshold, n1,r,prior_analysis ,location ,scale,dff , alternative)
    t2_TNE(t,lower,r,alternative)-true_rate
  }
  N1.TN <-  stats::uniroot(TN_root,lower = lower,upper =  upper)$root
  t  <-  t2_BF01_bound(threshold,  N1.TN,r,prior_analysis ,location ,scale,dff , alternative)
  FNE <- if (de_an_prior == 1)
    t2_FNE(t,N1.TN,r,prior_analysis ,location ,scale,dff , alternative ) else
      t2_FNE(t,N1.TN,r,prior_design ,location_d ,scale_d,dff_d , alternative )
  if (FNE <= false_rate) return(N1.TN)

  FN.root <- function(n1) {
    t <- t2_BF01_bound(threshold,  n1,r,prior_analysis ,location ,scale,dff , alternative)
    FNE <- if (de_an_prior == 1)
      t2_FNE(t,n1,r,prior_analysis ,location ,scale,dff , alternative ) else
        t2_FNE(t,n1,r,prior_design ,location_d ,scale_d,dff_d , alternative )
    FNE -false_rate
  }
  N1.FN <- stats::uniroot(FN.root, lower = N1.TN, upper = upper)$root
  return( N1.FN )
}

# probability table
t2_Table <- function(threshold,r,true_rate,prior_analysis,location,scale,dff, alternative,
                     prior_design,location_d,scale_d,dff_d, de_an_prior,N1,N2, mode_bf ,false_rate ,type_rate){

  bound01 = as.numeric(0)
  bound10 = as.numeric(0)

  if (mode_bf == 1) {
    n1 <- switch(type_rate,
                 "positive" = {ceiling(t2_N_finder(threshold, r, true_rate, prior_analysis, location, scale, dff,
                                                   alternative, prior_design, location_d, scale_d, dff_d,
                                                   de_an_prior, false_rate))},
                 "negative" = {ceiling(t2_N_01_finder(threshold, r, true_rate, prior_analysis, location, scale, dff,
                                                      alternative, prior_design, location_d, scale_d, dff_d,
                                                      de_an_prior, false_rate))}     )
    n2 <- n1 * r
  } else {
    n1 <- N1
    n2 <- N2
    r  <- n2 / n1
  }

  # t bounds:
  t10 <- t2_BF10_bound(threshold, n1,r,prior_analysis,location,scale,dff , alternative)
  t01 <- t2_BF01_bound(threshold, n1,r,prior_analysis,location,scale,dff , alternative)

  # max BF10 possible:
  max_BF <- 1 / t2_BF10(0,n1,r,prior_analysis ,location,scale,dff , alternative )
  BF_D   <- t10

  # FPE and TPE:
  FPE       <- t2_FPE(t10,n1,r, alternative)
  if (de_an_prior == 1) {
    TPE       <- t2_TPE(t10,n1,r,prior_analysis ,location ,scale,dff , alternative )
    TPR_prior <- prior_analysis
    TPR_loc   <- location
    TPR_scale <- scale
    TPR_dff   <- dff
  } else {
    TPE       <- t2_TPE(t10,n1,r,prior_design ,location_d ,scale_d,dff_d , alternative )
    TPR_prior <- prior_design
    TPR_loc   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }

  # FNE and TNE:
  if (any(alternative == "two.sided" & max_BF < threshold | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t2_FNE(t01, n1,r,TPR_prior, TPR_loc, TPR_scale, TPR_dff, alternative )
    TNE <- t2_TNE(t01 , n1,r,alternative)
  }

  # table:
  tab.names <- c(
    "TruePositve",
    "FalseNegative",
    "TrueNegative",
    "FalsePositive",
    "Required N1",
    "Required N2"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n1,n2, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}


# plots for showing the relationship between BF and t-values

t2_BF <- function(threshold, n1, r, target,
                  prior_analysis, location, scale, dff, alternative) {

  tt <- seq(-5, 5, 0.2)

  ## ---------- BF10 ----------
  BF10   <- t2_BF10(tt, n1, r, prior_analysis, location, scale, dff, alternative)
  t.BF10 <- t2_BF10_bound(threshold, n1, r, prior_analysis, location, scale, dff, alternative)

  df10 <- data.frame(t = tt, BF = BF10)

  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(threshold)~" when t = "~.(round(t.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(threshold)~" when t = "~.(round(t.BF10[1], 2))~
                  " or "~.(round(t.BF10[2], 2))))
  }

  x_breaks_10 <- sort(unique(c(-5, 5, round(t.BF10, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(t, BF)) +
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
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
  t.BF01 <- t2_BF01_bound(threshold, n1, r, prior_analysis, location, scale, dff, alternative)

  df01 <- data.frame(t = tt, BF = BF01)

  max.BF01   <- 1 / t2_BF10(0, n1, r, prior_analysis, location, scale, dff, "two.sided")
  impossible <- (alternative == "two.sided") &&
    (max.BF01 < threshold || identical(t.BF01, "bound cannot be found"))

  if (impossible) {

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(linewidth = 1.2, color = "black") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = c(-5, 5)) +
      ggplot2::labs(
        x = "t-value",
        y = bquote("BF"['01'] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[01]~"="~.(threshold)))
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
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when t = "~.(round(t.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when t = "~.(round(t.BF01[1], 2))~
                    " or "~.(round(t.BF01[2], 2))))
    }

    x_breaks_01 <- sort(unique(c(-5, 5, round(t.BF01, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(linewidth = 1.2, color = "black") +
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



Power_t2 <- function(threshold, prior_analysis, location, scale, dff, alternative,
                     prior_design, location_d, scale_d, dff_d,
                     de_an_prior, n1, r) {

  Total_ <- n1 + n1 * r
  smin <- 4
  smax <- Total_ * 1.2
  sdf  <- seq(smin, smax, length.out = 31)
  sn1  <- sdf / (1 + r)

  # Initialize vectors
  TPE <- FPE <- TNE <- FNE <- numeric(length(sdf))

  for (i in seq_along(sdf)) {

    t10 <- t2_BF10_bound(threshold, sn1[i], r, prior_analysis, location, scale, dff, alternative)
    t01 <- t2_BF01_bound(threshold, sn1[i], r, prior_analysis, location, scale, dff, alternative)

    TPE[i] <- if (de_an_prior == 1) {
      t2_TPE(t10, sn1[i], r, prior_analysis, location, scale, dff, alternative)
    } else {
      t2_TPE(t10, sn1[i], r, prior_design, location_d, scale_d, dff_d, alternative)
    }

    FPE[i] <- t2_FPE(t10, sn1[i], r, alternative)

    FNE[i] <- if (de_an_prior == 1) {
      t2_FNE(t01, sn1[i], r, prior_analysis, location, scale, dff, alternative)
    } else {
      t2_FNE(t01, sn1[i], r, prior_design, location_d, scale_d, dff_d, alternative)
    }

    TNE[i] <- t2_TNE(t01, sn1[i], r, alternative)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  p2 <- ggplot2::ggplot(df2,
                        ggplot2::aes(x = SampleSize,
                                     y = Probability,
                                     color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  ## ---------- Combine ----------

  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}


# ---- twosample_e.r ----

t2e_BF10i <-function(t,n1,r,prior_analysis ,location,scale,dff , alternative,ROPE ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  bound_h1  <- switch(alternative,
                      "greater" = c(a = ROPE, b = Inf),
                      "less" = c(a = -Inf, b = ROPE),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )
  bound_h0  <- switch(alternative,
                      "greater" = c(a = 0, b = ROPE),
                      "less" = c(a = ROPE, b = 0),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )
  normalizationh1 <- norm_h1(alternative, prior_analysis, bound_h1, location, scale, dff)


  normalizationh0 <- norm_h0(prior_analysis, bound_h0, location, scale, dff)


  int  <- function(delta){
    stats::dt(t,df,ncp = delta *constant)* te_prior(delta,location,scale,dff,prior_analysis)/normalizationh1}

  error = 1e-4

  if (alternative == "two.sided"){
    lh1 = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value
  }else{
    lh1 = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value}


  int  <- function(delta){
    stats::dt(t,df,ncp = delta *constant)* te_prior(delta,location,scale,dff,prior_analysis)/normalizationh0}

  lh0 = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(lh1/lh0)
}

t2e_BF10 <-function(t,n1,r,prior_analysis,location,scale,dff , alternative,ROPE ){
  sapply(t, function(ti) t2e_BF10i(ti,n1,r,prior_analysis ,location,scale,dff , alternative,ROPE ))
}

t2e_BF10_bound <-function(threshold, n1,r,prior_analysis,location,scale,dff , alternative,ROPE){

  y <- numeric(0)
  Bound_finding <-function(t){
    t2e_BF10(t,n1,r,prior_analysis,location,scale,dff , alternative,ROPE )- threshold
  }

  switch(alternative,
         "two.sided" ={
           x <- tryCatch(stats::uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
           y <- tryCatch(stats::uniroot(Bound_finding, lower =  0, upper = 20)$root, error = function(e) NA)
         },
         "greater"={
           x <- tryCatch(stats::uniroot(Bound_finding, lower = 0, upper = 20)$root, error = function(e) NA)
         },
         "less" = {
           x <- tryCatch(stats::uniroot(Bound_finding, lower = -20, upper = 0)$root, error = function(e) NA)
         })

  results <- c(x, y)
  results <- results[!is.na(results)]
  if (length(results) == 0) return("bound cannot be found")
  BF.vals  <- t2e_BF10(results,n1,r,prior_analysis,location,scale,dff , alternative,ROPE )
  BF.close <- which(round(BF.vals, 2) == round(threshold, 2))
  if (length(BF.close) == 0) return("bound cannot be found")
  return(results[BF.close])
}

# finding the t that correspond to BF10=D
t2e_BF01_bound <-function(threshold, n1,r,prior_analysis,location,scale,dff , alternative,ROPE){
  t2e_BF10_bound(1/threshold, n1,r,prior_analysis,location,scale,dff , alternative,ROPE)

}

t2e_TPE <-function(t,n1,r,prior_analysis ,location,scale,dff , alternative ,ROPE){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  if (prior_analysis =="Point"){
    x = switch(alternative,
               "two.sided" = {pnct(min(t),df,ncp= location*constant,lower = T)+ pnct(max(t),df,ncp=location*constant,lower = F)},
               "less"  = {pnct(t,df,ncp = location *constant,lower  = T)},
               "greater"  = {pnct(t,df,ncp = location *constant,lower  = F)}
    )
    return(x)
  }

  bound_h1  <- switch(alternative,
                      "greater" = c(a = ROPE, b = Inf),
                      "less" = c(a = -Inf, b = ROPE),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )

  normalizationh1 <- norm_h1(alternative, prior_analysis, bound_h1, location, scale, dff)



  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(alternative,
                  "two.sided" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = FALSE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 + pro2
                  },
                  "greater" = pnct(t, df, ncp = ncp, lower = FALSE),
                  "less" = pnct(t, df, ncp = ncp, lower = TRUE)
    )

    pro * te_prior(delta,location, scale, dff, prior_analysis) / normalizationh1
  }


  error = 1e-4

  if (alternative == "two.sided"){
    x = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value
  }else{
    x = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value

  }
  return(x)

}

t2e_FNE <-function(t,n1,r,prior_analysis ,location,scale,dff , alternative ,ROPE){

  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  if (prior_analysis =="Point"){
    x = switch(alternative,
               "two.sided" = {pnct(max(t),df,ncp= location*constant,lower = T)- pnct(min(t),df,ncp=location*constant,lower = T)},
               "less"  = {pnct(t,df,ncp = location *constant,lower  = F)},
               "greater"  = {pnct(t,df,ncp = location *constant,lower  = T)}
    )
    return(x)
  }
  bound_h1  <- switch(alternative,
                      "greater" = c(a = ROPE, b = Inf),
                      "less" = c(a = -Inf, b = ROPE),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )

  normalizationh1 <- norm_h1(alternative, prior_analysis, bound_h1, location, scale, dff)


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(alternative,
                  "two.sided" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = TRUE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 - pro2
                  },
                  "greater" = pnct(t, df, ncp = ncp, lower = TRUE),
                  "less" = pnct(t, df, ncp = ncp, lower = FALSE)
    )

    pro * te_prior(delta, location,scale, dff, prior_analysis) / normalizationh1
  }

  error = 1e-4
  if (alternative == "two.sided"){
    x = stats::integrate(int,lower = -Inf,upper = bound_h1[1], rel.tol=error,stop.on.error = F)$value+stats::integrate(int,lower =  bound_h1[2],upper = Inf, rel.tol=error,stop.on.error = F)$value
  }else{
    x = stats::integrate(int,lower = bound_h1[1],upper = bound_h1[2], rel.tol=error,stop.on.error = F)$value

  }
  return(x)

}


t2e_TNE <-function(t,n1,r,prior_analysis ,location,scale,dff , alternative ,ROPE){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))


  bound_h0  <- switch(alternative,
                      "greater" = c(a = 0, b = ROPE),
                      "less" = c(a = ROPE, b = 0),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )

  normalizationh0 <- norm_h0(prior_analysis, bound_h0, location, scale, dff)


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(alternative,
                  "two.sided" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = TRUE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 - pro2
                  },
                  "greater" = pnct(t, df, ncp = ncp, lower = TRUE),
                  "less" = pnct(t, df, ncp = ncp, lower = FALSE)
    )

    pro * te_prior(delta,location, scale, dff, prior_analysis) / normalizationh0
  }
  error = 1e-4
  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value
  return(x)

}

t2e_FPE <-function(t,n1,r,prior_analysis ,location,scale,dff , alternative ,ROPE){
  if (any(t == "bound cannot be found") || length(t) == 0) return(0)

  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))

  bound_h0  <- switch(alternative,
                      "greater" = c(a = 0, b = ROPE),
                      "less" = c(a = ROPE, b = 0),
                      "two.sided" = c(a = ROPE[1], b = ROPE[2])
  )

  normalizationh0 <- norm_h0(prior_analysis, bound_h0, location, scale, dff)


  int <- function(delta) {
    ncp <- delta * constant

    pro <- switch(alternative,
                  "two.sided" = {
                    pro1 <- pnct(max(t), df, ncp = ncp, lower = FALSE)
                    pro2 <- pnct(min(t), df, ncp = ncp, lower = TRUE)
                    pro1 + pro2
                  },
                  "greater" = pnct(t, df, ncp = ncp, lower = FALSE),
                  "less" = pnct(t, df, ncp = ncp, lower = TRUE)
    )

    pro * te_prior(delta,location, scale, dff, prior_analysis) / normalizationh0
  }
  error = 1e-4
  x = stats::integrate(int,lower = bound_h0[1],upper = bound_h0[2], rel.tol=error,stop.on.error = F)$value

  return(x)

}


t2e_N_finder<-function(threshold,r,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                       prior_design,location_d,scale_d,dff_d, de_an_prior,false_rate ){

  lower <- 2
  t2 <-t2e_BF10_bound(threshold, lower,r,prior_analysis,location,scale,dff , alternative,ROPE)
  p2 <- if (de_an_prior == 1)
    t2e_TPE (t2,lower,r,prior_analysis ,location,scale,dff , alternative,ROPE) else
      t2e_TPE (t2,lower,r,prior_design ,location_d,scale_d,dff_d , alternative,ROPE)
  if (p2 > true_rate) return(lower)

  Power_root <- function(n1) {

    t <- t2e_BF10_bound(threshold, n1,r,prior_analysis,location,scale,dff , alternative,ROPE)

    pro <- if (de_an_prior == 1) {
      t2e_TPE (t,n1,r,prior_analysis ,location,scale,dff , alternative,ROPE )
    } else {
      t2e_TPE (t,n1,r,prior_design ,location_d,scale_d,dff_d , alternative,ROPE)
    }

    true_rate - pro
  }
  N1.power <- robust_uniroot(Power_root, lower = 2)
  t <- t2e_BF10_bound(threshold, N1.power,r,prior_analysis,location,scale,dff , alternative,ROPE)
  FPE <-t2e_FPE(t,N1.power,r,prior_analysis,location,scale,dff , alternative ,ROPE)

  if (FPE <= false_rate) return(N1.power + 1)

  alpha.root <- function(n1) {
    t <- t2e_BF10_bound(threshold, n1,r,prior_analysis,location,scale,dff , alternative,ROPE)
    pro <- t2e_FPE(t,n1,r,prior_analysis,location ,scale,dff , alternative ,ROPE)
    return(pro - false_rate)
  }
  N1.alpha <- robust_uniroot(alpha.root , lower = N1.power)
  return(N1.alpha)
}
t2e_N_01_finder<-function(threshold,r,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                          prior_design,location_d,scale_d,dff_d, de_an_prior,false_rate ){

  lower <- 10
  t2 <-t2e_BF01_bound(threshold, lower,r,mode,location,scale,dff , alternative,ROPE)
  TNE_lo <- t2e_TNE(t2,lower,r,prior_analysis ,location,scale,dff , alternative,ROPE)
  if (TNE_lo > true_rate) return(lower)
  FNE_lo <- if (de_an_prior == 1)
    t2e_FNE (t2,lower,r,prior_analysis,location ,scale,dff , alternative,ROPE ) else
      t2e_FNE (t2,lower,r,prior_design ,location_d,scale_d,dff_d , alternative,ROPE )
  if (TNE_lo > true_rate&FNE_lo<false_rate) return(lower)

  TN_root <- function(n1) {

    t <- t2e_BF01_bound(threshold, n1,r,prior_analysis,location,scale,dff , alternative,ROPE)

    pro <- t2e_TNE(t,n1,r,prior_analysis ,location,scale,dff , alternative,ROPE)

    true_rate - pro
  }
  N1.TN <- robust_uniroot(TN_root, lower = 2)
  t <- t2e_BF01_bound(threshold, N1.TN,r,prior_analysis,location,scale,dff , alternative,ROPE)
  FNE <-t2e_FNE(t,N1.TN,r,prior_analysis ,location,scale,dff , alternative ,ROPE)

  if (FNE <= false_rate) return(N1.TN + 1)

  FN.root <- function(n1) {
    t <- t2e_BF01_bound(threshold, n1,r,prior_analysis,location,scale,dff , alternative,ROPE)
    pro <- if (de_an_prior == 1) {
      t2e_FNE (t,n1,r,prior_analysis ,location,scale,dff , alternative,ROPE )
    } else {
      t2e_FNE (t,n1,r,prior_design ,location_d,scale_d,dff_d , alternative,ROPE)
    }
    return(pro - false_rate)
  }
  N1.FN <- robust_uniroot(FN.root , lower = N1.TN)
  return(N1.FN)
}

t2e_table<-function(threshold,r,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                    prior_design,location_d,scale_d,dff_d, de_an_prior,mode_bf,N1,N2,false_rate ,type_rate){
  bound01 = as.numeric(0)
  bound10 = as.numeric(0)

  if (mode_bf == 1){

    n1 = switch(type_rate,
                "positive" = ceiling(t2e_N_finder(threshold,r,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                                                  prior_design,location_d,scale_d,dff_d, de_an_prior,false_rate )),
                "negative" = ceiling(t2e_N_01_finder(threshold,r,true_rate,prior_analysis,location,scale,dff, alternative,ROPE ,
                                                     prior_design,location_d,scale_d,dff_d, de_an_prior,false_rate ) ))
    n2 = n1*r
  } else {
    n1 = N1
    n2 = N2
    r= n2/n1
  }
  # t bounds:
  t10 <- t2e_BF10_bound(threshold, n1,r,prior_analysis,location,scale,dff , alternative,ROPE)
  t01 <- t2e_BF01_bound(threshold, n1,r,prior_analysis,location,scale,dff , alternative,ROPE)

  # max BF10 possible:
  max_BF <- 1 / t2e_BF10i(0,n1,r,prior_analysis ,location,scale,dff , alternative,ROPE )
  BF_D   <- t10

  # FPE and TPE:
  FPE       <- t2e_FPE(t10,n1,r,prior_analysis,location ,scale,dff , alternative ,ROPE)
  if (de_an_prior == 1) {
    TPE       <- t2e_TPE(t10,n1,r,prior_analysis ,location,scale,dff , alternative,ROPE )
    TPR_prior <- prior_analysis
    TPR_location   <- location
    TPR_scale <- scale
    TPR_dff   <- dff
  } else {
    TPE       <- t2e_TPE(t10,n1,r,prior_design ,location_d,scale_d,dff_d , alternative,ROPE )
    TPR_prior <- prior_design
    TPR_location   <- location_d
    TPR_scale <- scale_d
    TPR_dff   <- dff_d
  }

  # FNE and TNE:
  if (any(alternative == "two.sided" & max_BF < threshold | BF_D == "bound cannot be found")) {
    FNE <- 0
    TNE <- 0
  } else {
    FNE <- t2e_FNE(t01,n1,r,TPR_prior ,TPR_location,TPR_scale,TPR_dff , alternative ,ROPE)
    TNE <- t2e_TNE(t01,n1,r,prior_analysis,location ,scale,dff , alternative ,ROPE)
  }

  # table:
  tab.names <- c(
    "TruePositve",
    "FalseNegative",
    "TrueNegative",
    "FalsePositive",
    "Required N1",
    "Required N2"
  )
  table <- data.frame(TPE, FNE, TNE, FPE, n1,n2, check.names = FALSE, row.names = NULL)
  colnames(table) <- tab.names
  table
}

t2e_BF <- function(threshold, n1, r,
                   prior_analysis, location, scale, dff, alternative, ROPE) {

  tt <- seq(-5, 5, 0.2)
  xlim_range <- c(-5, 5)  # plot limits

  ## ---------- BF10 ----------
  BF10   <- t2e_BF10(tt, n1, r, prior_analysis, location, scale, dff, alternative, ROPE)
  t.BF10 <- t2e_BF10_bound(threshold, n1, r, prior_analysis, location, scale, dff, alternative, ROPE)

  df10 <- data.frame(t = tt, BF = BF10)

  main.bf10 <- if (length(t.BF10) == 1) {
    bquote(bold("BF"[10]~"="~.(threshold)~" when t = "~.(round(t.BF10, 2))))
  } else {
    bquote(bold("BF"[10]~"="~.(threshold)~" when t = "~.(round(t.BF10[1], 2))~
                  " or "~.(round(t.BF10[2], 2))))
  }

  # Keep only BF10 bounds inside plot
  t.BF10_plot <- t.BF10[t.BF10 >= xlim_range[1] & t.BF10 <= xlim_range[2]]

  x_breaks_10 <- sort(unique(c(xlim_range, round(t.BF10_plot, 2))))

  p1 <- ggplot2::ggplot(df10, ggplot2::aes(t, BF)) +
    ggplot2::geom_line(linewidth = 1.2, color = "black") +
    ggplot2::geom_vline(xintercept = t.BF10_plot, linetype = "dashed") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(limits = xlim_range, breaks = x_breaks_10) +
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
  t.BF01 <- t2e_BF01_bound(threshold, n1, r, prior_analysis, location, scale, dff, alternative, ROPE)

  df01 <- data.frame(t = tt, BF = BF01)

  max.BF01   <- 1 / t2e_BF10i(0, n1, r, prior_analysis, location, scale, dff, alternative, ROPE)
  impossible <- (max.BF01 < threshold || identical(t.BF01, "bound cannot be found"))

  if (impossible) {
    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(linewidth = 1.2, color = "black") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = xlim_range, breaks = xlim_range) +
      ggplot2::labs(
        x = "t-value",
        y = bquote("BF"['01'] * " (log scale)"),
        title = bquote(bold("It is impossible to have BF"[01]~"="~.(threshold)))
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
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when t = "~.(round(t.BF01, 2))))
    } else {
      bquote(bold("BF"[0][1]~"="~.(threshold)~" when t = "~.(round(t.BF01[1], 2))~
                    " or "~.(round(t.BF01[2], 2))))
    }

    # Keep only BF01 bounds inside plot
    t.BF01_plot <- t.BF01[t.BF01 >= xlim_range[1] & t.BF01 <= xlim_range[2]]

    x_breaks_01 <- sort(unique(c(xlim_range, round(t.BF01_plot, 2))))

    p2 <- ggplot2::ggplot(df01, ggplot2::aes(t, BF)) +
      ggplot2::geom_line(linewidth = 1.2, color = "black") +
      ggplot2::geom_vline(xintercept = t.BF01_plot, linetype = "dashed") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(limits = xlim_range, breaks = x_breaks_01) +
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
Power_t2e <- function(threshold, prior_analysis, location, scale, dff, alternative,
                      prior_design, location_d, scale_d, dff_d,
                      de_an_prior, n1, r, ROPE) {

  Total_ <- n1 + n1 * r
  smin   <- 4
  smax   <- Total_ * 1.2
  sdf    <- seq(smin, smax, length.out = 31)
  sn1    <- sdf / (1 + r)

  TPE <- FPE <- TNE <- FNE <- numeric(length(sdf))

  for (i in seq_along(sdf)) {

    t10 <- t2e_BF10_bound(threshold, sn1[i], r, prior_analysis, location, scale, dff, alternative, ROPE)
    t01 <- t2e_BF01_bound(threshold, sn1[i], r, prior_analysis, location, scale, dff, alternative, ROPE)

    TPE[i] <- if (de_an_prior == 1) {
      t2e_TPE(t10, sn1[i], r, prior_analysis, location, scale, dff, alternative, ROPE)
    } else {
      t2e_TPE(t10, sn1[i], r, prior_design, location_d, scale_d, dff_d, alternative, ROPE)
    }

    FPE[i] <- t2e_FPE(t10, sn1[i], r, prior_analysis, location, scale, dff, alternative, ROPE)
    TNE[i] <- t2e_TNE(t01, sn1[i], r, prior_analysis, location, scale, dff, alternative, ROPE)

    FNE[i] <- if (de_an_prior == 1) {
      t2e_FNE(t01, sn1[i], r, prior_analysis, location, scale, dff, alternative, ROPE)
    } else {
      t2e_FNE(t01, sn1[i], r, prior_design, location_d, scale_d, dff_d, alternative, ROPE)
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
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[10]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  ## ---------- BF01 plot ----------
  p2 <- ggplot2::ggplot(df2,
                        ggplot2::aes(x = SampleSize, y = Probability, color = Type)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = type_colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "Total sample size",
      y = "Probability",
      title = bquote(bold("Power curve for BF"[0][1]~">"~.(threshold)))
    ) +
    clean_theme +
    legend_theme

  ## ---------- Combine ----------
  print(patchwork::wrap_plots(p1, p2, ncol = 2))
}

