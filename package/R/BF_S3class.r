##### T-tests #####
#' @export
#' @method print BFpower_t
print.BFpower_t <- function (x,...) {
  # Box symbols (Unicode escapes)
  UL <- "\u2554"; UR <- "\u2557"
  LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  digits <- 3

  # Subscript-safe labels
  H0 <- "H\u2080"
  H1 <- "H\u2081"

  # Symbols
  delta <- "\u03b4"
  neq   <- "\u2260"
  leq   <- "\u2264"
  geq   <- "\u2265"
  union <- "\u222a"

  # Header with centered text
  header_text <- ifelse(x$mode_bf == 1, "SAMPLE SIZE CALCULATION", "POWER CALCULATION")
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypothesis section
  cat("Hypotheses\n")
  cat(line)

  if (is.null(x$ROPE)) {
    if (x$alternative == "!=") {
      cat("  ", H0, " (Null)         : ", delta, " = 0\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " ", neq, " 0\n\n", sep = "")
    }
    if (x$alternative == ">") {
      cat("  ", H0, " (Null)         : ", delta, " = 0\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " > 0\n\n", sep = "")
    }
    if (x$alternative == "<") {
      cat("  ", H0, " (Null)         : ", delta, " = 0\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " < 0\n\n", sep = "")
    }
  } else {
    if (x$alternative == "!=") {
      cat("  ", H0, " (Null)         : ",
          min(x$ROPE), " ", leq, " ", delta, " ", leq, " ", max(x$ROPE), "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ",
          delta, " < ", min(x$ROPE), " ", union, " ",
          delta, " > ", max(x$ROPE), "\n\n", sep = "")
    }
    if (x$alternative == ">") {
      cat("  ", H0, " (Null)         : 0 ", leq, " ", delta, " ", leq, " ", x$ROPE, "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " > ", x$ROPE, "\n\n", sep = "")
    }
    if (x$alternative == "<") {
      cat("  ", H0, " (Null)         : ", x$ROPE, " ", leq, " ", delta, " ", leq, " 0\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " < ", x$ROPE, "\n\n", sep = "")
    }
  }

  cat(line)
  cat("Analysis/design prior : \n")
  cat(line)

  # Prior under H0
  cat("Analysis/design prior under ", H0, ":\n", sep = "")
  if (is.null(x$ROPE)) {
    cat("  Point prior (location = 0)\n")
  } else {
    cat("  ", x$analysis_h1$prior,
        " (location = ", x$analysis_h1$prior,
        ", scale = ", x$analysis_h1$scale,
        if (!is.null(x$analysis_h1$dff)) paste0(", df = ", x$analysis_h1$dff),
        ")\n", sep = "")
  }

  # Prior under H1
  cat("Analysis prior under ", H1, ":\n", sep = "")
  cat( "  ",x$analysis_h1$prior,
       " (location = ", x$analysis_h1$prior,
       ", scale = ", x$analysis_h1$scale,
       if (!is.null(x$analysis_h1$dff)) paste0(", df = ", x$analysis_h1$dff),
       ")\n", sep = "")

  if (!is.null(x$design_h1$prior)) {
    cat("Design prior under ", H1, ":\n", sep = "")
    cat( "  ",x$design_h1$prior,
         " (location = ", x$design_h1$prior,
         ", scale = ", x$design_h1$scale,
         if (!is.null(x$design_h1$dff)) paste0(", df = ", x$design_h1$dff),
         ")\n", sep = "")
  }

  cat("\n")
  cat(line)
  cat("Threshold of compelling evidence: ", x$threshold, "\n\n", sep = "")
  cat(line)

  # Results
  cat("Results:\n")
  cat(line)
  cat(paste0("  True Positive Rate           = ", round(x$results[1], digits), "\n"))
  cat(paste0("  False Positive Rate          = ", round(x$results[4], digits), "\n"))
  cat(paste0("  True Negative Rate           = ", round(x$results[3], digits), "\n"))
  cat(paste0("  False Negative Rate          = ", round(x$results[2], digits), "\n"))

  # Required sample sizes (unlist for CRAN)
  if (x$type == "One-sample t-test") {
    cat("  Required Sample size         = ",
        round(unlist(x$results[5]), 1), "  <---\n", sep = "")
  } else {
    cat("  Required Sample size Group 1 = ",
        round(unlist(x$results[5]), 1), "  <---\n", sep = "")
    cat("  Required Sample size Group 2 = ",
        round(unlist(x$results[6]), 1), "  <---\n", sep = "")
  }
}


#' @export
#' @method plot BFpower_t
#' @keywords internal
plot.BFpower_t<-function (x,...){
  if(x$analysis_h1$prior!="t-distribution"){x$analysis_h1$dff=0}

  if (is.null(x$design_h1$prior)){
    de_an_prior=1
    x$design_h1$location=0
    x$design_h1$scale=0
    x$design_h1$dff=0
  }else{
    de_an_prior=0
    if(x$design_h1$prior!="t-distribution"){x$design_h1$dff=0}

  }
  if (is.null(x$ROPE)){
    t1_prior_plot(x$threshold, 0, x$analysis_h1$prior, x$analysis_h1$location, x$analysis_h1$scale, x$analysis_h1$dff,
                  x$alternative,
                  x$design_h1$prior, x$design_h1$location, x$design_h1$scale, x$design_h1$dff, de_an_prior)

  }else {

    t1e_prior_plot(x$analysis_h1$prior,x$analysis_h1$location, x$analysis_h1$scale, x$analysis_h1$dff,
                   x$alternative, x$ROPE,
                   de_an_prior, x$design_h1$prior, x$design_h1$scale, x$design_h1$dff, x$design_h1$location)

  }

  if(x$plot_power){
    if (is.null(x$ROPE)){
      if (x$type=="One-sample t-test"){
        suppressWarnings(Power_t1(x$threshold, x$analysis_h1$prior, x$analysis_h1$location, x$analysis_h1$scale, x$analysis_h1$dff,
                                  x$alternative,x$design_h1$prior, x$design_h1$location, x$design_h1$scale, x$design_h1$dff,de_an_prior, unlist(x$results[5])))


      }else{
       suppressWarnings(Power_t2(x$threshold,x$analysis_h1$prior, x$analysis_h1$location, x$analysis_h1$scale, x$analysis_h1$dff,
                                  x$alternative,
                                  x$design_h1$prior, x$design_h1$location, x$design_h1$scale, x$design_h1$dff,de_an_prior,x$results[1,5],x$results[1,6]/x$results[1,5]))

      }



    }else {
      if (x$type=="One-sample t-test"){
       suppressWarnings(Power_t1e(x$threshold,x$analysis_h1$prior, x$analysis_h1$location, x$analysis_h1$scale, x$analysis_h1$dff,
                                   x$alternative,x$design_h1$prior, x$design_h1$location, x$design_h1$scale, x$design_h1$dff,de_an_prior,unlist(x$results[5]) ,x$ROPE))

      } else{
        suppressWarnings(Power_t2e(x$threshold,x$analysis_h1$prior, x$analysis_h1$location, x$analysis_h1$scale, x$analysis_h1$dff,
                                   x$alternative,
                                   x$design_h1$prior, x$design_h1$location, x$design_h1$scale, x$design_h1$dff,de_an_prior,x$results[1,5],x$results[1,6]/x$results[1,5],x$ROPE))
      }
    }
  }


  if(x$plot_rel){
    if (is.null(x$ROPE)){

      if (x$type=="One-sample t-test"){
        suppressWarnings(bf10_t1(x$threshold, unlist(x$result[5])-1, 0,x$analysis_h1$prior , x$analysis_h1$location , x$analysis_h1$scale , x$analysis_h1$dff , x$alternative))

      }else{
        suppressWarnings(t2_BF(x$threshold ,x$result[1,5],x$results[1,6]/x$results[1,5], 0,x$analysis_h1$prior , x$analysis_h1$location , x$analysis_h1$scale , x$analysis_h1$dff , x$alternative ))
      }


    }else{
      if (x$type=="One-sample t-test"){

        suppressWarnings(te1_BF(x$threshold,unlist(x$result[5])-1,x$analysis_h1$prior , x$analysis_h1$location , x$analysis_h1$scale , x$analysis_h1$dff , x$alternative,x$ROPE))
      }else{
        suppressWarnings(t2e_BF (x$threshold ,x$result[1,5],x$results[1,6]/x$results[1,5],x$analysis_h1$prior , x$analysis_h1$location , x$analysis_h1$scale , x$analysis_h1$dff , x$alternative,x$ROPE))
      }

    }

  }
}


#' @export
#' @method print BFvalue_t

print.BFvalue_t <- function (x,...) {

  # Unicode box-drawing characters
  UL <- "\u2554"; UR <- "\u2557"
  LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")

  digits <- 3

  # Hypothesis symbols
  H0 <- "H\u2080"
  H1 <- "H\u2081"
  delta <- "\u03b4"
  neq <- "\u2260"
  leq <- "\u2264"
  geq <- "\u2265"
  union <- "\u222a"

  # Bayes Factor subscripts
  BF10 <- "BF\u2081\u2080"
  BF01 <- "BF\u2080\u2081"

  # Header
  header_text <- "BAYES FACTOR CALCULATION"
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypotheses
  cat("Hypotheses\n")
  cat(line)

  if (is.null(x$ROPE)) {

    if (x$alternative == "!=") {
      cat("  ", H0, " (Null)         : ", delta, " = 0\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " ", neq, " 0\n\n", sep = "")
    }

    if (x$alternative == ">") {
      cat("  ", H0, " (Null)         : ", delta, " = 0\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " > 0\n\n", sep = "")
    }

    if (x$alternative == "<") {
      cat("  ", H0, " (Null)         : ", delta, " = 0\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " < 0\n\n", sep = "")
    }

  } else {

    if (x$alternative == "!=") {
      cat("  ", H0, " (Null)         : ",
          min(x$ROPE), " ", leq, " ", delta, " ", leq, " ", max(x$ROPE), "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ",
          delta, " < ", min(x$ROPE), " ", union, " ", delta, " > ", max(x$ROPE), "\n\n", sep = "")
    }

    if (x$alternative == ">") {
      cat("  ", H0, " (Null)         : 0 ", leq, " ", delta, " ", leq, " ", x$ROPE, "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " > ", x$ROPE, "\n\n", sep = "")
    }

    if (x$alternative == "<") {
      cat("  ", H0, " (Null)         : ", x$ROPE, " ", leq, " ", delta, " ", leq, " 0\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", delta, " < ", x$ROPE, "\n\n", sep = "")
    }
  }

  cat(line)

  # Analysis prior
  cat("Analysis prior : \n")
  cat(line)

  cat("Analysis prior under ", H0, ":\n", sep = "")
  if (is.null(x$ROPE)) {
    cat("  Point prior (location = 0)\n")
  } else {
    cat("  ", x$analysis_h1$prior,
        " (location = ", x$analysis_h1$location,
        ", scale = ", x$analysis_h1$scale,
        if (!is.null(x$analysis_h1$dff)) paste0(", df = ", x$analysis_h1$dff),
        ")\n", sep = "")
  }

  cat("Analysis prior under ", H1, ":\n", sep = "")
  cat("  ", x$analysis_h1$prior,
      " (location = ", x$analysis_h1$location,
      ", scale = ", x$analysis_h1$scale,
      if (!is.null(x$analysis_h1$dff)) paste0(", df = ", x$analysis_h1$dff),
      ")\n", sep = "")

  cat("\n")
  cat(line)
  cat("Results\n")
  cat(line)

  cat(
    "t(df=", x$df, ") = ", x$tval,
    ", ", BF10, " = ", x$bf10,
    ", ", BF01, " = ", 1 / x$bf10,
    "\n", sep = ""
  )

  if (x$type == "Indepedent-samples t-test (equal variance)") {
    cat("N1 = ", x$N1, ", N2 = ", x$N2, "\n", sep = "")
  }

}




##### correlation #####
#' @export
#' @method print BFpower_r
print.BFpower_r <- function(x,...) {
  # box and line symbols
  UL <- "\u2554"; UR <- "\u2557"
  LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  digits <- 3

  # hypothesis symbols
  H0 <- "H\u2080"
  H1 <- "H\u2081"
  rho <- "\u03c1"
  neq <- "\u2260"
  leq <- "\u2264"
  geq <- "\u2265"
  union <- "\u222a"

  # Header text with centered padding
  header_text <- ifelse(x$mode_bf == 1, "SAMPLE SIZE CALCULATION", "POWER CALCULATION")
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)

  # Print box header
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypotheses
  cat("Hypotheses\n")
  cat(line)

  if (is.null(x$ROPE)) {
    cat("  Null value        : ", rho, "\u2080 = ", x$h0, " \n", sep = "")
    if (x$alternative == "!=") {
      cat("  ", H0, " (Null)         : ", rho, " = ", rho, "\u2080\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " ", neq, " ", rho, "\u2080\n\n", sep = "")
    }
    if (x$alternative == ">") {
      cat("  ", H0, " (Null)         : ", rho, " = ", rho, "\u2080\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " > ", rho, "\u2080\n\n", sep = "")
    }
    if (x$alternative == "<") {
      cat("  ", H0, " (Null)         : ", rho, " = ", rho, "\u2080\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " < ", rho, "\u2080\n\n", sep = "")
    }
  } else {
    if (x$alternative == "!=") {
      cat("  ", H0, " (Null)         : ", x$h0 + min(x$ROPE), " ", leq, " ", rho, " ", leq, " ", x$h0 + max(x$ROPE), "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " < ", x$h0 + min(x$ROPE), " ", union, " ", rho, " > ", x$h0 + max(x$ROPE), "\n\n", sep = "")
    }
    if (x$alternative == ">") {
      cat("  ", H0, " (Null)         : ", x$h0, " ", leq, " ", rho, " ", leq, " ", x$h0 + x$ROPE, "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " > ", x$h0 + x$ROPE, "\n\n", sep = "")
    }
    if (x$alternative == "<") {
      cat("  ", H0, " (Null)         : ", x$h0 + x$ROPE, " ", leq, " ", rho, " ", leq, " ", x$h0, "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " < ", x$h0 + x$ROPE, "\n\n", sep = "")
    }
  }

  cat(line)
  cat("Analysis/design prior:\n")
  cat(line)

  # prior under H0
  cat("Analysis/design prior under ", H0, ":\n", sep = "")
  if (is.null(x$ROPE)) {
    cat("  Point prior (location = ", x$h0, ")\n", sep = "")
  } else {
    if (x$analysis_h1$prior == "NLP") cat("  Normal-moment (location =", x$h0, ", scale =", x$analysis_h1$scale, ")\n")
    if (x$analysis_h1$prior == "d_beta") cat("  Default Stretched Beta (k =", x$analysis_h1$k, ")\n")
    if (x$analysis_h1$prior == "beta") cat("  Stretched Beta (alpha =", x$analysis_h1$alpha, ", beta =", x$analysis_h1$beta, ")\n")
  }

  # prior under H1
  cat("Analysis/design prior under ", H1, ":\n", sep = "")
  if (is.null(x$design_h1$prior)) {
    if (x$analysis_h1$prior == "NLP") cat("  Normal-moment (location =", x$h0, ", scale =", x$analysis_h1$scale, ")\n")
    if (x$analysis_h1$prior == "d_beta") cat("  Default Stretched Beta (k =", x$analysis_h1$k, ")\n")
    if (x$analysis_h1$prior == "beta") cat("  Stretched Beta (alpha =", x$analysis_h1$alpha, ", beta =", x$analysis_h1$beta, ")\n")
  } else {
    if (x$analysis_h1$prior == "NLP") cat("  Normal-moment (location =", x$h0, ", scale =", x$analysis_h1$scale, ")\n")
    if (x$analysis_h1$prior == "d_beta") cat("  Default Stretched Beta (k =", x$analysis_h1$k, ")\n")
    if (x$analysis_h1$prior == "beta") cat("  Stretched Beta (alpha =", x$analysis_h1$alpha, ", beta =", x$analysis_h1$beta, ")\n")
    cat("Design prior under ", H1, ":\n", sep = "")
    if (x$design_h1$prior == "NLP") cat("  Normal-moment (location =", x$design_h1$location, ", scale =", x$analysis_h1$scale, ")\n")
    if (x$design_h1$prior == "d_beta") cat("  Default Stretched Beta (k =", x$design_h1$k, ")\n")
    if (x$design_h1$prior == "beta") cat("  Stretched Beta (alpha =", x$design_h1$alpha, ", beta =", x$design_h1$beta, ")\n")
    if (x$design_h1$prior == "Point") cat("  Point (location =", x$design_h1$location, ")\n")
  }

  # Threshold
  cat("\n")
  cat(line)
  cat("Threshold of compelling evidence: ", x$D, "\n\n", sep = "")
  cat(line)

  # Results
  cat("Results:\n")
  cat(line)
  cat("  True Positive Rate           = ", round(unlist(x$results[1]), digits = digits), "\n", sep = "")
  cat("  False Positive Rate          = ", round(unlist(x$results[4]), digits = digits), "\n", sep = "")
  cat("  True Negative Rate           = ", round(unlist(x$results[3]), digits = digits), "\n", sep = "")
  cat("  False Negative Rate          = ", round(unlist(x$results[2]), digits = digits), "\n", sep = "")
  cat("  Required Sample size         = ", round(unlist(x$results[5]), digits = 1), "  <---\n", sep = "")
}

#' @export
#' @method plot BFpower_r
plot.BFpower_r<-function (x,...){

  if (is.null(x$design_h1$prior)){
    de_an_prior=1

  }else{
    de_an_prior=0


  }
  if (is.null(x$ROPE)){
    r_prior_plot(x$analysis_h1$k, x$analysis_h1$alpha, x$analysis_h1$beta,
                 x$h0,x$h0,x$analysis_h1$scale,
                 1,x$analysis_h1$prior,de_an_prior,
                 x$design_h1$k, x$design_h1$alpha, x$design_h1$beta,
                 x$design_h1$location,x$design_h1$scale,1,
                 x$design_h1$prior,x$alternative)
  }else {

    re_prior_plot(x$analysis_h1$k, x$analysis_h1$alpha, x$analysis_h1$beta,
                  x$h0,x$h0,x$analysis_h1$scale,
                  1,x$analysis_h1$prior,de_an_prior,
                  x$design_h1$k, x$design_h1$alpha, x$design_h1$beta,
                  x$design_h1$location,x$design_h1$scale,1,
                  x$design_h1$prior,x$alternative,x$ROPE)



  }

  if(x$plot_power){
    if (is.null(x$ROPE)){
      Power_r(x$threshold,x$analysis_h1$k, x$analysis_h1$alpha,
              x$analysis_h1$beta,x$h0,x$alternative,
              x$h0,x$analysis_h1$scale,x$analysis_h1$dff,x$analysis_h1$prior,
              x$design_h1$k, x$design_h1$alpha, x$design_h1$beta,
              x$design_h1$location,x$design_h1$scale,1,x$design_h1$prior,
              de_an_prior,x$results[1,5])
    }else {
      Power_re(x$threshold,x$analysis_h1$k, x$analysis_h1$alpha,
               x$analysis_h1$beta,x$h0,x$alternative,
               x$h0,x$analysis_h1$scale,x$analysis_h1$dff,x$analysis_h1$prior,
               x$design_h1$k, x$design_h1$alpha, x$design_h1$beta,
               x$design_h1$location,x$design_h1$scale,1,x$design_h1$prior,
               de_an_prior,x$results[1,5],x$ROPE)
    }
  }


  if(x$plot_rel){
    if (is.null(x$ROPE)){
      r_bf10_p(x$threshold,x$results[1,5],x$analysis_h1$k,x$analysis_h1$alpha, x$analysis_h1$beta,
               x$h0,
               x$alternative,x$h0,x$analysis_h1$scale,1,x$analysis_h1$prior)
    }else{
      re_bf10_p(x$threshold,x$results[1,5],x$analysis_h1$k,x$analysis_h1$alpha, x$analysis_h1$beta,
                x$h0,
                x$alternative,x$h0,x$analysis_h1$scale,1,x$analysis_h1$prior,x$ROPE)
    }

  }
}

#' @export
#' @method print BFvalue_r
print.BFvalue_r <- function(x,...) {
  # box and line symbols
  UL <- "\u2554"; UR <- "\u2557"
  LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  digits <- 3

  # hypothesis symbols
  H0 <- "H\u2080"
  H1 <- "H\u2081"
  rho <- "\u03c1"
  neq <- "\u2260"
  union <- "\u222a"

  # Centered header
  header_text <- "BAYES FACTOR CALCULATION"
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)

  # Print box header
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypotheses
  cat("Hypotheses\n")
  cat(line)

  if (is.null(x$ROPE)) {
    cat("  Null value        : ", rho, "\u2080 = ", x$h0, "\n", sep = "")
    if (x$alternative == "!=") {
      cat("  ", H0, " (Null)         : ", rho, " = ", rho, "\u2080\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " ", neq, " ", rho, "\u2080\n\n", sep = "")
    }
    if (x$alternative == ">") {
      cat("  ", H0, " (Null)         : ", rho, " = ", rho, "\u2080\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " > ", rho, "\u2080\n\n", sep = "")
    }
    if (x$alternative == "<") {
      cat("  ", H0, " (Null)         : ", rho, " = ", rho, "\u2080\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " < ", rho, "\u2080\n\n", sep = "")
    }
  } else {
    if (x$alternative == "!=") {
      cat("  ", H0, " (Null)         : ", x$h0 + min(x$ROPE), " <= ", rho, " <= ", x$h0 + max(x$ROPE), "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " < ", x$h0 + min(x$ROPE), " ", union, " ", rho, " > ", x$h0 + max(x$ROPE), "\n\n", sep = "")
    }
    if (x$alternative == ">") {
      cat("  ", H0, " (Null)         : ", x$h0, " <= ", rho, " <= ", x$h0 + x$ROPE, "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " > ", x$h0 + x$ROPE, "\n\n", sep = "")
    }
    if (x$alternative == "<") {
      cat("  ", H0, " (Null)         : ", x$h0 + x$ROPE, " <= ", rho, " <= ", x$h0, "\n", sep = "")
      cat("  ", H1, " (Alternative)  : ", rho, " < ", x$h0 + x$ROPE, "\n\n", sep = "")
    }
  }

  cat(line)
  cat("Analysis prior:\n")
  cat(line)

  # prior under H0
  cat("Analysis under ", H0, ":\n", sep = "")
  if (is.null(x$ROPE)) {
    cat("  Point prior (location = ", x$h0, ")\n", sep = "")
  } else {
    if (x$analysis_h1$prior == "Moment") cat("  Normal-moment (location =", x$h0, ", scale =", x$analysis_h1$scale, ")\n")
    if (x$analysis_h1$prior == "d_beta") cat("  Default Stretched Beta (k =", x$analysis_h1$k, ")\n")
    if (x$analysis_h1$prior == "beta") cat("  Stretched Beta (alpha =", x$analysis_h1$alpha, ", beta =", x$analysis_h1$beta, ")\n")
  }

  # prior under H1
  cat("Analysis prior under ", H1, ":\n", sep = "")
  if (x$analysis_h1$prior == "Moment") cat("  Normal-moment (location =", x$h0, ", scale =", x$analysis_h1$scale, ")\n")
  if (x$analysis_h1$prior == "d_beta") cat("  Default Stretched Beta (k =", x$analysis_h1$k, ")\n")
  if (x$analysis_h1$prior == "beta") cat("  Stretched Beta (alpha =", x$analysis_h1$alpha, ", beta =", x$analysis_h1$beta, ")\n")

  # Threshold
  cat("\n")

  # Results
  cat(line)
  cat("Results\n")
  cat(line)
  cat("  r (n = ", x$n, ") = ", x$r, ", BF\u2081\u2080 = ", x$bf10, ", BF\u2080\u2081 = ", 1 / x$bf10, "\n", sep = "")
}



##### ANOVA/regression #####

#' @export
#' @method print BFpower_f
print.BFpower_f <- function(x,...) {
  # Box and line symbols
  UL <- "\u2554"; UR <- "\u2557"
  LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  digits <- 3

  # Header
  header_text <- ifelse(x$mode_bf == 1, "SAMPLE SIZE CALCULATION", "POWER CALCULATION")
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)

  # Print box header
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypotheses
  f2 <- "f\u00B2"  # Unicode superscript 2
  cat("Hypotheses\n")
  cat(line)
  if (is.null(x$ROPE)) {
    cat("  H\u2080 (Null)         : ", f2, " = 0\n", sep = "")
    cat("  H\u2081 (Alternative)  : ", f2, " > 0\n\n", sep = "")
  } else {
    cat("  H\u2080 (Null)         : 0 <= ", f2, " <= ", x$ROPE, "\n", sep = "")
    cat("  H\u2081 (Alternative)  : ", f2, " > ", x$ROPE, "\n\n", sep = "")
  }

  cat(line)
  cat("Number of predictors:\n")
  cat(line)
  cat("  p (reduced model) : ", x$p, "\n", sep = "")
  cat("  k (full model)    : ", x$k, "\n\n", sep = "")
  cat(line)

  # Analysis/design prior
  cat("Analysis/design prior:\n")
  cat(line)

  # Prior under H0
  cat("Analysis/design prior under H\u2080:\n")
  if (is.null(x$ROPE)) {
    cat("  Point prior (location = 0)\n")
  } else {
    if (x$analysis_h1$prior == "effectsize") {
      cat("  Effect size prior (scale =", x$analysis_h1$rscale,
          ", ", f2, " = ", x$analysis_h1$f_m, ", df =", x$analysis_h1$dff, ")\n")
    }
    if (x$analysis_h1$prior == "Moment") {
      cat("  Moment prior (", f2, " =", x$analysis_h1$f_m, ", df =", x$analysis_h1$dff, ")\n")
    }
  }

  # Prior under h1
  cat("Analysis/design prior under H\u2081:\n")
  if (x$analysis_h1$prior == "effectsize") {
    cat("  Effect size prior (scale =", x$analysis_h1$rscale,
        ", ", f2, " =", x$analysis_h1$f_m, ", df =", x$analysis_h1$dff, ")\n")
  }
  if (x$analysis_h1$prior == "Moment") {
    cat("  Moment prior (", f2, " =", x$analysis_h1$f_m, ", df =", x$analysis_h1$dff, ")\n")
  }

  # Design prior under h1 (if available)
  if (!is.null(x$design_h1$prior)) {
    cat("Design prior under H\u2081:\n")
    if (x$design_h1$prior == "effectsize") {
      cat("  Effect size prior (scale =", x$design_h1$rscale,
          ", ", f2, " =", x$design_h1$f_m, ", df =", x$design_h1$dff, ")\n")
    }
    if (x$design_h1$prior == "Moment") {
      cat("  Moment prior (", f2, " =", x$design_h1$f_m, ", df =", x$design_h1$dff, ")\n")
    }
    if (x$design_h1$prior == "Point") {
      cat("  Point prior (", f2, " =", x$design_h1$f_m, ")\n")
    }
  }

  # Threshold
  cat("\n")
  cat(line)
  cat("Threshold of compelling evidence: ", x$D, "\n\n", sep = "")
  cat(line)

  # Results
  cat("Results:\n")
  cat(line)
  cat("  True Positive Rate   = ", round(unlist(x$results[1]), digits = digits), "\n", sep = "")
  cat("  False Positive Rate  = ", round(unlist(x$results[4]), digits = digits), "\n", sep = "")
  cat("  True Negative Rate   = ", round(unlist(x$results[3]), digits = digits), "\n", sep = "")
  cat("  False Negative Rate  = ", round(unlist(x$results[2]), digits = digits), "\n", sep = "")
  cat("  Required Sample size = ", round(unlist(x$results[5]), 1), "  <---\n", sep = "")
}



#' @export
#' @method plot BFpower_f
plot.BFpower_f<-function (x,...){
  q= x$k-x$p
  if (is.null(x$design_h1$prior)){
    de_an_prior=1

  }else{
    de_an_prior=0

  }


  if (is.null(x$ROPE)){

    suppressWarnings(prior_plot_f(q,x$analysis_h1$dff,x$analysis_h1$rscale,x$analysis_h1$f_m,
                                  x$analysis_h1$prior,
                                  x$design_h1$dff,x$design_h1$rscale,x$design_h1$f_m,
                                  x$design_h1$prior,de_an_prior))

  }else {

    suppressWarnings(prior_plot_fe(q,x$analysis_h1$dff,x$analysis_h1$rscale,x$analysis_h1$f_m,
                                   x$analysis_h1$prior,
                                   x$design_h1$dff,x$design_h1$rscale,x$design_h1$f_m,
                                   x$design_h1$prior,de_an_prior,x$ROPE))

  }

  if(x$plot_power){
    if (is.null(x$ROPE)){

      suppressWarnings( Power_f(x$threshold,x$k,x$p,
                                x$analysis_h1$dff,x$analysis_h1$rscale,x$analysis_h1$f_m,x$analysis_h1$prior,
                                1,1,x$design_h1$dff,x$design_h1$rscale,x$design_h1$f_m
                                ,x$design_h1$prior,de_an_prior,unlist(x$result[1,5])))


    }else {

      suppressWarnings( Power_fe(x$threshold,x$k,x$p,
                                 x$analysis_h1$dff,x$analysis_h1$rscale,x$analysis_h1$f_m,x$analysis_h1$prior,
                                 1,1,x$design_h1$dff,x$design_h1$rscale,x$design_h1$f_m
                                 ,x$design_h1$prior,de_an_prior,unlist(x$result[1,5]),x$ROPE))
    }
  }


  if(x$plot_rel){
    if (is.null(x$ROPE)){
      suppressWarnings(bf10_f(x$threshold,unlist(x$result[1,5]),x$k,x$p,x$analysis_h1$dff,
                              x$analysis_h1$rscale,x$analysis_h1$f_m,x$analysis_h1$prior))
    }else{
      suppressWarnings(bf10_fe(x$threshold,unlist(x$result[1,5]),x$k,x$p,x$analysis_h1$dff,
                               x$analysis_h1$rscale,x$analysis_h1$f_m,x$analysis_h1$prior,x$ROPE))
    }

  }
}




#' @export
#' @method print BFvalue_f
print.BFvalue_f <- function(x,...) {
  # Box and line symbols
  UL <- "\u2554"; UR <- "\u2557"
  LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  digits <- 3

  # Unicode symbols
  sub1 <- "\u2081"  # subscript 1
  sub2 <- "\u2082"  # subscript 2
  sub0 <- "\u2080"  # subscript 0
  f2 <- "\u00B2"   # superscript 2

  # Bayes Factor subscripts
  BF10 <- paste0("BF", sub1, sub0)
  BF01 <- paste0("BF", sub0, sub1)

  # Header
  header_text <- "BAYES FACTOR CALCULATION"
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypotheses
  cat("Hypotheses\n")
  cat(line)
  if (is.null(x$ROPE)) {
    cat("  H", sub0, " (Null)         : f", f2, " = 0\n", sep = "")
    cat("  H", sub1, " (Alternative)  : f", f2, " > 0\n\n", sep = "")
  } else {
    cat("  H", sub0, " (Null)         : 0 <= f", f2, " <= ", x$ROPE, "\n", sep = "")
    cat("  H", sub1, " (Alternative)  : f", f2, " > ", x$ROPE, "\n\n", sep = "")
  }

  cat(line)
  cat("Analysis prior:\n")
  cat(line)

  # Prior under H0
  cat("Analysis prior under H", sub0, ":\n", sep = "")
  if (is.null(x$ROPE)) {
    cat("  Point prior (location = 0)\n")
  } else {
    if (x$analysis_h1$prior == "effectsize") {
      cat("  Effect size prior (scale =", x$analysis_h1$rscale,
          ", f =", x$analysis_h1$f_m, ", df =", x$analysis_h1$dff, ")\n")
    }
    if (x$analysis_h1$prior == "Moment") {
      cat("  Moment prior (f =", x$analysis_h1$f_m, ", df =", x$analysis_h1$dff, ")\n")
    }
  }

  # Prior under H1
  cat("Analysis prior under H", sub1, ":\n", sep = "")
  if (x$analysis_h1$prior == "effectsize") {
    cat("  Effect size prior (scale =", x$analysis_h1$rscale,
        ", f =", x$analysis_h1$f_m, ", df =", x$analysis_h1$dff, ")\n")
  }
  if (x$analysis_h1$prior == "Moment") {
    cat("  Moment prior (f =", x$analysis_h1$f_m, ", df =", x$analysis_h1$dff, ")\n")
  }

  # Threshold
  cat("\n")
  cat(line)

  # Results
  cat("Results\n")
  cat(line)
  cat("  F(df", sub1, " =", x$df1, ", df", sub2, " =", x$df2, ") = ", x$fval,
      ", ", BF10, " = ", x$bf10, ", ", BF01, " = ", 1 / x$bf10, "\n", sep = "")
}


#### one-proportion ####
#' @export
#' @method print BFpower_bin
print.BFpower_bin <- function(x,...) {
  # Box and line symbols
  UL <- "\u2554"; UR <- "\u2557"
  LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  digits <- 3

  # Unicode symbols
  sub0 <- "\u2080"  # subscript 0
  sub1 <- "\u2081"  # subscript 1
  union <- "\u222a" #
  neq <- "\u2260"   #
  theta <- "\u03b8" #

  # Header
  header_text <- ifelse(x$mode_bf == 1, "SAMPLE SIZE CALCULATION", "POWER CALCULATION")
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypotheses
  cat("Hypotheses\n")
  cat(line)

  if (is.null(x$ROPE)) {
    cat("  Null value        : ", theta, sub0, " = ", x$h0, "\n", sep = "")
    if (x$alternative == "!=") {
      cat("  H", sub0, " (Null)         : ", theta, " = ", theta, sub0, "\n", sep = "")
      cat("  H", sub1, " (Alternative)  : ", theta, " ", neq, " ", theta, sub0, "\n\n", sep = "")
    }
    if (x$alternative == ">") {
      cat("  H", sub0, " (Null)         : ", theta, " = ", theta, sub0, "\n", sep = "")
      cat("  H", sub1, " (Alternative)  : ", theta, " > ", theta, sub0, "\n\n", sep = "")
    }
    if (x$alternative == "<") {
      cat("  H", sub0, " (Null)         : ", theta, " = ", theta, sub0, "\n", sep = "")
      cat("  H", sub1, " (Alternative)  : ", theta, " < ", theta, sub0, "\n\n", sep = "")
    }
  } else {
    if (x$alternative == "!=") {
      cat("  H", sub0, " (Null)         : ", x$h0 + min(x$ROPE), " <= ", theta, " <= ", x$h0 + max(x$ROPE), "\n", sep = "")
      cat("  H", sub1, " (Alternative)  : ", theta, " < ", x$h0 + min(x$ROPE), " ", union, " ", theta, " > ", x$h0 + max(x$ROPE), "\n\n", sep = "")
    }
    if (x$alternative == ">") {
      cat("  H", sub0, " (Null)         : ", x$h0, " <= ", theta, " <= ", x$h0 + x$ROPE, "\n", sep = "")
      cat("  H", sub1, " (Alternative)  : ", theta, " > ", x$h0 + x$ROPE, "\n\n", sep = "")
    }
    if (x$alternative == "<") {
      cat("  H", sub0, " (Null)         : ", x$h0 + x$ROPE, " <= ", theta, " <= ", x$h0, "\n", sep = "")
      cat("  H", sub1, " (Alternative)  : ", theta, " < ", x$h0 + x$ROPE, "\n\n", sep = "")
    }
  }

  cat(line)
  cat("Analysis/design prior:\n")
  cat(line)

  # Prior under H0
  cat("Analysis/design prior under H", sub0, ":\n", sep = "")
  if (is.null(x$ROPE)) {
    cat("  Point prior (location = ", x$h0, ")\n", sep = "")
  } else {
    if (x$analysis_h1$prior == "Moment") {
      cat("  Normal-moment (location =", x$h0, ", scale =", x$analysis_h1$scale, ")\n")
    }
    if (x$analysis_h1$prior == "beta") {
      cat("  Stretched Beta (alpha =", x$analysis_h1$alpha, ", beta =", x$analysis_h1$beta, ")\n")
    }
  }

  # Prior under H1
  cat("Analysis/design prior under H", sub1, ":\n", sep = "")
  if (x$analysis_h1$prior == "Moment") {
    cat("  Normal-moment (location =", x$h0, ", scale =", x$analysis_h1$scale, ")\n")
  }
  if (x$analysis_h1$prior == "beta") {
    cat("  Beta (alpha =", x$analysis_h1$alpha, ", beta =", x$analysis_h1$beta, ")\n")
  }

  # Design prior under H1
  if (!is.null(x$design_h1$prior)) {
    cat("Design prior under H", sub1, ":\n", sep = "")
    if (x$design_h1$prior == "Moment") {
      cat("  Normal-moment (location =", x$design_h1$location, ", scale =", x$analysis_h1$scale, ")\n")
    }
    if (x$design_h1$prior == "beta") {
      cat("  Beta (alpha =", x$design_h1$alpha, ", beta =", x$design_h1$beta, ")\n")
    }
    if (x$design_h1$prior == "Point") {
      cat("  Point (location =", x$design_h1$location, ")\n")
    }
  }

  # Threshold
  cat("\n")
  cat(line)
  cat("Threshold of compelling evidence: ", x$threshold, "\n\n", sep = "")
  cat(line)

  # Results
  cat("Results:\n")
  cat(line)
  cat("  True Positive Rate  = ", round(unlist(x$results[1]), digits), "\n", sep = "")
  cat("  False Positive Rate = ", round(unlist(x$results[4]), digits), "\n", sep = "")
  cat("  True Negative Rate  = ", round(unlist(x$results[3]), digits), "\n", sep = "")
  cat("  False Negative Rate = ", round(unlist(x$results[2]), digits), "\n", sep = "")
  cat("  Required Sample size= ", round(unlist(x$results[5]), 1), "  <---\n", sep = "")
}





#' @export
#' @method plot BFpower_bin
plot.BFpower_bin<-function (x,...){

  if (is.null(x$design_h1$prior)){
    de_an_prior=1

  }else{
    de_an_prior=0


  }
  if (is.null(x$ROPE)){
    bin_prior_plot(x$h0,x$analysis_h1$alpha,x$analysis_h1$beta,x$h0,x$analysis_h1$scale,x$analysis_h1$prior,
                   x$design_h1$alpha,x$design_h1$beta,x$design_h1$location,
                   x$design_h1$scale,x$design_h1$prior,
                   x$alternative,de_an_prior)
  }else {

    bin_e_prior_plot(x$h0,x$analysis_h1$alpha,x$analysis_h1$beta,x$h0,x$analysis_h1$scale,x$analysis_h1$prior,
                     x$design_h1$alpha,x$design_h1$beta,x$design_h1$location,
                     x$design_h1$scale,x$design_h1$prior,
                     x$alternative,de_an_prior,x$ROPE)
  }

  if(x$plot_power){
    if (is.null(x$ROPE)){
      Power_bin(x$threshold,x$h0,x$analysis_h1$alpha,x$analysis_h1$beta,
                x$h0,x$analysis_h1$scale,x$analysis_h1$prior,
                x$alternative,
                x$design_h1$alpha_d,x$design_h1$beta_d,x$design_h1$location_d,
                x$design_h1$scale_d,x$design_h1$prior_d, de_an_prior,
                unlist(x$results[1,5]))
    }else {
      Power_e_bin(x$threshold,x$h0,x$analysis_h1$alpha,x$analysis_h1$beta,
                  x$h0,x$analysis_h1$scale,x$analysis_h1$prior,
                  x$alternative,
                  x$design_h1$alpha_d,x$design_h1$beta_d,x$design_h1$location_d,
                  x$design_h1$scale_d,x$design_h1$prior_d, de_an_prior,
                  unlist(x$results[1,5]),x$ROPE)
    }
  }


  if(x$plot_rel){
    if (is.null(x$ROPE)){
      bin_bf10(x$threshold,unlist(x$results[1,5]),x$analysis_h1$alpha,
               x$analysis_h1$beta,x$h0,x$analysis_h1$scale,x$analysis_h1$prior,
               x$alternative)
    }else{
      bin_e_bf10(x$threshold,unlist(x$results[1,5]),x$analysis_h1$alpha,
                 x$analysis_h1$beta,x$h0,x$analysis_h1$scale,x$analysis_h1$prior,
                 x$alternative,x$ROPE)
    }

  }
}

#' @export
#' @method print BFvalue_bin
print.BFvalue_bin <- function(x,...) {
  # Box and line symbols
  UL <- "\u2554"; UR <- "\u2557"
  LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  digits <- 3

  # Unicode symbols
  sub0 <- "\u2080"
  sub1 <- "\u2081"
  theta <- "\u03b8"
  union <- "\u222a"
  neq <- "\u2260"

  # Header
  header_text <- "BAYES FACTOR CALCULATION"
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypotheses
  cat("Hypotheses\n")
  cat(line)

  if (is.null(x$ROPE)) {
    cat("  Null value        : ", theta, sub0, " = ", x$h0, "\n", sep = "")
    if (x$alternative == "!=") cat("  H", sub0, " (Null)         : ", theta, " = ", theta, sub0, "\n  H", sub1, " (Alternative)  : ", theta, " ", neq, " ", theta, sub0, "\n\n", sep = "")
    if (x$alternative == ">")  cat("  H", sub0, " (Null)         : ", theta, " = ", theta, sub0, "\n  H", sub1, " (Alternative)  : ", theta, " > ", theta, sub0, "\n\n", sep = "")
    if (x$alternative == "<")  cat("  H", sub0, " (Null)         : ", theta, " = ", theta, sub0, "\n  H", sub1, " (Alternative)  : ", theta, " < ", theta, sub0, "\n\n", sep = "")
  } else {
    if (x$alternative == "!=") {
      cat("  H", sub0, " (Null)         : ", x$h0 + min(x$ROPE), " <= ", theta, " <= ", x$h0 + max(x$ROPE), "\n", sep = "")
      cat("  H", sub1, " (Alternative)  : ", theta, " < ", x$h0 + min(x$ROPE), " ", union, " ", theta, " > ", x$h0 + max(x$ROPE), "\n\n", sep = "")
    }
    if (x$alternative == ">")  cat("  H", sub0, " (Null)         : ", x$h0, " <= ", theta, " <= ", x$h0 + x$ROPE, "\n  H", sub1, " (Alternative)  : ", theta, " > ", x$h0 + x$ROPE, "\n\n", sep = "")
    if (x$alternative == "<")  cat("  H", sub0, " (Null)         : ", x$h0 + x$ROPE, " <= ", theta, " <= ", x$h0, "\n  H", sub1, " (Alternative)  : ", theta, " < ", x$h0 + x$ROPE, "\n\n", sep = "")
  }

  cat(line)
  cat("Analysis prior:\n")
  cat(line)

  # Prior under H0
  cat("Analysis prior under H", sub0, ":\n", sep = "")
  if (is.null(x$ROPE)) {
    cat("  Point prior (location = ", x$h0, ")\n", sep = "")
  } else {
    if (x$analysis_h1$prior == "Moment") cat("  Normal-moment (location =", x$h0, ", scale =", x$analysis_h1$scale, ")\n")
    if (x$analysis_h1$prior == "beta")   cat("  Stretched Beta (alpha =", x$analysis_h1$alpha, ", beta =", x$analysis_h1$beta, ")\n")
  }

  # Prior under H1
  cat("Analysis prior under H", sub1, ":\n", sep = "")
  if (x$analysis_h1$prior == "Moment") cat("  Normal-moment (location =", x$h0, ", scale =", x$analysis_h1$scale, ")\n")
  if (x$analysis_h1$prior == "beta")   cat("  Beta (alpha =", x$analysis_h1$alpha, ", beta =", x$analysis_h1$beta, ")\n")

  # Results
  cat(line)
  cat("Results\n")
  cat(line)
  cat("  n = ", x$n, ", x = ", x$x, ", BF", sub1, sub0, " = ", x$bf10, ", BF", sub0, sub1, " = ", 1 / x$bf10, "\n", sep = "")
}


#### two-proportions ####
#' @export
#' @method print BFpower_2p
print.BFpower_2p <- function(x,...) {
  # Box and line symbols
  UL <- "\u2554"; UR <- "\u2557"; LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  digits <- 3

  # Unicode symbols
  theta <- "\u03B8"
  sub0 <- "\u2080"; sub1 <- "\u2081"; sub2 <- "\u2082"
  neq <- "\u2260"

  theta0 <- paste0(theta, sub0)
  theta1 <- paste0(theta, sub1)
  theta2 <- paste0(theta, sub2)

  # Header
  header_text <- ifelse(x$mode_bf == 1, "SAMPLE SIZE CALCULATION", "POWER CALCULATION")
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypotheses
  cat("Hypotheses\n")
  cat(line)
  cat("  H", sub0, " (Null)         : ", theta1, " = ", theta2, " = ", theta0, "\n", sep = "")
  cat("  H", sub1, " (Alternative)  : ", theta1, " ", neq, " ", theta2, "\n\n", sep = "")

  cat(line)
  cat("Analysis/design prior:\n")
  cat(line)

  # Prior under H0
  cat("Analysis/design prior under H", sub0, "\n", sep = "")
  cat(" ", theta0, " ~ Beta( alpha =", x$analysis_h0$a,
      ", beta = ", x$analysis_h0$b, ")\n", sep = "")

  # Prior under H1
  cat("Analysis prior under H", sub1, "\n", sep = "")
  cat(" ", theta1, " ~ Beta( alpha =", x$analysis_h1_theta_1$a,
      ", beta = ", x$analysis_h1_theta_1$b, ")\n", sep = "")
  cat(" ", theta2, " ~ Beta( alpha =", x$analysis_h1_theta_2$a,
      ", beta = ", x$analysis_h1_theta_2$b, ")\n\n", sep = "")

  # Design prior under H1
  cat("Design prior under H", sub1, "\n", sep = "")
  if (x$design_h1_theta_1$prior == "same") {
    cat(" ", theta1, " ~ Beta( alpha =", x$analysis_h1_theta_1$a,
        ", beta = ", x$analysis_h1_theta_1$b, ")\n", sep = "")
    cat(" ", theta2, " ~ Beta( alpha =", x$analysis_h1_theta_2$a,
        ", beta = ", x$analysis_h1_theta_2$b, ")\n", sep = "")
  } else if (x$design_h1_theta_1$prior == "beta") {
    cat(" ", theta1, " ~ Beta( alpha =", x$design_h1_theta_1$a,
        ", beta = ", x$design_h1_theta_1$b, ")\n", sep = "")
    cat(" ", theta2, " ~ Beta( alpha =", x$design_h1_theta_2$a,
        ", beta = ", x$design_h1_theta_2$b, ")\n", sep = "")
  } else if (x$design_h1_theta_1$prior == "Point") {
    cat(" ", theta1, " = ", x$design_h1_theta_1$p, "\n", sep = "")
    cat(" ", theta2, " = ", x$design_h1_theta_2$p, "\n", sep = "")
  }

  # Threshold
  cat("\n")
  cat(line)
  cat("Threshold of compelling evidence: ", x$ROPE, "\n\n", sep = "")
  cat(line)

  # Results
  cat("Results:\n")
  cat(line)
  cat(paste0("  True Positive Rate                    = ", round(x$results[1], digits = digits), "\n"))
  cat(paste0("  False Positive Rate                   = ", round(x$results[4], digits = digits), "\n"))
  cat(paste0("  True Negative Rate                    = ", round(x$results[3], digits = digits), "\n"))
  cat(paste0("  False Negative Rate                   = ", round(x$results[2], digits = digits), "\n"))
  cat(paste0("  Required Sample size per group        = ", round(x$results[5], digits = 1), "  <---\n"))
}





#' @export
#' @method plot BFpower_2p
plot.BFpower_2p<-function (x,...){

  p2_prior_plot(x$analysis_h0$a, x$analysis_h0$b, 1, 1, 0, "same", 0)
  p2_prior_plot(x$analysis_h1_theta_1$a, x$analysis_h1_theta_1$b,
                x$design_h1_theta_1$a, x$design_h1_theta_1$b,
                x$design_h1_theta_1$p, x$design_h1_theta_1$prior, 1)
  p2_prior_plot(x$analysis_h1_theta_2$a, x$analysis_h1_theta_2$b,
                x$design_h1_theta_2$a, x$design_h1_theta_2$b,
                x$design_h1_theta_2$p, x$design_h1_theta_2$prior, 2)


  if(x$plot_power){
    Power_p2(
      x$threshold, unlist(x$results[1,5]), x$analysis_h0$a, x$analysis_h0$b,
      x$analysis_h1_theta_1$a, x$analysis_h1_theta_1$b,
      x$analysis_h1_theta_2$a,x$analysis_h1_theta_2$b,
      unlist(x$results[1,6]) / unlist(x$results[1,5]),
      x$design_h1_theta_1$prior, x$design_h1_theta_1$a, x$design_h1_theta_1$b,
      x$design_h1_theta_1$p,
      x$design_h1_theta_2$prior, x$design_h1_theta_2$a, x$design_h1_theta_2$b,
      x$design_h1_theta_2$p
    )
  }


  if(x$plot_rel){


    print(heatmap_p2(x$grid, x$threshold))
  }
}

#' @export
#' @method print BFvalue_2p
print.BFvalue_2p <- function(x,...) {
  # Box and line symbols
  UL <- "\u2554"; UR <- "\u2557"; LL <- "\u255a"; LR <- "\u255d"
  HL <- "\u2550"; VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  digits <- 3

  # Unicode symbols
  theta <- "\u03B8"
  sub0 <- "\u2080"; sub1 <- "\u2081"; sub2 <- "\u2082"
  neq <- "\u2260"

  theta0 <- paste0(theta, sub0)
  theta1 <- paste0(theta, sub1)
  theta2 <- paste0(theta, sub2)

  # Header
  header_text <- "BAYES FACTOR CALCULATION"
  pad_total <- 50 - nchar(header_text)
  pad_left <- floor(pad_total / 2)
  pad_right <- ceiling(pad_total / 2)

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, strrep(" ", pad_left), header_text, strrep(" ", pad_right), VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Type of test
  cat(x$type, "\n\n", sep = "")
  cat(line)

  # Hypotheses
  cat("Hypotheses\n")
  cat(line)
  cat("  H", sub0, " (Null)         : ", theta1, " = ", theta2, " = ", theta0, "\n", sep = "")
  cat("  H", sub1, " (Alternative)  : ", theta1, " ", neq, " ", theta2, "\n\n", sep = "")

  cat(line)
  cat("Analysis prior:\n")
  cat(line)

  # Prior under H0
  cat("Analysis prior under H", sub0, "\n", sep = "")
  cat(" ", theta0, " ~ Beta( alpha =", x$analysis_h0$a,
      ", beta =", x$analysis_h0$b, ")\n", sep = "")

  # Prior under H1
  cat("Analysis prior under H", sub1, "\n", sep = "")
  cat(" ", theta1, " ~ Beta( alpha =", x$analysis_h1_theta_1$a,
      ", beta =", x$analysis_h1_theta_1$b, ")\n", sep = "")
  cat(" ", theta2, " ~ Beta( alpha =", x$analysis_h1_theta_2$a,
      ", beta =", x$analysis_h1_theta_2$b, ")\n", sep = "")

  # Results
  cat(line)
  cat("Results:\n")
  cat(line)
  cat("  n", sub1, " = ", x$n1, ", x", sub1, " = ", x$x1, ", n", sub2, " = ", x$n2, ", x", sub2, " = ", x$x2, "\n", sep = "")
  cat("  BF", sub1, sub0, " = ", x$bf10, ", BF", sub0, sub1, " = ", 1 / x$bf10, "\n", sep = "")
}


