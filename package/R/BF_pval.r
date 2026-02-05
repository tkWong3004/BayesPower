########## t-test ##########
t.pval <- function(tval, n1, n2 = NULL, alternative, ROPE = NULL, type = "One-sample t-test") {

  # Degrees of freedom and scaling constant
  if (type == "One-sample t-test") {
    df <- n1 - 1
    constant <- sqrt(n1)
  } else {
    df <- n1 + n2 - 2
    constant <- sqrt(n1 * n2 / (n1 + n2))
  }

  # No ROPE: standard p-values
  if (is.null(ROPE)) {
    p <- switch(alternative,
                "<" = stats::pt(tval, df),
                ">" = stats::pt(tval, df, lower.tail = FALSE),
                "!=" = 2 * stats::pt(abs(tval), df, lower.tail = FALSE))
  } else {
    # ROPE specified
    if (alternative %in% c("<", ">")) {
      # ROPE must be length 1 for one-sided tests
      ncp <- ROPE * constant
      p <- switch(alternative,
                  "<" = stats::pt(tval, df, ncp = ncp, lower.tail = FALSE),
                  ">" = stats::pt(tval, df, ncp = ncp, lower.tail = TRUE))
    } else if (alternative == "!=") {
      # ROPE must be length 2 for two-sided
      if (length(ROPE) != 2) stop("For alternative '!=', ROPE must be of length 2.")
      ncp_lower <- ROPE[1] * constant
      ncp_upper <- ROPE[2] * constant
      p <- max(
        stats::pt(tval, df, ncp = ncp_lower, lower.tail = FALSE),
        stats::pt(tval, df, ncp = ncp_upper, lower.tail = TRUE)
      )
    } else {
      stop("Invalid alternative. Must be '<', '>', or '!='.")
    }
  }

  return(p)
}


########## cor ##########

r.pval <- function(r, n,h0, alternative , ROPE = NULL) {
   Z.obs <- r_mean(r)
   Z.sd  <- r_sd(n)
   Z.h0  <- r_mean(h0)
  if (is.null(ROPE)) {
    p <- switch(alternative,
                "<" = stats::pnorm(Z.obs, mean=Z.h0,sd=Z.sd),
                ">" = stats::pnorm(Z.obs, mean=Z.h0,sd=Z.sd,lower=F),
                "!=" = min(stats::pnorm(Z.obs, mean=Z.h0,sd=Z.sd,lower=F),stats::pnorm(Z.obs, mean=Z.h0,sd=Z.sd,lower=T))*2)
  } else {
    ROPE = h0+ROPE

    Z.h1  <- r_mean(ROPE)

    p <- switch(alternative,
                "<" = stats::pnorm(Z.obs, mean=Z.h1,sd=Z.sd,lower.tail = F),
                ">" = stats::pnorm(Z.obs, mean=Z.h1,sd=Z.sd,lower.tail = T),
                "!=" = max(stats::pnorm(Z.obs, mean=max(Z.h1),sd=Z.sd,lower.tail = T),
                           stats::pnorm(Z.obs, mean=min(Z.h1),sd=Z.sd,lower.tail = F)))
  }

  return(p)
}

########## f.test ##########
f.pval <- function(fval, df1,df2,ROPE=NULL) {
  m=df1+df2
  p<-if(is.null(ROPE))  stats::pf(fval, df1, df2, lower.tail = FALSE) else{
     stats::pf(fval, df1, df2, ncp=m*ROPE,lower.tail = T)
  }

  return(p)
}

########## ONE-proportion ##########
bin.pval <-function(x,n,h0,alternative,ROPE=NULL){
  alternative <- switch(alternative,
                        "!=" = "two.sided",
                        "<"  = "less",
                        ">"  = "greater")
  if (is.null(ROPE)) {
    p <- stats::binom.test(x, n, p = h0,
                    alternative = alternative,
                    conf.level = 0.95)$p.value

  } else {
    ROPE = h0+ROPE

    p <- switch(alternative,
                "less" = stats::binom.test(x, n, p = min(ROPE),
                                 alternative = "greater",
                                 conf.level = 0.95)$p.value,
                "greater" = stats::binom.test(x, n, p = max(ROPE),
                                 alternative = "less",
                                 conf.level = 0.95)$p.value,
                "two.sided" = {
                  max(stats::binom.test(x, n, p = min(ROPE),
                                 alternative = "greater",
                                 conf.level = 0.95)$p.value,
                      stats::binom.test(x, n, p = max(ROPE),
                                 alternative = "less",
                                 conf.level = 0.95)$p.value)

                })
  }
  return(p)
}
########## two-proportion ##########

p2.pval <-function(n1,x1,n2,x2){
  tab <- matrix(
    c(x1, n1 - x1,
      x2, n2 - x2),
    nrow = 2,
    byrow = TRUE
  )
  stats::fisher.test(tab)$p.value

}
