
show_t1_code <- function(x) {

  args <- c(
    "hypothesis","e","D",
    "target","alpha",
    "model","location","scale","dff",
    "model_d","location_d","scale_d","dff_d",
    "direct",
    "pc","rela"
  )
  # NOTE: N and mode_bf handled separately

  code_lines <- sapply(args, function(arg) {

    val <- x[[arg]]

    ## OMITTING RULES
    if (!is.null(x$interval) && x$interval == 1 && arg == "e") return(NULL)

    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
        arg %in% c("model_d","location_d","scale_d","dff_d")) return(NULL)

    if (is.null(val)) return(NULL)

    ## SPECIAL RENAMING
    arg_print <- arg
    if (arg == "target") arg_print <- "true_rate"
    if (arg == "alpha") arg_print <- "false_rate"

    if (arg == "direct") {
      if (val != "h0") return(NULL)
      arg_print <- "type_rate"
      val <- "\"negative\""
    }

    if (arg == "pc") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_power"
      val <- TRUE
    }

    if (arg == "rela") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_rel"
      val <- TRUE
    }


    ## HYPOTHESIS > ALTERNATIVE RULE
    if (arg == "hypothesis") {

      arg_print <- "alternative"

      val <- switch(val,
                    "<"  = "less",
                    "!=" = "two.sided",
                    ">"  = "greater",
                    stop("Invalid hypothesis")
      )

      val <- shQuote(val)
    }

    ## D > threshold
    if (arg == "D") arg_print <- "threshold"

    ## e > ROPE
    if (arg == "e") arg_print <- "ROPE"


    ## model > prior_analysis
    if (arg == "model") arg_print <- "prior_analysis"

    ## model_d > prior_design
    if (arg == "model_d") arg_print <- "prior_design"


    ## VALUE FORMATTING
    if (is.character(val) && !grepl("^\"", val)) val <- shQuote(val)
    else if (is.vector(val) && length(val) > 1)
      val <- paste0("c(", paste(val, collapse = ", "), ")")

    glue::glue("  {arg_print} = {val},")
  })

  code_lines <- code_lines[!sapply(code_lines, is.null)]

  ## -----------------------------------
  ## N LOGIC (the corrected part)
  ## -----------------------------------
  if ( x$mode_bf != 1) {
    # Print N always (even if NULL)
    Nval <- if (is.null(x$N)) "NULL" else x$N
    code_lines <- c(code_lines, glue::glue("  N = {Nval}"))
  }
  # If mode_bf == 1 > skip printing N entirely

  ## Remove trailing comma
  if (length(code_lines) > 0)
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])

  paste0(
    "BFpower.ttest.OneSample(\n",
    paste(code_lines, collapse = "\n"),
    "\n)"
  )
}



show_t2_code <- function(x) {

  args <- c(
    "hypothesis","e","D",
    "target","alpha",
    "model","location","scale","dff",
    "model_d","location_d","scale_d","dff_d",
    "direct",
    "pc","rela",
    "r"    # keep r inside main arg loop
  )

  code_lines <- sapply(args, function(arg) {

    val <- x[[arg]]

    ## -----------------------------
    ## OMITTING RULES
    ## -----------------------------
    if (!is.null(x$interval) && x$interval == 1 && arg == "e") return(NULL)

    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
        arg %in% c("model_d","location_d","scale_d","dff_d")) return(NULL)

    ## NEW RULE:
    ## If mode_bf != 1 > DO NOT print r
    if (x$mode_bf != 1 && arg == "r") return(NULL)

    if (is.null(val)) return(NULL)

    ## -----------------------------
    ## SPECIAL RENAMING
    ## -----------------------------
    arg_print <- arg
    if (arg == "target") arg_print <- "true_rate"
    if (arg == "alpha") arg_print <- "false_rate"

    if (arg == "direct") {
      if (val != "h0") return(NULL)
      arg_print <- "type_rate"
      val <- "\"negative\""
    }

    if (arg == "pc") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_power"
      val <- TRUE
    }

    if (arg == "rela") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_rel"
      val <- TRUE
    }



    ## HYPOTHESIS > ALTERNATIVE RULE
    if (arg == "hypothesis") {

      arg_print <- "alternative"

      val <- switch(val,
                    "<"  = "less",
                    "!=" = "two.sided",
                    ">"  = "greater",
                    stop("Invalid hypothesis")
      )

      val <- shQuote(val)
    }

    ## D > threshold
    if (arg == "D") arg_print <- "threshold"

    ## e > ROPE
    if (arg == "e") arg_print <- "ROPE"


    ## model > prior_analysis
    if (arg == "model") arg_print <- "prior_analysis"

    ## model_d > prior_design
    if (arg == "model_d") arg_print <- "prior_design"


    ## -----------------------------
    ## VALUE FORMATTING
    ## -----------------------------
    if (is.character(val) && !grepl("^\"", val)) {
      val <- shQuote(val)
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")
    }

    glue::glue("  {arg_print} = {val},")
  })

  code_lines <- code_lines[!sapply(code_lines, is.null)]

  ## -----------------------------
  ## N1 / N2 LOGIC (parallel to t1 code)
  ## -----------------------------
  if (x$mode_bf != 1) {

    N1val <- if (is.null(x$N1)) "NULL" else x$N1
    N2val <- if (is.null(x$N2)) "NULL" else x$N2

    code_lines <- c(
      code_lines,
      glue::glue("  N1 = {N1val},"),
      glue::glue("  N2 = {N2val}")
    )
  }

  ## -----------------------------
  ## Remove trailing comma
  ## -----------------------------
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  }

  paste0(
    "BFpower.ttest.TwoSample(\n",
    paste(code_lines, collapse = "\n"),
    "\n)"
  )
}



show_cor_code <- function(x) {

  args <- c(
    "hypothesis","h0","e","D","target","FP",
    "model","k","alpha","beta","scale",
    "model_d","alpha_d","beta_d","location_d","k_d","scale_d","dff_d",
    "N","direct",
    "pc","rela"
  )

  code_lines <- sapply(args, function(arg) {

    val <- x[[arg]]

    ## ---------------------------------------------------------
    ## OMISSION RULES
    ## ---------------------------------------------------------

    # hide e if interval==1 ("no equivalence")
    if (!is.null(x$interval) && x$interval == 1 && arg == "e") return(NULL)

    # hide N unless mode_bf == 0
    if (arg == "N" && x$mode_bf != 0) return(NULL)

    # hide target/FP (true/false rates) when mode_bf==0 (sample-size mode)
    if (x$mode_bf == 0 && arg %in% c("target","FP")) return(NULL)

    # hide all design-prior args when de_an_prior == 1
    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
        arg %in% c("model_d","alpha_d","beta_d","location_d","k_d","scale_d","dff_d"))
      return(NULL)

    # hide internal dff and dff_d always
    if (arg %in% c("dff_d")) return(NULL)

    # analysis prior model-specific rules
    if (!is.null(x$model)) {
      if (x$model == "d_beta" &&
          arg %in% c("alpha","beta","scale")) return(NULL)
      if (x$model == "beta" &&
          arg %in% c("k","scale")) return(NULL)
      if (x$model == "NLP" &&
          arg %in% c("k","alpha","beta")) return(NULL)
    }

    # design prior model-specific rules
    if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$model_d)) {

      if (x$model_d == "d_beta" &&
          arg %in% c("alpha_d","beta_d","location_d","scale_d")) return(NULL)
      if (x$model_d == "beta" &&
          arg %in% c("k_d","location_d","scale_d")) return(NULL)
      if (x$model_d == "NLP" &&
          arg %in% c("k_d","alpha_d","beta_d")) return(NULL)
      if (x$model_d == "Point" &&
          arg %in% c("alpha_d","beta_d","k_d","scale_d")) return(NULL)
    }

    # skip null fields generally
    if (is.null(val)) return(NULL)

    ## ---------------------------------------------------------
    ## RENAME ARGUMENTS
    ## ---------------------------------------------------------

    arg_print <- arg

    # target > true_rate, FP > false_rate
    if (arg == "target") arg_print <- "true_rate"
    if (arg == "FP")     arg_print <- "false_rate"

    # direct > positive only when val=="h0" > positive="negative"
    if (arg == "direct") {
      if (val == "h0") {
        arg_print <- "type_rate"
        val <- "negative"
      } else {
        return(NULL) # skip h1 / positive
      }
    }

    # pc > plot_power only if TRUE
    if (arg == "pc") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_power"
      val <- TRUE
    }

    # rel > plot_rel only if TRUE
    if (arg == "rela") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_rel"
      val <- TRUE
    }

    ## HYPOTHESIS > ALTERNATIVE RULE
    if (arg == "hypothesis") {

      arg_print <- "alternative"

      val <- switch(val,
                    "<"  = "less",
                    "!=" = "two.sided",
                    ">"  = "greater",
                    stop("Invalid hypothesis")
      )

      val <- shQuote(val)
    }

    ## D > threshold
    if (arg == "D") arg_print <- "threshold"

    ## e > ROPE
    if (arg == "e") arg_print <- "ROPE"


    ## model > prior_analysis
    if (arg == "model") arg_print <- "prior_analysis"

    ## model_d > prior_design
    if (arg == "model_d") arg_print <- "prior_design"



    ## ---------------------------------------------------------
    ## VALUE FORMATTING
    ## ---------------------------------------------------------

    if (is.character(val) && !grepl("^\"", val))
      val <- shQuote(val)

    if (is.vector(val) && length(val) > 1)
      val <- paste0("c(", paste(val, collapse = ", "), ")")

    glue::glue("  {arg_print} = {val},")
  })

  # remove NULL entries
  code_lines <- code_lines[!sapply(code_lines, is.null)]

  ## remove trailing comma
  if (length(code_lines) > 0)
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])

  paste0(
    "BFpower.cor(\n",
    paste(code_lines, collapse = "\n"),
    "\n)"
  )
}




show_f_code <- function(x) {
  args <- c("inter","e", "D", "target","alpha", "p", "k",
            "model", "dff", "rscale", "f_m",
            "model_d","dff_d", "rscale_d", "f_m_d",
            "de_an_prior","N","mode_bf","direct",
            "pc","rela") # include plotting flags

  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]
    arg_print <- arg

    # ----- RULES: inter / e -----
    if (arg %in% c("inter","e")) {
      # never print inter
      if (arg == "inter") return(NULL)
      # print e only if inter != 1/"1"
      if (!is.null(x$inter) && (x$inter == 1 || x$inter == "1")) return(NULL)
    }

    # ----- RULES: target / FP -----
    if (arg == "target") {
      arg_print <- "true_rate"
      if (!is.null(x$mode_bf) && x$mode_bf != 1) return(NULL)
    }
    if (arg == "alpha") {
      arg_print <- "false_rate"
      if (!is.null(x$mode_bf) && x$mode_bf != 1) return(NULL)
    }



    # ----- RULES: direct > positive only if val=="h0" -----
    if (arg == "direct") {
      if (!is.null(val) && val == "h0") {
        arg_print <- "type_rate"
        val <- "negative"  # <-- no quotes!
      } else {
        return(NULL)  # skip h1
      }
    }

    # ----- RULES: plotting -----
    if (arg == "pc") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_power"
      val <- TRUE
    }
    if (arg == "rela") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_rel"
      val <- TRUE
    }

    # ----- RULES: de_an_prior -----
    if (arg == "de_an_prior") return(NULL) # never printed
    if (!is.null(x$de_an_prior) && (x$de_an_prior == 1 || x$de_an_prior == "1") &&
        arg %in% c("model_d","dff_d","rscale_d","f_m_d")) {
      return(NULL)
    }
    if (!is.null(x$de_an_prior) && (x$de_an_prior == 0 || x$de_an_prior == "0") &&
        !is.null(x$model_d)) {
      if (x$model_d == "Moment" && arg %in% c("rscale_d")) return(NULL)
      if (x$model_d == "Point" && arg %in% c("dff_d","rscale_d")) return(NULL)
    }

    # ----- RULES: model -----
    if (!is.null(x$model) && x$model != "effectsize" && arg == "rscale") return(NULL)

    # ----- RULES: mode_bf / N -----
    if (arg == "mode_bf") return(NULL)
    if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg == "N") return(NULL)


    ## D > threshold
    if (arg == "D") arg_print <- "threshold"

    ## e > ROPE
    if (arg == "e") arg_print <- "ROPE"


    ## model > prior_analysis
    if (arg == "model") arg_print <- "prior_analysis"

    ## model_d > prior_design
    if (arg == "model_d") arg_print <- "prior_design"



    # Skip NULL args
    if (is.null(val)) return(NULL)





    # Format values
    if (is.character(val)) {
      val <- shQuote(val, type = "cmd")
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")
    }

    glue::glue("  {arg_print} = {val},")
  })

  # Clean up commas
  code_lines <- code_lines[!sapply(code_lines, is.null)]
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  }

  # Final output
  paste0("BFpower.f.test(\n", paste(code_lines, collapse = "\n"), "\n)")
}





show_bin_code <- function(x) {
  args <- c("hypothesis","interval", "D", "target", "FP","h0","location",
            "model","alpha", "beta", "scale",
            "model_d","alpha_d", "beta_d", "location_d", "scale_d",
            "de_an_prior",
            "N", "mode_bf", "e", "direct","rela","pc")

  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]
    arg_print <- arg  # default printed name



    # ----- RULES: plotting -----
    if (arg == "pc") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_power"
      val <- TRUE
    }
    if (arg == "rela") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_rel"
      val <- TRUE
    }
    # ----- RULES: always skip location -----
    if (arg == "location") return(NULL)

    # ----- RULES: target / FP -----
    if (arg == "target") {
      arg_print <- "true_rate"
      if (!is.null(x$mode_bf) && x$mode_bf != 1) return(NULL)
    }
    if (arg == "FP") {
      arg_print <- "false_rate"
      if (!is.null(x$mode_bf) && x$mode_bf != 1) return(NULL)
    }

    # ----- RULES: direct > positive only if val=="h0" -----
    if (arg == "direct") {
      if (!is.null(val) && val == "h0") {
        arg_print <- "type_rate"
        val <- "negative"  # <-- no quotes!
      } else {
        return(NULL)  # skip h1
      }
    }
    # ----- RULES: model -----
    if (!is.null(x$model) && x$model == "beta" && arg == "scale") val <- NULL
    if (!is.null(x$model) && x$model != "beta" && arg %in% c("alpha","beta")) val <- NULL

    # ----- RULES: de_an_prior -----
    if (!is.null(x$de_an_prior)) {
      if (x$de_an_prior == 1 && arg %in% c("model_d","alpha_d","beta_d","location_d","scale_d")) val <- NULL
      if (arg == "de_an_prior") val <- NULL
    }

    # ----- RULES: model_d when de_an_prior = 0 -----
    if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$model_d)) {
      if (x$model_d == "beta" && arg %in% c("scale_d","location_d")) val <- val
      else if (x$model_d == "Moment" && arg %in% c("alpha_d","beta_d")) val <- NULL
      else if (x$model_d == "Point" && arg %in% c("alpha_d","beta_d","scale_d")) val <- NULL
    }

    # ----- RULES: mode_bf -----
    if (!is.null(x$mode_bf)) {
      if (arg == "mode_bf") val <- NULL
      if (x$mode_bf == 1 && arg == "N") val <- NULL
    }

    # ----- RULES: interval / e -----
    if (arg == "interval") return(NULL)  # interval is never printed
    if (arg == "e" && !is.null(x$interval) && x$interval == "1") val <- NULL



    ## HYPOTHESIS > ALTERNATIVE RULE
    if (arg == "hypothesis") {

      arg_print <- "alternative"

      val <- switch(val,
                    "<"  = "less",
                    "!=" = "two.sided",
                    ">"  = "greater",
                    stop("Invalid hypothesis")
      )
    }


    ## D > threshold
    if (arg == "D") arg_print <- "threshold"

    ## e > ROPE
    if (arg == "e") arg_print <- "ROPE"


    ## model > prior_analysis
    if (arg == "model") arg_print <- "prior_analysis"

    ## model_d > prior_design
    if (arg == "model_d") arg_print <- "prior_design"




    # Skip NULL args
    if (is.null(val)) return(NULL)

    # Format values
    if (is.character(val)) {
      val <- shQuote(val, type = "cmd")  # straight quotes
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")
    }

    glue::glue("  {arg_print} = {val},")
  })

  code_lines <- code_lines[!sapply(code_lines, is.null)]
  if (length(code_lines) > 0) code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  paste0("BFpower.bin(\n", paste(code_lines, collapse = "\n"), "\n)")
}


show_props_code <- function(x) {

  args <- c(
    "D", "target", "a0", "b0",
    "model1","a1", "b1", "a2", "b2",
    "model2", "a1d", "b1d", "dp1", "a2d", "b2d", "dp2",
    "mode_bf", "n1", "n2", "direct","pc","rela"
  )

  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # --- RULES: renaming ---
    # target → true_rate
    if (arg == "target") {
      arg_print <- "true_rate"
    } else {
      arg_print <- arg
    }

    # mode_bf: never printed
    if (arg == "mode_bf") return(NULL)

    # suppress n1/n2 if mode_bf == 1
    if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg %in% c("n1","n2")) {
      return(NULL)
    }

    # model1 rules
    if (!is.null(x$model1)) {
      if (x$model1 == "same" && arg %in% c("model1","a1d","b1d","dp1")) return(NULL)
      if (x$model1 == "beta" && arg == "dp1") val <- NULL
      if (x$model1 == "Point" && arg %in% c("a1d","b1d")) val <- NULL
    }

    # model2 rules
    if (!is.null(x$model2)) {
      if (x$model2 == "same" && arg %in% c("model2","a2d","b2d","dp2")) return(NULL)
      if (x$model2 == "beta" && arg == "dp2") val <- NULL
      if (x$model2 == "Point" && arg %in% c("a2d","b2d")) val <- NULL
    }

    # ----- RULES: plotting -----
    if (arg == "pc") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_power"
      val <- TRUE
    }
    if (arg == "rela") {
      if (!isTRUE(val)) return(NULL)
      arg_print <- "plot_rel"
      val <- TRUE
    }

    # ----- RULES: direct -----
    if (arg == "direct") {
      if (!is.null(val)) {
        if (val == "h1") return(NULL)        # ignore
        if (val == "h0") {                   # print type_rate = "negative"
          arg_print <- "type_rate"
          val <- "negative"                  # keep as plain string
        }
      }
    }


    # Skip NULL
    if (is.null(val)) return(NULL)

    # ----- Special renaming -----
    ## D → threshold
    if (arg == "D") arg_print <- "threshold"

    ## model1 → prior_design_1
    if (arg == "model1") arg_print <- "prior_design_1"

    ## model2 → prior_design_2
    if (arg == "model2") arg_print <- "prior_design_2"

    # ----- Value formatting -----
    # Format all character values
    if (is.character(val)) {
      val <- shQuote(val, type = "cmd")  # <- this will quote "negative" properly
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")
    }


    glue::glue("  {arg_print} = {val},")
  })

  # Remove NULL entries
  code_lines <- code_lines[!sapply(code_lines, is.null)]

  # Remove trailing comma
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  }

  paste0("BFpower.props(\n", paste(code_lines, collapse = "\n"), "\n)")
}
