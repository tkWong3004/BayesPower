show_t1_code <- function(x) {

  args <- c(
    "alternative","ROPE","threshold",
    "true_rate","false_rate","type_rate",
    "prior_analysis","location","scale","dff",
    "prior_design","location_d","scale_d","dff_d"
  )

  code_lines <- sapply(args, function(arg) {

    val <- x[[arg]]

    ## -----------------------------
    ## OMITTING RULES
    ## -----------------------------
    if (!is.null(x$interval) && x$interval == 1 && arg == "ROPE") return(NULL)
    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
        arg %in% c("prior_design","location_d","scale_d","dff_d")) return(NULL)
    if (is.null(val)) return(NULL)
    # analysis prior model-specific rules
    if (!is.null(x$prior_analysis)) {
      if (x$prior_analysis == "Normal" && arg %in% c("dff")) return(NULL)
      if (x$prior_analysis == "Moment" && arg %in% c("dff")) return(NULL)
    }

    # design prior model-specific rules
    if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$prior_design)) {
      if (x$prior_design == "Normal" && arg %in% c("dff_d")) return(NULL)
      if (x$prior_design == "Moment" && arg %in% c("dff_d")) return(NULL)
      if (x$prior_design == "Point"  && arg %in% c("dff_d","scale_d")) return(NULL)
    }


    ## -----------------------------
    ## SPECIAL RENAMING
    ## -----------------------------
    arg_print <- arg
    if (arg == "type_rate") { if (val == "positive") return(NULL)}

    ## -----------------------------
    ## VALUE FORMATTING
    ## -----------------------------
    if (is.character(val)) {
      val <- paste0('"', gsub('^"|"$', '', val), '"')
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")
    }

    glue::glue("  {arg_print} = {val},")
  })

  code_lines <- code_lines[!sapply(code_lines, is.null)]

  ## -----------------------------
  ## N LOGIC
  ## -----------------------------
  if (x$mode_bf != 1) {
    Nval <- if (is.null(x$N)) "NULL" else x$N
    code_lines <- c(code_lines, glue::glue("  N = {Nval}"))
  }

  ## Remove trailing comma
  if (length(code_lines) > 0)
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])

  ## -----------------------------
  ## Determine plot() arguments based on "pc" and "rela"
  ## -----------------------------
  plot_args <- c()

  if (isTRUE(x$pc)) {
    plot_args <- c(plot_args, "plot_power = TRUE")
  }

  if (isTRUE(x$rela)) {
    plot_args <- c(plot_args, "plot_rel = TRUE")
  }

  plot_args_text <- if (length(plot_args)) {
    paste0(", ", paste(plot_args, collapse = ", "))
  } else {
    ""
  }

  ## -----------------------------
  ## Final code
  ## -----------------------------
  code <- paste0(
    "results <- BFpower.ttest.OneSample(\n",
    paste(code_lines, collapse = "\n"),
    "\n)\n",
    "print(results)\n",
    "plot(results",
    plot_args_text,
    ")"
  )

  code
}

show_t2_code <- function(x) {

  args <- c(
    "alternative","ROPE","threshold",
    "true_rate","false_rate","type_rate",
    "prior_analysis","location","scale","dff",
    "prior_design","location_d","scale_d","dff_d",
    "r"    # keep r inside main arg loop
  )

  code_lines <- sapply(args, function(arg) {

    val <- x[[arg]]

    ## -----------------------------
    ## OMITTING RULES
    ## -----------------------------
    if (!is.null(x$interval) && x$interval == 1 && arg == "ROPE") return(NULL)
    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
        arg %in% c("prior_design","location_d","scale_d","dff_d")) return(NULL)
    if (x$mode_bf != 1 && arg == "r") return(NULL)
    if (arg == "type_rate") { if (val == "positive") return(NULL)}
    if (!is.null(x$prior_analysis)) {
      if (x$prior_analysis == "Normal" && arg %in% c("dff")) return(NULL)
      if (x$prior_analysis == "Moment" && arg %in% c("dff")) return(NULL)
    }

    # design prior model-specific rules
    if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$prior_design)) {
      if (x$prior_design == "Normal" && arg %in% c("dff_d")) return(NULL)
      if (x$prior_design == "Moment" && arg %in% c("dff_d")) return(NULL)
      if (x$prior_design == "Point"  && arg %in% c("dff_d","scale_d")) return(NULL)
    }
    if (is.null(val)) return(NULL)

    ## -----------------------------
    ## SPECIAL RENAMING
    ## -----------------------------
    arg_print <- arg


    ## -----------------------------
    ## VALUE FORMATTING
    ## -----------------------------
    if (is.character(val)) {
      val <- shQuote(val, type = "cmd")  # always double quotes
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")
    }

    glue::glue("  {arg_print} = {val},")
  })

  code_lines <- code_lines[!sapply(code_lines, is.null)]

  ## -----------------------------
  ## N1 / N2 LOGIC
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

  ## Remove trailing comma
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  }

  ## -----------------------------
  ## Determine plot() arguments based on "pc" and "rela"
  ## -----------------------------
  plot_args <- c()

  if (isTRUE(x$pc)) {
    plot_args <- c(plot_args, "plot_power = TRUE")
  }

  if (isTRUE(x$rela)) {
    plot_args <- c(plot_args, "plot_rel = TRUE")
  }

  plot_args_text <- if (length(plot_args)) {
    paste0(", ", paste(plot_args, collapse = ", "))
  } else {
    ""
  }
  ## -----------------------------
  ## Final code
  ## -----------------------------
  code <- paste0(
    "results <- BFpower.ttest.TwoSample(\n",
    paste(code_lines, collapse = "\n"),
    "\n)\n",
    "print(results)\n",
    "plot(results",
    plot_args_text,
    ")"
  )

  code
}

show_cor_code <- function(x) {

  args <- c(
    "alternative","h0","ROPE","threshold","true_rate","false_rate",
    "prior_analysis","k","alpha","beta","scale",
    "prior_design","alpha_d","beta_d","location_d","k_d","scale_d","dff_d",
    "N","type_rate"
  )

  code_lines <- sapply(args, function(arg) {

    val <- x[[arg]]

    ## -----------------------------
    ## OMITTING RULES
    ## -----------------------------
    if (!is.null(x$interval) && x$interval == 1 && arg == "ROPE") return(NULL)
    if (arg == "N" && x$mode_bf != 0) return(NULL)
    if (x$mode_bf == 0 && arg %in% c("true_rate","false_rate")) return(NULL)
    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
        arg %in% c("prior_design","alpha_d","beta_d","location_d","k_d","scale_d","dff_d"))
      return(NULL)
    if (arg %in% c("dff_d")) return(NULL)

    # analysis prior model-specific rules
    if (!is.null(x$prior_analysis)) {
      if (x$prior_analysis == "d_beta" && arg %in% c("alpha","beta","scale")) return(NULL)
      if (x$prior_analysis == "beta"   && arg %in% c("k","scale")) return(NULL)
      if (x$prior_analysis == "Moment" && arg %in% c("k","alpha","beta")) return(NULL)
    }

    # design prior model-specific rules
    if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$prior_design)) {
      if (x$prior_design == "d_beta" && arg %in% c("alpha_d","beta_d","location_d","scale_d")) return(NULL)
      if (x$prior_design == "beta"   && arg %in% c("k_d","location_d","scale_d")) return(NULL)
      if (x$prior_design == "Moment" && arg %in% c("k_d","alpha_d","beta_d")) return(NULL)
      if (x$prior_design == "Point"  && arg %in% c("alpha_d","beta_d","k_d","scale_d")) return(NULL)
    }

    if (is.null(val)) return(NULL)

    ## -----------------------------
    ## RENAME ARGUMENTS
    ## -----------------------------
    arg_print <- arg
    if (arg == "true_rate") arg_print <- "true_rate"
    if (arg == "false_rate")     arg_print <- "false_rate"
    if (arg == "type_rate") { if (val == "positive") return(NULL)}
    if (arg == "alternative") arg_print <- "alternative"
    if (arg == "threshold") arg_print <- "threshold"
    if (arg == "ROPE") arg_print <- "ROPE"
    if (arg == "prior_analysis") arg_print <- "prior_analysis"
    if (arg == "prior_design") arg_print <- "prior_design"

    ## -----------------------------
    ## VALUE FORMATTING
    ## -----------------------------
    if (is.character(val)) val <- shQuote(val, type = "cmd")
    if (is.vector(val) && length(val) > 1) val <- paste0("c(", paste(val, collapse = ", "), ")")

    glue::glue("  {arg_print} = {val},")
  })

  code_lines <- code_lines[!sapply(code_lines, is.null)]

  ## Remove trailing comma
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  }

  ## -----------------------------
  ## Determine plot() arguments based on "pc" and "rela"
  ## -----------------------------
  plot_args <- c()

  if (isTRUE(x$pc)) {
    plot_args <- c(plot_args, "plot_power = TRUE")
  }

  if (isTRUE(x$rela)) {
    plot_args <- c(plot_args, "plot_rel = TRUE")
  }

  plot_args_text <- if (length(plot_args)) {
    paste0(", ", paste(plot_args, collapse = ", "))
  } else {
    ""
  }

  ## -----------------------------
  ## Final code
  ## -----------------------------
  code <- paste0(
    "results <- BFpower.cor(\n",
    paste(code_lines, collapse = "\n"),
    "\n)\n",
    "print(results)\n",
    "plot(results",
    plot_args_text,
    ")"
  )

  code
}

show_f_code <- function(x) {
  args <- c("inter","ROPE", "threshold", "true_rate","false_rate", "p", "k",
            "prior_analysis", "dff", "rscale", "f_m",
            "prior_design","dff_d", "rscale_d", "f_m_d",
            "de_an_prior","N","type_rate")

  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]
    arg_print <- arg

    # ----- RULES: inter / ROPE -----
    if (arg == "inter") return(NULL)
    if (arg == "ROPE" && !is.null(x$inter) && (x$inter == 1 || x$inter == "1")) return(NULL)

    # ----- RULES: target / alpha -----
    if (arg == "true_rate") {
      arg_print <- "true_rate"
      if (!is.null(x$mode_bf) && x$mode_bf != 1) return(NULL)
    }
    if (arg == "false_rate") {
      arg_print <- "false_rate"
      if (!is.null(x$mode_bf) && x$mode_bf != 1) return(NULL)
    }

    # ----- RULES: type_rate -----
    if (arg == "type_rate") { if (val == "positive") return(NULL)}

    # ----- RULES: de_an_prior / prior_design -----
    if (arg == "de_an_prior") return(NULL)
    if (!is.null(x$de_an_prior) && (x$de_an_prior == 1 || x$de_an_prior == "1") &&
        arg %in% c("prior_design","dff_d","rscale_d","f_m_d")) return(NULL)
    if (!is.null(x$de_an_prior) && (x$de_an_prior == 0 || x$de_an_prior == "0") &&
        !is.null(x$prior_design)) {
      if (x$prior_design == "Moment" && arg %in% c("rscale_d")) return(NULL)
      if (x$prior_design == "Point"  && arg %in% c("dff_d","rscale_d")) return(NULL)
    }

    # ----- RULES: prior_analysis / rscale -----
    if (!is.null(x$prior_analysis) && x$prior_analysis != "effectsize" && arg == "rscale") return(NULL)

    # ----- RULES: mode_bf / N -----
    if (arg == "mode_bf") return(NULL)
    if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg == "N") return(NULL)

    # ----- Rename arguments -----
    if (arg == "threshold") arg_print <- "threshold"
    if (arg == "ROPE") arg_print <- "ROPE"
    if (arg == "prior_analysis") arg_print <- "prior_analysis"
    if (arg == "prior_design") arg_print <- "prior_design"

    # Skip NULL args
    if (is.null(val)) return(NULL)

    # Format values
    if (is.character(val)) val <- shQuote(val, type = "cmd")
    if (is.vector(val) && length(val) > 1) val <- paste0("c(", paste(val, collapse = ", "), ")")

    glue::glue("  {arg_print} = {val},")
  })

  # Clean up commas
  code_lines <- code_lines[!sapply(code_lines, is.null)]
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  }

  # ----- Determine plot() arguments based on pc/rela -----
  plot_args <- c()

  if (isTRUE(x$pc)) {
    plot_args <- c(plot_args, "plot_power = TRUE")
  }

  if (isTRUE(x$rela)) {
    plot_args <- c(plot_args, "plot_rel = TRUE")
  }

  plot_args_text <- if (length(plot_args)) {
    paste0(", ", paste(plot_args, collapse = ", "))
  } else {
    ""
  }

  # ----- Final code -----
  code <- paste0(
    "results <- BFpower.f.test(\n",
    paste(code_lines, collapse = "\n"),
    "\n)\n",
    "print(results)\n",
    "plot(results",
    plot_args_text,
    ")"
  )

  code
}

show_bin_code <- function(x) {
  args <- c("alternative", "threshold", "true_rate", "false_rate","h0",
            "prior_analysis","alpha", "beta", "scale",
            "prior_design","alpha_d", "beta_d","location_d","scale_d",
            "de_an_prior", "N","mode_bf","ROPE","type_rate") # removed pc/rela/location/interval from args

  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]
    arg_print <- arg  # default printed name

    # ----- RULES: skip location/interval -----
    if (arg %in% c("location","interval")) return(NULL)

    # ----- RULES: target / FP -----
    if (arg == "true_rate") {
      arg_print <- "true_rate"
      if (!is.null(x$mode_bf) && x$mode_bf != 1) return(NULL)
    }
    if (arg == "false_rate") {
      arg_print <- "false_rate"
      if (!is.null(x$mode_bf) && x$mode_bf != 1) return(NULL)
    }

    # ----- RULES
    if (arg == "type_rate") { if (val == "positive") return(NULL)}

    # ----- RULES: prior_analysis / prior_design / de_an_prior -----
    if (!is.null(x$prior_analysis)) {
      if (x$prior_analysis == "beta" && arg == "scale") val <- NULL
      if (x$prior_analysis != "beta" && arg %in% c("alpha","beta")) val <- NULL
    }
    if (!is.null(x$de_an_prior)) {
      if (x$de_an_prior == 1 && arg %in% c("prior_design","alpha_d","beta_d","location_d","scale_d")) val <- NULL
      if (arg == "de_an_prior") val <- NULL
    }
    if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$prior_design)) {
      if (x$prior_design == "Moment" && arg %in% c("alpha_d","beta_d")) val <- NULL
      if (x$prior_design == "Point" && arg %in% c("alpha_d","beta_d","scale_d")) val <- NULL
    }

    # ----- RULES: mode_bf / N -----
    if (!is.null(x$mode_bf)) {
      if (arg == "mode_bf") val <- NULL
      if (x$mode_bf == 1 && arg == "N") val <- NULL
    }

    # ----- RULES: ROPE -----
    if (arg == "ROPE" && !is.null(x$interval) && x$interval == "1") val <- NULL

    # ----- RENAME ARGUMENTS -----
    if (arg == "alternative") arg_print <- "alternative"
    if (arg == "threshold") arg_print <- "threshold"
    if (arg == "ROPE") arg_print <- "ROPE"
    if (arg == "prior_analysis") arg_print <- "prior_analysis"
    if (arg == "prior_design") arg_print <- "prior_design"

    # Skip NULL args
    if (is.null(val)) return(NULL)

    # Format values
    if (is.character(val)) val <- shQuote(val, type = "cmd")
    if (is.vector(val) && length(val) > 1) val <- paste0("c(", paste(val, collapse = ", "), ")")

    glue::glue("  {arg_print} = {val},")
  })

  # Clean up commas
  code_lines <- code_lines[!sapply(code_lines, is.null)]
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  }

  # ----- Determine plot() arguments based on pc/rela -----
  plot_args <- c()

  if (isTRUE(x$pc)) {
    plot_args <- c(plot_args, "plot_power = TRUE")
  }

  if (isTRUE(x$rela)) {
    plot_args <- c(plot_args, "plot_rel = TRUE")
  }

  plot_args_text <- if (length(plot_args)) {
    paste0(", ", paste(plot_args, collapse = ", "))
  } else {
    ""
  }

  # ----- Final code -----
  code <- paste0(
    "results <- BFpower.bin(\n",
    paste(code_lines, collapse = "\n"),
    "\n)\n",
    "print(results)\n",
    "plot(results",
    plot_args_text,
    ")"
  )

  code
}

show_props_code <- function(x) {

  args <- c(
    "threshold", "true_rate", "a0", "b0",
    "prior_design_1","a1", "b1", "a2", "b2",
    "prior_design_2", "a1d", "b1d", "dp1", "a2d", "b2d", "dp2",
    "mode_bf", "N1", "N2", "type_rate"   # removed pc/rela
  )

  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]
    arg_print <- arg

    # ----- RENAME -----
    if (arg == "true_rate") arg_print <- "true_rate"

    # mode_bf never printed
    if (arg == "mode_bf") return(NULL)

    # suppress n1/n2 if mode_bf == 1
    if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg %in% c("n1","n2"))
      return(NULL)

    # ----- prior prior_design_1 rules -----
    if (!is.null(x$prior_design_1)) {
      if (x$prior_design_1 == "same" && arg %in% c("prior_design_1","a1d","b1d","dp1"))
        return(NULL)
      if (x$prior_design_1 == "beta" && arg == "dp1")
        val <- NULL
      if (x$prior_design_1 == "Point" && arg %in% c("a1d","b1d"))
        val <- NULL
    }

    # ----- prior prior_design_2 rules -----
    if (!is.null(x$prior_design_2)) {
      if (x$prior_design_2 == "same" && arg %in% c("prior_design_2","a2d","b2d","dp2"))
        return(NULL)
      if (x$prior_design_2 == "beta" && arg == "dp2")
        val <- NULL
      if (x$prior_design_2 == "Point" && arg %in% c("a2d","b2d"))
        val <- NULL
    }

    if (arg == "type_rate") { if (val == "positive") return(NULL)}

    # Skip NULL values
    if (is.null(val)) return(NULL)

    # ----- Special renaming -----
    if (arg == "threshold") arg_print <- "threshold"
    if (arg == "prior_design_1") arg_print <- "prior_design_1"
    if (arg == "prior_design_2") arg_print <- "prior_design_2"

    # ----- Value formatting -----
    if (is.character(val))
      val <- shQuote(val, type = "cmd")

    if (is.vector(val) && length(val) > 1)
      val <- paste0("c(", paste(val, collapse = ", "), ")")

    glue::glue("  {arg_print} = {val},")
  })

  # Remove NULL entries
  code_lines <- code_lines[!sapply(code_lines, is.null)]

  # Remove trailing comma
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <-
      sub(",$", "", code_lines[length(code_lines)])
  }

  # ----- Determine plot() arguments -----
  plot_args <- c()

  if (isTRUE(x$pc)) {
    plot_args <- c(plot_args, "plot_power = TRUE")
  }

  if (isTRUE(x$rela)) {
    plot_args <- c(plot_args, "plot_rel = TRUE")
  }

  plot_args_text <- if (length(plot_args)) {
    paste0(", ", paste(plot_args, collapse = ", "))
  } else {
    ""
  }

  # ----- Final code -----
  code <- paste0(
    "results <- BFpower.props(\n",
    paste(code_lines, collapse = "\n"),
    "\n)\n",
    "print(results)\n",
    "plot(results",
    plot_args_text,
    ")"
  )
  code
}
