
show_t1_code <- function(x) {

  # List of all arguments for bp_t.test_one_sample
  args <- c( "hypothesis","e","interval", "D", "target", "alpha",
            "model", "location", "scale", "dff",
            "model_d", "location_d", "scale_d", "dff_d",
            "de_an_prior",
            "N", "mode_bf",  "direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # If de_an_prior == 1, set specific args to NULL
    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 && arg %in% c("model_d","location_d","scale_d","dff_d")) {
      val <- "NULL"
    }

    # If mode_bf == 1, set N to NULL
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg == "N") {
      val <- "NULL"
    }

    # If interval == "1", set e to NULL
    else if (!is.null(x$interval) && x$interval == "1" && arg == "e") {
      val <- "NULL"
    }

    # Handle normal formatting
    else if (is.character(val)) {
      val <- shQuote(val, type = "cmd")  # use straight quotes
    } else if (is.null(val)) {
      val <- "NULL"
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")  # format vectors
    }

    glue::glue("  {arg} = {val},")
  })

  # Remove trailing comma from last line
  code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])

  # Combine into multi-line function call
  code <- paste0("bp_t.test_one_sample(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}




show_t2_code <- function(x) {

  # List of all arguments for two-sample function
  args <- c("hypothesis","e","interval","D","target","alpha",
            "model","location","scale","dff",
            "model_d","location_d","scale_d","dff_d",
            "de_an_prior",
            "N1","N2","r","mode_bf","direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # If de_an_prior == 1, set *_d args to NULL
    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
        arg %in% c("model_d","location_d","scale_d","dff_d")) {
      val <- "NULL"
    }

    # If mode_bf == 1, set N1 and N2 to NULL
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 &&
             arg %in% c("N1","N2")) {
      val <- "NULL"
    }

    # If interval == "1", set e to NULL
    else if (!is.null(x$interval) && x$interval == "1" && arg == "e") {
      val <- "NULL"
    }

    # Normal formatting rules
    else if (is.character(val)) {
      val <- shQuote(val, type = "cmd")       # wrap in straight quotes
    } else if (is.null(val)) {
      val <- "NULL"
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")
    }

    glue::glue("  {arg} = {val},")
  })

  # Remove trailing comma from last line
  code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])

  # Combine into multi-line function call
  code <- paste0("bp_t.test_two_sample(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}
show_cor_code <- function(x) {

  # List of all arguments for correlation function
  args <- c("hypothesis","h0","e","interval", "D","target","FP",
            "model","k","alpha","beta","scale",
            "model_d","alpha_d","beta_d","location_d","k_d","scale_d",
            "de_an_prior",
            "N","mode_bf","direct")

  # Build code lines dynamically
  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # ----- RULES for model -----
    if (!is.null(x$model) && x$model == "d_beta" &&
        arg %in% c("alpha","beta","scale","dff")) {
      val <- "NULL"
    }
    else if (!is.null(x$model) && x$model == "beta" &&
             arg %in% c("k","scale","dff")) {
      val <- "NULL"
    }
    else if (!is.null(x$model) && x$model == "NLP" &&
             arg %in% c("k","alpha","beta")) {
      val <- "NULL"
    }

    # ----- RULES for de_an_prior + model_d -----
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
             arg %in% c("model_d","alpha_d","beta_d","location_d","k_d","scale_d","dff_d")) {
      val <- "NULL"
    }
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$model_d)) {
      if (x$model_d == "d_beta" &&
          arg %in% c("alpha_d","beta_d","location_d","scale_d","dff_d")) {
        val <- "NULL"
      }
      else if (x$model_d == "beta" &&
               arg %in% c("k_d","location_d","scale_d","dff_d")) {
        val <- "NULL"
      }
      else if (x$model_d == "NLP" &&
               arg %in% c("k_d","alpha_d","beta_d")) {
        val <- "NULL"
      }
      else if (x$model_d == "Point" &&
               arg %in% c("alpha_d","beta_d","k_d","scale_d","dff_d")) {
        val <- "NULL"
      }
    }

    # ----- RULES for mode_bf -----
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg == "N") {
      val <- "NULL"
    }

    # ----- RULES for interval -----
    else if (!is.null(x$interval) && x$interval == "1" && arg == "e") {
      val <- "NULL"
    }

    # ----- DEFAULT formatting -----
    else if (is.character(val)) {
      val <- shQuote(val, type = "cmd")       # wrap strings in quotes
    } else if (is.null(val)) {
      val <- "NULL"
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")  # format vectors
    }

    glue::glue("  {arg} = {val},")
  })

  # Remove trailing comma from last line
  code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])

  # Combine into multi-line function call
  code <- paste0("bp_cor(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}

show_f_code <- function(x) {

  # List of all arguments for bpf
  args <- c("interval","e", "D", "target","FP", "p", "k",
            "model", "dff", "rscale", "f_m",
            "model_d","dff_d", "rscale_d", "f_m_d",
            "de_an_prior",
            "N","mode_bf", "direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # ----- RULES -----
    # de_an_prior = 1 → NULL out all *_d args
    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
        arg %in% c("model_d","dff_d","rscale_d","f_m_d")) {
      val <- "NULL"
    }

    else if (!is.null(x$model) && x$model != "effectsize" &&
             arg == "rscale") {
      val <- "NULL"
    }
    # mode_bf == 1 → N = NULL
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 &&
             arg == "N") {
      val <- "NULL"
    }
    # interval == "1" → e = NULL
    else if (!is.null(x$interval) && x$interval != "1" &&
             arg == "e") {
      val <- "NULL"
    }

    # ----- Extra RULES when de_an_prior == 0 -----
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$model_d)) {
      # model_d = "effectsize" → rscale_d = NULL
      if (x$model_d != "effectsize" && arg == "rscale_d") {
        val <- "NULL"
      }
      # model_d = "Point" → dff_d, rscale_d = NULL
      else if (x$model_d == "Point" && arg %in% c("dff_d","rscale_d")) {
        val <- "NULL"
      }
    }

    # ----- DEFAULT formatting -----
    else if (is.character(val)) {
      val <- shQuote(val, type = "cmd")       # quote strings
    } else if (is.null(val)) {
      val <- "NULL"
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")  # format vectors
    }

    glue::glue("  {arg} = {val},")
  })

  # Remove trailing comma from last line
  code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])

  # Combine into multi-line function call
  code <- paste0("bp_f(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}


show_bin_code <- function(x) {

  # List of all arguments for bpbin
  args <- c("hypothesis","interval", "D", "target", "FP","location",
            "model","alpha", "beta", "scale",
            "model_d","alpha_d", "beta_d", "location_d", "scale_d",
            "de_an_prior",
            "N", "mode_bf", "e", "direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # ----- RULES -----
    # Model-specific
    if (!is.null(x$model) && x$model == "beta" && arg == "scale") {
      val <- "NULL"
    }
    else if (!is.null(x$model) && x$model != "beta" &&
             arg %in% c("alpha","beta")) {
      val <- "NULL"
    }

    # de_an_prior = 1 → NULL out all *_d
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
             arg %in% c("model_d","alpha_d","beta_d","location_d","scale_d")) {
      val <- "NULL"
    }

    # de_an_prior = 0 → model_d-specific rules
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$model_d)) {
      # model_d = "beta"
      if (x$model_d == "beta" && arg %in% c("scale_d","location_d")) {
        val <- "NULL"
      }
      # model_d = "Moment"
      else if (x$model_d == "Moment" && arg %in% c("alpha_d","beta_d","location_d")) {
        val <- "NULL"
      }
      # model_d = "Point"
      else if (x$model_d == "Point" && arg %in% c("alpha_d","beta_d","scale_d")) {
        val <- "NULL"
      }
    }

    # mode_bf == 1 → N = NULL
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg == "N") {
      val <- "NULL"
    }

    # interval == "1" → e = NULL
    else if (!is.null(x$interval) && x$interval == "1" && arg == "e") {
      val <- "NULL"
    }

    # ----- DEFAULT formatting -----
    else if (is.character(val)) {
      val <- shQuote(val, type = "cmd")       # quote strings
    } else if (is.null(val)) {
      val <- "NULL"
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")  # format vectors
    }

    glue::glue("  {arg} = {val},")
  })

  # Remove trailing comma from last line
  code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])

  # Combine into multi-line function call
  code <- paste0("bp_bin(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}

show_props_code <- function(x) {

  # List of all arguments for bpprops
  args <- c("D", "target", "a0", "b0",
            "model1","a1", "b1", "a2", "b2",
            "model2", "a1d", "b1d", "dp1", "a2d", "b2d", "dp2",
            "mode_bf", "n1", "n2", "direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # ----- RULES -----
    # If model1 == model2 → NULL out all *_d terms
    if (!is.null(x$model1) && !is.null(x$model2) &&
        x$model1 == x$model2 &&
        arg %in% c("a1d","b1d","dp1","a2d","b2d","dp2")) {
      val <- "NULL"
    }

    # If mode_bf == 1 → n1, n2 = NULL
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 &&
             arg %in% c("n1","n2")) {
      val <- "NULL"
    }

    # If model1 != model2 → apply model-specific rules
    else if (!is.null(x$model1) && !is.null(x$model2) && x$model1 != x$model2) {
      # model1 = beta → dp1 = NULL
      if (x$model1 == "beta" && arg == "dp1") {
        val <- "NULL"
      }
      # model1 = Point → a1d, b1d = NULL
      else if (x$model1 == "Point" && arg %in% c("a1d","b1d")) {
        val <- "NULL"
      }
      # model2 = beta → dp2 = NULL
      else if (x$model2 == "beta" && arg == "dp2") {
        val <- "NULL"
      }
      # model2 = Point → a2d, b2d = NULL
      else if (x$model2 == "Point" && arg %in% c("a2d","b2d")) {
        val <- "NULL"
      }
    }

    # ----- DEFAULT formatting -----
    else if (is.character(val)) {
      val <- shQuote(val, type = "cmd")       # quote strings
    } else if (is.null(val)) {
      val <- "NULL"
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")  # format vectors
    }

    glue::glue("  {arg} = {val},")
  })

  # Remove trailing comma from last line
  code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])

  # Combine into multi-line function call
  code <- paste0("bp_props(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}

