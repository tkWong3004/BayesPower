
show_t1_code <- function(x) {

  # List of all arguments for bp_t.test_one_sample
  args <- c("hypothesis","e","interval", "D", "target", "alpha",
            "model", "location", "scale", "dff",
            "model_d", "location_d", "scale_d", "dff_d",
            "de_an_prior",
            "N", "mode_bf", "direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # If de_an_prior == 1, set design prior args to NULL
    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 && arg %in% c("model_d","location_d","scale_d","dff_d")) {
      val <- NULL
    }

    # If mode_bf == 1, set N to NULL
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg == "N") {
      val <- NULL
    }

    # If interval == "1", set e to NULL
    else if (!is.null(x$interval) && x$interval == "1" && arg == "e") {
      val <- NULL
    }

    # Skip arguments that are NULL
    if (is.null(val)) return(NULL)

    # Format values
    if (is.character(val)) {
      val <- shQuote(val, type = "cmd")  # use straight quotes
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")
    }

    glue::glue("  {arg} = {val},")
  })

  # Remove NULL entries (those that were skipped)
  code_lines <- code_lines[!sapply(code_lines, is.null)]

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
      val <- NULL
    }

    # If mode_bf == 1, set N1 and N2 to NULL
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 &&
             arg %in% c("N1","N2")) {
      val <- NULL
    }

    # If interval == "1", set e to NULL
    else if (!is.null(x$interval) && x$interval == "1" && arg == "e") {
      val <- NULL
    }

    # Skip NULL arguments entirely
    if (is.null(val)) return(NULL)

    # Normal formatting rules
    if (is.character(val)) {
      val <- shQuote(val, type = "cmd")       # wrap in straight quotes
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")
    }

    glue::glue("  {arg} = {val},")
  })

  # Remove NULL entries
  code_lines <- code_lines[!sapply(code_lines, is.null)]

  # Remove trailing comma from last line
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  }

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
            "model_d","alpha_d","beta_d","location_d","k_d","scale_d","dff_d",
            "de_an_prior",
            "N","mode_bf","direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # ----- RULES for model -----
    if (!is.null(x$model) && x$model == "d_beta" &&
        arg %in% c("alpha","beta","scale","dff")) {
      val <- NULL
    }
    else if (!is.null(x$model) && x$model == "beta" &&
             arg %in% c("k","scale","dff")) {
      val <- NULL
    }
    else if (!is.null(x$model) && x$model == "NLP" &&
             arg %in% c("k","alpha","beta")) {
      val <- NULL
    }

    # ----- RULES for de_an_prior + model_d -----
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
             arg %in% c("model_d","alpha_d","beta_d","location_d","k_d","scale_d","dff_d")) {
      val <- NULL
    }
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$model_d)) {
      if (x$model_d == "d_beta" &&
          arg %in% c("alpha_d","beta_d","location_d","scale_d","dff_d")) {
        val <- NULL
      }
      else if (x$model_d == "beta" &&
               arg %in% c("k_d","location_d","scale_d","dff_d")) {
        val <- NULL
      }
      else if (x$model_d == "NLP" &&
               arg %in% c("k_d","alpha_d","beta_d")) {
        val <- NULL
      }
      else if (x$model_d == "Point" &&
               arg %in% c("alpha_d","beta_d","k_d","scale_d","dff_d")) {
        val <- NULL
      }
    }

    # ----- RULES for mode_bf -----
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg == "N") {
      val <- NULL
    }

    # ----- RULES for interval -----
    else if (!is.null(x$interval) && x$interval == "1" && arg == "e") {
      val <- NULL
    }

    # Skip NULL arguments entirely
    if (is.null(val)) return(NULL)

    # ----- DEFAULT formatting -----
    if (is.character(val)) {
      val <- shQuote(val, type = "cmd")       # wrap strings in quotes
    } else if (is.vector(val) && length(val) > 1) {
      val <- paste0("c(", paste(val, collapse = ", "), ")")  # format vectors
    }

    glue::glue("  {arg} = {val},")
  })

  # Remove NULL entries
  code_lines <- code_lines[!sapply(code_lines, is.null)]

  # Remove trailing comma from last line
  if (length(code_lines) > 0) {
    code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  }

  # Combine into multi-line function call
  code <- paste0("bp_cor(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}


show_f_code <- function(x) {
  args <- c("interval","e", "D", "target","FP", "p", "k",
            "model", "dff", "rscale", "f_m",
            "model_d","dff_d", "rscale_d", "f_m_d",
            "de_an_prior",
            "N","mode_bf", "direct")

  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    # ----- RULES -----
    if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
        arg %in% c("model_d","dff_d","rscale_d","f_m_d")) val <- "NULL"
    else if (!is.null(x$model) && x$model != "effectsize" && arg == "rscale") val <- "NULL"
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg == "N") val <- "NULL"
    else if (!is.null(x$interval) && x$interval != "1" && arg == "e") val <- "NULL"
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$model_d)) {
      if (x$model_d != "effectsize" && arg == "rscale_d") val <- "NULL"
      else if (x$model_d == "Point" && arg %in% c("dff_d","rscale_d")) val <- "NULL"
    }

    # Skip NULL args
    if (is.null(val)) return(NULL)

    if (is.character(val)) val <- shQuote(val, type = "cmd")
    else if (is.vector(val) && length(val) > 1) val <- paste0("c(", paste(val, collapse = ", "), ")")

    glue::glue("  {arg} = {val},")
  })

  code_lines <- code_lines[!sapply(code_lines, is.null)]
  code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  paste0("bp_f(\n", paste(code_lines, collapse = "\n"), "\n)")
}

show_bin_code <- function(x) {
  args <- c("hypothesis","interval", "D", "target", "FP","location",
            "model","alpha", "beta", "scale",
            "model_d","alpha_d", "beta_d", "location_d", "scale_d",
            "de_an_prior",
            "N", "mode_bf", "e", "direct")

  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    if (!is.null(x$model) && x$model == "beta" && arg == "scale") val <- "NULL"
    else if (!is.null(x$model) && x$model != "beta" && arg %in% c("alpha","beta")) val <- "NULL"
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 1 &&
             arg %in% c("model_d","alpha_d","beta_d","location_d","scale_d")) val <- "NULL"
    else if (!is.null(x$de_an_prior) && x$de_an_prior == 0 && !is.null(x$model_d)) {
      if (x$model_d == "beta" && arg %in% c("scale_d","location_d")) val <- "NULL"
      else if (x$model_d == "Moment" && arg %in% c("alpha_d","beta_d","location_d")) val <- "NULL"
      else if (x$model_d == "Point" && arg %in% c("alpha_d","beta_d","scale_d")) val <- "NULL"
    }
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg == "N") val <- "NULL"
    else if (!is.null(x$interval) && x$interval == "1" && arg == "e") val <- "NULL"

    # Skip NULL args
    if (is.null(val)) return(NULL)

    if (is.character(val)) val <- shQuote(val, type = "cmd")
    else if (is.vector(val) && length(val) > 1) val <- paste0("c(", paste(val, collapse = ", "), ")")

    glue::glue("  {arg} = {val},")
  })

  code_lines <- code_lines[!sapply(code_lines, is.null)]
  code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  paste0("bp_bin(\n", paste(code_lines, collapse = "\n"), "\n)")
}

show_props_code <- function(x) {
  args <- c("D", "target", "a0", "b0",
            "model1","a1", "b1", "a2", "b2",
            "model2", "a1d", "b1d", "dp1", "a2d", "b2d", "dp2",
            "mode_bf", "n1", "n2", "direct")

  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    if (!is.null(x$model1) && !is.null(x$model2) && x$model1 == x$model2 &&
        arg %in% c("a1d","b1d","dp1","a2d","b2d","dp2")) val <- "NULL"
    else if (!is.null(x$mode_bf) && x$mode_bf == 1 && arg %in% c("n1","n2")) val <- "NULL"
    else if (!is.null(x$model1) && !is.null(x$model2) && x$model1 != x$model2) {
      if (x$model1 == "beta" && arg == "dp1") val <- "NULL"
      else if (x$model1 == "Point" && arg %in% c("a1d","b1d")) val <- "NULL"
      else if (x$model2 == "beta" && arg == "dp2") val <- "NULL"
      else if (x$model2 == "Point" && arg %in% c("a2d","b2d")) val <- "NULL"
    }

    # Skip NULL args
    if (is.null(val)) return(NULL)

    if (is.character(val)) val <- shQuote(val, type = "cmd")
    else if (is.vector(val) && length(val) > 1) val <- paste0("c(", paste(val, collapse = ", "), ")")

    glue::glue("  {arg} = {val},")
  })

  code_lines <- code_lines[!sapply(code_lines, is.null)]
  code_lines[length(code_lines)] <- sub(",$", "", code_lines[length(code_lines)])
  paste0("bp_props(\n", paste(code_lines, collapse = "\n"), "\n)")
}


