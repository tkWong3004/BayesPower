
show_t1_code <- function(x) {

  # List of all arguments for bayespower_t.test_one_sample
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
  code <- paste0("bayespower_t.test_one_sample(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}




show_t2_code <- function(x) {

  # List of all arguments for two-sample function
  args <- c("D","r","target","model","location","scale","dff","hypothesis",
            "model_d","location_d","scale_d","dff_d","de_an_prior",
            "N1","N2","mode_bf","alpha","direct","e","interval")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    if (is.character(val)) {
      val <- dQuote(val)       # quote strings
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
  code <- paste0("bayespower_t.test_two_sample(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}
show_cor_code <- function(x) {

  # List of all arguments for correlation function
  args <- c("interval", "D","target","model","k","alpha","beta","h0","location",
            "scale","dff","hypothesis","model_d","location_d","k_d","alpha_d",
            "beta_d","scale_d","dff_d","de_an_prior","N","mode_bf","FP","e","direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    if (is.character(val)) {
      val <- dQuote(val)       # quote strings
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
  code <- paste0("bayespower_cor(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}
show_f_code <- function(x) {

  # List of all arguments for bayespower_f
  args <- c("interval", "D", "target", "p", "k", "dff", "rscale", "f_m", "model",
            "dff_d", "rscale_d", "f_m_d", "model_d", "de_an_prior", "N",
            "mode_bf", "FP", "e", "direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    if (is.character(val)) {
      val <- dQuote(val)       # quote strings
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
  code <- paste0("bayespower_f(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}
show_bin_code <- function(x) {

  # List of all arguments for bayespower_bin
  args <- c("interval", "D", "target", "alpha", "beta", "location", "scale", "model",
            "hypothesis", "alpha_d", "beta_d", "location_d", "scale_d", "model_d",
            "de_an_prior", "N", "mode_bf", "FP", "e", "direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    if (is.character(val)) {
      val <- dQuote(val)       # quote strings
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
  code <- paste0("bayespower_bin(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}
show_props_code <- function(x) {

  # List of all arguments for bayespower_props
  args <- c("D", "target", "a0", "b0", "a1", "b1", "a2", "b2", "r",
            "model1", "a1d", "b1d", "dp1", "model2", "a2d", "b2d", "dp2",
            "mode_bf", "n1", "n2", "direct")

  # Build code lines dynamically
  code_lines <- sapply(args, function(arg) {
    val <- x[[arg]]

    if (is.character(val)) {
      val <- dQuote(val)       # quote strings
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
  code <- paste0("bayespower_props(\n",
                 paste(code_lines, collapse = "\n"),
                 "\n)")

  return(code)
}
