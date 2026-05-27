############################################################
## Comprehensive test script for:
## 1) BF10.cor()
## 2) BFpower.cor()
############################################################

## Load your package first, e.g.
## devtools::load_all()
## or
## library(BayesPower)

############################################################
## Helper
############################################################


############################################################
## testthat setup
############################################################

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Package 'testthat' is required. Install it with install.packages('testthat').", call. = FALSE)
}

testthat::local_edition(3)

run_test <- function(description, expr, expect_error = grepl("^ERROR", description)) {
  testthat::test_that(description, {
    if (isTRUE(expect_error)) {
      testthat::expect_error(eval(expr))
    } else {
      testthat::expect_no_error(eval(expr))
    }
  })
}




############################################################
## PART A: BF10.cor()
############################################################

############################################################
## A1. Point-null BF10 tests: ROPE = NULL
############################################################

for (alt in c("two.sided", "greater", "less")) {
  for (prior in c("d_beta", "beta", "Moment")) {

    args <- list(
      r = 0.3930924,
      n = 46,
      h0 = 0,
      alternative = alt,
      prior_analysis = prior
    )

    if (prior == "d_beta") {
      args$k <- 1
    }

    if (prior == "beta") {
      args$alpha <- 1
      args$beta <- 1
    }

    if (prior == "Moment") {
      args$scale <- 0.707
    }

    run_test(
      paste("BF10.cor | point-null |", alt, "|", prior),
      as.call(c(quote(BF10.cor), args))
    )
  }
}


############################################################
## A2. Interval-null BF10 tests: ROPE supplied
############################################################

rope_by_alt <- list(
  "two.sided" = c(-0.2, 0.2),
  "greater" = 0.2,
  "less" = -0.2
)

for (alt in c("two.sided", "greater", "less")) {
  for (prior in c("d_beta", "beta", "Moment")) {

    args <- list(
      r = 0.3930924,
      n = 46,
      h0 = 0,
      ROPE = rope_by_alt[[alt]],
      alternative = alt,
      prior_analysis = prior
    )

    if (prior == "d_beta") {
      args$k <- 1
    }

    if (prior == "beta") {
      args$alpha <- 1
      args$beta <- 1
    }

    if (prior == "Moment") {
      args$scale <- 0.707
    }

    run_test(
      paste("BF10.cor | interval-null |", alt, "|", prior),
      as.call(c(quote(BF10.cor), args))
    )
  }
}


############################################################
## A3. Non-zero h0 tests
############################################################

run_test(
  "BF10.cor | point-null | h0 = 0.1 | greater | d_beta",
  quote(
    BF10.cor(
      r = 0.3930924,
      n = 46,
      h0 = 0.1,
      alternative = "greater",
      prior_analysis = "d_beta",
      k = 1
    )
  )
)

run_test(
  "BF10.cor | interval-null | h0 = 0.1 | two.sided | d_beta",
  quote(
    BF10.cor(
      r = 0.3930924,
      n = 46,
      h0 = 0.1,
      ROPE = c(-0.2, 0.2),
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  )
)

run_test(
  "BF10.cor | interval-null | h0 = 0.1 | greater | beta",
  quote(
    BF10.cor(
      r = 0.3930924,
      n = 46,
      h0 = 0.1,
      ROPE = 0.2,
      alternative = "greater",
      prior_analysis = "beta",
      alpha = 1,
      beta = 1
    )
  )
)

run_test(
  "BF10.cor | interval-null | h0 = 0.1 | less | Moment",
  quote(
    BF10.cor(
      r = -0.25,
      n = 46,
      h0 = 0.1,
      ROPE = -0.2,
      alternative = "less",
      prior_analysis = "Moment",
      scale = 0.707
    )
  )
)




############################################################
## A3b. Additional non-zero h0 tests across alternatives/priors
############################################################

for (h0_val in c(-0.3, 0.1, 0.3)) {
  for (alt in c("two.sided", "greater", "less")) {
    for (prior in c("d_beta", "beta", "Moment")) {

      args <- list(
        r = 0.3930924,
        n = 46,
        h0 = h0_val,
        alternative = alt,
        prior_analysis = prior
      )

      if (prior == "d_beta") {
        args$k <- 1
      }

      if (prior == "beta") {
        args$alpha <- 1
        args$beta <- 1
      }

      if (prior == "Moment") {
        args$scale <- 0.707
      }

      run_test(
        paste("BF10.cor | point-null | h0 =", h0_val, "|", alt, "|", prior),
        as.call(c(quote(BF10.cor), args))
      )
    }
  }
}



############################################################
## A3c. Additional non-zero h0 BF10.cor tests with ROPE
############################################################

for (h0_val in c(-0.3, 0.1, 0.3)) {
  for (alt in c("two.sided", "greater", "less")) {
    for (prior in c("d_beta", "beta", "Moment")) {

      args <- list(
        r = 0.3930924,
        n = 46,
        h0 = h0_val,
        ROPE = rope_by_alt[[alt]],
        alternative = alt,
        prior_analysis = prior
      )

      if (prior == "d_beta") {
        args$k <- 1
      }

      if (prior == "beta") {
        args$alpha <- 1
        args$beta <- 1
      }

      if (prior == "Moment") {
        args$scale <- 0.707
      }

      run_test(
        paste("BF10.cor | interval-null | h0 =", h0_val, "|", alt, "|", prior),
        as.call(c(quote(BF10.cor), args))
      )
    }
  }
}
############################################################
## A4. BF10.cor object structure test
############################################################

test_that("BF10.cor returns the expected object structure", {
  bf10_cor_obj <- BF10.cor(
    r = 0.3930924,
    n = 46,
    prior_analysis = "d_beta",
    k = 1,
    h0 = 0,
    alternative = "two.sided"
  )

  expect_s3_class(bf10_cor_obj, "BFvalue")

  expect_true("type" %in% names(bf10_cor_obj))
  expect_true("bf10" %in% names(bf10_cor_obj))
  expect_true("h0" %in% names(bf10_cor_obj))
  expect_true("r" %in% names(bf10_cor_obj))
  expect_true("n" %in% names(bf10_cor_obj))
  expect_true("analysis_h1" %in% names(bf10_cor_obj))
  expect_true("alternative" %in% names(bf10_cor_obj))
  expect_true("ROPE" %in% names(bf10_cor_obj))
  expect_true("p.value" %in% names(bf10_cor_obj))
})


############################################################
## A5. BF10.cor validation/error tests
############################################################

run_test(
  "ERROR BF10.cor | invalid r greater than 1",
  quote(
    BF10.cor(
      r = 1.2,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid r less than -1",
  quote(
    BF10.cor(
      r = -1.2,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid r NA",
  quote(
    BF10.cor(
      r = NA,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid n <= 3",
  quote(
    BF10.cor(
      r = 0.3,
      n = 3,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid h0 too large",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0.9,
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid h0 too small",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = -0.9,
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid alternative",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      alternative = "wrong",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid prior_analysis",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "Normal",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | missing k for d_beta",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "d_beta"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid k",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | missing alpha/beta for beta prior",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "beta"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid alpha",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "beta",
      alpha = 0,
      beta = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid beta",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "beta",
      alpha = 1,
      beta = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | missing scale for Moment",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "Moment"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | invalid scale for Moment",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      alternative = "two.sided",
      prior_analysis = "Moment",
      scale = 0
    )
  ),
  expect_error = TRUE
)


############################################################
## A6. BF10.cor ROPE validation/error tests
############################################################

run_test(
  "ERROR BF10.cor | two.sided ROPE length 1",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      ROPE = 0.2,
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | two.sided ROPE both positive",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      ROPE = c(0.1, 0.2),
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | two.sided ROPE both negative",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      ROPE = c(-0.3, -0.2),
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | two.sided ROPE reversed signs",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      ROPE = c(0.2, -0.2),
      alternative = "two.sided",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | greater ROPE negative",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      ROPE = -0.2,
      alternative = "greater",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.cor | less ROPE positive",
  quote(
    BF10.cor(
      r = 0.3,
      n = 46,
      h0 = 0,
      ROPE = 0.2,
      alternative = "less",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)


############################################################
## PART B: BFpower.cor()
############################################################

############################################################
## B1. Sample size determination mode, point-null
## N = NULL
############################################################

for (alt in c("two.sided", "greater", "less")) {
  for (prior in c("d_beta", "beta", "Moment")) {

    args <- list(
      alternative = alt,
      h0 = 0,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = prior
    )

    if (prior == "d_beta") {
      args$k <- 1
    }

    if (prior == "beta") {
      args$alpha <- 1
      args$beta <- 1
    }

    if (prior == "Moment") {
      args$scale <- 0.707
    }

    run_test(
      paste("BFpower.cor | sample size | point-null |", alt, "|", prior),
      as.call(c(quote(BFpower.cor), args))
    )
  }
}


############################################################
## B2. Sample size determination mode, interval-null
############################################################

for (alt in c("two.sided", "greater", "less")) {
  for (prior in c("d_beta", "beta", "Moment")) {

    args <- list(
      alternative = alt,
      h0 = 0,
      ROPE = rope_by_alt[[alt]],
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = prior
    )

    if (prior == "d_beta") {
      args$k <- 1
    }

    if (prior == "beta") {
      args$alpha <- 1
      args$beta <- 1
    }

    if (prior == "Moment") {
      args$scale <- 0.707
    }

    run_test(
      paste("BFpower.cor | sample size | interval-null |", alt, "|", prior),
      as.call(c(quote(BFpower.cor), args))
    )
  }
}


############################################################
## B3. Fixed-sample power mode, point-null
############################################################

for (alt in c("two.sided", "greater", "less")) {
  for (prior in c("d_beta", "beta", "Moment")) {

    args <- list(
      alternative = alt,
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = prior
    )

    if (prior == "d_beta") {
      args$k <- 1
    }

    if (prior == "beta") {
      args$alpha <- 1
      args$beta <- 1
    }

    if (prior == "Moment") {
      args$scale <- 0.707
    }

    run_test(
      paste("BFpower.cor | fixed N | point-null |", alt, "|", prior),
      as.call(c(quote(BFpower.cor), args))
    )
  }
}


############################################################
## B4. Fixed-sample power mode, interval-null
############################################################

for (alt in c("two.sided", "greater", "less")) {
  for (prior in c("d_beta", "beta", "Moment")) {

    args <- list(
      alternative = alt,
      h0 = 0,
      ROPE = rope_by_alt[[alt]],
      N = 46,
      threshold = 3,
      prior_analysis = prior
    )

    if (prior == "d_beta") {
      args$k <- 1
    }

    if (prior == "beta") {
      args$alpha <- 1
      args$beta <- 1
    }

    if (prior == "Moment") {
      args$scale <- 0.707
    }

    run_test(
      paste("BFpower.cor | fixed N | interval-null |", alt, "|", prior),
      as.call(c(quote(BFpower.cor), args))
    )
  }
}


############################################################
## B5. Non-zero h0 power tests
############################################################

run_test(
  "BFpower.cor | fixed N | h0 = 0.1 | two.sided | interval-null",
  quote(
    BFpower.cor(
      alternative = "two.sided",
      h0 = 0.1,
      ROPE = c(-0.2, 0.2),
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1
    )
  )
)

run_test(
  "BFpower.cor | sample size | h0 = 0.1 | greater | interval-null",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0.1,
      ROPE = 0.2,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "Point",
      location_d = 0.3
    )
  )
)

run_test(
  "BFpower.cor | sample size | h0 = 0.1 | less | interval-null",
  quote(
    BFpower.cor(
      alternative = "less",
      h0 = 0.1,
      ROPE = -0.2,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "Moment",
      scale = 0.707,
      prior_design = "Point",
      location_d = -0.2
    )
  )
)



############################################################
## B5b. Additional non-zero h0 BFpower.cor tests
############################################################

for (h0_val in c(-0.3, 0.1, 0.3)) {
  for (alt in c("two.sided", "greater", "less")) {
    for (prior in c("d_beta", "beta", "Moment")) {

      args <- list(
        alternative = alt,
        h0 = h0_val,
        N = 46,
        threshold = 3,
        prior_analysis = prior
      )

      if (prior == "d_beta") {
        args$k <- 1
      }

      if (prior == "beta") {
        args$alpha <- 1
        args$beta <- 1
      }

      if (prior == "Moment") {
        args$scale <- 0.707
      }

      run_test(
        paste("BFpower.cor | fixed N | h0 =", h0_val, "|", alt, "|", prior),
        as.call(c(quote(BFpower.cor), args))
      )
    }
  }
}

############################################################
## B5c. Additional non-zero h0 BFpower.cor tests with ROPE
############################################################

for (h0_val in c(-0.3, 0.1, 0.3)) {
  for (alt in c("two.sided", "greater", "less")) {
    for (prior in c("d_beta", "beta", "Moment")) {

      args <- list(
        alternative = alt,
        h0 = h0_val,
        ROPE = rope_by_alt[[alt]],
        N = 46,
        threshold = 3,
        prior_analysis = prior
      )

      if (prior == "d_beta") {
        args$k <- 1
      }

      if (prior == "beta") {
        args$alpha <- 1
        args$beta <- 1
      }

      if (prior == "Moment") {
        args$scale <- 0.707
      }

      run_test(
        paste("BFpower.cor | fixed N | interval-null | h0 =", h0_val, "|", alt, "|", prior),
        as.call(c(quote(BFpower.cor), args))
      )
    }
  }
}
############################################################
## B6. Type-rate tests
############################################################

for (tr in c("positive", "negative")) {
  run_test(
    paste("BFpower.cor | sample size | type_rate =", tr),
    quote(
      BFpower.cor(
        alternative = "greater",
        h0 = 0,
        threshold = 3,
        true_rate = 0.8,
        false_rate = 0.05,
        type_rate = tr,
        prior_analysis = "d_beta",
        k = 1,
        prior_design = "Point",
        location_d = 0.3
      )
    )
  )
}


############################################################
## B7. Design prior tests in fixed-sample mode
############################################################

run_test(
  "BFpower.cor | fixed N | design prior d_beta",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "d_beta",
      k_d = 1
    )
  )
)

run_test(
  "BFpower.cor | fixed N | design prior beta",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "beta",
      alpha_d = 1,
      beta_d = 1
    )
  )
)

run_test(
  "BFpower.cor | fixed N | design prior Moment",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "Moment",
      scale_d = 0.707,location_d=0
    )
  )
)

run_test(
  "BFpower.cor | fixed N | design prior Point",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "Point",
      location_d = 0.3
    )
  )
)


############################################################
## B8. Design prior tests in sample-size mode
############################################################

run_test(
  "BFpower.cor | sample size | design prior d_beta",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "d_beta",
      k_d = 1
    )
  )
)

run_test(
  "BFpower.cor | sample size | design prior beta",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "beta",
      alpha_d = 1,
      beta_d = 1
    )
  )
)

run_test(
  "BFpower.cor | sample size | design prior Moment",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "Moment",
      scale_d = 0.707,location_d=.1
    )
  )
)

run_test(
  "BFpower.cor | sample size | design prior Point",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "Point",
      location_d = 0.3
    )
  )
)


############################################################
## B9. BFpower.cor object structure test
############################################################

test_that("BFpower.cor returns the expected object structure", {
  bfpower_cor_obj <- BFpower.cor(
    alternative = "greater",
    h0 = 0,
    threshold = 3,
    true_rate = 0.8,
    false_rate = 0.05,
    type_rate = "positive",
    prior_analysis = "d_beta",
    k = 1,
    prior_design = "Point",
    location_d = 0.3
  )

  expect_s3_class(bfpower_cor_obj, "BFpower")

  expect_true("type" %in% names(bfpower_cor_obj))
  expect_true("alternative" %in% names(bfpower_cor_obj))
  expect_true("h0" %in% names(bfpower_cor_obj))
  expect_true("ROPE" %in% names(bfpower_cor_obj))
  expect_true("analysis_h1" %in% names(bfpower_cor_obj))
  expect_true("design_h1" %in% names(bfpower_cor_obj))
  expect_true("results" %in% names(bfpower_cor_obj))
  expect_true("threshold" %in% names(bfpower_cor_obj))
  expect_true("setting" %in% names(bfpower_cor_obj))
})


############################################################
## B10. BFpower.cor validation/error tests
############################################################

run_test(
  "ERROR BFpower.cor | invalid alternative",
  quote(
    BFpower.cor(
      alternative = "wrong",
      h0 = 0,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid h0",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0.9,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid threshold",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      threshold = 0.5,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid true_rate",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      threshold = 3,
      true_rate = 0.5,
      false_rate = 0.05,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid false_rate",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.2,
      type_rate = "positive",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid type_rate",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "wrong",
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid N",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 3,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid prior_analysis",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "Normal",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid analysis k",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid analysis alpha",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "beta",
      alpha = 0,
      beta = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid analysis beta",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "beta",
      alpha = 1,
      beta = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid analysis scale",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "Moment",
      scale = 0
    )
  ),
  expect_error = TRUE
)


############################################################
## B11. BFpower.cor ROPE validation/error tests
############################################################

run_test(
  "ERROR BFpower.cor | two.sided ROPE length 1",
  quote(
    BFpower.cor(
      alternative = "two.sided",
      h0 = 0,
      ROPE = 0.2,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | two.sided ROPE both positive",
  quote(
    BFpower.cor(
      alternative = "two.sided",
      h0 = 0,
      ROPE = c(0.1, 0.2),
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | two.sided ROPE both negative",
  quote(
    BFpower.cor(
      alternative = "two.sided",
      h0 = 0,
      ROPE = c(-0.3, -0.2),
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | two.sided ROPE reversed signs",
  quote(
    BFpower.cor(
      alternative = "two.sided",
      h0 = 0,
      ROPE = c(0.2, -0.2),
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | greater ROPE negative",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      ROPE = -0.2,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | less ROPE positive",
  quote(
    BFpower.cor(
      alternative = "less",
      h0 = 0,
      ROPE = 0.2,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1
    )
  ),
  expect_error = TRUE
)


############################################################
## B12. BFpower.cor design prior validation/error tests
############################################################

run_test(
  "ERROR BFpower.cor | invalid prior_design",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "Normal",
      location_d = 0.3
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid design k_d",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "d_beta",
      k_d = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid design alpha_d",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "beta",
      alpha_d = 0,
      beta_d = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid design beta_d",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "beta",
      alpha_d = 1,
      beta_d = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid design scale_d",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "Moment",
      scale_d = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.cor | invalid design Point location_d outside correlation range",
  quote(
    BFpower.cor(
      alternative = "greater",
      h0 = 0,
      N = 46,
      threshold = 3,
      prior_analysis = "d_beta",
      k = 1,
      prior_design = "Point",
      location_d = 1.2
    )
  ),
  expect_error = TRUE
)


############################################################
## B13. Boundary-value tests
############################################################

run_test(
  "BF10.cor | boundary r close to 1",
  quote(
    BF10.cor(
      r = 0.999,
      n = 46,
      h0 = 0,
      alternative = "greater",
      prior_analysis = "d_beta",
      k = 1
    )
  )
)

run_test(
  "BF10.cor | boundary r close to -1",
  quote(
    BF10.cor(
      r = -0.999,
      n = 46,
      h0 = 0,
      alternative = "less",
      prior_analysis = "d_beta",
      k = 1
    )
  )
)

run_test(
  "BF10.cor | h0 upper allowed boundary",
  quote(
    BF10.cor(
      r = 0.5,
      n = 46,
      h0 = 0.8,
      alternative = "less",
      prior_analysis = "d_beta",
      k = 1
    )
  )
)

run_test(
  "BF10.cor | h0 lower allowed boundary",
  quote(
    BF10.cor(
      r = -0.5,
      n = 46,
      h0 = -0.8,
      alternative = "greater",
      prior_analysis = "d_beta",
      k = 1
    )
  )
)


############################################################
## End
############################################################

