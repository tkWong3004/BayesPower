############################################################
## Comprehensive test script for:
## 1) BF10.f.test()
## 2) BFpower.f.test()
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
## PART A: BF10.f.test()
############################################################

############################################################
## A1. Point-null BF10 tests
############################################################

run_test(
  "BF10.f.test | point-null | effectsize prior",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  )
)

run_test(
  "BF10.f.test | point-null | Moment prior",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 3,
      f_m = 0.1,
      prior_analysis = "Moment"
    )
  )
)


############################################################
## A2. Interval-null BF10 tests
############################################################

run_test(
  "BF10.f.test | interval-null | effectsize prior",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize",
      ROPE = 0.2
    )
  )
)

run_test(
  "BF10.f.test | interval-null | Moment prior",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 3,
      f_m = 0.1,
      prior_analysis = "Moment",
      ROPE = 0.2
    )
  )
)


############################################################
## A3. Boundary valid BF10 tests
############################################################

run_test(
  "BF10.f.test | fval = 0 valid boundary",
  quote(
    BF10.f.test(
      fval = 0,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  )
)

run_test(
  "BF10.f.test | very small positive ROPE",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize",
      ROPE = 0.001
    )
  )
)


############################################################
## A4. BF10.f.test object structure
############################################################

test_that("BF10.f.test returns the expected object structure", {
  bf10_f_obj <- BF10.f.test(
    fval = 4.5,
    df1 = 2,
    df2 = 12,
    dff = 12,
    rscale = 0.707,
    f_m = 0.1,
    prior_analysis = "effectsize"
  )

  expect_s3_class(bf10_f_obj, "BFvalue")

  expect_true("type" %in% names(bf10_f_obj))
  expect_true("bf10" %in% names(bf10_f_obj))
  expect_true("fval" %in% names(bf10_f_obj))
  expect_true("df1" %in% names(bf10_f_obj))
  expect_true("df2" %in% names(bf10_f_obj))
  expect_true("analysis_h1" %in% names(bf10_f_obj))
  expect_true("ROPE" %in% names(bf10_f_obj))
  expect_true("p.value" %in% names(bf10_f_obj))
})


############################################################
## A5. BF10.f.test validation/error tests
############################################################

run_test(
  "ERROR BF10.f.test | invalid fval negative",
  quote(
    BF10.f.test(
      fval = -1,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid fval NA",
  quote(
    BF10.f.test(
      fval = NA,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid df1",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 0,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid df2",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 0,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | missing prior_analysis",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid prior_analysis",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "Normal"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | missing rscale for effectsize",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid rscale",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | missing dff",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid dff",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 0,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | Moment dff < 3",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 2,
      f_m = 0.1,
      prior_analysis = "Moment"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | missing f_m",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid f_m",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0,
      prior_analysis = "effectsize"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid ROPE zero",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize",
      ROPE = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid ROPE negative",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize",
      ROPE = -0.2
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.f.test | invalid ROPE vector length 2",
  quote(
    BF10.f.test(
      fval = 4.5,
      df1 = 2,
      df2 = 12,
      dff = 12,
      rscale = 0.707,
      f_m = 0.1,
      prior_analysis = "effectsize",
      ROPE = c(0.1, 0.2)
    )
  ),
  expect_error = TRUE
)


############################################################
## PART B: BFpower.f.test()
############################################################

############################################################
## B1. Sample size determination mode: point-null
############################################################

for (prior in c("effectsize", "Moment")) {

  args <- list(
    threshold = 3,
    true_rate = 0.8,
    false_rate = 0.05,
    p = 3,
    k = 4,
    prior_analysis = prior,
    dff = if (prior == "Moment") 3 else 12,
    f_m = 0.1,
    type_rate = "positive"
  )

  if (prior == "effectsize") {
    args$rscale <- 0.18
  }

  run_test(
    paste("BFpower.f.test | sample size | point-null |", prior),
    as.call(c(quote(BFpower.f.test), args))
  )
}


############################################################
## B2. Sample size determination mode: interval-null
############################################################

for (prior in c("effectsize", "Moment")) {

  args <- list(
    threshold = 3,
    true_rate = 0.8,
    false_rate = 0.05,
    p = 3,
    k = 4,
    prior_analysis = prior,
    dff = if (prior == "Moment") 3 else 12,
    f_m = 0.1,
    type_rate = "positive",
    ROPE = 0.2
  )

  if (prior == "effectsize") {
    args$rscale <- 0.18
  }

  run_test(
    paste("BFpower.f.test | sample size | interval-null |", prior),
    as.call(c(quote(BFpower.f.test), args))
  )
}


############################################################
## B3. Fixed-N power mode: point-null
############################################################

for (prior in c("effectsize", "Moment")) {

  args <- list(
    threshold = 3,
    p = 3,
    k = 4,
    prior_analysis = prior,
    dff = if (prior == "Moment") 3 else 12,
    f_m = 0.1,
    N = 50
  )

  if (prior == "effectsize") {
    args$rscale <- 0.18
  }

  run_test(
    paste("BFpower.f.test | fixed N | point-null |", prior),
    as.call(c(quote(BFpower.f.test), args))
  )
}


############################################################
## B4. Fixed-N power mode: interval-null
############################################################

for (prior in c("effectsize", "Moment")) {

  args <- list(
    threshold = 3,
    p = 3,
    k = 4,
    prior_analysis = prior,
    dff = if (prior == "Moment") 3 else 12,
    f_m = 0.1,
    N = 50,
    ROPE = 0.2
  )

  if (prior == "effectsize") {
    args$rscale <- 0.18
  }

  run_test(
    paste("BFpower.f.test | fixed N | interval-null |", prior),
    as.call(c(quote(BFpower.f.test), args))
  )
}


############################################################
## B5. Design prior tests in fixed-N mode
############################################################

run_test(
  "BFpower.f.test | fixed N | design prior effectsize",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "effectsize",
      dff_d = 12,
      rscale_d = 0.18,
      f_m_d = 0.1,
      N = 50
    )
  )
)

run_test(
  "BFpower.f.test | fixed N | design prior Moment",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "Moment",
      dff_d = 3,
      f_m_d = 0.1,
      N = 50
    )
  )
)

run_test(
  "BFpower.f.test | fixed N | design prior Point",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "Point",
      f_m_d = 0.1,
      N = 50
    )
  )
)


############################################################
## B6. Design prior tests in sample-size mode
############################################################

run_test(
  "BFpower.f.test | sample size | design prior effectsize",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "effectsize",
      dff_d = 12,
      rscale_d = 0.18,
      f_m_d = 0.1
    )
  )
)

run_test(
  "BFpower.f.test | sample size | design prior Moment",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "Moment",
      dff_d = 3,
      f_m_d = 0.1
    )
  )
)

run_test(
  "BFpower.f.test | sample size | design prior Point",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "Point",
      f_m_d = 0.1
    )
  )
)


############################################################
## B7. Type-rate tests
############################################################

for (type_rate in c("positive", "negative")) {
  run_test(
    paste("BFpower.f.test | sample size | type_rate =", type_rate),
    quote(
      BFpower.f.test(
        threshold = 3,
        true_rate = 0.8,
        false_rate = 0.05,
        p = 3,
        k = 4,
        prior_analysis = "effectsize",
        dff = 12,
        rscale = 0.18,
        f_m = 0.1,
        prior_design = "Point",
        f_m_d = 0.1,
        type_rate = type_rate
      )
    )
  )
}


############################################################
## B8. Object structure test
############################################################

test_that("BFpower.f.test returns the expected object structure", {
  bfpower_f_obj <- BFpower.f.test(
    threshold = 3,
    p = 3,
    k = 4,
    prior_analysis = "effectsize",
    dff = 12,
    rscale = 0.18,
    f_m = 0.1,
    prior_design = "Point",
    f_m_d = 0.1,
    N = 50
  )

  expect_s3_class(bfpower_f_obj, "BFpower")

  expect_true("type" %in% names(bfpower_f_obj))
  expect_true("k" %in% names(bfpower_f_obj))
  expect_true("p" %in% names(bfpower_f_obj))
  expect_true("ROPE" %in% names(bfpower_f_obj))
  expect_true("analysis_h1" %in% names(bfpower_f_obj))
  expect_true("design_h1" %in% names(bfpower_f_obj))
  expect_true("results" %in% names(bfpower_f_obj))
  expect_true("threshold" %in% names(bfpower_f_obj))
  expect_true("setting" %in% names(bfpower_f_obj))
})


############################################################
## B9. Validation/error tests
############################################################

run_test(
  "ERROR BFpower.f.test | invalid p missing",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid k missing",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | k <= p",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 4,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid prior_analysis",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "Normal",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | missing rscale for effectsize",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid rscale",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid dff",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 0,
      rscale = 0.18,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | Moment dff < 3",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "Moment",
      dff = 2,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid f_m",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid type_rate",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      type_rate = "wrong"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid true_rate",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.5,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid false_rate",
  quote(
    BFpower.f.test(
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.2,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid threshold in sample-size mode",
  quote(
    BFpower.f.test(
      threshold = .9,
      true_rate = 0.8,
      false_rate = 0.05,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1
    )
  ),
  expect_error = TRUE
)


############################################################
## B10. Design prior validation tests
############################################################

run_test(
  "ERROR BFpower.f.test | invalid prior_design",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "Normal",
      f_m_d = 0.1,
      N = 50
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid design rscale_d",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "effectsize",
      dff_d = 12,
      rscale_d = 0,
      f_m_d = 0.1,
      N = 50
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid design dff_d",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "effectsize",
      dff_d = 0,
      rscale_d = 0.18,
      f_m_d = 0.1,
      N = 50
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid design f_m_d",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "effectsize",
      dff_d = 12,
      rscale_d = 0.18,
      f_m_d = 0,
      N = 50
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | Moment design dff_d < 3",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "Moment",
      dff_d = 2,
      f_m_d = 0.1,
      N = 50
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | Point design invalid f_m_d",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      prior_design = "Point",
      f_m_d = 0,
      N = 50
    )
  ),
  expect_error = TRUE
)


############################################################
## B11. ROPE validation tests
############################################################

run_test(
  "ERROR BFpower.f.test | invalid ROPE zero",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      ROPE = 0,
      N = 50
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid ROPE negative",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      ROPE = -0.2,
      N = 50
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.f.test | invalid ROPE vector",
  quote(
    BFpower.f.test(
      threshold = 3,
      p = 3,
      k = 4,
      prior_analysis = "effectsize",
      dff = 12,
      rscale = 0.18,
      f_m = 0.1,
      ROPE = c(0.1, 0.2),
      N = 50
    )
  ),
  expect_error = TRUE
)


############################################################
## End
############################################################


