############################################################
## Comprehensive test script for:
## 1) BF10.ttest.TwoSample()
## 2) BFpower.ttest.TwoSample()
############################################################

## Load your package first, e.g.
## devtools::load_all()
## or
## library(BayesPower)


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
## PART A: BF10.ttest.TwoSample()
############################################################

############################################################
## A1. Point-null BF10 tests: ROPE = NULL
############################################################

for (alt in c("two.sided", "less", "greater")) {
  for (prior in c("t-distribution", "Normal", "Moment")) {

    args <- list(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = prior,
      location = 0,
      scale = 0.707,
      alternative = alt
    )

    if (prior == "t-distribution") {
      args$dff <- 1
    }

    run_test(
      paste("BF10 TwoSample | point-null |", alt, "|", prior),
      as.call(c(quote(BF10.ttest.TwoSample), args))
    )
  }
}


############################################################
## A2. Interval-null BF10 tests: ROPE supplied
############################################################

rope_by_alt <- list(
  "two.sided" = c(-0.36, 0.36),
  "less" = -0.36,
  "greater" = 0.36
)

for (alt in c("two.sided", "less", "greater")) {
  for (prior in c("t-distribution", "Normal", "Moment")) {

    args <- list(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = prior,
      location = 0,
      scale = 0.707,
      alternative = alt,
      ROPE = rope_by_alt[[alt]]
    )

    if (prior == "t-distribution") {
      args$dff <- 1
    }

    run_test(
      paste("BF10 TwoSample | interval-null |", alt, "|", prior),
      as.call(c(quote(BF10.ttest.TwoSample), args))
    )
  }
}


############################################################
## A3. BF10 object structure test
############################################################

test_that("BF10.ttest.TwoSample returns the expected object structure", {
  bf10_obj <- BF10.ttest.TwoSample(
    tval = -1.148,
    N1 = 53,
    N2 = 48,
    prior_analysis = "t-distribution",
    location = 0,
    scale = 0.707,
    dff = 1,
    alternative = "two.sided",
    ROPE = c(-0.36, 0.36)
  )

  expect_s3_class(bf10_obj, "BFvalue")

  expect_true("type" %in% names(bf10_obj))
  expect_true("bf10" %in% names(bf10_obj))
  expect_true("tval" %in% names(bf10_obj))
  expect_true("df" %in% names(bf10_obj))
  expect_true("analysis_h1" %in% names(bf10_obj))
  expect_true("alternative" %in% names(bf10_obj))
  expect_true("ROPE" %in% names(bf10_obj))
  expect_true("N1" %in% names(bf10_obj))
  expect_true("N2" %in% names(bf10_obj))
  expect_true("d" %in% names(bf10_obj))
  expect_true("p.value" %in% names(bf10_obj))
})


############################################################
## A4. BF10 validation/error tests
############################################################

run_test(
  "ERROR BF10 | invalid tval",
  quote(
    BF10.ttest.TwoSample(
      tval = NA,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | invalid N1",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 2,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | invalid N2",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 2,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | invalid prior_analysis",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Uniform",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | invalid location",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = NA,
      scale = 0.707,
      alternative = "two.sided"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | invalid scale",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0,
      alternative = "two.sided"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | missing dff for t-distribution",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | invalid dff",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 0,
      alternative = "two.sided"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | invalid alternative",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "wrong"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | two.sided ROPE length 1",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = 0.36
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | two.sided ROPE both positive",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = c(0.10, 0.36)
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | two.sided ROPE both negative",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = c(-0.36, -0.10)
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | two.sided ROPE reversed signs",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = c(0.36, -0.36)
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | greater ROPE negative",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "greater",
      ROPE = -0.36
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BF10 | less ROPE positive",
  quote(
    BF10.ttest.TwoSample(
      tval = -1.148,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "less",
      ROPE = 0.36
    )
  ),
  expect_error = TRUE
)


############################################################
## PART B: BFpower.ttest.TwoSample()
############################################################

############################################################
## B1. Sample size determination mode
## N1 = NULL, N2 = NULL, r supplied
############################################################

for (alt in c("two.sided", "less", "greater")) {
  for (prior in c("t-distribution", "Normal", "Moment")) {

    args <- list(
      alternative = alt,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = prior,
      location = 0,
      scale = 0.707
    )

    if (prior == "t-distribution") {
      args$dff <- 1
    }

    run_test(
      paste("BFpower TwoSample | sample size | point-null |", alt, "|", prior),
      as.call(c(quote(BFpower.ttest.TwoSample), args))
    )
  }
}


############################################################
## B2. Sample size determination mode with ROPE
############################################################

for (alt in c("two.sided", "less", "greater")) {
  for (prior in c("t-distribution", "Normal", "Moment")) {

    args <- list(
      alternative = alt,
      ROPE = rope_by_alt[[alt]],
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = prior,
      location = 0,
      scale = 0.707
    )

    if (prior == "t-distribution") {
      args$dff <- 1
    }

    run_test(
      paste("BFpower TwoSample | sample size | interval-null |", alt, "|", prior),
      as.call(c(quote(BFpower.ttest.TwoSample), args))
    )
  }
}


############################################################
## B3. Sample size determination with different allocation ratios
############################################################

for (ratio in c(0.5, 1, 2, 3)) {
  run_test(
    paste("BFpower TwoSample | sample size | ratio r =", ratio),
    quote(
      BFpower.ttest.TwoSample(
        alternative = "two.sided",
        threshold = 3,
        true_rate = 0.8,
        false_rate = 0.05,
        r = ratio,
        prior_analysis = "Normal",
        location = 0,
        scale = 0.707
      )
    )
  )
}


############################################################
## B4. Fixed-sample power mode
## N1 and N2 supplied
############################################################

for (alt in c("two.sided", "less", "greater")) {
  for (prior in c("t-distribution", "Normal", "Moment")) {

    args <- list(
      alternative = alt,
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = prior,
      location = 0,
      scale = 0.707
    )

    if (prior == "t-distribution") {
      args$dff <- 1
    }

    run_test(
      paste("BFpower TwoSample | fixed N | point-null |", alt, "|", prior),
      as.call(c(quote(BFpower.ttest.TwoSample), args))
    )
  }
}


############################################################
## B5. Fixed-sample power mode with ROPE
############################################################

for (alt in c("two.sided", "less", "greater")) {
  for (prior in c("t-distribution", "Normal", "Moment")) {

    args <- list(
      alternative = alt,
      ROPE = rope_by_alt[[alt]],
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = prior,
      location = 0,
      scale = 0.707
    )

    if (prior == "t-distribution") {
      args$dff <- 1
    }

    run_test(
      paste("BFpower TwoSample | fixed N | interval-null |", alt, "|", prior),
      as.call(c(quote(BFpower.ttest.TwoSample), args))
    )
  }
}


############################################################
## B6. Type-rate tests
############################################################

for (type_rate in c("positive", "negative")) {
  run_test(
    paste("BFpower TwoSample | sample size | type_rate =", type_rate),
    quote(
      BFpower.ttest.TwoSample(
        alternative = "two.sided",
        threshold = 3,
        true_rate = 0.8,
        false_rate = 0.05,
        type_rate = type_rate,
        r = 1,
        prior_analysis = "Normal",
        location = 0,
        scale = 0.707
      )
    )
  )
}


############################################################
## B7. Design prior tests in fixed-sample mode
############################################################

run_test(
  "BFpower TwoSample | design prior Normal",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1,
      prior_design = "Normal",
      location_d = 0,
      scale_d = 0.707
    )
  )
)

run_test(
  "BFpower TwoSample | design prior Moment",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "Moment",
      location_d = 0,
      scale_d = 0.707
    )
  )
)

run_test(
  "BFpower TwoSample | design prior t-distribution",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "t-distribution",
      location_d = 0,
      scale_d = 0.707,
      dff_d = 1
    )
  )
)

run_test(
  "BFpower TwoSample | design prior Point",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "Point",
      location_d = 0.5
    )
  )
)


############################################################
## B8. Design prior tests in sample-size mode
############################################################

run_test(
  "BFpower TwoSample | sample size | design prior Normal",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "Normal",
      location_d = 0,
      scale_d = 0.707
    )
  )
)

run_test(
  "BFpower TwoSample | sample size | design prior Moment",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "Moment",
      location_d = 0,
      scale_d = 0.707
    )
  )
)

run_test(
  "BFpower TwoSample | sample size | design prior t-distribution",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "t-distribution",
      location_d = 0,
      scale_d = 0.707,
      dff_d = 1
    )
  )
)

run_test(
  "BFpower TwoSample | sample size | design prior Point",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "Point",
      location_d = 0.5
    )
  )
)


############################################################
## B9. BFpower object structure test
############################################################

test_that("BFpower.ttest.TwoSample returns the expected object structure", {
  bfpower_obj <- BFpower.ttest.TwoSample(
    alternative = "two.sided",
    threshold = 3,
    N1 = 53,
    N2 = 48,
    prior_analysis = "t-distribution",
    location = 0,
    scale = 0.707,
    dff = 1,
    ROPE = c(-0.36, 0.36)
  )

  expect_s3_class(bfpower_obj, "BFpower")

  expect_true("type" %in% names(bfpower_obj))
  expect_true("alternative" %in% names(bfpower_obj))
  expect_true("ROPE" %in% names(bfpower_obj))
  expect_true("analysis_h1" %in% names(bfpower_obj))
  expect_true("design_h1" %in% names(bfpower_obj))
  expect_true("results" %in% names(bfpower_obj))
  expect_true("threshold" %in% names(bfpower_obj))
  expect_true("setting" %in% names(bfpower_obj))
})

############################################################
## B10. BFpower validation/error tests
############################################################

run_test(
  "ERROR BFpower | invalid alternative",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "wrong",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid prior_analysis",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "Uniform",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid scale",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | missing dff for t-distribution",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid dff",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid threshold",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 0.5,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid true_rate",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.5,
      false_rate = 0.05,
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid false_rate",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.2,
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid type_rate",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      type_rate = "wrong",
      r = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid r",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      r = 0,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid N1",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 1,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid N2",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 53,
      N2 = 1,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)


############################################################
## B11. BFpower ROPE validation/error tests
############################################################

run_test(
  "ERROR BFpower | two.sided ROPE length 1",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      ROPE = 0.36,
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | two.sided ROPE both positive",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      ROPE = c(0.10, 0.36),
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | two.sided ROPE both negative",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      ROPE = c(-0.36, -0.10),
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | two.sided ROPE reversed signs",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      ROPE = c(0.36, -0.36),
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | greater ROPE negative",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "greater",
      ROPE = -0.36,
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | less ROPE positive",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "less",
      ROPE = 0.36,
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  ),
  expect_error = TRUE
)


############################################################
## B12. Design prior validation/error tests
############################################################

run_test(
  "ERROR BFpower | invalid prior_design",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "Uniform",
      location_d = 0,
      scale_d = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid scale_d",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "Normal",
      location_d = 0,
      scale_d = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | missing dff_d for design t-distribution",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "t-distribution",
      location_d = 0,
      scale_d = 0.707
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower | invalid dff_d for design t-distribution",
  quote(
    BFpower.ttest.TwoSample(
      alternative = "two.sided",
      threshold = 3,
      N1 = 53,
      N2 = 48,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "t-distribution",
      location_d = 0,
      scale_d = 0.707,
      dff_d = 0
    )
  ),
  expect_error = TRUE
)


############################################################
## End
############################################################

