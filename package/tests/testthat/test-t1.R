############################################################
## Exhaustive test script for BF10.ttest.OneSample()
############################################################

## Assumption:
## Your package / source file is already loaded, e.g.
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
## 1. POINT-NULL TESTS
## alternative = "two.sided"
############################################################

run_test(
  "BF10.ttest.OneSample | two.sided | point-null | t-distribution",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1,
      alternative = "two.sided"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | two.sided | point-null | Normal",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | two.sided | point-null | Moment",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  )
)


############################################################
## 2. INTERVAL-NULL TESTS
## alternative = "two.sided"
############################################################

run_test(
  "BF10.ttest.OneSample | two.sided | interval-null | t-distribution",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1,
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2)
    )
  )
)

run_test(
  "BF10.ttest.OneSample | two.sided | interval-null | Normal",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2)
    )
  )
)

run_test(
  "BF10.ttest.OneSample | two.sided | interval-null | Moment",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2)
    )
  )
)


############################################################
## 3. POINT-NULL TESTS
## alternative = "less"
############################################################

run_test(
  "BF10.ttest.OneSample | less | point-null | t-distribution",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1,
      alternative = "less"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | point-null | Normal",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "less"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | point-null | Moment",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707,
      alternative = "less"
    )
  )
)


############################################################
## 4. INTERVAL-NULL TESTS
## alternative = "less"
############################################################

run_test(
  "BF10.ttest.OneSample | less | interval-null | t-distribution",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1,
      alternative = "less",
      ROPE = -0.2
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | interval-null | Normal",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "less",
      ROPE = -0.2
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | interval-null | Moment",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707,
      alternative = "less",
      ROPE = -0.2
    )
  )
)


############################################################
## 5. POINT-NULL TESTS
## alternative = "greater"
############################################################

run_test(
  "BF10.ttest.OneSample | greater | point-null | t-distribution",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1,
      alternative = "greater"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | point-null | Normal",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "greater"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | point-null | Moment",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707,
      alternative = "greater"
    )
  )
)


############################################################
## 6. INTERVAL-NULL TESTS
## alternative = "greater"
############################################################

run_test(
  "BF10.ttest.OneSample | greater | interval-null | t-distribution",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1,
      alternative = "greater",
      ROPE = 0.2
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | interval-null | Normal",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "greater",
      ROPE = 0.2
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | interval-null | Moment",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707,
      alternative = "greater",
      ROPE = 0.2
    )
  )
)


############################################################
## 7. ADDITIONAL SCALE = 0.71 TESTS
## Point-null hypothesis
############################################################

run_test(
  "BF10.ttest.OneSample | two.sided | point-null | t-distribution | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.71,
      dff = 1,
      alternative = "two.sided"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | two.sided | point-null | Normal | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.71,
      alternative = "two.sided"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | two.sided | point-null | Moment | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.71,
      alternative = "two.sided"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | point-null | t-distribution | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.71,
      dff = 1,
      alternative = "greater"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | point-null | Normal | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.71,
      alternative = "greater"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | point-null | Moment | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.71,
      alternative = "greater"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | point-null | t-distribution | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.71,
      dff = 1,
      alternative = "less"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | point-null | Normal | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.71,
      alternative = "less"
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | point-null | Moment | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.71,
      alternative = "less"
    )
  )
)


############################################################
## 8. ADDITIONAL SCALE = 0.71 TESTS
## Interval-null hypothesis
############################################################

run_test(
  "BF10.ttest.OneSample | two.sided | interval-null | t-distribution | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.71,
      dff = 1,
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2)
    )
  )
)

run_test(
  "BF10.ttest.OneSample | two.sided | interval-null | Normal | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.71,
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2)
    )
  )
)

run_test(
  "BF10.ttest.OneSample | two.sided | interval-null | Moment | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.71,
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2)
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | interval-null | t-distribution | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.71,
      dff = 1,
      alternative = "greater",
      ROPE = 0.2
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | interval-null | Normal | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.71,
      alternative = "greater",
      ROPE = 0.2
    )
  )
)

run_test(
  "BF10.ttest.OneSample | greater | interval-null | Moment | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.71,
      alternative = "greater",
      ROPE = 0.2
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | interval-null | t-distribution | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.71,
      dff = 1,
      alternative = "less",
      ROPE = -0.2
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | interval-null | Normal | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.71,
      alternative = "less",
      ROPE = -0.2
    )
  )
)

run_test(
  "BF10.ttest.OneSample | less | interval-null | Moment | scale 0.71",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.71,
      alternative = "less",
      ROPE = -0.2
    )
  )
)


############################################################
## 9. VALIDATION TESTS
## These should give errors
############################################################

run_test(
  "ERROR TEST | BF10.ttest.OneSample | invalid alternative",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "invalid"
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | invalid prior_analysis",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Uniform",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | invalid scale",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0,
      alternative = "two.sided"
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | missing dff for t-distribution",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | invalid dff",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 0,
      alternative = "two.sided"
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | invalid tval",
  quote(
    BF10.ttest.OneSample(
      tval = NA,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | invalid df",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 0,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided"
    )
  )
)


############################################################
## 10. ROPE VALIDATION TESTS
## These should give errors
############################################################

run_test(
  "ERROR TEST | BF10.ttest.OneSample | two.sided ROPE length 1",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = 0.2
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | two.sided ROPE both positive",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = c(0.1, 0.2)
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | two.sided ROPE both negative",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = c(-0.3, -0.2)
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | two.sided ROPE reversed signs",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "two.sided",
      ROPE = c(0.2, -0.2)
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | greater ROPE negative",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "greater",
      ROPE = -0.2
    )
  )
)

run_test(
  "ERROR TEST | BF10.ttest.OneSample | less ROPE positive",
  quote(
    BF10.ttest.OneSample(
      tval = 2,
      df = 50,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      alternative = "less",
      ROPE = 0.2
    )
  )
)


############################################################
## 11. OBJECT STRUCTURE TEST
############################################################

testthat::test_that("BF10.ttest.OneSample returns the expected object structure", {
  obj <- BF10.ttest.OneSample(
    tval = 2,
    df = 50,
    prior_analysis = "t-distribution",
    location = 0,
    scale = 0.707,
    dff = 1,
    alternative = "two.sided"
  )

  testthat::expect_s3_class(obj, "BFvalue")

  testthat::expect_true("type" %in% names(obj))
  testthat::expect_true("bf10" %in% names(obj))
  testthat::expect_true("tval" %in% names(obj))
  testthat::expect_true("df" %in% names(obj))
  testthat::expect_true("analysis_h1" %in% names(obj))
  testthat::expect_true("alternative" %in% names(obj))
  testthat::expect_true("ROPE" %in% names(obj))
  testthat::expect_true("p.value" %in% names(obj))
})


############################################################
## End of test script
############################################################


############################################################
## 1. SAMPLE SIZE DETERMINATION
## Point-null hypothesis: ROPE = NULL
############################################################

## two.sided + t-distribution
run_test(
  "Sample size determination | two.sided | point-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

## two.sided + Normal
run_test(
  "Sample size determination | two.sided | point-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## two.sided + Moment
run_test(
  "Sample size determination | two.sided | point-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 2. SAMPLE SIZE DETERMINATION
## Interval-null hypothesis: ROPE is supplied
############################################################

## two.sided + t-distribution + ROPE
run_test(
  "Sample size determination | two.sided | interval-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2),
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

## two.sided + Normal + ROPE
run_test(
  "Sample size determination | two.sided | interval-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2),
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## two.sided + Moment + ROPE
run_test(
  "Sample size determination | two.sided | interval-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2),
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 3. ONE-SIDED SAMPLE SIZE DETERMINATION
## alternative = "less"
############################################################

## less + t-distribution
run_test(
  "Sample size determination | less | point-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

## less + Normal
run_test(
  "Sample size determination | less | point-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## less + Moment
run_test(
  "Sample size determination | less | point-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)

## less + ROPE + t-distribution
run_test(
  "Sample size determination | less | interval-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      ROPE = -0.2,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

## less + ROPE + Normal
run_test(
  "Sample size determination | less | interval-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      ROPE = -0.2,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## less + ROPE + Moment
run_test(
  "Sample size determination | less | interval-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      ROPE = -0.2,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 4. ONE-SIDED SAMPLE SIZE DETERMINATION
## alternative = "greater"
############################################################

## greater + t-distribution
run_test(
  "Sample size determination | greater | point-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

## greater + Normal
run_test(
  "Sample size determination | greater | point-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## greater + Moment
run_test(
  "Sample size determination | greater | point-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)

## greater + ROPE + t-distribution
run_test(
  "Sample size determination | greater | interval-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      ROPE = 0.2,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

## greater + ROPE + Normal
run_test(
  "Sample size determination | greater | interval-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      ROPE = 0.2,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## greater + ROPE + Moment
run_test(
  "Sample size determination | greater | interval-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      ROPE = 0.2,
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 5. FIXED-N POWER CALCULATION
## Point-null hypothesis
############################################################

## two.sided fixed N + t-distribution
run_test(
  "Fixed N | two.sided | point-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      N = 50,
      threshold = 3,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

## two.sided fixed N + Normal
run_test(
  "Fixed N | two.sided | point-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## two.sided fixed N + Moment
run_test(
  "Fixed N | two.sided | point-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      N = 50,
      threshold = 3,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 6. FIXED-N POWER CALCULATION
## Interval-null hypothesis
############################################################

## two.sided fixed N + t-distribution + ROPE
run_test(
  "Fixed N | two.sided | interval-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2),
      N = 50,
      threshold = 3,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

## two.sided fixed N + Normal + ROPE
run_test(
  "Fixed N | two.sided | interval-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2),
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## two.sided fixed N + Moment + ROPE
run_test(
  "Fixed N | two.sided | interval-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = c(-0.2, 0.2),
      N = 50,
      threshold = 3,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 7. FIXED-N ONE-SIDED TESTS
############################################################

## less fixed N
run_test(
  "Fixed N | less | point-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      N = 50,
      threshold = 3,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

run_test(
  "Fixed N | less | point-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

run_test(
  "Fixed N | less | point-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      N = 50,
      threshold = 3,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)

## less fixed N with ROPE
run_test(
  "Fixed N | less | interval-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      ROPE = -0.2,
      N = 50,
      threshold = 3,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

run_test(
  "Fixed N | less | interval-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      ROPE = -0.2,
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

run_test(
  "Fixed N | less | interval-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      ROPE = -0.2,
      N = 50,
      threshold = 3,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)

## greater fixed N
run_test(
  "Fixed N | greater | point-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      N = 50,
      threshold = 3,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

run_test(
  "Fixed N | greater | point-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

run_test(
  "Fixed N | greater | point-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      N = 50,
      threshold = 3,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)

## greater fixed N with ROPE
run_test(
  "Fixed N | greater | interval-null | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      ROPE = 0.2,
      N = 50,
      threshold = 3,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 1
    )
  )
)

run_test(
  "Fixed N | greater | interval-null | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      ROPE = 0.2,
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

run_test(
  "Fixed N | greater | interval-null | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      ROPE = 0.2,
      N = 50,
      threshold = 3,
      prior_analysis = "Moment",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 8. DESIGN PRIOR TESTS
############################################################

## Design prior = Normal
run_test(
  "Design prior | Normal",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      N = 50,
      threshold = 3,
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

## Design prior = Moment
run_test(
  "Design prior | Moment",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "Moment",
      location_d = 0,
      scale_d = 0.707
    )
  )
)

## Design prior = t-distribution
run_test(
  "Design prior | t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      N = 50,
      threshold = 3,
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

## Design prior = Point
run_test(
  "Design prior | Point",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707,
      prior_design = "Point",
      location_d = 0.5
    )
  )
)


############################################################
## 9. TYPE RATE TESTS
############################################################

## positive rate
run_test(
  "Sample size determination | type_rate = positive",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      type_rate = "positive",
      true_rate = 0.8,
      false_rate = 0.05,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## negative rate
run_test(
  "Sample size determination | type_rate = negative",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      type_rate = "negative",
      true_rate = 0.8,
      false_rate = 0.05,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 10. VALIDATION TESTS
## These should give errors
############################################################

## invalid alternative
run_test(
  "ERROR TEST | invalid alternative",
  quote(
    BFpower.ttest.OneSample(
      alternative = "invalid",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## invalid prior_analysis
run_test(
  "ERROR TEST | invalid prior_analysis",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Uniform",
      location = 0,
      scale = 0.707
    )
  )
)

## invalid scale
run_test(
  "ERROR TEST | scale <= 0",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0
    )
  )
)

## missing dff for t-distribution
run_test(
  "ERROR TEST | missing dff for t-distribution",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707
    )
  )
)

## invalid dff
run_test(
  "ERROR TEST | invalid dff",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "t-distribution",
      location = 0,
      scale = 0.707,
      dff = 0
    )
  )
)

## invalid threshold
run_test(
  "ERROR TEST | threshold < 1",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 0.5,
      true_rate = 0.8,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## invalid true_rate
run_test(
  "ERROR TEST | true_rate outside valid range",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.5,
      false_rate = 0.05,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## invalid false_rate
run_test(
  "ERROR TEST | false_rate outside valid range",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      threshold = 3,
      true_rate = 0.8,
      false_rate = 0.2,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## invalid N
run_test(
  "ERROR TEST | invalid N",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      N = 0,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 11. ROPE VALIDATION TESTS
## These should give errors
############################################################

## two.sided ROPE not length 2
run_test(
  "ERROR TEST | two.sided ROPE length 1",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = 0.2,
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## two.sided ROPE both positive
run_test(
  "ERROR TEST | two.sided ROPE both positive",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = c(0.1, 0.2),
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## two.sided ROPE both negative
run_test(
  "ERROR TEST | two.sided ROPE both negative",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = c(-0.3, -0.2),
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## two.sided ROPE reversed
run_test(
  "ERROR TEST | two.sided ROPE reversed signs",
  quote(
    BFpower.ttest.OneSample(
      alternative = "two.sided",
      ROPE = c(0.2, -0.2),
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## greater ROPE negative
run_test(
  "ERROR TEST | greater ROPE negative",
  quote(
    BFpower.ttest.OneSample(
      alternative = "greater",
      ROPE = -0.2,
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)

## less ROPE positive
run_test(
  "ERROR TEST | less ROPE positive",
  quote(
    BFpower.ttest.OneSample(
      alternative = "less",
      ROPE = 0.2,
      N = 50,
      threshold = 3,
      prior_analysis = "Normal",
      location = 0,
      scale = 0.707
    )
  )
)


############################################################
## 12. OBJECT STRUCTURE TESTS
############################################################
test_that("BFpower.ttest.OneSample returns the expected object structure", {
  obj <- BFpower.ttest.OneSample(
    alternative = "two.sided",
    N = 50,
    threshold = 3,
    prior_analysis = "t-distribution",
    location = 0,
    scale = 0.707,
    dff = 1
  )

  expect_s3_class(obj, "BFpower")

  expect_true("type" %in% names(obj))
  expect_true("alternative" %in% names(obj))
  expect_true("analysis_h1" %in% names(obj))
  expect_true("results" %in% names(obj))
  expect_true("threshold" %in% names(obj))
  expect_true("setting" %in% names(obj))
})


############################################################
## End of test script
############################################################

