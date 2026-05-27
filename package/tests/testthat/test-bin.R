############################################################
## Comprehensive tests for:
## 1) BF10.bin.test()
## 2) BFpower.bin()
############################################################

## Load package/functions first:
## devtools::load_all()
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
## PART A: BF10.bin.test()
############################################################

############################################################
## A1. Valid point-null tests
############################################################

for (alt in c("two.sided", "greater", "less")) {
  run_test(
    paste("BF10.bin.test | point-null | beta |", alt),
    substitute(
      BF10.bin.test(
        x = 42,
        n = 52,
        h0 = 0.5,
        prior_analysis = "beta",
        alpha = 1,
        beta = 1,
        alternative = ALT
      ),
      list(ALT = alt)
    )
  )

  run_test(
    paste("BF10.bin.test | point-null | Moment |", alt),
    substitute(
      BF10.bin.test(
        x = 42,
        n = 52,
        h0 = 0.5,
        prior_analysis = "Moment",
        scale = 0.707,
        alternative = ALT
      ),
      list(ALT = alt)
    )
  )
}


############################################################
## A2. Valid interval-null tests
############################################################

rope_by_alt <- list(
  "two.sided" = c(-0.2, 0.2),
  "greater" = 0.2,
  "less" = -0.2
)

for (alt in c("two.sided", "greater", "less")) {
  run_test(
    paste("BF10.bin.test | interval-null | beta |", alt),
    substitute(
      BF10.bin.test(
        x = 42,
        n = 52,
        h0 = 0.5,
        ROPE = ROPE_VAL,
        prior_analysis = "beta",
        alpha = 1,
        beta = 1,
        alternative = ALT
      ),
      list(ALT = alt, ROPE_VAL = rope_by_alt[[alt]])
    )
  )

  run_test(
    paste("BF10.bin.test | interval-null | Moment |", alt),
    substitute(
      BF10.bin.test(
        x = 42,
        n = 52,
        h0 = 0.5,
        ROPE = ROPE_VAL,
        prior_analysis = "Moment",
        scale = 0.707,
        alternative = ALT
      ),
      list(ALT = alt, ROPE_VAL = rope_by_alt[[alt]])
    )
  )
}


############################################################
## A3. BF10 object structure
############################################################

bf10_obj <- BF10.bin.test(
  x = 42,
  n = 52,
  h0 = 0.5,
  prior_analysis = "beta",
  alpha = 1,
  beta = 1,
  alternative = "greater"
)

test_that("BF10.bin.test returns the expected object structure", {
  expect_s3_class(bf10_obj, "BFvalue")
  expect_equal(bf10_obj$type, "One-proportion")

  expect_named(
    bf10_obj,
    c("type", "bf10", "h0", "x", "n", "analysis_h1", "alternative", "ROPE", "p.value"),
    ignore.order = TRUE
  )
})


############################################################
## A4. BF10 input-validation tests
############################################################

run_test(
  "ERROR BF10.bin.test | invalid n = 0",
  quote(BF10.bin.test(x = 1, n = 0, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | invalid n non-integer",
  quote(BF10.bin.test(x = 1, n = 10.5, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | invalid x negative",
  quote(BF10.bin.test(x = -1, n = 10, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | invalid x non-integer",
  quote(BF10.bin.test(x = 1.5, n = 10, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | x > n",
  quote(BF10.bin.test(x = 11, n = 10, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | h0 too small",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.05, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | h0 too large",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.95, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | invalid alternative",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "wrong")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | missing prior_analysis",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)


run_test(
  "ERROR BF10.bin.test | invalid prior_analysis",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, prior_analysis = "d_beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | beta prior missing alpha",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, prior_analysis = "beta", beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | beta prior invalid alpha",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, prior_analysis = "beta", alpha = 0, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | beta prior missing beta",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, prior_analysis = "beta", alpha = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | beta prior invalid beta",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 0, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | Moment prior missing scale",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, prior_analysis = "Moment", alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | Moment prior invalid scale",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, prior_analysis = "Moment", scale = 0, alternative = "greater")),
  expect_error = TRUE
)


############################################################
## A5. BF10 ROPE validation tests
############################################################

run_test(
  "ERROR BF10.bin.test | two.sided ROPE length 1",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, ROPE = 0.2, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "two.sided")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | two.sided ROPE both positive",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, ROPE = c(0.1, 0.2), prior_analysis = "beta", alpha = 1, beta = 1, alternative = "two.sided")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | two.sided ROPE both negative",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, ROPE = c(-0.3, -0.2), prior_analysis = "beta", alpha = 1, beta = 1, alternative = "two.sided")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | two.sided ROPE reversed signs",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, ROPE = c(0.2, -0.2), prior_analysis = "beta", alpha = 1, beta = 1, alternative = "two.sided")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | two.sided ROPE outside width bound",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, ROPE = c(-0.6, 0.2), prior_analysis = "beta", alpha = 1, beta = 1, alternative = "two.sided")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | two.sided h0 + ROPE outside [0,1]",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.9, ROPE = c(-0.2, 0.2), prior_analysis = "beta", alpha = 1, beta = 1, alternative = "two.sided")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | greater ROPE negative",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, ROPE = -0.2, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | greater ROPE too large",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, ROPE = 0.6, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "greater")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | less ROPE positive",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, ROPE = 0.2, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "less")),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.bin.test | less ROPE too small",
  quote(BF10.bin.test(x = 5, n = 10, h0 = 0.5, ROPE = -0.6, prior_analysis = "beta", alpha = 1, beta = 1, alternative = "less")),
  expect_error = TRUE
)


############################################################
## PART B: BFpower.bin()
############################################################

############################################################
## B1. Sample-size mode, point-null
############################################################

for (alt in c("two.sided", "greater", "less")) {
  run_test(
    paste("BFpower.bin | sample-size | point-null | beta |", alt),
    substitute(
      BFpower.bin(
        alternative = ALT,
        threshold = 3,
        true_rate = 0.8,
        false_rate = 0.05,
        h0 = 0.5,
        prior_analysis = "beta",
        alpha = 1,
        beta = 1
      ),
      list(ALT = alt)
    )
  )

  run_test(
    paste("BFpower.bin | sample-size | point-null | Moment |", alt),
    substitute(
      BFpower.bin(
        alternative = ALT,
        threshold = 3,
        true_rate = 0.8,
        false_rate = 0.05,
        h0 = 0.5,
        prior_analysis = "Moment",
        scale = 0.707
      ),
      list(ALT = alt)
    )
  )
}


############################################################
## B2. Sample-size mode, interval-null
############################################################

for (alt in c("two.sided", "greater", "less")) {
  run_test(
    paste("BFpower.bin | sample-size | interval-null | beta |", alt),
    substitute(
      BFpower.bin(
        alternative = ALT,
        threshold = 3,
        true_rate = 0.8,
        false_rate = 0.05,
        h0 = 0.5,
        ROPE = ROPE_VAL,
        prior_analysis = "beta",
        alpha = 1,
        beta = 1
      ),
      list(ALT = alt, ROPE_VAL = rope_by_alt[[alt]])
    )
  )

  run_test(
    paste("BFpower.bin | sample-size | interval-null | Moment |", alt),
    substitute(
      BFpower.bin(
        alternative = ALT,
        threshold = 3,
        true_rate = 0.8,
        false_rate = 0.05,
        h0 = 0.5,
        ROPE = ROPE_VAL,
        prior_analysis = "Moment",
        scale = 0.707
      ),
      list(ALT = alt, ROPE_VAL = rope_by_alt[[alt]])
    )
  )
}


############################################################
## B3. Fixed-N mode, point-null
############################################################

for (alt in c("two.sided", "greater", "less")) {
  run_test(
    paste("BFpower.bin | fixed N | point-null | beta |", alt),
    substitute(
      BFpower.bin(
        alternative = ALT,
        threshold = 3,
        h0 = 0.5,
        N = 52,
        prior_analysis = "beta",
        alpha = 1,
        beta = 1
      ),
      list(ALT = alt)
    )
  )

  run_test(
    paste("BFpower.bin | fixed N | point-null | Moment |", alt),
    substitute(
      BFpower.bin(
        alternative = ALT,
        threshold = 3,
        h0 = 0.5,
        N = 52,
        prior_analysis = "Moment",
        scale = 0.707
      ),
      list(ALT = alt)
    )
  )
}


############################################################
## B4. Fixed-N mode, interval-null
############################################################

for (alt in c("two.sided", "greater", "less")) {
  run_test(
    paste("BFpower.bin | fixed N | interval-null | beta |", alt),
    substitute(
      BFpower.bin(
        alternative = ALT,
        threshold = 3,
        h0 = 0.5,
        N = 52,
        ROPE = ROPE_VAL,
        prior_analysis = "beta",
        alpha = 1,
        beta = 1
      ),
      list(ALT = alt, ROPE_VAL = rope_by_alt[[alt]])
    )
  )

  run_test(
    paste("BFpower.bin | fixed N | interval-null | Moment |", alt),
    substitute(
      BFpower.bin(
        alternative = ALT,
        threshold = 3,
        h0 = 0.5,
        N = 52,
        ROPE = ROPE_VAL,
        prior_analysis = "Moment",
        scale = 0.707
      ),
      list(ALT = alt, ROPE_VAL = rope_by_alt[[alt]])
    )
  )
}


############################################################
## B5. Design-prior tests
############################################################

run_test(
  "BFpower.bin | fixed N | design prior beta",
  quote(
    BFpower.bin(
      alternative = "greater",
      threshold = 3,
      h0 = 0.5,
      N = 52,
      prior_analysis = "beta",
      alpha = 1,
      beta = 1,
      prior_design = "beta",
      alpha_d = 1,
      beta_d = 1
    )
  )
)

run_test(
  "BFpower.bin | fixed N | design prior Moment",
  quote(
    BFpower.bin(
      alternative = "greater",
      threshold = 3,
      h0 = 0.5,
      N = 52,
      prior_analysis = "beta",
      alpha = 1,
      beta = 1,
      prior_design = "Moment",
      location_d = 0.7,
      scale_d = 0.707
    )
  )
)

run_test(
  "BFpower.bin | fixed N | design prior Point | greater",
  quote(
    BFpower.bin(
      alternative = "greater",
      threshold = 3,
      h0 = 0.5,
      N = 52,
      prior_analysis = "beta",
      alpha = 1,
      beta = 1,
      prior_design = "Point",
      location_d = 0.7
    )
  )
)

run_test(
  "BFpower.bin | fixed N | design prior Point | less",
  quote(
    BFpower.bin(
      alternative = "less",
      threshold = 3,
      h0 = 0.5,
      N = 52,
      prior_analysis = "beta",
      alpha = 1,
      beta = 1,
      prior_design = "Point",
      location_d = 0.3
    )
  )
)

run_test(
  "BFpower.bin | fixed N | design prior Point | two.sided",
  quote(
    BFpower.bin(
      alternative = "two.sided",
      threshold = 3,
      h0 = 0.5,
      N = 52,
      prior_analysis = "beta",
      alpha = 1,
      beta = 1,
      prior_design = "Point",
      location_d = 0.7
    )
  )
)


############################################################
## B6. type_rate tests
############################################################

for (tr in c("positive", "negative")) {
  run_test(
    paste("BFpower.bin | sample-size | type_rate =", tr),
    substitute(
      BFpower.bin(
        alternative = "greater",
        threshold = 3,
        true_rate = 0.8,
        false_rate = 0.05,
        h0 = 0.5,
        prior_analysis = "beta",
        alpha = 1,
        beta = 1,
        prior_design = "Point",
        location_d = 0.7,
        type_rate = TR
      ),
      list(TR = tr)
    )
  )
}


############################################################
## B7. BFpower object structure
############################################################
test_that("BFpower.bin returns the expected object structure", {
  bfpower_obj <- BFpower.bin(
    alternative = "greater",
    threshold = 3,
    h0 = 0.5,
    N = 52,
    prior_analysis = "beta",
    alpha = 1,
    beta = 1,
    prior_design = "Point",
    location_d = 0.7
  )

  expect_s3_class(bfpower_obj, "BFpower")
  expect_equal(bfpower_obj$type, "One-proportion")

  expect_true("alternative" %in% names(bfpower_obj))
  expect_true("h0" %in% names(bfpower_obj))
  expect_true("ROPE" %in% names(bfpower_obj))
  expect_true("analysis_h1" %in% names(bfpower_obj))
  expect_true("design_h1" %in% names(bfpower_obj))
  expect_true("results" %in% names(bfpower_obj))
  expect_true("threshold" %in% names(bfpower_obj))
  expect_true("setting" %in% names(bfpower_obj))
})


############################################################
## B8. BFpower validation tests
############################################################

run_test(
  "ERROR BFpower.bin | invalid h0 too small",
  quote(BFpower.bin(alternative = "greater", threshold = 3, true_rate = 0.8, false_rate = 0.05, h0 = 0.05, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | invalid h0 too large",
  quote(BFpower.bin(alternative = "greater", threshold = 3, true_rate = 0.8, false_rate = 0.05, h0 = 0.95, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | invalid N zero",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 0, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | invalid N non-integer",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 10.5, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | invalid alternative",
  quote(BFpower.bin(alternative = "wrong", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | invalid threshold",
  quote(BFpower.bin(alternative = "greater", threshold = 0.5, true_rate = 0.8, false_rate = 0.05, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "BFpower.bin | threshold = 1 allowed",
  quote(BFpower.bin(alternative = "greater", threshold = 1, true_rate = 0.8, false_rate = 0.05, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1))
)

run_test(
  "ERROR BFpower.bin | invalid true_rate",
  quote(BFpower.bin(alternative = "greater", threshold = 3, true_rate = 0.5, false_rate = 0.05, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | invalid false_rate",
  quote(BFpower.bin(alternative = "greater", threshold = 3, true_rate = 0.8, false_rate = 0.2, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | invalid type_rate",
  quote(BFpower.bin(alternative = "greater", threshold = 3, true_rate = 0.8, false_rate = 0.05, h0 = 0.5, prior_analysis = "beta", alpha = 1, beta = 1, type_rate = "wrong")),
  expect_error = TRUE
)


############################################################
## B9. BFpower design prior validation tests
############################################################

run_test(
  "ERROR BFpower.bin | invalid prior_design",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "wrong")),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | beta design missing alpha_d",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "beta", beta_d = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | beta design invalid alpha_d",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "beta", alpha_d = 0, beta_d = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | beta design missing beta_d",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "beta", alpha_d = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | beta design invalid beta_d",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "beta", alpha_d = 1, beta_d = 0)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Moment design missing location_d",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Moment", scale_d = 0.707)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Moment design invalid location_d <= 0",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Moment", location_d = 0, scale_d = 0.707)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Moment design invalid location_d >= 1",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Moment", location_d = 1, scale_d = 0.707)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Moment design missing scale_d",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Moment", location_d = 0.7)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Moment design invalid scale_d",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Moment", location_d = 0.7, scale_d = 0)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Point design missing location_d",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Point")),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Point greater location_d <= h0",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Point", location_d = 0.5)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Point greater location_d >= 1",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Point", location_d = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Point less location_d >= h0",
  quote(BFpower.bin(alternative = "less", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Point", location_d = 0.5)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Point less location_d <= 0",
  quote(BFpower.bin(alternative = "less", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Point", location_d = 0)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | Point two.sided location_d == h0",
  quote(BFpower.bin(alternative = "two.sided", threshold = 3, h0 = 0.5, N = 52, prior_analysis = "beta", alpha = 1, beta = 1, prior_design = "Point", location_d = 0.5)),
  expect_error = TRUE
)


############################################################
## B10. BFpower ROPE validation tests
############################################################

run_test(
  "ERROR BFpower.bin | two.sided ROPE length 1",
  quote(BFpower.bin(alternative = "two.sided", threshold = 3, h0 = 0.5, N = 52, ROPE = 0.2, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | two.sided ROPE both positive",
  quote(BFpower.bin(alternative = "two.sided", threshold = 3, h0 = 0.5, N = 52, ROPE = c(0.1, 0.2), prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | greater ROPE negative",
  quote(BFpower.bin(alternative = "greater", threshold = 3, h0 = 0.5, N = 52, ROPE = -0.2, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.bin | less ROPE positive",
  quote(BFpower.bin(alternative = "less", threshold = 3, h0 = 0.5, N = 52, ROPE = 0.2, prior_analysis = "beta", alpha = 1, beta = 1)),
  expect_error = TRUE
)


############################################################
## End
############################################################


