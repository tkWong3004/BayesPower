############################################################
## Comprehensive tests for:
## 1) BF10.props()
## 2) BFpower.props()
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
## PART A: BF10.props()
############################################################

run_test(
  "BF10.props | valid balanced data",
  quote(
    BF10.props(
      a0 = 1, b0 = 1,
      a1 = 1, b1 = 1,
      a2 = 1, b2 = 1,
      N1 = 493, N2 = 488,
      x1 = 155, x2 = 150
    )
  )
)

run_test(
  "BF10.props | valid zero successes",
  quote(
    BF10.props(
      a0 = 1, b0 = 1,
      a1 = 1, b1 = 1,
      a2 = 1, b2 = 1,
      N1 = 20, N2 = 20,
      x1 = 0, x2 = 0
    )
  )
)

run_test(
  "BF10.props | valid all successes",
  quote(
    BF10.props(
      a0 = 1, b0 = 1,
      a1 = 1, b1 = 1,
      a2 = 1, b2 = 1,
      N1 = 20, N2 = 20,
      x1 = 20, x2 = 20
    )
  )
)

run_test(
  "BF10.props | valid extreme difference",
  quote(
    BF10.props(
      a0 = 1, b0 = 1,
      a1 = 1, b1 = 1,
      a2 = 1, b2 = 1,
      N1 = 20, N2 = 20,
      x1 = 0, x2 = 20
    )
  )
)


############################################################
## A2. BF10.props object structure
############################################################
test_that("BF10.props returns the expected object structure", {
  bf10_props_obj <- BF10.props(
    a0 = 1, b0 = 1,
    a1 = 1, b1 = 1,
    a2 = 1, b2 = 1,
    N1 = 493, N2 = 488,
    x1 = 155, x2 = 150
  )

  expect_s3_class(bf10_props_obj, "BFvalue")
  expect_equal(bf10_props_obj$type, "Two-proportions")

  expect_true("bf10" %in% names(bf10_props_obj))
  expect_true("N1" %in% names(bf10_props_obj))
  expect_true("N2" %in% names(bf10_props_obj))
  expect_true("x1" %in% names(bf10_props_obj))
  expect_true("x2" %in% names(bf10_props_obj))
  expect_true("analysis_h0" %in% names(bf10_props_obj))
  expect_true("analysis_h1_theta_1" %in% names(bf10_props_obj))
  expect_true("analysis_h1_theta_2" %in% names(bf10_props_obj))
  expect_true("OddsRatio" %in% names(bf10_props_obj))
  expect_true("p.value" %in% names(bf10_props_obj))
})


############################################################
## A3. BF10.props validation tests
############################################################

run_test(
  "ERROR BF10.props | invalid a0 = 0",
  quote(BF10.props(a0 = 0, b0 = 1, a1 = 1, b1 = 1, a2 = 1, b2 = 1, N1 = 10, N2 = 10, x1 = 5, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | invalid b0 = 0",
  quote(BF10.props(a0 = 1, b0 = 0, a1 = 1, b1 = 1, a2 = 1, b2 = 1, N1 = 10, N2 = 10, x1 = 5, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | invalid a1 = -1",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = -1, b1 = 1, a2 = 1, b2 = 1, N1 = 10, N2 = 10, x1 = 5, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | invalid b1 = Inf",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = 1, b1 = Inf, a2 = 1, b2 = 1, N1 = 10, N2 = 10, x1 = 5, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | invalid a2 = NA",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = 1, b1 = 1, a2 = NA, b2 = 1, N1 = 10, N2 = 10, x1 = 5, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | invalid b2 = 0",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = 1, b1 = 1, a2 = 1, b2 = 0, N1 = 10, N2 = 10, x1 = 5, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | invalid N1 = 0",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = 1, b1 = 1, a2 = 1, b2 = 1, N1 = 0, N2 = 10, x1 = 0, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | invalid N2 non-integer",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = 1, b1 = 1, a2 = 1, b2 = 1, N1 = 10, N2 = 10.5, x1 = 5, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | invalid x1 negative",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = 1, b1 = 1, a2 = 1, b2 = 1, N1 = 10, N2 = 10, x1 = -1, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | invalid x2 non-integer",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = 1, b1 = 1, a2 = 1, b2 = 1, N1 = 10, N2 = 10, x1 = 5, x2 = 5.5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | x1 > N1",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = 1, b1 = 1, a2 = 1, b2 = 1, N1 = 10, N2 = 10, x1 = 11, x2 = 5)),
  expect_error = TRUE
)

run_test(
  "ERROR BF10.props | x2 > N2",
  quote(BF10.props(a0 = 1, b0 = 1, a1 = 1, b1 = 1, a2 = 1, b2 = 1, N1 = 10, N2 = 10, x1 = 5, x2 = 11)),
  expect_error = TRUE
)


############################################################
## PART B: BFpower.props()
############################################################

############################################################
## B1. Sample-size mode
############################################################

run_test(
  "BFpower.props | sample-size | same design priors",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 156, b1 = 339,
      a2 = 151, b2 = 339,
      prior_design_1 = "same",
      prior_design_2 = "same"
    )
  )
)

run_test(
  "BFpower.props | sample-size | beta design priors",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 156, b1 = 339,
      a2 = 151, b2 = 339,
      prior_design_1 = "beta",
      a1d = 2, b1d = 3,
      prior_design_2 = "beta",
      a2d = 3, b2d = 4
    )
  )
)

run_test(
  "BFpower.props | sample-size | Point design priors",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 156, b1 = 339,
      a2 = 151, b2 = 339,
      prior_design_1 = "Point",
      dp1 = 0.35,
      prior_design_2 = "Point",
      dp2 = 0.25
    )
  )
)

run_test(
  "BFpower.props | sample-size | mixed design priors",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 156, b1 = 339,
      a2 = 151, b2 = 339,
      prior_design_1 = "same",
      prior_design_2 = "Point",
      dp2 = 0.25
    )
  )
)


############################################################
## B2. Fixed-sample mode
############################################################

run_test(
  "BFpower.props | fixed sample | same design priors",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 156, b1 = 339,
      a2 = 151, b2 = 339,
      N1 = 493,
      N2 = 488,
      prior_design_1 = "same",
      prior_design_2 = "same"
    )
  )
)

run_test(
  "BFpower.props | fixed sample | beta design priors",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 156, b1 = 339,
      a2 = 151, b2 = 339,
      N1 = 493,
      N2 = 488,
      prior_design_1 = "beta",
      a1d = 2, b1d = 3,
      prior_design_2 = "beta",
      a2d = 3, b2d = 4
    )
  )
)

run_test(
  "BFpower.props | fixed sample | Point design priors",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 156, b1 = 339,
      a2 = 151, b2 = 339,
      N1 = 493,
      N2 = 488,
      prior_design_1 = "Point",
      dp1 = 0.35,
      prior_design_2 = "Point",
      dp2 = 0.25
    )
  )
)


############################################################
## B3. type_rate tests
############################################################

for (tr in c("positive", "negative")) {
  run_test(
    paste("BFpower.props | sample-size | type_rate =", tr),
    substitute(
      BFpower.props(
        threshold = 3,
        true_rate = 0.8,
        a0 = 1, b0 = 1,
        a1 = 156, b1 = 339,
        a2 = 151, b2 = 339,
        prior_design_1 = "Point",
        dp1 = 0.35,
        prior_design_2 = "Point",
        dp2 = 0.25,
        type_rate = TR
      ),
      list(TR = tr)
    )
  )
}


############################################################
## B4. BFpower.props object structure
############################################################
test_that("BFpower.props returns the expected object structure", {
  bfpower_props_obj <- BFpower.props(
    threshold = 3,
    true_rate = 0.8,
    a0 = 1, b0 = 1,
    a1 = 156, b1 = 339,
    a2 = 151, b2 = 339,
    N1 = 493,
    N2 = 488,
    prior_design_1 = "Point",
    dp1 = 0.2,
    prior_design_2 = "Point",
    dp2 = 0.1
  )

  expect_s3_class(bfpower_props_obj, "BFpower")
  expect_equal(bfpower_props_obj$type, "Two-proportions")

  expect_true("analysis_h0" %in% names(bfpower_props_obj))
  expect_true("analysis_h1_theta_1" %in% names(bfpower_props_obj))
  expect_true("analysis_h1_theta_2" %in% names(bfpower_props_obj))
  expect_true("design_h1_theta_1" %in% names(bfpower_props_obj))
  expect_true("design_h1_theta_2" %in% names(bfpower_props_obj))
  expect_true("results" %in% names(bfpower_props_obj))
  expect_true("grid" %in% names(bfpower_props_obj))
  expect_true("threshold" %in% names(bfpower_props_obj))
  expect_true("mode_bf" %in% names(bfpower_props_obj))
})


############################################################
## B5. BFpower.props validation tests
############################################################

run_test(
  "ERROR BFpower.props | only N1 supplied",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      N1 = 100
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | only N2 supplied",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      N2 = 100
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid N1 = 0",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      N1 = 0,
      N2 = 100
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid N2 non-integer",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      N1 = 100,
      N2 = 100.5
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid threshold",
  quote(
    BFpower.props(
      threshold = 0.5,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2
    )
  ),
  expect_error = TRUE
)

run_test(
  "BFpower.props | threshold = 1 allowed",
  quote(
    BFpower.props(
      threshold = 1,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2
    )
  )
)

run_test(
  "ERROR BFpower.props | invalid true_rate too small",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.5,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid type_rate",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      type_rate = "wrong"
    )
  ),
  expect_error = TRUE
)


############################################################
## B6. BFpower.props prior validation tests
############################################################

run_test(
  "ERROR BFpower.props | invalid a0",
  quote(BFpower.props(threshold = 3, true_rate = 0.8, a0 = 0, b0 = 1, a1 = 2, b1 = 2, a2 = 2, b2 = 2)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid b0",
  quote(BFpower.props(threshold = 3, true_rate = 0.8, a0 = 1, b0 = 0, a1 = 2, b1 = 2, a2 = 2, b2 = 2)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid a1",
  quote(BFpower.props(threshold = 3, true_rate = 0.8, a0 = 1, b0 = 1, a1 = 0, b1 = 2, a2 = 2, b2 = 2)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid b1",
  quote(BFpower.props(threshold = 3, true_rate = 0.8, a0 = 1, b0 = 1, a1 = 2, b1 = 0, a2 = 2, b2 = 2)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid a2",
  quote(BFpower.props(threshold = 3, true_rate = 0.8, a0 = 1, b0 = 1, a1 = 2, b1 = 2, a2 = 0, b2 = 2)),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid b2",
  quote(BFpower.props(threshold = 3, true_rate = 0.8, a0 = 1, b0 = 1, a1 = 2, b1 = 2, a2 = 2, b2 = 0)),
  expect_error = TRUE
)


############################################################
## B7. BFpower.props design-prior validation tests
############################################################

run_test(
  "ERROR BFpower.props | invalid prior_design_1",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_1 = "wrong"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | invalid prior_design_2",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_2 = "wrong"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | beta design 1 missing a1d",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_1 = "beta",
      b1d = 2
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | beta design 1 invalid a1d",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_1 = "beta",
      a1d = 0,
      b1d = 2
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | beta design 1 invalid b1d",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_1 = "beta",
      a1d = 2,
      b1d = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | beta design 2 invalid a2d",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_2 = "beta",
      a2d = 0,
      b2d = 2
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | beta design 2 invalid b2d",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_2 = "beta",
      a2d = 2,
      b2d = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | Point design 1 missing dp1",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_1 = "Point"
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | Point design 1 dp1 <= 0",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_1 = "Point",
      dp1 = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | Point design 1 dp1 >= 1",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_1 = "Point",
      dp1 = 1
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | Point design 2 dp2 <= 0",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_2 = "Point",
      dp2 = 0
    )
  ),
  expect_error = TRUE
)

run_test(
  "ERROR BFpower.props | Point design 2 dp2 >= 1",
  quote(
    BFpower.props(
      threshold = 3,
      true_rate = 0.8,
      a0 = 1, b0 = 1,
      a1 = 2, b1 = 2,
      a2 = 2, b2 = 2,
      prior_design_2 = "Point",
      dp2 = 1
    )
  ),
  expect_error = TRUE
)


############################################################
## End
############################################################


