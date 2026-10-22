test_that("msr returns a converged fit with the expected structure", {
  d <- msr_simulate(n = 3000, seed = 11)
  fit <- msr(d$target, d$synthetic, d$X)
  expect_s3_class(fit, "msr")
  expect_true(fit$converged)
  expect_equal(dim(fit$coefficients), c(3L, 3L))
  expect_equal(colnames(fit$coefficients), c("target", "S1", "S2"))
  expect_equal(fit$n_obs, sum(!is.na(d$target)))
  expect_equal(fit$df, fit$n_obs - 3 - 2)
  tst <- msr_test(fit, "exposure")
  expect_named(tst, c("estimate", "se", "t", "df", "p"))
  expect_equal(unname(tst["estimate"]), unname(coef(fit)["exposure"]))
  expect_equal(unname(sqrt(diag(vcov(fit)))["exposure"]), unname(tst["se"]))
})

test_that("with one complete synthetic phenotype msr reproduces SurrogateRegression::FitBNR", {
  skip_if_not_installed("SurrogateRegression")
  d <- msr_simulate(n = 2000, p_assayed = 1, seed = 3)
  X <- cbind(1, d$X)
  ref <- SurrogateRegression::FitBNR(t = d$target, s = d$synthetic[, "S1"], X = X)@Regression.tab
  ref <- ref[ref$Outcome == "Target", ]
  fit <- msr(d$target, d$synthetic[, "S1"], X, intercept = FALSE)
  expect_equal(unname(coef(fit)), ref$Point, tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(vcov(fit)))), ref$SE, tolerance = 1e-4)
})

test_that("a synthetic phenotype available only where the target is gives the OLS estimate", {
  d <- msr_simulate(n = 3000, seed = 5)
  s <- d$synthetic[, "S1"]; s[is.na(d$target)] <- NA
  fit <- msr(d$target, s, d$X)
  o <- !is.na(d$target)
  ols <- stats::lm.fit(cbind(1, d$X[o, ]), d$target[o])$coefficients
  expect_equal(unname(coef(fit)), unname(ols), tolerance = 1e-6)
  expect_equal(msr_expected_ratio(fit), 1, tolerance = 1e-8)
})

test_that("the estimate is close to the truth and more precise than complete case", {
  d <- msr_simulate(n = 20000, beta = 0.1, seed = 7)
  fit <- msr(d$target, d$synthetic, d$X)
  tst <- msr_test(fit, "exposure")
  expect_lt(abs(tst[["estimate"]] - d$beta), 4 * tst[["se"]])
  o <- !is.na(d$target)
  cc <- summary(stats::lm(d$target[o] ~ d$X[o, ]))$coefficients[3, 2]
  expect_lt(tst[["se"]], cc)
  expect_lt(msr_expected_ratio(fit), 1)
  expect_equal(sqrt(msr_expected_ratio(fit)), tst[["se"]] / cc, tolerance = 0.05)
})

test_that("bad input is rejected", {
  d <- msr_simulate(n = 500, seed = 1)
  expect_error(msr(d$target, d$synthetic[-1, ], d$X), "one row per")
  X <- d$X; X[1, 1] <- NA
  expect_error(msr(d$target, d$synthetic, X), "missing values")
  expect_error(msr(d$target, d$synthetic, cbind(d$X, d$X[, 1])), "full rank")
  expect_error(msr_test(msr(d$target, d$synthetic, d$X), "nope"), "unknown term")
})
