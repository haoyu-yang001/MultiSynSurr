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

test_that("synthetic phenotypes never observed together do not make the fit depend on the start", {
  d <- msr_simulate(n = 4000, seed = 9)
  half <- seq_along(d$target) <= 2000
  S <- cbind(d$synthetic, S3 = d$synthetic[, "S2"] + stats::rnorm(length(d$target)))
  S[half, "S2"] <- NA; S[!half, "S3"] <- NA       # S2 and S3 on disjoint halves
  st0 <- diag(4)
  st1 <- st0; st1[3, 4] <- st1[4, 3] <- 0.5
  f0 <- msr(d$target, S, d$X, start = st0)
  f1 <- msr(d$target, S, d$X, start = st1)
  expect_true(f0$converged && f1$converged)
  expect_equal(f0$vcov, f1$vcov, tolerance = 1e-5)
  expect_equal(msr_test(f0, "exposure"), msr_test(f1, "exposure"), tolerance = 1e-5)
  ## before the finite-sample correction the pair is independent given the target
  f2 <- msr(d$target, S, d$X, df_correct = FALSE)
  Om <- f2$Omega
  expect_equal(Om["S2", "S3"], Om["target", "S2"] * Om["target", "S3"] / Om["target", "target"])
})

test_that("robust covariance reduces to HC1 when the target and the synthetic phenotype are never observed together", {
  set.seed(3)
  n <- 300; X <- cbind(1, x = rnorm(n), e = rnorm(n))
  y <- drop(X %*% c(0, 0.3, 0)) + rnorm(n) * exp(0.4 * X[, 3])
  s <- y + rnorm(n); y[sample(n, 150)] <- NA
  s0 <- s; s0[!is.na(y)] <- NA                    # target block decouples: GLS = OLS on the labelled
  f0 <- msr(y, s0, X, intercept = FALSE, robust = TRUE)
  o <- !is.na(y); l <- lm.fit(X[o, ], y[o]); Xi <- chol2inv(qr.R(l$qr))
  hc1 <- Xi %*% crossprod(X[o, ] * l$residuals) %*% Xi * sum(o) / (sum(o) - 3)
  expect_equal(unname(f0$vcov_robust), unname(hc1), tolerance = 1e-8)
  f <- msr(y, s, X, intercept = FALSE, robust = TRUE)
  expect_equal(dim(f$vcov_robust), c(3L, 3L))
  r <- msr_test(f, "e", robust = TRUE)
  expect_true(r[["se"]] > 0 && r[["p"]] > 0 && r[["p"]] <= 1)
  expect_error(msr_test(msr(y, s, X, intercept = FALSE), robust = TRUE))
})
