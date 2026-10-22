## MultiSynSurr: SynSurr with several synthetic phenotypes, each available on its
## own subset of participants. The target Y and the synthetic phenotypes
## S_1, ..., S_K are regressed jointly on the same design X,
##   (Y_i, S_i1, ..., S_iK) = X_i' B + e_i,   e_i ~ N(0, Omega),
## with every outcome missing where it was not measured or not built. Given Omega
## the coefficients are the GLS solution; Omega maximises the profile likelihood.
## Everything is computed from per-missingness-pattern sufficient statistics, so
## the cost does not grow with n.

#' Fit MultiSynSurr
#'
#' Regresses a partially measured target and one or more synthetic phenotypes
#' jointly on a common design, as separate outcomes of a multivariate normal
#' regression in which every outcome may be missing (MultiSynSurr). Each
#' participant contributes the outcomes she or he has. With one synthetic
#' phenotype this is the bivariate model of synthetic surrogate (SynSurr)
#' analysis.
#'
#' The coefficients are the generalised least-squares solution given the
#' residual covariance \eqn{\Omega}, which is estimated by profile maximum
#' likelihood (Cholesky parametrisation, BFGS). With `df_correct = TRUE` the
#' maximum-likelihood \eqn{\Omega} is corrected for the degrees of freedom used
#' by the design, factor by factor in the ordering of the outcomes by
#' availability; the coefficient estimates are unchanged. Wald tests of the
#' target's coefficients are referred to a t distribution with
#' \eqn{n_{obs} - p - K_{obs}} degrees of freedom, where \eqn{n_{obs}} is the
#' number of participants with the target, \eqn{p} the number of design columns
#' and \eqn{K_{obs}} the number of synthetic phenotypes available to them.
#'
#' Validity does not depend on how accurate the synthetic phenotypes are, but
#' each must be built in the same way for every participant, without using any
#' variable whose association is tested, and the model must hold with the same
#' parameters where the target is observed and where it is missing (missing at
#' random). Do not merge synthetic phenotypes built from different surrogates
#' into one column: pass them as separate columns of `synthetic`.
#'
#' @param target Numeric vector, the partially measured outcome; `NA` where not
#'   measured.
#' @param synthetic Numeric vector, matrix or data frame of synthetic phenotypes
#'   (one column each), `NA` where not available. Rows correspond to `target`.
#' @param X Numeric design matrix or data frame (covariates and the tested
#'   variables), without missing values.
#' @param intercept Add an intercept column if `X` has no constant column.
#' @param df_correct Apply the finite-sample correction to \eqn{\Omega}.
#' @param start Optional starting value for \eqn{\Omega}; by default the
#'   per-outcome least-squares residual variances.
#' @param reltol Relative convergence tolerance of the optimiser.
#' @return An object of class `"msr"`: a list with the coefficient matrix
#'   `coefficients` (one column per outcome, the target first), their covariance
#'   `vcov`, the residual covariance `Omega` and correlation `cor`, the degrees
#'   of freedom `df` of the target's Wald tests, `n` (participants with any
#'   outcome), `n_obs` (participants with the target), the missingness
#'   `patterns`, `converged` and `loglik`.
#' @seealso [msr_test()], [msr_expected_ratio()], [msr_simulate()].
#' @references McCaw ZR, Gao J, Lin X, Gronsbell J (2024). Synthetic surrogates
#'   improve power for genome-wide association studies of partially missing
#'   phenotypes in population biobanks. Nature Genetics 56:1527-1536.
#' @examples
#' d <- msr_simulate(n = 3000, seed = 1)
#' fit <- msr(d$target, d$synthetic, d$X)
#' fit
#' msr_test(fit, "exposure")
#' @export
msr <- function(target, synthetic, X, intercept = TRUE, df_correct = TRUE,
                start = NULL, reltol = 1e-10) {
  cl <- match.call()
  target <- as.numeric(target)
  n <- length(target)
  synthetic <- as.matrix(synthetic)
  if (!is.numeric(synthetic)) stop("'synthetic' must be numeric")
  if (nrow(synthetic) != n) stop("'synthetic' must have one row per element of 'target'")
  if (is.null(colnames(synthetic))) colnames(synthetic) <- paste0("S", seq_len(ncol(synthetic)))
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("'X' must be numeric; build factors into columns, e.g. with model.matrix()")
  if (nrow(X) != n) stop("'X' must have one row per element of 'target'")
  if (anyNA(X)) stop("'X' has missing values")
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  const <- apply(X, 2, function(v) all(v == v[1]))
  if (intercept && !any(const)) X <- cbind(`(Intercept)` = 1, X)
  W <- cbind(target = target, synthetic)
  if (sum(!is.na(target)) <= ncol(X)) stop("too few participants with the target for the design")
  qr_obs <- qr(X[!is.na(target), , drop = FALSE])
  if (qr_obs$rank < ncol(X)) stop("'X' is not of full rank among the participants with the target")

  st <- msr_stats(W, X)
  q <- ncol(W)
  if (is.null(start)) {
    v <- vapply(seq_len(q), function(k) {
      o <- !is.na(W[, k])
      if (sum(o) <= ncol(X)) stop("too few participants with outcome '", colnames(W)[k], "'")
      stats::var(stats::lm.fit(X[o, , drop = FALSE], W[o, k])$residuals)
    }, numeric(1))
    start <- diag(v, q)
  }
  prof <- function(th) {
    Om <- msr_om(th, q)
    B <- tryCatch(msr_gls(st, Om)$B, error = function(e) NULL)
    if (is.null(B)) return(1e100)
    -msr_loglik(st, B, Om)
  }
  op <- stats::optim(msr_th(start), prof, method = "BFGS",
                     control = list(maxit = 1000, reltol = reltol))
  Om <- msr_om(op$par, q)
  if (df_correct) Om <- msr_df_correct(Om, st, ncol(X))
  gl <- msr_gls(st, Om)
  V <- chol2inv(gl$chol)
  dimnames(gl$B) <- list(colnames(X), colnames(W))
  lab <- paste(rep(colnames(W), each = ncol(X)), colnames(X), sep = ":")
  dimnames(V) <- list(lab, lab)
  dimnames(Om) <- list(colnames(W), colnames(W))
  n_obs <- sum(vapply(st, function(g) if (g$o[1]) g$n else 0, numeric(1)))
  k_obs <- sum(Reduce(`|`, lapply(st, function(g) if (g$o[1]) g$o[-1] else rep(FALSE, q - 1))))
  structure(list(coefficients = gl$B, vcov = V, Omega = Om, cor = stats::cov2cor(Om),
                 df = n_obs - ncol(X) - k_obs, n = sum(vapply(st, `[[`, numeric(1), "n")),
                 n_obs = n_obs, patterns = lapply(st, function(g) list(o = g$o, n = g$n)),
                 converged = op$convergence == 0, loglik = -op$value, call = cl),
            class = "msr")
}

#' Wald test of a target coefficient
#'
#' @param fit An `"msr"` fit.
#' @param term Name or index of a design column; by default the last one.
#' @return A named numeric vector: estimate, standard error, t statistic,
#'   degrees of freedom and two-sided p-value (t reference).
#' @examples
#' d <- msr_simulate(n = 3000, seed = 1)
#' msr_test(msr(d$target, d$synthetic, d$X), "exposure")
#' @export
msr_test <- function(fit, term = nrow(fit$coefficients)) {
  if (is.character(term)) term <- match(term, rownames(fit$coefficients))
  if (is.na(term) || term < 1 || term > nrow(fit$coefficients)) stop("unknown term")
  b <- fit$coefficients[term, 1]; se <- sqrt(fit$vcov[term, term])
  t <- b / se
  c(estimate = unname(b), se = se, t = unname(t), df = fit$df, p = 2 * stats::pt(-abs(t), fit$df))
}

#' Predicted variance ratio to complete-case analysis
#'
#' The variance of the target's coefficient relative to complete-case least
#' squares, at the fitted \eqn{\Omega} and missingness patterns, assuming every
#' pattern has the same design distribution:
#' \eqn{n_{obs} [(\sum_m n_m M_m)^{-1}]_{TT} / \Omega_{TT}}, where \eqn{M_m} is
#' the inverse of the block of \eqn{\Omega} observed in pattern \eqn{m}. With one
#' synthetic phenotype available to everyone with the target this is
#' \eqn{1 - \pi\rho^2}. Its square root predicts the ratio of standard errors.
#'
#' @param fit An `"msr"` fit.
#' @return The predicted variance ratio (below 1 is a gain).
#' @examples
#' d <- msr_simulate(n = 3000, seed = 1)
#' sqrt(msr_expected_ratio(msr(d$target, d$synthetic, d$X)))
#' @export
msr_expected_ratio <- function(fit) {
  Om <- fit$Omega; q <- nrow(Om)
  P <- matrix(0, q, q); n_obs <- 0
  for (g in fit$patterns) {
    P <- P + g$n * msr_prec(Om, g$o)
    if (g$o[1]) n_obs <- n_obs + g$n
  }
  n_obs * solve(P)[1, 1] / Om[1, 1]
}

#' Simulate data for MultiSynSurr
#'
#' A cohort in which a fraction of participants is assayed and has a first
#' synthetic phenotype, the target is measured on part of the assayed, and
#' everyone has a second synthetic phenotype; useful for examples and checks.
#'
#' @param n Number of participants.
#' @param p_assayed Fraction with the first synthetic phenotype.
#' @param p_measured Fraction of the assayed with the target measured.
#' @param beta Effect of the exposure on the target.
#' @param alpha Effects of the exposure on the two synthetic phenotypes.
#' @param rho Residual correlations: target with each synthetic phenotype, and
#'   between the two.
#' @param seed Optional random seed.
#' @return A list with `target` (with `NA`), `synthetic` (two columns, `S1` with
#'   `NA`), `X` (intercept, a covariate and `exposure`) and the true `beta`.
#' @examples
#' d <- msr_simulate(n = 1000, seed = 1)
#' colMeans(!is.na(cbind(d$target, d$synthetic)))
#' @export
msr_simulate <- function(n = 5000, p_assayed = 0.3, p_measured = 0.25, beta = 0.1,
                         alpha = c(0.5 * beta, 0), rho = c(t1 = 0.6, t2 = 0.3, s12 = 0.18),
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Om <- diag(3)
  Om[1, 2] <- Om[2, 1] <- rho[[1]]; Om[1, 3] <- Om[3, 1] <- rho[[2]]; Om[2, 3] <- Om[3, 2] <- rho[[3]]
  e <- matrix(stats::rnorm(n * 3), n) %*% chol(Om)
  x <- stats::rnorm(n); E <- 0.3 * x + stats::rnorm(n)
  y <- beta * E + 0.3 * x + e[, 1]
  s1 <- alpha[1] * E + 0.2 * x + e[, 2]
  s2 <- alpha[2] * E + 0.3 * x + e[, 3]
  assayed <- stats::runif(n) < p_assayed
  s1[!assayed] <- NA
  y[!assayed | stats::runif(n) >= p_measured] <- NA
  list(target = y, synthetic = cbind(S1 = s1, S2 = s2), X = cbind(covariate = x, exposure = E),
       beta = beta)
}

## ---- internals ---------------------------------------------------------------------

## Sufficient statistics by missingness pattern. W: n x q outcomes (column 1 the
## target), NA where missing; X: n x p design. Rows with no outcome are dropped.
msr_stats <- function(W, X) {
  R <- !is.na(W)
  keep <- rowSums(R) > 0
  W <- W[keep, , drop = FALSE]; X <- X[keep, , drop = FALSE]; R <- R[keep, , drop = FALSE]
  code <- drop(R %*% 2^(seq_len(ncol(W)) - 1))
  W0 <- W; W0[!R] <- 0
  lapply(split(seq_len(nrow(W)), code), function(i) {
    Xi <- X[i, , drop = FALSE]; Wi <- W0[i, , drop = FALSE]
    list(o = R[i[1], ], n = length(i), XtX = crossprod(Xi), XtW = crossprod(Xi, Wi),
         WtW = crossprod(Wi))
  })
}

## Omega restricted to the observed outcomes, inverted and embedded back (q x q).
msr_prec <- function(Om, o) {
  M <- matrix(0, nrow(Om), nrow(Om))
  M[o, o] <- chol2inv(chol(Om[o, o, drop = FALSE]))
  M
}

## GLS for B given Omega; returns B (p x q) and the information for vec(B).
msr_gls <- function(st, Om) {
  p <- nrow(st[[1]]$XtX); q <- nrow(Om)
  info <- matrix(0, p * q, p * q); rhs <- numeric(p * q)
  for (g in st) {
    M <- msr_prec(Om, g$o)
    info <- info + kronecker(M, g$XtX)
    rhs  <- rhs + as.vector(g$XtW %*% M)
  }
  ch <- chol(info)
  list(B = matrix(backsolve(ch, forwardsolve(t(ch), rhs)), p, q), info = info, chol = ch)
}

## Observed-data log-likelihood (up to a constant) at B, Omega.
msr_loglik <- function(st, B, Om) {
  ll <- 0
  for (g in st) {
    o <- g$o
    Q <- g$WtW - crossprod(B, g$XtW) - crossprod(g$XtW, B) + crossprod(B, g$XtX %*% B)
    ch <- chol(Om[o, o, drop = FALSE])
    ll <- ll - 0.5 * (g$n * 2 * sum(log(diag(ch))) + sum(chol2inv(ch) * Q[o, o, drop = FALSE]))
  }
  ll
}

## Omega from its Cholesky factor's free parameters (log-diagonal), and back.
msr_om <- function(th, q) {
  L <- matrix(0, q, q); L[lower.tri(L, TRUE)] <- th
  diag(L) <- exp(diag(L)); tcrossprod(L)
}
msr_th <- function(Om) {
  L <- t(chol(Om)); diag(L) <- log(diag(L)); L[lower.tri(L, TRUE)]
}

## Finite-sample correction of the MLE of Omega, factor by factor in the ordering
## of the outcomes by availability (synthetic phenotypes from the most to the
## least available, the target last): each outcome's variance given the earlier
## ones is multiplied by n_j / (n_j - p - m_j), with n_j the participants with it,
## p the design columns and m_j the earlier outcomes observed together with it;
## the regressions on the earlier outcomes are held fixed, so B stays at the MLE.
msr_df_correct <- function(Om, st, p) {
  q <- nrow(Om)
  n_k <- Reduce(`+`, lapply(st, function(g) g$n * g$o))
  co  <- Reduce(`|`, lapply(st, function(g) outer(g$o, g$o)))
  S   <- seq_len(q)[-1]
  ord <- c(S[order(-n_k[S])], 1L)
  out <- matrix(0, q, q)
  for (i in seq_along(ord)) {
    j <- ord[i]; prev <- ord[seq_len(i - 1)]
    m_j <- sum(co[j, prev])
    den <- n_k[j] - p - m_j
    if (den <= 0) stop("too few participants with outcome ", j, " for the design")
    f <- n_k[j] / den
    if (length(prev)) {
      b <- solve(Om[prev, prev, drop = FALSE], Om[prev, j])
      cond <- Om[j, j] - sum(Om[prev, j] * b)
      cr <- drop(out[prev, prev, drop = FALSE] %*% b)
      out[prev, j] <- out[j, prev] <- cr
      out[j, j] <- cond * f + sum(cr * b)
    } else {
      out[j, j] <- Om[j, j] * f
    }
  }
  out
}
