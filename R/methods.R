#' Methods for MultiSynSurr fits
#'
#' `coef()` and `vcov()` return the coefficients of one outcome (by default the
#' target) and their covariance; `summary()` tabulates the target's coefficients
#' with t tests.
#'
#' @param object,x An `"msr"` fit.
#' @param outcome Name or index of the outcome: 1 or `"target"` for the target,
#'   or a synthetic phenotype.
#' @param ... Unused.
#' @name msr-methods
NULL

#' @rdname msr-methods
#' @export
coef.msr <- function(object, outcome = 1, ...) {
  object$coefficients[, outcome]
}

#' @rdname msr-methods
#' @export
vcov.msr <- function(object, outcome = 1, ...) {
  if (is.character(outcome)) outcome <- match(outcome, colnames(object$coefficients))
  p <- nrow(object$coefficients)
  i <- (outcome - 1) * p + seq_len(p)
  V <- object$vcov[i, i, drop = FALSE]
  dimnames(V) <- list(rownames(object$coefficients), rownames(object$coefficients))
  V
}

#' @rdname msr-methods
#' @export
summary.msr <- function(object, ...) {
  b <- object$coefficients[, 1]
  se <- sqrt(diag(vcov.msr(object)))
  t <- b / se
  tab <- cbind(Estimate = b, `Std. Error` = se, `t value` = t,
               `Pr(>|t|)` = 2 * stats::pt(-abs(t), object$df))
  structure(list(coefficients = tab, df = object$df, n = object$n, n_obs = object$n_obs,
                 cor = object$cor, call = object$call, converged = object$converged),
            class = "summary.msr")
}

#' @rdname msr-methods
#' @export
print.summary.msr <- function(x, ...) {
  cat("MultiSynSurr fit\nCall: ", deparse(x$call), "\n\n", sep = "")
  cat(sprintf("Participants: %d with any outcome, %d with the target; t reference on %d df\n",
              x$n, x$n_obs, x$df))
  r <- x$cor[1, -1]
  cat("Residual correlation of the target with ",
      paste(sprintf("%s %.3f", names(r), r), collapse = ", "), "\n\n", sep = "")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  if (!x$converged) cat("\nWarning: the optimiser did not converge.\n")
  invisible(x)
}

#' @rdname msr-methods
#' @export
print.msr <- function(x, ...) {
  print(summary.msr(x))
  invisible(x)
}
