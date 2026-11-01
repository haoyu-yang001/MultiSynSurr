# MultiSynSurr 0.1.2

* `msr(robust = TRUE)` adds `vcov_robust`, a sandwich covariance of the target's
  coefficients (generalised least-squares information as bread, participants' score
  contributions as meat, HC1 scaling), and `msr_test(robust = TRUE)` uses it. It
  stays valid when the residual covariance depends on the design; without
  co-observed synthetic phenotypes it is the HC1 estimator of least squares.
* Documentation: the finite-sample correction leaves the coefficients unchanged for
  monotone missingness patterns and changes them negligibly otherwise.

# MultiSynSurr 0.1.1

* `msr()`: the covariance of two synthetic phenotypes that no participant has
  together is not identified by the likelihood. It is now set to the value that
  makes them independent given the target before the finite-sample correction,
  so the corrected covariance, the standard errors and the tests no longer depend
  on where the optimiser left it. Fits in which every pair of outcomes is
  observed together are unchanged.

# MultiSynSurr 0.1.0

* First release.
