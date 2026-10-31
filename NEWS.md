# MultiSynSurr 0.1.1

* `msr()`: the covariance of two synthetic phenotypes that no participant has
  together is not identified by the likelihood. It is now set to the value that
  makes them independent given the target before the finite-sample correction,
  so the corrected covariance, the standard errors and the tests no longer depend
  on where the optimiser left it. Fits in which every pair of outcomes is
  observed together are unchanged.

# MultiSynSurr 0.1.0

* First release.
