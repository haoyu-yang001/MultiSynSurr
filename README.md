# MultiSynSurr

Association testing for a partially measured outcome that borrows strength from
one or more **synthetic phenotypes**: predictions of the outcome from different
sets of surrogates, each available on its own subset of participants.

MultiSynSurr regresses the outcome and the synthetic phenotypes jointly on a
common design, as separate outcomes of a multivariate normal model with
missingness in every outcome. It extends synthetic surrogate (SynSurr) analysis
([McCaw et al. 2024, *Nature Genetics*](https://doi.org/10.1038/s41588-024-01793-9))
from one synthetic phenotype to several. Coefficients are estimated by
generalised least squares given the residual covariance, which is estimated by
profile maximum likelihood with a finite-sample correction; Wald tests use a t
reference. Inference on the outcome does not depend on how accurate the
synthetic phenotypes are.

A typical use is a cohort whose metabolomics was assembled from several
sub-studies: a core metabolite panel measured on most assayed participants gives
an accurate synthetic phenotype for them, and questionnaire data give a weaker
one for everyone.

## Installation

```r
# install.packages("remotes")
remotes::install_github("haoyu-yang001/MultiSynSurr")
```

## Example

```r
library(MultiSynSurr)
d <- msr_simulate(n = 10000, seed = 1)    # target with NA, two synthetic phenotypes, design
fit <- msr(d$target, d$synthetic, d$X)   # MultiSynSurr
fit                                      # target coefficients, t tests
msr_test(fit, "exposure")                # Wald test of one coefficient
sqrt(msr_expected_ratio(fit))            # predicted SE ratio to complete case
```

* `target`: numeric vector, `NA` where not measured.
* `synthetic`: vector or matrix of synthetic phenotypes, one column each, `NA`
  where not available. Keep synthetic phenotypes from different surrogates as
  separate columns; do not merge them into one.
* `X`: numeric design matrix of covariates and tested variables, without missing
  values (an intercept is added if there is none).

Each synthetic phenotype must be built the same way for every participant and
must not use a variable whose association is tested.

## Citation

Yang H, et al. Synthetic phenotypes from a core metabolite panel improve power
for diet–metabolome association studies with partially measured metabolomes.
(In preparation.)

## License

GPL (>= 3)
