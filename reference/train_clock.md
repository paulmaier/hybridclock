# Train the hybrid clock on a reference set

Fits one Cubist model (Quinlan 1993; Kuhn and Quinlan 2023) per target,
from the summaries of each simulated sample to its true history, and
gives each estimate an interval by split-conformal calibration (Vovk et
al. 2005; Lei et al. 2018). The simulated runs are divided into three
sets so that no run informs more than one step: one to fit the models,
one to calibrate the intervals, and one to test the finished estimator.

## Usage

``` r
train_clock(
  reference,
  targets = c("hybrid_age_both", "hybrid_age_rare", "backcrosses", "arriving",
    "arrival_rate"),
  level = 0.9,
  split = c(fit = 0.6, calibrate = 0.2, test = 0.2),
  max_rows = 40000,
  committees = 20,
  bins = 8,
  intervals = c("normalized", "binned", "cqr"),
  cores = 1,
  seed = 1
)
```

## Arguments

- reference:

  Output of
  [`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md).

- targets:

  Which parts of the history to estimate. See
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)
  for definitions. `arriving` and `arrival_rate` describe new
  hybridization: whether pure animals of the rarer species are still
  arriving, and if so what share of breeders they make up per
  generation.

- level:

  Interval coverage.

- split:

  Shares of runs used to fit, calibrate and test.

- max_rows:

  Most samples used to fit each model.

- committees:

  Number of Cubist committees.

- bins:

  Number of ranges of the estimate within which intervals are
  calibrated.

- intervals:

  How intervals are made: `"binned"` (split-conformal residuals within
  ranges of the estimate and by sample type) or `"cqr"` (conformalized
  quantile regression: a quantile forest gives each sample its own
  interval from its summaries, widened on calibration runs to the stated
  coverage; Romano et al. 2019), or `"normalized"` (intervals centred on
  the estimate, scaled by a forest's prediction of each sample's typical
  error, and calibrated on held-out runs; Papadopoulos et al. 2011).

- cores:

  Number of processes (targets are fitted in parallel).

- seed:

  Random seed.

## Value

An object of class `hybridclock_model`: the fitted models, the interval
calibration, the features used, the reference settings and a
`validation` table of accuracy on the test runs.

## Details

New hybridization is estimated in two parts. `arriving` is the
probability that pure animals of the rarer species are still arriving
(over the three generations before sampling), from a probability forest
(Wright and Ziegler 2017) trained on equal numbers of histories with and
without arrivals, so the prior odds are even. Its probabilities are then
recalibrated on the calibration runs by logistic regression on the
forest's log-odds, the estimated hybrid age and their interaction (Platt
1999), because the same summaries mean different things in young and old
swarms: in a young swarm, pure founders still form new hybrids without
any arrivals. `arrival_rate` is the rate if they are, from a Cubist
model trained only on histories with arrivals.

Ages and rates are modeled on a log scale (`log(age)`,
`log(1 + backcrosses)`, `log(arrival_rate)`). Fitting and calibration
sets are balanced so that every range of the target is equally common;
otherwise the most common ages in the reference set pull estimates
toward them. Intervals are calibrated separately within ranges of the
estimate (Mondrian conformal prediction), so that young and old
estimates each get intervals of the right width. Coverage is guaranteed
on average over histories like those in the reference set, which is why
its settings should match the study system.

## References

Kuhn M, Quinlan R (2023) Cubist: Rule- and instance-based regression
modeling. R package.

Lei J, G'Sell M, Rinaldo A, Tibshirani RJ, Wasserman L (2018)
Distribution-free predictive inference for regression. Journal of the
American Statistical Association 113:1094-1111.

Quinlan JR (1993) Combining instance-based and model-based learning.
Proceedings of the Tenth International Conference on Machine Learning,
236-243.

Vovk V, Gammerman A, Shafer G (2005) Algorithmic Learning in a Random
World. Springer.

Platt J (1999) Probabilistic outputs for support vector machines and
comparisons to regularized likelihood methods. In: Advances in Large
Margin Classifiers, MIT Press, 61-74.

Wright MN, Ziegler A (2017) ranger: A fast implementation of random
forests for high dimensional data in C++ and R. Journal of Statistical
Software 77:1-17.

## Examples

``` r
# \donttest{
ref <- build_reference(n_sims = 150, generations = 12, n_markers = 200, breeders = c(20, 100))
clock <- train_clock(ref, targets = "hybrid_age_both", committees = 5)
clock
#> hybridclock model: hybrid_age_both 
#>   reference: 200 markers, 0 shared loci, samples of 20 animals, breeders 20 to 100 
#>   90% intervals; accuracy on held-out runs:
#>           target typical_error coverage
#>  hybrid_age_both          0.12     0.95
#>   individual classes are used when at least 90% of hybrids fit an early class
# }
```
