# Predict histories from sample summaries

Predict histories from sample summaries

## Usage

``` r
# S3 method for class 'hybridclock_model'
predict(object, newdata, ...)
```

## Arguments

- object:

  Output of
  [`train_clock()`](https://paulmaier.github.io/hybridclock/reference/train_clock.md).

- newdata:

  A data.frame (or matrix) of summaries, one row per sample, as from
  [`sample_summaries()`](https://paulmaier.github.io/hybridclock/reference/sample_summaries.md)
  or
  [`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md).

- ...:

  Not used.

## Value

A data.frame with an estimate and interval (`_lower`, `_upper`) for each
target.

## Examples

``` r
# \donttest{
ref <- build_reference(n_sims = 150, generations = 12, n_markers = 200, breeders = c(20, 100))
clock <- train_clock(ref, targets = "hybrid_age_both", committees = 5)
head(predict(clock, ref))
#>   hybrid_age_both hybrid_age_both_lower hybrid_age_both_upper
#> 1        1.002518              1.000000              1.011505
#> 2        1.660544              1.425282              1.950226
#> 3        2.202945              1.623903              3.036716
#> 4        3.406930              2.320812              5.103216
#> 5        5.008823              3.200990              8.024186
#> 6        4.179676              2.602793              6.880985
# }
```
