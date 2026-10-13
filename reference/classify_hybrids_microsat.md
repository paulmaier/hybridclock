# Classify recent hybrids from microsatellites

The multi-allelic version of
[`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md):
a parent with species-2 ancestry `a` passes on allele `x` with
probability `a * freq2(x) + (1 - a) * freq1(x)`. Each animal is tested
against the pure, F1, F2, first- and second-backcross classes and
against the best fit under any pair of parental ancestries; if no class
fits (likelihood-ratio test at `alpha`), it is called `"later"`.

## Usage

``` r
classify_hybrids_microsat(alleles, ms, alpha = 0.01, error = 0.01)
```

## Arguments

- alleles:

  Allele matrix as in
  [`microsat_markers()`](https://paulmaier.github.io/hybridclock/reference/microsat_markers.md).

- ms:

  Output of
  [`microsat_markers()`](https://paulmaier.github.io/hybridclock/reference/microsat_markers.md).

- alpha:

  Significance level of the fit test. Many later-generation animals (for
  example backcrosses of backcrosses) have genotypes close to an early
  class and cannot be told apart one at a time. A larger `alpha` calls
  more of them "later" at the cost of rejecting some true early hybrids:
  in a bighorn sheep pedigree with 187 microsatellites, `alpha = 0.01`
  kept 91% of early calls right and called 20% of later-generation sheep
  "later"; `alpha = 0.1` gave 85% and 44%. Population ages from
  [`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md)
  do not depend on this choice.

- error:

  Probability that a called allele is replaced by a random allele of
  that locus (genotyping error, null alleles).

## Value

A data.frame with `class`, `probability` and `fit_p`, one row per
animal.

## Examples

``` r
set.seed(1)
a <- cbind(L1.a1 = c(rep(100, 10), rep(104, 10)), L1.a2 = c(rep(100, 10), rep(104, 10)))
g <- rep(c("sp1", "sp2"), each = 10)
ms <- microsat_markers(a, g, "sp1", "sp2")
classify_hybrids_microsat(rbind(a, c(100, 104)), ms)
#>    class probability fit_p
#> 1  pure1   0.3860411     1
#> 2  pure1   0.3860411     1
#> 3  pure1   0.3860411     1
#> 4  pure1   0.3860411     1
#> 5  pure1   0.3860411     1
#> 6  pure1   0.3860411     1
#> 7  pure1   0.3860411     1
#> 8  pure1   0.3860411     1
#> 9  pure1   0.3860411     1
#> 10 pure1   0.3860411     1
#> 11 pure2   0.3860411     1
#> 12 pure2   0.3860411     1
#> 13 pure2   0.3860411     1
#> 14 pure2   0.3860411     1
#> 15 pure2   0.3860411     1
#> 16 pure2   0.3860411     1
#> 17 pure2   0.3860411     1
#> 18 pure2   0.3860411     1
#> 19 pure2   0.3860411     1
#> 20 pure2   0.3860411     1
#> 21    F1   0.3035949     1
```
