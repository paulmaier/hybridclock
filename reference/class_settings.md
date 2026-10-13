# Settings for the classes that choose between the two tiers

The share of hybrids in each early class is part of the summaries the
clock learns from, and it decides whether a population is dated by its
individual classes or by the clock (see
[`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md)).
Those classes must be made the same way in the simulated reference and
in the data, so the settings are given once, to
[`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md),
and are then used for both.

## Usage

``` r
class_settings(
  error = 0.005,
  het_dropout = 0,
  effective_markers = NULL,
  markers = NULL,
  alpha = 0.01
)
```

## Arguments

- error:

  Genotyping error rate assumed by the classes.

- het_dropout:

  Share of heterozygotes called homozygous (see
  [`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md));
  for RAD or GBS data where F1s show clearly less than full
  heterozygosity.

- effective_markers:

  Optional effective number of independent markers (see
  [`effective_markers()`](https://paulmaier.github.io/hybridclock/reference/effective_markers.md)),
  for many linked markers.

- markers:

  Optional output of
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md)
  (columns `freq1`, `freq2`) for markers that are not fixed; with it,
  animals are classified from their genotype counts at the markers.

- alpha:

  Significance level of the fit test for "later".

## Value

A list of class `hybridclock_class_settings`.

## Details

With markers that are not fixed, class probabilities are averaged over
markers using the reference allele frequencies, which needs only a few
summary moments of the frequencies and keeps the classification fast
enough to run on every simulated sample.

## Examples

``` r
class_settings(het_dropout = 0.2)
#> class settings: error 0.005  het_dropout 0.2  effective markers none  
```
