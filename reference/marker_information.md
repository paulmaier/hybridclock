# Information carried by a set of markers

The sum over markers of the squared allele-frequency difference between
the two reference groups. A fixed marker contributes 1, so the total is
roughly the number of fixed markers the set is worth. In simulations,
accuracy of hybrid ages fell gradually with it: typical errors were
about 20% at 300, 26% at 126 and 36-39% at 7-12, with honest intervals
throughout.

## Usage

``` r
marker_information(markers)
```

## Arguments

- markers:

  Output of
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md)
  (columns `freq1`, `freq2`).

## Value

A number.

## Examples

``` r
ex <- hybridclock_example()
marker_information(diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2"))
#> [1] 314
```
