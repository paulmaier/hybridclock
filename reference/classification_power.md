# Expected accuracy of individual classification

Simulates animals of each class from the reference allele frequencies of
the diagnostic markers, adds genotyping error, and classifies them with
[`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md).
The result shows how often each class can be called correctly with these
markers, and what it is mistaken for.

## Usage

``` r
classification_power(markers, n = 1000, error = 0.005, alpha = 0.01, seed = 1)
```

## Arguments

- markers:

  Output of
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md)
  (columns `freq1`, `freq2`).

- n:

  Animals simulated per class.

- error:

  Genotyping error rate used to simulate and classify.

- alpha:

  Passed to
  [`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md).

- seed:

  Random seed.

## Value

A matrix of proportions: true classes in rows, called classes in
columns. The diagonal is the expected accuracy.

## Examples

``` r
ex <- hybridclock_example()
mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
round(classification_power(mk[1:15, ], n = 200), 2)
#>               called
#> true           pure1 pure2 F1   F2 BC1_species1 BC1_species2 BC2_species1
#>   pure1         0.96  0.00  0 0.00         0.00         0.00         0.04
#>   pure2         0.00  0.96  0 0.00         0.00         0.00         0.00
#>   F1            0.00  0.00  1 0.00         0.00         0.00         0.00
#>   F2            0.00  0.00  0 0.94         0.04         0.01         0.00
#>   BC1_species1  0.00  0.00  0 0.00         0.88         0.00         0.12
#>   BC1_species2  0.00  0.00  0 0.03         0.00         0.83         0.00
#>   BC2_species1  0.01  0.00  0 0.00         0.18         0.00         0.80
#>   BC2_species2  0.00  0.03  0 0.00         0.00         0.14         0.00
#>               called
#> true           BC2_species2 later
#>   pure1                0.00     0
#>   pure2                0.04     0
#>   F1                   0.00     0
#>   F2                   0.00     0
#>   BC1_species1         0.00     0
#>   BC1_species2         0.14     0
#>   BC2_species1         0.00     0
#>   BC2_species2         0.84     0
```
