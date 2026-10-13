# Hybrid index and interspecies heterozygosity

For each animal, the hybrid index S is the share of species-2 alleles at
the diagnostic markers, and the interspecies heterozygosity H is the
share of markers carrying one allele from each species. With markers
that are fixed or nearly fixed, these simple counts match the
maximum-likelihood estimates of Fitzpatrick (2012) closely.

## Usage

``` r
hybrid_scores(geno, markers, method = c("count", "likelihood"))
```

## Arguments

- geno:

  Genotype matrix as in
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md).

- markers:

  Output of
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md).

- method:

  `"count"` (simple counts; for fixed or nearly fixed markers) or
  `"likelihood"` (maximum-likelihood estimates that use each marker's
  allele frequencies in the two reference groups, for markers that are
  not fixed; Fitzpatrick 2012).

## Value

A data.frame with `S`, `H` and `n_markers` (markers genotyped).

## References

Fitzpatrick BM (2012) Estimating ancestry and heterozygosity of hybrids
using molecular markers. BMC Evolutionary Biology 12:131.

## Examples

``` r
ex <- hybridclock_example()
mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
head(hybrid_scores(ex$genotypes, mk))
#>            S H n_markers
#> species1_1 0 0       314
#> species1_2 0 0       314
#> species1_3 0 0       314
#> species1_4 0 0       314
#> species1_5 0 0       314
#> species1_6 0 0       314
```
