# Count species-2 alleles at diagnostic markers

Count species-2 alleles at diagnostic markers

## Usage

``` r
orient_genotypes(geno, markers)
```

## Arguments

- geno:

  Genotype matrix as in
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md).

- markers:

  Output of
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md).

## Value

An integer matrix (animals x markers) counting species-2 alleles.

## Examples

``` r
ex <- hybridclock_example()
mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
dim(orient_genotypes(ex$genotypes, mk))
#> [1] 100 314
```
