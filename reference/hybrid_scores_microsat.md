# Hybrid index and heterozygosity from microsatellites

Maximum-likelihood hybrid index (S) and interspecies heterozygosity (H)
for each animal from multi-allelic markers, following Fitzpatrick
(2012): at each locus an animal carries two species-1 copies, one of
each, or two species-2 copies, with probabilities `1 - S - H/2`, `H` and
`S - H/2`, and each copy's allele is drawn from that species'
frequencies. Alleles not seen in either reference are treated as
missing.

## Usage

``` r
hybrid_scores_microsat(alleles, ms)
```

## Arguments

- alleles:

  Allele matrix as in
  [`microsat_markers()`](https://paulmaier.github.io/hybridclock/reference/microsat_markers.md).

- ms:

  Output of
  [`microsat_markers()`](https://paulmaier.github.io/hybridclock/reference/microsat_markers.md).

## Value

A data.frame with `S`, `H` and `n_markers` (loci genotyped).

## Examples

``` r
set.seed(1)
a <- cbind(L1.a1 = c(rep(100, 10), rep(104, 10)), L1.a2 = c(rep(100, 10), rep(104, 10)))
g <- rep(c("sp1", "sp2"), each = 10)
ms <- microsat_markers(a, g, "sp1", "sp2")
hybrid_scores_microsat(rbind(a, c(100, 104)), ms)
#>      S H n_markers
#> 1  0.0 0         1
#> 2  0.0 0         1
#> 3  0.0 0         1
#> 4  0.0 0         1
#> 5  0.0 0         1
#> 6  0.0 0         1
#> 7  0.0 0         1
#> 8  0.0 0         1
#> 9  0.0 0         1
#> 10 0.0 0         1
#> 11 1.0 0         1
#> 12 1.0 0         1
#> 13 1.0 0         1
#> 14 1.0 0         1
#> 15 1.0 0         1
#> 16 1.0 0         1
#> 17 1.0 0         1
#> 18 1.0 0         1
#> 19 1.0 0         1
#> 20 1.0 0         1
#> 21 0.5 1         1
```
