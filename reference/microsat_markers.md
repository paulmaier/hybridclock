# Microsatellite markers

Allele frequencies of multi-allelic markers (such as microsatellites) in
two reference groups, for
[`hybrid_scores_microsat()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores_microsat.md)
and for simulating reference histories (`microsat` in
[`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)
and
[`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md)).
Microsatellites rarely have alleles fixed between species, so each
animal's hybrid index and interspecies heterozygosity are estimated by
likelihood from these frequencies.

## Usage

``` r
microsat_markers(alleles, groups, ref1, ref2, prior = 0.5)
```

## Arguments

- alleles:

  Integer matrix of allele sizes, two columns per locus (named
  `locus.a1`, `locus.a2`), animals in rows; `NA` for missing.

- groups:

  Group label of each animal.

- ref1, ref2:

  Labels of the two reference groups.

- prior:

  Pseudo-count added to every allele seen in either group, so that an
  allele absent from one reference sample is rare, not impossible, in
  that species.

## Value

A list with `loci` (locus names), `alleles` (per locus, the allele
sizes), and `freq1`, `freq2` (per locus, their frequencies in each
reference group).

## Examples

``` r
set.seed(1)
a <- cbind(L1.a1 = c(rep(100, 10), rep(104, 10)), L1.a2 = c(rep(100, 10), rep(104, 10)))
ms <- microsat_markers(a, rep(c("sp1", "sp2"), each = 10), "sp1", "sp2")
ms$freq1
#> $L1
#> [1] 0.97619048 0.02380952
#> 
```
