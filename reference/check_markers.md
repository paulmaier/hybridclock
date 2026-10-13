# Check that markers are diagnostic in the parents that actually bred

Markers chosen in reference samples can still be shared by the animals
that founded a hybrid population, for example when a release strain or
cross parent differs from the reference sample. First-generation hybrids
expose such markers: at a truly diagnostic marker every F1 is
heterozygous. This function finds animals that look like F1s (hybrid
index near 0.5 and high interspecies heterozygosity) and reports, for
each marker, how often they are heterozygous there.

## Usage

``` r
check_markers(
  geno,
  markers,
  scores,
  S_range = c(0.4, 0.6),
  min_H = 0.7,
  flag = 0.7
)
```

## Arguments

- geno:

  Genotype matrix as in
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md).

- markers:

  Output of
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md).

- scores:

  Output of
  [`hybrid_scores()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores.md)
  for the same animals.

- S_range:

  Range of hybrid index for likely F1s.

- min_H:

  Smallest interspecies heterozygosity for likely F1s.

- flag:

  Markers heterozygous in fewer than this share of likely F1s are
  flagged.

## Value

A data.frame with one row per marker: `name`, `n_f1` (likely F1s
genotyped), `het_f1` (their share of heterozygotes) and `flagged`.
Attribute `n_f1` gives the number of likely F1s; with fewer than about
5, the check says little.

## Examples

``` r
ex <- hybridclock_example()
mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
chk <- check_markers(ex$genotypes, mk, hybrid_scores(ex$genotypes, mk))
attr(chk, "n_f1")
#> [1] 0
```
